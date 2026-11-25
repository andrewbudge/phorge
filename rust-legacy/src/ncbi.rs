//! Minimal async client for the NCBI E-utilities used by phorge.
//!
//! This is the network boundary, so it exposes a typed [`NcbiError`]
//! (thiserror) for the failures callers may want to distinguish — notably the
//! NCBI habit of returning HTTP 200 with an error payload in the body. Command
//! code layers `anyhow` context on top of these.

use governor::clock::DefaultClock;
use governor::state::{InMemoryState, NotKeyed};
use governor::{Quota, RateLimiter};
use serde_json::Value;
use std::num::NonZeroU32;
use std::time::Duration;
use thiserror::Error;
use tracing::warn;

const EUTILS_BASE: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";

#[derive(Debug, Error)]
pub enum NcbiError {
    #[error("network error contacting NCBI: {0}")]
    Http(#[from] reqwest::Error),
    /// NCBI replied with HTTP 200 but an error message in the JSON body.
    #[error("NCBI returned an error: {0}")]
    Api(String),
    /// The JSON parsed but did not contain the fields we expected.
    #[error("unexpected NCBI response shape: {0}")]
    Shape(String),
}

impl NcbiError {
    /// True for transient failures worth retrying: rate-limit (429), timeouts,
    /// connection drops, and mid-stream body/decode failures (NCBI sometimes
    /// truncates a large response). Deterministic `Api`/`Shape` errors are never
    /// retried — the same request would fail the same way.
    fn is_retryable(&self) -> bool {
        match self {
            NcbiError::Http(e) => {
                e.is_timeout()
                    || e.is_connect()
                    || e.is_body()
                    || e.is_decode()
                    || e.status() == Some(reqwest::StatusCode::TOO_MANY_REQUESTS)
            }
            NcbiError::Api(_) | NcbiError::Shape(_) => false,
        }
    }
}

type DirectLimiter = RateLimiter<NotKeyed, InMemoryState, DefaultClock>;

/// A handle to a result set parked on NCBI's history server, returned by
/// [`EutilsClient::esearch_history`]. esummary pages are pulled against it
/// rather than shipping huge UID lists back and forth.
pub struct SearchHandle {
    pub count: usize,
    pub web_env: String,
    pub query_key: String,
}

/// A candidate taxon returned when resolving a scientific/common name to a
/// TaxID via `esearch db=taxonomy`. `division` (e.g. "Vertebrates", "Bacteria")
/// is carried so the caller can tell cross-kingdom homonyms apart in the confirm
/// prompt.
pub struct TaxonMatch {
    pub taxid: u64,
    pub name: String,
    pub rank: String,
    pub division: String,
}

pub struct EutilsClient {
    http: reqwest::Client,
    limiter: DirectLimiter,
    api_key: Option<String>,
    email: String,
    tool: String,
}

impl EutilsClient {
    /// Build a client. Rate limit follows NCBI policy: 3 req/s without an API
    /// key, 10 req/s with one. `email` is sent on every request per NCBI's ToS.
    pub fn new(api_key: Option<String>, email: String) -> Result<Self, NcbiError> {
        let rps = if api_key.is_some() { 10 } else { 3 };
        let quota = Quota::per_second(NonZeroU32::new(rps).expect("rps is non-zero"));
        let limiter = RateLimiter::direct(quota);
        let http = reqwest::Client::builder()
            .user_agent(format!("phorge/{} ({})", env!("CARGO_PKG_VERSION"), email))
            .build()?;
        Ok(Self {
            http,
            limiter,
            api_key,
            email,
            tool: "phorge".to_string(),
        })
    }

    /// Params attached to every request: tool + email identify us to NCBI,
    /// retmode=json keeps us off the XML path, api_key (if any) unlocks 10 req/s.
    fn common_params(&self) -> Vec<(&'static str, String)> {
        let mut params = vec![
            ("tool", self.tool.clone()),
            ("email", self.email.clone()),
            ("retmode", "json".to_string()),
        ];
        if let Some(key) = &self.api_key {
            params.push(("api_key", key.clone()));
        }
        params
    }

    /// Run one request-and-body-read as a single retryable unit, under the rate
    /// limiter, with exponential backoff. Crucially this wraps the *whole*
    /// operation — send, status check, and body decode — so a response NCBI drops
    /// mid-stream ("end of file before message length reached") is retried, not
    /// just a failed send. `op` is a closure so each attempt issues a fresh
    /// request. Only transient errors (see [`NcbiError::is_retryable`]) retry;
    /// deterministic API/shape errors fail fast.
    async fn retrying<T, Fut>(&self, op: impl Fn() -> Fut) -> Result<T, NcbiError>
    where
        Fut: std::future::Future<Output = Result<T, NcbiError>>,
    {
        const MAX_RETRIES: u32 = 5;
        let mut attempt = 0;
        loop {
            attempt += 1;
            self.limiter.until_ready().await;
            match op().await {
                Ok(value) => return Ok(value),
                Err(e) if attempt <= MAX_RETRIES && e.is_retryable() => {
                    let backoff = Duration::from_millis(400 * 2u64.pow(attempt - 1));
                    warn!(
                        attempt,
                        backoff_ms = backoff.as_millis() as u64,
                        error = %e,
                        "NCBI request failed transiently; backing off"
                    );
                    tokio::time::sleep(backoff).await;
                }
                Err(e) => return Err(e),
            }
        }
    }

    /// Rate-limited GET that returns the parsed JSON body. Send + body decode are
    /// retried as a unit; HTTP-200-with-error-body is handled by each caller,
    /// which knows where the error field lives.
    async fn get(&self, endpoint: &str, extra: Vec<(&str, String)>) -> Result<Value, NcbiError> {
        let mut params = self.common_params();
        params.extend(extra);
        let url = format!("{EUTILS_BASE}/{endpoint}");
        self.retrying(|| async {
            let resp = self.http.get(&url).query(&params).send().await?;
            let resp = resp.error_for_status()?;
            Ok(resp.json::<Value>().await?)
        })
        .await
    }

    /// esearch with `usehistory=y`; parks the full result set on NCBI's history
    /// server and returns its size + handle.
    pub async fn esearch_history(&self, db: &str, term: &str) -> Result<SearchHandle, NcbiError> {
        let body = self
            .get(
                "esearch.fcgi",
                vec![
                    ("db", db.to_string()),
                    ("term", term.to_string()),
                    ("usehistory", "y".to_string()),
                ],
            )
            .await?;

        let result = body
            .get("esearchresult")
            .ok_or_else(|| NcbiError::Shape("missing 'esearchresult'".to_string()))?;

        // esearch reports usage errors as an ERROR field inside a 200 response.
        if let Some(err) = result.get("ERROR").and_then(|v| v.as_str()) {
            return Err(NcbiError::Api(err.to_string()));
        }

        let count = result
            .get("count")
            .and_then(|v| v.as_str())
            .and_then(|s| s.parse::<usize>().ok())
            .ok_or_else(|| NcbiError::Shape("missing or non-numeric 'count'".to_string()))?;
        let web_env = result
            .get("webenv")
            .and_then(|v| v.as_str())
            .ok_or_else(|| NcbiError::Shape("missing 'webenv'".to_string()))?
            .to_string();
        let query_key = result
            .get("querykey")
            .and_then(|v| v.as_str())
            .ok_or_else(|| NcbiError::Shape("missing 'querykey'".to_string()))?
            .to_string();

        Ok(SearchHandle {
            count,
            web_env,
            query_key,
        })
    }

    /// Fetch one page of esummary docsums from a history handle. Returns the raw
    /// JSON; docsum parsing lives in the command so the client stays generic.
    pub async fn esummary_page(
        &self,
        db: &str,
        handle: &SearchHandle,
        retstart: usize,
        retmax: usize,
    ) -> Result<Value, NcbiError> {
        self.get(
            "esummary.fcgi",
            vec![
                ("db", db.to_string()),
                ("query_key", handle.query_key.clone()),
                ("WebEnv", handle.web_env.clone()),
                ("retstart", retstart.to_string()),
                ("retmax", retmax.to_string()),
            ],
        )
        .await
    }

    /// efetch a chunk of records as FASTA text. Unlike the JSON endpoints this
    /// returns plain text, and the id list is POSTed as form data — accession
    /// chunks (up to ~500) are too long to ride on a GET query string. Callers
    /// pass explicit accessions (never a WebEnv handle), so a resumed download
    /// never depends on NCBI's history-server TTL.
    pub async fn efetch_fasta(&self, db: &str, ids: &[&str]) -> Result<String, NcbiError> {
        // efetch fasta is text, not JSON, so we don't reuse common_params (which
        // forces retmode=json).
        let mut params: Vec<(&str, String)> = vec![
            ("tool", self.tool.clone()),
            ("email", self.email.clone()),
            ("db", db.to_string()),
            ("rettype", "fasta".to_string()),
            ("retmode", "text".to_string()),
            ("id", ids.join(",")),
        ];
        if let Some(key) = &self.api_key {
            params.push(("api_key", key.clone()));
        }
        let url = format!("{EUTILS_BASE}/efetch.fcgi");
        self.retrying(|| async {
            let resp = self.http.post(&url).form(&params).send().await?;
            let resp = resp.error_for_status()?;
            let body = resp.text().await?;
            // efetch reports bad requests as a 200 carrying a plain-text/HTML error
            // rather than FASTA; anything without a record marker is an API error.
            if !body.contains('>') {
                let snippet: String = body.trim().chars().take(200).collect();
                return Err(NcbiError::Api(format!(
                    "efetch returned no FASTA records: {snippet}"
                )));
            }
            Ok(body)
        })
        .await
    }

    /// Resolve a TaxID to its scientific name via taxonomy esummary
    /// (e.g. 89829 -> "Leptophlebiidae").
    pub async fn taxonomy_name(&self, taxid: u64) -> Result<String, NcbiError> {
        let body = self
            .get(
                "esummary.fcgi",
                vec![("db", "taxonomy".to_string()), ("id", taxid.to_string())],
            )
            .await?;
        body.pointer(&format!("/result/{taxid}/scientificname"))
            .and_then(|v| v.as_str())
            .map(|s| s.to_string())
            .ok_or_else(|| NcbiError::Shape(format!("no scientificname for taxid {taxid}")))
    }

    /// Resolve a taxon name to its candidate TaxID(s) via `esearch db=taxonomy`.
    /// Usually one match; more than one means a homonym the caller must
    /// disambiguate. An empty vec means no such taxon. A single esummary follows
    /// to attach name/rank/division to each candidate.
    pub async fn resolve_taxon_name(&self, name: &str) -> Result<Vec<TaxonMatch>, NcbiError> {
        let body = self
            .get(
                "esearch.fcgi",
                vec![("db", "taxonomy".to_string()), ("term", name.to_string())],
            )
            .await?;
        let result = body
            .get("esearchresult")
            .ok_or_else(|| NcbiError::Shape("missing 'esearchresult'".to_string()))?;
        if let Some(err) = result.get("ERROR").and_then(|v| v.as_str()) {
            return Err(NcbiError::Api(err.to_string()));
        }
        // For the taxonomy DB the UIDs in idlist are the TaxIDs themselves.
        let ids: Vec<String> = result
            .get("idlist")
            .and_then(|v| v.as_array())
            .map(|a| {
                a.iter()
                    .filter_map(|v| v.as_str().map(String::from))
                    .collect()
            })
            .unwrap_or_default();
        if ids.is_empty() {
            return Ok(Vec::new());
        }

        let summary = self
            .get(
                "esummary.fcgi",
                vec![("db", "taxonomy".to_string()), ("id", ids.join(","))],
            )
            .await?;
        let result = summary
            .get("result")
            .ok_or_else(|| NcbiError::Shape("esummary missing 'result'".to_string()))?;

        let mut matches = Vec::with_capacity(ids.len());
        for id in &ids {
            let Some(doc) = result.get(id) else { continue };
            let field = |k: &str| {
                doc.get(k)
                    .and_then(|v| v.as_str())
                    .unwrap_or_default()
                    .to_string()
            };
            matches.push(TaxonMatch {
                taxid: id.parse().unwrap_or(0),
                name: field("scientificname"),
                rank: doc
                    .get("rank")
                    .and_then(|v| v.as_str())
                    .unwrap_or("no rank")
                    .to_string(),
                division: field("division"),
            });
        }
        Ok(matches)
    }
}
