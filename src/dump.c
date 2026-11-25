#define _POSIX_C_SOURCE 200809L

#include "cli_help.h"
#include "cli_util.h"
#include "dataset.h"
#include "version.h"

#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define ES_URL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
#define EF_URL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

typedef struct {
    char **data;
    size_t len;
    size_t cap;
} StrVec;

static void vec_push(StrVec *v, char *s) {
    if (v->len == v->cap) {
        size_t cap = v->cap ? v->cap * 2 : 64;
        char **p = realloc(v->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        v->data = p;
        v->cap = cap;
    }
    v->data[v->len++] = s;
}

static void vec_free(StrVec *v) {
    for (size_t i = 0; i < v->len; i++) free(v->data[i]);
    free(v->data);
    v->data = NULL;
    v->len = v->cap = 0;
}

static void pace(bool has_key) {
    struct timespec ts = {0, has_key ? 110000000L : 340000000L};
    nanosleep(&ts, NULL);
}

static void curl_base(char **args, size_t *n, const char *url) {
    *n = 0;
    args[(*n)++] = (char *)"curl";
    args[(*n)++] = (char *)"-sS";
    args[(*n)++] = (char *)"-f";
    args[(*n)++] = (char *)"--retry";
    args[(*n)++] = (char *)"3";
    args[(*n)++] = (char *)"--retry-delay";
    args[(*n)++] = (char *)"1";
    args[(*n)++] = (char *)"--connect-timeout";
    args[(*n)++] = (char *)"15";
    args[(*n)++] = (char *)"-G";
    args[(*n)++] = (char *)url;
}

static void add_field(char **args, size_t *n, const char *name, const char *value,
                      StrVec *owned) {
    args[(*n)++] = (char *)"--data-urlencode";
    char *s = xmalloc(strlen(name) + strlen(value) + 2);
    sprintf(s, "%s=%s", name, value);
    args[(*n)++] = s;
    vec_push(owned, s);
}

static void esearch(const char *db, const char *term, size_t retstart, size_t retmax,
                    const char *key, const char *out_path, const char *err_path) {
    char *args[32];
    size_t n = 0;
    StrVec owned = {0};
    char start[32], max[32];
    snprintf(start, sizeof(start), "%zu", retstart);
    snprintf(max, sizeof(max), "%zu", retmax);
    curl_base(args, &n, ES_URL);
    add_field(args, &n, "db", db, &owned);
    add_field(args, &n, "term", term, &owned);
    add_field(args, &n, "retmode", "json", &owned);
    add_field(args, &n, "retstart", start, &owned);
    add_field(args, &n, "retmax", max, &owned);
    if (key) add_field(args, &n, "api_key", key, &owned);
    args[n] = NULL;
    int rc = cli_spawn_wait("curl", args, out_path, err_path, true);
    vec_free(&owned);
    if (rc != 0) xdie("dump: ESearch failed (see %s)", err_path);
}

static void efetch(const char *db, const char *ids, const char *format,
                   const char *key, const char *out_path, const char *err_path) {
    char *args[32];
    size_t n = 0;
    StrVec owned = {0};
    curl_base(args, &n, EF_URL);
    add_field(args, &n, "db", db, &owned);
    add_field(args, &n, "id", ids, &owned);
    add_field(args, &n, "retmode", "text", &owned);
    add_field(args, &n, "rettype", !strcmp(format, "gb") ? "gbwithparts" : "fasta", &owned);
    if (key) add_field(args, &n, "api_key", key, &owned);
    args[n] = NULL;
    int rc = cli_spawn_wait("curl", args, out_path, err_path, true);
    vec_free(&owned);
    if (rc != 0) xdie("dump: EFetch failed (see %s)", err_path);
}

static void parse_ids(const char *raw, StrVec *ids) {
    const char *p = strchr(raw ? raw : "", '[');
    if (!p) return;
    while ((p = strchr(p, '"')) != NULL) {
        p++;
        const char *end = strchr(p, '"');
        if (!end) break;
        if (end > p) vec_push(ids, strndup(p, (size_t)(end - p)));
        p = end + 1;
    }
}

static void read_ids_file(const char *path, StrVec *ids) {
    FILE *f = fopen(path, "rb");
    if (!f) xdie("dump: could not read %s: %s", path, strerror(errno));
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    while ((len = getline(&line, &cap, f)) != -1) {
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = '\0';
        char *s = line;
        while (*s && (*s == ' ' || *s == '\t')) s++;
        if (!*s || *s == '#') continue;
        char *end = s + strlen(s);
        while (end > s && (end[-1] == ' ' || end[-1] == '\t')) *--end = '\0';
        vec_push(ids, xstrdup(s));
    }
    free(line);
    fclose(f);
}

static char *join_ids(const StrVec *ids, size_t start, size_t count) {
    size_t cap = 256;
    char *out = xmalloc(cap);
    size_t len = 0;
    for (size_t i = 0; i < count; i++) {
        const char *id = ids->data[start + i];
        size_t n = strlen(id);
        if (len + n + 2 > cap) {
            while (len + n + 1 > cap) cap *= 2;
            char *p = realloc(out, cap);
            if (!p) xdie("out of memory");
            out = p;
        }
        if (len) out[len++] = ',';
        memcpy(out + len, id, n);
        len += n;
    }
    out[len] = '\0';
    return out;
}

static void json_string(FILE *f, const char *s) {
    fputc('"', f);
    for (const char *p = s; *p; p++) {
        unsigned char c = (unsigned char)*p;
        if (c == '"' || c == '\\') { fputc('\\', f); fputc((int)c, f); }
        else if (c == '\n') fputs("\\n", f);
        else if (c == '\t') fputs("\\t", f);
        else if (c < 0x20) fprintf(f, "\\u%04x", c);
        else fputc((int)c, f);
    }
    fputc('"', f);
}

static char *read_query_file(const char *path) {
    char *text = ds_read_text(path);
    if (!text) xdie("dump: could not read %s: %s", path, strerror(errno));
    for (char *p = text; *p; p++)
        if (*p == '\r' || *p == '\n') *p = ' ';
    char *start = text;
    while (*start && (*start == ' ' || *start == '\t')) start++;
    char *end = start + strlen(start);
    while (end > start && (end[-1] == ' ' || end[-1] == '\t')) *--end = '\0';
    if (!*start) xdie("dump: query file %s is empty", path);
    char *out = xstrdup(start);
    free(text);
    return out;
}

static void emit_ids(const StrVec *ids, const char *path) {
    FILE *f;
    if (!strcmp(path, "-")) {
        f = stdout;
    } else {
        f = fopen(path, "wb");
        if (!f) xdie("dump: could not create %s: %s", path, strerror(errno));
    }
    for (size_t i = 0; i < ids->len; i++) fprintf(f, "%s\n", ids->data[i]);
    if (f != stdout) fclose(f);
}

static void write_provenance(const char *path, const char *db, const char *query,
                             size_t count, size_t fetched, const char *format,
                             size_t batch, bool has_key) {
    FILE *f = fopen(path, "wb");
    if (!f) xdie("dump: could not create %s: %s", path, strerror(errno));
    char stamp[32] = "unknown";
    time_t now = time(NULL);
    struct tm utc;
    if (gmtime_r(&now, &utc)) strftime(stamp, sizeof(stamp), "%Y-%m-%dT%H:%M:%SZ", &utc);
    fputs("{\n  \"tool\": \"phorge " PHORGE_VERSION "\",\n  \"db\": ", f);
    json_string(f, db);
    fputs(",\n  \"query\": ", f);
    json_string(f, query ? query : "");
    fprintf(f, ",\n  \"count\": %zu,\n  \"fetched\": %zu,\n  \"format\": ", count, fetched);
    json_string(f, format);
    fprintf(f, ",\n  \"batch\": %zu,\n  \"date_utc\": ", batch);
    json_string(f, stamp);
    fprintf(f, ",\n  \"api_key_used\": %s\n}\n", has_key ? "true" : "false");
    fclose(f);
}

int cmd_dump(int argc, char **argv) {
    const char *query = NULL, *output = NULL, *ids_file = NULL, *ids_out = NULL;
    const char *db = "nuccore", *format = "fasta";
    char *owned_query = NULL;
    size_t limit = 0, batch = 200;
    bool dry_run = false, force = false;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-q") || !strcmp(a, "--query")) query = next_value(argc, argv, &i);
        else if (!strncmp(a, "--query=", 8)) query = a + 8;
        else if (!strcmp(a, "--query-file")) {
            if (query && query[0] != '@') xdie("dump: choose either --query or --query-file");
            owned_query = read_query_file(next_value(argc, argv, &i));
            query = owned_query;
        }
        else if (!strncmp(a, "--query-file=", 13)) {
            if (query && query[0] != '@') xdie("dump: choose either --query or --query-file");
            owned_query = read_query_file(a + 13);
            query = owned_query;
        }
        else if (!strcmp(a, "-o") || !strcmp(a, "--output")) output = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) output = a + 9;
        else if (!strcmp(a, "--ids")) ids_file = next_value(argc, argv, &i);
        else if (!strcmp(a, "--ids-out")) ids_out = next_value(argc, argv, &i);
        else if (!strncmp(a, "--ids-out=", 10)) ids_out = a + 10;
        else if (!strcmp(a, "--db")) db = next_value(argc, argv, &i);
        else if (!strcmp(a, "--format")) format = next_value(argc, argv, &i);
        else if (!strcmp(a, "-n") || !strcmp(a, "--limit")) limit = (size_t)strtoull(next_value(argc, argv, &i), NULL, 10);
        else if (!strncmp(a, "--limit=", 8)) limit = (size_t)strtoull(a + 8, NULL, 10);
        else if (!strcmp(a, "--batch")) batch = (size_t)strtoull(next_value(argc, argv, &i), NULL, 10);
        else if (!strcmp(a, "--dry-run")) dry_run = true;
        else if (!strcmp(a, "--force")) force = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("dump"); free(owned_query); return 0; }
        else if (a[0] == '-') xdie("dump: unrecognized argument: %s", a);
        else xdie("dump: unexpected argument: %s", a);
    }
    if (query && query[0] == '@' && !owned_query) {
        owned_query = read_query_file(query + 1);
        query = owned_query;
    }
    if (!output && !ids_out) xdie("dump: requires --output or --ids-out");
    if (query && ids_file) xdie("dump: choose either --query or --ids");
    if (!query && !ids_file) xdie("dump: missing --query (or --ids)");
    if (strcmp(format, "fasta") && strcmp(format, "gb"))
        xdie("dump: --format must be fasta or gb");
    if (!batch) xdie("dump: --batch must be positive");

    const char *key = getenv("NCBI_API_KEY");
    bool has_key = key && *key;
    char *curl_ver[] = {(char *)"curl", (char *)"--version", NULL};
    if (cli_spawn_wait("curl", curl_ver, "/dev/null", "/dev/null", false) != 0)
        xdie("dump: curl not found or failed; install curl and retry");

    const char *base = output ? output
                     : (ids_out && strcmp(ids_out, "-")) ? ids_out : "phorge_dump";
    char *err_path = xmalloc(strlen(base) + 24);
    char *tmp_path = xmalloc(strlen(base) + 24);
    char *part_path = xmalloc(strlen(base) + 16);
    char *prov_path = xmalloc(strlen(base) + 24);
    sprintf(err_path, "%s.curl.log", base);
    sprintf(tmp_path, "%s.curl.tmp", base);
    sprintf(part_path, "%s.partial", base);
    sprintf(prov_path, "%s.dump.json", base);
    FILE *log = fopen(err_path, "wb");
    if (!log) xdie("dump: could not create %s: %s", err_path, strerror(errno));
    fclose(log);

    StrVec ids = {0};
    size_t count = 0;
    if (ids_file) {
        read_ids_file(ids_file, &ids);
        count = ids.len;
    } else {
        esearch(db, query, 0, 1, has_key ? key : NULL, tmp_path, err_path);
        char *text = ds_read_text(tmp_path);
        if (!text) xdie("dump: could not read ESearch response");
        DsKVList kvs = {0};
        ds_parse_json_scalars(text, &kvs);
        const char *cs = ds_kv_get(&kvs, "esearchresult.count");
        if (!cs) xdie("dump: ESearch response has no count (see %s)", err_path);
        count = strtoull(cs, NULL, 10);
        ds_kv_free(&kvs);
        free(text);
        pace(has_key);
    }
    if (dry_run && !ids_out) {
        printf("%zu\n", count);
        fprintf(stderr, "dump: %zu record(s) for %s (dry run)\n", count, ids_file ? ids_file : query);
        vec_free(&ids);
        free(err_path); free(tmp_path); free(part_path); free(prov_path); free(owned_query);
        return 0;
    }
    if (count == 0) xdie("dump: no records match the query");
    if (!ids_file) {
        size_t want = limit && limit < count ? limit : count;
        size_t page = want < 10000 ? want : 10000;
        for (size_t start = 0; start < want; start += page) {
            size_t n = start + page < want ? page : want - start;
            esearch(db, query, start, n, has_key ? key : NULL, tmp_path, err_path);
            char *text = ds_read_text(tmp_path);
            if (!text) xdie("dump: could not read ESearch response");
            DsKVList kvs = {0};
            ds_parse_json_scalars(text, &kvs);
            parse_ids(ds_kv_get(&kvs, "esearchresult.idlist"), &ids);
            ds_kv_free(&kvs);
            free(text);
            pace(has_key);
            fprintf(stderr, "dump: collected %zu/%zu UID(s)\r", ids.len, want);
        }
        fprintf(stderr, "\n");
    }
    size_t fetched = ids.len;
    if (limit && limit < fetched) {
        for (size_t i = limit; i < ids.len; i++) free(ids.data[i]);
        ids.len = limit;
        fetched = limit;
    }
    if (!fetched) xdie("dump: no UIDs to fetch");
    if (ids_out) {
        emit_ids(&ids, ids_out);
        fprintf(stderr, "dump: %zu UID(s) -> %s\n", fetched,
                !strcmp(ids_out, "-") ? "stdout" : ids_out);
        vec_free(&ids);
        free(err_path); free(tmp_path); free(part_path); free(prov_path); free(owned_query);
        return 0;
    }

    if (!force && access(output, F_OK) == 0)
        xdie("dump: %s exists (use --force)", output);
    FILE *partial = fopen(part_path, "wb");
    if (!partial) xdie("dump: could not create %s: %s", part_path, strerror(errno));

    for (size_t start = 0; start < fetched; start += batch) {
        size_t n = start + batch < fetched ? batch : fetched - start;
        char *joined = join_ids(&ids, start, n);
        efetch(db, joined, format, has_key ? key : NULL, tmp_path, err_path);
        free(joined);
        char *text = ds_read_text(tmp_path);
        if (!text) xdie("dump: could not read EFetch response");
        size_t len = strlen(text);
        size_t nonws = 0;
        while (nonws < len && isspace((unsigned char)text[nonws])) nonws++;
        if (nonws == len)
            fprintf(stderr, "dump: warning: batch %zu-%zu returned no sequence data\n", start + 1, start + n);
        else if (fwrite(text, 1, len, partial) != len)
            xdie("dump: write failed");
        free(text);
        fprintf(stderr, "dump: fetched %zu/%zu record(s)\r", start + n, fetched);
        pace(has_key);
    }
    fprintf(stderr, "\n");
    if (fclose(partial) != 0) xdie("dump: could not close %s", part_path);
    if (rename(part_path, output) != 0) xdie("dump: could not rename %s: %s", part_path, strerror(errno));

    write_provenance(prov_path, db, query, count, fetched, format, batch, has_key);
    fprintf(stderr, "dump: %zu record(s) -> %s (provenance: %s)\n", fetched, output, prov_path);

    vec_free(&ids);
    free(err_path); free(tmp_path); free(part_path); free(prov_path); free(owned_query);
    return 0;
}
