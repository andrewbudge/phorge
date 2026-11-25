#define _POSIX_C_SOURCE 200809L

#include "cli_help.h"
#include "cli_util.h"
#include "dataset.h"
#include "phorge_c.h"

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static bool is_refseq(const char *acc) {
    static const char *const prefixes[] = {"NM_", "NR_", "NP_", "XM_", "XR_", "XP_",
                                           "NG_", "NC_", "NX_", "NT_", "NW_", "NZ_", "WP_"};
    for (size_t i = 0; i < sizeof(prefixes) / sizeof(*prefixes); i++)
        if (!strncmp(acc, prefixes[i], strlen(prefixes[i]))) return true;
    return false;
}

static uint64_t parse_u64_or_zero(const char *s) {
    if (!s || !*s) return 0;
    char *end;
    unsigned long long v = strtoull(s, &end, 10);
    return end == s ? 0 : (uint64_t)v;
}

static int cmp_accession(const void *a, const void *b) {
    return strcmp(((const Accession *)a)->accession, ((const Accession *)b)->accession);
}

static void push_record(const DsFastaRecord *rec, const DsMetaList *meta,
                        const char *query_taxid, TaxonGroup group, AccessionList *rows) {
    char *prot = ds_header_tag(rec->header, "protein_accession");
    char *gene = ds_header_tag(rec->header, "GeneID");
    char *organism = ds_header_tag(rec->header, "organism");
    char *name = ds_header_tag(rec->header, "name");
    char *first = xstrdup(rec->header);
    char *space = strchr(first, ' ');
    if (space) *space = '\0';
    const char *acc = prot ? prot : first;
    const char *key = gene ? gene : acc;

    const char *taxid_s = ds_meta_lookup(meta, key, DS_META_TAXID);
    const char *taxname = ds_meta_lookup(meta, key, DS_META_TAXNAME);
    const char *desc = ds_meta_lookup(meta, key, DS_META_DESC);
    char *fallback = NULL;
    if (!taxname) taxname = organism ? organism : "";
    if (!desc) {
        if (name) desc = name;
        else if (strchr(rec->header, '[')) {
            fallback = strndup(rec->header, (size_t)(strchr(rec->header, '[') - rec->header));
            desc = fallback;
        } else desc = rec->header;
    }

    Accession a = {0};
    a.accession = xstrdup(acc);
    a.taxon_name = xstrdup(taxname);
    a.taxid = parse_u64_or_zero(taxid_s);
    a.length = (uint64_t)strlen(rec->seq);
    a.gene_annotation = xstrdup(desc);
    a.refseq = is_refseq(acc) ? 1 : 0;
    a.source_db = xstrdup(a.refseq ? "RefSeq" : "GenBank");
    a.query_taxid = query_taxid ? parse_u64_or_zero(query_taxid) : a.taxid;
    a.taxon_group = group;
    a.taxonomic_outlier = 0;

    if (rows->len == rows->cap) {
        size_t cap = rows->cap ? rows->cap * 2 : 16;
        Accession *p = realloc(rows->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        rows->data = p;
        rows->cap = cap;
    }
    rows->data[rows->len++] = a;
    free(prot); free(gene); free(organism); free(name); free(first); free(fallback);
}

int cmd_report(int argc, char **argv) {
    const char *input = NULL, *output = NULL;
    const char *query_taxid = NULL, *group = "ingroup";

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-o") || !strcmp(a, "--output")) output = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) output = a + 9;
        else if (!strcmp(a, "--query-taxid")) query_taxid = next_value(argc, argv, &i);
        else if (!strcmp(a, "--group")) group = next_value(argc, argv, &i);
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("report"); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else if (input) xdie("unexpected argument: %s", a);
        else input = a;
    }
    if (!input) xdie("report: missing dataset directory or zip");
    TaxonGroup tg = TAXON_GROUP_INGROUP;
    if (!strcmp(group, "outgroup")) tg = TAXON_GROUP_OUTGROUP;
    else if (strcmp(group, "ingroup")) xdie("report: --group must be ingroup or outgroup");

    DsMetaList meta = {0};
    DsFastaList fastas = {0};
    ds_load_reports(input, &meta);
    ds_load_fasta(input, NULL, &fastas);
    if (fastas.len == 0) xdie("report: no FASTA files found in %s", input);

    AccessionList rows = {0};
    for (size_t i = 0; i < fastas.len; i++)
        push_record(&fastas.data[i], &meta, query_taxid, tg, &rows);

    /* One row per accession; keep the longest sequence. */
    AccessionList uniq = {0};
    for (size_t i = 0; i < rows.len; i++) {
        Accession *found = NULL;
        for (size_t j = 0; j < uniq.len; j++)
            if (!strcmp(uniq.data[j].accession, rows.data[i].accession)) { found = &uniq.data[j]; break; }
        if (!found) {
            if (uniq.len == uniq.cap) {
                size_t cap = uniq.cap ? uniq.cap * 2 : 16;
                Accession *p = realloc(uniq.data, cap * sizeof(*p));
                if (!p) xdie("out of memory");
                uniq.data = p;
                uniq.cap = cap;
            }
            uniq.data[uniq.len++] = rows.data[i];
        } else if (rows.data[i].length > found->length) {
            free(found->accession); free(found->taxon_name);
            free(found->gene_annotation); free(found->source_db);
            *found = rows.data[i];
        } else {
            free(rows.data[i].accession); free(rows.data[i].taxon_name);
            free(rows.data[i].gene_annotation); free(rows.data[i].source_db);
        }
    }
    free(rows.data);
    rows = uniq;

    qsort(rows.data, rows.len, sizeof(*rows.data), cmp_accession);
    char *tsv = render_query_tsv(&rows);
    FILE *out = stdout;
    if (output) {
        out = fopen(output, "wb");
        if (!out) xdie("Could not create %s: %s", output, strerror(errno));
    }
    fputs(tsv, out);
    if (output) fclose(out);
    else fprintf(stderr, "report: %zu accessions from %zu FASTA file(s)\n", rows.len, fastas.len);
    free(tsv);

    accession_list_free(&rows);
    ds_fasta_free(&fastas);
    ds_meta_free(&meta);
    return 0;
}
