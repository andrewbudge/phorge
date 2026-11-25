#define _POSIX_C_SOURCE 200809L

#include "cli_help.h"
#include "cli_util.h"
#include "dataset.h"

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

typedef struct {
    char *group;
    char *taxon;
    char *symbol;
    char *gene_id;
    char *taxid;
    char *taxname;
    char *acc;
    char *header;
    char *seq;
} OutRec;

typedef struct {
    OutRec *data;
    size_t len;
    size_t cap;
} OutList;

static char *str_or_empty(const char *s) { return xstrdup(s ? s : ""); }

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n);
    if (!q) xdie("out of memory");
    return q;
}

static void out_push(OutList *list, const OutRec *r) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        OutRec *p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    OutRec copy = {
        str_or_empty(r->group), str_or_empty(r->taxon), str_or_empty(r->symbol),
        str_or_empty(r->gene_id), str_or_empty(r->taxid), str_or_empty(r->taxname),
        str_or_empty(r->acc), str_or_empty(r->header), str_or_empty(r->seq)
    };
    list->data[list->len++] = copy;
}

static void out_free_rec(OutRec *r) {
    free(r->group); free(r->taxon); free(r->symbol); free(r->gene_id);
    free(r->taxid); free(r->taxname); free(r->acc); free(r->header); free(r->seq);
}

static long out_find(const OutList *list, const char *group, const char *taxon) {
    for (size_t i = 0; i < list->len; i++)
        if (!strcmp(list->data[i].group, group) && !strcmp(list->data[i].taxon, taxon))
            return (long)i;
    return -1;
}

static int cmp_out(const void *a, const void *b) {
    const OutRec *x = a, *y = b;
    int c = strcmp(x->group, y->group);
    if (c) return c;
    c = strcmp(x->symbol, y->symbol);
    if (c) return c;
    return strcmp(x->header, y->header);
}

static char *first_token_dup(const char *s) {
    char *p = xstrdup(s ? s : "");
    char *space = strchr(p, ' ');
    if (space) *space = '\0';
    return p;
}

static char *clean_field(const char *s) {
    if (!s) return xstrdup("");
    char *out = xstrdup(s);
    for (char *p = out; *p; p++)
        if (*p == '\t' || *p == '\n' || *p == '\r') *p = ' ';
    return out;
}

static const char *preferred_fasta(const DsFastaList *fastas) {
    static const char *const order[] = {"cds.fna", "rna.fna", "gene.fna", "protein.faa", NULL};
    for (size_t j = 0; order[j]; j++) {
        for (size_t i = 0; i < fastas->len; i++)
            if (!strcmp(fastas->data[i].source, order[j])) return order[j];
    }
    return NULL;
}

static char *default_output(const char *input) {
    const char *slash = strrchr(input, '/');
    const char *base = slash ? slash + 1 : input;
    char *stem, *ext;
    stem_and_ext(base, &stem, &ext);
    free(ext);
    char *out = xmalloc(strlen(stem) + 12);
    sprintf(out, "%s_orthologs", stem);
    free(stem);
    return out;
}

int cmd_orthologs(int argc, char **argv) {
    const char *input = NULL, *outdir = NULL, *which = NULL;
    bool keep_all = false, id_taxid = false, force = false;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-o") || !strcmp(a, "--output")) outdir = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) outdir = a + 9;
        else if (!strcmp(a, "--fasta")) which = next_value(argc, argv, &i);
        else if (!strncmp(a, "--fasta=", 8)) which = a + 8;
        else if (!strcmp(a, "--keep-all")) keep_all = true;
        else if (!strcmp(a, "--id")) {
            const char *v = next_value(argc, argv, &i);
            if (!strcmp(v, "taxid")) id_taxid = true;
            else if (strcmp(v, "name")) xdie("orthologs: --id must be taxid or name");
        } else if (!strcmp(a, "--force")) force = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("orthologs"); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else if (input) xdie("unexpected argument: %s", a);
        else input = a;
    }
    if (!input) xdie("orthologs: missing dataset zip or directory");
    char *outdir_owned = NULL;
    if (!outdir) {
        outdir_owned = default_output(input);
        outdir = outdir_owned;
    }
    char *manifest_path = cli_join_path(outdir, "orthologs.tsv");
    if (!force && access(manifest_path, F_OK) == 0)
        xdie("orthologs: %s exists (use --force)", manifest_path);
    free(manifest_path);

    DsMetaList meta = {0};
    DsFastaList fastas = {0};
    ds_load_reports(input, &meta);
    ds_load_fasta(input, NULL, &fastas);
    if (fastas.len == 0) xdie("orthologs: no FASTA records found in %s", input);
    if (!which) which = preferred_fasta(&fastas);

    OutList rows = {0};
    size_t nselected = 0;
    for (size_t i = 0; i < fastas.len; i++) {
        const DsFastaRecord *rec = &fastas.data[i];
        if (which && strcmp(rec->source, which)) continue;
        nselected++;

        char *gene = ds_header_tag(rec->header, "GeneID");
        char *prot = ds_header_tag(rec->header, "protein_accession");
        char *acc = prot ? prot : first_token_dup(rec->header);
        const char *key = gene ? gene : acc;
        const DsMeta *m = ds_meta_find(&meta, key);
        if (!m && gene) m = ds_meta_find(&meta, acc);

        const char *group = (m && m->group_id && *m->group_id) ? m->group_id
                           : (gene && *gene) ? gene : acc;
        const char *symbol = (m && m->symbol && *m->symbol) ? m->symbol
                            : (gene && *gene) ? gene : group;
        const char *taxid = (m && m->taxid && *m->taxid) ? m->taxid : "";
        const char *taxname = (m && m->taxname && *m->taxname) ? m->taxname : "";
        char *organism = ds_header_tag(rec->header, "organism");
        const char *name = taxname[0] ? taxname : (organism ? organism : acc);
        const char *taxon = taxid[0] ? taxid : name;

        char *header;
        if (id_taxid) header = xstrdup(taxid[0] ? taxid : name);
        else header = xstrdup(name);
        if (keep_all) {
            char *with_acc = xmalloc(strlen(header) + strlen(acc) + 2);
            sprintf(with_acc, "%s %s", header, acc);
            free(header);
            header = with_acc;
        }

        OutRec r = {0};
        r.group = xstrdup(group);
        r.taxon = xstrdup(taxon);
        r.symbol = xstrdup(symbol);
        r.gene_id = xstrdup(gene ? gene : "");
        r.taxid = xstrdup(taxid);
        r.taxname = xstrdup(taxname);
        r.acc = xstrdup(acc);
        r.header = header;
        r.seq = xstrdup(rec->seq);

        if (!keep_all) {
            long idx = out_find(&rows, group, taxon);
            if (idx >= 0) {
                OutRec *prev = &rows.data[idx];
                if (strlen(r.seq) > strlen(prev->seq)) {
                    out_free_rec(prev);
                    *prev = r;
                } else {
                    out_free_rec(&r);
                }
                free(gene); free(prot); free(acc); free(organism);
                continue;
            }
        }
        out_push(&rows, &r);
        out_free_rec(&r);
        free(gene); free(prot); free(acc); free(organism);
    }
    if (rows.len == 0) xdie("orthologs: no records to write");

    qsort(rows.data, rows.len, sizeof(*rows.data), cmp_out);
    if (mkdir_p(outdir) != 0) xdie("Could not create %s", outdir);

    char *manifest = cli_join_path(outdir, "orthologs.tsv");
    FILE *mf = fopen(manifest, "wb");
    if (!mf) xdie("Could not create %s: %s", manifest, strerror(errno));
    fputs("group\tsymbol\tgene_id\ttaxid\ttaxon_name\taccession\toutput_file\n", mf);

    typedef struct { char *group; char *file; size_t count; } OutFile;
    OutFile *files = xmalloc(rows.len * sizeof(*files));
    size_t nfiles = 0;
    for (size_t i = 0; i < rows.len; i++) {
        OutRec *r = &rows.data[i];
        size_t fi = nfiles;
        for (size_t j = 0; j < nfiles; j++)
            if (!strcmp(files[j].group, r->group)) { fi = j; break; }
        if (fi < nfiles) {
            char *path = cli_join_path(outdir, files[fi].file);
            FILE *out = fopen(path, "ab");
            if (!out) xdie("Could not append %s: %s", path, strerror(errno));
            fprintf(out, ">%s\n%s\n", r->header, r->seq);
            fclose(out);
            free(path);
            files[fi].count++;
            char *c_group = clean_field(r->group), *c_sym = clean_field(r->symbol);
            char *c_gid = clean_field(r->gene_id), *c_tid = clean_field(r->taxid);
            char *c_name = clean_field(r->taxname), *c_acc = clean_field(r->acc);
            fprintf(mf, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", c_group, c_sym, c_gid,
                    c_tid, c_name, c_acc, files[fi].file);
            free(c_group); free(c_sym); free(c_gid); free(c_tid); free(c_name); free(c_acc);
            continue;
        }
        char *base = ds_sanitize_name(r->symbol);
        char *gid = ds_sanitize_name(r->group);
        char *stem = xmalloc(strlen(base) + strlen(gid) + 2);
        sprintf(stem, "%s_%s", base, gid);
        free(base);
        free(gid);
        char *file = NULL;
        for (unsigned suffix = 1;; suffix++) {
            size_t cap = strlen(stem) + 24;
            file = xrealloc(file, cap);
            if (suffix == 1) snprintf(file, cap, "%s.fa", stem);
            else snprintf(file, cap, "%s_%u.fa", stem, suffix);
            bool taken = false;
            for (size_t j = 0; j < nfiles; j++)
                if (!strcmp(files[j].file, file)) { taken = true; break; }
            if (!taken) break;
        }
        free(stem);
        files[nfiles].group = xstrdup(r->group);
        files[nfiles].file = xstrdup(file);
        files[nfiles].count = 0;
        nfiles++;
        char *path = cli_join_path(outdir, file);
        FILE *out = fopen(path, files[fi].count ? "ab" : "wb");
        if (!out) xdie("Could not create %s: %s", path, strerror(errno));
        fprintf(out, ">%s\n%s\n", r->header, r->seq);
        fclose(out);
        free(path);
        files[fi].count++;
        free(file);

        char *c_group = clean_field(r->group), *c_sym = clean_field(r->symbol);
        char *c_gid = clean_field(r->gene_id), *c_tid = clean_field(r->taxid);
        char *c_name = clean_field(r->taxname), *c_acc = clean_field(r->acc);
        fprintf(mf, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", c_group, c_sym, c_gid,
                c_tid, c_name, c_acc, files[fi].file);
        free(c_group); free(c_sym); free(c_gid); free(c_tid); free(c_name); free(c_acc);
    }
    fclose(mf);

    /* Deterministic stdout listing. */
    for (size_t i = 0; i < nfiles; i++)
        printf("%s\t%zu\n", files[i].file, files[i].count);
    fprintf(stderr, "orthologs: %zu group(s), %zu sequence(s) -> %s\n",
            nfiles, rows.len, outdir);

    for (size_t i = 0; i < nfiles; i++) { free(files[i].group); free(files[i].file); }
    free(files);
    free(manifest);
    for (size_t i = 0; i < rows.len; i++) out_free_rec(&rows.data[i]);
    free(rows.data);
    ds_fasta_free(&fastas);
    ds_meta_free(&meta);
    free(outdir_owned);
    return 0;
}
