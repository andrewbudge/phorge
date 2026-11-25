#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_help.h"
#include "cli_util.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

typedef struct {
    void **ptrs;
    size_t len;
    size_t cap;
} Arena;

static void *arena_alloc(Arena *a, size_t n) {
    void *p = calloc(n ? n : 1, 1);
    if (!p) xdie("out of memory");
    if (a->len == a->cap) {
        size_t cap = a->cap ? a->cap * 2 : 64;
        void **q = realloc(a->ptrs, cap * sizeof(*q));
        if (!q) xdie("out of memory");
        a->ptrs = q;
        a->cap = cap;
    }
    a->ptrs[a->len++] = p;
    return p;
}

static char *arena_dup(Arena *a, const char *s) {
    size_t n = strlen(s) + 1;
    char *p = arena_alloc(a, n);
    memcpy(p, s, n);
    return p;
}

static void arena_free(Arena *a) {
    for (size_t i = 0; i < a->len; i++) free(a->ptrs[i]);
    free(a->ptrs);
    a->ptrs = NULL;
    a->len = a->cap = 0;
}

typedef struct {
    char **keys;
    size_t cap;
    size_t len;
} XSet;

static uint64_t xhash(const char *s) {
    uint64_t h = 1469598103934665603ULL;
    for (; *s; s++) {
        h ^= (unsigned char)*s;
        h *= 1099511628211ULL;
    }
    return h;
}

static void xset_grow(XSet *s) {
    size_t cap = s->cap ? s->cap * 2 : 64;
    char **keys = calloc(cap, sizeof(*keys));
    if (!keys) xdie("out of memory");
    for (size_t i = 0; i < s->cap; i++) {
        if (!s->keys[i]) continue;
        size_t j = (size_t)(xhash(s->keys[i]) & (cap - 1));
        while (keys[j]) j = (j + 1) & (cap - 1);
        keys[j] = s->keys[i];
    }
    free(s->keys);
    s->keys = keys;
    s->cap = cap;
}

static bool xset_contains(const XSet *s, const char *key) {
    if (!s->cap) return false;
    size_t j = (size_t)(xhash(key) & (s->cap - 1));
    while (s->keys[j]) {
        if (strcmp(s->keys[j], key) == 0) return true;
        j = (j + 1) & (s->cap - 1);
    }
    return false;
}

static void xset_insert(XSet *s, const char *key) {
    if (!s->cap || (s->len + 1) * 10 >= s->cap * 7) xset_grow(s);
    size_t j = (size_t)(xhash(key) & (s->cap - 1));
    while (s->keys[j]) {
        if (strcmp(s->keys[j], key) == 0) return;
        j = (j + 1) & (s->cap - 1);
    }
    s->keys[j] = xstrdup(key);
    s->len++;
}

static void xset_free(XSet *s) {
    for (size_t i = 0; i < s->cap; i++) free(s->keys[i]);
    free(s->keys);
    s->keys = NULL;
    s->cap = s->len = 0;
}

typedef struct {
    char *key;
    size_t value;
} XEntry;

typedef struct {
    XEntry *items;
    size_t cap;
    size_t len;
} XMap;

static void xmap_grow(XMap *m) {
    size_t cap = m->cap ? m->cap * 2 : 64;
    XEntry *items = calloc(cap, sizeof(*items));
    if (!items) xdie("out of memory");
    for (size_t i = 0; i < m->cap; i++) {
        if (!m->items[i].key) continue;
        size_t j = (size_t)(xhash(m->items[i].key) & (cap - 1));
        while (items[j].key) j = (j + 1) & (cap - 1);
        items[j] = m->items[i];
    }
    free(m->items);
    m->items = items;
    m->cap = cap;
}

static void xmap_set(XMap *m, const char *key, size_t value) {
    if (!m->cap || (m->len + 1) * 10 >= m->cap * 7) xmap_grow(m);
    size_t j = (size_t)(xhash(key) & (m->cap - 1));
    while (m->items[j].key) {
        if (strcmp(m->items[j].key, key) == 0) {
            m->items[j].value = value;
            return;
        }
        j = (j + 1) & (m->cap - 1);
    }
    m->items[j].key = xstrdup(key);
    m->items[j].value = value;
    m->len++;
}

static const size_t *xmap_get(const XMap *m, const char *key) {
    if (!m->cap) return NULL;
    size_t j = (size_t)(xhash(key) & (m->cap - 1));
    while (m->items[j].key) {
        if (strcmp(m->items[j].key, key) == 0) return &m->items[j].value;
        j = (j + 1) & (m->cap - 1);
    }
    return NULL;
}

static void xmap_free(XMap *m) {
    for (size_t i = 0; i < m->cap; i++) free(m->items[i].key);
    free(m->items);
    m->items = NULL;
    m->cap = m->len = 0;
}

/* Growable char buffer; caller frees buf_take() results. */
typedef struct {
    char *data;
    size_t len;
    size_t cap;
} XBuf;

static void xbuf_reserve(XBuf *b, size_t need) {
    if (b->cap - b->len >= need) return;
    size_t cap = b->cap ? b->cap : 64;
    while (cap - b->len < need) cap *= 2;
    char *p = realloc(b->data, cap);
    if (!p) xdie("out of memory");
    b->data = p;
    b->cap = cap;
}

static void xbuf_puts(XBuf *b, const char *s) {
    size_t n = strlen(s);
    xbuf_reserve(b, n + 1);
    memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = '\0';
}

static void xbuf_putc(XBuf *b, char c) {
    xbuf_reserve(b, 2);
    b->data[b->len++] = c;
    b->data[b->len] = '\0';
}

static void xbuf_repeat(XBuf *b, const char *s, size_t count) {
    size_t n = strlen(s);
    if (!n || !count) return;
    xbuf_reserve(b, n * count + 1);
    if (n == 1) {
        memset(b->data + b->len, (unsigned char)*s, count);
        b->len += count;
    } else {
        for (size_t i = 0; i < count; i++) {
            memcpy(b->data + b->len, s, n);
            b->len += n;
        }
    }
    b->data[b->len] = '\0';
}

typedef struct {
    const char *path;
    SequenceList seqs;
    size_t length;
    bool aa;
} Gene;

typedef struct {
    size_t taxon;
    size_t count;
} Coverage;

static bool detect_aa(const char *seq) {
    return strpbrk(seq, "EFILPQ") != NULL;
}

static void lower_match(const char *s, char *out, size_t cap) {
    size_t n = 0;
    for (; *s && n + 1 < cap; s++, n++) out[n] = (char)tolower((unsigned char)*s);
    out[n] = '\0';
}

static void write_provenance(const char *dest, const char *taxa_file,
                             const StringList *taxa, const Gene *genes, size_t ngenes,
                             const ssize_t **matches) {
    FILE *f = dest ? fopen(dest, "wb") : stdout;
    if (!f) xdie("Failed to create provenance log file");
    XBuf row = {0};
    xbuf_puts(&row, taxa_file);
    for (size_t g = 0; g < ngenes; g++) {
        xbuf_putc(&row, '\t');
        xbuf_puts(&row, basename_of(genes[g].path));
    }
    xbuf_putc(&row, '\n');
    fwrite(row.data, 1, row.len, f);
    row.len = 0;
    for (size_t t = 0; t < taxa->len; t++) {
        xbuf_puts(&row, taxa->data[t]);
        for (size_t g = 0; g < ngenes; g++) {
            xbuf_putc(&row, '\t');
            ssize_t h = matches[g][t];
            xbuf_puts(&row, h >= 0 ? genes[g].seqs.data[h].header : "MISSING");
        }
        xbuf_putc(&row, '\n');
        fwrite(row.data, 1, row.len, f);
        row.len = 0;
    }
    free(row.data);
    if (dest) fclose(f);
}

static char *resolve_log_path(Arena *arena, const char *path) {
    struct stat st;
    if (stat(path, &st) == 0 && S_ISDIR(st.st_mode)) {
        size_t n = strlen(path) + strlen("/phorge_concat.log") + 1;
        char *p = arena_alloc(arena, n);
        snprintf(p, n, "%s/phorge_concat.log", path);
        return p;
    }
    return arena_dup(arena, path);
}

int cmd_concat(int argc, char **argv) {
    Arena arena = {0};
    const char *alias = NULL;
    const char *format = "FASTA";
    const char *missing = NULL;
    const char *partitions = "raxml";
    const char *log = NULL;
    bool dry_run = false;
    const char **files = xmalloc((size_t)argc * sizeof(*files));
    size_t nfiles = 0;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-a") || !strcmp(a, "--alias")) alias = next_value(argc, argv, &i);
        else if (!strncmp(a, "--alias=", 8)) alias = a + 8;
        else if (!strcmp(a, "-f") || !strcmp(a, "--format")) format = next_value(argc, argv, &i);
        else if (!strncmp(a, "--format=", 9)) format = a + 9;
        else if (!strcmp(a, "-m") || !strcmp(a, "--missing")) missing = next_value(argc, argv, &i);
        else if (!strncmp(a, "--missing=", 10)) missing = a + 10;
        else if (!strcmp(a, "-p") || !strcmp(a, "--partitions")) partitions = next_value(argc, argv, &i);
        else if (!strncmp(a, "--partitions=", 13)) partitions = a + 13;
        else if (!strcmp(a, "-l") || !strcmp(a, "--log")) log = next_value(argc, argv, &i);
        else if (!strncmp(a, "--log=", 6)) log = a + 6;
        else if (!strcmp(a, "--dry-run")) dry_run = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("concat"); free(files); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else files[nfiles++] = a;
    }
    if (nfiles == 0) xdie("concat: at least one alignment file is required");

    if (alias) {
        struct stat st;
        if (stat(alias, &st) == 0 && S_ISDIR(st.st_mode)) {
            xdie("Error: '%s' is a directory. -a/--alias expects a plain text file with one taxon name per line.", alias);
        }
        TextInput alias_input;
        if (text_input_open(alias, &alias_input) == 0) {
            char *buf = NULL;
            size_t bcap = 0;
            ssize_t n = text_input_read_line(&alias_input, &buf, &bcap);
            if (n > 0) {
                char *p = buf;
                while (*p && isspace((unsigned char)*p)) p++;
                if (*p == '>') {
                    free(buf);
                    text_input_close(&alias_input);
                    xdie("Error: '%s' looks like a FASTA file. -a/--alias expects a plain text taxa list (one name per line, no '>' headers).", alias);
                }
            }
            free(buf);
            text_input_close(&alias_input);
        }
    }

    Gene *genes = xmalloc(nfiles * sizeof(*genes));
    for (size_t g = 0; g < nfiles; g++) {
        genes[g].path = files[g];
        genes[g].seqs.data = NULL;
        genes[g].seqs.len = genes[g].seqs.cap = 0;
        genes[g].length = 0;
        if (parse_fasta(files[g], true, &genes[g].seqs, &genes[g].length) != 0) {
            xdie("Failed to parse fasta file: %s", phorge_last_error());
        }
    }

    bool smart = alias != NULL;
    StringList taxa = {0};
    if (smart) {
        if (load_taxa_list(alias, &taxa) != 0) xdie("Failed to load taxa list: %s", phorge_last_error());
    } else {
        XSet seen = {0};
        for (size_t g = 0; g < nfiles; g++) {
            for (size_t h = 0; h < genes[g].seqs.len; h++) {
                const char *hdr = genes[g].seqs.data[h].header;
                if (!xset_contains(&seen, hdr)) {
                    xset_insert(&seen, hdr);
                    char *copy = xstrdup(hdr);
                    if (taxa.len == taxa.cap) {
                        size_t cap = taxa.cap ? taxa.cap * 2 : 16;
                        char **p = realloc(taxa.data, cap * sizeof(*p));
                        if (!p) xdie("out of memory");
                        taxa.data = p;
                        taxa.cap = cap;
                    }
                    taxa.data[taxa.len++] = copy;
                }
            }
        }
        xset_free(&seen);
    }
    size_t ntaxa = taxa.len;

    ssize_t **matches = xmalloc(nfiles * sizeof(*matches));
    for (size_t g = 0; g < nfiles; g++) {
        matches[g] = xmalloc((ntaxa ? ntaxa : 1) * sizeof(**matches));
        for (size_t t = 0; t < ntaxa; t++) matches[g][t] = -1;
        if (smart) {
            size_t *order = arena_alloc(&arena, ntaxa * sizeof(*order));
            for (size_t t = 0; t < ntaxa; t++) order[t] = t;
            for (size_t i = 1; i < ntaxa; i++) { /* stable insertion sort by descending length */
                size_t cur = order[i];
                size_t j = i;
                while (j > 0) {
                    size_t prev = order[j - 1];
                    size_t lp = strlen(taxa.data[prev]), lc = strlen(taxa.data[cur]);
                    if (lp > lc || (lp == lc && prev < cur)) break;
                    order[j] = prev;
                    j--;
                }
                order[j] = cur;
            }
            char **tlower = arena_alloc(&arena, ntaxa * sizeof(*tlower));
            for (size_t t = 0; t < ntaxa; t++) tlower[t] = str_lower(taxa.data[t]);
            char **hlower = xmalloc((genes[g].seqs.len ? genes[g].seqs.len : 1) * sizeof(*hlower));
            for (size_t h = 0; h < genes[g].seqs.len; h++) hlower[h] = str_lower(genes[g].seqs.data[h].header);
            XSet claimed = {0};
            for (size_t oi = 0; oi < ntaxa; oi++) {
                size_t t = order[oi];
                for (size_t h = 0; h < genes[g].seqs.len; h++) {
                    const char *hdr = genes[g].seqs.data[h].header;
                    if (xset_contains(&claimed, hdr)) continue;
                    if (strstr(hlower[h], tlower[t]) != NULL) {
                        xset_insert(&claimed, hdr);
                        matches[g][t] = (ssize_t)h;
                        break;
                    }
                }
            }
            for (size_t h = 0; h < genes[g].seqs.len; h++) free(hlower[h]);
            free(hlower);
            for (size_t t = 0; t < ntaxa; t++) free(tlower[t]);
            xset_free(&claimed);
        } else {
            XMap by_header = {0};
            for (size_t h = 0; h < genes[g].seqs.len; h++) {
                xmap_set(&by_header, genes[g].seqs.data[h].header, h);
            }
            for (size_t t = 0; t < ntaxa; t++) {
                const size_t *idx = xmap_get(&by_header, taxa.data[t]);
                matches[g][t] = idx ? (ssize_t)*idx : -1;
            }
            xmap_free(&by_header);
        }
    }

    if (smart && !dry_run && !log) {
        xdie("Error: -a/--alias requires -l/--log to write the provenance TSV.\n"
             "  To preview matching without writing output, use --dry-run.");
    }

    bool *aa = arena_alloc(&arena, nfiles * sizeof(*aa));
    bool has_dna = false, has_aa = false;
    for (size_t g = 0; g < nfiles; g++) {
        ssize_t first = -1;
        for (size_t t = 0; t < ntaxa && first < 0; t++) {
            if (matches[g][t] >= 0) first = matches[g][t];
        }
        aa[g] = first >= 0 ? detect_aa(genes[g].seqs.data[first].seq) : false;
        if (aa[g]) has_aa = true; else has_dna = true;
    }
    bool mixed = has_dna && has_aa;

    if (dry_run) {
        fprintf(stderr, "=== Dry Run: Matching Summary ===\n\n");
        fprintf(stderr, "%-40s %8s %8s %8s\n", "Gene", "Seqs", "Matched", "Missing");
        fprintf(stderr, "%s\n", "--------------------------------------------------------------------");
        for (size_t g = 0; g < nfiles; g++) {
            size_t matched_count = 0;
            for (size_t t = 0; t < ntaxa; t++) if (matches[g][t] >= 0) matched_count++;
            fprintf(stderr, "%-40s %8zu %8zu %8zu\n", basename_of(genes[g].path),
                    genes[g].seqs.len, matched_count, ntaxa - matched_count);
        }
        fprintf(stderr, "\n%-40s %8s\n", "Taxon", "Genes");
        fprintf(stderr, "%s\n", "--------------------------------------------------");
        Coverage *cov = arena_alloc(&arena, (ntaxa ? ntaxa : 1) * sizeof(*cov));
        for (size_t t = 0; t < ntaxa; t++) {
            cov[t].taxon = t;
            size_t count = 0;
            for (size_t g = 0; g < nfiles; g++) if (matches[g][t] >= 0) count++;
            cov[t].count = count;
        }
        for (size_t i = 1; i < ntaxa; i++) { /* stable descending by count */
            Coverage cur = cov[i];
            size_t j = i;
            while (j > 0 && (cov[j - 1].count < cur.count ||
                             (cov[j - 1].count == cur.count && cov[j - 1].taxon > cur.taxon))) {
                cov[j] = cov[j - 1];
                j--;
            }
            cov[j] = cur;
        }
        for (size_t i = 0; i < ntaxa; i++) {
            fprintf(stderr, "%-40s %8zu/%zu\n", taxa.data[cov[i].taxon], cov[i].count, nfiles);
        }
        size_t full = 0, any = 0, none = 0;
        for (size_t i = 0; i < ntaxa; i++) {
            if (cov[i].count == nfiles) full++;
            if (cov[i].count > 0) any++;
            if (cov[i].count == 0) none++;
        }
        fprintf(stderr, "\n=== Summary ===\n");
        fprintf(stderr, "Taxa: %zu total, %zu with all genes, %zu with some, %zu with none\n",
                ntaxa, full, any, none);
        fprintf(stderr, "Genes: %zu\n", nfiles);
        if (smart) {
            char *dest = log ? resolve_log_path(&arena, log) : NULL;
            write_provenance(dest, alias, &taxa, genes, nfiles, (const ssize_t **)matches);
            if (dest) fprintf(stderr, "\nTentative provenance log written to: %s\n", dest);
            else fprintf(stderr, "\n=== Tentative Provenance TSV ===\n");
        }
        goto cleanup;
    }

    const char *mc;
    if (missing) mc = missing;
    else if (mixed) mc = "?";
    else mc = has_aa ? "X" : "N";

    char fmt[64], pfmt[64];
    lower_match(format, fmt, sizeof(fmt));
    if (!strcmp(fmt, "nexus") || !strcmp(fmt, "n") || !strcmp(fmt, "nex")) {
        const char *dtype = mixed ? "MIXED" : has_aa ? "PROTEIN" : "DNA";
        const char *nmissing = missing ? missing : (mixed ? "?" : (has_aa ? "X" : "N"));
        size_t total = 0;
        for (size_t g = 0; g < nfiles; g++) total += genes[g].length;
        printf("#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=%zu NCHAR=%zu;\n", ntaxa, total);
        printf("  FORMAT DATATYPE=%s MISSING=%s GAP=-;\n  MATRIX\n", dtype, nmissing);
        XBuf matrix = {0};
        for (size_t t = 0; t < ntaxa; t++) {
            matrix.len = 0;
            for (size_t g = 0; g < nfiles; g++) {
                ssize_t h = matches[g][t];
                if (h >= 0) xbuf_puts(&matrix, genes[g].seqs.data[h].seq);
                else xbuf_repeat(&matrix, mc, genes[g].length);
            }
            printf("  %s    %s\n", taxa.data[t], matrix.data ? matrix.data : "");
        }
        free(matrix.data);
        printf(";\nEND;\nBEGIN SETS;\n");
        size_t pos = 1;
        for (size_t g = 0; g < nfiles; g++) {
            printf("  CHARSET %s = %zu-%zu;\n", basename_of(genes[g].path), pos, pos + genes[g].length - 1);
            pos += genes[g].length;
        }
        printf("END;\n");
    } else {
        XBuf row = {0};
        for (size_t t = 0; t < ntaxa; t++) {
            row.len = 0;
            for (size_t g = 0; g < nfiles; g++) {
                ssize_t h = matches[g][t];
                if (h >= 0) xbuf_puts(&row, genes[g].seqs.data[h].seq);
                else xbuf_repeat(&row, mc, genes[g].length);
            }
            printf(">%s\n%s\n", taxa.data[t], row.data ? row.data : "");
        }
        free(row.data);
    }

    lower_match(partitions, pfmt, sizeof(pfmt));
    bool pnexus = !strcmp(pfmt, "nexus") || !strcmp(pfmt, "n") || !strcmp(pfmt, "nex");
    size_t pos = 1;
    for (size_t g = 0; g < nfiles; g++) {
        size_t end = pos + genes[g].length - 1;
        if (pnexus) fprintf(stderr, "CHARSET %s = %zu-%zu;\n", basename_of(genes[g].path), pos, end);
        else fprintf(stderr, "%s, %s = %zu-%zu\n", aa[g] ? "WAG" : "DNA", basename_of(genes[g].path), pos, end);
        pos += genes[g].length;
    }

    if (smart) {
        char *dest = resolve_log_path(&arena, log);
        fprintf(stderr, "Provenance log written to: %s\n", dest);
        write_provenance(dest, alias, &taxa, genes, nfiles, (const ssize_t **)matches);
    }

cleanup:
    for (size_t g = 0; g < nfiles; g++) {
        sequence_list_free(&genes[g].seqs);
        free(matches[g]);
    }
    free(genes);
    free(matches);
    string_list_free(&taxa);
    free(files);
    arena_free(&arena);
    return 0;
}
