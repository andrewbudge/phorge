#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_help.h"
#include "cli_util.h"

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char **cells;
    size_t ncols;
} Row;

typedef struct {
    Row *data;
    size_t len;
    size_t cap;
} Rows;

static char *fmt_string(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);
    if (n < 0) xdie("formatting failed");
    char *s = malloc((size_t)n + 1);
    if (!s) xdie("out of memory");
    va_start(ap, fmt);
    vsnprintf(s, (size_t)n + 1, fmt, ap);
    va_end(ap);
    return s;
}

static void rows_add(Rows *rows, const char *const *cells, size_t ncols) {
    if (rows->len == rows->cap) {
        size_t cap = rows->cap ? rows->cap * 2 : 16;
        Row *p = realloc(rows->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        rows->data = p;
        rows->cap = cap;
    }
    Row *row = &rows->data[rows->len++];
    row->cells = malloc(ncols * sizeof(*row->cells));
    if (!row->cells) xdie("out of memory");
    for (size_t i = 0; i < ncols; i++) row->cells[i] = xstrdup(cells[i]);
    row->ncols = ncols;
}

static void rows_free(Rows *rows) {
    for (size_t i = 0; i < rows->len; i++) {
        for (size_t j = 0; j < rows->data[i].ncols; j++) free(rows->data[i].cells[j]);
        free(rows->data[i].cells);
    }
    free(rows->data);
    rows->data = NULL;
    rows->len = rows->cap = 0;
}

/* Rust's f64 `format!("{:.1}")` uses uppercase NaN. */
static char *fmt_pct(double v) {
    if (isnan(v)) return xstrdup("NaN");
    if (isinf(v)) return xstrdup(v > 0 ? "inf" : "-inf");
    return fmt_string("%.1f", v);
}

static void print_table(const Rows *rows, bool pretty) {
    if (rows->len == 0) return;
    if (!pretty) {
        for (size_t i = 0; i < rows->len; i++) {
            const Row *r = &rows->data[i];
            for (size_t j = 0; j < r->ncols; j++) {
                if (j) putchar('\t');
                fputs(r->cells[j], stdout);
            }
            putchar('\n');
        }
        return;
    }
    size_t ncols = rows->data[0].ncols;
    size_t *widths = calloc(ncols, sizeof(*widths));
    if (!widths) xdie("out of memory");
    for (size_t i = 0; i < rows->len; i++) {
        const Row *r = &rows->data[i];
        for (size_t j = 0; j < r->ncols && j < ncols; j++) {
            size_t w = strlen(r->cells[j]);
            if (w > widths[j]) widths[j] = w;
        }
    }
    for (size_t i = 0; i < rows->len; i++) {
        const Row *r = &rows->data[i];
        for (size_t j = 0; j < r->ncols; j++) {
            if (j) fputs("  ", stdout);
            printf("%-*s", (int)widths[j], r->cells[j]);
        }
        putchar('\n');
    }
    free(widths);
}

/* ---------------------------- getheaders -------------------------------- */
typedef struct {
    size_t gc;
    size_t missing;
    size_t total;
    size_t length;
} SeqStats;

static SeqStats calc_seq_stats(const char *seq) {
    SeqStats st = {0};
    for (const unsigned char *p = (const unsigned char *)seq; *p; p++) {
        st.total++;
        if (*p == 'G' || *p == 'C') st.gc++;
        if (phorge_is_missing_char(*p)) st.missing++;
    }
    st.length = strlen(seq);
    return st;
}

int cmd_stats(int argc, char **argv) {
    bool detailed = false, pretty = false;
    const char **files = calloc((size_t)argc + 1, sizeof(*files));
    if (!files) xdie("out of memory");
    size_t nfiles = 0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-d") || !strcmp(a, "--detailed")) {
            detailed = true;
        } else if (!strcmp(a, "-p") || !strcmp(a, "--pretty")) {
            pretty = true;
        } else if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            print_help("stats");
            free(files);
            return 0;
        } else if (a[0] == '-') {
            xdie("unrecognized argument: %s", a);
        } else {
            files[nfiles++] = a;
        }
    }

    Rows rows = {0};
    if (detailed) {
        const char *header[] = {"file", "header", "length", "gc_pct", "missing_pct"};
        rows_add(&rows, header, 5);
        for (size_t f = 0; f < nfiles; f++) {
            SequenceList seqs = {0};
            size_t length = 0;
            if (parse_fasta(files[f], false, &seqs, &length) != 0) {
                xdie("Error: failed to parse '%s': %s", files[f], phorge_last_error());
            }
            bool dna = is_dna(&seqs);
            char *basename = phorge_basename_dup(files[f]);
            for (size_t i = 0; i < seqs.len; i++) {
                SeqStats st = calc_seq_stats(seqs.data[i].seq);
                char *missing = fmt_pct(st.total ? (double)st.missing / (double)st.total * 100.0 : NAN);
                char *gc = dna ? fmt_pct((double)st.gc / (double)(st.total - st.missing) * 100.0)
                               : xstrdup("NA");
                const char *cells[] = {
                    basename,
                    seqs.data[i].header,
                    fmt_string("%zu", st.length),
                    gc,
                    missing,
                };
                rows_add(&rows, cells, 5);
                free(missing);
                free(gc);
                free((void *)cells[2]);
            }
            free(basename);
            sequence_list_free(&seqs);
        }
    } else {
        const char *header[] = {"file", "sequences", "length", "type", "gc_pct",
                                "missing_pct", "variable", "variable_pct",
                                "informative", "informative_pct"};
        rows_add(&rows, header, 10);
        for (size_t f = 0; f < nfiles; f++) {
            SequenceList seqs = {0};
            size_t length = 0;
            if (parse_fasta(files[f], false, &seqs, &length) != 0) {
                xdie("Error: failed to parse '%s': %s", files[f], phorge_last_error());
            }
            size_t num_sequences = seqs.len;
            bool dna = is_dna(&seqs);
            size_t gc_count = 0, missing_count = 0, total_chars = 0;
            for (size_t i = 0; i < seqs.len; i++) {
                SeqStats st = calc_seq_stats(seqs.data[i].seq);
                gc_count += st.gc;
                missing_count += st.missing;
                total_chars += st.total;
            }
            char *missing_pct = fmt_pct(total_chars ? (double)missing_count / (double)total_chars * 100.0 : NAN);
            char *gc_str = dna ? fmt_pct((double)gc_count / (double)(total_chars - missing_count) * 100.0)
                               : xstrdup("NA");
            const char *data_type = dna ? "DNA" : "AA";
            bool all_equal = true;
            for (size_t i = 0; i < seqs.len; i++) {
                if (strlen(seqs.data[i].seq) != length) { all_equal = false; break; }
            }
            char *basename = phorge_basename_dup(files[f]);
            char *length_str = all_equal ? fmt_string("%zu", length) : xstrdup("NA");
            char *variable = xstrdup("NA"), *variable_pct = xstrdup("NA");
            char *informative = xstrdup("NA"), *informative_pct = xstrdup("NA");
            if (all_equal) {
                size_t nvar = 0, ninfo = 0;
                unsigned long counts[256];
                size_t *lens = malloc(seqs.len * sizeof(*lens));
                if (!lens) xdie("out of memory");
                for (size_t s = 0; s < seqs.len; s++) lens[s] = strlen(seqs.data[s].seq);
                for (size_t i = 0; i < length; i++) {
                    memset(counts, 0, sizeof(counts));
                    for (size_t s = 0; s < seqs.len; s++) {
                        if (i >= lens[s]) continue;
                        unsigned char ch = (unsigned char)seqs.data[s].seq[i];
                        if (!phorge_is_missing_char(ch)) counts[ch]++;
                    }
                    size_t bases_present = 0, bases_twice = 0;
                    for (size_t c = 0; c < 256; c++) {
                        if (counts[c]) bases_present++;
                        if (counts[c] >= 2) bases_twice++;
                    }
                    if (bases_present > 1) nvar++;
                    if (bases_twice >= 2) ninfo++;
                }
                free(variable);
                variable = fmt_string("%zu", nvar);
                free(variable_pct);
                variable_pct = fmt_pct(length ? (double)nvar / (double)length * 100.0 : NAN);
                free(informative);
                informative = fmt_string("%zu", ninfo);
                free(informative_pct);
                informative_pct = fmt_pct(length ? (double)ninfo / (double)length * 100.0 : NAN);
                free(lens);
            }
            const char *cells[] = {
                basename,
                fmt_string("%zu", num_sequences),
                length_str,
                data_type,
                gc_str,
                missing_pct,
                variable,
                variable_pct,
                informative,
                informative_pct,
            };
            rows_add(&rows, cells, 10);
            for (size_t j = 0; j < 10; j++) {
                if (j == 0 || j == 3) continue; /* borrowed basename / static string */
                free((void *)cells[j]);
            }
            free(basename);
            sequence_list_free(&seqs);
        }
    }
    print_table(&rows, pretty);
    rows_free(&rows);
    free(files);
    return 0;
}
