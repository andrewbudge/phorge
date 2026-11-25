#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_util.h"

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cli_help.h"

typedef struct {
    void **data;
    size_t len;
    size_t cap;
} PtrList;

static void ptr_push(PtrList *list, void *p) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        void **data = realloc(list->data, cap * sizeof(*data));
        if (!data) xdie("out of memory");
        list->data = data;
        list->cap = cap;
    }
    list->data[list->len++] = p;
}

/* All uses are short numeric cells ("n/m", "x.x%"); fixed buffer is ample. */
static char *makef(const char *fmt, ...) {
    char buf[128];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    return xstrdup(buf);
}

static char *pct_cell(size_t numerator, size_t denominator) {
    if (denominator == 0) return xstrdup(numerator == 0 ? "NaN%" : "inf%");
    return makef("%.1f%%", (double)numerator / (double)denominator * 100.0);
}

typedef struct {
    PtrList rows; /* each row: array of ncols malloc'd cell strings */
    size_t ncols;
} Table;

static void table_add(Table *table, ...) {
    char **cells = xmalloc(table->ncols * sizeof(*cells));
    va_list ap;
    va_start(ap, table);
    for (size_t j = 0; j < table->ncols; j++) cells[j] = va_arg(ap, char *);
    va_end(ap);
    ptr_push(&table->rows, cells);
}

static void table_free(Table *table) {
    for (size_t i = 0; i < table->rows.len; i++) {
        char **cells = table->rows.data[i];
        for (size_t j = 0; j < table->ncols; j++) free(cells[j]);
        free(cells);
    }
    free(table->rows.data);
}

static void print_table(const Table *table, bool pretty) {
    if (!pretty) {
        for (size_t i = 0; i < table->rows.len; i++) {
            char **cells = table->rows.data[i];
            for (size_t j = 0; j < table->ncols; j++)
                printf("%s%s", j ? "\t" : "", cells[j]);
            putchar('\n');
        }
        return;
    }
    size_t *widths = xmalloc(table->ncols * sizeof(*widths));
    for (size_t j = 0; j < table->ncols; j++) widths[j] = 0;
    for (size_t i = 0; i < table->rows.len; i++) {
        char **cells = table->rows.data[i];
        for (size_t j = 0; j < table->ncols; j++) {
            size_t n = strlen(cells[j]);
            if (n > widths[j]) widths[j] = n;
        }
    }
    for (size_t i = 0; i < table->rows.len; i++) {
        char **cells = table->rows.data[i];
        for (size_t j = 0; j < table->ncols; j++) {
            if (j) fputs("  ", stdout);
            printf("%-*s", (int)widths[j], cells[j]);
        }
        putchar('\n');
    }
    free(widths);
}

static void strings_push(StringList *list, char *s) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        char **data = realloc(list->data, cap * sizeof(*data));
        if (!data) xdie("out of memory");
        list->data = data;
        list->cap = cap;
    }
    list->data[list->len++] = s;
}

static size_t count_fields(const char *s) {
    size_t n = 1;
    for (; *s; s++) if (*s == '\t') n++;
    return n;
}

/* Splits in place; returned pointers point into `line`, which must stay alive. */
static size_t split_fields(char *line, char ***out) {
    size_t cap = 8, n = 0;
    char **fields = xmalloc(cap * sizeof(*fields));
    char *start = line;
    for (char *p = line;; p++) {
        if (*p == '\t' || *p == '\0') {
            char saved = *p;
            *p = '\0';
            if (n == cap) {
                cap *= 2;
                char **q = realloc(fields, cap * sizeof(*fields));
                if (!q) xdie("out of memory");
                fields = q;
            }
            fields[n++] = start;
            if (saved == '\0') break;
            start = p + 1;
        }
    }
    *out = fields;
    return n;
}

static void read_lines(const char *path, StringList *lines) {
    FILE *f = fopen(path, "rb");
    if (!f) xdie("Error: could not open '%s': %s", path, strerror(errno));
    char *buf = NULL;
    size_t cap = 0;
    ssize_t n;
    while ((n = getline(&buf, &cap, f)) != -1) {
        while (n > 0 && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) buf[--n] = '\0';
        strings_push(lines, xstrdup(buf));
    }
    free(buf);
    fclose(f);
    if (lines->len == 0) xdie("Error: '%s' is empty", path);
}

static void coverage_by_taxa(const StringList *lines, size_t total_loci, bool pretty) {
    Table table = {.ncols = 4};
    table_add(&table, xstrdup("taxa"), xstrdup("loci_present"),
              xstrdup("loci_missing"), xstrdup("pct_missing"));
    for (size_t i = 1; i < lines->len; i++) {
        char **fields;
        size_t nf = split_fields(lines->data[i], &fields);
        /* Rust counts every field after the first, even when the row has more
           columns than the header; keep the same (wrapping) arithmetic. */
        size_t present = 0;
        for (size_t j = 1; j < nf; j++)
            if (strcmp(fields[j], "MISSING") != 0) present++;
        size_t missing = total_loci - present;
        table_add(&table, xstrdup(fields[0]),
                  makef("%zu/%zu", present, total_loci),
                  makef("%zu/%zu", missing, total_loci),
                  pct_cell(missing, total_loci));
        free(fields);
    }
    print_table(&table, pretty);
    table_free(&table);
}

static void coverage_by_loci(const StringList *lines, const char *header,
                             size_t total_loci, bool pretty) {
    size_t *counts = calloc(total_loci ? total_loci : 1, sizeof(*counts));
    if (!counts) xdie("out of memory");
    size_t total_taxa = 0;
    for (size_t i = 1; i < lines->len; i++) {
        total_taxa++;
        char **fields;
        size_t nf = split_fields(lines->data[i], &fields);
        for (size_t j = 1; j < nf && j <= total_loci; j++)
            if (strcmp(fields[j], "MISSING") != 0) counts[j - 1]++;
        free(fields);
    }

    Table table = {.ncols = 3};
    table_add(&table, xstrdup("loci"), xstrdup("appearance_count"), xstrdup("missing_pct"));
    char *header_copy = xstrdup(header);
    char **names;
    size_t nnames = split_fields(header_copy, &names);
    for (size_t i = 0; i < total_loci && i + 1 < nnames; i++) {
        size_t missing = total_taxa - counts[i];
        table_add(&table, xstrdup(names[i + 1]),
                  makef("%zu/%zu", counts[i], total_taxa),
                  pct_cell(missing, total_taxa));
    }
    free(names);
    free(header_copy);
    print_table(&table, pretty);
    table_free(&table);
    free(counts);
}

int cmd_coverage(int argc, char **argv) {
    bool taxa = false, loci = false, pretty = false;
    const char *path = NULL;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-t") || !strcmp(a, "--taxa")) taxa = true;
        else if (!strcmp(a, "-l") || !strcmp(a, "--loci")) loci = true;
        else if (!strcmp(a, "-p") || !strcmp(a, "--pretty")) pretty = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("coverage"); return 0; }
        else if (a[0] == '-' && a[1] && a[1] != '-' && strchr(a, '=') == NULL &&
                 strspn(a + 1, "tlp") == strlen(a + 1)) {
            for (const char *p = a + 1; *p; p++) {
                if (*p == 't') taxa = true;
                else if (*p == 'l') loci = true;
                else pretty = true;
            }
        } else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else if (!path) path = a;
        else xdie("unexpected argument: %s", a);
    }
    if (!path) xdie("Error: missing provenance TSV");
    if (!taxa && !loci) xdie("Error: choose a coverage view: -t/--taxa or -l/--loci");

    StringList lines = {0};
    read_lines(path, &lines);
    size_t total_loci = count_fields(lines.data[0]) - 1;

    if (taxa) coverage_by_taxa(&lines, total_loci, pretty);
    else coverage_by_loci(&lines, lines.data[0], total_loci, pretty);

    string_list_free(&lines);
    return 0;
}
