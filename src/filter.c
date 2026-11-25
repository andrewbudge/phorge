#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_help.h"
#include "cli_util.h"

#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct {
    char *key;
    size_t value;
} MapEntry;

typedef struct {
    MapEntry *items;
    size_t cap;
    size_t len;
} CountMap;

static void map_grow(CountMap *map) {
    size_t cap = map->cap ? map->cap * 2 : 64;
    MapEntry *items = calloc(cap, sizeof(*items));
    if (!items) xdie("out of memory");
    for (size_t i = 0; i < map->cap; i++) {
        if (!map->items[i].key) continue;
        size_t j = (size_t)(cli_hash_str(map->items[i].key) & (cap - 1));
        while (items[j].key) j = (j + 1) & (cap - 1);
        items[j] = map->items[i];
    }
    free(map->items);
    map->items = items;
    map->cap = cap;
}

static void map_set(CountMap *map, const char *key, size_t value) {
    if (!map->cap || (map->len + 1) * 10 >= map->cap * 7) map_grow(map);
    size_t j = (size_t)(cli_hash_str(key) & (map->cap - 1));
    while (map->items[j].key) {
        if (strcmp(map->items[j].key, key) == 0) {
            map->items[j].value = value;
            return;
        }
        j = (j + 1) & (map->cap - 1);
    }
    map->items[j].key = xstrdup(key);
    map->items[j].value = value;
    map->len++;
}

static const size_t *map_get(const CountMap *map, const char *key) {
    if (!map->cap) return NULL;
    size_t j = (size_t)(cli_hash_str(key) & (map->cap - 1));
    while (map->items[j].key) {
        if (strcmp(map->items[j].key, key) == 0) return &map->items[j].value;
        j = (j + 1) & (map->cap - 1);
    }
    return NULL;
}

static void map_free(CountMap *map) {
    for (size_t i = 0; i < map->cap; i++) free(map->items[i].key);
    free(map->items);
    map->items = NULL;
    map->cap = map->len = 0;
}

/* ---------------------------- table rows -------------------------------- */
/* ------------------------------- filter --------------------------------- */

int cmd_filter(int argc, char **argv) {
    bool have_threshold = false;
    double threshold = 0.0;
    const char *log_path = NULL;
    size_t min_loci = 0;
    bool have_min_loci = false;
    const char *input = NULL;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-m") || !strcmp(a, "--max-missing") || !strcmp(a, "--max_missing")) {
            if (++i >= argc) xdie("missing value for %s", a);
            threshold = strtod(argv[i], NULL);
            have_threshold = true;
        } else if (!strcmp(a, "-l") || !strcmp(a, "--log")) {
            if (++i >= argc) xdie("missing value for %s", a);
            log_path = argv[i];
        } else if (!strcmp(a, "-n") || !strcmp(a, "--min-loci") || !strcmp(a, "--min_loci")) {
            if (++i >= argc) xdie("missing value for %s", a);
            min_loci = strtoull(argv[i], NULL, 10);
            have_min_loci = true;
        } else if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            print_help("filter");
            return 0;
        } else if (a[0] == '-') {
            xdie("unrecognized argument: %s", a);
        } else {
            if (input) xdie("unexpected argument: %s", a);
            input = a;
        }
    }
    if (!input) xdie("filter: missing input FASTA");

    SequenceList seqs = {0};
    size_t length = 0;
    if (parse_fasta(input, true, &seqs, &length) != 0) {
        xdie("Could not read alignment (file not found or sequences are not the same length): %s",
            phorge_last_error());
    }

    CountMap loci_counts = {0};
    if (log_path) {
        FILE *f = fopen(log_path, "rb");
        if (!f) xdie("Could not open provenance TSV");
        char *line = NULL;
        size_t cap = 0;
        ssize_t nread;
        bool first = true;
        while ((nread = getline(&line, &cap, f)) != -1) {
            cli_strip_line_ending(line, &nread);
            if (first) { first = false; continue; }
            char *cols[64];
            size_t nc = 0;
            char *start = line;
            for (char *p = line;; p++) {
                if (*p == '\t' || *p == '\0') {
                    bool end = *p == '\0';
                    *p = '\0';
                    if (nc < 64) cols[nc] = start;
                    nc++;
                    if (end) break;
                    start = p + 1;
                }
            }
            if (nc < 2) continue;
            size_t count = 0;
            for (size_t i = 1; i < nc; i++) {
                if (strcmp(cols[i], "MISSING") != 0) count++;
            }
            map_set(&loci_counts, cols[0], count);
        }
        free(line);
        fclose(f);
    }

    size_t kept = 0;
    for (size_t i = 0; i < seqs.len; i++) {
        bool keep = true;
        const char *header = seqs.data[i].header;
        const char *seq = seqs.data[i].seq;
        if (have_threshold) {
            size_t missing = 0;
            for (const unsigned char *p = (const unsigned char *)seq; *p; p++) {
                if (phorge_is_missing_char(*p)) missing++;
            }
            double fraction = seq[0] ? (double)missing / (double)strlen(seq) : 0.0;
            if (fraction > threshold) keep = false;
        }
        if (have_min_loci) {
            const size_t *count = map_get(&loci_counts, header);
            if (count && *count < min_loci) keep = false;
        }
        if (keep) {
            printf(">%s\n%s\n", header, seq);
            kept++;
        }
    }
    fprintf(stderr, "Total taxa: %zu\nKept taxa: %zu\nDropped taxa: %zu\n",
            seqs.len, kept, seqs.len - kept);
    map_free(&loci_counts);
    sequence_list_free(&seqs);
    return 0;
}
