#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_help.h"
#include "cli_util.h"

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

/* -------------------------------- curate --------------------------------- */

typedef struct {
    double gappiness;
    bool informative;
    bool constant;
} SiteStats;

typedef struct {
    double g;
    double count;
} GapPoint;

static double smart_gap_threshold(const double *gappiness, size_t n) {
    if (n == 0) return 1.0;
    size_t *counts = calloc(10001, sizeof(*counts));
    if (!counts) xdie("out of memory");
    for (size_t i = 0; i < n; i++) {
        uint32_t key = (uint32_t)round(gappiness[i] * 10000.0);
        if (key > 10000) key = 10000;
        counts[key]++;
    }
    GapPoint *arr = xmalloc(n * sizeof(*arr));
    size_t m = 0;
    for (int k = 10000; k >= 0; k--) {
        if (counts[k]) {
            arr[m].g = (double)k / 10000.0;
            arr[m].count = (double)counts[k];
            m++;
        }
    }
    free(counts);
    if (m < 2) {
        double result = m ? arr[0].g : 1.0;
        free(arr);
        return result;
    }
    double total = (double)n;
    double *cumsum = xmalloc(m * sizeof(*cumsum));
    double running = 0.0;
    for (size_t i = 0; i < m; i++) {
        running += arr[i].count / total;
        cumsum[i] = running;
    }
    size_t nslopes = 0;
    double *slopes = xmalloc(m * sizeof(*slopes));
    for (size_t i = 1; i < m; i++) {
        double delta = arr[i].g - arr[i - 1].g;
        if (fabs(delta) > DBL_EPSILON) {
            slopes[nslopes++] = fabs((cumsum[i] - cumsum[i - 1]) / delta);
        }
    }
    size_t usable = nslopes / 2;
    double result;
    if (usable < 2) {
        result = usable == 0 ? 1.0 : arr[0].g;
    } else {
        size_t best = 0;
        double best_diff = fabs(slopes[1] - slopes[0]);
        for (size_t i = 1; i + 1 < usable; i++) {
            double diff = fabs(slopes[i + 1] - slopes[i]);
            if (diff >= best_diff) { /* Rust max_by returns the last maximum */
                best_diff = diff;
                best = i;
            }
        }
        result = arr[best].g;
    }
    free(cumsum);
    free(slopes);
    free(arr);
    return result;
}

int cmd_curate(int argc, char **argv) {
    const char *keep = "ps";
    double gap_threshold = 0.9;
    const char *extension = "_curated";
    const char *outdir = NULL;
    const char **files = xmalloc((size_t)argc * sizeof(*files));
    size_t nfiles = 0;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-k") || !strcmp(a, "--keep")) keep = next_value(argc, argv, &i);
        else if (!strncmp(a, "--keep=", 7)) keep = a + 7;
        else if (!strcmp(a, "--gap-threshold")) gap_threshold = strtod(next_value(argc, argv, &i), NULL);
        else if (!strncmp(a, "--gap-threshold=", 16)) gap_threshold = strtod(a + 16, NULL);
        else if (!strcmp(a, "-e") || !strcmp(a, "--extension")) extension = next_value(argc, argv, &i);
        else if (!strncmp(a, "--extension=", 12)) extension = a + 12;
        else if (!strcmp(a, "-o") || !strcmp(a, "--output")) outdir = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) outdir = a + 9;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("curate"); free(files); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else files[nfiles++] = a;
    }
    if (nfiles == 0) xdie("curate: at least one input FASTA is required");
    if (outdir && mkdir_p(outdir) != 0) xdie("Could not create output directory");

    bool keep_p = strchr(keep, 'p') != NULL;
    bool keep_c = strchr(keep, 'c') != NULL;
    bool keep_s = strchr(keep, 's') != NULL;
    bool keep_g = strchr(keep, 'g') != NULL;
    bool kpi_only = keep_p && !keep_c;

    for (size_t f = 0; f < nfiles; f++) {
        SequenceList seqs = {0};
        size_t aln_len = 0;
        if (parse_fasta(files[f], true, &seqs, &aln_len) != 0) {
            xdie("Could not read alignment (file not found or sequences are not the same length): %s",
                 phorge_last_error());
        }
        size_t ntaxa = seqs.len;
        if (ntaxa == 0) xdie("empty alignment: %s", files[f]);

        bool unique[256] = {false};
        int ucount = 0;
        const char *first = seqs.data[0].seq;
        for (size_t i = 0; i < strlen(first) && ucount < 100; i++) {
            unsigned char b = (unsigned char)toupper((unsigned char)first[i]);
            if (first[i] == '-' || first[i] == '*') continue;
            if (!unique[b]) { unique[b] = true; ucount++; }
        }
        const char *gap_chars = ucount > 5 ? "-?*Xx" : "-?*XxNn";
#ifndef CURATE_STRCHR
        /* Build the gap-character membership table once. The 10M-character
           hot loop then does one array load instead of a strchr() call. */
        bool is_gap[256] = {false};
        for (const char *p = gap_chars; *p; p++) is_gap[(unsigned char)*p] = true;
#endif

        SiteStats *stats = xmalloc(aln_len * sizeof(*stats));
        for (size_t col = 0; col < aln_len; col++) {
            size_t freq[256] = {0};
            size_t gap_count = 0;
            for (size_t s = 0; s < ntaxa; s++) {
                unsigned char b = (unsigned char)seqs.data[s].seq[col];
#ifndef CURATE_STRCHR
                if (is_gap[b]) gap_count++;
#else
                /* A/B baseline kept for bench/run_simd_bench.sh. */
                if (strchr(gap_chars, b) != NULL) gap_count++;
#endif
                else freq[b]++;
            }
            double gappiness = round((double)gap_count / (double)ntaxa * 10000.0) / 10000.0;
            size_t qualifying = 0, distinct = 0;
            for (size_t c = 0; c < 256; c++) {
                if (freq[c]) distinct++;
                if (freq[c] >= 2) qualifying++;
            }
            stats[col].gappiness = gappiness;
            stats[col].informative = qualifying >= 2;
            stats[col].constant = distinct == 1 && qualifying == 1;
        }

        double *vals = xmalloc(aln_len * sizeof(*vals));
        for (size_t i = 0; i < aln_len; i++) vals[i] = stats[i].gappiness;
        double smart_threshold = keep_s ? smart_gap_threshold(vals, aln_len) : gap_threshold;
        free(vals);

        size_t *kept_cols = xmalloc(aln_len * sizeof(*kept_cols));
        size_t nkept = 0;
        for (size_t col = 0; col < aln_len; col++) {
            const SiteStats *st = &stats[col];
            bool gap_ok;
            if (keep_s) gap_ok = kpi_only ? st->gappiness <= smart_threshold
                                         : st->gappiness < smart_threshold;
            else if (keep_g) gap_ok = kpi_only ? st->gappiness <= gap_threshold
                                              : st->gappiness < gap_threshold;
            else gap_ok = true;
            bool class_ok = (keep_p && keep_c) ? (st->informative || st->constant)
                          : keep_p ? st->informative : true;
            if (gap_ok && class_ok) kept_cols[nkept++] = col;
        }

        char *stem, *ext;
        stem_and_ext(files[f], &stem, &ext);
        FILE *writer = stdout;
        if (outdir) {
            char path[1024];
            snprintf(path, sizeof(path), "%s/%s%s.%s", outdir, stem, extension, ext);
            writer = fopen(path, "wb");
            if (!writer) xdie("Could not create output file: %s", path);
        }
        char *trimmed = xmalloc(nkept + 1);
        for (size_t s = 0; s < ntaxa; s++) {
            const char *seq = seqs.data[s].seq;
            for (size_t k = 0; k < nkept; k++) trimmed[k] = seq[kept_cols[k]];
            trimmed[nkept] = '\0';
            fprintf(writer, ">%s\n%s\n", seqs.data[s].header, trimmed);
        }
        if (outdir) {
            fclose(writer);
            fprintf(stderr, "%s: %zu → %zu sites (%zu removed)\n",
                    stem, aln_len, nkept, aln_len - nkept);
        } else {
            fprintf(stderr, "Sites in:      %zu\nSites kept:    %zu\nSites removed: %zu\n",
                    aln_len, nkept, aln_len - nkept);
        }
        free(trimmed);
        free(stem);
        free(ext);
        free(stats);
        free(kept_cols);
        sequence_list_free(&seqs);
    }
    free(files);
    if (outdir && nfiles > 1) fprintf(stderr, "Done. Curated %zu files.\n", nfiles);
    return 0;
}

