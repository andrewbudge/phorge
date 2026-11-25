#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_util.h"

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

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

typedef struct {
    char *key;
    void *value;
} MapEntry;

typedef struct {
    MapEntry *items;
    size_t cap;
    size_t len;
} StrMap;

static uint64_t hash_str(const char *s) {
    uint64_t h = 1469598103934665603ULL;
    for (; *s; s++) {
        h ^= (unsigned char)*s;
        h *= 1099511628211ULL;
    }
    return h;
}

static void map_grow(StrMap *m) {
    size_t cap = m->cap ? m->cap * 2 : 64;
    MapEntry *items = calloc(cap, sizeof(*items));
    if (!items) xdie("out of memory");
    for (size_t i = 0; i < m->cap; i++) {
        if (!m->items[i].key) continue;
        size_t j = (size_t)(hash_str(m->items[i].key) & (cap - 1));
        while (items[j].key) j = (j + 1) & (cap - 1);
        items[j] = m->items[i];
    }
    free(m->items);
    m->items = items;
    m->cap = cap;
}

static void map_set(StrMap *m, const char *key, void *value) {
    if (!m->cap || (m->len + 1) * 10 >= m->cap * 7) map_grow(m);
    size_t j = (size_t)(hash_str(key) & (m->cap - 1));
    while (m->items[j].key) {
        if (!strcmp(m->items[j].key, key)) { m->items[j].value = value; return; }
        j = (j + 1) & (m->cap - 1);
    }
    m->items[j].key = xstrdup(key);
    m->items[j].value = value;
    m->len++;
}

static void *map_get(const StrMap *m, const char *key) {
    if (!m->cap) return NULL;
    size_t j = (size_t)(hash_str(key) & (m->cap - 1));
    while (m->items[j].key) {
        if (!strcmp(m->items[j].key, key)) return m->items[j].value;
        j = (j + 1) & (m->cap - 1);
    }
    return NULL;
}

static void map_free(StrMap *m) {
    for (size_t i = 0; i < m->cap; i++) free(m->items[i].key);
    free(m->items);
    m->items = NULL;
    m->cap = m->len = 0;
}

static char *substr_dup(const char *s, size_t n) {
    char *out = xmalloc(n + 1);
    memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

static char *first_token(const char *s) {
    while (*s && isspace((unsigned char)*s)) s++;
    const char *start = s;
    while (*s && !isspace((unsigned char)*s)) s++;
    return substr_dup(start, (size_t)(s - start));
}

static char *replace_spaces(const char *s) {
    char *out = xstrdup(s);
    for (char *p = out; *p; p++) if (*p == ' ') *p = '_';
    return out;
}

static char *join_path(const char *a, const char *b) {
    size_t n = strlen(a), m = strlen(b);
    bool slash = n && a[n - 1] == '/';
    char *out = xmalloc(n + (slash ? 0 : 1) + m + 1);
    if (slash) snprintf(out, n + m + 1, "%s%s", a, b);
    else snprintf(out, n + m + 2, "%s/%s", a, b);
    return out;
}

static char *dir_name(const char *path) {
    const char *slash = strrchr(path, '/');
    if (!slash) return xstrdup("");
    if (slash == path) return xstrdup("/");
    return substr_dup(path, (size_t)(slash - path));
}

static bool has_fasta_ext(const char *name) {
    const char *dot = strrchr(name, '.');
    if (!dot) return false;
    const char *ext = dot + 1;
    return !strcmp(ext, "fasta") || !strcmp(ext, "fa") ||
           !strcmp(ext, "fna") || !strcmp(ext, "fas");
}

static void dir_files(const char *dir, StringList *out) {
    DIR *d = opendir(dir);
    if (!d) xdie("Error: could not read directory '%s': %s", dir, strerror(errno));
    struct dirent *entry;
    while ((entry = readdir(d)) != NULL) {
        if (!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, "..")) continue;
        if (!has_fasta_ext(entry->d_name)) continue;
        strings_push(out, join_path(dir, entry->d_name));
    }
    closedir(d);
}

static void collect_targets(const PtrList *args, StringList *out) {
    for (size_t i = 0; i < args->len; i++) {
        const char *path = args->data[i];
        struct stat st;
        if (stat(path, &st) == 0 && S_ISDIR(st.st_mode)) dir_files(path, out);
        else if (stat(path, &st) == 0 && S_ISREG(st.st_mode)) strings_push(out, xstrdup(path));
        else fprintf(stderr, "Warning: '%s' is not a file or directory, skipping.\n", path);
    }
}

static void warn_refs_leaked(const PtrList *refs, const StringList *targets) {
    if (!refs->len) return;
    PtrList ref_dirs = {0};
    for (size_t i = 0; i < refs->len; i++) {
        struct stat st;
        if (stat(refs->data[i], &st) == 0 && S_ISREG(st.st_mode))
            ptr_push(&ref_dirs, dir_name(refs->data[i]));
    }
    PtrList leaked = {0};
    for (size_t i = 0; i < targets->len; i++) {
        char *dir = dir_name(targets->data[i]);
        for (size_t j = 0; j < ref_dirs.len; j++) {
            if (!strcmp(dir, ref_dirs.data[j])) { ptr_push(&leaked, targets->data[i]); break; }
        }
        free(dir);
    }
    if (leaked.len) {
        fprintf(stderr, "Warning: %zu target file(s) sit in a --refs directory: ", leaked.len);
        for (size_t i = 0; i < leaked.len; i++)
            fprintf(stderr, "%s%s", i ? ", " : "", (char *)leaked.data[i]);
        fputc('\n', stderr);
        fputs("  --refs takes ONE value; a glob like `--refs refs/*.fasta` leaves the rest as targets.\n", stderr);
        const char *suggest = ref_dirs.len ? ref_dirs.data[0] : ".";
        fprintf(stderr, "  Pass the directory instead: --refs %s\n", *suggest ? suggest : ".");
    }
    free(leaked.data);
    for (size_t i = 0; i < ref_dirs.len; i++) free(ref_dirs.data[i]);
    free(ref_dirs.data);
}

static size_t pool_references(const char *reference, const PtrList *refs, const char *path) {
    FILE *w = fopen(path, "wb");
    if (!w) xdie("Could not create pooled reference file");
    size_t counter = 0;

    if (reference) {
        SequenceList seqs = {0};
        if (parse_fasta(reference, false, &seqs, NULL) != 0)
            xdie("Failed to read reference FASTA: %s", phorge_last_error());
        for (size_t i = 0; i < seqs.len; i++) {
            char *gene = first_token(seqs.data[i].header);
            fprintf(w, ">%s::%zu\n%s\n", gene, counter, seqs.data[i].seq);
            free(gene);
            counter++;
        }
        sequence_list_free(&seqs);
    } else {
        StringList files = {0};
        collect_targets(refs, &files);
        for (size_t i = 0; i < files.len; i++) {
            char *stem, *ext;
            stem_and_ext(files.data[i], &stem, &ext);
            char *gene = replace_spaces(stem);
            SequenceList seqs = {0};
            if (parse_fasta(files.data[i], false, &seqs, NULL) != 0)
                xdie("Failed to read reference FASTA: %s", phorge_last_error());
            for (size_t j = 0; j < seqs.len; j++) {
                fprintf(w, ">%s::%zu\n%s\n", gene, counter, seqs.data[j].seq);
                counter++;
            }
            sequence_list_free(&seqs);
            free(gene);
            free(stem);
            free(ext);
        }
        string_list_free(&files);
    }
    fclose(w);
    return counter;
}

typedef struct {
    char *header;
    char *seq;
    char *filename;
} Target;

static void pool_targets(const StringList *files, const char *path,
                         StrMap *lookup, PtrList *records) {
    FILE *w = fopen(path, "wb");
    if (!w) xdie("Could not create pooled targets file");
    for (size_t i = 0; i < files->len; i++) {
        const char *filename = basename_of(files->data[i]);
        char *stem, *ext;
        stem_and_ext(files->data[i], &stem, &ext);
        char *organism = replace_spaces(stem);
        SequenceList seqs = {0};
        if (parse_fasta(files->data[i], false, &seqs, NULL) != 0)
            xdie("Failed to read FASTA: %s", phorge_last_error());
        for (size_t j = 0; j < seqs.len; j++) {
            char *seq_id = first_token(seqs.data[j].header);
            size_t keylen = strlen(organism) + 2 + strlen(seq_id) + 1;
            char *key = xmalloc(keylen);
            snprintf(key, keylen, "%s::%s", organism, seq_id);
            fprintf(w, ">%s\n%s\n", key, seqs.data[j].seq);

            Target *t = xmalloc(sizeof(*t));
            t->header = xstrdup(seqs.data[j].header);
            t->seq = xstrdup(seqs.data[j].seq);
            t->filename = xstrdup(filename);
            ptr_push(records, t);
            map_set(lookup, key, t);
            free(key);
            free(seq_id);
        }
        sequence_list_free(&seqs);
        free(organism);
        free(stem);
        free(ext);
    }
    fclose(w);
}

typedef struct {
    char *gene;
    char *target;
    double identity;
    double coverage;
    size_t tstart;
    size_t tend;
    bool minus;
} Hit;

static double parse_f64_or(const char *s, double fallback) {
    if (!*s) return fallback;
    char *end;
    double v = strtod(s, &end);
    return end == s ? fallback : v;
}

static size_t parse_size_or(const char *s, size_t fallback) {
    if (!*s || s[0] == '-') return fallback;
    errno = 0;
    char *end;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno || *end) return fallback;
    return (size_t)v;
}

static void split_fields(char *line, char ***out, size_t *count) {
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
    *count = n;
}

static void parse_hits(const char *path, PtrList *hits) {
    FILE *f = fopen(path, "rb");
    if (!f) xdie("Could not open MMseqs2 output");
    char *line = NULL;
    size_t cap = 0;
    ssize_t nread;
    while ((nread = getline(&line, &cap, f)) != -1) {
        while (nread > 0 && (line[nread - 1] == '\n' || line[nread - 1] == '\r'))
            line[--nread] = '\0';
        char **fields;
        size_t nf;
        split_fields(line, &fields, &nf);
        if (nf >= 8) {
            const char *sep = strstr(fields[0], "::");
            Hit *h = xmalloc(sizeof(*h));
            h->gene = substr_dup(fields[0], sep ? (size_t)(sep - fields[0]) : strlen(fields[0]));
            h->target = xstrdup(fields[1]);
            h->identity = parse_f64_or(fields[2], 0.0);
            h->coverage = parse_f64_or(fields[3], 0.0);
            size_t qstart = parse_size_or(fields[4], 1);
            size_t qend = parse_size_or(fields[5], 1);
            size_t tstart = parse_size_or(fields[6], 1);
            size_t tend = parse_size_or(fields[7], 1);
            h->tstart = tstart < tend ? tstart : tend;
            h->tend = tstart > tend ? tstart : tend;
            h->minus = qstart > qend;
            ptr_push(hits, h);
        }
        free(fields);
    }
    free(line);
    fclose(f);
}

static char complement_char(unsigned char c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'U': return 'A'; case 'R': return 'Y'; case 'Y': return 'R';
        case 'K': return 'M'; case 'M': return 'K'; case 'B': return 'V';
        case 'V': return 'B'; case 'D': return 'H'; case 'H': return 'D';
        default: return (char)c;
    }
}

static char *revcomp_slice(const char *s, size_t n) {
    char *out = xmalloc(n + 1);
    out[n] = '\0';
    for (size_t i = 0; i < n; i++) out[i] = complement_char((unsigned char)s[n - 1 - i]);
    return out;
}

static void remove_tree(const char *path) {
    DIR *d = opendir(path);
    if (d) {
        struct dirent *entry;
        while ((entry = readdir(d)) != NULL) {
            if (!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, "..")) continue;
            char *child = join_path(path, entry->d_name);
            struct stat st;
            if (lstat(child, &st) == 0 && S_ISDIR(st.st_mode)) remove_tree(child);
            else unlink(child);
            free(child);
        }
        closedir(d);
    }
    rmdir(path);
}

static bool log_contains(const char *path, const char *needle) {
    FILE *f = fopen(path, "rb");
    if (!f) return false;
    char *buf = NULL;
    size_t cap = 0;
    ssize_t n = getline(&buf, &cap, f);
    bool found = false;
    while (n != -1) {
        if (strstr(buf, needle)) { found = true; break; }
        n = getline(&buf, &cap, f);
    }
    free(buf);
    fclose(f);
    return found;
}

static char *num_d(double v) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%g", v);
    return xstrdup(buf);
}

static char *num_z(size_t v) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%zu", v);
    return xstrdup(buf);
}

int cmd_extract(int argc, char **argv) {
    const char *reference = NULL, *output = NULL, *max_mem = NULL;
    double min_identity = 0.7, coverage = 0.5, sensitivity = 5.7;
    size_t flank = 0, max_seqs = 300, kmer_length = 13;
    bool keep_intermediates = false;
    PtrList refs = {0}, target_args = {0};
    bool after_double_dash = false;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (after_double_dash) { ptr_push(&target_args, xstrdup(a)); continue; }
        if (!strcmp(a, "--")) { after_double_dash = true; continue; }
        if (!strcmp(a, "-r") || !strcmp(a, "--reference")) reference = next_value(argc, argv, &i);
        else if (!strncmp(a, "--reference=", 12)) reference = a + 12;
        else if (!strcmp(a, "--refs")) ptr_push(&refs, xstrdup(next_value(argc, argv, &i)));
        else if (!strncmp(a, "--refs=", 7)) ptr_push(&refs, xstrdup(a + 7));
        else if (!strcmp(a, "-o") || !strcmp(a, "--output")) output = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) output = a + 9;
        else if (!strcmp(a, "-m") || !strcmp(a, "--min-identity")) min_identity = strtod(next_value(argc, argv, &i), NULL);
        else if (!strncmp(a, "--min-identity=", 15)) min_identity = strtod(a + 15, NULL);
        else if (!strcmp(a, "-c") || !strcmp(a, "--coverage")) coverage = strtod(next_value(argc, argv, &i), NULL);
        else if (!strncmp(a, "--coverage=", 11)) coverage = strtod(a + 11, NULL);
        else if (!strcmp(a, "--flank")) flank = (size_t)strtoull(next_value(argc, argv, &i), NULL, 10);
        else if (!strncmp(a, "--flank=", 8)) flank = (size_t)strtoull(a + 8, NULL, 10);
        else if (!strcmp(a, "-s") || !strcmp(a, "--sensitivity")) sensitivity = strtod(next_value(argc, argv, &i), NULL);
        else if (!strncmp(a, "--sensitivity=", 14)) sensitivity = strtod(a + 14, NULL);
        else if (!strcmp(a, "--max-seqs")) max_seqs = (size_t)strtoull(next_value(argc, argv, &i), NULL, 10);
        else if (!strncmp(a, "--max-seqs=", 11)) max_seqs = (size_t)strtoull(a + 11, NULL, 10);
        else if (!strcmp(a, "--max-memory-limit")) max_mem = next_value(argc, argv, &i);
        else if (!strncmp(a, "--max-memory-limit=", 19)) max_mem = a + 19;
        else if (!strcmp(a, "--kmer-length")) kmer_length = (size_t)strtoull(next_value(argc, argv, &i), NULL, 10);
        else if (!strncmp(a, "--kmer-length=", 14)) kmer_length = (size_t)strtoull(a + 14, NULL, 10);
        else if (!strcmp(a, "--keep-intermediates")) keep_intermediates = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("extract"); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else ptr_push(&target_args, xstrdup(a));
    }

    if (!reference && !refs.len) xdie("Error: choose either --reference or --refs");
    if (reference && refs.len) xdie("Error: --reference and --refs are mutually exclusive");
    if (!output) xdie("Error: --output is required");
    if (!target_args.len) xdie("Error: no target FASTA files found.");

    char *version_args[] = {(char *)"mmseqs", (char *)"--version", NULL};
    if (cli_spawn_wait("mmseqs", version_args, "/dev/null", "/dev/null", false) < 0)
        xdie("Error: mmseqs not found. Make sure it is installed and in your PATH.");

    StringList target_files = {0};
    collect_targets(&target_args, &target_files);
    if (target_files.len == 0) xdie("Error: no target FASTA files found.");
    warn_refs_leaked(&refs, &target_files);

    if (mkdir_p(output) != 0) xdie("Could not create output directory");
    const char *tmp_base = getenv("TMPDIR");
    if (!tmp_base || !*tmp_base) tmp_base = "/tmp";
    size_t tmp_len = strlen(tmp_base) + 40;
    char *tmp_dir = xmalloc(tmp_len);
    snprintf(tmp_dir, tmp_len, "%s/phorge_extract_%ld", tmp_base, (long)getpid());
    if (mkdir_p(tmp_dir) != 0) xdie("Could not create temp directory");

    char *ref_path = join_path(tmp_dir, "pooled_reference.fasta");
    char *pool_path = join_path(tmp_dir, "pooled_targets.fasta");
    char *results_path = join_path(tmp_dir, "results.tsv");
    char *mmseqs_tmp = join_path(tmp_dir, "mmseqs_tmp");
    char *log_path = join_path(output, "mmseqs.log");

    size_t n_refs = pool_references(reference, &refs, ref_path);
    if (n_refs == 0) xdie("Error: no reference sequences found.");
    fprintf(stderr, "Pooled %zu reference sequence(s).\n", n_refs);
    fprintf(stderr, "Pooling %zu target files...\n", target_files.len);

    StrMap lookup = {0};
    PtrList records = {0};
    pool_targets(&target_files, pool_path, &lookup, &records);

    FILE *logf = fopen(log_path, "wb");
    if (!logf) xdie("Could not create mmseqs.log");
    fclose(logf);

    char *id_s = num_d(min_identity), *cov_s = num_d(coverage), *sens_s = num_d(sensitivity);
    char *maxseq_s = num_z(max_seqs), *kmer_s = num_z(kmer_length);
    fprintf(stderr, "Running MMseqs2 easy-search (min identity %s, min coverage %s)...\n", id_s, cov_s);

    char *args[64];
    size_t k = 0;
    args[k++] = (char *)"mmseqs";
    args[k++] = (char *)"easy-search";
    args[k++] = ref_path;
    args[k++] = pool_path;
    args[k++] = results_path;
    args[k++] = mmseqs_tmp;
    args[k++] = "--search-type"; args[k++] = "3";
    args[k++] = "-s"; args[k++] = sens_s;
    args[k++] = "--min-seq-id"; args[k++] = id_s;
    args[k++] = "-c"; args[k++] = cov_s;
    args[k++] = "--cov-mode"; args[k++] = "2";
    args[k++] = "--max-seqs"; args[k++] = maxseq_s;
    args[k++] = "--format-output";
    args[k++] = "query,target,fident,qcov,qstart,qend,tstart,tend";
    if (max_mem) { args[k++] = "--split-memory-limit"; args[k++] = (char *)max_mem; }
    args[k++] = "-k"; args[k++] = kmer_s;
    args[k++] = NULL;

    int status = cli_spawn_wait("mmseqs", args, log_path, log_path, false);
    if (status < 0) xdie("Failed to run mmseqs");
    if (status != 0) {
        fprintf(stderr, "Error: mmseqs easy-search failed. See %s\n", log_path);
        if (log_contains(log_path, "Cannot fit databases into"))
            fprintf(stderr,
                    "\nThe search ran out of memory building its k-mer index. "
                    "Splitting cannot help: the index table is sized by k, not by "
                    "the database. Retry with a shorter k than the current %zu, "
                    "e.g. --kmer-length 11.\n", kmer_length);
        exit(1);
    }

    fprintf(stderr, "Parsing results...\n");
    PtrList hits = {0};
    parse_hits(results_path, &hits);

    StrMap writers = {0};
    for (size_t i = 0; i < hits.len; i++) {
        Hit *h = hits.data[i];
        Target *t = map_get(&lookup, h->target);
        if (!t) {
            fprintf(stderr, "Warning: '%s' not found in lookup, skipping.\n", h->target);
            continue;
        }
        size_t seq_len = strlen(t->seq);
        size_t start = 0;
        if (h->tstart >= 1) {
            start = h->tstart - 1;
            start = start > flank ? start - flank : 0;
        }
        size_t end = seq_len;
        if (flank <= SIZE_MAX - h->tend) {
            size_t candidate = h->tend + flank;
            if (candidate < end) end = candidate;
        }
        if (end < start) end = start;
        char *extracted = h->minus ? revcomp_slice(t->seq + start, end - start)
                                   : substr_dup(t->seq + start, end - start);

        FILE **slot = map_get(&writers, h->gene);
        if (!slot) {
            size_t name_len = strlen(h->gene) + sizeof(".fasta");
            char *name = xmalloc(name_len);
            snprintf(name, name_len, "%s.fasta", h->gene);
            char *out_path = join_path(output, name);
            FILE *w = fopen(out_path, "wb");
            if (!w) xdie("Could not create output file");
            slot = xmalloc(sizeof(*slot));
            *slot = w;
            map_set(&writers, h->gene, slot);
            free(name);
            free(out_path);
        }
        fprintf(*slot,
                ">%s [gene=%s ident=%.3f cov=%.3f strand=%c src=%s %zu-%zu]\n%s\n",
                t->header, h->gene, h->identity, h->coverage, h->minus ? '-' : '+',
                t->filename, start, end, extracted);
        free(extracted);
    }
    fprintf(stderr, "Done. Extracted %zu gene(s) from %zu hits.\n", writers.len, hits.len);

    for (size_t i = 0; i < writers.cap; i++) {
        if (!writers.items[i].key) continue;
        FILE **slot = writers.items[i].value;
        fclose(*slot);
        free(slot);
    }
    map_free(&writers);
    for (size_t i = 0; i < hits.len; i++) {
        free(((Hit *)hits.data[i])->gene);
        free(((Hit *)hits.data[i])->target);
        free(hits.data[i]);
    }
    free(hits.data);
    map_free(&lookup);
    for (size_t i = 0; i < records.len; i++) {
        Target *t = records.data[i];
        free(t->header);
        free(t->seq);
        free(t->filename);
        free(t);
    }
    free(records.data);
    string_list_free(&target_files);
    for (size_t i = 0; i < refs.len; i++) free(refs.data[i]);
    free(refs.data);
    for (size_t i = 0; i < target_args.len; i++) free(target_args.data[i]);
    free(target_args.data);
    free(ref_path); free(pool_path); free(results_path); free(mmseqs_tmp); free(log_path);
    free(id_s); free(cov_s); free(sens_s); free(maxseq_s); free(kmer_s);

    if (!keep_intermediates) remove_tree(tmp_dir);
    else fprintf(stderr, "Intermediates kept at: %s\n", tmp_dir);
    free(tmp_dir);
    return 0;
}
