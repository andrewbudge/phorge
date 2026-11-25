#define _POSIX_C_SOURCE 200809L

#include "cli_help.h"
#include "cli_util.h"

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

typedef struct {
    char *out;
    char *rel;
    unsigned long long bytes;
} Harvested;

typedef struct {
    char **data;
    size_t len;
    size_t cap;
} PathList;

static void path_push(PathList *list, const char *rel) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        char **p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len++] = xstrdup(rel);
}

static int cmp_strptr(const void *a, const void *b) {
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

static bool ends_with_ci(const char *s, const char *suffix) {
    size_t n = strlen(s), m = strlen(suffix);
    if (m > n) return false;
    for (size_t i = 0; i < m; i++) {
        if (tolower((unsigned char)s[n - m + i]) != tolower((unsigned char)suffix[i]))
            return false;
    }
    return true;
}

static bool is_fasta(const char *name) {
    if (ends_with_ci(name, ".gz")) {
        size_t n = strlen(name) - 3;
        char *stem = xmalloc(n + 1);
        memcpy(stem, name, n);
        stem[n] = '\0';
        bool ok = is_fasta(stem);
        free(stem);
        return ok;
    }
    static const char *const exts[] = {".fna", ".faa", ".fasta", ".fa", ".ffn", ".fas"};
    for (size_t i = 0; i < sizeof(exts) / sizeof(*exts); i++)
        if (ends_with_ci(name, exts[i])) return true;
    return false;
}

static bool same_real_path(const char *a, const char *b) {
    char *ra = realpath(a, NULL), *rb = realpath(b, NULL);
    bool same = ra && rb && !strcmp(ra, rb);
    free(ra);
    free(rb);
    return same;
}

static void walk(const char *root, const char *rel, const char *skip, PathList *out) {
    char *dir = rel[0] ? cli_join_path(root, rel) : xstrdup(root);
    DIR *d = opendir(dir);
    if (!d) xdie("Could not open directory %s", dir);

    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, "..")) continue;
        char *full = cli_join_path(dir, e->d_name);
        char *child = rel[0] ? cli_join_path(rel, e->d_name) : xstrdup(e->d_name);
        struct stat st;
        if (lstat(full, &st) != 0) xdie("Could not stat %s", full);
        if (S_ISDIR(st.st_mode)) {
            if (!same_real_path(full, skip)) walk(root, child, skip, out);
        } else if (S_ISREG(st.st_mode) && is_fasta(e->d_name)) {
            path_push(out, child);
        }
        free(child);
        free(full);
    }
    closedir(d);
    free(dir);
}

static char *output_name(const char *rel, bool decompress) {
    const char *slash = strrchr(rel, '/');
    const char *base = slash ? slash + 1 : rel;
    if (decompress && ends_with_ci(base, ".gz")) {
        size_t n = strlen(base) - 3;
        char *name = xmalloc(n + 1);
        memcpy(name, base, n);
        name[n] = '\0';
        return name;
    }
    return xstrdup(base);
}

static bool name_used(Harvested *items, size_t len, const char *name) {
    for (size_t i = 0; i < len; i++)
        if (!strcmp(items[i].out, name)) return true;
    return false;
}

static char *unique_name(Harvested *items, size_t len, const char *rel,
                         const char *candidate) {
    if (!name_used(items, len, candidate)) return xstrdup(candidate);

    const char *slash = strrchr(rel, '/');
    char *parent = NULL;
    if (slash && slash > rel) {
        const char *start = rel;
        for (const char *q = rel; q < slash; q++)
            if (*q == '/') start = q + 1;
        size_t plen = (size_t)(slash - start);
        parent = xmalloc(plen + 1);
        memcpy(parent, start, plen);
        parent[plen] = '\0';
    }
    char *prefixed = parent ? xmalloc(strlen(parent) + 2 + strlen(candidate) + 1) : NULL;
    if (prefixed) sprintf(prefixed, "%s__%s", parent, candidate);
    free(parent);

    const char *base = prefixed ? prefixed : candidate;
    for (unsigned n = 2;; n++) {
        if (!name_used(items, len, base)) return prefixed ? prefixed : xstrdup(base);
        char *next = xmalloc(strlen(base) + 16);
        sprintf(next, "%s_%u", base, n);
        if (prefixed) free(prefixed);
        prefixed = next;
        base = next;
    }
}

static unsigned long long copy_fasta(const char *src, const char *dst, bool gz, bool force) {
    if (!force) {
        FILE *existing = fopen(dst, "rb");
        if (existing) {
            fclose(existing);
            xdie("refusing to overwrite %s (use --force)", dst);
        }
    }
    FILE *out = fopen(dst, "wb");
    if (!out) xdie("Could not create %s: %s", dst, strerror(errno));
    unsigned long long total = 0;
    char buf[65536];

    if (gz) {
        gzFile in = gzopen(src, "rb");
        if (!in) {
            fclose(out);
            xdie("Could not open %s", src);
        }
        int n;
        while ((n = gzread(in, buf, sizeof(buf))) > 0) {
            if (fwrite(buf, 1, (size_t)n, out) != (size_t)n) {
                gzclose(in);
                fclose(out);
                xdie("Could not write %s", dst);
            }
            total += (unsigned long long)n;
        }
        if (n < 0) {
            int errnum;
            const char *msg = gzerror(in, &errnum);
            gzclose(in);
            fclose(out);
            xdie("Could not read %s: %s", src, msg ? msg : "gzip error");
        }
        gzclose(in);
    } else {
        FILE *in = fopen(src, "rb");
        if (!in) {
            fclose(out);
            xdie("Could not open %s: %s", src, strerror(errno));
        }
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), in)) > 0) {
            if (fwrite(buf, 1, n, out) != n) {
                fclose(in);
                fclose(out);
                xdie("Could not write %s", dst);
            }
            total += n;
        }
        if (ferror(in)) {
            fclose(in);
            fclose(out);
            xdie("Could not read %s", src);
        }
        fclose(in);
    }
    if (fclose(out) != 0) xdie("Could not close %s", dst);
    return total;
}

int cmd_harvest(int argc, char **argv) {
    const char *root = NULL, *outdir = NULL;
    bool decompress = false, force = false;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-o") || !strcmp(a, "--output")) outdir = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) outdir = a + 9;
        else if (!strcmp(a, "--decompress")) decompress = true;
        else if (!strcmp(a, "--force")) force = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("harvest"); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else if (root) xdie("unexpected argument: %s", a);
        else root = a;
    }
    if (!root) xdie("harvest: missing dataset directory");

    char *outdir_owned = NULL;
    if (!outdir) {
        outdir_owned = cli_join_path(root, "fasta");
        outdir = outdir_owned;
    }
    if (mkdir_p(outdir) != 0) xdie("Could not create %s", outdir);
    char *skip = realpath(outdir, NULL);

    PathList paths = {0};
    walk(root, "", skip ? skip : outdir, &paths);
    qsort(paths.data, paths.len, sizeof(*paths.data), cmp_strptr);

    Harvested *items = xmalloc((paths.len ? paths.len : 1) * sizeof(*items));
    size_t nitems = 0;
    for (size_t i = 0; i < paths.len; i++) {
        const char *rel = paths.data[i];
        char *candidate = output_name(rel, decompress);
        char *name = unique_name(items, nitems, rel, candidate);
        free(candidate);
        char *dst = cli_join_path(outdir, name);
        char *src = cli_join_path(root, rel);
        bool gz = ends_with_ci(rel, ".gz");
        unsigned long long bytes = copy_fasta(src, dst, gz && decompress, force);
        free(src);
        items[nitems].out = name;
        items[nitems].rel = xstrdup(rel);
        items[nitems].bytes = bytes;
        nitems++;
        free(dst);
    }

    for (size_t i = 0; i < nitems; i++)
        printf("%s\t%s\t%llu\n", items[i].out, items[i].rel, items[i].bytes);
    fprintf(stderr, "harvest: %zu FASTA files -> %s\n", nitems, outdir);
    if (nitems == 0) {
        fprintf(stderr, "harvest: no FASTA files found under %s\n", root);
        free(skip);
        free(outdir_owned);
        return 1;
    }

    for (size_t i = 0; i < nitems; i++) {
        free(items[i].out);
        free(items[i].rel);
    }
    free(items);
    for (size_t i = 0; i < paths.len; i++) free(paths.data[i]);
    free(paths.data);
    free(skip);
    free(outdir_owned);
    return 0;
}
