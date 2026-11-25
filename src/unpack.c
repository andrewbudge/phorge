#define _POSIX_C_SOURCE 200809L

#include "cli_help.h"
#include "cli_util.h"

#include <ctype.h>
#include <errno.h>
#include <minizip/unzip.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

typedef struct {
    char *name;
    unsigned long long uncompressed;
    unsigned long long compressed;
} Entry;

typedef struct {
    Entry *data;
    size_t len;
    size_t cap;
} Entries;

static void entries_push(Entries *list, char *name, unsigned long long usize,
                         unsigned long long csize) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        Entry *p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len].name = name;
    list->data[list->len].uncompressed = usize;
    list->data[list->len].compressed = csize;
    list->len++;
}

static int cmp_entry(const void *a, const void *b) {
    return strcmp(((const Entry *)a)->name, ((const Entry *)b)->name);
}

/* Reject absolute paths, drive letters, backslashes and '..' components.
   Returns a normalized relative path or NULL. */
static char *safe_relpath(const char *raw) {
    if (!raw || !*raw) return NULL;
    if (raw[0] == '/' || raw[0] == '\\') return NULL;
    if (isalpha((unsigned char)raw[0]) && raw[1] == ':') return NULL;
    if (strchr(raw, '\\')) return NULL;

    char *out = xmalloc(strlen(raw) + 1);
    size_t o = 0;
    const char *p = raw;
    while (*p) {
        const char *slash = strchr(p, '/');
        size_t n = slash ? (size_t)(slash - p) : strlen(p);
        if (n == 0 || (n == 1 && p[0] == '.')) {
            /* skip */
        } else if (n == 2 && p[0] == '.' && p[1] == '.') {
            free(out);
            return NULL;
        } else {
            memcpy(out + o, p, n);
            o += n;
            out[o++] = '/';
        }
        if (!slash) break;
        p = slash + 1;
    }
    if (o == 0) {
        free(out);
        return NULL;
    }
    out[o - 1] = '\0';
    return out;
}

static void make_parents(const char *dest, const char *rel) {
    char *path = xstrdup(dest);
    const char *p = rel;
    for (;;) {
        const char *slash = strchr(p, '/');
        if (!slash) break;
        size_t n = (size_t)(slash - p);
        char *part = xmalloc(n + 1);
        memcpy(part, p, n);
        part[n] = '\0';
        char *next = cli_join_path(path, part);
        free(part);
        if (mkdir(next, 0777) != 0 && errno != EEXIST)
            xdie("mkdir %s: %s", next, strerror(errno));
        free(path);
        path = next;
        p = slash + 1;
    }
    free(path);
}

static void extract_entry(unzFile zip, const char *dest, const char *rel, bool force) {
    make_parents(dest, rel);
    char *out = cli_join_path(dest, rel);

    if (!force) {
        FILE *existing = fopen(out, "rb");
        if (existing) {
            fclose(existing);
            xdie("refusing to overwrite %s (use --force)", out);
        }
    }

    FILE *f = fopen(out, "wb");
    if (!f) xdie("Could not create %s: %s", out, strerror(errno));
    if (unzOpenCurrentFile(zip) != UNZ_OK) {
        fclose(f);
        unlink(out);
        xdie("Could not open zip entry '%s' (encrypted?)", rel);
    }

    char buf[65536];
    int n;
    while ((n = unzReadCurrentFile(zip, buf, sizeof(buf))) > 0) {
        if (fwrite(buf, 1, (size_t)n, f) != (size_t)n) {
            unzCloseCurrentFile(zip);
            fclose(f);
            xdie("Could not write %s", out);
        }
    }
    if (n < 0) {
        unzCloseCurrentFile(zip);
        fclose(f);
        xdie("Could not read zip entry '%s'", rel);
    }
    if (unzCloseCurrentFile(zip) == UNZ_CRCERROR) {
        fclose(f);
        xdie("CRC error in zip entry '%s'", rel);
    }
    if (fclose(f) != 0) xdie("Could not close %s", out);
    free(out);
}

static char *default_dest(const char *zip_path) {
    const char *slash = strrchr(zip_path, '/');
    const char *base = slash ? slash + 1 : zip_path;
    char *stem, *ext;
    stem_and_ext(base, &stem, &ext);
    free(ext);
    char *out;
    if (slash) {
        size_t n = (size_t)(slash - zip_path);
        out = xmalloc(n + 1 + strlen(stem) + 1);
        memcpy(out, zip_path, n);
        out[n] = '/';
        strcpy(out + n + 1, stem);
    } else {
        out = xstrdup(stem);
    }
    free(stem);
    return out;
}

int cmd_unpack(int argc, char **argv) {
    const char *zip_path = NULL, *dest = NULL;
    bool list = false, force = false, prune = false;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-o") || !strcmp(a, "--output")) dest = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) dest = a + 9;
        else if (!strcmp(a, "--list")) list = true;
        else if (!strcmp(a, "--force")) force = true;
        else if (!strcmp(a, "--prune-zip")) prune = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("unpack"); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else if (zip_path) xdie("unexpected argument: %s", a);
        else zip_path = a;
    }
    if (!zip_path) xdie("unpack: missing ZIP file");

    char *dest_owned = NULL;
    if (!list) {
        dest_owned = dest ? xstrdup(dest) : default_dest(zip_path);
        dest = dest_owned;
        if (mkdir_p(dest) != 0) xdie("Could not create %s", dest);
    }

    unzFile zip = unzOpen(zip_path);
    if (!zip) xdie("Could not open zip '%s'", zip_path);

    Entries entries = {0};
    int ret = unzGoToFirstFile(zip);
    while (ret == UNZ_OK) {
        unz_file_info64 info;
        if (unzGetCurrentFileInfo64(zip, &info, NULL, 0, NULL, 0, NULL, 0) != UNZ_OK)
            xdie("Could not read zip directory");
        char *raw = xmalloc((size_t)info.size_filename + 1);
        if (unzGetCurrentFileInfo64(zip, &info, raw, (uLong)info.size_filename + 1,
                                    NULL, 0, NULL, 0) != UNZ_OK)
            xdie("Could not read zip entry name");
        raw[info.size_filename] = '\0';

        char *rel = safe_relpath(raw);
        if (!rel) {
            fprintf(stderr, "unpack: skipping unsafe path '%s'\n", raw);
            free(raw);
            ret = unzGoToNextFile(zip);
            continue;
        }
        entries_push(&entries, rel, info.uncompressed_size, info.compressed_size);

        if (!list) {
            size_t rawlen = strlen(raw);
            if (rawlen == 0 || raw[rawlen - 1] != '/') extract_entry(zip, dest, rel, force);
        }
        free(raw);
        ret = unzGoToNextFile(zip);
    }
    if (ret != UNZ_OK && ret != UNZ_END_OF_LIST_OF_FILE) xdie("error reading zip directory");
    unzClose(zip);

    qsort(entries.data, entries.len, sizeof(*entries.data), cmp_entry);
    for (size_t i = 0; i < entries.len; i++)
        printf("%s\t%llu\t%llu\n", entries.data[i].name,
               entries.data[i].uncompressed, entries.data[i].compressed);

    unsigned long long total = 0;
    for (size_t i = 0; i < entries.len; i++) total += entries.data[i].uncompressed;
    fprintf(stderr, "unpack: %zu entries, %llu bytes %s %s\n", entries.len, total,
            list ? "in" : "->", list ? zip_path : dest);

    if (prune && !list) {
        if (unlink(zip_path) != 0) xdie("Could not remove %s: %s", zip_path, strerror(errno));
        fprintf(stderr, "unpack: removed %s\n", zip_path);
    }

    for (size_t i = 0; i < entries.len; i++) free(entries.data[i].name);
    free(entries.data);
    free(dest_owned);
    return 0;
}
