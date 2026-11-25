#define _POSIX_C_SOURCE 200809L

#include "dataset.h"
#include "cli_util.h"

#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <minizip/unzip.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>

void ds_kv_free(DsKVList *list) {
    for (size_t i = 0; i < list->len; i++) {
        free(list->data[i].key);
        free(list->data[i].val);
    }
    free(list->data);
    list->data = NULL;
    list->len = list->cap = 0;
}

const char *ds_kv_get(const DsKVList *list, const char *key) {
    for (size_t i = 0; i < list->len; i++)
        if (!strcmp(list->data[i].key, key)) return list->data[i].val;
    return NULL;
}

static void kv_push(DsKVList *list, char *key, char *val) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        DsKV *p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len].key = key;
    list->data[list->len].val = val;
    list->len++;
}

static void skip_ws(const char **pp) {
    while (**pp && isspace((unsigned char)**pp)) (*pp)++;
}

static char *json_string(const char **pp) {
    const char *p = *pp;
    if (*p != '"') return NULL;
    p++;
    char *out = xmalloc(strlen(p) + 1);
    size_t n = 0;
    while (*p && *p != '"') {
        if (*p == '\\') {
            p++;
            switch (*p) {
                case 'n': out[n++] = '\n'; break;
                case 'r': out[n++] = '\r'; break;
                case 't': out[n++] = '\t'; break;
                case 'b': out[n++] = '\b'; break;
                case 'f': out[n++] = '\f'; break;
                case 'u':
                    for (int i = 0; i < 4 && p[1]; i++) p++;
                    out[n++] = '?';
                    break;
                default: out[n++] = *p; break;
            }
        } else {
            out[n++] = *p;
        }
        p++;
    }
    if (*p == '"') p++;
    out[n] = '\0';
    *pp = p;
    return out;
}

/* Skip a JSON value, including nested objects/arrays, and return its raw text. */
static char *skip_value_raw(const char **pp) {
    skip_ws(pp);
    const char *start = *pp;
    const char *p = *pp;
    if (*p == '"') {
        free(json_string(pp));
        return NULL;
    }
    if (*p == '{' || *p == '[') {
        char open = *p, close = open == '{' ? '}' : ']';
        int depth = 0;
        bool in_string = false, escaped = false;
        for (; *p; p++) {
            if (in_string) {
                if (escaped) escaped = false;
                else if (*p == '\\') escaped = true;
                else if (*p == '"') in_string = false;
                continue;
            }
            if (*p == '"') in_string = true;
            else if (*p == open) depth++;
            else if (*p == close && --depth == 0) { p++; break; }
        }
        *pp = p;
        return strndup(start, (size_t)(p - start));
    }
    while (*p && *p != ',' && *p != '}' && *p != ']') p++;
    *pp = p;
    return strndup(start, (size_t)(p - start));
}

static char *trim_copy(const char *s, size_t n) {
    while (n && isspace((unsigned char)*s)) { s++; n--; }
    while (n && isspace((unsigned char)s[n - 1])) n--;
    char *out = xmalloc(n + 1);
    memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

static char *join_dotted(const char *prefix, const char *key) {
    if (!prefix) return xstrdup(key);
    size_t n = strlen(prefix), m = strlen(key);
    char *out = xmalloc(n + 1 + m + 1);
    memcpy(out, prefix, n);
    out[n] = '.';
    memcpy(out + n + 1, key, m);
    out[n + 1 + m] = '\0';
    return out;
}

static void collect_scalars(const char **pp, const char *prefix, DsKVList *out) {
    skip_ws(pp);
    if (**pp != '{') { free(skip_value_raw(pp)); return; }
    (*pp)++;
    for (;;) {
        skip_ws(pp);
        if (**pp == ',') { (*pp)++; continue; }
        if (**pp == '}' || !**pp) { if (**pp == '}') (*pp)++; break; }
        if (**pp != '"') { free(skip_value_raw(pp)); continue; }
        char *key = json_string(pp);
        skip_ws(pp);
        if (**pp == ':') (*pp)++;
        skip_ws(pp);
        char *full = join_dotted(prefix, key);
        free(key);
        if (**pp == '{') {
            collect_scalars(pp, full, out);
        } else if (**pp == '"') {
            kv_push(out, full, json_string(pp));
            full = NULL;
        } else {
            char *raw = skip_value_raw(pp);
            kv_push(out, full, raw);
            full = NULL;
        }
        free(full);
    }
}

void ds_parse_json_scalars(const char *json, DsKVList *out) {
    const char *p = json;
    collect_scalars(&p, NULL, out);
}

bool ds_is_zip(const char *source) {
    struct stat st;
    if (stat(source, &st) != 0) return false;
    return !S_ISDIR(st.st_mode);
}

char *ds_read_text(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    size_t cap = 4096, len = 0;
    char *buf = xmalloc(cap);
    size_t n;
    while ((n = fread(buf + len, 1, cap - len - 1, f)) > 0) {
        len += n;
        if (len + 1 >= cap) {
            cap *= 2;
            char *p = realloc(buf, cap);
            if (!p) xdie("out of memory");
            buf = p;
        }
    }
    fclose(f);
    buf[len] = '\0';
    return buf;
}

static char *read_zip_entry(unzFile zip) {
    if (unzOpenCurrentFile(zip) != UNZ_OK) return NULL;
    size_t cap = 4096, len = 0;
    char *buf = xmalloc(cap);
    for (;;) {
        if (len + 1 >= cap) {
            cap *= 2;
            char *p = realloc(buf, cap);
            if (!p) xdie("out of memory");
            buf = p;
        }
        int n = unzReadCurrentFile(zip, buf + len, (unsigned)(cap - len - 1));
        if (n > 0) len += (size_t)n;
        else if (n == 0) break;
        else {
            unzCloseCurrentFile(zip);
            free(buf);
            return NULL;
        }
    }
    unzCloseCurrentFile(zip);
    buf[len] = '\0';
    return buf;
}

static void meta_push(DsMetaList *list, DsMeta m) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        DsMeta *p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len++] = m;
}

static const char *nonempty(const char *s) { return s && *s ? s : NULL; }

static void push_meta_row(DsMetaList *out, const char *key, const DsKVList *kvs,
                          const char *group_id, const char *group_method) {
    DsMeta m = {0};
    m.gene_id = xstrdup(key ? key : "");
    const char *v;
    if ((v = nonempty(ds_kv_get(kvs, "accession"))) || (v = nonempty(ds_kv_get(kvs, "proteinAccession"))))
        m.accession = xstrdup(v);
    if (nonempty(group_id)) m.group_id = xstrdup(group_id);
    if (nonempty(group_method)) m.group_method = xstrdup(group_method);
    if ((v = nonempty(ds_kv_get(kvs, "taxId"))) || (v = nonempty(ds_kv_get(kvs, "organism.taxId"))))
        m.taxid = xstrdup(v);
    if ((v = nonempty(ds_kv_get(kvs, "taxname"))) || (v = nonempty(ds_kv_get(kvs, "organism.organismName"))))
        m.taxname = xstrdup(v);
    if ((v = nonempty(ds_kv_get(kvs, "symbol"))))
        m.symbol = xstrdup(v);
    if ((v = nonempty(ds_kv_get(kvs, "description"))) || (v = nonempty(ds_kv_get(kvs, "proteinName"))))
        m.description = xstrdup(v);
    meta_push(out, m);
}

static void parse_group_array(const char *raw, DsMetaList *out, const char *key,
                              const DsKVList *kvs) {
    const char *p = strchr(raw, '[');
    if (!p) { push_meta_row(out, key, kvs, NULL, NULL); return; }
    p++;
    size_t ngroups = 0;
    for (;;) {
        while (*p && *p != '{' && *p != ']') p++;
        if (*p == ']' || !*p) break;
        const char *start = p;
        int depth = 0;
        bool in_string = false, escaped = false;
        for (; *p; p++) {
            if (in_string) {
                if (escaped) escaped = false;
                else if (*p == '\\') escaped = true;
                else if (*p == '"') in_string = false;
                continue;
            }
            if (*p == '"') in_string = true;
            else if (*p == '{') depth++;
            else if (*p == '}' && --depth == 0) { p++; break; }
        }
        char *obj = strndup(start, (size_t)(p - start));
        DsKVList okvs = {0};
        ds_parse_json_scalars(obj, &okvs);
        const char *id = ds_kv_get(&okvs, "id");
        const char *method = ds_kv_get(&okvs, "method");
        push_meta_row(out, key, kvs, id, method);
        ngroups++;
        ds_kv_free(&okvs);
        free(obj);
    }
    if (ngroups == 0) push_meta_row(out, key, kvs, NULL, NULL);
}

static void parse_report_line(const char *line, DsMetaList *out) {
    DsKVList kvs = {0};
    ds_parse_json_scalars(line, &kvs);
    const char *gene = nonempty(ds_kv_get(&kvs, "geneId"));
    const char *acc = nonempty(ds_kv_get(&kvs, "accession"));
    const char *prot = nonempty(ds_kv_get(&kvs, "proteinAccession"));
    if (gene || acc || prot) {
        const char *key = gene ? gene : acc ? acc : prot;
        const char *groups = ds_kv_get(&kvs, "geneGroups");
        if (groups) parse_group_array(groups, out, key, &kvs);
        else push_meta_row(out, key, &kvs, NULL, NULL);
    }
    ds_kv_free(&kvs);
}

static void parse_report_text(const char *text, DsMetaList *out) {
    const char *p = text;
    while (*p) {
        const char *line = p;
        while (*p && *p != '\n') p++;
        if (*p == '\n') { ((char *)p)[0] = '\0'; p++; }
        if (*line) parse_report_line(line, out);
    }
}

static bool report_entry(const char *name) {
    const char *base = strrchr(name, '/');
    base = base ? base + 1 : name;
    return !strcmp(base, "data_report.jsonl") || !strcmp(base, "annotation_report.jsonl");
}

void ds_load_reports(const char *source, DsMetaList *out) {
    static const char *const names[] = {"data_report.jsonl", "annotation_report.jsonl"};
    if (!ds_is_zip(source)) {
        char *ncbi = cli_join_path(source, "ncbi_dataset");
        char *data = cli_join_path(ncbi, "data");
        for (size_t j = 0; j < sizeof(names) / sizeof(*names); j++) {
            char *path = cli_join_path(data, names[j]);
            char *text = ds_read_text(path);
            if (text) {
                parse_report_text(text, out);
                free(text);
            }
            free(path);
        }
        free(ncbi);
        free(data);
        return;
    }
    unzFile zip = unzOpen(source);
    if (!zip) xdie("could not open zip '%s'", source);
    int ret = unzGoToFirstFile(zip);
    while (ret == UNZ_OK) {
        unz_file_info64 info;
        unzGetCurrentFileInfo64(zip, &info, NULL, 0, NULL, 0, NULL, 0);
        char *name = xmalloc((size_t)info.size_filename + 1);
        unzGetCurrentFileInfo64(zip, &info, name, (uLong)info.size_filename + 1,
                                NULL, 0, NULL, 0);
        name[info.size_filename] = '\0';
        if (report_entry(name)) {
            char *text = read_zip_entry(zip);
            if (text) {
                parse_report_text(text, out);
                free(text);
            }
        }
        free(name);
        ret = unzGoToNextFile(zip);
    }
    unzClose(zip);
}

const char *ds_meta_lookup(const DsMetaList *list, const char *key, int field) {
    for (size_t i = 0; i < list->len; i++) {
        if (strcmp(list->data[i].gene_id, key) && strcmp(list->data[i].accession ? list->data[i].accession : "", key))
            continue;
        const char *v = field == DS_META_TAXID ? list->data[i].taxid
                      : field == DS_META_TAXNAME ? list->data[i].taxname
                      : field == DS_META_DESC ? list->data[i].description
                                              : list->data[i].symbol;
        if (nonempty(v)) return v;
    }
    return NULL;
}

const DsMeta *ds_meta_find(const DsMetaList *list, const char *key) {
    for (size_t i = 0; i < list->len; i++) {
        if (!strcmp(list->data[i].gene_id, key)) return &list->data[i];
        if (list->data[i].accession && !strcmp(list->data[i].accession, key))
            return &list->data[i];
    }
    return NULL;
}

void ds_meta_free(DsMetaList *list) {
    for (size_t i = 0; i < list->len; i++) {
        free(list->data[i].gene_id); free(list->data[i].accession);
        free(list->data[i].group_id); free(list->data[i].group_method);
        free(list->data[i].taxid); free(list->data[i].taxname);
        free(list->data[i].symbol); free(list->data[i].description);
    }
    free(list->data);
    list->data = NULL;
    list->len = list->cap = 0;
}

char *ds_header_tag(const char *header, const char *key) {
    char needle[128];
    snprintf(needle, sizeof(needle), "[%s=", key);
    const char *p = strstr(header, needle);
    if (!p) return NULL;
    p += strlen(needle);
    const char *end = strchr(p, ']');
    if (!end) return NULL;
    return trim_copy(p, (size_t)(end - p));
}

static char *trim_in_place(char *s) {
    while (*s && isspace((unsigned char)*s)) s++;
    char *end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) *--end = '\0';
    return s;
}

static void fasta_push(DsFastaList *list, const char *header, const char *seq,
                       const char *source) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        DsFastaRecord *p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len].header = xstrdup(header);
    list->data[list->len].seq = xstrdup(seq);
    list->data[list->len].source = xstrdup(source);
    list->len++;
}

static void parse_fasta_text(const char *text, const char *source, DsFastaList *out) {
    const char *p = text;
    char *header = NULL;
    char *seq = NULL;
    size_t seq_len = 0, seq_cap = 0;
    while (*p) {
        const char *line = p;
        while (*p && *p != '\n') p++;
        if (*p == '\n') { ((char *)p)[0] = '\0'; p++; }
        char *trimmed = trim_in_place((char *)line);
        if (trimmed[0] == '>') {
            if (header) fasta_push(out, header, seq ? seq : "", source);
            free(header);
            free(seq);
            header = trimmed[1] ? xstrdup(trimmed + 1) : NULL;
            seq = NULL;
            seq_len = seq_cap = 0;
        } else if (trimmed[0] && header) {
            size_t n = strlen(trimmed);
            if (seq_len + n + 1 > seq_cap) {
                seq_cap = (seq_len + n + 1) * 2;
                char *q = realloc(seq, seq_cap);
                if (!q) xdie("out of memory");
                seq = q;
            }
            for (size_t i = 0; i < n; i++)
                seq[seq_len++] = (trimmed[i] >= 'a' && trimmed[i] <= 'z')
                               ? (char)(trimmed[i] - 'a' + 'A') : trimmed[i];
            seq[seq_len] = '\0';
        }
    }
    if (header) fasta_push(out, header, seq ? seq : "", source);
    free(header);
    free(seq);
}

static bool is_fasta_name(const char *name) {
    static const char *const exts[] = {".fna", ".faa", ".fasta", ".fa", ".ffn", ".fas"};
    for (size_t i = 0; i < sizeof(exts) / sizeof(*exts); i++) {
        size_t n = strlen(name), m = strlen(exts[i]);
        if (m <= n && !strcasecmp(name + n - m, exts[i])) return true;
    }
    return false;
}

static void collect_fasta_files(const char *dir, DsFastaList *unused,
                                void (*cb)(const char *, const char *, void *), void *ctx) {
    (void)unused;
    struct dirent *e;
    DIR *d = opendir(dir);
    if (!d) return;
    while ((e = readdir(d)) != NULL) {
        if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, "..")) continue;
        char *full = cli_join_path(dir, e->d_name);
        struct stat st;
        if (lstat(full, &st) == 0) {
            if (S_ISDIR(st.st_mode)) collect_fasta_files(full, NULL, cb, ctx);
            else if (S_ISREG(st.st_mode) && is_fasta_name(e->d_name)) cb(full, e->d_name, ctx);
        }
        free(full);
    }
    closedir(d);
}

typedef struct {
    const char *which;
    DsFastaList *out;
} FastaCtx;

static void load_fasta_file(const char *path, const char *base, void *ctxp) {
    FastaCtx *ctx = ctxp;
    if (ctx->which && strcmp(ctx->which, base)) return;
    char *text = ds_read_text(path);
    if (!text) xdie("could not read %s", path);
    parse_fasta_text(text, base, ctx->out);
    free(text);
}

void ds_load_fasta(const char *source, const char *which, DsFastaList *out) {
    if (!ds_is_zip(source)) {
        char *ncbi = cli_join_path(source, "ncbi_dataset");
        char *data = cli_join_path(ncbi, "data");
        FastaCtx ctx = {which, out};
        collect_fasta_files(data, NULL, load_fasta_file, &ctx);
        free(data);
        free(ncbi);
        return;
    }
    unzFile zip = unzOpen(source);
    if (!zip) xdie("could not open zip '%s'", source);
    int ret = unzGoToFirstFile(zip);
    while (ret == UNZ_OK) {
        unz_file_info64 info;
        unzGetCurrentFileInfo64(zip, &info, NULL, 0, NULL, 0, NULL, 0);
        char *name = xmalloc((size_t)info.size_filename + 1);
        unzGetCurrentFileInfo64(zip, &info, name, (uLong)info.size_filename + 1,
                                NULL, 0, NULL, 0);
        name[info.size_filename] = '\0';
        const char *base = strrchr(name, '/');
        base = base ? base + 1 : name;
        if (is_fasta_name(base) && (!which || !strcmp(which, base))) {
            char *text = read_zip_entry(zip);
            if (text) {
                parse_fasta_text(text, base, out);
                free(text);
            }
        }
        free(name);
        ret = unzGoToNextFile(zip);
    }
    unzClose(zip);
}

char *ds_sanitize_name(const char *s) {
    if (!s || !*s) return xstrdup("NA");
    char *out = xmalloc(strlen(s) + 1);
    size_t n = 0;
    for (const char *p = s; *p; p++) {
        unsigned char c = (unsigned char)*p;
        if (isalnum(c) || c == '-' || c == '_' || c == '.') out[n++] = (char)c;
        else out[n++] = '_';
    }
    out[n] = '\0';
    return out;
}

void ds_fasta_free(DsFastaList *list) {
    for (size_t i = 0; i < list->len; i++) {
        free(list->data[i].header);
        free(list->data[i].seq);
        free(list->data[i].source);
    }
    free(list->data);
    list->data = NULL;
    list->len = list->cap = 0;
}
