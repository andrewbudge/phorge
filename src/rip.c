#define _POSIX_C_SOURCE 200809L

#include "cli_help.h"
#include "cli_util.h"
#include "phorge_c.h"

#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct {
    char **data;
    size_t len;
    size_t cap;
} StrList;

static void str_push(StrList *list, const char *s) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 8;
        char **p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len++] = xstrdup(s);
}

static void str_free(StrList *list) {
    for (size_t i = 0; i < list->len; i++) free(list->data[i]);
    free(list->data);
    list->data = NULL;
    list->len = list->cap = 0;
}

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n);
    if (!q) xdie("out of memory");
    return q;
}

static char *str_ndup(const char *s, size_t n) {
    char *p = xmalloc(n + 1);
    memcpy(p, s, n);
    p[n] = '\0';
    return p;
}

static char *trim_ws_copy(const char *s) {
    while (*s && isspace((unsigned char)*s)) s++;
    const char *end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) end--;
    return str_ndup(s, (size_t)(end - s));
}

static char *trim_quotes(char *s) {
    size_t n = strlen(s);
    size_t start = 0, end = n;
    while (start < end && s[start] == '"') start++;
    while (end > start && s[end - 1] == '"') end--;
    if (start) memmove(s, s + start, end - start);
    s[end - start] = '\0';
    return s;
}

typedef struct {
    size_t start;
    size_t end;
} Loc;

typedef struct {
    char *type;
    Loc *locs;
    size_t nlocs;
    bool complement;
    char *gene;
    char *product;
} Feature;

typedef struct {
    char *accession;
    char *organism;
    char *seq;
    Feature *features;
    size_t nfeatures;
    size_t cap_features;
} GBRecord;

static void feature_push(GBRecord *rec, Feature f) {
    if (rec->nfeatures == rec->cap_features) {
        size_t cap = rec->cap_features ? rec->cap_features * 2 : 16;
        Feature *p = realloc(rec->features, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        rec->features = p;
        rec->cap_features = cap;
    }
    rec->features[rec->nfeatures++] = f;
}

static void feature_free(Feature *f) {
    free(f->type);
    free(f->locs);
    free(f->gene);
    free(f->product);
}

static void record_free(GBRecord *rec) {
    free(rec->accession);
    free(rec->organism);
    free(rec->seq);
    for (size_t i = 0; i < rec->nfeatures; i++) feature_free(&rec->features[i]);
    free(rec->features);
    memset(rec, 0, sizeof(*rec));
}

static bool valid_feature_type(const char *s) {
    static const char *const types[] = {"rRNA", "CDS", "gene", "tRNA", "misc_RNA"};
    for (size_t i = 0; i < sizeof(types) / sizeof(*types); i++)
        if (!strcmp(s, types[i])) return true;
    return false;
}

static void parse_ranges(char *s, Feature *f) {
    char *token = strtok(s, ",");
    while (token) {
        while (*token && isspace((unsigned char)*token)) token++;
        size_t n = strlen(token);
        while (n && isspace((unsigned char)token[n - 1])) n--;
        char *clean = str_ndup(token, n);
        for (char *p = clean; *p; p++)
            if (*p == '<' || *p == '>') *p = ' ';
        char *a = clean;
        while (*a == ' ') a++;
        char *dots = strstr(a, "..");
        size_t start = 0, end = 0;
        if (dots) {
            *dots = '\0';
            char *b = dots + 2;
            start = strtoull(a, NULL, 10);
            end = strtoull(b, NULL, 10);
        } else {
            start = end = strtoull(a, NULL, 10);
        }
        if (start || end) {
            f->locs = xrealloc(f->locs, (f->nlocs + 1) * sizeof(*f->locs));
            f->locs[f->nlocs].start = start;
            f->locs[f->nlocs].end = end;
            f->nlocs++;
        }
        free(clean);
        token = strtok(NULL, ",");
    }
}

static void parse_location(const char *loc, Feature *f) {
    char *s = xstrdup(loc);
    if (!strncmp(s, "complement(", 11) && s[strlen(s) - 1] == ')') {
        f->complement = true;
        size_t n = strlen(s);
        memmove(s, s + 11, n - 12);
        s[n - 12] = '\0';
    }
    if (!strncmp(s, "join(", 5) && s[strlen(s) - 1] == ')') {
        size_t n = strlen(s);
        memmove(s, s + 5, n - 6);
        s[n - 6] = '\0';
    }
    parse_ranges(s, f);
    free(s);
}

static void parse_feature_block(const StrList *lines, GBRecord *rec) {
    if (lines->len == 0) return;
    char **fields = NULL;
    size_t nf = 0;
    char *first = xstrdup(lines->data[0]);
    char *tok = strtok(first, " \t");
    while (tok) {
        fields = xrealloc(fields, (nf + 1) * sizeof(*fields));
        fields[nf++] = tok;
        tok = strtok(NULL, " \t");
    }
    if (nf < 2 || !valid_feature_type(fields[0])) {
        free(fields);
        free(first);
        return;
    }

    Feature f = {0};
    f.type = xstrdup(fields[0]);
    parse_location(fields[1], &f);
    free(fields);
    free(first);

    char *key = NULL, *val = NULL;
    for (size_t i = 1; i < lines->len; i++) {
        char *line = xstrdup(lines->data[i]);
        char *trimmed = line;
        while (*trimmed && isspace((unsigned char)*trimmed)) trimmed++;
        if (*trimmed == '/') {
            if (key) {
                trim_quotes(val);
                if (!strcmp(key, "gene")) { free(f.gene); f.gene = xstrdup(val); }
                else if (!strcmp(key, "product")) { free(f.product); f.product = xstrdup(val); }
                free(key);
                free(val);
                key = val = NULL;
            }
            trimmed++;
            char *eq = strchr(trimmed, '=');
            if (eq) {
                *eq = '\0';
                key = xstrdup(trimmed);
                const char *raw = eq + 1;
                while (*raw && isspace((unsigned char)*raw)) raw++;
                if (*raw == '"') raw++;
                val = xstrdup(raw);
            } else {
                key = xstrdup(trimmed);
                val = xstrdup("");
            }
        } else if (key) {
            char *end = val + strlen(val);
            size_t need = strlen(trimmed) + 1;
            size_t used = (size_t)(end - val);
            val = xrealloc(val, used + need);
            memcpy(val + used, trimmed, need);
        }
        free(line);
    }
    if (key) {
        trim_quotes(val);
        if (!strcmp(key, "gene")) f.gene = xstrdup(val);
        else if (!strcmp(key, "product")) f.product = xstrdup(val);
        free(key);
        free(val);
    }
    feature_push(rec, f);
}

static void parse_genbank(const char *path, GBRecord **out, size_t *nout) {
    FILE *fp = fopen(path, "rb");
    if (!fp) xdie("could not open %s: %s", path, strerror(errno));

    GBRecord *records = NULL;
    size_t nrecords = 0;
    GBRecord rec = {0};
    size_t seq_len = 0, seq_cap = 0;
    bool in_features = false, in_origin = false;
    StrList feature_lines = {0};
    char *line = NULL;
    size_t linecap = 0;
    ssize_t nread;

    while ((nread = getline(&line, &linecap, fp)) != -1) {
        while (nread > 0 && (line[nread - 1] == '\n' || line[nread - 1] == '\r'))
            line[--nread] = '\0';

        if (!strncmp(line, "//", 2)) {
            if (feature_lines.len) parse_feature_block(&feature_lines, &rec);
            str_free(&feature_lines);
            if (rec.seq || rec.nfeatures) {
                if (!rec.accession) {
                    const char *slash = strrchr(path, '/');
                    const char *base = slash ? slash + 1 : path;
                    char *dot = strrchr(base, '.');
                    rec.accession = dot ? str_ndup(base, (size_t)(dot - base)) : xstrdup(base);
                }
                records = xrealloc(records, (nrecords + 1) * sizeof(*records));
                records[nrecords++] = rec;
                memset(&rec, 0, sizeof(rec));
                seq_len = seq_cap = 0;
            }
            in_features = in_origin = false;
            continue;
        }
        if (!strncmp(line, "VERSION", 7)) {
            char *p = line + 7;
            while (*p && isspace((unsigned char)*p)) p++;
            char *end = p;
            while (*end && !isspace((unsigned char)*end)) end++;
            if (end > p) {
                free(rec.accession);
                rec.accession = str_ndup(p, (size_t)(end - p));
            }
            continue;
        }
        if (!rec.organism && strstr(line, "ORGANISM")) {
            const char *p = strstr(line, "ORGANISM") + 8;
            while (*p && isspace((unsigned char)*p)) p++;
            rec.organism = xstrdup(p);
            continue;
        }
        if (!strncmp(line, "FEATURES", 8)) { in_features = true; in_origin = false; continue; }
        if (!strncmp(line, "ORIGIN", 6)) {
            if (feature_lines.len) parse_feature_block(&feature_lines, &rec);
            str_free(&feature_lines);
            in_features = false;
            in_origin = true;
            continue;
        }
        if (in_features) {
            if ((int)strlen(line) > 5 && line[5] != ' ') {
                if (feature_lines.len) parse_feature_block(&feature_lines, &rec);
                str_free(&feature_lines);
                str_push(&feature_lines, line);
            } else if (feature_lines.len) {
                str_push(&feature_lines, line);
            }
        }
        if (in_origin) {
            size_t need = seq_len + (size_t)nread + 1;
            if (need > seq_cap) {
                seq_cap = seq_cap ? seq_cap : 4096;
                while (seq_cap < need) seq_cap *= 2;
                rec.seq = xrealloc(rec.seq, seq_cap);
            }
            for (ssize_t i = 0; i < nread; i++)
                if (isalpha((unsigned char)line[i]))
                    rec.seq[seq_len++] = (char)toupper((unsigned char)line[i]);
            rec.seq[seq_len] = '\0';
        }
    }
    if (feature_lines.len) parse_feature_block(&feature_lines, &rec);
    if (rec.seq || rec.nfeatures) {
        if (!rec.accession) rec.accession = xstrdup("");
        records = xrealloc(records, (nrecords + 1) * sizeof(*records));
        records[nrecords++] = rec;
    }
    str_free(&feature_lines);
    free(line);
    fclose(fp);
    *out = records;
    *nout = nrecords;
}

typedef struct {
    char *name;
    char **syns;
    size_t nsyns;
} Target;

typedef struct {
    Target *data;
    size_t len;
    size_t cap;
} Targets;

static const char *const SYN_12S[] = {"12S", "12S rRNA", "s-rRNA", "12S ribosomal RNA", "rrnS",
                                      "rrn12", "small subunit ribosomal RNA",
                                      "12S_ribosomal_RNA", "small_subunit_ribosomal_RNA"};
static const char *const SYN_16S[] = {"16S", "16S rRNA", "l-rRNA", "16S ribosomal RNA", "rrnL",
                                      "rrn16", "large subunit ribosomal RNA",
                                      "16S_ribosomal_RNA", "large_subunit_ribosomal_RNA"};
static const char *const SYN_COI[] = {"COI", "COXI", "COX1", "CO1", "cytochrome oxidase subunit 1",
                                      "cytochrome c oxidase subunit I", "cytochrome oxidase subunit I",
                                      "cytochrome c oxidase subunit 1"};
static const char *const SYN_COII[] = {"COII", "COXII", "COX2", "CO2", "cytochrome oxidase subunit 2",
                                       "cytochrome c oxidase subunit II",
                                       "cytochrome c oxidase subunit 2"};
static const char *const SYN_COIII[] = {"COIII", "COXIII", "COX3", "CO3", "cytochrome oxidase subunit 3",
                                        "cytochrome c oxidase subunit III",
                                        "cytochrome oxidase subunit III",
                                        "cytochrome c oxidase subunit 3"};
static const char *const SYN_CYTB[] = {"CYTB", "Cyt b", "cytochrome b", "cob", "CYB",
                                      "cyt_b", "cytochrome_b"};
static const char *const SYN_ND1[] = {"ND1", "NAD1", "NADH1", "NADH dehydrogenase subunit 1"};
static const char *const SYN_ND2[] = {"ND2", "NAD2", "NADH2", "NADH dehydrogenase subunit 2"};
static const char *const SYN_ND3[] = {"ND3", "NAD3", "NADH3", "NADH dehydrogenase subunit 3"};
static const char *const SYN_ND4[] = {"ND4", "NAD4", "NADH4", "NADH dehydrogenase subunit 4"};
static const char *const SYN_ND4L[] = {"ND4L", "NAD4L", "NADH4L", "NADH dehydrogenase subunit 4L"};
static const char *const SYN_ND5[] = {"ND5", "NAD5", "NADH5", "NADH dehydrogenase subunit 5"};
static const char *const SYN_ND6[] = {"ND6", "NAD6", "NADH6", "NADH dehydrogenase subunit 6"};
static const char *const SYN_ATP6[] = {"ATP6", "ATPase6", "ATPase 6", "ATP synthase F0 subunit 6"};
static const char *const SYN_ATP8[] = {"ATP8", "ATPase8", "ATPase 8", "ATP synthase F0 subunit 8"};

typedef struct { const char *name; const char *const *syns; size_t n; } SynEntry;
static const SynEntry SYNONYMS[] = {
    {"12S", SYN_12S, sizeof(SYN_12S) / sizeof(*SYN_12S)},
    {"16S", SYN_16S, sizeof(SYN_16S) / sizeof(*SYN_16S)},
    {"COI", SYN_COI, sizeof(SYN_COI) / sizeof(*SYN_COI)},
    {"COII", SYN_COII, sizeof(SYN_COII) / sizeof(*SYN_COII)},
    {"COIII", SYN_COIII, sizeof(SYN_COIII) / sizeof(*SYN_COIII)},
    {"CYTB", SYN_CYTB, sizeof(SYN_CYTB) / sizeof(*SYN_CYTB)},
    {"ND1", SYN_ND1, sizeof(SYN_ND1) / sizeof(*SYN_ND1)},
    {"ND2", SYN_ND2, sizeof(SYN_ND2) / sizeof(*SYN_ND2)},
    {"ND3", SYN_ND3, sizeof(SYN_ND3) / sizeof(*SYN_ND3)},
    {"ND4", SYN_ND4, sizeof(SYN_ND4) / sizeof(*SYN_ND4)},
    {"ND4L", SYN_ND4L, sizeof(SYN_ND4L) / sizeof(*SYN_ND4L)},
    {"ND5", SYN_ND5, sizeof(SYN_ND5) / sizeof(*SYN_ND5)},
    {"ND6", SYN_ND6, sizeof(SYN_ND6) / sizeof(*SYN_ND6)},
    {"ATP6", SYN_ATP6, sizeof(SYN_ATP6) / sizeof(*SYN_ATP6)},
    {"ATP8", SYN_ATP8, sizeof(SYN_ATP8) / sizeof(*SYN_ATP8)},
};

static void target_push(Targets *list, Target t) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 8;
        Target *p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len++] = t;
}

static char *ascii_upper(const char *s) {
    char *out = xstrdup(s);
    for (char *p = out; *p; p++) *p = (char)toupper((unsigned char)*p);
    return out;
}

static char *ascii_lower(const char *s) {
    char *out = xstrdup(s);
    for (char *p = out; *p; p++) *p = (char)tolower((unsigned char)*p);
    return out;
}

static void parse_targets(const char *genes, Targets *out) {
    char *copy = xstrdup(genes);
    char *tok = strtok(copy, ",");
    while (tok) {
        char *name = trim_ws_copy(tok);
        char *upper = ascii_upper(name);
        bool found = false;
        for (size_t i = 0; i < sizeof(SYNONYMS) / sizeof(*SYNONYMS); i++) {
            if (strcmp(upper, SYNONYMS[i].name)) continue;
            Target t = {0};
            t.name = xstrdup(SYNONYMS[i].name);
            t.nsyns = SYNONYMS[i].n;
            t.syns = xmalloc(t.nsyns * sizeof(*t.syns));
            for (size_t j = 0; j < t.nsyns; j++) t.syns[j] = xstrdup(SYNONYMS[i].syns[j]);
            target_push(out, t);
            found = true;
            break;
        }
        if (!found) {
            Target t = {0};
            t.name = upper;
            t.nsyns = 3;
            t.syns = xmalloc(3 * sizeof(*t.syns));
            t.syns[0] = xstrdup(name);
            t.syns[1] = xstrdup(upper);
            t.syns[2] = ascii_lower(name);
            target_push(out, t);
        } else free(upper);
        free(name);
        tok = strtok(NULL, ",");
    }
    free(copy);
}

static bool is_word_char(unsigned char c) {
    return isalnum(c) || c == '_';
}

/* Substring match with word boundaries: "subunit III" must not match "subunit I". */
static bool contains_gene_ci(const char *haystack, const char *needle) {
    char *h = ascii_lower(haystack), *n = ascii_lower(needle);
    size_t hn = strlen(h), nn = strlen(n);
    bool found = false;
    for (size_t i = 0; i + nn <= hn; i++) {
        if (strncmp(h + i, n, nn)) continue;
        if (i && is_word_char((unsigned char)h[i - 1])) continue;
        if (i + nn < hn && is_word_char((unsigned char)h[i + nn])) continue;
        found = true;
        break;
    }
    free(h);
    free(n);
    return found;
}

static const char *match_gene(const char *annotation, const Targets *targets) {
    char *ann = ascii_lower(annotation);
    const char *result = NULL;
    for (size_t i = 0; i < targets->len && !result; i++)
        for (size_t j = 0; j < targets->data[i].nsyns; j++) {
            char *syn = ascii_lower(targets->data[i].syns[j]);
            if (!strcmp(ann, syn)) { result = targets->data[i].name; free(syn); break; }
            free(syn);
        }
    if (!result)
        for (size_t i = 0; i < targets->len && !result; i++)
            for (size_t j = 0; j < targets->data[i].nsyns; j++)
                if (contains_gene_ci(ann, targets->data[i].syns[j])) {
                    result = targets->data[i].name;
                    break;
                }
    free(ann);
    return result;
}

static void build_canonical_targets(Targets *out) {
    for (size_t i = 0; i < sizeof(SYNONYMS) / sizeof(*SYNONYMS); i++) {
        Target t = {0};
        t.name = xstrdup(SYNONYMS[i].name);
        t.nsyns = SYNONYMS[i].n;
        t.syns = xmalloc(t.nsyns * sizeof(*t.syns));
        for (size_t j = 0; j < t.nsyns; j++) t.syns[j] = xstrdup(SYNONYMS[i].syns[j]);
        target_push(out, t);
    }
}

static bool target_index(const Targets *targets, const char *name, size_t *idx) {
    for (size_t i = 0; i < targets->len; i++)
        if (!strcmp(targets->data[i].name, name)) { *idx = i; return true; }
    return false;
}

static char *clean_name(const char *s) {
    char *out = xstrdup(s);
    for (char *p = out; *p; p++)
        if (*p == ' ' || *p == '/' || *p == '\\' || *p == ':') *p = '_';
    return out;
}

static bool is_trna(const Feature *f, const char *gene_name) {
    if (!strcmp(f->type, "tRNA")) return true;
    char *lower = ascii_lower(gene_name);
    bool ret = !strncmp(lower, "trn", 3) || strstr(lower, "trna") != NULL;
    free(lower);
    return ret;
}

static char *extract_sequence(const GBRecord *rec, const Feature *f) {
    const char *seq = rec->seq ? rec->seq : "";
    size_t slen = strlen(seq);
    size_t cap = 64, n = 0;
    char *out = xmalloc(cap);
    for (size_t i = 0; i < f->nlocs; i++) {
        size_t start = f->locs[i].start ? f->locs[i].start - 1 : 0;
        size_t end = f->locs[i].end < slen ? f->locs[i].end : slen;
        if (start < slen && end > start) {
            size_t need = end - start;
            while (n + need + 1 > cap) { cap *= 2; out = xrealloc(out, cap); }
            memcpy(out + n, seq + start, need);
            n += need;
            out[n] = '\0';
        }
    }
    if (f->complement) {
        char *rc = reverse_complement(out);
        free(out);
        return rc;
    }
    return out;
}

static char *location_string(const Feature *f) {
    if (f->nlocs == 0) return xstrdup("1-1");
    char *out = xmalloc(64);
    snprintf(out, 64, "%zu-%zu", f->locs[0].start, f->locs[f->nlocs - 1].end);
    return out;
}

typedef struct { char *base; int count; } NameCount;

static char *unique_filename(NameCount **usedp, size_t *nused, size_t *cap,
                             const char *base) {
    NameCount *used = *usedp;
    for (size_t i = 0; i < *nused; i++) {
        if (strcmp(used[i].base, base)) continue;
        used[i].count++;
        if (used[i].count == 1) return xstrdup(base);
        const char *dot = strrchr(base, '.');
        char *stem = dot ? str_ndup(base, (size_t)(dot - base)) : xstrdup(base);
        char *ext = dot ? xstrdup(dot) : xstrdup("");
        char *out = xmalloc(strlen(stem) + strlen(ext) + 24);
        sprintf(out, "%s_%d%s", stem, used[i].count - 1, ext);
        free(stem);
        free(ext);
        return out;
    }
    if (*nused == *cap) {
        *cap = *cap ? *cap * 2 : 16;
        NameCount *p = realloc(used, *cap * sizeof(*p));
        if (!p) xdie("out of memory");
        used = p;
        *usedp = p;
    }
    used[*nused].base = xstrdup(base);
    used[*nused].count = 1;
    (*nused)++;
    return xstrdup(base);
}

typedef struct {
    char *accession;
    char *gene;
    size_t length;
    char *status;
    char *filename;
} Result;

typedef struct {
    Result *data;
    size_t len;
    size_t cap;
} Results;

static void result_push(Results *list, const char *acc, const char *gene,
                        size_t length, const char *status, const char *file) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 32;
        Result *p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    Result r = {xstrdup(acc), xstrdup(gene), length, xstrdup(status), xstrdup(file)};
    list->data[list->len++] = r;
}

static int cmp_result(const void *a, const void *b) {
    const Result *x = a, *y = b;
    int c = strcmp(x->accession, y->accession);
    return c ? c : strcmp(x->gene, y->gene);
}

static void append_fasta(FILE *f, const char *header, const char *seq) {
    fprintf(f, ">%s\n", header);
    size_t len = strlen(seq);
    for (size_t i = 0; i < len; i += 60) {
        size_t n = len - i < 60 ? len - i : 60;
        fwrite(seq + i, 1, n, f);
        fputc('\n', f);
    }
}

static void write_fasta(const char *path, const char *header, const char *seq) {
    FILE *f = fopen(path, "wb");
    if (!f) xdie("Could not create %s: %s", path, strerror(errno));
    append_fasta(f, header, seq);
    fclose(f);
}

typedef struct {
    char *gene;
    char *seq;
    char *loc;
    char *header;
    bool cds;
    bool used;
} Hit;

typedef struct {
    Hit *data;
    size_t len;
    size_t cap;
} Hits;

static void hit_push(Hits *list, char *gene, char *seq, char *loc, char *header, bool cds) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        Hit *p = realloc(list->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        list->data = p;
        list->cap = cap;
    }
    Hit h = {gene, seq, loc, header, cds, false};
    list->data[list->len++] = h;
}

static void hits_free(Hits *list) {
    for (size_t i = 0; i < list->len; i++) {
        free(list->data[i].gene);
        free(list->data[i].seq);
        free(list->data[i].loc);
        free(list->data[i].header);
    }
    free(list->data);
    list->data = NULL;
    list->len = list->cap = 0;
}

typedef struct {
    char *name;
    FILE *f;
} GeneOut;

typedef struct {
    GeneOut *data;
    size_t len;
    size_t cap;
} GeneOuts;

static FILE *gene_out(GeneOuts *outs, const char *outdir, const char *gene) {
    for (size_t i = 0; i < outs->len; i++)
        if (!strcmp(outs->data[i].name, gene)) return outs->data[i].f;
    if (outs->len == outs->cap) {
        size_t cap = outs->cap ? outs->cap * 2 : 16;
        GeneOut *p = realloc(outs->data, cap * sizeof(*p));
        if (!p) xdie("out of memory");
        outs->data = p;
        outs->cap = cap;
    }
    char *name = xstrdup(gene);
    size_t file_len = strlen(gene) + sizeof(".fasta");
    char *file = xmalloc(file_len);
    snprintf(file, file_len, "%s.fasta", gene);
    char *path = cli_join_path(outdir, file);
    FILE *f = fopen(path, "wb");
    if (!f) xdie("Could not create %s: %s", path, strerror(errno));
    free(path);
    free(file);
    outs->data[outs->len].name = name;
    outs->data[outs->len].f = f;
    outs->len++;
    return f;
}

static void gene_outs_close(GeneOuts *outs) {
    for (size_t i = 0; i < outs->len; i++) {
        if (fclose(outs->data[i].f) != 0) xdie("Could not close %s.fasta", outs->data[i].name);
        free(outs->data[i].name);
    }
    free(outs->data);
    outs->data = NULL;
    outs->len = outs->cap = 0;
}

int cmd_rip(int argc, char **argv) {
    const char *genes_arg = NULL, *outdir = NULL;
    bool extract_all = false, skip_trna = false, aggregate = false;
    bool verbose = false, quiet = false;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-g") || !strcmp(a, "--genes")) genes_arg = next_value(argc, argv, &i);
        else if (!strncmp(a, "--genes=", 8)) genes_arg = a + 8;
        else if (!strcmp(a, "--all")) extract_all = true;
        else if (!strcmp(a, "--no-trna")) skip_trna = true;
        else if (!strcmp(a, "-a") || !strcmp(a, "--aggregate"))
            aggregate = true;
        else if (!strcmp(a, "-o") || !strcmp(a, "--output")) outdir = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) outdir = a + 9;
        else if (!strcmp(a, "-v") || !strcmp(a, "--verbose")) verbose = true;
        else if (!strcmp(a, "-q") || !strcmp(a, "--quiet")) quiet = true;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("rip"); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else {
            /* accept trailing positional inputs, like SeqRipper */
            continue;
        }
    }
    if (!extract_all && !genes_arg) xdie("rip: either -g/--genes or --all is required");
    if (!outdir) xdie("rip: -o/--output is required");

    Targets targets = {0};
    if (genes_arg) parse_targets(genes_arg, &targets);
    Targets canonical = {0};
    if (aggregate) build_canonical_targets(&canonical);
    if (mkdir_p(outdir) != 0) xdie("Could not create %s", outdir);

    GeneOuts outs = {0};
    NameCount *used = NULL;
    size_t nused = 0, ucap = 0;
    Results results = {0};

    for (int i = 1; i < argc; i++) {
        const char *path = argv[i];
        if (path[0] == '-') continue;
        /* skip option values */
        if ((!strcmp(argv[i - 1], "-g") || !strcmp(argv[i - 1], "--genes") ||
             !strcmp(argv[i - 1], "-o") || !strcmp(argv[i - 1], "--output"))) continue;

        GBRecord *records = NULL;
        size_t nrecords = 0;
        if (access(path, F_OK) != 0) {
            if (!quiet) fprintf(stderr, "Warning: file not found: %s\n", path);
            const char *slash = strrchr(path, '/');
            result_push(&results, slash ? slash + 1 : path, "-", 0, "error", "-");
            continue;
        }
        parse_genbank(path, &records, &nrecords);

        for (size_t r = 0; r < nrecords; r++) {
            GBRecord *rec = &records[r];
            if (verbose) fprintf(stderr, "  Processing: %s (%s)\n", rec->accession, rec->organism);
            bool *found = calloc(targets.len ? targets.len : 1, sizeof(*found));

            if (aggregate) {
                Hits hits = {0};
                for (size_t fi = 0; fi < rec->nfeatures; fi++) {
                    Feature *f = &rec->features[fi];
                    const char *name = f->gene ? f->gene : f->product;
                    if (!name) continue;
                    if (skip_trna && is_trna(f, name)) continue;

                    const char *m = match_gene(name, &canonical);
                    char *gene;
                    if (m) gene = xstrdup(m);
                    else if (genes_arg) {
                        const char *um = match_gene(name, &targets);
                        gene = um ? xstrdup(um) : clean_name(name);
                    } else gene = clean_name(name);

                    size_t t = 0;
                    if (genes_arg && !target_index(&targets, gene, &t)) {
                        free(gene);
                        continue;
                    }
                    char *seq = extract_sequence(rec, f);
                    if (!seq || !*seq) { free(seq); free(gene); continue; }
                    char *loc = location_string(f);
                    char *header = xmalloc(strlen(rec->accession) + strlen(loc) +
                                           strlen(rec->organism) + strlen(gene) + 16);
                    sprintf(header, "%s:%s %s [gene=%s]", rec->accession, loc,
                            rec->organism, gene);
                    hit_push(&hits, gene, seq, loc, header, !strcmp(f->type, "CDS"));
                }

                for (size_t hi = 0; hi < hits.len; hi++) {
                    Hit *h = &hits.data[hi];
                    if (h->used) continue;
                    size_t best = hi;
                    for (size_t hj = hi + 1; hj < hits.len; hj++) {
                        Hit *c = &hits.data[hj];
                        if (c->used || strcmp(c->gene, h->gene)) continue;
                        if ((c->cds && !hits.data[best].cds) ||
                            (c->cds == hits.data[best].cds &&
                             strlen(c->seq) > strlen(hits.data[best].seq)))
                            best = hj;
                    }
                    for (size_t hj = hi; hj < hits.len; hj++) {
                        Hit *c = &hits.data[hj];
                        if (!strcmp(c->gene, h->gene)) c->used = true;
                    }
                    Hit *w = &hits.data[best];
                    size_t flen = strlen(w->gene) + sizeof(".fasta");
                    char *filename = xmalloc(flen);
                    snprintf(filename, flen, "%s.fasta", w->gene);
                    append_fasta(gene_out(&outs, outdir, w->gene), w->header, w->seq);
                    result_push(&results, rec->accession, w->gene, strlen(w->seq),
                                "extracted", filename);
                    size_t t = 0;
                    if (genes_arg && target_index(&targets, w->gene, &t)) found[t] = true;
                    if (verbose) fprintf(stderr, "    %s -> %s (%zu bp)\n",
                                         w->gene, filename, strlen(w->seq));
                    free(filename);
                }
                hits_free(&hits);
            } else {
                for (size_t fi = 0; fi < rec->nfeatures; fi++) {
                    Feature *f = &rec->features[fi];
                    const char *name = f->gene ? f->gene : f->product;
                    if (!name) continue;
                    if (skip_trna && is_trna(f, name)) continue;

                    char *target_name = NULL;
                    size_t target_idx = 0;
                    if (extract_all) {
                        target_name = clean_name(name);
                    } else {
                        const char *m = match_gene(name, &targets);
                        if (!m) continue;
                        target_name = xstrdup(m);
                        for (size_t t = 0; t < targets.len; t++)
                            if (!strcmp(targets.data[t].name, m)) { target_idx = t; break; }
                        if (found[target_idx]) { free(target_name); continue; }
                    }

                    char *seq = extract_sequence(rec, f);
                    if (!seq || !*seq) { free(seq); free(target_name); continue; }

                    char *base = xmalloc(strlen(rec->accession) + strlen(target_name) + 8);
                    sprintf(base, "%s_%s.fasta", rec->accession, target_name);
                    char *filename = unique_filename(&used, &nused, &ucap, base);
                    free(base);
                    char *path_out = cli_join_path(outdir, filename);
                    char *loc = location_string(f);
                    char *header = xmalloc(strlen(rec->accession) + strlen(loc) +
                                           strlen(rec->organism) + strlen(target_name) + 16);
                    sprintf(header, "%s:%s %s [gene=%s]", rec->accession, loc,
                            rec->organism, target_name);
                    write_fasta(path_out, header, seq);
                    result_push(&results, rec->accession, target_name, strlen(seq),
                                "extracted", filename);
                    if (verbose) fprintf(stderr, "    %s -> %s (%zu bp)\n", target_name, filename, strlen(seq));
                    if (!extract_all) found[target_idx] = true;
                    free(header); free(loc); free(path_out); free(filename);
                    free(seq); free(target_name);
                }
            }

            if (!extract_all) {
                for (size_t t = 0; t < targets.len; t++)
                    if (!found[t])
                        result_push(&results, rec->accession, targets.data[t].name, 0,
                                    "missing", "-");
            }
            free(found);
        }
        for (size_t r = 0; r < nrecords; r++) record_free(&records[r]);
        free(records);
    }
    gene_outs_close(&outs);

    qsort(results.data, results.len, sizeof(*results.data), cmp_result);
    char *log_path = cli_join_path(outdir, "seqripper_log.tsv");
    FILE *logf = fopen(log_path, "wb");
    if (!logf) xdie("Could not create %s: %s", log_path, strerror(errno));
    fputs("Accession\tGene\tLength_bp\tStatus\tFilename\n", logf);
    for (size_t i = 0; i < results.len; i++)
        fprintf(logf, "%s\t%s\t%zu\t%s\t%s\n", results.data[i].accession,
                results.data[i].gene, results.data[i].length,
                results.data[i].status, results.data[i].filename);
    fclose(logf);

    if (!quiet) {
        size_t extracted = 0, errors = 0;
        for (size_t i = 0; i < results.len; i++) {
            if (!strcmp(results.data[i].status, "extracted")) extracted++;
            if (!strcmp(results.data[i].status, "error")) errors++;
        }
        fprintf(stderr, "SUMMARY:\nSequences extracted: %zu\n", extracted);
        if (errors) fprintf(stderr, "Errors: %zu\n", errors);
        fprintf(stderr, "Output written to: %s\n", outdir);
    }

    free(log_path);
    for (size_t i = 0; i < results.len; i++) {
        free(results.data[i].accession); free(results.data[i].gene);
        free(results.data[i].status); free(results.data[i].filename);
    }
    free(results.data);
    for (size_t i = 0; i < nused; i++) free(used[i].base);
    free(used);
    for (size_t i = 0; i < targets.len; i++) {
        free(targets.data[i].name);
        for (size_t j = 0; j < targets.data[i].nsyns; j++) free(targets.data[i].syns[j]);
        free(targets.data[i].syns);
    }
    free(targets.data);
    for (size_t i = 0; i < canonical.len; i++) {
        free(canonical.data[i].name);
        for (size_t j = 0; j < canonical.data[i].nsyns; j++) free(canonical.data[i].syns[j]);
        free(canonical.data[i].syns);
    }
    free(canonical.data);
    return 0;
}
