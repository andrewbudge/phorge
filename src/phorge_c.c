#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

static char error_buf[1024];

static void set_error(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(error_buf, sizeof(error_buf), fmt, ap);
    va_end(ap);
}

const char *phorge_last_error(void) { return error_buf; }

typedef struct {
    char *data;
    size_t len;
    size_t cap;
} Buf;

static void buf_free(Buf *b) {
    free(b->data);
    b->data = NULL;
    b->len = b->cap = 0;
}

static int buf_reserve(Buf *b, size_t need) {
    if (b->cap - b->len >= need) return 0;
    size_t cap = b->cap ? b->cap : 64;
    while (cap - b->len < need) {
        if (cap > (size_t)-1 / 2) {
            set_error("buffer too large");
            return -1;
        }
        cap *= 2;
    }
    char *p = realloc(b->data, cap);
    if (!p) {
        set_error("out of memory");
        return -1;
    }
    b->data = p;
    b->cap = cap;
    return 0;
}

static int buf_putc(Buf *b, char c) {
    if (buf_reserve(b, 2) != 0) return -1;
    b->data[b->len++] = c;
    b->data[b->len] = '\0';
    return 0;
}

static int buf_puts(Buf *b, const char *s) {
    size_t n = strlen(s);
    if (buf_reserve(b, n + 1) != 0) return -1;
    memcpy(b->data + b->len, s, n);
    b->len += n;
    b->data[b->len] = '\0';
    return 0;
}

static int buf_printf(Buf *b, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);
    if (n < 0 || buf_reserve(b, (size_t)n + 1) != 0) return -1;
    va_start(ap, fmt);
    vsnprintf(b->data + b->len, (size_t)n + 1, fmt, ap);
    va_end(ap);
    b->len += (size_t)n;
    return 0;
}

static char *buf_take(Buf *b) {
    char *p = b->data;
    b->data = NULL;
    b->len = b->cap = 0;
    return p;
}

static char *str_dup(const char *s) {
    char *p = malloc(strlen(s) + 1);
    if (!p) {
        set_error("out of memory");
        return NULL;
    }
    strcpy(p, s);
    return p;
}

/* Rust BufRead::lines(): strip '\n', and strip a preceding '\r' only when a
   newline was actually present. */
static void strip_line_ending(char *s, ssize_t *len) {
    if (*len > 0 && s[*len - 1] == '\n') {
        (*len)--;
        if (*len > 0 && s[*len - 1] == '\r') (*len)--;
    }
    s[*len] = '\0';
}

/* Same rule, for rows carved out of an in-memory buffer: the '\n' has already
   been located but not yet overwritten. */
static void strip_crlf(char *s, ssize_t *len, int had_nl) {
    if (had_nl && *len > 0 && s[*len - 1] == '\r') (*len)--;
    s[*len] = '\0';
}

static char *trim_ascii(char *s) {
    while (*s && isspace((unsigned char)*s)) s++;
    char *end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) end--;
    *end = '\0';
    return s;
}

static void uppercase_ascii(char *s) {
    for (; *s; s++) {
        if (*s >= 'a' && *s <= 'z') *s = (char)(*s - 'a' + 'A');
    }
}

/* ---------------------------- text input -------------------------------- */

static int line_buf_reserve(char **buf, size_t *cap, size_t need) {
    if (*cap >= need) return 0;
    size_t newcap = *cap ? *cap : 256;
    while (newcap < need) {
        if (newcap > (size_t)-1 / 2) return -1;
        newcap *= 2;
    }
    char *p = realloc(*buf, newcap);
    if (!p) return -1;
    *buf = p;
    *cap = newcap;
    return 0;
}

int text_input_open(const char *path, TextInput *in) {
    in->stream = NULL;
    in->use_zlib = 0;

    FILE *f = fopen(path, "rb");
    if (!f) {
        set_error("Could not open %s: %s", path, strerror(errno));
        return -1;
    }

    unsigned char magic[2];
    size_t n = fread(magic, 1, sizeof(magic), f);
    if (n == 2 && magic[0] == 0x1f && magic[1] == 0x8b) {
        fclose(f);
        gzFile gz = gzopen(path, "rb");
        if (!gz) {
            set_error("Could not open %s: %s", path, strerror(errno));
            return -1;
        }
        in->stream = gz;
        in->use_zlib = 1;
        return 0;
    }

    /* Put the sniffed bytes back so the plain-file path is byte-identical. */
    if (n > 0 && fseek(f, 0, SEEK_SET) != 0) {
        for (size_t i = n; i > 0; i--) {
            if (ungetc(magic[i - 1], f) == EOF) {
                set_error("Error reading file: %s", strerror(errno));
                fclose(f);
                return -1;
            }
        }
        clearerr(f);
    }
    in->stream = f;
    return 0;
}

int text_input_stdin(TextInput *in) {
    int fd = dup(STDIN_FILENO);
    if (fd < 0) {
        set_error("Could not read stdin: %s", strerror(errno));
        return -1;
    }
    gzFile gz = gzdopen(fd, "rb");
    if (!gz) {
        set_error("Could not read stdin: %s", strerror(errno));
        close(fd);
        return -1;
    }
    /* gzdopen() accepts both plain and gzip streams without sniffing. */
    in->stream = gz;
    in->use_zlib = 1;
    return 0;
}

ssize_t text_input_read_line(TextInput *in, char **line, size_t *cap) {
    if (!in->use_zlib) return getline(line, cap, (FILE *)in->stream);

    gzFile gz = (gzFile)in->stream;
    size_t len = 0;
    if (line_buf_reserve(line, cap, 1) != 0) {
        set_error("out of memory");
        return -1;
    }
    (*line)[0] = '\0';

    for (;;) {
        size_t avail = *cap - len;
        if (avail > (size_t)INT_MAX) avail = (size_t)INT_MAX;
        char *dest = *line + len;
        char *got = gzgets(gz, dest, (int)avail);
        if (!got) {
            int errnum;
            const char *msg = gzerror(gz, &errnum);
            if (gzeof(gz) && (errnum == Z_OK || errnum == Z_STREAM_END)) {
                return len ? (ssize_t)len : -1;
            }
            set_error("Error reading file: %s", msg && *msg ? msg : "gzip error");
            return -1;
        }

        size_t n = strlen(dest);
        len += n;
        if (n && dest[n - 1] == '\n') break;
        if (gzeof(gz)) {
            if (text_input_error(in)) {
                set_error("Error reading file: truncated gzip stream");
                return -1;
            }
            return len ? (ssize_t)len : -1;
        }
        if (n + 1 < avail) {
            /* gzgets() stopped early without a newline: only a NUL can do that. */
            set_error("Error reading file: NUL byte in input");
            return -1;
        }
        if (line_buf_reserve(line, cap, len + 2) != 0) {
            set_error("out of memory");
            return -1;
        }
    }
    return (ssize_t)len;
}

int text_input_error(TextInput *in) {
    if (!in->stream) return 0;
    if (!in->use_zlib) return ferror((FILE *)in->stream) != 0;
    int errnum;
    gzerror((gzFile)in->stream, &errnum);
    return errnum != Z_OK && errnum != Z_STREAM_END;
}

void text_input_close(TextInput *in) {
    if (!in->stream) return;
    if (in->use_zlib) gzclose((gzFile)in->stream);
    else fclose((FILE *)in->stream);
    in->stream = NULL;
}

/* ------------------------------- FASTA ---------------------------------- */

int sequence_list_push(SequenceList *list, char *header, char *seq) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        Sequence *p = realloc(list->data, cap * sizeof(*p));
        if (!p) {
            set_error("out of memory");
            return -1;
        }
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len].header = header;
    list->data[list->len].seq = seq;
    list->len++;
    return 0;
}

void sequence_list_free(SequenceList *list) {
    for (size_t i = 0; i < list->len; i++) {
        free(list->data[i].header);
        free(list->data[i].seq);
    }
    free(list->data);
    list->data = NULL;
    list->len = list->cap = 0;
}

int parse_fasta(const char *filename, int validate_equal, SequenceList *out, size_t *out_length) {
    TextInput input;
    if (text_input_open(filename, &input) != 0) return -1;

    out->data = NULL;
    out->len = out->cap = 0;
    if (out_length) *out_length = 0;

    char *line = NULL;
    size_t line_cap = 0;
    ssize_t nread;
    char *header = NULL; /* NULL represents Rust's empty current_header */
    Buf seq = {0};
    size_t expected_length = 0;
    int have_expected = 0;
    int rc = -1;

    while ((nread = text_input_read_line(&input, &line, &line_cap)) != -1) {
        strip_line_ending(line, &nread);
        char *trimmed = trim_ascii(line);

        if (trimmed[0] == '>') {
            if (header != NULL) {
                char *seq_copy = str_dup(seq.data ? seq.data : "");
                if (!seq_copy) goto done;
                uppercase_ascii(seq_copy);
                if (validate_equal) {
                    size_t len = strlen(seq_copy);
                    if (!have_expected) {
                        expected_length = len;
                        have_expected = 1;
                    } else if (len != expected_length) {
                        set_error("Error: Sequence length mismatch in %s : %s", filename, header);
                        free(seq_copy);
                        goto done;
                    }
                }
                if (sequence_list_push(out, header, seq_copy) != 0) {
                    free(seq_copy);
                    goto done;
                }
                header = NULL;
            }
            char *new_header = trimmed + 1;
            header = new_header[0] ? str_dup(new_header) : NULL;
            if (new_header[0] && !header) goto done;
            seq.len = 0;
            if (seq.data) seq.data[0] = '\0';
        } else if (trimmed[0] != '\0') {
            if (buf_puts(&seq, trimmed) != 0) goto done;
        }
    }
    if (text_input_error(&input)) {
        if (!input.use_zlib) set_error("Error reading file: %s", strerror(errno));
        goto done;
    }

    if (header != NULL) {
        char *seq_copy = str_dup(seq.data ? seq.data : "");
        if (!seq_copy) goto done;
        uppercase_ascii(seq_copy);
        if (validate_equal && have_expected && strlen(seq_copy) != expected_length) {
            set_error("Error: Sequence length mismatch in %s : %s", filename, header);
            free(seq_copy);
            goto done;
        }
        if (sequence_list_push(out, header, seq_copy) != 0) {
            free(seq_copy);
            goto done;
        }
        header = NULL;
    }

    if (out_length) *out_length = out->len ? strlen(out->data[0].seq) : 0;
    rc = 0;

done:
    free(header);
    buf_free(&seq);
    free(line);
    text_input_close(&input);
    if (rc != 0) sequence_list_free(out);
    return rc;
}

/* ---------------------------- sequence ops ------------------------------ */

int is_dna(const SequenceList *sequences) {
    static const char *allowed = "ATCGN-RYSWKMBDHV?";
    for (size_t i = 0; i < sequences->len; i++) {
        const char *s = sequences->data[i].seq;
        for (; *s; s++) {
            if (!strchr(allowed, *s)) return 0;
        }
    }
    return 1;
}

char *reverse_complement(const char *seq) {
    size_t n = strlen(seq);
    char *out = malloc(n + 1);
    if (!out) {
        set_error("out of memory");
        return NULL;
    }
    out[n] = '\0';
    for (size_t i = 0; i < n; i++) {
        unsigned char c = (unsigned char)seq[n - 1 - i];
        char comp;
        switch (c) {
            case 'A': comp = 'T'; break;
            case 'T': comp = 'A'; break;
            case 'C': comp = 'G'; break;
            case 'G': comp = 'C'; break;
            case 'U': comp = 'A'; break;
            case 'R': comp = 'Y'; break;
            case 'Y': comp = 'R'; break;
            case 'K': comp = 'M'; break;
            case 'M': comp = 'K'; break;
            case 'B': comp = 'V'; break;
            case 'V': comp = 'B'; break;
            case 'D': comp = 'H'; break;
            case 'H': comp = 'D'; break;
            default: comp = (char)c; break;
        }
        out[i] = comp;
    }
    return out;
}

/* ---------------------------- taxa list --------------------------------- */

static int string_list_push(StringList *list, char *s) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        char **p = realloc(list->data, cap * sizeof(*p));
        if (!p) {
            set_error("out of memory");
            return -1;
        }
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len++] = s;
    return 0;
}

void string_list_free(StringList *list) {
    for (size_t i = 0; i < list->len; i++) free(list->data[i]);
    free(list->data);
    list->data = NULL;
    list->len = list->cap = 0;
}

int load_taxa_list(const char *filename, StringList *out) {
    TextInput input;
    if (text_input_open(filename, &input) != 0) return -1;
    out->data = NULL;
    out->len = out->cap = 0;
    char *line = NULL;
    size_t cap = 0;
    ssize_t nread;
    int rc = -1;
    while ((nread = text_input_read_line(&input, &line, &cap)) != -1) {
        strip_line_ending(line, &nread);
        char *trimmed = trim_ascii(line);
        if (*trimmed) {
            char *copy = str_dup(trimmed);
            if (!copy) goto done;
            if (string_list_push(out, copy) != 0) {
                free(copy);
                goto done;
            }
        }
    }
    if (text_input_error(&input)) {
        if (!input.use_zlib) set_error("Error reading file: %s", strerror(errno));
        goto done;
    }
    rc = 0;
done:
    free(line);
    text_input_close(&input);
    if (rc != 0) string_list_free(out);
    return rc;
}

/* ------------------------------ query TSV ------------------------------- */

static const char *TSV_COLUMNS[10] = {
    "accession", "taxon_name", "taxid", "length", "gene_annotation",
    "refseq", "source_db", "query_taxid", "taxon_group", "taxonomic_outlier"
};

enum {
    COL_ACCESSION, COL_TAXON_NAME, COL_TAXID, COL_LENGTH, COL_GENE_ANNOTATION,
    COL_REFSEQ, COL_SOURCE_DB, COL_QUERY_TAXID, COL_TAXON_GROUP, COL_OUTLIER,
    NCOLS
};

static char *read_whole_file(const char *path) {
    TextInput input;
    if (text_input_open(path, &input) != 0) {
        set_error("reading %s: %s", path, strerror(errno));
        return NULL;
    }
    Buf b = {0};
    char *line = NULL;
    size_t cap = 0;
    ssize_t n;
    while ((n = text_input_read_line(&input, &line, &cap)) != -1) {
        if (buf_reserve(&b, (size_t)n + 1) != 0) goto fail;
        memcpy(b.data + b.len, line, (size_t)n);
        b.len += (size_t)n;
        b.data[b.len] = '\0';
    }
    if (text_input_error(&input)) {
        set_error("reading %s: failed", path);
        goto fail;
    }
    free(line);
    text_input_close(&input);
    if (!b.data) {
        if (buf_putc(&b, '\0') != 0) {
            buf_free(&b);
            return NULL;
        }
        b.len = 0;
    }
    return buf_take(&b);
fail:
    free(line);
    buf_free(&b);
    text_input_close(&input);
    return NULL;
}

/* Split a mutable, NUL-terminated TSV row in place. Returns column count. */
static size_t split_tsv(char *line, char **cols, size_t max_cols) {
    size_t n = 0;
    char *start = line;
    for (char *p = line;; p++) {
        if (*p == '\t') {
            *p = '\0';
            if (n < max_cols) cols[n] = start;
            n++;
            start = p + 1;
        } else if (*p == '\0') {
            if (n < max_cols) cols[n] = start;
            n++;
            break;
        }
    }
    return n;
}

static int accession_push(AccessionList *list, Accession a) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        Accession *p = realloc(list->data, cap * sizeof(*p));
        if (!p) {
            set_error("out of memory");
            return -1;
        }
        list->data = p;
        list->cap = cap;
    }
    list->data[list->len++] = a;
    return 0;
}

void accession_list_free(AccessionList *list) {
    for (size_t i = 0; i < list->len; i++) {
        free(list->data[i].accession);
        free(list->data[i].taxon_name);
        free(list->data[i].gene_annotation);
        free(list->data[i].source_db);
    }
    free(list->data);
    list->data = NULL;
    list->len = list->cap = 0;
}

static int parse_u64_field(const char *v, const char *path, size_t line_no,
                           const char *name, uint64_t *out) {
    if (!*v) {
        *out = 0;
        return 0;
    }
    const char *p = v;
    if (*p == '+') p++;
    if (!*p) goto bad;
    for (const char *q = p; *q; q++) {
        if (*q < '0' || *q > '9') goto bad;
    }
    errno = 0;
    char *end = NULL;
    unsigned long long x = strtoull(p, &end, 10);
    if (errno == ERANGE || !end || *end) goto bad;
    *out = (uint64_t)x;
    return 0;
bad:
    set_error("%s:%zu: bad %s value \"%s\"", path, line_no, name, v);
    return -1;
}

int read_query_tsv(const char *path, AccessionList *out) {
    out->data = NULL;
    out->len = out->cap = 0;
    char *content = read_whole_file(path);
    if (!content) return -1;

    char *cursor = content;
    size_t remaining = strlen(content);
    char *header_line = NULL;
    if (remaining > 0) {
        size_t raw = strcspn(cursor, "\n");
        int has_nl = cursor[raw] == '\n';
        ssize_t hl = (ssize_t)raw;
        strip_crlf(cursor, &hl, has_nl);
        header_line = cursor;
        cursor += raw + (has_nl ? 1 : 0);
        remaining -= raw + (has_nl ? 1 : 0);
    }
    if (!header_line || !*header_line) {
        set_error("%s: empty query file (no header row)", path);
        free(content);
        return -1;
    }

    char *hcols[NCOLS];
    size_t nh = split_tsv(header_line, hcols, NCOLS);
    int indices[NCOLS];
    for (int i = 0; i < NCOLS; i++) indices[i] = -1;
    for (size_t i = 0; i < nh; i++) {
        for (int j = 0; j < NCOLS; j++) {
            if (strcmp(hcols[i], TSV_COLUMNS[j]) == 0) indices[j] = (int)i;
        }
    }
    if (indices[COL_ACCESSION] < 0) {
        set_error("%s: query TSV missing required 'accession' column", path);
        free(content);
        return -1;
    }

    int rc = -1;
    size_t row = 0;
    while (remaining > 0) {
        size_t raw = strcspn(cursor, "\n");
        int has_nl = cursor[raw] == '\n';
        ssize_t ll = (ssize_t)raw;
        strip_crlf(cursor, &ll, has_nl);
        char *next = cursor + raw + (has_nl ? 1 : 0);
        size_t next_remaining = remaining - raw - (has_nl ? 1 : 0);

        if (cursor[0] != '\0') {
            char *cols[NCOLS * 4];
            size_t nc = split_tsv(cursor, cols, NCOLS * 4);
            const char *empty = "";
#define GET(name_idx) \
    ((indices[name_idx] >= 0 && (size_t)indices[name_idx] < nc) \
         ? cols[indices[name_idx]] : empty)
            size_t line_no = row + 2;
            Accession a = {0};
            const char *acc = GET(COL_ACCESSION);
            const char *taxon = GET(COL_TAXON_NAME);
            const char *gene = GET(COL_GENE_ANNOTATION);
            const char *src = GET(COL_SOURCE_DB);
            a.accession = str_dup(acc);
            a.taxon_name = str_dup(taxon);
            a.gene_annotation = str_dup(gene);
            a.source_db = str_dup(src);
            if (!a.accession || !a.taxon_name || !a.gene_annotation || !a.source_db) {
                free(a.accession); free(a.taxon_name); free(a.gene_annotation); free(a.source_db);
                set_error("out of memory");
                goto done;
            }
            if (parse_u64_field(GET(COL_TAXID), path, line_no, "taxid", &a.taxid) != 0 ||
                parse_u64_field(GET(COL_LENGTH), path, line_no, "length", &a.length) != 0 ||
                parse_u64_field(GET(COL_QUERY_TAXID), path, line_no, "query_taxid", &a.query_taxid) != 0) {
                free(a.accession); free(a.taxon_name); free(a.gene_annotation); free(a.source_db);
                goto done;
            }
            a.refseq = strcmp(GET(COL_REFSEQ), "true") == 0;
            a.taxonomic_outlier = strcmp(GET(COL_OUTLIER), "true") == 0;
            const char *group = GET(COL_TAXON_GROUP);
            if (!*group || strcmp(group, "ingroup") == 0) {
                a.taxon_group = TAXON_GROUP_INGROUP;
            } else if (strcmp(group, "outgroup") == 0) {
                a.taxon_group = TAXON_GROUP_OUTGROUP;
            } else {
                set_error("%s:%zu: unknown taxon_group \"%s\" (expected ingroup/outgroup)",
                          path, line_no, group);
                free(a.accession); free(a.taxon_name); free(a.gene_annotation); free(a.source_db);
                goto done;
            }
            if (accession_push(out, a) != 0) {
                free(a.accession); free(a.taxon_name); free(a.gene_annotation); free(a.source_db);
                goto done;
            }
#undef GET
        }
        row++;
        cursor = next;
        remaining = next_remaining;
    }
    rc = 0;
done:
    free(content);
    if (rc != 0) accession_list_free(out);
    return rc;
}

static void append_clean(Buf *b, const char *s) {
    for (; *s; s++) {
        char c = (*s == '\t' || *s == '\n' || *s == '\r') ? ' ' : *s;
        if (buf_putc(b, c) != 0) return;
    }
}

char *render_query_tsv(const AccessionList *records) {
    static const char *header = "accession\ttaxon_name\ttaxid\tlength\tgene_annotation\t"
                                "refseq\tsource_db\tquery_taxid\ttaxon_group\ttaxonomic_outlier";
    Buf b = {0};
    if (buf_puts(&b, header) != 0 || buf_putc(&b, '\n') != 0) {
        buf_free(&b);
        return NULL;
    }
    for (size_t i = 0; i < records->len; i++) {
        const Accession *r = &records->data[i];
        append_clean(&b, r->accession);
        buf_putc(&b, '\t');
        append_clean(&b, r->taxon_name);
        buf_putc(&b, '\t');
        buf_printf(&b, "%llu", (unsigned long long)r->taxid);
        buf_putc(&b, '\t');
        buf_printf(&b, "%llu", (unsigned long long)r->length);
        buf_putc(&b, '\t');
        append_clean(&b, r->gene_annotation);
        buf_putc(&b, '\t');
        buf_puts(&b, r->refseq ? "true" : "false");
        buf_putc(&b, '\t');
        append_clean(&b, r->source_db);
        buf_putc(&b, '\t');
        buf_printf(&b, "%llu", (unsigned long long)r->query_taxid);
        buf_putc(&b, '\t');
        buf_puts(&b, r->taxon_group == TAXON_GROUP_OUTGROUP ? "outgroup" : "ingroup");
        buf_putc(&b, '\t');
        buf_puts(&b, r->taxonomic_outlier ? "true" : "false");
        buf_putc(&b, '\n');
    }
    return buf_take(&b);
}

/* ------------------------------ CLI helpers ----------------------------- */

int phorge_is_missing_char(unsigned char c) {
    return c == '-' || c == 'N' || c == 'X' || c == '?';
}

char *phorge_basename_dup(const char *path) {
    const char *slash = strrchr(path, '/');
    return str_dup(slash ? slash + 1 : path);
}
