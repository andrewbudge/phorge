#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_help.h"
#include "cli_util.h"

#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static void lower_match(const char *s, char *out, size_t cap) {
    size_t n = 0;
    for (; *s && n + 1 < cap; s++, n++) out[n] = (char)tolower((unsigned char)*s);
    out[n] = '\0';
}

static char *trim_ws(char *s) {
    while (*s && isspace((unsigned char)*s)) s++;
    char *end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) end--;
    *end = '\0';
    return s;
}

static void strip_nl(char *s, ssize_t *len) {
    if (*len > 0 && s[*len - 1] == '\n') {
        (*len)--;
        if (*len > 0 && s[*len - 1] == '\r') (*len)--;
    }
    s[*len] = '\0';
}

static char *uppercase_dup(const char *s) {
    char *out = xstrdup(s);
    for (char *p = out; *p; p++) {
        if (*p >= 'a' && *p <= 'z') *p = (char)(*p - 'a' + 'A');
    }
    return out;
}

static void convert_push(SequenceList *seqs, const char *header, const char *seq) {
    char *h = xstrdup(header);
    char *s = uppercase_dup(seq);
    if (sequence_list_push(seqs, h, s) != 0) xdie("out of memory");
}

/* Returns 1 if a record was added. */
static int push_from_fields(SequenceList *seqs, char *line, size_t min_fields) {
    char *fields[8];
    size_t n = 0;
    char *p = line;
    for (;;) {
        while (isspace((unsigned char)*p)) p++;
        char *start = p;
        while (*p && !isspace((unsigned char)*p)) p++;
        if (n < 8) fields[n] = start;
        n++;
        if (!*p) break;
        *p++ = '\0';
    }
    if (n < min_fields) return 0;
    convert_push(seqs, fields[0], fields[1]);
    return 1;
}

static void write_convert_fasta(const SequenceList *seqs) {
    for (size_t i = 0; i < seqs->len; i++) {
        printf(">%s\n%s\n", seqs->data[i].header, seqs->data[i].seq);
    }
}

static void write_convert_nexus(const SequenceList *seqs) {
    bool dna = is_dna(seqs) != 0;
    const char *datatype = dna ? "DNA" : "PROTEIN";
    const char *missing = dna ? "N" : "X";
    size_t length = seqs->len ? strlen(seqs->data[0].seq) : 0;
    printf("#NEXUS\nBEGIN DATA;\n");
    printf("  DIMENSIONS NTAX=%zu NCHAR=%zu;\n", seqs->len, length);
    printf("  FORMAT DATATYPE=%s MISSING=%s GAP=-;\n  MATRIX\n", datatype, missing);
    for (size_t i = 0; i < seqs->len; i++) {
        printf("  %s    %s\n", seqs->data[i].header, seqs->data[i].seq);
    }
    printf(";\nEND;\n");
}

static void write_convert_relaxed_phylip(const SequenceList *seqs) {
    size_t length = seqs->len ? strlen(seqs->data[0].seq) : 0;
    printf("%zu %zu\n", seqs->len, length);
    for (size_t i = 0; i < seqs->len; i++) {
        printf("%s    %s\n", seqs->data[i].header, seqs->data[i].seq);
    }
}

static void write_convert_strict_phylip(const SequenceList *seqs) {
    size_t length = seqs->len ? strlen(seqs->data[0].seq) : 0;
    printf("%zu %zu\n", seqs->len, length);
    for (size_t i = 0; i < seqs->len; i++) {
        char label[11];
        size_t n = strlen(seqs->data[i].header);
        if (n > 10) n = 10;
        memcpy(label, seqs->data[i].header, n);
        label[n] = '\0';
        printf("%-10s%s\n", label, seqs->data[i].seq);
    }
}

static bool starts_with_ci(const char *s, const char *prefix) {
    size_t n = strlen(prefix);
    for (size_t i = 0; i < n; i++) {
        if (!s[i] || tolower((unsigned char)s[i]) != tolower((unsigned char)prefix[i])) return false;
    }
    return true;
}

int cmd_convert(int argc, char **argv) {
    const char *input = NULL;
    const char *format = NULL;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-o") || !strcmp(a, "--output_format")) format = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output_format=", 16)) format = a + 16;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("convert"); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else if (!input) input = a;
        else xdie("unexpected argument: %s", a);
    }
    if (!input || !format) {
        fputs("error: the following required arguments were not provided:\n"
              "  --output_format <OUTPUT_FORMAT>\n"
              "\n"
              "Usage: phorge convert --output_format <OUTPUT_FORMAT> <INPUT_FILE>\n"
              "\n"
              "For more information, try '--help'.\n", stderr);
        exit(2);
    }

    TextInput input_stream;
    if (text_input_open(input, &input_stream) != 0) xdie("Failed to open file");
    char *line = NULL;
    size_t cap = 0;
    ssize_t nread = text_input_read_line(&input_stream, &line, &cap);
    if (nread == -1) {
        if (text_input_error(&input_stream) && input_stream.use_zlib)
            xdie("Failed to read input: %s", phorge_last_error());
        xdie("empty file");
    }
    strip_nl(line, &nread);
    char *first = trim_ws(line);

    SequenceList seqs = {0};
    if (first[0] == '>') {
        text_input_close(&input_stream);
        if (parse_fasta(input, false, &seqs, NULL) != 0) xdie("Failed to parse fasta file: %s", phorge_last_error());
    } else if (first[0] == '#') {
        bool in_matrix = false;
        while ((nread = text_input_read_line(&input_stream, &line, &cap)) != -1) {
            strip_nl(line, &nread);
            char *trimmed = trim_ws(line);
            if (!*trimmed) continue;
            if (starts_with_ci(trimmed, "MATRIX")) { in_matrix = true; continue; }
            if (!in_matrix) continue;
            if (trimmed[0] == ';') break;
            push_from_fields(&seqs, trimmed, 2);
        }
        if (text_input_error(&input_stream) && input_stream.use_zlib)
            xdie("Failed to read input: %s", phorge_last_error());
        text_input_close(&input_stream);
    } else if (isdigit((unsigned char)first[0])) {
        while ((nread = text_input_read_line(&input_stream, &line, &cap)) != -1) {
            strip_nl(line, &nread);
            char *trimmed = trim_ws(line);
            if (!*trimmed) continue;
            push_from_fields(&seqs, trimmed, 2);
        }
        if (text_input_error(&input_stream) && input_stream.use_zlib)
            xdie("Failed to read input: %s", phorge_last_error());
        text_input_close(&input_stream);
    } else {
        text_input_close(&input_stream);
        free(line);
        xdie("Error: could not detect input format");
    }
    free(line);

    char fmt[32];
    lower_match(format, fmt, sizeof(fmt));
    if (!strcmp(fmt, "f") || !strcmp(fmt, "fasta")) write_convert_fasta(&seqs);
    else if (!strcmp(fmt, "n") || !strcmp(fmt, "nexus") || !strcmp(fmt, "nex")) write_convert_nexus(&seqs);
    else if (!strcmp(fmt, "rp") || !strcmp(fmt, "phylip")) write_convert_relaxed_phylip(&seqs);
    else if (!strcmp(fmt, "sp")) write_convert_strict_phylip(&seqs);
    else xdie("Error: unknown output format '%s'. Use: f, n, sp, or rp", format);

    sequence_list_free(&seqs);
    return 0;
}

