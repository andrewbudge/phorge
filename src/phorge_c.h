#ifndef PHORGE_C_H
#define PHORGE_C_H

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>

typedef struct {
    char *header;
    char *seq;
} Sequence;

typedef struct {
    Sequence *data;
    size_t len;
    size_t cap;
} SequenceList;

typedef struct {
    char **data;
    size_t len;
    size_t cap;
} StringList;

typedef enum {
    TAXON_GROUP_INGROUP = 0,
    TAXON_GROUP_OUTGROUP = 1
} TaxonGroup;

typedef struct {
    char *accession;
    char *taxon_name;
    char *gene_annotation;
    char *source_db;
    uint64_t taxid;
    uint64_t length;
    uint64_t query_taxid;
    int refseq;
    int taxonomic_outlier;
    TaxonGroup taxon_group;
} Accession;

typedef struct {
    Accession *data;
    size_t len;
    size_t cap;
} AccessionList;

/* Human-readable detail for the most recent failure, valid until next call. */
const char *phorge_last_error(void);

/* Transparent text input: plain files use getline(); gzip files use zlib.
   text_input_open() sniffs the gzip magic bytes. */
typedef struct {
    void *stream;   /* FILE* or gzFile; opaque to callers */
    int use_zlib;
} TextInput;

/* 0 on success, -1 on failure (see phorge_last_error). */
int text_input_open(const char *path, TextInput *in);
/* Wraps stdin with zlib, which transparently accepts plain or gzip data. */
int text_input_stdin(TextInput *in);
/* getline()-compatible: returns line length including '\n', or -1. */
ssize_t text_input_read_line(TextInput *in, char **line, size_t *cap);
int text_input_error(TextInput *in);
void text_input_close(TextInput *in);

/* Returns 0 on success, -1 on failure (see phorge_last_error). */
int parse_fasta(const char *filename, int validate_equal, SequenceList *out, size_t *out_length);
/* Takes ownership of both strings. */
int sequence_list_push(SequenceList *list, char *header, char *seq);
void sequence_list_free(SequenceList *list);

int is_dna(const SequenceList *sequences);
/* Caller frees the returned malloc'd string. */
char *reverse_complement(const char *seq);

int load_taxa_list(const char *filename, StringList *out);
void string_list_free(StringList *list);

int read_query_tsv(const char *path, AccessionList *out);
void accession_list_free(AccessionList *list);
/* Caller frees the returned malloc'd string. */
char *render_query_tsv(const AccessionList *records);

/* CLI helpers. */
int phorge_is_missing_char(unsigned char c);
char *phorge_basename_dup(const char *path);

#endif /* PHORGE_C_H */
