#ifndef DATASET_H
#define DATASET_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* Shared reader for NCBI Datasets packages: a zip path or an unpacked dir. */

typedef struct {
    char *key;
    char *val;
} DsKV;

typedef struct {
    DsKV *data;
    size_t len;
    size_t cap;
} DsKVList;

void ds_kv_free(DsKVList *list);
const char *ds_kv_get(const DsKVList *list, const char *key);
/* Parse JSON object scalars into dotted paths; array values are kept raw. */
void ds_parse_json_scalars(const char *json, DsKVList *out);

typedef struct {
    char *gene_id;
    char *accession;
    char *group_id;
    char *group_method;
    char *taxid;
    char *taxname;
    char *symbol;
    char *description;
} DsMeta;

typedef struct {
    DsMeta *data;
    size_t len;
    size_t cap;
} DsMetaList;

void ds_load_reports(const char *source, DsMetaList *out);
void ds_meta_free(DsMetaList *list);

enum { DS_META_TAXID = 0, DS_META_TAXNAME = 1, DS_META_DESC = 2, DS_META_SYMBOL = 3 };
const char *ds_meta_lookup(const DsMetaList *list, const char *key, int field);
/* Returns the first row whose gene_id or accession matches key. */
const DsMeta *ds_meta_find(const DsMetaList *list, const char *key);

typedef struct {
    char *header;   /* text after '>' */
    char *seq;      /* trimmed, concatenated, uppercased */
    char *source;   /* basename of the FASTA file it came from */
} DsFastaRecord;

typedef struct {
    DsFastaRecord *data;
    size_t len;
    size_t cap;
} DsFastaList;

/* Loads all FASTA records, or only the file whose basename equals `which`. */
void ds_load_fasta(const char *source, const char *which, DsFastaList *out);
void ds_fasta_free(DsFastaList *list);

char *ds_read_text(const char *path);
bool ds_is_zip(const char *source);
char *ds_header_tag(const char *header, const char *key);
char *ds_sanitize_name(const char *s);

#endif /* DATASET_H */
