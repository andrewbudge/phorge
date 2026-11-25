#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

static int failures = 0;
static char tmpdir[256];

#define CHECK(cond)                                                                     \
    do {                                                                                \
        if (!(cond)) {                                                                  \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);            \
            failures++;                                                                 \
        }                                                                               \
    } while (0)

static void write_file(const char *name, const char *content) {
    char path[512];
    snprintf(path, sizeof(path), "%s/%s", tmpdir, name);
    FILE *f = fopen(path, "wb");
    if (!f) {
        perror("fopen");
        exit(2);
    }
    fputs(content, f);
    fclose(f);
}

static const char *path_of(const char *name) {
    static char path[512];
    snprintf(path, sizeof(path), "%s/%s", tmpdir, name);
    return path;
}

static void write_gzip_file(const char *name, const char *content) {
    gzFile gz = gzopen(path_of(name), "wb");
    if (!gz) {
        perror("gzopen");
        exit(2);
    }
    int n = (int)strlen(content);
    if (gzwrite(gz, content, (unsigned)n) != n) {
        perror("gzwrite");
        exit(2);
    }
    gzclose(gz);
}

static void test_fasta_basics(void) {
    write_file("basic.fa", ">seq1 description\nacgt\nNN\n>seq2\nTTTT\n");
    SequenceList seqs = {0};
    size_t length = 0;
    CHECK(parse_fasta(path_of("basic.fa"), false, &seqs, &length) == 0);
    CHECK(seqs.len == 2);
    CHECK(length == 6);
    CHECK(strcmp(seqs.data[0].header, "seq1 description") == 0);
    CHECK(strcmp(seqs.data[0].seq, "ACGTNN") == 0);
    CHECK(strcmp(seqs.data[1].header, "seq2") == 0);
    CHECK(strcmp(seqs.data[1].seq, "TTTT") == 0);
    CHECK(is_dna(&seqs) == 1);
    sequence_list_free(&seqs);
}

static void test_fasta_crlf_and_no_final_newline(void) {
    write_file("crlf.fa", ">a\r\nac\r\ngt\r\n>b\r\nTT");
    SequenceList seqs = {0};
    size_t length = 0;
    CHECK(parse_fasta(path_of("crlf.fa"), false, &seqs, &length) == 0);
    CHECK(seqs.len == 2);
    CHECK(strcmp(seqs.data[0].seq, "ACGT") == 0);
    CHECK(strcmp(seqs.data[1].seq, "TT") == 0);
    sequence_list_free(&seqs);
}

static void test_fasta_validate_equal(void) {
    write_file("equal.fa", ">a\nACGT\n>b\nTTTT\n");
    SequenceList seqs = {0};
    size_t length = 0;
    CHECK(parse_fasta(path_of("equal.fa"), true, &seqs, &length) == 0);
    CHECK(seqs.len == 2 && length == 4);
    sequence_list_free(&seqs);

    write_file("unequal.fa", ">a\nACGT\n>b\nTTT\n");
    CHECK(parse_fasta(path_of("unequal.fa"), true, &seqs, &length) == -1);
    CHECK(strstr(phorge_last_error(), "Sequence length mismatch") != NULL);
    sequence_list_free(&seqs);
}

static void test_fasta_edge_cases(void) {
    write_file("empty.fa", "");
    SequenceList seqs = {0};
    size_t length = 0;
    CHECK(parse_fasta(path_of("empty.fa"), false, &seqs, &length) == 0);
    CHECK(seqs.len == 0 && length == 0);
    sequence_list_free(&seqs);

    write_file("header_only.fa", ">only\n");
    CHECK(parse_fasta(path_of("header_only.fa"), false, &seqs, &length) == 0);
    CHECK(seqs.len == 1 && strcmp(seqs.data[0].seq, "") == 0);
    sequence_list_free(&seqs);

    write_file("blank_header.fa", ">\nACGT\n");
    CHECK(parse_fasta(path_of("blank_header.fa"), false, &seqs, &length) == 0);
    CHECK(seqs.len == 0); /* Rust treats an empty header as "no record" */
    sequence_list_free(&seqs);

    CHECK(parse_fasta(path_of("does-not-exist.fa"), false, &seqs, &length) == -1);
    CHECK(strstr(phorge_last_error(), "Could not open") != NULL);
}

static void test_reverse_complement(void) {
    char *rc = reverse_complement("ATCGGGCTA");
    CHECK(strcmp(rc, "TAGCCCGAT") == 0);
    free(rc);
    rc = reverse_complement("RYKMBVDH");
    CHECK(strcmp(rc, "DHBVKMRY") == 0);
    free(rc);
    rc = reverse_complement("SWN");
    CHECK(strcmp(rc, "NWS") == 0);
    free(rc);
    rc = reverse_complement("AT-CG");
    CHECK(strcmp(rc, "CG-AT") == 0);
    free(rc);
    rc = reverse_complement("AT?CG");
    CHECK(strcmp(rc, "CG?AT") == 0);
    free(rc);
    const char *iupac = "ATCGNRYSWKMBDHV-?";
    rc = reverse_complement(iupac);
    CHECK(strlen(rc) == strlen(iupac));
    char *again = reverse_complement(rc);
    CHECK(strcmp(again, iupac) == 0);
    free(rc);
    free(again);
}

static void test_is_dna(void) {
    write_file("dna.fa", ">a\nATCGNRYSWKMBDHV-?\n");
    write_file("protein.fa", ">a\nMKFP\n");
    SequenceList seqs = {0};
    size_t length = 0;
    CHECK(parse_fasta(path_of("dna.fa"), false, &seqs, &length) == 0);
    CHECK(is_dna(&seqs) == 1);
    sequence_list_free(&seqs);
    CHECK(parse_fasta(path_of("protein.fa"), false, &seqs, &length) == 0);
    CHECK(is_dna(&seqs) == 0);
    sequence_list_free(&seqs);
}

static void test_taxa_list(void) {
    write_file("taxa.txt", "  Homo sapiens  \n\nMus musculus\n\t\n");
    StringList taxa = {0};
    CHECK(load_taxa_list(path_of("taxa.txt"), &taxa) == 0);
    CHECK(taxa.len == 2);
    CHECK(strcmp(taxa.data[0], "Homo sapiens") == 0);
    CHECK(strcmp(taxa.data[1], "Mus musculus") == 0);
    string_list_free(&taxa);
}

static void test_tsv_round_trip(void) {
    AccessionList records = {0};
    write_file("q.tsv",
               "accession\ttaxon_name\ttaxid\tlength\tgene_annotation\trefseq\tsource_db\t"
               "query_taxid\ttaxon_group\ttaxonomic_outlier\n"
               "AB123.1\tMus musculus\t10090\t658\tcytochrome c oxidase subunit I\tfalse\t"
               "GenBank\t10088\toutgroup\tfalse\n");
    CHECK(read_query_tsv(path_of("q.tsv"), &records) == 0);
    CHECK(records.len == 1);
    CHECK(strcmp(records.data[0].accession, "AB123.1") == 0);
    CHECK(records.data[0].taxid == 10090);
    CHECK(records.data[0].length == 658);
    CHECK(records.data[0].taxon_group == TAXON_GROUP_OUTGROUP);
    CHECK(records.data[0].refseq == 0);
    char *rendered = render_query_tsv(&records);
    CHECK(strcmp(rendered,
                 "accession\ttaxon_name\ttaxid\tlength\tgene_annotation\trefseq\tsource_db\t"
                 "query_taxid\ttaxon_group\ttaxonomic_outlier\n"
                 "AB123.1\tMus musculus\t10090\t658\tcytochrome c oxidase subunit I\tfalse\t"
                 "GenBank\t10088\toutgroup\tfalse\n") == 0);
    free(rendered);
    accession_list_free(&records);
}

static void test_tsv_partial_and_reordered(void) {
    AccessionList records = {0};
    write_file("partial.tsv", "taxid\taccession\ttaxon_group\n10090\tAB999.1\toutgroup\n");
    CHECK(read_query_tsv(path_of("partial.tsv"), &records) == 0);
    CHECK(records.len == 1);
    CHECK(strcmp(records.data[0].accession, "AB999.1") == 0);
    CHECK(records.data[0].taxid == 10090);
    CHECK(records.data[0].length == 0);
    CHECK(records.data[0].taxon_group == TAXON_GROUP_OUTGROUP);
    CHECK(records.data[0].taxon_name[0] == '\0');
    accession_list_free(&records);
}

static void test_tsv_errors(void) {
    AccessionList records = {0};
    write_file("noacc.tsv", "taxid\tlength\n10090\t658\n");
    CHECK(read_query_tsv(path_of("noacc.tsv"), &records) == -1);
    CHECK(strstr(phorge_last_error(), "missing required 'accession'") != NULL);

    write_file("empty.tsv", "");
    CHECK(read_query_tsv(path_of("empty.tsv"), &records) == -1);
    CHECK(strstr(phorge_last_error(), "empty query file") != NULL);

    write_file("bad.tsv", "accession\ttaxid\nA1\tnot-a-number\n");
    CHECK(read_query_tsv(path_of("bad.tsv"), &records) == -1);
    CHECK(strstr(phorge_last_error(), "bad taxid value") != NULL);
}

static void test_fasta_gzip(void) {
    write_gzip_file("gz_basic.fa.gz", ">seq1 description\nacgt\nNN\n>seq2\nTTTT\n");
    SequenceList seqs = {0};
    size_t length = 0;
    CHECK(parse_fasta(path_of("gz_basic.fa.gz"), false, &seqs, &length) == 0);
    CHECK(seqs.len == 2 && length == 6);
    CHECK(strcmp(seqs.data[0].header, "seq1 description") == 0);
    CHECK(strcmp(seqs.data[0].seq, "ACGTNN") == 0);
    CHECK(strcmp(seqs.data[1].seq, "TTTT") == 0);
    sequence_list_free(&seqs);

    write_gzip_file("gz_crlf.fa.gz", ">a\r\nac\r\ngt\r\n>b\r\nTT");
    CHECK(parse_fasta(path_of("gz_crlf.fa.gz"), false, &seqs, &length) == 0);
    CHECK(seqs.len == 2 && strcmp(seqs.data[0].seq, "ACGT") == 0);
    CHECK(strcmp(seqs.data[1].seq, "TT") == 0);
    sequence_list_free(&seqs);

    write_gzip_file("gz_empty.fa.gz", "");
    CHECK(parse_fasta(path_of("gz_empty.fa.gz"), false, &seqs, &length) == 0);
    CHECK(seqs.len == 0 && length == 0);
    sequence_list_free(&seqs);

    /* Long lines exercise gzgets() buffer growth. */
    char *long_content = malloc(12000);
    CHECK(long_content != NULL);
    if (!long_content) return;
    long_content[0] = '\0';
    strcat(long_content, ">long\n");
    for (int i = 0; i < 10000; i++) long_content[6 + i] = (i % 2) ? 'a' : 'c';
    long_content[10006] = '\0';
    strcat(long_content, "\n");
    write_gzip_file("gz_long.fa.gz", long_content);
    free(long_content);
    CHECK(parse_fasta(path_of("gz_long.fa.gz"), true, &seqs, &length) == 0);
    CHECK(seqs.len == 1 && length == 10000);
    sequence_list_free(&seqs);

    write_gzip_file("gz_unequal.fa.gz", ">a\nACGT\n>b\nTTT\n");
    CHECK(parse_fasta(path_of("gz_unequal.fa.gz"), true, &seqs, &length) == -1);
    CHECK(strstr(phorge_last_error(), "Sequence length mismatch") != NULL);
    sequence_list_free(&seqs);
}

static void test_gzip_errors(void) {
    SequenceList seqs = {0};
    size_t length = 0;

    /* Magic bytes without a valid deflate stream. */
    FILE *bad = fopen(path_of("bad_magic.fa.gz"), "wb");
    CHECK(bad != NULL);
    if (bad) {
        const unsigned char junk[] = {0x1f, 0x8b, 'n', 'o', 't', 'g', 'z', 0};
        fwrite(junk, 1, sizeof(junk), bad);
        fclose(bad);
    }
    CHECK(parse_fasta(path_of("bad_magic.fa.gz"), false, &seqs, &length) == -1);
    sequence_list_free(&seqs);

    write_gzip_file("truncated.fa.gz", ">a\nACGTACGT\n");
    FILE *f = fopen(path_of("truncated.fa.gz"), "r+b");
    CHECK(f != NULL);
    if (f) {
        fseek(f, 0, SEEK_END);
        long size = ftell(f);
        CHECK(size > 5);
        CHECK(ftruncate(fileno(f), (off_t)(size - 5)) == 0);
        fclose(f);
    }
    CHECK(parse_fasta(path_of("truncated.fa.gz"), false, &seqs, &length) == -1);
    sequence_list_free(&seqs);
}

static void test_gzip_metadata(void) {
    write_gzip_file("gz_taxa.txt.gz", "  Homo sapiens  \n\nMus musculus\n");
    StringList taxa = {0};
    CHECK(load_taxa_list(path_of("gz_taxa.txt.gz"), &taxa) == 0);
    CHECK(taxa.len == 2 && strcmp(taxa.data[0], "Homo sapiens") == 0);
    string_list_free(&taxa);

    write_gzip_file("gz_q.tsv.gz",
                    "accession\ttaxon_name\ttaxid\tlength\tgene_annotation\trefseq\t"
                    "source_db\tquery_taxid\ttaxon_group\ttaxonomic_outlier\n"
                    "AB123.1\tMus musculus\t10090\t658\tcytochrome c oxidase subunit I\tfalse\t"
                    "GenBank\t10088\toutgroup\tfalse\n");
    AccessionList records = {0};
    CHECK(read_query_tsv(path_of("gz_q.tsv.gz"), &records) == 0);
    CHECK(records.len == 1 && records.data[0].taxid == 10090);
    accession_list_free(&records);
}

int main(void) {
    snprintf(tmpdir, sizeof(tmpdir), "/tmp/opencode/phorge-tests-%ld", (long)getpid());
    if (mkdir(tmpdir, 0700) != 0 && errno != EEXIST) {
        perror("mkdir");
        return 2;
    }
    test_fasta_basics();
    test_fasta_crlf_and_no_final_newline();
    test_fasta_validate_equal();
    test_fasta_edge_cases();
    test_reverse_complement();
    test_is_dna();
    test_taxa_list();
    test_tsv_round_trip();
    test_tsv_partial_and_reordered();
    test_tsv_errors();
    test_fasta_gzip();
    test_gzip_errors();
    test_gzip_metadata();

    if (failures) {
        fprintf(stderr, "%d test(s) failed\n", failures);
        return 1;
    }
    puts("All C parsing tests passed");
    return 0;
}
