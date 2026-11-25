#define _POSIX_C_SOURCE 200809L

#include "cli_help.h"
#include "commands.h"
#include "version.h"

#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
    /* Default C buffering: line-buffered on a TTY, block-buffered when
       redirected. Output bytes are identical to Rust's LineWriter. */
    if (argc < 2) {
        print_usage();
        return 1;
    }
    if (!strcmp(argv[1], "getheaders") || !strcmp(argv[1], "ghd")) return cmd_getheaders(argc - 1, argv + 1);
    if (!strcmp(argv[1], "stats")) return cmd_stats(argc - 1, argv + 1);
    if (!strcmp(argv[1], "filter")) return cmd_filter(argc - 1, argv + 1);
    if (!strcmp(argv[1], "curate")) return cmd_curate(argc - 1, argv + 1);
    if (!strcmp(argv[1], "concat") || !strcmp(argv[1], "liger")) return cmd_concat(argc - 1, argv + 1);
    if (!strcmp(argv[1], "convert")) return cmd_convert(argc - 1, argv + 1);
    if (!strcmp(argv[1], "coverage")) return cmd_coverage(argc - 1, argv + 1);
    if (!strcmp(argv[1], "align") || !strcmp(argv[1], "aln")) return cmd_align(argc - 1, argv + 1);
    if (!strcmp(argv[1], "extract")) return cmd_extract(argc - 1, argv + 1);
    if (!strcmp(argv[1], "unpack")) return cmd_unpack(argc - 1, argv + 1);
    if (!strcmp(argv[1], "harvest")) return cmd_harvest(argc - 1, argv + 1);
    if (!strcmp(argv[1], "report")) return cmd_report(argc - 1, argv + 1);
    if (!strcmp(argv[1], "orthologs") || !strcmp(argv[1], "groups")) return cmd_orthologs(argc - 1, argv + 1);
    if (!strcmp(argv[1], "rip")) return cmd_rip(argc - 1, argv + 1);
    if (!strcmp(argv[1], "--version") || !strcmp(argv[1], "-V")) {
        printf("phorge %s\n", PHORGE_VERSION);
        return 0;
    }
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help") || !strcmp(argv[1], "help")) {
        print_help(NULL);
        return 0;
    }
    fprintf(stderr, "phorge: unknown command '%s'\n"
                    "Run 'phorge --help' for the command list.\n", argv[1]);
    return 1;
}
