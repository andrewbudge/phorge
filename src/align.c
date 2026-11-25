#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_util.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cli_help.h"

typedef struct {
    void **data;
    size_t len;
    size_t cap;
} PtrList;

static void ptr_push(PtrList *list, void *p) {
    if (list->len == list->cap) {
        size_t cap = list->cap ? list->cap * 2 : 16;
        void **data = realloc(list->data, cap * sizeof(*data));
        if (!data) xdie("out of memory");
        list->data = data;
        list->cap = cap;
    }
    list->data[list->len++] = p;
}

static char *join_path(const char *dir, const char *name) {
    size_t n = strlen(dir), m = strlen(name);
    bool slash = n && dir[n - 1] == '/';
    char *out = xmalloc(n + (slash ? 0 : 1) + m + 1);
    if (slash) snprintf(out, n + m + 1, "%s%s", dir, name);
    else snprintf(out, n + 1 + m + 1, "%s/%s", dir, name);
    return out;
}

int cmd_align(int argc, char **argv) {
    const char *program = "mafft";
    const char *extension = "_aln";
    const char *output = NULL;
    PtrList inputs = {0};
    PtrList passthrough = {0};
    bool after_double_dash = false;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (after_double_dash) { ptr_push(&passthrough, (void *)a); continue; }
        if (!strcmp(a, "--")) { after_double_dash = true; continue; }
        if (!strcmp(a, "-p") || !strcmp(a, "--program")) program = next_value(argc, argv, &i);
        else if (!strncmp(a, "--program=", 10)) program = a + 10;
        else if (!strcmp(a, "-e") || !strcmp(a, "--extension")) extension = next_value(argc, argv, &i);
        else if (!strncmp(a, "--extension=", 12)) extension = a + 12;
        else if (!strcmp(a, "-o") || !strcmp(a, "--output")) output = next_value(argc, argv, &i);
        else if (!strncmp(a, "--output=", 9)) output = a + 9;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { print_help("align"); return 0; }
        else if (a[0] == '-') xdie("unrecognized argument: %s", a);
        else ptr_push(&inputs, (void *)a);
    }
    if (!output) xdie("align: --output/-o is required");
    if (inputs.len == 0) xdie("align: at least one input file is required");

    char *version_args[] = {(char *)program, (char *)"--version", NULL};
    if (cli_spawn_wait(program, version_args, "/dev/null", "/dev/null", false) < 0)
        xdie("Error: '%s' not found. Make sure it is installed and in your PATH.", program);

    if (mkdir_p(output) != 0) xdie("Could not create output directory");
    char *log_path = join_path(output, "align.log");
    FILE *logf = fopen(log_path, "wb");
    if (!logf) xdie("Could not create log file");
    fclose(logf);

    char *pname = str_lower(basename_of(program));
    bool is_muscle = strstr(pname, "muscle") != NULL;
    free(pname);

    for (size_t i = 0; i < inputs.len; i++) {
        const char *input = inputs.data[i];
        char *stem, *ext;
        stem_and_ext(input, &stem, &ext);
        size_t name_len = strlen(stem) + strlen(extension) + 1 + strlen(ext) + 1;
        char *name = xmalloc(name_len);
        snprintf(name, name_len, "%s%s.%s", stem, extension, ext);
        char *out_path = join_path(output, name);

        fprintf(stderr, "Aligning %s...", stem);

        size_t nargs = 1 + (passthrough.len ? passthrough.len : 1) + 1 + (is_muscle ? 2 : 0) + 1;
        char **args = xmalloc(nargs * sizeof(*args));
        size_t k = 0;
        args[k++] = (char *)program;
        if (passthrough.len) {
            for (size_t j = 0; j < passthrough.len; j++) args[k++] = passthrough.data[j];
        } else {
            args[k++] = (char *)(is_muscle ? "-align" : "--auto");
        }
        args[k++] = (char *)input;
        if (is_muscle) { args[k++] = "-output"; args[k++] = out_path; }
        args[k] = NULL;

        int status = cli_spawn_wait(program, args, is_muscle ? NULL : out_path, log_path, true);
        if (status < 0) xdie("Failed to run aligner");
        if (status != 0) {
            fprintf(stderr, "FAILED (see %s)\n", log_path);
            exit(1);
        }
        fputs("done\n", stderr);

        free(args);
        free(out_path);
        free(name);
        free(stem);
        free(ext);
    }

    fprintf(stderr, "Done. Aligned %zu files.\n", inputs.len);
    free(log_path);
    free(inputs.data);
    free(passthrough.data);
    return 0;
}
