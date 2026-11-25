#define _POSIX_C_SOURCE 200809L

#include "phorge_c.h"
#include "cli_help.h"
#include "cli_util.h"

#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --------------------------- string hash set ----------------------------- */

typedef struct {
    char **keys;
    size_t cap;
    size_t len;
} StringSet;

static void set_grow(StringSet *set) {
    size_t cap = set->cap ? set->cap * 2 : 64;
    char **keys = calloc(cap, sizeof(*keys));
    if (!keys) xdie("out of memory");
    for (size_t i = 0; i < set->cap; i++) {
        if (!set->keys[i]) continue;
        size_t j = (size_t)(cli_hash_str(set->keys[i]) & (cap - 1));
        while (keys[j]) j = (j + 1) & (cap - 1);
        keys[j] = set->keys[i];
    }
    free(set->keys);
    set->keys = keys;
    set->cap = cap;
}

/* Returns 1 if s was newly inserted, 0 if already present. */
static int set_insert(StringSet *set, const char *s) {
    if (!set->cap || (set->len + 1) * 10 >= set->cap * 7) set_grow(set);
    size_t j = (size_t)(cli_hash_str(s) & (set->cap - 1));
    while (set->keys[j]) {
        if (strcmp(set->keys[j], s) == 0) return 0;
        j = (j + 1) & (set->cap - 1);
    }
    set->keys[j] = xstrdup(s);
    set->len++;
    return 1;
}

static void set_free(StringSet *set) {
    for (size_t i = 0; i < set->cap; i++) free(set->keys[i]);
    free(set->keys);
    set->keys = NULL;
    set->cap = set->len = 0;
}

int cmd_getheaders(int argc, char **argv) {
    bool unique = false;
    const char **files = calloc((size_t)argc + 1, sizeof(*files));
    if (!files) xdie("out of memory");
    size_t nfiles = 0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-u") || !strcmp(a, "--unique")) {
            unique = true;
        } else if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            print_help("getheaders");
            free(files);
            return 0;
        } else if (a[0] == '-') {
            xdie("unrecognized argument: %s", a);
        } else {
            files[nfiles++] = a;
        }
    }

    StringSet seen = {0};
    char *line = NULL;
    size_t cap = 0;
    ssize_t nread;

    TextInput *streams = calloc(nfiles ? nfiles : 1, sizeof(*streams));
    if (!streams) xdie("out of memory");
    size_t nstreams = 0;
    if (nfiles == 0) {
        if (text_input_stdin(&streams[nstreams]) != 0) xdie("%s", phorge_last_error());
        nstreams++;
    } else {
        for (size_t i = 0; i < nfiles; i++) {
            if (text_input_open(files[i], &streams[nstreams]) != 0)
                xdie("Error: could not open '%s': %s", files[i], strerror(errno));
            nstreams++;
        }
    }

    for (size_t s = 0; s < nstreams; s++) {
        while ((nread = text_input_read_line(&streams[s], &line, &cap)) != -1) {
            cli_strip_line_ending(line, &nread);
            if (line[0] != '>') continue;
            const char *header = line + 1;
            if (!unique || set_insert(&seen, header)) puts(header);
        }
        if (text_input_error(&streams[s])) {
            if (!streams[s].use_zlib) xdie("Error reading file: %s", strerror(errno));
            xdie("%s", phorge_last_error());
        }
    }
    for (size_t i = 0; i < nstreams; i++) text_input_close(&streams[i]);
    free(streams);
    free(line);
    free(files);
    set_free(&seen);
    return 0;
}

