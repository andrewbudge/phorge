#ifndef CLI_UTIL_H
#define CLI_UTIL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>

/* Shared CLI plumbing. */
void xdie(const char *fmt, ...);
uint64_t cli_hash_str(const char *s);
void cli_strip_line_ending(char *s, ssize_t *len);
void *xmalloc(size_t n);
char *xstrdup(const char *s);
char *str_lower(const char *s);
const char *next_value(int argc, char **argv, int *i);
int mkdir_p(const char *path);
char *cli_join_path(const char *dir, const char *name);
void stem_and_ext(const char *path, char **stem, char **ext);
const char *basename_of(const char *path);

/* Run a child with stdout/stderr redirected to the given paths (NULL = inherit).
   Returns 0 on success, 1 if the child failed, -1 if it could not be spawned. */
int cli_spawn_wait(const char *program, char **args, const char *stdout_path,
                   const char *stderr_path, bool append_stderr);

#endif /* CLI_UTIL_H */
