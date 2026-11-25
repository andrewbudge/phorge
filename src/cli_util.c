#define _POSIX_C_SOURCE 200809L

#include "cli_util.h"

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <spawn.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

extern char **environ;

void xdie(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

uint64_t cli_hash_str(const char *s) {
    uint64_t h = 1469598103934665603ULL;
    for (; *s; s++) {
        h ^= (unsigned char)*s;
        h *= 1099511628211ULL;
    }
    return h;
}

void cli_strip_line_ending(char *s, ssize_t *len) {
    if (*len > 0 && s[*len - 1] == '\n') {
        (*len)--;
        if (*len > 0 && s[*len - 1] == '\r') (*len)--;
    }
    s[*len] = '\0';
}

void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) xdie("out of memory");
    return p;
}

char *xstrdup(const char *s) {
    size_t n = strlen(s) + 1;
    char *p = xmalloc(n);
    memcpy(p, s, n);
    return p;
}

char *str_lower(const char *s) {
    char *out = xmalloc(strlen(s) + 1);
    size_t i = 0;
    for (; *s; i++, s++) out[i] = (char)tolower((unsigned char)*s);
    out[i] = '\0';
    return out;
}

const char *next_value(int argc, char **argv, int *i) {
    if (*i + 1 >= argc) xdie("missing value for %s", argv[*i]);
    return argv[++(*i)];
}

int mkdir_p(const char *path) {
    char tmp[1024];
    size_t n = 0;
    for (const char *p = path;; p++) {
        if (*p == '/' || *p == '\0') {
            if (n > 0 && n < sizeof(tmp)) {
                tmp[n] = '\0';
                if (mkdir(tmp, 0777) != 0 && errno != EEXIST) return -1;
            }
            if (*p == '\0') break;
        }
        if (n + 1 < sizeof(tmp)) tmp[n++] = *p;
    }
    return 0;
}

char *cli_join_path(const char *dir, const char *name) {
    size_t n = strlen(dir), m = strlen(name);
    bool slash = n && dir[n - 1] == '/';
    char *out = xmalloc(n + (slash ? 0 : 1) + m + 1);
    if (slash) snprintf(out, n + m + 1, "%s%s", dir, name);
    else snprintf(out, n + 1 + m + 1, "%s/%s", dir, name);
    return out;
}

void stem_and_ext(const char *path, char **stem, char **ext) {
    const char *slash = strrchr(path, '/');
    const char *base = slash ? slash + 1 : path;
    const char *dot = strrchr(base, '.');
    if (dot && dot != base) {
        *stem = xmalloc((size_t)(dot - base) + 1);
        memcpy(*stem, base, (size_t)(dot - base));
        (*stem)[dot - base] = '\0';
        *ext = xstrdup(dot + 1);
    } else {
        *stem = xstrdup(base);
        *ext = xstrdup("fasta");
    }
}

const char *basename_of(const char *path) {
    const char *slash = strrchr(path, '/');
    return slash ? slash + 1 : path;
}

int cli_spawn_wait(const char *program, char **args, const char *stdout_path,
                   const char *stderr_path, bool append_stderr) {
    posix_spawn_file_actions_t actions;
    int rc = posix_spawn_file_actions_init(&actions);
    if (rc != 0) xdie("posix_spawn_file_actions_init: %s", strerror(rc));
    if (stdout_path && stderr_path && !strcmp(stdout_path, stderr_path)) {
        /* One open description shared by both streams, like Rust's try_clone. */
        int fd = open(stdout_path, O_WRONLY | O_CREAT | O_TRUNC, 0666);
        if (fd < 0) xdie("Could not open %s: %s", stdout_path, strerror(errno));
        rc = posix_spawn_file_actions_adddup2(&actions, fd, STDOUT_FILENO);
        if (rc != 0) xdie("posix_spawn_file_actions_adddup2: %s", strerror(rc));
        rc = posix_spawn_file_actions_adddup2(&actions, fd, STDERR_FILENO);
        if (rc != 0) xdie("posix_spawn_file_actions_adddup2: %s", strerror(rc));
        rc = posix_spawn_file_actions_addclose(&actions, fd);
        if (rc != 0) xdie("posix_spawn_file_actions_addclose: %s", strerror(rc));
    } else {
        if (stdout_path) {
            rc = posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, stdout_path,
                                                  O_WRONLY | O_CREAT | O_TRUNC, 0666);
            if (rc != 0) xdie("posix_spawn_file_actions_addopen: %s", strerror(rc));
        }
        if (stderr_path) {
            int flags = O_WRONLY | O_CREAT | (append_stderr ? O_APPEND : O_TRUNC);
            rc = posix_spawn_file_actions_addopen(&actions, STDERR_FILENO, stderr_path, flags, 0666);
            if (rc != 0) xdie("posix_spawn_file_actions_addopen: %s", strerror(rc));
        }
    }
    pid_t pid;
    rc = posix_spawnp(&pid, program, &actions, NULL, args, environ);
    posix_spawn_file_actions_destroy(&actions);
    if (rc != 0) return -1;

    int status;
    while (waitpid(pid, &status, 0) == -1) {
        if (errno != EINTR) return -1;
    }
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 0;
    return 1;
}
