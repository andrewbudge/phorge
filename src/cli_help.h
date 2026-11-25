#ifndef CLI_HELP_H
#define CLI_HELP_H

void print_help(const char *cmd); /* NULL = command list, to stdout */
void print_usage(void);           /* command list, to stderr */

#endif /* CLI_HELP_H */
