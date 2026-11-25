#ifndef COMMANDS_H
#define COMMANDS_H

int cmd_getheaders(int argc, char **argv);
int cmd_stats(int argc, char **argv);
int cmd_filter(int argc, char **argv);
int cmd_curate(int argc, char **argv);
int cmd_convert(int argc, char **argv);
int cmd_concat(int argc, char **argv);
int cmd_coverage(int argc, char **argv);
int cmd_align(int argc, char **argv);
int cmd_extract(int argc, char **argv);
int cmd_dump(int argc, char **argv);
int cmd_unpack(int argc, char **argv);
int cmd_harvest(int argc, char **argv);
int cmd_report(int argc, char **argv);
int cmd_orthologs(int argc, char **argv);
int cmd_rip(int argc, char **argv);

#endif /* COMMANDS_H */
