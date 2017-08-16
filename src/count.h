#ifndef COUNT_H
#define COUNT_H

#ifndef DEF_H
#include "def.h"
#endif

int count_lines(char* fname);
int count_entries(char* fname);
int count_fr(int *fr, int frame, int readframe, int nread, int twin, int nbuffer, int freq);
int count_time(int flag, FILE *lfp);

#endif
