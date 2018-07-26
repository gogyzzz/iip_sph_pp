#ifndef IIP_TIME_H
#define IIP_TIME_H

#include "iip_type.h"

#if OS_UNIX
void stopwatch(int flag);
#else

LARGE_INTEGER getFILETIMEoffset();
int clock_gettime(struct timeval *tv);
void stopwatch(int flag);

#endif

long long get_micro_sec();
#endif
