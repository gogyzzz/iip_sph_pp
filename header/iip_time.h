/*
 * ===========================================================
 *           Copyright (c) 2018, __IIPLAB__
 *                All rights reserved.
 *
 * This Source Code Form is subject to the terms of
 * the Mozilla Public License, v. 2.0.
 * If a copy of the MPL was not distributed with this file,
 *  You can obtain one at http://mozilla.org/MPL/2.0/.
 * ===========================================================
 */
#ifndef IIP_TIME_H
#define IIP_TIME_H

#include "iip_type.h"

/* Usage of stopwatch()
 ****************************
 * stopwatch(0);
 *
 *  do something you want to test
 *
 * stopwatch(1);
 ****************************
 *
 * stopwatch() function
 * print elapsed time between
 * stopwatch(0) and stopwatch(1)
 *
 * */

#if OS_UNIX
void stopwatch(int flag);
#elif OS_WIN
LARGE_INTEGER get_filetime_offset();
int clock_gettime(struct timeval *tv);
void stopwatch(int flag);
#endif

/* return micro sec of currnet time, used for srand()*/
long long get_micro_sec();
#endif
