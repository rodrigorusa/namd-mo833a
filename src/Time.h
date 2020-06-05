/*
 * This file is just instrumentation to capture paramount iterations time
 *
 */
#ifndef TIME_H
#define TIME_H

#include <sys/time.h>

extern double T_START_MAIN;
extern double T_INIT;
extern double T_LAST_PARAMOUNT;
extern double T_FINALIZE;
extern double T_PARAMOUNT_TOTAL;
extern int MAX_PI;

// Function to get time in seconds
inline static double mysecond() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

#endif