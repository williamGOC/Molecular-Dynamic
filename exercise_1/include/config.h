#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

// Size of the simulation box
#ifndef L
#define L 54.77225575
#endif

// 
#ifndef NC
#define NC 30
#endif

// 
#ifndef R_CUT
#define R_CUT 2.5
#endif

// 
#ifndef RHO 
#define RHO 0.3
#endif

// 
#define EMPTY -1
#define N 900

#endif // __CONFIG_H__