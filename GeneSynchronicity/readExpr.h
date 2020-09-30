// -------------------------------------------------------------------------
// readExpr.h -   Header file for reading in expression values
//
// written by Sharlee Climer, March 2019
//
// ------------------------------------------------------------------------

#ifndef _READEXPR_H
#define _READEXPR_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>

#include "timer.h"

float** ReadAndNormalize(char*,int,int);

const int NUMHEADROWS = 1; // number of header rows before data values
const int NUMHEADCOLS = 1; // number of header columns before data values

const int QUIET = 0;  // set to one to eliminate output to screen
const int VERBOSE = 1;  // set to one to display maximum output to screen

const int STRNGLENGTH = 200; // maximum string length
const double TOL = 0.00001; // tolerance

// following used for error checking, can be modified as needed
const float MAX_EXP = 100000.0; // maximum expression level allowed in input
const float MIN_EXP = -100000.0; // minimum expression level allowed in input
const int MAX_NUM_INDIVIDUALS = 100000; // maximum number of individuals
const int MAX_NUM_GENES = 10000000; // maximum number of genes


inline void warning(const char* p) { fprintf(stderr,"Warning: %s \n",p); }
inline void fatal(const char* string) {fprintf(stderr,"\nFatal: %s\n\n",string); exit(1); }

#endif
