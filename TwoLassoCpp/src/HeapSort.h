/* ========================================================================== */
/*                                                                            */
/*   HeapSort.h                                                               */
/*   (c) 2010 Alan Leanrcic                                                   */
/*                                                                            */
/*   Experiment in sorting XTResid in order                                   */
/*   Code modified from Numerical Recipes in C.                               */
/*   There is probably an extended C++ library function better for this.      */
/*                                                                            */
/* ========================================================================== */


#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <Rdefines.h>
  #include <Rinternals.h>
  #include <stdio.h>
  #define RMATH 0
#endif

int *HeapSortAll(int N, double *Doubles);
int ReverseString(int N, int *ReverseScores);
