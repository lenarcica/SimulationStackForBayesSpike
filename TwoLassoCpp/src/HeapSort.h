/* ========================================================================== */
/*                                                                            */
/*   HeapSort.h                                                            */
/*   (c) 2010 Author                                                          */
/*                                                                            */
/*   Experiment in sorting XTResid in order                                                             */
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
