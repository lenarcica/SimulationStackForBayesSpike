/* ========================================================================== */
/*                                                                            */
/*   HorseshoeC.h                                                               */
/*   (c) 2010 Alan Lenarcic                                                          */
/*                                                                            */
/*   h File for HorseShoe Integration                                                              */
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
#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif

int GetAllLambdaj(SEXP SlZVals, SEXP SdlStandard, SEXP SBetaj, SEXP StauSq,
  SEXP SFillerDens, SEXP SUnifDraw, SEXP SOutLambdaj, SEXP Sexpn2lz2);
double GetLambdaj(int LengthZ2, double *dlStandard, double Betaj, 
     double tauSq, double *lZ2, double *FillerDens, double GetDraw, double *expn2lz2);
int AllGetManyLambdaj(SEXP SlZVals, SEXP SdlStandard, SEXP SOneBetaj, SEXP StauSq,
  SEXP SFillerDens, SEXP SortedUnifDraw, SEXP SOutLambdaj, SEXP Sexpn2lz2);
int GetManyLambdaj(int LengthZ2, double *dlStandard, double Betaj, 
     double tauSq, double *lZ2, double *FillerDens, 
     int CountDraws, double *SortedDraw, double *Output, double *expn2lz2);

       

