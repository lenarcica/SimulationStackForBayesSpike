/* ========================================================================== */
/*                                                                            */
/*   HorseshoeC.h                                                               */
/*   (c) 2010 Alan Lenarcic                                                   */
/*                                                                            */
/*   h File for HorseShoe Integration                                         */
/*   Used to Integrate Lambdaj Densities                                      */
/*     This attempt to implement Horseshoe is not very good and it would  be  */
/*     recommend users look into users of TwoSimR5 for more professional      */
/*     versions of this estimator.                                            */
/*  In this implementation a n attempt to ingtegrate to do the penalty doesn't*/
/*    quite work as well as hoped.                                            */
/* ========================================================================== */
/******************************************************************************/
//// LICENSE INFO: C CODE
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  A copy of the GNU General Public License is available at
//  https://www.R-project.org/Licenses/
//
/******************************************************************************/


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

       

