/* ========================================================================== */
/*                                                                            */
/*   IncidentalOp2009.cc                                                       */
/*   (c) 2010 Alan Lenarcic                                                   */
/*                                                                            */
/*   C code not used often, mostly to free memory and set matrix positions.   */
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

int ceiling (double Inpt);
int GetSign(double Value);
int SetMeI ( int kSLen, long double *GMat );
int SetMeXIn (int kLen, int NLen, int OnL, 
              int* ActiveOns, int *ActiveSigns, double *XInputMat , 
              long double *XOutputMat);
int SetMeXIn (int kLen, int NLen, int OnL, int* ActiveOns,  int*ActiveSigns,
              long double *XInputMat , 
              long double *XOutputMat);
              
#define RFP()  R_FlushConsole();  R_ProcessEvents()

void FFree(void **MyPointer);
void FFree(double **MyPointer);
void FFree(long double **MyPointer);
void FFree(int **MyPointer);
void FFree(short int **MyPointer);
int EqualToZero( long double A);              
 int EqualToZero( double A); 
int IntMin(int A1, int A2);                  


