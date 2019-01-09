/* ========================================================================== */
/*                                                                            */
/*   MyMatrixOp2009.h                                                         */
/*   (c) 2010 Alan Lenarcic                                                   */
/*                                                                            */
/*   C code not used often, mostly to print Matrix pieces and compare BLAS    */
/*   to less performative C implementations.                                  */
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
 
 ///////////////////////////////////////////////////////
 // MyMatrixOp.h
 //
 //   Self C coded Functions For Matrix Operation 
 //    Helpful for Running Matrices Etc.
 //
 //
 void PrintRpMatrix(double **pMat, int NR, int NC);
 void PrintRMatrix(double *Mat, int NR, int NC);
 void PrintRMatrix(long double *Mat, int NR, int NC);
 void PrintRMatrix(int *Mat, int NR, int NC);
 void PrintRMatrix(short int *Mat, int NR, int NC);
 void PrintRMatrixDenseAll(double *Mat, int NR, int NC);
 int InvertMatrix(int MatLeng, long double *InputMatrix, long double *OutputMatrix);
 int InvertMatrix(int MatLeng, double *InputMatrix, double *OutputMatrix);
 
 void PrintVector(long double *Vec, int Len);
 void PrintVector(double *Vec, int Len);
 void PrintVector(int *Vec, int Len);
 void PrintVector(short int *Vec, int Len);
 void PrintVectorAll(int *Vec, int Len);
 void PrintVectorAll(double *Vec, int Len);
 void PrintVectorAll(long double *Vec, int Len);
void PrintRMatrixElim(double *Mat, int NR, int NC, int ElimR, int ElimC);
void PrintVectorElim(double *Vec, int Len, int Elim);
void PrintVectorTimes(double *Vec, int Len, double Factor);


int MatTimesVec(int kOnLen, int NLen, double *XXA, double *InputV, 
       double *OutPutV); 
int MatTimesVec(int kOnLen, int NLen, long double *XXA, long double *InputV, long double *OutPutV);
int tMatTimesVec(int kOnLen, int NLen, double *XX, long double *InputV, long double *OutPutV);
int tMatTimesVec(int kOnLen, int NLen, double *XX, double *InputV, long double *OutPutV);
int tMatTimesVec(int kOnLen, int NLen, double *XX, double *InputV, double *OutPutV);
int tMatTimesVec(int kOnLen, int NLen,  long double *XX,  long double *InputV, double *OutPutV);
int tMatTimesVec(int kOnLen, int NLen,  long double *XX,  long double *InputV, long double *OutPutV);
int Cholesky(int Nlen, long double *InputMatrix, long double *OutputMatrix);

void SetToZero(int Nlen, long double *Matrix);
int MatTimesVecVec(int kOnLen, int NLen, long double *XXA, long double *WeightV, 
                  long double *InputV, long double *OutPutV);
int MatTimesVecVec(int kOnLen, int NLen, double *XXA, double *WeightV, 
                  double *InputV, double *WIV, double *OutPutV);                  
int SqMat (long double *RetMat, int n1, int n2, double *InpMat) ; 
int SqMat (double *RetMat, int n1, int n2, double *InpMat) ;
int pSqMat(double**RetMat, int n1, int n2, double *InpMat); 
int SqMat1 (long double *RetMat, int n1, int n2, double *InpMat);
int MatTMat (long double *RetMat, int n1, int n2, double *NMat, int m1, int m2, double * MMat); 
int ReSquareLU(double *RetMat, int kLen);
int ReSquareUL(double *RetMat, int kLen);

int CopyMatrix(int n1, int n2, long double *InputMatrix, long double *OutputMatrix);        
int ReduceFormatToOne(int KappaT, int Onk2, long double *OnDiag, long double *OffDiag, 
                    long double *OutW11);
   int OnToOff(int KappaT, int Onk2, double *InptSS,
                    long double *OnDiag, long double *OffDiag);                 

//// Lapack Based functions
int SolveSymPosMatrix(int dimk, double *Input, double *Output);
