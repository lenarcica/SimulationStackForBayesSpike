

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


