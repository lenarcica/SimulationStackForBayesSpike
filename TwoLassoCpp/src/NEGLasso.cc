/* ========================================================================== */
/*                                                                            */
/*   NEGLasso.c                                                               */
/*   (c) 2011 Alan Lenarcic                                                          */
/*                                                                            */
/*   Description                                                              */
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

#ifndef NEGLASSOH
  #include "NEGLasso.h"
  #define NEGLASSOH 1
#endif


SEXP NEGLassoMakeTable(SEXP OutInvPsiValues, SEXP PsiValues,
  SEXP logInvPsiNearZeroValues, SEXP BetaValues,
  SEXP SOnLambda, SEXP SOnGammaSq, 
  SEXP SIntegral, SEXP SInvSum, SEXP SSecondInvSum, SEXP VerboseInt
  )
{
  //double dPsi = REAL(PsiValues)[1] - REAL(PsiValues)[2];
  
  int Betaii;
  int Count = 0;
  for (Betaii = 0; Betaii < length(BetaValues); Betaii++) {
    REAL(OutInvPsiValues)[Betaii] =  
      NEGLassoMakePSiValue(PsiValues, logInvPsiNearZeroValues,
        REAL(BetaValues)[Betaii],
        REAL(SOnLambda)[0], REAL(SOnGammaSq)[0],
        REAL(SIntegral)+ Betaii, REAL(SInvSum) + Betaii,
        REAL(SSecondInvSum) + Betaii);
    if (REAL(VerboseInt)[0] > 0 && Count == REAL(VerboseInt)[0]) {
      Rprintf("Finished %d, Beta = %f, OutInvPsi = %f\n",
        Betaii, REAL(BetaValues)[Betaii], 
        REAL(OutInvPsiValues)[Betaii]);
      R_FlushConsole(); Count = 0;
    }
    Count++;
  }
  return(OutInvPsiValues);
}


double NEGLassoMakePSiValue( 
  SEXP PsiValues, SEXP logInvPsiNearZeroValues, 
  double BetaValue,
  double OnLambda, double OnGammaSq, 
  double *pIntegral, double *pInvSum, double *pSecondInvSum) {

int lPsi = length(PsiValues) -1;
int ii;
double Integral = 0.0;
double InvSum = 0.0;
//const double NearZero = .0000001;
double BetaValueSq = BetaValue*BetaValue;
for (ii = 0; ii < lPsi; ii++) {
   Integral +=  exp(-.5 * BetaValueSq /REAL(PsiValues)[ii]) * 
     pow( REAL(PsiValues)[ii], -.5) * 
     pow( OnGammaSq + REAL(PsiValues)[ii], -OnLambda-1.0) *
     (REAL(PsiValues)[ii+1] - REAL(PsiValues)[ii-1]);
}  
Integral +=  exp( -.5 *  BetaValueSq / REAL(PsiValues)[0]) *
  pow(  REAL(PsiValues)[0], -.5) * 
  pow( OnGammaSq + REAL(PsiValues)[0], -OnLambda-1.0) *
  ( REAL(PsiValues)[1] - REAL(PsiValues)[0]);
Integral +=  exp( - .5 * BetaValueSq / REAL(PsiValues)[lPsi] ) *
  pow( REAL(PsiValues)[lPsi], -.5) * 
  pow( OnGammaSq + REAL(PsiValues)[lPsi], -OnLambda-1.0)  *
  ( REAL(PsiValues)[lPsi] - REAL(PsiValues)[lPsi-1]);
Integral = Integral * .5;
Integral += 2 * sqrt( REAL(PsiValues)[0]) * 
  exp( -.5 * BetaValueSq / REAL(PsiValues)[0]) *
  pow( OnGammaSq + REAL(PsiValues)[0], -OnLambda-1.0);

  
for (ii = 0; ii < lPsi; ii++) {
  InvSum += exp(-.5 * BetaValueSq /REAL(PsiValues)[ii])  *
    pow(REAL(PsiValues)[ii], -1.5) * 
    pow( OnGammaSq + REAL(PsiValues)[ii], -OnLambda-1.0) *
    (REAL(PsiValues)[ii+1] - REAL(PsiValues)[ii-1]);
}
InvSum +=  exp( -.5 * BetaValueSq / REAL(PsiValues)[0] ) *
  pow(  REAL(PsiValues)[0], -1.5) * 
  pow( OnGammaSq + REAL(PsiValues)[0], -OnLambda-1.0) *
  (REAL(PsiValues)[1] - REAL(PsiValues)[0]);
InvSum +=  exp( -.5 *  BetaValueSq / REAL(PsiValues)[lPsi]) *
  pow( REAL(PsiValues)[lPsi], -1.5) * 
  pow( OnGammaSq + REAL(PsiValues)[lPsi], -OnLambda-1.0) *
  (REAL(PsiValues)[lPsi] - REAL(PsiValues)[lPsi-1]);
InvSum = InvSum * .5;

pInvSum[0] = InvSum;
pIntegral[0] = Integral;
double SecondInvSum = 0;
lPsi = length(logInvPsiNearZeroValues)-1;
double logBetaValueSq = log(BetaValueSq);
for (ii = 1; ii < lPsi; ii++) {
  SecondInvSum += exp(-.5 * exp(REAL(logInvPsiNearZeroValues)[ii] + logBetaValueSq)
   + .5*REAL(logInvPsiNearZeroValues)[ii] )  *
   pow( OnGammaSq + (exp(-REAL(logInvPsiNearZeroValues)[ii])), -OnLambda-1.0);
  //Rprintf(" BetaValueSq = %f; OnLambda = %f; OnGammaSq = %f; IPNZ = %f; Out = %f \n",
  //  BetaValueSq, OnLambda, OnGammaSq, REAL(InvPsiNearZeroValues)[ii],
  //  exp(-.5 * REAL(InvPsiNearZeroValues)[ii] * BetaValueSq)  *
  //  pow(REAL(InvPsiNearZeroValues)[ii], -.5) * 
  //  pow( OnGammaSq + 1.0/REAL(InvPsiNearZeroValues)[ii], -OnLambda-1.0));
  //R_FlushConsole();
  //if (ii > 10) {error("Quit Now Temp");}
}
ii = 0;
  SecondInvSum+= .5 * exp(-.5 * exp(REAL(logInvPsiNearZeroValues)[ii] + logBetaValueSq)
    + REAL(logInvPsiNearZeroValues)[ii]* .5)  *
    pow( OnGammaSq + exp(-REAL(logInvPsiNearZeroValues)[ii]), -OnLambda-1.0);
  //Rprintf(" BetaValueSq = %f; OnLambda = %f; OnGammaSq = %f; IPNZ = %f; Out = %f \n",
  //  BetaValueSq, OnLambda, OnGammaSq, REAL(InvPsiNearZeroValues)[ii],
  //  exp(-.5 * REAL(InvPsiNearZeroValues)[ii] * BetaValueSq)  *
  //  pow(REAL(InvPsiNearZeroValues)[ii], -.5) * 
  //  pow( OnGammaSq + 1.0/REAL(InvPsiNearZeroValues)[ii], -OnLambda-1.0));
  //R_FlushConsole();
ii = lPsi;
  SecondInvSum+= .5 * exp(-.5 * exp(REAL(logInvPsiNearZeroValues)[lPsi] + BetaValueSq) +
    REAL(logInvPsiNearZeroValues)[lPsi] *.5) * 
    pow( OnGammaSq + exp(-REAL(logInvPsiNearZeroValues)[lPsi]), -OnLambda-1.0);    
  //Rprintf(" BetaValueSq = %f; OnLambda = %f; OnGammaSq = %f; IPNZ = %f; Out = %f \n",
  //  BetaValueSq, OnLambda, OnGammaSq, REAL(InvPsiNearZeroValues)[ii],
  //  exp(-.5 * REAL(InvPsiNearZeroValues)[ii] * BetaValueSq)  *
  //  pow( REAL(InvPsiNearZeroValues)[ii], -.5) * 
  //  pow( OnGammaSq + 1.0/(REAL(InvPsiNearZeroValues)[ii]), -OnLambda-1.0));
  //R_FlushConsole();
SecondInvSum *= ( REAL(logInvPsiNearZeroValues)[1] - REAL(logInvPsiNearZeroValues)[0]);
pSecondInvSum[0] = SecondInvSum;
return( (InvSum+SecondInvSum) / Integral);  
}
