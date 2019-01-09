/* ========================================================================== */
/*                                                                            */
/*   NEGLasso.c                                                               */
/*   (c) 2009-2019 Alan Lenarcic                                              */
/*       Work with Edoardo Airoldi Lab                                        */
/*                                                                            */
/*       Work accompanies effort by Lenarcic and Valdar on BayesSpike         */
/*       This is a companion algorithm using Coordinate Descent and EM        */
/*       to generate Model Inclusion Probability estimates and credibility    */
/*       without using Gibbs Sampling Integration.                            */
/*                                                                            */
/*   Code to implement NEG penalized selector                                 */
/*   This was a competing estimator was implemented after being strongly      */
/*   suggested as similar in approach to 2Lasso.              However,        */
/*  implementation required use of complex functions and did not seem to     */
/*  converge.                                                                 */
/*                                                                            */
/*   Description                                                              */
/*                                                                            */
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
//  On NEG Code:
//  You probably don't wnant to use this code, it was an experimental attempt
//  that probably did not choose the right implementation.
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
