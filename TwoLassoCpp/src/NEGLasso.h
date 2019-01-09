/* ========================================================================== */
/*                                                                            */
/*   NEGLasso.h                                                               */
/*   (c) 2009-2019 Alan Lenarcic                                              */
/*       Work with Edoardo Airoldi Lab                                        */
/*                                                                            */
/*       Work accompanies effort by Lenarcic and Valdar on BayesSpike         */
/*       This is a companion algorithm using Coordinate Descent and EM        */
/*       to generate Model Inclusion Probability estimates and credibility    */
/*       without using Gibbs Sampling Integration.                            */
/*                                                                            */
/*   Code to implement NEG penalized selector                                 */
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


double NEGLassoMakePSiValue( 
  SEXP PsiValues, SEXP logInvPsiNearZeroValues, 
  double BetaValue,
  double OnLambda, double OnGammaSq, 
  double *pIntegral, double *pInvSum, double *pSecondInvSum);
SEXP NEGLassoMakeTable(SEXP OutInvPsiValues, SEXP PsiValues,
  SEXP logInvPsiNearZeroValues, SEXP BetaValues,
  SEXP SOnLambda, SEXP SOnGammaSq, 
  SEXP SIntegral, SEXP SInvSum, SEXP SSecondInvSum, SEXP VerboseInt
  );
  

