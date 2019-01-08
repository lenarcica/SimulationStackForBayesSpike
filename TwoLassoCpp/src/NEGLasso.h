/* ========================================================================== */
/*                                                                            */
/*   NEGLasso.h                                                               */
/*   (c) 2011 Alan Lenarcic                                                   */
/*                                                                            */
/*   Code to implement NEG penalized selector                                 */
/*                                                                            */
/* ========================================================================== */



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
  

