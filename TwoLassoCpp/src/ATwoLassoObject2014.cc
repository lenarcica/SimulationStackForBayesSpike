/* ========================================================================== */
/*                                                                            */
/*   ATwoLassoObject2014.cc                                                   */
/*   (c) 2009-2019 Alan Lenarcic                                              */
/*       Work with Edoardo Airoldi Lab                                        */
/*                                                                            */
/*       Work accompanies effort by Lenarcic and Valdar on BayesSpike         */
/*       This is a companion algorithm using Coordinate Descent and EM        */
/*       to generate Model Inclusion Probability estimates and credibility    */
/*       without using Gibbs Sampling Integration.                            */
/*                                                                            */
/*   TwoLasso code that takes advantage of SEXPS in RCpp format               */
/*                                                                            */
/*   This is the main ".cc" code for the main 2Lasso algorithm.               */
/*   Class Object "TwoLassoSexp" is declared and it is the easiest R console  */
/*   accesibile object using RCpp modules interface.                          */
/*   Calling "TwoLassoRegression()" should generate TwoLassoSexp and launch   */
/*   algorithm.                                                               */
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


#ifndef TWOLASSO2011
  #include "TwoLassoObject2014.h"
  #define TWOLASSO2011 0
#endif

#ifndef RDYNLOADH
  #include <R_ext/Rdynload.h>
  #define RDYNLOADH 0
#endif

  int SeekHPD(double *QuantilesList, int LenQuantile,  double Epsilon,
    double A, double B, double PiA, double LambdaA, double LambdaD,  double NegCons,
    int SkipGo,
    double *pOutLeft, double *pOutRight, double *PutMax);  

SEXP GiveUnpackedSEXP(double *SpM, int dM, int Verb, const char *NameMat) {
  SEXP sOut = R_NilValue;
  if (SpM == NULL) { return(R_NilValue);}
  Rf_protect(sOut = Rf_allocMatrix(REALSXP, dM, dM));
  if (Rf_isNull(sOut)) {
    Rf_error("Error Processing to print unpacked matrix %s\n", NameMat);
  }
  int i = 0; int j = 0; int Oo = 0; int Oij = 0;  int Oji = 0;
  for (j = 0; j < dM-1; j++) {
    Oji = j+j*dM;  Oij = Oji;
    REAL(sOut)[Oij] = SpM[Oo];
    Oij++;  Oo++;  Oji += dM;
    for (i = j +1; i < dM; i++) {
      REAL(sOut)[Oij] = SpM[Oo]; 
      REAL(sOut)[Oji] = SpM[Oo];
      Oo++;  Oji += dM; Oij++;
    }
  }
  REAL(sOut)[dM * dM -1] = SpM[ dM * (dM+1)/2 -1];
  Rf_unprotect(1);
  return(sOut);
}

void PrintPackedDouble(double *SpM, int dM, int Verb, const char *NameMat) {
  SEXP sOut = GiveUnpackedSEXP(SpM, dM, Verb, NameMat); 
  PrintRMatrix(REAL(sOut), dM, dM);
  return;
}

void DeleteSTLS(SEXP sTLS) {
 //Rprintf("TwoLassoSexp:: You Probably want to know, we're deleting TLS!");
 //R_FlushConsole();
 if (Rf_isNull(sTLS)) { return; }
 TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);


 if (sTLS != NULL) {
   //Rprtf("TwoLassoSexp: Setting PrintFlag to Delete\n"); R_FlushConsole();
   Rprintf("TwoLassoSexp: Deleting Through DeleteSTLS \n"); 
   TLS->Verbose = 4;
   //HM->PrintFlag = 4; 
   delete(TLS);
 }
 R_SetExternalPtrAddr(sTLS, NULL);
}

#ifndef HowToGet
  #define HowToGet(GetFunction, SGetFunction, SSGetFunction)    SEXP SSGetFunction(SEXP sTLS) { \
      if (Rf_isNull(sTLS)) { return(R_NilValue); }     \
      TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);  \
      if (TLS == NULL) {                 \
        Rprintf("HowToGet Error! TCS not cast to TLS\n");   \
      }                                  \
      return(TLS->GetFunction());           \
    }                         \
  extern "C" { \
    SEXP SGetFunction(SEXP sTLS) { \
      return(SSGetFunction(sTLS)); \
    }                                \
  }     
#endif
  HowToGet(getBetaCVInputs, SgetBetaCVInputs, SSgetBetaCVInputs)                                                 
  HowToGet(getSigmaVectorInputs, SgetSigmaVectorInputs, SSgetSigmaVectorInputs)    
  HowToGet( getPiAVectorInputs, SgetPiAVectorInputs, SSgetPiAVectorInputs)    

SEXP  sWhatIsSuccessFlag(SEXP sTLS) {
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no SuccessFlag!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("sTLS:WhatIsSuccessFlag: did not get TLS from TCSPointer, assume fail!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -1; Rf_unprotect(1); return(sOut);
  }
  INTEGER(sOut)[0] = TLS->get_SuccessFlag();
  Rf_unprotect(1);
  return(sOut);
}

SEXP sSetupConfidenceQuantiles(SEXP sTLS, SEXP iConfidenceQuantiles){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no SuccessFlag!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("sSetupConfidenceQuantiles: Error: TLS was NULL from TCS pointer!\n");
    R_FlushConsole(); INTEGER(sOut)[0] = -1;  Rf_unprotect(1); return(sOut);
  }
  TLS->set_ConfidenceQuantiles(iConfidenceQuantiles);
  INTEGER(sOut)[0] = 1;
  Rf_unprotect(1);
  return(sOut);
}


SEXP sgetHPDMatrix(SEXP sTLS){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no SuccessFlag!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("sgetHPDMatrix: Error, TLS from TCS is NUL!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -1;  Rf_unprotect(1); return(sOut);
  }
  SEXP sOut2 = R_NilValue;
  Rf_protect(sOut2 = TLS->get_HPDMatrix());
  if (Rf_isNull(sOut2)) {
    Rprintf("sgetHPDMatrix: get_ConfidenceMatrix returned a NULL!\n");
    R_FlushConsole(); Rf_unprotect(2); return(sOut);
  }
  Rf_unprotect(2);
  return(sOut2);
}
SEXP sgetConfidenceMatrix(SEXP sTLS){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no SuccessFlag!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("sgetConfidenceMatrix: Error, TLS from TCS is NUL!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -1;  Rf_unprotect(1); return(sOut);
  }
  SEXP sOut2 = R_NilValue;
  Rf_protect(sOut2 = TLS->get_ConfidenceMatrix());
  if (Rf_isNull(sOut2)) {
    Rprintf("sgetConfidenceMatrix: get_ConfidenceMatrix returned a NULL!\n");
    R_FlushConsole(); Rf_unprotect(2); return(sOut);
  }
  Rf_unprotect(2);
  return(sOut2);
}
SEXP sgetConfidenceQuantiles(SEXP sTLS){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no SuccessFlag!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("getConfidenceQuantiles: Error, TCS pointer not cast to TLS!\n"); R_FlushConsole();
    Rf_unprotect(1); return(sOut);
  }
  SEXP sOut2 = R_NilValue;
  Rf_protect(sOut2 = TLS->get_ConfidenceQuantiles());
  if (Rf_isNull(sOut2)) {
    Rprintf("sgetUnshrunkConfidenceMatrix: sOut2 is NULL return error!\n"); R_FlushConsole();
    Rf_unprotect(2); return(sOut);
  }
  Rf_unprotect(2);
  return(sOut2);

}

SEXP sgetUnshrunkConfidenceMatrix(SEXP sTLS){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  INTEGER(sOut)[0] = -666;
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no UnshrunkConfidenceMatrix!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("sgetUnshrunkConfidenceMatrix: failed to get TLS from pointer!\n");
    INTEGER(sOut)[0] = -666;
    Rf_unprotect(1);
    return(sOut);
  }
  SEXP sOut2 = R_NilValue;
  Rf_protect(sOut2 = TLS->get_UnshrunkConfidenceMatrix());
  if (Rf_isNull(sOut2)) {
    Rprintf("sgetUnshrunkConfidenceMatrix: sOut2 is NULL return error!\n"); R_FlushConsole();
    Rf_unprotect(2); return(sOut);
  }
  Rf_unprotect(2);
  return(sOut2);

}
SEXP sgetLambdaIndexForConfidenceIntervals(SEXP sTLS){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no LambdaIndexForConfidenceIntervals!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("sgetLambdaIndexForConfidenceIntervals: failed to get TLS from pointer!\n");
    INTEGER(sOut)[0] = -666; Rf_unprotect(1);
    return(sOut);
  }
  INTEGER(sOut)[0] = TLS->get_LambdaIndexForConfidenceIntervals();
  Rf_unprotect(1);
  return(sOut);
}
SEXP sgetVerbose(SEXP sTLS){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no sgetVerbose!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("sgetVerbose: failed to get TLS from pointer!\n");
    INTEGER(sOut)[0] = -666; Rf_unprotect(1);
    return(sOut);
  }
  INTEGER(sOut)[0] = TLS->get_Verbose();
  Rf_unprotect(1);
  return(sOut);
}


SEXP ssetLambdaIndexForConfidenceIntervals(SEXP sTLS, SEXP sIn){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  INTEGER(sOut)[0] = -666;
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no LambdaIndexForConfidenceIntervals!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  int iIn = -66;
  if (Rf_isInteger(sIn) && Rf_length(sIn) >= 1) { iIn = INTEGER(sIn)[0]; }
  if (Rf_isReal(sIn) && Rf_length(sIn) >= 1) { iIn = (int) REAL(sIn)[0]; }
  if (sIn < 0) {
    Rprintf("ssetLambdaIndexForConfidenceIntervals, sIn was not good!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;
    Rf_unprotect(1); return(sOut);
  }
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  INTEGER(sOut)[0] = -666;
  TLS->set_LambdaIndexForConfidenceIntervals(iIn);
  INTEGER(sOut)[0] = 1;
  Rf_unprotect(1);
  return(sOut);
}

SEXP ssetVerbose(SEXP sTLS, SEXP siVerbose){
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(INTSXP,1));
  INTEGER(sOut)[0] = -666;
  if (Rf_isNull(sTLS)) { 
    Rprintf("sTLS is Null there is no ssetVerbose!\n"); R_FlushConsole();
    INTEGER(sOut)[0] = -666;  Rf_unprotect(1);
    return(sOut);
  }
  int iVerbose = -66;
  if (Rf_isInteger(siVerbose) && Rf_length(siVerbose) >= 1) { iVerbose = INTEGER(siVerbose)[0]; }
  if (Rf_isReal(siVerbose) && Rf_length(siVerbose) >= 1) { iVerbose = (int) REAL(siVerbose)[0]; }

  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  INTEGER(sOut)[0] = -666;
  TLS->set_Verbose(iVerbose);
  INTEGER(sOut)[0] = 1;
  Rf_unprotect(1);
  return(sOut);
}



void SSetCVInputs(SEXP sTLS, SEXP rSPiAVectorInputs,  SEXP rSSigmaVectorInputs,
    SEXP rSLambdaAKInputs,  SEXP rSLambdaDKInputs, SEXP rSStartBetaMatrix,
    SEXP rSRecordBetaCV)  {
    if (Rf_isNull(sTLS)) { return; }
    TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
    TLS->SetupInputs(rSPiAVectorInputs, rSSigmaVectorInputs,
      rSLambdaAKInputs,  rSLambdaDKInputs, rSStartBetaMatrix,
      rSRecordBetaCV );
    return;
}
extern "C" {
  void SetCVInputs(SEXP sTLS, SEXP rSPiAVectorInputs,  SEXP rSSigmaVectorInputs,
  SEXP rSLambdaAKInputs,  SEXP rSLambdaDKInputs, SEXP rSStartBetaMatrix,
  SEXP rSRecordBetaCV) {
    SSetCVInputs(sTLS, rSPiAVectorInputs, rSSigmaVectorInputs,
      rSLambdaAKInputs, rSLambdaDKInputs, rSStartBetaMatrix,
      rSRecordBetaCV);
  }
}


SEXP TwoLassoCreateTwoLassoSexp(SEXP rSn, SEXP rSp, SEXP rSyys, 
    SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaA, SEXP rSLambdaD, SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
    SEXP rSVerbose,
    SEXP rGroupSexp, SEXP RunFlag) {
 int Verbose = 0;
 if (isInteger(rSVerbose)) {
   Verbose = INTEGER(rSVerbose)[0];
 }  else if (Rf_isReal(rSVerbose)) {
   Verbose = (int) REAL(rSVerbose)[0];
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoCreateTwoLassoSexp: Starting by loading object:\n"); R_FlushConsole();
 }
 TwoLassoSexp *TLS = new TwoLassoSexp( rSn, rSp, rSyys, 
     rSxxs,  rSXtY,
     rSXtX,  rStt1, rStt2, rSLambdaA,  rSLambdaD,  rSOrderSeq, 
     rSBBOn1,  rSOnBeta,  rRecBBOn1,  rRecOnBeta,
     rSOnGammas,  rSPiA,  rSSigma,
     rSm,  rSSigmaPrior,
     rSInitKKs,  rSInverseGammaConstant, 
     rSCauchyEpsilon,  rSMaxCauchy, 
     rSTDFNu,  rSiiWeights,  rSL2ShrinkagePrior, rSL2ShrinkageRecords,
     rSVerbose);
 TLS->BackXTYResid = (double *) Calloc(TLS->p, double); 
 if (Verbose > 0) {
   Rprintf("TwoLassoCreateTwoLassoSexp: Loaded Object:\n"); R_FlushConsole();
 }
 if (TLS == NULL || TLS->SuccessFlag < 0) {
   Rprintf("TwoLassoRegression: Sorry, some error in the Load");
   REAL(rSn)[0] = -999;
   R_FlushConsole();
   return(R_NilValue);
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoCreateTwoLassoSexp: Setup Groups?\n"); R_FlushConsole();
 }
 if (FALSE && !Rf_isNull(rGroupSexp)) {
   TLS->SetupGroups(rGroupSexp);
 }
 if (!Rf_isNull(RunFlag) && GetFirstInteger(RunFlag) >= 1) {
   if (Verbose > 0) {
    Rprintf("TwoLassoCreateTwoLassoSexp: Run regress\n"); R_FlushConsole();
   }
   TLS->RunTwoLassoRegression(1);  
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoCreateTwoLassoSexp: Make External\n"); R_FlushConsole();
 }
 if (TLS == NULL || TLS->SuccessFlag < 0) {
   Rprintf("TwoLassoRegression: Sorry, some error in the algorithm");
   REAL(rSn)[0] = -999;
   R_FlushConsole();
   return(R_NilValue);
 }
 SEXP sTLS =  R_MakeExternalPtr((void*) TLS, R_NilValue, R_NilValue);          
    
  // The Following Command registers to R garbage collector to delete the HM
  //   whenever it sees that references have disappeared.  
  // I don't think we need weak references here however
 R_RegisterCFinalizer(sTLS, DeleteSTLS);   
 return(sTLS);       
}

#define CheckCrossValidate(AThing, AText) \
  if (Rf_isNull(AThing)) {                                         \
    Rprintf("ERROR: CheckCrossValidate, we got %s as NULL Input!\n",\
      AText);  R_FlushConsole();                                    \
  }  else if (Rf_length(AThing) <= 0) { \
    Rprintf("ERROR: CheckCrossValidate, we got %s has 0 length!\n",  \
      AText); R_FlushConsole();                                     \
  }
  
void CheckTwoLassoCrossValidate(SEXP rSVerbose, SEXP rSn, 
    SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaAK, SEXP rSLambdaDK, 
    SEXP rSPiAVectorInputs, SEXP rSSigmaVectorInputs,
    SEXP rSLambdaAInputs, SEXP rSLambdaDInputs,   
    SEXP rSBetaStatStartMatrix,
    SEXP rBetaCVInputs,
    SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
    SEXP rGroupSexp) {
  CheckCrossValidate(rSp, "rSp")
  CheckCrossValidate(rSyys, "rSyys")
  CheckCrossValidate(rSxxs, "rSxxs")
  CheckCrossValidate(rSXtY, "rSXtY")
  CheckCrossValidate(rSXtX, "rSXtX")
  CheckCrossValidate(rStt1, "rStt1")
  CheckCrossValidate(rStt2, "rStt2")
  CheckCrossValidate(rSLambdaAK, "rSLambdaAK")    
  CheckCrossValidate(rSLambdaDK, "rSLambdaDK")
  CheckCrossValidate(rSOrderSeq, "rSOrderSeq")
  CheckCrossValidate(rSBBOn1, "rSBBOn1")
  CheckCrossValidate(rSOnBeta, "rSOnBeta")  
  CheckCrossValidate(rRecBBOn1, "rRecBBOn1")
  CheckCrossValidate(rRecOnBeta, "rRecOnBeta")
  CheckCrossValidate(rSOnGammas, "rSOnGammas")
  CheckCrossValidate(rSPiA, "rSPiA")  
  CheckCrossValidate(rSSigma, "rSSigma")
  CheckCrossValidate(rSm, "rSm")
  CheckCrossValidate(rSSigmaPrior, "rSSigmaPrior")
  CheckCrossValidate(rSInitKKs, "rSInitKKs")   
  CheckCrossValidate(rSInverseGammaConstant, "rSInverseGammaConstant")
  CheckCrossValidate(rSCauchyEpsilon, "rSCauchyEpsilon")
  CheckCrossValidate(rSMaxCauchy, "rSMaxCauchy")
  CheckCrossValidate(rSTDFNu, "rSTDFNu") 
  CheckCrossValidate(rSiiWeights, "rSiiWeights")
  CheckCrossValidate(rSMaxCauchy, "rSMaxCauchy")
  CheckCrossValidate(rSL2ShrinkagePrior, "rSL2ShrinkagePrior")     
  CheckCrossValidate(rSL2ShrinkageRecords, "rSL2ShrinkageRecords")   
  CheckCrossValidate(rSVerbose, "rSVerbose")  
  CheckCrossValidate(rGroupSexp, "rGroupSexp")  
  CheckCrossValidate(rSPiAVectorInputs, "rSPiAVectorInputs")
  CheckCrossValidate(rSSigmaVectorInputs, "rSSigmaVectorInputs")     
  CheckCrossValidate(rSLambdaAInputs, "rSLambdaAInputs")   
  CheckCrossValidate(rSLambdaDInputs, "rSLambdaDInputs")  
  CheckCrossValidate(rSBetaStatStartMatrix, "rSBetaStatStartMatrix")   
  CheckCrossValidate(rBetaCVInputs, "rBetaCVInputs")  
  
  if ( !Rf_isNull(rSVerbose) && (Rf_isInteger(rSVerbose) ||
    Rf_isReal(rSVerbose)) && GetFirstInteger(rSVerbose) >= 1) {
    if (Rf_isNull(rSLambdaDInputs) || rSLambdaDInputs == NULL) {
      Rprintf("Weird, rSLambdaDInputs is NULL!\n");  R_FlushConsole();
    } else {
      Rprintf("rSLambdaDInputs is not NULL, that's good, length=%d\n",
        Rf_length(rSLambdaDInputs));
      R_FlushConsole();
    }
    //if (Rf_isNull(rSLambdaAInputs) || rSLambdaAInputs == NULL) {
    //  Rprintf("Weird, rSLambdaAInputs is NULL!\n");  R_FlushConsole();
    //} else {
    //  Rprintf("rSLambdaAInputs is not NULL, that's good, length=%d\n",
    //    Rf_length(rSLambdaAInputs));
    //  R_FlushConsole();
    //}
  }

  return;
}
SEXP TwoLassoCrossValidate(SEXP rSVerbose, SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaAK, SEXP rSLambdaDK, 
    SEXP rSPiAVectorInputs, SEXP rSSigmaVectorInputs,
    SEXP rSLambdaAInputs, SEXP rSLambdaDInputs,   
    SEXP rSBetaStatStartMatrix,
    SEXP rBetaCVInputs,
    SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
    SEXP rGroupSexp) {
 int Verbose = 0;
 if (Rf_isNull(rSVerbose)) {
   Rprintf("TwoLassoCrossValidate: Issue, rSVerbose is NULL!\n"); 
   R_FlushConsole();
   Verbose = 10;
 } else if (isInteger(rSVerbose)) {
   Verbose = INTEGER(rSVerbose)[0];
 }  else if (Rf_isReal(rSVerbose)) {
   Verbose = (int) REAL(rSVerbose)[0];
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoCreateTwoLassoSexp: Starting by loading object:\n"); R_FlushConsole();
 }
 if (Verbose > 4) {
   if (Rf_isNull(rSL2ShrinkagePrior)) {
     Rprintf("L2ShrinkagePrior is NULL!\n"); R_FlushConsole();
   } else {
     Rprintf("Length L2ShrinakgePrior = %d. \n", 
       Rf_length(rSL2ShrinkagePrior) );
   }
   if (Rf_isNull(rSL2ShrinkageRecords)) {
     Rprintf("rSL2ShrinkageRecords is NULL!\n"); R_FlushConsole();
   } else {
     Rprintf("Length rSL2ShrinkageRecords = %d. \n",
       Rf_length(rSL2ShrinkageRecords));
   }

   R_FlushConsole();
 }
 
 if (Verbose >= 1) {
   Rprintf("We are running CheckTwoLassoCrossValidate! \n"); R_FlushConsole();
 }
 /*
 CheckTwoLassoCrossValidate(rSVerbose, rSn, 
    rSp, rSyys, rSxxs, rSXtY,
    rSXtX, rStt1, rStt2,
    rSLambdaAK, rSLambdaDK, 
    rSPiAVectorInputs, rSSigmaVectorInputs,
    rSLambdaAInputs, rSLambdaDInputs,   
    rSBetaStatStartMatrix,
    rBetaCVInputs,
    rSOrderSeq, 
    rSBBOn1, rSOnBeta, rRecBBOn1, rRecOnBeta,
    rSOnGammas, rSPiA, rSSigma,
    rSm, rSSigmaPrior,
    rSInitKKs, rSInverseGammaConstant, 
    rSCauchyEpsilon, rSMaxCauchy, 
    rSTDFNu, rSiiWeights, 
    rSL2ShrinkagePrior, rSL2ShrinkageRecords,
    rGroupSexp);
 */
 if (Rf_isNull(rSn)) {
   Rprintf("TwoLassoCreateTwoLassoSexp:: Error rSn is NULL!"); R_FlushConsole();
 }
 TwoLassoSexp *TLS = new TwoLassoSexp( rSn, rSp, rSyys, 
     rSxxs,  rSXtY,
     rSXtX,  rStt1, rStt2, rSLambdaAK,  rSLambdaDK,  rSOrderSeq, 
     rSBBOn1,  rSOnBeta,  rRecBBOn1,  rRecOnBeta,
     rSOnGammas,  rSPiA,  rSSigma,
     rSm,  rSSigmaPrior,
     rSInitKKs,  rSInverseGammaConstant, 
     rSCauchyEpsilon,  rSMaxCauchy, 
     rSTDFNu,  rSiiWeights,  
     rSL2ShrinkagePrior, rSL2ShrinkageRecords, rSVerbose);
 if (Verbose > 2) {
   if (Rf_isNull(RGET_DIM(rBetaCVInputs))) {
     Rf_error("Before SetupInput: rBetaCVInputs has 0 dim!\n");
   }
   Rprintf("Before SetupInputs, we got a dim of rBetaCVInputs = %d, %d\n",
     INTEGER(RGET_DIM(rBetaCVInputs))[0],
      INTEGER(RGET_DIM(rBetaCVInputs))[1]); R_FlushConsole();
 }
  
 //if (rSLambdaDInputs == NULL || Rf_isNull(rSLambdaDInputs)) {
 //  Rprintf("Before Setup Inputs, OHNO, SLambdaDInputs is NULL!\n"); R_FlushConsole();
 //} else {
 //  Rprintf(" WEIRD, SLambdaDInputs is not null, length = %d, right before Setup!\n",
 //    Rf_length(rSLambdaDInputs)); R_FlushConsole();
 //}     
 TLS->SetupInputs(rSPiAVectorInputs, rSSigmaVectorInputs,
   rSLambdaAInputs, rSLambdaDInputs, rSBetaStatStartMatrix, rBetaCVInputs);
 if (GetFirstInteger(rSPiAVectorInputs) < 0) {
   TLS->FixKa = abs((int) GetFirstInteger(rSPiAVectorInputs));
 }
 if (Verbose >= 1) {
   Rprintf("TwoLassoCreateTwoLassoSexp: Allocating backup items Rf_length %d\n",
     TLS->p); R_FlushConsole();
 }
 if (Verbose >= 1) {
   Rprintf("TwoLassoCreateTwoLassoSexp: Setup Groups?\n"); R_FlushConsole();
 }
 if (FALSE && !Rf_isNull(rGroupSexp)) {
   TLS->SetupGroups(rGroupSexp);
 }
 

 if (TLS == NULL || TLS->SuccessFlag < 0) {
   Rprintf("TwoLassoRegression: (ATwoLassoObject20011.cc) Sorry, some error in the Load");
   R_FlushConsole();
   REAL(rSn)[0] = -999;
   R_FlushConsole();
   return(R_NilValue);
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoCrossValidate: (ATwoLassoObject20011.cc): Loaded Object.\n"); R_FlushConsole();
 }
 
 if (Verbose > 0 || TLS->p >= 50000) {
   Rprintf("TwoLassoCrossValidate: (ATwoLassoObject20011.cc): Time to run TLS->RunCrossValidate()\n"); R_FlushConsole();
 }
 TLS->RunCrossValidate();
 
 if (TLS->SuccessFlag < 0) {
   Rprintf("TwoLassoCrossValidate: we got a failed Success flag somehow! \n"); R_FlushConsole();
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoCrossValidate: Making External Pointer for export"); R_FlushConsole();
 }
 if (TLS == NULL || TLS->SuccessFlag < 0) {
   Rprintf("TwoLassoRegression: Sorry, some error in the algorithm");
   REAL(rSn)[0] = -999;
   R_FlushConsole();
   return(R_NilValue);
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoCrossValidate: Making External Pointer! \n"); R_FlushConsole();
 }
 SEXP sTLS =  R_MakeExternalPtr((void*) TLS, R_NilValue, R_NilValue);          
    
  // The Following Command registers to R garbage collector to delete the HM
  //   whenever it sees that references have disappeared.  
  // I don't think we need weak references here however
 if (Verbose > 0) {
   Rprintf("TwoLassoCrossValidate: Registering Finalizer \n"); R_FlushConsole();
 }
 R_RegisterCFinalizer(sTLS, DeleteSTLS);   
 if (Verbose > 0) {
   Rprintf("TwoLassoCrossValidate: Returnb sTLS\n"); R_FlushConsole();
 }
 return(sTLS);       
}

int TwoLassoSexp::VerbLessTestCrossValidate() {
  double CBeta; int ACount = 0;
  if (BackOrigBeta == NULL) {
    Rprintf("\n*** TwoLassoCpp::TestCrossValidate, Error, BackOrigBeta is NULL!  ACount = %d\n", ACount); R_FlushConsole();
    Rf_error("TwoLassoCpp:: BackOrigBeta NULL");
  }
  CBeta = 0.0;
  int ii = 0;
  for (ii = 0; ii < p; ii++) {
    CBeta += this->BackOrigBeta[ii];
  }
  if (this->SXtY == NULL) {
     Rprintf("*** SXtY is actually a NULL value\n");   R_FlushConsole();
  } else if (Rf_isNull(this->SXtY)) {
    Rprintf("*** SXtY is A NIL Object!\n"); R_FlushConsole();
  } else if (!Rf_isReal(this->SXtY)) {
    Rprintf("*** SXtY is Not Real!\n"); R_FlushConsole();
  } else if (!Rf_isNull(this->SXtY) && Rf_length(this->SXtY) != p) {
     Rprintf("*** SXtY not a Nil Value, but length SXtY not length p! is %d\n", Rf_length(this->SXtY)); R_FlushConsole();
  } else if (Rf_length(this->SXtY) == p) {
  } else {
     Rprintf("*** SXtY is a Not good at all!\n");     R_FlushConsole();
  }                      
  if (this->SXtY == NULL || Rf_isNull(this->SXtY) || Rf_length(this->SXtY) != p ||
   !Rf_isReal(this->SXtY)) {
    if (CDO->XTY == NULL) {
      Rf_error("*** Error SXtY Fails!\n"); R_FlushConsole();
    }
  } else {
  CBeta = 0.0;
  for (ii = 0; ii < p; ii++) {
    CBeta += REAL(this->SXtY)[ii];
  }
  }

  if (p != CDO->kLen) {
    Rf_error("*** TestCrossValidate, no way, p = %d, CDO->kLen=%d!\n", p, CDO->kLen);
  }
  if (CDO == NULL) {
    Rprintf("***   Woah CDO is Not Allocated!\n"); R_FlushConsole();
    Rf_error("*** TestCrossValidate: CDO is not allocated!\n");
  }
  if (CDO->XTY == NULL) {
    Rprintf("***  Woah, CDO->XTY is NULL!\n"); R_FlushConsole();
    Rf_error("  CDO->XTY not Allocated!\n");
  }
  CBeta =0.0;
  for (ii = 0; ii < p; ii++) {
    CBeta += CDO->XTY[ii];
    CDO->XTY[ii] = CDO->XTY[ii];
  }
  if (this->BackXTYResid == NULL) {
    Rprintf(" -- Fail, BackXTYResid was not Allocated!\n"); R_FlushConsole();
    Rf_error(" -- Fail BackXTYResid was not allocated. \n");
  }
  CBeta = 0.0;
  for (ii = 0; ii < p; ii++) { CBeta += this->BackXTYResid[ii];
    this->BackXTYResid[ii] = this->BackXTYResid[ii]; 
  }
  if (this->BackOrigBeta == NULL) {
    Rprintf(" -- It's NULL!\n"); R_FlushConsole(); 
  } else {
    Rprintf(" -- It's not NULL -- "); R_FlushConsole();
    CBeta = 0.0;
    for (ii = 0; ii < p; ii++ ) {
      CBeta += this->BackOrigBeta[ii];
      this->BackOrigBeta[ii] = this->BackOrigBeta[ii];
    }
    Rprintf("  Pass with total = %f \n", CBeta); R_FlushConsole();
  }        

  CBeta = 0.0;
  double MinPiA = REAL(SPiAVectorInputs)[0];
  double MaxPiA = REAL(SPiAVectorInputs)[0];
  for (ii = 0; ii < Rf_length(SPiAVectorInputs); ii++) {
    CBeta += REAL(SPiAVectorInputs)[ii];
    if (R_isnancpp(REAL(SPiAVectorInputs)[ii])) {
       Rf_error("Error TestCrossValidate, SPiA Input %d is NAN!\n",
         ii);
    }
    if (REAL(SPiAVectorInputs)[ii] < 0 || REAL(SPiAVectorInputs)[ii] > 1.0) {
       Rf_error("Error TestCrossValidate, SPiA Input %d is %f!\n",
         ii, REAL(SPiAVectorInputs)[ii]);
    }
    if (MinPiA > REAL(SPiAVectorInputs)[ii]) {
       MinPiA = REAL(SPiAVectorInputs)[ii];
    }
    if (MaxPiA < REAL(SPiAVectorInputs)[ii]) {
       MaxPiA = REAL(SPiAVectorInputs)[ii];
    }    
  }
  
  double MinSigma = REAL(SSigmaVectorInputs)[0];
  double MaxSigma = REAL(SSigmaVectorInputs)[0];
  for (ii = 0; ii < Rf_length(SSigmaVectorInputs); ii++) {
    CBeta += REAL(SSigmaVectorInputs)[ii];
    if (R_isnancpp(REAL(SSigmaVectorInputs)[ii])) {
       Rf_error("Error TestCrossValidate, Sigma Input %d is NAN!\n",
         ii);
    }
    if (REAL(SSigmaVectorInputs)[ii] < 0 || !R_finite(REAL(SSigmaVectorInputs)[ii])) {
       Rf_error("Error TestCrossValidate, Sigma Input %d is %f!\n",
         ii, REAL(SSigmaVectorInputs)[ii]);
    }
    if (MinSigma > REAL(SSigmaVectorInputs)[ii]) {
       MinSigma = REAL(SSigmaVectorInputs)[ii];
    }
    if (MaxSigma < REAL(SSigmaVectorInputs)[ii]) {
       MaxSigma = REAL(SSigmaVectorInputs)[ii];
    }    
  }

      
  double MinLambdaA = REAL(SLambdaAKInputs)[0];
  double MaxLambdaA = REAL(SLambdaAKInputs)[0];

  for (ii = 0; ii < Rf_length(SLambdaAKInputs); ii++) {
    CBeta += REAL(SLambdaAKInputs)[ii];
    if (R_isnancpp(REAL(SLambdaAKInputs)[ii])) {
       Rf_error("Error TestCrossValidate, LambdaA Input %d is NAN!\n",
         ii);
    }
    if (REAL(SLambdaAKInputs)[ii] < 0 || !R_finite(REAL(SLambdaAKInputs)[ii])) {
       Rf_error("Error TestCrossValidate, LambdaA Input %d is %f!\n",
         ii, REAL(SLambdaAKInputs)[ii]);
    }
    if (MinLambdaA > REAL(SLambdaAKInputs)[ii]) {
       MinLambdaA = REAL(SLambdaAKInputs)[ii];
    }
    if (MaxLambdaA < REAL(SLambdaAKInputs)[ii]) {
       MaxLambdaA = REAL(SLambdaAKInputs)[ii];
    }    
  }
    
  double MinLambdaD = REAL(SLambdaDKInputs)[0];
  double MaxLambdaD = REAL(SLambdaDKInputs)[0];
  for (ii = 0; ii < Rf_length(SLambdaDKInputs); ii++) {
    CBeta += REAL(SLambdaDKInputs)[ii];
    if (R_isnancpp(REAL(SLambdaDKInputs)[ii])) {
       Rf_error("Error TestCrossValidate, LambdaD Input %d is NAN!\n",
         ii);
    }
    if (REAL(SLambdaDKInputs)[ii] < 0 || !R_finite(REAL(SLambdaDKInputs)[ii])) {
       Rf_error("Error TestCrossValidate, LambdaD Input %d is %f!\n",
         ii, REAL(SLambdaDKInputs)[ii]);
    }
    if (MinLambdaD > REAL(SLambdaDKInputs)[ii]) {
       MinLambdaD = REAL(SLambdaDKInputs)[ii];
    }
    if (MaxLambdaD < REAL(SLambdaDKInputs)[ii]) {
       MaxLambdaD = REAL(SLambdaDKInputs)[ii];
    }    
  }
//  CDO->TestCDO();
  return(1);
}

int TwoLassoSexp::TestCrossValidate() {
  Rprintf("*** TwoLassoCpp::TestCrossValidate(n=%d,p=%d): Start. \n", n, p); R_FlushConsole();
  double CBeta; int ACount = 0;
  Rprintf("***   -- Test BackOrigBeta "); R_FlushConsole();
  if (BackOrigBeta == NULL) {
    Rprintf("\n*** TwoLassoCpp::TestCrossValidate, Error, BackOrigBeta is NULL!  ACount = %d\n", ACount); R_FlushConsole();
    Rf_error("TwoLassoCpp:: BackOrigBeta NULL");
  }
  CBeta = 0.0;
  int ii = 0;
  for (ii = 0; ii < p; ii++) {
    CBeta += this->BackOrigBeta[ii];
  }
  Rprintf("  -- PASS CBeta = %f. \n", CBeta); R_FlushConsole();
  Rprintf("***        Test SXtY. \n"); R_FlushConsole();
  if (this->SXtY == NULL) {
     Rprintf("*** SXtY is actually a NULL value\n");   R_FlushConsole();
  } else if (Rf_isNull(this->SXtY)) {
    Rprintf("*** SXtY is A NIL Object!\n"); R_FlushConsole();
  } else if (!Rf_isReal(this->SXtY)) {
    Rprintf("*** SXtY is Not Real!\n"); R_FlushConsole();
  } else if (!Rf_isNull(this->SXtY) && Rf_length(this->SXtY) != p) {
     Rprintf("*** SXtY not a Nil Value, but length SXtY not length p! is %d\n", Rf_length(this->SXtY)); R_FlushConsole();
  } else if (Rf_length(this->SXtY) == p) {
  } else {
     Rprintf("*** SXtY is a Not good at all!\n");     R_FlushConsole();
  }                      
  if (this->SXtY == NULL || Rf_isNull(this->SXtY) || Rf_length(this->SXtY) != p ||
   !Rf_isReal(this->SXtY)) {
    if (CDO->XTY == NULL) {
      Rf_error("*** Error SXtY Fails!\n"); R_FlushConsole();
    }
  } else {
  CBeta = 0.0;
  for (ii = 0; ii < p; ii++) {
    CBeta += REAL(this->SXtY)[ii];
  }
  Rprintf(" -- PASS, SXtY total = %f\n", CBeta); R_FlushConsole();
  }
  if (CDO == NULL) {
    Rprintf("***   Woah CDO is Not Allocated!\n"); R_FlushConsole();
    Rf_error("*** TestCrossValidate: CDO is not allocated!\n");
  }
  Rprintf(" -- Test CDO  \n"); R_FlushConsole();
  Rprintf(" -- Test CDO->XTY : "); R_FlushConsole();
  if (CDO->XTY == NULL) {
    Rprintf("***  Woah, CDO->XTY is NULL!\n"); R_FlushConsole();
    Rf_error("  CDO->XTY not Allocated!\n");
  }
  CBeta =0.0;
  for (ii = 0; ii < p; ii++) {
    CBeta += CDO->XTY[ii];
    CDO->XTY[ii] = CDO->XTY[ii];
  }
  Rprintf(" -- CDO->XTY Pass, total = %f\n", CBeta); R_FlushConsole(); 
  Rprintf("*** Test this->BackXTYResid ");
  if (SPiA != NULL && !Rf_isNull(SPiA)) {
    Rprintf(" --  Test SPiA -- ");  R_FlushConsole();
    CBeta = 0.0;  
    for (ii = 0; ii < Rf_length(SPiA); ii++) {
      CBeta += REAL(SPiA)[ii] / Rf_length(SPiA);
      REAL(SPiA)[ii] = REAL(SPiA)[ii];
    }
    Rprintf(" -- PASS Av = %f\n", CBeta); R_FlushConsole();
  }
  if (SSigma != NULL && !Rf_isNull(SSigma)) {
    Rprintf(" --  Test SSigma -- ");  R_FlushConsole();
    CBeta = 0.0;  
    for (ii = 0; ii < Rf_length(SSigma); ii++) {
      CBeta += REAL(SSigma)[ii] / Rf_length(SSigma);
      REAL(SSigma)[ii] = REAL(SSigma)[ii];
    }
    Rprintf(" -- PASS Av = %f \n", CBeta); R_FlushConsole();
  }
  if (this->BackXTYResid == NULL) {
    Rprintf(" -- Fail, BackXTYResid was not Allocated!\n"); R_FlushConsole();
    Rf_error(" -- Fail BackXTYResid was not allocated. \n");
  }
  CBeta = 0.0;
  for (ii = 0; ii < p; ii++) { CBeta += this->BackXTYResid[ii];
    this->BackXTYResid[ii] = this->BackXTYResid[ii]; 
  }
  Rprintf(" -- PASS, CBeta = %f \n", CBeta); R_FlushConsole();
  Rprintf("*** Test BackOrigBeta  -- ");
  if (this->BackOrigBeta == NULL) {
    Rprintf(" -- It's NULL!\n"); R_FlushConsole(); 
  } else {
    Rprintf(" -- It's not NULL -- "); R_FlushConsole();
    CBeta = 0.0;
    for (ii = 0; ii < p; ii++ ) {
      CBeta += this->BackOrigBeta[ii];
      this->BackOrigBeta[ii] = this->BackOrigBeta[ii];
    }
    Rprintf("  Pass with total = %f \n", CBeta); R_FlushConsole();
  }        
  Rprintf("*** CDO is about to be tested \n"); R_FlushConsole();
  CDO->TestCDO();
  Rprintf("*** CDO Has been Tested \n"); R_FlushConsole();
  Rprintf("*** Length SPiAVectorInputs is %d \n", Rf_length(SPiAVectorInputs)); R_FlushConsole();
  CBeta = 0.0;
  double MinPiA = REAL(SPiAVectorInputs)[0];
  double MaxPiA = REAL(SPiAVectorInputs)[0];
  for (ii = 0; ii < Rf_length(SPiAVectorInputs); ii++) {
    CBeta += REAL(SPiAVectorInputs)[ii];
    if (R_isnancpp(REAL(SPiAVectorInputs)[ii])) {
       Rf_error("Error TestCrossValidate, SPiA Input %d is NAN!\n",
         ii);
    }
    if (REAL(SPiAVectorInputs)[ii] < 0 || REAL(SPiAVectorInputs)[ii] > 1.0) {
       Rf_error("Error TestCrossValidate, SPiA Input %d is %f!\n",
         ii, REAL(SPiAVectorInputs)[ii]);
    }
    if (MinPiA > REAL(SPiAVectorInputs)[ii]) {
       MinPiA = REAL(SPiAVectorInputs)[ii];
    }
    if (MaxPiA < REAL(SPiAVectorInputs)[ii]) {
       MaxPiA = REAL(SPiAVectorInputs)[ii];
    }    
  }
  Rprintf("*** TestCrossValidate, PiA ranges from (%.9f, %.6f) \n",
    MinPiA, MaxPiA); R_FlushConsole();
  
  double MinSigma = REAL(SSigmaVectorInputs)[0];
  double MaxSigma = REAL(SSigmaVectorInputs)[0];
  Rprintf("*** TestCrossValidate, testing Sigma \n");
  for (ii = 0; ii < Rf_length(SSigmaVectorInputs); ii++) {
    CBeta += REAL(SSigmaVectorInputs)[ii];
    if (R_isnancpp(REAL(SSigmaVectorInputs)[ii])) {
       Rf_error("Error TestCrossValidate, Sigma Input %d is NAN!\n",
         ii);
    }
    if (REAL(SSigmaVectorInputs)[ii] < 0 || !R_finite(REAL(SSigmaVectorInputs)[ii])) {
       Rf_error("Error TestCrossValidate, Sigma Input %d is %f!\n",
         ii, REAL(SSigmaVectorInputs)[ii]);
    }
    if (MinSigma > REAL(SSigmaVectorInputs)[ii]) {
       MinSigma = REAL(SSigmaVectorInputs)[ii];
    }
    if (MaxSigma < REAL(SSigmaVectorInputs)[ii]) {
       MaxSigma = REAL(SSigmaVectorInputs)[ii];
    }    
  }
  Rprintf("*** TestCrossValidate Sigma ranges from (%f, %f), CBeta = %f\n",
    MinSigma, MaxSigma, CBeta); R_FlushConsole();
      
  double MinLambdaA = REAL(SLambdaAKInputs)[0];
  double MaxLambdaA = REAL(SLambdaAKInputs)[0];
  Rprintf("*** TestCrossValidate, testing LambdaA, length %d\n", 
    Rf_length(SLambdaAKInputs));
  for (ii = 0; ii < Rf_length(SLambdaAKInputs); ii++) {
    CBeta += REAL(SLambdaAKInputs)[ii];
    if (R_isnancpp(REAL(SLambdaAKInputs)[ii])) {
       Rf_error("Error TestCrossValidate, LambdaA Input %d is NAN!\n",
         ii);
    }
    if (REAL(SLambdaAKInputs)[ii] < 0 || !R_finite(REAL(SLambdaAKInputs)[ii])) {
       Rf_error("Error TestCrossValidate, LambdaA Input %d is %f!\n",
         ii, REAL(SLambdaAKInputs)[ii]);
    }
    if (MinLambdaA > REAL(SLambdaAKInputs)[ii]) {
       MinLambdaA = REAL(SLambdaAKInputs)[ii];
    }
    if (MaxLambdaA < REAL(SLambdaAKInputs)[ii]) {
       MaxLambdaA = REAL(SLambdaAKInputs)[ii];
    }    
  }
  Rprintf("*** TestCrossValidate LambdaA ranges from (%f, %f), CBeta = %f\n",
    MinLambdaA, MaxLambdaA, CBeta); R_FlushConsole();
    
  double MinLambdaD = REAL(SLambdaDKInputs)[0];
  double MaxLambdaD = REAL(SLambdaDKInputs)[0];
  Rprintf("*** TestCrossValidate, testing LambdaD, length %d\n", 
    Rf_length(SLambdaDKInputs));
  for (ii = 0; ii < Rf_length(SLambdaDKInputs); ii++) {
    CBeta += REAL(SLambdaDKInputs)[ii];
    if (R_isnancpp(REAL(SLambdaDKInputs)[ii])) {
       Rf_error("Error TestCrossValidate, LambdaD Input %d is NAN!\n",
         ii);
    }
    if (REAL(SLambdaDKInputs)[ii] < 0 || !R_finite(REAL(SLambdaDKInputs)[ii])) {
       Rf_error("Error TestCrossValidate, LambdaD Input %d is %f!\n",
         ii, REAL(SLambdaDKInputs)[ii]);
    }
    if (MinLambdaD > REAL(SLambdaDKInputs)[ii]) {
       MinLambdaD = REAL(SLambdaDKInputs)[ii];
    }
    if (MaxLambdaD < REAL(SLambdaDKInputs)[ii]) {
       MaxLambdaD = REAL(SLambdaDKInputs)[ii];
    }    
  }
  Rprintf("*** TestCrossValidate LambdaD ranges from (%f, %f), CBeta = %f\n",
    MinLambdaD, MaxLambdaD, CBeta); R_FlushConsole();
  Rprintf("*** Test RSRecordBetaCV  ");
  if (RSRecordBetaCV == NULL || Rf_isNull(RSRecordBetaCV->asSexp())) {
    Rprintf(" -- FAIL RSRecordBetaCV is NULL!\n "); R_FlushConsole();
  } else if (Rf_length(RSRecordBetaCV->asSexp()) <= 0 ) {
    Rprintf(" -- FAIL RSRecordBetaCV is Zero Length!\n "); R_FlushConsole();
  } else {
    CBeta = 0.0;
    for (ii = 0; ii < Rf_length(RSRecordBetaCV->asSexp()); ii++) {
      CBeta += REAL(RSRecordBetaCV->asSexp())[ii] / Rf_length(RSRecordBetaCV->asSexp());      
    }
    Rprintf(" -- Average is %f for length %d \n", CBeta,
        Rf_length(RSRecordBetaCV->asSexp())); R_FlushConsole();
  }

  Rprintf("*** Test RSRecordBetaCV  ");
  if (RSRecordBetaCV == NULL || Rf_isNull(RSRecordBetaCV->asSexp())) {
    Rprintf(" -- FAIL RSRecordBetaCV is NULL!\n "); R_FlushConsole();
  } else if (Rf_length(RSRecordBetaCV->asSexp()) <= 0 ) {
    Rprintf(" -- FAIL RSRecordBetaCV is Zero Length!\n "); R_FlushConsole();
  } else {
    CBeta = 0.0;
    for (ii = 0; ii < Rf_length(RSRecordBetaCV->asSexp()); ii++) {
      CBeta += REAL(RSRecordBetaCV->asSexp())[ii] / Rf_length(RSRecordBetaCV->asSexp());

    }
    Rprintf(" -- Average is %f for length %d \n", CBeta,
      Rf_length(RSRecordBetaCV->asSexp())); R_FlushConsole();
  }  
  
  Rprintf("*** TestCrossValidate, Now Test CDO. \n"); R_FlushConsole();
  CDO->TestCDO();
  Rprintf("*** ATwoLassoObject2011.cc: TestCrossValidate is concluding. \n"); R_FlushConsole();   
  Rprintf("*** TestCrossValidate concluded. \n"); R_FlushConsole();
  return(1);
}
int TwoLassoSexp::RunCrossValidate() {
 int One = 1; double ZeroD = 0.0; 
 if (Verbose >= 1 || p >= LARGEDATA) {
   Rprintf("**************************************************************\n");
   Rprintf("*** In ATwoLassoObject20011.cc \n");
   Rprintf("*** TwoLassoCpp:RunCrossValidate(n=%d,p=%d) Starting.",  
     this->n, this->p); R_FlushConsole();                                     
   Rprintf("***  We have length SPiAVectorInputs = %d.\n", Rf_length(SPiAVectorInputs));
 }
 if (this->BackOrigBeta == NULL) {
   Rf_error("RunCrossValidate: BackOrigBeta is NULL!\n");
 }
 if (this->p >= LARGEDATA) {
   Rprintf("*** TwoLassoCpp::RunCrossValidate(n=%d,p=%d), Run Initial TestCrossValidate.\n", n,p); R_FlushConsole();
   TestCrossValidate();
   Rprintf("*** TwoLassoCpp::RunCrossValidate(n=%d,p=%d), We successfully TestCrossValidate. \n", n,p); R_FlushConsole();
   Rprintf("*** Onto run RunCrossValidate. \n"); R_FlushConsole();
 }
 F77_CALL(dscal)(&this->p, &ZeroD, this->BackOrigBeta, &One);
 double Stabm1 = 0.0; double Stabm2 = 0.0;  double StabSigmaBar = 0.0;
 if (Verbose >= 2) {
   Rprintf("*** ATwoLassoObject20011.cc: TwoLassoCreate: TwoLassoSexp: Copying into XTY %d, StabSigmaBar = %f\n",  
     this->p, StabSigmaBar); R_FlushConsole();
   Rprintf("*** What is going on?");  R_FlushConsole();
   if (this->SXtY == NULL) {
     Rprintf("*** SXtY is actually a NULL value\n");   R_FlushConsole();
   } else if (!Rf_isNull(this->SXtY)) {
     Rprintf("*** SXtY not a Nil Value, but we'll investigate Rf_length \n"); R_FlushConsole();
     Rprintf("*** SXtY is not null, Rf_length is %d\n", 
       Rf_length(this->SXtY));   R_FlushConsole();
   } else {
     Rprintf("*** SXtY is a Nil\n");     R_FlushConsole();
   }
   if (Verbose > 5) {
   Rprintf("*** Was SXtY null or not? \n"); R_FlushConsole();
   Rprintf("*** CDO's XTYflag is %d\n", this->CDO->XTYPassByPointer);  
     R_FlushConsole();
   }
   if (this->CDO->XTY == NULL) {
     Rprintf("*** this->CDO->XTY is NULL though, so we're going to break.\n");
     R_FlushConsole();
   }
   if (this->BackXTYResid == NULL) {
     Rprintf("*** this->BackXTYResid is going to be null.\n"); R_FlushConsole();
   }
   if (Verbose > 5) {
     for (int jj = 0; jj < this->p; jj++) {
       Rprintf("*** this->CDO->XTY[jj] = %f.\n", this->CDO->XTY[jj]); 
       R_FlushConsole();
     } 
   }
 }
 if (this->BackXTYResid == NULL) {
   Rf_error("Error: No copy into BackXTYResid because it is already NULL!\n");
 }
 F77_CALL(dcopy)(&this->p, (const double*) this->CDO->XTY, &One, 
   this->BackXTYResid, &One);
  

 this->SetupOrigBetaFromCDO();
 if (Verbose >= 1 || (Verbose >= 0 && p >= LARGEDATA)) {
   Rprintf("*** RunCrossValidate(n=%d,p=%d) finished SetupOrigBetaFromCDO\n", n, p);
 }
 jti = 0;  int copylen = 0;
 if (m1 > 0 && m2 > 0) {
   Stabm1 = m1;  Stabm2 = m2;
 }
 if (SigmaBar > 0) {
   StabSigmaBar = SigmaBar;
 }
 if (Verbose >= 2 || (Verbose >= 1 && p >= LARGEDATA)) {
   Rprintf("*** RunCrossValidate(n=%d,p=%d), Logit = %d, we have length inputs: \n",
     this->n, this->p, GLMCDO == NULL ? 0 : 1); R_FlushConsole();
   Rprintf("***  -- RCV(n=%d,p=%d),  PiAVec=[%d], SigVec=[%d], LamAK=[%d], LamDK=[%d].\n",
     n,p,Rf_length(this->SPiAVectorInputs),  Rf_length(this->SSigmaVectorInputs),
     Rf_length(this->SLambdaAKInputs), Rf_length(this->SLambdaDKInputs));
   R_FlushConsole();
 }
 if (Rf_length(this->SPiAVectorInputs) <= 0) {
   Rf_error("*** TwoLassoCpp:: Error SPiAVectorInputs is NULL \n");
 }
 double LastOldPiA = OnPiA;
 if (Verbose >= 5) {
   Rprintf("*** TwoLassoCpp:: LastOldPiA = %f\n", LastOldPiA);
 }
 for (jti = 0; jti < Rf_length(this->SPiAVectorInputs); jti++) {
   if (Verbose == 2 || (Verbose >= 0 && Verbose <= 2 && p >= LARGEDATA)) {
     Rprintf("*** RCV(jti=%d/%d), Logit=%d, Pia=%.7f, Sig=%f, LA=%f, LD=%f, OnKappaS=%d/%d/%d: START : ",
       jti, Rf_length(this->SPiAVectorInputs), GLMCDO == NULL ? 0 : 1, REAL(this->SPiAVectorInputs)[jti],
       REAL(this->SSigmaVectorInputs)[jti], 
       REAL(this->SLambdaAKInputs)[jti * Rf_length(this->SLambdaAK)],
       REAL(this->SLambdaDKInputs)[jti * Rf_length(this->SLambdaDK)],
       CDO->OnKappaS, CDO->OnKappaMem,p);
     R_FlushConsole();
   }

   if (Verbose >= 2) {
     Rprintf("****************************************************\n");
     R_FlushConsole();
     Rprintf("*** TwoLassoCpp::RunCrossValidate()  in ATwoLassoObject2011.cc jti=%d/%d\n", 
       jti,  Rf_length(this->SPiAVectorInputs));
     R_FlushConsole();
     Rprintf("*** TwoLassoCpp::RunCrossValidate() jti=%d/%d, Pia = %f, Sig=%f.\n",
       jti, Rf_length(this->SPiAVectorInputs),  
       REAL(this->SPiAVectorInputs)[jti], REAL(this->SSigmaVectorInputs)[jti]);
     R_FlushConsole();
   }
   if (this->CDO->MemoryFailFlag == 1) {
     Rprintf("*** ATwoLassoObject2011.cc::RunCrossValidate, before jti=%d/%d we had a memory fail, correcting. \n",
       jti, Rf_length(this->SPiAVectorInputs));
     if (REAL(SOnBeta)[0] == -666.0) {
       for (int kti = 0; kti < p; kti++) {
         REAL(SOnBeta)[kti] = 0.0;
       }
       this->CDO->MakeXTResid();
     }
     Rprintf("*** ATwoLassoObject2011.cc::RunCrossValidate,  before jti=%d/%d we will test integrity for fail. \n",
       jti, Rf_length(this->SPiAVectorInputs));
     this->TestCrossValidate();
     this->CDO->MemoryFailFlag = 0;
     this->CDO->SuccessFlag = 1;
     this->CDO->MaximumAllocation += 50; 
   }
   if (REAL(this->SPiAVectorInputs)[jti]>= 0.0) {
     this->OnPiA = REAL(this->SPiAVectorInputs)[jti];
     this->FixKa = -100;
   } else {
     if (REAL(this->SPiAVectorInputs)[jti] <= -1.0) {
       this->OnPiA = fabs((double) REAL(this->SPiAVectorInputs)[jti]) /
          ((double) this->p);
       this->FixKa = abs((int) REAL(this->SPiAVectorInputs)[jti]);
     } else {
       this->OnPiA = fabs((double) REAL(this->SPiAVectorInputs)[jti]);
       this->FixKa = (int) fabs( REAL(this->SPiAVectorInputs)[jti] * this->p);     
     }
   }
   if (RSSigmaVectorInputs != NULL &&
     !Rf_isNull(this->RSSigmaVectorInputs->asSexp()) && Rf_length(this->RSSigmaVectorInputs->asSexp()) > jti) {
     this->OnSigma = REAL(this->RSSigmaVectorInputs->asSexp())[jti];
   } 
   if (!Rf_isNull(SSigma) && Rf_isReal(SSigma)) {
     for (int itt = 0; itt < Rf_length(SSigma); itt++) {
       REAL(SSigma)[itt] = this->OnSigma;
     }
   }
   copylen = Rf_length(this->SLambdaAK);
   if (Rf_length(this->SLambdaAKInputs) >= (jti+1) * copylen) {
     F77_CALL(dcopy)(&copylen, REAL(this->SLambdaAKInputs) + jti * Rf_length(this->SLambdaAK),
       &One, REAL(this->SLambdaAK), &One);
   }
   copylen = Rf_length(this->SLambdaDK);

   if (Rf_length(this->SLambdaDKInputs) >= (jti+1) * copylen) {
     F77_CALL(dcopy)(&copylen, REAL(this->SLambdaDKInputs) + jti * copylen,
       &One, REAL(this->SLambdaDK), &One);
   }
   if (Verbose >= 3) {
     Rprintf("*** RunCrossValidate() jti=%d/%d \n", jti, Rf_length(this->SPiAVectorInputs));
     R_FlushConsole();
     Rprintf("***    SLambdaAK = ");  PrintVector(REAL(this->SLambdaAK), 
       Rf_length(this->SLambdaAK));  R_FlushConsole();
     Rprintf("\n***    SLambdaDK = ");  PrintVector(REAL(this->SLambdaDK), 
       Rf_length(this->SLambdaDK));  Rprintf("\n***\n"); R_FlushConsole();     
   }
   // If we are using ``loose PiAPrior'' for cross validation
   //  Consider changing  m1/m2 according to the PiAVectorInput;
   if (m1 > 0 && m2 > 0) {
      m1 = Stabm1 * OnPiA;
      m2 = Stabm2  * (1.0 - OnPiA);
   }

   // If we are using ``loose SigmaPrior'' for cross validation
   //  Consider changing SigmaBar according to the SigmaVectorInput;
   if (this->SigmaBar > 0  && this->SigmaDf > 0) {
     this->SigmaBar = this->OnSigma;
   } 
    
   this->OnLambdaA = REAL(this->SLambdaAK)[0];
   this->OnLambdaD = REAL(this->SLambdaDK)[0];
   this->tt1 = 0; this->tt2 = 0;
   if (RStt1 != NULL && !Rf_isNull(RStt1->asSexp())) {
     if (Rf_isInteger(RStt1->asSexp())){ INTEGER(RStt1->asSexp())[0] = 0; }
     if (Rf_isReal(RStt1->asSexp())){ REAL(RStt1->asSexp())[0] = 0; }   
   }
   if (RStt2 != NULL && !Rf_isNull(RStt2->asSexp())) {
     if (Rf_isInteger(RStt2->asSexp())){ INTEGER(RStt2->asSexp())[0] = 0; }
     if (Rf_isReal(RStt2->asSexp())){ REAL(RStt2->asSexp())[0] = 0; }   
   }
   
   if (Verbose >= 2) {
     Rprintf("*** TwoLassoCpp->RunCrossValidate, Inside C++, jti = %d/%d, pi = %f, sig = %f\n",
       jti, Rf_length(this->SPiAVectorInputs), this->OnPiA, this->OnSigma);
     R_FlushConsole();
   }
   if (this->Verbose >= 3) {
     Rprintf("*** TwoLassoCpp->RunCrossValidate: Fitting jti = %d/%d, PiA=%f, Sig=%f, LA = %f, LD = %f \n",
       jti, Rf_length(this->SPiAVectorInputs), 
       this->OnPiA, this->OnSigma, this->SLambdaAK, this->SLambdaDK);
     R_FlushConsole();
   }
   /*
   if (jti == 0 || jti % 5 == 0 || LastOldPiA != OnPiA) {
   for (int ii = 0; ii < p; ii++) {
     CDO->OnBeta[ii] = 0.0;
   }
   for (int ii = 0; ii < p; ii++) {
     REAL(this->RSOnBeta->asSexp())[ii];
   }
   for (int ii = 0; ii < p; ii++) {
     REAL(SBBOn1)[ii] = OnPiA;
   }
   for (int ii = 0; ii < p; ii++) {
     CDO->XTResid[ii] = CDO->XTY[ii];
   }
   double NewPen = OnSigma*(OnPiA* OnLambdaA + (1.0-OnPiA) * OnLambdaD);
   for (int ii = 0; ii < p; ii++) {
	   CDO->OnGammas[ii] = NewPen;
   }
   if ((p >= LARGEDATA && this->Verbose >= 2)|| this->Verbose >= 3) {
     Rprintf("*** TwoLassoCpp->RunCrossValidate, jti=%d/%d, we'll run one Coordinatedescent, every Gamma is %f, CDOFactor = %f\n",
       jti, Rf_length(this->SPiAVectorInputs), CDO->OnGammas[0], CDO->FactorGammaDenominator); R_FlushConsole();
   }
   //Rprintf("***  TwoLassoCpp->RunCrossValidate, jti=%d/%d, Test Before initial RunAlgorithm \n", jti, Rf_length(this->SPiAVectorInputs));
   //R_FlushConsole();
   //CDO->TestCDO();
   CDO->RunAlgorithm(MaxCauchy, CauchyEpsilon);
   if ((p >= LARGEDATA && this->Verbose >= 2) || this->Verbose >= 3) {
     Rprintf("*** TwoLassoCpp->RunCrossValidate, jti=%d/%d, we'll one CoordinatedescentRan, every Gamma is %f, CDOFactor=%f\n",
     jti, Rf_length(this->SPiAVectorInputs), CDO->OnGammas[0], CDO->FactorGammaDenominator); R_FlushConsole();
     Rprintf("*** TestCrossValidate! \n");
     int MyTest = 0;
     MyTest = CDO->TestCDO();
     Rprintf("*** CDO->Test has Concluded with MyTest = %d, good?\n", MyTest);
     R_FlushConsole();
     if (MyTest < 0) {
       Rf_error("*** TwoLassoCpp->RunCrossValidate, MyTest = %d signifies Error for TestCDO!\n", MyTest); 
     }
   } 
   }
   */
   LastOldPiA = OnPiA;
     //Rprintf("TwoLassoCrossValidate, Running now with LambdA = ");
     // PrintVector(REAL(this->SLambdaA), Rf_length(this->SLambdaA));
     // Rprintf("\n"); R_FlushConsole();
     //Rprintf("TwoLassoCrossValidate, Running now with LambdD = ");
     // PrintVector(REAL(this->SLambdaD), Rf_length(this->SLambdaD));
     // Rprintf("\n");
   if (Verbose >= 2) {
     Rprintf("*** RunCrossValidate: Calling RunTwoLassoRegression jti=%d/%d\n",
       jti,  Rf_length(this->SPiAVectorInputs)); R_FlushConsole();
   }
   if ((p >= LARGEDATA && Verbose >= -1)) {
     int CNZ = 0;
     for (int ii = 0; ii < p; ii++) {
       if (CDO->OnBeta[ii] != 0.0) {
         CNZ++;
       }
       if (!R_finite(CDO->OnBeta[ii]) || R_isnancpp(CDO->OnBeta[ii])) {
         Rf_error("*** RunCrossValidate(jti=%d/%d), Oh no, OnBeta[%d] = %f!\n", jti, Rf_length(this->SPiAVectorInputs),
           ii, CDO->OnBeta[ii]);
       }
     }
     double AvXtY = 0.0;  double MaxXtY = 0.0;
     for (int ii = 0; ii < p; ii++) {
       if (fabs(CDO->XTY[ii]) > MaxXtY) {
         MaxXtY = fabs(CDO->XTY[ii]);
       }
       AvXtY += fabs(CDO->XTY[ii]) * (1.0 / p);
     }
     //Rprintf("***  RCV(jti=%d/%d): Running now RunTwoLassoRegression() OnKappaS=%d/%d, piA=%.8f, Sig=%f, LA=%f, LD=%f, Count NonZero = %d, AvXtY = %f, MaxXtY = %f\n",
     //  jti, Rf_length(this->SPiAVectorInputs), CDO->OnKappaS, CDO->OnKappaMem,
     //   this->OnPiA, this->OnSigma, this->OnLambdaA, this->OnLambdaD, CNZ, AvXtY, MaxXtY);
     //R_FlushConsole();
     //Rprintf("**** Let's TestCDO !\n");     R_FlushConsole();
     //int AC = TestCDO();
     //if (AC < 0) {
     //  Rf_error("RCV Fail!\n");
     //}
   }
   if (CVQuitTime >= 0 && CVQuitTime == jti) {
     Rprintf("**********************************************************************************************\n");
     Rprintf("*** TwoLassoCpp::RunCrossValidate() jti = %d/%d, we've been ordered to quit CVQuitTime = %d\n",
       jti, Rf_length(this->SPiAVectorInputs), CVQuitTime);
     Rprintf("***\n");
     Rprintf("***              RunCrossValidate, triggering an error. \n"); R_FlushConsole();
     Rprintf("**********************************************************************************************\n");
     Rf_error("*** TwoLassoCpp::RunCrossValidate, CVQuitTime Issue. \n");
   }
   ///Rprintf("*** RCV(jti=%d/%d): TestCDO Before RunTwoLassoRegression \n", jti, Rf_length(this->SPiAVectorInputs));
   //R_FlushConsole();
   //CDO->TestCDO();
   if (jti == -40 && p >= LARGEDATA) {
     this->Verbose = this->Verbose + 3;
     this->CDO->PrintFlag += 1;
     Rprintf("***************************************************************************");
     Rprintf("**  Hey, Test on %d, let's see a memory leak. \n", jti);
     CDO->TestCDO();
     R_FlushConsole(); 
   }

   this->RunTwoLassoRegression(0);  
   //Rprintf("*** RCV(jti=%d/%d): TestCDO After RunTwoLassoRegression \n", jti, Rf_length(this->SPiAVectorInputs));
   //R_FlushConsole();
   //CDO->TestCDO();
   if (jti == -40) {
     Rprintf("***************************************************************************");
     Rprintf("**  Hey, Test on %d, End test to see a memory leak. \n", jti);
     CDO->TestCDO();
     R_FlushConsole(); 
     this->Verbose = this->Verbose - 3;
     this->CDO->PrintFlag -= 1;     
   }   
   int InGoOn = 0;
   if (Verbose == 2 || (Verbose >= 0 && Verbose <= 2 && p >= LARGEDATA)) {
     for (int iti = 0; iti < CDO->OnKappaS; iti++) {
       if (CDO->OnBeta[CDO->kXFinder[iti]] != 0.0) { InGoOn++; }
     }
     Rprintf(": FINISHED : (jti=%d/%d): %d are active, PiA=%.8f, Sig=%f, LA=%f, LD=%f, OKs=%d/%d/%d\n",
       jti, Rf_length(this->SPiAVectorInputs), InGoOn, OnPiA,
       OnSigma, 
       REAL(this->SLambdaAKInputs)[jti * Rf_length(this->SLambdaAK)],
       REAL(this->SLambdaDKInputs)[jti * Rf_length(this->SLambdaDK)],
       CDO->OnKappaS, CDO->OnKappaMem,p);
     R_FlushConsole();
   }
   if (p >= LARGEDATA && Verbose >= 3) {
     InGoOn = 0;
     for (int iti = 0; iti < CDO->OnKappaS; iti++) {
       if (CDO->OnBeta[CDO->kXFinder[iti]] != 0.0) { InGoOn++; }
     }
     Rprintf("**** RCV(jti=%d/%d): DoLogit = %d, completed RunTwoLassoRegression(), %d are on, now OnKappaS=%d/%d/%d, piA=%.8f, Sig=%f, LA=%f, LD=%f\n",
       jti, Rf_length(this->SPiAVectorInputs), GLMCDO == NULL ? 0 : 1, InGoOn, CDO->OnKappaS, CDO->OnKappaMem,
       this->p,
       this->OnPiA, this->OnSigma, this->OnLambdaA, this->OnLambdaD);
     R_FlushConsole();
   }   
   if (this->CDO->MemoryFailFlag == 1) {
      Rprintf("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\n");
      Rprintf("MMM ATwoLassoObject2011.cc::RunCrossValidate() -- Memory failed for jti=%d/%d, piA=%.6f, sig=%.4f, OnA=%.4f, OnD=%.4f, A=%.4f, D=%.4f\n",
       jti, Rf_length(this->SPiAVectorInputs), this->OnPiA, this->OnSigma, this->OnLambdaA, this->OnLambdaD,
       REAL(this->SLambdaAK)[ Rf_length(this->SLambdaDK)-1], 
       REAL(this->SLambdaDK)[ Rf_length(this->SLambdaDK)-1]); R_FlushConsole();   
      int CountOn = 0;
      for (int iit = 0; iit < p; iit++) {
        if (CDO->OnBeta[iit] != 0.0) {
          CountOn++;
        }
      }
      Rprintf("MMM The number of On Coefficients is %d/%d, OnKS = %d, OnKappaMem=%d at Memory Fail. \n",
        CountOn, p, CDO->OnKappaS, CDO->OnKappaMem);
      Rprintf("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\n");
      R_FlushConsole();
       for (int kti = 0; kti < p; kti++) {
         REAL(this->SOnBeta)[kti] = -666;
       } 
   } else if (this->SuccessFlag < 0) {
     Rprintf("*** ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n"); R_FlushConsole();
     Rprintf("*** RunCrossValidate: In Middle of ");
     Rprintf("jti=%/d%d, piA=%f, sig=%f, A=%f, D=%f\n",
       jti, Rf_length(this->SPiAVectorInputs), this->OnPiA, this->OnSigma,
       REAL(this->SLambdaAK)[ Rf_length(this->SLambdaAK)-1],
       REAL(this->SLambdaDK)[ Rf_length(this->SLambdaDK)-1]);
     Rprintf("*** Inspect for this error. \n");
     R_FlushConsole();
     return(-1);
   }
   if (Verbose >= 3 || (this->Verbose >= 2 && p >= LARGEDATA)) {
     int NCountOn = 0;
     for (int iti = 0; iti < CDO->OnKappaS; iti++) {
       if (CDO->kXFinder[iti] >= 0 && CDO->kXFinder[iti] <= p &&
         CDO->OnBeta[CDO->kXFinder[iti]] != 0.0) {
         NCountOn++;
       }
     }
      Rprintf("***** ATwoLassoObject2011.cc::RunCrossValidate(), Successful RunTwoLassoRegression jti=%d/%d, %d Are NonZero, OnKappaS=%d/%d/%d, piA=%.8f, sig=%f, A=%f, D=%f\n",
       jti,  Rf_length(this->SPiAVectorInputs), NCountOn, CDO->OnKappaS, CDO->OnKappaMem, p,
       this->OnPiA, this->OnSigma,
       REAL(this->SLambdaAK)[ Rf_length(this->SLambdaAK)-1],
       REAL(this->SLambdaDK)[ Rf_length(this->SLambdaDK)-1]);
      R_FlushConsole();
   }
   if (GLMCDO != NULL) {
       if (RSRecordBeta0CV == NULL) {
         RSRecordBeta0CV = new AObject(Rf_allocVector(REALSXP, Rf_length(RSPiAVectorInputs->asSexp())));
         for (int idi = 0; idi < Rf_length(RSRecordBeta0CV->asSexp()); idi++) {
           REAL(RSRecordBeta0CV->asSexp())[idi] = 0.0;
         }      
       }
       if (RSRecordBeta0CV == NULL || !Rf_isReal(RSRecordBetaCV->asSexp())) {
         Rf_error("TwoLassoObject2011.cc: Why is RSRecordBeta0CV not a good vector. \n"); 
       } else if (Rf_length(RSRecordBeta0CV->asSexp()) <= jti) {
         Rf_error("TwoLassoObject2011.cc: Why is RSRecordBeta0CV length %d but jti = %d!\n", 
           Rf_length(RSRecordBeta0CV->asSexp()), jti);
       }
       REAL(RSRecordBeta0CV->asSexp())[jti] = GLMCDO->pBeta0[0];
   }
   if (Rf_length(this->RSRecordBetaCV->asSexp()) >=
     (jti+1) * Rf_length(SOnBeta)) {
     One = 1;
     if (Verbose >= 3 || (Verbose >= 2 && p >= LARGEDATA)) {
       int NNonZeroOn = 0;
       for (int itt = 0; itt < p; itt++) {
         if (REAL(this->SOnBeta)[itt] != 0.0) {
           NNonZeroOn++;
         }
       }
       Rprintf("*** ATwoLassoObject2011.cc::RunCrossValidate(), Now about to copy OnBeta length %d/%d for RunTwoLassoRegression jti=%d/%d\n",
         NNonZeroOn, Rf_length(SOnBeta), jti,  Rf_length(this->SPiAVectorInputs));
       R_FlushConsole();
     }
     F77_CALL(dcopy)(&this->p, REAL(this->SOnBeta), &One,
       REAL(this->RSRecordBetaCV->asSexp()) + 
       jti * Rf_length(SOnBeta), &One);

   } else {
     Rprintf("HEYAHEYAHEYAHEYAHEYAHEYAHEYAHEYAHEYAHEYAHEYAHEYAHEYAHEYA\n");
     Rprintf("***  ATwoLassoObject20011.cc -- warning! \n");    R_FlushConsole();
     Rprintf("*** We're doing TwoLassoCpp->RunCrossValidate. \n"); R_FlushConsole();
     Rprintf("***  However, length(RSRecordBetaCV) = %d, length(OnBeta)=%d, jti=%d/%d\n",
       Rf_length(this->RSRecordBetaCV->asSexp()), Rf_length(this->SOnBeta),
       jti, Rf_length(this->SPiAVectorInputs));
     Rprintf("***  This does not work good for saving anything!\n"); R_FlushConsole();
     Rprintf("***  Go Back and prepare a better RecordBetaCV!\n"); R_FlushConsole();
     Rprintf("***  ERROR, Heya!\n"); R_FlushConsole();
     return(-1);
   }
   if (this->Verbose >= 2) {
     Rprintf("*** ATwoLassoObject2011: TwoLassoCpp->RunCrossValidate, finished Run Regression jti=%d/%d\n",
       jti, Rf_length(this->SPiAVectorInputs) ); R_FlushConsole();
     One = 1;
     Rprintf("***    Sum Beta = %f, Beta =",
       F77_CALL(dasum)(&this->p, REAL(this->SOnBeta), &One));
     PrintVector(REAL(this->SOnBeta), this->p); Rprintf("\n"); R_FlushConsole();
   }
   if (p > 5000) {
     for (int ii3 = 0; ii3 < p; ii3++) {
       REAL(this->SOnBeta)[ii3] = 0.0;  CDO->OnBeta[ii3] = 0.0;
     }
     this->CDO->MakeXTResid();
   } else if (!Rf_isNull(this->RSStartBetaMatrix->asSexp()) && p < 5000) {
     if (this->Verbose >= 3 || (this->Verbose >= 2 && p >= LARGEDATA)) {
       Rprintf("*** ATwoLassoObject2011::TwoLassoCpp->RunCrossValidate(), jti=%d/%d, RStartBetaMatrix is length %d to go\n",
         jti, Rf_length(this->SPiAVectorInputs), Rf_length(this->RSStartBetaMatrix->asSexp()));
         R_FlushConsole();
     }
     int TTSSB = Rf_length(this->RSStartBetaMatrix->asSexp()) / this->p;
     int TLen =  Rf_length(this->RSStartBetaMatrix->asSexp());
     if (TTSSB <= 1) {
       Rf_error("** Error TwoLassoCpp: TTSSB = %d, p = %d, bad! \n", TTSSB, this->p);
       R_FlushConsole();
     }
     int iStart = jti % TTSSB;
     if (iStart == TTSSB) { iStart = 0; }
     if (iStart >= TTSSB) { iStart = 0; }
     if (iStart < 0) { iStart = 0; }
     if (iStart*this->p + this->p > TLen) {
       Rf_error("Hey, we ain't copying from beta start, iStart=%d, jti=%d, because length is %d, TLen = %d, p=%d!\n",
         iStart, jti, TTSSB, TLen, this->p);
     }
     if (iStart*this->p + this->p > Rf_length(this->RSStartBetaMatrix->asSexp())) {
       Rf_error("ATwoLassoObject2011::TwoLassoCpp->RunCrossValidate(jti=%d/%d, iStart=%d, length(StartBetaMatrix) = %d, this don't work to copy from StartBetaMatrix\n",
         jti, Rf_length(this->SPiAVectorInputs), iStart, Rf_length(this->RSStartBetaMatrix->asSexp()));
     }
     One = 1;
     F77_CALL(dcopy)(&this->p, 
       REAL(this->RSStartBetaMatrix->asSexp()) + iStart * this->p, 
       &One, REAL(this->SOnBeta), &One);
     if (this->Verbose >= 3 || (this->Verbose >= 2 && p >= LARGEDATA) ) {
         Rprintf("*** ATwoLassoObject2011::TwoLassoCpp->RunCrossValidate(), jti=%d/%d, successful copy of RStartBetaMatrix p=%d from length %d\n",
         jti, Rf_length(this->SPiAVectorInputs), p, Rf_length(this->RSStartBetaMatrix->asSexp()));
         R_FlushConsole();
     }
     this->CDO->MakeXTResid();
     //Rprintf("*** ATwoLassoObject2011::TwoLassoCpp->RunCrossValidate(), jti=%d/%d, we copied a new Beta in, After MakeXTResid is?\n",
     //  jti, Rf_length(SPiAVectorInputs)); R_FlushConsole();
     //CDO->TestCDO();
   } else if (this->BackOrigBeta != NULL)  {
     if (this->Verbose >= 3 || (this->Verbose >= 2 && p >= LARGEDATA) ) {
       Rprintf("*** ATwoLassoObject2011::TwoLassoCpp->RunCrossValidate(), jti=%d/%d, will just copy BackOrigBeta. \n",
         jti, Rf_length(this->SPiAVectorInputs), p, Rf_length(this->RSStartBetaMatrix->asSexp()));
       R_FlushConsole();
     }
     F77_CALL(dcopy)(&this->p, this->BackOrigBeta, &One, 
       REAL(this->SOnBeta), &One);
     F77_CALL(dcopy)(&this->p, this->CDO->OnBeta, &One, 
       this->BackOrigBeta, &One); 
     this->CDO->MakeXTResid();     
     // Rprintf("*** ATwoLassoObject2011::TwoLassoCpp->RunCrossValidate(), jti=%d/%d, we copied a new Beta in, after BackOrigBeta move, After XTResid is?\n",
     //  jti, Rf_length(SPiAVectorInputs)); R_FlushConsole();
     //CDO->TestCDO(); 
   }
   if (Verbose >= 2) {
     Rprintf("***       ATwoLassoObject20011::RunCrossValidate - finished jti = %d \n", jti);
     R_FlushConsole();
   }
   //if (jti <= 5) {
   //  Rprintf("**** ATwoLassoObject20011::RunCrossValidate at End of a loop, jti = %d/%d, TestCDO \n", jti, Rf_length(this->SPiAVectorInputs));
   //  CDO->TestCDO();
   //}
 }
 if (Verbose > 0) {
   Rprintf("*** ATwoLassoObject20011: TwoLassoCpp->RunCrossValidate: Finished Loop\n"); R_FlushConsole();
 }
 if (this->SuccessFlag < 0) {
   Rprintf("*** ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n"); R_FlushConsole();
   Rprintf("*** RunCrossValidate: Error, SuccessFlag is less than zero!");  R_FlushConsole();
   return(-1);
 }
 m1 = Stabm1;  m2 = Stabm2;
 return(1);
}

SEXP TwoLassoReRegress(SEXP sTLS, SEXP sRunFlag, SEXP newMVec, 
  SEXP newSigmaPrior, SEXP newVerbose) {
  int RunFlag = GetFirstInteger(sRunFlag);
  TwoLassoSexp *TLS = (TwoLassoSexp*) R_ExternalPtrAddr(sTLS);
  if (TLS == NULL) {
    Rprintf("TwoLassoReRegression: Sorry, initial TLS is Bad!");
    R_FlushConsole();
    return(R_NilValue);
  }

  if (TLS->Verbose > 0) {
    Rprintf("TwoLassoReRegress: We're in Business!");
  }
  if (TLS->get_SuccessFlag() <= 0) {
    Rprintf("TwoLassoReRegress: No fricken way, successflag <= 0, won't start! \n");
    SEXP sOut = R_NilValue; Rf_protect(sOut = Rf_allocVector(INTSXP, 1));
    INTEGER(sOut)[0] = -1; Rf_unprotect(1);
    return(sOut);
  }
  int FVerbose = -25;
  if (!Rf_isNull(newVerbose) && 
      (Rf_isReal(newVerbose) || isInteger(newVerbose))) {
    if (Rf_length(newVerbose) <= 0) {
      Rf_error("TwoLassoReRegress: you gave newVerbose that is less than= 0 in length!\n");
    }
    if (Rf_isInteger(newVerbose)) {
      FVerbose = INTEGER(newVerbose)[0];
    } else if (Rf_isReal(newVerbose)) {
      FVerbose = (int) REAL(newVerbose)[0];
    } else {
      FVerbose = -25;
    }
    if (FVerbose > -20) {
      if (FVerbose >= 1) {
        Rprintf("TwoLassoReRegress: Assigning a new Verbose=%d! \n", FVerbose); R_FlushConsole();
      }
      TLS->Verbose = FVerbose;
    }
  }
  if (TLS->Verbose >= 1) {
    Rprintf("TwoLassoReRegress: Assigning times back to zero! \n"); R_FlushConsole();
  }
  if (TLS->RStt1 != NULL  && !Rf_isNull(TLS->RStt1->asSexp())) {
    if (Rf_length(TLS->RStt1->asSexp()) <= 0) {
      Rf_error("TwoLassoReRegress: uh oh, Stt1 has <= 0 length!\n");
    }
    if (Rf_isReal(TLS->RStt1->asSexp())) {
      REAL(TLS->RStt1->asSexp())[0] = 0.0;
    } else if (Rf_isInteger(TLS->RStt1->asSexp())) {
      INTEGER(TLS->RStt1->asSexp())[0] = 0;
    }
  }
  if (TLS->RStt2 != NULL) {
    if (Rf_length(TLS->RStt2->asSexp()) <= 0) {
      Rf_error("TwoLassoReRegress: uh oh, Stt2 has <= 0 length!\n");
    }
    if (Rf_isReal(TLS->RStt2->asSexp())) {
      REAL(TLS->RStt2->asSexp())[0] = 0.0;
    } else if (Rf_isInteger(TLS->RStt2->asSexp())) {
      INTEGER(TLS->RStt2->asSexp())[0] = 0;
    }
  }
  if (TLS->Verbose >= 1) {
    Rprintf("TwoLassoReRegress: asigning tt1,tt2 = 0\n."); R_FlushConsole();
  }
  
  TLS->tt1 = 0; TLS->tt2 = 0;
  if (!Rf_isReal(TLS->SPiA)) {
    Rprintf("Warning: TwoLassoReRegress:: SPiA is not a real! \n"); R_FlushConsole();
  }

  if (Rf_isNull(TLS->SSigma) || Rf_length(TLS->SSigma) <= 0 ) {
    Rf_error("TwoLassoReRegress, SSigma of <= 0 length!\n");
  }
  if (!Rf_isReal(TLS->SSigma)) {
    Rprintf("Warning: TwoLassoReRegress:: SSigma is not a real! \n"); R_FlushConsole();
  }
  if (TLS->tt1 < 0 || TLS->tt1 >= Rf_length(TLS->SPiA)) {
    Rf_error("TwoLassoReRegress: Cannot go to piA because tt1=%d, lenPiA=%d\n",
      TLS->tt1, Rf_length(TLS->SPiA)); 
  }
  if (TLS->tt1 < 0 || TLS->tt1 >= Rf_length(TLS->SSigma)) {
    Rf_error("TwoLassoReRegress: Cannot go to OnSigma because tt1=%d, lenSSigma=%d\n",
      TLS->tt1, Rf_length(TLS->SSigma)); 
  }
  if (Rf_isNull(TLS->SPiA) || Rf_length(TLS->SPiA) <= 0 ) {
    Rf_error("TwoLassoReRegress, SSigma of <= 0 length!\n");
  }
  if (!Rf_isReal(TLS->SPiA)) {
    Rprintf("Warning: Setting SPiA to tt1 position of SPiA \n"); R_FlushConsole();
  }
  TLS->OnPiA = REAL(TLS->SPiA)[TLS->tt1];
  TLS->OnSigma = REAL(TLS->SSigma)[TLS->tt1];
  if (TLS->Verbose > 0) {
    Rprintf(" TwoLassoReRegress: Summarizing Sigma/Pi Vectors:\n"); R_FlushConsole();
    Rprintf(" For New ReRegress:: SigmaVector[%d] = ", Rf_length(TLS->SSigma));  R_FlushConsole();
    PrintVector(REAL(TLS->SSigma), Rf_length(TLS->SSigma));
    Rprintf(" For New ReRegress:: PiAVector[%d] = ", Rf_length(TLS->SPiA));  R_FlushConsole();
    PrintVector(REAL(TLS->SPiA), Rf_length(TLS->SPiA));  R_FlushConsole();
  }
  if (newMVec != NULL && !Rf_isNull(newMVec) && !Rf_isReal(newMVec)) {
    Rprintf("Warning: TwoLassoReRegress:: newMVec is not null but also not a real! \n"); R_FlushConsole();
  }
  if (newMVec != NULL && !Rf_isNull(newMVec) && Rf_length(newMVec) == 2) {
    if (TLS->Verbose > 0) {
      Rprintf("TwoLassoReREgress:  Adding new m1, m2\n"); R_FlushConsole();
    }
    TLS->m1 = REAL(newMVec)[0];  TLS->m2 = REAL(newMVec)[1];
  }
  if (!Rf_isNull(newSigmaPrior) && Rf_length(newSigmaPrior) == 2) {
    if (!Rf_isReal(newSigmaPrior)) {
      Rprintf("TwoLassoReRegress:  uh oh, newSigmaPrior is not NULL, but also not real!\n");
      R_FlushConsole();
    }
    TLS->SigmaBar = REAL(newSigmaPrior)[0];  
    TLS->SigmaDf = REAL(newSigmaPrior)[1];
  }
  //  Note:  One really should not adjust OnBeta if one is doing 
  //    the TwoLasso that keeps the ReRegress
  //  This could potentially blow out the system
  //   Also, there's no point to it anyway
  //    However, it's sort of cool to see that you could do it
  if (!Rf_isReal(TLS->SLambdaAK)) {
    Rprintf("TwoLassoReRegress: uh oh, SLambdaAK is not a real!\n"); R_FlushConsole();
  }
  if (!Rf_isReal(TLS->SLambdaDK)) {
    Rprintf("TwoLassoReRegress: uh oh, SLambdaDK is not a real!\n"); R_FlushConsole();
  }
  if (TLS->SRecOnBeta != NULL && !Rf_isNull(TLS->SRecOnBeta)  &&
    !Rf_isReal(TLS->SRecOnBeta)) {
    Rprintf("TwoLassoReRegress: uh oh, SRecOnBeta is not null but no real either! \n"); R_FlushConsole(); 
  }
  if (TLS->Verbose >= 1) {
    Rprintf("TwoLassoReRegress: Setup done, now onto RunFlag1 or RunFlag3 \n"); R_FlushConsole();
  }
  if (RunFlag == 1 || RunFlag == 3) {
    if (TLS->Verbose >= 2) {
      Rprintf("TwoLassoReRegress: RunFlag = %d, in 1or3!\n"); R_FlushConsole();
    }
    TLS->tt1 = 0; TLS->tt2 = 0;
    if (TLS->RStt1 != NULL && !Rf_isNull(TLS->RStt1->asSexp())) {
      if (Rf_isInteger(TLS->RStt1->asSexp())) {
        INTEGER(TLS->RStt1->asSexp())[0] = 0;
      }
      if (Rf_isReal(TLS->RStt1->asSexp())) {
        REAL(TLS->RStt1->asSexp())[0] = 0;
      }
    }
    if (TLS->RStt2 != NULL && !Rf_isNull(TLS->RStt2->asSexp())) {
      if (Rf_isInteger(TLS->RStt2->asSexp())) {
        INTEGER(TLS->RStt2->asSexp())[0] = 0;
      }
      if (Rf_isReal(TLS->RStt2->asSexp())) {
        REAL(TLS->RStt2->asSexp())[0] = 0;
      }
    }
    if (TLS->Verbose >= 2) {
      Rprintf("TwoLassoReRegress, getting first members of OnLambdaA/D.\n"); R_FlushConsole();
    }
    if (Rf_isNull(TLS->SLambdaAK) || Rf_length(TLS->SLambdaAK) <= 0  ||
      !Rf_isReal(TLS->SLambdaAK)) {
      Rf_error("TwoLassoReRegress: SLambdaAK is problematic!\n");  
    }
    if (Rf_isNull(TLS->SLambdaDK) || Rf_length(TLS->SLambdaDK) <= 0  ||
      !Rf_isReal(TLS->SLambdaDK)) {
      Rf_error("TwoLassoReRegress: SLambdaAK is problematic!\n");  
    }
    TLS->OnLambdaA = REAL(TLS->SLambdaAK)[0];
    TLS->OnLambdaD = REAL(TLS->SLambdaDK)[0];
    TLS->OnPiA = REAL(TLS->SPiA)[0];  TLS->OnSigma = REAL(TLS->SSigma)[0];

    if (TLS->Verbose >= 2) {
      Rprintf("TwoLassoReRegress: Running first UpdateBBOn1\n"); R_FlushConsole();
    }
    TLS->UpdateBBOn1();
    if (RunFlag == 3 && TLS->BackXTYResid != NULL && 
      TLS->SRecOnBeta != NULL && !Rf_isNull(TLS->SRecOnBeta) && Rf_length(TLS->SRecOnBeta) >= TLS->p) {
      int One = 1;
      F77_CALL(dcopy)(&(TLS->p), REAL(TLS->SRecOnBeta), &One,
        REAL(TLS->SOnBeta), &One);
      F77_CALL(dcopy)(&(TLS->p), TLS->BackXTYResid, &One,
        TLS->CDO->XTResid, &One);         
    } else {
      TLS->CDO->MakeXTResid();
    }
    if (TLS->Verbose >= 2) {
      Rprintf("TwoLassoReRegress: About to RunTowLassoRegression.\n"); R_FlushConsole();
    }
    TLS->RunTwoLassoRegression(1);
  } else if (RunFlag > 3) {
    if (TLS->Verbose >= 2) {
      Rprintf("TwoLassoReRegress: Execute RunFlag = %d version!\n", RunFlag); R_FlushConsole();
    }
    TLS->tt1 = 0; TLS->tt2 = 0;
    if (TLS->RStt1 != NULL && !Rf_isNull(TLS->RStt1->asSexp())) {
      if (Rf_isInteger(TLS->RStt1->asSexp())) {
        INTEGER(TLS->RStt1->asSexp())[0] = 0;
      }
      if (Rf_isReal(TLS->RStt1->asSexp())) {
        REAL(TLS->RStt1->asSexp())[0] = 0;
      }
    }
    if (TLS->RStt2 != NULL && !Rf_isNull(TLS->RStt2->asSexp())) {
      if (Rf_isInteger(TLS->RStt2->asSexp())) {
        INTEGER(TLS->RStt2->asSexp())[0] = 0;
      }
      if (Rf_isReal(TLS->RStt2->asSexp())) {
        REAL(TLS->RStt2->asSexp())[0] = 0;
      }
    }
    TLS->SetupCDO(1);     
    if (TLS->Verbose >= 2 || (TLS->p >= LARGEDATA && TLS->Verbose >= 1)) {
      Rprintf("TwoLassoReRegress: Execute AddAllNewNonzeroCoords. \n"); R_FlushConsole();
    }                                  
    TLS->CDO->AddAllNewNonzeroCoords();
    if (TLS->Verbose >= 2 || (TLS->p >= LARGEDATA && TLS->Verbose >= 1)) {
      Rprintf("TwoLassoReRegress: Execute MakeXTResid. \n"); R_FlushConsole();
    } 
    TLS->CDO->MakeXTResid();	
    if (TLS->Verbose >= 2 || (TLS->p >= LARGEDATA && TLS->Verbose >= 1)) {
      Rprintf("TwoLassoReRegress: Execute NowRunTwoLassoRegression. \n"); R_FlushConsole();
    }   
    TLS->RunTwoLassoRegression(1); 
  }
  
  if (TLS->Verbose >= 1) {
    Rprintf("TwoLassoReRegress: Finished  Return sTLS!\n"); R_FlushConsole();
  }
  //REAL(TLS->SOnBeta)[10] = 0;      Why would I ever do this?
  return(sTLS);
}
    
////////////////////////////////////////////////////////////////////////////////
SEXP TwoLassoRegression(SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaA, SEXP rSLambdaD, SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords, SEXP rSVerbose,
    SEXP rGroupSexp) {    
 int Verbose = 0;
 if (isInteger(rSVerbose)) {
   Verbose = INTEGER(rSVerbose)[0];
 }  else if (Rf_isReal(rSVerbose)) {
   Verbose = (int) REAL(rSVerbose)[0];
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoRegression: Starting by loading object:\n"); R_FlushConsole();
 }
 TwoLassoSexp *TLS = new TwoLassoSexp(rSn, rSp,  rSyys, 
   rSxxs,  rSXtY, rSXtX,  rStt1, rStt2, rSLambdaA,  rSLambdaD,  rSOrderSeq, 
   rSBBOn1,  rSOnBeta,  rRecBBOn1,  rRecOnBeta,
   rSOnGammas,  rSPiA,  rSSigma,
   rSm,  rSSigmaPrior,
   rSInitKKs,  rSInverseGammaConstant, 
   rSCauchyEpsilon,  rSMaxCauchy, 
   rSTDFNu,  rSiiWeights,  
   rSL2ShrinkagePrior, rSL2ShrinkageRecords, rSVerbose);

 if (Verbose > 0) {
   Rprintf("TwoLassoRegression: We've got an object\n"); R_FlushConsole();
 }
 if (TLS == NULL || TLS->SuccessFlag < 0) {
   Rprintf("TwoLassoRegression: Sorry, some error in the Load");
   REAL(rSn)[0] = -999;
   R_FlushConsole();
   return(R_NilValue);
 }
 if (FALSE && !Rf_isNull(rGroupSexp)) {
   TLS->SetupGroups(rGroupSexp);
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoRegression: Going to Run Now!\n"); R_FlushConsole();
 }
 TLS->RunTwoLassoRegression(1);  
 if (TLS == NULL || TLS->SuccessFlag < 0) {
   Rprintf("TwoLassoRegression: Sorry, some error in the algorithm");
   REAL(rSn)[0] = -999;
   R_FlushConsole();
   return(R_NilValue);
 }
 if (Verbose > 0) {
   Rprintf("TwoLassoRegression: Starting to Delete!\n"); R_FlushConsole();
 } 
 delete(TLS);  
 if (Verbose > 0) {
   Rprintf("TwoLassoRegression: All Done, returning!\n"); R_FlushConsole();
 } 
 return(rSOnBeta);
} 

int TwoLassoSexp::UpdateL2Shrinkage() {
  double TargetL2 = 0;       int One = 1;
  int ii;
  double TotalActive = 0;
  if (CDO->kXFinder != NULL && CDO->OnKappaS < .5 *p) {
    for (ii = 0; ii < CDO->OnKappaS; ii++) {
      TargetL2 = REAL(SOnBeta)[CDO->kXFinder[ii]] * REAL(SOnBeta)[CDO->kXFinder[ii]];
    }
    for (ii = 0; ii < CDO->OnKappaS; ii++) {
      TotalActive += REAL(SBBOn1)[CDO->kXFinder[ii]];
    }
  } else {
    TargetL2 = F77_CALL(ddot)(&p, REAL(SOnBeta), &One, REAL(SOnBeta), &One);
    TotalActive = F77_CALL(dasum)(&p, REAL(SBBOn1), &One);  
  }
  if (Rf_isNull(SL2ShrinkagePrior) || !Rf_isReal(SL2ShrinkagePrior)) {
    Rprintf("UpdateL2Shrinkage: Cannot do since L2Shrinkage Prior is NULL!\n");
    R_FlushConsole();
    Rf_error("UpdateL2Shrinkage:  An Error \n");
  }
  TargetL2 = (TargetL2 + REAL(SL2ShrinkagePrior)[1])/
    (REAL(SL2ShrinkagePrior)[0] + TotalActive);
  if (TargetL2 > 0) {
   CDO->DiagShrinkage = REAL(SSigma)[0] / TargetL2;
  } else {
    CDO->DiagShrinkage = 0;
  }
  if (!Rf_isNull(SL2ShrinkageRecords) && Rf_length(SL2ShrinkageRecords) > tt1+1) {
    REAL(SL2ShrinkageRecords)[tt1] = CDO->DiagShrinkage;
  }
  return(1);
}

////////////////////////////////////////////////////////////////////////////////

int TwoLassoSexp::RunTwoLassoRegression(int DoHistory) {

  double MovedBeta = 0.0;
  if (Rf_isNull(Stt1)) { Rf_error("RunTwoLassoRegression: Stt1 is NIL");}
  if (Rf_isNull(Stt2)) { Rf_error("RunTwoLassoRegression: Stt2 is NIL");}
  if (Rf_isNull(SLambdaAK)) { Rf_error("RunTwoLassoRegression: SLambdaAK is still NIL");}
  if (Rf_isNull(SLambdaDK)) { Rf_error("RunTwoLassoRegression: SLambdaDK is still NIL");}
  if (!Rf_isReal(SLambdaAK)) {Rf_error("RunTwoLassoRegression: Error SLambdaAK is not real!");}
  if (Rf_length(SLambdaDK) != Rf_length(SLambdaAK)) {
    Rf_error("RunTwoLassoRegression: Error: SLambdaDK is short compared to AK!\n");}
  if (Rf_length(SLambdaAK) <= 0) {
    Rf_error("RunTwoLassoRegression: SLambdaAK has negligible length!");
  }
  if (!Rf_isReal(SLambdaDK)) { Rf_error("RunTwoLassoRegression: Error SLambdaDK is not real!");}
  if (OrderSeq == NULL) { Rf_error("RunTwoLassoRegression: OrderSeq is still NULL pointer");}  
  if (Rf_isNull(SSigma)) { Rf_error("RunTwoLassoRegression: SSigma is still NIL");}  
  if (Rf_isNull(SSigma)) { Rf_error("RunTwoLassoRegression: SSigma is still NIL");}
  //if (Verbose > 1) {
  //  Rprintf("TwoLassoRegression: tt1 = %d, off to start!\n", tt1); R_FlushConsole();
  //}
  if(Verbose >= 2) {
    Rprintf("############################################################\n");
    R_FlushConsole();
    Rprintf("#### TwoLassoSexp: Starting RunTwoLassoRegression \n"); R_FlushConsole();
    Rprintf("####   tt1 = %d(+1=%d)/%d, tt2=%d/%d, Pia=%.4f. Sig=%.4f, A[0]=%.3f, D[0]=%.3f\n",
      tt1, 
      (GetFirstInteger(Stt1)), Rf_length(SLambdaAK),
      tt2, OrderSeq[tt1], 
       this->OnPiA, this->OnSigma,
      REAL(SLambdaAK)[0], REAL(SLambdaDK)[0]); R_FlushConsole();
    Rprintf("####   Begin.\n"); R_FlushConsole();
  }
  //tt1 = (GetFirstInteger(Stt1))-1;
  for (; tt1 < Rf_length(SLambdaAK); tt1++) {
    if (Rf_isReal(Stt1)) { REAL(Stt1)[0] = (double) tt1+1; }
    if (isInteger(Stt1)) { INTEGER(Stt1)[0] = tt1+1; }
    if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = 1; }
    if (isInteger(Stt2)) { INTEGER(Stt2)[0] = (int) 1; }
    if(Verbose >= 2) {
      Rprintf("#### RunTwoLassoRegression: tt1 = %d/%d\n", tt1,
        Rf_length(SLambdaAK)); R_FlushConsole();
    }
    if (Rf_length(SLambdaAK) > tt1) {
      OnLambdaA = REAL(SLambdaAK)[tt1];
    }
    if (Rf_length(SLambdaDK) > tt1) {
      OnLambdaD = REAL(SLambdaDK)[tt1];
    }
    if(Verbose >=1) {
      Rprintf("#### RunTwoLassoRegression(tt1=%d/%d): LA = %.4f, LD = %.4f, OS = %d\n", 
        tt1, Rf_length(SLambdaAK), 
        OnLambdaA, OnLambdaD, OrderSeq[tt1]); R_FlushConsole();
    }
    if (Rf_length(SOrderSeq) <= tt1) {
      Rf_error("#### Error: RunTwoLassoRegression: SOrderSeq is two short\n");
    }                                                           
    MovedBeta = 0.0;
    if (!Rf_isNull(SSigma) && Rf_isReal(SSigma) && Rf_length(SSigma) >= 1 && 
      Rf_length(SSigma) > tt1) {
      OnSigma = REAL(SSigma)[tt1];
    }
    if (R_isnancpp(CDO->OnGammas[0])) {      
      Rf_error("Error: RunTwoLassoRegression: CDO->OnGammas went nan!\n");
    }
    for (tt2 = 0; tt2 < abs(OrderSeq[tt1]); tt2++) {
      if (Rf_isNull(Stt2)) {
      } else if (Rf_length(Stt2) == 1) {
        if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = tt2+1; }
        if (isInteger(Stt2)) { INTEGER(Stt2)[0] = (int) tt2+1; }
      } else if (Rf_length(Stt2) > tt1) {
        if (Rf_isReal(Stt2)) { REAL(Stt2)[tt1] = tt2+1; }
        if (isInteger(Stt2)) { INTEGER(Stt2)[tt1] = (int) tt2+1; }      
      }
      if(Verbose > 3) {
        Rprintf("####  RunTwoLassoRegression(tt1=%%d/%d, tt2 = %d/%d)\n", 
          tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
        R_FlushConsole();
      }
      if (OrderSeq[tt1] >= 0) {
        UpdateBeta();
        if(Verbose > 3) {
          Rprintf("####  RunTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): Finished UpdateBeta, Onto SigmaSq\n",
             tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
          R_FlushConsole();
        }

        if (CDO->MemoryFailFlag == 1 || SuccessFlag < 0) {
          Rprintf("#### RunTwoLassoRegression(tt1=%d/%d), after UpdateBeta(1), Memory Fail! \n", tt1, tt2);
          Rprintf("####    OnLambdaA=%f, OnLambdAD=%f, OnPiA=%f, OnSigma=%f, OnKappaS=%d/%d/%d",
            OnLambdaA, OnLambdaD, OnPiA, OnSigma, CDO->OnKappaS, CDO->OnKappaMem, p); R_FlushConsole();
          return(-1);
        }
        if (TDFNu > 0 && CDO->TDFNoise > 0 && CDO->TDFSigma > 0) {
          CDO->UpdateTDFNoise();
        }
        if (SigmaBar > 0) {
          UpdateSigmaSq();
        }
        if(Verbose > 3) {
          Rprintf("####  RunTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): Finished UpdateSigma, onto BBOn1\n", 
            tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
          R_FlushConsole();
        }             
        if (GLMCDO != NULL) {
          if (Verbose >= 1) {
            Rprintf("#### RunGLMTwoLassoRegression(tt1=%d, tt2=%d), LambdaA(%f), LambdaD(%f), Sigma=%f, Now UpdateWeightLambdaALambdaDSeq()\n",
             tt1, tt2, OnLambdaA, OnLambdaD, OnSigma); R_FlushConsole();
          }
          UpdateWeightLambdaALambdaDSeq();
          if (Verbose >= 1) {
            Rprintf("#### RunGLMTwoLassoRegression(tt1=%d, tt2=%d), After UpdateWeightLambdaALambdaDSeq(), LambdaA(%f), LambdaD(%f), Sigma=%f\n",
             tt1, tt2, OnLambdaA, OnLambdaD, OnSigma); R_FlushConsole();
          }
        }
        UpdateBBOn1();
        if(Verbose > 3) {
          Rprintf("####  RunTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): Finished UpdateBBOn1\n", 
            tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
          R_FlushConsole();
        }   
        if (m1 > 0) {
          if(Verbose > 3) {
            Rprintf("####  RunTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): UpdateHatPiA\n", 
              tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
            R_FlushConsole();
          }
          UpdateHatPiA();
        }
        if (SL2ShrinkagePrior != NULL && !Rf_isNull(SL2ShrinkagePrior)) {
          UpdateL2Shrinkage();
        }    
        
      } else {
        UpdateBBOn1();
        if (m1 > 0) {
          UpdateHatPiA();
        }  
        UpdateBeta();
        if (CDO->MemoryFailFlag == 1 || SuccessFlag < 0) {
          Rprintf("#### RunTwoLassoRegression(tt1=%d/%d), after UpdateBeta(2), Memory Fail! \n", tt1, tt2);
          Rprintf("####    OnLambdaA=%f, OnLambdAD=%f, OnPiA=%f, OnSigma=%f, OnKappaS=%d/%d/%d",
            OnLambdaA, OnLambdaD, OnPiA, OnSigma, CDO->OnKappaS, CDO->OnKappaMem, p); R_FlushConsole();
          return(-1);
        }
        if (GLMCDO != NULL) {
          UpdateWeightLambdaALambdaDSeq();
        }
        if(Verbose > 3) {
          Rprintf("####  RunTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): Finished UpdateBeta, onto SigmaBar\n", 
            tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
          R_FlushConsole();
        }
        if (TDFNu > 0 && CDO->TDFNoise > 0 && CDO->TDFSigma > 0) {
          CDO->UpdateTDFNoise();
        }
        if (GLMCDO == NULL && SigmaBar > 0.0) {
          UpdateSigmaSq();
        }
        if (SL2ShrinkagePrior != NULL && !Rf_isNull(SL2ShrinkagePrior)) {
          UpdateL2Shrinkage();
        }
      }
      if(Verbose > 3) {
        if (GLMCDO == NULL) {
        Rprintf("####  RunTwoLassoRegression(tt1=%d, tt2=%d): RecentMove = %f, Cep = %f, CDO->TotMove = %f\n",
          tt1, tt2, CDO->RecentMove, CauchyEpsilon, CDO->TotMove); 
        } else {
        Rprintf("####  RunTwoLassoRegression(tt1=%d, tt2=%d): RecentMove = %f, Cep = %f, CDO->GLMTotMove = %f\n",
          tt1, tt2, CDO->RecentMove, CauchyEpsilon, GLMCDO->GLMTotMove); 
        }
        R_FlushConsole();
      }
      if ((GLMCDO != NULL && GLMCDO->GLMTotMove < CauchyEpsilon) ||
        (GLMCDO == NULL && CDO->TotMove < CauchyEpsilon)) {
        tt2  = abs(OrderSeq[tt1]);
      } else {
        MovedBeta += CDO->TotMove;
      }
      if(Verbose > 3) {
        Rprintf("####  RunTwoLassoRegression: tt1=%d,tt2=%d,Moved= %f\n",
          tt1, tt2, MovedBeta); 
        R_FlushConsole();
      }
    }
    if(Verbose > 3) {
      Rprintf("####  RunTwoLassoRegression(tt1=%d,tt2=%d) Finished Loop, Record?\n",
          tt1, tt2); 
        R_FlushConsole();
    }
    if (DoHistory > 0) {
      RecordHistory();
    }
    if (BackXTYResid != NULL) {
      int One = 1;
      F77_CALL(dcopy)(&p, CDO->XTResid, &One, BackXTYResid, &One);
    }
    if (BackOrigBeta != NULL) {
      int One = 1;
      F77_CALL(dcopy)(&p, CDO->OnBeta, &One, BackOrigBeta, &One);
    }
  }
  return(1);
}

int TwoLassoSexp::RunGLMTwoLassoRegression(int DoHistory) {
  if (GLMCDO == NULL) {
    Rf_error("RunGLMTwoLassoRegression: Cannot run this if GLMCDO is NULL!\n");
  }
  double MovedBeta = 0.0;
  if (Rf_isNull(Stt1)) { Rf_error("RunGLMTwoLassoRegression: Stt1 is NIL");}
  if (Rf_isNull(Stt2)) { Rf_error("RunGLMTwoLassoRegression: Stt2 is NIL");}
  if (Rf_isNull(SLambdaAK)) { Rf_error("RunGLMTwoLassoRegression: SLambdaAK is still NIL");}
  if (Rf_isNull(SLambdaDK)) { Rf_error("RunGLMTwoLassoRegression: SLambdaDK is still NIL");}
  if (!Rf_isReal(SLambdaAK)) {Rf_error("RunGLMTwoLassoRegression: Error SLambdaAK is not real!");}
  if (Rf_length(SLambdaDK) != Rf_length(SLambdaAK)) {
    Rf_error("RunGLMTwoLassoRegression: Error: SLambdaDK is short compared to AK!\n");}
  if (Rf_length(SLambdaAK) <= 0) {
    Rf_error("RunGLMTwoLassoRegression: SLambdaAK has negligible length!");
  }
  if (!Rf_isReal(SLambdaDK)) { Rf_error("RunGLMTwoLassoRegression: Error SLambdaDK is not real!");}
  if (OrderSeq == NULL) { Rf_error("RunGLMTwoLassoRegression: OrderSeq is still NULL pointer");}  
  if (Rf_isNull(SSigma)) { Rf_error("RunGLMTwoLassoRegression: SSigma is still NIL");}  
  if (Rf_isNull(SSigma)) { Rf_error("RunGLMTwoLassoRegression: SSigma is still NIL");}
  //if (Verbose > 1) {
  //  Rprintf("TwoLassoRegression: tt1 = %d, off to start!\n", tt1); R_FlushConsole();
  //}
  if(Verbose >= 2) {
    Rprintf("############################################################\n");
    R_FlushConsole();
    Rprintf("#### TwoLassoSep: Starting RunGLMTwoLassoRegression \n"); R_FlushConsole();
    Rprintf("####   tt1 = %d(+1=%d)/%d, tt2=%d/%d, Pia=%.4f. Sig=%.4f, A[0]=%.3f, D[0]=%.3f\n",
      tt1, 
      (GetFirstInteger(Stt1)), Rf_length(SLambdaAK),
      tt2, OrderSeq[tt1], 
       this->OnPiA, this->OnSigma,
      REAL(SLambdaAK)[0], REAL(SLambdaDK)[0]); R_FlushConsole();
    Rprintf("####   Begin.\n"); R_FlushConsole();
  }
  //tt1 = (GetFirstInteger(Stt1))-1;
  for (; tt1 < Rf_length(SLambdaAK); tt1++) {
    if (Rf_isReal(Stt1)) { REAL(Stt1)[0] = (double) tt1+1; }
    if (isInteger(Stt1)) { INTEGER(Stt1)[0] = tt1+1; }
    if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = 1; }
    if (isInteger(Stt2)) { INTEGER(Stt2)[0] = (int) 1; }
    if(Verbose >= 2) {
      Rprintf("#### RunGLMTwoLassoRegression: tt1 = %d/%d\n", tt1,
        Rf_length(SLambdaAK)); R_FlushConsole();
    }
    if (Rf_length(SLambdaAK) > tt1) {
      OnLambdaA = REAL(SLambdaAK)[tt1];
      OnLambdaD = REAL(SLambdaDK)[tt1];
    }
    if(Verbose >=1) {
      Rprintf("#### ATwoLassoObject2014.cc::RunGLMTwoLassoRegression() GLMCDO(tt1=%d/%d): LA = %.4f, LD = %.4f, OS = %d\n", 
        tt1, Rf_length(SLambdaAK), 
        OnLambdaA, OnLambdaD, OrderSeq[tt1]); R_FlushConsole();
    }
    if (Rf_length(SOrderSeq) <= tt1) {
      Rf_error("#### Error: GLMCDO: SOrderSeq is two short\n");
    }                                                           
    MovedBeta = 0.0;
    if (Rf_length(SSigma) > tt1) {
      OnSigma = REAL(SSigma)[tt1];
    }
    if (R_isnancpp(GLMCDO->OnGammas[0])) {      
      Rf_error("Error: ATwoLassoObject2014.cc::RunGLMTwoLassoRegression: GLMCDO->OnGammas went nan!\n");
    }
    for (tt2 = 0; tt2 < abs(OrderSeq[tt1]); tt2++) {
      if (Rf_length(Stt2) == 1) {
        if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = tt2+1; }
        if (isInteger(Stt2)) { INTEGER(Stt2)[0] = (int) tt2+1; }
      } else if (Rf_length(Stt2) > tt1) {
        if (Rf_isReal(Stt2)) { REAL(Stt2)[tt1] = tt2+1; }
        if (isInteger(Stt2)) { INTEGER(Stt2)[tt1] = (int) tt2+1; }      
      }
      if(Verbose > 3) {
        Rprintf("####  RunGLMTwoLassoRegression(tt1=%%d/%d, tt2 = %d/%d)\n", 
          tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
        R_FlushConsole();
      }
      if (OrderSeq[tt1] >= 0) {
        UpdateBeta();
        if (CDO->SuccessFlag < 0) {
          Rprintf("#### RunGLMTwoLassoRegression(tt1=%d/%d), updateBeta(3) Failed \n", tt1, tt2);
          Rprintf("#### OnLambdaA=%f, OnLambdaD=%f, what went wrong?\n", OnLambdaA, OnLambdaD);
          R_FlushConsole();
        }
        if (CDO->MemoryFailFlag == 1 || SuccessFlag < 0) {
          Rprintf("#### RunGLMTwoLassoRegression(tt1=%d/%d), after UpdateBeta(3), Memory Fail! \n", tt1, tt2);
          Rprintf("####    OnLambdaA=%f, OnLambdAD=%f, OnPiA=%f, OnSigma=%f, OnKappaS=%d/%d/%d",
            OnLambdaA, OnLambdaD, OnPiA, OnSigma, CDO->OnKappaS, CDO->OnKappaMem, p); R_FlushConsole();
          return(-1);
        }
        
        if (Verbose > 1) {
          Rprintf("### RunGLMTwoLasso(tt1=%d/%d), LambdaALambdaD will be updated. \n",tt1,tt2); R_FlushConsole();
        }
        UpdateWeightLambdaALambdaDSeq();
        if (Verbose > 1) {
          Rprintf("### RunGLMTwoLasso(tt1=%d/%d), Completed LambdaA(%f), LambdaD(%f) updated. \n",tt1,tt2,
            OnLambdaA, OnLambdaD); R_FlushConsole();
        }        
        if(Verbose > 3) {
          Rprintf("####  RunGLMTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): Finished UpdateBeta, Onto SigmaSq\n",
             tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
          R_FlushConsole();
        }

        UpdateBBOn1();
        if(Verbose > 3) {
          Rprintf("####  RunGLMTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): Finished UpdateBBOn1\n", 
            tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
          R_FlushConsole();
        }
        if (TDFNu > 0 && CDO->TDFNoise > 0 && CDO->TDFSigma > 0) {
          CDO->UpdateTDFNoise();
        }
        if (GLMCDO == NULL && SigmaBar > 0.0) {
          UpdateSigmaSq();
        }
        if(Verbose > 3) {
          Rprintf("####  RunGLMTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): Finished UpdateSigma, onto BBOn1\n", 
            tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
          R_FlushConsole();
        }
        if (m1 > 0) {
          if(Verbose > 3) {
            Rprintf("####  RunGLMTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): UpdateHatPiA\n", 
              tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
            R_FlushConsole();
          }
          UpdateHatPiA();
        }
        if (SL2ShrinkagePrior != NULL && !Rf_isNull(SL2ShrinkagePrior)) {
          UpdateL2Shrinkage();
        }    
        
      } else {
        UpdateBBOn1();
        if (Verbose > 1) {
          Rprintf("### RunGLMTwoLasso(tt1=%d/%d), LambdaA(%f), LambdaD(%f), OnSigma=%f will be updated. \n",tt1,tt2, OnLambdaA,
            OnLambdaD, OnSigma); R_FlushConsole();
        }        
        UpdateWeightLambdaALambdaDSeq();
        if (Verbose > 1) {
          Rprintf("### RunTwoLasso(tt1=%d/%d), Completed LambdaA(%f), LambdaD(%f), OnSigma=%f update. \n",
            tt1,tt2, OnLambdaA, OnLambdaD, OnSigma); R_FlushConsole();
        }

        if (TDFNu > 0 && CDO->TDFNoise > 0 && CDO->TDFSigma > 0) {
          CDO->UpdateTDFNoise();
        }
        if (GLMCDO == NULL && SigmaBar > 0) {
          UpdateSigmaSq();
        }
        if (m1 > 0) {
          UpdateHatPiA();
        }  
        UpdateBeta();
        if (CDO->SuccessFlag < 0) {
          Rprintf("#### RunGLMTwoLassoRegression(tt1=%d/%d), updateBeta(3c) Failed \n", tt1, tt2);
          Rprintf("#### OnLambdaA=%f, OnLambdaD=%f, what went wrong?\n", OnLambdaA, OnLambdaD);
          R_FlushConsole();
        }
        if (CDO->MemoryFailFlag == 1 || SuccessFlag < 0) {
          Rprintf("#### RunGLMTwoLassoRegression(tt1=%d/%d), after UpdateBeta(3), Memory Fail! \n", tt1, tt2);
          Rprintf("####    OnLambdaA=%f, OnLambdAD=%f, OnPiA=%f, OnSigma=%f, OnKappaS=%d/%d/%d",
            OnLambdaA, OnLambdaD, OnPiA, OnSigma, CDO->OnKappaS, CDO->OnKappaMem, p); R_FlushConsole();
          return(-1);
        }

        if(Verbose > 3) {
          Rprintf("####  RunGLMTwoLassoRegression(tt1=%d/%d, tt2=%d/%d): Finished UpdateBeta, onto SigmaBar\n", 
            tt1, Rf_length(SLambdaAK), tt2, abs(OrderSeq[tt1])); 
          R_FlushConsole();
        }
        if (SL2ShrinkagePrior != NULL && !Rf_isNull(SL2ShrinkagePrior)) {
          UpdateL2Shrinkage();
        }
      }
      if(Verbose > 3) {
        Rprintf("####  RunGLMTwoLassoRegression: RecentMove = %f, Cep = %f\n",
          GLMCDO->RecentMove, CauchyEpsilon); 
        R_FlushConsole();
      }
      if ((GLMCDO != NULL && GLMCDO->GLMTotMove < CauchyEpsilon) ||
       (GLMCDO == NULL && CDO->TotMove < CauchyEpsilon)) {
        tt2  = abs(OrderSeq[tt1]);
      } else {
        MovedBeta += GLMCDO->GLMTotMove;
      }
      if(Verbose > 3) {
        Rprintf("####  RunGLMTwoLassoRegression: t1%d,t2%d,Moved= %f\n",
          tt1, tt2, MovedBeta); 
        R_FlushConsole();
      }
      Rprintf("#### RunGLMTwoLassoRegression:GLM, checking Integrity. \n"); R_FlushConsole();
      if (GLMCDO != NULL) {
        if (fabs(GLMCDO->OnBeta[0]) > 100000000.0) {
          Rprintf("RunGLMTwoLassoRegression:GLM, We're going to say that we got OnBeta[0] = %f, and we're broke\n",
            GLMCDO->OnBeta[0]); R_FlushConsole();
          Rf_error("  RunGLMTwoLassoRegression:GLM Error. \n");
        }
        GLMCDO->DataGLMIntegrity(0);
      }
    }
    if(Verbose > 3) {
      Rprintf("####  RunGLMTwoLassoRegression: t1%d,t2%d Finished Loop, Record?\n",
          tt1, tt2); 
        R_FlushConsole();
    }
    if (DoHistory > 0) {
      RecordHistory();
    }
    if (BackXTYResid != NULL) {
      int One = 1;
      F77_CALL(dcopy)(&p, GLMCDO->XTResid, &One, BackXTYResid, &One);
    }
    if (BackOrigBeta != NULL) {
      int One = 1;
      F77_CALL(dcopy)(&p, GLMCDO->OnBeta, &One, BackOrigBeta, &One);
    }
  }
  return(1);
}

int TwoLassoSexp::UpdateBBOn1() {
  if(Verbose > 4) {
    Rprintf("#### ATwoLassoObject2011.cc:: UpdateBBOn1(tt1=%d, tt2=%d):", tt1, tt2); 
    Rprintf("  Updating BBOn1\n",
      tt2); R_FlushConsole();
  }
  if (NumGroups >= 1) {
    return(GroupTwoUpdateBBOn1());
  }
	int ii;
	double OnlambdaDiff = (OnLambdaD - OnLambdaA) / 2.0;
	if (Verbose > 4) {
    Rprintf("####    Update BBOn1: OnlambdaDiff = %f, OA = %f, OD = %f\n", 
      OnlambdaDiff, OnLambdaA, OnLambdaD);
    R_FlushConsole();
  }
	int SuccFlag = 0;
	if (SBBOn1 == NULL || Rf_isNull(SBBOn1) || Rf_length(SBBOn1) < p ||
    !Rf_isReal(SBBOn1)) {
    Rf_error("UpdateBBOn1: Can't work because SBBOn1 is defective!");
  }
	if (FixKa > 0) {
		SuccFlag = FixKCalculateBBOn();
		if (SuccFlag < 0) {
			Rprintf("####  ATwoLassoObject2011.cc:: UpdateBBOn1(), we got a return FixKCalculateBBOn is angry \n");
			R_FlushConsole(); return(-1);
		}
		return(SuccFlag);
	}
	if (Rf_length(SBBOn1) < p) {
    Rprintf("#### ATwoLassoObject2011.cc::UpdateBBOn1() SBBOn1 only has Rf_length %d, for p = %d", 
      Rf_length(SBBOn1), p); R_FlushConsole();
    Rf_error("No UpdateBBOn1 this cannot work!");
  }
  
  double LogitLambdaDiff = log( (1.0 - OnPiA)) -log(OnPiA) + log(OnLambdaD) -log(OnLambdaA);
  double LambdaDiff = (OnLambdaD - OnLambdaA);
  double ZeroProb = 1.0 / ( 1.0+ exp(LogitLambdaDiff) );
  if (CDO->XTXFlag == 2) {  
    for (ii = 0; ii < p; ii++) {
      REAL(SBBOn1)[ii] = ZeroProb;
    }
    int Onii = 0;
    for (ii = 0; ii < CDO->OnKappaS; ii++) {
      Onii = CDO->kXFinder[ii];
      REAL(SBBOn1)[Onii] = fabs(REAL(SOnBeta)[Onii]) * LambdaDiff;
      if ( REAL(SBBOn1)[Onii] - LogitLambdaDiff > 40 ) {
        REAL(SBBOn1)[Onii] = 1.0;
      } else if ( REAL(SBBOn1)[Onii] - LogitLambdaDiff < -15) {
        REAL(SBBOn1)[Onii] = exp ( REAL(SBBOn1)[Onii] - LogitLambdaDiff );
      } else {
        REAL(SBBOn1)[Onii] =  exp ( REAL(SBBOn1)[Onii] - LogitLambdaDiff );
        REAL(SBBOn1)[Onii] = REAL(SBBOn1)[Onii] / (REAL(SBBOn1)[Onii] + 1.0);
      }
    }
  } else {
    for (ii = 0; ii < p; ii++) {
      if (REAL(SOnBeta)[ii] == 0) {
        REAL(SBBOn1)[ii] = ZeroProb;
      } else {
        REAL(SBBOn1)[ii] = fabs(REAL(SOnBeta)[ii]) * LambdaDiff;
        if ( REAL(SBBOn1)[ii] - LogitLambdaDiff > 40 ) {
          REAL(SBBOn1)[ii] = 1.0;
        } else if ( REAL(SBBOn1)[ii] - LogitLambdaDiff < -15) {
          REAL(SBBOn1)[ii] = exp ( REAL(SBBOn1)[ii] - LogitLambdaDiff );
        } else {
          REAL(SBBOn1)[ii] =  exp ( REAL(SBBOn1)[ii] - LogitLambdaDiff );
          REAL(SBBOn1)[ii] = REAL(SBBOn1)[ii] / (REAL(SBBOn1)[ii] + 1.0);
        }
      }
    }
  }  

  if (NoShrinkColumns != NULL) {
   for (int iti = 0; iti < Rf_length(NoShrinkColumns->asSexp()); iti++) {
     REAL(SBBOn1)[INTEGER(NoShrinkColumns->asSexp())[iti]] = 1.0;
   }
  }
  if (Verbose > 4) {
      Rprintf("####  UpdateBBOn1(): Returning!\n",
        OnLambdaA, OnLambdaD); R_FlushConsole();
  }
	return(1);
}

int TwoLassoSexp::GLMRecordsNow(int tt1) {
    if (Verbose > 2) {
      Rprintf("ATwoLassoObject2011.cc:: GLMRecordsNow(tt1=%d, tt2=%d):", tt1, tt2); 
      Rprintf(" Implementing the Records Now, tt1 = %d \n", tt1);
      R_FlushConsole();
    }
    int jj=0;  int jjtt1 = tt1 *p ;
    if (RSRecOnBeta != NULL) {
       for (jj = 0; jj < p; jj++) {
         REAL(RSRecOnBeta->asSexp())[jjtt1] = GLMCDO->OnBeta[jj];
         jjtt1++;
       } 
    }
    if (RBeta0Records != NULL && tt1 >= 0 && !Rf_isNull(RBeta0Records->asSexp()) &&
      Rf_length(RBeta0Records->asSexp()) > tt1) {
      REAL(RBeta0Records->asSexp())[tt1] = GLMCDO->pBeta0[0];
    }
    if (RSRecBBOn1 != NULL) {
      jjtt1 = tt1 * p;
      for (jj=0; jj < p; jj++) {
        REAL(RSRecBBOn1->asSexp())[jjtt1] = REAL(SBBOn1)[jj];
        jjtt1++; 
      }
    }
    if (Verbose > 2) {
      Rprintf("ATwoLassoObject2011.cc:: GLMRecordsNow(tt1=%d, tt2=%d):", tt1, tt2); 
      Rprintf(" Done Recording \n"); R_FlushConsole();
    }
    return(1);
}
int TwoLassoSexp::OneRoundUpdateBeta() {
 if(Verbose > 4) {
    Rprintf("ATwoLassoObject2011.cc:: UpdateOneRoundBeta(tt1=%d, tt2=%d):", tt1, tt2); 
    Rprintf("  One Round tt1 = %d, tt2 = %d, Updating Beta\n",
      tt1, tt2); R_FlushConsole();
  }
  int SuccessFlag = 0;
  if (GLMCDO != NULL) {
    if (Verbose >= 3) {
      Rprintf("ATwoLassoObject2011.cc:: UpdateOneRoundBeta(tt1=%d, tt2=%d):", tt1, tt2); 
      Rprintf(" Setting TotMove and SetupGLMCDO \n"); R_FlushConsole();
    }
    GLMCDO->GLMTotMove = 0.0;
    this->SetupCDO(1);
    GLMCDO->MakeXTResid();
    if (Verbose >= 3) {
      Rprintf("ATwoLassoObject2011.cc:: UpdateOneRoundBeta(tt1=%d, tt2=%d):", tt1, tt2); 
      Rprintf(" About to Run GLMLasso \n"); R_FlushConsole();
    }
    SuccessFlag = GLMCDO->RunGLMLasso();
    if (GLMCDO->MemoryFailFlag == 1) {
      Rprintf("ATwoLAssoObject2011.cc:: GLMCDO is Memory Fail!\n");
      return(-1);
    }
    if (Verbose >= 4) {
      Rprintf("ATwoLassoObject2011.cc:: UpdateOneRoundBeta(tt1=%d, tt2=%d):", tt1, tt2);  
      Rprintf(" Finished Running GLMLasso returns SuccessFlag = %d\n",
       SuccessFlag); R_FlushConsole();
    }
    ////////////////////////////////////////////////////////
    //  E step - after CDO now take data			                             
    this->AfterCDO();	
    return(1);
  }
  
	CDO->TotMove = 0.0;
  if (CDO->ModifiedXTResid == 1) { 
    SetupCDO(1);  
    if (Verbose >= 2 || (p >= LARGEDATA && Verbose >= 1)) {
      Rprintf("ATwoLassoObject2011.cc:: UpdateOneRoundBeta(tt1=%d, tt2=%d): One Round now MakeXTresid. \n", tt1, tt2); R_FlushConsole();
    } 
		CDO->MakeXTResid();	
    if (Verbose >= 2 || (p >= LARGEDATA && Verbose >= 1)) {
      Rprintf("ATwoLassoObject2011.cc:: UpdateOneRoundBeta(tt1=%d, tt2=%d): One Round Finished XTResid. \n", tt1,tt2); R_FlushConsole();
    } 
	} else {
    SetupCDO(0);
    if (Verbose >= 2 || (p >= LARGEDATA && Verbose >= 1)) {
      Rprintf("ATwoLassoObject2011.cc:: UpdateOneRoundBeta(tt1=%d, tt2=%d): One Round Finished SetupCDO(0). \n", tt1,tt2); R_FlushConsole();
    } 
  }
  if (tt1 == 0 && Verbose > 4) {	
     Rprintf("ATwoLassoObject2011.cc:: UpdateOneRoundBeta(tt1=%d, tt2=%d): One Round PrintDiagonostic.\n", tt1, tt2);
     R_FlushConsole();
     TwoLassoPrintDiagnosticA();	
  }	
  
  SuccessFlag = CDO->RunAlgorithm(1, CauchyEpsilon);
  if (CDO->MemoryFailFlag == 1) {
    Rprintf("ATwoLassoObject2011.cc:OneRoundUpdateBeta(tt1=%d, tt2=%d), MemoryFails. A=%f, D=%f, OnPiA=%f\n",
      tt1, tt2, OnLambdaA, OnLambdaD, OnPiA);
    Rprintf("####                  OnLoop=%d, OnCoord=%d, TotMove=%f\n",
      CDO->OnLoop, CDO->OnCoord, CDO->TotMove);
    return(-1);
  }
	if (Verbose > 3 || (p >= LARGEDATA && Verbose >= 2)) {
		Rprintf((char*)"ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d) After CDO, OnKappaS=%d/%d/%d\n", 
      tt1,tt2, CDO->OnKappaS, CDO->OnKappaMem,p);   
		R_FlushConsole(); R_ProcessEvents(); 
	}  
  if (Verbose >= 4 || (p >= LARGEDATA && Verbose >= 3)) {
    Rprintf("ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d):  ", tt1, tt2);
    Rprintf(" Finish Coordinate Descent, OnLoop = %d, OnCoord = %d, TotMove=%f\n", 
      CDO->OnLoop, CDO->OnCoord, CDO->TotMove);
    R_FlushConsole();
  }
	if (SuccessFlag < 0) {
    Rprintf("ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d):", tt1, tt2);
		Rprintf(" Returned a Fail \n");
		R_FlushConsole(); R_ProcessEvents(); error("Oh No");
	}
  ////////////////////////////////////////////////////////
  //  E step - after CDO now take data			                             
  this->AfterCDO();	
  return(1);
}

int TwoLassoSexp::OneCoordUpdateBeta(int DoCoord) {
 if(Verbose > 4) {
    Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2,DoCoord);
    Rprintf("RunTwoLassoRegression: One Round tt1 = %d, tt2 = %d, DoCoord = %d, Updating Beta\n",
      tt1, tt2, DoCoord); R_FlushConsole();
  }
  if (DoCoord < 0 || DoCoord >= p) {
      Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
    Rf_error("OneCoordUpdateBeta: No DoCoord is %d\n", DoCoord); 
  }
  int SuccessFlag = 0;
  if (GLMCDO != NULL) {
    if (Verbose >= 3) {
      Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
      Rprintf(" Setting TotMove and SetupGLMCDO \n"); R_FlushConsole();
    }
    GLMCDO->GLMTotMove = 0.0; GLMCDO->TotMove = 0.0;
    this->SetupCDO(1);
    GLMCDO->MakeXTResid();
    GLMCDO->OnCoord = DoCoord;
    GLMCDO->UpDateCoord();
    if (Verbose >= 3) {
      Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
      Rprintf("  About to Run GLMLasso \n"); R_FlushConsole();
    }
    SuccessFlag = GLMCDO->RunGLMLasso();
    if (GLMCDO->MemoryFailFlag == 1) {
      Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
      Rprintf("####  Memory produced a failure on RunGLMLasso = %d\n",
        GLMCDO->MemoryFailFlag); R_FlushConsole();
    }
    if (Verbose >= 4) {
      Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
      Rprintf(" Finished Running GLMLasso returns SuccessFlag = %d\n",
       SuccessFlag); R_FlushConsole();
    }
    ////////////////////////////////////////////////////////
    //  E step - after CDO now take data			                             
    this->AfterCDO();	
    return(1);
  }
  
	CDO->TotMove = 0.0;
  CDO->OnCoord = DoCoord;
  if (CDO->ModifiedXTResid == 1) { 
    SetupCDO(1);  
    CDO->OnCoord = DoCoord;
    if (Verbose >= 3 || (p >= LARGEDATA && Verbose >= 2)) {
      Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
      Rprintf(" One Coordinate %d now MakeXTresid. \n", DoCoord); R_FlushConsole();
    } 
		CDO->MakeXTResid();	
    if (Verbose >= 3 || (p >= LARGEDATA && Verbose >= 2)) {
      Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
      Rprintf(" One Coordinate %d Finished XTResid. \n", DoCoord); R_FlushConsole();
    } 
	} else {
    SetupCDO(0);
    if (Verbose >= 3 || (p >= LARGEDATA &&Verbose >= 2)) {
      Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
      Rprintf(" One Coordinate %d Finished SetupCDO(0). \n", DoCoord); R_FlushConsole();
    } 
  }
  if (tt1 == 0 && Verbose > 4) {	
     Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
     Rprintf(" One Round PrintDiagonostic.\n");
     R_FlushConsole();
     TwoLassoPrintDiagnosticA();	
  }	
  SuccessFlag = CDO->UpDateCoord();
	if (Verbose > 3 || (p >= LARGEDATA && Verbose >= 2)) {
		Rprintf((char*)"TwoLassoSexp:: After CDO\n");   
		R_FlushConsole(); R_ProcessEvents(); 
	}  
  if (Verbose >= 4 || (p >= LARGEDATA && Verbose >= 2)) {
    Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
    Rprintf(" Finish Coordinate Descent, OnLoop = %f, OnCoord = %d, TotMove=%f\n", 
      CDO->OnLoop, CDO->OnCoord, CDO->TotMove);
    R_FlushConsole();
  }
	if (SuccessFlag < 0) {
    Rprintf("ATwoLassoObject2011.cc::OneCoordUpdateBeta(tt1=%d, tt2=%d, OnCoord=%d):", tt1, tt2, DoCoord);
		Rprintf(" UpdateBeta Returned a Fail \n");
		R_FlushConsole(); R_ProcessEvents(); error("Oh No");
	}
  ////////////////////////////////////////////////////////
  //  E step - after CDO now take data			                             
  this->AfterCDO();	
  return(1);
}


int TwoLassoSexp::UpdateBeta() {
  if (Verbose > 4) {
    Rprintf("****************************UpdateBeta**********************************************\n");
    Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d, OnCoord=%d): OnKappaS=%d/%d/%d", tt1, tt2, CDO->OnCoord,
      CDO->OnKappaS, CDO->OnKappaMem, p);
    Rprintf(" On The Start. \n",
      tt1, tt2); R_FlushConsole();
  }
  int SuccessFlag = 0;
  if (GLMCDO != NULL) {
    if (Verbose >= 3) {
      Rprintf("ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d, OnCoord=%d): OnKappaS=%d/%d/%d", tt1, tt2, CDO->OnCoord,
        CDO->OnKappaS, CDO->OnKappaMem, p);
      Rprintf(" Setting TotMove and SetupGLMCDO \n"); R_FlushConsole();
    }
    GLMCDO->GLMTotMove = 0.0; GLMCDO->TotMove = 0.0;
    this->SetupCDO(1);
    GLMCDO->MakeXTResid();
    if (Verbose >= 3) {
      Rprintf("**** UpdateBeta(tt1=%d,tt2=%d): About to Run GLMLasso \n", tt1, tt2); R_FlushConsole();
    }
    SuccessFlag = GLMCDO->RunGLMLasso();
    if (SuccessFlag < 0 || GLMCDO->SuccessFlag < 0) {
      Rprintf("ATwoLassoObject2011.cc::UpdateBeta(tt1=%d,tt2=%d), SuccessFlag error = %d\n", tt1, tt2, SuccessFlag);
      Rprintf(" -----------------------------------------------------------------------------------------------\n");
      Rprintf(" --  This is an Error, OnLambdaA=%d, OnLambdaD=%f, OnCoord=%d, OnKappaS=%d/%d/%d. \n", OnLambdaA,
        OnLambdaD, GLMCDO->OnCoord, GLMCDO->OnKappaS, GLMCDO->OnKappaMem, p);
      Rprintf(" -- TotMove was %f \n", GLMCDO->TotMove); R_FlushConsole();
      Rprintf(" --  OnSigma is %f \n", OnSigma);
      Rprintf(" -- Thiw as an error we return -666");
      Rprintf(" --  UpdateBeta(ATwoLassoObject2011.cc) Bad Error. \n"); R_FlushConsole();
      Rprintf(" -----------------------------------------------------------------------------------------------\n"); R_FlushConsole();
      return(-666);
    }
    if (Verbose >= 4) {
      Rprintf("**** UpdateBeta: Finished Running GLMLasso returns SuccessFlag = %d\n",
       SuccessFlag); R_FlushConsole();
    }
    ////////////////////////////////////////////////////////
    //  E step - after CDO now take data			                             
    this->AfterCDO();	
    return(1);
  }
  
	CDO->TotMove = 0.0;
  if (CDO->ModifiedXTResid == 1) { 
    SetupCDO(1);  
    if (Verbose >= 3 || (p >= LARGEDATA && Verbose >= 2)) {
      Rprintf("**** UpdateBeta(tt1=%d, tt2=%d): now MakeXTresid. \n", tt1, tt2); R_FlushConsole();
    } 
		CDO->MakeXTResid();	
    if (Verbose >= 3 || (p >= LARGEDATA && Verbose >= 2)) {
      Rprintf("**** UpdateBeta(tt1=%d, tt2=%d): Finished XTResid. \n", tt1, tt2); R_FlushConsole();
    } 
	} else {
    SetupCDO(0);
    if (Verbose >= 3 || (p >= LARGEDATA && Verbose >= 2)) {
      int CNZ = 0;
      double MBBOn1 = 0.0;
      double MOnGamma = 0.0;
      for (int iii = 0; iii < p; iii++) {
        if (CDO->OnBeta[iii] != 0.0) {
          CNZ++;
        }
        if (R_isnancpp(CDO->OnBeta[iii])) {
          Rf_error("**** UpdateBeta(tt1=%d, tt2=%d), SetupCDO, CDO->OnBeta[iii=%d] is a NAN! \n", tt1, tt2, iii);
        }
        if (!R_finite(CDO->OnBeta[iii])) {
          Rf_error("**** UpdateBeta(tt1=%d, tt2=%d), SetupCDO, CDO->OnBeta[iii=%d] is %f! \n", tt1, tt2, iii,
            CDO->OnBeta[iii]);
        }
        MBBOn1 += REAL(SBBOn1)[iii] * (1.0/p);
        if (CDO->OnGammas != NULL) {
          MOnGamma += CDO->OnGammas[iii] * (1.0/p);
        }
      }
      Rprintf("**** UpdateBeta(tt1=%d, tt2=%d): LocAAAA Finished SetupCDO(0). OnKappaS=%d/%d/%d, NonZero = %d.  Mean BBOn1 = %f, Mean OnGamma = %f.  PiA=%.8f\n", tt1, tt2,
        CDO->OnKappaS, CDO->OnKappaMem, p, CNZ, MBBOn1, MOnGamma, OnPiA); R_FlushConsole();
    } 
  }
  if ((tt1 == 0 && Verbose > 4)  || jti == -40) {	
     Rprintf("**** UpdateBeta(tt1=%d, tt2=%d): PrintDiagonostic.\n", tt1, tt2);
     R_FlushConsole();
     TwoLassoPrintDiagnosticA();	
  }	
  if (Verbose > 4 || (p >= LARGEDATA && Verbose >= 3) || (p >= LARGEDATA && jti == -40)) {
    Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d): ", tt1, tt2);
    Rprintf("  about to run coordinate descent: tt1=%d, MaxCauchy=%d, CauchyEpsilon = %f, OnLambdaA=%f, OnLambdaD=%f\n",
      tt1, MaxCauchy, CauchyEpsilon, OnLambdaA, OnLambdaD);
    R_FlushConsole();
    Rprintf("**** p = %d, CDO->OnKappaS=%d, CDO->OnKappaMem = %d, good luck!\n", p, CDO->OnKappaS, CDO->OnKappaMem);
    R_FlushConsole();
    Rprintf("**** Let's check Max Size XTX!\n"); R_FlushConsole();
    if (CDO->OnKappaS==0 && CDO->OnKappaMem == p) {
      Rprintf("**** -Issue- This is absolutely ludicrous, CDO->OnKappaS = 0, but CDO->OnKappaMem = p=%d, nothing is active!\n", CDO->OnKappaMem);
    }
    int cOnS = 0;
    for (int ii = 0; ii < p; ii++) {
      if (CDO->OnBeta[ii] != 0.0) {
        cOnS++;
      }
    }
    Rprintf("**** Count of Nonzero OnBeta = %d\n", cOnS);  R_FlushConsole();
    Rprintf("**** Attempting to trigger last element of pXTX, OnKappaS = %d, OnKappaMem = %d\n", CDO->OnKappaS, CDO->OnKappaMem); R_FlushConsole();
    if (CDO->pXTX == NULL) {
      Rf_error("Error: pXTX is NULL!\n");
    }
    if (CDO->OnKappaMem <= 0) {
      Rprintf("** ATwoLassoObject2014.cc: UpdateBeta, error, OnKappaMem = %d, not good!\n", CDO->OnKappaMem); R_FlushConsole();
      Rf_error("** OnKappaMem based not set right error. \n");
    }
    if (CDO->pXTX[CDO->OnKappaS-1] == NULL)  {
      Rprintf("** ATwoLassoObject2014.cc: UpdateBeta at pXTX[OnKappaS-1] for OnKappaS = %d, The pXTX is null!\n", CDO->OnKappaS);
      R_FlushConsole();
      Rf_error("** OnKappaMem has an error. \n");
    }
    if (CDO->OnKappaS >= 1) {
      if (CDO->pXTX[CDO->OnKappaS-1] == NULL) {
        Rf_error("ATwoLassoObject2014.cc:: Error, OnKappaS = %d, but pXTX[%d-1] = NULL!\n",
          CDO->OnKappaS, CDO->OnKappaS);
      }
      Rprintf("****  pXtX[(OnKappaS=%d)-1]+((p=%d)-1)) = %f\n",
        CDO->OnKappaS, p, *(CDO->pXTX[CDO->OnKappaS-1]+p-1)); R_FlushConsole();
      Rprintf("****    CDO->OnKappaS=%d, CDO->kXFinder[OnKappaS-1=%d]=%d\n",
      CDO->OnKappaS, CDO->OnKappaS-1, CDO->kXFinder[CDO->OnKappaS-1]); R_FlushConsole();
    } else {
      Rprintf("**** Doesn't matter, OnKappaS is zero anyway. \n");
      Rprintf("****    CDO->OnKappaS=%d. \n",
      CDO->OnKappaS); R_FlushConsole();
    }
    if (CDO->OnKappaMem <= 0) {
      Rf_error("  Can't tell you what pXtX is because OnKappaMem = %d!\n", CDO->OnKappaMem);
    }
    Rprintf("***     *(pXtX[(OnKappaMem=%d)-1] + ((p=%d)-1)) = %f\n",
      CDO->OnKappaMem, p,  *(CDO->pXTX[CDO->OnKappaMem-1] + p-1)); R_FlushConsole();
    Rprintf("****    XtResid[%d-1=%d] = %f\n", CDO->kLen, CDO->kLen-1, CDO->XTResid[CDO->kLen-1]); R_FlushConsole();
    Rprintf("****    CDO->OnGammas[0] = %f,\n", CDO->OnGammas[0]); R_FlushConsole();
    Rprintf("****    and CDO->OnGammas[%d] = %f\n", CDO->kLen-1, CDO->OnGammas[CDO->kLen-1]); R_FlushConsole();

    //CDO->PrintFlag =3;
    if (CDO->OnRecordBetas != NULL) {
      Rprintf("**** Note that OnRecordBetas is not NULL!\n"); R_FlushConsole();
    } else {
      Rprintf("**** -Issue- But OnRecordBetas is NULL\n"); R_FlushConsole();
    }
    CDO->TriggerReallocFlag = 1;  CDO->SuccessDOP = 0;
    CDO->TriggerInvestigateFlag = 1;
    Rprintf("****  Run One  loop of Algorithm settting checkers!\n"); R_FlushConsole();
      SuccessFlag = CDO->RunAlgorithm(0, .01);
    //CDO->PrintFlag = Verbose-2;
    Rprintf("****  Algorithm was run once, now more times. SuccessDOP = %d \n", CDO->SuccessDOP); R_FlushConsole();
    CDO->TriggerReallocFlag = 1;  CDO->SuccessDOP = 0;  CDO->TriggerInvestigateFlag = 0;
    if (CDO->MemoryFailFlag == 1) {
      Rprintf("****  ATwoLassoObject2011.cc::UpdateBeta() In examine ERROR ERROR ERROR ERROR\n");
      Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta() Algorithm Was run with a Memory Failure!\n");
      Rprintf("**** CDO->OnLoop=%d, OnKappaS=%d/%d/%d, A=%f, D=%f\n",
        CDO->OnLoop, CDO->OnKappaS, CDO->OnKappaMem, p, OnLambdaA, OnLambdaD);
      R_FlushConsole();
      return(-1);
    }
      SuccessFlag = CDO->RunAlgorithmDiagonostic(MaxCauchy, CauchyEpsilon);   
    if (Verbose >= 4 || (p >= LARGEDATA && Verbose >= 3)) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d): ", tt1, tt2);
      Rprintf("**** Just Finish CDO->RunAlgorithm, OnLoop = %d, OnCoord=%d\n",
        CDO->OnLoop, CDO->OnCoord); R_FlushConsole();
    } 
  } else{
    if (Verbose >= 3) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d: LAAAB About to Run CDO Algorithm(%d,%f)\n", tt1, tt2, MaxCauchy, CauchyEpsilon);
      R_FlushConsole();
    }
     SuccessFlag = CDO->RunAlgorithm(MaxCauchy, CauchyEpsilon);
    if (Verbose >= 3) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d: LAAAC Completed Run CDO Algorithm(%d,%f)\n", tt1, tt2, MaxCauchy, CauchyEpsilon);
      R_FlushConsole();
    }
  }
  if (SuccessFlag < 0 || CDO->MemoryFailFlag == 1) {
    Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta, CDO->RunAlgorithm returned a memory fail. \n"); R_FlushConsole();
    return(-1);
  }
	if (Verbose > 3 || (p >= LARGEDATA && Verbose >= 3)) {
    Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d): ", tt1, tt2);
		Rprintf((char*)"TwoLassoSexp:: After CDO\n");   
		R_FlushConsole(); R_ProcessEvents(); 
	}  
  if (Verbose >= 4 || (p >= LARGEDATA && Verbose >= 3)) {
    Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d): ", tt1, tt2);
    Rprintf("TwoLassoSexp: Finish Coordinate Descent, OnLoop = %f, OnCoord = %d, TotMove=%f\n", 
      CDO->OnLoop, CDO->OnCoord, CDO->TotMove);
    R_FlushConsole();
  }
	if (SuccessFlag < 0) {
    Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d): ERROR ISSUE ", tt1, tt2);
		Rprintf(" -- UpdateBeta Returned a Fail \n");
		R_FlushConsole(); R_ProcessEvents(); Rf_error(" - SuccessFlag = %d \n", SuccessFlag);
	}
  if (jti == -40) {
    Rprintf("**** ATwoLassoObject2011.cc::Update(tt1=%d, tt2=%d): LAAAD Enact After CDO\n", tt1, tt2);
    R_FlushConsole();
  }
  ////////////////////////////////////////////////////////
  //  E step - after CDO now take data			                             
  this->AfterCDO();	
  if (Verbose >= 4 || (p >= LARGEDATA && Verbose >= 3)) {
    Rprintf("**** ATwoLassoObject2011.cc::UpdateBeta(tt1=%d, tt2=%d):  ALL Completed \ns", tt1, tt2);  R_FlushConsole();
    Rprintf("******************************************UpdateBeta()********************\n"); 
    R_FlushConsole();
  }
  if (jti == -40) {
    Rprintf("**** ATwoLassoObject2011.cc::Update(tt1=%d, tt2=%d): LAAAE Enact Completed AfterCDO, All UpdateBeta completed.\n", tt1, tt2);
    R_FlushConsole();
  }
  return(1);
}
/////////////////////////////////////////////////////////////////////////////
//  UpdateSigmaSq()
// We're really putting into sigmaNoiseSq 1 / E[ 1 / sigmaNoiseSq | Beta ]
// AKA We don't need sigmaNoiseSq itself, but the expected inverse.
int TwoLassoSexp::UpdateSigmaSq() {
  if(Verbose > 4) {
    Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): Updating SigmaSq\n",
      tt1, tt2); R_FlushConsole();
  }
  if (SigmaBar <= 0) {
    return(0);
  }
	//if (Verbose > 3) {
  //  Rprintf((char*)"UpdateSigmaSq: this->SumCurResids = %.4f \n", 
	//    (double) this->OnLARS->SumCurResids);
	//  R_FlushConsole();
  //}

  int ii;
  NonZero = 0; NonZeroG = 0; CurSBBar= 0.0; CountGroupSq = 0.0;
  if (NumGroups <= 0) {
  for (ii = 0; ii < p; ii++) {
	  CurSBBar += REAL(SBBOn1)[ii]; 
    if (CDO->OnBeta[ii] != 0.0) {
      NonZero++;
    }
  }
  } else if (NumGroups > 0) {
    for (ii = 0; ii < FirstRandomIndex; ii++) {
	    CurSBBar += REAL(SBBOn1)[ii]; 
      if (CDO->OnBeta[ii] != 0.0) {
        NonZero++;
      }
    }
    /*
    int St = FirstRandomIndex;
    for (ii = 0; ii < NumGroups; ii++) {
      for (int jj = St; jj < EndGroupLocations[ii]+1; jj++) {
        if (CDO->OnBeta[jj] > 0.0) {
          NonZeroG++;
          jj = EndGroupLocations[ii]+1;
          break;
        }
      }
      St = EndGroupLocations[ii]+1;
    } */
    for (ii = 0; ii < NumGroups; ii++) {
      if (REAL(GroupsBBOn1->asSexp())[ii] > 0.0 &&
        REAL(GroupsBBOn1->asSexp())[ii] <= 1.0) {
        NonZeroG += (double) REAL(GroupsBBOn1->asSexp())[ii];
      }
    }
    int St = FirstRandomIndex;
    double Adder = 0.0;
    for (ii = 0; ii < NumGroups; ii++) {
       for (int jj = St; jj < EndGroupLocations[ii]+1; jj++) {
         if (CDO->OnBeta[jj] != 0.0 && CDO->XLC != NULL && CDO->XLC[jj] >= 0 &&
           CDO->XLC[jj] <= p) {
           if (CDO->XLC[jj] >= CDO->OnKappaMem) {
             Rf_error("Error GroupsBBOns1, jj = %d, CDO->XLC[%d] = %d, OnKappaMem = %d!\n",
               jj, jj, CDO->OnKappaMem);
           }
           if (CDO->pXTX[CDO->XLC[jj]] == NULL) {
             Rf_error("Error, GroupsBBOn1, jj = %d, XLC[%d]=%d, but pXTX[%d] is NULL!\n",
              jj, jj, CDO->XLC[jj], CDO->XLC[jj]); 
           }
           if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[CDO->XLC[jj]] <= 0) {
             CDO->ReweightCoordinate(jj);
           }
           Adder = CDO->XTResid[jj] + *(CDO->pXTX[CDO->XLC[jj]] + jj) * CDO->OnBeta[jj];
           if (R_finite(Adder) && Adder != 0.0) {
             Adder = *(CDO->pXTX[CDO->XLC[jj]] + jj) * CDO->OnBeta[jj] / Adder;
             CountGroupSq += Adder * Adder; 
           }
         } else if (CDO->OnBeta[jj] != 0.0) {
           Adder = 0.0;
           if (CDO->iiWeights != NULL) {
           for (int iti = 0; iti < n; iti++) {
             Adder += CDO->XX[jj * n + iti] * CDO->XX[jj*n + iti] * CDO->iiWeights[iti];
           }
           } else {
           for (int iti = 0; iti < n; iti++) {
             Adder += CDO->XX[jj * n + iti] * CDO->XX[jj*n + iti];
           }
           }
           Adder = CDO->XTResid[jj] + Adder * CDO->OnBeta[jj];
           if (R_finite(Adder) && Adder != 0.0) {
             Adder = CDO->OnBeta[jj] / Adder;
             CountGroupSq += Adder * Adder;
           }
         }
       }
       St = EndGroupLocations[ii]+1;
    }
    if (R_isnancpp(CountGroupSq)) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdatesigmaSq(tt1=%d, tt2=%d) Error, CountGroupSq is NAN!\n");
      R_FlushConsole();
    }
    if (!R_finite(CountGroupSq)) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d, tt2=%d) Error CountGroupSq is not finite is %f \n", CountGroupSq);
      R_FlushConsole();
    }
  }
  

  if (Verbose > 3) {
    Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d):", tt1, tt2);
    Rprintf(" CurSBBar = %f \n", CurSBBar); R_FlushConsole();
  }
  double TotalResids = 0.0;
  double TempSumYYSq;  double SumRDSquared;
  int One = 1;   int Zero = 0; double OneD = 1; double ZeroD = 0.0;
  double SumWeightsii = 0;
  double OtherSumRDS = 0; int ii1 = 0; int jj1 = 0;
  double OnYResid =0;  int MarkMe = 0;
  if (MarkMe == -1  || OnYResid == -1 || ii1 == -1 || jj1 == -1) {
  }
     
  ii1 = 1; jj1 = 1;OtherSumRDS++; MarkMe = 0; OnYResid = 0;
  if (TDFNu > 0 && this->CDO != NULL) {
    if (Verbose >= 3) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d):  ", tt1, tt2);
      Rprintf(" TDFNu = %f, CDO is not null doing p =%d > n=%d version. \n",
        TDFNu, p, n); R_FlushConsole();
    }
    if (iiWeights == NULL) {
      Rprintf("ERROR ERROR ERROR ERROR On iiWeights \n"); R_FlushConsole();
      Rprintf("**** ERROR ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): ", tt1, tt2);
      Rprintf(" Error in UpdateSigmaSq: TDFNu = %f, but iiWeights is NULL\n", iiWeights);
      Rf_error("UpdatesigmaSq, TDFNu=%f!", TDFNu);
    }
    TempSumYYSq = 0;
    for (ii = 0; ii < n; ii++) {
      TempSumYYSq += CDO->YY[ii] * CDO->YY[ii] * iiWeights[ii];
    }
    SumRDSquared =   TempSumYYSq;
    //if (CDO->XTXFlag == 2) {
    //  for (ii = 0; ii < CDO->OnKappaS; ii++) {
    //    SumRDSquared -= CDO->OnBeta[CDO->kXFinder[ii]] * 
    //      (CDO->XTResid[CDO->kXFinder[ii]] + 
    //       CDO->XTY[CDO->kXFinder[ii]]);
    //  }
    //} else {
    SumRDSquared -= F77_CALL(ddot)(&p, CDO->OnBeta, &One, CDO->XTResid, &One);
    SumRDSquared -= F77_CALL(ddot)(&p, CDO->OnBeta, &One, CDO->XTY, &One);
    //}
    SumWeightsii = 0; 
    for (ii = 0; ii < n; ii++) {
      SumWeightsii += iiWeights[ii];
    }
     
    if (Verbose >= 3) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d):  ", tt1, tt2);
      Rprintf(" OtherSumRDS = %f, SigmaBar=%f, SigmaDF = %f, SumWeightsii = %f, SigmaDf = %f\n",
        OtherSumRDS, SigmaBar, SigmaDf, SumWeightsii, SigmaDf); R_FlushConsole();
    }
    //this->OnSigma = (SumRDSquared + SigmaBar * SigmaDf) / (SumWeightsii + SigmaDf);  
    
    if (p < 100) {
     this->OnSigma = (SumRDSquared + SigmaBar * SigmaDf) / 
       (this->n + SigmaDf);
    } else {
     if (Verbose >= 3) {
       Rprintf("**** ATwoLassoObject2011.cc:::UpdateSigmaSq(%d,%d): GroupLargeVersion::: CountGroupSq=%f, CurSBBar=%f, NonZero=%f, NonZeroG=%f!\n",
         CountGroupSq, CurSBBar, NonZero, NonZeroG); R_FlushConsole();
     }
     if (CountGroupSq + CurSBBar > NonZero+NonZeroG) {
       if (CountGroupSq + CurSBBar > this->n + SigmaDf) {
         this->OnSigma = (SumRDSquared + SigmaBar * SigmaDf) / ((double) SigmaDf);
       }  else {
         this->OnSigma = (SumRDSquared + SigmaBar * SigmaDf) / 
           ( (double) this->n + (double) SigmaDf- (double) CountGroupSq - CurSBBar);      
       }
     } else if (NonZero+NonZeroG > this->n + SigmaDf) {
       this->OnSigma = (SumRDSquared+SigmaBar*SigmaDf) /
         (SigmaDf);
     } else {
       this->OnSigma = (SumRDSquared+SigmaBar*SigmaDf) /
         (this->n + SigmaDf - NonZero - NonZeroG);
     }
    }
    CDO->TDFSigma = this->OnSigma;
    if (Rf_isNull(SSigma) || Rf_isString(SSigma)) {
      if (Rf_isNull(SSigma)) {
        Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): ", tt1, tt2);
        Rprintf(" Weird: SSigma is null here!\n");  R_FlushConsole();
      }
      //if (Rf_isChar(SSigma)) {
      //  Rprintf("UpdateSigmasq: weird: SSigma is char here!\n");
      //}
      if (Rf_isString(SSigma)) {    
        Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): ", tt1, tt2);
        Rprintf("  Weirder:  SSigma is string here!\n");
      }
      R_FlushConsole();
      if (RSSigma != NULL) { DDelete(RSSigma, "RSSigma"); }
      RSSigma = new AObject(Rf_allocVector(REALSXP,1));
      SSigma = RSSigma->asSexp();
    }
    if (RSSigma == NULL || Rf_isNull(SSigma) || !Rf_isReal(SSigma) ||
      Rf_length(SSigma) <= 0) {
    } else if (Rf_length(SSigma) > tt1+1) {
      REAL(SSigma)[tt1] = this->OnSigma;
    } else if (Rf_length(SSigma) > tt1+2) {
      (REAL(SSigma))[tt1+1] = this->OnSigma;
    } else if (Rf_length(SSigma) == 1) {
      REAL(SSigma)[0] = this->OnSigma;
    }
    if (Verbose >= 3) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): ", tt1, tt2);
      Rprintf(" Finished OtherSumRDS version, Sigma = %f, OtherSumRDS = %f!\n",
        this->OnSigma, OtherSumRDS); R_FlushConsole();
    }
     return(1);
  } else if (this->CDO != NULL) {
	  if (Verbose > 3) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): CDO Not Null, no TDFNoise, ", tt1, tt2);
	    Rprintf((char*) "****                     Trying to UpdateSigmaSq, with SumYYSq = %.4f, CDO != NULL\n", (double) SumYYSq);
	    R_FlushConsole(); R_ProcessEvents();
    }
    TotalResids = SumYYSq;
    if (Verbose > 3) {
      Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): ", tt1, tt2);
      Rprintf(" SumYYSq = %f \n", SumYYSq); R_FlushConsole();
    }
    if (CDO->XTXFlag == 1) {
	    if (Verbose > 3) {
        Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): ", tt1, tt2);
	      Rprintf((char*) "Trying to UpdateSigmaSq, XTXFlag == 1, SumYYSq = %f\n", 
          (double) SumYYSq);
	      R_FlushConsole(); R_ProcessEvents();	  
      }   
      for (ii = 0; ii < p; ii++) {
        if (CDO->OnBeta != NULL && CDO->OnBeta[ii] != -999.00 && 
          CDO->OnBeta[ii] != 0.0) {
          TotalResids -= CDO->OnBeta[ii] * (CDO->XTResid[ii] + CDO->XTY[ii]);
        }
      }
      if (Verbose >= 3) {
        Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d): ", tt1, tt2);
        Rprintf(" After update we have TotalResids = %f \n", TotalResids); R_FlushConsole();
      }
    } else  if (CDO->XTXFlag == 2) {
	    if (Verbose > 3) {
        Rprintf("**** ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d):  ", tt1, tt2);
	      Rprintf((char*) ": Trying to UpdateSigmaSq, XTXFlag == 2, so version 2\n" );
	      R_FlushConsole(); R_ProcessEvents();	     
      }
      for (ii = 0; ii < CDO->OnKappaS; ii++) {
        if (CDO->OnBeta != NULL && CDO->OnBeta[CDO->kXFinder[ii]] != -999 &&
          CDO->OnBeta[CDO->kXFinder[ii]] != 0.0) {
          TotalResids -= CDO->OnBeta[CDO->kXFinder[ii]] * 
            (CDO->XTResid[CDO->kXFinder[ii]] + CDO->XTY[CDO->kXFinder[ii]]);
        }
      }
      if (Verbose >= 3) {
        Rprintf("****      with OnKappaS = %d, we got TotalResids = %f \n",
          CDO->OnKappaS, TotalResids); R_FlushConsole();
      }
    }     
    if (TotalResids < 0) {
      Rprintf("**** ERROR ATwoLassoObject2011.cc::UpdateSigmaSq(tt1=%d,tt2=%d):  ", tt1, tt2);
      Rprintf("  ERROR EstimateSigmaSq, CDO Error total resids ");
      Rprintf(" is less than zero? = %.4f\n", TotalResids);
      Rprintf("\n****E TwoLassoSexp:: Current CDO->OnBeta\n");
       PrintVector(CDO->OnBeta, p);
      Rprintf("\n****E TwoLassoSexp:: Current CDO->OnGamma = %.4f, CDO->OnGammas\n", 
        CDO->OnGamma);
      Rprintf("\n****E TwoLassoSexp:: OnGammas = "); PrintVector(CDO->OnGammas, p);
      Rprintf("\n****E TwoLassoSexp:: Current OnBetas\n");
        Rprintf("\n****E OnBetas = "); PrintVector(REAL(SOnBeta), p); R_FlushConsole();
      Rprintf("\n****E TwoLassoSexp:: CurrentXTResid\n");
        Rprintf("\n****E XTResid = "); PrintVector(CDO->XTResid, p);
      Rprintf("\n****E SumYYSq = %f \n and XtY = ", SumYYSq); R_FlushConsole();
      PrintVector(CDO->XTY, p); R_FlushConsole();
      double *ErrorFind = Calloc(p, double);
      int jj2;
      for (jj2 = 0; jj2 < p; jj2++) { ErrorFind[jj2] = CDO->XTResid[jj2]; }
      R_FlushConsole(); R_ProcessEvents();
      SetupSumYYSq();
      if (Sxxs != NULL && Syys != NULL && CDO->XTY != NULL) {
        tMatTimesVec(p, n,  REAL(Sxxs), REAL(Syys), CDO->XTY); 
      }
      Rprintf("\n****E After that Error: Making XTResid for Residuals \n");
        CDO->MakeXTResid();
      for (jj2 = 0; jj2 < p; jj2++) { ErrorFind[jj2] -= CDO->XTResid[jj2]; }
      Rprintf("****E Pre-PostCalculation Errors are: \n"); R_FlushConsole();
      PrintVector(ErrorFind, p);  Free(ErrorFind); ErrorFind = NULL;
      Rprintf("****E TwoLassoSexp:: We're going to try this one more time\n"); 
      R_FlushConsole();
      TotalResids = SumYYSq;
      if (CDO->XTXFlag == 1) {
  	    if (Verbose > 3) {
  	      Rprintf((char*) "****E Trying to updateSigmaSq, XTXFlag == 1\n");
  	      R_FlushConsole(); R_ProcessEvents();	  
        }   
        for (ii = 0; ii < p; ii++) {
          if (CDO->OnBeta != NULL && CDO->OnBeta[ii] != -999 && 
            CDO->OnBeta[ii] != 0.0) {
            TotalResids -= CDO->OnBeta[ii] * (CDO->XTResid[ii] + CDO->XTY[ii]);
          }
        }
     } else  if (CDO->XTXFlag == 2) {
	     if (Verbose > 3) {
	       Rprintf((char*) "****E Trying to updateSigmaSq, XTXFlag == 2\n", 
           (double) SumYYSq);
	       R_FlushConsole(); R_ProcessEvents();	     
       }
       for (ii = 0; ii < CDO->OnKappaS; ii++) {
         if (CDO->OnBeta != NULL && CDO->OnBeta[CDO->kXFinder[ii]] != -999 &&
           CDO->OnBeta[CDO->kXFinder[ii]] != 0.0) {
           TotalResids -= CDO->OnBeta[CDO->kXFinder[ii]] * 
             (CDO->XTResid[CDO->kXFinder[ii]] + CDO->XTY[CDO->kXFinder[ii]]);
         }
       }
     }    
	   if (TotalResids < 0) {
		   Rprintf("****EE TwoLassoSexp:: EstimateSigmaSq, CDO Error total resids is ");
       Rprintf(" less than zero? = %.4f\n", TotalResids);
		   R_FlushConsole(); R_ProcessEvents();      
		   Rprintf((char*) "****EE SumYYSq = %.5f, Printing XTY\n", SumYYSq); 
		   PrintVector(CDO->XTY, p);
		   Rprintf((char*) "****EE Printing XTResid\n"); 
		   PrintVector(CDO->XTResid, p); R_FlushConsole();       
		     
		   Rprintf((char*) "****EE OnBeta: \n");
		   PrintVector(CDO->OnBeta, p);
		   Rprintf((char*) "****EE OnKappaS = %d, the kXFinder: \n", CDO->OnKappaS);
		   Rprintf((char*) "****EE TotalResid Pretend Calc = %.4f\n", TotalResids);
		   PrintVector(CDO->kXFinder, CDO->OnKappaS);
		   int jj; double Resid;
		   TotalResids = 0;
		   if (Sxxs != NULL) {
		     for (ii = 0; ii < n; ii++) {
			     Resid = REAL(Syys)[ii];
			     for (jj = 0; jj < CDO->OnKappaS; jj++) {
				     Resid -= REAL(Sxxs)[ CDO->kXFinder[jj] * n + ii] * 
				       CDO->OnBeta[CDO->kXFinder[jj]];
			      }
			      TotalResids += Resid*Resid;
		      }
	     } else {
		     Rprintf("****EE Needed RealXXS to not be NULL, failed \n");
         Rprintf("EEEE Error Return UpdateSigmaSq() Error \n");
		     R_FlushConsole(); R_ProcessEvents();
		     SuccessFlag = -1; return(-1);
	     }
		   Rprintf((char*) "****EE Calculated it Manually and it is %.4f\n", TotalResids);
       Rprintf("EEEE Error Return UpdateSigmaSq() Error \n");
		   R_FlushConsole(); R_ProcessEvents();
		   SuccessFlag = -1; return(-1);
	   }
   } else if (Verbose >= 3) {
     Rprintf("**** ATwoLassoObject2011.cc:::UpdateSigmaSq(%d,%d):  TotalResids=%f, passes the test!\n", tt1, tt2, TotalResids);
     R_FlushConsole();
   }

   if (p < 100) {
     this->OnSigma = (TotalResids + SigmaBar * SigmaDf) / 
       (this->n + SigmaDf);
   } else {
     if (Verbose >= 3) {
       Rprintf("**** ATwoLassoObject2011.cc:::UpdateSigmaSq(%d,%d): GroupLargeVersion::: CountGroupSq=%f, CurSBBar=%f, NonZero=%f, NonZeroG=%f!\n",
         CountGroupSq, CurSBBar, NonZero, NonZeroG); R_FlushConsole();
     }
     if (CountGroupSq + CurSBBar > NonZero+NonZeroG) {
       if (CountGroupSq + CurSBBar > this->n + SigmaDf) {
         this->OnSigma = (TotalResids + SigmaBar * SigmaDf) / ((double) SigmaDf);
       }  else {
         this->OnSigma = (TotalResids + SigmaBar * SigmaDf) / 
           ( (double) this->n + (double) SigmaDf- (double) CountGroupSq - CurSBBar);      
       }
     } else if (NonZero+NonZeroG > this->n + SigmaDf) {
       this->OnSigma = (TotalResids+SigmaBar*SigmaDf) /
         (SigmaDf);
     } else {
       this->OnSigma = (TotalResids+SigmaBar*SigmaDf) /
         (this->n + SigmaDf - NonZero - NonZeroG);
     }
   }
   if (Verbose >= 3) {
     Rprintf("**** UpdateSigma(%d,%d): OnSigma = %f, SigmaBar = %f, SigmaDf = %f, TotalResids = %f, n = %d\n",
       tt1, tt2, OnSigma, SigmaBar, SigmaDf, TotalResids, n); R_FlushConsole();
   }
   if (this->OnSigma < 0 && SigmaBar > 0) { this->OnSigma = SigmaBar; }
   if (R_isnancpp(this->OnSigma)) {
     Rprintf("**** ATwoLassoObject2011.cc:::OnSigma, not a good day, OnSigma was set to NAN!\n"); R_FlushConsole();
   }
   if (!R_finite(this->OnSigma)) {
     Rprintf("**** ATwoLassoObject2011.cc:::OnSigma, not a good day, OnSigma was set to non finite!\n"); R_FlushConsole();
   }
   if (RSSigma != NULL && !Rf_isNull(SSigma) && Rf_isReal(SSigma)) {
     if (Verbose >= 3) {
       Rprintf("**** UpdateSigma(%d,%d), OnSigma=%f, SSigma is not null with length %d, tt1 = %d, insertion time.  \n",
         tt1, tt2, this->OnSigma, Rf_length(SSigma), tt1); R_FlushConsole();
     }
     int WeGo = 0;
     if (Rf_length(SSigma) > tt1 +1) {
       if (Verbose >= 3) {
         Rprintf("**** UpdateSigma(%d,%d): OnSigma=%f we have tt1 (tt1+1) will put in SSigma, length=%d\n", tt1, tt2, OnSigma, Rf_length(SSigma)); R_FlushConsole();
         Rprintf("**** Putting into tt1=%d SSigma, this->OnSigma = %f, NonZero=%d, NonZeroG=%f, TotalResids=%f, CurSBBar = %f, CountGroupSq = %f\n", tt1, this->OnSigma,
           NonZero, NonZeroG, TotalResids, CurSBBar, CountGroupSq); R_FlushConsole();
       }
       WeGo = 1;
       REAL(SSigma)[tt1] = this->OnSigma;
     }
     if (Rf_length(SSigma) > tt1+2) {
       if (Verbose >= 3) {
         Rprintf("**** UpdateSigma(%d,%d): OnSigma=%f we have tt1 (tt1+2) will put in SSigma. length=%d\n", tt1, tt2, OnSigma, Rf_length(SSigma)); R_FlushConsole();
         Rprintf("**** Putting into tt1=%d SSigma, this->OnSigma = %f, NonZero=%d, NonZeroG=%f, TotalResids=%f, CurSBBar = %f, CountGroupSq = %f\n", tt1, this->OnSigma,
           NonZero, NonZeroG, TotalResids, CurSBBar, CountGroupSq); R_FlushConsole();
       }
       WeGo = 1;
       REAL(SSigma)[tt1+1] = OnSigma;
     }
     
     if (!Rf_isNull(SSigma) && Rf_length(SSigma) == 1) {
       REAL(SSigma)[0] = this->OnSigma;
       WeGo = 1; 
     } 
     if (WeGo == 0) {
       Rprintf("**** UpdateSigma(%d,%d)  Woah, we could not put anything into SSigma!\n"); R_FlushConsole();
     }
   } else {
     if (Verbose >= 3) {
       Rprintf("**** ERROR!!!! ATwoLassoObject2011.cc:::UpdateSigma(%d,%d) problem with RSSigma to account for!\n"); R_FlushConsole();
     }
     if (RSSigma == NULL) {
       Rprintf("****    RSSigma is really NULL\n");
     }
     if (Rf_isNull(SSigma)) {
       Rprintf("****   SSigma is a NULL Object\n");
     }
     if (!Rf_isReal(SSigma)) {
       Rprintf("****  SSigma is not a real Object!\n"); R_FlushConsole();
     }
   } 
	 if (Verbose > 3) {
		 Rprintf((char*) "**** Finished XTX Calculating SigmaNoiseSq, Sigma = %f\n", OnSigma);
     Rprintf("********************************UpdatesigmaSq()******************\n");
		 R_FlushConsole(); R_ProcessEvents();	     
   }     
 } else {
    Rprintf("**** TwoLassoSexp:: Estimate SigmaSq, you've got to be kidding!\n");
    R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
 }
 return(1);	
}


/* ========================================================================== */
/*                                                                            */
/*   GroupTwoLasso functions                                                            */
/*   (c) 2012 Author                                                          */
/*                                                                            */
/*   Description                                                              */
/*    With Random Effects, BBOn Calculations must be maximized differently.                                                                       */
/* ========================================================================== */


int TwoLassoSexp::GroupTwoUpdateBBOn1() {
  if(Verbose > 4) {
    Rprintf("###################################################################\n");
    Rprintf("####  ATwoLassoObject2011.cc::GroupTwoUpdateBBOn1(tt1=%d,tt=%d):  \n", tt1, tt2);
    Rprintf("####          Updating BBOn1\n",
      tt2); R_FlushConsole();
  }
	int ii;
	double OnlambdaDiff = (OnLambdaD - OnLambdaA) / 2.0;
	if (Verbose > 4) {
    Rprintf("####          OnlambdaDiff = %f, OA = %f, OD = %f\n", 
      OnlambdaDiff, OnLambdaA, OnLambdaD);
    R_FlushConsole();
  }
	int SuccFlag = 0;
	if (SBBOn1 == NULL || Rf_isNull(SBBOn1) || Rf_length(SBBOn1) < NumGroups+FirstRandomIndex ||
    !Rf_isReal(SBBOn1)) {
    Rf_error("#### ERROR ERROR ERROR Can't work because SBBOn1 is defective!");
  }
	if (FixKa > 0) {
	  Rprintf("####          GroupTwoLassoSexp::UpdateBBOn1: No, FixKa = %d, does not work!\n",
      FixKa); R_FlushConsole();
    Rf_error("#### ERROR ERROR  GroupTwoLassoSexp:: No way!\n"); 
		SuccFlag = FixKCalculateBBOn();
		if (SuccFlag < 0) {
			Rprintf("#### ERROR ERROR UpdateBBOn1, we got a return FixKCalculateBBOn is angry \n");
			R_FlushConsole(); return(-1);
		}
		return(SuccFlag);
	}
	if (Rf_length(SBBOn1) < NumGroups+FirstRandomIndex) {
    Rprintf("#### ERROR ERROR SBBOn1 only has Rf_length %d, for p = %d", 
      Rf_length(SBBOn1), p); R_FlushConsole();
    Rf_error("This ain't going to work!");
  }
  double OnPiA = this->OnPiA;
  double OnRPiA = this->OnPiA;
  if (this->OnRPiA > 0.0 && this->OnRPiA < 1.0) {
    OnRPiA = this->OnRPiA;
  }
  
  double LogitLambdaDiff = log( (1.0 - OnPiA)) -log(OnPiA) + log(OnLambdaD) -log(OnLambdaA);
  double LambdaDiff = (OnLambdaD - OnLambdaA);
  double ZeroProb = 1.0 / ( 1.0+ exp(LogitLambdaDiff) );
  if (CDO->XTXFlag == 2) {  
    if (FirstRandomIndex > 0) {
    for (ii = 0; ii < FirstRandomIndex; ii++) {
      //if (REAL(SOnBeta)[Onii] == 0.0) {
        REAL(SBBOn1)[ii] = ZeroProb;
      /*} else {
        REAL(SBBOn1)[ii] = LogitLambdaDiff+ LambdaDiff*fabs(REAL(SBBOn1)[ii]);
        if (REAL(SBBOn1)[ii] <= -10) {
          REAL(SBBOn1)[ii] = exp(REAL(SBBOn1)[ii]);
        } else if (REAL(SBBOn1)[ii] >= 10) {
          REAL(SBBOn1)[ii] = 1.0-exp(-1.0*REAL(SBBOn1)[ii]);
        } else {
          REAL(SBBOn1)[ii] = exp(REAL(SBBOn1)[ii])/(1.0+exp(REAL(SBBOn1)[ii]));
        }
      } */
    }
    int Onii = 0;
    if (CDO->OnKappaS > 0) {
    for (ii = 0; ii < CDO->OnKappaS; ii++) {
      Onii = CDO->kXFinder[ii];
      if (Onii < FirstRandomIndex) {
        REAL(SBBOn1)[Onii] = fabs(REAL(SOnBeta)[Onii]) * LambdaDiff;
        if ( REAL(SBBOn1)[Onii] - LogitLambdaDiff > 40 ) {
          REAL(SBBOn1)[Onii] = 1.0;
        } else if ( REAL(SBBOn1)[Onii] - LogitLambdaDiff < -15) {
          REAL(SBBOn1)[Onii] = exp ( REAL(SBBOn1)[Onii] - LogitLambdaDiff );
        } else {
          REAL(SBBOn1)[Onii] =  exp ( REAL(SBBOn1)[Onii] - LogitLambdaDiff );
          REAL(SBBOn1)[Onii] = REAL(SBBOn1)[Onii] / (REAL(SBBOn1)[Onii] + 1.0);
        }
      }
    }
    }
    }
  } else {
    if (FirstRandomIndex > 0) {
    for (ii = 0; ii < FirstRandomIndex; ii++) {
      if (REAL(SOnBeta)[ii] == 0) {
        REAL(SBBOn1)[ii] = ZeroProb;
      } else {
        REAL(SBBOn1)[ii] = fabs(REAL(SOnBeta)[ii]) * LambdaDiff;
        if ( REAL(SBBOn1)[ii] - LogitLambdaDiff > 40 ) {
          REAL(SBBOn1)[ii] = 1.0;
        } else if ( REAL(SBBOn1)[ii] - LogitLambdaDiff < -15) {
          REAL(SBBOn1)[ii] = exp ( REAL(SBBOn1)[ii] - LogitLambdaDiff );
        } else {
          REAL(SBBOn1)[ii] =  exp ( REAL(SBBOn1)[ii] - LogitLambdaDiff );
          REAL(SBBOn1)[ii] = REAL(SBBOn1)[ii] / (REAL(SBBOn1)[ii] + 1.0);
        }
      }
    }
    }
  }
  // int exp(- lambda_i sum(|x_i|)) / lambdai    =  int exp(-a) / a da
  //  = Ei(LambdaA * sum(|x_i|) )  - Ei(LambdaD * sum(|x_i|) ) 
  //
  //  Have lambda_k ~  B_k Expo(1/LambdaA) + (1-B_k) Expo(1/LambdaD)
  //   Or density:  
  //     piA exp(-lambda_k/LambdaA)/LambdaA 
  //     + (1-piA)exp(-lambda_k/LambdaD)/LambdaD
  //
  //  If Beta_j(k) ~  Laplace(lambda_k) or
  //          -lambda_k/2 * exp(-|beta_j|*lambda_k)
  // 
  //  Then Posterior lambda_k| Beta_j propto
  //   -(lambda_k/2)^Jexp(- lambda_k* sum|beta_j|)  * 
  //     [piA exp(-lambda_k/LambdaA)/LambdaA 
  //     + (1-piA)exp(-lambda_k/LambdaD)/LambdaD]
  //
  //  We get then that E[lambda_k] is
  //    (piA C_1 + (1-piA)C_2 )/(Z)
  //
  //  Where C_1, C_2 are the integrating constants and Z is partition integration
  //
  //  int e^(-ax) dx = 1/a
  //  int x^J e^(-ax) dx =  int (ax)^J e^(-ax) adx / a^(J+1) = Gamma(J+1)/ a^(J+1)
  //
  //  So Z = piA (J+1) / ( sum |beta_j| + 1/lambdaA)^(J+1) + (1-piA) (J+1)  / ( sum |beta_j| + 1/lambdaD)^(J+1)
  //
  //  And   C_1 =  (J+2)/ (sum|beta_j| +1/lambdaA)^(J+2)
  //        C_2 =  (J+2)/ (sum|beta_j| +1/lambdaD)^(J+2)
  //
  //
  //  So true we can consider posterior of B:
  //    E[B_k|betas J(k)] = C_A / Z
  //  where C_A =   piA (J+1) / (sum|beta_j| +1/lambdaA)^(J+2)
  //
  double SumBeta = 0.0;    int CLen;
  double C1, C2, Z = 0.0;
  double InvLambdaA = 1.0 / OnLambdaA;
  double InvLambdaD = 1.0 / OnLambdaD;
  int One = 1;
  double lPiA = log(OnRPiA) - log(1-OnRPiA) + log(OnLambdaD)-log(OnLambdaA);
  double KBB = 0.0;
  for (ii = 0; ii < NumGroups; ii++) {
    SumBeta = 0.0;
    if (ii == 0) {
      CLen = EndGroupLocations[0] - FirstRandomIndex +1;
      if (CLen <= 0) {
        Rf_error("Error: We have for ii = 0 that CLen =%d!\n", CLen);
      }
      SumBeta = F77_CALL(dasum)(&CLen, REAL(SOnBeta) + FirstRandomIndex, &One);
    } else {
      CLen = EndGroupLocations[ii] - EndGroupLocations[ii-1];
      if (CLen <= 0) {
        Rf_error("Error: We have for ii = %d that CLen =%d!\n", ii, CLen);
      }
      SumBeta = F77_CALL(dasum)(&CLen, REAL(SOnBeta) + EndGroupLocations[ii-1]+1, &One);
    }
    C1 = - (CLen+1)* log(SumBeta+InvLambdaA)+ lPiA;
    C2 = - (CLen+1)* log(SumBeta+InvLambdaD);
    if (ii >= get_NumGroups()) {
      Rf_error("No, we will not fill GroupLambdaEstimates into ii = %d NumGroups=%d!\n",
        ii, get_NumGroups());
    }
    if (Rf_length(GroupsBBOn1->asSexp()) <= ii) {
      Rf_error("Error:UpdateGroupBBOn1: no, CLen = %f \n");
    }
    if (C1 <= -10 && C2 >= 10) {
      REAL(GroupsBBOn1->asSexp())[ii] = exp(C1-C2);  
      KBB = C2; Z = 1.0;
      GroupLambdaEstimates[ii] = CLen/(SumBeta+InvLambdaD);
    } else if (C2 <= -10 && C1 >= 10) {
      REAL(GroupsBBOn1->asSexp())[ii] = 1.0-exp(C2-C1);  
      KBB = C2; Z = 1.0;
      GroupLambdaEstimates[ii] = CLen/(SumBeta+InvLambdaA);  
    } else if (C2 >= C1) {
      KBB = C2;  Z = exp(C1-C2) + 1.0;
      REAL(GroupsBBOn1->asSexp())[ii] =  exp(C1-C2)/(exp(C1-C2)+1.0);     
      C1 = - (CLen+2)* log(SumBeta+InvLambdaA)+ lPiA - C2;
      C2 = - (CLen+2)* log(SumBeta+InvLambdaD) - C2; 
      GroupLambdaEstimates[ii] = CLen*(exp(C1)+exp(C2))/Z;  
    } else if (C2 >= C1) {
      KBB = C1;  Z = 1.0 + exp(C2-C1);
      REAL(GroupsBBOn1->asSexp())[ii] = exp(C1-C2)/(exp(C1-C2)+1.0);   
      C2 = - (CLen+2)* log(SumBeta+InvLambdaD) - C1; 
      C1 = - (CLen+2)* log(SumBeta+InvLambdaA)+ lPiA - C1;
      GroupLambdaEstimates[ii] = CLen*(exp(C1)+exp(C2))/Z;   
    }
    if (Verbose >= 4) {
      Rprintf("#### GroupBBOn1 ii=%d/%d: we get CLen=%d and SumBeta = %f.  Z=%f, C1=%f, C2=%f, KBB=%f, lPiA=%f\n",
        ii, NumGroups, CLen, SumBeta, Z, C1, C2, KBB, lPiA);R_FlushConsole();
    }
    if (FALSE) {
    if (SumBeta > 0.0) {
       // OnLambdaA < OnLambdaD
       // 1/OnLambdaA > 1/OnLambdaD
       // 
       C1 = log(SumBeta + InvLambdaA) - log(SumBeta + InvLambdaD);
       if (C1 * (CLen+1) + lPiA < -10) {
         Rprintf("Group Lasso, very negative!  How did we did we get here?");
          GroupLambdaEstimates[ii] = 1 / (SumBeta + InvLambdaD);   
          REAL(GroupsBBOn1->asSexp())[ii] = exp(log(SumBeta + InvLambdaA) * (CLen+1)  +lPiA);     
       } else if (C1 * (CLen+1) + lPiA > 10) {
          GroupLambdaEstimates[ii] = 1 / (SumBeta + InvLambdaA);
          REAL(GroupsBBOn1->asSexp())[ii] = 1.0;
       } else {
          GroupLambdaEstimates[ii] = 
           (CLen+1) * (exp(-C1*(CLen+1.0) - log(SumBeta+InvLambdaA) + lPiA) +
            exp( -log(SumBeta+InvLambdaD) )) /
            (exp(-C1*(CLen+1.0) + lPiA) +1.0);
          REAL(GroupsBBOn1->asSexp())[ii] = exp( -C1 * (CLen+1.0) + lPiA -
            log(exp(-C1 * (CLen+1.0) + lPiA)+ 1.0));
       } 
    } else {
       C1 = log(OnLambdaA) - log(OnLambdaD);
       if (C1 * (CLen+1) + lPiA < -10) {
         //Rprintf("Group Lasso, very negative!  How did we did we get here?\n");
          GroupLambdaEstimates[ii] = OnLambdaD;
          REAL(GroupsBBOn1->asSexp())[ii] = exp(CLen*C1 + log(OnPiA) - log(1.0-OnPiA));    
       } else if (C1 * (CLen+1) + lPiA> 10) {
         //Rprintf("Group Lasso, very positive!  Why do you like OnLambdaA!?\n");
          GroupLambdaEstimates[ii] = OnLambdaD;
          REAL(GroupsBBOn1->asSexp())[ii] = 1.0;
       } else {
          GroupLambdaEstimates[ii] = (CLen+1) * (
            exp((CLen+1.0) * C1) * OnPiA/(1.0-OnPiA) +
             1.0) /
           (exp(C1 * (CLen+1.0) + lPiA ) + 1.0);
          REAL(GroupsBBOn1->asSexp())[ii] = exp( C1 * CLen + log(OnPiA) - log(1.0-OnPiA) -
             log(exp(C1 * CLen+log(OnPiA)- log(1.0-OnPiA))+1.0));
       } 
    }
    }
  }
  if (Verbose >= 4) {
    Rprintf("####      About to Update BBOn1 \n"); R_FlushConsole();
  }
  if (Rf_length(SBBOn1) == p) {
    int St = FirstRandomIndex;  int CNN = St;
    for (int ii = 0; ii < NumGroups; ii++) {
       for (CNN = St; CNN < EndGroupLocations[ii]+1; CNN++) {
         if (CNN >= p || CNN < 0 || CNN >= EndGroupLocations[ii]+1) {
           Rf_error("UpdateGroupBBOn1, no, CNN=%d, p=%d\n", CNN, p); 
         }
         REAL(SBBOn1)[CNN] = REAL(GroupsBBOn1->asSexp())[ii];
         if (CNN < 0 || CNN >= p) {
           Rf_error("UpdateGroupBBOn1, no, CNN=%d, p=%d\n", CNN,p);
         }
         CDO->OnGammas[CNN] = GroupLambdaEstimates[ii];
       }
       St = EndGroupLocations[ii]+1;
    }
  }
  if (tt1 <= Rf_length(SLambdaDK)-2) {
    int ACopyLength = Rf_length(GroupsBBOn1->asSexp());
    int AOne = 1;
    if (BackGroupsBBOn1 != NULL && !Rf_isNull(BackGroupsBBOn1->asSexp()) &&
      Rf_length(BackGroupsBBOn1->asSexp()) >= ACopyLength) {
      F77_CALL(dcopy)(&ACopyLength, REAL(GroupsBBOn1->asSexp()), &AOne,
        REAL(BackGroupsBBOn1->asSexp()), &AOne);
    }
  }
  if (tt1 <= Rf_length(SLambdaDK)-1 && RsGroupLambdaRecord != NULL &&
    !(Rf_length(RsGroupLambdaRecord->asSexp()) <= 0)) {
    int AOne2 = 1;
    if (tt1 >= 0 && NumGroups*tt1+NumGroups <= Rf_length(RsGroupLambdaRecord->asSexp())) {
      F77_CALL(dcopy)(&NumGroups, GroupLambdaEstimates, &AOne2,
      REAL(RsGroupLambdaRecord->asSexp()) + tt1 * NumGroups, &AOne2);
    }
  }
  if (NoShrinkColumns != NULL && !Rf_isNull(NoShrinkColumns->asSexp()) && 
    Rf_length(NoShrinkColumns->asSexp()) >= 1) {
    if (Verbose >= 4) {
      Rprintf("#### Group Update BBOn1 updating NoShrinkColumns\n"); R_FlushConsole();
    }
    for (int iti = 0; iti < Rf_length(NoShrinkColumns->asSexp()); iti++) {
      if (Rf_length(SBBOn1) > INTEGER(NoShrinkColumns->asSexp())[iti]) {
        REAL(SBBOn1)[INTEGER(NoShrinkColumns->asSexp())[iti]] = 1.0;
      }
    }
  }
  if (Verbose > 4) {
      Rprintf("#### ATwoLassoObject2011.cc::UpdateGroupBBOn1() Returning!\n",
        OnLambdaA, OnLambdaD); R_FlushConsole();
      Rprintf("###############################################################\n"); R_FlushConsole();
  }
	return(1);
}



int TwoLassoSexp::GroupTwoFixKCalculateBBOn() {
  Rprintf("########################################################################\n"); R_FlushConsole();
  Rprintf("#### ATwoLassoObject2011.cc: GroupTwoFixKCalculateBBOn1() GroupTwoLassoSexp::");
  Rprintf("####    Hey, no, we will not do this, this is confusing for group variables!\n");
  R_FlushConsole(); return(-999);
  if (FixKa > p || FixKa < 0) {
	  Rprintf("FixKCalculateBBOn:  Cannot make this specious Kcount = %d, p=%d \n", 
      FixKa, p);
  }
  if (Verbose > 3) {
    Rprintf(" We're Running FixKCalculateBBOn with FixKa = %d, p=%d\n", 
      FixKa, p); R_FlushConsole();
  }
  int ii;
	int SuccFlag = 0;
  if (SuccFlag < 0) {
    Rprintf("GroupTwoFixKCalculate, shouldn't have happened.\n");
  }
	double Solvedmu = 0;    
  double CurrentCount = 0; 
	double LogitOnPiA;
	if (OnPiA > 0.0 && OnPiA < 1.0) {
	  LogitOnPiA = log(OnPiA / (1-OnPiA));
  } else if (OnPiA == 1.0) {
	   LogitOnPiA = 20;
  } else if (OnPiA == 0.0) {
	  LogitOnPiA = -20;
  } else {
	  LogitOnPiA = 0;
  }
	double OnCurrentPiACount = AofX(LogitOnPiA);
	double ApparentMax = 0;
	double ApparentMin  = 0;
	double CountAtMax = 0;
	double CountAtMin = 0;
	int BreakItYouBoughtIt = 0;
 
	if (Verbose > 4) {
		Rprintf(" \n \n#### FixKCalculateBBOn(): tt1 = %d, tt2 = %d, ", tt1, tt2);
    Rprintf(" OnlambdaA = %.4e, OnlambdaD = %.4e \n",
			(double) OnLambdaA, (double) OnLambdaD);
		Rprintf("#### OnGammas  = "); PrintVector(CDO->OnGammas, p);
		Rprintf("\n#### OnBetas = "); PrintVector(REAL(SOnBeta), p);
		Rprintf("\n####  Now to Run Calculations \n");
		R_FlushConsole(); R_ProcessEvents();
  }
	if (OnLambdaA == OnLambdaD) {
		OnPiA = ((double)FixKa / ((double)p));
		if (OnPiA <= 0 || OnPiA >= 1) {
			Rprintf("#### 2L:FixKCalculateBBOn(), we have a ppiuse = %.4e Error, ", OnPiA);
      Rprintf("#### 2L:quitting \n", OnPiA);
			R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return(-1);
    } 
		LogitOnPiA = log(OnPiA / (1-OnPiA));
		Solvedmu = LogitOnPiA;
		CurrentCount = AofX(LogitOnPiA);
		if (Verbose > 4) {
			Rprintf("#### 2L: FixKCalculateBBOn(): tt1 = %d, tt2 = %d, OnlambdaA = %.4e, ",
        tt1, tt2, (double) OnLambdaA);
      Rprintf(" OnLambdaD = %.4e \n", (double) OnLambdaD);
			Rprintf("#### 2L: Started at Default, ppiuse = %.4e, LogitOnppiuse = %.4e, ",
			  (double) OnPiA, (double) LogitOnPiA);
      Rprintf("#### 2L: CurrentCount = %.4e \n", (double) CurrentCount);
		  R_FlushConsole();
    }
  } else {       
	  if (Verbose > 4) {
			Rprintf("#### 2L:FixKCalculateBBOn(): Starting to Solve: ppiuse = %.4e,", OnPiA);
      Rprintf("#### 2L: OnCurrentPiACount = %.4e \n",  OnCurrentPiACount);
			R_FlushConsole(); R_ProcessEvents();
		}
		Solvedmu = SolveFormu(FixKa, SUFFEPSILON_FD, MAXITERS_FD, LogitOnPiA);   
	    CurrentCount = AofX(Solvedmu);
		if (Verbose > 4) {
		  Rprintf("##############################################\n");
		  R_FlushConsole();
			Rprintf("#### 2L:FixKCalculateBBOn(): Initial Solvedmu = %.4e, ", Solvedmu);
      Rprintf("#### CurrentCount = %.4e \n", (double) CurrentCount);
		  R_FlushConsole(); R_ProcessEvents();
		}	        
		if (ISNAN(Solvedmu) ||  !R_FINITE(Solvedmu) ||
		  (fabs( ((double) FixKa )- CurrentCount ) > SUFFEPSILON_FD * 10 )  ) {			            
	    ApparentMax = 60;
	    ApparentMin  = -60;
	    CountAtMax = AofX(ApparentMax);
	    CountAtMin = AofX(ApparentMin);
	    if (OnCurrentPiACount < CountAtMax && OnCurrentPiACount > (double) FixKa) {
		    ApparentMax = LogitOnPiA;
		    CountAtMax = OnCurrentPiACount;
      }  
	    if (OnCurrentPiACount > CountAtMin && OnCurrentPiACount < (double) FixKa) {
		    ApparentMax = LogitOnPiA;
		    CountAtMax = OnCurrentPiACount;
      } 
      if (CountAtMax < FixKa) {
        if (Verbose > 1) {
          Rprintf("####   CountAtMax=%d < FixKa = %d.\n", CountAtMax, FixKa); R_FlushConsole();
	        Rprintf("#### 2L: FindBBOn: Apparently FixKa is already less than Max.\n");
          Rprintf("####   Probabilities to big to tweak.\n");  R_FlushConsole();
	        R_FlushConsole(); R_FlushConsole();
	      }
	      LogitOnPiA = ApparentMax; OnCurrentPiACount = CountAtMax;
      } else if (CountAtMin > FixKa) {
        if (Verbose > 1) {
          Rprintf("#### 2L: FindBBOn: CountAtMin=%d > FixKa = %d.\n", CountAtMin, FixKa); R_FlushConsole();
	        Rprintf("#### 2L:FindBBOn: Apparently FixKa is already less than Max.\n");
          Rprintf("####  Probabilities to big to tweak.\n");
	        R_FlushConsole(); R_FlushConsole();
	      }
	      LogitOnPiA = ApparentMin; OnCurrentPiACount = CountAtMin;	      
      } else {	
	      BreakItYouBoughtIt = 0;	
	      OnCurrentPiACount = AofX(LogitOnPiA);            
			  for (ii = 0; ii < MAXITERS_FD; ii++) {			       
				  if ( fabs( ((double) FixKa )- OnCurrentPiACount ) < 
            SUFFEPSILON_FD * 10) {
					  Solvedmu = LogitOnPiA; 
            CurrentCount = OnCurrentPiACount;
					  BreakItYouBoughtIt = 1;
					  break;    
				  }
				  if ( ISNAN(Solvedmu) || !R_FINITE(Solvedmu) ||
				    fabs(((double) FixKa )- CurrentCount ) >
				    fabs(((double) FixKa) - OnCurrentPiACount) ) {
					  if (FixKa > OnCurrentPiACount) {
						  LogitOnPiA = (ApparentMax + LogitOnPiA) / 2.0 ;
						  OnCurrentPiACount = AofX(LogitOnPiA);
						  if (OnCurrentPiACount < FixKa) {
							  ApparentMin = LogitOnPiA;
							  CountAtMin = OnCurrentPiACount;
						  } else if (OnCurrentPiACount > FixKa) {
							  ApparentMax = LogitOnPiA;
							  CountAtMax = OnCurrentPiACount;							       
						  }
						  Solvedmu = LogitOnPiA; 			       
					  } else {
						  LogitOnPiA = (ApparentMin + LogitOnPiA) / 2.0 ;
						  OnCurrentPiACount = AofX(LogitOnPiA);
						  if (OnCurrentPiACount < FixKa) {
							  ApparentMin = LogitOnPiA;
							  CountAtMin = OnCurrentPiACount;
						  } else if (OnCurrentPiACount > FixKa) {
							  ApparentMax = LogitOnPiA;
							  CountAtMax = OnCurrentPiACount;							       
						  }						       
						  Solvedmu = LogitOnPiA	;
					  }         
					  Solvedmu = SolveFormu(FixKa,SUFFEPSILON_FD, 
              MAXITERS_FD, LogitOnPiA);   
				  } else {
            if (CurrentCount < FixKa) {
							ApparentMin = Solvedmu;
							CountAtMin = CurrentCount;
						} else if (CurrentCount > FixKa) {
							ApparentMax = Solvedmu;
							CountAtMax = CurrentCount;							       
						}
					  LogitOnPiA = Solvedmu; OnCurrentPiACount = CurrentCount;
					  Solvedmu = SolveFormu(FixKa,SUFFEPSILON_FD, 
              MAXITERS_FD, LogitOnPiA);   
			    }
			    CurrentCount = AofX(Solvedmu);
			    if (!ISNAN(Solvedmu) && R_FINITE(Solvedmu) && 
			      fabs( ((double) FixKa )- CurrentCount ) < SUFFEPSILON_FD * 10) {
				    LogitOnPiA = Solvedmu;
				    OnCurrentPiACount = CurrentCount;
				    BreakItYouBoughtIt = 2;
				    break;
			    }
			    if (LogitOnPiA > 50) {
				    Solvedmu = LogitOnPiA;
				    CurrentCount = OnCurrentPiACount;
				    BreakItYouBoughtIt = 3;
				    break;   
			    } else if (LogitOnPiA < -50) {
				    Solvedmu = LogitOnPiA;
				    CurrentCount = OnCurrentPiACount;
				    BreakItYouBoughtIt = 4;
				    break;   				           			           
			    }
		    }
		    if (ISNAN(LogitOnPiA) || !R_FINITE(LogitOnPiA) )  {
			    Rprintf("#### 2L: FixKCalculateBBOn ISNAN(LogitOnPiA) failed, \n");
          Rprintf("####    we've got FixKa = %d, p = %d, ", FixKa, p);
          Rprintf("####    LogitOnPiA = %.4e, piLogitOnPiA = %.4e \n",
			      LogitOnPiA, exp(LogitOnPiA) / ( 1 + exp(LogitOnPiA)) );
			    Rprintf("####    ii = %d \n", ii); R_FlushConsole();
			    Rprintf("####    BreakItYouBoughtIt = %d\n", BreakItYouBoughtIt);			                  
			    R_FlushConsole(); R_ProcessEvents();
			    SuccessFlag = -1; return(-1);
		    }		           
		    if (fabs( ((double) FixKa )- OnCurrentPiACount ) > SUFFEPSILON_FD * 4000
		      && LogitOnPiA < 50 && LogitOnPiA > -50) {
			    Rprintf("#### 2L: FixKCalculateBBOn we failed, we've got FixKa = %d, ", FixKa);
          Rprintf("####    p = %d, LogitOnPiA = %.4e, piLogitOnPiA = %.4e \n",
			      p, (double) LogitOnPiA,(double) exp(LogitOnPiA) /
            ( 1 + exp(LogitOnPiA)) ); R_FlushConsole();
			    Rprintf("#### 2L: OnCurrentPiACount = %.4e, CurrentCount = %.4e \n",
			      OnCurrentPiACount, CurrentCount);
			    Rprintf("#### 2L: Solvedmu = %.4e, piSolvedmu = %.4e \n", 
			      Solvedmu, exp(Solvedmu) / ( 1.0 + exp(Solvedmu)));
			    Rprintf("#### 2L: ApparentMax = %.4e, CountAtMax = %.4e \n", 
			      ApparentMax, CountAtMax);	R_FlushConsole();
			    Rprintf("#### 2L: ApparentMin = %.4e, CountAtMin = %.4e \n", 
			      ApparentMin, CountAtMin);	R_FlushConsole();		               		               
			    Rprintf("#### 2L: ii = %d \n", ii);  R_FlushConsole();
			    Rprintf("#### 2L: BreakItYouBoughtIt = %d\n", BreakItYouBoughtIt);
			    R_FlushConsole(); R_ProcessEvents();
			    SuccessFlag = -1; return(-1);
		    }
	    }	           
		} else {
			LogitOnPiA = Solvedmu;
			OnCurrentPiACount = CurrentCount;
		}
		if (ISNAN(LogitOnPiA) || !R_FINITE(LogitOnPiA)) {
		  Rprintf("###################################################### \n");
		  Rprintf("#### 2L  Error on Logit \n"); R_FlushConsole();
			Rprintf("#### 2L: FixKCalculateBBOn: LogitOnPiA ISNAN!!, ");
			R_FlushConsole();
      Rprintf("#### 2L: LogitOnPiA = %.4e, FixKa = %d \n",  LogitOnPiA, FixKa);
			SuccessFlag = -1; return(-1);
		}		       
		if (LogitOnPiA > 50) {
			OnPiA = 1.0;
		} else if (LogitOnPiA < -50) {
			OnPiA= exp(LogitOnPiA);
		} else { 
			OnPiA = exp(LogitOnPiA) / ( 1.0 + exp(LogitOnPiA));
		}
  }
	if (Verbose > 3) {
	  Rprintf("#### 2L: tt1=%d, tt2=%d, FixKCalculateBBOn settles upon ", tt1,tt2);
    Rprintf("#### 2L: Solvedmu = %.4e, LogitOnppiuse = %.4e, ppiuse = %.4e, ",
      Solvedmu, LogitOnPiA, OnPiA);   R_FlushConsole();
    Rprintf("#### 2L: CurrentCount = %.4e\n", CurrentCount);
	  R_FlushConsole(); R_ProcessEvents();
  }
  if (Rf_length(SPiA) >= tt1+1)  {
    REAL(SPiA)[tt1] = OnPiA;
  }
  if (Rf_length(SPiA) >= tt1+2) {
    REAL(SPiA)[tt1+1] = OnPiA;
  }
  if (Verbose > 3) {
    Rprintf("#### 2L: FixBBOn: We reached near the end, now Inser %f\n",
      LogitOnPiA); R_FlushConsole();
  }
	SuccFlag = InsertBBOn1ofX(LogitOnPiA);                   
  if (Verbose > 3) {
    Rprintf("#### 2L: FixBBOn: Really at the end\n"); R_FlushConsole();
  }
	return(1);
}

	  
///////////////////////////////////////////////////////////////////////////////
//   TwoLassoSexp::GroupSetupCDO(int SetBeta)
//
//   Setup CoordinateDescentObject for a run (Grouped Version)
//
//    if (SetBeta == 0) 
//     then the residuals are already okay in the CDO (XTYResid specifically)
//    Else, we need to recalulate residuals according to moved Betas
int TwoLassoSexp::GroupTwoSetupCDO(int SetBeta) {
  if(Verbose >= 4) {
    Rprintf("----------------------------------------------------------------------\n");
    Rprintf("----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d):  \n", tt1, tt2);
    Rprintf("----      tt2 = %d, Running SetupCDO\n",
      tt2); R_FlushConsole();
  }
  if (tt1 == 0 && Verbose >= 1) {
    if (SetBeta == 0) {
      Rprintf("----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d) Start, SetBeta = %d, so we are keeping current Beta (default)\n",
        tt1, tt2, SetBeta); R_FlushConsole();
    } else {
      Rprintf("----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d) Start, SetBeta = %d, so we are refilling Beta\n",
        tt1, tt2, SetBeta); R_FlushConsole();
    }
  }
  if (CDO == NULL) {
    Rprintf((char*)"----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d), CDO is NULL\n", tt1, tt2);
    R_FlushConsole(); R_ProcessEvents(); 
    SuccessFlag = -1; error("No CDO"); }
  if (Verbose > 0 && SetBeta > 0) {
	  Rprintf((char*) "----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d), Verbose = %d, SetBeta = %d so we are refilling Beta in CDO\n",  tt1, 
      tt2, Verbose, SetBeta);
	  R_FlushConsole();
  }
  if (NumGroups <= 0) {
    Rf_error("---- ATwoLassoObject2011.cc:GroupTwoSetupCDO(tt1=%d, tt2=%d), why, NumGroups = %d?\n",
      tt1, tt2, NumGroups);
  }
  /*
  if (EndGroupLocations[NumGroups-1] != CDO->kLen-1) {
    Rprintf("TwoLassoSexp::GroupTwoSetupCDO: big Error, kLen=%d, NumGroups = %d, EndGroups[%d-1]=%d\n",
      CDO->kLen, NumGroups, NumGroups, EndGroupLocations[NumGroups-1]);
    int Acn = 0;
    Rprintf("EndGroupLoctions =\n( ");
    for (int jj1 = 0; jj1 < NumGroups; jj1++) {
      if (Acn == 10) {Rprintf("\n  ");  R_FlushConsole(); Acn = 0; }
      if (jj1 == NumGroups-1) {
        Rprintf("%d)\n", EndGroupLocations[jj1]); R_FlushConsole();
      } else {
        Rprintf("%d, ", EndGroupLocations[jj1]);
      }
      Acn++;
    }
    Rf_error("TwoLassoSexp: Can't go beyond this!\n");
  }        */
  int ii;
  int AStart, jj;
	AStart = FirstRandomIndex;
  if (InverseGammaConstant != 1.0 && InverseGammaConstant > 0.0) {
    if (tt1 == 0 && Verbose >= 2) {
      Rprintf("---- ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d) InverseGammaConstant is %f, nonzero game.\n", 
        tt1, tt2, InverseGammaConstant); 
        R_FlushConsole();
    }
    if (FirstRandomIndex > 0) {
	    for (ii = 0; ii < FirstRandomIndex; ii++) {
	      CDO->OnGammas[ii] = (REAL(SBBOn1)[ii] * OnLambdaA + 
	        (1.0- REAL(SBBOn1)[ii]) *OnLambdaD)  * OnSigma * 
	        InverseGammaConstant;
	    }
	  }

    if (AStart < 0) { AStart = 0; }
	  for  (ii = 0; ii < NumGroups; ii++) {
	    for (jj = 0; jj < EndGroupLocations[ii] - AStart+1; jj++) {
 	      CDO->OnGammas[AStart + jj] = GroupLambdaEstimates[ii]  * OnSigma * 
	        InverseGammaConstant;  
      }
      AStart = EndGroupLocations[ii]+1;  
	  }
  } else {
    if (tt1 == 0 && Verbose >= 2) {
      Rprintf("----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d): ", tt1, tt2);
      Rprintf(" InverseGammaConstant is %f, or zero game!, FirstRandom = %d\n", InverseGammaConstant, FirstRandomIndex); R_FlushConsole();
    }
    if (FirstRandomIndex > 0) {
    for (ii = 0; ii < FirstRandomIndex; ii++) {
	    CDO->OnGammas[ii] = (REAL(SBBOn1)[ii] * OnLambdaA + 
	      (1.0- REAL(SBBOn1)[ii]) *OnLambdaD) * OnSigma;
	  }
	  }
	  AStart = FirstRandomIndex;
    if (AStart < 0 || AStart >= CDO->kLen) { AStart = 0; }
	  for  (ii = 0; ii < NumGroups; ii++) {
	    for (jj = 0; jj < EndGroupLocations[ii] - AStart+1; jj++) {
        if (AStart+jj > EndGroupLocations[ii] || AStart+jj >= CDO->kLen) {
          Rprintf("--- Error in GroupTwoSetupCDO(%d,%d), AStart=%d, jj=%d, CDO->kLen=%d, ii=%d, EndGroup=%d\n",
            tt1, tt2, AStart, jj, CDO->kLen, ii, EndGroupLocations[ii]); R_FlushConsole();
          Rf_error("Error this is a quit!\n");
        }
 	      CDO->OnGammas[AStart + jj] = GroupLambdaEstimates[ii]  * OnSigma * 
	        InverseGammaConstant;              
      }
      AStart = EndGroupLocations[ii]+1;  
    }
    if (tt1 == 0 && Verbose >= 2) {
      Rprintf("----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d) ", tt1, tt2);
      Rprintf(" Well We've setup the whole GroupTwoLassoCpp, AStart=%d. ii = %d\n", AStart, ii); R_FlushConsole();
    }
  }
  if (SetBeta == 1 && CDO->DOnBeta == 0) {
    int One = 1;
    if (Verbose >= -10) {
      Rprintf("----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d) ", tt1, tt2);
      Rprintf("  SetBeta=1, DOnBeta = 0 now copying.\n"); R_FlushConsole();
    }
    F77_CALL(dcopy)(&p, REAL(SOnBeta), &One, CDO->OnBeta, &One);  
    if (Verbose >= -10) {
      Rprintf("----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d) ", tt1, tt2);
      Rprintf("  Now Add AllNewNonzeroCoords.() \n"); R_FlushConsole();
    }
    CDO->PrintFlag = 4;
    CDO->AddAllNewNonzeroCoords();
    CDO->PrintFlag = Verbose -3;
    CDO->ModifiedXTResid = 1;
  }
  if (Verbose > 5) {
    Rprintf("----  ATwoLassoObject2011.cc::GroupTwoSetupCDO(tt1=%d,tt2=%d)  ", tt1, tt2);
	  Rprintf(" Let's show CDO StatusReport(). \n"); R_FlushConsole(); 
	  CDO->StatusReport();
  }

  return(1);
}

///////////////////////////////////////////////////////////////////////////////
//   TwoLassoSexp::SetupCDO(int SetBeta)
//
//   Setup CoordinateDescentObject for a run
//
//    if (SetBeta == 0) 
//     then the residuals are already okay in the CDO (XTYResid specifically)
//    Else, we need to recalulate residuals according to moved Betas
int TwoLassoSexp::SetupCDO(int SetBeta) {
  if(Verbose > 4) {
    Rprintf("----  ATwoLassoObject2014.cc::SetupCDO(tt1=%d,tt2=%d)  ", tt1, tt2);
    Rprintf("  tt2 = %d, Running SetupCDO\n",
      tt2); R_FlushConsole();
  }
  if (NumGroups >= 0) {
    if (Verbose >= 2) {
      Rprintf("----  ATwoLassoObject2014.cc::SetupCDO(tt1=%d,tt2=%d) ", tt1, tt2);
      Rprintf(" SetupCDO: We will send to a setup version: GroupTwoSetupCDO().\n"); R_FlushConsole();
    }
    int AMark =  GroupTwoSetupCDO(SetBeta);
    if (Verbose >= 2) {
      Rprintf("----  ATwoLassoObject2014.cc::SetupCDO(tt1=%d,tt2=%d) ", tt1, tt2);
      Rprintf("SetupCDO: We have completed GroupTwoSetupCDO() and AMark = %d\n", AMark);
      R_FlushConsole();
    }
    return(AMark);
  }
  if (CDO == NULL) {
    Rprintf("----  ATwoLassoObject2014.cc::SetupCDO(tt1=%d,tt2=%d) ", tt1, tt2);
    Rprintf((char*)"   SetupCDO, CDO is NULL\n");
    R_FlushConsole(); R_ProcessEvents(); 
    SuccessFlag = -1; Rf_error("No CDO"); }
  if (Verbose > 0 && SetBeta > 0) {
    Rprintf("----  ATwoLassoObject2014.cc::SetupCDO(tt1=%d,tt2=%d) ", tt1, tt2);
	  Rprintf((char*) "TwoLassoSexp::SetupCDO, Verbose = %d, SetBeta = %d \n",  Verbose, SetBeta);
	  R_FlushConsole();
  }
  int ii;

  if (InverseGammaConstant != 1.0 && InverseGammaConstant > 0.0) {
	  for (ii = 0; ii < p; ii++) {
	    CDO->OnGammas[ii] = (REAL(SBBOn1)[ii] * OnLambdaA + 
	      (1.0- REAL(SBBOn1)[ii]) *OnLambdaD)  * OnSigma * 
	      InverseGammaConstant;
	  }

  } else {
    for (ii = 0; ii < p; ii++) {
	    CDO->OnGammas[ii] = (REAL(SBBOn1)[ii] * OnLambdaA + 
	      (1.0- REAL(SBBOn1)[ii]) *OnLambdaD) * OnSigma;
	  }
  }
  if (SetBeta == 1 && GLMCDO== NULL && CDO->DOnBeta == 0) {
    int One = 1;
    F77_CALL(dcopy)(&p, REAL(SOnBeta), &One, CDO->OnBeta, &One);  
    CDO->ModifiedXTResid = 1;
  }
  if (Verbose > 5) {
    Rprintf("----  ATwoLassoObject2014.cc::SetupCDO(tt1=%d,tt2=%d) ", tt1, tt2);
	  Rprintf(" SetupBeta: Let's show CDO StatusUpdate() \n"); R_FlushConsole(); 
	  CDO->StatusReport();
  }

  return(1);
}

///////////////////////////////////////////////////////////////////////////////
// AfterCDO, 
//   Copy things that might need to be copied, or print things to print
//
//
int TwoLassoSexp::AfterCDO() {
  if (CDO == NULL) {Rprintf((char*)"TwoLassoSexp::AfterCDO, CDO is NULL\n");
    R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1); }
  //int ii;
  if (Verbose > 2) {
	  Rprintf((char*)"........................................................................\n");
    Rprintf("... ATwoLassoObject2014.cc():: AfterCDO(tt1=%d,tt2=%d) After CDO \n", tt1, tt2);  R_FlushConsole();
  } 
  int One = 1;
  if (CDO->OnBeta == NULL) {
    Rf_error("... ATwoLassoObject2014.cc():: No Way (tt1=%d,tt=%d, CDO->OnBeta is NULL!\n", tt1, tt2); R_FlushConsole();
  }
  if (CDO->kLen != p) {
    Rf_error("... ATwoLassoObject2014.cc():: No  Way (tt1=%d,tt=%d, CDO->OnBeta is NULL! p = %d, kLen=%d\n", tt1, tt2, p, CDO->kLen); R_FlushConsole();
  }
  if (GLMCDO != NULL) {
    // Beta is all pointers in GLMCDO.  No copies.
    //F77_CALL(dcopy)(&p, GLMCDO->OnBeta, &One, REAL(SOnBeta), &One);     
  } else if (CDO != NULL && CDO->OnBeta != NULL && !Rf_isNull(SOnBeta) && Rf_isReal(SOnBeta) &&
    Rf_length(SOnBeta) == p && CDO->DOnBeta == 0) {    
    //if (p >= LARGEDATA) {
    //  Rprintf("AfterCDO: copying OnBeta into Our Object \n"); R_FlushConsole();
    //}
    F77_CALL(dcopy)(&p, CDO->OnBeta, &One, REAL(SOnBeta), &One);         
  }
  if (Verbose > 2) {
	  Rprintf((char*) ".... AfterCDO(tt1=%d,tt2-=%d) The CDO OnGammas Are : \n",tt1,tt2); R_FlushConsole();
	  PrintVector(CDO->OnGammas, p);
    Rprintf((char*) ".... AfterCDO(tt1=%d,tt2=%d), final OnBeta Are : \n", tt1, tt2); R_FlushConsole();
    PrintVector(REAL(SOnBeta), p);
  }   
  return(1);
}
 
    
///////////////////////////////////////////////////////////////////
//  UpdateMaxPi()
//    For those who need to find maximum pi_A, this asses the result
//   using the current BBOn1 estimates.
//
//
//
int TwoLassoSexp::UpdateHatPiA() {
  if (m1 <= 0) {
    return(-1);
  }
  double SupBHat = 0;
  int ii;
  if (NumGroups <= 0 || FirstRandomIndex >= p || OnRPiA <= 0.0 || OnRPiA >= 1.0){
  if (m1 <= -65) {
    return(2);
  }
  for (ii = 0; ii < p; ii++) {
	    SupBHat += REAL(SBBOn1)[ii];
  }
  OnPiA = (m1 +1.0+ SupBHat) / ( m1 + m2 + p-2.0);
  if (OnPiA < 0) { 
	  OnPiA = (m1) / ( m1 + m2 + p);
  } else if (OnPiA > 1) {
	  OnPiA = (m1 + p) / ( m1 + m2 + p);
  }
  if (Rf_length(SPiA) >= tt1) {
    REAL(SPiA)[tt1] = OnPiA;
  }
  return(1);
  }

  if (FirstRandomIndex > 0) {
  for (ii = 0; ii < FirstRandomIndex; ii++) {
	    SupBHat += REAL(SBBOn1)[ii];
  }
  if (m1 >= -2) {
    OnPiA = (m1 - 1.0+ SupBHat) / ( m1 + m2 + FirstRandomIndex-2.0);
  if (OnPiA < 0.0) { 
	  OnPiA = (m1) / ( m1 + m2 + p);
  } else if (OnPiA > 1.0) {
	  OnPiA = (m1 + p) / ( m1 + m2 + p);
  }
  }
  if (Rf_length(SPiA) >= tt1) {
    REAL(SPiA)[tt1] = OnPiA;
  }
  }
  if (GroupsBBOn1 == NULL || Rf_isNull(GroupsBBOn1->asSexp()) || Rf_length(GroupsBBOn1->asSexp()) <= 0) {
    Rf_error("UpdateHatBBOn1: How do we update PiA if GroupsBBOn1 is bad!\n");
  }
  SupBHat = 0.0;
  for (ii = 0; ii < NumGroups; ii++) {
    SupBHat += REAL(GroupsBBOn1->asSexp())[ii];
  }
  if (m3 >= 0.0 && m4 >= 0.0) {
    OnRPiA = (m3-1.0 + SupBHat) / ( m3 + m4 + NumGroups - 2.0);
    if (OnRPiA < 0.0) {
      OnRPiA = m3 / (m3 + m4 + NumGroups);
    } else if (OnRPiA > 1.0) {
      OnRPiA = (m3 + NumGroups) / (m3 + m4 + NumGroups);
    }
  } else {
    if (m1 >= 0) {
    OnRPiA = (m1-1.0 + SupBHat) / ( m1 + m2 + NumGroups - 2.0);
    if (OnRPiA < 0.0) {
      OnRPiA = m1 / (m1 + m2 + NumGroups);
    } else if (OnRPiA > 1.0) {
      OnRPiA = (m1 + NumGroups) / (m1 + m2 + NumGroups);
    }  
    }
  }
  return(1);
 }	 
/*********************************
// Don't Need until I active OnLARS
//////////////////////////////////////////////////////////
//   Set DiagXTX
//     DiagXTX is a vector containing only the diagonals of 
//     the GRAM matrix XTX.  Since this may be used for weighting
//     Procedures independent of GramXTX, it may be useful to calculate
//     simply the diagonal alone, instead of full matrix, unless of course
//     the matrix itself is already calculated.
int TwoLassoSexp::SetDiagXTX() {
  
  if (this->OnLARS == NULL) {
    Rprintf("ConvergentLARS:: SetDiagXTX, OnLARS NULL, why bother???\n"); 
    SuccessFlag = -1; R_FlushConsole(); R_ProcessEvents(); return(-1);
  }
  if (DiagXTX != NULL) { FFree(&DiagXTX); }
	DiagXTX = (double *) Calloc( this->p+1, double);
	if (DiagXTX == NULL) {
		Rprintf((char*) "Error DiagXTX is Null, cannot apportion Space \n");
		R_FlushConsole();
		return(-1);
  }
  int ii;
	if (SXtX != NULL)	 {
		ii = 0;
		for (ii = 0; ii < p; ii++) {
		      DiagXTX[ii] = REAL(SXtX)[ (ii) * p + ii ];	
    }
    return(1);
  }
  
  int cn1, kk;
  int One = 1;  double OneD = 1.0; int Zero = 0; double ZeroD = 0.0;
  if (iiWeights == NULL) {
    for (ii = 0; ii < p; ii++ ) {
      DiagXTX[ii] = F77_CALL(ddot)(&n, REAL(Sxxs) + ii * n, &One, 
        REAL(Sxxs) + ii * n, &One);
    }
  } else {
    if (p < n) {
      for (ii = 0; ii < p; ii++) {
        DiagXTX[ii] = 0; cn1 = n * ii;
        for (kk = 0; kk < n; kk++) {
          DiagXTX[ii] += REAL(Sxxs)[cn1] * REAL(Sxxs)[cn1] * iiWeights[kk];
          cn1++;
        }
      }
    } else {
      F77_CALL(dscal)(&p, &ZeroD, DiagXTX, &One);
      for (ii = 0; ii < n; ii++) {
        F77_CALL(dsbmv)("U", &p, &Zero, iiWeights + ii, REAL(Sxxs) + ii * n, &One,
          REAL(Sxxs) + ii * n, &One, &OneD, DiagXTX + ii, &One);          
      }
    }
  }

  return(1);
}
***********/


int TwoLassoSexp::RefreshTD2Lasso() {
  int ii;  double ResidOnii;  int One = 1; 
  double InvSigmaSq = 1.0/ OnSigma;
  int jj;
  if (CDO->YY == NULL || CDO->XX == NULL) {
    error("RefreshTD2Lasso: I can't do that I have NULL values\n");
  }
  if (CDO->XTXFlag == 2 && CDO->OnKappaS > 0 && 
    CDO->OnKappaS < n /2  && CDO->OnKappaS < p / 2) {
    F77_CALL(dcopy)(&n, CDO->YY, &One, iiWeights, &One);
    for (jj = 0; jj < CDO->OnKappaS; jj++) {
      if (CDO->OnBeta[CDO->kXFinder[jj]] != 0.0) {
        F77_CALL(daxpy)(&n, &CDO->OnBeta[CDO->kXFinder[jj]], 
          CDO->XX + n * CDO->kXFinder[jj], &One, 
          iiWeights, &One);
      }
    }
    for (ii = 0; ii < n; ii++) {
      iiWeights[ii] = (TDFNu + 1) / 
        (TDFNu + iiWeights[ii] * iiWeights[ii] * InvSigmaSq);
    } 
  } else {
    for (ii = 0; ii < n; ii++) {
      //Rprintf("RTD2Lasso: Updating Coordinate[ii=%d] \n", ii); R_FlushConsole();
      ResidOnii = CDO->YY[ii] - F77_CALL(ddot)(&p, CDO->XX + ii, &n, CDO->OnBeta, &One);
      //Rprintf("RTD2Lasso: Updated Coordinate[ii=%d]  to ResidOnii = %.4f\n", ii, ResidOnii); R_FlushConsole();
      iiWeights[ii] = (TDFNu +1)  / ( TDFNu + ResidOnii * ResidOnii * InvSigmaSq); 
      //Rprintf("RTD2Lasso: Updated Coordinate[ii=%d]  to iiWeights[%d] = %f\n",
      //    ii, ii, iiWeights[ii]); R_FlushConsole();      
    }
  }
  if (Verbose > 3) {
    Rprintf("RefreshTD2Lasso: Got to End of updating all iiWeights\n");
  }  
  //Rprintf("RefreshTD2Lasso: Going to Update Weights now on CDO I Promise!!!\n"); 
  SuccessFlag = CDO->RefreshWeights(iiWeights);
  if (SuccessFlag < 0) {
    Rprintf("RefreshTD2Lasso: Refresh Weights Failed: Abort\n");
    R_FlushConsole();
    error("RefreshTD2Lasso: error");
  }
  return(SuccessFlag);
}


///////////////////////////////////////////////////////////////////////////////
//   InitFixKa:
//       Useful for ``Fermi-Dirac'' Approximation Two-Lasso
//        This sets up PiRecVec to record changes to Pi necessary
//        to define a fixed KActive set as  algorithm progresses.                                                       
//
int TwoLassoSexp::InitFixKa(int rFixKa) {
	FixKa = rFixKa;
	OnPiA = ((double)FixKa )/((double) p);
  return(1);	  
}


////////////////////////////////////////////////////////////////////////////////
//      AofX:
//          Used in ``Fermi-Dirac'' Approximation Two-Lasso
//          Sums exponentiated Logit probabilities, 
//           necessary for robust computation
double TwoLassoSexp::AofX(double Onmu) {
	double logLambdasRat = log( OnLambdaA / OnLambdaD);
	double LambdaDiff = OnLambdaD-OnLambdaA;
	double RetTotal = 0;
	int ii; 
	for (ii = 0; ii < p; ii++) {
    if (-logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu < - 35) {
		  RetTotal += 1;
    } else if (-logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu > 35)  {
	    RetTotal += exp(- log( 1.0 + 
	      exp( -logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu)) );
    } else {
      RetTotal += exp(- log(1.0 + 
        exp( -logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu)) );		
    }
  }
  return(RetTotal);
}

////////////////////////////////////////////////////////////////////////////////
//      InsertBBOn1ofX:
//          Used in ``Fermi-Dirac'' Approximation Two-Lasso
//          Defines new BBOn1s after declaration of chemical fugacity Onmu. 
int TwoLassoSexp::InsertBBOn1ofX(double Onmu) {
	double logLambdasRat = log( OnLambdaA / OnLambdaD);
	double LambdaDiff = OnLambdaD-OnLambdaA;
	int ii; 
	for (ii = 0; ii < p; ii++) {
	  if (-logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu < - 35) {
		   REAL(SBBOn1)[ii] = 1;
    } else if (-logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu > 35)  {
	     REAL(SBBOn1)[ii] = exp(-log(1 + 
	       exp(-logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu)) );
    }  else {
       REAL(SBBOn1)[ii] = exp(-log( 1 +
         exp(-logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu)) );		
    }
  }
  return(1);
}

////////////////////////////////////////////////////////////////////////////////
//      DAofX:
//          Used in ``Fermi-Dirac'' Approximation Two-Lasso
//          Calculates local derivative of AofX as function of fugacity Onmu
double TwoLassoSexp::DAofX(double Onmu) {
	double logLambdasRat = log(OnLambdaA / OnLambdaD);
	double LambdaDiff = OnLambdaD-OnLambdaA;
	double RetTotal = 0;
	double LDenom = 0;
	int ii; 
	for (ii = 0; ii < p; ii++) {
    if (-logLambdasRat - LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu < - 35) {
		  RetTotal += exp(-logLambdasRat 
        -LambdaDiff * fabs((REAL(SOnBeta))[ii]) - Onmu);
    } else if ( -logLambdasRat - LambdaDiff * 
      fabs(REAL(SOnBeta)[ii]) - Onmu > 35)  {
      LDenom = - (log( 1 + exp(-logLambdasRat - LambdaDiff * 
        fabs(REAL(SOnBeta)[ii]) - Onmu)));	
      RetTotal += exp(LDenom);	
    } else {
      LDenom = - (log( 1 + exp(-logLambdasRat - 
        LambdaDiff * fabs(REAL(SOnBeta)[ii]) - Onmu)));	
      RetTotal += exp(-logLambdasRat - LambdaDiff * 
        fabs(REAL(SOnBeta)[ii]) - Onmu + LDenom + LDenom);		
    }
 }
 return(RetTotal);
}

////////////////////////////////////////////////////////////////////////////////
//      SolveFormu:
//          Used in ``Fermi-Dirac'' Approximation Two-Lasso
//          Solves for best fugacity Onmu that fixes Sum OnBBOn1 = kActive
double TwoLassoSexp::SolveFormu(double GoalTotal, 
         double SuffEpsilon, int MaxIters, double StartOnmu) {
	double Onmu = StartOnmu;
	double Prevmu = StartOnmu;
	double OnTotal = AofX(Onmu);
	double PrevTotal = OnTotal;
	
	double OnDerA;
	double ChangeMe;
	int ii; 
	if (Verbose > 5) {
    Rprintf("...............................................................\n");
		Rprintf("..... ATwoLassoObject.cc::SolveFormu(tt1=%d, tt2=%d): ", tt1, tt2);
    Rprintf(" Starting with StartOnmu = %.4e, OnTotal = %.4e, ",
      StartOnmu, OnTotal);
    Rprintf(" GoalTotal Will be %.4e, PrevTotal = %f, Prevmu = %f\n", GoalTotal,
      PrevTotal, Prevmu);
		R_FlushConsole(); R_ProcessEvents();
  }
	for (ii = 0; ii < MaxIters; ii++) {
		if (fabs(OnTotal - GoalTotal) < SuffEpsilon)  {
			return(Onmu);
    }
		OnDerA = DAofX(Onmu);
		ChangeMe = - (OnTotal-GoalTotal) / OnDerA;
		 
	  if (Verbose > 5) {
		  Rprintf(".... SolveFormu(): ii = %d, Onmu = %.4e, OnTotal = %.4e, ",
		    (int) ii, (double) Onmu, (double) OnTotal);
      Rprintf(" OnDerA = %.4e, ChangeMe = %.4e\n",
		    (double) OnDerA, (double) ChangeMe);
		  R_FlushConsole(); R_ProcessEvents();
    }
		if ( ChangeMe < 0 && OnTotal < GoalTotal) {
      if (Verbose > 5) {
		    Rprintf(".... SolveFormu: ii = %d, ChangeMe = %.4e < 0, but OnTotal = %.4e",
		      (int) ii, (double) ChangeMe , (double) OnTotal);
        Rprintf(" < GoalTotal %.4e, will do %.4e\n",
		      (double) GoalTotal, -.1 *  (OnTotal-GoalTotal) / fabs(OnDerA));
		    R_FlushConsole(); R_ProcessEvents();
      }
			ChangeMe = -.1 *  (OnTotal-GoalTotal) / fabs(OnDerA);
    } else if (ChangeMe > 0 && (OnTotal > GoalTotal)) {
      if (Verbose > 5) {
		    Rprintf("SolveFormu: ii = %d, ChangeMe = %.4e > 0, but OnTotal = %.4e",
		      (int) ii, (double) ChangeMe , (double) OnTotal);
        Rprintf(" > GoalTotal %.4e, will do %.4e\n",
		      (double) GoalTotal, -.1 *  (OnTotal-GoalTotal) / fabs(OnDerA));
		    R_FlushConsole(); R_ProcessEvents();
      }	         
	    ChangeMe = -.1 *  (OnTotal-GoalTotal) / fabs(OnDerA);
    }
    Prevmu = Onmu;
		Onmu = Onmu + ChangeMe;
		PrevTotal = OnTotal;
		OnTotal = AofX(Onmu);		
		if (!R_FINITE(Onmu) || ISNAN(Onmu)) {
			return(Onmu);
    } 
  }
	return(Onmu);
}
int TwoLassoSexp::FitBBOnWithGroups() {
  //VECTOR_ELT(sCodaChains,0);
 int NumElts = Rf_length(sGroupsSexp);
 SEXP GroupIterList;
 int ii;
 double piAFactor = log(OnPiA) - log(1.0-OnPiA);
 double LambdaFactor = log(OnLambdaA) - log(OnLambdaD);
 double BBHatEst = 0;
 double LDf = OnLambdaD-OnLambdaA;
 int jj;
 for (ii = 0; ii < NumElts; ii++) {
   GroupIterList = VECTOR_ELT(sGroupsSexp,ii);
   BBHatEst = piAFactor;
   if (isInteger(GroupIterList)) {
     for (jj = 0; jj < Rf_length(GroupIterList); jj++) {
       BBHatEst += LambdaFactor + 
         LDf * fabs(REAL(SOnBeta)[INTEGER(GroupIterList)[jj]]);
     }
   }  else if (Rf_isReal(GroupIterList)) {
     for (jj = 0; jj < Rf_length(GroupIterList);jj++) {
       BBHatEst += LambdaFactor + 
         LDf * fabs(REAL(SOnBeta)[(int) REAL(GroupIterList)[jj]]);
     }      
   }
   if (BBHatEst > 15) {BBHatEst = 1.0;}
   else if (BBHatEst < -15) { BBHatEst = exp(BBHatEst);}
   else { BBHatEst = exp(BBHatEst)/(1.0 + exp(BBHatEst));}
   if (isInteger(GroupIterList)) {
     for (jj = 0; jj < Rf_length(GroupIterList);jj++) {
       REAL(SBBOn1)[INTEGER(GroupIterList)[jj]] = BBHatEst;
     }
   }  else if (Rf_isReal(GroupIterList)) {
     for (jj = 0; jj < Rf_length(GroupIterList);jj++) {
       REAL(SBBOn1)[(int) REAL(GroupIterList)[jj]] = BBHatEst;
     }     
   }
 }
 return(1);
}

int TwoLassoSexp::FixKCalculateBBOn() {
  if (FixKa > p || FixKa < 0) {
    Rprintf("----- ATwoLassoObject2014.cc:: FixKCalculateBBOn(tt1=%d,tt2=%d): \n",tt1, tt2);
	  Rprintf("----- FixKCalculateBBOn:  Cannot make this specious Kcount = %d, p=%d \n", 
      FixKa, p);
  }
  if (Verbose > 3) {
    Rprintf("----- ATwoLassoObject2014.cc:: FixKCalculateBBOn(tt1=%d,tt2=%d): \n",tt1, tt2);
    Rprintf("-----    We're Running FixKCalculateBBOn with FixKa = %d, p=%d\n", 
      FixKa, p); R_FlushConsole();
  }
  int ii;
	int SuccFlag = 0;
  if (SuccFlag  == -1) {
  }
	double Solvedmu = 0;    
  double CurrentCount = 0; 
	double LogitOnPiA;
	if (OnPiA > 0.0 && OnPiA < 1.0) {
	  LogitOnPiA = log(OnPiA / (1-OnPiA));
  } else if (OnPiA == 1.0) {
	   LogitOnPiA = 20;
  } else if (OnPiA == 0.0) {
	  LogitOnPiA = -20;
  } else {
	  LogitOnPiA = 0;
  }
	double OnCurrentPiACount = AofX(LogitOnPiA);
	double ApparentMax = 0;
	double ApparentMin  = 0;
	double CountAtMax = 0;
	double CountAtMin = 0;
	int BreakItYouBoughtIt = 0;
 
	if (Verbose > 4) {
		Rprintf(" \n \n#### FixKCalculateBBOn(): tt1 = %d, tt2 = %d, ", tt1, tt2);
    Rprintf(" OnlambdaA = %.4e, OnlambdaD = %.4e \n",
			(double) OnLambdaA, (double) OnLambdaD);
		Rprintf("#### OnGammas  = "); PrintVector(CDO->OnGammas, p);
		Rprintf("\n#### OnBetas = "); PrintVector(REAL(SOnBeta), p);
		Rprintf("\n####  Now to Run Calculations \n");
		R_FlushConsole(); R_ProcessEvents();
  }
	if (OnLambdaA == OnLambdaD) {
		OnPiA = ((double)FixKa / ((double)p));
		if (OnPiA <= 0 || OnPiA >= 1) {
			Rprintf("#### 2L:FixKCalculateBBOn(), we have a ppiuse = %.4e Error, ", OnPiA);
      Rprintf("#### 2L:quitting \n", OnPiA);
			R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return(-1);
    } 
		LogitOnPiA = log(OnPiA / (1-OnPiA));
		Solvedmu = LogitOnPiA;
		CurrentCount = AofX(LogitOnPiA);
		if (Verbose > 4) {
			Rprintf("#### 2L: FixKCalculateBBOn(): tt1 = %d, tt2 = %d, OnlambdaA = %.4e, ",
        tt1, tt2, (double) OnLambdaA);
      Rprintf(" OnLambdaD = %.4e \n", (double) OnLambdaD);
			Rprintf("#### 2L: Started at Default, ppiuse = %.4e, LogitOnppiuse = %.4e, ",
			  (double) OnPiA, (double) LogitOnPiA);
      Rprintf("#### 2L: CurrentCount = %.4e \n", (double) CurrentCount);
		  R_FlushConsole();
    }
  } else {       
	  if (Verbose > 4) {
			Rprintf("#### 2L:FixKCalculateBBOn(): Starting to Solve: ppiuse = %.4e,", OnPiA);
      Rprintf("#### 2L: OnCurrentPiACount = %.4e \n",  OnCurrentPiACount);
			R_FlushConsole(); R_ProcessEvents();
		}
		Solvedmu = SolveFormu(FixKa, SUFFEPSILON_FD, MAXITERS_FD, LogitOnPiA);   
	    CurrentCount = AofX(Solvedmu);
		if (Verbose > 4) {
		  Rprintf("##############################################\n");
		  R_FlushConsole();
			Rprintf("#### 2L:FixKCalculateBBOn(): Initial Solvedmu = %.4e, ", Solvedmu);
      Rprintf("#### CurrentCount = %.4e \n", (double) CurrentCount);
		  R_FlushConsole(); R_ProcessEvents();
		}	        
		if (ISNAN(Solvedmu) ||  !R_FINITE(Solvedmu) ||
		  (fabs( ((double) FixKa )- CurrentCount ) > SUFFEPSILON_FD * 10 )  ) {			            
	    ApparentMax = 60;
	    ApparentMin  = -60;
	    CountAtMax = AofX(ApparentMax);
	    CountAtMin = AofX(ApparentMin);
	    if (OnCurrentPiACount < CountAtMax && OnCurrentPiACount > (double) FixKa) {
		    ApparentMax = LogitOnPiA;
		    CountAtMax = OnCurrentPiACount;
      }  
	    if (OnCurrentPiACount > CountAtMin && OnCurrentPiACount < (double) FixKa) {
		    ApparentMax = LogitOnPiA;
		    CountAtMax = OnCurrentPiACount;
      } 
      if (CountAtMax < FixKa) {
        if (Verbose > 1) {
          Rprintf("####   CountAtMax=%d < FixKa = %d.\n", CountAtMax, FixKa); R_FlushConsole();
	        Rprintf("#### 2L: FindBBOn: Apparently FixKa is already less than Max.\n");
          Rprintf("####   Probabilities to big to tweak.\n");  R_FlushConsole();
	        R_FlushConsole(); R_FlushConsole();
	      }
	      LogitOnPiA = ApparentMax; OnCurrentPiACount = CountAtMax;
      } else if (CountAtMin > FixKa) {
        if (Verbose > 1) {
          Rprintf("#### 2L: FindBBOn: CountAtMin=%d > FixKa = %d.\n", CountAtMin, FixKa); R_FlushConsole();
	        Rprintf("#### 2L:FindBBOn: Apparently FixKa is already less than Max.\n");
          Rprintf("####  Probabilities to big to tweak.\n");
	        R_FlushConsole(); R_FlushConsole();
	      }
	      LogitOnPiA = ApparentMin; OnCurrentPiACount = CountAtMin;	      
      } else {	
	      BreakItYouBoughtIt = 0;	
	      OnCurrentPiACount = AofX(LogitOnPiA);            
			  for (ii = 0; ii < MAXITERS_FD; ii++) {			       
				  if ( fabs( ((double) FixKa )- OnCurrentPiACount ) < 
            SUFFEPSILON_FD * 10) {
					  Solvedmu = LogitOnPiA; 
            CurrentCount = OnCurrentPiACount;
					  BreakItYouBoughtIt = 1;
					  break;    
				  }
				  if ( ISNAN(Solvedmu) || !R_FINITE(Solvedmu) ||
				    fabs(((double) FixKa )- CurrentCount ) >
				    fabs(((double) FixKa) - OnCurrentPiACount) ) {
					  if (FixKa > OnCurrentPiACount) {
						  LogitOnPiA = (ApparentMax + LogitOnPiA) / 2.0 ;
						  OnCurrentPiACount = AofX(LogitOnPiA);
						  if (OnCurrentPiACount < FixKa) {
							  ApparentMin = LogitOnPiA;
							  CountAtMin = OnCurrentPiACount;
						  } else if (OnCurrentPiACount > FixKa) {
							  ApparentMax = LogitOnPiA;
							  CountAtMax = OnCurrentPiACount;							       
						  }
						  Solvedmu = LogitOnPiA; 			       
					  } else {
						  LogitOnPiA = (ApparentMin + LogitOnPiA) / 2.0 ;
						  OnCurrentPiACount = AofX(LogitOnPiA);
						  if (OnCurrentPiACount < FixKa) {
							  ApparentMin = LogitOnPiA;
							  CountAtMin = OnCurrentPiACount;
						  } else if (OnCurrentPiACount > FixKa) {
							  ApparentMax = LogitOnPiA;
							  CountAtMax = OnCurrentPiACount;							       
						  }						       
						  Solvedmu = LogitOnPiA	;
					  }         
					  Solvedmu = SolveFormu(FixKa,SUFFEPSILON_FD, 
              MAXITERS_FD, LogitOnPiA);   
				  } else {
            if (CurrentCount < FixKa) {
							ApparentMin = Solvedmu;
							CountAtMin = CurrentCount;
						} else if (CurrentCount > FixKa) {
							ApparentMax = Solvedmu;
							CountAtMax = CurrentCount;							       
						}
					  LogitOnPiA = Solvedmu; OnCurrentPiACount = CurrentCount;
					  Solvedmu = SolveFormu(FixKa,SUFFEPSILON_FD, 
              MAXITERS_FD, LogitOnPiA);   
			    }
			    CurrentCount = AofX(Solvedmu);
			    if (!ISNAN(Solvedmu) && R_FINITE(Solvedmu) && 
			      fabs( ((double) FixKa )- CurrentCount ) < SUFFEPSILON_FD * 10) {
				    LogitOnPiA = Solvedmu;
				    OnCurrentPiACount = CurrentCount;
				    BreakItYouBoughtIt = 2;
				    break;
			    }
			    if (LogitOnPiA > 50) {
				    Solvedmu = LogitOnPiA;
				    CurrentCount = OnCurrentPiACount;
				    BreakItYouBoughtIt = 3;
				    break;   
			    } else if (LogitOnPiA < -50) {
				    Solvedmu = LogitOnPiA;
				    CurrentCount = OnCurrentPiACount;
				    BreakItYouBoughtIt = 4;
				    break;   				           			           
			    }
		    }
		    if (ISNAN(LogitOnPiA) || !R_FINITE(LogitOnPiA) )  {
			    Rprintf("#### 2L: FixKCalculateBBOn ISNAN(LogitOnPiA) failed, \n");
          Rprintf("####    we've got FixKa = %d, p = %d, ", FixKa, p);
          Rprintf("####    LogitOnPiA = %.4e, piLogitOnPiA = %.4e \n",
			      LogitOnPiA, exp(LogitOnPiA) / ( 1 + exp(LogitOnPiA)) );
			    Rprintf("####    ii = %d \n", ii); R_FlushConsole();
			    Rprintf("####    BreakItYouBoughtIt = %d\n", BreakItYouBoughtIt);			                  
			    R_FlushConsole(); R_ProcessEvents();
			    SuccessFlag = -1; return(-1);
		    }		           
		    if (fabs( ((double) FixKa )- OnCurrentPiACount ) > SUFFEPSILON_FD * 4000
		      && LogitOnPiA < 50 && LogitOnPiA > -50) {
			    Rprintf("#### 2L: FixKCalculateBBOn we failed, we've got FixKa = %d, ", FixKa);
          Rprintf("####    p = %d, LogitOnPiA = %.4e, piLogitOnPiA = %.4e \n",
			      p, (double) LogitOnPiA,(double) exp(LogitOnPiA) /
            ( 1 + exp(LogitOnPiA)) ); R_FlushConsole();
			    Rprintf("#### 2L: OnCurrentPiACount = %.4e, CurrentCount = %.4e \n",
			      OnCurrentPiACount, CurrentCount);
			    Rprintf("#### 2L: Solvedmu = %.4e, piSolvedmu = %.4e \n", 
			      Solvedmu, exp(Solvedmu) / ( 1.0 + exp(Solvedmu)));
			    Rprintf("#### 2L: ApparentMax = %.4e, CountAtMax = %.4e \n", 
			      ApparentMax, CountAtMax);	R_FlushConsole();
			    Rprintf("#### 2L: ApparentMin = %.4e, CountAtMin = %.4e \n", 
			      ApparentMin, CountAtMin);	R_FlushConsole();		               		               
			    Rprintf("#### 2L: ii = %d \n", ii);  R_FlushConsole();
			    Rprintf("#### 2L: BreakItYouBoughtIt = %d\n", BreakItYouBoughtIt);
			    R_FlushConsole(); R_ProcessEvents();
			    SuccessFlag = -1; return(-1);
		    }
	    }	           
		} else {
			LogitOnPiA = Solvedmu;
			OnCurrentPiACount = CurrentCount;
		}
		if (ISNAN(LogitOnPiA) || !R_FINITE(LogitOnPiA)) {
		  Rprintf("###################################################### \n");
		  Rprintf("#### 2L  Error on Logit \n"); R_FlushConsole();
			Rprintf("#### 2L: FixKCalculateBBOn: LogitOnPiA ISNAN!!, ");
			R_FlushConsole();
      Rprintf("#### 2L: LogitOnPiA = %.4e, FixKa = %d \n",  LogitOnPiA, FixKa);
			SuccessFlag = -1; return(-1);
		}		       
		if (LogitOnPiA > 50) {
			OnPiA = 1.0;
		} else if (LogitOnPiA < -50) {
			OnPiA= exp(LogitOnPiA);
		} else { 
			OnPiA = exp(LogitOnPiA) / ( 1.0 + exp(LogitOnPiA));
		}
  }
	if (Verbose > 3) {
	  Rprintf("#### 2L: tt1=%d, tt2=%d, FixKCalculateBBOn settles upon ", tt1,tt2);
    Rprintf("#### 2L: Solvedmu = %.4e, LogitOnppiuse = %.4e, ppiuse = %.4e, ",
      Solvedmu, LogitOnPiA, OnPiA);   R_FlushConsole();
    Rprintf("#### 2L: CurrentCount = %.4e\n", CurrentCount);
	  R_FlushConsole(); R_ProcessEvents();
  }
  if (Rf_length(SPiA) >= tt1+1)  {
    REAL(SPiA)[tt1] = OnPiA;
  }
  if (Rf_length(SPiA) >= tt1+2) {
    REAL(SPiA)[tt1+1] = OnPiA;
  }
  if (Verbose > 3) {
    Rprintf("#### 2L: FixBBOn: We reached near the end, now Inser %f\n",
      LogitOnPiA); R_FlushConsole();
  }
	SuccFlag = InsertBBOn1ofX(LogitOnPiA);                   
  if (Verbose > 3) {
    Rprintf("#### 2L: FixBBOn: Really at the end\n"); R_FlushConsole();
  }
	return(1);
}



 
  
#define TestVector(RVector, SVector, CharVector)                               \
  NewRVector = NULL;                                                           \
  if (DoP >= 2 || Verbose > 2) {                                               \
    Rprintf("TestVector: Testing %s\n", CharVector); R_FlushConsole();         \
  }                                                                            \
  if (RVector != NULL ) {                                                      \
    if (Rf_isReal(RVector->asSexp())) {  \
      NewRVector = new AObject(Rf_allocVector(REALSXP,           \
        Rf_length(RVector->asSexp())));                                \
      for (int ii = 0; ii < Rf_length(NewRVector->asSexp()); ii++) {   \
        REAL(NewRVector->asSexp())[ii] = REAL(RVector->asSexp())[ii];  \
      }                                                                \
    } else  if (Rf_isReal(RVector->asSexp())) {                        \
      NewRVector = new AObject(Rf_allocVector(INTSXP,              \
      Rf_length(RVector->asSexp())));                                          \
      for (int ii = 0; ii < Rf_length(NewRVector->asSexp()); ii++) {           \
        INTEGER(NewRVector->asSexp())[ii] = INTEGER(RVector->asSexp())[ii];    \
      }                                                                \
    } else if (Rf_isNull(RVector->asSexp())) {                         \
      Rprintf(" -- Hey, %s is NIL Sexp? \n", CharVector); R_FlushConsole(); \
      NewRVector = new AObject(R_NilValue);                           \
    }                                                                       \
    delete(RVector);  SVector = NULL;                                       \
    RVector = NewRVector;  SVector = NewRVector->asSexp();                  \
  }                                                                         \
  NewRVector = NULL
  
#define TestRVector(RVector, CharVector)                                       \
  NewRVector = NULL;                                                           \
  if (DoP >= 2 || Verbose > 2) {                                               \
    Rprintf("TestVector: Testing %s\n", CharVector); R_FlushConsole();         \
  }                                                                            \
  if (RVector != NULL ) {                                                      \
    if (Rf_isReal(RVector->asSexp())) {  \
      NewRVector = new AObject(Rf_allocVector(REALSXP,           \
        Rf_length(RVector->asSexp())));                                \
      for (int ii = 0; ii < Rf_length(NewRVector->asSexp()); ii++) {   \
        REAL(NewRVector->asSexp())[ii] = REAL(RVector->asSexp())[ii];  \
      }                                                                \
    } else  if (Rf_isReal(RVector->asSexp())) {                        \
      NewRVector = new AObject(Rf_allocVector(INTSXP,              \
      Rf_length(RVector->asSexp())));                                          \
      for (int ii = 0; ii < Rf_length(NewRVector->asSexp()); ii++) {   \
        INTEGER(NewRVector->asSexp())[ii] = INTEGER(RVector->asSexp())[ii];  \
      }                                                                \
    } else if (Rf_isNull(RVector->asSexp())) {                         \
      Rprintf(" -- Hey, %s is NIL Sexp? \n", CharVector); R_FlushConsole(); \
      NewRVector = new AObject(R_NilValue);                           \
    }                                                                       \
    delete(RVector);                                                        \
    RVector = NewRVector;                                                   \
  }                                                                         \
  NewRVector = NULL
  
  
#define TestMatrix(RMatrix, CharMatrix)                                     \
  NewRMatrix = NULL;  NewRDim = NULL;                                       \
  if (DoP >=2 || Verbose >= 2)  {                                           \
    Rprintf("CVDataIntegrity(jti=%d), now testing matrix %s \n",            \
      jti, CharMatrix);                                                     \
    R_FlushConsole();                                                       \
  }                                                                         \
  if (RMatrix == NULL || Rf_isNull(RMatrix->asSexp())) {                    \
  }  else if (Rf_isNull(Rf_getAttrib(RMatrix->asSexp(), R_DimSymbol))) {    \
    Rprintf("  -- CVData: Actually %s is not a matrix. \n", CharMatrix);    \
    R_FlushConsole();                                                       \
  } else {                                                                  \
    NewRDim = Rf_getAttrib(RMatrix->asSexp(), R_DimSymbol);                 \
    if (Rf_isInteger(RMatrix->asSexp())) {                                  \
       NewRMatrix = new AObject(Rf_allocMatrix(INTSXP,                \
         INTEGER(NewRDim)[0],                                               \
         INTEGER(NewRDim)[1]));                                             \
       for (int ii = 0; ii < Rf_length(RMatrix->asSexp()); ii++) {          \
         INTEGER(NewRMatrix->asSexp())[ii] = INTEGER(RMatrix->asSexp())[ii];\
       }                                                                    \
       delete(RMatrix); RMatrix = NULL;                                     \
       RMatrix = NewRMatrix;                                                \
    } else if (Rf_isReal(RMatrix->asSexp())) {                              \
       NewRMatrix = new AObject(Rf_allocMatrix(REALSXP,               \
         INTEGER(NewRDim)[0],                                               \
         INTEGER(NewRDim)[1]));                                             \
       for (int ii = 0; ii < Rf_length(RMatrix->asSexp()); ii++) {          \
         REAL(NewRMatrix->asSexp())[ii] = REAL(RMatrix->asSexp())[ii];      \
       }                                                                    \
       delete(RMatrix); RMatrix = NULL;                                     \
       RMatrix = NewRMatrix;                                                \
    }                                                                       \
  }                                                                         \
  NewRMatrix = NULL



int TwoLassoSexp::CVDataIntegrity(int DoP) {
   AObject *NewRVector = NULL;
   AObject *NewRMatrix = NULL;
   SEXP NewRDim = R_NilValue;
   if (DoP >= 1) {
     Rprintf("CVDataIntegrity(jti=%d): Testing CVDataIntegrity. \n", jti); R_FlushConsole();
   }
   TestVector(RListFlags, ListFlags, "ListFlags");
   TestVector(RSPiAVectorInputs, SPiAVectorInputs, "SPiAVectorInputs");
   TestVector(RSLambdaAKInputs, SLambdaAKInputs, "SLambdaAKInputs");
   TestVector(RSSigmaVectorInputs, SSigmaVectorInputs, "RSSigmaVectorInputs");
   TestVector(RSLambdaDKInputs, SLambdaDKInputs, "RSLambdaDKInputs");
   if (iiWeights != NULL) {
     TestVector(RSiiWeights, SiiWeights, "RSiiWeights");
   }
   TestRVector(RSRecordBeta0CV, "RSRecordBeta0CV");
   TestMatrix(RSStartBetaMatrix, "RSStartBetMatrix");
   if (RSStartBetaMatrix != NULL) {
     SStartBetaMatrix = RSStartBetaMatrix->asSexp();
   }
   TestMatrix(RSRecordBetaCV, "RSRecordBetaCV");
   if (RSRecordBetaCV != NULL) {
     SRecordBetaCV = RSRecordBetaCV->asSexp();
   }

   if (GLMCDO != NULL) {
     if (Verbose >= 2 || DoP >= 2) {
       Rprintf("CVCheckIntegirty for GLMCDO Object. \n"); R_FlushConsole();
     }
     TestRVector(RBetasPrev, "RBetasPrev"); 
     GLMCDO->pBetasPrev = REAL(RBetasPrev->asSexp());
     if (RBeta0Records != NULL) {
     if (Verbose >= -1) {
       Rprintf("CVCheckIntegrity: Beta0Records about to analyze"); R_FlushConsole();
       Rprintf(", length %d\n", Rf_length(RBeta0Records->asSexp()));  R_FlushConsole();
     }
     NewRVector = new AObject(Rf_allocVector(REALSXP, Rf_length(RBeta0Records->asSexp())));
     for (int iti = 0; iti < Rf_length(RBeta0Records->asSexp()); iti++) {
       if (iti != 0 && iti % 20 == 0) {
         Rprintf("\n"); R_FlushConsole();
       }
       Rprintf(" - %d", iti); R_FlushConsole();
       REAL(NewRVector->asSexp())[iti] = REAL(RBeta0Records->asSexp())[iti];
     }
     if (Verbose >= -1) {
       Rprintf(" -- Completed Filling NewRVector, now delete !\n"); R_FlushConsole();
     }
     delete(RBeta0Records); RBeta0Records = NULL;
     RBeta0Records = NewRVector;
     if (Verbose >= -1) {
       Rprintf(" -- Now setting RBeta0Records to vector. \n"); R_FlushConsole();
     }
     GLMCDO->pBeta0Records = REAL(RBeta0Records->asSexp());
      TestRVector(RBeta0Records, "RBeta0Records"); 
     GLMCDO->pBeta0Records = REAL(RBeta0Records->asSexp());
     if (Verbose >= -1) {
       Rprintf("CVCheckIntegrity: Completed Beta0Records\n"); R_FlushConsole();
     }
     }
     TestRVector(RProbWeights, "RProbWeights");
     GLMCDO->nProbWeights = REAL(RProbWeights->asSexp());
     TestRVector(RLogitCurrentProb, "RLogitCurrentProb");
     GLMCDO->nlogitCurrentProb = REAL(RLogitCurrentProb->asSexp());
     TestRVector(RLogitPrevProb, "RLogitPrevProb");
     GLMCDO->nlogitCurrentProb = REAL(RLogitPrevProb->asSexp());
     TestRVector(RSOnGammas, "RSOnGammas");
     SOnGammas = RSOnGammas->asSexp();
     GLMCDO->OnGammas = REAL(RSOnGammas->asSexp());     
     TestRVector(RSBBOn1, "RSBBOn1");
     SBBOn1 = RSBBOn1->asSexp();
     TestRVector(RSOnBeta, "RSOnBeta");
     SOnBeta = RSOnBeta->asSexp();
     GLMCDO->OnBeta = REAL(RSOnBeta->asSexp());    
     if (RSRecOnBeta != NULL && !Rf_isNull(RSRecOnBeta->asSexp())) {
       if (Rf_isNull(Rf_getAttrib(RSRecOnBeta->asSexp(), R_DimSymbol))) {
         TestRVector(RSRecOnBeta, "RSRecOnBeta");
         SRecOnBeta = RSRecOnBeta->asSexp();
         GLMCDO->pRecordOnBetas = REAL(RSRecOnBeta->asSexp());
       } else {
         TestMatrix(RSRecOnBeta, "RSRecOnBeta");
         SRecOnBeta = RSRecOnBeta->asSexp();
         GLMCDO->pRecordOnBetas = REAL(RSRecOnBeta->asSexp());
       }
     }
     if (RSiiWeights != NULL) {
       TestRVector(RSiiWeights, "RSiiWeights");
       SiiWeights = RSiiWeights->asSexp();
       iiWeights = REAL(RSiiWeights->asSexp());
       GLMCDO->iiWeights = REAL(RSiiWeights->asSexp());
     }
     if (GLMCDO->iWeightedXtX != NULL) {
       if (DoP >= 2 || Verbose >= 2) {
         Rprintf("CVCheckIntegrity: checking iWeightedXtX\n"); R_FlushConsole();
       }
       int *NewiWeighted = (int*) Calloc(p, int);
       for (int iti = 0; iti < p; iti++) {
         NewiWeighted[iti] = GLMCDO->iWeightedXtX[iti];
       }
       Free(GLMCDO->iWeightedXtX);  GLMCDO->iWeightedXtX = NewiWeighted;
       
     }
   }
   
   if (Verbose >= 3 || DoP >= 2) {
     Rprintf(" -------------------------------------------------------------------- \n");
     Rprintf(" -- CVDataIntegrity(jti=%d, tt1=%d) - PASS   \n", jti, tt1);     
     Rprintf(" ---------------------------------------------------------------------\n");
     R_FlushConsole();
   }
   if (Verbose >= 1 || DoP >= 1) {
     Rprintf("  -- CVDataIntegrity: We have tested all integrity. \n"); R_FlushConsole();
   }
   return(1);
}
  

///////////////////////////////////////////////////////////////////////////////
// A PrintDiagnostic of current Info
int TwoLassoSexp::TwoLassoPrintDiagnosticA() {
	if (CDO != NULL) {
		Rprintf((char*) "TwoLasso:: Just Prior to Run CDO \n");
		Rprintf((char*) "TwoLasso:: tt1 = %d, tt2 = %d; sigmaNoiseSq  = %.8f\n",
			tt1, tt2, OnSigma);
    Rprintf("   Thus SSigma is: "); PrintVector(REAL(SSigma), Rf_length(SSigma));
    Rprintf("\n"); R_FlushConsole();
  }  	
  if (SXtX != NULL && !Rf_isNull(SXtX)) {		 
    Rprintf((char*) "  TwoLassoSXtX = \n"); R_FlushConsole();
    PrintRMatrix(REAL(SXtX), p, p);
    Rprintf("\n  CDOSXtX = ");  R_FlushConsole();
    PrintRpMatrix(CDO->pXTX, p,CDO->OnKappaMem); Rprintf("\n"); R_FlushConsole();
  }
  if (CDO->XTY == NULL) {
    Rprintf("Weird: XTY is NULL!\n"); R_FlushConsole();
  } else {
    Rprintf((char*) "TwoLasso::  CDO->XTY = \n");
    PrintVector(CDO->XTY, p);  Rprintf("\n"); R_FlushConsole();
  }
  
  if (CDO->OnBeta == NULL) {
    Rprintf("CDO->OnBeta is NULL!\n"); R_FlushConsole();
  } else {
    Rprintf((char*) "TwoLasso::  CDO->OnBeta = \n");
    PrintVector(CDO->OnBeta, p);
  }
  if (CDO->XTResid == NULL) {
    Rprintf("CDO->XTResid is NULL!\n"); R_FlushConsole();
  } else {
    Rprintf((char*) "TwoLasso::  CDO->XTResid = \n");
    PrintVector(CDO->XTResid, p);  Rprintf("\n"); R_FlushConsole();
  }		
  Rprintf((char*) "TwoLasso::  n = %d, p = %d \n", n, p); R_FlushConsole();
		
  if (CDO->OnGammas != NULL) {
    Rprintf((char*) "TwoLasso:: OnGammas = \n");
    PrintVector(CDO->OnGammas, p);
  }
  return(1);
}

int TwoLassoSexp::SetupGroups(SEXP GroupSexp) {
  if (GroupSexp == NULL || Rf_isNull(GroupSexp) || Rf_length(GroupSexp) <= 0) {
    return(1);
  }
  Rprintf("- ATwoLassoObject2014.cc::SetupGroups() -- Issue -- ");
  Rprintf("  Note: Setup Groups is Not Activated!\n");  R_FlushConsole();
  return(1);
}
int TwoLassoSexp::SetupTDFNu(SEXP rSTDFNu, SEXP rSiiWeights) {
  TDFNu = GetFirstDouble(rSTDFNu);  STDFNu = rSTDFNu;
  if (Verbose >= 1) {
    Rprintf("TwoLassoSexp::SetupTDFNu, TDFNu = %f \n", TDFNu); R_FlushConsole();
  }
  if (TDFNu > 0) {
    if (Rf_isNull(rSiiWeights) || Rf_length(rSiiWeights) != n) {
      RSiiWeights = new AObject(Rf_allocVector(REALSXP, n));
      SiiWeights = RSiiWeights->asSexp();
      iiWeights = REAL(SiiWeights);
      if (iiWeights == NULL) {
        Rprintf("TDFNu weights fail! \n"); R_FlushConsole();
      }
      for (int ii = 0; ii < n; ii++) {
        iiWeights[ii] = 1.0;
      }
    } else {
      SiiWeights = rSiiWeights;
      if (RSiiWeights != NULL) { 
        if (Verbose >= 3) {
          Rprintf("SetupTDFNu: Deleting Old RSiiWeights \n"); R_FlushConsole();
        }
        delete(RSiiWeights); 
      }
      RSiiWeights = new AObject(rSiiWeights);
      SiiWeights = RSiiWeights->asSexp();  iiWeights = REAL(SiiWeights);
    }
  } else {
    if (CDO != NULL) {
      CDO->TDFNoise = TDFNu;
      CDO->TDFSigma = -1.0;
    }
  }
  if (CDO != NULL) {
    CDO->SetupWeights(iiWeights);
    CDO->TDFNoise = TDFNu;  CDO->TDFSigma = OnSigma;
  }
  return(1);
}

int TwoLassoSexp::OnLARSCDODiagnostic() {
/*
	if (this->OnLARS != NULL) {
		Rprintf((char*)"\n  Now OnLARS->OnBetas \n");
		PrintVector(this->OnLARS->OnBetas, p);
	  R_FlushConsole(); Rprintf((char*)"\n");
		Rprintf((char*) "\n  OnLARS->OnNus = c( ");
		for (ii = 0; ii < this->p-1; ii++) {
			Rprintf((char*) "%.8f, ", (double)(double) this->OnLARS->OnBetas[ii]);
		}
		Rprintf((char*) " %.8f ) \n", this->OnLARS->OnBetas[p-1]);	          
		R_FlushConsole();
	}
  */ // If We implement OnLARS Version
	if (this->CDO != NULL) {
		Rprintf((char*)"\n  CDO->OnBeta\n");
	  R_FlushConsole(); PrintVector(CDO->OnBeta, p);
	}          
	Rprintf((char*)"  Now SOnBeta \n");
	R_FlushConsole(); PrintVector(REAL(SOnBeta), p);
	Rprintf((char*)"\n  Now BBOn1 \n");
  R_FlushConsole(); PrintVector(REAL(SBBOn1), p);  R_FlushConsole();
	return(1);
}

int TwoLassoSexp::RecordHistory() {
   int One = 1;  int MyPlace = tt1 * p;
   if (Verbose > 4) {
     Rprintf("#### ATwoLassoObject2014.cc::RecordHistory(tt1=%d) ", tt1);
     Rprintf(": Starting!\n"); R_FlushConsole();
   }
   if (SRecBBOn1 != NULL && 
     !Rf_isNull(SRecBBOn1) && Rf_length(SRecBBOn1) >= p * (tt1+1)) {
      F77_CALL(dcopy)(&p, REAL(SBBOn1), &One,
        REAL(SRecBBOn1) + MyPlace, &One);
   } 
   if (SRecOnBeta != NULL && 
     !Rf_isNull(SRecOnBeta) && Rf_length(SRecOnBeta) >= p * (tt1+1) ) {
      F77_CALL(dcopy)(&p, REAL(SOnBeta), &One,
        REAL(SRecOnBeta) + MyPlace, &One);
   }
   if (Verbose > 4) {
     Rprintf("#### RecordHistory: Finished\n"); R_FlushConsole();
   }
   return(1);
}


SEXP TwoLassoSexp::get_HPDQuantiles() {
  if (ConfidenceQuantiles == NULL) {
    return(R_NilValue);
  }
  SetHPDQuantiles();
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(REALSXP, Rf_length(HPDQuantiles->asSexp())));
  int One = 1;
  int Len = Rf_length(HPDQuantiles->asSexp());
  F77_CALL(dcopy)(&Len, REAL(HPDQuantiles->asSexp()), &One,
    REAL(sOut), &One);
  Rf_unprotect(1);
  return(sOut);
}
SEXP TwoLassoSexp::get_ConfidenceQuantiles() {
  if (ConfidenceQuantiles == NULL) {
    return(R_NilValue);
  }
  SEXP sOut = R_NilValue;
  Rf_protect(sOut = Rf_allocVector(REALSXP, Rf_length(ConfidenceQuantiles->asSexp())));
  int One = 1;
  int Len = Rf_length(ConfidenceQuantiles->asSexp());
  F77_CALL(dcopy)(&Len, REAL(ConfidenceQuantiles->asSexp()), &One,
    REAL(sOut), &One);
  Rf_unprotect(1);
  return(sOut);
}
void TwoLassoSexp::set_ConfidenceQuantiles(SEXP iConfidenceQuantiles) {
  if (Verbose >= 1) {
    Rprintf("ATwoLassoObject2014.cc::set_ConfidenceQuantiles() \n"); R_FlushConsole();
    Rprintf("** Setting Confidence Quantiles.\n");
    if (ConfidenceQuantiles != NULL) {
      Rprintf("** -- ISSUE -- Current ConfidenceQuantiles is null, will reset!\n");
    }
  }
  DDelete(ConfidenceQuantiles, "ConfidenceQuantiles");
  if (Rf_isNull(iConfidenceQuantiles)) {
    return;
  }
  if (!Rf_isReal(iConfidenceQuantiles)  || Rf_length(iConfidenceQuantiles) <= 0) {
    Rprintf("** set_ConfidenceQuantiles() ; Hey, y");
    return;
  }
  SEXP sSave = R_NilValue;
  int Len =  Rf_length(iConfidenceQuantiles);
  if (Verbose >= 2) {
    Rprintf("** set_ConfidenceQuantiles(): about to allocate vector of length %d\n",
      Len); R_FlushConsole();
  }
  Rf_protect(sSave = Rf_allocVector(REALSXP, Len));
  int One = 1;
  F77_CALL(dcopy)(&Len, REAL(iConfidenceQuantiles), &One,
    REAL(sSave), &One);  
  if (Verbose >= 2) {
    Rprintf("** set_ConfidenceQuantiles(): About to create RObject for sSave");
    R_FlushConsole();
  }
  ConfidenceQuantiles = new AObject(sSave); 
  Rf_unprotect(1);
  if (Verbose >= 2) {
    Rprintf("** set_ConfidenceQuantiles(): All the way finished\n");
    R_FlushConsole();
  }
  return;
}

SEXP TwoLassoSexp::get_HPDMatrix() {
  int One = 1;  int GoodFlag;
  int Alln= p * Rf_length(ConfidenceQuantiles->asSexp());
  if (Verbose >= 2) {
    Rprintf("-- ATwoLassoObject2014.cc:: get_HPDMatrix(): Starting!\n"); R_FlushConsole();
  }
  if (HPDMatrix == NULL || Rf_isNull(HPDMatrix->asSexp()) ||
    (!Rf_isNull(HPDMatrix->asSexp()) && 
    F77_CALL(dasum)(&Alln, REAL(HPDMatrix->asSexp()), &One) == 0.0)
    ) {
    if (Verbose >= 3) {
      Rprintf("--  get_HPDMatrix(): have to GenreteConfidenceMatrix\n");
      R_FlushConsole();
    }
    GoodFlag = GenerateConfidenceMatrix();
    if (GoodFlag < 0) {
      Rprintf("--   get_HPDMatrix(): Failure to generate good confidence matrix!\n");
      R_FlushConsole();
    }
  }
  if (HPDMatrix == NULL) {
    return(R_NilValue);
  }              
  return(HPDMatrix->asSexp());
}
SEXP TwoLassoSexp::get_ConfidenceMatrix() {
  int One = 1;  int GoodFlag;
  if (ConfidenceQuantiles == NULL) {
    Rf_error("get_ConfidenceMatrix: Error, there are no Confidence Quantiles. \n");
  }
  int Alln= p * Rf_length(ConfidenceQuantiles->asSexp());
  if (Verbose >= 2) {
    Rprintf("-- ATwoLassoObject2011.cc:: get_ConfidenceMatrix(): Starting!\n"); R_FlushConsole();
  }
  if (ConfidenceMatrix == NULL || Rf_isNull(ConfidenceMatrix->asSexp()) ||
    (!Rf_isNull(ConfidenceMatrix->asSexp()) && 
    F77_CALL(dasum)(&Alln, REAL(ConfidenceMatrix->asSexp()), &One) == 0.0)
    ) {
    if (Verbose >= 3) {
      Rprintf("--  get_ConfidenceMatrix(): have to GenreteConfidenceMatrix\n");
      R_FlushConsole();
    }
    GoodFlag = GenerateConfidenceMatrix();
    if (GoodFlag < 0) {
      Rprintf("-- get_ConfidenceMatrix() : Failure to generate good confidence matrix!\n");
      R_FlushConsole();
    }
  }
  if (ConfidenceMatrix == NULL) {
    return(R_NilValue);
  }
  return(ConfidenceMatrix->asSexp());
}
SEXP TwoLassoSexp::get_UnshrunkConfidenceMatrix() {
  int One = 1;  int GoodFlag;
  int AllIn = p *  Rf_length(ConfidenceQuantiles->asSexp());
 
  if (Verbose >= 2) {
    Rprintf("-- ATwoLassoObject2011.cc:: get_UnshrunkConfidenceMatrix() : Starting!\n"); R_FlushConsole();
  }
  if (UnshrunkConfidenceMatrix == NULL || Rf_isNull(UnshrunkConfidenceMatrix->asSexp()) ||
    (!Rf_isNull(ConfidenceMatrix->asSexp()) && 
    F77_CALL(dasum)(&AllIn, REAL(UnshrunkConfidenceMatrix->asSexp()), &One) == 0.0)
    ) {
    if (Verbose >= 3) {
      Rprintf("--   get_UnshrunkConfidenceMatrix() : have to GenreteConfidenceMatrix\n");
      R_FlushConsole();
    }
    GoodFlag = GenerateConfidenceMatrix();
    if (GoodFlag < 0) {
      Rprintf("--  get_UnshrunkConfidenceMatrix() : Failure to generate good confidence matrix!\n");
      R_FlushConsole();
    }
  }
  if (UnshrunkConfidenceMatrix == NULL) {
    return(R_NilValue);
  }
  return(UnshrunkConfidenceMatrix->asSexp());
}
SEXP TwoLassoSexp::SetHPDQuantiles() {
 if (ConfidenceQuantiles==NULL) {
   Rf_error("SetHPDQuantiles: Sorry, we feed this off of ConfidenceQuantiles.  Cannot Run function.\n");
 }
 if (HPDQuantiles!=NULL) {
   DDelete(HPDQuantiles, "HPDQuantiles");
 }
 int Find50 = 0;
 for (int ii = 0; ii < Rf_length(ConfidenceQuantiles->asSexp()); ii++) {
   if (REAL(ConfidenceQuantiles->asSexp())[ii] == .5) {
     Find50 = 1;  break;
   }
 }
 const double Epsilon = .0000001;
 int CountFs = 0;  //double OnFind;
 for (int ii = 0; ii < Rf_length(ConfidenceQuantiles->asSexp())-1; ii++) { 
   if (REAL(ConfidenceQuantiles->asSexp())[ii] != .50) {
     for (int jj = ii+1; jj < Rf_length(ConfidenceQuantiles->asSexp()); jj++) {
       if (fabs(1.0-REAL(ConfidenceQuantiles->asSexp())[ii] - 
         REAL(ConfidenceQuantiles->asSexp())[jj]) <  Epsilon) {
         CountFs++;   break;
       }
     }
   }
 }
 HPDQuantiles = new AObject(Rf_allocVector(REALSXP, Find50+CountFs));
 int on = 0;
 if (Find50 == 1) {
   REAL(HPDQuantiles->asSexp())[0] = 1.5;  on++;
 }
 for (int ii = 0; ii < Rf_length(ConfidenceQuantiles->asSexp())-1; ii++) { 
   if (REAL(ConfidenceQuantiles->asSexp())[ii] != .50) {
     for (int jj = ii+1; jj < Rf_length(ConfidenceQuantiles->asSexp()); jj++) {
       if (fabs(1.0-REAL(ConfidenceQuantiles->asSexp())[ii] - 
         REAL(ConfidenceQuantiles->asSexp())[jj]) <  Epsilon) {
         REAL(HPDQuantiles->asSexp())[on] = fabs(
           REAL(ConfidenceQuantiles->asSexp())[ii] - 
         REAL(ConfidenceQuantiles->asSexp())[jj]
          ); on++;  break;
       }
     }
   }
 }
 //int One = 1;
 int WantCols = Find50 + 2*CountFs;
 if (Verbose >= 2) {
   Rprintf("--- SetHPDQuantiles() Check to allocate HPDMatrix, WantCols=%d, Find50=%f, CountFs=%d\n",
     WantCols, Find50, CountFs); R_FlushConsole();
 }
 int Len = WantCols*p; 
 if (HPDMatrix == NULL) {
   HPDMatrix = new AObject(Rf_allocMatrix(REALSXP, p, Find50+2*CountFs));
   for (int ii = 0; ii < Len; ii++) {
     REAL(HPDMatrix->asSexp())[ii] = 0.0; }
 } else if ( Rf_isNull(HPDMatrix->asSexp()) ||  !Rf_isReal(HPDMatrix->asSexp()) ||  
   Rf_length(HPDMatrix->asSexp()) != p * (Find50+2*CountFs) ) {
   DDelete(HPDMatrix, "HPDMatrix");
   HPDMatrix = new AObject(Rf_allocMatrix(REALSXP, p, WantCols));    
   for (int ii = 0; ii < Len; ii++) {
   REAL(HPDMatrix->asSexp())[ii] = 0.0; }
 }  else {
   SEXP MyDim = RGET_DIM(HPDMatrix->asSexp());
   if (Rf_isNull(MyDim) || !Rf_isInteger(MyDim) || Rf_length(MyDim) != 2 ||
     INTEGER(MyDim)[1] != p || INTEGER(MyDim)[1] != WantCols) {
     DDelete(HPDMatrix, "HPDMatrix");  
     HPDMatrix = new AObject(Rf_allocMatrix(REALSXP, p, WantCols)); 
     for (int ii = 0; ii < Len; ii++) {
       REAL(HPDMatrix->asSexp())[ii] = 0.0; }       
   }
 }
 if (Verbose >= 2) {
   Rprintf("----  SetHPDQuantiles() Allocate Rownames, colnames for HPDMatrix\n");
 }
  char ALabel[100];
  SEXP sColNames = R_NilValue, sRowNames = R_NilValue, sList = R_NilValue;
  Rf_protect(sColNames = Rf_allocVector(STRSXP, WantCols));
  on = 0;  int on2  = 0;
  if (Find50 == 1) {
    for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
      ALabel[jtj2] = '\0';
    }
    sprintf(ALabel, "Maximum");
    SET_STRING_ELT(sColNames, on2, Rf_mkChar(ALabel));
    on++;  on2++;
  } 
  for (; on < Rf_length(HPDQuantiles->asSexp()); on++) {
    for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
      ALabel[jtj2] = '\0';
    }
    sprintf(ALabel, "HPDLow:%f", REAL(HPDQuantiles->asSexp())[on]);
    SET_STRING_ELT(sColNames, on2, Rf_mkChar(ALabel));
    on2++;
    for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
      ALabel[jtj2] = '\0';
    }
    sprintf(ALabel, "HPDHigh:%f", REAL(HPDQuantiles->asSexp())[on]);
    SET_STRING_ELT(sColNames, on2, Rf_mkChar(ALabel));
    on2++;
  } 
  Rf_protect(sRowNames = Rf_allocVector(STRSXP, p));  
  for (int jtj = 0; jtj < p; jtj++) {
    for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
      ALabel[jtj2] = '\0';
    }
    sprintf(ALabel, "Beta:%d", jtj);
    SET_STRING_ELT(sRowNames, jtj, Rf_mkChar(ALabel));
  } 
  Rf_protect(sList = Rf_allocVector(VECSXP, 2));                 
  if (Verbose >= 3) {
    Rprintf("--  SetHPDQuantiles() : I believe we allocated string list.\n");
  }
  if (Verbose >= 3){
    Rprintf("--  SetHPDQuantiles(): fill List.\n"); R_FlushConsole();
  }
  SET_VECTOR_ELT(sList, 0, sRowNames);
  SET_VECTOR_ELT(sList,1, sColNames); 
  Rf_setAttrib(HPDMatrix->asSexp(), R_DimNamesSymbol, sList); 
  Rf_unprotect(3);    
 return(HPDQuantiles->asSexp());
}
int TwoLassoSexp::GenerateConfidenceMatrix() {
  SEXP sOut = R_NilValue;  int One = 1;
  int jj = 0;  int ii = 0;

  if (Verbose >= 1) {
    Rprintf("--- ATwoLassoObject2011.cc::GenerateConfidenceMatrix(): Starting. \n");
    R_FlushConsole();
  }
  if (ConfidenceQuantiles == NULL) {
    Rprintf("-- GenerateConfidenceMatrix(): Cannot GenerateConfidenceMatrix without Quantiles supplied.\n");
    R_FlushConsole();
    return(-1);
  }
  if (LambdaIndexForConfidenceIntervals < Rf_length(SLambdaAK)-1 && tt1 >
    LambdaIndexForConfidenceIntervals && RSRecOnBeta == NULL) {
    Rprintf("-- GenerateConfidenceMatrix(): Can't work, you gave LambdaIndexForConfidenceIntervals = %d,",
      LambdaIndexForConfidenceIntervals);
    Rprintf("\n  But we are on tt1 = %d, with LambdaAK = %d and gone past!  But RecOnBeta is NULL! \n",
      tt1, Rf_length(SLambdaAK));  
    Rf_error("--GenerateConfidenceMatrix() : Not completing.\n");
  }     
  if (Verbose >= 3) {
      Rprintf("GenerateConfidenceMatrix: Going to RecOnBeta in old occasion.\n",
        LambdaIndexForConfidenceIntervals);
      R_FlushConsole();
  }
  if (Verbose >= 3) {
    Rprintf("--  GenerateConfidenceMatrix(): Allocating matrix\n");
    R_FlushConsole();
  }
  char ALabel[100];
  SEXP sColNames = R_NilValue, sRowNames = R_NilValue, sList = R_NilValue;
  Rf_protect(sColNames = Rf_allocVector(STRSXP, Rf_length(ConfidenceQuantiles->asSexp())));  
  for (int jtj = 0; jtj < Rf_length(ConfidenceQuantiles->asSexp()); jtj++) {
    for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
      ALabel[jtj2] = '\0';
    }
    sprintf(ALabel, "%f", REAL(ConfidenceQuantiles->asSexp())[jtj]);
    SET_STRING_ELT(sColNames, jtj, Rf_mkChar(ALabel));
  } 
  Rf_protect(sRowNames = Rf_allocVector(STRSXP, p));  
  for (int jtj = 0; jtj < p; jtj++) {
    for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
      ALabel[jtj2] = '\0';
    }
    sprintf(ALabel, "Beta:%d", jtj);
    SET_STRING_ELT(sRowNames, jtj, Rf_mkChar(ALabel));
  }       
  Rf_protect(sList = Rf_allocVector(VECSXP, 2));                 
  if (Verbose >= 3) {
    Rprintf("--   GenerateConfidenceMatrix(): I believe we allocated string list.\n");
  }
  if (Verbose >= 3){
    Rprintf("--   GenerateConfidenceMatrix(): fill List.\n"); R_FlushConsole();
  }
  SET_VECTOR_ELT(sList, 0, sRowNames);
  SET_VECTOR_ELT(sList,1, sColNames);
  if (Verbose >= 3){
    Rprintf("--   GenerateConfidenceMatrix(): We filled Row and Colnames try to set names.\n"); R_FlushConsole();
  } 



  if ((ConfidenceMatrix == NULL || Rf_isNull(ConfidenceMatrix->asSexp())) 
    && ConfidenceQuantiles != NULL ) {
    if (Verbose >= 3) {
      Rprintf("--   GenerateConfidenceMatrix(): Allocating ConfidenceMatrix for first time.\n");
      R_FlushConsole();
    }
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, p,
       Rf_length(ConfidenceQuantiles->asSexp())));
    Rf_setAttrib(sOut, R_DimNamesSymbol, sList);
    ConfidenceMatrix = new AObject(sOut);
    Rf_unprotect(1);
  } else if (Rf_length(ConfidenceMatrix->asSexp()) != 
    p * Rf_length(ConfidenceQuantiles->asSexp())) {
    if (Verbose >= 3) {
      Rprintf("--   GenerateConfidenceMatrix(): Bad length on previous confidence matrix =%d not %d.\n",
        Rf_length(ConfidenceMatrix->asSexp()),
        p* Rf_length(ConfidenceQuantiles->asSexp()));
      R_FlushConsole();
    }
    DDelete(ConfidenceMatrix, "ConfidenceMatrix");  
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, p,
       Rf_length(ConfidenceQuantiles->asSexp())));
    Rf_setAttrib(sOut, R_DimNamesSymbol, sList);
    ConfidenceMatrix = new AObject(sOut);
    Rf_unprotect(1);
  }
  Rf_unprotect(3);
  
  Rf_protect(sColNames = Rf_allocVector(STRSXP, Rf_length(ConfidenceQuantiles->asSexp())));  
  for (int jtj = 0; jtj < Rf_length(ConfidenceQuantiles->asSexp()); jtj++) {
    for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
      ALabel[jtj2] = '\0';
    }  
    sprintf(ALabel, "%f", REAL(ConfidenceQuantiles->asSexp())[jtj]);
    SET_STRING_ELT(sColNames, jtj, Rf_mkChar(ALabel));
  } 
  Rf_protect(sRowNames = Rf_allocVector(STRSXP, p));  
  for (int jtj = 0; jtj < p; jtj++) {
    for (int jtj2 = 0; jtj2 < 100-1; jtj2++) {
      ALabel[jtj2] = '\0';
    }  
    sprintf(ALabel, "Beta:%d", jtj+1);
    SET_STRING_ELT(sRowNames, jtj, Rf_mkChar(ALabel));
  }       
  Rf_protect(sList = Rf_allocVector(VECSXP, 2));                 
  if (Verbose >= 3) {
    Rprintf("--   GenerateConfidenceMatrix(): I believe we allocated string list.\n");
  }
  if (Verbose >= 3){
    Rprintf("--   GenerateConfidenceMatrix(): fill List.\n"); R_FlushConsole();
  }
  SET_VECTOR_ELT(sList, 0, sRowNames);
  SET_VECTOR_ELT(sList,1, sColNames);
  if (Verbose >= 3){
    Rprintf("--   GenerateConfidenceMatrix(): We filled Row and Colnames try to set names.\n"); R_FlushConsole();
  } 

  if ((UnshrunkConfidenceMatrix == NULL || Rf_isNull(UnshrunkConfidenceMatrix->asSexp())) 
    && ConfidenceQuantiles != NULL ) {
    if (Verbose >= 3) {
      Rprintf("--    GenerateConfidenceMatrix(): Allocating UnshrunkConfidenceMatrix for first time.\n");
      R_FlushConsole();
    }
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, p,
       Rf_length(ConfidenceQuantiles->asSexp())));
    Rf_setAttrib(sOut, R_DimNamesSymbol, sList);
    UnshrunkConfidenceMatrix = new AObject(sOut);
    Rf_unprotect(1);
  } else if (Rf_length(UnshrunkConfidenceMatrix->asSexp()) != 
    p * Rf_length(ConfidenceQuantiles->asSexp())) {
    if (Verbose >= 3) {
      Rprintf("--    GenerateConfidenceMatrix(): Bad length on previous confidence matrix =%d not %d.\n",
        Rf_length(UnshrunkConfidenceMatrix->asSexp()),
        p* Rf_length(ConfidenceQuantiles->asSexp()));
      R_FlushConsole();
    }
    DDelete(UnshrunkConfidenceMatrix, "ConfidenceMatrix");  
    Rf_protect(sOut = Rf_allocMatrix(REALSXP, p,
       Rf_length(ConfidenceQuantiles->asSexp())));
    Rf_setAttrib(sOut, R_DimNamesSymbol, sList);
    UnshrunkConfidenceMatrix = new AObject(sOut);
    Rf_unprotect(1);
  }
  Rf_unprotect(3);
  if (Verbose >= 3) {
    Rprintf("--    GenerateConfidenceMatrix(): Allocated Unshrunk and attached all of the Lists.\n"); R_FlushConsole();
  }
  
  if (XjXj_CI_Vector != NULL) { Free(XjXj_CI_Vector); XjXj_CI_Vector = NULL;}
  if (C_CI_Vector != NULL) { Free(C_CI_Vector); C_CI_Vector = NULL;}
  if (Sig_CI_Vector != NULL) { Free(Sig_CI_Vector); Sig_CI_Vector = NULL;}
  C_CI_Vector = new AObject(Rf_allocVector(REALSXP,p));
  XjXj_CI_Vector = new AObject(Rf_allocVector(REALSXP,p));
  Sig_CI_Vector = new AObject(Rf_allocVector(REALSXP, p));
  SEXP sHPDQ = SetHPDQuantiles();
  if (Rf_isReal(sHPDQ)) {
  }
  int HPDLen = Rf_length(HPDMatrix->asSexp());
  for (int ii = 0; ii < HPDLen; ii++) {
    REAL(HPDMatrix->asSexp())[ii] = 0.0; 
  }
  
  if (Verbose >= 3) {
    Rprintf("--    GenerateConfidenceMatrix(): Done creating Matrices\n");
    R_FlushConsole();
  } 
  //Rprintf("GenerateConfidenceMatrix: Experiment, Quit NOW!\n");
  // R_FlushConsole();
  // return(1);  
  sOut = ConfidenceMatrix->asSexp();
  
  if (Verbose >= 3) {
    Rprintf("--     GenerateConfidenceMatrix(): Allocating WorkXtResid.\n");
    R_FlushConsole();
  }
  double *WorkXtResid = NULL;
  WorkXtResid = (double *) Calloc(p+2, double);
 
  if (WorkXtResid == NULL) {
    Rf_error("GenerateConfidenceMatrix: Could not reallocate WorkXtResid!\n");
  }
  One = 1;
  // Throw XtResid = X^T(Y-X*Beta) into a working vector:
  F77_CALL(dcopy)(&p, CDO->XTY, &One, WorkXtResid, &One); 
  double ACons = -1.0;
  if  (INTEGER(Stt1)[0] == LambdaIndexForConfidenceIntervals ||
    (INTEGER(Stt1)[0] == Rf_length(SLambdaDK)  && 
    LambdaIndexForConfidenceIntervals >= Rf_length(SLambdaDK) -1)) {
    if (Verbose >= 3) {
      Rprintf("--     GenerateConfidenceMatrix(): Using LambdaIndex = %d, length LambdaAK = %d\n",
      LambdaIndexForConfidenceIntervals, Rf_length(SLambdaDK));
      R_FlushConsole();
    }
    for (jj = 0; jj < p; jj++) { 
      if (CDO->OnBeta[jj] != 0.0) {
        ACons = -1.0 * CDO->OnBeta[jj];
        if (CDO->XTXFlag >= 2) {
          if (CDO->XLC[jj] < 0) {
            if (WorkXtResid != NULL) { Free(WorkXtResid); WorkXtResid = NULL; }
            Rf_error("GenerateConfidenceMatrix: Error, though Beta[%d] = %f, XLC[%d] = %d!\n",
              jj, CDO->OnBeta[jj], jj, CDO->XLC[jj]);
          }
          if (jj < 0 || CDO->XLC[jj]  < 0) {
            Rprintf("** Error ATwoLassoObject2014.cc::GenerateConfidenceMatrix, XLC[jj=%d] = %d, not good! daxpy would have errored.\n",
             jj, CDO->XLC[jj]); R_FlushConsole();
            Rf_error(" GenerateConfidenceMatrix error. \n");
          }
          if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[CDO->XLC[jj]] <= 0) {
            CDO->ReweightCoordinate(jj);
          }
          F77_CALL(daxpy)(&p, &ACons, CDO->pXTX[CDO->XLC[jj]],  &One,
            WorkXtResid, &One);
        } else {
          if (jj < 0) {
            Rprintf("** Error ATwoLassoObject2014.cc::GenerateConfidenceMatrix, jj=%d, XTXFlag = %d, not good!\n",
             jj, CDO->XTXFlag); R_FlushConsole();
            Rf_error(" GenerateConfidenceMatrix error. \n");
          }
          if (CDO->pXTX[jj] == NULL) {
            Rprintf("** Error ATwoLassoObject2014.cc:: pXTX[jj=%d] is NULL!\n", jj);
            Rf_error("ATwoLassoObject2014.cc:: GenerateConfidenceMatrix error. \n");
          }
          if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[jj] <= 0) {
            CDO->ReweightCoordinate(jj);
          }
          F77_CALL(daxpy)(&p, &ACons, CDO->pXTX[jj], &One,
            WorkXtResid, &One);
        }
      }
    }
  } else if (RSRecOnBeta == NULL || Rf_length(RSRecOnBeta->asSexp()) <
     (LambdaIndexForConfidenceIntervals+1)*p) {
    if (WorkXtResid != NULL) { Free(WorkXtResid); WorkXtResid = NULL; }
    Rf_error("Hey: GenreateConfidenceMatrix: You did not record Beta!\n");
  } else {
    if (Verbose >= 3) {
      Rprintf("--     GenerateConfidenceMatrix(): Going to RecOnBeta in old occasion.\n",
      LambdaIndexForConfidenceIntervals);
      R_FlushConsole();
    }
    for (jj = 0; jj < p; jj++) {   
      if (REAL(RSRecOnBeta->asSexp())[LambdaIndexForConfidenceIntervals * p + jj] != 0.0) {  
        ACons = -1.0 * REAL(RSRecOnBeta->asSexp())[LambdaIndexForConfidenceIntervals * p + jj];
        if (CDO->XTXFlag == 2) {
          if (CDO->XLC[jj] < 0) {
            if (WorkXtResid != NULL) { Free(WorkXtResid); WorkXtResid = NULL; }   
            Rf_error("GenerateConfidenceMatrix: Backrecord %d, Error, though Recorded Beta[%d] = %f, XLC[%d] = %d!\n",
              LambdaIndexForConfidenceIntervals, jj, CDO->OnBeta[jj], jj, CDO->XLC[jj]);
          }
          if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[CDO->XLC[jj]] <= 0) {
            CDO->ReweightCoordinate(jj);
          }
          F77_CALL(daxpy)(&p, &ACons, CDO->pXTX[CDO->XLC[jj]], &One, WorkXtResid, &One);
        } else {
          if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[jj] <= 0) {
            CDO->ReweightCoordinate(jj);
          }
          F77_CALL(daxpy)(&p, &ACons, CDO->pXTX[jj], &One, WorkXtResid, &One);
        }
      }
    }
  }
  if (Verbose >= 3) {
    Rprintf("--    GenerateConfidenceMatrix(): About to calculate Sum YSq\n");
    R_FlushConsole();
  }
  //   (Y-XBeta)^T (Y-XBeta) = Y^T %*% Y - Y^TX Beta  - Beta^T X^T(Y-XBeta)
  double SumYSq = 0.0;
  for (ii = 0; ii < n; ii++) {
    SumYSq += CDO->YY[ii]  *CDO->YY[ii]; 
  }
  double BetaTXResid = 0.0;  double XTYBeta = 0.0;
  if (CDO->XTY == NULL) {
    Rprintf("--   GenerateConfidenceMatrix(): CDO XTY does not exist!\n");
    Rf_error("Thats not good!\n");
  }
  if  (INTEGER(Stt1)[0] == LambdaIndexForConfidenceIntervals ||
    (INTEGER(Stt1)[0] == Rf_length(SLambdaDK)  && 
    LambdaIndexForConfidenceIntervals >= Rf_length(SLambdaDK) -1)) {
    BetaTXResid = F77_CALL(ddot)(&p, CDO->OnBeta, &One, WorkXtResid, &One);
    XTYBeta = F77_CALL(ddot)(&p, CDO->OnBeta, &One, CDO->XTY, &One);
  } else if (!(RSRecOnBeta == NULL) &&
    !(Rf_isNull(RSRecOnBeta->asSexp())) &&
    Rf_length(RSRecOnBeta->asSexp()) >= (LambdaIndexForConfidenceIntervals+1)*p) {
    BetaTXResid = F77_CALL(ddot)(&p, 
      REAL(RSRecOnBeta->asSexp()) + 
      LambdaIndexForConfidenceIntervals * p, &One, WorkXtResid, &One);
    XTYBeta = F77_CALL(ddot)(&p,
      REAL(RSRecOnBeta->asSexp()) +  
      LambdaIndexForConfidenceIntervals * p, &One, CDO->XTY, &One);  
  }  else {
    Rprintf("ERROR, GetConfidenceIntervals: UhOh, sad error, LambdaIndex=%d\n",
      LambdaIndexForConfidenceIntervals); R_FlushConsole();
    if (RSRecOnBeta == NULL) {
      Rprintf("--   Yeah, RSRecOnBeta is just NULL!\n"); R_FlushConsole();
    } else if (Rf_isNull(RSRecOnBeta->asSexp())) {
      Rprintf("--   RSRecOnBeta is NULL Sexp holder!\n");  R_FlushConsole();
    } else {
      Rprintf("--   Length RSRecOnBeta = %d \n",
        Rf_length(RSRecOnBeta->asSexp())); R_FlushConsole();
    }
    if (WorkXtResid != NULL) { Free(WorkXtResid); WorkXtResid = NULL; }
    Rf_error("Bad Error, check it out!\n");
  }
  if (SumYSq - BetaTXResid - XTYBeta <= 0.0) {
    Rprintf("--    GenerateConfidenceMatrix() Error, SumYSq = %f, BetaTXResid=%f, XTYBeta = %f, sub = %f\n",
      SumYSq, BetaTXResid, XTYBeta, SumYSq - BetaTXResid - XTYBeta);   R_FlushConsole();
    if (WorkXtResid != NULL) { Free(WorkXtResid); WorkXtResid = NULL; }
    Rf_error("See if you can spot the problem!\n");
  }
  // int_{-infty}^b 1/sqrt(2piA sigma^2)exp (Betaj Xj (Y-XnjBetanj) - Betaj^2 Xj^2)/sigmahat^2
 /*
  if  (INTEGER(Stt1)[0] == LambdaIndexForConfidenceIntervals ||
    INTEGER(Stt1)[0] == Rf_length(SLambdaDK)  && 
    LambdaIndexForConfidenceIntervals >= Rf_length(SLambdaDK) -1) {
    for (jj = 0; jj < p; jj++) {
      if (CDO->XTXFlag == 2) {
        WorkXtResid[jj] += CDO->XTX[CDO->XLC[jj] * p + jj] * CDO->OnBeta[jj];
      }
    }
    BetaTXResid = F77_CALL(ddot)(&p, CDO->OnBeta, &One, WorkXtResid, &One);
    XTYBeta = F77_CALL(ddot)(&p, CDO->OnBeta, &One, CDO->XTY, &One);
  } else {
    BetaTXResid = F77_CALL(ddot)(&p, 
      REAL(RSRecOnBeta->asSexp()) + 
      LambdaIndexForConfidenceIntervals * p, &One, WorkXtResid, &One);
    XTYBeta = F77_CALL(ddot)(&p,
      REAL(RSRecOnBeta->asSexp()) +  
      LambdaIndexForConfidenceIntervals * p, &One, CDO->XTY, &One);  
  }
  */
  double SigOnPoint = 0.0;
  double OnPointBeta = 0.0;
  double A;  double B;  double XjXj;
  
  int kActive = 0;
  double *pBetas;  double iSigOnPoint = 1.0;    double IntegratedEffect = 0.0;
  double onFind;
  if (INTEGER(Stt1)[0] == LambdaIndexForConfidenceIntervals ||
      (INTEGER(Stt1)[0] == Rf_length(SLambdaDK)  && 
      LambdaIndexForConfidenceIntervals >= Rf_length(SLambdaDK) -1)) {
     pBetas = CDO->OnBeta;  
  } else {
    pBetas = REAL(RSRecOnBeta->asSexp()) +  LambdaIndexForConfidenceIntervals *p;
  }
  for (int iii = 0; iii < p; iii++) {
    if (pBetas[iii] != 0.0) { kActive++; }
  }
  
  // Allocate XnXn which is a matrix of the active non-zero factors.
  double *XnXn = (double *) Calloc(kActive *(kActive+1)/2+1, double);
  double *XnXj = (double *) Calloc(kActive, double);
  double *WORK = (double *) Calloc(kActive, double);
  int *IPIV = (int *) Calloc(kActive, double);
  int Info1 = 0;   double CPart = 0.0;  double ZeroD = 0.0;  double OneD = 1.0;
  double PutInA;

  if (Verbose >= 3) {
    Rprintf("--    GenerateConfidenceMatrix(): about to start calculating confidence intervals!\n");
    R_FlushConsole();
  }
  // (Y-X Beta + xjbetaj)^2 = YtY - Beta Xt(Y-XBeta) - YtXBeta + betajxj^T(Y-XBeta) + xj^2betaj^2 
  // 
  int CountOn = 0;
  for (jj = 0; jj < p; jj++) {
    if (fabs(pBetas[jj]) >= .0001) {
      CountOn++;
    }
  }
  double UsePiA = 0.0;
  if (!Rf_isNull(SPiA) && Rf_length(SPiA) >= Rf_length(SLambdaDK) &&
    LambdaIndexForConfidenceIntervals >= Rf_length(SLambdaDK)-1) {
    UsePiA = REAL(SPiA)[Rf_length(SLambdaDK)-1];
  } else if (!Rf_isNull(SPiA) && Rf_length(SPiA) >= LambdaIndexForConfidenceIntervals +1) {
    UsePiA = REAL(SPiA)[LambdaIndexForConfidenceIntervals];
  } else {
    UsePiA = OnPiA;
  }
  double UseLambdaA = 0.0;  double UseLambdaD = 0.0;
  if (!Rf_isNull(SLambdaAK) &&
    LambdaIndexForConfidenceIntervals >= Rf_length(SLambdaAK)-1) {
    UseLambdaA = REAL(SLambdaAK)[Rf_length(SLambdaAK)-1];
  } else if (!Rf_isNull(SLambdaAK) && LambdaIndexForConfidenceIntervals  < Rf_length(SLambdaAK)-1) {
    UseLambdaA = REAL(SLambdaAK)[LambdaIndexForConfidenceIntervals];
  } else {
    UseLambdaA = OnLambdaA;
  }
  

  if (!Rf_isNull(SLambdaDK) &&
    LambdaIndexForConfidenceIntervals >= Rf_length(SLambdaDK)-1) {
    UseLambdaD = REAL(SLambdaDK)[Rf_length(SLambdaDK)-1];
  } else if (!Rf_isNull(SLambdaDK) && LambdaIndexForConfidenceIntervals  < Rf_length(SLambdaDK)-1) {
    UseLambdaD = REAL(SLambdaDK)[LambdaIndexForConfidenceIntervals];
  } else {
    UseLambdaD = OnLambdaD;
  }

  if (Verbose >= 3) {
    Rprintf(" --  GenerateConfidenceMatrix: UseLambdaA=%f, UseLambdaD = %f\n", UseLambdaA, UseLambdaD);
    R_FlushConsole();
  }
  int q = 0;
  int Ooo = 0; int Oo2 = 0;
  int akk = 0;
  
  char UpOrO;
  
  UpOrO = 'L';
  // XjXj is vector of the diagonal of X^T X.
  // XnXn is a sub matrix of X^TX where n is the active (non zero set) 
  //  Of Two Lasso Beta estimate.
  // XnXj is correlation of acitve set with all Xj.
  int TryBad = 0;  
  int IDoneLast = -1; 
  int OnOo2 = 0;

  for (jj = 0; jj < p; jj++) {
    if (Verbose >= 5) {
      Rprintf(" -- On Coord %d, XLC[%d] = %d: ", jj, jj, CDO->XLC[jj]); R_FlushConsole();
      Rprintf("  Note: iWeighted[%d] = %d\n", CDO->iiWeights && CDO->iWeightedXtX != NULL  && CDO->XLC[jj] >= 0 ?
        CDO->iWeightedXtX[CDO->XLC[jj]] : -999); R_FlushConsole();
    }
    TryBad = 0;
    while(TryBad <= 10) {
      OnPointBeta = pBetas[jj];
      if (CDO->XTXFlag >= 2) {
        if (CDO->XLC[jj] >= 0) {
          if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL &&  
            CDO->iWeightedXtX[CDO->XLC[jj]] <= 0) {
            if (Verbose >= 2) {
               Rprintf(" -- Conducting a manual n*jj = %d * %d \n", n, jj); R_FlushConsole();
            }
            XjXj = 0.0; double *Ap = CDO->XX + n * jj;
            for (int ii = 0; ii < n; ii++) {
              XjXj+=Ap[ii]*Ap[ii] * CDO->iiWeights[ii]; 
            }
          } else {
            XjXj = *(CDO->pXTX[CDO->XLC[jj]]+jj);
          }
        } else {
          XjXj = 0.0; double *Ap = CDO->XX + n * jj;
          if (CDO->iiWeights == NULL) {
          for (int ii = 0; ii < n; ii++) {
            XjXj+=Ap[ii]*Ap[ii]; 
          }
          } else {
          for (int ii = 0; ii < n; ii++) {
            XjXj+=Ap[ii]*Ap[ii] * CDO->iiWeights[ii]; 
          }          
          }
        }
        Ooo = 0; Oo2 = 0;
        if (OnPointBeta != 0.0 || IDoneLast <= 0) {
          for (akk = 0;  akk < CDO->OnKappaS; akk++) {
            if (pBetas[CDO->kXFinder[akk]] != 0.0 &&  
              CDO->kXFinder[akk] != jj) {
              if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[akk] <= 0) {
                if (Verbose >= 3) {
                  Rprintf(" -- ConfidenceMatrix: Looks like we will Reweight Coordinate %d/%d  is %d \n", akk, CDO->OnKappaS,
                    CDO->kXFinder[akk]); R_FlushConsole();
                }
                CDO->ReweightCoordinate(CDO->kXFinder[akk]);
              }
              XnXj[Oo2] = *(CDO->pXTX[akk] + jj); Oo2++;
              for (int ajj = akk; ajj < CDO->OnKappaS; ajj++) {
                if (pBetas[CDO->kXFinder[ajj]] != 0.0 && CDO->kXFinder[ajj] != jj) {
                  XnXn[Ooo] = *(CDO->pXTX[akk] + CDO->kXFinder[ajj]); Ooo++;
                }
              }
            } 
          }
          if (TryBad > 0) {
            Ooo = 0;
            for (akk = 0; akk < CDO->OnKappaS; akk++) {
              XnXn[Ooo] = powl(1.05, TryBad) * XnXn[Ooo];
              Ooo+=CDO->OnKappaS-akk;
            }
          }
        } else {
          for (akk = 0;  akk < CDO->OnKappaS; akk++) {
            if (pBetas[CDO->kXFinder[akk]] != 0.0 &&  
              CDO->kXFinder[akk] != jj) {
              if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[akk] <= 0) {
                if (Verbose >= 3) {
                  Rprintf(" -- ConfidenceMatrix: Looks in other Coordinate we will Reweight Coordinate %d/%d  is %d \n", akk, CDO->OnKappaS,
                    CDO->kXFinder[akk]); R_FlushConsole();
                }
                CDO->ReweightCoordinate(CDO->kXFinder[akk]);
              }
              XnXj[Oo2] = *(CDO->pXTX[akk] + jj); Oo2++;
            }
          }         
        }
      } else {
        XjXj = *(CDO->pXTX[jj] + jj);
        Ooo = 0; Oo2 = 0;
        if (OnPointBeta != 0.0 || IDoneLast <= 0) {
          for (akk = 0; akk < p; akk++) {
            if (pBetas[akk] != 0.0 && akk != jj) {
              if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[akk] <= 0) {
                if (Verbose >= 3) {
                  Rprintf(" -- ConfidenceMatrix: Reweight 3 Coordinate %d/%d  is %d \n", akk, CDO->OnKappaS,
                    CDO->kXFinder[akk]); R_FlushConsole();
                }
                CDO->ReweightCoordinate(CDO->kXFinder[akk]);
              }
              XnXj[Oo2] = *(CDO->pXTX[akk] + jj); Oo2++;
              for (int ajj = akk; ajj < p; ajj++) {
                if (pBetas[ajj] != 0.0 &&  ajj != jj) {
                  XnXn[Ooo] = *(CDO->pXTX[akk] + ajj);  Ooo++;
                }
              }
            }
          }
          if (TryBad >= 1) {
            Ooo = 0;
            for (akk = 0; akk < p; akk++) {
              if (pBetas[akk] != 0.0 && akk != jj) {
                  XnXn[Ooo] = powl(1.05, TryBad)*XnXn[Ooo]; 
                  Ooo += Oo2-akk;
              }
            }
          }
        } else {
          for (akk = 0; akk < p; akk++) {
            if (pBetas[akk] != 0.0 && akk != jj) {
              if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[akk] <= 0) {
                if (Verbose >= 3) {
                  Rprintf(" -- ConfidenceMatrix: V 4 will Reweight Coordinate %d/%d  is %d \n", akk, CDO->OnKappaS,
                    CDO->kXFinder[akk]); R_FlushConsole();
                }
                CDO->ReweightCoordinate(CDO->kXFinder[akk]);
              }
              XnXj[Oo2] = *(CDO->pXTX[ akk ]+ jj); Oo2++;
            }
          }
        }
        
      }
      
      //Rprintf("XnXn: We have filled it %d/%d but Oo2 = %d, and we have.\n",
      //  jj, p, Oo2);
      //PrintPackedDouble(XnXn, Oo2, 0, "XnXn");    R_FlushConsole();
      //Rprintf("Beta: "); R_FlushConsole(); 
      //PrintVector(pBetas, p); Rprintf("\n"); R_FlushConsole();
      
      if (OnPointBeta != 0 || IDoneLast <= 0) {
        F77_CALL(dsptrf)(&UpOrO, &Oo2, XnXn, IPIV, &Info1);
      }
      if (Info1 == 0) {
        TryBad = 20;
      } else {
        TryBad++;
      }
    }

      // Invert XaXa

      if (Info1 == 0) {
        if (OnPointBeta != 0 || IDoneLast <= 0) {
          F77_CALL(dsptri)(&UpOrO, &Oo2, XnXn, IPIV, WORK, &Info1);
          if (OnPointBeta == 0.0  && kActive == Oo2) {
            IDoneLast = 1;
          } else {
            IDoneLast = -1;
          }
        }
        //Rprintf("XnXn: We have Inverted it %d/%d and we have Oo2 = %d. kActive %d\n",
        //  jj, p, Oo2, kActive);   R_FlushConsole();
        //PrintPackedDouble(XnXn, Oo2, 0, "XnXn");    R_FlushConsole();
        if (Info1 == 0) {
          //DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
          F77_CALL(dspmv)(&UpOrO, &Oo2, &OneD, XnXn,
            XnXj, &One, &ZeroD, WORK, &One);
          CPart = F77_CALL(ddot)(&Oo2, WORK, &One, XnXj, &One);
        } else {
          Rprintf("--    Failure to invert XnXn after dsptrf, jj = %d\n", jj); R_FlushConsole();
          CPart = 0.0;
        }
      } else {
        Rprintf("--      Failure to dsptrf XnXn j = %d\n", jj); R_FlushConsole();
        CPart = 0.0;
      }
       
      SigOnPoint = SumYSq - BetaTXResid - XTYBeta + OnPointBeta * WorkXtResid[jj] +
        OnPointBeta * OnPointBeta * XjXj;
      if (fabs(pBetas[jj]) >= .0001) {
        SigOnPoint = SigOnPoint / (n - CountOn+1);
      } else {
        SigOnPoint = SigOnPoint / (n - CountOn);
      }
      if (Sig_CI_Vector != NULL) {
        REAL(Sig_CI_Vector->asSexp())[jj] = SigOnPoint;
      }
      if (C_CI_Vector != NULL) {
        REAL(C_CI_Vector->asSexp())[jj] = CPart;
      }
      // Old CPart: CPart = Xa %*% solve(Xa^T Xa) %*% Xa^T
      // CPart = sqrt(SigOnPoint) * CPart;
      
      if (SigOnPoint <= 0.0) {
        Rprintf("--  GenerateConfidence[jj=%d]: really unfortunate Sigma miscalculation \n",jj); R_FlushConsole();
        Rprintf("--     We have that  n=%d, p = %d: SigOnPoint = %f, SumYSq = %f \n",
            n, p, SigOnPoint, SumYSq); R_FlushConsole();
        Rprintf("--     But BetaXTResid = %f, XTYBeta = %f, Beta = %f, OnPointBet = %f\n",
          BetaTXResid, XTYBeta, OnPointBeta);
        Rprintf("--     With WorkXtResid[jj=%d] = %f \n",jj, WorkXtResid);
        Rprintf("--  GenerateConfidence[jj=%d]: No way, SigOnPoint = %f, SumYSq = %f, BetaXTResid = %f, XTYBeta = %f, OnPB=%f, WorkXtR=%f, XjXj = %f\n",
          jj, SigOnPoint, SumYSq, BetaTXResid, XTYBeta, OnPointBeta, WorkXtResid[jj],
          XjXj); R_FlushConsole();
        DidIFail = 1;
        if (WorkXtResid != NULL) { Free(WorkXtResid); WorkXtResid = NULL; }
        if (XnXj != NULL) { Free(XnXj); XnXj = NULL; }
        if (WORK != NULL) { Free(WORK); WORK = NULL; }
        if (IPIV != NULL) { Free(IPIV); IPIV = NULL; }
        if (XnXn != NULL) { Free(XnXn); XnXn = NULL; }
        Rf_error("GenerateConfidence: Error, please diagnose!\n");
      }
      iSigOnPoint = 1.0 / SigOnPoint;
      PutInA = 0.0;   OnOo2 = 0;
      for (akk = 0; akk < CDO->OnKappaS; akk++) {
        if (pBetas[CDO->kXFinder[akk]] != 0.0 && CDO->kXFinder[akk] != jj) {
          PutInA += WorkXtResid[CDO->kXFinder[akk]] * WORK[OnOo2]; OnOo2++;
        }
      }
      if (OnPointBeta != 0.0) {
        PutInA += OnPointBeta * F77_CALL(ddot)(&Oo2, WORK, &One, XnXj, &One);
      }
      A = iSigOnPoint * (WorkXtResid[jj] + XjXj * OnPointBeta - 
        PutInA);
        
      // Correction for OnPoint Beta;
      //if (OnPointBeta != 0.0) {
      //  A += OnPointBeta * CPart * iSigOnPoint; 
      //}
      if (Verbose >= 4) {
        Rprintf("--     NewA: j=%d/%d: A=%f. OnPointBeta=%f, CPart=%f, PutInA = %f, XjXj=%f, iSigOnPoint=%f, WorkXtResid[%d]=%f, Oo2=%d, XnXj = ",
          jj, p, A, OnPointBeta, CPart, PutInA, XjXj, iSigOnPoint, jj, WorkXtResid[jj], Oo2); 
        PrintVector(XnXj, Oo2);
        Rprintf("\n");
        R_FlushConsole();
      }
      //
      //  A should be 1/sigma^2 * (  Xj^T(Y-X_a hatBeta_a) - 
      //   (Y-X_ahatBeta_a)^T %*% X_a solve(X_a^T X_a) X_a^T X_j
      
      //  (Y-XahatBeta_a - XjBetaj + XjBetaj ) %*% 
      //
      //  Betaj^T X_j^T X_a solve(X_a^T X_a) X_a^T X_j
  
      B = XjXj * iSigOnPoint;
      if (XjXj_CI_Vector != NULL && Rf_length(XjXj_CI_Vector->asSexp()) > jj) {
        REAL(XjXj_CI_Vector->asSexp())[jj] = XjXj;
      }
      // Previous code, why is CPart squared?  I don't know, I don't think it matters!
      //if (XjXj >= CPart*CPart) {
      //  B = B - CPart * CPart * iSigOnPoint;
      // } 
      // New Version:
      if (XjXj >= CPart) {
        B = B - CPart * iSigOnPoint;
      }     
      if (Verbose >= 4) {
        Rprintf("--      GenerateConfidence: we have jj=%d, OnBeta=%f, SigOnPoint = %f, A=%f, B=%f\n",
          jj, pBetas[jj], SigOnPoint, A, B); R_FlushConsole();
      }
      IntegratedEffect =  AllInt(A, B, UsePiA, UseLambdaA, UseLambdaD);
      if (Verbose >= 6) {
        Rprintf("--       got IntegratedEffect = \n",
          IntegratedEffect); R_FlushConsole();
      }
      onFind = 0.0;
      for (q = 0; q < Rf_length(ConfidenceQuantiles->asSexp()); q++) {
         onFind= SeekQuantile(REAL(ConfidenceQuantiles->asSexp())[q], 
           onFind, .0001,  A, B, UsePiA, UseLambdaA, UseLambdaD, 
           IntegratedEffect);
         if (Verbose >= 6) {
           Rprintf("--       Found Quantile %d or %f  at position %f \n", q, 
             REAL(ConfidenceQuantiles->asSexp())[q], onFind);
           R_FlushConsole();
         }
         REAL(ConfidenceMatrix->asSexp())[q * p + jj]  = onFind; 
      }
      // Now without Shrinkage!
      //  Shrinkage penalty devalues confidence near zero this undoes this.
      //
      IntegratedEffect =  AllInt(A, B, UsePiA, UseLambdaA, UseLambdaA);
      if (Verbose >= 6) {
        Rprintf("--     got IntegratedEffect = \n",
          IntegratedEffect); R_FlushConsole();
      }
      onFind = 0.0;
      for (q = 0; q < Rf_length(ConfidenceQuantiles->asSexp()); q++) {
         onFind= SeekQuantile(REAL(ConfidenceQuantiles->asSexp())[q], 
           onFind, .0001,  A, B, UsePiA, UseLambdaA, UseLambdaA, 
           IntegratedEffect);
         if (Verbose >= 6) {
           Rprintf("-- Found Quantile %d or %f  at position %f \n", q,
             REAL(ConfidenceQuantiles->asSexp())[q], onFind);
           R_FlushConsole();
         }
         REAL(UnshrunkConfidenceMatrix->asSexp())[q * p + jj]  = onFind; 
      }
      if (HPDQuantiles != NULL && HPDMatrix != NULL) {
        if (REAL(HPDQuantiles->asSexp())[0] == 1.5) {
          SeekHPD(REAL(HPDQuantiles->asSexp())+1,  
            Rf_length(HPDQuantiles->asSexp()) -1,  .000001,
             A,  B, UsePiA,  UseLambdaA, UseLambdaD, 0.0,
            2*p,
            REAL(HPDMatrix->asSexp())+p+jj, 
            REAL(HPDMatrix->asSexp())+2*p+jj, 
             REAL(HPDMatrix->asSexp())+jj);
        } else {
          SeekHPD(REAL(HPDQuantiles->asSexp()),  
            Rf_length(HPDQuantiles->asSexp()),  .000001,
            A,  B, UsePiA,  UseLambdaA, UseLambdaD, 0.0,
            2*p,
            REAL(HPDMatrix->asSexp())+jj, 
            REAL(HPDMatrix->asSexp())+p+jj, 
            NULL);        
        }
      }
      
    }
  if (Verbose >= 4) {
    Rprintf("GenerateConfidence: All done! Freeing \n"); R_FlushConsole();
  }
  if (WorkXtResid != NULL) {
    //Rprintf("Free WorkXtResid\n"); R_FlushConsole();
    Free(WorkXtResid);
    WorkXtResid = NULL;
  }
  if (XnXn != NULL) {
    //Rprintf("Free XnXn \n"); R_FlushConsole();
    Free(XnXn); XnXn = NULL;
  }
  //Rprintf("XnXj: Free \n"); R_FlushConsole();
  if (XnXj != NULL) { Free(XnXj); XnXj = NULL; }
  //Rprintf("WORK: Free \n"); R_FlushConsole();
  if (WORK != NULL) {Free(WORK); WORK = NULL; }
  //Rprintf("Free IPIV \n"); R_FlushConsole();
  if (IPIV != NULL) { Free(IPIV); IPIV = NULL; }
  //Rprintf("Well we should return now. \n"); R_FlushConsole();
  return(1);
}

double TwoLassoSexp::TestAllInt(double A, double B, double PiA, 
  double LambdaA, double LambdaD) {
  return(AllInt(A,B,PiA, LambdaA, LambdaD));
}

double TwoLassoSexp::TestSuperInt(double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons, double b) {
  return(SuperInt(A,B,PiA, LambdaA, LambdaD, NegCons, b));
}

// Sd = 1/ sqrt(B);
// -(beta - MeanBeta)^2/(2Sd^2) = a beta - b beta^2/2 + Cons
//  MeanBeta / Sd^2 = a;  MeanBeta  = a / b
//  Cons = MeanBeta^2 / (2 Sd^2) = .5 * B * MeanBeta^2;
// int_-infty^infty
// piA*LambdaA/sqrt(2*pi * Sig)*exp( A beta - B beta^2/2 - LambdaA|beta|) dbeta +
// (1-piA)*LambdaD/sqrt(2*pi * Sig)*exp( A beta - B beta^2 - LambdaD|beta|) dbeta
double AllInt(double A, double B, double PiA, double LambdaA, double LambdaD) {
  double PartA1, PartA2, PartB1, PartB2;
  if (B <= 0.0) {
    // Note B cannot be zero if A is nonzero.  Thus we assume A is also zero and integrate easily.
    // int .5 exp(-LambdaA|beta|) * piA * LambdaA + .5 (1-piA) * LambdaD* exp(-LambdaD|beta|) = 
    // return(log( PiA / LambdaA + (1.0-PiA) / LambdaD));
    return(log( PiA + (1.0-PiA) ));
  }
  double Sd = 1/sqrt(B);
  double MyCons = -log(2.0)+M_LN_SQRT_2PI - .5 * log(B);
  PartA1 =  (A + LambdaA) / B;
  PartA1 = Rf_pnorm5(0, PartA1, Sd,  1, 1) +  PartA1 * PartA1 * B * .5;
  PartA2 =  (A - LambdaA) / B;
  PartA2 = Rf_pnorm5(0, PartA2, Sd,  0, 1) +  PartA2 * PartA2 * B * .5;
  PartB1 =  (A + LambdaD) / B;
  PartB1 = Rf_pnorm5(0, PartB1, Sd,  1, 1) +  PartB1 * PartB1 * B * .5;
  PartB2 =  (A - LambdaD) / B;
  PartB2 = Rf_pnorm5(0, PartB2, Sd,  0, 1) +  PartB2 * PartB2 * B * .5; 
  if (PartA1 >= PartB1 + 8.0 && PartA1 >= PartB2 + 8.0) {
    if (PartA1 >= PartA2 + 8.0) {
      return(MyCons + PartA1 + log( PiA * LambdaA * (1.0 + exp(PartA2-PartA1)) + 
        (1.0-PiA) * LambdaD * (exp(PartB1-PartA1) + exp(PartB2-PartA1))));
    } else if (PartA2 >= PartA1 + 8.0) {
      return(MyCons +PartA2 + log( PiA * LambdaA* (1.0 + exp(PartA1-PartA2)) + 
        (1.0-PiA) * LambdaD*(exp(PartB1-PartA2) + exp(PartB2-PartA2))));    
    } else {
      return(MyCons +PartA1 + log( PiA * LambdaA* (1.0 + exp(PartA2-PartA1)) + 
        (1.0-PiA) * LambdaD * (exp(PartB1-PartA1) + exp(PartB2-PartA1))));       
    }
  } else if (PartA2 >= PartB1 + 8.0 && PartA2 >= PartB2 + 8.0) {
    return(MyCons +PartA2 + log( PiA * LambdaA* (1.0 + exp(PartA1-PartA2)) + 
        (1.0-PiA) * LambdaD * (exp(PartB1-PartA2) + exp(PartB2-PartA2)))); 
  } else if (PartB1 >= PartA1 + 8.0 && PartB1 >= PartA2 + 8.0) {
    if (PartB1 >= PartB2 + 8.0) {
      return(MyCons +PartB1 + log( (1.0-PiA) * LambdaD * (1.0 + exp(PartB2-PartB1)) + 
        PiA * LambdaA * (exp(PartA1-PartB1) + exp(PartA2-PartB1))));
    } else if (PartB2 >= PartB1 + 8.0) {
      return(-log(2.0) +PartB2 + log( (1.0-PiA) * LambdaD * (1.0 + exp(PartB1-PartB2)) + 
        (PiA) * LambdaA *  (exp(PartA1-PartB2) + exp(PartA2-PartB2))));    
    } else {
      return(MyCons +PartB1 + log( (1.0-PiA) *LambdaD* (1.0 + exp(PartB2-PartB1)) + 
        PiA * LambdaA * (exp(PartA1-PartB1) + exp(PartA2-PartB1))));      
    }  
  } else if (PartB2 >= PartA1 + 8.0 && PartB2 >= PartA2 + 8.0) {
    return(MyCons +PartB2 + log( (1.0-PiA) * LambdaD * (1.0 + exp(PartB1-PartB2)) + 
        (PiA) * LambdaA * (exp(PartA1-PartB2) + exp(PartA2-PartB2))));    
  }
  return(MyCons +log( PiA * LambdaA* (exp(PartA1) + exp(PartA2)) + (1.0-PiA) * LambdaD * (exp(PartB1) + exp(PartB2)))); 
  // pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

}

//   int_-infty^B
// LambdaA*piA/sqrt(2*pi * Sig)*exp( A beta - B beta^2 - LambdaA|beta|) dbeta +
// LambdaD*(1-piA)/sqrt(2*pi * Sig)*exp( A beta - B beta^2 - LambdaD|beta|) dbeta
double SuperInt(double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons, double b) {
  double PartA1, PartA2, PartB1, PartB2;
  double Sd = 1/sqrt(B);
  double MyCons = -log(2.0)+M_LN_SQRT_2PI - .5 * log(B);
  if (B <= 0.0) {
    // Note B cannot be zero (never less than zero) if A is nonzero.  Thus we assume A is also zero and integrate easily.
    // int exp(-LambdaA|beta|) * piA + (1-piA) * exp(-LambdaD|beta|) = 
    if (b <= 0.0) {
      return( PiA * exp(-LambdaA *fabs(b)-NegCons+ MyCons) +
        (1.0-PiA) * exp(-LambdaD *fabs(b)-NegCons+MyCons));
    }

    PartA1 = (PiA/LambdaA +(1.0-PiA)/LambdaD) * exp(-NegCons+MyCons);
    PartA2 = (PiA * (1.0-exp(-LambdaA * b)) + 
      (1.0-PiA) * (1.0-exp(-LambdaD*b))) * exp(-NegCons+MyCons);
    if (PiA == 1.0 || PiA == 0.0) {
      return(PartA1 + PartA2);
    } 
    return(PartA1 + PartA2);
  }  
  PartA2 =  (A - LambdaA) / B;
  if (b <= 0.0) {
    PartA1 =  (A + LambdaA) / B;
    PartA1 = Rf_pnorm5(b, PartA1, Sd,  1, 1) +  PartA1 * PartA1 * B * .5-NegCons+MyCons;
    PartB1 =  (A + LambdaD) / B;
    PartB1 = Rf_pnorm5(b, PartB1, Sd,  1, 1) +  PartB1 * PartB1 * B * .5-NegCons+MyCons;
    if (PiA == 1.0) {
      return(exp(PartA1));
    } else if (PiA == 0.0) {
      return(exp(PartB1));
    }
    if (log(PiA) - log(1.0-PiA) + log(LambdaA) + PartA1 >= log(LambdaD)+5+ PartB1) {
      return( exp(PartA1) * PiA * LambdaA * ( 
        1.0 + exp(PartB1-PartA1 + log(1.0-PiA) - log(PiA) - log(LambdaA) + log(LambdaD) )));
    } else if (log(PiA) - log(1.0-PiA) + log(LambdaA) + PartA1 <= log(LambdaD)-5+ PartB1) {
      return( exp(PartB1) * (1.0-PiA) * LambdaD * ( 
        1.0 + exp(PartA1-PartB1 - log(1.0-PiA) + log(PiA) + log(LambdaA) - log(LambdaD) )));    
    }
    return(PiA * LambdaA * exp(PartA1) + (1.0-PiA) * LambdaD* exp(PartB1));
  }
  PartA1 =  (A + LambdaA) / B;
  PartA1 = Rf_pnorm5(0, PartA1, Sd,  1, 1) +  PartA1 * PartA1 * B * .5-NegCons+MyCons;
  PartB1 =  (A + LambdaD) / B;
  PartB1 = Rf_pnorm5(0, PartB1, Sd,  1, 1) +  PartB1 * PartB1 * B * .5-NegCons+MyCons;
  PartA2 =  (A - LambdaA) / B;
  double PzA0 = Rf_pnorm5(0, PartA2,Sd,0,1);
  double PzAb = Rf_pnorm5(b, PartA2,Sd,0,1);  
  PartA2 = PzA0 + log(1.0-exp(PzAb-PzA0)) + PartA2 * PartA2 * B * .5 - NegCons+MyCons;
  
  PartB2 =  (A - LambdaD) / B;
  double PzB0 = Rf_pnorm5(0, PartB2,Sd,0,1);
  double PzBb = Rf_pnorm5(b, PartB2,Sd,0,1);  
  PartB2 = PzB0 + log(1.0-exp(PzBb-PzB0)) + PartB2 * PartB2 * B * .5 - NegCons+MyCons;
  
  return(PiA * LambdaA * (exp(PartA1) + exp(PartA2)) + LambdaD * (1.0-PiA) * 
    (exp(PartB1) + exp(PartB2)) );
}
///////////////////////////////////////////////////////////////////////
// SeekQuantile
//
//  Starting at StartAt, within a tolerance of Epislon, this seeks a quantile
//  For Semi-Gaussian multimodal density:  exp(-Ax^2+2Bx)( LambdaA PiA exp(-LambdaA |x|) + LambdaD(1-PiA) exp{-LambdaD|x|} )
//
double SeekQuantile(double QuantileSeek, double StartAt, double Epsilon,
  double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons) {
  double PastQuantile = -999.0;
  double OnQuantile = 0.0;
  double AoA, BoB;
  if (!R_finite(StartAt)) {
    StartAt = -10.0; PastQuantile = 0.0;  OnQuantile = -10.0;
    AoA = SuperInt(A, B, PiA, LambdaA, LambdaD, NegCons, OnQuantile);
    BoB = SuperInt(A, B, PiA, LambdaA, LambdaD, NegCons, PastQuantile); 
  } else {
    PastQuantile = StartAt + 2.0;  OnQuantile = StartAt;
    AoA = SuperInt(A, B, PiA, LambdaA, LambdaD, NegCons, OnQuantile);
    BoB = SuperInt(A, B, PiA, LambdaA, LambdaD, NegCons, PastQuantile); 
  }
  if (fabs(AoA - QuantileSeek) <= Epsilon) {
    return(OnQuantile);
  } else if (fabs(BoB - QuantileSeek) <= Epsilon) {
    return(PastQuantile);
  }
  //Rprintf("SeekQuantile: Start, we want %f, have A at %f:%f, and B at %f:%f\n",
  //  QuantileSeek, OnQuantile, AoA, PastQuantile, BoB); R_FlushConsole();
  int countBads = 0;
  double Movediffs = 1.0;
  if (Movediffs == 0.0) {
  }
  double NewPoint, NewMove;
  while(fabs(AoA - QuantileSeek) > Epsilon && countBads <= 400) {
    //Rprintf("SeekQuantile: On countBads=%d, we want %f, have A at %f:%f, and B at %f:%f\n",
    //  countBads, QuantileSeek, OnQuantile, AoA, PastQuantile, BoB); R_FlushConsole();
    if (BoB < QuantileSeek) {
      Movediffs = 1.0;
      OnQuantile = PastQuantile;
      AoA = BoB;
      PastQuantile += 1.0;
      BoB = SuperInt(A, B, PiA, LambdaA, LambdaD, NegCons, PastQuantile); 
    } else if (AoA > QuantileSeek) {
      Movediffs = 1.0;
      PastQuantile = OnQuantile;
      BoB = AoA;
      OnQuantile -= 1.0;
      AoA = SuperInt(A, B, PiA, LambdaA, LambdaD, NegCons, OnQuantile);     
    } else {
      NewPoint = (OnQuantile + PastQuantile) * .5;
      NewMove =  SuperInt(A, B, PiA, LambdaA, LambdaD, NegCons, NewPoint);
      if (fabs(NewMove-QuantileSeek) <= Epsilon) {
        return(NewPoint);
      } else if (NewMove > QuantileSeek) {
        BoB = NewMove;  PastQuantile = NewPoint;
      } else {
        AoA = NewMove;  OnQuantile = NewPoint;
      } 
    }
    countBads++;
  }
  if (fabs(AoA - QuantileSeek) <= Epsilon) {
    return(OnQuantile);
  }
  // Rprintf("We tried to seek, but we moved countBads = %d times and didn't find it!\n",
  //   countBads);
  return(-999);
}

double GetVal(double AtX, double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons) {
  if (AtX == 0.0) {
     return( log(PiA * LambdaA * exp(-NegCons) +
        (1.0-PiA) * LambdaD * exp(-NegCons) ) - log(2.0) );
  }
  double MyA = A * AtX - .5 * B * AtX * AtX -NegCons;
  return(MyA - LambdaA * fabs(AtX)-log(2.0)+
    log( PiA * LambdaA +
      (1-PiA) * LambdaD * exp( - (LambdaD-LambdaA) * fabs(AtX)
        ) ));
}

SEXP TwoLassoSexp::TestGetVal(SEXP sAtX, SEXP sA, SEXP sB, SEXP sPiA, 
  SEXP sLambdaA, SEXP sLambdaD, SEXP sNegCons) {
  SEXP sOuts = R_NilValue;
  Rf_protect(sOuts = Rf_allocVector(REALSXP, Rf_length(sAtX)));
  double AtX = REAL(sAtX)[0];
  if (AtX <= 0.0) {
  }
  double A = REAL(sA)[0];  double B = REAL(sB)[0];  double PiA = REAL(sPiA)[0];
  double LambdaA = REAL(sLambdaA)[0];  double LambdaD = REAL(sLambdaD)[0];
  double NegCons = REAL(sNegCons)[0];
  for (int ii = 0; ii < Rf_length(sAtX); ii++) {
    REAL(sOuts)[ii] = GetVal(REAL(sAtX)[ii], A, B, PiA, LambdaA, LambdaD, NegCons);
  }
  Rf_unprotect(1); return(sOuts);
}

double SeekLevel(int GoLeft1GoRight0, double LevelSeek, double UpBound, double DownBound, double Epsilon,
  double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons) {
  double MidBound; double vMidBound; 
  const int TooManyMax = 200;
  int Count = 0;
  if (UpBound <= DownBound) {
    MidBound = UpBound;
    UpBound = DownBound;
    DownBound = MidBound;
  }
  double vUpBound = GetVal(UpBound, A, B, PiA, LambdaA, LambdaD, NegCons);
  double vDownBound = GetVal(DownBound, A, B, PiA, LambdaA, LambdaD, NegCons); 
  double OtherBound = 0.0;
  if (fabs(vUpBound-LevelSeek) <= Epsilon) {
    return(UpBound);
  }
  if (fabs(vDownBound-LevelSeek) <= Epsilon) {
    return(DownBound);
  }
  if (GoLeft1GoRight0 == 1) {
    if (vUpBound+Epsilon < LevelSeek) {
      OtherBound =  GetVal(UpBound, A, B, PiA, LambdaA, LambdaD, NegCons);
      Rprintf("**************************************************************\n");
      Rprintf("*****: Error  trying to get OtherBound \n");
      Rprintf(" Error: Seek Level: hey, vUpBound = %f < LevelSeek = %f on supply left! Epsilon=%f\n", vUpBound, LevelSeek, Epsilon);
      Rprintf("      : The vUpBound = %f +Epsilon = %f  does not satisfy less than LevelSeek = %f\n", vUpBound, Epsilon, LevelSeek);
      Rprintf("  A=%f; B=%f; PiA = %f; LambdaA = %f; LambdaD = %f; NegCons=%f; UpBound = %f\n",
        A, B, PiA, LambdaA, LambdaD, NegCons, UpBound); R_FlushConsole();
      Rprintf("Error: SeekLevel, NegCons=%f, UpBound=%f, vUpBound=%f can't find LevelSeek=%f, OtherBound = %f\n",
        NegCons, UpBound, vUpBound, LevelSeek, OtherBound);
      return(R_NaReal);
    }
    if (DownBound >= 0.0) {
      vMidBound = GetVal(0.0, A,B, PiA, LambdaA, LambdaD, NegCons);
      if (vMidBound > LevelSeek) {
        DownBound = -UpBound;
        vDownBound = GetVal(DownBound, A, B, PiA, LambdaA, LambdaD, NegCons);
        UpBound =  0.0;
        vUpBound = vMidBound;
      }
    }
    if (vDownBound > LevelSeek) {
      if (UpBound == DownBound) {
        DownBound = UpBound-1;
        vDownBound = GetVal(DownBound, A, B, PiA, LambdaA, LambdaD, NegCons); 
      }
      while(vDownBound > LevelSeek && Count < TooManyMax) {
        DownBound = DownBound-(UpBound-DownBound);
        vDownBound = GetVal(DownBound, A, B, PiA, LambdaA, LambdaD, NegCons); 
        Count ++;
      }
    }  
    MidBound = (UpBound+DownBound)*.5;
    vMidBound = GetVal(MidBound, A, B, PiA, LambdaA, LambdaD, NegCons);  
    Count=0;
    while (fabs(vMidBound-LevelSeek) > Epsilon && Count < TooManyMax)  {
      if (fabs(vMidBound-LevelSeek) <= Epsilon) {
        return(MidBound);
      } else if (vMidBound > LevelSeek) {
        UpBound = MidBound;
        vUpBound = vMidBound;
      } else {
        DownBound = MidBound;
        vDownBound = vMidBound;
      }
      Count++; 
      MidBound = (UpBound+DownBound)*.5;
      vMidBound = GetVal(MidBound, A, B, PiA, LambdaA, LambdaD, NegCons);  
    }
    if (Count >= TooManyMax) {
      Rprintf("** Error: we didn't get close for LevelSeek = %f go left.\n", 
        LevelSeek);   R_FlushConsole();
      return(R_NaReal);
    }
  }
  
  if (GoLeft1GoRight0 == 0) {
    if (vDownBound+Epsilon < LevelSeek) {
      Rprintf("************************************************************************************\n");
      Rprintf("**  Error on Go right, vDownBound = %f (at %f), vUpBound = %f (at %f), Epsilon=%f, LevelSeek = %f\n",
        vDownBound, DownBound, vUpBound, UpBound, Epsilon, LevelSeek);
      Rprintf("** Well this is bad and we should not continue!!\n "); R_FlushConsole();
      Rprintf("** Error: SeekLevel, NegCons=%f, DownBound=%f, vDownBound=%f can't find LevelSeek=%f\n",
        NegCons, DownBound, vDownBound, LevelSeek);  R_FlushConsole();
      Rprintf("A=%f; B=%f; PiA=%f; LambdaA=%f; LambdaD=%f; NegCons=%f\n",
        A, B, PiA, LambdaA, LambdaD, NegCons); R_FlushConsole();
      Rprintf("Well lets see why you did this \n");
      return(R_NaReal);
    }
    if (UpBound < 0.0) {
      vMidBound = GetVal(0.0, A, B, PiA, LambdaA, LambdaD, NegCons);
      if (vMidBound > LevelSeek) {
        UpBound = -DownBound;
        vUpBound = GetVal(UpBound, A, B, PiA, LambdaA, LambdaD, NegCons);
        DownBound = 0.0;
        vDownBound = vMidBound;
      }
      //Rprintf("We did the wacky, vMidBound=%f at 0.0, go up.  We have DownBound=0.0, UpBound=%f, vUpBound=%f, seek %f\n",
      //  vMidBound, UpBound, vUpBound, LevelSeek); R_FlushConsole();
    }
    if (UpBound < DownBound) {
      Rprintf("Error: shouldn't be here, UpBound=%f, DownBound = %f\n",
        UpBound, DownBound); R_FlushConsole();
      return(R_NaReal);
    }
    if (vUpBound > LevelSeek) {
      if (UpBound == DownBound) {
        UpBound = DownBound+1;
        vUpBound = GetVal(UpBound, A, B, PiA, LambdaA, LambdaD, NegCons); 
      }
      while(vUpBound > LevelSeek && Count < TooManyMax) {
        UpBound = UpBound+(UpBound-DownBound);
        vUpBound = GetVal(UpBound, A, B, PiA, LambdaA, LambdaD, NegCons); 
        Count++;
      }
    } 
    vDownBound = GetVal(DownBound, A, B, PiA, LambdaA, LambdaD, NegCons);
    MidBound = (UpBound+DownBound)*.5;
    vMidBound = GetVal(MidBound, A, B, PiA, LambdaA, LambdaD, NegCons);  
    vUpBound = GetVal(UpBound, A, B, PiA, LambdaA, LambdaD, NegCons);
    if (vUpBound > LevelSeek) {
      Rprintf("hey, we did a lot of crazyness at an end, vUpBound = %f, LevelSeek = %f\n",
        vUpBound, LevelSeek);
      Rprintf("UpBound = %f?\n", UpBound);
      Rprintf("  and vDownBound = %f, and DownBound = %f, vMidBound = Rf, MidBound = %f\n",
        vDownBound, DownBound, vMidBound, MidBound); R_FlushConsole();
    }
    if (vDownBound < LevelSeek) {
      Rprintf("hey, we did a lot of crazyness at an end, vDownBound = %f, DownBound = %f, vUpBound = %f, LevelSeek = %f\n",
        vDownBound, DownBound, vUpBound, LevelSeek);
      Rprintf("UpBound = %f?\n", UpBound);
      Rprintf("  and vDownBound = %f, and DownBound = %f, vMidBound = Rf, MidBound = %f\n",
        vDownBound, DownBound, vMidBound, MidBound); R_FlushConsole();
    }
    Count = 0;
    while (fabs(vMidBound-LevelSeek) >= Epsilon && Count < TooManyMax)  {   
      if (fabs(vMidBound-LevelSeek) < Epsilon) {
        return(MidBound);
      } else if (vMidBound > LevelSeek) {
        DownBound = MidBound;
        vDownBound = vMidBound;
      } else {
        UpBound = MidBound;
        vUpBound = vMidBound;
      }
      Count++;
      MidBound = (UpBound+DownBound)*.5;
      vMidBound = GetVal(MidBound, A, B, PiA, LambdaA, LambdaD, NegCons); 
    }  
    if (Count >= TooManyMax) {
      Rprintf("** Count Max Out, we have Count = %d/%d but vMidBound = %f (at %f), vUpBound = %f (at %f), vDownBound = %f (at %f). \n",
        Count, TooManyMax, vMidBound, MidBound, vUpBound, UpBound, vDownBound, DownBound); R_FlushConsole();
      Rprintf("** Error: we didn't get close for LevelSeek = %f go right.\n", 
        LevelSeek);
      return(R_NaReal);
    }
    MidBound = (UpBound+DownBound)*.5;
    vMidBound = GetVal(MidBound, A, B, PiA, LambdaA, LambdaD, NegCons);  
  }
  return(.5 * (UpBound+DownBound));
}

double TwoLassoSexp::TestSeekLevel(SEXP sGoLeft1GoRight0, SEXP sLevelSeek,
    SEXP sUpBound, SEXP sDownBound, SEXP sEpsilon,
    SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD, SEXP sNegCons) {
  int GoLeft1GoRight0 = 0;
  if (Rf_isReal(sGoLeft1GoRight0)) {
    GoLeft1GoRight0 = REAL(sGoLeft1GoRight0)[0];
  } else {
    GoLeft1GoRight0 = INTEGER(sGoLeft1GoRight0)[0];
  }
  return(SeekLevel(GoLeft1GoRight0, REAL(sLevelSeek)[0], REAL(sUpBound)[0], 
    REAL(sDownBound)[0], REAL(sEpsilon)[0], REAL(sA)[0], REAL(sB)[0],
    REAL(sPiA)[0], REAL(sLambdaA)[0], REAL(sLambdaD)[0], REAL(sNegCons)[0])); 
}
    
double LDet(double xx, double A, double B, double PiA, double LambdaA, double LambdaD) {
  double ZCons = (1.0-PiA) / PiA * LambdaD / LambdaA;
  if (xx >= 0.0) {
  return(A - B * xx  - LambdaD + ( (LambdaD-LambdaA)  ) /
  (1.0 + ZCons * exp(-(LambdaD-LambdaA)*xx) ));
  }
  return(A - B * xx  + LambdaD - ( (LambdaD-LambdaA) * 1.0 ) /
  (1.0 + ZCons * exp((LambdaD-LambdaA)*xx)));
} 
SEXP TwoLassoSexp::TestLDet(SEXP sxx, SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD) {
 SEXP sOut = R_NilValue;
 Rf_protect(sOut = Rf_allocVector(REALSXP, Rf_length(sxx)));
 for (int ii = 0; ii < Rf_length(sxx); ii++) {
   REAL(sOut)[ii] = LDet(REAL(sxx)[ii], REAL(sA)[0], REAL(sB)[0], 
     REAL(sPiA)[0], REAL(sLambdaA)[0], REAL(sLambdaD)[0]);
 }
 Rf_unprotect(1);
 return(sOut);
}
  
double SecDet(double xx, double A, double B, double PiA, double LambdaA, double LambdaD) {
  double Lss;
  double ZCons = (1.0-PiA) * LambdaD / (PiA * LambdaA);
  double MySign = 1.0;
  if (MySign == 0.0) {
    Rprintf("SecDet, really?\n");
  }
  if (xx >= 0.0) {
    Lss = (LambdaD-LambdaA) * xx;
  } else {
    Lss = - (LambdaD-LambdaA) * xx;
    MySign = -1.0;
  }
  double Go = 1.0 + ZCons * exp(-Lss);
  return(-B + (LambdaD-LambdaA)*(LambdaD-LambdaA) * exp(-Lss) * ZCons / (Go*Go));
} 

SEXP TwoLassoSexp::TestSecDet(SEXP sxx, SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD) {
 SEXP sOut = R_NilValue;
 Rf_protect(sOut = Rf_allocVector(REALSXP, Rf_length(sxx)));
 for (int ii = 0; ii < Rf_length(sxx); ii++) {
   REAL(sOut)[ii] = SecDet(REAL(sxx)[ii], REAL(sA)[0], REAL(sB)[0], 
     REAL(sPiA)[0], REAL(sLambdaA)[0], REAL(sLambdaD)[0]);
 }
 Rf_unprotect(1);
 return(sOut);
}

double MaxPosterior(double A, double B, double PiA, double LambdaA, double LambdaD, 
  double Epsilon, double NegCons ) {
  double At0 = log(PiA * LambdaA + (1-PiA) * LambdaD)-NegCons;
  if (At0 < -1.0) {
  }
  double MaxOfA = 0.0;
  double AtMaxOfA;
  if (A > 0) {
    MaxOfA = (A - LambdaA) / B;
    if (MaxOfA < 0) { MaxOfA = 0.0;}

  }
  if (A < 0) {
    MaxOfA = (A + LambdaA) / B;
    if (MaxOfA > 0) { MaxOfA = 0.0;}
  }
  AtMaxOfA = GetVal(MaxOfA, A, B, PiA, LambdaA, LambdaD, NegCons);

  int tt;  const int Maxtt = 100;
  double dDet =  LDet(AtMaxOfA, A, B, PiA, LambdaA, LambdaD);
  double sDet = 0.0;
  double OnMax = MaxOfA;
  if (OnMax <= -1.0) {
  }
  double NewMax;  double ValAtNewMax = 0.0;
  if (ValAtNewMax >= 1.0) {
  }
  double OlddDet = dDet;
  if (fabs(dDet) <= Epsilon || OlddDet <= Epsilon) {
    return(OnMax);
  }
  for (tt = 0; tt < Maxtt; tt++) {
    sDet =  SecDet(OnMax,  A,  B,  PiA,  LambdaA,  LambdaD);
    OlddDet = dDet;
    NewMax = OnMax - dDet / sDet;
    dDet =  LDet(NewMax, A, B, PiA, LambdaA, LambdaD);
    //Rprintf("tt=%d/%d: we moved, sDet=%f, OlddDet = %f, dDet=%f, OnMax=%f, ValAtNewMax = %f, NewMax=%f\n ",
    //    tt, Maxtt, sDet, OlddDet, dDet, OnMax, ValAtNewMax, NewMax, tt); R_FlushConsole();
    if (fabs(dDet) <= Epsilon) {
      ValAtNewMax = GetVal(NewMax, A,B, PiA, LambdaA, LambdaD, NegCons);
      //Rprintf("ValAtNewMax = %f, NewMax=%f, steps is tt=%d/%d \n",
      //  ValAtNewMax, NewMax, tt, Maxtt); R_FlushConsole();
      return(NewMax);
    }
    OnMax = NewMax;
  }
  //Rprintf("Well, the derviatives didn't move, MaxOfA = %f, AtMaxOfA = %f\n",
  //  MaxOfA, AtMaxOfA);
  return(MaxOfA);
}
double TwoLassoSexp::TestMaxPosterior(SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD, 
  SEXP sEpsilon, SEXP sNegCons )  {
  return( MaxPosterior(REAL(sA)[0], REAL(sB)[0], REAL(sPiA)[0], 
    REAL(sLambdaA)[0], REAL(sLambdaD)[0], 
    REAL(sEpsilon)[0], REAL(sNegCons)[0] ));
}

///////////////////////////////////////////////////////////////////////////////
//  SeekHPD
//    Seek Highest Posterior Density Interval
//
//  Given A, B, PiA, LambdaA, LambdaD, One searches for quantiles of 
//   exp(- Ax^2 + 2 Bx) * { LambdaA PiA exp(-LambdaA|x|) + LambfdaD(1-PiA) exp(-LambdaD|x|)}
//   regions.
// 
//
//
int SeekHPD(double *QuantilesList, int LenQuantile,  double Epsilon,
  double A, double B, double PiA, double LambdaA, double LambdaD,  double NegCons,
  int SkipGo,
  double *pOutLeft, double *pOutRight, double *PutMax) {
  double FoundMax   = MaxPosterior( A, B, PiA, LambdaA, LambdaD, Epsilon, 0.0);
  double TopMax = GetVal(FoundMax, A, B, PiA, LambdaA, LambdaD, 0.0);
  
  double OtherMax = GetVal(0.0, A, B, PiA, LambdaA, LambdaD, 0.0);
  if (OtherMax > TopMax) {
    TopMax = OtherMax;  FoundMax = 0.0;
  }
  // TopMax is the maximum location, which could be at zero or at MaxPosterior's Newton Raphson solution.
  
  // Integrate entire density to know quantiles.  IntTo is log f(x)
  double IntTo = AllInt( A,  B,  PiA, LambdaA, LambdaD);  
  NegCons = IntTo;
  if (SkipGo <= 0) { SkipGo = 1;}
  double AboveLevel = TopMax-IntTo;
  //Rprintf("SeekHPD: seek "); PrintVector(QuantilesList, LenQuantile); 
  //Rprintf("\n"); R_FlushConsole();
  //Rprintf("SeekHPD: We get FoundMax = %f, TopMax = %f\n", FoundMax, TopMax);
  //R_FlushConsole();

  //Rprintf("SeekHPD: We have A=%f, B=%f, PiA=%f, IntTo = %f \n",
  //  A, B, PiA, IntTo);  R_FlushConsole();
  int OnGo = 0;

  double UpLeft = FoundMax, UpRight = FoundMax;
  double DownLeft = FoundMax, DownRight = FoundMax;
  if (PutMax != NULL) {
    PutMax[0] = FoundMax;
  }
  // R_NaReal  R_isnancpp 
  DownLeft = SeekLevel(1, AboveLevel-1.0, UpLeft, UpLeft-1.0, Epsilon,
    A, B, PiA,  LambdaA,  LambdaD,  NegCons);
  if (R_isnancpp(DownLeft)) {
    Rprintf("**************************************************************\n");
    Rprintf("** Error on SeekLevel to make DownLeft                        \n");
    Rprintf("##\n   GoLeftRight=%d, AboveLevel=%f, UpLeft=%f, Epsilon=%f, A=%f, B=%f",
      1, AboveLevel, UpLeft, Epsilon, A, B);
    Rprintf("\n     PiA=%f, LambdaA=%f, LambdaD = %f, NegCons=%f \n",
      PiA, LambdaA, LambdaD, NegCons); R_FlushConsole(); 
    Rprintf("\n*** Return Error \n");  R_FlushConsole();
    Rprintf("\n* TopMax = %f, OtherMax = %f, FoundMax = %f\n", TopMax, OtherMax, FoundMax);
    Rf_error("Bad Seek Level!\n");
  }
  DownRight = SeekLevel(0, AboveLevel-1.0, UpRight+1.0, UpRight, Epsilon,
    A, B, PiA,  LambdaA,  LambdaD,  NegCons); 
  if (R_isnancpp(DownRight)) {
    Rprintf("**************************************************************\n");
    Rprintf("** Error on SeekLevel to make DownRight                        \n");
    Rprintf("##\n   GoLeftRight=%d, AboveLevel=%f, UpRight=%f, Epsilon=%f, A=%f, B=%f",
      0, AboveLevel, UpRight, Epsilon, A, B);
    Rprintf("\n     PiA=%f, LambdaA=%f, LambdaD = %f, NegCons=%f \n",
      PiA, LambdaA, LambdaD, NegCons); R_FlushConsole(); 
    Rprintf("\n*** Return Error \n");  R_FlushConsole();

    Rf_error("Bad Seek Level!\n");
  }

   
  double vUpLeft = TopMax - NegCons;
  double vUpRight = TopMax - NegCons;   
  double vDownLeft = GetVal(DownLeft, A, B, PiA, LambdaA, LambdaD, NegCons); 
  double vDownRight = GetVal(DownRight, A, B, PiA, LambdaA, LambdaD, NegCons); 
  int Count = 0;
  const int TooManyMax = 200;
  
  if (vUpLeft < 0) {
  } else if (vUpRight < 0) {
  }

  //Rprintf("SeekHPD: OnGo=%d: NegCons=%f, got to -1.0 from Ups(%f,%f), levels of %f, %f\n",
  //  OnGo, NegCons, UpLeft, UpRight, 
  //  vUpLeft, vUpRight);
  //R_FlushConsole();  
  double NewLeft, NewRight;
  double vNewLeft, vNewRight;
  double UpLevel = AboveLevel; double UpInt = 0.0;

  double DownLevel=0.0, DownInt=0.0;
  DownInt = 0.0;  
  DownLevel = AboveLevel - 1.0;
  double MidLevel=0.0, MidInt = 0.0;

  if (fabs(vDownLeft - DownLevel) > Epsilon) {
     NewLeft = SeekLevel(1,  DownLevel,  UpLeft, DownLeft,  Epsilon,
       A,  B,  PiA,  LambdaA,  LambdaD,  NegCons);
    if (R_isnancpp(NewLeft)) {
      Rprintf("**************************************************************\n");
      Rprintf("** Error on SeekLevel to make NewLeft                        \n");
      Rprintf("##\n   GoLeftRight=%d, DownLevel=%f, UpLeft=%f, DownLeft = %f, Epsilon=%f, A=%f, B=%f",
        0, DownLevel, UpLeft, DownLeft, Epsilon, A, B);
      Rprintf("\n     PiA=%f, LambdaA=%f, LambdaD = %f, NegCons=%f \n",
        PiA, LambdaA, LambdaD, NegCons); R_FlushConsole(); 
      Rprintf("\n*** Return Error \n");  R_FlushConsole();
      Rf_error("Bad Seek Level!\n");
    }
     vNewLeft =  GetVal(DownLeft, A, B, PiA, LambdaA, LambdaD, NegCons); 
     if (fabs(vNewLeft - DownLevel) > 2*Epsilon) {
       Rprintf("**************************************************************\n");
       Rprintf("** Hey: Seek level %f failed on go Left from DownLeft = %f,  FoundMax=%f, A=%f, B=%f, NegCons=%f, TopMax=%f\n",
         DownLevel, DownLeft, FoundMax, A, B, NegCons, TopMax);
       Rprintf("** LambdaA=%f, LambdaD=%f, Epsilon=%f, PiA=%f, IntTo=%f\n",
         LambdaA, LambdaD, Epsilon, PiA, IntTo);R_FlushConsole();
       Rprintf("We sought DownLevel = %f, but got vNewRight = %f\n", DownLevel, vNewRight);
       Rprintf("FoundMax=%f, and supposedly its level = AboveLevel = %f\n", FoundMax, AboveLevel);
       Rf_error("** Well Search don't work. \n");
     }
     DownLeft = NewLeft;  vDownLeft = vNewLeft;
  } else {
  } 
  if (fabs(vDownRight - DownLevel) > Epsilon) {
     NewRight = SeekLevel(0,  DownLevel,  DownRight, UpRight,  Epsilon,
       A,  B,  PiA,  LambdaA,  LambdaD,  NegCons);
     if (R_isnancpp(NewRight)) {
       Rprintf("**************************************************************\n");
       Rprintf("** Error on SeekLevel to make NewRight                        \n");
       Rprintf("##\n   GoLeftRight=%d, DownLevel=%f, DownRight=%f, UpRight = %f, Epsilon=%f, A=%f, B=%f",
         0, DownLevel, DownRight, UpRight, Epsilon, A, B);
       Rprintf("\n     PiA=%f, LambdaA=%f, LambdaD = %f, NegCons=%f \n",
         PiA, LambdaA, LambdaD, NegCons); R_FlushConsole(); 
       Rprintf("\n*** Return Error \n");  R_FlushConsole();
       Rf_error("Bad Seek Level!\n");
     }
     vNewRight =  GetVal(NewRight, A, B, PiA, LambdaA, LambdaD, NegCons); 
     if (fabs(vNewRight - DownLevel) > 2*Epsilon) {
       Rprintf("**************************************************************\n");
       Rprintf("** Hey: Seek level %f failed on go Right from DownRight = %f,  FoundMax=%f, A=%f, B=%f, NegCons=%f, TopLevel=%f\n",
         DownLevel, DownRight, FoundMax, A, B, NegCons, TopMax);
       Rprintf("We sought DownLevel = %f, but got vNewRight = %f\n", DownLevel, vNewRight);
       Rprintf("FoundMax=%f;  AboveLevel = %f; ## supposedly its level\n", FoundMax, AboveLevel);
       Rprintf("** LambdaA=%f; LambdaD=%f; Epsilon=%f; PiA=%f; IntTo = %f;\n",
         LambdaA, LambdaD, Epsilon, PiA, IntTo); R_FlushConsole();
       Rf_error("** Well Search don't work. \n");
     }
     DownRight = NewRight;  vDownRight = vNewRight;
  } else {
  }
  
  //DownInt = (SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, DownRight) - 
  //  SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, DownLeft)); 
  /*Rprintf("***\n***\n");
  Rprintf("*** - Okay before beginning we have DownLeft=%f, UpLeft=%f, UpRight=%f, DownRight = %f for DownLevel = %f\n",
    DownLeft, UpLeft, UpRight, DownRight, DownLevel);
  Rprintf("*** - DownInt = %f, first seek is %f, Wonder what the start of the algorithm brings\n",
    DownInt, QuantilesList[0]);
  Rprintf("***\n"); R_FlushConsole();
  */
 for (OnGo = 0; OnGo < LenQuantile; OnGo++) {
  DownInt = SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, DownRight) - 
    SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, DownLeft); 
  DownLevel = GetVal(DownLeft, A, B, PiA, LambdaA, LambdaD, NegCons); 
  //Rprintf("Well new OnGo=%d/%d seek %f, starting with DownInt=%f, left-right = (%f,%f)\n",
  //  OnGo, LenQuantile, QuantilesList[OnGo], DownInt, DownLeft, DownRight); R_FlushConsole();
  Count = 0;
  while (DownInt < QuantilesList[OnGo]) {
    UpLevel = DownLevel;  UpInt = DownInt;
    vUpRight = DownLevel;  vUpLeft = DownLevel;
    UpRight = DownRight;  UpLeft =  DownLeft;
    DownLevel -= 1.0;
     NewRight = SeekLevel(0,  DownLevel,  DownRight+1.0, UpRight,  .5*Epsilon,
       A,  B,  PiA,  LambdaA,  LambdaD,  NegCons);
     if (R_isnancpp(NewRight)) {
      Rprintf("**************************************************************\n");
      Rprintf("** Error on SeekLevel to make NewRight, OnGo=%d, looking for Quant %d\n", OnGo, QuantilesList[OnGo]);
      Rprintf("##\n   0=%d, DownLevel=%f, DownRight+1.0=%f, UpRight = %f, Epsilon=%f, A=%f, B=%f",
        0, DownLevel, DownRight+1.0, UpRight, Epsilon, A, B);
      Rprintf("\n     PiA=%f, LambdaA=%f, LambdaD = %f, NegCons=%f \n",
        PiA, LambdaA, LambdaD, NegCons); R_FlushConsole(); 
      Rprintf("\n*** Return Error \n");  R_FlushConsole();
      Rf_error("Bad Seek Level!\n");
     }
     vNewRight =  GetVal(NewRight, A, B, PiA, LambdaA, LambdaD, NegCons); 
     NewLeft = SeekLevel(1,  DownLevel,  UpLeft, DownLeft-1.0,  .5*Epsilon,
       A,  B,  PiA,  LambdaA,  LambdaD,  NegCons);
      if (R_isnancpp(NewLeft)) {
      Rprintf("**************************************************************\n");
      Rprintf("** Error on SeekLevel to make NewLeft, OnGo=%d, looking for Quant %d\n", OnGo, QuantilesList[OnGo]);
      Rprintf("##\n   0=%d, DownLevel=%f, UpLeft=%f, DownLeft-1.0 = %f, Epsilon=%f, A=%f, B=%f",
        0, DownLevel, UpLeft, DownLeft-1.0, Epsilon, A, B);
      Rprintf("\n     PiA=%f, LambdaA=%f, LambdaD = %f, NegCons=%f \n",
        PiA, LambdaA, LambdaD, NegCons); R_FlushConsole(); 
      Rprintf("\n*** Return Error \n");  R_FlushConsole();
      }
    
     vNewLeft =  GetVal(NewLeft, A, B, PiA, LambdaA, LambdaD, NegCons); 
     if (fabs(vNewRight - DownLevel) > 2*Epsilon) {
       Rprintf("**************************************************************\n");
       Rprintf("** Hey: While Seek level %f failed on go Left from DownRight = %f,  FoundMax=%f, A=%f, B=%f, NegCons=%f\n",
         DownLevel, DownRight, FoundMax, A, B, NegCons);
       Rprintf("** LambdaA=%f, LambdaD=%f, Epsilon=%f, PiA=%f, DownInt=%f, IntTo = %f\n",
         LambdaA, LambdaD, Epsilon, PiA, DownInt, IntTo); R_FlushConsole();
       Rprintf("** Turns out vNewRight = %f, DownLevel=%f, diff = %f for NewRight = %f\n",
         vNewRight, DownLevel, vNewRight - DownLevel, NewRight);
       Rf_error("** Well Search don't work. \n");
      } 
      if (fabs(vNewLeft - DownLevel) > 2*Epsilon) {
       Rprintf("**************************************************************\n");
       Rprintf("** Hey: While Seek level %f failed on go Left from DownLeft = %f,  FoundMax=%f, A=%f, B=%f, NegCons=%f\n",
         DownLevel, DownLeft, FoundMax, A, B, NegCons);
       Rprintf("** LambdaA=%f; LambdaD=%f; Epsilon=%f; PiA=%f; DownInt = %f\n",
         LambdaA, LambdaD, Epsilon, PiA, DownInt);R_FlushConsole();
       Rprintf("** Turns out vNewRight = %f, DownLevel=%f, diff = %f for NewRight = %f\n",
         vNewLeft, DownLevel, vNewLeft - DownLevel, NewLeft);
       Rf_error("** Well Search don't work. \n");
      } 
     DownLeft = NewLeft;  DownRight = NewRight;
     vDownLeft = vNewLeft;  vDownRight = vNewRight; 
     DownInt = (SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, DownRight) - 
       SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, DownLeft)); 
      Count++;
     if (Count >= TooManyMax) {
       Rprintf("****************************************************************\n");
       Rprintf("** Error on moving Down below.  We have Count = %d >= %d toomany!\n",
         Count, TooManyMax); R_FlushConsole();
       Rprintf("## LengthQuantilesList = %d\n", LenQuantile);
       Rprintf("Epsilon=%f;  A=%f; B=%f; PiA = %f; LambdaA = %f; LambdaD = %f; NegCons=%f\n",
         Epsilon, A, B, PiA, LambdaA, LambdaD, NegCons);
       Rprintf("SkipGo = %d\n", SkipGo);
       Rprintf("QuantilesList = "); PrintVector(QuantilesList, LenQuantile);
       Rprintf("NewRight = %f; vNewRight = %f; NewLeft = %f;  vNewLeft = %f\n", 
         NewRight, vNewRight, NewLeft, vNewLeft);
       Rprintf("DownLeft = %f, vDownLeft=%f, DownRight = %f, vDownRight = %f\n", DownLeft, vDownLeft, DownRight, vDownRight);
       Rprintf("  UpLeft=%f, UpRight = %f, If we redo UpInt we get %f \n", UpLeft, UpRight,
         (SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, UpRight) - 
          SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, UpLeft)) 
       );
       Rprintf("  DownLeft=%f, DownRight = %f, If we redo DownInt we get %f \n", DownLeft, DownRight,
         (SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, DownRight) - 
          SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, DownLeft)) 
       );
       R_FlushConsole();
       R_FlushConsole(); Rprintf("\n");
       Rprintf("OutLeft = "); PrintVector(pOutLeft, LenQuantile);
       Rprintf("\n"); R_FlushConsole();
       Rprintf("OutRight = "); PrintVector(pOutRight, LenQuantile);   
       Rprintf("\n"); R_FlushConsole(); 
       Rprintf("PutMax = "); PrintVector(PutMax, LenQuantile);        
       Rprintf("\n");  R_FlushConsole();
       Rprintf("** OnGo = %d/%d, Epsilon = %f, MidInt=%f, DownInt=%f, UpInt=%f, DownLevel = %f, UpLevel = %f\n",
          OnGo, LenQuantile, Epsilon, MidInt, DownInt, UpInt, DownLevel, UpLevel);
       Rf_error("Toomany loop.\n");
     }
   }
   /*Rprintf("****: We pass on DownInt=%f (level %f) at (%f,%f), UpInt=%f (level=%f) at (%f,%f) seek QuantilesList[OnGo=%d] = %f\n",
     DownInt, DownLevel, DownLeft, DownRight, UpInt, UpLevel, 
       UpLeft, UpRight, OnGo, QuantilesList[OnGo]); R_FlushConsole();
   */
   MidLevel = (DownLevel + UpLevel)  * .5;
   if (DownInt == UpInt) {
     Rf_error("Hey: This can't be, UpInt=DownInt in the middle!\n");
   }
   //Rprintf("***: We will seek MidLevel = %f, with  for OnGo = %d \n",
   //  MidLevel, OnGo); R_FlushConsole();
   Count = 0;
   if (fabs(UpInt - QuantilesList[OnGo]) <= Epsilon) {
     Rprintf("This is weird, OnGo=%d, at beginning, UpInt = %f, QuantilesList[%d]=%f, UpLeft = %f, UpRight = %f\n",
       OnGo, UpInt, OnGo, QuantilesList[OnGo], UpLeft, UpRight);
     Rprintf("  and DownInt = %f, DownLevel = %f, UpLevel=%f, Down = %f,%f, what's up? (no real error though) \n",
       DownInt, DownLevel, UpLevel, DownLeft, DownRight); R_FlushConsole();
     pOutLeft[OnGo*SkipGo] = UpLeft;  pOutRight[OnGo*SkipGo] = UpRight;
   }
   if (DownLevel+Epsilon >= UpLevel) {
     if (DownInt >= QuantilesList[OnGo] && UpInt <= QuantilesList[OnGo]) {
       pOutLeft[OnGo * SkipGo] = DownLeft;  pOutRight[OnGo * SkipGo] = DownRight;
     } else {
       Rprintf("Hey, this ain't good, before the loop, DownLevel=%f (int %f), UpLevel=%f (int %f)\n",
         DownLevel, DownInt, UpLevel, UpInt);
       Rprintf("MidInt = %f, DownInt=%f, UpInt=%f, QW is %\n",
         MidInt, DownInt, UpInt, QuantilesList[OnGo]);
       Rprintf("Our QuatilesWant is %f\n", QuantilesList[OnGo]);   R_FlushConsole();
       Rprintf("  At (%f,%f) and (%f,%f) but if Epsilon=%f we have violation!\n",
         DownLeft, DownRight, UpLeft, UpRight, Epsilon);
       Rf_error("Seek this!\n"); 
     }
   }
   do  {
     MidLevel = .5 * (UpLevel + DownLevel);
     NewRight = SeekLevel(0,  MidLevel,  DownRight, UpRight,  Epsilon,
       A,  B,  PiA,  LambdaA,  LambdaD,  NegCons);
      if (R_isnancpp(NewRight)) {
      Rprintf("**************************************************************\n");
      Rprintf("** Error on SeekLevel in do to make NewRight, OnGo=%d, looking for Quant %d, MidLevel=%f\n", OnGo, QuantilesList[OnGo], MidLevel);
      Rprintf("##\n   0=%d, DownLevel=%f, DownRight=%f, UpRight = %f, Epsilon=%f, A=%f, B=%f",
        0, DownLevel, DownRight, UpRight, Epsilon, A, B);
      Rprintf("\n     PiA=%f, LambdaA=%f, LambdaD = %f, NegCons=%f \n",
        PiA, LambdaA, LambdaD, NegCons); R_FlushConsole(); 
      Rprintf("\n*** Return Error \n");  R_FlushConsole();
      }
      
     vNewRight =  GetVal(NewRight, A, B, PiA, LambdaA, LambdaD, NegCons); 
     NewLeft = SeekLevel(1,  MidLevel,  UpLeft, DownLeft,  Epsilon,
       A,  B,  PiA,  LambdaA,  LambdaD,  NegCons);
     if (R_isnancpp(NewLeft)) {
      Rprintf("**************************************************************\n");
      Rprintf("** Error on SeekLevel in do to make NewLeft, OnGo=%d, looking for Quant %d, MidLevel=%f\n", OnGo, QuantilesList[OnGo], MidLevel);
      Rprintf("##\n   0=%d, MidLevel=%f, UpLeft=%f, DownLeft = %f, Epsilon=%f, A=%f, B=%f",
        1, MidLevel, UpLeft, DownLeft, Epsilon, A, B);
      Rprintf("\n     PiA=%f, LambdaA=%f, LambdaD = %f, NegCons=%f \n",
        PiA, LambdaA, LambdaD, NegCons); R_FlushConsole(); 
      Rprintf("\n*** Return Error \n");  R_FlushConsole();
     }
     vNewLeft =  GetVal(NewLeft, A, B, PiA, LambdaA, LambdaD, NegCons); 
     if (fabs(vNewRight - MidLevel) > 10*Epsilon) {
       Rprintf("**************************************************************\n");
       Rprintf("** Hey: While Seek Midlevel %f failed on go Right from DownRight = %f and UpRight=%f,  FoundMax=%f, A=%f, B=%f, NegCons=%f\n",
         MidLevel, DownRight, UpRight, FoundMax, A, B, NegCons);
       Rprintf("** LambdaA=%f, LambdaD=%f, Epsilon=%f, PiA=%f, DownInt=%f\n",
         LambdaA, LambdaD, Epsilon, PiA, DownInt); R_FlushConsole();
       Rf_error("** Well Search don't work. \n");
      } 
      if (fabs(vNewLeft - MidLevel) > 10*Epsilon) {
       Rprintf("**************************************************************\n");
       Rprintf("** Hey: While Seek level %f failed on go Left from UpLeft = %f to DownLeft=%f,  FoundMax=%f, A=%f, B=%f, NegCons=%f\n",
         MidLevel, UpLeft, DownLeft, FoundMax, A, B, NegCons);
       Rprintf("** LambdaA=%f, LambdaD=%f, Epsilon=%f, PiA=%f, DownInt = %f\n",
         LambdaA, LambdaD, Epsilon, PiA, DownInt);R_FlushConsole();
       Rf_error("** Well Search don't work. \n");
      } 
      if (NewRight < NewLeft) {
        Rf_error("Error at this point because NewRight = %f < NewLeft = %f we'll fail!\n",
          NewRight, NewLeft);
      } 
      MidInt = (SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, NewRight) - 
       SuperInt(A, B, PiA, LambdaA, LambdaD, IntTo, NewLeft)); 
      if (MidInt < 0) {
        Rf_error("Error, we integrated from %f to %f but got MidInt = %f!\n",
          NewLeft, NewRight, MidInt); 
      }
      /*Rprintf("We get MidLevel=%f at (%f,%f) between %f at (%f,%f) and %f at (%f,%f) has Int = %f but Seek %f\n",
        MidLevel, NewLeft, NewRight, DownLevel, DownLeft, DownRight, UpLevel, UpLeft, UpRight, MidInt,
        QuantilesList[OnGo]);
      R_FlushConsole();*/
      if (MidInt < UpInt) {
        Rprintf("*************************************************************\n");
        Rprintf("** Error detection code  MidInt < UpInt\n"); R_FlushConsole();
        Rprintf("**  We have DownLevel = %f (%f,%f) Int is %f, UpLevel = %f (%f,%f), UpInt = %f\n",
          DownLevel, DownLeft, DownRight, DownInt, UpLevel, UpLeft, UpRight, UpInt); R_FlushConsole();
        Rf_error("hey, we are in the split level, but MidInt %f < UpInt=%f!\n", MidInt, UpInt);
      }
      if (fabs(MidInt-QuantilesList[OnGo]) <= Epsilon) {
        pOutLeft[OnGo * SkipGo] = NewLeft;  pOutRight[OnGo * SkipGo] = NewRight;
        UpInt = MidInt;  UpLevel = MidLevel;
        UpLeft = NewLeft;  UpRight = NewRight; 
        break;
      } else if (MidInt > QuantilesList[OnGo]) {
        DownInt = MidInt;  DownLevel = MidLevel;
        DownLeft = NewLeft;  DownRight = NewRight;
      } else {
        UpInt = MidInt;  UpLevel = MidLevel;
        UpLeft = NewLeft;  UpRight = NewRight;
      }  
      if (DownLevel + Epsilon >= UpLevel) {
        if (DownInt > MidInt && UpInt < MidInt) {
          pOutLeft[OnGo * SkipGo] = DownLeft;  pOutRight[OnGo * SkipGo] = DownRight;
          break;
        }
      } else if (UpLevel == DownLevel) {
        Rf_error("Error: OnGo = %d, we're here but UpLevel = %f = DownLevel = %f but we shouldn't!\n",
          OnGo, UpLevel, DownLevel);
      }  
      if (UpInt == DownInt) {
        Rf_error("Error: OnGo = %d, we're here but UpInt = %f = DownInt = %f, but MidInt = %f, but we shouldn't,\n",
          OnGo, UpInt, DownInt, MidInt);
      }
      Count++;
      if (Count >= TooManyMax) {
       Rprintf("** Error on moving Down to meet Up.  We have Count = %d >= %d toomany!\n",
         Count, TooManyMax); R_FlushConsole();
        Rprintf("** OnGo = %d/%d, Epsilon = %f, MidInt=%f, DownInt=%f, UpInt=%f, DownLevel = %f, UpLevel = %f\n",
          OnGo, Epsilon, MidInt, DownInt, UpInt, DownLevel, UpLevel);
        Rf_error("Error Go Home\n");
     }
   } while (fabs(UpInt- QuantilesList[OnGo]) > Epsilon  && DownLevel+Epsilon < UpLevel);  
   if (fabs(MidInt-QuantilesList[OnGo]) < Epsilon) { 
     pOutLeft[OnGo * SkipGo] = NewLeft;  pOutRight[OnGo * SkipGo] = NewRight;
   } else   if (UpLevel < DownLevel+Epsilon) {
     if (DownInt >= QuantilesList[OnGo] && UpInt <= QuantilesList[OnGo]) {
       pOutLeft[OnGo*SkipGo] = DownLeft;  pOutRight[OnGo*SkipGo] = DownRight;
     } else {
       Rf_error("Hey: UpLevel=%f, DownLevel = %f, UpInt=%f, DownInt=%f, bad result!\n",
         UpLevel, DownLevel, UpInt, DownInt);
     }
   } else if (MidInt > QuantilesList[OnGo]) {
     pOutLeft[OnGo*SkipGo] = DownLeft;  pOutRight[OnGo*SkipGo] = NewRight;
   }
   if (pOutLeft[OnGo*SkipGo] == 0.0 && pOutRight[OnGo*SkipGo] == 0.0  &&
     (fabs(NewLeft) > Epsilon || fabs(NewRight) > Epsilon)) {
    Rprintf("*************************************************************\n");
    Rprintf("** Hey OnGo=%d/%d, we have double zeros. \n", OnGo, LenQuantile);
    Rprintf("**  DownLevel = %f, int %f, (%f,%f) \n",
      DownLevel, DownInt, DownLeft, DownRight);
    Rprintf("**  UpLevel = %f, int %f, (%f,%f) \n",
      UpLevel, UpInt, UpLeft, UpRight);      
    Rprintf("**  MidLevel = %f, MidInt %f, (%f,%f)",
      MidLevel, MidInt, NewLeft, NewRight);
    Rprintf("** We were seeking Quant %f which is %d/%d\n", 
      QuantilesList[OnGo], OnGo, LenQuantile);
    Rprintf("**  Epsilon = %f, What did we do? \n", Epsilon); R_FlushConsole();
    Rprintf("** \n"); R_FlushConsole(); 
    if (UpLevel < DownLevel + Epsilon) {
      Rprintf("But we do have UpLevel < DownLevel+Epsilon!, (Up-Down)/Epsilon=%f\n",
        (UpLevel-DownLevel)/Epsilon);  R_FlushConsole();
    } else {
      Rprintf("But we do not have UpLevel < DownLevel+Epsilon!, (Up-Down)/Epsilon=%f\n",
        (UpLevel-DownLevel)/Epsilon);   R_FlushConsole(); 
    }
    } else {
      //Rprintf("Well OnGo=%d/%d we go give (%f,%f)\n",
      //  OnGo, LenQuantile, pOutLeft[OnGo], pOutRight[OnGo]); R_FlushConsole();
    }
  }
 
  //Rprintf("Okay we're at the end of SeekHPD algorithm\n"); R_FlushConsole();
  return(1);
}

SEXP TwoLassoSexp::TestSeekHPD(SEXP sQuantilesList, SEXP sEpsilon,
  SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD,  SEXP sNegCons) {
  SEXP Outs = R_NilValue;
  Rf_protect(Outs = Rf_allocMatrix(REALSXP, Rf_length(sQuantilesList), 2));
  SEXP sMax = R_NilValue;
  Rf_protect(sMax = Rf_allocVector(REALSXP, 1));
  SEXP sList = R_NilValue;
  Rf_protect(sList = Rf_allocVector(VECSXP, 2));
  for (int ii = 0; ii < Rf_length(Outs); ii++) {
    REAL(Outs)[ii] = 0.0;
  }
  int Go;
  Go = SeekHPD(REAL(sQuantilesList), Rf_length(sQuantilesList),  
    REAL(sEpsilon)[0],
    REAL(sA)[0], REAL(sB)[0], REAL(sPiA)[0], REAL(sLambdaA)[0], 
    REAL(sLambdaD)[0],  REAL(sNegCons)[0],
    1,
    REAL(Outs), REAL(Outs)+Rf_length(sQuantilesList), REAL(sMax));
  if (Go < 0 ) {
    Rprintf("TestSeekHPD, Go Returns %d \n", Go);  R_FlushConsole();
  }
  SET_VECTOR_ELT(sList, 0, sMax);
  SET_VECTOR_ELT(sList, 1, Outs);  
  Rf_unprotect(3);
  return(sList);
}
extern "C" {
RCPP_MODULE(modTwoLassoSexpCL) {
  using namespace Rcpp;
  
  class_<TwoLassoSexp>( "TwoLassoSexp" )
  .constructor< SEXP, SEXP, SEXP, SEXP, SEXP>()
  .field("Verbose", &TwoLassoSexp::Verbose)
  .field_readonly("Syys", &TwoLassoSexp::Syys )
  .field_readonly("Sxxs", &TwoLassoSexp::Sxxs )
  .field_readonly("SXtX", &TwoLassoSexp::SXtX )
  .field_readonly("SXtY", &TwoLassoSexp::SXtY )
  .field_readonly("SumYYSq", &TwoLassoSexp::SumYYSq )
  .field_readonly("DidIFail", &TwoLassoSexp::DidIFail )
  .property("X", &TwoLassoSexp::get_X, "X")
  .property("Y", &TwoLassoSexp::get_Y, "Y")
  .property("XTXFlag", &TwoLassoSexp::get_XTXFlag, "XTXFlag")
  .property("XtXFlag", &TwoLassoSexp::get_XTXFlag, "XTXFlag")
  .property("XtX", &TwoLassoSexp::get_XTX, "XtX sum of squares matrix (or allocated part) ")
  .property("XTX", &TwoLassoSexp::get_XTX, "XtX sum of squares matrix (or allocated part)")
  .method("XTXjj", &TwoLassoSexp::get_XTXjj, "XtX sum of squares matrix (or allocated part)")
  .property("XTY", &TwoLassoSexp::get_XTY, "XtY sum of leverage of X versus Y")
  .property("XtY", &TwoLassoSexp::get_XTY, "XtY sum of leverage of X versus Y")
  .property("XLC", &TwoLassoSexp::get_XLC, "location of X in kXFinder allocated vector/XTX")
  .property("XtResid", &TwoLassoSexp::get_XTResid, "t(X) %*% (Y - X %*% Beta) residual vector")
  .property("XTResid", &TwoLassoSexp::get_XTResid, "   alias for XtResid") 
  .property("XtYResid", &TwoLassoSexp::get_XTResid, "   alias for XtResid") 
  .property("LambdaAK", &TwoLassoSexp::get_LambdaAK, "Vector of Lambda_A values to use for active set shrinkage") 
  .property("SLambdaAK", &TwoLassoSexp::get_LambdaAK, "  alias for LambdaAK") 
  .property("LambdaDK", &TwoLassoSexp::get_LambdaDK, "Vector of Lambda_D values to use for active set shrinkage") 
  .property("SLambdaDK", &TwoLassoSexp::get_LambdaDK, "  alias for LambdaDK")
  .property("LamdaAK", &TwoLassoSexp::get_LambdaAK, "Vector of Lambda_A values to use for active set shrinkage") 
  .property("SLamdaAK", &TwoLassoSexp::get_LambdaAK, "  alias for LambdaAK") 
  .property("LamdaDK", &TwoLassoSexp::get_LambdaDK, "Vector of Lambda_D values to use for active set shrinkage") 
  .property("SLamdaDK", &TwoLassoSexp::get_LambdaDK, "  alias for LambdaDK")
  .property("L2ShrinkagePrior", &TwoLassoSexp::get_SL2ShrinkagePrior, "L2 Shrinkage Prior")
  .property("L2ShrinkageRecords", &TwoLassoSexp::get_SL2ShrinkageRecords, "L2 Shrinkage as found")
  .property("SigmaPrior", &TwoLassoSexp::get_SSigmaPrior, 
    &TwoLassoSexp::set_SSigmaPrior, "SigmaPrior")
  .property("PiAVectorInputs", &TwoLassoSexp::getPiAVectorInputs,
    "PiA Vector Inputs for CV")
  .property("SigmaVectorInputs", &TwoLassoSexp::getSigmaVectorInputs,
    "Sigma Vector Inputs for CV")
  .property("LambdaAKVectorInputs", &TwoLassoSexp::getLambdaAKInputs,
    "LambdaAK Inputs Matrix for CV")
  .property("LambdaDKInputs", &TwoLassoSexp::getLambdaDKInputs,
    "LambdaDK Inputs Matrix for CV")
  .property("StartBetaMatrix", &TwoLassoSexp::get_StartBetaMatrix, "StartBetaMatrix for CrossValidation")
  .property("FixKa", &TwoLassoSexp::get_FixKa, &TwoLassoSexp::set_FixKa,
    "FixKa that fixes total amount of values on")
  .property("MaxCauchy", &TwoLassoSexp::get_MaxCauchy,  
    &TwoLassoSexp::set_MaxCauchy, "Cauchy Convergence maximum number of iterations")
  .property("CauchyEpsilon", &TwoLassoSexp::get_CauchyEpsilon, 
    &TwoLassoSexp::set_CauchyEpsilon, "Get Cauchy Move Epsilon")
  .property("tt1", &TwoLassoSexp::get_tt1, &TwoLassoSexp:: set_tt1, "Get tt1")
  .property("tt2", &TwoLassoSexp::get_tt2, &TwoLassoSexp:: set_tt2, "Get tt2")  
  .property("Stt1", &TwoLassoSexp::get_Stt1, &TwoLassoSexp:: set_Stt1, "Get Stt1 -- pointer SEXP ")
  .property("Stt2", &TwoLassoSexp::get_Stt2, &TwoLassoSexp:: set_Stt2, "Get Stt2 -- pointer SEXP ")  
  .property("NoShrinkColumns", &TwoLassoSexp::get_NoShrinkColumns, &TwoLassoSexp::set_NoShrinkColumns, "Columns that shouldn't Shrink?")
  .property("ListFlags", &TwoLassoSexp::get_ListFlags, "List of Flags")
  .property("XjXjCIVector", &TwoLassoSexp::get_XjXj_CI_Vector, "Vector Important to calculating confidence intervals")
  .property("CCIVector", &TwoLassoSexp::get_C_CI_Vector, "Vector important to calculating confidence intervals, residual relation ship between coordinate j to all other coords in model.")
  .property("SigCIVector", &TwoLassoSexp::get_Sig_CI_Vector, "Sigma calculated on one coordinate j during calculation of confidence intervals.")  
  .property("PiAPrior", &TwoLassoSexp::get_SPiAPrior, 
    &TwoLassoSexp::set_SPiAPrior, "PiAPrior")
  .property("RPiAPrior", &TwoLassoSexp::get_RPiAPrior, 
    &TwoLassoSexp::set_RPiAPrior, "Random Coefficients PiAPrior")
  .property("InitKappaMem", &TwoLassoSexp::get_InitKappaMem,
    "Initial Memory Count that Powered CDO")
  .property("InitKKs", &TwoLassoSexp::get_InitKKs,
    "Initial Memory Count that 2Lasso would like to give CDO->InitKappaMem")
  .property("ForceXTXFlag", &TwoLassoSexp::get_InitKKs,
    "Flag to force XTX method")
  .method("ClearOnGammas", &TwoLassoSexp::ClearOnGammas, "Clear the OnGammas Field")
  .method("ReUpdateOnGammas", &TwoLassoSexp::ReUpdateOnGammas, "Reupdate the OnGammas Field based upon BBOn1, Lambdas")
  .method("RefreshOnGammas", &TwoLassoSexp::ReUpdateOnGammas, "Reupdate the OnGammas Field based upon BBOn1, Lambdas")
  //.property("LambdaAK", &TwoLassoSexp::get_LambdaAK, "Sequence of LambdaA values")
  //.property("LambdaDK", &TwoLassoSexp::get_LambdaDK, "Sequence of LambdaD values")
  .method("MakeXTResid", &TwoLassoSexp::MakeXTResid, "Re update XTResid")
  .method("UpdateXTResid", &TwoLassoSexp::MakeXTResid, "Re update XTResid")
  .method("ShowSFunction", &TwoLassoSexp::ShowSFunction, "Show S Function Ran")
  .method("SFunction", &TwoLassoSexp::ShowSFunction, "Show S Function Ran")
  .method("UpdateCoord", &TwoLassoSexp::run_UpDateCoord, "Run One UpdateCoord")  
  .method("RunUpdateCoord", &TwoLassoSexp::run_UpDateCoord, "Run One UpdateCoord") 
  .method("SupplyNewWeights", &TwoLassoSexp::SupplyNewWeights, "Refresh the Weights") 
  .method("RefreshWeights", &TwoLassoSexp::RefreshWeights, "RefreshWeights")   
  .property("InvSXX", &TwoLassoSexp::get_InvSXX, "Inverse XX from CDO")
  .property("DiagShrinkage", &TwoLassoSexp::get_DiagShrinkage, 
    &TwoLassoSexp::set_DiagShrinkage, "L2 Shrinkage in Coordinate Descent")
  .property("OnLoop", &TwoLassoSexp::get_OnLoop,  &TwoLassoSexp::set_OnLoop, "On Loop for Coordinate Descent")
  .property("OnCoord", &TwoLassoSexp::get_OnCoord, &TwoLassoSexp::set_OnCoord, "On Coordinate for Coordinate Descent")
  .property("TotMove", &TwoLassoSexp::get_TotMove, "Total Moved in last Coordinate Descent")
  .property("OnEp", &TwoLassoSexp::get_OnEp, "Total Moved in last Coordinate Descent Loop")  
  .property("XtResidModify", &TwoLassoSexp::get_XTResidModify, "XTResid Modify for each coordinate.")
  .property("BBOn1", &TwoLassoSexp::get_BBOn1, "Current B-Hat vector in algorithm")
  .property("OnBBOn1", &TwoLassoSexp::get_BBOn1, "Current B-Hat vector in algorithm")
  .property("CDOOnGammas", &TwoLassoSexp::get_CDOOnGammas, "Current On Gammas Value for the p coefficients in OnGammas")
  .property("OnGammas", &TwoLassoSexp::get_OnGammas, "Current On Gammas Value for the p coefficients")
  .property("OrderSeq", &TwoLassoSexp::get_OrderSeq, "Sequence of number of EM iterations to use per LambdaA/LambdaD value")
  .property("dfTNosie", &TwoLassoSexp::get_dfTNoise, "number of degrees of freedom from hypothesized T distribution for noise")
  .property("dfTNoise", &TwoLassoSexp::get_dfTNoise, "number of degrees of freedom from hypothesized T distribution for noise")
  .property("TDFNu", &TwoLassoSexp::get_dfTNoise, "number of degrees of freedom from hypothesized T distribution for noise")
  .property("RecOnBeta", &TwoLassoSexp::get_SRecOnBeta, "Records of OnBeta progress")
  .property("RecOnBetas", &TwoLassoSexp::get_SRecOnBeta, "Records of OnBeta progress")
  .property("RecBBOn1", &TwoLassoSexp::get_SRecBBOn1, "Records of BBOn1 progress") 
  .property("PiA", &TwoLassoSexp::get_PiA, &TwoLassoSexp::set_PiA, "PiA - potentially vector of setups")
  .property("Z", &TwoLassoSexp::get_Z, &TwoLassoSexp::set_Z, "Z binomial Logit 1 or zero ")
  .property("zzs", &TwoLassoSexp::get_Z, &TwoLassoSexp::set_Z, "Z binomial Logit 1 or zero ")
  .property("GLMZ", &TwoLassoSexp::get_GLMZ, "Current Z Vector in GLMCDO")
  .property("PiAVector", &TwoLassoSexp::get_PiA, &TwoLassoSexp::set_PiA, "PiA - potentially vector of setups")
  .property("OnPiA", &TwoLassoSexp::get_OnPiA, &TwoLassoSexp::set_OnPiA, "PiA")
  .property("OnRPiA", &TwoLassoSexp::get_OnRPiA, &TwoLassoSexp::set_OnRPiA, "RPiA")
  .property("RPiA", &TwoLassoSexp::get_OnRPiA, &TwoLassoSexp::set_OnRPiA, "RPiA")
  .property("NumGroups", &TwoLassoSexp::get_NumGroups, "Number Of Groups of Parameters (Group Option On)")
  .property("EndGroupList", &TwoLassoSexp::get_EndGroupList, "List End locations for Groups")
  .property("tauEndList", &TwoLassoSexp::get_EndGroupList, "List End locations for Groups")
  .property("EndGroupLocations", &TwoLassoSexp::get_EndGroupList, "List End locations for Groups")
  .property("FirstRandomIndex", &TwoLassoSexp::get_FirstRandomIndex, "Index of first position of Random Groups Beta")
  .property("iFirstRandom", &TwoLassoSexp::get_FirstRandomIndex, "Index of first position of Random Groups Beta")
  .property("GroupLambdaEstimates", &TwoLassoSexp::get_GroupLambdaEstimates, "Group Parameter Lambda")
  .property("RandBBOn1", &TwoLassoSexp::get_RandBBOn1, "Values of BBOn1 for random Coefficients")     
  .method("CDORunAlgorithm", &TwoLassoSexp::CDORunAlgorithm, "Run the Coordinate Descent Algorithm")  
  .method("SetupGroupTwoLasso", &TwoLassoSexp::SetupGroupTwoLasso, "Setup For Group 2Lasso Estimation")
  .property("iiWeights", &TwoLassoSexp::get_Weights,  "Weights")
  .property("GLMTotMove", &TwoLassoSexp::get_GLMTotMove, "GLM Total Move")
  .property("Weights", &TwoLassoSexp::get_Weights,  "Weights")
  .property("RSiiWeights", &TwoLassoSexp::get_Weights,  "Weights")
  .property("CDOWeights", &TwoLassoSexp::get_CDOWeights,  "CDO Weights")
  .property("Sigma", &TwoLassoSexp::get_Sigma, &TwoLassoSexp::set_Sigma, "Sigma")
  .property("SigmaSq", &TwoLassoSexp::get_Sigma, &TwoLassoSexp::set_Sigma, "Sigma")
  .property("sigma", &TwoLassoSexp::get_Sigma, &TwoLassoSexp::set_Sigma, "Sigma")
  .property("Lambda2", &TwoLassoSexp::get_OnLambda2, &TwoLassoSexp::set_OnLambda2, "Lambda 2, L2shrinkage parameter")
  .property("OnLambda2", &TwoLassoSexp::get_OnLambda2, &TwoLassoSexp::set_OnLambda2, "Lambda 2, L2shrinkage parameter")
  .property("RealSSigma", &TwoLassoSexp::get_RealSSigma, &TwoLassoSexp::set_RealSSigma, "Sigma")  
  .property("sigmaSqNoise", &TwoLassoSexp::get_Sigma, &TwoLassoSexp::set_Sigma, "Sigma")
  .property("SigmaNoiseSq", &TwoLassoSexp::get_Sigma, &TwoLassoSexp::set_Sigma, "Sigma")
  .property("SigmaVector", &TwoLassoSexp::get_Sigma, &TwoLassoSexp::set_Sigma, "Sigma - potentially vector of setups")
  .field_readonly("NonZero", &TwoLassoSexp::NonZero) 
  .field_readonly("NonZeroG", &TwoLassoSexp::NonZeroG) 
  .field_readonly("CountGroupSq", &TwoLassoSexp::CountGroupSq) 
  .field_readonly("CurSBBar", &TwoLassoSexp::CurSBBar) 
  .method("DataGLMIntegrity", &TwoLassoSexp::CDODataGLMIntegrity, "Check for integrity on GLMCDO")
  .method("DataIntegrity", &TwoLassoSexp::CDODataIntegrity, "Check for integity of CDO")
  .property("OnSigma", &TwoLassoSexp::get_OnSigma, &TwoLassoSexp::set_OnSigma, "Sigma")
  .property("sigmaNoiseSq", &TwoLassoSexp::get_Sigma, &TwoLassoSexp::set_Sigma, "Sigma")
  .property("OnKappaS", &TwoLassoSexp::get_OnKappaS, "Number of active variables in kXFinder")
  .property("OnKappaMem", &TwoLassoSexp::get_OnKappaMem, "Number of columns allocated for XtX/kXFinder")
  .property("MaximumAllocation", &TwoLassoSexp::get_MaximumAllocation, &TwoLassoSexp::set_MaximumAllocation, "Maximum of number columns to allocate for XtX/kXFinder")
  .property("MemoryFailFlag", &TwoLassoSexp::get_MemoryFailFlag, "Memory Fails")
  .method("ResetMemoryFail", &TwoLassoSexp::ResetMemoryFail, "Reset MemoryFail Flag")
  .method("DoubleMemory", &TwoLassoSexp::DoubleMemory, "Run Double Memory")
  .property("kXFinder", &TwoLassoSexp::get_kXFinder, "kXFinder")
  .property("AllkXFinder", &TwoLassoSexp::get_AllkXFinder, "All OnKappaMem elements of kXFinder")
  .property("Beta", &TwoLassoSexp::get_Beta, &TwoLassoSexp::set_Beta, "Beta")
  .property("GLMBeta", &TwoLassoSexp::get_GLMBeta, "Beta of the GLM Object")
  .property("ReturnBetas", &TwoLassoSexp::get_Beta,  "Beta")
  .property("OnBeta", &TwoLassoSexp::get_Beta, &TwoLassoSexp::set_Beta, "Beta")
  .property("OnBetas", &TwoLassoSexp::get_Beta, &TwoLassoSexp::set_Beta, "Beta")
  .property("n", &TwoLassoSexp::get_n,  "n number of samples")
  .property("nSample", &TwoLassoSexp::get_n, "n number of samples")
  .property("p", &TwoLassoSexp::get_p, "p")
  .property("pCoef", &TwoLassoSexp::get_p,  "p")
  .property("Verbose", &TwoLassoSexp::get_Verbose, &TwoLassoSexp::set_Verbose, "Verbose for TwoLasso Object")
  .property("PrintFlag", &TwoLassoSexp::get_PrintFlag, &TwoLassoSexp::set_PrintFlag, "Print Flag Verbosity for CDO Object")
  .property("CDOOnBeta", &TwoLassoSexp::get_CDOBeta, "Get the OnBeta")
  .property("CDOBeta", &TwoLassoSexp::get_CDOBeta, "Get the OnBeta")
  .property("kLen", &TwoLassoSexp::get_p, "p number of coefficients")
  .property("OnLambdaA", &TwoLassoSexp::get_OnLambdaA, "Current Lambda_A value in algorithm")
  .property("OnLambdaD", &TwoLassoSexp::get_OnLambdaD, "Current Lambda_D value in algorithm")
  .property("sRecBBOn1", &TwoLassoSexp::get_sRecBBOn1, "S-expression pointer for recorded BBOn1")
  .property("sRecOnBeta", &TwoLassoSexp::get_sRecOnBeta, "S-expression pointer for recorded OnBeta")
  .property("SuccessFlag", &TwoLassoSexp::get_SuccessFlag, "Get SuccessFlag!")
  .property("RecordBetaCVFinish", &TwoLassoSexp::get_RecordBetaCV, "Betas fit in Cross Validate")  
  .property("RecordBetaCV", &TwoLassoSexp::get_RecordBetaCV, "Betas fit in Cross Validate")
  .property("RecordBeta0CV", &TwoLassoSexp::get_RSRecordBeta0CV, "Beta0 fits if GLM Binomial is performed")
  .property("InverseGammaConstant", &TwoLassoSexp::get_InverseGammaConstant,
    &TwoLassoSexp::set_InverseGammaConstant, "Invese tweak to Gammas effect on CDO")
  .property("CDOInverseGammaConstant", &TwoLassoSexp::get_CDOInverseGammaConstant,
    &TwoLassoSexp::set_InverseGammaConstant, "Inverse tweak to Gammas effect on CDO")
  .property("FactorGammaDenominator", &TwoLassoSexp::get_CDOInverseGammaConstant,
    &TwoLassoSexp::set_InverseGammaConstant, "Inverse tweak to Gammas effect on CDO")   
  .property("LambdaIndexForConfidenceIntervals", &TwoLassoSexp::get_LambdaIndexForConfidenceIntervals,
    &TwoLassoSexp::set_LambdaIndexForConfidenceIntervals, "Lambda Index For Confidence Intervals")
  .property("ConfidenceQuantiles", &TwoLassoSexp::get_ConfidenceQuantiles,
    &TwoLassoSexp::set_ConfidenceQuantiles, "What Quantiles to seek for Credibility?")
  .property("CredibilityQuantiles", &TwoLassoSexp::get_ConfidenceQuantiles,
    &TwoLassoSexp::set_ConfidenceQuantiles, "What Quantiles to seek for Credibility?")
  .property("LogitNoise", &TwoLassoSexp::get_LogitNoise, &TwoLassoSexp::set_LogitNoise,
    "Noise during Logit for relaxation of ideals")
  .property("HPDMatrix", &TwoLassoSexp::get_HPDMatrix, "Get HPD format intervals")
  .property("HPDQuantiles", &TwoLassoSexp::get_HPDQuantiles, "Get HPD format intervals")
  .property("ConfidenceMatrix", &TwoLassoSexp::get_ConfidenceMatrix,
    "Estimates for credibility  Quantiles")
  .property("CredibilityMatrix", &TwoLassoSexp::get_ConfidenceMatrix,
    "Estimates for credibility  Quantiles")
  .property("UnshrunkConfidenceMatrix", &TwoLassoSexp::get_UnshrunkConfidenceMatrix,
    "Estimates for credibility  Quantiles without near zero shrinking")
   .property("UnshrinkCredibilityMatrix", &TwoLassoSexp::get_UnshrunkConfidenceMatrix,
    "Estimates for credibility  Quantiles without near zero shrinking")
   .property("FailureFlag", &TwoLassoSexp::get_FailureFlag,
      &TwoLassoSexp::set_FailureFlag,
    "Set Failure Flag")   
  .property("PrintFlagGLMB", &TwoLassoSexp::get_PrintFlagGLMB,
    &TwoLassoSexp::set_PrintFlagGLMB,
    "Print Flag GLMB")
  .property("jti", &TwoLassoSexp::get_jti,
    "jti for Cross Validation")
  .field("TargetMinimumBeta", &TwoLassoSexp::TargetMinimumBeta, 
    "Minimum Beta To set SLambdaAK")
  .method("UpdateWeightLambdaALambdaDSeq", &TwoLassoSexp::UpdateWeightLambdaALambdaDSeq,
    "Update Weight of LambdaA, LambdaD")
  .property("Beta0", &TwoLassoSexp::get_Beta0, &TwoLassoSexp::SetBeta0, 
    "Beta0 intercept value")
  .property("CalculateBeta0", &TwoLassoSexp::get_CalculateBeta0,
    "Beta0 intercept after calculate")
  .property("Beta0Prev", &TwoLassoSexp::get_Beta0Prev,
    "Last move Beta0Prev intercept value")
  .property("RBetasPrev", &TwoLassoSexp::get_RBetasPrev,
    "Recording in GLMBinomial of previous Beta")
  .property("BetasPrev", &TwoLassoSexp::get_RBetasPrev,
    "Recording in GLMBinomial of previous Beta")
  .property("RProbWeights", &TwoLassoSexp::get_RProbWeights,
    "Probability Weights for GLM Binomial regression")
  .property("RBeta0Records", &TwoLassoSexp::get_RBeta0Records,
    "Beta 0 records of intercept value")
  .property("RLogitCurrentProb", &TwoLassoSexp::get_RLogitCurrentProb,
    "Logit of current probability")
  .property("nlogitCurrentProb", &TwoLassoSexp::get_nlogitCurrentProb,
    "Logit of current probability in GLM")
  .property("RLogitPrevProb", &TwoLassoSexp::get_RLogitPrevProb,
    "Logit of previous moves probability")
  .property("ProbWeights", &TwoLassoSexp::get_RProbWeights,
    "Probability Weights for GLM Binomial regression")
  .property("Beta0Records", &TwoLassoSexp::get_RBeta0Records,
    "Beta 0 records of intercept value")
  .property("LogitCurrentProb", &TwoLassoSexp::get_RLogitCurrentProb,
    "Logit of current probability")
  .property("LogitPrevProb", &TwoLassoSexp::get_RLogitPrevProb,
    "Logit of previous moves probability")
  .property("GiveCDOIter", &TwoLassoSexp::get_GiveCDOIter,
    "Loop Iteration that GLMBinom object is on.")
  .property("RsGroupLambdaRecord", &TwoLassoSexp::get_RsGroupLambdaRecord,
    "Get records of old Group Lambdas")
  .property("GroupLambdaRecord", &TwoLassoSexp::get_RsGroupLambdaRecord,
    "Get records of old Group Lambdas")
  .method("TestCDO", &TwoLassoSexp::TestCDO, "Test Coordinate Descent Object integrity")
  .method("TestAllInt", &TwoLassoSexp::TestAllInt,
      "Test Integration Method")
  .method("TestSuperInt", &TwoLassoSexp::TestSuperInt,
      "Test Integration Method")
  .method("TestSeekQuantile", &TwoLassoSexp::TestSeekQuantile,
      "Test Find Integration Method")
  .field("CVQuitTime", &TwoLassoSexp::CVQuitTime, "CVQuitTime")
  .property("GroupsBBOn1", &TwoLassoSexp::get_GroupsBBOn1, "Group BBOn1 at final fit")
  .property("BackGroupsBBOn1", &TwoLassoSexp::get_BackGroupsBBOn1, "Group BBOn1 at final fit")
  .property("y01", &TwoLassoSexp::get_y01, "Integer y zero or 1")
  .method("RefreshGLMWeights", &TwoLassoSexp::RefreshGLMWeights, "Refresh GLM Weights")
  .method("RefreshCurrentLogitProb", &TwoLassoSexp::RefreshCurrentLogitProb, "Refresh Current Logit Prob for GLM")
  .method("RefreshBBOn1", &TwoLassoSexp::RefreshBBOn1, "Set BBOn1 values to PiA")
  .method("ARefreshBBOn1", &TwoLassoSexp::RefreshBBOn1, "Set BBOn1 values to PiA")
  .method("DoRefreshBBOn1", &TwoLassoSexp::RefreshBBOn1, "Set BBOn1 values to PiA")  
  .method("SetupTDFNu", &TwoLassoSexp::SetupTDFNu, "Set noise to t distributed and weights")
  .method("SetupRecords", &TwoLassoSexp::SetupRecords, "Setup Recording TwoSexp performance")
  .method("SetupL2Shrinkage", &TwoLassoSexp::SetupL2Shrinkage, "Setup L2 Shrinkage")
  .method("Setuptt", &TwoLassoSexp::Setuptt, "Sets the tt times, tt1 and tt2")
  .method("SetupSumYYSq", &TwoLassoSexp::SetupSumYYSq, "Calculates YSq if it must be calculated")
  .method("GiveYYSq", &TwoLassoSexp::GiveYYSq, "Give us YSq if it must be given")
  .method("UpdateBeta", &TwoLassoSexp::UpdateBeta, "Updates the Beta using CDO algorithm")
  .method("OneRoundUpdateBeta", &TwoLassoSexp::OneRoundUpdateBeta, "Updates the Beta One round using CDO algorithm")
  .method("OneCoordUpdateBeta", &TwoLassoSexp::OneCoordUpdateBeta, "Updates the Beta One Coordinate using CDO algorithm")
  //.method("RefreshBBOn1", &TwoLassoSexp::RefreshBBOn1, "Refresh BBOn1 vector")
  .method("UpdateBBOn1", &TwoLassoSexp::UpdateBBOn1, "Updates BBOn1 vector")
  .method("SetupPiSigma", &TwoLassoSexp::SetupPiSigma, "Set in PiA, Sigma that are NULL")
  .method("UpdateSigmaSq", &TwoLassoSexp::UpdateSigmaSq, "Update Sigma^2 estimate")
  .property("SigmaBar", &TwoLassoSexp::get_SigmaBar, "SigmaBar")
  .property("SigmaDf", &TwoLassoSexp::get_SigmaDf, "SigmaDf")  
  .property("SigmaDF", &TwoLassoSexp::get_SigmaDf, "SigmaDf")
  .method("UpdateTDFNoise", &TwoLassoSexp::UpdateTDFNoise, "Update T Noise based discounting in Coordinate Descent Object")
  .property("TDFSigma", &TwoLassoSexp::get_TDFSigma, " TDFSigma inside the CDO.  ")
  .property("CDOTDFNoise", &TwoLassoSexp::get_CDOTDFNoise, " CDOTDFNoise inside the CDO.  ")  
  .method("UpdateHatPiA", &TwoLassoSexp::UpdateHatPiA, "Gives PiA estimate")      
  .method("SetupCauchy", &TwoLassoSexp::SetupCauchy, "Sets Cauchy Epsilon and Max iters")
  .method("SetupInputs", &TwoLassoSexp::SetupInputs, "Supply Input LambdaDK/LambdaAK/StartBeta/PiA/Sigma for Cross Validation")
  .method("RunTwoLassoRegression", &TwoLassoSexp::RunTwoLassoRegression, "Run Two Lasso Regression Algorithm")
  .method("RunGLMLasso", &TwoLassoSexp::RunGLMLasso, "Run GLMLasso on GLMBinom object")
  .method("SetupCDO", &TwoLassoSexp::RunSetupCDO, "Setup CDO")
  .method("AfterCDO", &TwoLassoSexp::AfterCDO, "After CDO configure")
  //.field("PrintFlagGLMB", &TwoLassoSexp::get_PrintFlagGLMB, &TwoLassoSexp::set_PrintFlagGLMB, "Print Flag for GLMCDO")
  .method("RunCrossValidate", &TwoLassoSexp::RunCrossValidate, "Run Cross Validate Algorithm")
  .method("TestGetVal", &TwoLassoSexp::TestGetVal, "Test Get a Validation")
  .method("CVDataIntegrity", &TwoLassoSexp::CVDataIntegrity, "Test Memory Elements for CV")
  .method("TestSeekLevel", &TwoLassoSexp::TestSeekLevel, "Run Cross Validate Algorithm")  
  .method("TestSeekHPD", &TwoLassoSexp::TestSeekHPD, "Run Cross Validate Algorithm")    
  .method("TestMaxPosterior", &TwoLassoSexp::TestMaxPosterior, "Seek point of max posterior")
  .method("TestSecDet", &TwoLassoSexp::TestSecDet, "Calculate Second Derivative at some points")  
  .method("TestLDet", &TwoLassoSexp::TestLDet, "Calculate First Derviative at some points")      
  .method("SetupLambda", &TwoLassoSexp::SetupLambda, "Supply LambdaAK/LambdaDK ") 
  .method("ReweightCoordinate", &TwoLassoSexp::TReweightCoordinate, "Run a Reweighting of a coordinate in CDO")
  .property("iWeightedXtX", &TwoLassoSexp::get_iWeightedXtX, "iWeightedXtX"); 
     
}}

  

