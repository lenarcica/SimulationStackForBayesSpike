/* ========================================================================== */
/*                                                                            */
/*   TwoLassoObject2014.h                                                     */
/*   (c) 2011 Alan Lenarcic                                                   */
/*                                                                            */
/*   TwoLasso code that takes advantage of SEXPS                              */
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



#ifndef RCPPH
  #include <Rcpp.h>
  #include <cmath>
  #include "Acpp.h"
  #define RCPPH 0
#endif


#ifndef ACPPH
  #include "Acpp.h"
  #define ACPPH 0
#endif 

#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif

#ifndef LARSOBJECT2009DD
  #include "LarsObject2009.h"
  #define LARSOBJECT2009DD 0
#endif 

#ifndef ANINT
  #define ANINT(x,ii) ( (int) Rf_isInteger((x)) ? INTEGER((x))[(ii)] : (int) REAL((x))[(ii)] )
#endif
#ifndef AREAL
  #define AREAL(x,ii) ( (int) Rf_isReal((x)) ? REAL((x))[(ii)] : (double) INTEGER((x))[(ii)] )
#endif

#ifndef get_Attrib
  #define get_Attrib(X,Y)  Rf_getAttrib((X), (Y))
#endif 

#ifndef RGET_DIM
  #define RGET_DIM(X) Rf_getAttrib(X, R_DimSymbol)
#endif

#ifndef Rf_isNull
  #define Rf_isNull(X)  Rf_isNull(X)
#endif
#ifndef isInteger
  #define isInteger(X) Rf_isInteger(X)
#endif
//#ifndef ifReal
//  #define Rf_isReal(X) Rf_isReal(X)
//#endif
//#ifndef length
//  #define length(X) Rf_length(X)
//#endif
#ifndef error
  #define error(X) Rf_error(X)
#endif

#ifndef LARGEDATA
  #define LARGEDATA 999
#endif

#ifndef CoordinateDescentLassoDD
  #include "CoordinateDescentLasso2014.h"
  #define CoordinateDescentLassoDD 0
#endif

#ifndef MyMatrixOp2009DD
  #include "MyMatrixOp2009.h"
  #define MyMatrixOp2009DD 0
#endif


#ifndef MAXINITKKs
  #define MAXINITKKs 1000
#endif 

#define SUFFEPSILON_FD .0000001
#define MAXITERS_FD  140
#define KLENCOPY 200

#ifndef GetFirstInteger
  #define GetFirstInteger(sO) isInteger(sO) ? INTEGER(sO)[0] : (int) REAL(sO)[0]
#endif

#ifndef GetFirstDouble
  #define GetFirstDouble(sO) Rf_isReal(sO) ? REAL(sO)[0] :  \
    isInteger(sO) ? (double) INTEGER(sO)[0] \
    : 0 
#endif



#ifndef DDelete
  #define DDelete(X, charS)  if (Verbose > 2) {  \
    Rprintf("deleting  %s \n ", charS); R_FlushConsole(); } \
    if ( (X) != NULL ) {  delete((X));  (X) = NULL; }
#endif

#ifndef HEAPSORTH
  #define HEAPSORTH 1
  #include "HeapSort.h"
#endif
#define getInputsVector( RSNameInput, getNameI) \
  SEXP getNameI() {      \
    if (RSNameInput == NULL) {  \
      return(R_NilValue);        \
    }                             \
    return(RSNameInput->asSexp()); \
  }                                 \
  
inline SEXP GiveADoubleVectorOut(double *XXX, int L, int Verbose, const char *charS) {
  if (Verbose > 2)  {
    Rprintf(" Making a Double Vector for %s \n", charS); R_FlushConsole(); 
  } 
  SEXP sOut = R_NilValue; int One = 1;  int Len = L;
  if (XXX == NULL) { return(R_NilValue); }
    Rf_protect(sOut = Rf_allocVector(REALSXP,Len)); 
    F77_CALL(dcopy)(&Len, XXX, &One, REAL(sOut), &One); 
  Rf_unprotect(1); 
  return(sOut);
}

double SeekQuantile(double QuantileSeek, double StartAt, double Epsilon,
  double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons);

class TwoLassoSexp {
  private:
  /////////////////////////////////////////////////////////////////////////////
  // RListFlags - 10 length vector of scalar inputs to TwoLassoSexp() for Rcpp
  //   loadin. 
  AObject *RListFlags; SEXP ListFlags;
  
  // Related to Credibility measure: sets which LambdaA in LambdaAK sequence
  //  to choose as the variance for credibility measures
  int LambdaIndexForConfidenceIntervals;
    
  public: 
  
  ///////////////////////////////////////////////////////////////
  //  These are Cross validation input vectors
  //
  //    If these are filled, then cross validation will cycle through
  //     The values.
  SEXP SPiAVectorInputs;  SEXP SSigmaVectorInputs;
  SEXP SLambdaAKInputs;  SEXP SLambdaDKInputs;
  SEXP SRecordBetaCV; SEXP SStartBetaMatrix;
  AObject* RSPiAVectorInputs;  AObject* RSSigmaVectorInputs;
  AObject* RSLambdaAKInputs;  AObject*  RSLambdaDKInputs;
  AObject* RSStartBetaMatrix;
  AObject*  RSRecordBetaCV;
  AObject* RSRecordBeta0CV;
  
  /////////////////////////////////////////////////////////////////////////////
  //  The following go into GLMBin object if we are doing logit regression.
  double Beta0; double Beta0Prev; 
  AObject* RBetasPrev;
  AObject* RBeta0Records;
  AObject* RProbWeights; 
  AObject* RLogitCurrentProb;
  AObject* RLogitPrevProb;
  int CVDataIntegrity(int DoP);

  int *GiveCDOIter;

  
  ////////////////////////////////////////////////////////////////////////////
  //  These are related to Credibility measures
  //
  AObject *XjXj_CI_Vector;
  AObject *C_CI_Vector; AObject *Sig_CI_Vector;
  // A_CI, C_CI, Sig_CI store information about local variance for each
  //   coordinate Beta_j for j in 1:p related to confidence variance
  AObject* ConfidenceMatrix;
  AObject* HPDMatrix;  AObject *HPDQuantiles;
  AObject* UnshrunkConfidenceMatrix;
  AObject* ConfidenceQuantiles;
  AObject *NoShrinkColumns;
  int DidIFail;
  

  ////////////////////////////////////////////////////////////////////
  //   Setuptt(SEXP rStt1, SEXP rStt2)
  //
  //   Attempts to set up iteration time as SEXPs
  void Setuptt(SEXP rStt1, SEXP rStt2) {
    DDelete(RStt1, "RStt1"); DDelete(RStt2, "RStt2");
    Stt1 = R_NilValue; Stt2 = R_NilValue;
    RStt1 = new AObject(Rf_allocVector(INTSXP,1));
    Stt1 = RStt1->asSexp();
    RStt2 = new AObject(Rf_allocVector(INTSXP,1));
    Stt2 = RStt2->asSexp();
    if (Rf_isNull(Stt1)  || Rf_length(Stt1) <= 0) {
      Rf_error("Setuptt: Damn weird Error, Stt1 is NULL!");
    }
    if (Rf_isNull(Stt2) || Rf_length(Stt2) <= 0) {
      Rf_error("Setuptt: Error, Stt2 is NULL!");
    }
    if (Rf_length(Stt1) <= 0) { Rf_error("Setuptt: Stt1 length <= 0\n");}
    if (Rf_length(Stt2) <= 0) { Rf_error("Setuptt: Stt2 length <= 0\n"); }
    if (Verbose >= 1) {
      Rprintf("Setuptt, you supply Stt1=%d, Stt2=%d!\n",
        GetFirstInteger(rStt1), GetFirstInteger(rStt2)); R_FlushConsole();
    }
    int TEntry =   ((int) GetFirstInteger(rStt1))-1;    
    if (Rf_isInteger(rStt1)) {
      if (TEntry+1 != INTEGER(rStt1)[0]) {
        Rprintf("Setuptt, we have TEntry = %d, GetfirstStt1-1=%d, but INT(Stt1)=%d!\n",
          TEntry, ((int) GetFirstInteger(rStt1))-1, INTEGER(rStt1)[0]);
        R_FlushConsole();
      }
    } else {
      if (TEntry+1 != (int) REAL(rStt1)[0]) {
        Rprintf("Setuptt, we have TEntry = %d, GetfirstStt1-1=%d, but REAL(Stt1)=%d!\n",
          TEntry, ((int) GetFirstInteger(rStt1))-1, (int) REAL(rStt1)[0]);
        R_FlushConsole();
      }    
    }
    tt1 = ((int) TEntry);   
    TEntry =   ((int) GetFirstInteger(rStt2))-1;    
    if (Rf_isInteger(rStt2)) {
      if (TEntry+1 != INTEGER(rStt2)[0]) {
        Rprintf("Setuptt, we have TEntry = %d, GetfirstStt2-1=%d, but INT(Stt2)=%d!\n",
          TEntry, ((int) GetFirstInteger(rStt2))-1, INTEGER(rStt2)[0]);
        R_FlushConsole();
      }
    } else {
      if (TEntry+1 != (int) REAL(rStt2)[0]) {
        Rprintf("Setuptt, we have TEntry = %d, GetfirstStt2-1=%d, but REAL(Stt2)=%d!\n",
          TEntry, ((int) GetFirstInteger(rStt2))-1, (int) REAL(rStt2)[0]);
        R_FlushConsole();
      }    
    }
    tt2 = ((int) TEntry);
    if (Rf_isInteger(Stt1)) { INTEGER(Stt1)[0] = tt1; }
    if (Rf_isReal(Stt1)) { REAL(Stt1)[0] = (double) tt1;  }
    if (Rf_isInteger(Stt2)) { INTEGER(Stt2)[0] = tt2; }
    if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = (double) tt2;  }
    if (Verbose >= 3) {
    if (Rf_isReal(Stt1)) {
      Rprintf("Setuptt:  Oh doh, Stt1 is REAL = %d\n", (int) REAL(Stt1)[0]); R_FlushConsole();
    }
    if (Rf_isReal(Stt2)) {
      Rprintf("Setuptt:  Oh doh, Stt2 is REAL = %d\n", (int) REAL(Stt2)[0]); R_FlushConsole();
    }
    if (Rf_isInteger(Stt1)) {
      Rprintf("Setuptt:  Oh doh, Stt1 is INTEGER = %d\n", INTEGER(Stt1)[0]); R_FlushConsole();
    }
    if (Rf_isInteger(Stt2)) {
      Rprintf("Setuptt:  Oh doh, Stt2 is INTEGER = %d\n", INTEGER(Stt2)[0]); R_FlushConsole();
    }    
    }
    if (Verbose >= 1) {
      if (Rf_isInteger(Stt1) && Rf_isInteger(Stt2)) {      
        Rprintf("Setuptt, After Stt1,Stt2 TEntry set, we have tt1 = %d, tt2 = %d, Stt1=%d, Stt2=%d.\n", 
          tt1, tt2, INTEGER(Stt1)[0], INTEGER(Stt2)[0] ); R_FlushConsole();
      } else if (Rf_isReal(Stt1) && Rf_isReal(Stt2)) {
        Rprintf("Setuptt, After Stt1,Stt2 TEntry set, we have tt1 = %d, tt2 = %d, Stt1=%d, Stt2=%d.\n", 
          tt1, tt2, (int) REAL(Stt1)[0], (int) REAL(Stt2)[0] ); R_FlushConsole();
      } else {
        Rprintf("Setuptt, After Stt1,Stt2 TEntry set, we have tt1 = %d, tt2 = %d, Stt1, Stt2 are bad.\n", 
          tt1, tt2, (int) REAL(Stt1)[0], (int) REAL(Stt2)[0] ); R_FlushConsole();    
      }
    }
    if (tt1 < 0) { tt1 = 0; }
    if (tt2 < 0) { tt2 = 0; }
    if (!Rf_isNull(SLambdaAK) && Rf_length(SLambdaAK) >= 1 && Rf_length(SLambdaAK) <= tt1) {
      tt1 = Rf_length(SLambdaAK)-1;
    }
    if (!Rf_isNull(SLambdaAK) && Rf_length(SLambdaAK) >= 1 && Rf_length(SLambdaAK) <= tt2) {
      tt2 = Rf_length(SLambdaAK)-1;
    }
    if (Verbose >= 1) {   
      if (Rf_isInteger(Stt1) && Rf_isInteger(Stt2)) {
        Rprintf("Setuptt, After SLambdaAK max set we have tt1 = %d, tt2 = %d, Stt1=%d, Stt2=%d\n", 
          tt1, tt2, (int) INTEGER(Stt1)[0], INTEGER(Stt2)[0] );
        R_FlushConsole();
      } else if (Rf_isReal(Stt1) && Rf_isReal(Stt2)) {
        Rprintf("Setuptt, After SLambdaAK max set we have tt1 = %d, tt2 = %d, (REAL) Stt1=%d, Stt2=%d\n", 
          tt1, tt2, (int) REAL(Stt1)[0], (int) REAL(Stt2)[0] ); R_FlushConsole();     
      } else {
        Rprintf("Setuptt, After SLambdaAK: tt1 = %d, tt2 = %d, weird Stt1 is something different from Stt2\n",
          tt1, tt2); R_FlushConsole();
      }
    } 
    if (Rf_isInteger(Stt1)) { INTEGER(Stt1)[0] = tt1+1; }
    if (Rf_isReal(Stt1)) { REAL(Stt1)[0] = (double) (tt1+1); }  
    if (Rf_isInteger(Stt2)) { INTEGER(Stt2)[0] = tt2+1; }
    if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = (double) (tt2+1); }  
    if (Verbose >= 1) {
      Rprintf("After Setuptt, At End tt1 = %d, tt2 = %d, Stt1=%d, Stt2=%d\n", 
        tt1, tt2, (int) (GetFirstInteger(Stt1)), (int) (GetFirstInteger(Stt2)) );R_FlushConsole();
    } 
  }


  SEXP TestSeekHPD(SEXP sQuantilesList, SEXP sEpsilon,
  SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD, SEXP sNegCons);
  double TestMaxPosterior(SEXP sA, SEXP sB, SEXP sPiA, 
    SEXP sLambdaA, SEXP sLambdaD, 
    SEXP sEpsilon, SEXP sNegCons );
  SEXP TestSecDet(SEXP sxx, SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD);
  SEXP TestLDet(SEXP sxx, SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD);
  double TestSeekLevel(SEXP sGoLeft1GoRight0, SEXP sLevelSeek, SEXP sUpBound, SEXP sDownBound, SEXP sEpsilon,
    SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD, SEXP sNegCons);
  SEXP TestGetVal(SEXP sAtX, SEXP sA, SEXP sB, SEXP sPiA, SEXP sLambdaA, SEXP sLambdaD, SEXP NegCons);
  
  //////////////////////////////////////////////////////////
  // getters for iteration time
  //
  //   tt1 is which LambdaAK value to choose
  //     tt2 is how many EM iterations on this LambdaAK value to perform.
  int get_tt1() { 
    if (Rf_isNull(Stt1)  || Rf_length(Stt1) <= 0) {
      Rf_error("Setuptt: Error, Stt1 is NULL!");
    }  
    return(tt1); }
  int get_tt2() { 
    if (Rf_isNull(Stt2) || Rf_length(Stt2) <= 0) {
      Rf_error("Setuptt: Error, Stt2 is NULL!");
    }  
    return(tt2); }
  int get_Stt1() { 
    if (Rf_isNull(Stt1)  || Rf_length(Stt1) <= 0) {
      Rf_error("Setuptt: Error, Stt1 is NULL!");
    }  
    return((int) GetFirstInteger(Stt1)); }
  int get_Stt2() { 
    if (Rf_isNull(Stt2) || Rf_length(Stt2) <= 0) {
      Rf_error("Setuptt: Error, Stt2 is NULL!");
    }  
    return((int) GetFirstInteger(Stt2)); }
  int jti;
  int get_jti() { return(jti); }
  void set_tt1(int ans) {
    if (Rf_isNull(Stt1) || Rf_length(Stt1) <= 0) {
      Rf_error("Set_tt1, run Setuptt First!\n");
    }
    if (Rf_isNull(SLambdaAK)) {
      Rf_error("Set_tt1, run SetupLambdaFirst!\n");
    }
    if (ans < 0) {
      Rf_error("Set_tt1, can't set less than zero!\n");
    }
    if (ans >= Rf_length(SLambdaAK)) {
      Rf_error("Set tt2, LambdaAK len=%d, can't set to %d\n",
        Rf_length(SLambdaAK), ans); 
    }
    tt1 = ans;  if (Rf_isReal(Stt1)) { REAL(Stt1)[0] = (double) tt1+1; }
    if (Rf_isInteger(Stt1)) { INTEGER(Stt1)[0] = tt1+1; }
    if (!Rf_isNull(SLambdaAK) && 
      Rf_length(SLambdaAK) > tt1 && Rf_isReal(SLambdaAK)) {
      OnLambdaA = REAL(SLambdaAK)[tt1];  OnLambdaD = REAL(SLambdaDK)[tt1];
    }
    if (!Rf_isNull(SSigma) && Rf_isReal(SSigma) && Rf_length(SSigma) > tt1 &&
      REAL(SSigma)[tt1] > 0.0) {
      OnSigma = REAL(SSigma)[tt1];
    }
  }
  void set_tt2(int ans) {
    if (Rf_isNull(Stt2) || Rf_length(Stt2) <= 0) {
      Rf_error("Set_tt2, run Setuptt First!\n");
    }
    if (Rf_isNull(SLambdaAK)) {
      Rf_error("Set_tt2, run SetupLambdaFirst!\n");
    }
    if (ans < 0) {
      Rf_error("Set_tt2, can't set less than zero!\n");
    }
    if (ans >= Rf_length(SLambdaAK)) {
      Rf_error("Set tt2, LambdaAK len=%d, can't set to %d\n",
        Rf_length(SLambdaAK), ans); 
    }
    tt2 = ans;  if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = (double) tt2+1; }
    if (Rf_isInteger(Stt2)) { INTEGER(Stt2)[0] = tt2+1; }
    if (CDO != NULL) {
      CDO->OnLoop = 0;  CDO->OnCoord = 0;
    }
  }
  void set_Stt1(int ans) {
    if (Rf_isNull(Stt1) || Rf_length(Stt1) <= 0) {
      Rf_error("Set_Stt1, run Setuptt First!\n");
    }
    if (Rf_isNull(SLambdaAK)) {
      Rf_error("Set_Stt1, run SetupLambdaFirst!\n");
    }
    if (ans < 1) {
      Rf_error("Set_Stt1, can't set less than one!\n");
    }
    if (ans > Rf_length(SLambdaAK)) {
      Rf_error("Set Stt1, LambdaAK len=%d, can't set to %d\n",
        Rf_length(SLambdaAK), ans); 
    }
    tt1 = ans-1;  if (Rf_isReal(Stt1)) { REAL(Stt1)[0] = (double) tt1+1; }
    if (Rf_isInteger(Stt1)) { INTEGER(Stt1)[0] = tt1+1; }
    OnLambdaA = REAL(SLambdaAK)[tt1];  OnLambdaD = REAL(SLambdaDK)[tt1];
  }
  void set_Stt2(int ans) {
    if (Rf_isNull(Stt2) || Rf_length(Stt2) <= 0) {
      Rf_error("Set_Stt2, run Setuptt First!\n");
    }
    if (Rf_isNull(SLambdaAK)) {
      Rf_error("Set_Stt2, run SetupLambdaFirst!\n");
    }
    if (ans < 1) {
      Rf_error("Set_Stt2, can't set less than one!\n");
    }
    if (ans > Rf_length(SLambdaAK)) {
      Rf_error("Set Stt2, LambdaAK len=%d, can't set to %d\n",
        Rf_length(SLambdaAK), ans); 
    }
    tt2 = ans-1;  if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = (double) tt2+1; }
    if (Rf_isInteger(Stt2)) { INTEGER(Stt2)[0] = tt2+1; }
  }
  
  ///////////////////////////////////////////////////////////////////
  // get_LambdaAK, get_LambdaDK, get_OrderSeq
  //
  //  These are 2Lasso sequences of LambdaA, LambdaD parameters and target
  //   number of EM iterations
  //
  SEXP get_LambdaAK() {
    if (Rf_isNull(SLambdaAK)) { return(R_NilValue); }
    SEXP sOn = R_NilValue;  int One = 1; int Len = Rf_length(SLambdaAK);
    Rf_protect(sOn = Rf_allocVector(REALSXP, Len));
    F77_CALL(dcopy)(&Len, REAL(SLambdaAK), &One, REAL(sOn), &One);
    Rf_unprotect(1); return(sOn);
  }
  SEXP get_LambdaDK() {
    if (Rf_isNull(SLambdaDK)) { return(R_NilValue); }
    SEXP sOn = R_NilValue;  int One = 1; int Len = Rf_length(SLambdaDK);
    Rf_protect(sOn = Rf_allocVector(REALSXP, Len));
    F77_CALL(dcopy)(&Len, REAL(SLambdaDK), &One, REAL(sOn), &One);
    Rf_unprotect(1); return(sOn);
  }        
  
  // Sets up RLambdaA/K vectors and OrderSeq for TwoLasso
  int SetupLambda(SEXP rSLambdaAK, SEXP rSLambdaDK, SEXP rSOrderSeq
  ) {
    if (Rf_length(rSLambdaAK) != Rf_length(rSLambdaDK) ||
      Rf_length(rSLambdaAK) != Rf_length(rSOrderSeq)) {
      Rf_error("SetupLambda: Error, OrderSeq, LambdaAK, SLambdaDK not right Rf_length !");  
    }
    DDelete(RSLambdaAK, "RSLambdaAK");
    RSLambdaAK = new AObject(rSLambdaAK);
    SLambdaAK = RSLambdaAK->asSexp();

    DDelete(RSLambdaDK, "RSLambdaDK");
    RSLambdaDK = new AObject(rSLambdaDK);
    SLambdaDK = RSLambdaDK->asSexp();   
    
    DDelete(RSOrderSeq, "RSOrderSeq");
    if (Rf_isReal(rSOrderSeq)) {
      if (Rf_length(rSOrderSeq) < Rf_length(rSLambdaAK)) {
        Rf_error("Sorry, OrderSeq is bad length, fail!");
      }
      Rf_protect(SOrderSeq = Rf_allocVector(INTSXP, Rf_length(rSOrderSeq)));
      RSOrderSeq = new AObject(SOrderSeq); Rf_unprotect(1); 
      SOrderSeq = RSOrderSeq->asSexp();
      for (int tt = 0; tt < Rf_length(SOrderSeq); tt++) {
        INTEGER(SOrderSeq)[tt] = (int) (REAL(rSOrderSeq)[tt]);
      }
      OrderSeq = INTEGER(SOrderSeq);
    } else if (Rf_isInteger(rSOrderSeq)) {
      RSOrderSeq = new AObject(rSOrderSeq);
      SOrderSeq = RSOrderSeq->asSexp();
      OrderSeq = INTEGER(SOrderSeq);
    } else {
      OrderSeq = NULL;  SOrderSeq = R_NilValue;
      Rf_error("SetupLambda:  Failed with improper OrderSeq");
    } 
    OnLambdaA = REAL(RSLambdaAK->asSexp())[0];
    OnLambdaD = REAL(RSLambdaDK->asSexp())[0];
    return(1);
  }
  
  //////////////////////////////////////////////////////////////////
  //  SetupInputs
  //
  //    For Cross Validation, these input vectors are vectors of values
  //    to consider and fit on a dataset.  Potentially, PiA/Sigma/LambdaA/LambdaD
  //    are all slightly dependent on each other, so these sequences are not
  //    Presumed to represent a 4d space
  void SetupInputs(   SEXP rSPiAVectorInputs,  SEXP rSSigmaVectorInputs,
    SEXP rSLambdaAKInputs,  SEXP rSLambdaDKInputs,
    SEXP rSStartBetaMatrix, SEXP rSRecordBetaCV) {
    int ii;
    int MaxLen = 0;   jti = 0;
    if (rSLambdaDKInputs == NULL || Rf_isNull(rSLambdaDKInputs)) {
      Rprintf("NOTE: SetupInputs: RSLambdaDKInputs null here!\n"); R_FlushConsole();
    }
   
    DDelete(RSPiAVectorInputs, "RSPiAVectorInputs");
    SPiAVectorInputs = R_NilValue;
    if (!Rf_isNull(rSPiAVectorInputs) && Rf_length(rSPiAVectorInputs) > 0) { 
    RSPiAVectorInputs = new AObject(rSPiAVectorInputs);
    SPiAVectorInputs  = RSPiAVectorInputs->asSexp();
    MaxLen = Rf_length(SPiAVectorInputs);
    }

    DDelete(RSSigmaVectorInputs, "RSSigmaVectorInputs");
    SSigmaVectorInputs = R_NilValue;
    if (!Rf_isNull(rSSigmaVectorInputs) && Rf_length(rSSigmaVectorInputs) > 0) {
      RSSigmaVectorInputs = new AObject(rSSigmaVectorInputs);
      SSigmaVectorInputs  = RSSigmaVectorInputs->asSexp();
      if (Rf_length(SSigmaVectorInputs) > MaxLen) {
        MaxLen = Rf_length(SSigmaVectorInputs);
      }    
    }
    if (Rf_length(SSigmaVectorInputs) != Rf_length(SPiAVectorInputs)) {
      Rprintf("Error: SSigmaVI has Rf_length = %d, PiA = Rf_length %d \n",
        Rf_length(SSigmaVectorInputs), Rf_length(SPiAVectorInputs));
    }
    
    DDelete(RSLambdaAKInputs, "RSLambdaAKInputs");
    SLambdaAKInputs = R_NilValue;
    if (!Rf_isNull(rSLambdaAKInputs) && Rf_length(rSLambdaAKInputs) > 0 &&
      ( (Rf_length(rSLambdaAKInputs) == Rf_length(SLambdaAK)) ||
      ( !Rf_isNull(RGET_DIM(rSLambdaAKInputs)) &&
        INTEGER(RGET_DIM(rSLambdaAKInputs))[0]  == Rf_length(SLambdaAK))  
      ) ) {
    RSLambdaAKInputs = new AObject(rSLambdaAKInputs);
    SLambdaAKInputs  = RSLambdaAKInputs->asSexp();
    } else {
      Rprintf("Uh Oh, before SLambdaAKInputs setup.");R_FlushConsole();
      if (Rf_isNull(rSLambdaDKInputs) || rSLambdaAKInputs == NULL) {
        Rprintf("SLambdaAKInputs was NULL!\n");R_FlushConsole();
      } else if (Rf_isNull(RGET_DIM(rSLambdaAKInputs))) {
        Rprintf("Doh, rSLambdaAKInputs has NULL dimension!\n");
      } else {
        Rprintf("SLambdaAKInputs not nULL, but length is %d. \n",
          Rf_length(rSLambdaAKInputs)); R_FlushConsole();
      }
      Rprintf("rSLambdaAKInputs:  Error Maybe? I'm not sure this ");
      Rprintf("is a permissible result.\n");    R_FlushConsole();
      SLambdaDKInputs = R_NilValue;  RSLambdaDKInputs = NULL;    
    }
    
    DDelete(RSLambdaDKInputs, "RSLambdaDKInputs");
    SLambdaDKInputs = R_NilValue;
    if (!Rf_isNull(rSLambdaDKInputs) && Rf_length(rSLambdaDKInputs) > 0  &&
      ( (Rf_length(rSLambdaDKInputs) == Rf_length(SLambdaDK)) ||
      ( !Rf_isNull(RGET_DIM(rSLambdaDKInputs)) &&
        INTEGER(RGET_DIM(rSLambdaDKInputs))[0]  == Rf_length(SLambdaDK))  
      ) ) {
      RSLambdaDKInputs = new AObject(rSLambdaDKInputs);
      SLambdaDKInputs  = RSLambdaDKInputs->asSexp();
    } else {
      Rprintf("Uh Oh, before SLambdaDKInputs setup.");R_FlushConsole();
      if (Rf_isNull(rSLambdaDKInputs) || rSLambdaDKInputs == NULL) {
        Rprintf("SLambdaDKInputs was NULL!\n");R_FlushConsole();
      } else if (Rf_isNull(RGET_DIM(rSLambdaDKInputs))) {
        Rprintf("Doh, rSLambdaDKInputs has NULL dimension!\n");
      } else {
        Rprintf("SLambdaDKInputs not nULL, but length is %d. \n",
          Rf_length(rSLambdaDKInputs)); R_FlushConsole();
      }
      Rprintf("rSLambdaDKInputs: Error Maybe? I'm not sure ");
      Rprintf("this is a permissible result.");   R_FlushConsole();
      SLambdaDKInputs = R_NilValue;  RSLambdaDKInputs = NULL;
    }
  
    DDelete(RSStartBetaMatrix, "RSStartBetaMatrix");
    SStartBetaMatrix = R_NilValue;
    if (rSStartBetaMatrix != NULL &&
      !Rf_isNull(rSStartBetaMatrix) && Rf_length(rSStartBetaMatrix) > 0 &&
      ( Rf_length(rSStartBetaMatrix) == p ||
      ( !Rf_isNull(RGET_DIM(rSStartBetaMatrix)) &&
        INTEGER(RGET_DIM(rSStartBetaMatrix))[0]  == p)  
      )) {
      RSStartBetaMatrix = new AObject(rSStartBetaMatrix);
      SStartBetaMatrix  = RSStartBetaMatrix->asSexp(); 
    }
    
    DDelete(RSRecordBetaCV, "RSRecordBetaCV"); 
    DDelete(RSRecordBeta0CV, "RSRecordBeta0CV");
    if (rSRecordBetaCV == NULL || Rf_isNull(rSRecordBetaCV)) {
      Rf_error("Error Setup Inputs: You didn't give RecordBetaCV! "); 
    }
    if (Verbose >= 1) {
      Rprintf("SetupInputs now to test for Dimension of BetaCV"); R_FlushConsole();
    }
    if (Rf_isNull(RGET_DIM(rSRecordBetaCV))) {
      Rf_error("SetupInputs: rSRecordBetaCV has no dimension!");
    }
    if (Rf_length(RGET_DIM(rSRecordBetaCV)) < 2) {
      Rf_error("SetupInputs: rSRecordBetaCV has short Dimension!");
    }
    if (Rf_isNull(RGET_DIM(rSRecordBetaCV)) ||
      INTEGER(RGET_DIM(rSRecordBetaCV))[0] != p ||
      INTEGER(RGET_DIM(rSRecordBetaCV))[1] < MaxLen) {
      Rprintf("Error Setup Record Beta: You didn't give right dim to");
      Rprintf(" RecordBetaCV! (%d,%d) but MaxLen = %d\n",
        INTEGER(RGET_DIM(rSRecordBetaCV))[0],
        INTEGER(RGET_DIM(rSRecordBetaCV))[1], MaxLen); 
      Rprintf(" Length PiAInputs = %d, Length Sigma = %d \n",
        Rf_length(SPiAVectorInputs), 
        Rf_length(SSigmaVectorInputs));
      Rf_error(" So we're going to error. ");
    }  
    if (Verbose >= 2) {
      Rprintf("SetupInputs: We are inserting RSRecordBetaCV\n");
    }      
    RSRecordBetaCV = new AObject(rSRecordBetaCV);
    RSRecordBeta0CV = NULL;
    SRecordBetaCV  = RSRecordBetaCV->asSexp();  
    if (Verbose >= 2) {
      Rprintf("SetupInputs: RecordBetaCV is good, dim=%d,%d. \n",
        INTEGER(RGET_DIM(SRecordBetaCV))[0],
        INTEGER(RGET_DIM(SRecordBetaCV))[1]); R_FlushConsole();      
    }
    int LambdaALen = 0; int LambdaDLen = 0; int LambdaTLen = 0;  
    if (Verbose >= 1) {
      Rprintf("SetupInputs: Decide what to do with SLambdaDK/AK\n");
      R_FlushConsole();
    }
    if (Rf_isNull(SLambdaDK) || Rf_isNull(SLambdaAK)) {
       Rprintf("SetupInputs:: Setting up LambdaDK/ LambdaAK in Setup Inputs as they are null\n");
            int LambdaTLen; int LambdaALen; int LambdaDLen; 
      if (Rf_isNull(RGET_DIM(SLambdaAKInputs))) {
        LambdaALen = Rf_length(SLambdaAKInputs);
      } else {
        LambdaALen = INTEGER(RGET_DIM(SLambdaAKInputs))[0];
      }
      if (Rf_isNull(RGET_DIM(SLambdaDKInputs))) {
        LambdaDLen = Rf_length(SLambdaDKInputs);
      }  else {
        LambdaDLen = INTEGER(RGET_DIM(SLambdaDKInputs))[0];
      }
      if (LambdaALen > LambdaDLen) { 
        LambdaTLen = LambdaALen;
      } else { LambdaTLen = LambdaDLen; }
      if (Verbose >= 3) {
        Rprintf("SetupInputs, we have LambdaALen =%d,DLen=%d, TLen=%d\n",
          LambdaALen, LambdaDLen, LambdaTLen); R_FlushConsole();
      }
      Rf_protect(SLambdaAK = Rf_allocVector(REALSXP, LambdaTLen));
      DDelete(RSLambdaAK, "RSLambdaAK");
      RSLambdaAK = new AObject(SLambdaAK); Rf_unprotect(1);
      SLambdaAK = RSLambdaAK->asSexp();
      Rf_protect(SLambdaDK = Rf_allocVector(REALSXP, LambdaTLen));
      DDelete(RSLambdaDK, "RSLambdaDK");
      RSLambdaDK = new AObject(SLambdaDK); Rf_unprotect(1);
      SLambdaDK = RSLambdaDK->asSexp();    
      for (ii = 0; ii < LambdaALen; ii++) {
        REAL(SLambdaAK)[ii] = REAL(SLambdaAKInputs)[ii];
      }  
      if (LambdaALen < LambdaTLen) {
        for (ii = LambdaALen; ii < LambdaTLen; ii++) {
          REAL(SLambdaAK)[ii] = REAL(SLambdaAKInputs)[LambdaALen-1];
        }
      }
      for (ii = 0; ii < LambdaDLen; ii++) {
        REAL(SLambdaDK)[ii] = REAL(SLambdaDKInputs)[ii];
      }  
      if (LambdaDLen < LambdaTLen) {
        for (ii = LambdaDLen; ii < LambdaTLen; ii++) {
          REAL(SLambdaDK)[ii] = REAL(SLambdaDKInputs)[LambdaDLen-1];
        }
      } 
    } else {     
      LambdaALen = Rf_length(SLambdaAK);
      LambdaDLen = Rf_length(SLambdaDK);
      LambdaTLen = LambdaALen;  if (LambdaDLen > LambdaTLen) {
        LambdaTLen = LambdaDLen;
      }     
      if (Verbose >=3) {
        Rprintf("SetupInputs: From LambdaAK/DK have ALen=%d, DLen=%d, TLen=%d\n",
          LambdaALen, LambdaDLen, LambdaTLen); R_FlushConsole();
      }
      if (Rf_isNull(SLambdaAKInputs) || !Rf_isReal(SLambdaAKInputs)) {
        Rf_error("SetupInputs: have LambdaALen=%d, but LambdaAKInputs are NULL\n",
          LambdaALen); 
      } else if (Rf_length(SLambdaAKInputs) < LambdaALen) {
        Rf_error("SetupInputs: LambdaALen = %d, but Inputs has length %d\n",
          LambdaALen, Rf_length(SLambdaAKInputs));
      }
      for (ii = 0; ii < LambdaALen; ii++) {
        REAL(SLambdaAK)[ii] = REAL(SLambdaAKInputs)[ii];
      }  
      if (LambdaALen < LambdaTLen) {
        for (ii = LambdaALen; ii < LambdaTLen; ii++) {
          REAL(SLambdaAK)[ii] = REAL(SLambdaAKInputs)[LambdaALen-1];
        }
      }
      if (Rf_isNull(SLambdaDKInputs) || !Rf_isReal(SLambdaDKInputs)) {
        Rf_error("SetupInputs, have LambdaDLen=%d, but SLambdaDKInputs are NULL\n", LambdaDLen); 
      } else if (Rf_length(SLambdaDKInputs) < LambdaALen) {
        Rf_error("SetupInputs, LambdaDLen = %d, but Inputs has length %d\n",
          LambdaDLen, Rf_length(SLambdaDKInputs));
      }
      for (ii = 0; ii < LambdaDLen; ii++) {
        REAL(SLambdaDK)[ii] = REAL(SLambdaDKInputs)[ii];
      }  
      if (LambdaDLen < LambdaTLen) {
        for (ii = LambdaDLen; ii < LambdaTLen; ii++) {
          REAL(SLambdaDK)[ii] = REAL(SLambdaDKInputs)[LambdaDLen-1];
        }
      }     
    }
    if (Verbose >= 3) {
      Rprintf("SetupInputs: Allocate BackXTYResid and BackOrigBeta\n");
      R_FlushConsole();
    }
    BackXTYResid = (double *) Calloc(this->p+1, double); 
    BackOrigBeta = (double *) Calloc(this->p+1, double); 
    int One = 1;
    if (SOnBeta == NULL || Rf_isNull(SOnBeta) || !Rf_isReal(SOnBeta)) {
      Rf_error("SetupInputs: SOnBeta is Null, can't copy.");
    }
    if (Rf_length(SOnBeta) < p) {
      Rf_error("SetupInputs: SOnBeta only has length %d, not p=%d\n",
        Rf_length(SOnBeta), p); 
    }
    F77_CALL(dcopy)(&p, CDO->XTResid, &One, BackXTYResid, &One);
    F77_CALL(dcopy)(&p, REAL(SOnBeta), &One, BackOrigBeta, &One);    
  }
  getInputsVector( RSPiAVectorInputs , getPiAVectorInputs )
  getInputsVector( RSSigmaVectorInputs , getSigmaVectorInputs )
  getInputsVector( RSLambdaAKInputs , getLambdaAKInputs )
  getInputsVector( RSLambdaDKInputs , getLambdaDKInputs )
  getInputsVector( RSRecordBetaCV , getBetaCVInputs )
        

  int n,p;
  
  int TestCrossValidate(); 
  int VerbLessTestCrossValidate();
   int TestCDO() {
    int AT=0;
    if (CDO != NULL) {
      AT = CDO->TestCDO();
    }
    return(AT);
  }
  int CVQuitTime; 
  int SuccessFlag;
  int get_SuccessFlag() { return(this->SuccessFlag); }  
  SEXP Syys; SEXP Sxxs;
  AObject *RSyys;  AObject *RSxxs;
  int *zzs;
  SEXP SXtY; SEXP SXtX;
  AObject *RSXtY;  AObject *RSXtX;
  SEXP SLambdaAK;  SEXP SLambdaDK;
  AObject *RSLambdaAK;  AObject *RSLambdaDK;
  SEXP SOrderSeq;  int *OrderSeq;
  AObject *RSOrderSeq; 
  
  SEXP SSigmaPrior;  AObject *RSSigmaPrior;  
  void set_SSigmaPrior(SEXP rSSigmaPrior) {
    if (Rf_length(rSSigmaPrior) != 2) {
      Rf_error("set_SSigmaPrior: faulty input");
    }
    DDelete(RSSigmaPrior, "RSSigmaPrior");
    RSSigmaPrior = new AObject(rSSigmaPrior);
    SSigmaPrior = RSSigmaPrior->asSexp();
    SigmaBar = REAL(SSigmaPrior)[0];
    SigmaDf = REAL(SSigmaPrior)[1];
  }
  //int get_Verbose() { return(Verbose); }
  int get_PrintFlag() { if (CDO != NULL) { return(CDO->PrintFlag);} 
    Rprintf("PrintFlag: CDO not setup!\n"); R_FlushConsole();
    return(-1);
  }
  //void set_Verbose(int iVerbose) {
  //  Verbose = iVerbose;
  //}
  void set_PrintFlag(int iPrintFlag) {
    if (CDO != NULL) {
      CDO->PrintFlag = iPrintFlag; return;
    }
    Rprintf("set PrintFlag: CDO not setup!\n"); R_FlushConsole();
  }
  SEXP get_SSigmaPrior() {
    if (RSSigmaPrior != NULL) {
      return(RSSigmaPrior->asSexp()); 
    }
    Rprintf("RSSigma Prior is not set yet, here is SigmaBar, SigmaDf anyway.\n");
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, 2));
    REAL(sOut)[0] = SigmaBar;  REAL(sOut)[1] = SigmaDf;
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_X() { return(Sxxs);}
  SEXP get_Y() { return(Syys);}
  SEXP get_Z() {
    if (zzs == NULL)  { return(R_NilValue); }
    SEXP sOut = R_NilValue;  Rf_protect(sOut = Rf_allocVector(INTSXP, n));
    for (int ii = 0; ii < n; ii++) { INTEGER(sOut)[ii] = zzs[ii]; }
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_GLMZ() {
    if (GLMCDO == NULL)  { Rf_error("Error: GLM Is Not Setup. \n"); }
    return(GLMCDO->get_Z());
  }
  void set_Z(SEXP iZ) {
    if (Rf_isNull(iZ) || Rf_length(iZ) != n) {
      Rf_error("Error: set zzs, can't do unless length = %d, not %d\n",
        n, Rf_length(iZ)); 
    }
    if (zzs != NULL) { Free(zzs); zzs = NULL; }
    zzs = (int *) Calloc(n, int);
    if (Rf_isReal(iZ)) {
      for (int ii = 0; ii < n; ii++) {
        zzs[ii] = (int)  REAL(iZ)[ii];
      }
    } else if (Rf_isInteger(iZ)) {
      for (int ii = 0; ii < n; ii++) {
        zzs[ii] = (int)  INTEGER(iZ)[ii];
      }    
    }
  }
  
  
  SEXP SBBOn1; SEXP SOnBeta;  
  AObject *RSBBOn1;  AObject *RSOnBeta;
  SEXP SRecBBOn1; SEXP SRecOnBeta;
  AObject *RSRecBBOn1; AObject *RSRecOnBeta;
  int RecordFlag;
  SEXP SL2ShrinkagePrior;
  SEXP SL2ShrinkageRecords;
  AObject *RSL2ShrinkagePrior;  AObject *RSL2ShrinkageRecords;
  
  SEXP Stt1, Stt2;
  AObject *RStt1;  AObject *RStt2;
  
  SEXP SOnGammas;
  AObject *RSOnGammas;
  SEXP SPiA, SSigma;
  AObject *RSPiA;  AObject *RSSigma;
  double OnPiA, OnSigma; double OnRPiA;
  
  double m1, m2, SigmaBar, SigmaDf;
  double m3, m4;
  int ForceXTXFlag;
  
  void set_SPiAPrior(SEXP rSPiAPrior) {
    if (Rf_length(rSPiAPrior) != 2) {
      Rf_error("set_SPiAPrior: faulty input");
    }
    m1 = REAL(rSPiAPrior)[0]; m2 = REAL(rSPiAPrior)[1];
  }
  SEXP get_SPiAPrior() { SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, 2));
    REAL(sOn)[0] = m1;  REAL(sOn)[1] = m2;  Rf_unprotect(1);
    return(sOn);  
  }
  void set_RPiAPrior(SEXP rSRPiAPrior) {
    if (Rf_length(rSRPiAPrior) != 2) {
      Rf_error("set_RPiAPrior: faulty input");
    }
    m3 = REAL(rSRPiAPrior)[0]; m4 = REAL(rSRPiAPrior)[1];
  }
  SEXP get_RPiAPrior() {
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, 2));
    if (m3 > 0 && m4 > 0) {
      REAL(sOn)[0] = m3;  REAL(sOn)[1] = m4;  
    } else {
      Rprintf("Note RPiA Prior must use m1 and m2");
      REAL(sOn)[0] = m1;  REAL(sOn)[1] = m2;  
    }
    Rf_unprotect(1);
    return(sOn);  
  }
  	
  CoordinateDescentObject *CDO;
  GLMBinomObject *GLMCDO;
  //int ChoiceCDOOnLARSFlag;
  //LARSObject  *OnLARS;
  
  double OnLambdaA, OnLambdaD;
  int tt1, tt2;
  double InverseGammaConstant;
  int InitKKs;
  
  int Verbose;
  
  int MaxCauchy; double CauchyEpsilon;
  void SetupCauchy(SEXP rsMaxCauchy, SEXP rsCauchyEpsilon) {
    MaxCauchy = GetFirstInteger(rsMaxCauchy);
    CauchyEpsilon = REAL(rsCauchyEpsilon)[0];
    this->CDO->SetEpsilonAndLoops(CauchyEpsilon, MaxCauchy);
  }

  int RunTwoLassoRegression(int DoHistory);
  int RunGLMTwoLassoRegression(int DoHistory);
  double get_InverseGammaConstant() {
    return(InverseGammaConstant);
  }
  double get_CDOInverseGammaConstant() {
    return(CDO->FactorGammaDenominator);
  }
  void set_InverseGammaConstant(double rInverseGammaConstant) {
    InverseGammaConstant = rInverseGammaConstant;
    if (CDO != NULL) {
      CDO->FactorGammaDenominator = InverseGammaConstant;
    }
    return;
  }
  int UpdateBeta();
  int UpdateBBOn1();
  int OneRoundUpdateBeta();
  int OneCoordUpdateBeta(int DoCoord);

  int UpdateSigmaSq();  
  int NonZero; double NonZeroG;  double CurSBBar, CountGroupSq;
  int UpdateHatPiA();
  
  double *BackXTYResid; double *BackOrigBeta; 
  void SetupOrigBeta(double *OrigBeta, double *OnXTYResid) {
    if (BackXTYResid == NULL) {
      BackXTYResid = (double *) Calloc( p+1, double );
    }
    if (BackOrigBeta == NULL) {
      BackOrigBeta = (double *) Calloc( p+1, double );    
    }
    int One = 1;  
    F77_CALL(dcopy)(&p, OrigBeta, &One, CDO->OnBeta, &One);
    F77_CALL(dcopy)(&p, OnXTYResid, &One, BackXTYResid, &One);    
  } 
  void SetupOrigBetaFromCDO() {
    if (BackXTYResid == NULL) {
      BackXTYResid = (double *) Calloc( p+1, double );
    }
    if (BackOrigBeta == NULL) {
      BackOrigBeta = (double *) Calloc( p+1, double );    
    }
    int One = 1;  
  
    F77_CALL(dcopy)(&p, CDO->OnBeta, &One, BackOrigBeta, &One);
    F77_CALL(dcopy)(&p, CDO->XTResid, &One, BackXTYResid, &One);
    
  } 
  
   int FirstRandomIndex;
   int NumGroups;
   double *GroupLambdaEstimates;
   int *EndGroupLocations; 
   AObject *RsGroupLambdaRecord; 
   AObject *RsGroupLambdaEstimates;
   AObject *GroupsBBOn1;
   AObject *BackGroupsBBOn1;
   int UpdateGroupBBOn1();
   int FailureFlag;
   int get_FailureFlag () {
      return(FailureFlag);
   }
   void set_FailureFlag(int iFailureFlag) {
     FailureFlag = iFailureFlag; }
   
  int RunTwoLasso(int Starttt1);
  
  int SetupCDO(int SetBeta); int AfterCDO();
  
  SEXP sGroupsSexp;
  AObject *RsGroupsSexp;
  int SetupGroups(SEXP rsGroupsSexp);
  int RecordHistory();
  
  int *RunCoords;
  

  
  double SumYYSq;
  void GiveYYSq(double NewYYSq) {
    if ((Syys != NULL) && !Rf_isNull(Syys)) {
      Rf_error("You gave yys, so I'm not going to trust you");
    }
    SumYYSq = NewYYSq;
  }
  int SetupSumYYSq() {
   int One = 1;
   if (Syys == NULL || Rf_isNull(Syys)) {Rf_error("Cannot SumYYSq if Syys is NULL!");}
   if (Rf_length(Syys) < n) { Rf_error("Cannot, because length Syys = %d, n = %d\n", 
     Rf_length(Syys), n); }
   SumYYSq = 0.0;
   if (iiWeights == NULL) {
     for (One = 0; One < n; One++) {
       SumYYSq += REAL(Syys)[One] * REAL(Syys)[One];
     }
     return(1);
   } else {
     for (One = 0; One < n; One++) {
       SumYYSq += REAL(Syys)[One] * REAL(Syys)[One] * iiWeights[One];
     }
     return(1);
   }
   //SumYYSq = F77_CALL(ddot)(&n, REAL(Syys), &One, REAL(Syys), &One);
   return(1);
  }
  SEXP get_OrderSeq() {
    if (OrderSeq == NULL || Rf_isNull(SLambdaAK) || Rf_length(SLambdaAK) <= 0) {
      return(R_NilValue);
    }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, Rf_length(SLambdaAK)));
    for (int ii = 0; ii < Rf_length(SLambdaAK); ii++) {
      INTEGER(sOn)[ii] = OrderSeq[ii];
    }
    Rf_unprotect(1);  
    return(sOn);
  }
  //int ForceUpdateBBOn1;
  //void set_ForceUpdateBBOn1(int AGo) {
  //  if (AGo == 1) { ForceUpdateBBOn1 = 1; }
  //  if (Ago == 0) { ForceUpdateBBOn1 = 0; }
  //}
  //int get_ForceUpdateBBOn1() { return(ForceUpdateBBOn1); }
  int get_ForceXTXFlag() { return(ForceXTXFlag); }
  int RefreshBBOn1() {
    if (Rf_isNull(SBBOn1)) {
      Rf_error("Error, can't RefreshBBOn1 because it is NULL!\n");
    }
    if (NumGroups <= 0 || FirstRandomIndex >= p || FirstRandomIndex < 0) {
    for (int ii = 0; ii < p; ii ++)  {
      REAL(SBBOn1)[ii] = OnPiA;
    }
    return(1);   
    }
    if (NumGroups >= 0 && FirstRandomIndex > 0) {
      for (int ii = 0; ii < FirstRandomIndex; ii++)  {
        REAL(SBBOn1)[ii] = OnPiA;
      }     
    } 
      double ROnPiA = OnPiA;
      if (this->OnRPiA > 0.0 && this->OnRPiA < 1.0) {
        ROnPiA = this->OnRPiA;
      } 
      for (int ii = FirstRandomIndex; ii < p; ii++) {
        REAL(SBBOn1)[ii] = ROnPiA;
      }
      for (int ii = 0; ii < NumGroups; ii++) {
        REAL(GroupsBBOn1->asSexp())[ii] =  ROnPiA;
      }
    //ReUpdateOnGammas();
    return(1);
  }
  //SEXP get_LambdaAK() { return(SLambdaAK);}
  //SEXP get_LambdaDK() { return(SLambdaDK);}
  SEXP get_CDOOnGammas() {
    if (CDO == NULL) { Rprintf("Error: get_CDOOnGammas: CDO is NULL!\n"); R_FlushConsole();}
    if (CDO->OnGammas == NULL) { Rprintf("Error: get_CDOOnGammas: CDO->OnGammas = NULL!\n"); R_FlushConsole();}
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, p));
    int One = 1;
    F77_CALL(dcopy)(&p, CDO->OnGammas, &One, REAL(sOn), &One);
    Rf_unprotect(1); return(sOn);
  }
  SEXP get_OnGammas() {
    if (CDO == NULL) { return(R_NilValue); }
    if (!Rf_isNull(SOnGammas)) {
      if (Verbose >=2) {
        Rprintf("Returning SOnGammas!\n"); R_FlushConsole();
      }
      if (REAL(SOnGammas) != CDO->OnGammas) {
        Rprintf("UhOhUhOhUhOhUhOhUhOh"); 
        Rprintf("Note SOnGammas pointer is different from CDO->OnGammas pointer!\n");
        R_FlushConsole();
      }
      return(SOnGammas);
    }
    if (CDO->OnGammas == NULL) { 
      if (Verbose >= 3) {   
        Rprintf("CDO->OnGammas is Null, get_OnGammas no show!\n"); R_FlushConsole();
      }
      return(R_NilValue); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, p));
    int One = 1;
    if (Verbose >= 3) {
      Rprintf("get_OnGammas, using dcopy to copy p=%d, One=1, OnGammas to sOn[0]= %f.\n",
        p, One, CDO->OnGammas[0]); R_FlushConsole();
    }
    F77_CALL(dcopy)(&p, CDO->OnGammas, &One, REAL(sOn), &One);
    Rf_unprotect(1); return(sOn);  
  }
  void ClearOnGammas() {
    if (CDO == NULL) { return; }
    if (!Rf_isNull(SOnGammas)) {
      for (int ii = 0; ii < Rf_length(SOnGammas); ii++) {
        REAL(SOnGammas)[ii] = 0.0;
      }
      return;
    }
    if (CDO->OnGammas == NULL) { 
      Rprintf("WARNING WARNINGS: Can't clear OnGammas, OnGammas not setup, still NULL!\n"); 
      return;
    }
    for (int ii = 0; ii < p; ii++) {
      CDO->OnGammas[ii] = 0.0;
    }
    return;   
  }
  int ReUpdateOnGammas() {
    if (CDO == NULL) { 
      Rf_error("Can't ReupdateOnGammas because CDO is NULL!\n"); 
    }
    if (CDO->OnGammas == NULL) {
      Rf_error("Can't ReupdateOnGammas because CDO->OnGammas is NULL\n"); R_FlushConsole();
    }
    if (Rf_isNull(SBBOn1)) {
      Rf_error("No SBBOn1 can't ReUpdateOnGammas\n");
    }
    if (Verbose >= 2) {
      Rprintf("-- ReUpdateOnGammas, OnLambdaA = %f, OnLambdaD = %f, p=%d, OnSigma=%f\n",
        OnLambdaA, OnLambdaD, p, OnSigma); R_FlushConsole();
    }
    if (NumGroups <= 0) {
    if (InverseGammaConstant != 1.0  && InverseGammaConstant > 0.0) {
	    for (int ii = 0; ii < p; ii++) {
	      CDO->OnGammas[ii] = (REAL(SBBOn1)[ii] * OnLambdaA + 
	        (1.0- REAL(SBBOn1)[ii]) *OnLambdaD)  * OnSigma * 
	        InverseGammaConstant;
	    }
    } else {
      for (int ii = 0; ii < p; ii++) {
	      CDO->OnGammas[ii] = (REAL(SBBOn1)[ii] * OnLambdaA + 
	        (1.0- REAL(SBBOn1)[ii]) *OnLambdaD) * OnSigma;
  	  }
    }
    } else {
    if (FirstRandomIndex > 0) {
    if (InverseGammaConstant != 1.0  && InverseGammaConstant > 0.0) {
	    for (int ii = 0; ii < FirstRandomIndex; ii++) {
	      CDO->OnGammas[ii] = (REAL(SBBOn1)[ii] * OnLambdaA + 
	        (1.0- REAL(SBBOn1)[ii]) *OnLambdaD)  * OnSigma * 
	        InverseGammaConstant;
	    }
    } else {
      for (int ii = 0; ii < FirstRandomIndex; ii++) {
	      CDO->OnGammas[ii] = (REAL(SBBOn1)[ii] * OnLambdaA + 
	        (1.0- REAL(SBBOn1)[ii]) *OnLambdaD) * OnSigma;
  	  }
    }} 
    int St = FirstRandomIndex;
      if (InverseGammaConstant != 1.0  && InverseGammaConstant > 0.0) {
	      for (int ii = 0; ii < NumGroups; ii++) {
          GroupLambdaEstimates[ii] = REAL(GroupsBBOn1->asSexp())[ii]  * OnLambdaA +
           (1.0-REAL(GroupsBBOn1->asSexp())[ii]) * OnLambdaD;
            for (int jj = St; jj < EndGroupLocations[ii]+1; jj++) {
	            CDO->OnGammas[ii] = (GroupLambdaEstimates[ii])  * OnSigma * 
	               InverseGammaConstant;
            }
          St = EndGroupLocations[ii]+1;
	      }
      } else {
	      for (int ii = 0; ii < NumGroups; ii++) {
          GroupLambdaEstimates[ii] = REAL(GroupsBBOn1->asSexp())[ii]  * OnLambdaA +
           (1.0-REAL(GroupsBBOn1->asSexp())[ii]) * OnLambdaD;
            for (int jj = St; jj < EndGroupLocations[ii]+1; jj++) {
	            CDO->OnGammas[jj] = (GroupLambdaEstimates[ii])  * OnSigma;
            }
          St = EndGroupLocations[ii]+1;
	      }
      }    
  }
    if (Verbose >= 4) {
      Rprintf("---- ReUpdateOnGammas: SBBOn1 = ");  PrintVector(REAL(SBBOn1), p);
      Rprintf("\n---- ReUpdateOnGammas: After Update, OnGammas = ");
        PrintVector(CDO->OnGammas, p); Rprintf("\n"); R_FlushConsole();
    }
    if (REAL(SOnGammas) == CDO->OnGammas) { return(1); }
    if (!Rf_isNull(SOnGammas)) {
      if (REAL(SOnGammas) != CDO->OnGammas) {
        Rprintf("UhOhUhOhUhOhUhOhUhOh"); 
        Rprintf("Note SOnGammas pointer is different from CDO->OnGammas pointer!\n");
        R_FlushConsole();
      }
      if (Verbose >= 4) {
        Rprintf("---- ReUpdateOnGammas, checking CDO->OnGammas, SOnGammas\n");
        for (int ii = 0; ii < p; ii++) {
          if (REAL(SOnGammas)[ii] != CDO->OnGammas[ii]) {
            Rprintf("----   SOnGammas[%d]=%.5f != CDO->OnGammas[%d]=%.5f\n",
              ii, REAL(SOnGammas)[ii], ii, CDO->OnGammas[ii]); R_FlushConsole();
          }
        }
      }
      for (int ii = 0; ii < p; ii++) {
        REAL(SOnGammas)[ii] = CDO->OnGammas[ii];
      }
    }
    if (NoShrinkColumns != NULL) {
      for (int ii = 0; ii < Rf_length(NoShrinkColumns->asSexp()); ii++) {
        REAL(SOnGammas)[INTEGER(NoShrinkColumns->asSexp())[ii]] = 0.0;
      }
    }
    return(1);
  }
  SEXP get_dfTNoise() { SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, 1));
    REAL(sOn)[0] = TDFNu;  Rf_unprotect(1); return(sOn);
  }
  

  /////////////////////////////////////////////////////////////
  //  TDFNu is a way of using t-distributed noise of 
  //   degrees of freedom TDFNu if TDFNu > 0
  //
  double TDFNu; double *iiWeights; SEXP STDFNu; 
  
  SEXP SiiWeights; AObject *RSiiWeights;
  
  int SetupTDFNu(SEXP rSTDFNu, SEXP rSiiWeights);  
  int SetupTD2Lasso(SEXP rSTDFNu, SEXP rSiiWeights);
  int RefreshTD2Lasso();
  
  
  ///////////////////////////////////////////////////////////////
  //  FixKa is set by a negative m1
  //    It is a method of fixing the exact number of coefficients
  //    to allow nonzero in the model
  int FixKa;  
  SEXP get_FixKa() { SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, 1));
    INTEGER(sOn)[0] = FixKa; Rf_unprotect(1); return(sOn); }
  int InitFixKa(int rFixKa);
  void set_FixKa(SEXP iFixKa) {
    FixKa = GetFirstInteger(iFixKa);  
    InitFixKa(FixKa);
  }
  int InsertBBOn1ofX(double Onmu);
  int FixKCalculateBBOn();
  double AofX(double Onmu);	
  double DAofX(double Onmu);
  double SolveFormu( double GoalTotal, 
         double SuffEpsilon, int MaxIters, double StartOnmu);	
  
  int TwoLassoPrintDiagnosticA();
  int OnLARSCDODiagnostic();
  int XTXFlag;
  
  int UpdateL2Shrinkage();
  int FitBBOnWithGroups();
  
  int RunCrossValidate();
  double TargetMinimumBeta;
  
  void SetupRecord(int RecordFlag) {
    if (Verbose >= 1) {
      Rprintf("**** SetupRecords: Version RecordFlag = %d\n", RecordFlag);
    }
    if (RSRecBBOn1 != NULL) { RSRecBBOn1 = NULL;  SRecBBOn1 = R_NilValue; }
    if (RSRecOnBeta != NULL) { RSRecOnBeta = NULL; SRecOnBeta = R_NilValue; }
    if (RSLambdaAK == NULL) {
      Rprintf("SetupRecords: Can't do if RSLambdaAK is NULL!\n"); 
      R_FlushConsole();
    }
    if (RecordFlag <= 0 || RecordFlag >= 3) {
      return;
    }
    if (RecordFlag == 1) {
      RSRecBBOn1 =  new AObject(
        Rf_allocMatrix(REALSXP, p, Rf_length(RSLambdaAK->asSexp())));
      RSRecOnBeta =  new AObject(
        Rf_allocMatrix(REALSXP, p, Rf_length(RSLambdaAK->asSexp()))); 
      SRecBBOn1 = RSRecBBOn1->asSexp();  SRecOnBeta = RSRecOnBeta->asSexp();
      return;       
    }
    if (OrderSeq == NULL || RSLambdaAK == NULL || Rf_isNull(RSLambdaAK->asSexp())) {
      Rf_error("SetupRecords: Can't do if OrderSeq is NULL\n");
    }
    int MOrderSeq = 0;
    for (int ii = 0; ii < Rf_length(RSLambdaAK->asSexp()); ii++) {
      if (OrderSeq[ii] > MOrderSeq) { MOrderSeq = OrderSeq[ii]; }
    }
    RSRecBBOn1 = new AObject(Rf_alloc3DArray(REALSXP,
      p, Rf_length(RSLambdaAK->asSexp()), MOrderSeq));
    RSRecOnBeta = new AObject(Rf_alloc3DArray(REALSXP,
      p, Rf_length(RSLambdaAK->asSexp()), MOrderSeq));      
    SRecBBOn1 = RSRecBBOn1->asSexp();  SRecOnBeta = RSRecOnBeta->asSexp();
  }
  void SetupRecords(SEXP rRecBBOn1, SEXP rRecOnBeta) {
    if (!Rf_isNull(rRecBBOn1) && Rf_length(rRecBBOn1) > 1) {
       if (RSRecBBOn1 != NULL) {
         if (Verbose >= 1) {
           Rprintf("SetupRecords: Deleting RSRecBBOn1 before anything.\n");
           R_FlushConsole();
         }
         delete(RSRecBBOn1); RSRecBBOn1 = NULL;  SRecBBOn1 = R_NilValue;
       }
       RSRecBBOn1 = new AObject(rRecBBOn1);
       SRecBBOn1 = RSRecBBOn1->asSexp();
    }
    if (!Rf_isNull(rRecOnBeta) && Rf_length(rRecOnBeta) > 1) {
       if (RSRecOnBeta != NULL) {
         if (Verbose >= 1) {
           Rprintf("SetupRecords: Deleting RSRecOnBeta before anything.\n");
           R_FlushConsole();
         }
         delete(RSRecOnBeta); RSRecOnBeta = NULL;  SRecOnBeta = R_NilValue;
       }
       RSRecOnBeta = new AObject(rRecOnBeta);
       SRecOnBeta = RSRecOnBeta->asSexp();
    }  
    return;
  } 
  
  SEXP get_StartBetaMatrix() {
    if (RSStartBetaMatrix == NULL) {
      return(R_NilValue);
    }
    return(RSStartBetaMatrix->asSexp());
  }
//  TwoLassoSexp<SEXP irSn, SEXP irSp> { 
  TwoLassoSexp () {
    RSyys = NULL; RSxxs = NULL; RSXtY = NULL;
    RSXtX = NULL; RStt1 = NULL; RStt2 = NULL; zzs = NULL;
    RSLambdaAK= NULL; RSLambdaDK = NULL; RSOrderSeq = NULL;
    RSBBOn1 = NULL; RSOnBeta = NULL; RSRecBBOn1 = NULL; RSRecOnBeta = NULL;
    RSOnGammas = NULL; RSPiA = NULL; RSSigma = NULL; 
    RSSigmaPrior = NULL;  BackOrigBeta = NULL;
    RSiiWeights = NULL;  SiiWeights = R_NilValue; RSL2ShrinkagePrior = NULL;  
    RSL2ShrinkageRecords = NULL;  RBeta0Records = NULL;
    RSStartBetaMatrix = NULL;     RecordFlag = -1;
    SStartBetaMatrix = R_NilValue; // ForceUpdateBBOn1 = 1;
    InverseGammaConstant = 1.0;  ConfidenceMatrix = NULL; 
    HPDMatrix = NULL; HPDQuantiles=NULL;
    UnshrunkConfidenceMatrix = NULL; XjXj_CI_Vector = NULL;  C_CI_Vector=NULL;
    ConfidenceQuantiles = NULL; FailureFlag = 0; Sig_CI_Vector= NULL;
    CDO = NULL; GLMCDO = NULL; DidIFail = 0;
    NonZero = 0; NonZeroG = 0.0; CurSBBar= 0.0; CountGroupSq = 0.0;
    
    BackOrigBeta = NULL;  BackXTYResid = NULL; TargetMinimumBeta = 1.0;
    
    FirstRandomIndex = -1; NumGroups = -1;
    GroupLambdaEstimates = NULL;  EndGroupLocations = NULL; 
    RsGroupLambdaRecord = NULL; 
    RsGroupLambdaEstimates = NULL;   GroupsBBOn1 = NULL; BackGroupsBBOn1=NULL; 

    RBetasPrev = NULL; RBeta0Records = NULL;
    RSRecordBetaCV = NULL; RSRecordBeta0CV = NULL;  
    RProbWeights = NULL; RLogitCurrentProb = NULL;
    RLogitPrevProb = NULL;    GiveCDOIter = NULL;
    
    RStt1 = new AObject(Rf_allocVector(INTSXP,1));
    Stt1= RStt1->asSexp();  INTEGER(Stt1)[0] = 0;
    RStt2 = new AObject(Rf_allocVector(INTSXP,1));
    Stt2= RStt2->asSexp();  INTEGER(Stt2)[0] = 0;    
   // RSn = new AObject(irSn); sSn = RSn->asSexp();
   // Rsp = new AObject(irSp); sSp = RSp->asSexp();
  }
  void set_Verbose(int iVerbose) {
    Verbose = iVerbose;
  } 
  int get_Verbose() { return (Verbose); }
  //int get_p() { return(p); }
  void set_Y(SEXP rSyys) {
    if (Rf_isNull(Sxxs)) {
      Rf_error("set_Y, sorry, no setting Y without a real X");
    }
    int One = 1;
    DDelete(RSyys, "RSyys");  
    RSyys = new AObject(rSyys);
    Syys = RSyys->asSexp();

      char Trans = 'T';
      double OneD = 1.0; double ZeroD = 0.0;    
    if (!Rf_isNull(SXtY)  && Rf_length(SXtY) == p) {

      F77_CALL(dgemv)( &Trans, &n, &p, &OneD, 
        REAL(Sxxs), &n, REAL(Syys), &One, &ZeroD, REAL(SXtY), &One);
    } else if (CDO != NULL && CDO->XTY != NULL) {
      F77_CALL(dgemv)(&Trans, &n, &p, &OneD, 
         REAL(Sxxs), &n, REAL(Syys), &One, &ZeroD, CDO->XTY, &One); 
      F77_CALL(dscal)(&p, &ZeroD, CDO->OnBeta, &One);
      F77_CALL(dcopy)(&p, CDO->XTY, &One, CDO->XTResid, &One);       
    }   
    if (n <= 0) {
      Rf_error("setY:  Please supply n first to a constructor! \n"); R_FlushConsole();
    }
    if (Rf_length(rSyys) != n) {
      Rf_error("setY:  Y is not Rf_length %d \n", n);
    }

  }
  //SEXP get_Y() { return(RSyys->asSexp()); }

  void set_X(SEXP sX2) {
    Rf_error("Sorry, We are not allowing reset of X Matrix");
  }
  TwoLassoSexp(SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaAK, SEXP rSLambdaDK, SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
    SEXP rSVerbose) {
    RSyys = NULL; RSxxs = NULL; RSXtY = NULL; zzs = NULL;
    RSXtX = NULL; RStt1 = NULL; RStt2 = NULL; 
    Syys = R_NilValue;  Sxxs = R_NilValue;   
    Stt1 = R_NilValue;  Stt2 = R_NilValue;
    RSSigma = NULL;  SSigma = R_NilValue;
    NonZero = 0; NonZeroG = 0; CurSBBar= 0.0; CountGroupSq = 0.0;
    RSLambdaAK = NULL; RSLambdaDK = NULL; RSOrderSeq = NULL; OrderSeq = NULL;
    RSBBOn1 = NULL; RSOnBeta = NULL; RSRecBBOn1 = NULL; RSRecOnBeta = NULL;
    RBeta0Records = NULL;
    RSOnGammas = NULL; RSPiA = NULL; RSSigma = NULL; RSiiWeights = NULL; SiiWeights = R_NilValue;
    iiWeights = NULL;   RListFlags = NULL; ListFlags = R_NilValue;
    RSSigmaPrior = NULL;  RsGroupsSexp = NULL; sGroupsSexp =  R_NilValue;
    RSL2ShrinkagePrior = NULL;  RSL2ShrinkageRecords = NULL;
    m1=-1;m2=-1;  m3 = -1;  m4 = -1; ConfidenceQuantiles = NULL; ConfidenceMatrix = NULL;
    HPDMatrix=NULL;  HPDQuantiles=NULL; 
    UnshrunkConfidenceMatrix = NULL;   RecordFlag = -1; XjXj_CI_Vector=NULL; C_CI_Vector=NULL;
    InverseGammaConstant = 1.0; DidIFail = 0;   SigmaBar = -1.0;  SigmaDf = -1.0;
    CDO = NULL; GLMCDO = NULL;   RListFlags = NULL;  ListFlags = R_NilValue;
    CVQuitTime = -1;  jti = 0;  NoShrinkColumns = NULL; TargetMinimumBeta = 1.0;
    
    
    ForceXTXFlag = -1;
    
    BackOrigBeta = NULL; BackXTYResid = NULL;
    
    RBetasPrev = NULL; RBeta0Records = NULL;
    RSRecordBetaCV = NULL; RSRecordBeta0CV = NULL;  
    RProbWeights = NULL; RLogitCurrentProb = NULL;
    RLogitPrevProb = NULL;   GiveCDOIter = NULL;
    BackXTYResid = NULL; BackOrigBeta = NULL;  SXtY = R_NilValue;
    HPDMatrix=NULL;   HPDQuantiles = NULL;
    ConfidenceQuantiles = NULL;  ConfidenceMatrix = NULL;
    UnshrunkConfidenceMatrix = NULL;
    LambdaIndexForConfidenceIntervals = 0; FailureFlag = 0;
    
    SPiAVectorInputs = R_NilValue;    SSigmaVectorInputs = R_NilValue;
    SLambdaAKInputs = R_NilValue;    SLambdaDKInputs = R_NilValue;
    SRecordBetaCV = R_NilValue; BackOrigBeta = NULL; 
    RSPiAVectorInputs = NULL;   RSSigmaVectorInputs  = NULL;
    RSLambdaAKInputs = NULL;   RSLambdaDKInputs = NULL;
    RSRecordBetaCV = NULL; RSRecordBeta0CV = NULL;  
    RSStartBetaMatrix = NULL;  RSSigma = NULL; SSigma=R_NilValue;
    SStartBetaMatrix = R_NilValue; //  ForceUpdateBBOn1 = 1;
    
    RunCoords = NULL; SumYYSq = 0;
    
    SL2ShrinkagePrior = R_NilValue; RunCoords = NULL;
    SL2ShrinkageRecords = R_NilValue; RSL2ShrinkageRecords = NULL;
    m1=-1;m2=-1;  m3 = -1; m4 = -1;
    CDO = NULL; GLMCDO = NULL; FailureFlag = 0;
      
    RBetasPrev = NULL; RBeta0Records = NULL;
    RProbWeights = NULL; RLogitCurrentProb = NULL;
    RLogitPrevProb = NULL;  GiveCDOIter = NULL;
    
    InverseGammaConstant = 1.0;
    
    ConfidenceQuantiles = NULL;  ConfidenceMatrix = NULL;
    HPDMatrix=NULL;  HPDQuantiles = NULL;
    UnshrunkConfidenceMatrix = NULL; XjXj_CI_Vector=NULL; C_CI_Vector=NULL;
    LambdaIndexForConfidenceIntervals = 0;
    
    Beta0 = 0.0; Beta0Prev = 0.0; 
    RBetasPrev = NULL; RBeta0Records = NULL;
    RProbWeights = NULL; RLogitCurrentProb = NULL;
    RLogitPrevProb = NULL;  GiveCDOIter = NULL;
    //double *rnProbWeights, double *rnlogitCurrentProb, double *rnlogitPrevProb,
    //double *rZ, 
    
    SPiAVectorInputs = R_NilValue;    SSigmaVectorInputs = R_NilValue;
    SLambdaAKInputs = R_NilValue;    SLambdaDKInputs = R_NilValue;
    SRecordBetaCV = R_NilValue;
    RSPiAVectorInputs = NULL;   RSSigmaVectorInputs  = NULL;
    RSLambdaAKInputs = NULL;   RSLambdaDKInputs = NULL;
    RSRecordBetaCV = NULL; RSRecordBeta0CV = NULL; 
    RSStartBetaMatrix = NULL;
    SStartBetaMatrix = R_NilValue;  Stt1 = R_NilValue;   RSStartBetaMatrix = NULL;
    SSigma = R_NilValue;    // ForceUpdateBBOn1 = 1;
    
    FirstRandomIndex = -1; NumGroups = -1;
    GroupLambdaEstimates = NULL; EndGroupLocations = NULL; 
    RsGroupLambdaRecord = NULL; 
    RsGroupLambdaEstimates= NULL; GroupsBBOn1 = NULL;  BackGroupsBBOn1=NULL;

    
    SuccessFlag = 1;  XTXFlag = 0;
    if (Rf_isNull(rSVerbose)) {
      Rprintf("rSVerbose is NULL: Bad\n"); R_FlushConsole();
      Verbose = 5;
    } else {
      Verbose = GetFirstInteger(rSVerbose);
    }
    if (Verbose >= 0) {
      Rprintf("TwoLassoSexp: Beginning Large Input Constructor;\n"); R_FlushConsole();
    }
    
    //if (rStt1 == NULL || Rf_isNull(rStt1) || Rf_length(rStt1) <= 0) {
    //  Rf_error("TwoLassoSexp: Error, Stt1 is NULL!"); R_FlushConsole();
    //}
    //if (rStt2 == NULL || Rf_isNull(rStt2) || Rf_length(rStt2) <= 0) {
    //   Rf_error("TwoLassoSexp: Error, Stt2 is NULL!"); R_FlushConsole();
    // }
    RStt1 = new AObject(Rf_allocVector(INTSXP,1));
    Stt1 = RStt1->asSexp();
    RStt2 = new AObject(Rf_allocVector(INTSXP,2));
    Stt2 = RStt2->asSexp();
    if (!Rf_isNull(rStt1) && (Rf_isReal(rStt1) || Rf_isInteger(rStt1))) {
      INTEGER(Stt1)[0] = GetFirstInteger(rStt1);
    }
    if (!Rf_isNull(rStt2) && (Rf_isReal(rStt2) || Rf_isInteger(rStt2))) {
      INTEGER(Stt2)[0] = GetFirstInteger(rStt2);    
    }

    if (Rf_length(Stt1) >= 1) { tt1 = GetFirstInteger(Stt1); }    
    if (Rf_length(Stt2) >= 1) { tt2 = GetFirstInteger(Stt2); } 
    if (tt1 < 0) { tt1 = 0;  
      if (Rf_isInteger(Stt1)) { 
        INTEGER(Stt1)[0] = 0;
      } else if (Rf_isReal(Stt1)) { REAL(Stt1)[0] = 0.0;     }
    }
    if (tt2 < 0) { tt2 = 0;  
      if (Rf_isInteger(Stt2)) { 
        INTEGER(Stt1)[0] = 0;
      } else if (Rf_isReal(Stt2)) { REAL(Stt2)[0] = 0.0;     }
    }
    SumYYSq = 0;
    if (rSn == NULL || Rf_isNull(rSn) || Rf_length(rSn) <= 0) {
      Rf_error("TwoLassoSexp::new rSn is NULL: Very Very Bad\n"); R_FlushConsole();
      Verbose = 5;
    }
    if (Verbose >= 1) {
      Rprintf("About to Input --n-- sample size.\n"); R_FlushConsole();
    }
    n = GetFirstInteger(rSn);
    if (rSp == NULL || Rf_isNull(rSp) || Rf_length(rSp) <= 0) {
      Rf_error("TwoLassoSexp::new rSp is NULL: Very Very Bad\n"); R_FlushConsole();
      Verbose = 5;
    }
    p = GetFirstInteger(rSp);
   
    if (Verbose >= 1) {
      Rprintf("TwoLassoSexp:new: n = %d, p = %d\n",n,p); R_FlushConsole();
    }
    if (n < 0) {
      XTXFlag = 2;  n = abs(n);
    } else {
      XTXFlag = 1;  n = abs(n);
    }
    if (Verbose >= 1) {
      Rprintf("TwoLassSexp:new: XTXFlag = %d, n=%d, p=%d\n", XTXFlag, n, p);
      R_FlushConsole();
    }
    if (Rf_isNull(rSMaxCauchy)) {
      Rf_error("TwoLassoSexp:: rSMaxCauchy is NULL: What?\n"); R_FlushConsole();
      Verbose = 5;
    }     
    MaxCauchy = GetFirstInteger(rSMaxCauchy);
    if (Rf_isNull(rSCauchyEpsilon)) {
      Rf_error("TwoLassoSexp:: rSCauchyEpsilon is NULL: What?\n"); R_FlushConsole();
      Verbose = 5;
    }
    CauchyEpsilon = GetFirstDouble(rSCauchyEpsilon); 
    if (Verbose >= 1) {
      Rprintf("TwoLassoSexp::new, CauchyEpsilon = %f \n", CauchyEpsilon);
      R_FlushConsole();
    } 
    if (Verbose >= 1) {
      Rprintf("TwoLassoSexp::new, MaxCauchy = %d, and CauchyEpsilon = %f\n",
        MaxCauchy, CauchyEpsilon); R_FlushConsole();
    }


    if (Verbose > 4) {
      Rprintf("TwoLassoSexp: Starting with InitKKs, IGC\n"); R_FlushConsole();
    } 
    if (Rf_isNull(rSInitKKs) || Rf_length(rSInitKKs) <= 0 || 
      !(Rf_isReal(rSInitKKs) || Rf_isInteger(rSInitKKs))) {
      Rf_error("TwoLassoSexp:new: InitKKs is in Error \n");   
    }   
    InitKKs = GetFirstInteger(rSInitKKs);
    InverseGammaConstant = GetFirstDouble(rSInverseGammaConstant);
    if (Verbose > 4) {
      Rprintf("TwoLassoSexp: InitKKs = %d, InverseGammaConstant = %f;\n",
        InitKKs, InverseGammaConstant); R_FlushConsole();
    }
    if (Verbose > 4) {
      Rprintf("Reading in SLambdaA = "); R_FlushConsole();
      PrintVector(REAL(rSLambdaAK), Rf_length(rSLambdaAK)); Rprintf("\n");
      Rprintf("Reading in SLambdaD = "); R_FlushConsole();
      PrintVector(REAL(rSLambdaDK), Rf_length(rSLambdaDK)); Rprintf("\n");
    } 
    
    RSLambdaAK = new AObject(rSLambdaAK);
    SLambdaAK = RSLambdaAK->asSexp();  
    RSLambdaDK = new AObject(rSLambdaDK);
    SLambdaDK =  RSLambdaDK->asSexp(); 
    RSOrderSeq = new AObject(rSOrderSeq);
    SOrderSeq = RSOrderSeq->asSexp();
    RSPiA = new AObject(rSPiA);   
    SPiA = RSPiA->asSexp(); 
    DDelete(RSSigma, "RSSigma"); 
    RSSigma = new AObject(rSSigma);  
    SSigma = RSSigma->asSexp();
    
    int ii = 0;
    
    OnSigma = GetFirstDouble(rSSigma); OnPiA = GetFirstDouble(rSPiA);
    if (Verbose > 4) {
      Rprintf("TwoLassoSexp: OnSigma = %f, OnPiA = %f;\n",
        OnSigma, OnPiA); R_FlushConsole();
      Rprintf("  SpiA = "); PrintVector(REAL(SPiA), (int) Rf_length(SPiA)); R_FlushConsole();
      Rprintf("  Sigma = "); PrintVector(REAL(SSigma), (int) Rf_length(SSigma)); R_FlushConsole();
    }
    
    if (Rf_isNull(rSOnBeta) || Rf_length(rSOnBeta) <= 0) {
      Rf_error("TwoLassoSexp: SOnBeta is Bad NULL!\n"); R_FlushConsole();
    }
    if (Rf_length(rSOnBeta) != p) {
      Rf_error("TwoLassoSexp: SOnBeta has length %d does not equal p = %d\n",
        Rf_length(rSOnBeta), p); 
    }
    RSOnBeta = new AObject(rSOnBeta);    
    SOnBeta = RSOnBeta->asSexp();  

    if (Rf_isNull(rSBBOn1) || Rf_length(rSBBOn1) <= 0) {
      Rf_error("TwoLassoSexp: SBBOn1 is Bad NULL!\n"); R_FlushConsole();
    }
    if (Rf_length(rSBBOn1) != p) {
      Rf_error("TwoLassoSexp: rSBBOn1 has length %d does not equal p = %d\n",
        Rf_length(rSBBOn1), p); 
    }    
    RSBBOn1 = new AObject(rSBBOn1);
    SBBOn1 = RSBBOn1->asSexp();

    if (Rf_isNull(rSOnGammas) || Rf_length(rSOnGammas) <= 0) {
      Rf_error("TwoLassoSexp: SOnGammas is Bad NULL!\n"); R_FlushConsole();
    }
    if (Rf_length(rSOnGammas) != p) {
      Rf_error("TwoLassoSexp: SOnGammas has length %d does not equal p = %d\n",
        Rf_length(rSOnGammas), p); 
    }     
    RSOnGammas = new AObject(rSOnGammas);
    SOnGammas = RSOnGammas->asSexp();
    if (Verbose > 4) {
      Rprintf("Reading in SOnBeta = "); R_FlushConsole();
      PrintVector(REAL(rSOnBeta), (int) Rf_length(rSOnBeta)); Rprintf("\n");
      Rprintf("Reading in SBBOn1 = "); R_FlushConsole();
      PrintVector(REAL(rSBBOn1), (int) Rf_length(rSBBOn1)); Rprintf("\n");
      Rprintf("Reading in SOnGammas = "); R_FlushConsole();
      PrintVector(REAL(rSOnGammas), (int) Rf_length(rSOnGammas)); Rprintf("\n");
      R_FlushConsole();
    }
    if (Rf_length(SOnBeta) != Rf_length(SBBOn1) || Rf_length(SOnBeta) != 
      Rf_length(SOnGammas)) {
      Rf_error("SOnBeta[%d], SBBOn1[%d], SOnBeta[%d] are of wrong Rf_length",
        Rf_length(SOnBeta), Rf_length(SBBOn1), Rf_length(SOnBeta));  
    }
    
    SRecOnBeta = R_NilValue; SRecBBOn1 = R_NilValue;
    if (rRecOnBeta != NULL && !Rf_isNull(rRecOnBeta) && Rf_length(rRecOnBeta) >= p) {
      RSRecOnBeta = new AObject(rRecOnBeta);
      SRecOnBeta = RSRecOnBeta->asSexp();
    } else if (RSRecOnBeta != NULL) {
      delete(RSRecOnBeta); RSRecOnBeta = NULL;  SRecOnBeta = R_NilValue;
    }
    if (rRecBBOn1 != NULL && !Rf_isNull(rRecBBOn1) && Rf_length(rRecBBOn1) >= p) {
      RSRecBBOn1 = new AObject(rRecBBOn1);
      SRecBBOn1 = RSRecBBOn1->asSexp();
    }  else if (RSRecBBOn1 != NULL) {
      delete(RSRecBBOn1); RSRecBBOn1 = NULL;  SRecBBOn1 = R_NilValue;
    }
        
    OnLambdaA = REAL(SLambdaAK)[0]; OnLambdaD = REAL(SLambdaDK)[0];
    if (Rf_length(SLambdaAK) <= 0 || Rf_length(SLambdaAK) != Rf_length(SLambdaDK)) {
      Rf_error("SLambdaA[%d] and SLambdaD[%d] not of the right Rf_lengths!",
        Rf_length(SLambdaAK), Rf_length(SLambdaDK));
    }

    if (Verbose > 4) {
      Rprintf("TwoLassoSexp:  Loading OrderSeq\n");
    }    
    if (Rf_length(rSOrderSeq) != Rf_length(SLambdaAK)) {
      Rprintf("Note Rf_length(rSOrderSeq) = %d, but Rf_length(Slambda) = %d\n",
        Rf_length(rSOrderSeq), Rf_length(SLambdaAK));R_FlushConsole();
      Rf_error("rSOrderSeq not same Rf_length as SLambdaA\n");
    }
    OrderSeq = (int *) Calloc(Rf_length(rSOrderSeq), int);
    if (OrderSeq == NULL) {
      Rf_error("Could Not Allocate OrderSeq");
    }
    if (Rf_isReal(rSOrderSeq)) {
      for (ii = 0; ii < Rf_length(rSOrderSeq); ii++) {
        OrderSeq[ii] = (int) REAL(rSOrderSeq)[ii];
      }
    } else if (isInteger(rSOrderSeq)) {
      for (ii = 0; ii < Rf_length(rSOrderSeq); ii++) {
        OrderSeq[ii] = (int) INTEGER(rSOrderSeq)[ii];
      }      
    } else {
      Rf_error("rSOrderSeq Is of the Wrong Type!");
    }
    if (Verbose > 4) {
      Rprintf("TwoLassoSexp: Into the rSm\n"); R_FlushConsole();
      if (Rf_isNull(rSm)) {
        Rprintf("   But rSm is Null! \n");
      } else if (Rf_isReal(rSm) && Rf_length(rSm) >= 2) {
        Rprintf("  rSm[0] = %f, rSm[1] = %f \n",
          REAL(rSm)[0], REAL(rSm)[1]); R_FlushConsole();
      } else {
        Rprintf("  rSm is not good enough\n");
      }
    }
    if (!Rf_isNull(rSm) && Rf_length(rSm) == 2) {
      m1 = REAL(rSm)[0]; m2 = REAL(rSm)[1];
    } else {m1 = 0; m2 = 0;}
    
    if (Verbose > 4) {
      Rprintf("Trying to read from rSSigmaPrior\n");
      if (Rf_isNull(rSSigmaPrior)) {
        Rprintf("Note that Prior is Null!\n");
      }
      R_FlushConsole();
    }
    if (!Rf_isNull(rSSigmaPrior) && Rf_length(rSSigmaPrior) == 2 && 
      Rf_isReal(rSSigmaPrior)) {
      SigmaBar = REAL(rSSigmaPrior)[0]; SigmaDf = REAL(rSSigmaPrior)[1];
    }
    if (m1 < 0 && !(m1 < -1.4 && m1 > -1.55)) {
      FixKa = (int) fabs(m1);
    } else {
      FixKa = -1;
    }
     
     
    if (Verbose > 4) {
      Rprintf("Trying  to look into rSXtX\n");
      if (Rf_isNull(rSXtX)) {
        Rprintf("Note that rSXtX is Null!\n");
      }
      if (Rf_isNull(RGET_DIM(rSXtX))) {
        Rprintf(" Note that rSXtX has Null dimension\n");
      }
      R_FlushConsole();
    } 
    if (Verbose > 4) {
      Rprintf("Trying  to look into rSxxs\n");
      if (Rf_isNull(rSxxs)) {
        Rprintf("Note that rSxxs is Null!\n");
      }
      if (Rf_isNull(RGET_DIM(rSxxs))) {
        Rprintf(" Note that rSxxs has Null dimension\n");
      } else {
        Rprintf(" rSxxs has dimension [%d, %d], but n = %d, p = %d\n",
          INTEGER(RGET_DIM(rSxxs))[0],INTEGER(RGET_DIM(rSxxs))[1],n,p); 
        R_FlushConsole(); 
      }
      R_FlushConsole();
    } 
    int BadDim = 0;
    if (Rf_isNull(rSXtX)) { BadDim = 1; 
    } else if (!Rf_isReal(rSXtX)) {  
      BadDim = 2; 
    } else if (Rf_isNull(RGET_DIM(rSXtX))) { BadDim = 1;
    } else if (INTEGER(RGET_DIM(rSXtX))[0] != p) { BadDim = 1;
    } else if (INTEGER(RGET_DIM(rSXtX))[1] != p) { BadDim = 1;
    }    
    if (BadDim >= 1) {
      if (Verbose > 4) {
        Rprintf("We've got a BadDim for rSXtX\n"); R_FlushConsole();
        if (BadDim == 2) { Rprintf("rSXtX is not NULL and not Real!\n"); }
      }

      if ( Rf_isNull(rSxxs) || Rf_isNull(RGET_DIM(rSxxs)) || 
        INTEGER(RGET_DIM(rSxxs))[0] != n ||
        INTEGER(RGET_DIM(rSxxs))[1] != p ||
        Rf_length(rSyys) != n) {
        if (Rf_isNull(RGET_DIM(rSxxs))) {
          Rprintf("rSxxs has no dimensions, Doh!! \n");
          Rf_error(" Neither does rSXtX \n");
        }
        if (!Rf_isNull(rSXtX) && Rf_isReal(rSXtX) && !Rf_isNull(RGET_DIM(rSXtX))) {
          Rprintf("rSxxs[%d,%d] or rSyys[%d]  or rSXtX[%d,%d] or  ",
            INTEGER(RGET_DIM(rSxxs))[0], INTEGER(RGET_DIM(rSxxs))[1],
            Rf_length(rSyys), INTEGER(RGET_DIM(rSXtX))[0],
            INTEGER(RGET_DIM(rSXtX))[1]);    R_FlushConsole();

        } else {
          Rprintf("rSxxs[%d,%d] or rSyys[%d]  or rSXtX[NULL] or  ",
            INTEGER(RGET_DIM(rSxxs))[0], INTEGER(RGET_DIM(rSxxs))[1],
            Rf_length(rSyys));  R_FlushConsole();     
        }
          Rprintf("  or rSXtY[%d] are wrong dimension", Rf_length(rSXtY));
          Rf_error(" It Don't work");
      } 
      SXtX = NULL;  SXtY = NULL;  
      RSyys = new AObject(rSyys);
      Syys = RSyys->asSexp();  
      RSxxs = new AObject(rSxxs);
      Sxxs = RSxxs->asSexp();
      if (Verbose > 0) {
        Rprintf("TwoLassoSexp: Creating CDO with XX, XY;\n"); R_FlushConsole();
      }
      CDO = new CoordinateDescentObject(n, p, REAL(Sxxs), REAL(Syys),
        REAL(SOnBeta), REAL(SOnGammas)[0], 
        NULL, // Do not record CoordinateDescent OnBetas 
        REAL(SOnGammas), InitKKs, InverseGammaConstant, Verbose - 3);
    } else {
      if (Rf_isNull(rSXtX) || Rf_length(rSXtX) <= 0) {
        Rprintf("TwoLassoSexp: This ain't going to sork setting Sexp XtX");
        Rprintf("  Supplied rSXtX is NULL! \n"); R_FlushConsole();
        Rf_error("TwoLassoSexp Error! \n"); 
      }
      if (Rf_isNull(rSXtY) || Rf_length(rSXtY) <= 0) {
        Rprintf("TwoLassoSexp: This ain't going to sork setting Sexp XtY");
        Rprintf("  Supplied rSXtY is NULL! \n"); R_FlushConsole();
        Rf_error("TwoLassoSexp Error! \n"); 
      }
      RSXtX = new AObject(rSXtX);
      RSXtY = new AObject(rSXtY);  
      SXtX = RSXtX->asSexp();  SXtY = RSXtY->asSexp();
      Syys = R_NilValue; Sxxs = R_NilValue;
      if (Verbose > 0) {
        Rprintf("TwoLassoSexp: Creating CDO with XtX, XtY;\n"); R_FlushConsole();
      }
      CDO = new CoordinateDescentObject(-n, p, REAL(SXtX), REAL(SXtY),
        REAL(SOnBeta), REAL(SOnGammas)[0], NULL, 
        REAL(SOnGammas), Verbose - 3);
      if (Rf_length(Syys) != n) { Syys = NULL;  }
      if (INTEGER(RGET_DIM(Sxxs))[0] != n ||
        INTEGER(RGET_DIM(Sxxs))[1] != p) { Sxxs = NULL; } 
    }
    if (Verbose > 4) {
      Rprintf("TwoLassoSexp: Setting Epsilon = %f, Cauchy = %f\n",
        CauchyEpsilon, (double) MaxCauchy);
      R_FlushConsole();
    }
    CDO->SetEpsilonAndLoops(CauchyEpsilon, (double) MaxCauchy);
    
    CDO->FactorGammaDenominator = InverseGammaConstant;

    if (Verbose > 4) {
      Rprintf("Loading in TDFNu \n"); R_FlushConsole();
      if (Rf_isNull(rSTDFNu)) {
        Rprintf(" Note that TDFNu is NULL! \n"); R_FlushConsole();
      }
    }
    iiWeights = NULL; SiiWeights = NULL; TDFNu = 0;
    TDFNu = GetFirstDouble(rSTDFNu);  STDFNu = rSTDFNu;
    if (Verbose > 4) {
      Rprintf("TDFNu = %f", TDFNu);
      if (TDFNu > 0) {
        if (Rf_isNull(rSiiWeights)) {
          Rprintf(", but rSiiWeights is NULL!\n"); R_FlushConsole();
        } else if (Rf_length(rSiiWeights) != n) {
          Rprintf(", but Rf_length rSiiWeights = %d\n",
            Rf_length(rSiiWeights)); R_FlushConsole();
        } else {
          Rprintf(", and Rf_length rSiiWeights = %d\n", n); R_FlushConsole();
        }
      }
    }
    if (TDFNu > 0) {
      if (Rf_isNull(rSiiWeights) || Rf_length(rSiiWeights) != n) {
        iiWeights = (double*) Calloc(n, double);
        if (iiWeights == NULL) {
          Rprintf("TDFNu weights fail! \n"); R_FlushConsole();
        }
        SiiWeights = R_NilValue; RSiiWeights = NULL;
      } else {
        DDelete(RSiiWeights, "RSiiWeights");
        RSiiWeights = new AObject(rSiiWeights);
        SiiWeights = RSiiWeights->asSexp();
        iiWeights = REAL(rSiiWeights);
      }
    } else if (TDFNu <= 0) {
      if (!Rf_isNull(rSiiWeights) && Rf_length(rSiiWeights) == n) {
        DDelete(RSiiWeights, "SiiWeights");
        RSiiWeights = new AObject(rSiiWeights);
        SiiWeights = RSiiWeights->asSexp();
        iiWeights = REAL(SiiWeights);
      }
    }
    if (Verbose > 4) {
      Rprintf("TwoLassoSexp: Reading rRecBBOn1, and rRecOnBeta\n");
      R_FlushConsole();
      //Rprintf("rRecBBOn1: = "); PrintRMatrix(rRecBBOn1)
    }
    if (rRecBBOn1 != NULL && !Rf_isNull(rRecBBOn1) && Rf_isReal(rRecBBOn1) &&
       Rf_length(rRecBBOn1) >= p) {
       if (RSRecOnBeta != NULL) {
         if (Verbose >= 2) {
           Rprintf("TwoLassoSexp: Delete RSRecOnBeta to start anew \n"); 
           R_FlushConsole();
         }
         delete(RSRecOnBeta); RSRecOnBeta = NULL; 
       }
       RSRecBBOn1 = new AObject(rRecBBOn1);
       SRecBBOn1 = RSRecBBOn1->asSexp();
    } else {
      if (Verbose >= 2) {
        Rprintf("TwoLassoSexp: Setup, Deleting RSRecBBOn1 for null input.\n"); R_FlushConsole();
      }
      if (RSRecBBOn1 != NULL) {
      delete(RSRecBBOn1);   SRecBBOn1 = R_NilValue;     RSRecBBOn1 = NULL;
      }
    }
    if (rRecOnBeta != NULL && !Rf_isNull(rRecOnBeta) && Rf_length(rRecOnBeta) >= p) {
       if (RSRecOnBeta != NULL) {
         if (Verbose >= 2) {
           Rprintf("TwoLassoSexp: Setup Delete RSRecOnBeta to start anew.\n");
           R_FlushConsole();
         }
         delete(RSRecOnBeta); RSRecOnBeta = NULL; 
       }
       RSRecOnBeta = new AObject(rRecOnBeta);
       SRecOnBeta = RSRecOnBeta->asSexp();
    }  else {
      if (Verbose >= 2) {
        Rprintf("TwoLassoSexp: Setup, Deleting RSRecOnBeta \n"); R_FlushConsole();
      }
      if (RSRecOnBeta != NULL) {
       delete(RSRecOnBeta);   SRecOnBeta = R_NilValue;  RSRecOnBeta = NULL;
      }
    }
    if (Verbose > 0) {
      Rprintf("TwoLassoSexp: End Constructor!\n"); R_FlushConsole();
    } 
    if (Verbose > 3) {
      CDO->PrintFlag = Verbose -3;
    }
    if (!Rf_isNull(rSL2ShrinkagePrior) && rSL2ShrinkagePrior != NULL &&
      REAL(rSL2ShrinkagePrior)[0] >= -3) {
      RSL2ShrinkagePrior = new AObject(rSL2ShrinkagePrior);
      SL2ShrinkagePrior = RSL2ShrinkagePrior->asSexp();
      if (!Rf_isNull(rSL2ShrinkageRecords) &&
        Rf_length(rSL2ShrinkageRecords) > 2 &&
        REAL(rSL2ShrinkageRecords)[0] >= 0) {
        RSL2ShrinkageRecords = new AObject(rSL2ShrinkageRecords);
        SL2ShrinkageRecords = RSL2ShrinkageRecords->asSexp();  
      } else {
        SL2ShrinkageRecords = NULL;
      }
    } else {
      SL2ShrinkagePrior = NULL;  SL2ShrinkageRecords = NULL;
    }
    
    if (Verbose > 0) {
      Rprintf("TwoLassoSexp: Extra super End Constructor?\n");
      R_FlushConsole();
    }
    if (p > 5000 && CDO != NULL) {
      RunCoords = HeapSortAll(p, CDO->XTY);
      ReverseString(p, RunCoords);
      Rprintf("We Just Ran HeapSortAll");
      CDO->SetupRunOrder(RunCoords);
    }
    SetupSumYYSq();
  }
  int get_OnKappaS() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_OnKappaS(), CDO is not set up!\n");
    }
    return(CDO->OnKappaS);
  } 
  int get_OnKappaMem() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_OnKappaMem(), CDO is not set up!\n");
    }
    return(CDO->OnKappaMem);
  } 
  SEXP get_kXFinder() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_kXfinder(), CDO is not set up!\n");
    }
    if (CDO->OnKappaS == 0) { 
      if(Verbose > 0) {
        Rprintf("TwoLassoSexp: get_kXFinder, OnKappaS == 0\n");
        R_FlushConsole();
      }
      return(R_NilValue); 
    }
    if (CDO->OnKappaS < 0 || CDO->OnKappaS > p) {
      Rprintf("TwoLassoSexp: get_kXFinder, uh, oh, OnKappaS = %d, not good. \n", CDO->OnKappaS); 
    }
    SEXP sOn = R_NilValue; 
    Rf_protect(sOn = Rf_allocVector(INTSXP, CDO->OnKappaS));
    for (int ii = 0; ii < CDO->OnKappaS; ii++) {
      INTEGER(sOn)[ii] = CDO->kXFinder[ii];
    }
    Rf_unprotect(1);  return(sOn);
  }
  SEXP get_AllkXFinder() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_kXfinder(), CDO is not set up!\n");
    }
    if (CDO->OnKappaS == 0) { 
      if(Verbose > 0) {
        Rprintf("TwoLassoSexp: get_kXFinder, OnKappaS == 0\n");
        R_FlushConsole();
      }
      return(R_NilValue); 
    }
    if (CDO->OnKappaS < 0 || CDO->OnKappaS > p) {
      Rprintf("TwoLassoSexp: get_kXFinder, uh, oh, OnKappaS = %d, not good. \n", CDO->OnKappaS); 
    }
    SEXP sOn = R_NilValue; 
    Rf_protect(sOn = Rf_allocVector(INTSXP, CDO->OnKappaMem));
    for (int ii = 0; ii < CDO->OnKappaMem; ii++) {
      INTEGER(sOn)[ii] = CDO->kXFinder[ii];
    }
    Rf_unprotect(1);  return(sOn);
  }
  SEXP get_XTResidModify() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_XTresid(), CDO is not set up!\n");
    }
    if (CDO->XTResid == NULL) { Rf_error("TwoLassoSexp: get_XTResid, it is NULL!\n"); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, p));
    for (int ii = 0; ii < p; ii++) { REAL(sOn)[ii] = CDO->XTResid[ii]; }
    if (CDO->XTXFlag == 2) {
      for (int ii = 0; ii < p; ii++) {
        if (CDO->OnBeta[ii] != 0.0 && 
          CDO->XLC[ii] >= 0) {
          if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[CDO->XLC[ii]] <= 0) {
            CDO->ReweightCoordinate(ii);
          }
          REAL(sOn)[ii] += *(CDO->pXTX[CDO->XLC[ii]] + ii) * CDO->OnBeta[ii];
        }
      }
    } else {
      for (int ii = 0; ii < p; ii++) {
        if (CDO->OnBeta[ii] != 0.0) {
          if (CDO->iiWeights != NULL && CDO->iWeightedXtX != NULL && CDO->iWeightedXtX[ii] <= 0) {
            CDO->ReweightCoordinate(ii);
          }
          REAL(sOn)[ii] += *(CDO->pXTX[ii] + ii) * CDO->OnBeta[ii];
        }
      }    
    }
      
    Rf_unprotect(1); return(sOn);  
  }
  int TReweightCoordinate(int TryC) {
    if (CDO == NULL) {
      Rf_error("ReweightCoordiante: TwoLassoObject2014.h, there is no CDO!\n");
    }
    if (CDO->iiWeights == NULL) {
      Rf_error("ReweightCoordiante:  TwoLassoObject2014.h, there is no iiWeights!\n");
    }
    if (TryC < 0 || TryC >= p) {
      Rf_error("ReweightCoordiante:  TwoLassoObject2014.h, TryC = %d not good if p=%d!!\n", TryC, p);
    }
    if (CDO->XTXFlag >= 2 && CDO->XLC[TryC] < 0) {
      Rf_error("ReweightCoordiante:  TwoLassoObject2014.h, TryC = %d but XLC[TryC=%d]=%d not in table yet!!!\n", TryC, TryC, CDO->XLC[TryC]);    
    }
    return(CDO->ReweightCoordinate(TryC));
  }
  SEXP get_RecordBetaCV() {
    if (RSRecordBetaCV != NULL) {
      return(RSRecordBetaCV->asSexp());
    }
    return(R_NilValue);
  }
  SEXP get_RSRecordBeta0CV() {
    if (GLMCDO == NULL) {
      Rprintf("get_RSRecordBetaCV: Why if GLMCDO is Null?\n"); R_FlushConsole();
      return(R_NilValue);
    }
    if (RSRecordBeta0CV == NULL) {
      return(R_NilValue);
    }
    return(RSRecordBeta0CV->asSexp());
  }
  SEXP get_XTResid() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_XTresid(), CDO is not set up!\n");
    } 
    if (GLMCDO != NULL) {
      if (GLMCDO->XTResid == NULL) {
        Rf_error("TwoLassoSexp: XTResid in GLMCDO is not setup!\n");
      }
      return(GiveADoubleVectorOut(GLMCDO->XTResid, p, Verbose, "GLMCDO->XTResid"));
    }
    if (CDO->XTResid == NULL) { Rf_error("TwoLassoSexp: get_XTResid, it is NULL!\n"); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, p));
    for (int ii = 0; ii < p; ii++) { REAL(sOn)[ii] = CDO->XTResid[ii]; }
    Rf_unprotect(1); return(sOn);
  }
  int get_XTXFlag() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_XTX(), CDO is not set up!\n");
    }
    return(CDO->XTXFlag);
  }
  SEXP get_XLC() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_XTX(), CDO is not set up!\n");
    }
    if (CDO->XTXFlag < 2) {
      Rprintf("get_XLC, can't, because XTXFlag < 2 this procedure wasn't used!\n");
      return(R_NilValue);
    }
    if (CDO->XLC == NULL) { Rf_error("TwoLassoSexp: get_XLC, it is NULL!\n"); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(INTSXP, p));
    for (int ii = 0; ii < p; ii++) {
      INTEGER(sOn)[ii] = CDO->XLC[ii];
    }
    Rf_unprotect(1); return(sOn);
  }
  int CDORunAlgorithm() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSExp: CDORunAlgorithm, CDO not here. \n");
    }
    return(CDO->RunAlgorithm(MaxCauchy, CauchyEpsilon));
  }
  int get_MaximumAllocation() {
    return(CDO->MaximumAllocation);
  }
  int get_MemoryFailFlag() {
    return(CDO->MemoryFailFlag);
  }
  int ResetMemoryFail(int AddMem) {
    if (this->CDO->MemoryFailFlag == 0) {
      Rprintf("TwoLassoObject2014.h::ResetMemoryFail() No, MemoryFailFlag is not a fail.\n");
      R_FlushConsole();
      return(-1);
    }
    for (int kti = 0; kti < p; kti++) {
         REAL(SOnBeta)[kti] = 0.0;
    }
    this->CDO->MakeXTResid();
    
     this->CDO->MemoryFailFlag = 0;
     this->CDO->SuccessFlag = 1;
    return(1);
  }
  void set_MaximumAllocation(int sIn) {
    if (sIn < 0) {
      Rprintf("Set MaximumAllocation, no, sIn=%d!\n", sIn);
      return;
    }
    if (sIn < CDO->MaximumAllocation && sIn < CDO->OnKappaMem) {
      Rprintf("Set MaximumAllocation, no set %d, CDO->OnKappaMem is allready %d, InitKappaMem=%d\n",
        sIn, CDO->OnKappaMem, CDO->InitKappaMem);
      R_FlushConsole();
      return;
    }
    CDO->MaximumAllocation = sIn;
  }
  int setupGLMNeeds() {
     Beta0 = 0.0; Beta0Prev = 0.0; 
     if (Verbose >= 1) {
       Rprintf("Run Setup GLMNeeds. \n"); R_FlushConsole();
     }
     DDelete(RBetasPrev, "RBetasPrev");
     DDelete(RBeta0Records, "RBeta0Records");
     DDelete(RProbWeights, "RProbWeights");
     DDelete(RLogitCurrentProb, "RLogitCurrentProb");
     DDelete(RLogitPrevProb, "RLogitPrevProb");
     if (GiveCDOIter != NULL)  {Free(GiveCDOIter); }
     RBetasPrev = new AObject(Rf_allocVector(REALSXP, p));
     if ((RecordFlag == 1 || RecordFlag == 2) && MaxCauchy >= 1) {
       RBeta0Records = new AObject(Rf_allocVector(REALSXP, MaxCauchy));
     }
     RProbWeights = new AObject(Rf_allocVector(REALSXP, n));
     for (int ii = 0; ii < n; ii++) {
       REAL(RProbWeights->asSexp())[ii] = 1.0;
     }
     RLogitCurrentProb = new AObject(Rf_allocVector(REALSXP, n));
     for (int ii = 0; ii < n; ii++) {
       REAL(RLogitCurrentProb->asSexp())[ii] = 0.0;
     }
     RLogitPrevProb = new AObject(Rf_allocVector(REALSXP, n));  
     for (int ii = 0; ii < n; ii++) {
       REAL(RLogitCurrentProb->asSexp())[ii] = 0.0;
     }
     GiveCDOIter = (int*)Calloc(3, int);
     return(1);
  }
  SEXP get_nlogitCurrentProb() {
    if (GLMCDO == NULL) {Rf_error("GLMCDO Not Setup. \n"); }
    if (GLMCDO->nlogitCurrentProb == NULL) {
      Rf_error("GLMCDO, nlogitCurrentProb is NULL. \n");
    }
    SEXP sOut = R_NilValue; Rf_protect(sOut = Rf_allocVector(REALSXP, n));
    for (int ii = 0; ii < n; ii++) {
      REAL(sOut)[ii] = GLMCDO->nlogitCurrentProb[ii];
    }
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_XTX() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_XTX(), CDO is not set up!\n");
    }
    if (Verbose >= 2) {
      Rprintf("** TwoLassoObject2014.h::get_XTX, we have pXTX in CDO length \n", p); R_FlushConsole();
    }
    if (CDO->pXTX == NULL) { Rf_error("TwoLassoSexp: get_XTX, pXTX it is NULL!\n"); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocMatrix(REALSXP, p, CDO->OnKappaS));
    int One = 1;  
    for (int jj = 0; jj < CDO->OnKappaS; jj++) {
      F77_CALL(dcopy)(&p, CDO->pXTX[jj], &One, REAL(sOn) + jj * p, &One);
    }   
    Rf_unprotect(1); return(sOn);
  }
  SEXP get_XTXjj(int jj) {
     if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_XTX(), CDO is not set up!\n");
    }
    if (Verbose >= 2) {
      Rprintf("** TwoLassoObject2014.h::get_XTX, we have pXTX in CDO length \n", p); R_FlushConsole();
    }
    if (CDO->pXTX == NULL) { Rf_error("TwoLassoSexp: get_XTX, pXTX it is NULL!\n"); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, p));
    int One = 1;           
    if (CDO->XTXFlag == 2) {
      if (jj < 0 || jj >= CDO->OnKappaS) {
        Rf_error("get_XTXjj, no jj = %d, but OnKappaS = %d\n", jj, CDO->OnKappaS);
      }
    } else {
      if (jj < 0 || jj >= p) {
        Rf_error("get_XTXjj, no jj = %d, but p = %d \n", jj, p); 
      }
    }
    F77_CALL(dcopy)(&p, CDO->pXTX[jj], &One, REAL(sOn), &One);   
    Rf_unprotect(1); return(sOn);
  }
  SEXP get_XTY() {
    if (CDO == NULL) {
      Rf_error("TwoLassoSexp: get_XTX(), CDO is not set up!\n");
    }
    if (GLMCDO != NULL) {
       return(GiveADoubleVectorOut(GLMCDO->XTY, p, Verbose, "GLMCDO->XTY"));
    }
    if (CDO->XTY == NULL) { Rf_error("TwoLassoSexp: get_XTY, it is NULL!\n"); }
    SEXP sOn = R_NilValue;
    Rf_protect(sOn = Rf_allocVector(REALSXP, p));
    int One =1 ;  int Cp = p;
    F77_CALL(dcopy)(&Cp, CDO->XTY, &One, REAL(sOn), &One);
    Rf_unprotect(1); return(sOn);
  }
  TwoLassoSexp(SEXP rListFlags, SEXP rSyys, SEXP rSxxs,
    SEXP rSOnBeta, SEXP sListLambdas) {
    RSyys = NULL; Syys = R_NilValue; RSxxs = NULL; Sxxs = R_NilValue;; 
    RSXtY = NULL; SXtY = R_NilValue; 
    RSXtX = NULL; SXtX = R_NilValue; RStt1 = NULL; 
    Stt1 = R_NilValue; RStt2 = NULL; Stt2 = R_NilValue;
     TDFNu = -1;   zzs = NULL;
    RSLambdaAK = NULL; RSLambdaDK = NULL; 
    SLambdaAK = R_NilValue; SLambdaDK = R_NilValue;
    RSOrderSeq = NULL; SOrderSeq = R_NilValue;
    RSBBOn1 = NULL; SBBOn1 = R_NilValue;   SigmaBar = -1.0;  SigmaDf = -1.0;
    RSOnBeta = NULL; SOnBeta = R_NilValue; RSSigma = NULL;  SSigma=R_NilValue; 
    NonZero = 0; NonZeroG = 0.0; CurSBBar= 0.0; CountGroupSq = 0.0;
    RSRecBBOn1 = NULL; SRecBBOn1 = R_NilValue;  RBeta0Records = NULL;
    RSRecOnBeta = NULL;  SRecOnBeta = R_NilValue;
    RSOnGammas = NULL; SOnGammas = R_NilValue; jti = 0;
    RSPiA = NULL; SPiA = R_NilValue;   RecordFlag = -1;
    RSSigma = NULL; SSigma = R_NilValue;   OnSigma = 1.0;   TargetMinimumBeta = 1.0;
    RSiiWeights = NULL; SiiWeights = R_NilValue; iiWeights = NULL; 
    RSSigmaPrior = NULL; SSigmaPrior = R_NilValue;
    RsGroupsSexp = NULL; sGroupsSexp = R_NilValue; GroupLambdaEstimates = NULL;
    RSiiWeights = NULL;  SiiWeights = R_NilValue; iiWeights = NULL; RSL2ShrinkagePrior = NULL;
    RSSigma = NULL;  DidIFail = 0;   RListFlags = NULL;  ListFlags = R_NilValue;
      SL2ShrinkagePrior = R_NilValue; RunCoords = NULL;  ForceXTXFlag = -1;
      SL2ShrinkageRecords = R_NilValue; RSL2ShrinkageRecords = NULL;
      m1=-1;m2=-1; BackOrigBeta = NULL;   m3 = 0-1; m4 = -1;
    CDO = NULL; GLMCDO = NULL; FailureFlag = 0;
    CVQuitTime = -1;  NoShrinkColumns = NULL;
      
    RBetasPrev = NULL; RBeta0Records = NULL;  OnSigma = 1.0;
    RProbWeights = NULL; RLogitCurrentProb = NULL;
    RLogitPrevProb = NULL;  GiveCDOIter = NULL;
    
    InverseGammaConstant = 1.0;
    
    ConfidenceQuantiles = NULL;  ConfidenceMatrix = NULL;
    HPDMatrix=NULL;  HPDQuantiles=NULL;
    UnshrunkConfidenceMatrix = NULL; XjXj_CI_Vector=NULL;
    Sig_CI_Vector = NULL; C_CI_Vector=NULL;
    LambdaIndexForConfidenceIntervals = 0;
    
    Beta0 = 0.0; Beta0Prev = 0.0; 
    RBetasPrev = NULL; RBeta0Records = NULL;
    RProbWeights = NULL; RLogitCurrentProb = NULL;
    RLogitPrevProb = NULL;  GiveCDOIter = NULL;
    //double *rnProbWeights, double *rnlogitCurrentProb, double *rnlogitPrevProb,
    //double *rZ, 
    
    
    BackOrigBeta = NULL;  BackXTYResid = NULL;
    SPiAVectorInputs = R_NilValue;    SSigmaVectorInputs = R_NilValue;
    SLambdaAKInputs = R_NilValue;    SLambdaDKInputs = R_NilValue;
    SRecordBetaCV = R_NilValue;
    RSPiAVectorInputs = NULL;   RSSigmaVectorInputs  = NULL;
    RSLambdaAKInputs = NULL;   RSLambdaDKInputs = NULL;
    RSRecordBeta0CV = NULL;
    RSRecordBetaCV = NULL; RSStartBetaMatrix = NULL;
    SStartBetaMatrix = R_NilValue;  Stt1 = R_NilValue;
    SSigma = R_NilValue;    // ForceUpdateBBOn1 = 1;
    
    FirstRandomIndex = -1; NumGroups = -1;
    GroupLambdaEstimates = NULL; EndGroupLocations = NULL; 
    RsGroupLambdaRecord = NULL; 
    RsGroupLambdaEstimates= NULL; GroupsBBOn1 = NULL;  BackGroupsBBOn1=NULL;

    Verbose = ANINT(rListFlags, 0);
    if (Verbose > 0) {
      Rprintf("TwoLassoCpp:: Starting Cpp Open \n"); R_FlushConsole();
    }
        
    if (Verbose > 0)  {
      Rprintf("-- About to create new AObject out of ListFlags. \n");
    }
    RListFlags = new AObject(rListFlags);
    ListFlags = RListFlags->asSexp();
    if (Verbose > 0)  {
      Rprintf("-- Successfully created AObject out of ListFlags. \n");
    }
        
    if (Rf_length(ListFlags) < 11) {
      Rf_error("Sorry: ListFlags supplied is not at least 11 length, is length %d\n",
        Rf_length(ListFlags));
    }

    n = ANINT(ListFlags, 1);
    if (Rf_isReal(ListFlags)) {  OnPiA = AREAL(ListFlags,2); }
    if (!Rf_isReal(ListFlags) && !Rf_isInteger(ListFlags)) {
      Rprintf("ERRORERRROERRROERRROERRROERRROERRR\n");
      Rprintf("ERROR: TwoLassoObject2014.h, ListFlags supplied is not real or Integer!\n");
    }
    FixKa = ANINT(ListFlags, (3));
    double OnGamma1 = AREAL(ListFlags,(4)); 
    if (Rf_isInteger(ListFlags)) {
      OnSigma = 1.0;
    } else if (!R_isnancpp(AREAL(ListFlags, (5))) && R_finite(AREAL(ListFlags, (5)))) {
      OnSigma = AREAL(ListFlags,(5));
    } else {
      OnSigma = 1.0;
    }
    if (Rf_isReal(ListFlags)) {
      InitKKs = (int) REAL(ListFlags)[6];
    } else if (Rf_isInteger(ListFlags)) {
      InitKKs = (int) INTEGER(ListFlags)[6];
    }
    if (Rf_isReal(ListFlags)) {
      CauchyEpsilon = REAL(ListFlags)[8];
    } else if (Rf_isInteger(ListFlags)) { CauchyEpsilon = 0.01; }
    if (Rf_isReal(ListFlags)) {
      MaxCauchy = (int) REAL(ListFlags)[9];
    } else if (Rf_isInteger(ListFlags)) { MaxCauchy = INTEGER(ListFlags)[9]; }
    
    if (Rf_length(ListFlags) >= 12) {
      ForceXTXFlag = ANINT(ListFlags, 11);
    }
    
    if (InitKKs <= 0) {
      Rprintf("TwoLassoCpp:: Issue, ListFlags produces InitKKs = %d!  Not good!\n",
        InitKKs); R_FlushConsole();
      Rprintf("Here Is ListFlags: (");  
      for (int jj = 0; jj < Rf_length(ListFlags); jj++) {
        if (Rf_isReal(ListFlags)) {
          if (jj < Rf_length(ListFlags)-1) {
            Rprintf("%.5f, ", REAL(ListFlags)[jj]);
          } else { Rprintf("%.5f);\n", REAL(ListFlags)[jj]); }
        }  else {
          if (jj < Rf_length(ListFlags)-1) {
            Rprintf("%d, ", INTEGER(ListFlags)[jj]);
          } else { Rprintf("%d);\n", INTEGER(ListFlags)[jj]); }
        }
      }
      Rf_error("Go Inspect this error!\n");
    }          
    if (InitKKs >= MAXINITKKs) {
      InitKKs = MAXINITKKs;  // At least force COO to load data at 1000 size first.
    }                                            
    if (Rf_isNull(sListLambdas) || Rf_length(sListLambdas) != 3) {
      Rprintf("Hey, Setup TwoLassoSexp: sListLambdas should be length 3!");
      SuccessFlag = 0;
      Rf_error("TwoLassoSexp: Must give sListLambdas first!\n");
    }
      SEXP sOne = VECTOR_ELT(sListLambdas,0);
      SEXP sTwo = VECTOR_ELT(sListLambdas,1);
      SEXP sThree = VECTOR_ELT(sListLambdas,2);
      SetupLambda(sOne, sTwo, sThree);
 
    RSPiA = new AObject(Rf_allocVector(REALSXP, 1));
    SPiA = RSPiA->asSexp();
    REAL(SPiA)[0] = OnPiA;
    RSSigma = new AObject(Rf_allocVector(REALSXP, 1));
    SSigma = RSSigma->asSexp();
    REAL(SSigma)[0] = OnSigma;
        
    SuccessFlag = -1;
  
    BackXTYResid = NULL; BackOrigBeta = NULL; SXtY = NULL;
    

    RunCoords = NULL; SumYYSq = -1.0;
    //RSn = new AObject(rSn);  Sn = RSn->asSexp();
    if (n == 0) { Rf_error("TwoLassoCpp:: Error, n must be supplied nonzero \n"); }
    if (Rf_isNull(rSxxs)) { Rf_error("TwoLassoCpp:: Error, rSxxs input must be NonNull!");
    } else if (Rf_isNull(RGET_DIM(rSxxs))) { Rf_error("TwoLassoCpp:: Error rSxxs must have dim!"); 
    } else if (Rf_length(RGET_DIM(rSxxs)) != 2) { Rf_error("TwoLassoCpp:: Error, rSxxs must have 2 dim "); 
    } else if (INTEGER(RGET_DIM(rSxxs))[1]  <= 0) { Rf_error("TwoLassoCpp:: Error, rSxxs has less than 0 dim"); }
    if (Rf_isNull(rSyys)) { Rf_error("TwoLassoCpp:: Error, rSyys input must be NonNull!");
    } else if (Rf_length(rSyys) <= 0) { Rf_error("TwoLassoCpp:: Error rSyys must have length > 0!"); }    
    
    p = INTEGER(RGET_DIM(rSxxs))[1];

   
 
    if (Rf_isNull(rSOnBeta)) {
      Rf_error("Sorry: rSOnBeta input is NULL!");
    }  else if (!Rf_isReal(rSOnBeta)) {
      Rf_error("Sorry rSOnBeta, is not REAL!");
    } else if (!Rf_length(rSOnBeta) == p) {    
      Rf_error("Sorry rSOnBeta not right length!");
    }
    RSOnBeta = new AObject(rSOnBeta);    
    SOnBeta = RSOnBeta->asSexp(); 
    if (Rf_isNull(SOnBeta) || !Rf_isReal(SOnBeta)) {
      Rf_error("TwoLassoCpp: Error: SOnBeta is NULL after protection !");
    }
    
    if (Verbose > 0) {
      Rprintf("TwoLassoCpp:: Allocating within BBOn1, SOnGammas \n");
      R_FlushConsole();
    }         
    if (p <= 0){
      SuccessFlag = -1;
      Rf_error("TwoLassoCpp: Cannot run Not if p = %d \n", p);
    }
    
    RecordFlag = -1;

    if (Rf_isInteger(ListFlags)) {
      RecordFlag = INTEGER(ListFlags)[7];
    } else {
      RecordFlag = (int) REAL(ListFlags)[7];
    }
    if (Verbose >= 2 && RecordFlag >= 0) {
      Rprintf("TwoLassoCpp:: RecordFlag = %d we will Setup Record\n"); R_FlushConsole();
    }
    if (RecordFlag == 1 || RecordFlag == 2) {
      SetupRecord(RecordFlag);
    }
    
    if (Verbose >= 2) {
      Rprintf("Allocating RSBBOn1 activation indicator vector. \n"); R_FlushConsole();
    }
    SEXP aProtect = R_NilValue;
    Rf_protect(aProtect=Rf_allocVector(REALSXP,p));
    if (!Rf_isNull(aProtect)) {
      RSBBOn1 = new AObject(aProtect); 
      Rf_unprotect(1);
    } else {
      Rf_unprotect(1);
      Rf_error("Well we have an a Protect RSBBOn1 Error!");
    }
    SBBOn1 = RSBBOn1->asSexp();
    for (int iti = 0; iti < p; iti++) { REAL(SBBOn1)[0] = OnPiA; }
    
    if (Verbose >= 2) {
      Rprintf("TwoLassoSexP:: Allocating a Protect SOnGammas within Constructor\n");
      R_FlushConsole();
    }
    Rf_protect(aProtect = Rf_allocVector(REALSXP, p));
    if (!Rf_isNull(aProtect)) {
      RSOnGammas = new AObject(aProtect);
      SOnGammas = RSOnGammas->asSexp();
      if (Verbose >= 2) {
        Rprintf("TwoLassoSexp:: Successfuly Allocated SOnGammas within Constructor\n");
        R_FlushConsole();
      }
      Rf_unprotect(1);
    } else {
      Rf_unprotect(1);
      Rf_error("Well we have a protect RSOnGammas Error!\n");
    }
    for (int iti = 0; iti < p; iti++) { REAL(SOnGammas)[iti] = OnGamma1;}

    double InverseGammaConstant = 1.0;
    if (rSxxs == NULL) {
      Rf_error("Doh!: rSxxs is literally a NULL pointer!");
    } else if (Rf_isNull(rSxxs)) {
      Rf_error("Doh!: rSxxs is a NULL type SEXP");
    } else if (rSyys == NULL) {
      Rf_error("Doh!: rSyys is a NULL pointer!");
    } else if (Rf_isNull(rSyys)) {
      Rf_error("Doh!: rSyys is a NULL type SEXP");
    } else if (Rf_isNull(RGET_DIM(rSxxs))) {
      Rf_error("Doh!: rSxxs had no dimension!");
    }
    if (Verbose > 0) {
      Rprintf("TwoLassoObject2014.h::TwoLassoSexp() n = %d, p=%d, length(rSxxs) = %d, wonder what we'll do?", 
        n,p, Rf_length(rSxxs)); R_FlushConsole();
    }
    if (n == 0) {
      Rf_error("No fucking way we're dealing with n = 0!"); R_FlushConsole();
    }
      if (Rf_isNull(rSyys)) {
        Rf_error("No, Syys is NULL!");
      }
      if (Rf_isNull(rSxxs)) {
        Rf_error("No, Sxxs is NULL!");
      }
      if (Rf_isNull(RGET_DIM(rSxxs))) {
        Rf_error("No, rSxxs has no dimension!");
      }
    if (Rf_isNull(SOnBeta) || SOnBeta == NULL || !Rf_isReal(SOnBeta)) {
        Rf_error("Before XtX/X attempt error, SOnBeta is NULL!");
    }
    if (Rf_isNull(SOnGammas) || SOnGammas == NULL || !Rf_isReal(SOnGammas)) {
      Rf_error("Before XtX/X attempt error, SOnGammas is NULL");
    }             
    if (n < 0) {
      if (Verbose > 1) {
        Rprintf("TwoLassoObject2014.h::TwoLassoSexp() n = %d <0, we're doing an XtX attempt now\n", n); R_FlushConsole();
      }
      n = -n;

      if (INTEGER(RGET_DIM(rSxxs))[0] != Rf_length(rSyys) || 
        INTEGER(RGET_DIM(rSxxs))[0] != p) {
        Rf_error("TwoLassoCpp:: Error: rSxxs/rSyys should have p dim for negative n"); 
        R_FlushConsole();
      }
      if (Rf_isNull(rSxxs)) {
        Rf_error("TwoLassoCpp: Setup Error (XTX, XTY version), rSyys is NULL Can't setup! SXtY! \n");
      }
      RSXtX = new AObject(rSxxs);  SXtX = RSXtX->asSexp();
      if (Rf_isNull(rSyys)) {
        Rf_error("TwoLassoCpp: Setup Error (XTX, XTY version), rSyys is NULL Can't setup! SXtY! \n");
      }
      RSXtY = new AObject(rSyys);  SXtY = RSXtY->asSexp();
      if (Rf_isNull(SXtY) || Rf_length(SXtY) <= 0) {
        Rf_error("XtX attempt error, SXtY is null after preserve");
      }
      if (Rf_isNull(SXtX) || !Rf_isReal(SXtX) || Rf_length(SXtX) <= 0) {
        Rf_error("XtX attempt error, SXtX is null after preserve");
      }
      
        CDO = new CoordinateDescentObject(-n, p, REAL(SXtX), REAL(SXtY),
        REAL(SOnBeta), REAL(SOnGammas)[0], NULL, 
        REAL(SOnGammas), Verbose - 2);
      if (CDO == NULL || CDO->SuccessFlag < 0) {
        Rprintf("Creation of CDO (Based upon SXtX) does not work!"); R_FlushConsole();
        if (CDO != NULL) { delete(CDO); }
        Rf_error(" Leaving !");
      }
      XTXFlag = 2;
    } else {
      n = n;
      if (Verbose > 1) {
        Rprintf("\n n = %d > 0, we're doing an X, not XtX  attempt\n", n); R_FlushConsole();
      }
      XTXFlag = 1;
      if (Rf_length(rSyys) != n) {  
        Rf_error("TwoLassoCpp:: Error, rSyys[%d] does not have n=%d Rf_length \n",
          Rf_length(rSyys), n);
      }
      if (INTEGER(RGET_DIM(rSxxs))[0] != n) {
        Rf_error("TwoLassoCpp:: Error, rSxxs[%d,%d] does not have n=%d Rf_length \n",
          INTEGER(RGET_DIM(rSxxs))[0], INTEGER(RGET_DIM(rSxxs))[1], n);
      }
      if (Rf_isNull(rSxxs)) {
        Rf_error("TwoLassoCpp:: Error, rSxxs ");
      } else if (Rf_isNull(RGET_DIM(rSxxs))) {
        Rf_error("TwoLassoCpp:: Error, rSxxs is not");
      }
      if (Rf_isNull(rSyys)) {
        Rf_error("TwoLassoCpp:: rSyys is null!");
      }
      RSxxs = new AObject(rSxxs);
      Sxxs = RSxxs->asSexp();
      int ADoLogit = ANINT(ListFlags, 10);
      if (ADoLogit >= 1 && Rf_isInteger(rSyys)) {
        int IntCount = 0;
        for (int ii = 0; ii > Rf_length(rSyys); ii++) {
          if (INTEGER(rSyys)[ii] != 0 && INTEGER(rSyys)[ii] != 1 &&
            INTEGER(rSyys)[ii] != -1) {
            IntCount = 1; break;  
          }
        }
        if (IntCount == 1) {
          Rf_error("Hey: Setup TwoLassoCpp ADoLogit = %d, but you have IntCount = %d, bad!\n", 
            ADoLogit, IntCount);
          SEXP nSyys = R_NilValue;
          Rf_protect(nSyys = Rf_allocVector(REALSXP, n));
          for (int iti = 0; iti < n; iti++) {
            REAL(nSyys)[iti] = (double) INTEGER(rSyys)[iti];
          }
          RSyys = new AObject(nSyys);  Syys = RSyys->asSexp();
          zzs = NULL;
          Rf_unprotect(1);
        } else {
          zzs = (int *) Calloc(n+1, int);
          SEXP nSyys = R_NilValue;
          Rf_protect(nSyys = Rf_allocVector(REALSXP, n));
          for (int iti = 0; iti < n; iti++) {
            if (INTEGER(rSyys)[iti] == 1) {
              zzs[iti] = 1;  REAL(nSyys)[iti] = 1.0;
            }  else {
              zzs[iti] = 0;  REAL(nSyys)[iti] = -1.0;
            }
          }
          RSyys = new AObject(nSyys);  Syys = RSyys->asSexp();
          if (Verbose >= -3) {
            Rprintf("TwoLassoCpp: We are doing GLM mode because with REAL(nSyys) \n"); R_FlushConsole();
            Rprintf("zzs = "); PrintVector(zzs, n);
            Rprintf("\n    but first nSyys = ");  PrintVector(REAL(nSyys), n); R_FlushConsole();
            Rprintf("\n"); R_FlushConsole();
            Rprintf("TwoLassoCpp: Next step is to create CDO. \n"); R_FlushConsole();
          }
          Rf_unprotect(1);
        }
      } else if (ADoLogit >= 1) {
        int IntCount = 0;
        for (int ii = 0; ii > Rf_length(rSyys); ii++) {
          if (REAL(rSyys)[ii] != 0.0 && REAL(rSyys)[ii] != 1.0 &&
            REAL(rSyys)[ii] != -1.0) {
            IntCount = 1; break;  
          }
        }
        if (IntCount == 1) {
          Rf_error("Hey: Setup TwoLassoCpp ADoLogit = %d, but you have IntCount = %d, bad!\n", 
            ADoLogit, IntCount);
          RSyys = new AObject(rSyys);  Syys = RSyys->asSexp();
          zzs = NULL;
        } else {
          zzs = (int *) Calloc(n+1, int);
          SEXP nSyys = R_NilValue;
          Rf_protect(nSyys = Rf_allocVector(REALSXP, n));
          for (int iti = 0; iti < n; iti++) {
            if (REAL(rSyys)[iti] == 1.0) {
              zzs[iti] = 1;  REAL(nSyys)[iti] = 1.0;
            }  else {
              zzs[iti] = 0;  REAL(nSyys)[iti] = -1.0;
            }
          }
          RSyys = new AObject(nSyys);  Syys = RSyys->asSexp();
          if (Verbose >= -3) {
            Rprintf("TwoLassoCpp: We are doing GLM mode because with REAL(nSyys) \n"); R_FlushConsole();
            Rprintf("zzs = "); PrintVector(zzs, n);
            Rprintf("\n   but first nSyys = ");  PrintVector(REAL(nSyys), n); R_FlushConsole();
            Rprintf("\n"); R_FlushConsole();
          }
          Rf_unprotect(1);
        }      
      }  else {
        zzs = NULL;  RSyys = new AObject(rSyys);
        Syys = RSyys->asSexp();
      }
      if (Verbose > 1) {
        if (!Rf_isNull(Sxxs)) {
          Rprintf("Sxxs is not null, ");
          if (Rf_isReal(Sxxs)) {
            Rprintf(" and it is real ");
            if (Rf_length(RGET_DIM(Sxxs)) == 2) {
              Rprintf(" with dimension %d, %d\n",
                INTEGER(RGET_DIM(Sxxs))[0], INTEGER(RGET_DIM(Sxxs))[1]);R_FlushConsole();
            } else { Rf_error("Error: it has no dimensions! ");}
          } else { Rf_error("Error: Sxxs is not real"); }
        } else {  Rf_error("Error, Sxxs is a null!"); }
        if (!Rf_isNull(Syys)) {
          Rprintf("Syys is not null, ");
          if (Rf_isReal(Syys)) {
            Rprintf(" and it is real.\n ");
           } else { Rf_error("Error: Syys is not real");}
        } else {  Rf_error("Error, Syys is a null!"); }        
        if (!Rf_isNull(SOnBeta)) {
          Rprintf("SOnBeta is not null, ");
          if (Rf_isReal(SOnBeta)) {
            Rprintf(" and it is real.\n ");
           } else { Rf_error("Error: SOnBeta is not real");}
        } else {  Rf_error("Error, SOnBeta is a null!"); }    
        if (!Rf_isNull(SOnGammas)) {
          Rprintf("SOnGammas is not null, ");
          if (Rf_isReal(SOnGammas)) {
            Rprintf(" and it is real. \n");
           } else { Rf_error("Error: SOnGammas is not real");}
        } else {  Rf_error("Error, SOnGammas is a null!"); }          
      }
      if (InitKKs <= 0) {
        InitKKs = 10;
      }
      if (InitKKs > p) {
        InitKKs = p ;
      }
      if (p <= 0) {
        Rf_error("TwoLassoCpp Start: not happening if p < 0!  CoordinateDescent Fail!\n");
        return;
      }
      if (InitKKs <= 0) {
         Rf_error("TwoLassoObject2014.h::TwoLassoSexp(): Allocate, why is InitKKs = %d? \n", InitKKs); 
      }
      if (Verbose > 1){
        Rprintf("TwoLassoObject2014.h::TwoLassoSexp(): Lets see if CoordinateDescent #2 Start will work! \n"); R_FlushConsole();
      }
     // CoordinateDescentObject(int rNLen, int rkLen, double *rXXX, double *rYYY,
     //        double *rOnBeta, double rOnGamma, double *rRecordOnBetas, 
     //        double *rOnGammas, int InitKKs, double InverseGammaConstant)
     if (zzs == NULL) {
      CDO = new CoordinateDescentObject((int) n, (int) p, (double*) REAL(Sxxs), REAL(Syys),
        REAL(SOnBeta), (double) REAL(SOnGammas)[0], 
        (double*) NULL, // Do not record CoordinateDescent OnBetas 
        REAL(SOnGammas), (int) InitKKs, (double) InverseGammaConstant, Verbose - 3);
       GLMCDO = NULL;  
     } else {
      if (Verbose >= -1) {
        Rprintf("TwoLassoObject2014::TwoLassoSexp() Setup GLMNeeds. \n"); R_FlushConsole();
      }
      setupGLMNeeds();   // RunSetupGLMNeeds
      double *pBeta0Records = NULL;
      if (Verbose >= -1) {
        Rprintf("TwoLassoObject2014::TwoLassoSexp() Setup Beta0Records. \n"); R_FlushConsole();
      }
      if (RBeta0Records == NULL && (RecordFlag == 1 || RecordFlag == 2) && MaxCauchy >= 1) {
       RBeta0Records = new AObject(Rf_allocVector(REALSXP, MaxCauchy));
     }
      if (RBeta0Records != NULL) { pBeta0Records = REAL(RBeta0Records->asSexp()); 
        pBeta0Records[0] = 0.0; }
      double *RecOnBetas = NULL;
      if (RSRecOnBeta != NULL  && Rf_length(RSRecOnBeta->asSexp()) >
        p * Rf_length(RSLambdaAK->asSexp()) ) { 
          RecOnBetas = REAL(RSRecOnBeta->asSexp()); 
          RecOnBetas[0] = Rf_length(RSLambdaAK->asSexp());
        if (RBeta0Records != NULL && Rf_length(RBeta0Records->asSexp()) < ceiling(Rf_length(RSRecOnBeta->asSexp()) / (1.0*p)) +1) {
          delete(RBeta0Records);  RBeta0Records = NULL;
          RBeta0Records = new AObject(
            Rf_allocVector(REALSXP, 
              ceiling(Rf_length(RSRecOnBeta->asSexp()) / (1.0*p)) +1));
          pBeta0Records = REAL(RBeta0Records->asSexp());
          pBeta0Records[0] = 0.0;
        }
      }
        int InMaximumAllocation = -1;
        if (Verbose >= 2) {
          Rprintf("TwoLassoObject: Allocating a GLM Version: \n"); R_FlushConsole();
        }
      if (Verbose >= -1) {
        Rprintf("TwoLassoObject2014::TwoLassoSexp() Allocate GLMBinomObject(). \n"); R_FlushConsole();
      }
      GLMCDO =  new GLMBinomObject((int)n, (int) p, REAL(Sxxs), zzs, 
        REAL(SOnBeta), &Beta0, &Beta0Prev, 
        REAL(RBetasPrev->asSexp()), REAL(SOnGammas)[0], 
        RecOnBetas, pBeta0Records,
        REAL(SOnGammas), (int) InitKKs, (double) InverseGammaConstant, 
        REAL(RProbWeights->asSexp()), REAL(RLogitCurrentProb->asSexp()),
        REAL(RLogitPrevProb->asSexp()),
        REAL(Syys), Verbose - 3, CauchyEpsilon, MaxCauchy,   //Verbose -3 for the 1
        GiveCDOIter, InMaximumAllocation);
       CDO = (CoordinateDescentObject*) GLMCDO;
     }
      if (Verbose > -1) {
        Rprintf("TwoLassoObject, Cpp version CDO: We are done from construction\n");
        R_FlushConsole();
      }
      if (CDO == NULL || CDO->SuccessFlag < 0) {
        Rprintf("Creation of CDO based upon Sxxs Approach does not work!"); R_FlushConsole();
        if (CDO != NULL) { delete(CDO); }                             
        Rf_error(" Leaving !");
      }      
      if (Verbose > 1) {
        Rprintf("TwoLassoCpp, Rcpp version, setting up SetupSumYYSq. \n"); R_FlushConsole();
      }
      SetupSumYYSq();   
    }
    if (REAL(SOnGammas) != CDO->OnGammas) {
      Rprintf("UhOhUhOhUhOhUhOhUhOh"); 
      Rprintf("Note SOnGammas pointer is different from CDO->OnGammas pointer!\n");
      R_FlushConsole();
    }
    
    RStt1 = new AObject(Rf_allocVector(INTSXP,1));
    Stt1 = RStt1->asSexp();  INTEGER(Stt1)[0] = 0;  tt1 = 0;
    RStt2 = new AObject(Rf_allocVector(INTSXP,1));
    Stt2 = RStt1->asSexp();  INTEGER(Stt2)[0] = 0;  tt2 = 0;

    SuccessFlag = 1;
    if (Verbose > 1) {
      Rprintf("TwoLassoCpp, Rcpp version, Complete and Return\n"); R_FlushConsole();
    }       
    return;
  }
  void RunSetupCDO(int AIn) {
    if (AIn < 0 || AIn >= 2) {
      Rprintf("SetupCDO: Only 0 or 1 please!\n");
      AIn = 1;
    }
    if (GLMCDO != NULL) {
      SetupCDO(AIn);
    } else {
      SetupCDO(AIn);
    }
    return;
  }
  int RunGLMLasso() {
    if (GLMCDO == NULL) {
      Rf_error("Error: GLMCDO is Null cannot run GLMLasso!\n");
    }
    return(GLMCDO->RunGLMLasso());
  }
  int CDODataIntegrity() {
    if (CDO != NULL) {
      CDO->DataIntegrity(1);
    }
    return(1);
  }
  int CDODataGLMIntegrity() {
    if (GLMCDO != NULL) {
      GLMCDO->DataGLMIntegrity(1);
    }
    return(1);
  }
  void set_PrintFlagGLMB(int OnIn) {
    if (GLMCDO == NULL) {
      Rf_error("Error: Set PrintFlagGLMB: No GLMCDO is NULL!\n");
    }
    GLMCDO->PrintFlagGLMB = OnIn;
  }
  int get_PrintFlagGLMB(int OnIn) {
    if (GLMCDO == NULL) {
      Rf_error("Error: Set PrintFlagGLMB: No GLMCDO is NULL!\n");
    }
    return(GLMCDO->PrintFlagGLMB);
  }
  void SetupTDFNoise(SEXP rSTDFNu, SEXP rSiiWeights) {
    if (CDO == NULL) {
      Rf_error("SetupTDFNoise, CDO not prepared.\n");
    }
    if (RSiiWeights != NULL) {
      DDelete(RSiiWeights, "RSiiWeights"); iiWeights= NULL; SiiWeights = R_NilValue;
    } else {
      if (iiWeights == CDO->iiWeights && CDO->iiWeights != NULL) {
        if (Verbose >= 2) {
          Rprintf("SetupTDFNoise: Free iiWeights, CDO->iiWeights"); R_FlushConsole();
        }
        Free(iiWeights); iiWeights = NULL;  CDO->iiWeights = NULL;
      } else if (iiWeights == CDO->iiWeights) {
      }                                              
      if (CDO->iiWeights != NULL) {
        Free(CDO->iiWeights);  CDO->iiWeights = NULL;
      }
      if (iiWeights != NULL) {
        Free(iiWeights);  iiWeights = NULL;
      }
      SiiWeights = R_NilValue;
      RSiiWeights = NULL;
    } 
    if (rSiiWeights == NULL || Rf_isNull(rSiiWeights) || !Rf_isReal(rSiiWeights)) {
      if (Verbose >= 3) {
        Rprintf("SetupTDFNoise: Sorry, rSiiWeights not good form.\n");
      }
      rSiiWeights = R_NilValue;
    } else if (Rf_length(rSiiWeights) != n) {
      if (Verbose >= 3) {
        Rprintf("SetupTDFNoise: Sorry, rSiiWeights has length %d, not good",
          Rf_length(rSiiWeights)); R_FlushConsole();
      }
      rSiiWeights = R_NilValue;
    }
    TDFNu = GetFirstDouble(rSTDFNu);  STDFNu = rSTDFNu;
    if (Verbose > 4) {
      Rprintf("TDFNu = %f", TDFNu);
      if (TDFNu > 0) {
        if (Rf_isNull(rSiiWeights)) {
          Rprintf(", but rSiiWeights is NULL!\n"); R_FlushConsole();
        } else if (Rf_length(rSiiWeights) != n) {
          Rprintf(", but Rf_length rSiiWeights = %d\n",
            Rf_length(rSiiWeights)); R_FlushConsole();
        } else {
          Rprintf(", and Rf_length rSiiWeights = %d\n", n); R_FlushConsole();
        }
      }
    }
    if (TDFNu > 0) {
      if (Rf_isNull(rSiiWeights) || Rf_length(rSiiWeights) != n) {
        RSiiWeights = new AObject(Rf_allocVector(REALSXP,n));
        SiiWeights = RSiiWeights->asSexp();
        iiWeights = REAL(SiiWeights);
        if (iiWeights == NULL) {
          Rprintf("TDFNu weights fail! \n"); R_FlushConsole();
        }
      } else {
        RSiiWeights = new AObject(rSiiWeights);
        SiiWeights = RSiiWeights->asSexp();
        iiWeights = REAL(SiiWeights);
      }
      CDO->TDFNoise = TDFNu;
      CDO->TDFSigma = OnSigma;
      CDO->SetupWeights(iiWeights); CDO->iiWeightsHere = 0;
    } else if (TDFNu <= 0) {
      if (!Rf_isNull(rSiiWeights) && Rf_length(rSiiWeights) == n) {
        iiWeights = REAL(rSiiWeights);
        RSiiWeights = new AObject(rSiiWeights);
        SiiWeights = RSiiWeights->asSexp();
      }
    }
    return;
  }
  int DoubleMemory() {
    if (CDO == NULL) {
      Rf_error("** TwoLassoObject2014.h::DoubleMemory() No CDO!\n");
    }
    
    if (CDO->pXTX == NULL) {
      Rf_error("** TwoLassoObject2014.h::DoubleMemory() No pXTX!\n");
    }
    if (CDO->OnKappaMem >= p) {
      Rprintf("** TwoLassoObject2014.h::DoubleMemory() OnKappaMem =%d, p=%d no more\n",
        CDO->OnKappaMem, p);
    }
    return(CDO->DoubleMemory());
  }
  SEXP get_CDOBeta() {
    if (CDO == NULL) {
      Rf_error("get_CDOBeta: Error, CDO is not setup!\n");
    }
    if (CDO->OnBeta == NULL) {
      Rf_error("get_CDOBeta: Error, CDO OnBeta not setup!\n");
    }
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(REALSXP, p));
    int One = 1;
    F77_CALL(dcopy)(&p, CDO->OnBeta, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  SEXP get_Beta() {
    return(RSOnBeta->asSexp());
  }
  void set_Beta(SEXP rsOnBeta) {
    if (Rf_isNull(rsOnBeta) || Rf_length(rsOnBeta) != p) {
      Rf_error("Set Beta: Sorry, must supply real vector of length %d", p);
    }
    int One = 1;
    if (Rf_length(rsOnBeta) != p) {
      Rf_error("Set Beta: sorry, must supply vector of length p=%d, not %d!\n",
        p, Rf_length(rsOnBeta));
    }
    if (!Rf_isReal(rsOnBeta) && Rf_isInteger(rsOnBeta)) {
      for (int ii = 0; ii < Rf_length(rsOnBeta); ii++) {
        REAL(RSOnBeta->asSexp())[ii] = REAL(rsOnBeta)[ii];
      }
    } else {
      F77_CALL(dcopy)(&p, REAL(rsOnBeta), &One, REAL(RSOnBeta->asSexp()), &One);
    }
    if (REAL(RSOnBeta->asSexp()) != CDO->OnBeta) {
      F77_CALL(dcopy)(&p, REAL(RSOnBeta->asSexp()), &One, CDO->OnBeta, &One);
    }
    CDO->MakeXTResid();
    UpdateBBOn1();
  }
  SEXP get_InvSXX() {
    if (CDO == NULL) {
      Rf_error("Can't get InvSXX if CDO is NULL!\n");
    }
    if (CDO->InvSXX == NULL) {
      Rf_error("Can't get InvSXX if InvSXX is NULL!\n");
    }
    SEXP sOut = R_NilValue; Rf_protect(sOut = Rf_allocVector(REALSXP, p));
    int One = 1; F77_CALL(dcopy)(&p, CDO->InvSXX, &One, REAL(sOut), &One);
    Rf_unprotect(1); return(sOut);
  }
  double get_DiagShrinkage() {
    if (CDO == NULL) {
      Rf_error("Can't get DiagShrinkage if CDO is NULL!\n");
    }
    return(CDO->DiagShrinkage);
  }
  void set_DiagShrinkage(double iDiagShrinkage) {
    if (CDO == NULL) {
      Rf_error("Can't get DiagShrinkage if CDO is NULL!\n");
    }
    if (iDiagShrinkage >= 0.0) {
      CDO->DiagShrinkage = iDiagShrinkage;
    }
    return;
  }
  
  double get_TotMove() {
    if (CDO == NULL) {
      Rf_error("Can't get TotMove if CDO is NULL!\n");
    }
    return(CDO->TotMove);
  }
  double get_OnEp() {
    if (CDO == NULL) {
      Rf_error("Can't get OnEp if CDO is NULL!\n"); R_FlushConsole();
    }
    return(CDO->OnEp);
  }
  double ShowSFunction() {
    if (CDO==NULL) {
      Rf_error("Can't ShowSFunction if CDO is NULL!\n");
    }
    double AReturn = CDO->ShowSFunction(CDO->OnCoord);
    return(AReturn);
  }
  int run_UpDateCoord() {
    if (CDO == NULL) {
      Rf_error("Can't Update Coord if CDO is NULL!\n");
    }
    if (Verbose >= 3) {
      Rprintf("About to manually run UpdateCoord on Coordinate %d. \n", CDO->OnCoord);
    }
    int Out = 0;
    if (CDO->OnCoord < 0 || CDO->OnCoord >= p) {
      Rf_error("Can't update Coord since OnCoord = %d! and p=%d!\n",
       CDO->OnCoord, p);
    }
    Out = CDO->UpDateCoord();
    if (REAL(RSOnBeta->asSexp()) != CDO->OnBeta) {
      REAL(RSOnBeta->asSexp())[CDO->OnCoord] = CDO->OnBeta[CDO->OnCoord];
    }
    if (Verbose >= 3) {
      Rprintf("We just ran Update Coord on Coordinate %d.  Value is %f \n",
        CDO->OnCoord, *(CDO->OnBeta + CDO->OnCoord)); R_FlushConsole();
    }
    return(Out);
  }
  void set_OnCoord(int iOnCoord) {
    if (CDO == NULL) {
      Rf_error("Can't get CDO if OnCoord is NULL.\n");
    }
    if (iOnCoord < 0 || iOnCoord >= p) {
      Rf_error(" Can't set OnCoord to %d, if p = %d\n", iOnCoord, p);
    }
    CDO->OnCoord = iOnCoord; 
  }
  int get_OnCoord() { if (CDO == NULL) {
      Rf_error("Can't get OnCoord if CDO is NULL!\n"); 
    }
    return(CDO->OnCoord);
  }
  void set_OnLoop(int iOnLoop) {
    if (CDO == NULL) {
      Rf_error("Can't get CDO if OnCoord is NULL.\n");
    }
    if (iOnLoop < 0 || iOnLoop >= MaxCauchy) {
      Rf_error(" Can't set OnLoop to %d, if MaxCauchy = %d\n", iOnLoop, MaxCauchy);
    }
    CDO->OnLoop = iOnLoop; 
  }
  int get_OnLoop() { if (CDO == NULL) {
      Rf_error("Can't get OnCoord if CDO is NULL!\n"); 
    }
    return(CDO->OnLoop);
  }
  int MakeXTResid() {
    if (CDO == NULL) {
      Rf_error("MakeXTResid: CDO must be set up. \n");
    }
    return(CDO->MakeXTResid());
  }
  void SetupPiSigma(SEXP rSPiA, SEXP rSSigma) {
    DDelete(RSPiA, "RSPiA"); DDelete(RSSigma, "RSSigma");
    RSPiA = new AObject(rSPiA);   
    SPiA = RSPiA->asSexp();  OnPiA = REAL(SPiA)[0];
    RSSigma = new AObject(rSSigma);  
    SSigma = RSSigma->asSexp(); OnSigma = REAL(SSigma)[0];
  }
  void set_PiA(SEXP rSPiA) {
    DDelete(RSPiA, "RSPiA"); 
    RSPiA = new AObject(rSPiA);   
    SPiA = RSPiA->asSexp();  OnPiA = REAL(SPiA)[0];  
  }
  double get_OnPiA() {
    return(OnPiA);
  }
  double get_OnRPiA() { return(OnRPiA); }
  void set_OnRPiA(double iOnRPiA) {
     OnRPiA = iOnRPiA;
  }
  void set_OnPiA(double iOnPiA) {
    if (iOnPiA <= 0.0 || iOnPiA >= 1.0) {
      Rprintf("Error: cannot set OnPiA outside of 0-1: %f \n",
        iOnPiA); R_FlushConsole();     return;
    }
    if (!Rf_isNull(SPiA) && Rf_length(SPiA) >= 1) {
      if (tt1 < 0) {
        REAL(SPiA)[0] = iOnPiA;
      } else if (tt1 >= Rf_length(SPiA)) {
        REAL(SPiA)[Rf_length(SPiA) - 1] = iOnPiA;
      } else {
        REAL(SPiA)[tt1] = iOnPiA;
      }
    }
    OnPiA = iOnPiA;  return;
  }
  int get_InitKKs() {
    return(this->InitKKs);
  }
  int get_InitKappaMem() {
    if (CDO == NULL) {
      Rf_error("get_InitKappaMem: CDO is NULL!\n");
    }
    return(CDO->InitKappaMem);
  }
  SEXP get_PiA() {
    SEXP sOn = R_NilValue;
    if (!Rf_isNull(SPiA) && Rf_length(SPiA) >= 1) {
      Rf_protect(sOn = Rf_allocVector(REALSXP, Rf_length(SPiA)));
      int One = 1; int Len = Rf_length(SPiA);
      F77_CALL(dcopy)(&Len, REAL(SPiA), &One, REAL(sOn), &One);
      Rf_unprotect(1); return(sOn);
    }
    Rf_protect(sOn = Rf_allocVector(REALSXP, 1));
    REAL(sOn)[0] = REAL(SPiA)[0]; Rf_unprotect(1); return(sOn);
  }
  void set_Sigma(SEXP rSSigma) {
    DDelete(RSSigma, "RSSigma");
    RSSigma = new AObject(rSSigma);  
    SSigma = RSSigma->asSexp(); OnSigma = REAL(SSigma)[0];
  }
  double get_OnLambda2() {
    if (CDO == NULL) { Rprintf("No, get_OnLambda2, CDO is NULL \n"); return(-666.0);}
    return(CDO->DiagShrinkage);
  }
  double get_OnSigma() { 
    if (R_isnancpp(OnSigma)) {
      Rprintf("ERRORERRORERRORERROR: OnSigma is NOW NAN!\n"); R_FlushConsole();
      return(-666);
    }
    if (!R_finite(OnSigma)) {
      Rprintf("ERRORERRORERROR: Onsigma is not Finite!\n"); R_FlushConsole();
    }
    return(OnSigma); 
  }
  void set_OnLambda2(double iIn) { 
     if (CDO == NULL) { Rf_error("No, get_OnLambda2, CDO is NULL \n");}
    if (iIn < 0.0) {
      Rf_error("set_OnLambda2, can't set less than 0.0\n");
    }
    double OldDiagShrinkage = CDO->DiagShrinkage;
    CDO->DiagShrinkage = iIn;
    if (CDO->InvSXX != NULL) {
    for (int ii = 0; ii < p; ii++) {
      CDO->InvSXX[ii] = 1.0/(1.0/CDO->InvSXX[ii] + CDO->DiagShrinkage - OldDiagShrinkage); 
    }
    }
  }
  void set_OnSigma(double iOnSigma) {
    if (R_isnancpp(iOnSigma)) {
      Rprintf("ERRORERRORERROR:  set_OnSigma: You supplied nan iOnSigma!\n"); R_FlushConsole();
    }
    if (iOnSigma <= 0.0) {
      Rprintf("Error: Can't set OnSigma to negative number: %f!\n", iOnSigma); R_FlushConsole();
    }
    if (!Rf_isNull(SSigma)) {
      if (tt1 < 0 || tt1 >= Rf_length(SSigma)) {
        if (Rf_length(SSigma) == 1) {
          REAL(SSigma)[0] = iOnSigma;
        } else if (Rf_length(SSigma) > 1) {
          REAL(SSigma)[Rf_length(SSigma)-1] = iOnSigma;
        }
      } else {
        REAL(SSigma)[tt1] = iOnSigma;
      }
    }
    OnSigma = iOnSigma;
  }
  int get_n() { return(n); }
  int get_p() { return(p); }
  SEXP get_Sigma() {
    SEXP sOn = R_NilValue;
    if (!Rf_isNull(SSigma) && Rf_length(SSigma) >= 1) {
      if (!Rf_isReal(SSigma)) {
        Rprintf("ERROR SSigma Is Not REAL!\n"); R_FlushConsole();
      }
      if (R_isnancpp(REAL(SSigma)[0])) {
        Rprintf("Warning SSigma, the first element is NAN!"); R_FlushConsole();
      }
      if (!R_finite(REAL(SSigma)[0])) {
        Rprintf("Warning SSigma, The First element is not finite!\n"); R_FlushConsole();
      }
      Rf_protect(sOn = Rf_allocVector(REALSXP, Rf_length(SSigma)));
      for (int ii = 0; ii < Rf_length(SSigma); ii++) {
        REAL(sOn)[ii] = REAL(SSigma)[ii];
      }
      Rf_unprotect(1); return(sOn);
    }
    Rf_protect(sOn = Rf_allocVector(REALSXP, 1));
    REAL(sOn)[0] = REAL(SSigma)[0]; Rf_unprotect(1); return(sOn);
  }
  SEXP get_RealSSigma() {
    if (SSigma == NULL || Rf_isNull(SSigma)) {
      return(R_NilValue);
    }
    return(SSigma);
  }
  void set_RealSSigma(SEXP iSSigma) {
    if (Rf_isNull(iSSigma) || Rf_length(iSSigma) <= 0 ||
      !Rf_isReal(iSSigma)) {
      Rprintf("set_RealSSigma: Nope, iSSigma is not qualified!\n"); R_FlushConsole();
    }
    if (RSSigma != NULL) {
      Rprintf("Deleting Current RSSigma\n"); R_FlushConsole();
      DDelete(RSSigma, "RSSigma");
    }
    RSSigma = new AObject(iSSigma);
    SSigma = RSSigma->asSexp();
  }
  double get_OnLambdaA() { return(OnLambdaA);}
  double get_OnLambdaD() { return(OnLambdaD);}
  SEXP get_BBOn1() {
    if (!Rf_isNull(SBBOn1)) { return(SBBOn1); }
    return(R_NilValue);
    //if (BBOn1 != NULL) {
    //  SEXP sOn = R_NilValue;  Rf_protect(sOn = Rf_allocVector(REALSXP, p));
    //  int One = 1;
    //  F77_CALL(dcopy)(&p, BBOn1, &One, REAL(sOn), &One);  
    //  Rf_unprotect(1); return(sOn);
    //} else {
    //  return(R_NilValue);
    //}
  }
  void SetupL2Shrinkage(SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords) {
      if (!Rf_isNull(rSL2ShrinkagePrior) && rSL2ShrinkagePrior != NULL &&
      REAL(rSL2ShrinkagePrior)[0] >= -3) {
      RSL2ShrinkagePrior = new AObject(rSL2ShrinkagePrior);
      SL2ShrinkagePrior = RSL2ShrinkagePrior->asSexp();
      if (!Rf_isNull(rSL2ShrinkageRecords) &&
        Rf_length(rSL2ShrinkageRecords) > 2 &&
        REAL(rSL2ShrinkageRecords)[0] >= 0) {
        RSL2ShrinkageRecords = new AObject(rSL2ShrinkageRecords);
        SL2ShrinkageRecords = RSL2ShrinkageRecords->asSexp();  
      } else {
        SL2ShrinkageRecords = NULL;
      }
    } else {
      SL2ShrinkagePrior = NULL;  SL2ShrinkageRecords = NULL;
    }
  }
  SEXP get_SL2ShrinkagePrior() {
    return(SL2ShrinkagePrior);
  }
  SEXP get_SL2ShrinkageRecords() {
    return(SL2ShrinkageRecords);
  }
  SEXP get_ListFlags() { return(ListFlags); }
  SEXP get_sRecBBOn1() {
    if (RSRecBBOn1 == NULL) { return(R_NilValue); }
    return(RSRecBBOn1->asSexp());
  }
  SEXP get_sRecOnBeta() {
    if (RSRecOnBeta == NULL) { return(R_NilValue); }
    return(RSRecOnBeta->asSexp());
  }
  

  SEXP get_SRecBBOn1() { return(SRecBBOn1); }
  SEXP get_SRecOnBeta() { return(SRecOnBeta); }
  SEXP get_y01() {
    if (GLMCDO == NULL) {
      Rf_error("GLMCDO No y01\n");
    }
    return(GLMCDO->get_y01());
  }
  int RefreshGLMWeights() {
    if (GLMCDO == NULL) {
      Rf_error("GLMCDO, Can't update Weights no GLMCDO\n");
    }
    Rprintf("Refreshing Manually GLMWeights. \n"); R_FlushConsole();
    return(GLMCDO->RefreshGLMWeights());
  }
  int RefreshCurrentLogitProb() {
    if (GLMCDO == NULL) {
      Rf_error("GLMCDO, Can't update Weights no GLMCDO\n");
    }
    Rprintf("Refreshing Manually RefreshCurrentLogitProb. \n"); R_FlushConsole();
    return(GLMCDO->RefreshCurrentLogitProb());
  }

~TwoLassoSexp() {
   //Verbose = 4;
   Rprintf("TwoLassoSexp(): Destructor begin\n"); R_FlushConsole();
   if (Verbose > 2) {
     Rprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
     Rprintf("~~~   TwoLassoSexp() Running Destructor\n"); R_FlushConsole();
     Rprintf("~~~    Should be Fun \n"); R_FlushConsole();
   }
   if (SXtX == NULL) {
     CDO->YY = NULL;  CDO->XX = NULL;
   } else {
   
   }
   DDelete(RSXtX, "RSXtX");
   DDelete(RSXtY, "RSXtY");
   DDelete(RSxxs, "RSxxs");
   DDelete(RSyys, "RSyys");
   if (zzs != NULL) {
     Free(zzs); zzs = NULL;
   }
   DDelete(RSBBOn1, "RSBBOn1");
   DDelete(RSOnBeta, "RSOnBeta");
   if (Verbose >= 1) {
     Rprintf(" ~~~TwoLassoSexp() Destructor, cleared RSBBOn1, RSOnBeta.\n");
     R_FlushConsole();
   }
   DDelete(RStt2, "RStt2"); DDelete(RStt1, "RStt1");
      
   DDelete(RSLambdaAK, "RSLambdaAK"); SLambdaAK = R_NilValue;
   DDelete(RSLambdaDK, "RSLambdaDK"); SLambdaDK = R_NilValue;
   if (Verbose >= 1) {
     Rprintf(" ~~~TwoLassoSexp() Destructor, Cleared: RSLambdaAK RSLambdaDK.\n");
     R_FlushConsole();
   }
   if (RSOrderSeq != NULL) {
     DDelete(RSOrderSeq, "RSOrderSeq"); SOrderSeq = R_NilValue;
     OrderSeq = NULL;
   } else if (!Rf_isNull(SOrderSeq)) {
      OrderSeq = NULL; SOrderSeq = R_NilValue;
   } else {
     if (Verbose > 0) { Rprintf("Deleting OrderSeq alone \n"); R_FlushConsole();}
     Free(OrderSeq);  OrderSeq = NULL;
   }
   if (Verbose >= 1) {
     Rprintf(" ~~~TwoLassoSexp() Destructor, Cleared: OrderSeq.\n");
     R_FlushConsole();
   }
   DDelete(RSL2ShrinkageRecords, "RSL2ShrinkageRecords");
   DDelete(RSL2ShrinkagePrior, "RSL2ShrinkagePrior");
   if (Verbose >= 1) {
     Rprintf(" ~~~TwoLassoSexp() Destructor, Cleared: L2Shrinkage.\n");
     R_FlushConsole();
   }
   DDelete(RSRecOnBeta, "RSRecOnBeta"); SRecOnBeta = R_NilValue;
   if (Verbose >= 1) {
     Rprintf(" ~~~ TwoLassoSexp() Destructor: about to clear RecBBOn1.\n"); R_FlushConsole();
     if (RSRecBBOn1 == NULL) {
       Rprintf(" ~~~ Note RSRecBBOn1 is null already, shouldn't be a problem.\n"); R_FlushConsole();
     } else {
       Rprintf(" ~~~ Note RSRecBBOn1 is Not Null!\n "); R_FlushConsole();
     }
   }
   DDelete(RSRecBBOn1, "RSRecBBOn1");  SRecBBOn1 = R_NilValue;
   if (Verbose >= 1) {
     Rprintf(" ~~~ Note We cleared RSRecBBOn1 \n"); R_FlushConsole();
   }
   if (Verbose >= 1) {
     Rprintf("What should we do with RSii Weights? \n"); R_FlushConsole();
     if (RSiiWeights == NULL) {
       Rprintf("Apparently RSii weights is NULL pointer \n"); R_FlushConsole();
     }
     if (RSiiWeights != NULL && Rf_isNull(SiiWeights)) {
       Rprintf("RsiiWeights is not null, but SiiWeights is, awkward? \n"); R_FlushConsole();
     }
     if (RSiiWeights == NULL && iiWeights == NULL) {
       Rprintf("iiWeights is NULL, RsiiWeights is null, we will do nothing about weights!\n"); R_FlushConsole();
     }
   }
   if (Verbose >= 1) {
     if (RSiiWeights != NULL) {
       Rprintf(" ~~~ Now we need to clear RSiiWeights. \n"); R_FlushConsole();
       if (RSiiWeights != NULL && (SiiWeights == NULL || Rf_isNull(SiiWeights)) ) {
          Rprintf(" ~~~ How could this be, Sii Weights Nil but RSiiWeights is not NULL?\n");
          R_FlushConsole();      
       }
     } else if (SiiWeights != NULL && !Rf_isNull(SiiWeights)) {
       Rprintf(" ~~~ We'll be clearing SiiWeights.\n");
       R_FlushConsole();
     } else { 
       if (iiWeights != NULL) {
         Rprintf(" ~~~ Free iiWeights manually. \n"); R_FlushConsole();
       } else {
         Rprintf(" ~~~ iiWeights is Null so we'll be doing nothing in weights. \n"); R_FlushConsole();
         iiWeights = NULL;
       }
     }
   }
   if (Verbose >= 2) {
     Rprintf(" ~~~ What fun things to do to iiWeights ?\n"); R_FlushConsole();
   }
   if (RSiiWeights != NULL) {
     if (Verbose >= 1) {
       Rprintf(" ~~~ We are truly deleting RSiiWeights \n"); R_FlushConsole();
     }
     DDelete(RSiiWeights, "RSiiWeights"); SiiWeights = R_NilValue;
     iiWeights = NULL;  
     if (CDO != NULL && CDO->iiWeightsHere >= 1 && CDO->iiWeights != NULL) {
       Free(CDO->iiWeights); CDO->iiWeights = NULL; CDO->iiWeightsHere = 0;
     }  else if (CDO != NULL && CDO->iiWeights != NULL && CDO->iiWeightsHere <= 0) {
       CDO->iiWeights = NULL; 
     }
   } else if (SiiWeights != NULL && !Rf_isNull(SiiWeights)) {
     if (Verbose >= 1) {
       Rprintf(" ~~~ Looks like we're doing something to SiiWeights.\n"); R_FlushConsole();
     }
     SiiWeights = R_NilValue;  iiWeights = NULL;
     if (Verbose >= 1) {
       Rprintf(" ~~~ Done doing stuff to iiWeights.\n"); R_FlushConsole();
     }
   } else {
     if (Verbose > 0) {
       Rprintf(" ~~~ Choosing to Delete iiWeights\n"); R_FlushConsole();}
     if (iiWeights != NULL) {
       Free(iiWeights);  iiWeights = NULL;
       if (Verbose >= 1) {
         Rprintf(" ~~~ Succesfully freed iiWeights. \n"); R_FlushConsole();
       }
     } else {
      if (Verbose >= 1) {
        Rprintf(" ~~~ We did nothing to iiWeights\n"); R_FlushConsole();
      }
     }
   }
   
   if (Verbose >= 1) {
     Rprintf(" ~~~ Attempting to delete GLM Vectors! \n"); R_FlushConsole();
   }
   
   DDelete(RBetasPrev, "RBetasPrev"); 
   DDelete(RBeta0Records, "RBeta0Records");
   DDelete(RProbWeights, "RProbWeights");
   DDelete(RLogitCurrentProb, "RLogitCurrentProb");
   DDelete(RLogitPrevProb, "RLogitPrevProb");
   DDelete(NoShrinkColumns, "NoShrinkColumns");

  if (GiveCDOIter != NULL) {
    if (Verbose >= 3) {
       Rprintf("  ~~~ Deleting GiveCDOIter \n"); R_FlushConsole();
    }
    Free(GiveCDOIter);  GiveCDOIter = NULL;
  }
   
  if (Verbose >= 1) {
    Rprintf(" ~~~ Thank you, something has been done to iiWeights. \n"); R_FlushConsole();
  }
  if (Verbose >= 2) {     
    Rprintf(" ~~~ Delete RSSigma? \n"); R_FlushConsole();
  }
  DDelete(RSSigma, "RSSigma"); SSigma = R_NilValue;
  if (Verbose >= 3) {
    Rprintf(" ~~~ Delete RsGroupsSexp? \n"); R_FlushConsole();
  }
  DDelete(RsGroupsSexp, "RsGroupsSexp");
  if (RSOnGammas != NULL) {
    DDelete(RSOnGammas, "RSOnGammas"); SOnGammas = R_NilValue;
  } else if (CDO->OnGammas != NULL) {
    if (CDO->OnGammas) {
      Rprintf("CDO->OnGammas is not Null, going to delete !\n"); R_FlushConsole();
    }
    Free(CDO->OnGammas);  CDO->OnGammas = NULL;
  }
  
  DDelete(RSPiA, "RSPiA");
  if (Verbose >= 3) {
    Rprintf("Okay, now to delete RListFlags. \n"); R_FlushConsole();
  }
  DDelete(RListFlags, "RListFlags");
   if (Verbose >= 1) {
     Rprintf(" ~~~TwoLassoSexp() Destructor, Cleared: RSPiA and RListFlags.\n");
     R_FlushConsole();
   }
  if (Verbose >= 1) {
    Rprintf(" ~~~ TwoLassoSexp: now to delete RSPiAVectorInputs \n"); R_FlushConsole();
  } 
  DDelete( RSPiAVectorInputs, "RSPiAVectorInputs");  
  DDelete( RSSigmaVectorInputs, "RSSigmaVectorInputs");
  DDelete( RSLambdaAKInputs, "RSLambdaAKInputs");  
  DDelete(  RSLambdaDKInputs, "RSLambdaDKInputs");
  DDelete(RSStartBetaMatrix, "RSStartBetaMatrix");
  if (RSRecordBetaCV == NULL) {
    if (Verbose > 1) {
      Rprintf("RSRecordBetaCV  is NULL, not deleting\n"); R_FlushConsole();
    }
  } else {
    DDelete(RSRecordBetaCV, "RSRecordBetaCV");
  } 
  if (RSRecordBeta0CV == NULL) {
    if (Verbose > 1) {
      Rprintf("RSRecordBeta0CV  is NULL, not deleting\n"); R_FlushConsole();
    }
  } else {
    DDelete(RSRecordBeta0CV, "RSRecordBeta0CV");
  } 
  
  if (Verbose >= 3) {
    Rprintf("Deleting ConfidenceQuantiles and Confidence Matrix!\n");
    R_FlushConsole();
  } 
  DDelete(ConfidenceQuantiles, "ConfidenceQuantiles");  
  if (Verbose >= 2) {
    Rprintf("Deleting HPDMatrix \n"); R_FlushConsole();
  }
  DDelete(HPDMatrix, "HPDMatrix");  DDelete(HPDQuantiles, "HPDQuantiles");
  if (Verbose >= 2) {
    Rprintf("Deleting: Confidence Matrix!\n");
    R_FlushConsole();
  } 
  DDelete(ConfidenceMatrix, "ConfidenceMatrix");
  if (Verbose >= 2) {
    Rprintf("Deleting UnshrunkConfidenceMatrix!\n");
    R_FlushConsole();
  } 
  DDelete(UnshrunkConfidenceMatrix, "UnshrunkConfidenceMatrix");
  if (Verbose >= 1) {
    Rprintf("Deleting XjXj_CI_Vector!\n");
    R_FlushConsole();
  } 
  DDelete(XjXj_CI_Vector, "XjXj_CI_Vector");
  DDelete(C_CI_Vector, "C_CI_Vector");
  DDelete(Sig_CI_Vector, "Sig_CI_Vector");
  if (Verbose >= 3) {
    Rprintf("Successfully deleted Confidence Quantiles\n"); R_FlushConsole();
  }
  LambdaIndexForConfidenceIntervals = 0;
  // DDelete(RStt1, "RStt1"); DDelete(RStt2, "RStt2");
  // DDelete(RSn, "RSn");  DDelete(RSp, "RSp");
   if (RunCoords != NULL) { Free(RunCoords); RunCoords = NULL; }
   
   if (Verbose >= 3) {
     Rprintf("Deleting Group Parameter and Estimates\n"); R_FlushConsole();
   }
   
   if (Verbose >= 3) { 
     Rprintf("Deleting GroupLambdaEstimates\n"); R_FlushConsole();
   }
   if (GroupLambdaEstimates != NULL) {
     if (Verbose >= 3) {
       Rprintf("GroupLambdaEstimates is not null, deleting!\n"); R_FlushConsole();
     }
     Free(GroupLambdaEstimates); GroupLambdaEstimates = NULL;
   }
   if (Verbose >= 3) { 
     Rprintf("Deleting EndGroupLocations\n"); R_FlushConsole();
   }
   if (EndGroupLocations != NULL) {
     Free(EndGroupLocations); EndGroupLocations = NULL;
   }
   if (Verbose >= 3) { 
     Rprintf("Deleting RsGroupLambdaEstimates\n"); R_FlushConsole();
   }
   if (RsGroupLambdaEstimates != NULL) {
     DDelete(RsGroupLambdaEstimates, "RsGroupLambdaEstimates"); RsGroupLambdaEstimates = NULL;
   }
   if (Verbose >= 3) { 
     Rprintf("Deleting RsGroupLambdaRecord\n"); R_FlushConsole();
   }
   if (RsGroupLambdaRecord != NULL) {
     DDelete(RsGroupLambdaRecord, "RsGroupLambdaRecord"); RsGroupLambdaRecord = NULL;
   }
   if (GroupsBBOn1 != NULL) {
     if (Verbose >= 3) {
       Rprintf("Deleting GroupsBBOn1\n"); R_FlushConsole();
     }
     DDelete(GroupsBBOn1, "GroupsBBOn1"); GroupsBBOn1 = NULL;
   }
   if (BackGroupsBBOn1 != NULL) {
     DDelete(BackGroupsBBOn1, "BackGroupsBBOn1"); BackGroupsBBOn1 = NULL;
   }
   
   if (Verbose >= 1) {
     Rprintf(" On to Delete CDO members and weights!\n"); R_FlushConsole();
   }
   SOnGammas = R_NilValue;  SOnBeta = R_NilValue; SBBOn1 = R_NilValue;
   CDO->OnGammas = NULL; CDO->OnBeta = NULL;
   if (SiiWeights == NULL && iiWeights != NULL) {
     Free(iiWeights);  iiWeights = NULL; 
   }
   SiiWeights = NULL;
   if (Verbose > 2) {
     Rprintf("~TwoLassoSexp() Finished Destructor\n"); R_FlushConsole();
   }
   if (BackXTYResid != NULL) {Free(BackXTYResid); BackXTYResid = NULL; }
   if (BackOrigBeta != NULL) {Free(BackOrigBeta); BackOrigBeta = NULL; }
    if (GLMCDO != NULL) {
       if (Verbose >= -1) {
         Rprintf("Deleting GLMCDO!\n"); R_FlushConsole(); 
       }
         delete(GLMCDO); GLMCDO = NULL; CDO = NULL;
    } else {
      if (CDO != NULL) {
        if (Verbose >= 1) {
           Rprintf("Deleting CDO \n"); R_FlushConsole(); 
        }
        delete(CDO); GLMCDO = NULL; CDO = NULL;
      } 
    }
  }
  SEXP get_CNoShrinkColumns() {
    if (NoShrinkColumns == NULL) {
      return(R_NilValue);
    }
    return(NoShrinkColumns->asSexp());
  }
  SEXP get_NoShrinkColumns() {
    if (NoShrinkColumns == NULL) {
      return(R_NilValue);
    }
    SEXP sOut = R_NilValue;
    Rf_protect(sOut = Rf_allocVector(INTSXP, Rf_length(NoShrinkColumns->asSexp())));
    for (int iti = 0; iti < Rf_length(sOut); iti++) {
      INTEGER(sOut)[iti] = INTEGER(NoShrinkColumns->asSexp())[iti]+1;
    }
    Rf_unprotect(1);
    return(sOut);
  }
  void set_NoShrinkColumns(SEXP iIn) {
    if (Verbose >= -1) {
      Rprintf("set_NoShrinkColumns: Running \n"); R_FlushConsole();
    }
    if (iIn == NULL || Rf_isNull(iIn) || Rf_length(iIn) <= 0) {
      if (NoShrinkColumns != NULL) {
        Rprintf("Deleting Current No ShrinkColumns\n"); R_FlushConsole();
        delete(NoShrinkColumns); NoShrinkColumns = NULL;
      }
      return;
    }
    NoShrinkColumns = new AObject(Rf_allocVector(INTSXP, Rf_length(iIn)));
    for (int iti = 0; iti < Rf_length(NoShrinkColumns->asSexp()); iti++) {
      INTEGER(NoShrinkColumns->asSexp())[iti] = -1;
    }
    int OnColumns = 0;
    for (int iti = 0 ; iti < Rf_length(iIn); iti++) {
      if (Rf_isInteger(iIn)) {
        if (INTEGER(iIn)[iti] <= 0 || INTEGER(iIn)[iti] > p) {
        } else {
          INTEGER(NoShrinkColumns->asSexp())[OnColumns] = INTEGER(iIn)[iti]-1;
          OnColumns++;
        } 
      } else if (Rf_isReal(iIn)) {
        if (REAL(iIn)[iti] <= 0.0 || REAL(iIn)[iti] > p) {
        } else {
          INTEGER(NoShrinkColumns->asSexp())[OnColumns] = (int) REAL(iIn)[iti]-1;
          OnColumns++;
        } 
      }
    }
    if (OnColumns < Rf_length(NoShrinkColumns->asSexp())) {
      AObject *AlterShrink = new AObject(Rf_allocVector(INTSXP, OnColumns));
      for (int iti = 0; iti < OnColumns; iti++)  {
        INTEGER(AlterShrink->asSexp())[iti] =INTEGER(NoShrinkColumns->asSexp())[iti];
      }
      delete(NoShrinkColumns); NoShrinkColumns = AlterShrink;
    }
    return;
  }
  SEXP get_ConfidenceQuantiles();
  SEXP get_ConfidenceMatrix();
  SEXP get_HPDMatrix();  SEXP get_HPDQuantiles(); SEXP SetHPDQuantiles();
  SEXP get_UnshrunkConfidenceMatrix();
  SEXP get_XjXj_CI_Vector() {
    if (XjXj_CI_Vector == NULL) { return(R_NilValue); }
    return(XjXj_CI_Vector->asSexp());
  }
  SEXP get_C_CI_Vector() {
    if (C_CI_Vector == NULL) { return(R_NilValue); }
    return(C_CI_Vector->asSexp());
  }
  SEXP get_Sig_CI_Vector() {
    if (Sig_CI_Vector == NULL) { return(R_NilValue); }
    return(Sig_CI_Vector->asSexp());
  }    
  void set_ConfidenceQuantiles(SEXP iConfidenceQuantiles);
  int GenerateConfidenceMatrix();
  int UpdateTDFNoise() {
    if (CDO != NULL) {
       return(CDO->UpdateTDFNoise());
    }
    return(-1);
  } 
  
  void set_LambdaIndexForConfidenceIntervals(int iIn) {
    if (iIn < 0 || iIn >= Rf_length(SLambdaDK)) {
      Rprintf("Error: Cannot set Lambda Index for ConfidenceIntervals to %d with LambdaDK length %d!\n",
        iIn, Rf_length(SLambdaDK));
      Rf_error("setLambdaIndexForConfidenceIntervals Error!\n");
    }
    LambdaIndexForConfidenceIntervals = iIn;
  }
  int get_LambdaIndexForConfidenceIntervals() { 
    return(LambdaIndexForConfidenceIntervals);
  }
  double TestAllInt(double A, double B, double PiA, double LambdaA, double LambdaD);
  double TestSuperInt(double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons, double b);
  double TestSeekQuantile(double QuantileSeek, double StartAt, double Epsilon,
    double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons){
    return(SeekQuantile(QuantileSeek, StartAt, Epsilon, A, B, PiA, LambdaA, LambdaD, NegCons));
  }
  int SetupGroupTwoLasso(SEXP iEndGroupList) {
  if (p <= 0) {
    Rprintf("SetupGroupTwoLasso: p=%d, you got problems.\n",p); R_FlushConsole();
  }
   if (Verbose >= 1) {
     Rprintf("SetupGroupTwoLasso: Start. \n"); R_FlushConsole();
   }
   int DoSSigma = 1;
   if (RSSigma == NULL) {
     Rprintf("   SetupGroupTwoLasso:  RSSigma is not Setup!\n"); R_FlushConsole(); DoSSigma = 0;    
   } else if (RSSigma != NULL && !Rf_isNull(SSigma) && RSSigma->asSexp() != SSigma) {
     Rprintf("   SetupGroupTwoLasso:  Uh Oh, SSigm and RSSigma do not co-point!\n"); DoSSigma=0; R_FlushConsole();
     DDelete(RSSigma, "RSSigma");  RSSigma = NULL; SSigma = R_NilValue;
   } else if (Rf_isNull(RSSigma->asSexp())) {
     Rprintf("   SetupGroupTwoLasso:  Uh Oh, SSigma is NULL object\n"); R_FlushConsole(); DoSSigma = 0;
     DDelete(RSSigma, "RSSigma");  SSigma = R_NilValue;
   } else if (!Rf_isReal(SSigma)) {
     Rprintf("   SetupGroupTwoLasso:  SSigma is Not Real!\n"); R_FlushConsole(); DoSSigma = 0;
     DDelete(RSSigma, "RSSigma");  SSigma = R_NilValue;
   } else if (RSSigma != NULL && Rf_length(SSigma)  < Rf_length(SLambdaAK)+1) {
     Rprintf("   SetupGroupTwoLasso: SSigma is not a good length, %d but SLambdaAK=%d\n", Rf_length(SSigma), Rf_length(SLambdaAK));
     R_FlushConsole();
     DDelete(RSSigma, "RSSigma");
     DoSSigma = 0;
   }
   if (DoSSigma == 0) {
     Rprintf("    SetupGroupTwoLasso:  We are setting up SSigma. \n"); R_FlushConsole();
     RSSigma = new AObject(Rf_allocVector(REALSXP, Rf_length(SLambdaAK)+1));
     SSigma = RSSigma->asSexp();
     for (int ii = 0; ii < Rf_length(SSigma); ii++) {
       REAL(SSigma)[ii] = OnSigma;
     }
   }
   
   if (Rf_isNull(iEndGroupList)) {
     Rf_error("GroupTwoLassoCpp:: Hey, iEndGroupList must be supplied!\n");
   } else if (Rf_length(iEndGroupList) <= 1) {
     Rf_error("GroupTwoLassoCpp:: Hey, iEndGroupList must have length > 1 for groups!\n");
   } else if (Rf_isInteger(iEndGroupList)) {
     FirstRandomIndex = INTEGER(iEndGroupList)[0];
     if (FirstRandomIndex < 0 || FirstRandomIndex >= p) {
       Rf_error("GroupTwoLassoCpp:: Hey, FirstRandomIndex = %d, not good!\n",
         FirstRandomIndex);
     }
     NumGroups = Rf_length(iEndGroupList)-1;
     EndGroupLocations = (int *) Calloc(NumGroups+1, int);
     if (INTEGER(iEndGroupList)[NumGroups] == p) {
     for (int jtjt = 0; jtjt < NumGroups; jtjt++) {
       if (INTEGER(iEndGroupList)[jtjt+1] < FirstRandomIndex ||
         INTEGER(iEndGroupList)[jtjt+1] < 
         INTEGER(iEndGroupList)[jtjt]) {
         Rf_error("GroupTwoLassoCpp:: Hey, iEndGroupList is not monotonic!\n");
       }
       EndGroupLocations[jtjt] = INTEGER(iEndGroupList)[jtjt+1]-1;
     }
     } else  if (INTEGER(iEndGroupList)[NumGroups] == p-1) {
     for (int jtjt = 0; jtjt < NumGroups; jtjt++) {
       if (INTEGER(iEndGroupList)[jtjt+1] < FirstRandomIndex ||
         INTEGER(iEndGroupList)[jtjt+1] < 
         INTEGER(iEndGroupList)[jtjt]) {
         Rf_error("GroupTwoLassoCpp:: Hey, iEndGroupList is not monotonic!\n");
       }
       EndGroupLocations[jtjt] = INTEGER(iEndGroupList)[jtjt+1];
     }     
     } else {
       Rprintf("GroupTwoLassoCpp: can't set EndGroupList because last iEndGroupList is %d\n",
         INTEGER(iEndGroupList)[NumGroups]);
       R_FlushConsole();
       Rf_error("Bad EndGroupList!\n");
     }
   } else if (Rf_isReal(iEndGroupList)) {
     FirstRandomIndex = (int) REAL(iEndGroupList)[0];
     if (FirstRandomIndex < 0  || FirstRandomIndex >= p) {
       Rf_error("GroupTwoLassoCpp:: Hey, FirstRandomIndex = %d, not good!\n",
         FirstRandomIndex);
     }
     NumGroups = Rf_length(iEndGroupList)-1;
     EndGroupLocations = (int *) Calloc(NumGroups+1, int);
     if (REAL(iEndGroupList)[NumGroups] == p) {
     for (int jtjt = 0; jtjt < NumGroups; jtjt++) {
       if (REAL(iEndGroupList)[jtjt+1] < FirstRandomIndex ||
         REAL(iEndGroupList)[jtjt+1] < 
         REAL(iEndGroupList)[jtjt]) {
         Rf_error("GroupTwoLassoCpp:: Hey, iEndGroupList is not monotonic!\n");
       }
       EndGroupLocations[jtjt] = (int) REAL(iEndGroupList)[jtjt+1]-1;
     }  
     } else if (REAL(iEndGroupList)[NumGroups] == p-1) {
     for (int jtjt = 0; jtjt < NumGroups; jtjt++) {
       if (REAL(iEndGroupList)[jtjt+1] < FirstRandomIndex ||
         REAL(iEndGroupList)[jtjt+1] < 
         REAL(iEndGroupList)[jtjt]) {
         Rf_error("GroupTwoLassoCpp:: Hey, iEndGroupList is not monotonic!\n");
       }
       EndGroupLocations[jtjt] = (int) REAL(iEndGroupList)[jtjt+1];
     }
     } else {
       Rf_error("GroupTwoLassoCpp: p = %d, but iEndGroupList[NumGroups=%d] = %d!\n",
         p, NumGroups, (int) REAL(iEndGroupList)[NumGroups]);
     }  
   }
   GroupLambdaEstimates = NULL;
   GroupLambdaEstimates = (double *) Calloc(NumGroups+1, double);
   GroupsBBOn1 = new AObject( 
     Rf_allocVector(REALSXP, NumGroups));
   BackGroupsBBOn1 = new AObject( 
     Rf_allocVector(REALSXP, NumGroups));
   if (Rf_isNull(SLambdaDK) || Rf_length(SLambdaDK) <= 0) {
     Rf_error("GroupTwoLassoCpp:: On Setup, SLambdaDK must already have length!\n");
   }
   RsGroupLambdaRecord = new AObject(
     Rf_allocMatrix(REALSXP, NumGroups, Rf_length(SLambdaDK))
   );
   return(1);
 }
 SEXP get_GroupsBBOn1() {
   if (GroupsBBOn1 == NULL) { return(R_NilValue); }
   return(GroupsBBOn1->asSexp());
 }
 SEXP get_RsGroupLambdaRecord() {
   if (RsGroupLambdaRecord == NULL) { return(R_NilValue); }
   return(RsGroupLambdaRecord->asSexp());
 }
SEXP get_BackGroupsBBOn1() {
   if (BackGroupsBBOn1 == NULL) { return(R_NilValue); }
   return(BackGroupsBBOn1->asSexp());
 }
 int get_NumGroups() { return(NumGroups); }
 SEXP get_EndGroupList() {
   if (EndGroupLocations == NULL || NumGroups <= 0) { return(R_NilValue); }
   SEXP sOut = R_NilValue;  Rf_protect(sOut = Rf_allocVector(INTSXP, NumGroups));
   for (int ii = 0; ii < NumGroups; ii++) {  INTEGER(sOut)[ii] = EndGroupLocations[ii]; }
   Rf_unprotect(1); return(sOut);
 }
 int get_FirstRandomIndex() { return(FirstRandomIndex); }
 SEXP get_GroupLambdaEstimates() {
   if (EndGroupLocations == NULL || NumGroups <= 0 || GroupLambdaEstimates  == NULL ||
     NumGroups <= 0) { Rprintf("Get_GroupLambdaEstimates: NULLS!\n"); R_FlushConsole();
     return(R_NilValue); }
   int One = 1;
   SEXP sOut = R_NilValue;
   Rf_protect(sOut = Rf_allocVector(REALSXP, NumGroups));
   F77_CALL(dcopy)(&NumGroups, GroupLambdaEstimates, &One, REAL(sOut), &One);
   Rf_unprotect(1); 
   return(sOut);
 }
 double get_SigmaBar() { return(SigmaBar); }
 double get_SigmaDf() { return(SigmaDf); }
 SEXP get_Weights() {
   if (GLMCDO != NULL) {
     return(GiveADoubleVectorOut(GLMCDO->nProbWeights, n, Verbose, "ProbWeights"));
   }
   if (RSiiWeights != NULL) {
     return(RSiiWeights->asSexp());
   }
   if (iiWeights == NULL) { return(R_NilValue); }
   SEXP sOut = R_NilValue;
   if (Verbose >= 2) {
     Rprintf("get_Weights: Have to pull from ii Weights. \n"); R_FlushConsole();
   }
   Rf_protect(sOut = Rf_allocVector(REALSXP, n));
   int One=1;
   F77_CALL(dcopy)(&n, iiWeights, &One, REAL(sOut), &One);
   Rf_unprotect(1); return(sOut);
 }
 SEXP get_RandBBOn1() {
   return(get_GroupsBBOn1());
   if (Rf_isNull(SBBOn1) || !Rf_isReal(SBBOn1)) {
     Rprintf("Can't get RandBBOn1 because SBBOn1 is bad.\n");
     R_FlushConsole();
     return(R_NilValue);
   }
   if (EndGroupLocations == NULL || NumGroups <= 0 || FirstRandomIndex < 0) { return(R_NilValue); }
   SEXP sOut = R_NilValue;  Rf_protect(sOut = Rf_allocVector(REALSXP, NumGroups));
   for (int ii = 0; ii < NumGroups; ii++) {  REAL(sOut)[ii] = REAL(SBBOn1)[FirstRandomIndex+ii];}
   Rf_unprotect(1); return(sOut);
 }

 SEXP get_iWeightedXtX() {
   if (CDO == NULL) { Rf_error("get_iWeightedXtX: CDO is NULL!\n"); }
   if (CDO->iWeightedXtX == NULL) { Rf_error("get_iWeightedXtX: iWeightedXtX is NULL!\n"); }
   if (CDO->iiWeights == NULL) { Rf_error("get)_iWeightedXtX: shouldn't be non null if iiWeights is NULL!\n"); }
   if (CDO->OnKappaS <= 0) { return(R_NilValue); }
   SEXP sOut = R_NilValue;
   Rf_protect(sOut = Rf_allocVector(INTSXP, CDO->OnKappaS));
   for (int ii = 0; ii < CDO->OnKappaS; ii++) { INTEGER(sOut)[ii] = CDO->iWeightedXtX[ii]; }
   Rf_unprotect(1); return(sOut);
 }
 int GroupTwoSetupCDO(int SetBeta);
 int GroupTwoUpdateBBOn1();
 int GroupTwoFixKCalculateBBOn();
 int get_MaxCauchy() { return(MaxCauchy); }
 void set_MaxCauchy(int iMaxCauchy) {
   if (iMaxCauchy <= 0) {
     Rf_error("set_MaxCauchy: Sorry, iMaxCauchy = %d, but must set to positive number \n",
      iMaxCauchy);
   }
   MaxCauchy = iMaxCauchy;
 }
 double get_CauchyEpsilon() { return(CauchyEpsilon); }
 void set_CauchyEpsilon(double iCauchyEpsilon) {
   if (iCauchyEpsilon <= 0) {
     Rf_error("set_CauchyEpsilon: Sorry, oCauchyEpsilon = %f, but must set to positive number \n",
      iCauchyEpsilon);
   }
   CauchyEpsilon = iCauchyEpsilon;
 }
 
 double get_Beta0() { return(Beta0); }
 double get_Beta0Prev() { return(Beta0Prev); }
 SEXP get_RBetasPrev() {
   if (RBetasPrev == NULL) { return(R_NilValue); }
   return(GiveADoubleVectorOut( (REAL(RBetasPrev->asSexp())), p, Verbose,"BetasPrev"));
 } 
 double get_CalculateBeta0() {
   if (GLMCDO == NULL) {
     Rf_error("get_CalculateBeta0: Won't work GLMCDO is null!\n"); R_FlushConsole();
   }
   return(GLMCDO->CalculateBeta0());
 }
 SEXP get_RBeta0Records() {
   if (RBeta0Records == NULL) { return(R_NilValue); }
   return(GiveADoubleVectorOut( (REAL(RBeta0Records->asSexp())), Rf_length(RBeta0Records->asSexp()), Verbose,"RBeta0Records"));
 } 
 SEXP get_RProbWeights() {
  if (RProbWeights == NULL) { return(R_NilValue); }
   return(GiveADoubleVectorOut( (REAL(RProbWeights->asSexp())), n, Verbose,"RProbWeights"));
 } 
 SEXP get_RLogitCurrentProb() {
   if (RLogitCurrentProb == NULL) { return(R_NilValue); }
   return(GiveADoubleVectorOut( (REAL(RLogitCurrentProb->asSexp())), n, Verbose,"RLogitCurrentProb"));
 } 
 SEXP get_RLogitPrevProb() {
   if (RLogitPrevProb == NULL) { return(R_NilValue); }
   return(GiveADoubleVectorOut( (REAL(RLogitPrevProb->asSexp())), n, Verbose, "RLogitPrevProb"));
 }
 int get_GiveCDOIter() {
   if (GiveCDOIter == NULL) { return(-1); }
   return(GiveCDOIter[0]);
 } 
  
 void SetBeta0(double InBeta0) {
   if (GLMCDO == NULL) {
     Rf_error("Error, can't setup Beta0 if GLMCDO is NULL!\n");
   }
   GLMCDO->SetBeta0(InBeta0);
   return;
 }
 int GLMRecordsNow(int tt1);
 int get_iiWeightsHere() {
   if (CDO == NULL) { Rf_error("iiWeights Here, no CDO\n"); }
   if (CDO->iiWeights == NULL) {
     Rprintf(" -- get_iiWeightsHere: There are no weights whatsoever. \n"); R_FlushConsole();
   }
   return(CDO->iiWeightsHere);
 }
 
 SEXP get_GLMBeta() {
   if (GLMCDO == NULL) {
     Rf_error("Error: there is no GLMBeta!");
   }
   SEXP sOn = R_NilValue;
   Rf_protect(sOn = Rf_allocVector(REALSXP, p));
   int One=1;
   F77_CALL(dcopy)(&p, GLMCDO->OnBeta, &One, REAL(sOn), &One);
   Rf_unprotect(1);
   return(sOn);
 }    
 double get_GLMTotMove() {
   if (GLMCDO == NULL) {
     Rf_error("GLMTotMove, no GLMB\n");
   }
   return(GLMCDO->GLMTotMove);
 }        
 int get_PrintFlagGLMB() {
   if (GLMCDO == NULL) {
     Rprintf("No GLMCDO! We are not a Logit structure.\n"); return(-1);
   }
   return(GLMCDO->PrintFlagGLMB);
 }
 int RefreshWeights() {
   if (CDO == NULL) {
     Rprintf("RefreshWeights: No, CDO is NULL\n"); return(-1);
   }
   if (CDO->iiWeights == NULL) {
     Rprintf("RefreshWeights: No, CDO->iiWeights is NULL\n"); return(-1);
   }
   return(CDO->RefreshWeights(CDO->iiWeights));
 }
 int SupplyNewWeights(SEXP sIn) {
   if (Rf_isNull(sIn) || Rf_length(sIn) != n || !Rf_isReal(sIn)) {
     Rprintf("SupplyNewWeights: No, there is no SIn\n");
   }
   if (RSiiWeights != NULL) {
     delete(RSiiWeights);  RSiiWeights = NULL;
     if (CDO != NULL && CDO->iiWeights != NULL) {
       if (CDO->iiWeightsHere == 1) {
         CDO->iiWeightsHere = 0; Free(CDO->iiWeights);  CDO->iiWeights = NULL;
       } else {
         CDO->iiWeightsHere = 0;
       }
     }
   }
   RSiiWeights = new AObject(Rf_allocVector(REALSXP, n));
   SiiWeights = RSiiWeights->asSexp();
   for (int ii = 0; ii < n; ii++) {
     REAL(RSiiWeights->asSexp())[ii] = REAL(sIn)[ii];
   }
   if (CDO->iiWeights == NULL) {
     CDO->SetupWeights(REAL(RSiiWeights->asSexp()));  CDO->iiWeightsHere = 0;
   }
   if (Verbose >= 1) {
     Rprintf("SupplyNewWeights: Supplied Setup, now Reweight\n");
     R_FlushConsole();
   }
   int AReturn = 0;
   AReturn = CDO->RefreshWeights(CDO->iiWeights);
   return(AReturn);
 }
 SEXP get_CDOWeights() {
   if (CDO == NULL) {
     Rf_error("get_CDOWeights: CDO is NULL \n");
   }
   if (CDO->iiWeights == NULL) {
     Rf_error("get_CDOWeights: CDO->iiWeights is NULL \n");
   }
   if (Verbose >= 1) {
     Rprintf("get_CDOWeights: well, iiWeights is going to be supplied. \n");
   }
   SEXP sOn = R_NilValue;
   Rf_protect(sOn = Rf_allocVector(REALSXP, n));
   for (int ii = 0; ii < n; ii++) {
     REAL(sOn)[ii] = CDO->iiWeights[ii];
   }
   Rf_unprotect(1);
   return(sOn);
 }
 double get_TDFSigma() {
   if (CDO == NULL) {
     Rf_error("get_CDOSigma: no, CDO is NULL!\n");
   }
   return(CDO->TDFSigma);
 }
  double get_CDOTDFSigma() {
   if (CDO == NULL) {
     Rf_error("get_CDOSigma: no, CDO is NULL!\n");
   }
   return(CDO->TDFSigma);
 }
  double get_CDOTDFNoise() {
   if (CDO == NULL) {
     Rf_error("get_CDONoise: no, CDO is NULL!\n");
   }
   return(CDO->TDFNoise);
 }
 int UpdateWeightLambdaALambdaDSeq() {
  if (Rf_length(SLambdaAK) >= 4) {
    return(3);
  }
  if (Verbose >= 1) {
    Rprintf("UpdateWeightLambdaALambdaDSeq running(tt1=%d, tt2=%d), p=%d, n=%d\n",
      tt1, tt2, p, n); R_FlushConsole();
  }
  if (CDO->iiWeights == NULL) {
    Rf_error("UpdateWeightLambdaALambdaDSeq: no weights. \n");
  } 
  if (GLMCDO == NULL) {
    Rf_error("No GLMCDO.  UpdateWeightLambdaALambdaDSeq. \n");
  }
  double MySD = 0.0;  double MinVar = -1.0;  double MyMean = 0.0;
  double SumWeights = 0.0;
  for (int ii = 0; ii < n; ii++) {
    SumWeights += CDO->iiWeights[ii];
  }
  if (SumWeights <= 0.0) {
    Rprintf("GLMCDO: UpdateWeightLambdaASeq, somehow SumWeights = %f\n", SumWeights);
    Rf_error("Error, no way we are doing this!\n");
  }
  if (Verbose >= 1) {
    Rprintf("GLMCDO, SumWeights = %f \n", SumWeights); R_FlushConsole();
  }
  int On = 0;  double Cons = 1.0 / SumWeights;
  for (int ii = 0; ii < p; ii++) {
    On = n * ii;
    MyMean = 0.0;  MySD = 0.0;
    for (int jj = 0; jj < n; jj++) {
      MyMean += CDO->iiWeights[jj] * REAL(Sxxs)[On+jj]*Cons;
    }
    for (int jj = 0; jj < n; jj++) {
      MySD += CDO->iiWeights[jj] * (REAL(Sxxs)[On+jj] - MyMean) * (REAL(Sxxs)[On+jj]-MyMean);
    }
    MySD = MySD;
    if (MinVar < 0) {
      MinVar = MySD;
    } else if (MySD < MinVar) {
      MinVar = MySD;
    }
  }
  int aMD = 0; if (SumWeights > aMD) { aMD =  SumWeights; }
  if (Rf_isNull(SLambdaDK)) {
    Rf_error("Error: SLambdaDK is not null. \n");
  }
  double SDD = sqrt(MinVar);
  

  
  if (Verbose >= 3 && OnSigma < 1.0) {
    Rprintf("TwoLassoObject2014.h::UpdateWeightLambdaALambdaDSeq: Note that OnSigma = %f, tt1=%d,tt2=%d\n", OnSigma, tt1, tt2);
    R_FlushConsole();
  }
  REAL(SLambdaDK)[0] = sqrt(3.7 / OnSigma) * SDD;
  REAL(SLambdaDK)[1] = aMD * REAL(SLambdaDK)[0];
  for (int iti = 2; iti < Rf_length(SLambdaDK); iti++) {
    if (Rf_length(SLambdaDK) > iti) {
      REAL(SLambdaDK)[iti] = REAL(SLambdaDK)[1];
    }
  }
  double SDA = SDD;
  if (SumWeights -1 >= 1) {
    SDA = SDD / sqrt(SumWeights-1.0);
  } else {
    SDA = SDD / sqrt(SumWeights);
  }
  double LambdaAS1 = exp( -TargetMinimumBeta * TargetMinimumBeta * SDA * SDA / (2.0*OnSigma)) *
    sqrt( 2* SDA*SDA / OnSigma);
  if (Verbose >= 1) {
    Rprintf("UpdateWeightLambdaALambdaDSeq Calculated SDD=%f, SumWeights=%f, SDA=%f, LambdaAS1 = %f\n", SDD, SumWeights, SDA, LambdaAS1); R_FlushConsole();
  }
  REAL(SLambdaAK)[0] = 1.1 * LambdaAS1;
  if (Rf_length(SLambdaAK) >= 2) {
  REAL(SLambdaAK)[1] = .05 * LambdaAS1;
  }
  if (Rf_length(SLambdaAK) > 2) {
    for (int iti = 2; iti < Rf_length(SLambdaAK); iti++) {
      REAL(SLambdaAK)[iti] = REAL(SLambdaAK)[1];
    }
  }
  if (REAL(SLambdaAK)[0] >= REAL(SLambdaDK)[0]) {
    double Swap = REAL(SLambdaAK)[0];
    REAL(SLambdaAK)[0] = REAL(SLambdaDK)[0];
    REAL(SLambdaDK)[0] = Swap;
  }
  if (tt1 >= 0 && tt1 < Rf_length(SLambdaAK)) {
    OnLambdaA = REAL(SLambdaAK)[tt1];
  }  
  if (tt1 >= 0 && tt1 < Rf_length(SLambdaDK)) {
    OnLambdaD = REAL(SLambdaDK)[tt1];
  } 
  if (Verbose >= 2) {
    Rprintf("Update LambdaAK(tt1=%d, tt2=%d), LambdaDK we have concluded. \n",tt1,tt2); R_FlushConsole();
    Rprintf("LambdaAK = "); PrintVector(REAL(SLambdaAK), Rf_length(SLambdaAK));
    Rprintf("\nLambdaDK = "); PrintVector(REAL(SLambdaDK), Rf_length(SLambdaDK));
    Rprintf("\n"); R_FlushConsole();
  }
  int UpdateLambdaAKNAN = 0;
  for (int ii = 0; ii < Rf_length(SLambdaAK); ii++) {
    if (R_isnancpp(REAL(SLambdaAK)[ii])) { UpdateLambdaAKNAN++; }
  }
  for (int ii = 0; ii < Rf_length(SLambdaDK); ii++) {
    if (R_isnancpp(REAL(SLambdaDK)[ii])) { UpdateLambdaAKNAN++; }
  }  
  if (R_isnancpp(OnLambdaA)  || R_isnancpp(OnLambdaD) || UpdateLambdaAKNAN > 0) {
    Rprintf("ERRORERRORERRORERRORERRORERRORERRORERRORERRORERRORERRRORERRORERRORERRORERRORERROR\n");
    Rprintf("ERROR: TwoLassoObject2014.h::UpdateWeightLambdaALambdaDSeq() --  \n");
    Rprintf("ERROR:   We Have Nans OA=%d, OD=%d, UpdateLambdaAKAN = %d \n",
      R_isnancpp(OnLambdaA) ? 1 : 0, R_isnancpp(OnLambdaD) ? 1 : 0,
      UpdateLambdaAKNAN); R_FlushConsole();
    Rprintf("ERROR: Note  SDA = %f, SDD=%f, SumWeights=%f, LambdaAS1 = %f, TargetMinimumBeta=%f, OnSigma=%f\n",
      SDA, SDD, SumWeights, LambdaAS1, TargetMinimumBeta, OnSigma); R_FlushConsole();
    Rprintf("ERROR: Check Out for this error. \n");
    Rprintf("ERROR: SLambdaAK: "); PrintVector(REAL(SLambdaAK), Rf_length(SLambdaAK));
    Rprintf("\nERROR: SLambdaDK: "); PrintVector(REAL(SLambdaDK), Rf_length(SLambdaDK));
    Rprintf("\nERROR: ");
    SuccessFlag = -1;
    Rf_error("ERROR: TwoLassoObject2014.h:: UpdateWeightLambdaALambdaDSeq Fail!\n");
  }
  
  if (FALSE && RSPiAVectorInputs != NULL) {
    Rprintf("TwoLassoObject2014.h::UpdateWeightLambdaALambdaDSeq: Checking CVDataIntegrity. \n"); R_FlushConsole();
    CVDataIntegrity(0);                                                                           R_FlushConsole();
  }
  return(1);
}
 void set_LogitNoise(double InNoise) {
   if (InNoise <= 0) {
     Rf_error("set_LogitNoise, no not setting InNoise = %f \n", InNoise);
   }
   if (GLMCDO == NULL) {
     Rf_error("set_LogitNoise, why, if this is not a Logit regression? \n");
   }
   OnSigma = InNoise;
   if (tt1 >= 0 && !Rf_isNull(SSigma) && Rf_length(SSigma) > tt1) {
     REAL(SSigma)[tt1] = OnSigma;
   }
   UpdateWeightLambdaALambdaDSeq();
 }
 double get_LogitNoise() {
   if (GLMCDO == NULL) {
     Rf_error("get_LogitNoise, why, if this is not a Logit regression? \n");
   }
   return(OnSigma);
 }
 //void set_PrintFlagGLMB(int iPrintFlagGLMB) {
 //  if (GLMCDO == NULL) {
 //    Rprintf("No GLMCDO! We are not a Logit structure.\n"); return;
 //  }
 //  GLMCDO->PrintFlagGLMB = iPrintFlagGLMB;
 //}
}; 

double AllInt(double A, double B, double PiA, double LambdaA, double LambdaD);
double SuperInt(double A, double B, double PiA, double LambdaA, double LambdaD, double NegCons, double b);


SEXP TwoLassoRegression(SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaA, SEXP rSLambdaD, SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
    SEXP rSVerbose, SEXP rGroupSexp);
SEXP TwoLassoReRegress(
  SEXP sTLS, SEXP sRunFlag, SEXP newMvec, SEXP newSigmaPrior,
  SEXP newVerbose);
SEXP TwoLassoCreateTwoLassoSexp(
  SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
  SEXP rSXtX, SEXP rStt1, SEXP rStt2,
  SEXP rSLambdaA, SEXP rSLambdaD, SEXP rSOrderSeq, 
  SEXP rSBBOn1, SEXP z, SEXP rRecBBOn1, SEXP rRecOnBeta,
  SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
  SEXP rSm, SEXP rSSigmaPrior,
  SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
  SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
  SEXP rSTDFNu, SEXP rSiiWeights, 
  SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
  SEXP rSVerbose,
  SEXP rGroupSexp, SEXP RunFlag);
void DeleteSTLS(SEXP sTLS);

SEXP TwoLassoCrossValidate(SEXP rSVerbose, SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaA, SEXP rSLambdaD, 
    SEXP rSPiAVectorInputs, SEXP rSSigmaVectorInputs,
    SEXP rSLambdaAInputs, SEXP rSLambdaDInputs, 
    SEXP rSBetaStartMatrix,
    SEXP rSBetaCVInputs,
    SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
    SEXP rGroupSexp);


SEXP sWhatIsSuccessFlag(SEXP sTLS);
SEXP sSetupConfidenceQuantiles(SEXP sTLS, SEXP iConfidenceQuantiles);
SEXP sgetConfidenceMatrix(SEXP sTLS); 
SEXP sgetHPDMatrix(SEXP sTLS);
SEXP sgetConfidenceQuantiles(SEXP sTLS); 
SEXP sgetUnshrunkConfidenceMatrix(SEXP sTLS);
SEXP sgetLambdaIndexForConfidenceIntervals(SEXP sTLS);
SEXP ssetLambdaIndexForConfidenceIntervals(SEXP sTLS, SEXP sIn);

SEXP ssetVerbose(SEXP sTLS, SEXP siVerbose);
SEXP sgetVerbose(SEXP sTLS);



