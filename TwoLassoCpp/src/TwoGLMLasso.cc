/* ========================================================================== */
/*                                                                            */
/*   TwoGLMLasso.cc                                                           */
/*   (c) 2010 Alan Lenarcic                                                   */
/*                                                                            */
/*   TwoLasso for Logistic Regression                                         */
/*       Work with Edoardo Airoldi Lab                                        */
/*                                                                            */
/*       Work accompanies effort by Lenarcic and Valdar on BayesSpike         */
/*       This is a companion algorithm using Coordinate Descent and EM        */
/*       to generate Model Inclusion Probability estimates and credibility    */
/*       without using Gibbs Sampling Integration.                            */
/*                                                                            */
/*                                                                            */
/*    This extends TwoLassoSEXP object to perform GLM-Binomial Coordinate     */
/*  Descent.  Performance of TwoLasso seems noticably less when dealing       */
/*  with binomial data, and there is significant information lost in 0/1      */
/*  versus continuous Ydata.                                                  */
/*  continuous Y data.                                                        */
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
/******************************************************************************/
#ifndef TWOGLMLASSODD
   #include "TwoGLMLasso.h"
   #define TWOGLMLASSODD 0
#endif


#ifndef  TwoLarsObject2009DD        
  #include "TwoLarsObject2009.h"
  #define TwoLarsObject2009DD 0
#endif

#ifndef CoordinateDescentLassoDD
  #include "CoordinateDescentLasso2014.h"
  #define CoordinateDescentLassoDD 0
#endif


#ifndef LARSOBJECT2009DD
  #include "LarsObject2009.h"
  #define LARSOBJECT2009DD 0
#endif 

#ifndef INCIDENTALOP2009DD
  #include "IncidentalOp2009.h"
  #define INCIDENTALOP2009DD 0 
#endif
#ifndef MyMatrixOp2009DD
  #include "MyMatrixOp2009.h"
  #define MyMatrixOp2009DD 0
#endif


#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <Rdefines.h>
  #include <Rinternals.h>
  #include <stdio.h>
  #define RMATH 0
#endif

//    TwoGLMLassoObject(int rn, int *ry01, double *rX, int rp, 
//              double *rBetas,  double *rpBeta0,  double *rOnGammas,
//              double *rBBHats, double *rpRecordBetas,
//              double *rpBeta0Records, double *rpRecordBBHats,
//              double rStlambdaD, double rStlambdaA,
//              double rsigmaNoiseSq, double rMultDConst, double rMultAConst,
//              int rStopAInt, double rppiuse, int rTotalRuns,
//              int rInitKKs, double InverseGammaConstant, int PrintFlag,
//              double *rLambdaDK, double *rLambdaAK, int *rOrderSeq,
//              double *rnProbWeights, double *rnlogitCurrentProb,
//              double *rz);

SEXP RunG2GLMLasso2(SEXP sXX, SEXP sy01, SEXP sBetas, SEXP sBeta0,
                   SEXP sOnGammas, SEXP sBBHats, SEXP sRecordBetas,
                   SEXP sBeta0prev, SEXP sBetasprev,
                   SEXP sBeta0Records, SEXP sRecordBBHats,
                   SEXP sLambdaDK, SEXP sLambdaAK,
                   SEXP sppiuse, SEXP sPiRecVec,
                   SEXP ssigmaNoiseSq, SEXP sSigmaRecVec,
                   SEXP sOrderSeq, SEXP sProbWeights, SEXP snLogitCurrentProb,
                   SEXP snlogitPrevProb,
                   SEXP sZ, SEXP sInverseGammaConstant, 
                   SEXP sInitKKs, SEXP sPrintFlag,
                   SEXP sMaxCDOEpsilon, SEXP sMaxCDORuns,
                   SEXP sIter, SEXP sFixKa,
                   SEXP sPriorS, SEXP sInMaximumAllocation) {
   int PrintFlag2 = 0;
   if (Rf_isInteger(sPrintFlag)) {PrintFlag2 = (int) INTEGER(sPrintFlag)[0];} else {
      PrintFlag2 = (int) REAL(sPrintFlag)[0];}
  if (PrintFlag2 > 0) {
     Rprintf("TwoLarsOperator: Starting the GLM2Lasso Loads \n\n"); R_FlushConsole();
  }                   
  TwoGLMLassoObject *TG2O = NULL;
  int n = INTEGER(RGET_DIM(sXX))[0]; int p =INTEGER(RGET_DIM(sXX))[1];
  int yLen = Rf_length(sy01);  
  if (PrintFlag2 > 1) {
     Rprintf("We get n = %d, p =%d, yLen = %d\n", n, p, yLen); R_FlushConsole();
  }
  
  if (yLen != n) { 
    Rprintf("yLen = %d, n = %d, not same length, reject!\n", yLen, n); 
    R_FlushConsole();
    return(sy01);
  }
  int BetasLen = Rf_length(sBetas);
  if (BetasLen < p) {
    Rprintf("BetasLen = %d, p = %d, not same length, reject!\n", BetasLen, p);
    R_FlushConsole(); return(sBetas);
  }
  if (PrintFlag2 > 2) {
    Rprintf("We've got Betas, here's sBetas\n"); R_FlushConsole();
    PrintVector(REAL(sBetas), p); R_FlushConsole();
  }
  int BBHatsLen = Rf_length(sBBHats);
  if (BBHatsLen < p) {
    Rprintf("BBHatsLen = %d, p = %d, not same length, reject!\n", BBHatsLen, p);
    R_FlushConsole(); return(sBetas);
  }
  if (Rf_length(sBetasprev) < p) {
   Rprintf("sBetasprev not long enough=%d,not %d\n", 
   Rf_length(sBetasprev), p);  R_FlushConsole();return(sBetasprev);
  }  
  double *pBetasPrev = REAL(sBetasprev);  
  double *pBeta0prev =REAL(sBeta0prev);
  if (Rf_length(sOnGammas) < p) {
    Rprintf(" sOnGammaslength = %d, too short to p = %d\n", Rf_length(sOnGammas),p);
    R_FlushConsole(); return(sOnGammas);
  }
  int OuterNumLoops = (int) Rf_length(sLambdaDK);
  if (Rf_length(sOrderSeq) < OuterNumLoops) {
    Rprintf(" OrderSeq[%d] is Too Small to %d", Rf_length(sOrderSeq), OuterNumLoops);
    R_FlushConsole(); return(sOrderSeq);
  }
  if (Rf_length(sLambdaAK) < OuterNumLoops) {
    Rprintf(" LambdaA[%d] is Too Small to %d", Rf_length(sLambdaAK), OuterNumLoops);
    R_FlushConsole(); return(sLambdaAK);
  }  
  double *pRecordBBHats = NULL;  double *pRecordBetas =NULL; 
  double *pBeta0Records=NULL;
  if (PrintFlag2 > 2) {
    Rprintf("Judging sRecordBBHats has length %d\n", Rf_length(sRecordBBHats));
    R_FlushConsole();
  }
  if (Rf_length(sSigmaRecVec) < OuterNumLoops) {
    Rprintf(" sSigmaRecVec is not Long Enough\n"); R_FlushConsole();
    R_FlushConsole();
    return(sSigmaRecVec);
  }
  if (Rf_length(sRecordBBHats) < p * OuterNumLoops) {
    if (Rf_length(sRecordBBHats) <= 1) {
    } else {
       Rprintf("sRecordBBHats not long enough, turning to NULL \n");
       R_FlushConsole();
    }
    pRecordBBHats = NULL;
  } else {pRecordBBHats=REAL(sRecordBBHats);}
  if (PrintFlag2 > 2) {
    Rprintf("And dim sRecordBetas =c(%d,%d)\n",
       INTEGER(RGET_DIM(sRecordBetas))[0],  INTEGER(RGET_DIM(sRecordBetas))[1]);
    R_FlushConsole();
  }
  if (Rf_length(sRecordBetas) < p * OuterNumLoops) {
    if (Rf_length(sRecordBetas) <= 1) {
    } else {
       Rprintf("sRecordBetas not long enough, turning to NULL \n");
       R_FlushConsole();
    }
    pRecordBetas = NULL;
  } else {pRecordBetas=REAL(sRecordBetas);}  
  if (PrintFlag2 > 2) {
    Rprintf("And length sBeta0Records =%d\n",
       Rf_length(sBeta0Records));
    R_FlushConsole();
  }
  if (Rf_length(sBeta0Records) < OuterNumLoops) {
    if (Rf_length(sBeta0Records) <= 1) {
    } else {
       Rprintf("sBeta0Records not long enough, turning to NULL \n");
       R_FlushConsole();
    }
    pBeta0Records = NULL;
  } else {
    if (sBeta0Records != NULL) {
       pBeta0Records=REAL(sBeta0Records);
    }
  }
  if (Rf_length(sPiRecVec) < OuterNumLoops) {
    Rprintf("Error, sPiRecVec is not Long enough\n"); R_FlushConsole();
    return(sPiRecVec);
  }
  if (PrintFlag2 > 2) {
    Rprintf("And length sBeta0Records =%d\n",
       Rf_length(sBeta0Records));
    R_FlushConsole();
  }
  if (Rf_length(snlogitPrevProb)< n){
     Rprintf("Error, Rf_length(snlogitPrevProb)=%d, smaller than n=%d\n",
      Rf_length(snlogitPrevProb), n); R_FlushConsole();  return(snlogitPrevProb);
  } 
  int *py01;  int Madepy01 = 0;
  if (PrintFlag2 > 2) {
    Rprintf("Check up on sy01 length =%d\n",
       Rf_length(sy01));
    R_FlushConsole();
  }
  if (!Rf_isInteger(sy01)) {
    if (PrintFlag2 > 2) {
      Rprintf("Uh Oh, sy01, not integer, got to make\n",
       Rf_length(sy01));
      R_FlushConsole();
    }
       Madepy01 = 1;
       py01 = (int*) Calloc(Rf_length(sy01)+1, int);
       for (int iti = 0; iti < Rf_length(sy01); iti++) {
         py01[iti] = (int) REAL(sy01)[iti];
       }
  } else {
      py01 = INTEGER(sy01);
  }   
  int PrintFlag;
  if (Rf_isInteger(sPrintFlag)) {
    PrintFlag = (int) INTEGER(sPrintFlag)[0];
  } else {
    PrintFlag = (int) REAL(sPrintFlag)[0];
  }              
  int InitKKs;
  if (PrintFlag2 > 2) {
    Rprintf("What is sInitKKs I wonder?\n",
       Rf_length(sy01));
    R_FlushConsole();
  }
  if (Rf_isInteger(sInitKKs)) {
    InitKKs = (int) INTEGER(sInitKKs)[0];
  } else {
    InitKKs = (int) REAL(sInitKKs)[0];
  }
  if (!Rf_isNumeric(sMaxCDOEpsilon)) {
   Rprintf("Error sMaxCDOEpsilon is not a Real\n"); 
   R_FlushConsole(); return(sMaxCDOEpsilon);
  }
  if (PrintFlag2 > 2) {
    Rprintf("sMaxCDOEpsilon is %.4f\n", REAL(sMaxCDOEpsilon)[0]);
    R_FlushConsole();
  }  
  int MaxCDORuns;
  if (PrintFlag2 > 2) {
    Rprintf("What is sMaxCDORuns I wonder?\n");
    R_FlushConsole();
  }

  if (sMaxCDORuns == NULL) {
   Rprintf("Error sMaxCDORuns is NULL Error \n"); R_FlushConsole();
   return(sMaxCDORuns);
  } else if (Rf_isReal(sMaxCDORuns)) {
   MaxCDORuns = (int) REAL(sMaxCDORuns)[0];
  }  else if (Rf_isInteger(sMaxCDORuns))  {
      MaxCDORuns = INTEGER(sMaxCDORuns)[0];
  }  else {
      Rprintf("MaxCDORuns is defective List, oopso \n");
      R_FlushConsole(); return(sMaxCDORuns);
  }
  int *pIter;  int MakepIter;
  if (PrintFlag2 > 2) {
    Rprintf("What is sIter,do I make?\n");
    R_FlushConsole();
  }
  if (sIter != NULL && Rf_isInteger(sIter)) {
     pIter = INTEGER(sIter);  pIter[0] = 0; MakepIter = 0;
  }  else {
    pIter = (int*) Calloc(2, int);  pIter[0] = 0; MakepIter = 1;
  }
  int InMaximumAllocation = -1;
  if (InMaximumAllocation != -1) {
    Rprintf("--- TwoLassoCpp:src: Hey Why is InMaximumAllocation not equal to -1?\n"); R_FlushConsole();
  }
  if (!Rf_isNull(sInMaximumAllocation)) {
    if (Rf_isReal(sInMaximumAllocation)) {
      InMaximumAllocation = REAL(sInMaximumAllocation)[0];
    } else if (Rf_isInteger(sInMaximumAllocation)) {
      InMaximumAllocation = INTEGER(sInMaximumAllocation)[0];
    }
  }
  if (sppiuse== NULL || !Rf_isReal(sppiuse)) {
      Rprintf("What is wrong with sppiuse? It's NULL?\n"); R_FlushConsole();
      return(sppiuse);
  }
  if (sOnGammas== NULL || !Rf_isReal(sOnGammas)) {
      Rprintf("What is wrong with sppiuse? It's NULL?\n"); R_FlushConsole();
      return(sOnGammas);
  } 
  if (sInverseGammaConstant== NULL || !Rf_isReal(sInverseGammaConstant)) {
      Rprintf("What is wrong with sInverseGammaConstant? It's NULL?\n"); R_FlushConsole();
      return(sInverseGammaConstant);
  }  
  if (PrintFlag2 > 2) {
    Rprintf("Well, we're going into make TG2O, wish us luck.\n");
    R_FlushConsole();
  }   
  int *pOrderSeq = NULL; int MakeOrderSeq = 0;
  if (sOrderSeq == NULL) {
     Rprintf("sOrderSeq Sucks, its null \n"); R_FlushConsole(); 
     return(sOrderSeq);
  } else if (Rf_isInteger(sOrderSeq)) {
    pOrderSeq =INTEGER(sOrderSeq); MakeOrderSeq = 0;
  } else if (Rf_isReal(sOrderSeq)){
    pOrderSeq = (int*)Calloc(OuterNumLoops, int);
    for (int itii = 0; itii < OuterNumLoops; itii++){
      pOrderSeq[itii] = (int) REAL(sOrderSeq)[itii];
    }
    MakeOrderSeq = 1;
  } else {
    Rprintf("sOrderSeq aint good, fail \n"); R_FlushConsole();
    return(sOrderSeq);
  }
  if (Rf_length(sBBHats) < p) {
    Rprintf("sBBHats is not long enough, quit\n"); R_FlushConsole();
    return(sBBHats);
  }
  int FixKa = -1;
  if (Rf_isInteger(sFixKa)) {
     FixKa = (int) INTEGER(sFixKa)[0];
  } else if (Rf_isNumeric(sFixKa)) {
     FixKa = (int) REAL(sFixKa)[0];
  }
  TG2O =new TwoGLMLassoObject((int)n, (int*) py01,
   (double*) REAL(sXX), (int) p, 
   (double*)REAL(sBetas),  REAL(sBeta0),  REAL(sOnGammas),
   REAL(sBBHats), pRecordBetas, pBeta0Records, pRecordBBHats,
   pBeta0prev, pBetasPrev,   
   REAL(sppiuse)[0], REAL(sPiRecVec), (int)OuterNumLoops,  REAL(ssigmaNoiseSq)[0],
   REAL(sSigmaRecVec),
    (int)InitKKs, (double)REAL(sInverseGammaConstant)[0],
    (int)PrintFlag,
// REAL(sLambdaDK)[0], REAL(sLambdaAK)[0],
//              double rsigmaNoiseSq, double rMultDConst, double rMultAConst,
  (double*)REAL(sLambdaDK), (double*)REAL(sLambdaAK), (int*)pOrderSeq,
  REAL(sProbWeights), REAL(snLogitCurrentProb),
  (double*)REAL(snlogitPrevProb),
  (double*)REAL(sZ),
  (double) REAL(sMaxCDOEpsilon)[0], (int) MaxCDORuns,
  (int*) pIter, -10); 
  if (Rf_length(sPriorS) ==2 && Rf_isReal(sPriorS) && 
      REAL(sPriorS)[0] > 0) {
      TG2O->SetupPiEst( REAL(sPriorS)[0], REAL(sPriorS)[1]);    
  }
  if (TG2O->IPLLF > 0 ) {
     Rprintf("TG2O created: now to run algorithm\n");
     R_FlushConsole();
  }
  if (FixKa > 0 && FixKa < n) {
    TG2O->SetFixKa(FixKa);
  }  
  TG2O->RunG2GLMLasso();
  if (Madepy01 == 1) {
     Free(py01); py01 = NULL;
  }
  if (MakepIter == 1) {
     REAL(sIter)[0] = pIter[0];
     TG2O->deletepIter();
     pIter = NULL;
  }
  if (MakeOrderSeq == 1){
     Free(pOrderSeq); pOrderSeq = NULL;  MakeOrderSeq = 0;
  }
  if (PrintFlag2 > 2) {
    Rprintf("About to delete TG2O\n"); R_FlushConsole();
  }
  delete(TG2O);
  if (PrintFlag2 > 2) {
    Rprintf("Done Deleting TG2O \n"); R_FlushConsole();
  }
  return(sBetas);              
} 
int  TwoGLMLassoObject::RunG2GLMLasso()
{
  double CountMoves =0.0; double CurrentMoves = 0.0;
  pIter[0] = 0;
  if (IPLLF > 0) {
     Rprintf("Starting RunG2GLMLasso\n");
     R_FlushConsole();
  }
  if (IPLLF > 0){
    Rprintf("CDOIter = %d, GiveCDOIter = %d, just a warning\n", 
       GiveCDOIter[0], GLMO->pIter[0]); R_FlushConsole();
  }
  for (tt1 = 0; tt1 < TotalRuns; tt1++) {
   if (IPLLF > 2) {
     Rprintf("Setting pIter\n"); R_FlushConsole();
   }
   pIter[0] = tt1 + 1;
   if (IPLLF > 2) {
     Rprintf("Finished Setting pIter, now taking OnLambdaA\n"); R_FlushConsole();
   }
   if (LambdaAK == NULL) {
     Rprintf("Error, you set LambdaAK to NULL\n"); R_FlushConsole();
     SuccessFlag = -1; return(-1);
   }
   if (LambdaDK == NULL) {
     Rprintf("Error, you set LambdaDK to NULL\n"); R_FlushConsole();
     SuccessFlag = -1; return(-1);
   } 
   if (OrderSeq == NULL) {
     Rprintf("Error, you set OrderSeq to NULL\n"); R_FlushConsole();
     SuccessFlag = -1; return(-1);
   }    
   OnlambdaA = LambdaAK[tt1];  OnlambdaD = LambdaDK[tt1];
   CountMoves = 0.0;
   //Rprintf(" On Iter %d, OnlambdaA = %f, OnlambdaD = %f\n",
   //  tt1, OnlambdaA, OnlambdaD); R_FlushConsole();
   if (IPLLF> 1) {
     Rprintf("G2GLMLasso iter %d, OnlambdaA = %.4f, OnlambdaD =%.4f\n",
          tt1, OnlambdaA, OnlambdaD);
     R_FlushConsole();
   }
   for (tt2=0; tt2 < abs(OrderSeq[tt1]); tt2++) {
     if (IPLLF>4) {
        Rprintf("G2GLMLasso iter %d, On tt2 = %d / %d\n",
          tt1, tt2, abs(OrderSeq[tt1]));
        R_FlushConsole();
     }
     if (OrderSeq[tt1] < 0) {
       if (IPLLF > 4) {
         Rprintf("     G2GLMLassso, tt1=%d, tt2=%d, Updating BBHats\n",tt1,tt2);
         R_FlushConsole();
       }
       UpdateBBHats();
       if (IPLLF > 4) {
          Rprintf("      G2GLMLassso, tt1=%d,tt2=%d, Updating Gammas\n", tt1, tt2);
          R_FlushConsole();
        }
        UpdateGammas();
        if (IPLLF > 4) {
          Rprintf("      G2GLMLassso, tt1=%d, tt2=%d, Running the GLMLasso\n", tt1,tt2);
          R_FlushConsole();
        }        
        CurrentMoves = GLMO->RunGLMLasso();
      } else {
        if (IPLLF > 4) {
          Rprintf("      G2GLMLassso, tt1=%d, tt2=%d, Running the GLMLasso\n", tt1,tt2);
          R_FlushConsole();
        }       
        CurrentMoves = GLMO->RunGLMLasso();
        if (IPLLF > 4) {
          Rprintf("      G2GLMLassso, tt1=%d, tt2=%d, Running the GLMLasso\n", tt1,tt2);
          R_FlushConsole();
        }    
        UpdateBBHats();
        if (IPLLF > 4) {
          Rprintf("      G2GLMLassso, tt1=%d, tt2=%d, UpdateGamma\n", tt1,tt2);
          R_FlushConsole();
        }    
        UpdateGammas();
      }
      
      CountMoves=CountMoves + CurrentMoves;
      if (IPLLF >= 2) {
        Rprintf("RunGLM2Lasso: CurrentMoves = %.8f, but MaxEpsilon = %.8f\n",
           CurrentMoves,GLMO->GMaxEpsilon());
        
      }
      if (CurrentMoves < GLMO->GMaxEpsilon()) {
         break;
      }
   }
   RecordsNow(tt1);
   if (CountMoves < GLMO->GMaxEpsilon()) {
       for (;tt1< TotalRuns; tt1++) {
           RecordsNow(tt1);
       }
       return(1);
   }
  }
  return(1);
}
int TwoGLMLassoObject::RecordsNow(int tt1) {
    if (IPLLF > 2) {
      Rprintf("Implementing the Records Now, tt1 = %d \n", tt1);
      R_FlushConsole();
    }
    int jj=0;  int jjtt1 = tt1 *p ;
    if (pRecordBetas != NULL) {
       for (jj = 0; jj < p; jj++) {
         pRecordBetas[jjtt1] = OnBetas[jj];
         jjtt1++;
       } 
    }
    if (pBeta0Records != NULL) {
      pBeta0Records[tt1] = GLMO->pBeta0[0];
    }
    if (pRecordBBHats != NULL) {
      jjtt1 = tt1 * p;
      for (jj=0; jj < p; jj++) {
        pRecordBBHats[jjtt1] = BBOn1[jj];
        jjtt1++; 
      }
    }
    if (IPLLF > 2) {
      Rprintf("Done Recording \n"); R_FlushConsole();
    }
    return(1);
}

///////////////////////////////////////////////////////////////////////////////
// TwoGLMLasso:  Updating BBHats based upon on Odds for both variables
int TwoGLMLassoObject::UpdateBBHats() {
  int jj;
  //double OnlambdaDiff = (OnlambdaD - OnlambdaA) / 2.0;
  int SuccFlag = 0;
	if (FixKa > 0) {
		SuccFlag = FixKCalculateBBOn();
	  if (SuccFlag < 0) {
			Rprintf(" CalculateBBOn1, we got a return FixKCalculateBBOn is angry \n");
			R_FlushConsole(); return(-1);
		}
		return(1);
  }
	double ZeroDud = OnlambdaA  * ppiuse/ (OnlambdaA *ppiuse + (1.0-ppiuse) * OnlambdaD);
  double OnOdds = 0.0;
  double lnFakerOdds = log ( OnlambdaA * ppiuse / ( (1.0-ppiuse) *OnlambdaD) );
  for (jj=0;jj < p; jj++) {
    if (OnBetas[jj] ==0) {
      BBOn1[jj] = ZeroDud;
    } else {
      OnOdds = lnFakerOdds  + (OnlambdaD -OnlambdaA) *
        fabs(OnBetas[jj]);
      if (OnOdds > 40) {
        BBOn1[jj] = 1.0;
      } else if (OnOdds < -40) {
        BBOn1[jj] = 0.0;
      } else {
        BBOn1[jj] = exp(OnOdds) / ( exp(OnOdds) + 1.0);
      }
    }    
  }
  return(1);
}

////////////////////////////////////////////////////////////////////////
//  In GLMO the OnGammas are the penalty weights for each coordinate j
//
int TwoGLMLassoObject::UpdateGammas() {
 int jj;
 for(jj=0; jj < p; jj++) {
    GLMO->OnGammas[jj] =  sigmaNoiseSq * ( BBOn1[jj] *OnlambdaA + 
        ( 1.0 - BBOn1[jj]) *OnlambdaD );   
 }
 return(1);
}

int TwoGLMLassoObject::SetupPiEst(double rPriorM1, double rPriorM2) {
   PriorM1 = rPriorM1; PriorM2 = rPriorM2;
   return(1);
}
int TwoGLMLassoObject::UpdatePiEst() {
  int StC = 0;  double SumBBOn1 = 0; int TotFree = p;
  if (X[0] == 1 && X[1] == 1) {
     StC = 1;  TotFree = p-1;
  }
  for(; StC < p; StC++) {
     SumBBOn1 += BBOn1[StC];
  }
  ppiuse = (PriorM1 + SumBBOn1 - 1) /
    (PriorM1 + PriorM2 + TotFree - 2);

  if (PiRecVec != NULL) {
     PiRecVec[tt1] = ppiuse;
  }
  return(1);
}
