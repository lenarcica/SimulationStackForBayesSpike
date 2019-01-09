/* ========================================================================== */
/*                                                                            */
/*   TwoGLMLasso.h                                                            */
/*   (c) 2009 2019 Alan Lenarcic                                              */
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
#ifndef RCPPH
  #include <Rcpp.h>
  #include <cmath>
  #define RCPPH 0
#endif
#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <Rdefines.h>
  #include <Rinternals.h>
  #include <stdio.h>
  #define RMATH 0
#endif

#ifndef RGET_DIM
  #define RGET_DIM(X) Rf_getAttrib(X, R_DimSymbol)
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

#ifndef SMALLESTInitKKs
   #define SMALLESTInitKKs 20
#endif
#ifndef SMALLESTkLen
 #define SMALLESTkLen 20
#endif
#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif

SEXP RunG2GLMLasso2(SEXP sXX, SEXP sy01, SEXP sBetas, SEXP sBeta0,
                   SEXP sOnGammas, SEXP sBBHats, SEXP sRecordBetas,
                   SEXP sBeta0prev, SEXP sBetasprev,
                   SEXP sBeta0Records, SEXP sRecordBBHats,
                   SEXP sLambdaDK, SEXP sLambdaAK,
                   SEXP sppiuse, SEXP sPiRecVec,
                   SEXP ssigmaNoiseSq,
                   SEXP sSigmaRecVec,
                   SEXP sOrderSeq, SEXP sProbWeights, SEXP snLogitCurrentProb,
                   SEXP snlogitPrevProb,
                   SEXP sZ, SEXP sInverseGammaConstant, 
                   SEXP sInitKKs, SEXP sPrintFlag,
                   SEXP sMaxCDOEpsilon, SEXP sMaxCDORuns,
                   SEXP sIter, SEXP s, SEXP sPriorS,
                   SEXP sInMaximumAllocation);
                   
class TwoGLMLassoObject:public  ConvergentLARS  {
public:
    int RunG2GLMLasso();

    int SetupPiEst(double rPriorM1, double rPriorM2);
    
    int SetFixKa(int r) {
      FixKa = (int) r;
      return(1);
    }
    void deletepIter() { if (pIter != NULL) {Free(pIter); pIter=NULL;}}   
private:
    int *y01;  // Vector of 0 or 1 responses (length n)
    double *X; // Matrix of Covariates (length n * p) 
    int n;     //Number of Samples
    int p;     // Number of Covariates
    double *Betas;  //  Covariate Loadings, length n
    double *pBeta0; //  Pointer to Intercept Beta0 value
    double *OnGammas; //  Weight vector of Gamma weights
    double *pBeta0Records; double *pRecordBetas; 
    double *pRecordBBHats;
    int UpdateBBHats();
    int UpdateGammas();
    int *pIter; //pIter, pointer to integer for what spot algorithm is on
    int RecordsNow(int tt1);
    int *GiveCDOIter;
    double PriorM1, PriorM2;
    int UpdatePiEst();
    
public :    
  GLMBinomObject  *GLMO;
    //CoordinateDescentObject *CDO;
  TwoGLMLassoObject(int rn, int *ry01, double *rX, int rp, 
    double *rBetas,  double *rpBeta0,  double *rOnGammas,
    double *rBBHats, double *rpRecordBetas,
    double *rpBeta0Records, double *rpRecordBBHats,
    double*rpBeta0prev,double*rpBetasPrev,
    double rppiuse, double *rPiRecVec,
    int rTotalRuns, double rsigmaSq,
    double *rSigmaRecVec, int rInitKKs, double rInverseGammaConstant, 
    int PrintFlag, double *rLambdaDK, double *rLambdaAK, 
    int *rOrderSeq, double *rnProbWeights, double *rnlogitCurrentProb,
    double *rnlogitPrevProb, double *rZ,
    double MaxCDOEpsilon, int MaxCDORuns, int *rpIter, int InMaximumAllocation)
        //ConvergentLARS((int)rn, (double *)ry01, (double *)rX, (int)p, 
        //               (double*)rInitBetas, 
        //              double rStlambdaD, double rStlambdaA,
        //              double rsigmaNoiseSq, double rMultDConst, double rMultAConst,
        //              int rStopAInt, double rppiuse, int rTotalRuns,
        //              int rInitKKs,
        //              double InverseGammaConstant, 
        //              int PrintFlag)
         {
  FixKa =-1;
  IPLLF = PrintFlag;
	OnLARS = NULL;  OnGammas = rOnGammas;
	LambdaDK = rLambdaDK;  RealXXS = NULL; yys = NULL;
	LambdaAK = rLambdaAK; OrderSeq = rOrderSeq; OnMultiCons = NULL;
	OnBetas = rBetas; BBOn1 = rBBHats; OnMultiCons = NULL;
	BBPrev = NULL; BetasKeep = NULL; 
  PriorM1 = 0; PriorM2 = 0;

  pBeta0Records = NULL;
  BBHatsKeep = NULL; 
  pRecordBetas = NULL;
        
  BBHatsNewStep = NULL;
	SigmaRecVec = NULL; PiRecVec = NULL; 
  CurrentPostPBj = NULL;  RecordPostPBj = NULL;  CurrentConfidenceInts = NULL;
  RecordConfidenceInts = NULL; CurGammaA = NULL; CurGammaD = NULL; DiagXTX=NULL;                  
	TDFNu = -1; iiWeights = NULL; Lambda3Seq = NULL; SigmaMultVec = NULL;
	if (rpRecordBBHats != NULL && rpRecordBBHats[0] >= 0) {
    pRecordBBHats = rpRecordBBHats; BBHatsKeep = rpRecordBBHats;
  } else {
    pRecordBBHats = NULL;  BBHatsKeep = NULL;
  }
       if (rpRecordBetas != NULL && rpRecordBetas[0] >= 0) {
           BetasKeep = rpRecordBetas;  pRecordBetas = rpRecordBetas;    
       }  else {
          BetasKeep  = NULL;  pRecordBetas=NULL;
       }    
       if (rpBeta0Records != NULL && rpBeta0Records[0] >= 0) { 
         pBeta0Records = rpBeta0Records;
       } else {
         pBeta0Records = NULL;
       }
       	    
      yys = rZ; 
	    pIter= rpIter;  // pIter is just a pointer to return what iteration
	                    // the algorithm stops upon
	    if (PrintFlag >= 1) {
        Rprintf("Declared pIter: pIter = %d", pIter[0]); R_FlushConsole();
      }
	    GiveCDOIter = (int*) Calloc(2, int);
	    n = rn; p = rp; y01 = ry01; X = rX;
	    Betas  = rBetas; pBeta0  = rpBeta0; 
      if (PrintFlag >= 0) { IPLLF = PrintFlag; } else {IPLLF = -1;} // No Printing
      SuccessFlag = 1;    // Flag triggers to -1 when TwoLasso Breaks
      ChoiceCDOOnLARSFlag = 2;   // This Will be Coordinate Descent Based TwoLasso
      if (IPLLF > 0) {
        Rprintf("G2Lasso Setup, about to Setting Up GLMO\n"); R_FlushConsole();
      }
      if (IPLLF > 0) {
        Rprintf("Note rGammas[0] =%.4f\n", rOnGammas[0]); R_FlushConsole();
        Rprintf("And rGammas[%d-1] = %4f\n", p,rOnGammas[p-1]); R_FlushConsole();
      }
      if (IPLLF > 0) {
        Rprintf("Note BBOn1[0] =%.4f\n", BBOn1[0]); R_FlushConsole();
        Rprintf("And BBOn1[%d-1] = %4f\n", p, BBOn1[p-1]); R_FlushConsole();
      }
      GLMO  = new GLMBinomObject(n, p, rX, (int *)ry01, 
                  rBetas, rpBeta0, rpBeta0prev, rpBetasPrev,
                  StlambdaA/2,
                  (double*) NULL, (double*) NULL,   //pRecordBetas, rpBeta0Records,                    
                  rOnGammas, rInitKKs, rInverseGammaConstant, 
                  rnProbWeights, rnlogitCurrentProb, rnlogitPrevProb,rZ,
                  PrintFlag-6, MaxCDOEpsilon, MaxCDORuns, GiveCDOIter,
                  InMaximumAllocation); 
      CDO = GLMO; 
      OnLARS = NULL;
      //->pBeta0prev = rpBeta0prev;    
       //////////////////  Now Insert Known Quantities
	     TotalRuns = rTotalRuns;
	     tt1 = 0; tt2 = 0;
	     SuccessFlag = 1;
	     OnLARS = NULL;  LambdaSeqF = 1;
	     FactorGammaDenominator = rInverseGammaConstant; NumEMConv = 0;
	     	                      
       //CDO = (CoordinateDescentObject*) GLMO; 
       this->CDO = NULL;
       
       OnLARS = NULL;  ppiuse = rppiuse; ChoiceCDOOnLARSFlag = 2;
	     TotalRuns = rTotalRuns;
	     SuccessFlag = 1;
	     OnLARS = NULL;
	     RealXXS = NULL;
	     OnNus = NULL;
	     OnMultiCons = NULL;
	     BBHatsNewStep = NULL;
	     SigmaRecVec = NULL; PiRecVec = NULL; CurrentPostPBj = NULL;
       RecordPostPBj = NULL;    CurrentConfidenceInts = NULL;
       RecordConfidenceInts = NULL; 
       NusKeep = NULL;
       
       CurGammaA = NULL; CurGammaD = NULL; DiagXTX=NULL;
       NusKeep = NULL; m1 = 0; m2 = 0;
        
       ChoiceCDOOnLARSFlag = 2; LambdaSeqF = 1;
       FactorGammaDenominator = rInverseGammaConstant;
       NumEMConv = 0;
          
       tt1 = 0; tt2 = 0;
           
	     //IPLLF = 2;

	    if (IPLLF > 0) {
		     Rprintf((char*)"Generating GLMTwoLasso\n");
		     R_FlushConsole();
	    }  
	    GAll = NULL;     

	 
      kLen =rp;
		  NusKeep = NULL;
	  
	    //BBPrev = (double *)Calloc(kLen+1, double);
	    sigmaNoiseSq = (double) rsigmaSq;
	    StlambdaD = 0.0;
	    StlambdaA = 0.0;
	    OnlambdaD = rLambdaDK[0];
	    OnlambdaA = rLambdaAK[0];
	    //MultDConst = (double) rMultDConst;
	    //MultAConst = (double) rMultAConst;
		  //StopAInt = rStopAInt;
		  ppiuse = (double) rppiuse;

     if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Survived CDO, now setup Gammas\n");
				     R_FlushConsole(); R_ProcessEvents();
		 }	    	
     SuccessFlag = GLMO->SetUpGammas(rOnGammas);
     if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Survived CDO, setup SigmaRecVec/PiRecVec\n");
				     R_FlushConsole(); R_ProcessEvents();
		 }	
                  
	   BBHatsNewStep = NULL;	
	   SigmaRecVec = NULL;
	   PiRecVec = NULL;    
  
	   CurrentPostPBj = NULL;		     	      
		 RecordPostPBj = NULL;		     	     


     OnNus = NULL;GAll = NULL;  RealXXS = NULL;
		 CurGammaA = NULL;
		 CurGammaD = NULL;
	  

      if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Deciding to make DiagXTX?\n");
				     R_FlushConsole(); R_ProcessEvents();
	    }			     	
  
	     if (IPLLF > 1) { 
				     Rprintf((char*)"new TwoLassoGLM: Made it to end of loading\n");
				     R_FlushConsole(); R_ProcessEvents();
	     }	
	     ChoiceCDOOnLARSFlag = 2; OnLARS = NULL;		
	     //IPLLF = 5;
	     //BBprev = (double*) Calloc(p+1, double);
	     //if (BBprev == NULL) {
       //  Rprintf("Error, BBprev set to NULL\n");
       //  SuccessLoad = -1;
       //}
       BBPrev = NULL;
       PiRecVec=rPiRecVec;
       SigmaRecVec = rSigmaRecVec;
   }
   
  public : 
   ~TwoGLMLassoObject() {
     //if(BetaPrev != NULL) {Free(BetaPrev); BetaPrev=NULL;}
     IPLLF = 0;
     if (IPLLF > 1) {
       Rprintf("G2LMO -> Deleting Items\n");
       R_FlushConsole();
     }
     if (PiRecVec != NULL) {PiRecVec = NULL;}
     if (IPLLF > 2) {
       Rprintf("G2LMO -> Deleting GiveCDOIter\n");
       R_FlushConsole();
     }
     if (GiveCDOIter != NULL) {
        GLMO->deletepIter();
        GiveCDOIter=NULL;
     }
     if (IPLLF > 2) {
       Rprintf("G2LMO-> Deleting pIter\n"); R_FlushConsole();
     }
     if (pIter != NULL) {
       Free(pIter);
     }
     if (IPLLF > 2) {
       Rprintf("G2MLO -> Deleting BBPrev\n");
       R_FlushConsole();
     }
     if (BBPrev != NULL) {
         Free(BBPrev);  BBPrev = NULL;
     }
     if (IPLLF > 2) {
      Rprintf("G2LMO -> Deleting GLMO\n");
       R_FlushConsole();
       GLMO->PrintFlagGLMB = 6;
     }
     delete(GLMO); GLMO = NULL;  CDO = NULL;
     if (IPLLF > 2) {
      Rprintf("G2LMO -> Finished Deleting\n");
       R_FlushConsole();
     }
     yys = NULL; RealXXS = NULL; BBOn1 = NULL; OnBetas=NULL;
     LambdaDK = NULL; LambdaAK = NULL; OrderSeq = NULL;
     SigmaRecVec = NULL; BetasKeep = NULL; 
     BBHatsKeep = NULL;  BetasKeep = NULL;
     this->CDO = NULL;
     if (this->CDO != NULL) {
        Rprintf("We're in for a load of shit since CDO is not NULL!\n");
        R_FlushConsole();
     }
     return;
   }


};


