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
  #define RCPPH 0
#endif

#ifndef RDYNLOADH
  #include <R_ext/Rdynload.h>
  #define RDYNLOADH 0
#endif

#undef isinf
#define isinf(x) (((x) > 999999999999.9) || ((x) < -999999999999.9))
#ifndef INCIDENTALOP2009DD
  #include "IncidentalOp2009.h"
  #define INCIDENTALOP2009DD 0 
#endif
#ifndef MyMatrixOp2009DD
  #include "MyMatrixOp2009.h"
  #define MyMatrixOp2009DD 0
#endif

#ifndef CoordinateDescentLassoDD
  #include "CoordinateDescentLasso2014.h"
  #define CoordinateDescentLassoDD 0
#endif
#ifndef LARSOBJECT2009DD
  #include "LarsObject2009.h"
  #define LARSOBJECT2009DD 0
#endif 
 
#ifndef  TwoLarsObject2009DD        
  #include "TwoLarsObject2009.h"
  #define TwoLarsObject2009DD 0
#endif

#ifndef CCEMRIDGE2009
   #include "CC_EMRIDGE2009.h"
   #define CCEMRIDGE2009 0
#endif

#ifndef TWOGLMLASSODD
   #include "TwoGLMLasso.h"
   #define TWOGLMLASSODD 0
#endif

#ifndef HorseShoeCDef
  #include "HorseShoeC.h"
  #define HorseShoeCDef 1
#endif 


#ifndef NEGLASSOH
  #include "NEGLasso.h" 
  #define NEGLASSOH 1
#endif
#ifndef TWOLASSO2011
  #include "TwoLassoObject2014.h"
  #define TWOLASSO2011 0
#endif


#ifndef Win32
  void R_ProcessEvents() {
     return;
  }
#endif

#undef isinf
#define isinf(x) (((x) > 999999999999.9) || ((x) < -999999999999.9))
//#define isinf
//#endif
/*





#ifndef DEFCONST 
   #define DEFCONST 2.5
#endif 





#ifndef CovMatrixEstimationDD
  #include "CovMatrixEstimation.h"
  #define CovMatrixEstimationDD 0
#endif

*/



extern "C" {

void CoordinateDescentLassoShell2009(int *pNLen, int *pkLen,  double *XXX, double *YYY,
  double *OnBeta, double *pOnGamma, double *OnRecordBetas,
  double *FinalBeta, int *pTotalLoops, double *pMaxEpsilon, double *OnGammas,
  int *pInitKKs, double *pInverseGammaConstant, double *WLSiiWeights,
  int *pPrintFlag) {
     int PrintFlag =0;
     if (pPrintFlag == NULL) { 
       Rprintf("pPrintFlag is NULL\n;"); R_FlushConsole(); return;
     }
     if (pPrintFlag[0] > 0) {PrintFlag = pPrintFlag[0];}
     if (PrintFlag > 0) {
       Rprintf("CoordinateDescentLassoShell2009: Entering\n");
     }
     if (pInitKKs != NULL && pInitKKs[0] > 0 && pInitKKs[0] > pkLen[0]) {
        Rprintf("CoordinateDescentLassoShell2009: ");
        Rprintf("pInitKKs[0] = %d, pkLen[0] = %d", pInitKKs[0], pkLen[0]);
        pInitKKs[0] = pkLen[0];
     }
     if (pInitKKs == NULL || pInitKKs[0] <= 0) {
       if (PrintFlag > 0) {
         Rprintf("CoordinateDescentLassoShell2009: Doing IK = %d <0 Version\n",
           pInitKKs[0]);
       }
       CoordinateDescentLasso(pNLen[0], pkLen[0],  XXX, YYY,
         OnBeta, pOnGamma[0], OnRecordBetas,
         FinalBeta, pTotalLoops, pMaxEpsilon[0], OnGammas, -1, 0, WLSiiWeights,
         PrintFlag);    
     }    else {
       if (PrintFlag > 0) {
         Rprintf("CoordinateDescentLassoShell2009: Doing IK = %d > 0 Version\n",
           pInitKKs[0]);
       }
       CoordinateDescentLasso(pNLen[0], pkLen[0],  XXX, YYY,
         OnBeta, pOnGamma[0], OnRecordBetas,
         FinalBeta, pTotalLoops, pMaxEpsilon[0], OnGammas, (int) pInitKKs[0],
         pInverseGammaConstant[0], WLSiiWeights,
         PrintFlag);
     }
}
   ////////////////////////////////////////////////////////////////////////////
   //  LarsLassAlgorithmBI
   //       Calls LarsLasso Algorithm  in LarsObject2009.cc
   void LarsLassoAlgorithmShell2( double * yys, int *pNLen, double *xxs, 
                      int *pkLen,  double *InitBetas, double *OldBetas, 
                      double *plambda, double *OutBetas, 
                      double *BetasKeep, int *JKeep, double *PenaltyKeep,
                      double *GammasKeep, double *gFKeep, int *pGFlag, 
                      int *pWeightLen,
                      double *Weights, int *pFinaltt) {
	    LarsLassoAlgorithmBI( yys, pNLen[0], xxs, pkLen[0], 
                      InitBetas, 
                      OldBetas, plambda[0],OutBetas, BetasKeep, JKeep, PenaltyKeep,
                      GammasKeep, gFKeep, pGFlag[0], pWeightLen[0],
                      Weights, pFinaltt);                  	                      
  }
   /////////////////////////////////////////////////////////////////////////////
   //    CalculateGFTAllShell
   //
   //
   //
   void CalculateGFTAllShell(double  *RetGamm, double *TFLV2, double *GFT0, 
          double *GFT1, int *pOnL, int *pTotk, 
          int *ActiveOns, double *OnBetas, double *BetaOlds, int *pNLen, double *OnResid, 
          double *wa, int *sjA, double *ua, double *pv1, double *pv2, double *pOverallMax)  {
	 //Rprintf( (char *) "About to Create Memory Buggers \n");
	 //R_FlushConsole();
	 double * PODOnBetas = (double *) Calloc(pTotk[0]+2, double );
	 double *PODBetaOlds = (double *) Calloc(pTotk[0]+2, double);
	 //Rprintf( (char *) "   To Create Memory OnResid/ua \n");
	 //R_FlushConsole();	 
	 double *PODOnResid = (double *) Calloc(pNLen[0] +2, double);
     double *PODua     = (double *) Calloc(pNLen[0] + 2, double);	 
	 //Rprintf( (char *) "   To Create Memory wa \n");
	 //R_FlushConsole();	 	 
	 double *PODwa = (double *) Calloc(pOnL[0] +2, double);
	 long  double OverallMax = pOverallMax[0];
	 int ii;
	 //Rprintf( (char *) "Filling in PODOnBetas/BetaOlds\n");
	 //R_FlushConsole();	 
	 for (ii = 0; ii < pTotk[0]; ii++) {
		   PODOnBetas[ii] = (double) OnBetas[ii];
		   PODBetaOlds[ii] = (double) BetaOlds[ii];
     }
     //Rprintf( (char *) "Filling in OnResid/ua\n");
	 //R_FlushConsole();	      
     for (ii = 0; ii < pNLen[0] ;ii++ ) {
	       PODOnResid[ii] = (double) OnResid[ii];
	       PODua[ii] = (double) ua[ii];
     }
     //Rprintf( (char *) "Filling in wa/ua\n");
	 //R_FlushConsole();	      
     for (ii = 0; ii < pOnL[0]; ii++) { PODwa[ii] = (double) wa[ii];}
	 //return;
	 RetGamm[0] = (double) CalculateGFTAll(TFLV2, GFT0, GFT1, pOnL[0],  pTotk[0],
          ActiveOns, PODOnBetas, PODBetaOlds, pNLen[0], PODOnResid, 
          PODwa, sjA, PODua, pv1[0], pv2[0], OverallMax);
     Free(PODOnBetas);
     Free(PODBetaOlds);
     Free(PODOnResid);  
     Free(PODwa);     
     Free(PODua);
     return;                	                      		
        	             
}

void  LarsConvergencyShell(double * yys, int *pNLen, double *xxs, int *pkLen, 
  double *InitBetas, 
  double *OldBetas, double *OutBetas, double *BetasKeep,
  int *pNumTotalRuns, int *pNumEMConv, double *pMultEMCons,
  double *pStlambdaD, double *pStlambdaA, 
  double *plambdaDMultC, double *plambdaAmultC,
  int *plambdaAmultStop, double *pppiuse,
  double *psigmaNoiseSq,
  double *BBHatsKeep, double *NusKeep, 
  double *LambdaDK, double *LambdaAK,
  double *RecordPostPBj,  double *pDConfidenceIntv,           
  double *RecordConfidenceInts, int *pLAKSeq, 
  int* OrderSeq, int* pNumCDOConv, double *pCDOEpsilon,
  int *pInitKKs, double *pInverseGammaConstant,
  int *pFixKa, double *PiRecVec, double *iiWeights, 
  double *pTDFNu, int *pPrintFlag,
  double *pSigmaVec,
  double *pLambda3, double *pLambda3Seq)   {
	int LCBFlag = 0;
  if (pInitKKs != NULL && pInitKKs[0] > 0 && pInitKKs[0] > pkLen[0]) {
    Rprintf("LarsConvergenyShell: ");
    Rprintf("pInitKKs[0] = %d, pkLen[0] = %d", pInitKKs[0], pkLen[0]);
    pInitKKs[0] = pkLen[0];
  }  
	//Rprintf((char*) "LarsConvergencyShell about start \n");
	//R_FlushConsole(); R_ProcessEvents(); 
	  
////////////////////////////////////////////////////////////////////////////
//////////  The Test operation tests out your variables to ensure
//////////   that you've assigned enough data space.   Uncomment for evaluation. 
//	  LCBFlag = LarsConvergencyBITest(yys, pNLen[0], xxs, pkLen[0], 
//                      InitBetas, OldBetas, OutBetas, BetasKeep,
//                      pNumTotalRuns[0], pNumEMConv[0], pMultEMCons[0],
//                      pStlambdaD[0], pStlambdaA[0], 
//                      plambdaDMultC[0], plambdaAmultC[0],
//                      plambdaAmultStop[0], pppiuse[0],
//                      psigmaNoiseSq[0], BBHatsKeep, NusKeep, LambdaDK, LambdaAK,
//                      RecordPostPBj, pDConfidenceIntv[0], RecordConfidenceInts,
//                      pLAKSeq[0], OrderSeq, pNumCDOConv[0], pCDOEpsilon[0],
//                        pInitKKs[0],
//                      pInverseGammaConstant[0],
//                      pFixKa[0], PiRecVec, iiWeights, pTDFNu[0],
//                      pPrintFlag[0]
//                      );    
                     // return; 	         
    double Lambda3;
    if (pLambda3 != NULL && pLambda3[0] > 0) {
       Lambda3 = pLambda3[0];
    }  else {
       Lambda3 = 0;
    }
    //Rprintf("A start pSigmaVec[0] = %f \n", pSigmaVec[0]); R_FlushConsole();
	  LCBFlag = LarsConvergencyBI(yys, pNLen[0], xxs, pkLen[0], 
      InitBetas, OldBetas, OutBetas, BetasKeep,
      pNumTotalRuns[0], pNumEMConv[0], pMultEMCons[0],
      pStlambdaD[0], pStlambdaA[0], 
      plambdaDMultC[0], plambdaAmultC[0],
      plambdaAmultStop[0], pppiuse[0],
      psigmaNoiseSq[0],
      BBHatsKeep, NusKeep, LambdaDK, LambdaAK,
      RecordPostPBj, pDConfidenceIntv[0], RecordConfidenceInts,
      pLAKSeq[0], OrderSeq, pNumCDOConv[0], pCDOEpsilon[0],
      pInitKKs[0],  pInverseGammaConstant[0], pFixKa[0], PiRecVec,
      iiWeights, pTDFNu[0], pPrintFlag[0],
      pSigmaVec, Lambda3, pLambda3Seq
                      );     
      if (LCBFlag < 0) {
	      Rprintf((char*)"LarsConvergency Fail, LCBFlag = %d\n", LCBFlag);
	      InitBetas[0] = -999;
	      R_FlushConsole();
      }                
	                      
   }
void  LarsConvergencyPIEstShell(double * yys, int *pNLen, double *xxs, int *pkLen, 
                      double *InitBetas, 
                      double *OldBetas, double *OutBetas, double *BetasKeep,
                      int *pNumTotalRuns, int *pNumEMConv, double *pMultEMCons,
                      double *pStlambdaD, double *pStlambdaA, 
                      double *plambdaDMultC, double *plambdaAmultC,
                      int *plambdaAmultStop,
                      double *BBHatsKeep, double *NusKeep, 
                      double *LambdaDK, double *LambdaAK, 
                      int *pLAKSeq, int *OrderSeq,
                      double *RecordPostPBj, double *pDConfidenceIntv, 
                      double *RecordConfidenceInts,
                      double *pm1, double *pm2, int *pNumEConv,
                      double *pSigmaEta, double *pSigmaBarEst,
                      double *PiRecVec, double *SigmaRecVec,
                      int *pNumCDOConv, double *pCDOEpsilon,
                      int *pInitKKs,
                      double *pInverseGammaConstant,
                      double *iiWeights, double *pTDFNu, int *pPrintFlag,
                      double *pSigmaMultVec) {
     if (pInitKKs != NULL && pInitKKs[0] > 0 && pInitKKs[0] > pkLen[0]) {
        Rprintf("LarsConvergencyPIEstShell: ");
        Rprintf("pInitKKs[0] = %d, pkLen[0] = %d", pInitKKs[0], pkLen[0]);
        pInitKKs[0] = pkLen[0];
     }
	  int LCBFlag = 0;  
	  //Rprintf((char*)"m1 = %.4f, m2 = %.4f, SigmaBarEst = %.4f \n", 
	  //    pm1[0], pm2[0], pSigmaBarEst[0]);
	  //    R_FlushConsole();                 
	  LCBFlag = LarsConvergencyPIRECBI(yys, pNLen[0], xxs, pkLen[0], 
                      InitBetas, OldBetas, OutBetas, BetasKeep,
                      pNumTotalRuns[0], pNumEMConv[0], pMultEMCons[0],
                      pStlambdaD[0], pStlambdaA[0], 
                      plambdaDMultC[0], plambdaAmultC[0],
                      plambdaAmultStop[0], BBHatsKeep, NusKeep, LambdaDK, LambdaAK,
                      pLAKSeq[0], OrderSeq,
                      RecordPostPBj, pDConfidenceIntv[0], RecordConfidenceInts,
                      pm1[0], pm2[0], (int) pNumEConv[0],
                      pSigmaEta[0], pSigmaBarEst[0],
                      PiRecVec, SigmaRecVec, pNumCDOConv[0], pCDOEpsilon[0],
                      pInitKKs[0],
                      pInverseGammaConstant[0],
                      iiWeights, pTDFNu[0], pPrintFlag[0],
                      pSigmaMultVec);     
      if (LCBFlag < 0) {
	      Rprintf((char*)"LarsConvergency Fail, LCBFlag = %d\n", LCBFlag);
	      InitBetas[0] = -999;
	      R_FlushConsole();
      }                
	                      
   }   

//void CoordinateDescentLassoForce2Shell(int *pNLen, int *pkLen,  double *XXX, double *YYY,
//             double *OnBeta, double *pOnGamma, double *OnRecordBetas,
//             double *FinalBeta, int *pTotalLoops, double *pMaxEpsilon, double *OnGammas, 
//             int *pInitKKs) {
//   CoordinateDescentLassoForce2(pNLen[0], pkLen[0],  XXX, YYY,
//              OnBeta, pOnGamma[0], OnRecordBetas,
//              FinalBeta, pTotalLoops, pMaxEpsilon[0], OnGammas, pInitKKs[0]);             	             
//}	

////////////////////////////////////////////////////////////////////////////////
//   EMRIDGE Shells
//      There are Multiple EM2Ridge shells, all with different speeds and 
//         accuracies for computation.  EM2Ridge is covered in Lenarcic 2009
//         thesis.  It did not seem to be competitive against the EM2Lasso
//         as a selection algorithm, and it is nearly identical to the
//         George 1994 Gibbs Sampling selection framework.
//
   void EMRIDGEShell (int *pNLen, int *pkLen, double *YData, double *XData, 
         int *pChangeNSteps, int *pConvergenceNSteps, double *pSigmaSqNoise,
         double *ppuse,
         double *psigma1sqStart, double *psigma2sqStart, 
         double *psigma1Mult, double *psigma2Mult, int *psigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut) {
     int RFlag = EMRIDGE(pNLen[0], pkLen[0], YData, XData, 
         pChangeNSteps[0], pConvergenceNSteps[0], pSigmaSqNoise[0],
         ppuse[0],
         psigma1sqStart[0], psigma2sqStart[0], 
         psigma1Mult[0], psigma2Mult[0], psigma2MultStop[0], 
         BStartBBs, BStartBetas,
         BetaPartBetas,
         BEndBetas, BBHatRec, BBHatOut);       
     if (RFlag < 0) {
	     Rprintf((char*)"EMRIDGE Fail, RFlag = %d \n", RFlag);
	     R_FlushConsole(); R_ProcessEvents();
	     BStartBetas[0] = -999;
     }
}


void EMRIDGEShell2 (int *pNLen, int *pkLen, double *YData, double *XData, 
         int *pChangeNSteps, int *pConvergenceNSteps, double *pSigmaSqNoise,
         double *ppuse,
         double *psigma1sqStart, double *psigma2sqStart, 
         double *psigma1Mult, double *psigma2Mult, int *psigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut, 
         double *plogp, double *plog1mp, double *pSigmaNoiseInv) {
	   //Rprintf("ConvergentEMShell2");
	   //Rprintf("plogp[0] = %.4f, plog1mp[0] = %.4f, pSigmaNoiseInv[0] = %.4f, BBHatRec[0] = %.4f\n",
	   //      (double) plogp[0],(double) plog1mp[0], (double) pSigmaNoiseInv[0], 
	   //      (double) BBHatRec[0]);
	   //R_FlushConsole();
	   //R_ProcessEvents(); 
     int RFlag = EMRIDGE (pNLen[0], pkLen[0], YData, XData, 
         pChangeNSteps[0], pConvergenceNSteps[0], pSigmaSqNoise[0],
         ppuse[0],
         psigma1sqStart[0], psigma2sqStart[0], 
         psigma1Mult[0], psigma2Mult[0], psigma2MultStop[0], 
         BStartBBs, BStartBetas,
         BetaPartBetas,
         BEndBetas, BBHatRec, BBHatOut, 
         plogp[0], plog1mp[0], pSigmaNoiseInv[0]);    
     if (RFlag < 0) {
	     Rprintf((char*)"EMRIDGE2 Fail, RFlag = %d \n", RFlag);
	     R_FlushConsole(); R_ProcessEvents();
	     BStartBetas[0] = -999;
     }
}



void EMDYMANICRIDGEShell2 (int *pNLen, int *pkLen, double *YData, double *XData, 
         int *pTotalRuns, int *pNumEMConv, double *pSigmaSqNoise,
         double *ppuse,
         double *pStTauSqD, double *pStTauSqA, 
         double *pTauDMult, double *pTauAMult, int *pTauAMultStop,
         double *TauASeq, double *TauDSeq, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut, 
         double *plogp, double *plog1mp, double *pSigmaNoiseInv, 
         int *pNumCDOConv, double *pCDOEpsilon, int *OrderSeq, int *CountDescentOperations) {
	//
	// Rprintf("EMDYMANICRIDGEShell2"); R_FlushConsole();
	//  Rprintf("plogp[0] = %.4f, plog1mp[0] = %.4f, pSigmaNoiseInv[0] = %.4f\n",
	//         (double) plogp[0],(double) plog1mp[0], (double) pSigmaNoiseInv[0]);
	//   R_FlushConsole();
	//   R_ProcessEvents(); 
	// Rprintf("ppuse[0] = %.4f\n", ppuse); R_FlushConsole(); R_ProcessEvents();
	// Rprintf("BBHatRec[0] = %f\n", BBHatRec[0]); R_FlushConsole();
	// if (BBHatRec != NULL && BBHatRec[0] >= 0) {
	//   Rprintf("BBHatRec[TotalRuns * kLen -1 = %d] = %f \n", 
	//       pTotalRuns[0] * pkLen[0] - 1, BBHatRec[ pTotalRuns[0] * pkLen[0] - 1] );
	//   R_FlushConsole(); R_ProcessEvents();	 
  //   }
  //   Rprintf("BBHatRec[0] = %f\n", BetaPartBetas[0]); R_FlushConsole();
	// if (BetaPartBetas != NULL && BetaPartBetas[0] >= 0) {
	//   Rprintf("BBHatRec[TotalRuns * kLen -1 = %d] = %f \n", 
	//        pTotalRuns[0] * pkLen[0] - 1, BetaPartBetas[ pTotalRuns[0] * pkLen[0] - 1] );
	//   R_FlushConsole(); R_ProcessEvents();	 
  //   }     
	// Rprintf("BBHatOut[0] = %f, BBHatOut[%d] = %f\n", BBHatOut[0], pkLen[0]-1, BBHatOut[pkLen[0]-1]);
	// R_FlushConsole(); R_ProcessEvents();
	// Rprintf("plogp[0] = %.4f, plog1mp = %.5f, pSigmaNoiseInv[0] = %f\n",
	//      plogp[0], plog1mp[0], pSigmaNoiseInv[0]); R_FlushConsole(); R_ProcessEvents();
	// Rprintf("pNumCDOConv[0] = %d, pCDOEpsilon[0] = %10f\n", pNumCDOConv[0], pCDOEpsilon[0]);
	// R_FlushConsole(); R_ProcessEvents();
	// Rprintf("OrderSeq[0] = %d and OrderSeq[%d] = %d \n", OrderSeq[0], 
	//                  pTotalRuns[0], OrderSeq[pTotalRuns[0]-1]);
	// R_FlushConsole(); R_ProcessEvents();
	// Rprintf("CountDescentOperations[0] = %d \n", CountDescentOperations[0]);
	// R_FlushConsole(); R_ProcessEvents();
	//  Rprintf("CountDescentOperations[TotalRuns-1 = %d] = %d \n",
	//     pTotalRuns[0] -1, CountDescentOperations[pTotalRuns[0]-1]);
	// R_FlushConsole(); R_ProcessEvents();
	//
     int RFlag =  EMDYNAMICRIDGE(pNLen[0], pkLen[0], YData, XData, 
         pTotalRuns[0], pNumEMConv[0], pSigmaSqNoise[0],
         ppuse[0],
         pStTauSqD[0], pStTauSqA[0], 
         pTauDMult[0], pTauAMult[0], pTauAMultStop[0], 
         TauASeq, TauDSeq,
         BStartBBs, BStartBetas,
         BetaPartBetas,
         BEndBetas, BBHatRec, BBHatOut, 
         plogp[0], plog1mp[0], pSigmaNoiseInv[0],
         pNumCDOConv[0], pCDOEpsilon[0], OrderSeq, CountDescentOperations);  
     //Rprintf("EMRIDGE2:  EMDYNAMICRIDGE finished returned RFlag = %d \n", RFlag); 
     //R_FlushConsole(); R_ProcessEvents();
     
     if (RFlag < 0) {
	     Rprintf((char*)"EMRIDGE2 Fail, RFlag = %d \n", RFlag);
	     R_FlushConsole(); R_ProcessEvents();
	     BStartBetas[0] = -999;
     }
     return;
}
//NLen = as.integer(NLen), kLen=as.integer(kLen),
//         YY = as.double(YY), XX = as.double(XX), 
//         ChangeNSteps = as.integer(n1), ConvergenceNSteps = as.integer(n2), 
//         ConvergeCloseEnough = as.double(.000001), 
//         SigmaSqNoise = as.double(SigmaSqNoise),
//         ppiuse = as.double(ppiuse), 
//         tauDsqStart =as.double(tauDsqStart), 
//         tauAsqStart = as.double(tauAsqStart), 
//         tauDMult = as.double(tauDMult), tauAMult = as.double(tauAMult), 
//         tauAMultStop = as.integer(tauAMultStop), 
//         BStartBBs = as.double(BStartBBs),
//         BStartBetas = as.double (BStartBetas),
//         ReturnBetas = as.double(vector("numeric", kLen)),
//         BBHatsFinal = as.double(vector("numeric", kLen)),
//         logp = as.double(logp), log1mp = as.double(log1mp), 
//         SigmaNoiseInv = as.double(SigmaNoiseInv)
void EMRIDGENoRecShell (int *pNLen, int *pkLen, double *YData, double *XData, 
         int *pChangeNSteps, int *pConvergenceNSteps, 
         double *pConvergeCloseEnough, double *pSigmaSqNoise,
         double *ppuse,
         double *ptauDsqStart, double *ptauAsqStart, 
         double *ptauDMult, double *ptauAMult, int *ptauAMultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BEndBetas, double *BBHatOut,
         double *plogp, double *plog1mp, double *pSigmaNoiseInv) {
     int RFlag = EMRIDGENoRec (pNLen[0], pkLen[0], YData, XData, 
         pChangeNSteps[0], pConvergenceNSteps[0], 
         pConvergeCloseEnough[0], pSigmaSqNoise[0],
         ppuse[0],
         ptauDsqStart[0], ptauAsqStart[0], 
         ptauDMult[0], ptauAMult[0], ptauAMultStop[0], 
         BStartBBs, BStartBetas,
         BEndBetas,  BBHatOut,
         plogp[0], plog1mp[0], pSigmaNoiseInv[0]); 
     if (RFlag < 0) {
	     Rprintf((char*)"EMRIDGENoRec Fail, RFlag = %d \n", RFlag);
	     R_FlushConsole(); R_ProcessEvents();
	     BStartBetas[0] = -999;
     }          
     //Rprintf((char*) "EMRIDGENoRecShell, finished, all done \n");
     //R_FlushConsole(); R_ProcessEvents();       
}

void SimCor(int *pLen, double *prho, double *RTVec) {
  double sqtRho = sqrt(1-prho[0] * prho[0]);
   for (int ii = 1; ii < pLen[0]; ii++) {
       RTVec[ii] = prho[0] * RTVec[ii-1] +  sqtRho * RTVec[ii];
   }
   return;
}

 /*
void CovEstimatorRunShell(int *pKappaT, int *pNData, double *YY, 
    int *pTotalSteps, int *pStartOnK2, double *ppiA, double *pStLambdaA,
    double *pStLambdaD, double *pLambdaAMult, double *pLambdaDMult, 
    int *pLambdaAStop, int *pTotalRuns, int *pNumEStep, double *LambdaASeq,
    double *LambdaDSeq, int *OrderSeq, int *pRecordFlag,
    double *Theta22PriorAlpha, double *Theta22PriorNus,
    double *OutTheta, double *OutThetaRecords, 
    double *OutBBOn, double *BBOnRecords, double *OutW, double *WRecords,
    int *pNumCDOConv, double *pCDOEpsilon,
    double *pStbb, double *ptausqB, double *plambdaB, double *pOutbb, double *OnbbRecs,
    double *pMaxTolerance) {	  
	      
	    //CETestInputsFunction(pKappaT[0], pNData[0], YY,
	    //               pTotalSteps[0], pStartOnK2[0],
	    //               ppiA[0], pStLambdaA[0], pStLambdaD[0],
	    //               pLambdaAMult[0], pLambdaDMult[0], 
	    //              pLambdaAStop[0],
	    //               pTotalRuns[0], pNumEStep[0],
	    //               LambdaASeq, LambdaDSeq,
	    //               OrderSeq, pRecordFlag[0], Theta22PriorAlpha, 
	    //               Theta22PriorNus, OutTheta, OutThetaRecords,
	    //               OutBBOn, BBOnRecords, OutW, WRecords,
	    //               pNumCDOConv[0], pCDOEpsilon[0],
      //                 pStbb[0], ptausqB[0], plambdaB[0], pOutbb, OnbbRecs,
      //                 pMaxTolerance[0]);	 
	    CovEstimatorRun(pKappaT[0], pNData[0], YY,
	                   pTotalSteps[0], pStartOnK2[0],
	                   ppiA[0], pStLambdaA[0], pStLambdaD[0],
	                   pLambdaAMult[0], pLambdaDMult[0], 
	                   pLambdaAStop[0],
	                   pTotalRuns[0], pNumEStep[0],
	                   LambdaASeq, LambdaDSeq,
	                   OrderSeq, pRecordFlag[0], Theta22PriorAlpha, 
	                   Theta22PriorNus, OutTheta, OutThetaRecords,
	                   OutBBOn, BBOnRecords, OutW, WRecords,
	                   pNumCDOConv[0], pCDOEpsilon[0],
                       pStbb[0], ptausqB[0], plambdaB[0], pOutbb, OnbbRecs,
                       pMaxTolerance[0]);
        return;
}  

void CovEstimatorRunShell2(int *pKappaT, int *pNData, double *YY, 
    int *pTotalSteps, int *pStartOnK2, double *ppiA, double *pStLambdaA,
    double *pStLambdaD, double *pLambdaAMult, double *pLambdaDMult, 
    int *pLambdaAStop, int *pTotalRuns, int *pNumEStep, double *LambdaASeq,
    double *LambdaDSeq, int *OrderSeq, int *pRecordFlag,
    double *Theta22PriorAlpha, double *Theta22PriorNus,
    double *OutTheta, double *OutThetaRecords, 
    double *OutBBOn, double *BBOnRecords, double *OutW, double *WRecords, 
    int *pNumCDOConv, double *pCDOEpsilon, 
    double *pStbb, double *ptausqB, double *plambdaB, double *pOutbb, double *OnbbRecs,
    double *pMaxTolerance, 
    double *AlphaMatrix, double *BBOnAlpha,
    double *pSigDif) {
	    // Rprintf((char*)"CovEstimatorRunShell2 \n");
	    //  Rprintf("pStbb[0] = %f, ptausqB[0] = %f, pOutbb[0] = %f, OnbbRecs[0] = %f\n",
	    //     pStbb[0], ptausqB[0], pOutbb[0], OnbbRecs[0]);
	    //     R_FlushConsole(); 
	    
	    //CETestInputsFunction(pKappaT[0], pNData[0], YY,
	    //               pTotalSteps[0], pStartOnK2[0],
	    //               ppiA[0], pStLambdaA[0], pStLambdaD[0],
	    //               pLambdaAMult[0], pLambdaDMult[0], 
	    //               pLambdaAStop[0],
	    //               pTotalRuns[0], pNumEStep[0],
	    //               LambdaASeq, LambdaDSeq,
	    //               OrderSeq, pRecordFlag[0], Theta22PriorAlpha, 
	    //               Theta22PriorNus, OutTheta, OutThetaRecords,
	    //               OutBBOn, BBOnRecords, OutW, WRecords,
	    //               pNumCDOConv[0], pCDOEpsilon[0],
      //                 pStbb[0], ptausqB[0], plambdaB[0], pOutbb, OnbbRecs,
      //                 pMaxTolerance[0]);	 
	    CovEstimatorRun(pKappaT[0], pNData[0], YY,
	                   pTotalSteps[0], pStartOnK2[0],
	                   ppiA[0], pStLambdaA[0], pStLambdaD[0],
	                   pLambdaAMult[0], pLambdaDMult[0], 
	                   pLambdaAStop[0],
	                   pTotalRuns[0], pNumEStep[0],
	                   LambdaASeq, LambdaDSeq,
	                   OrderSeq, pRecordFlag[0], Theta22PriorAlpha, 
	                   Theta22PriorNus, OutTheta, OutThetaRecords,
	                   OutBBOn, BBOnRecords, OutW, WRecords,
	                   pNumCDOConv[0], pCDOEpsilon[0],
                       pStbb[0], ptausqB[0], plambdaB[0], pOutbb, OnbbRecs,
                       pMaxTolerance[0]);
        return;
}  

*/     	             

SEXP GLMLassoShell(SEXP sX, SEXP sy01, SEXP sOnBeta, SEXP sBeta0,
   SEXP sBeta0prev, SEXP sBetasPrev, SEXP sOnGamma,
   SEXP sRecordOnBetas,SEXP sBeta0Records,
   SEXP sOnGammas,  SEXP sInitKKs, SEXP sInverseGammaConstant,
   SEXP snProbWeights, SEXP slogitCurrentProb, 
   SEXP snlogitPrevProb, SEXP sz, 
   SEXP sPrintFlag, SEXP sMaxCDOEpsilon, SEXP sMaxCDORuns,
   SEXP sIter){
   SEXP returner = 
    GLMLasso( sX,  sy01,  sOnBeta, sBeta0,  sBeta0prev,  sBetasPrev,  sOnGamma,
    sRecordOnBetas, sBeta0Records, sOnGammas, sInitKKs,
    sInverseGammaConstant, snProbWeights,  slogitCurrentProb, 
    snlogitPrevProb,  sz,  sPrintFlag, sMaxCDOEpsilon,
    sMaxCDORuns, sIter);
    return(returner);
   }
   //SEXP GLMLasso(SEXP sX, SEXP sy01, SEXP sOnBeta, double sBeta0,
   //SEXP sBeta0prev, SEXP sBetaprev, SEXP sOnGamma,
   //SEXP sRecordOnBetas,SEXP sBeta0Records,
   //SEXP sOnGammas,
   //SEXP sInitKKs, SEXP sInverseGammaConstant,
   //SEXP snProbWeights, SEXP slogitCurrentProb, 
   //SEXP snlogitPrevProb, SEXP rZ, SEXP sPrintFlag,
   //SEXP sMaxCDOEpsilon, SEXP sMaxCDORuns);     

SEXP RunG2GLMLasso(SEXP sXX, SEXP sy01, SEXP sBetas, SEXP sBeta0,
                   SEXP sOnGammas, SEXP sBBHats, SEXP sRecordBetas,
                   SEXP sBeta0prev, SEXP sBetasprev,
                   SEXP sBeta0Records, SEXP sRecordBBHats,
                   SEXP sLambdaDK, SEXP sLambdaAK,
                   SEXP sppiuse, SEXP sPiRecVec, SEXP ssigmaSq,
                   SEXP sSigmaRecVec,
                   SEXP sOrderSeq, SEXP sProbWeights, SEXP snLogitCurrentProb,
                   SEXP snlogitPrevProb,
                   SEXP sZ, SEXP sInverseGammaConstant, 
                   SEXP sInitKKs, SEXP sPrintFlag,
                   SEXP sMaxCDOEpsilon, SEXP sMaxCDORuns,
                   SEXP sIter, SEXP sFixKa, SEXP sPriorS, SEXP sInMaximumAllocation)   {
   int PrintFlag = 0;
   if (Rf_isInteger(sPrintFlag)) {PrintFlag = (int) INTEGER(sPrintFlag)[0];} else {
      PrintFlag = (int) REAL(sPrintFlag)[0];}
   if (PrintFlag > 0) {
     Rprintf("TwoLarsOperator: Into the Shell! \n\n"); R_FlushConsole();
   }
   
   SEXP returner = RunG2GLMLasso2(  sXX,   sy01,   sBetas,   sBeta0,
                     sOnGammas,   sBBHats,   sRecordBetas,
                     sBeta0prev,   sBetasprev,
                     sBeta0Records,   sRecordBBHats,
                     sLambdaDK,   sLambdaAK,
                     sppiuse,   sPiRecVec, ssigmaSq, sSigmaRecVec,
                     sOrderSeq,   sProbWeights,   snLogitCurrentProb,
                     snlogitPrevProb,
                     sZ,   sInverseGammaConstant, 
                     sInitKKs,   sPrintFlag,
                     sMaxCDOEpsilon,   sMaxCDORuns,
                     sIter, sFixKa, sPriorS, sInMaximumAllocation);
                
   return(returner);                
   }
   
SEXP OutGetAllLambdaj(SEXP SlZVals, SEXP SdlStandard, SEXP SBetaj, SEXP StauSq,
     SEXP SFillerDens, SEXP SUnifDraw, SEXP SOutLambdaj, SEXP Sexpn2lz2) {
  GetAllLambdaj(SlZVals, SdlStandard, SBetaj, StauSq,
       SFillerDens, SUnifDraw, SOutLambdaj, Sexpn2lz2);
  return(SOutLambdaj);
}

SEXP OutGetManyLambdaj(SEXP SlZVals, SEXP SdlStandard, SEXP SBetaj, SEXP StauSq,
     SEXP SFillerDens, SEXP SSortedDraws, SEXP SOutLambdaj, SEXP Sexpn2lz2) {
  AllGetManyLambdaj(SlZVals, SdlStandard, SBetaj, StauSq,
       SFillerDens, SSortedDraws, SOutLambdaj, Sexpn2lz2);
  return(SOutLambdaj);
}

SEXP TwoLassoRegressionShell(
  SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
  SEXP rSXtX, SEXP rStt1, SEXP rStt2, 
  SEXP rSLambdaA, SEXP rSLambdaD, SEXP rSOrderSeq, 
  SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
  SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
  SEXP rSm, SEXP rSSigmaPrior,
  SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
  SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
  SEXP rSTDFNu, SEXP rSiiWeights, 
  SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
  SEXP rSVerbose, SEXP rGroupSexp) {
  return( TwoLassoRegression( rSn,  rSp,  rSyys,  rSxxs,  rSXtY,
     rSXtX,  rStt1, rStt2, rSLambdaA,  rSLambdaD,  rSOrderSeq, 
     rSBBOn1,  rSOnBeta,  rRecBBOn1,  rRecOnBeta,
     rSOnGammas,  rSPiA,  rSSigma,
     rSm,  rSSigmaPrior,
     rSInitKKs,  rSInverseGammaConstant, 
     rSCauchyEpsilon,  rSMaxCauchy, 
     rSTDFNu,  rSiiWeights,  
     rSL2ShrinkagePrior, rSL2ShrinkageRecords, rSVerbose,
     rGroupSexp) );    
}
SEXP TwoLassoCreateTwoLassoSexpShell(SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaA, SEXP rSLambdaD, SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,SEXP rSVerbose,
    SEXP rGroupSexp, SEXP RunFlag) {
return(TwoLassoCreateTwoLassoSexp( rSn,  rSp,  rSyys,  rSxxs,  rSXtY,
     rSXtX,  rStt1, rStt2, rSLambdaA,  rSLambdaD,  rSOrderSeq, 
     rSBBOn1,  rSOnBeta,  rRecBBOn1,  rRecOnBeta,
     rSOnGammas,  rSPiA,  rSSigma,
     rSm,  rSSigmaPrior,
     rSInitKKs,  rSInverseGammaConstant, 
     rSCauchyEpsilon,  rSMaxCauchy, 
     rSTDFNu,  rSiiWeights, 
     rSL2ShrinkagePrior, rSL2ShrinkageRecords, rSVerbose,
     rGroupSexp,  RunFlag));    
}
SEXP TwoLassoReRegressShell(SEXP sTLS, SEXP sRunFlag, SEXP newMvec, SEXP newSigmaPrior,
  SEXP newVerbose) {
  return(TwoLassoReRegress(sTLS, sRunFlag, newMvec, newSigmaPrior, newVerbose));
}

SEXP AlterVector(SEXP VectorOb, SEXP IntLoc, SEXP ChangePart) {
  int ii;                            
  int rIntLoc = 0; if (Rf_isReal(IntLoc)) { rIntLoc = (int) REAL(IntLoc)[0];}
  if (Rf_isInteger(IntLoc)) { rIntLoc = INTEGER(IntLoc)[0]; }
  if (Rf_isReal(VectorOb)) {
    if (Rf_isReal(ChangePart)) {
      //Rprintf("Vector Ob and Change Part  Both are REAL\n"); R_FlushConsole();
      for (ii = 0; ii < Rf_length(ChangePart); ii++) {
        if (ii + rIntLoc > Rf_length(VectorOb)) {
          Rf_error("AlterVector, Rf_length of VectorOb only %d but ii=%d, GL = %d", 
            Rf_length(VectorOb), ii, rIntLoc );
        }
        REAL(VectorOb)[ii + (int) rIntLoc] = 
          REAL(ChangePart)[ii];
      }
    }  else if (Rf_isInteger(ChangePart)) {
      //Rprintf("Vector Ob REAL and Change Part are INTEGER\n"); R_FlushConsole();
      for (ii = 0; ii < Rf_length(ChangePart); ii++) {
        if (ii + rIntLoc > Rf_length(VectorOb)) {
          Rf_error("AlterVector, Rf_length of VectorOb only %d but ii=%d, GL = %d", 
            Rf_length(VectorOb), ii, rIntLoc );
        }
        REAL(VectorOb)[ii + (int) rIntLoc] = 
          (double) INTEGER(ChangePart)[ii];
      }     
    }
  } else if (Rf_isInteger(VectorOb)) {
    //Rprintf("Vector Ob INTEGER and Change Part is REAL\n"); R_FlushConsole();
    if (Rf_isReal(ChangePart)) {
      for (ii = 0; ii < Rf_length(ChangePart); ii++) {
        if (ii + rIntLoc > Rf_length(VectorOb)) {
          Rf_error("AlterVector, Rf_length of VectorOb only %d but ii=%d, GL = %d", 
            Rf_length(VectorOb), ii, rIntLoc );
        }
        INTEGER(VectorOb)[ii +(int)  rIntLoc ] = 
          (int) REAL(ChangePart)[ii];
      }
    }  else if (Rf_isInteger(ChangePart)) {
    //Rprintf("Vector Ob and ChangePart are INTEGER \n"); R_FlushConsole();
      for (ii = 0; ii < Rf_length(ChangePart); ii++) {
        if (ii + rIntLoc > Rf_length(VectorOb)) {
          Rf_error("AlterVector, Rf_length of VectorOb only %d but ii=%d, GL = %d", 
            Rf_length(VectorOb), ii, GetFirstInteger(IntLoc) );
        }
        INTEGER(VectorOb)[ii + GetFirstInteger(IntLoc)] = 
          INTEGER(ChangePart)[ii];
      }     
    }
  }
  return(VectorOb);
}

SEXP NEGLassoMakeTableShell(SEXP OutInvPsiValues, SEXP PsiValues,
  SEXP logInvPsiNearZeroValues, SEXP BetaValues,
  SEXP SOnLambda, SEXP SOnGammaSq, 
  SEXP SIntegral, SEXP SInvSum, SEXP SSecondInvSum, SEXP VerboseInt
  ) {
return(NEGLassoMakeTable(OutInvPsiValues, PsiValues,
  logInvPsiNearZeroValues, BetaValues,
  SOnLambda, SOnGammaSq, 
  SIntegral, SInvSum, SSecondInvSum, VerboseInt
  ));  
}

SEXP TwoLassoCrossValidateShell(SEXP rSVerbose, SEXP rSn, SEXP rSp, SEXP rSyys, SEXP rSxxs, SEXP rSXtY,
    SEXP rSXtX, SEXP rStt1, SEXP rStt2,
    SEXP rSLambdaA, SEXP rSLambdaD, 
    SEXP rSPiAVectorInputs, SEXP rSSigmaVectorInputs,
    SEXP rSLambdaAKInputs, SEXP rSLambdaDKInputs, 
    SEXP rSBetaStatStartMatrix,
    SEXP rSBetaCVVector,
    SEXP rSOrderSeq, 
    SEXP rSBBOn1, SEXP rSOnBeta, SEXP rRecBBOn1, SEXP rRecOnBeta,
    SEXP rSOnGammas, SEXP rSPiA, SEXP rSSigma,
    SEXP rSm, SEXP rSSigmaPrior,
    SEXP rSInitKKs, SEXP rSInverseGammaConstant, 
    SEXP rSCauchyEpsilon, SEXP rSMaxCauchy, 
    SEXP rSTDFNu, SEXP rSiiWeights, 
    SEXP rSL2ShrinkagePrior, SEXP rSL2ShrinkageRecords,
    SEXP rGroupSexp) {
  return(TwoLassoCrossValidate(rSVerbose, rSn, rSp, rSyys, rSxxs, rSXtY,
    rSXtX, rStt1, rStt2,
    rSLambdaA, rSLambdaD,
    rSPiAVectorInputs, rSSigmaVectorInputs,
    rSLambdaAKInputs, rSLambdaDKInputs, rSBetaStatStartMatrix,
    rSBetaCVVector,
    rSOrderSeq, 
    rSBBOn1, rSOnBeta, rRecBBOn1, rRecOnBeta,
    rSOnGammas, rSPiA, rSSigma,
    rSm, rSSigmaPrior,
    rSInitKKs, rSInverseGammaConstant, 
    rSCauchyEpsilon, rSMaxCauchy, 
    rSTDFNu, rSiiWeights, 
    rSL2ShrinkagePrior, rSL2ShrinkageRecords,
    rGroupSexp));    
    
}

SEXP WhatIsSuccessFlag(SEXP sTLS) {
  return(sWhatIsSuccessFlag(sTLS));
}

SEXP SetupConfidenceQuantiles(SEXP sTLS, SEXP iConfidenceQuantiles) {
  return(sSetupConfidenceQuantiles(sTLS, iConfidenceQuantiles));
}
SEXP getConfidenceQuantiles(SEXP sTLS) {
  return(sgetConfidenceQuantiles(sTLS));
}
SEXP getConfidenceMatrix(SEXP sTLS) {
  return(sgetConfidenceMatrix(sTLS));
}
SEXP getUnshrunkConfidenceMatrix(SEXP sTLS) {
  return(sgetUnshrunkConfidenceMatrix(sTLS));
}

SEXP getLambdaIndexForConfidenceIntervals(SEXP sTLS) {
  return(sgetLambdaIndexForConfidenceIntervals(sTLS));
}
SEXP setLambdaIndexForConfidenceIntervals(SEXP sTLS, SEXP sTwo) {
  return(ssetLambdaIndexForConfidenceIntervals(sTLS, sTwo));
}

SEXP getVerbose(SEXP sTLS) {
  return(sgetVerbose(sTLS));
}
SEXP setVerbose(SEXP sTLS, SEXP sTwo) {
  return(ssetVerbose(sTLS, sTwo));
}



/* ========================================================================== */
/*                                                                            */
/*   AcppExtern.c                                                             */
/*   (c) 2010 Alan Lenarcic                                                   */
/*                                                                            */
/*   Description                                                              */
/*                                                                            */
/* ========================================================================== */


#ifndef ACPPCONFIGUREME
  #include "AcppConfigureMe.h"
  #define ACPPCONFIGUREME 0
#endif


#define UpLo 'L'
#define TransWant 'N'
#ifndef RNeeds 
 #define RNeeds 0
 #include <R.h>
 #include <Rinternals.h>
 #include <Rdefines.h>
#endif


#ifndef RDYNLOADH
  #include <R_ext/Rdynload.h>
  #define RDYNLOADH 0
#endif


SEXP LockMeIn(SEXP aList, SEXP aOb) {
  SEXP gVI = R_NilValue;
  //Rprintf(" -- ACppExtern.c: Starting LockMeIn \n"); R_FlushConsole();
  //Rprintf(" -- aList is length %d\n", Rf_length(aList)); R_FlushConsole();
  if (Rf_isSymbol(aList)) {
    Rprintf("Acpp::AcppExtern.c:::LockMeIn:: Error aList is a symbol\n");
    if (Rf_isNull(aList)) {
       Rprintf("Acpp::AcppExtern.c:::LockMeIn:: Error aList is a NULL\n");
    }
    Rf_error("Acpp::AcppExtern.c:::LockMeIn:: This ain't good. \n");
  }
  if (Rf_length(aList) <= 0) {
    Rf_error("Acpp::AcppExtern:::LockMeIn:: not if a List is Length %d\n",
      Rf_length(aList));
  }
  if (Rf_isEnvironment(aList)) {
    Rf_error("Acpp::AcppExtern:::LockMeIn:: Not good because aList is an environment. \n");
  }
  for (int ii = 0; ii < Rf_length(aList); ii++) {
    gVI = VECTOR_ELT(aList, ii);
    if (Rf_isNull(gVI)) {
      //Rprintf(" -- ACppExtern.c::LockMeIn(): we now are setting into aList at postition %d \n", ii); R_FlushConsole();
      SET_VECTOR_ELT(aList, ii, aOb);
      SEXP iS = Rf_allocVector(INTSXP, 1);
      INTEGER(iS)[0] = ii;
      return(iS);
    }
  }
  SEXP iS = Rf_allocVector(INTSXP, 1);
  INTEGER(iS)[0] = -1;
  return(iS);
}
                  
SEXP GetMyPackageNamespace() {
  return(R_FindNamespace(Rf_mkString(PCKGNAME))); 
}
SEXP DoubleUpaList(SEXP aList) {
  SEXP newList = R_NilValue;
  Rf_protect(newList = Rf_allocVector(VECSXP, Rf_length(aList) *2));
  if (Rf_isNull(newList)) {
    Rf_error("DoubleUpaList: Error, we could not allocnewList");
  }
  int ii;
  for (ii = 0; ii < Rf_length(aList); ii++) {
    SET_VECTOR_ELT(newList, ii, VECTOR_ELT(aList, ii));
  }
  for (ii = Rf_length(aList); ii < Rf_length(newList); ii++) {
    SET_VECTOR_ELT(newList, ii, R_NilValue);
  }
  Rf_unprotect(1);
  return(newList);
}
SEXP KillFromaList(SEXP aList, SEXP skInt) {
  int kInt = -1;
  if (Rf_isNull(skInt) || Rf_length(skInt) <= 0) {
    Rf_error("Kill From a List, skInt is NULL length. \n");
  }
  if (Rf_isReal(skInt)) { kInt = (int) REAL(skInt)[0];
  } else if (Rf_isInteger(skInt)) { kInt = INTEGER(skInt)[0]; }
  if (kInt < 0) {
    Rf_error("Error:KillFromaList, skInt not supplied well. \n");
  }
  if (Rf_isNull(aList) || Rf_length(aList) <= 0) {
    Rf_error("Error:KillFromaList, aList is NULL length!\n ");
  }
  if (Rf_length(aList) <= kInt) {
    Rf_error("Error: KillFromaList, kInt=%d but length aList = %d\n",
      kInt, Rf_length(aList));
  }
  if (Rf_isSymbol(aList)) {
    Rprintf(" -- KillFromaList, this will fail because aList is a symbol.\n");
    R_FlushConsole();
    return(R_NilValue);
  }
  //Rprintf("--KillFromaList: Setting Now the aList[%d] to R_NilValue \n", kInt);
  R_FlushConsole();
  SET_VECTOR_ELT(aList, kInt, R_NilValue);
  return(skInt);
}

/*
R_CallMethodDef callMethods[] = {

  {NULL, NULL, 0}
};

void
R_init_Acpp(DllInfo *info)
{
R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
*/

R_CallMethodDef callMethods[] = {
//{"CoordinateDescentLassoShell2009", 
//  (DL_FUNC) &CoordinateDescentLassoShell2009,14},
{"RunG2GLMLasso", (DL_FUNC)&RunG2GLMLasso, 30},
//{"LarsConvergencyShell", (DL_FUNC)&LarsConvergencyShell, 39},
{"OutGetAllLambdaj", (DL_FUNC)&OutGetManyLambdaj, 8},
{"OutGetManyLambdaj", (DL_FUNC)&OutGetManyLambdaj, 8},
{"TwoLassoRegressionShell", 
  (DL_FUNC) &TwoLassoRegressionShell,30},
{"TwoLassoCreateTwoLassoSexpShell", 
  (DL_FUNC) &TwoLassoCreateTwoLassoSexpShell,31},
{"TwoLassoReRegressShell", 
  (DL_FUNC) &TwoLassoReRegressShell, 5},
{"NEGLassoMakeTableShell", 
  (DL_FUNC) &NEGLassoMakeTableShell,10},
{"AlterVector", (DL_FUNC) &AlterVector, 3},
{"TwoLassoCrossValidateShell", (DL_FUNC) &TwoLassoCrossValidateShell, 32},
{"WhatIsSuccessFlag", (DL_FUNC) &WhatIsSuccessFlag, 1},
{"SetupConfidenceQuantiles", (DL_FUNC) &SetupConfidenceQuantiles, 2},
{"getConfidenceQuantiles", (DL_FUNC) &getConfidenceQuantiles, 1},
{"getConfidenceMatrix", (DL_FUNC) &getConfidenceMatrix, 1},
{"getUnshrunkConfidenceMatrix", (DL_FUNC) &getUnshrunkConfidenceMatrix, 1},
{"getLambdaIndexForConfidenceIntervals", (DL_FUNC) &getLambdaIndexForConfidenceIntervals, 1},
{"setLambdaIndexForConfidenceIntervals", (DL_FUNC) &setLambdaIndexForConfidenceIntervals, 2},
{"getVerbose", (DL_FUNC) &getVerbose, 1},
{"setVerbose", (DL_FUNC) &setVerbose, 2},
{"LockMeIn", (DL_FUNC) &LockMeIn, 2},
{"DoubleUpaList", (DL_FUNC) &DoubleUpaList, 1},
{NULL, NULL,0}
};

void
R_init_TwoLasso(DllInfo *info)
{
R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

}

