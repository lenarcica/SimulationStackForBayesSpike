/* ========================================================================== */
/*                                                                            */
/*   CC_EMRIDGE2009.cc                                                        */
/*   (c) 2010 Alan Lenarcic                                                   */
/*                                                                            */
/*   An attempt at a "2 Sparse Ridge Estimator".   Bayesian estimators        */
/*   commonly compare two gaussians to compare inlusion.  This shows          */
/*   one way to try to implement an EM version of this comparison. Speed and  */
/*   memory issues tend to make this not a desireable approach.  Gibbs        */
/*   samplers may be more convenient with the challenges of Gaussian tails.   */
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

#ifndef CCEMRIDGE2009
   #include "CC_EMRIDGE2009.h"
   #define CCEMRIDGE2009 0
#endif
#ifndef MATHH
   #include "math.h"
#endif
#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif
#ifndef MyMatrixOp2009DD
  #include "MyMatrixOp2009.h"
  #define MyMatrixOp2009DD 0
#endif
#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <stdio.h>
  #define RMATH 0
#endif
#ifndef INCIDENTALOP2009DD
  #include "IncidentalOp2009.h"
  #define INCIDENTALOP2009DD 0 
#endif
#ifndef CoordinateDescentLassoDD
  #include "CoordinateDescentLasso.h"
  #define CoordinateDescentLassoDD 0
#endif

#ifndef IISNAN
   int isnan(double x) {return(R_isnancpp(x)); }
   #define IISNAN 0
#endif


int EMRIDGEStep (int NLen, int kLen, double *XTY, double *XTX, double sigmaSqNoise, 
    double tauDsq, double tauAsq, double TwotauDsq, double TwotauAsq,
    double sqrttauD2pi, double sqrttauA2pi, double puse,
    double *BHatVec, double *PutInMat, 
    double *PutInInvMat, double *OnBetas, double *OutBetas) {
	     if (isnan(TwotauDsq) || isinf(TwotauDsq)) {
		        Rprintf("TwotauDsq has gone nan\n");
		        R_FlushConsole();
		        return(-1);
	     }
	     if (isnan(TwotauAsq) || isinf(TwotauAsq)) {
		        Rprintf("TwotauASq has gone nan\n");
		        R_FlushConsole();
		        return(-1);
	     }	
	     if (isnan(sqrttauD2pi) || isinf(sqrttauD2pi)) {
		        Rprintf("sqrttau12pi has gone nan\n");
		        R_FlushConsole();
		        return(-1);
	     }
	     if (isnan(sqrttauA2pi) || isinf(sqrttauA2pi)) {
		        Rprintf("sqrttau12pi has gone nan\n");
		        R_FlushConsole();
		        return(-1);
	     }		     	
	double A1;
	double A2;
	double InvsqrttauA2pi = 1.0/sqrttauA2pi;
	double InvsqrttauD2pi = 1.0 / sqrttauD2pi;
	double InvTwotauASq = 1.0 / TwotauAsq;
	double InvTwotauDSq = 1.0 / TwotauDsq;
	double DNew;
	int ii,jj, cnt1 = 0;
	int InvertFlag;
	int FunFlag;
	for (ii = 0; ii < kLen; ii++) {
		for (jj =0; jj < kLen; jj++) {
			PutInMat[cnt1] = XTX[cnt1];
			cnt1++;
        }
    } 
    cnt1 = 0;
    for (ii = 0; ii < kLen; ii++) {
	     A1 = exp(-(double) OnBetas[ii] *(double)  OnBetas[ii] * InvTwotauDSq) * InvsqrttauD2pi;
	     A2 = exp(-(double) OnBetas[ii] * (double) OnBetas[ii] * InvTwotauASq) * InvsqrttauA2pi;
	     BHatVec[ii] = ((double)puse) * A2 / ( (double) puse * A2 + (1.0-(double) puse) * A1);
	     DNew = BHatVec[ii] / ((double) tauAsq) + (1-BHatVec[ii]) / (( double) tauDsq);
	     PutInMat[cnt1] += DNew * sigmaSqNoise;  
	     cnt1 = cnt1 + kLen + 1;       
    } 
     InvertFlag = InvertMatrix(kLen, PutInMat, PutInInvMat);
	 if (InvertFlag  < 0) {
		     Rprintf((char*)"EMRIDGEStep: Failure to Invert Matrix returning -1");
		     R_FlushConsole();
		     return(-1);
	 }
	 FunFlag = MatTimesVec(kLen, kLen, PutInInvMat, XTY, OutBetas);
	 return(FunFlag);  	   
}

int EMRIDGE(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, double SigmaSqNoise,
         double puse,
         double sigmaDsqStart, double sigmaAsqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut)   {
double *PutInInvMat; double *PutInMat; double *XTX;
double *OnBetas; double *OutBetas; double *BHatVec;
double *XTY;
OnBetas = NULL; OutBetas = NULL; BHatVec = NULL;
PutInMat = NULL; PutInInvMat = NULL; XTX = NULL;XTY = NULL;      	         
//Rprintf((char*)"Convergent EM: NLen = %d, kLen = %d, YData[0] = %f", NLen, kLen, YData[0]);
//Rprintf((char*)" puse = %f, ChangeNSteps = %d, ConvergenceNStpes = %d \n", 
//        puse, ChangeNSteps, ConvergenceNSteps);
//R_FlushConsole();
//Rprintf((char*)"SigmaSqNoise = %.4f, BStartBetas[0] = %f\n", SigmaSqNoise, BStartBetas[0]);
//R_FlushConsole();
XTY = (double *) Calloc(kLen+5, double);
if (XTY == NULL) {
	Rprintf((char*)"EMRIDGE: Not Enough Space for XTY \n");
	R_FlushConsole();	
	R_ProcessEvents();
}
XTX = (double *)Calloc((kLen+1) * (kLen+1), double);

SqMat (XTX, NLen, kLen, XData) ; 
tMatTimesVec(kLen, NLen, XData, YData, XTY);
PutInMat = (double *) Calloc((kLen+1)  * (kLen+1), double);
PutInInvMat = (double *)Calloc((kLen+1) * (kLen+1), double);

OnBetas = (double *)Calloc((kLen+1), double);
OutBetas = (double *)Calloc((kLen+1),double);     
BHatVec = (double *)Calloc( kLen+1, double);  
	         

double OnTauDsq = sigmaDsqStart;
double OnTauAsq = sigmaAsqStart;
double TwoOnTauDsq, TwoOnTauAsq, SqrtSigma12pi, SqrtSigma22pi;
int nnii, nnjj, cnt1;
int ii;
for (ii = 0; ii < kLen; ii++) { OnBetas[ii] = (double) BStartBetas[ii] ;
                                BHatVec[ii] = (double) BStartBBs[ii];}

  //Rprintf((char*)"What XTX is \n");
  int maxmax = kLen;
  if (maxmax > 6) { maxmax = 6; }
  nnii = 0;
  //for (ii = 0; ii < maxmax; ii++) {
  //	  Rprintf((char*)"%d [ ", (int) ii);
  //	  nnii = ii;
  //	  for (nnjj = 0; nnjj < maxmax; nnjj++) {
  //		  Rprintf((char*)"  %.4f,", (double) XTX[nnii]);
  //		  nnii += kLen;
  //    }
  //    Rprintf((char*)"]\n");
  //    R_FlushConsole();
  // }
  // Rprintf((char*)"\n, And XTY is: \n");
  // for (ii = 0; ii < maxmax; ii++) { 
  //	    Rprintf((char*)" %.4f,", (double) XTY[ii]);
  //	    R_FlushConsole();
  //  }
  //  Rprintf((char*)"\n");
  //  R_FlushConsole();
  
int ToyFlag = 0;
for (nnii = 0; nnii < ChangeNSteps; nnii++) {
	//Rprintf((char*)"EMConvergence:  On Values nnii = %d, OnTauDsq = %.4f, OnTauAsq = %.4f\n", 
	//           nnii, (double) OnTauDsq, (double) OnTauAsq);
	//R_FlushConsole();
	TwoOnTauDsq = 2.0 * OnTauDsq;
	TwoOnTauAsq = 2.0 * OnTauAsq;
	SqrtSigma12pi = sqrt( TwoOnTauDsq * M_PI);
	SqrtSigma22pi = sqrt( TwoOnTauAsq * M_PI);
	for (nnjj = 0; nnjj < ConvergenceNSteps; nnjj++) {
		ToyFlag =  EMRIDGEStep (NLen,  kLen, XTY, XTX, SigmaSqNoise, 
        OnTauDsq, OnTauAsq, TwoOnTauDsq, TwoOnTauAsq,
        SqrtSigma12pi, SqrtSigma22pi, puse,
        BHatVec, PutInMat, PutInInvMat, 
        OnBetas, OutBetas);
        if (ToyFlag >= 0) {
	         for (ii = 0; ii < kLen; ii++) {OnBetas[ii] = OutBetas[ii];}
        } else {
	        Rprintf((char*) "EMConvergence Error, Matrix not invertible\n");
	        R_FlushConsole();
	        return(-1);
        }
        //  Rprintf((char*)"What Solve(XTX+D) is \n");
        //  if (maxmax > 6) { maxmax = 6; }
        //  nnii = 0;
        // for (ii = 0; ii < maxmax; ii++) {
	    //    Rprintf((char*)"%d [ ", (int) ii);
	    //    for (ToyFlag = 0; ToyFlag < maxmax; ToyFlag++) {
	    // 	  Rprintf((char*)"  %.4f,", (double) PutInInvMat[ii  + kLen * ToyFlag]);
     	//    }
        //  Rprintf((char*)"]\n");
        //}
        //  R_FlushConsole();
        //Rprintf((char*)"\n,  Fit OnBetas is: \n");
        //  for (ii = 0; ii < maxmax; ii++) { 
	    //Rprintf((char*)" %.4f,", (double) OnBetas[ii]);
	    //R_FlushConsole();
        //}
        //Rprintf((char*)"\n");
        // R_FlushConsole();
        //Rprintf((char*)"\n,  Fit BHatVec is: \n");
        //  for (ii = 0; ii < maxmax; ii++) { 
	    //Rprintf((char*)" %.4f,", (double) BHatVec[ii]);
	    //R_FlushConsole();
        //}
        //Rprintf((char*)"\n");
        // R_FlushConsole();         
  
        cnt1 = (nnii * ConvergenceNSteps * kLen) + nnjj * kLen;
        for (ii = 0; ii < kLen; ii++) {
	        if (isnan(OnBetas[ii])) {
		           Rprintf("nnii = %d, nnjj = %d, OnBetas[%d] isnan\n",
		              nnii, nnjj, ii);
				   for (ii = 0; ii < kLen; ii++) {
                       Rprintf("BHatVec[%d] = %.4f, OnBetas[%d] = %.4f\n",
					      ii, (double) BHatVec[ii], ii, (double) OnBetas[ii]);					   
					   BBHatOut[ii ] = (double) BHatVec[ii];
					   BEndBetas[ii] = (double) OnBetas[ii];
					   
				   }
				   R_FlushConsole();		              
                   FFree(&BHatVec);
				   FFree(&OnBetas);
				   FFree(&PutInInvMat);
				   FFree(&PutInMat);
				   FFree(&OutBetas);
				   FFree(&XTY);
				   FFree(&XTX);
				  // Free(OnXTX);
				  // Free(OnXTY);		      		              
		        return(nnii);
	        }
	        if (isinf(OnBetas[ii])) {
		        Rprintf("nnii = %d, nnjj = %d, OnBetas[%d] isinf\n",
		              nnii, nnjj, ii);
				   for (ii = 0; ii < kLen; ii++) {
                       Rprintf("BHatVec[%d] = %.4f, OnBetas[%d] = %.4f\n",
					      ii, (double) BHatVec[ii], ii, (double) OnBetas[ii]);					   
					   BBHatOut[ii ] = (double) BHatVec[ii];
					   BEndBetas[ii] = (double) OnBetas[ii];
					   
				   }	
				   R_FlushConsole();	              
                   FFree(&BHatVec);
				   FFree(&OnBetas);
				   FFree(&PutInInvMat);
				   FFree(&PutInMat);
				   FFree(&OutBetas);
				   FFree(&XTY);
				   FFree(&XTX);
				  // Free(OnXTX);
				  // Free(OnXTY);		      		              
		        return(nnii);
	        }	        
	        BetaPartBetas[cnt1 + ii] = (double) OnBetas[ii];
	        BBHatRec[cnt1 + ii] = (double) BHatVec[ii];
	    }
    }
    //cnt1 = (nnii *ConvergenceNSteps) * kLen;
    //for (ii = 0; ii < kLen; ii++) {
	//     BetaPartBetas[cnt1 + ii] = (double) OnBetas[ii];
	//     BBHatRec[cnt1 + ii] = (double) BHatVec[ii];
	//}
    OnTauDsq = OnTauDsq *sigma1Mult;
    if (nnii < sigma2MultStop) { OnTauAsq= OnTauAsq * sigma2Mult; }		
}
   for (ii = 0; ii < kLen; ii++) {
	   BBHatOut[ii ] = (double) BHatVec[ii];
	   BEndBetas[ii] = (double) OnBetas[ii];
	   
   }
   FFree(&BHatVec);
   FFree(&OnBetas);
   FFree(&PutInInvMat);
   FFree(&PutInMat);
   FFree(&OutBetas);
   FFree(&XTY);
   FFree(&XTX);
    return(1);

}
///////////////////////////////////////////////////////////
//  EMRIDGEStep
//    Function that calculates the new expected B values
//    Given estimates for Betas.
//
int EMRIDGEStep (int NLen, int kLen, double *XTY, double *XTX, double sigmaSqNoise, 
    double tauDsq, double tauAsq, double TwotauDsq, double TwotauAsq,
    double sqrttauD2pi, double sqrttauA2pi, double puse,
    double *BHatVec, double *PutInMat, 
    double *PutInInvMat, double *OnBetas, double *OutBetas, 
    double logp, double log1mp, double SigmaNoiseInv) {
	     if (isnan(TwotauDsq) || isinf(TwotauDsq)) {
		        Rprintf("TwotauDSq has gone nan\n");
		        R_FlushConsole();
		        return(-1);
	     }
	     if (isnan(TwotauAsq) || isinf(TwotauAsq)) {
		        Rprintf("TwotauASq has gone nan\n");
		        R_FlushConsole();
		        return(-1);
	     }	
	     if (isnan(sqrttauD2pi) || isinf(sqrttauD2pi)) {
		        Rprintf("sqrttau12pi has gone nan\n");
		        R_FlushConsole();
		        return(-1);
	     }
	     if (isnan(sqrttauA2pi) || isinf(sqrttauA2pi)) {
		        Rprintf("sqrttau12pi has gone nan\n");
		        R_FlushConsole();
		        return(-1);
	     }		     	
	double A1;
	double A2;
	double InvsqrttauA2pi = 1.0/sqrttauA2pi;
	double InvsqrttauD2pi = 1.0 / sqrttauD2pi;
	double InvTwotauASq = 1.0 / TwotauAsq;
	double InvTwotauDSq = 1.0 / TwotauDsq;	    
	double DNew;
	int ii,jj, cnt1 = 0;
	int InvertFlag;
	int FunFlag;
	for (ii = 0; ii < kLen; ii++) {
		for (jj =0; jj < kLen; jj++) {
			PutInMat[cnt1] = XTX[cnt1];
			cnt1++;
        }
    } 
    double SigmaNoiseInvFlip = 1 / SigmaNoiseInv;
    cnt1 = 0;
    for (ii = 0; ii < kLen; ii++) {
	     A1 = exp(-(double) OnBetas[ii] *(double)  OnBetas[ii] * (InvTwotauDSq)) * InvsqrttauD2pi;
	     A2 = exp(-(double) OnBetas[ii] * (double) OnBetas[ii] * (InvTwotauASq)) * InvsqrttauA2pi;
	     BHatVec[ii] = ((double) puse * A2) / ( (double) puse * A2 + (1.0-(double) puse) * A1);
	     DNew = BHatVec[ii] / ((double) tauAsq) + (1-BHatVec[ii]) / (( double) tauDsq);
	     PutInMat[cnt1] += DNew * SigmaNoiseInvFlip;  
	     cnt1 = cnt1 + kLen + 1;       
    } 
     InvertFlag = InvertMatrix(kLen, PutInMat, PutInInvMat);
	 if (InvertFlag  < 0) {
		     Rprintf((char*)"EMRIDGEStep: Failure to Invert Matrix returning -1");
		     R_FlushConsole();
		     return(-1);
	 }
	 FunFlag = MatTimesVec(kLen, kLen, PutInInvMat, XTY, OutBetas);
	 return(FunFlag);  	   
}

int EMDYNAMICRIDGE(int NLen, int kLen, double *YData, double *XData, 
         int TotalRuns, int NumEMConv, double SigmaSqNoise,
         double puse,
         double StTauSqD, double StTauSqA, 
         double TauDMult, double TauAMult, int TauAMultStop, 
         double *TauSqASeq, double *TauSqDSeq,
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv,
         int NumCDOConv, double CDOEpsilon, int *OrderSeq,
         int *CountDescentOperations) {
	  int PrintFlag = 0;
	  int SFlag;
    if (kLen < TOOMUCHKAPPA || NumCDOConv < 0) {
	    if (PrintFlag > 0 ) {		        
		    Rprintf("kLen = %d, TMK = %d, Num CDOConv = %d \n", kLen, TOOMUCHKAPPA, NumCDOConv);
		    R_FlushConsole(); R_ProcessEvents();
	    }
	    if (BetaPartBetas == NULL || BetaPartBetas[0] < 0 ||
	        BBHatRec == NULL || BBHatRec[0] < 0 ) {
         if (PrintFlag > 0) {
	         if ( BetaPartBetas == NULL) { 
	             Rprintf(" BetaPartBetas is NULL \n"); R_FlushConsole(); R_ProcessEvents();
             } else  {
	             Rprintf(" BetaPartBetas[0] = %.4f \n", BetaPartBetas[0]); 
	             R_FlushConsole(); R_ProcessEvents();
             }
	         if ( BBHatRec == NULL) { 
	             Rprintf(" BBHatRec is NULL \n"); R_FlushConsole(); R_ProcessEvents();
             } else  {
	             Rprintf(" BBHatRec[0] = %.4f \n", BBHatRec[0]); 
	             R_FlushConsole(); R_ProcessEvents();
             }  
             Rprintf("  Obviously we're doing EMRIDGENoRec \n");           
        }
     // int EMRIDGENoRec(int NLen, int kLen, double *YData, double *XData, 
	 //        int ChangeNSteps, int ConvergenceNSteps, 
	 //        double ConvergeCloseEnough, double SigmaSqNoise,
	 //        double puse,
	 //        double sigma1sqStart, double sigma2sqStart, 
	 //        double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
	 //        double *BStartBBs, double *BStartBetas,
	 //       double *BEndBetas, double *BBHatOut,
	 //        double logp, double log1mp, double SigmaNoiseInv) ;		        
	       SFlag = EMRIDGENoRec( NLen, kLen, YData, XData, TotalRuns,
		       NumEMConv, CDOEpsilon, SigmaSqNoise, puse, StTauSqD, 
		       StTauSqA, TauDMult, TauAMult, TauAMultStop,
		       BStartBBs, BStartBetas, BEndBetas, 
		       BBHatOut, logp, log1mp, SigmaNoiseInv);		        	        
	    } else {
         if (PrintFlag > 0) {
	         if ( BetaPartBetas == NULL) { 
	             Rprintf(" BetaPartBetas is NULL \n"); R_FlushConsole(); R_ProcessEvents();
             } else  {
	             Rprintf(" BetaPartBetas[0] = %.4f \n", BetaPartBetas[0]); 
	             R_FlushConsole(); R_ProcessEvents();
	             Rprintf(" BetaPartBetas[ TotalRuns * kLen * NumEMConv-1= %d] = %.4f \n",
	                 TotalRuns * kLen * NumEMConv -1,
	                 BetaPartBetas[TotalRuns * kLen * NumEMConv -1]);
	             R_FlushConsole(); R_ProcessEvents();
             }
	         if ( BBHatRec == NULL) { 
	             Rprintf(" BBHatRec is NULL \n"); R_FlushConsole(); R_ProcessEvents();
             } else  {
	             Rprintf(" BBHatRec[0] = %.4f \n", BBHatRec[0]); 
	             R_FlushConsole(); R_ProcessEvents();
	             Rprintf(" BBHatRec[ TotalRuns * kLen * NumEMConv-1= %d] = %.4f \n",
	                 TotalRuns * kLen * NumEMConv -1,
	                 BBHatRec[TotalRuns * kLen * NumEMConv -1]);
	             R_FlushConsole(); R_ProcessEvents();	             
             }  
             Rprintf("  Obviously we're doing EMRIDGENoRec \n");           
        }		    	    
	       SFlag = EMRIDGE( NLen, kLen, YData, XData, TotalRuns,
	         NumEMConv, SigmaSqNoise, puse, StTauSqD, 
	         StTauSqA, TauDMult, TauAMult, TauAMultStop,
	         BStartBBs, BStartBetas, BetaPartBetas, BEndBetas, BBHatRec, 
	         BBHatOut, logp, log1mp, SigmaNoiseInv);				    
	    }
    } else {
	    //int EMCDORIDGE(int NLen, int kLen, double *YData, double * XData, int TotalRuns,
	    //   int NumEMSteps, double SigmaSqNoise, double piAuse, double StTauSqA, 
	    //   double StTauSqD, double TauAMult, double TauDMult, int TauAMultStop,
	    //   double *TauSqASeq, double *TauSqDSeq,
	    //   double *BStartBBs, double *BStartBetas, double *BetaPartBetas, 
	    //   double *BEndBetas, double *BBHatRec, 
	    //   double *BBHatOut, int NumCDOConv, double CDOEpsilon, int *OrderSeq);  
	    SFlag = EMCDORIDGE(NLen, kLen, YData, XData, TotalRuns,
	       NumEMConv, SigmaSqNoise, puse, StTauSqA, 
	       StTauSqD, TauAMult, TauDMult, TauAMultStop,
	       TauSqASeq, TauSqDSeq,
	       BStartBBs, BStartBetas, BetaPartBetas, BEndBetas, BBHatRec, 
	       BBHatOut, NumCDOConv, CDOEpsilon, OrderSeq,
	       CountDescentOperations);
    }
    //Rprintf("EMDynamicRIDGE, All Done \n"); R_FlushConsole(); R_ProcessEvents();
    return(SFlag);	         
}	         
/////////////////////////////////////////////////////////////////////////////
//  EMRIDGE
//      EMRIDGE where the part Betas are recorded in BetaPartBetas
//       And the part "B"'s are recorded in BBHatOut
int EMRIDGE(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, double SigmaSqNoise,
         double puse,
         double sigmaDsqStart, double sigmaAsqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv)   {
double *PutInInvMat; double *PutInMat; double *XTX;
double *OnBetas; double *OutBetas; double *BHatVec;
double *XTY;
OnBetas = NULL; OutBetas = NULL; BHatVec = NULL;
PutInMat = NULL; PutInInvMat = NULL; XTX = NULL;XTY = NULL;         
//Rprintf((char*)"Convergent EM: NLen = %d, kLen = %d, YData[0] = %f", NLen, kLen, YData[0]);
//Rprintf((char*)" puse = %f, ChangeNSteps = %d, ConvergenceNStpes = %d \n", 
//        puse, ChangeNSteps, ConvergenceNSteps);
//R_FlushConsole();
//Rprintf((char*)"SigmaSqNoise = %.4f, BStartBetas[0] = %f\n", SigmaSqNoise, BStartBetas[0]);
//R_FlushConsole();
XTY = (double *) Calloc(kLen+5, double);
if (XTY == NULL) {
	Rprintf((char*)"EMRIDGE: Not Enough Space for XTY \n");
	R_FlushConsole();	
	R_ProcessEvents();
}
XTX = (double *)Calloc((kLen+1) * (kLen+1), double);

SqMat (XTX, NLen, kLen, XData) ; 
tMatTimesVec(kLen, NLen, XData, YData, XTY);
PutInMat = (double *) Calloc((kLen+1)  * (kLen+1), double);
PutInInvMat = (double *)Calloc((kLen+1) * (kLen+1), double);





OnBetas = (double *)Calloc((kLen+1), double);
OutBetas = (double *)Calloc((kLen+1),double);     
BHatVec = (double *)Calloc( kLen+1, double);  
	         

double OnTauDsq = sigmaDsqStart;
double OnTauAsq = sigmaAsqStart;
double TwoOnTauDsq, TwoOnTauAsq, SqrtSigma12pi, SqrtSigma22pi;
int nnii, nnjj, cnt1;
int ii;
for (ii = 0; ii < kLen; ii++) { OnBetas[ii] = (double) BStartBetas[ii] ;
                                BHatVec[ii] = (double) BStartBBs[ii];}

  //Rprintf((char*)"What XTX is \n");
  int maxmax = kLen;
  if (maxmax > 6) { maxmax = 6; }
  nnii = 0;
  //for (ii = 0; ii < maxmax; ii++) {
  //	  Rprintf((char*)"%d [ ", (int) ii);
  //	  nnii = ii;
  //	  for (nnjj = 0; nnjj < maxmax; nnjj++) {
  //		  Rprintf((char*)"  %.4f,", (double) XTX[nnii]);
  //		  nnii += kLen;
  //    }
  //    Rprintf((char*)"]\n");
  //    R_FlushConsole();
  // }
  // Rprintf((char*)"\n, And XTY is: \n");
  // for (ii = 0; ii < maxmax; ii++) { 
  //	    Rprintf((char*)" %.4f,", (double) XTY[ii]);
  //	    R_FlushConsole();
  //  }
  //  Rprintf((char*)"\n");
  //  R_FlushConsole();
  
int ToyFlag = 0;
for (nnii = 0; nnii < ChangeNSteps; nnii++) {
	//Rprintf((char*)"EMConvergence:  On Values nnii = %d, OnTauDsq = %.4f, OnTauAsq = %.4f\n", 
	//           nnii, (double) OnTauDsq, (double) OnTauAsq);
	//R_FlushConsole();
	TwoOnTauDsq = 2.0 * OnTauDsq;
	TwoOnTauAsq = 2.0 * OnTauAsq;
	SqrtSigma12pi = sqrt( TwoOnTauDsq * M_PI);
	SqrtSigma22pi = sqrt( TwoOnTauAsq * M_PI);
	for (nnjj = 0; nnjj < ConvergenceNSteps; nnjj++) {
		ToyFlag =  EMRIDGEStep(NLen,  kLen, XTY, XTX, SigmaSqNoise, 
        OnTauDsq, OnTauAsq, TwoOnTauDsq, TwoOnTauAsq,
        SqrtSigma12pi, SqrtSigma22pi, puse,
        BHatVec, PutInMat, PutInInvMat, 
        OnBetas, OutBetas, logp, log1mp, SigmaNoiseInv);
        if (ToyFlag >= 0) {
	         for (ii = 0; ii < kLen; ii++) {OnBetas[ii] = OutBetas[ii];}
        } else {
	        Rprintf((char*) "EMConvergence Error, Matrix not invertible\n");
	        R_FlushConsole();
	        return(-1);
        }
        //  Rprintf((char*)"What Solve(XTX+D) is \n");
        //  if (maxmax > 6) { maxmax = 6; }
        //  nnii = 0;
        // for (ii = 0; ii < maxmax; ii++) {
	    //    Rprintf((char*)"%d [ ", (int) ii);
	    //    for (ToyFlag = 0; ToyFlag < maxmax; ToyFlag++) {
	    // 	  Rprintf((char*)"  %.4f,", (double) PutInInvMat[ii  + kLen * ToyFlag]);
     	//    }
        //  Rprintf((char*)"]\n");
        //}
        //  R_FlushConsole();
        //Rprintf((char*)"\n,  Fit OnBetas is: \n");
        //  for (ii = 0; ii < maxmax; ii++) { 
	    //Rprintf((char*)" %.4f,", (double) OnBetas[ii]);
	    //R_FlushConsole();
        //}
        //Rprintf((char*)"\n");
        // R_FlushConsole();
        //Rprintf((char*)"\n,  Fit BHatVec is: \n");
        //  for (ii = 0; ii < maxmax; ii++) { 
	    //Rprintf((char*)" %.4f,", (double) BHatVec[ii]);
	    //R_FlushConsole();
        //}
        //Rprintf((char*)"\n");
        // R_FlushConsole();         
  
        cnt1 = (nnii * ConvergenceNSteps * kLen) + nnjj * kLen;
        for (ii = 0; ii < kLen; ii++) {
            if (isnan(OnBetas[ii])) {
		        Rprintf("nnii = %d, nnjj = %d, OnBetas[%d] isnan\n",
		              nnii, nnjj, ii);
				   for (ii = 0; ii < kLen; ii++) {
                     Rprintf("BHatVec[%d] = %.4f, OnBetas[%d] = %.4f\n",
					      ii, (double) BHatVec[ii], ii, (double) OnBetas[ii]);					   
					   BBHatOut[ii ] = (double) BHatVec[ii];
					   BEndBetas[ii] = (double) OnBetas[ii];
					   
				   }		              
                   FFree(&BHatVec);
				   FFree(&OnBetas);
				   FFree(&PutInInvMat);
				   FFree(&PutInMat);
				   FFree(&OutBetas);
				   FFree(&XTY);
				   FFree(&XTX);
				  // Free(OnXTX);
				  // Free(OnXTY);		              
		        return(nnii);
	        }
	        if (isinf(OnBetas[ii])) {
		        Rprintf("nnii = %d, nnjj = %d, OnBetas[%d] isinf\n",
		              nnii, nnjj, ii);
		        Rprintf("Well Print you some stuff now\n");
				   for (ii = 0; ii < kLen; ii++) {
					   Rprintf("BHatVec[%d] = %.4f, OnBetas[%d] = %.4f\n",
					      ii, (double) BHatVec[ii], ii, (double) OnBetas[ii]);
					   BBHatOut[ii ] = (double) BHatVec[ii];
					   BEndBetas[ii] = (double) OnBetas[ii];
					   
				   }	
				   R_FlushConsole();	              
                   FFree(&BHatVec);
				   FFree(&OnBetas);
				   FFree(&PutInInvMat);
				   FFree(&PutInMat);
				   FFree(&OutBetas);
				   FFree(&XTY);
				   FFree(&XTX);
				  // Free(OnXTX);
				  // Free(OnXTY);		      		              
		        return(nnii);
	        }	   
	        BetaPartBetas[cnt1 + ii] = (double) OnBetas[ii];
	        BBHatRec[cnt1 + ii] = (double) BHatVec[ii];
	    }
    }
    //cnt1 = (nnii *ConvergenceNSteps) * kLen;
    //for (ii = 0; ii < kLen; ii++) {
	//     BetaPartBetas[cnt1 + ii] = (double) OnBetas[ii];
	//     BBHatRec[cnt1 + ii] = (double) BHatVec[ii];
	//}
    OnTauDsq = OnTauDsq *sigma1Mult;
    if (nnii < sigma2MultStop) { OnTauAsq= OnTauAsq * sigma2Mult; }		
}
   for (ii = 0; ii < kLen; ii++) {
	   BBHatOut[ii ] = (double) BHatVec[ii];
	   BEndBetas[ii] = (double) OnBetas[ii];
	   
   }
   FFree(&BHatVec);
   FFree(&OnBetas);
   FFree(&PutInInvMat);
   FFree(&PutInMat);
   FFree(&OutBetas);
   FFree(&XTY);
   FFree(&XTX);
    return(1);

}


////////////////////////////////////////////////////////////////////////
//  EMRIDGENoRec:  Run EMRIDGE algorithm, do not record
//
//
int EMRIDGENoRec(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, 
         double ConvergeCloseEnough,
         double SigmaSqNoise,
         double puse,
         double sigmaDsqStart, double sigmaAsqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BEndBetas, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv)   {
	         int TPFlag = 0;
	if (TPFlag > 0) {  
	   Rprintf((char*)"Convergent EM: NLen = %d, kLen = %d, YData[0] = %f", NLen, kLen, YData[0]);
	   Rprintf((char*)" puse = %f, ChangeNSteps = %d, ConvergenceNStpes = %d \n", 
	          puse, ChangeNSteps, ConvergenceNSteps);
	   R_FlushConsole();
	   Rprintf((char*)"SigmaSqNoise = %.4f, BStartBetas[0] = %f\n", SigmaSqNoise, BStartBetas[0]);
	   R_FlushConsole();
    }
double * XTY = NULL;
XTY = (double *) Calloc(kLen+5, double);
if (XTY == NULL) {
	Rprintf((char*)"EMRIDGE: Not Enough Space for XTY \n");
	R_FlushConsole();	
	R_ProcessEvents();
}
double *XTX = NULL;
XTX = (double *)Calloc((kLen+1) * (kLen+1), double);

SqMat (XTX, NLen, kLen, XData) ; 
tMatTimesVec(kLen, NLen, XData, YData, XTY);
double *PutInMat = NULL; double * PutInInvMat = NULL;
PutInMat = (double *) Calloc((kLen+1)  * (kLen+1), double);
PutInInvMat = (double *)Calloc((kLen+1) * (kLen+1), double);

double *OnBetas = NULL; double *OutBetas = NULL; 
double *OlderPBetas = NULL; double *BHatVec = NULL;
OnBetas = (double *)Calloc((kLen+1), double);
OutBetas = (double *)Calloc((kLen+1),double); 
OlderPBetas = (double *)Calloc((kLen+1), double);    
BHatVec = (double *)Calloc( kLen+1, double);  
	         
const int MaxPleasePrint = 150;
double OnTauDsq = sigmaDsqStart;
double OnTauAsq = sigmaAsqStart;
double TwoOnTauDsq, TwoOnTauAsq, SqrtSigma12pi, SqrtSigma22pi;
int nnii, nnjj, cnt1 = 0;
if (cnt1 < 0) {
}
int ii;
for (ii = 0; ii < kLen; ii++) { OnBetas[ii] = (double) BStartBetas[ii] ;
                                BHatVec[ii] = (double) BStartBBs[ii];}
for (ii = 0; ii < kLen; ii++) { OlderPBetas[ii] = (double) OnBetas[ii]; }
  int maxmax = kLen;
  if (maxmax > 6) { maxmax = 6; }   
  if (TPFlag > 3) {
    Rprintf((char*)"What XTX is \n"); R_FlushConsole();
      PrintRMatrix(XTX, kLen, kLen); Rprintf((char*) "\n"); R_FlushConsole();
  }
int ToyFlag = 0;
double TotalTT;
for (nnii = 0; nnii < ChangeNSteps; nnii++) {
	if (MaxPleasePrint < kLen && 
	    (double) round(nnii / MaxPleasePrint) == (double) nnii / (double)MaxPleasePrint)
	{
		  Rprintf((char*)"EMConvergence, on nnii = %d, kLen = %d\n", nnii, kLen);
		  R_FlushConsole();
		  R_ProcessEvents();
     }
    if (TPFlag > 1) {
	   Rprintf((char*)"EMConvergence:  On Values nnii = %d, OnTauDsq = %.4f, OnTauAsq = %.4f\n", 
	             nnii, (double) OnTauDsq, (double) OnTauAsq);
	   R_FlushConsole();
    }
	TwoOnTauDsq = 2.0 * OnTauDsq;
	TwoOnTauAsq = 2.0 * OnTauAsq;
	SqrtSigma12pi = sqrt( TwoOnTauDsq * M_PI);
	SqrtSigma22pi = sqrt( TwoOnTauAsq * M_PI);
	for (nnjj = 0; nnjj < ConvergenceNSteps; nnjj++) {
		ToyFlag =  EMRIDGEStep(NLen,  kLen, XTY, XTX, SigmaSqNoise, 
	        OnTauDsq, OnTauAsq, TwoOnTauDsq, TwoOnTauAsq,
	        SqrtSigma12pi, SqrtSigma22pi, puse,
	        BHatVec, PutInMat, PutInInvMat, 
	        OnBetas, OutBetas, logp, log1mp, SigmaNoiseInv);
	    //if (ISNAN(OnBetas[0])) {
		//    Rprintf((char*)"EMConvergence: We Have ISNAN OnBetas[0] = %.4f \n", (double) OnBetas[0]);
		//    Rprintf((char*)"nnii = %d, nnjj = %d, sigmaDsq = %.6f, sigmaAsq = %.6f\n",
		//       nnii, nnjj, (double) OnTauDsq, (double) OnTauAsq);
		//    return(-1);
	    // } 
        if (ToyFlag >= 0) {
	         TotalTT = 0;
	         for (ii = 0; ii < kLen; ii++) {
		         TotalTT = fabs((double) OutBetas[ii] - (double) OnBetas[ii]);
		         OnBetas[ii] = OutBetas[ii];
		     }
        } else {
	        Rprintf((char*) "EMConvergence Error, Matrix not invertible\n");
	        R_FlushConsole();
	        return(-1);
        }
        if (TotalTT < ConvergeCloseEnough) {
	        nnjj = ConvergenceNSteps;
	        break;    
        }

	    if (TPFlag > 4) {  
	          Rprintf((char*)"What Solve(XTX+D) is \n");
	          if (maxmax > 6) { maxmax = 6; }
	          nnii = 0;
	         for (ii = 0; ii < maxmax; ii++) {
		        Rprintf((char*)"%d [ ", (int) ii);
		        for (ToyFlag = 0; ToyFlag < maxmax; ToyFlag++) {
		     	  Rprintf((char*)"  %.4f,", (double) PutInInvMat[ii  + kLen * ToyFlag]);
	     	    }
	          Rprintf((char*)"]\n");
	        }
	          R_FlushConsole();
        }
        if (TPFlag > 2) {
		        Rprintf((char*)"\n,  Fit OnBetas is: \n");
		          for (ii = 0; ii < maxmax; ii++) { 
			    Rprintf((char*)" %.4f,", (double) OnBetas[ii]);
			    R_FlushConsole();
		        }
		        Rprintf((char*)"\n");
		         R_FlushConsole();
		        Rprintf((char*)"\n,  Fit BHatVec is: \n");
		          for (ii = 0; ii < maxmax; ii++) { 
			    Rprintf((char*)" %.4f,", (double) BHatVec[ii]);
			    R_FlushConsole();
		        }
		        Rprintf((char*)"\n");
		         R_FlushConsole();         
        }
        cnt1 = (nnii * ConvergenceNSteps * kLen) + nnjj * kLen;
    }

    OnTauDsq = OnTauDsq *sigma1Mult;
             TotalTT = 0;
             for (ii = 0; ii < kLen; ii++) {
		         TotalTT = fabs((double) OlderPBetas[ii] - (double) OnBetas[ii]);
		         OlderPBetas[ii] = OnBetas[ii];
		     }
    if (nnii < sigma2MultStop) { OnTauAsq= OnTauAsq * sigma2Mult; }		
}
   for (ii = 0; ii < kLen; ii++) {
	   BBHatOut[ii ] = (double) BHatVec[ii];
	   BEndBetas[ii] = (double) OnBetas[ii];
	   
   }
   if (TPFlag > 1) {
	   Rprintf((char*) "EMRIDGENoRec, end of Algorithm, Freeing \n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing BHatVec\n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   FFree(&BHatVec);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OnBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OnBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing PutInInvMat\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&PutInInvMat);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing PutInMat\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&PutInMat);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OutBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OutBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OlderPBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OlderPBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing XTY\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&XTY);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing XTX\n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   FFree(&XTX);
   if (TPFlag > 1) {
	   Rprintf((char*) "All Done I Survive\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
    return(1);

}


/////////////////////////////////////////////////////////////////
//  EMRIDGENoRecKiller
//  
//     Run the EMRidge algorithm do not record intermediate steps;
int EMRIDGENoRecKiller(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, 
         double ConvergeCloseEnough,
         double SigmaSqNoise,
         double puse,
         double sigmaDsqStart, double sigmaAsqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BEndBetas, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv)   {
int ii;	  
int TPFlag = 0; 
int ToyFlag = 0;
double TotalTT;
//Rprintf((char*)"Convergent EM: NLen = %d, kLen = %d, YData[0] = %f", NLen, kLen, YData[0]);
//Rprintf((char*)" puse = %f, ChangeNSteps = %d, ConvergenceNStpes = %d \n", 
//        puse, ChangeNSteps, ConvergenceNSteps);
//R_FlushConsole();
//Rprintf((char*)"SigmaSqNoise = %.4f, BStartBetas[0] = %f\n", SigmaSqNoise, BStartBetas[0]);
//R_FlushConsole();
double *XTY = NULL;
XTY = (double *) Calloc(kLen+5, double);
if (XTY == NULL) {
	Rprintf((char*)"EMRIDGE: Not Enough Space for XTY \n");
	R_FlushConsole();	
	R_ProcessEvents();
}
double *XTX = NULL;
XTX = (double *)Calloc((kLen+1) * (kLen+1), double);

SqMat (XTX, NLen, kLen, XData) ; 
tMatTimesVec(kLen, NLen, XData, YData, XTY);
double *PutInMat = NULL; double *PutInInvMat = NULL; double *OnBetas=NULL;
double * OutBetas = NULL; double *OlderPBetas = NULL;
double *BHatVec = NULL;
PutInMat = (double *) Calloc((kLen+1)  * (kLen+1), double);
PutInInvMat = (double *)Calloc((kLen+1) * (kLen+1), double);

OnBetas = (double *)Calloc((kLen+1), double);
OutBetas = (double *)Calloc((kLen+1),double); 
OlderPBetas = (double *)Calloc((kLen+1), double);    
BHatVec = (double *)Calloc( kLen+1, double);  
int *ActiveSet = NULL;
ActiveSet = (int *) Calloc(kLen+1, int);
int sKlen;
sKlen = kLen;
for (ii = 0; ii < kLen; ii++) { ActiveSet[ii] = ii; }

double *OnXTY = NULL; double *OnXTX = NULL;
OnXTY = (double *)Calloc( kLen +1, double );
OnXTX = (double *)Calloc( (kLen+1) * (kLen+1), double);         
MakeOnXTX(kLen, ActiveSet,  sKlen, XTX, OnXTX);
MakeOnXTY(kLen, ActiveSet,  sKlen, XTX, OnXTY);

double OnTauDsq = sigmaDsqStart;
double OnTauAsq = sigmaAsqStart;
double TwoOnTauDsq, TwoOnTauAsq, SqrtSigma12pi, SqrtSigma22pi;
int nnii, nnjj, cnt1 = 0;
for (ii = 0; ii < kLen; ii++) { OnBetas[ii] = (double) BStartBetas[ii] ;
                                BHatVec[ii] = (double) BStartBBs[ii];}
for (ii = 0; ii < kLen; ii++) { OlderPBetas[ii] = (double) OnBetas[ii]; }
   
  if (cnt1 == 0) {
  }
  //Rprintf((char*)"What XTX is \n");
  int maxmax = kLen;
  if (maxmax > 6) { maxmax = 6; }
  nnii = 0; nnjj = 0;
  //for (ii = 0; ii < maxmax; ii++) {
  //	  Rprintf((char*)"%d [ ", (int) ii);
  //	  nnii = ii;
  //	  for (nnjj = 0; nnjj < maxmax; nnjj++) {
  //		  Rprintf((char*)"  %.4f,", (double) XTX[nnii]);
  //		  nnii += kLen;
  //    }
  //    Rprintf((char*)"]\n");
  //    R_FlushConsole();
  // }
  // Rprintf((char*)"\n, And XTY is: \n");
  // for (ii = 0; ii < maxmax; ii++) { 
  //	    Rprintf((char*)" %.4f,", (double) XTY[ii]);
  //	    R_FlushConsole();
  //  }
  //  Rprintf((char*)"\n");
  //  R_FlushConsole();
  
for (nnii = 0; nnii < ChangeNSteps; nnii++) {

	//Rprintf((char*)"EMConvergence:  On Values nnii = %d, OnTauDsq = %.4f, OnTauAsq = %.4f\n", 
	//           nnii, (double) OnTauDsq, (double) OnTauAsq);
	//R_FlushConsole();
	TwoOnTauDsq = 2.0 * OnTauDsq;
	TwoOnTauAsq = 2.0 * OnTauAsq;
	SqrtSigma12pi = sqrt( TwoOnTauDsq * M_PI);
	SqrtSigma22pi = sqrt( TwoOnTauAsq * M_PI);
	for (nnjj = 0; nnjj < ConvergenceNSteps; nnjj++) {
		ToyFlag =  EMRIDGEStep (NLen,  kLen, XTY, XTX, SigmaSqNoise, 
        OnTauDsq, OnTauAsq, TwoOnTauDsq, TwoOnTauAsq,
        SqrtSigma12pi, SqrtSigma22pi, puse,
        BHatVec, PutInMat, PutInInvMat, 
        OnBetas, OutBetas, logp, log1mp, SigmaNoiseInv);
        if (ToyFlag >= 0) {
	         TotalTT = 0;
	         for (ii = 0; ii < kLen; ii++) {
		         TotalTT = fabs((double) OutBetas[ii] - (double) OnBetas[ii]);
		         OnBetas[ii] = OutBetas[ii];
		     }
        } else {
	        Rprintf((char*) "EMConvergence Error, Matrix not invertible\n");
	        R_FlushConsole();
	        return(-1);
        }
        if (TotalTT < ConvergeCloseEnough) {
	        nnjj = ConvergenceNSteps;
	        break;    
        }

	        
        //  Rprintf((char*)"What Solve(XTX+D) is \n");
        //  if (maxmax > 6) { maxmax = 6; }
        //  nnii = 0;
        // for (ii = 0; ii < maxmax; ii++) {
	    //    Rprintf((char*)"%d [ ", (int) ii);
	    //    for (ToyFlag = 0; ToyFlag < maxmax; ToyFlag++) {
	    // 	  Rprintf((char*)"  %.4f,", (double) PutInInvMat[ii  + kLen * ToyFlag]);
     	//    }
        //  Rprintf((char*)"]\n");
        //}
        //  R_FlushConsole();
        //Rprintf((char*)"\n,  Fit OnBetas is: \n");
        //  for (ii = 0; ii < maxmax; ii++) { 
	    //Rprintf((char*)" %.4f,", (double) OnBetas[ii]);
	    //R_FlushConsole();
        //}
        //Rprintf((char*)"\n");
        // R_FlushConsole();
        //Rprintf((char*)"\n,  Fit BHatVec is: \n");
        //  for (ii = 0; ii < maxmax; ii++) { 
	    //Rprintf((char*)" %.4f,", (double) BHatVec[ii]);
	    //R_FlushConsole();
        //}
        //Rprintf((char*)"\n");
        // R_FlushConsole();         
  
        cnt1 = (nnii * ConvergenceNSteps * kLen) + nnjj * kLen;
    }
    //cnt1 = (nnii *ConvergenceNSteps) * kLen;
    //for (ii = 0; ii < kLen; ii++) {
	//     BetaPartBetas[cnt1 + ii] = (double) OnBetas[ii];
	//     BBHatRec[cnt1 + ii] = (double) BHatVec[ii];
	//}
    OnTauDsq = OnTauDsq *sigma1Mult;
             TotalTT = 0;
             for (ii = 0; ii < kLen; ii++) {
		         TotalTT = fabs((double) OlderPBetas[ii] - (double) OnBetas[ii]);
		         OlderPBetas[ii] = OnBetas[ii];
		     }
    if (nnii < sigma2MultStop) { OnTauAsq= OnTauAsq * sigma2Mult; }		
}
   for (ii = 0; ii < kLen; ii++) {
	        if (isnan(OnBetas[ii])) {
		        Rprintf("nnii = %d, nnjj = %d, OnBetas[%d] isnan\n",
		              nnii, nnjj, ii);
				   for (ii = 0; ii < kLen; ii++) {
					   BBHatOut[ii ] = (double) BHatVec[ii];
					   BEndBetas[ii] = (double) OnBetas[ii];
					   
				   }		              
                   FFree(&BHatVec);
				   FFree(&OnBetas);
				   FFree(&PutInInvMat);
				   FFree(&PutInMat);
				   FFree(&OutBetas);
				   FFree(&XTY);
				   FFree(&XTX);
				   FFree(&OnXTX);
				   FFree(&OnXTY);		      
		        return(nnii);
	        }
	        if (isinf(OnBetas[ii])) {
		        Rprintf("nnii = %d, nnjj = %d, OnBetas[%d] isinf\n",
		              nnii, nnjj, ii);
				   for (ii = 0; ii < kLen; ii++) {
					   BBHatOut[ii ] = (double) BHatVec[ii];
					   BEndBetas[ii] = (double) OnBetas[ii];
					   
				   }		              
if (TPFlag > 1) {
	   Rprintf((char*) "EMRIDGENoRec, end of Algorithm, Freeing \n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing BHatVec\n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   FFree(&BHatVec);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OnBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OnBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing PutInInvMat\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&PutInInvMat);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing PutInMat\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&PutInMat);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OutBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OutBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OlderPBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OlderPBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing XTY\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&XTY);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing XTX\n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   FFree(&XTX);	      	              
		        return(nnii);
	        }	
	   BBHatOut[ii ] = (double) BHatVec[ii];
	   BEndBetas[ii] = (double) OnBetas[ii];	   
   }
if (TPFlag > 1) {
	   Rprintf((char*) "EMRIDGENoRec, end of Algorithm, Freeing \n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing BHatVec\n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   FFree(&BHatVec);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OnBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OnBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing PutInInvMat\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&PutInInvMat);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing PutInMat\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&PutInMat);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OutBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OutBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing OlderPBetas\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&OlderPBetas);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing XTY\n");
	   R_FlushConsole(); R_ProcessEvents();
   }   
   FFree(&XTY);
   if (TPFlag > 1) {
	   Rprintf((char*) "Freeing XTX\n");
	   R_FlushConsole(); R_ProcessEvents();
   }
   FFree(&XTX);   
    return(1);

}

/////////////////////////////////////////////////////////////////////
//  MakeOnXTX
//
//    If one wishes to use a deliberate subset of XTX matrix, 
//       this subset constructs it.
//
int MakeOnXTX(int kLen, int *ActiveLen, int sKlen, double *XTX, double *OnXTX) {
   int cn2 = 0, cn3 = 0, cn4 = 0;
   int ii,jj;
   for (ii = 0; ii < sKlen; ii++) {
	   cn2 = kLen * ActiveLen[ii];
	   cn3 = ii + sKlen  *ii;
	   cn4 = sKlen * ii + ii;
	   for (jj = ii; jj < sKlen; jj++) {
		   	OnXTX[cn3] = XTX[cn2 + ActiveLen[jj]];
		   	OnXTX[cn4] = XTX[cn2 + ActiveLen[jj]];
		    cn3+=sKlen;
		    cn4++;
      }
    }
 	return(1);	   	
}

////////////////////////////////////////////////////////////////////////
//  MakeOnXTY
//
//    If one wishes to use a deliberate subset of XTY matrix, this subset constructs it.
//
int MakeOnXTY(int kLen, int *ActiveLen, int sKlen, double *XTY, double *OnXTY) {
    int ii;
    for (ii = 0; ii < sKlen; ii++) {
	    OnXTY[ii] = XTY[ActiveLen[ii]];	
    }
    return(1);
}


int EMCDORIDGE(int NLen, int kLen, double *YData, double * XData, int TotalRuns,
	       int NumEMSteps, double SigmaSqNoise, double piAuse, double StTauSqA, 
	       double StTauSqD, double TauAMult, double TauDMult, int TauAMultStop,
	       double *TauSqASeq, double *TauSqDSeq,
	       double *BStartBBs, double *BStartBetas, double *BetaPartBetas, 
	       double *BEndBetas, double *BBHatRec, 
	       double *BBHatOut, int NumCDOConv, double CDOEpsilon, int *OrderSeq,
	       int *CountDescentOperations) {
int SFlag = 0;       	
EMRIDGEObject *EMO = new EMRIDGEObject(NLen, kLen, XData, YData, 
               SigmaSqNoise,  piAuse, 
               TotalRuns, NumEMSteps, 
               StTauSqA, StTauSqD, TauAMult, TauDMult, 
               TauAMultStop, BStartBBs, BStartBetas,
               TauSqASeq, TauSqDSeq,
               NumCDOConv, CDOEpsilon,
               OrderSeq);	
               
if (EMO == NULL) {
	Rprintf("EMCDORIDGE, Error: EMRIDGEObject EMO did not load: NULL \n");
	R_FlushConsole(); R_ProcessEvents();  return(-1);
}	
if (EMO->SuccessFlag < 0) {
	Rprintf("EMCDORIDGE, Error: EMRIDGEObject EMO->SuccessFlag =%d\n", 
	             EMO->SuccessFlag);
	R_FlushConsole(); R_ProcessEvents();  return(-1);
}	
EMO->PrintFlag = 0;
if (BBHatRec == NULL || BetaPartBetas == NULL || BBHatRec[0] < 0 ||
      BetaPartBetas[0] < 0 ) {
	    if (EMO->PrintFlag >= 1) { 
	      if  (BBHatRec == NULL || BetaPartBetas == NULL) {
		       Rprintf("No Rec: BBHat Rec is null or BetaPartBetas \n");
		       R_FlushConsole(); R_ProcessEvents();
           } else {
               Rprintf("No Rec: BBHatRec[0] == %.4f and BetaPartBetas[0] = %.4f \n",
                         BBHatRec[0], BetaPartBetas[0]);
		       R_FlushConsole(); R_ProcessEvents();	       
           }
        }
} else {
	    if (EMO->PrintFlag >= 1) { 
		       Rprintf("Creating Rec: BBHatRec[0] == %.4f and BetaPartBetas[0] = %.4f \n",
		             BBHatRec[0], BetaPartBetas[0]);
		       R_FlushConsole(); R_ProcessEvents();
        }	
      SFlag = EMO->SetupRecs();
}
if (SFlag < 0) {
	Rprintf("SetupRecs Did not Operate Error \n\n\n\n"); R_FlushConsole(); R_ProcessEvents();
	delete(EMO); return(-1);
}
if (CountDescentOperations != NULL && CountDescentOperations[0] >= 0) {
	 if (EMO->PrintFlag >= 1) { 
		       Rprintf("Creating CountDescOp: CountDescentOperations[0] = %.4f \n",
		             CountDescentOperations[0]);
		       R_FlushConsole(); R_ProcessEvents();
		       Rprintf("Creating CountDescOp: CountDescentOperations[TotalRuns-1=%d] = %.4f \n",
		             TotalRuns-1, CountDescentOperations[TotalRuns-1]);
		       R_FlushConsole(); R_ProcessEvents();
        }	
	SFlag = EMO->SetupCountDescentOperations();
} else {
	 Rprintf("NOT Creating CountDescOp: CountDescentOperations[0] = %.4f \n",
		             CountDescentOperations[0]);	
}
if (SFlag < 0) {
	Rprintf("SetupCountDescentOperations Did not Operate Error \n\n\n\n"); 
	  R_FlushConsole(); R_ProcessEvents();
	delete(EMO); return(-1);
}
SFlag = EMO->TwoRidgeSequence();
if (EMO->SuccessFlag < 0 || SFlag < 0) {
	Rprintf("EMCDORIDGE, ERROR ERROR ERROR SFlag = %d, and EMO->SuccessFlag = %d\n",
	    SFlag, EMO->SuccessFlag); R_FlushConsole(); R_ProcessEvents();
}	
if (EMO->PrintFlag > 0) {
	Rprintf("EMCDORIDGE: Finished The TwoRidgeSequence Moving to copy Betas\n",
	    SFlag, EMO->SuccessFlag); R_FlushConsole(); R_ProcessEvents();
}       


int ii;
  for (ii = 0; ii < kLen; ii++) {
	  BEndBetas[ii] = EMO->OnBetas[ii];
	  BBHatOut[ii] = EMO->OnBBs[ii];
  } 

int MaxBB;
if (EMO->PrintFlag > 0) {
	Rprintf("EMCDORIDGE: Finished The TwoRidgeSequence Moving to copy BBHatRecs\n",
	    SFlag, EMO->SuccessFlag); R_FlushConsole(); R_ProcessEvents();
}
 
if (BBHatRec == NULL || BetaPartBetas == NULL || BBHatRec[0] < 0 ||
      BetaPartBetas[0] < 0 ) {
} else {
	  MaxBB = kLen * TotalRuns;
	  for (ii = 0; ii < MaxBB; ii++) {
		    BBHatRec[ii] = EMO->RecBBs[ii];
		    BetaPartBetas[ii] = EMO->RecBetas[ii];		    
	  }   	
}

if (EMO->PrintFlag > 0) {
	Rprintf("EMCDORIDGE: Finished About to Record the CountDescentOperations\n"); 
	 R_FlushConsole(); R_ProcessEvents();
} 
if (CountDescentOperations != NULL && CountDescentOperations[0] >= 0) {
	for (ii = 0; ii < TotalRuns; ii++) {
		 CountDescentOperations[ii] = EMO->CountDescentOperations[ii];
   }
}

delete(EMO);
    //Rprintf("EMCDORidge, All Done \n"); R_FlushConsole(); R_ProcessEvents();
    return(SFlag);
}







int EMRIDGEObject::SetupPISeq(double rpiASt, double rm1, double rm2, double rSigmaSqEta,
                              double rSigmaSqBar) {
	    if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject::SetupPISeq()     Starting \n");		
		      R_FlushConsole();
		      R_ProcessEvents();    
	    } 	                              
	    if (rm1 < 0 || rm2 < 0) {
		    piSeq = NULL;
		    piAuse = rpiASt;
	    } else {
		    m1 = rm1; m2 = rm2;
		    piAuse = rpiASt;
			   if (piSeq != NULL) {
				   Rprintf("EMRIDGEObject::SetupPISeq, piSeq already not null!\n");
				   R_FlushConsole(); R_ProcessEvents();
			   } else {
			       piSeq = (double *) Calloc(TotalRuns +2, double);
			       if (piSeq == NULL) {
				   Rprintf("EMRIDGEObject::SetupPISeq, No room for piSeq\n");
				   R_FlushConsole(); R_ProcessEvents();	SuccessFlag = 1; return(-1);       
			       }
			    }
         }
	    if (rSigmaSqEta < 0 || rSigmaSqBar < 0) {
		    SigmaSqEta = -1;
		    SigmaSqBar = -1;
		    SigmaSqSeq = NULL;
	    } else {
		    SigmaSqEta = rSigmaSqEta; SigmaSqBar = rSigmaSqBar;
		    SigmaSqNoise = SigmaSqBar;
			   if (SigmaSqSeq != NULL) {
				   Rprintf("EMRIDGEObject::SetupPISeq, SigmaSqSeq already not null!\n");
				   R_FlushConsole(); R_ProcessEvents();
			   } else {
			       SigmaSqSeq = (double *) Calloc(TotalRuns +2, double);
			       if (piSeq == NULL) {
				   Rprintf("EMRIDGEObject::SetupPISeq, No room for SigmaSqSeq\n");
				   R_FlushConsole(); R_ProcessEvents();	SuccessFlag = 1; return(-1);       
			       }
			    }
         }         
         return(1);
}
int EMRIDGEObject::SetupOrderSeq(int *rOrderSeq) {
	    if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject::SetupOrderSeq()     Starting \n");		
		      R_FlushConsole();
		      R_ProcessEvents();    
	    } 	
   if (OrderSeq != NULL) {
	   Rprintf("EMRIDGEObject::SetupOrderSeq, OrderSeq already not null!\n");
	   R_FlushConsole(); R_ProcessEvents();
   } else {
       OrderSeq = (int *) Calloc(TotalRuns +2, int);
       if (OrderSeq == NULL) {
	   Rprintf("EMRIDGEObject::SetupOrderSeq, No room for OrderSeq\n");
	   R_FlushConsole(); R_ProcessEvents();	SuccessFlag = 1; return(-1);       
       }
    }
   int ii; 
   if (rOrderSeq == NULL) {
         if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject::SetupOrderSeq:  rOrderSeq was NULL\n");
		    R_FlushConsole(); R_ProcessEvents();	
	     }   
	   for (ii = 0; ii < TotalRuns; ii++) {
		   OrderSeq[ii] = NumEMSteps;
       }
   } else {
         if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject::SetupOrderSeq:  Copying in rOrderSeq\n");
		    R_FlushConsole(); R_ProcessEvents();		    
		      Rprintf("      That rOrderSeq is : \n");
		      PrintVector( rOrderSeq, TotalRuns );
		    R_FlushConsole(); R_ProcessEvents();	
	     } 	   
	   for (ii = 0; ii < TotalRuns; ii++) {
		   OrderSeq[ii] = rOrderSeq[ii];
       }
   }	
	    if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject::SetupOrderSeq:  Giving OrderSeq: \n");
		    R_FlushConsole(); R_ProcessEvents();
		    PrintVector(OrderSeq, TotalRuns);
		    R_FlushConsole(); R_ProcessEvents();
		    Rprintf("EMRIDGEObject:: SetupOrderSeq: OrderSeq[%d] =%d \n", 
		                        TotalRuns-1, OrderSeq[TotalRuns-1]);

		    R_FlushConsole(); R_ProcessEvents();
	    
	    }    
   return(1);
}
double EMRIDGEObject::EMCDONewNuj() {
		if (PrintFlag > 6) {
		   Rprintf((char*) "EMRIDGEObject: EMCDONewNuj(), OnCoord = %d\n", CDO->OnCoord); 
		   R_FlushConsole(); R_ProcessEvents();
	   	}	
	  int OnCoord = CDO->OnCoord;
	  //int OnCoordOnCoord = CDO->OnCoord * CDO->kLen + CDO->OnCoord;
	  double Nuj = (CDO->XTResid[OnCoord] - CDO->OnGammas[OnCoord] * CDO->OnBeta[OnCoord]) / 
	         (*(CDO->pXTX[OnCoord]+OnCoord) +  CDO->OnGammas[OnCoord] );
	  return(Nuj);
};

int EMRIDGEObject::EMCDOneStep() {
	  double NewNuj = EMCDONewNuj();
		if (PrintFlag > 6) {
		   Rprintf((char*) "EMRIDGEObject: EMCDOneStep(), CDO->OnCoord = %d\n", CDO->OnCoord); 
		   R_FlushConsole(); R_ProcessEvents();
	   	}

	  if (  (double) NewNuj == (double) 0.0 ||
	        fabs(NewNuj) <- CDOEpsilon / ((double)kLen) )  {
		CDO->RecentMove = 0;
		return(0);
      }
     CDO->RecentMove = NewNuj;
    
	 int ii;
	 for (ii = 0; ii < kLen; ii++) {
		 CDO->XTResid[ii] = ((double) CDO->XTResid[ii]) - 
		 ((double) *(CDO->pXTX[CDO->OnCoord]+ii) * (double) NewNuj);
	 }
	 CDO->OnBeta[CDO->OnCoord] = CDO->OnBeta[CDO->OnCoord] + NewNuj;
   return(1);
}
int EMRIDGEObject::CoordinateDescentRidge() {
	if (PrintFlag > 4) {
	    Rprintf("EMRIDGEObject:CoordinateDescentRidge()      About to Start\n");
	    R_FlushConsole(); R_ProcessEvents();
	}		
	double TotalEpForLoop = 0;
	int SFlag = 0;
    for (CDiiLoop = 0; CDiiLoop < NumCDOConv; CDiiLoop++) {
	    if (CountDescentOperations != NULL) {
		     (CountDescentOperations[Seqtt])++;
	     }
		if (PrintFlag > 5) {
		    Rprintf("EMRIDGEObject:CoordinateDescentRidge()      CDiiLoop = %d\n", 
					                CDiiLoop);
		    R_FlushConsole(); R_ProcessEvents();
		}
		 TotalEpForLoop = 0;	
		 for (CDOnjj = 0; CDOnjj < kLen; CDOnjj++) {
			if (PrintFlag > 6) {
			    Rprintf("EMRIDGEObject:CoordinateDescentRidge()           CDOnjj = %d\n", 
						                CDOnjj);
			    R_FlushConsole(); R_ProcessEvents();
		    }	
			 CDO->OnCoord = CDOnjj;
			 CDO->RecentMove = 0;
			 SFlag = EMCDOneStep();
			 if (SFlag < 0) {
				 Rprintf("EMRIDGEObject::CoordinateDescentRidge  Step Failed\n");
				 R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
	         }
			 TotalEpForLoop += fabs(CDO->RecentMove);
         }
         if (TotalEpForLoop < CDOEpsilon) {
	         return(1);
         }
    }
    return(1);
}
int EMRIDGEObject::InitiateOnBBs() {
			if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:InitiateOnBBs()      Starting\n");
			    R_FlushConsole(); R_ProcessEvents();
		    }	
		    if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:InitiateOnBBs()      piAuse is %.5e\n", piAuse);
			    R_FlushConsole(); R_ProcessEvents();
		    }	
	   double InvTauSqA = 1.0 / (TauSqASeq[0]);    
       double InvTauSqD = 1.0 / (TauSqDSeq[0]);    
       int ii;       
       for (ii = 0; ii < kLen; ii++) {
	       if (OnBetas[ii] > 0 && OnBBs[ii] <= 0) {
		       OnBBs[ii] = piAuse * exp( - .5 * (InvTauSqA-InvTauSqD) * OnBetas[ii] * OnBetas[ii] ) * 
		                sqrt( InvTauSqA * .5);
		       OnBBs[ii] = OnBBs[ii] / ( OnBBs[ii] + 
		            (1.0-piAuse) * 1.0 * 
		                sqrt( InvTauSqD * .5));
	       } else if (OnBBs[ii] <= 0) {
		       OnBBs[ii] = piAuse;
	       }
		       OnInvTauSq[ii] = OnBBs[ii] * InvTauSqA + ( 1.0 - OnBBs[ii]) * InvTauSqD; 
		             
       }
       CDO->SetUpGammas(OnInvTauSq);	
			if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:InitiateOnBBs()      Ending\n");
			    R_FlushConsole(); R_ProcessEvents();
		    }       
       return(1);	
}
int EMRIDGEObject::CalcOnBBs() {
		    if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:CalcOnBBs()        Starting\n");
			    R_FlushConsole(); R_ProcessEvents();
		    }	
	        if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:CalcOnBBs,  Before OnBBs is \n");
			       PrintVector(OnBBs, kLen);
			    Rprintf("EMRIDGEObject:CalcOnBBs()        Ending\n");
			    R_FlushConsole(); R_ProcessEvents();
		    }			    	
	   double InvTauSqA = 1.0 / (TauSqASeq[Seqtt]);    
       double InvTauSqD = 1.0 / (TauSqDSeq[Seqtt]);  
       //double sqrtInvTauSqAp5 = sqrt(InvTauSqA * .5);
       //double sqrtInvTauSqDp5 = sqrt(InvTauSqD * .5);
       double sqrtInvTauSqAdivD = sqrt( TauSqDSeq[Seqtt] / TauSqASeq[Seqtt] );
       int ii;  
       	     if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:CalcOnBBs InvTauSqA = %.5e, InvTauSqD = %.5e \n",
			        InvTauSqA, InvTauSqD );
			    Rprintf("EMRIDGEObject:  And sqrtInvTauSqAdivD = %.5e \n", sqrtInvTauSqAdivD);
			    Rprintf("EMRIDGEObject: and piAuse = %.5e \n", piAuse);
			    Rprintf("EMRIDGEObject: exp(-.5 * (InvTauSqD - InvTauSqA) * 1 ) = %.5e \n", 
			                exp(- .5 * (InvTauSqD - InvTauSqA) ) );
			    R_FlushConsole(); R_ProcessEvents();
		    }       
       for (ii = 0; ii < kLen; ii++) {
		       OnBBs[ii] = piAuse *
		                sqrtInvTauSqAdivD;
		       OnBBs[ii] = OnBBs[ii] / ( OnBBs[ii] + 
		            (1.0-piAuse) * exp( - .5 * (InvTauSqD - InvTauSqA) * OnBetas[ii] * OnBetas[ii] ));
		       OnInvTauSq[ii] = OnBBs[ii] * InvTauSqA + ( 1.0 - OnBBs[ii]) * InvTauSqD; 
       }	
            CDO->SetUpGammas(OnInvTauSq);	
		    if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:CalcOnBBs,  Immediately after OnBBs is \n");
			       PrintVector(OnBBs, kLen);
			    Rprintf("EMRIDGEObject:CalcOnBBs()        Ending\n");
			    R_FlushConsole(); R_ProcessEvents();
		    }	         
       return(1);
}
int EMRIDGEObject::TwoRidgeSequence() {
		if (PrintFlag > 0) {
		    Rprintf("EMRIDGEObject:TwoRidgeSequence() Starting \n");
		    R_FlushConsole(); R_ProcessEvents();
	    } 
    int SFlag = 0;
    int jj;
    
    int MaxEMStepOn = 0;
    InitiateOnBBs();
     if (PrintFlag > 2) {
	      Rprintf("EMRIDGEObject:After InitiateOnBBs, the OnBBs is \n");
	      PrintVector(OnBBs, kLen); R_FlushConsole();
     }
     int RecttOn,  iiRec;
    for (Seqtt = 0; Seqtt < TotalRuns; Seqtt++) {	    
	    if (OrderSeq[Seqtt] < 0) {
		    SFlag = CalcOnBBs();
		    MaxEMStepOn = -OrderSeq[Seqtt];
		    if (SFlag < 0) {
			    Rprintf("EMRIDGEObject::TwoRidgeSequence, CalcOnBBs() returned -1\n");
			    Rprintf("    Seqtt = %d/%d, At Begginning because OrderSeq[%d] = %d\n", 
			                Seqtt,TotalRuns, Seqtt, OrderSeq[Seqtt]);
			    R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		    }			    
	    } else {
		    MaxEMStepOn = OrderSeq[Seqtt];
	    }
	    if (CountDescentOperations != NULL) { 
	          (CountDescentOperations[Seqtt]) = 0;
        }
	    for (EMStept2 = 0; EMStept2 <  MaxEMStepOn; EMStept2++) {
		    if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:TwoRidgeSequence()    Starting  EMStept2 = %d/%d, Seqtt = %d/%d\n", 
				                EMStept2, MaxEMStepOn, Seqtt,TotalRuns);
			    R_FlushConsole(); R_ProcessEvents();
		    } 
		    SFlag = CDO->MakeXTResid();
		    if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:TwoRidgeSequence()      About to CDRidge\n", 
				                EMStept2, MaxEMStepOn, Seqtt,TotalRuns);
			    R_FlushConsole(); R_ProcessEvents();
		    }	
    
		    SFlag = CoordinateDescentRidge();
		    if (SFlag < 0) {
			    Rprintf("EMRIDGEObject::TwoRidgeSequence, CDORidge returned -1\n");
			    Rprintf("    Seqtt = %d/%d, EMStept2 = %d/%d\n", 
			                Seqtt,TotalRuns, EMStept2, MaxEMStepOn);
			    R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		    }
		    if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:TwoRidgeSequence()      CDRidge Finished Copying Beta\n", 
				                EMStept2, MaxEMStepOn, Seqtt,TotalRuns);
			    R_FlushConsole(); R_ProcessEvents();
		    }			    
		    for (jj = 0; jj < kLen; jj++) {
			    OnBetas[jj] = CDO->OnBeta[jj];
		    }
		    SFlag = CalcOnBBs();
		    if (SFlag < 0) {
			    Rprintf("EMRIDGEObject::TwoRidgeSequence, CalcOnBBs() returned -1\n");
			    Rprintf("    Seqtt = %d/%d, EMStept2 = %d/%d\n", 
			                Seqtt,TotalRuns, EMStept2, MaxEMStepOn);
			    R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		    }
		    if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:TwoRidgeSequence() Current OnBetas:\n", 
				                EMStept2, MaxEMStepOn, Seqtt,TotalRuns);
				           PrintVector(OnBetas, kLen);
			    R_FlushConsole(); R_ProcessEvents();
			    Rprintf("EMRIDGEObject:TwoRidgeSequence() Current OnBBs:\n", 
				                EMStept2, MaxEMStepOn, Seqtt,TotalRuns);
				           PrintVector(OnBBs, kLen);
			    R_FlushConsole(); R_ProcessEvents();			    
		    } 			    	
			if (PrintFlag > 3) {
			    Rprintf("EMRIDGEObject:TwoRidgeSequence() Finished  EMStept2 = %d/%d, Seqtt = %d/%d\n", 
				                EMStept2, MaxEMStepOn, Seqtt,TotalRuns);
			    R_FlushConsole(); R_ProcessEvents();
		    } 		
	    	    
	    }
	   	if (RecBBs != NULL) {
		   	   RecttOn = Seqtt * kLen;
		   	   for (iiRec = 0; iiRec < kLen; iiRec++) {
			   	   RecBBs[RecttOn] = OnBBs[iiRec]; RecttOn++;
		       }
	     } 
  	     if (RecBetas != NULL) {
		   	   RecttOn = Seqtt * kLen;
		   	   for (iiRec = 0; iiRec < kLen; iiRec++) {
			   	   RecBetas[RecttOn] = OnBetas[iiRec]; RecttOn++;
		       }
	     } 
		if (PrintFlag > 2) {
		    Rprintf("EMRIDGEObject:TwoRidgeSequence() Finished  Seqtt = %d/%d, EMStept2 = %d/%d\n", 
			                Seqtt,TotalRuns, EMStept2, MaxEMStepOn);
		    R_FlushConsole(); R_ProcessEvents();
	    } 	     	     
    }
		if (PrintFlag > 0) {
		    Rprintf("EMRIDGEObject:TwoRidgeSequence() Finished \n");
		    R_FlushConsole(); R_ProcessEvents();
	    }     	
	return(1);
}

