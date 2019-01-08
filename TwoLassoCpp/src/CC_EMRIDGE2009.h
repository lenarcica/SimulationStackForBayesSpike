#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <Rdefines.h>
  #include <Rinternals.h>
  #include <stdio.h>
  #define RMATH 0
#endif

#undef isinf
#define isinf(x) (((x) > 999999999999.9) || ((x) < -999999999999.9))
#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif


#ifndef CoordinateDescentLassoDD
  #include "CoordinateDescentLasso2014.h"
  #define CoordinateDescentLassoDD 0
#endif

#define TOOMUCHKAPPA 50

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
         int NumCDOConv, double CDOEpsilon, int *OrderSeq, int *CountDescentOperations);
         
int EMRIDGE(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, double SigmaSqNoise,
         double puse,
         double sigma1sqStart, double sigma2sqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut);
int EMRIDGE (int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, double SigmaSqNoise,
         double puse,
         double sigma1sqStart, double sigma2sqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut, double logp, double log1mp, 
         double SigmaNoiseInv, int *OrderSeq);    
int EMCDORIDGE(int NLen, int kLen, double *YData, double * XData, int TotalRuns,
	       int NumEMSteps, double SigmaSqNoise, double piAuse, double StTauSqA, 
	       double StTauSqD, double TauAMult, double TauDMult, int TauAMultStop,
	       double *TauSqASeq, double *TauSqDSeq,
	       double *BStartBBs, double *BStartBetas, double *BetaPartBetas, 
	       double *BEndBetas, double *BBHatRec, 
	       double *BBHatOut, int NumCDOConv, double CDOEpsilon, int *OrderSeq,
	       int *CountDescentOperations);           

int EMRIDGENoRec(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, 
         double ConvergeCloseEnough, double SigmaSqNoise,
         double puse,
         double sigma1sqStart, double sigma2sqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BEndBetas, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv) ;       
         
         
int EMRIDGE (int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, double SigmaSqNoise,
         double puse,
         double sigma1sqStart, double sigma2sqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut);
int EMRIDGEStep (int NLen, int kLen, double *XTY, double *XTX, double sigmaSqNoise, 
    double tauDsq, double tauAsq, double TwotauDsq, double TwotauAsq,
    double sqrttauD2pi, double sqrttauA2pi, double puse,
    double *BHatVec, double *PutInMat, 
    double *PutInInvMat, double *OnBetas, double *OutBetas);
int EMRIDGEStep (int NLen, int kLen, double *XTY, double *XTX, double sigmaSqNoise, 
    double tauDsq, double tauAsq, double TwotauDsq, double TwotauAsq,
    double sqrttauD2pi, double sqrttauA2pi, double puse,
    double *BHatVec, double *PutInMat, 
    double *PutInInvMat, double *OnBetas, double *OutBetas, 
    double logp, double log1mp, double SigmaNoiseInv);
int EMRIDGE(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, double SigmaSqNoise,
         double puse,
         double sigma1sqStart, double sigma2sqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv);
int EMRIDGENoRec(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, 
         double ConvergeCloseEnough,
         double SigmaSqNoise,
         double puse,
         double sigma1sqStart, double sigma2sqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BEndBetas, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv);
int EMRIDGENoRecKiller(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, 
         double ConvergeCloseEnough,
         double SigmaSqNoise,
         double puse,
         double sigma1sqStart, double sigma2sqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BEndBetas, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv);            
int MakeOnXTX(int kLen, int *ActiveLen, int sKlen, double *XTX, double *OnXTX);
int MakeOnXTY(int kLen, int *ActiveLen, int sKlen, double *XTY, double *OnXTY);
 
int EMDYNAMICRIDGE(int NLen, int kLen, double *YData, double *XData, 
         int ChangeNSteps, int ConvergenceNSteps, double SigmaSqNoise,
         double puse,
         double sigma1sqStart, double sigma2sqStart, 
         double sigma1Mult, double sigma2Mult, int sigma2MultStop, 
         double *BStartBBs, double *BStartBetas,
         double *BetaPartBetas,
         double *BEndBetas, double *BBHatRec, double *BBHatOut,
         double logp, double log1mp, double SigmaNoiseInv);        
         

class EMRIDGEObject {
	public:
	
int PrintFlag;
int SuccessFlag;

int TotalRuns;
int NLen; int kLen;

CoordinateDescentObject *CDO;	

double *XTX; double *XTY; double *xxs; double *yys;

double *piSeq;
double *SigmaSqSeq;


double *TauSqASeq; double *TauSqDSeq;


double *OnBBs;
double *OnBetas;
double *RecBBs; double *RecBetas;
int *CountDescentOperations;
double *OnInvTauSq;
int *OrderSeq;

int TauAMultStop;

double OnTauSqA; double OnTauSqD;
double StTauSqA; double StTauSqD;
double TauAMult; double TauDMult;
double SigmaSqNoise;

double piAuse;
double m1, m2;
double SigmaSqEta; double SigmaSqBar;

int Seqtt; int EMStept2; int CDiiLoop; int CDOnjj;


int EMCDOneStep();
int EMCDOUpdate();
double EMCDONewNuj();
int TauFlag;

int InitiateOnBBs();
int CalcOnBBs();
int CoordinateDescentRidge();
int TwoRidgeSequence();

int NumCDOConv; double CDOEpsilon;
int NumEMSteps;


int SetupOrderSeq(int *rOrderSeq);

int SetupPISeq(double rpiASt, double rm1, double rm2, double rSigmaSqEta,
                              double rSigmaSqBar);

EMRIDGEObject( int rNlen, int rkLen, double *rXXX, double *rYYY, 
               double rSigmaSqNoise, double rpiAuse, 
               int rTotalRuns, int rNumEMSteps, 
               double rStTauSqA, double rStTauSqD, double  rTauAMult, double rTauDMult, 
               int rTauAMultStop, double *BStartBBs, double *BStartBetas,
               double *rTauSqASeq, double *rTauSqDSeq,
               int rNumCDOConv, double  rCDOEpsilon,
               int *rOrderSeq)   {
	    XTX = NULL; XTY = NULL; xxs = NULL; yys = NULL; OnBBs = NULL; OnBetas = NULL;
	    RecBBs = NULL; RecBetas = NULL;  TauSqASeq = NULL; TauSqDSeq = NULL; 
	    OrderSeq = NULL; OnInvTauSq = NULL;
	    CDO= NULL;  piSeq = NULL; SigmaSqSeq = NULL;
	    CountDescentOperations = NULL;
	   
	    
	    TauFlag = 0;
	    PrintFlag = 0;
	    SuccessFlag = 0;   
	    TotalRuns = rTotalRuns; NumEMSteps = rNumEMSteps;
	    piAuse = rpiAuse;
	    
	    if (PrintFlag > 0) {
		    Rprintf("EMRIDGEObject:NEW  Starting \n");
		    R_FlushConsole(); R_ProcessEvents();
	    }
	    if (rCDOEpsilon < 0) { 
		    CDOEpsilon = .00000001;
	    } else { CDOEpsilon = rCDOEpsilon; }
	    if (rNumCDOConv < 0) {
		    NumCDOConv = 40;
	    } else { NumCDOConv = rNumCDOConv; }
	   
	   int ii, maxCopy;
	   if (rNlen < 0) {
		   NLen = -rNlen; kLen = rkLen;
		if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW  NLen < 0, =%d, Starting with XTX \n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }
		   XTX = (double *) Calloc(rkLen * (rkLen+1) + 1, double);
		   if (XTX == NULL) {
			   Rprintf("EMRIDGEObject: Load, rNlen = %d, No space for XTX \n", rNlen);
			   R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return;
	       }
	       maxCopy = rkLen * rkLen;
	       for (ii = 0; ii < maxCopy; ii++) { XTX[ii] = rXXX[ii]; }
	       
	       XTY = (double *) Calloc(rkLen + 2, double);
		   if (XTY == NULL) {
			   Rprintf("EMRIDGEObject: Load, rNlen = %d, No space for XTY \n", rNlen);
			   R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return;
	       }
	       maxCopy = rkLen;
	       for (ii = 0; ii < maxCopy; ii++) { XTY[ii] = rYYY[ii]; }	      		   
       } else {
		   NLen = rNlen; kLen = rkLen;
		if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW  NLen > 0, =%d, Creating the XTX \n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }
		   XTX = (double *) Calloc(rkLen * (rkLen+1) + 1, double);
		   if (XTX == NULL) {
			   Rprintf("EMRIDGEObject: Load, rNlen = %d, No space for XTX \n", rNlen);
			   R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return;
	       }
	       //int SqMat (double *RetMat, int n1, int n2, double *InpMat)
	       SqMat(XTX, NLen, kLen, rXXX);
	       XTY = (double *) Calloc(rkLen + 2, double);
		   if (XTY == NULL) {
			   Rprintf("EMRIDGEObject: Load, rNlen = %d, No space for XTY \n", rNlen);
			   R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return;
	       }
           //int tMatTimesVec(int kOnLen, int NLen,  double *XX,  double *InputV, double *OutPutV)      		   	       
           tMatTimesVec(kLen, NLen, rXXX, rYYY, XTY);
       } 
        if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW OnBBs and OnBetas Create \n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }
       OnBBs = (double *) Calloc( rkLen + 2, double);
       OnBetas = (double *) Calloc( rkLen + 2, double);
       OnInvTauSq = (double *) Calloc(rkLen +2, double);
       if (OnBBs == NULL || OnBetas == NULL || OnInvTauSq == NULL)  {
	       Rprintf("EMRIDGEObject: Error, no room for OnBBs, OnBetas, or OnInvTauSq \n");
	       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
       }
        if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW  Initiating by inserting StartBetas, StartBBs \n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }       
       maxCopy = rkLen;
       for (ii = 0; ii < rkLen; ii++ ) {
	       OnBBs[ii] = BStartBBs[ii]; OnBetas[ii] = BStartBetas[ii];
       }
    
         if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW  Setup the Tau\n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }       
       SetupTau(rTauSqASeq, rTauSqDSeq, rStTauSqA, rStTauSqD, rTauAMult,
                rTauDMult, rTauAMultStop);  
       
	   //CoordinateDescentObject(int rNLen, int rkLen, double *rXXX, double *rYYY,
       //      double *rOnBeta, double rOnGamma, double *rRecordOnBetas, double *rOnGammas,
       //      double InverseGammaConstant)
        if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW  To Create the CDO object \n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }        
            // CDO = new CoordinateDescentObject(-NLen, kLen, XTX, XTY,
            //            BStartBetas, TauSqASeq[0] / 2.0, (double*) NULL, (double*) OnInvTauSq,
            //            1.0);	
             CDO = new CoordinateDescentObject(-NLen, kLen, XTX, XTY,
                        BStartBetas, TauSqASeq[0] / 2.0, (double*) NULL, 
                        (double*) OnInvTauSq, PrintFlag);	                        
       if (CDO == NULL) {
	          Rprintf("EMRIDGEObject: Load, CDO did not load \n");
	          R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return;
	               
       }
        if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW  Initiating OnBBs \n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }        
       InitiateOnBBs();
        if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW  SetupOrderSeq\n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }        
       SetupOrderSeq(rOrderSeq);
        if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:NEW  MakingXTResid \n", rNlen);
		    R_FlushConsole(); R_ProcessEvents();
	    }     
	    if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject: Printing CDO->XTX \n"); R_FlushConsole();
		    PrintRpMatrix(CDO->pXTX, kLen, kLen);
		    Rprintf("And EMRIDGEObject: Printing CDO->XTY \n"); R_FlushConsole();
		    PrintVector(CDO->XTY, kLen); R_FlushConsole();
	    }   
       CDO->MakeXTResid();  
       return; 
    }
    int SetupRecs() {
	     if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:SetupRecs() Starting \n");
		    R_FlushConsole(); R_ProcessEvents();
	    } 
           if (RecBBs != NULL) {
		       Rprintf("EMRIDGESetup:: SetupRecs(), Issue RecBBs is not NULL\n");
		       R_FlushConsole(); R_ProcessEvents();// SuccesFlag = -1; return(-1);
	        } else {
		       RecBBs = (double *)Calloc((TotalRuns+2)* (kLen+2) + 2, double);
		       if (RecBBs == NULL) {
			       Rprintf("EMRIDGESetup:: SetupRecs(), no room for RecBBs \n");
			       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		       }
	        }
	        if (RecBetas != NULL) {
		       Rprintf("EMRIDGESetup:: SetupRecs(), Issue RecBetas is not NULL\n");
		       R_FlushConsole(); R_ProcessEvents(); // SuccesFlag = -1; return(-1);
	        } else {
		       RecBetas = (double *)Calloc((TotalRuns+2)*(kLen+2)+2, double);
		       if (RecBetas == NULL) {
			       Rprintf("EMRIDGESetup:: SetupRecs(), no room for RecBetas \n");
			       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		       }
	        }	  
	     return(1);
    }
    int SetupCountDescentOperations() {
		  if (CountDescentOperations != NULL) {
			  Rprintf("CountDescentOperations is ALREADY not NULL! \n");
			  R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
	      } else {
	             CountDescentOperations = (int *)Calloc((TotalRuns*2), int);
			       if (CountDescentOperations == NULL) {
				       Rprintf("EMRIDGESetup:: SetupCountDescentOperations(), no room for CountDescentOperations \n");
				       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
			       }  
	     }
	     int ii;
	     for (ii = 0; ii < TotalRuns; ii++) { CountDescentOperations[ii] = 0;}
	     return(1);
    }
    int SetupTau(double *rTauSqASeq, double *rTauSqDSeq, 
             double rStTauSqA, double rStTauSqD, double rTauAMult,
                double rTauDMult, int rTauAMultStop) { 
	    if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:SetupTau() Starting \n");
		    R_FlushConsole(); R_ProcessEvents();
	    } 	             
	        if (TauSqASeq != NULL) {
		       Rprintf("EMRIDGESetup:: SetupTau, Issue TauSqASeq is not NULL\n");
		       R_FlushConsole(); R_ProcessEvents();// SuccesFlag = -1; return(-1);
	        } else {
		       TauSqASeq = (double *)Calloc(TotalRuns+1, double);
		       if (TauSqASeq == NULL) {
			       Rprintf("EMRIDGESetup:: SetupTau, no room for TauSqASeq \n");
			       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		       }
	        }
	        if (TauSqDSeq != NULL) {
		       Rprintf("EMRIDGESetup:: SetupTau, Issue TauSqDSeq is not NULL\n");
		       R_FlushConsole(); R_ProcessEvents(); // SuccesFlag = -1; return(-1);
	        } else {
		       TauSqDSeq = (double *)Calloc(TotalRuns+1, double);
		       if (TauSqDSeq == NULL) {
			       Rprintf("EMRIDGESetup:: SetupTau, no room for TauSqDSeq \n");
			       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		       }
	        }	  
	        int ii;      
	        if ( rTauSqASeq == NULL || rTauSqASeq[0] <= 0  ||
	             rTauSqDSeq == NULL || rTauSqDSeq[0] <= 0) {
		       if (rStTauSqA <= 0 || rStTauSqD <= 0 || 
		           rTauAMult <= 0 || rTauDMult <= 0) {
			         Rprintf("EMRIDGEObject::SetupTau  Sorry Everything is NULL! \n");
			         R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		       }
		       TauFlag = 1;
		       StTauSqA = rStTauSqA;  TauSqASeq[0] = StTauSqA;
		       StTauSqD = rStTauSqD;  TauSqDSeq[0] = StTauSqD;	
		       if (rTauAMultStop >= 0) { 
		          TauAMultStop = rTauAMultStop;
	           } else {
		          TauAMultStop = TotalRuns+1;
	           }
		       TauAMult = rTauAMult; TauDMult = rTauDMult;
		       if (TotalRuns > 1) {
			       for (ii = 1; ii < TotalRuns; ii++) {
				      if (ii  < TauAMultStop) {
					      TauSqASeq[ii] = TauSqASeq[ii-1] * TauAMult;
				      }
				      TauSqDSeq[ii] = TauSqDSeq[ii-1] * TauDMult; 
			       }
		       }
		       return(1);
	        } 	 
	        TauFlag = 2;      
		    TauSqASeq[0] = rTauSqASeq[0];
		    TauSqDSeq[0] = rTauSqDSeq[0];
		    if (TotalRuns > 1) {
			    for (ii = 1; ii < TotalRuns; ii++) {
				    TauSqASeq[ii] = rTauSqASeq[ii];
				    TauSqDSeq[ii] = rTauSqDSeq[ii];
			    }
		    }
	    if (PrintFlag > 1) {
		    Rprintf("EMRIDGEObject:SetupTau() Giving TauSqASeq: \n");
		      PrintVector(TauSqASeq, TotalRuns);
		    R_FlushConsole(); R_ProcessEvents();
		    Rprintf("EMRIDGEOBject:SetupTau()  TauSqASeq[%d-1] = %.4f\n",
		                  TotalRuns, TauSqASeq[TotalRuns-1]);
		    R_FlushConsole(); R_ProcessEvents();
		    Rprintf("EMRIDGEObject:SetupTau() Giving TauSqDSeq: \n");
		      PrintVector(TauSqDSeq, TotalRuns);
		    R_FlushConsole(); R_ProcessEvents();	
		    Rprintf("EMRIDGEOBject:SetupTau()  TauSqDSeq[%d-1] = %.4f\n",
		                  TotalRuns, TauSqDSeq[TotalRuns-1]);
		    R_FlushConsole(); R_ProcessEvents();		    	    
	    } 			    
		    return(1); 
    }  
    ~EMRIDGEObject() {
	    if (PrintFlag > 0) {
		    Rprintf("EMRIDGEObject:DELETE Starting \n");
		    R_FlushConsole(); R_ProcessEvents();
	    } 	    
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete XTX \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (XTX != NULL) { Free(XTX); XTX = NULL; }
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete XTY \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (XTY != NULL) { Free(XTY); XTY = NULL; }	
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete xxs \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (xxs != NULL) { Free(xxs); xxs = NULL; }
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete yys \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (yys != NULL) { Free(yys); yys = NULL; }	   	 
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete TauSqASeq \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (TauSqASeq != NULL) { Free(TauSqASeq); TauSqASeq = NULL; }	        
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete TauSqDSeq \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (TauSqDSeq != NULL) { Free(TauSqDSeq); TauSqDSeq = NULL; }	
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete OnBBs \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (OnBBs != NULL) { Free(OnBBs); OnBBs = NULL; }	
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete OnBetas \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (OnBetas != NULL) { Free(OnBetas); OnBetas = NULL; }	
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete RecBBs \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (RecBBs != NULL) { Free(RecBBs); RecBBs = NULL; }	
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete RecBetas \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (RecBetas != NULL) { Free(RecBetas); RecBetas = NULL; }		   	   	   	   	   
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete OnInvTauSq \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (OnInvTauSq != NULL) { Free(OnInvTauSq); OnInvTauSq = NULL; }		   
       if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete OrderSeq \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (OrderSeq != NULL) { Free(OrderSeq); OrderSeq = NULL; }		   	   	   	   	   
       if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete SigmaSqSeq \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (SigmaSqSeq != NULL) { Free(SigmaSqSeq); SigmaSqSeq = NULL; }		   	   	   	   	   
       if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete piSeq \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (piSeq != NULL) { Free(piSeq); piSeq = NULL; }		   	   	   	   	   
       if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete CountDescentOperations \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (CountDescentOperations != NULL) { Free(CountDescentOperations); CountDescentOperations = NULL; }
	   	   	   	   	   	   
	   if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: delete CDO \n"); R_FlushConsole(); R_ProcessEvents(); };
	   if (CDO != NULL) { delete(CDO); CDO = NULL; }		   	 
       if (PrintFlag >= 1) { Rprintf("~EMRDIGEObject: Deleted Everything \n");
            R_FlushConsole(); R_ProcessEvents(); }
	   return;
    }
};


