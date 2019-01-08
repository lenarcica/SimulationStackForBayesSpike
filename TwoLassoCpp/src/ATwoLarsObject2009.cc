#ifndef RMATH

  #include <R.h>
  #include <Rdefines.h>
  #include <Rinternals.h>
  #include <Rmath.h>
  #include <stdio.h>
  #define RMATH 0
#endif

#ifndef RCPPH
  #include <Rcpp.h>
  #include <cmath>
  #define RCPPH 0
#endif

#ifndef LARSOBJECT2009DD
  #include "LarsObject2009.h"
  #define LARSOBJECT2009DD 0
#endif 

#ifndef TWOLARS2009DD
  #include "TwoLarsObject2009.h"
  #define TWOLARS2009DD 0
#endif 

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
#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif 

#ifndef RDYNLOADH
  #include <R_ext/Rdynload.h>
  #define RDYNLOADH 0
#endif

#ifndef TWOLASSO2011
  #include "TwoLassoObject2014.h"
  #define TWOLASSO2011 0
#endif


int LarsConvergencyBITest(double * yys, int NLen, double *xxs, int kLen,    //4
  double *InitBetas,                                                        //5
  double *OldBetas, double *OutBetas, double *BetasKeep,                    //8
  int NumTotalRuns, int NumEMConv, double MultEMCons,                       //11
  double StlambdaD, double StlambdaA,                                       //13
  double lambdaDMultC, double lambdaAmultC,                                 //15
  int lambdaAmultStop, double ppiuse,                                       //17
  double sigmaNoiseSq, double *BBHatsKeep,                                  //19
  double *NusKeep,                                                          //20
  double *LambdaDK, double *LambdaAK,                                       //22
  double *RecordPostPBj, double DConfidenceInt, double *RecordConfidenceInts,//25 
  int LAKSeq, int *OrderSeq, int NumCDOConv, double CDOEpsilon,             //29
  int InitKKs,                                                              //30
  double InverseGammaConstant, int FixKa, double *PiRecVec,                 //36
  double *iiWeights, double TDFNu, int PrintFlag,
  double Lambda3, double *Lambda3Seq) {

Rprintf("LarsConvergencyBITest:  A test of the input data \n"); 
R_FlushConsole(); R_ProcessEvents();
Rprintf("LarsConvergencyBITest: yys[0] = %.4f, xxs[0] = %.4f, NLen = %d, kLen = %d\n", (double) yys[0], (double) xxs[0], NLen, kLen); R_FlushConsole();       
if ( NLen < 0) { 
	Rprintf("  NLen < 0, and so yys[kLen -1 = %d] = %.4f, xxs[kLen * kLen -1 = %d] = %.4f \n", kLen, yys[kLen-1], kLen * kLen -1, xxs[kLen * kLen-1]);
} else {
	Rprintf("  NLen > 0, and so yys[NLen-1= %d] = %.4f, xxs[NLen * kLen -1 = %d] = %.4f \n", NLen-1, yys[NLen-1], NLen * kLen -1, xxs[NLen * kLen-1]);
}
R_FlushConsole(); R_ProcessEvents();
Rprintf(" InitBetas[0] = %.4f, InitBetas[kLen-1 = %d] = %.4f \n", InitBetas[0], kLen-1, InitBetas[kLen-1]); R_FlushConsole();
Rprintf(" OldBetas[0] = %.4f, OldBetas[kLen-1 = %d] = %.4f \n", OldBetas[0], kLen-1, OldBetas[kLen-1]); R_FlushConsole();   
Rprintf(" OutBetas[0] = %.4f, OutBetas[kLen-1 = %d] = %.4f \n", OutBetas[0], kLen-1, OutBetas[kLen-1]); R_FlushConsole(); 
Rprintf("NumTotalRuns = %d, NumEMConv = %d, MultEMCons = %.4f \n", NumTotalRuns, NumEMConv, MultEMCons); R_FlushConsole();
Rprintf(" StlambdaD = %.4f, StlambdaA = %.4f, lambdaDMultC = %.4f, lambdaAmultC = %.4f, lambdaAmultStop = %d\n", 
   StlambdaD, StlambdaA, lambdaDMultC, lambdaAmultStop); R_FlushConsole();
Rprintf( " ppiuse = %.4f \n", ppiuse); R_FlushConsole();
Rprintf( " sigmaNoiseSq = %.4f \n", sigmaNoiseSq); R_FlushConsole();
Rprintf(" BetasKeep[0] = %.4f, BetasKeep[NumTotalRuns * kLen - 1= %d ] = %.4f \n", BetasKeep[0], BetasKeep[NumTotalRuns * kLen - 1]);
R_FlushConsole();
Rprintf(" BBHatsKeep[0] = %.4f, BBHatsKeep[NumTotalRuns * kLen - 1= %d ] = %.4f \n", BBHatsKeep[0], BBHatsKeep[NumTotalRuns * kLen - 1]);
R_FlushConsole();
Rprintf(" BBHaNusKeeptsKeep[0] = %.4f, NusKeep[NumTotalRuns * kLen - 1= %d ] = %.4f \n", NusKeep[0], NusKeep[NumTotalRuns * kLen - 1]);
R_FlushConsole();
Rprintf("  LambdaDK = %.4f \n", LambdaDK); Rprintf(  "LambdaAK = %.4f \n", LambdaAK); R_FlushConsole();
Rprintf("  LAKSeq = %d \n", LAKSeq); R_FlushConsole();
if (LAKSeq == 0) {
	Rprintf("LAKSeq == 0, Is that an error ? \n");
} else if (LAKSeq == 1) {
	Rprintf("LAKSeq == 1, Nothing Printing \n"); 
	R_FlushConsole();
} else if (LAKSeq == 2) {
	Rprintf("LAKSeq == 2,  Printing: \n"); R_FlushConsole();
	Rprintf("           LambdaDK[NumTotalRuns-1 = %d] = %.4f \n", NumTotalRuns-1, LambdaDK[NumTotalRuns-1]);
	Rprintf("           LambdaAK[NumTotalRuns-1 = %d] = %.4f \n", NumTotalRuns-1, LambdaAK[NumTotalRuns-1]);	
} else if (LAKSeq ==3) {
	Rprintf("LAKSeq == 3,  Printing: \n"); R_FlushConsole();
	Rprintf("           OrderSeq[0] = %.4f \n", OrderSeq);
	Rprintf("           OrderSeq[NumTotalRuns-1 = %d] = %.4f \n", NumTotalRuns-1, OrderSeq[NumTotalRuns-1]);	
} else if (LAKSeq ==4) {
	Rprintf("LAKSeq == 4,  Printing: \n"); R_FlushConsole();
	Rprintf("           OrderSeq[0] = %.4f \n", OrderSeq);
	Rprintf("           OrderSeq[NumTotalRuns-1 = %d] = %.4f \n", NumTotalRuns-1, OrderSeq[NumTotalRuns-1]);	
	
    Rprintf("           LambdaDK[NumTotalRuns-1 = %d] = %.4f \n", NumTotalRuns-1, LambdaDK[NumTotalRuns-1]);
	Rprintf("           LambdaAK[NumTotalRuns-1 = %d] = %.4f \n", NumTotalRuns-1, LambdaAK[NumTotalRuns-1]);		
}
R_FlushConsole();
Rprintf("DConfidenceInt = %.4f \n");
if (DConfidenceInt != -1.0) {
	Rprintf("Error, DConfidence Int is greater than zero, are we doing this?\n");
	R_FlushConsole();
}
Rprintf(" NumCDOConv = %d, CDOEpsilon = %.4e, InverseGammaConstant = %.4e\n", NumCDOConv, CDOEpsilon, InverseGammaConstant);
Rprintf(" FixK = %d \n", FixKa); R_FlushConsole(); R_ProcessEvents();
Rprintf("PiRecVec[0] = %.4e \n", PiRecVec[0]);R_FlushConsole(); R_ProcessEvents();
Rprintf("PiRecVec[TotalRuns-1] = %.4e \n", PiRecVec[NumTotalRuns-1]);
R_FlushConsole(); R_ProcessEvents();
if (iiWeights == NULL) {
   Rprintf("iiWeights is NULL\n");
} else if  (iiWeights[0] < 0) { 
   Rprintf("iiWeights[0] = %f, no iiWeights\n", (double) iiWeights[0]);
}  else {
   Rprintf("iiWeights[0] = %f and going to end: \n", iiWeights[0]);
   R_FlushConsole(); R_ProcessEvents();
   Rprintf("iiWeights[NLen-1] = %f \n", iiWeights[NLen-1]);
   R_FlushConsole();
}
Rprintf("TDFNu = %.4f\n", TDFNu);   R_FlushConsole();
Rprintf("InitKKs = %d\n", InitKKs); R_FlushConsole();
R_FlushConsole();
if (Lambda3 > 0) {
  Rprintf("Lambda3 is greater 0, = %f \n", Lambda3);
  R_FlushConsole();
}
if (Lambda3Seq != NULL && Lambda3Seq[0] > 0) {
  Rprintf("Lambda3Seq[0] = %f, and Lambda3Seq[TotalRuns-1] = .4f",
     Lambda3Seq[0], Lambda3Seq[NumTotalRuns-1]);
  R_FlushConsole();
}
Rprintf("PrintFlag = %d\n", PrintFlag); R_FlushConsole();
return(1);
}	                     
////////////////////////////////////////////////////////////////////////////
//  LarsConvergencyBI
//
//    ConvergentLARS algorithm is the EMLARS algorithm applied to a sequence of 
//   Decreasing/increasing lambda2/lambda1 values.  
//   
//  This directs the program, loads a ConvergentLARS object, and then commands the 
//   running of the appropriate program.  
//              
int LarsConvergencyBI(double * yys, int NLen, double *xxs, int kLen, 
  double *InitBetas, double *OldBetas, double *OutBetas, double *BetasKeep,
  int NumTotalRuns, int NumEMConv, double MultEMCons,
  double StlambdaD, double StlambdaA, double lambdaDMultC, 
  double lambdaAmultC, int lambdaAmultStop, double ppiuse,
  double sigmaNoiseSq, 
  double *BBHatsKeep, 
  double *NusKeep, double *LambdaDK, double *LambdaAK,
  double *RecordPostPBj, double DConfidenceInt, double *RecordConfidenceInts, 
  int LAKSeq, int *OrderSeq, int NumCDOConv, double CDOEpsilon, int InitKKs,
  double InverseGammaConstant, 
  int FixKa, double *PiRecVec, 
  double *iiWeights, double TDFNu, 
  int PrintFlag, double *pSigmaVec, 
  double Lambda3, double *Lambda3Seq) {
int IPF = 0;
int SuccFlag = 0;
if (IPF > 0) {
	 Rprintf((char*)"LarsConvergencyBI:  Loading MyCLARS \n");
	 R_FlushConsole();
	 R_ProcessEvents();
}
ConvergentLARS *MyCLARS  = NULL;
int rNLen = NLen;
if (NLen < 0) {
	rNLen = NLen;
	NLen = abs(rNLen);
}
   
if (InitKKs > kLen) {
  Rprintf("LarsConvergencyBI, ISSUE, InitKKs = %d and kLen = %d\n",
          InitKKs, kLen);
  InitKKs = kLen -1;
  R_FlushConsole();
}

if (NumCDOConv > 0) {
  MyCLARS = new ConvergentLARS(rNLen, yys, xxs, kLen, InitBetas, 
    StlambdaD, StlambdaA, 
    sigmaNoiseSq,
    lambdaDMultC, lambdaAmultC, 
    lambdaAmultStop, ppiuse, NumTotalRuns,
    InitKKs,
    InverseGammaConstant, PrintFlag);
} else {
  MyCLARS = new ConvergentLARS(yys, rNLen, xxs, kLen, InitBetas, 
    StlambdaD, StlambdaA, 
    sigmaNoiseSq,
    lambdaDMultC, lambdaAmultC, 
    lambdaAmultStop, ppiuse, NumTotalRuns, InitKKs,
    InverseGammaConstant, PrintFlag);

}
if (MyCLARS->CDO->rWXVec != NULL) {
  Rprintf("MyCLARS->Freaky, the rWXVec is Nonnull and we haven't looked at weights yet!\n");
  R_FlushConsole();
  return(-1);
}
if (PrintFlag > 2 && kLen > 5000) {
    MyCLARS->CDO->PrintFlag = 1;
}
if (Lambda3 > 0) {
   MyCLARS->CDO->OnLambda3 = Lambda3;
}
if (Lambda3Seq != NULL && Lambda3Seq[0] >= 0) {
   MyCLARS->SetupLambda3Seq(Lambda3Seq);
   MyCLARS->CDO->OnLambda3 = Lambda3Seq[0];
}
int UsedpSigmaVec= 0;
if (pSigmaVec != NULL && pSigmaVec[0] > 0) {
  //Rprintf("Going to use pSigamVec = %f\n", pSigmaVec[0]);
  MyCLARS->SigmaRecVec = pSigmaVec;  MyCLARS->LSigmaRecVec = NumTotalRuns;
  //Rprintf(" So we've set MyCLARS->SigmaRecVec to \n");
  //PrintVector(MyCLARS->SigmaRecVec, NumTotalRuns); Rprintf("\n"); R_FlushConsole();
  UsedpSigmaVec = 1;
}
//Rprintf("MyCLARS->CDO->XTXFlag = %d, OnKappaS = %d\n",
//      MyCLARS->CDO->XTXFlag, MyCLARS->CDO->OnKappaS); R_FlushConsole();
//if (MyCLARS->CDO->rWXVec == NULL) {
//  Rprintf("MyCLARS->CDO, rWXVec is NULL as it should Be \n"); R_FlushConsole();
//}     else {
//   Rprintf("MyCLARS->CDO, rWXVec is Not NULL as it shouldn't Be, not right\n");
//   R_FlushConsole();
//}
//if (MyCLARS->CDO->iiWeights == NULL) {
//  Rprintf("MyCLARS->CDO, iiWeights is NULL as it should Be \n"); R_FlushConsole();
//}  else {
//   Rprintf("MyCLARS->CDO->iiWeights is Not NULL, as it shouldn't be, not right\n");
//   R_FlushConsole();
//}
if (TDFNu > 0 && iiWeights != NULL && iiWeights[0] > 0 && rNLen > 0) {
   if (MyCLARS->IPLLF > 0) {
       Rprintf("CDO-> Setting Up TD2Lasso for df = %f", TDFNu);
       R_FlushConsole();
  }
  MyCLARS->SetupTD2Lasso(TDFNu, iiWeights);
} else if (iiWeights != NULL && iiWeights[0] > 0 && rNLen > 0) {
   if (MyCLARS->IPLLF > 0) {
       Rprintf("CDO-> Setting Up Weights EM2Lasso", TDFNu);
       R_FlushConsole();
   }
  MyCLARS->CDO->SetupWeights(iiWeights);
}
    if (DConfidenceInt > 0) { MyCLARS->DConfidenceInt = DConfidenceInt; }
int ii, ttall = kLen * NumTotalRuns; 
//MyCLARS->IPLLF = 0;                     
if (MyCLARS->OnLARS != NULL) {MyCLARS->OnLARS->PrintFlag = 0;}
if (MyCLARS->CDO != NULL) {MyCLARS->CDO->PrintFlag = 0;}
if (ppiuse < 0 && FixKa > 0 && FixKa < kLen) {
	//Rprintf(" We're going to Initiate: FixKa = %d \n", FixKa); R_FlushConsole();
	MyCLARS->InitFixKa(FixKa);
	//MyCLARS->IPLLF = 0;
} else if (ppiuse < 0.0 || ppiuse > 1.0) {
	Rprintf("LarsConvergencyAlgorithm:: ppiuse = %.4e, FixKa =%d, is incorrect, going home \n",
	    ppiuse, FixKa); R_FlushConsole(); R_ProcessEvents();
	return(-1);
}
  // Rprintf((char*) " MyCLARS->InverseGammaConstant is %.4f\n", 
  //          (double) MyCLARS->FactorGammaDenominator );
if (IPF > 0) {
	 Rprintf((char*)"LarsConvergencyBI:  Loaded, now RunAlgorithm\n");
	 R_FlushConsole();
	 R_ProcessEvents();
}
if (LAKSeq == 2 || LAKSeq == 4) {
	    MyCLARS->SetupLambdaKSeq(LambdaDK, LambdaAK);
} 
MyCLARS->LambdaSeqF = LAKSeq;
if (MyCLARS->LambdaSeqF == 3 || MyCLARS->LambdaSeqF == 4) {
	MyCLARS->SetupOrderSeq( OrderSeq );
}
SuccFlag = MyCLARS->RunConvergentLARS(NumTotalRuns, NumEMConv, MultEMCons, NumCDOConv, CDOEpsilon);
if (IPF > 0) {
	Rprintf((char*) "Finished Algorithm, SuccFlag = %d\n", SuccFlag);
	R_FlushConsole();
	R_ProcessEvents();
}

//Retrieves Weights so they are not deleted;
if (TDFNu > 0 && iiWeights != NULL && iiWeights[0] > 0 && rNLen > 0) {
  MyCLARS->iiWeights = NULL;  MyCLARS->CDO->iiWeights = NULL;
} else if (iiWeights != NULL && iiWeights[0] > 0 && rNLen > 0) {
  MyCLARS->iiWeights = NULL;  MyCLARS->CDO->iiWeights = NULL;
}
if (SuccFlag < 0) {
	Rprintf((char*) "MyCLARS->ConvergentLARS:  Flag, SuccFlag = %d\n", SuccFlag);
	R_FlushConsole();
	R_ProcessEvents();
	for (ii = 0; ii < MyCLARS->Ontt1 * MyCLARS->kLen; ii++) {
	  BetasKeep[ii] = (double) MyCLARS->BetasKeep[ii];
	  NusKeep[ii] = (double) MyCLARS->NusKeep[ii];
	  BBHatsKeep[ii] = (double) MyCLARS->BBHatsKeep[ii];	
    }
	if (MyCLARS->Ontt1 * kLen < ttall) {
    for (ii = MyCLARS->Ontt1 * kLen; ii < ttall; ii++) { 
	   BetasKeep[ii] = (double) -999.0;
	   BBHatsKeep[ii] = (double) -999.0;
	   NusKeep[ii] = (double) -999.0;
    }	
    }
    for (ii = 0; ii < kLen; ii++) { OutBetas[ii] = (double) -999.0; }
    if (kLen > 3) {
	    OutBetas[1] = (double) MyCLARS->Ontt1;
	    OutBetas[2] = (double) MyCLARS->Ontt2;
    }
	for (ii = 0; ii < NumTotalRuns; ii++) { 
		 LambdaDK[ii] = (double) -999.0;
		 LambdaAK[ii] = (double) -999.0;	
	}
	if (MyCLARS->PiRecVec != NULL && PiRecVec != NULL && PiRecVec[0] > 0.0) {
		for (ii = 0; ii < NumTotalRuns; ii++) { 
			 PiRecVec[ii] = MyCLARS->PiRecVec[0];	
		}
    }
    if (MyCLARS->RecordPostPBj != NULL) {
		ttall = NumTotalRuns * kLen;
		for (ii = 0; ii <ttall; ii++) { RecordPostPBj[ii] = (double) MyCLARS->RecordPostPBj[ii]; }
	}
	if (MyCLARS->RecordConfidenceInts != NULL) {
		ttall = NumTotalRuns * kLen *2;
		for (ii = 0; ii <ttall; ii++) { RecordConfidenceInts[ii] = 
		                                  (double) MyCLARS->RecordConfidenceInts[ii]; }	
	}
  if (UsedpSigmaVec  > 0) {
    MyCLARS->SigmaRecVec = NULL;
  }	
	Rprintf((char*)"MyCLARS->ConvergentLARS: Error Apparent, About to delete MyCLARS\n");
	R_FlushConsole();
	delete(MyCLARS);
	return(-1);  
}

for (ii = 0; ii < ttall; ii++) { 
	   BetasKeep[ii] = (double) MyCLARS->BetasKeep[ii];
	   BBHatsKeep[ii] = (double) MyCLARS->BBHatsKeep[ii];
	   NusKeep[ii] = (double) MyCLARS->NusKeep[ii];
}
for (ii = 0; ii < kLen; ii++) { OutBetas[ii] = (double) MyCLARS->OnBetas[ii]; }
 	for (ii = 0; ii < NumTotalRuns; ii++) { 
		 LambdaDK[ii] = (double) MyCLARS->LambdaDK[ii];
		 LambdaAK[ii] = (double) MyCLARS->LambdaAK[ii];	
	}
    if (MyCLARS->RecordPostPBj != NULL) {
		ttall = NumTotalRuns * kLen;
		for (ii = 0; ii <ttall; ii++) { RecordPostPBj[ii] = (double) MyCLARS->RecordPostPBj[ii]; }
	}
	if (MyCLARS->RecordConfidenceInts != NULL) {
		ttall = NumTotalRuns * kLen *2;
		for (ii = 0; ii <ttall; ii++) { RecordConfidenceInts[ii] = 
		                                  (double) MyCLARS->RecordConfidenceInts[ii]; }	
	}
  if (UsedpSigmaVec  > 0) {
    MyCLARS->SigmaRecVec = NULL;
  }	  	
delete(MyCLARS);
if (IPF > 0) {
	 Rprintf((char*)"LarsConvergencyBI:  Done Freeing, going home \n");
	 R_FlushConsole();
	 R_ProcessEvents();
}   
return(1);                 
}

//////////////////////////////////////////////////////////////////////////
// LarsConvergencyPIRECBI
//
//    This directs the EMRIDGE algorithm in the case that PI and SIGMA are
//  unknown quantities.
//                   
int LarsConvergencyPIRECBI(double * yys, int NLen, double *xxs, int kLen, 
                      double *InitBetas, 
                      double *OldBetas, double *OutBetas, double *BetasKeep,
                      int NumTotalRuns, int NumEMConv, double MultEMCons,
                      double StlambdaD, double StlambdaA, 
                      double lambdaDMultC, double lambdaAmultC,
                      int lambdaAmultStop, double *BBHatsKeep, 
                      double *NusKeep,
                      double *LambdaDK, double *LambdaAK, int LAKSeq, int *OrderSeq,
                      double *RecordPostPBj, 
                      double DConfidenceInt, double *RecordConfidenceInts,
                      double m1, double m2, int NumEConv, 
                      double SigmaSqEta, double SigmaSqBar, 
                      double *PiRecVec, double *SigmaRecVec, int NumCDOConv, double CDOEpsilon,
                      int InitKKs,
                      double InverseGammaConstant, double *iiWeights, double TDFNu,
                      int PrintFlag,
                      double *pSigmaMultVec) {
int IPF = 0;
int SuccFlag = 0;
if (IPF > 0) {
	 Rprintf((char*)"LarsConvergencyPIREC:  Loading MyCLARS \n");
	 R_FlushConsole();
	 R_ProcessEvents();
}
double sigmaNoiseSq = SigmaSqBar;
double ppiuse = m1 / (m1 + m2);

ConvergentLARS *MyCLARS = NULL;
if (InitKKs > kLen) {
  Rprintf("LarsConvergencyPIRECBI, ISSUE, InitKKs = %d and kLen = %d\n",
          InitKKs, kLen);
  InitKKs = kLen -1;
  R_FlushConsole();
}
if (NumCDOConv > 0) {
    MyCLARS = new ConvergentLARS(NLen, yys, xxs, kLen, InitBetas, 
                        StlambdaD, StlambdaA, 
                        sigmaNoiseSq, lambdaDMultC, lambdaAmultC, 
                      lambdaAmultStop, ppiuse, NumTotalRuns,
                      InitKKs,
                      InverseGammaConstant, PrintFlag);
} else {
    MyCLARS = new ConvergentLARS(yys, NLen, xxs, kLen, InitBetas, 
                        StlambdaD, StlambdaA, 
                        sigmaNoiseSq, lambdaDMultC, lambdaAmultC, 
                      lambdaAmultStop, ppiuse, NumTotalRuns,
                      InitKKs,
                      InverseGammaConstant, PrintFlag);
}
if (TDFNu > 0 && iiWeights != NULL && iiWeights[0] >= 0 && NLen > 0) {
   if (MyCLARS->IPLLF > 0) {
       Rprintf("CDO-> Setting Up TD2Lasso for df = %f", TDFNu);
       R_FlushConsole();
   }
   MyCLARS->SetupTD2Lasso(TDFNu, iiWeights);
} else if (iiWeights!= NULL &&iiWeights[0] >= 0 && NLen > 0) {
   if (MyCLARS->IPLLF > 0) {
       Rprintf("CDO-> Setting Up Weights for ConvergentLars", TDFNu);
       R_FlushConsole();
   }
   MyCLARS->CDO->SetupWeights(iiWeights);
   MyCLARS->TDFNu = -1;
}
  if (DConfidenceInt > 0) { MyCLARS->DConfidenceInt = DConfidenceInt; }                      
if (pSigmaMultVec != NULL && pSigmaMultVec[0] > 0) {
   MyCLARS->SigmaMultVec = pSigmaMultVec;
}

int ii, ttall = kLen * NumTotalRuns; 
MyCLARS->IPLLF = PrintFlag;                     
if (MyCLARS->OnLARS != NULL) {MyCLARS->OnLARS->PrintFlag = 0;}
if (MyCLARS->CDO != NULL) {MyCLARS->CDO->PrintFlag = 0; }

if (IPF > 0) {
	 Rprintf((char*)"LarsConvergencyPIREC:  Loading LambadAK \n");
	 R_FlushConsole();
	 R_ProcessEvents();
}
       if (LAKSeq == 2 || LAKSeq == 4) {
	    MyCLARS->SetupLambdaKSeq(LambdaDK, LambdaAK);
       } 
        MyCLARS->LambdaSeqF = LAKSeq;
        if (MyCLARS->LambdaSeqF == 3 || MyCLARS->LambdaSeqF == 4) {
	         MyCLARS->SetupOrderSeq( OrderSeq );
        }
  if (IPF > 0) {
	 Rprintf((char*)"LarsConvergencyPIREC:  Loaded, now RunAlgorithm\n");
	 R_FlushConsole();
	 R_ProcessEvents();
  }       
SuccFlag = MyCLARS->RunConvergentPILARS(NumTotalRuns, NumEMConv, MultEMCons, 
           m1, m2, NumEConv, SigmaSqEta, SigmaSqBar, NumCDOConv, CDOEpsilon);
if (IPF > 0) {
	Rprintf((char*) "Finished Algorithm, SuccFlag = %d\n", SuccFlag);
	R_FlushConsole();
	R_ProcessEvents();
}
//Retrieves Weights so they are not deleted;
if (TDFNu > 0 && iiWeights != NULL && iiWeights[0] > 0) {
  MyCLARS->iiWeights = NULL;  MyCLARS->CDO->iiWeights = NULL;
} else if (iiWeights != NULL && iiWeights[0] > 0) {
  MyCLARS->iiWeights = NULL;  MyCLARS->CDO->iiWeights = NULL;
}
if (SuccFlag < 0) {
	Rprintf((char*) "MyCLARS->ConvergentLARS:  Flag, SuccFlag = %d, Suggesting Delete for PiRec\n", SuccFlag);
	R_FlushConsole();
	R_ProcessEvents();
	for (ii = 0; ii < MyCLARS->Ontt1 * MyCLARS->kLen; ii++) {
	  BetasKeep[ii] = (double) MyCLARS->BetasKeep[ii];
	  NusKeep[ii] = (double) MyCLARS->NusKeep[ii];
	  BBHatsKeep[ii] = (double) MyCLARS->BBHatsKeep[ii];	
    }
	if (MyCLARS->Ontt1 * kLen < ttall) {
	    for (ii = MyCLARS->Ontt1 * kLen; ii < ttall; ii++) { 
		   BetasKeep[ii] = (double) -999.0;
		   BBHatsKeep[ii] = (double) -999.0;
		   NusKeep[ii] = (double) -999.0;
	    }	
    }
    for (ii = 0; ii < kLen; ii++) { OutBetas[ii] = (double) -999.0; }
    if (kLen > 3) {
	    OutBetas[1] = (double) MyCLARS->Ontt1;
	    OutBetas[2] = (double) MyCLARS->Ontt2;
    }
	for (ii = 0; ii < NumTotalRuns; ii++) { 
		 LambdaDK[ii] = (double) -999.0;
		 LambdaAK[ii] = (double) -999.0;	
	}
    if (MyCLARS->RecordPostPBj != NULL) {
		ttall = NumTotalRuns * kLen;
		for (ii = 0; ii <ttall; ii++) { RecordPostPBj[ii] = (double) MyCLARS->RecordPostPBj[ii]; }
	}
	if (MyCLARS->RecordConfidenceInts != NULL) {
		ttall = NumTotalRuns * kLen *2;
		for (ii = 0; ii <ttall; ii++) { RecordConfidenceInts[ii] = 
		                                  (double) MyCLARS->RecordConfidenceInts[ii]; }	
	}

  if (pSigmaMultVec != NULL && pSigmaMultVec[0] > 0) {
   MyCLARS->SigmaMultVec = NULL;
  }
	Rprintf((char*)"MyCLARS->ConvergentLARS: Error Apparent, About to delete MyCLARS\n");
	R_FlushConsole();
	delete(MyCLARS);
	return(-1);  
}
//Rprintf((char*)"MyCLARS->PrintingBBHatsKeep");
//R_FlushConsole();
// for (ii = 0; ii < NumTotalRuns; ii++) {
//	 Rprintf((char*)"\nNumRun %d : ", ii);
//	 for (int kkk = 0; kkk < kLen; kkk++) {
//		 Rprintf((char*)"%.3f ", (double) MyCLARS->BBHatsKeep[ii * kLen + kkk]);
//     }
//     R_FlushConsole();
//     R_ProcessEvents();
// }    
for (ii = 0; ii < ttall; ii++) { 
	   BetasKeep[ii] = (double) MyCLARS->BetasKeep[ii];
	   BBHatsKeep[ii] = (double) MyCLARS->BBHatsKeep[ii];
	   NusKeep[ii] = (double) MyCLARS->NusKeep[ii];
}
for (ii = 0; ii < kLen; ii++) { OutBetas[ii] = (double) MyCLARS->OnBetas[ii]; }
 	for (ii = 0; ii < NumTotalRuns; ii++) { 
		 LambdaDK[ii] = (double) MyCLARS->LambdaDK[ii];
		 LambdaAK[ii] = (double) MyCLARS->LambdaAK[ii];	
		 PiRecVec[ii] = (double) MyCLARS->PiRecVec[ii];
		 SigmaRecVec[ii] = (double) MyCLARS->SigmaRecVec[ii];
	}
    if (MyCLARS->RecordPostPBj != NULL) {
		ttall = NumTotalRuns * kLen;
		for (ii = 0; ii <ttall; ii++) { RecordPostPBj[ii] = (double) MyCLARS->RecordPostPBj[ii]; }
	}
	if (MyCLARS->RecordConfidenceInts != NULL) {
		ttall = NumTotalRuns * kLen *2;
		for (ii = 0; ii <ttall; ii++) { RecordConfidenceInts[ii] = 
		                                  (double) MyCLARS->RecordConfidenceInts[ii]; }	
	}
	if (pSigmaMultVec != NULL && pSigmaMultVec[0] > 0) {
     MyCLARS->SigmaMultVec = NULL;
  }
delete(MyCLARS);
if (IPF > 0) {
	 Rprintf((char*)"LarsConvergencyBI:  Done Freeing, going home \n");
	 R_FlushConsole();
	 R_ProcessEvents();
}   
return(1);                 
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
// RunConvergentLARS: Actually runs the EMLARS algorithm on a sequence
//
//   Update 1:  This is algorithm is designed to run either ConvergentLARS
//    or ConvergentCDLasso or Coordinate Descent Lasso.  It seems that
//    usually Coordinate Descent is a much easier algorithm to run to 
//    inspect convergence for this algorithm.  
//
//
//
int ConvergentLARS::RunConvergentLARS(int NumTotalRuns, int NumEMConv, 
          double MultEMCons, int NumCDOConv, double CDOEpsilon) {
 int cnt1 = 0;
 int NumEMCWant = NumEMConv;
 int NumCDOConvWant =NumCDOConv;
 //this->IPLLF = 0;
 int ii; int ErrorLarsFlag = -1; 
 
this->RunConvergentAlgorithmInitializeStuff(NumCDOConv, NumTotalRuns);
      
for (tt1 = 0; tt1 < NumTotalRuns; tt1++) {	
		if (SuccessFlag <0 ) { 
			Rprintf((char*) "ConvergentLARS:: SuccessFlag < 0\n");
			R_FlushConsole(); R_ProcessEvents(); return(-1);
        }
    if (this->IPLLF >1) {
	    Rprintf((char*) "RunConvergentLars: Checking out LambdaSeqF = %d \n", LambdaSeqF);
	    R_FlushConsole();
    }         
		if (this->LambdaSeqF == 0 || this->LambdaSeqF == 1 || this->LambdaSeqF == 3) {
			LambdaDK[tt1] = this->OnlambdaD;
			LambdaAK[tt1] = this->OnlambdaA;
    }       
		this->Ontt1 = tt1+1;
		if (Lambda3Seq != NULL) {CDO->OnLambda3 = Lambda3Seq[tt1];}
		if (this->IPLLF > 1) {
			Rprintf((char*)"\n RunConvergentLars: OnRun: tt1 = %d/%d, OnlambdaA = %.4f, LambdaSeqF = %d\n", 
			    tt1, NumTotalRuns, (double) OnlambdaA, LambdaSeqF);
			R_FlushConsole();
    }
  
	if (tt1 > 2) {
		NumEMCWant = ceiling( pow(MultEMCons, tt1-1) * NumEMConv);
    NumCDOConvWant = ceiling( pow(MultEMCons, tt1-1) * NumCDOConv);
		if (NumEMCWant < 1) {NumEMCWant = 2;}
  }
  if (SigmaRecVec != NULL && LSigmaRecVec > tt1) {
      sigmaNoiseSq = SigmaRecVec[tt1];
      //Rprintf("C Line 607:  tt1 = %d, Here we are moving SigmaRecVec to %f\n", tt1,sigmaNoiseSq);
      //R_FlushConsole();
  }
  //Rprintf("For tt1 = %d, sigmaNoiseSq = %f\n", tt1, sigmaNoiseSq);  R_FlushConsole();
  if (this->LambdaSeqF >= 3 && OrderSeq != NULL) { NumEMCWant = abs(OrderSeq[tt1]); }
      if (this->IPLLF > 1 && this->LambdaSeqF >= 3 && OrderSeq!= NULL ) {
	        Rprintf((char*)"\n RunConvergentLars: OnRun: Going with OrderSeq[%d] = %d \n", tt1, OrderSeq[tt1]);
      }
      if (this->IPLLF >1) {
	       Rprintf((char*) "RunConvergentLars: Ready to check tt2 \n");
	       R_FlushConsole();
      }   
      if (this->IPLLF > 0) {
	        Rprintf((char*) "Beginning of tt1 = %d/%d, LA = %.4f, LD = %.4f, Now OnBetas = ", 
	           tt1, NumTotalRuns, (double) OnlambdaA, (double) OnlambdaD);
	        PrintVector(this->OnBetas, this->kLen); Rprintf((char*)"");
	        R_FlushConsole(); R_ProcessEvents();
	        Rprintf((char*) " and CDO->OnGammas = ");
	        PrintVector(this->CDO->OnGammas, this->kLen); R_FlushConsole();
	        Rprintf((char*) " and BBOn1 = "); PrintVector(this->BBOn1, this->kLen); R_FlushConsole();
      }
  for (tt2 = 0; tt2 < NumEMCWant; tt2++) {
	  if (this->OnlambdaD == this->OnlambdaA) { tt2 = NumEMCWant; }
	   if (tt1 == 0 && IPLLF > 3) {
		   Rprintf((char*) "RunConvergentLars:: About to run only tt1 = 0 run\n");
		   R_FlushConsole();
		   if (CDO != NULL) {
		     TwoLassoPrintDiagnosticA();
       }
     }
	  this->Ontt2 = tt2+1;
     if (IPLLF > 2) {
	      Rprintf((char*)"RunConvergentLars::Current OnlambdaD = %.4f, OnlambdaA = %.4f\n", 
	                                    (double) OnlambdaD,(double) OnlambdaA);
	      R_FlushConsole();
     }
	   if (IPLLF > 2) {
	      Rprintf((char*)"RunConvergentLars::  About To Run\n", 
	                                 (double) OnlambdaD,(double) OnlambdaA);
	      R_FlushConsole();
     }            

	   if (this->LambdaSeqF <= 2 || tt2 > 0 ||
	             (this->LambdaSeqF >= 3 && OrderSeq[tt1] >= 0)) {     
		    if (IPLLF >= 4) {
		         Rprintf((char*)"RunConvergentLars::About to run OnMultiCons First\n");
	           R_FlushConsole();		          
		    } 
		    this->ConvertCreateOnMultiCons();
		    if (IPLLF >= 4) {
		        Rprintf((char*)"RunConvergentLars::Finished OnMultiCons");
            Rprintf((char*)" decide OnLARS or Convergent Descent.\n");
	          R_FlushConsole();		        
			  } 	
		    if (IPLLF >= 4) {
		        Rprintf("RunConvergentLars:: This is when life gets hairy \n");
	          R_FlushConsole();  R_FlushConsole();			        
		    } 	        
        if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL) { 
           ErrorLarsFlag=RunAlgorithmRunLARS();

				    
				    /////////////////////////////////////////////////////////
				    //  Undoing weigthst from LARS to get true Betas:
				   this->ConvertOnNusOnBetas();  
       } else if (ChoiceCDOOnLARSFlag == 2 && CDO != NULL) {
			          if (IPLLF > 2) {
				        Rprintf((char*)"RunConvergentLARS:: Setting up CDO\n");   
				        R_FlushConsole(); R_ProcessEvents(); 
			          }		                 
			          if (tt1 == 0) {                
                         this->SetupCDO(1);
                } else {
	                      this->SetupCDO(0);
                }
			          if (IPLLF > 3) {
				        Rprintf((char*)"RunConvergentLARS:: PrintInvSXX\n");  
				        R_FlushConsole(); R_ProcessEvents();
				        PrintVector(CDO->InvSXX, kLen); 
				        R_FlushConsole(); R_ProcessEvents(); 
			          }		                       
			          if (IPLLF > 2) {
					        Rprintf((char*)"RunConvergentLARS:: Running CDO\n");   
					        R_FlushConsole(); R_ProcessEvents(); 
					        Rprintf((char*)"NumCDOConvWant = %.d, CDOEpsilon = %.6f\n",
					             NumCDOConvWant, CDOEpsilon); 
					        R_FlushConsole(); R_ProcessEvents();
			          }
				   if (IPLLF > 3) {
					   Rprintf((char*)"ConvergentLARS:: Settingup CDO, Making XTResid \n"); R_FlushConsole();
				   }
				   if (CDO->ModifiedXTResid == 1) {   
				       CDO->MakeXTResid();	
			     }		          
			          //if (tt1 == 0 && tt2 == 0) {CDO->PrintFlag = 5;}
			          //else {CDO->PrintFlag = 0;}  
      if (tt1 == 0 && this->IPLLF > 3) {	
        TwoLassoPrintDiagnosticA();		          

      }		            
                 
       //////////////////////////////////////////////////////////////
       //  Run Coordinate Descent LASSO for current Weight Settings    
       SuccessFlag = CDO->RunAlgorithm(NumCDOConvWant, CDOEpsilon);
			   if (IPLLF > 3) {
					    Rprintf((char*)"RunConvergentLARS:: After CDO\n");   
					    R_FlushConsole(); R_ProcessEvents(); 
			   }  
	       if (SuccessFlag < 0) {
		          Rprintf("RunConvergetnLARS:: RunAlgorithm Returned a Fail \n");
		          R_FlushConsole();
		          R_ProcessEvents(); return(-1);
	       }
	       
        ////////////////////////////////////////////////////////
        //  E step - after CDO now take data			                             
        this->AfterCDO();
         if (SuccessFlag < 0) {
		          Rprintf("RunConvergetnLARS:: AfterCDO Returned a Fail \n");
		          R_FlushConsole();
		          R_ProcessEvents(); return(-1);
	       }
			   if (IPLLF > 3) {
					    Rprintf((char*)"RunConvergentLARS:: CDO Ending\n");   
					    R_FlushConsole(); R_ProcessEvents(); 
			   }                     
      } 
         if (IPLLF > 3) {
					   Rprintf((char*)"RunConvergentLARS:: Calculate BBOn1\n");   
					   R_FlushConsole(); R_ProcessEvents(); 
			   }                   
	       this->CalculateBBOn1();	
         
         if (IPLLF > 3) {
					   Rprintf((char*)"RunConvergentLARS:: Done Calculating BBOn1\n");   
					   R_FlushConsole(); R_ProcessEvents(); 
			   } 	                
	       if (SuccessFlag < 0) {
		         Rprintf("RunConvergetnLARS:: CalculateBBOn1 Returned a Fail \n");
		         R_FlushConsole();
		         R_ProcessEvents(); return(-1);
	       }
	       if (TDFNu > 0) {
             if (IPLLF > 3) {
					       Rprintf("RunConvergentLARS:: TDFNu = %f > 0, Going to Update\n", TDFNu);   
					       R_FlushConsole(); R_ProcessEvents(); 
			       } 	                
             SuccessFlag = RefreshTD2Lasso();
             if (SuccessFlag < 0) {
                  Rprintf("RunConvLARS: TD2Lasso Refresh Failed\n"); R_FlushConsole();
                  R_ProcessEvents();
                  return(-1);
             }
         } else {
             if (IPLLF > 3) {
					        Rprintf((char*)"RunConvergentLARS:: TDFNu = %f <= 0, didn't do anything\n", TDFNu);   
					        R_FlushConsole(); R_ProcessEvents(); 
			       } 	 
         }
	    } else if (tt2 == 0 && this->LambdaSeqF >= 3 
                        && OrderSeq != NULL && OrderSeq[tt1] < 0) {
	           // Rprintf((char*) "Running Reverse Operation, OnLA = %.4f, OnLAD = %4f\n",
	           //     (double) OnlambdaA, (double) OnlambdaD); R_FlushConsole();
               if (IPLLF > 3) {
				        Rprintf((char*)"RunConvergentLARS:: Calculate BBOn1\n"); 
                R_FlushConsole();
				        if (OrderSeq != NULL) {
					       Rprintf((char*)"OrderSeq[tt1] = %d\n", 
				                       OrderSeq[tt1]);
			            }   
				        R_FlushConsole(); R_ProcessEvents(); 
			          }  	           
	                this->CalculateBBOn1();	
	                if (SuccessFlag < 0) {
		                Rprintf("RunConvergetnLARS:: CalculateBBOn1 Returned a Fail \n");
		                R_FlushConsole();
		                R_ProcessEvents(); return(-1);
	                }	
                  if (IPLLF > 3) {
					          Rprintf((char*)"RunConvergentLARS:: Done Calculating BBOn1\n");   
					          R_FlushConsole(); R_ProcessEvents(); 
			            } 	                    
	                if (TDFNu > 0) {
                    if (IPLLF > 3) {
					            Rprintf((char*)"RunConvergentLARS:: Refreshing TD2Lasso, TDFNu = %f\n", TDFNu);   
					            R_FlushConsole(); R_ProcessEvents(); 
			              } 		                
                       RefreshTD2Lasso();
                       if (SuccessFlag < 0) {
                         Rprintf("RunConvLARS: TD2Lasso Refresh Failed\n"); R_FlushConsole();
                         R_ProcessEvents();
                         return(-1);
                       }                       
                  }                                
                if (IPLLF > 3) {
				          Rprintf((char*)"RunConvergentLARS:: Calculate Creating OnMultiCons\n"); 
				          R_FlushConsole(); R_ProcessEvents(); 
			          }  	                	            
	                this->ConvertCreateOnMultiCons();
	                if (SuccessFlag < 0) {
		                Rprintf("RunConvergetnLARS:: ConvertCreateOnMultiCons Returned a Fail \n");
		                R_FlushConsole();
		                R_ProcessEvents(); return(-1);
	                }
	               if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL) {               
					          this->OnLARS->lambda = 2 * sigmaNoiseSq;
					           if (IPLLF > 3) {
						        this->OnLARS->PrintWeights();   
					          }		        
					          this->OnLARS->LARSInit();
					          ErrorLarsFlag = this->OnLARS->LarsLassoAlgorithm2();
						        if (ErrorLarsFlag < 0) {
							        Rprintf((char*)"ErrorLarsFlag =%d\n", ErrorLarsFlag);
							        this->OnBetas[0] = -999;
							        R_FlushConsole();
							        return(ErrorLarsFlag);
						        }
						        this->ConvertOnNusOnBetas();  
                   }  else if (ChoiceCDOOnLARSFlag == 2 && CDO != NULL) {
                      if (IPLLF > 3) {
				        Rprintf((char*)"RunConvergentLARS:: Setting up CDO\n");   
				        R_FlushConsole(); R_ProcessEvents(); 
			          }
			          if (tt1 == 0) {                
                         this->SetupCDO(1);
                } else {
	                      this->SetupCDO(0);
                }
			if (IPLLF > 3) {
						   Rprintf((char*)"ConvergentLARS:: Settingup CDO, Making XTResid \n"); R_FlushConsole();
			}   
			if (CDO->ModifiedXTResid == 1) {
					       CDO->MakeXTResid();  
			}                    
			if (IPLLF > 3) {
				        Rprintf((char*)"RunConvergentLARS:: Running CDO\n");   
				        R_FlushConsole(); R_ProcessEvents(); 
				        Rprintf((char*)"NumCDOConvWant = %d, CDOEpsilon = %.6f\n",
				             NumCDOConvWant, CDOEpsilon); R_FlushConsole(); R_ProcessEvents();
			}		       
 
      if (tt1 == 0 && this->IPLLF > 3) {
        TwoLassoPrintDiagnosticA();
      }	
      if (IPLLF > 2) {
          Rprintf("RunConvergentLARS:: Now Run CDO Algorithm\n"); R_FlushConsole();
          R_ProcessEvents();
      }				                         
      CDO->RunAlgorithm(NumCDOConvWant, CDOEpsilon);
			if (IPLLF > 3) {
				        Rprintf((char*)"RunConvergentLARS:: After CDO\n");   
				        R_FlushConsole(); R_ProcessEvents(); 
			}                     
      this->AfterCDO();
    } else {
             SuccessFlag = -1;
             Rprintf((char*) "Is CDO Non NULL or OnLARS\n"); R_FlushConsole();
             R_ProcessEvents(); return(-1);
    }
    if (IPLLF > 3) {
	       Rprintf((char*)"RunConvergentLARS:: Calculate BBOn1 after OrderSeqRun\n");
	       R_FlushConsole(); R_ProcessEvents();
    } 
         this->CalculateBBOn1();				                       	            
	            
    }	 else {
	      Rprintf((char*) "RunConvergentLARS:Not Running Anything\n"); 
        R_FlushConsole();
	      Rprintf((char*) "The LambdaSeqF == %d\n", LambdaSeqF);
	      if (OrderSeq == NULL) {
		         Rprintf((char*) "OrderSeq is NULL\n");
	      }
	      Rprintf((char*) "OrderSeq[tt1 = %d] is %d\n", tt1, OrderSeq[tt1]);
	      R_FlushConsole(); R_ProcessEvents();
	         	                    
    }
	  if (IPLLF > 3) {
		      Rprintf((char*)"RunConvergentLARS: Making it to end of a run\n");
		      R_FlushConsole();
		      if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL) {
			      Rprintf((char*) "Printing OnLARS, Weights\n");
		          this->OnLARS->PrintWeights();   
	          } else if (ChoiceCDOOnLARSFlag == 2 && this->CDO != NULL) {
		          Rprintf((char*) "Printing CDO->OnGammas\n"); R_FlushConsole();
		          if (this->CDO->OnGammas != NULL) {
		              PrintVector(this->CDO->OnGammas, this->kLen);
	              }
              } else {
	              Rprintf((char*)"RunConvergentLARS:: End Run, should print something but can't\n");
	              R_FlushConsole(); R_ProcessEvents();
	              SuccessFlag = -1; return(-1);
              }
	        }
	  
        if (this->IPLLF > 2) {
		          Rprintf("End of tt1=%d, tt2 = %d/%d\n",tt1,  tt2, NumTotalRuns);
		          R_FlushConsole();
		          AllLARSCDODiagnositc();
		    }
    }
        //return(1);
		    if (LambdaSeqF == 1 || LambdaSeqF == 3) {
		       OnlambdaD = OnlambdaD * MultDConst;
			     if (this->IPLLF > 1 ) {
				      Rprintf("StopAInt = %d, tt1 = %d, MultAConst = %.4f, OnlambdaA = %.4f\n",
				        StopAInt, tt1, (double) MultAConst, (double) OnlambdaA );
			     }
			     if (StopAInt > tt1) {
			         OnlambdaA = OnlambdaA * MultAConst;
			     }
		    }  else if (LambdaSeqF == 2|| LambdaSeqF == 4) {
			     if (tt1 >= NumTotalRuns -1) { 
			     } else { 
			        OnlambdaD = this->LambdaDK[tt1+1];
			        OnlambdaA = this->LambdaAK[tt1+1];
		       }
		    }
		        
		    // If keeping a running tally, BetasKeep holds information    
		    if (this->IPLLF > 2) {
			    Rprintf((char*) "RunConvergentLARS:: Filling in Keep Variables\n");
			    R_FlushConsole(); R_ProcessEvents();
		    }
	      for (ii = 0; ii < this->kLen; ii++) {
		        this->NusKeep[cnt1] = (double) this->OnNus[ii];
		        this->BetasKeep[cnt1] = (double) this->OnBetas[ii];
		        this->BBHatsKeep[cnt1] = (double) this->BBOn1[ii];
		        cnt1++;
	      }
		   if (this->IPLLF > 2) {
			    Rprintf((char*) "RunConvergentLARS:: Filling in Confidence Intervals\n");
			    R_FlushConsole(); R_ProcessEvents();
		  }
       //  Code attempts to get confidence intervals for this content 	        
	     if (this->DConfidenceInt >  0 && this->DConfidenceInt < 1.0 ) {	        
		        this->CalculateCurrentPostPBj();
		        this->SetCopyCurrentPostPBj(tt1);
		        this->CalculateCurrentConfidenceInts(this->DConfidenceInt);
		        this->SetCopyCurrentConfidenceInts(tt1);
	     }

   }
   return(1);	
}


/////////////////////////////////////////////////////////////////////////////////////
// RunConvergentPILARS
//    Running algorithm that simulatneously estimates PI and Sigma with prior information
//
//
int ConvergentLARS::RunConvergentPILARS(int NumTotalRuns, int NumEMConv, double MultEMCons, 
    double rm1, double rm2, int NumEConv, double SigmaEta, double SigmaBarSq, int NumCDOConv
	, double CDOEpsilon)  {
	int cnt1 = 0;
	int NumEMCWant = NumEMConv;
  int NumCDOConvWant = NumCDOConv;

  int ii;
  int ErrorLarsFlag = -1;
     
  this->m1 = (double) rm1; this->m2 = (double) rm2;
  this->ppiuse = this->m1 / (this->m1 + this->m2);
  this->sigmaNoiseSq = SigmaBarSq;
  this->BBHatsNewStep = NULL;
	this->BBHatsNewStep = (double *)Calloc( this->kLen+2, double);
	if (this->BBHatsNewStep == NULL) {
		  Rprintf((char*)"LarsConv Cannot Apportion BBHatsNewStep \n");
		  R_FlushConsole();
		  return(-1);
	}	
  this->SigmaRecVec = NULL;
	this->SigmaRecVec = (double *)Calloc( NumTotalRuns+2, double);
	if (this->SigmaRecVec == NULL) {
	   Rprintf((char*)"LarsConv Cannot Apportion SigmaRecVec \n");
     R_FlushConsole();
     return(-1);
	}	
  this->PiRecVec = NULL;
	this->PiRecVec = (double *)Calloc( NumTotalRuns+2, double);
	if (this->PiRecVec == NULL) {
	     Rprintf((char*)"LarsConv Cannot Apportion PiRecVec \n");
	     R_FlushConsole();
		   return(-1);
	}		     	     
    
	int iijjkk = 0;
  // if (LambdaDK != NULL) {
  //	    Free(LambdaDK);
  //  } 
  //  if (LambdaAK != NULL) {
  //	    Free(LambdaDK);
  //  }
  if (LambdaDK == NULL) { 
    LambdaDK = (double *) Calloc(NumTotalRuns, double);
     if (LambdaDK == NULL) {
	    Rprintf((char*)"RunCLARS: LambdaDK = NULL\n");
	    R_FlushConsole(); SuccessFlag = -1; return(-1);
    }
    LambdaDK[0] = -999;
  }
  if (LambdaAK == NULL) {
    LambdaAK = (double *) Calloc(NumTotalRuns, double);
    if (LambdaAK == NULL) {
	    Rprintf((char*)"RunCLARS: LambdaAK = NULL\n");
	    R_FlushConsole(); SuccessFlag = -1; return(-1);
    } 
    LambdaAK[0] = -999;   
  }
  
    if (this->IPLLF > 4) {
	    Rprintf((char*) "RunConvergentLARS, RealSet, here's the YY\n");
	    int dii, djj;
	    if (Nlen < 20) {
	    for (dii = 0; dii < Nlen; dii++) {
		    Rprintf((char*) "%d=%.4f,", dii, (double) this->yys[dii]);
	    }
        }
	    if (Nlen < 20) {
	    Rprintf((char*) "\n And here is the XX\n");
	    for (dii = 0; dii < Nlen; dii++) {
		    Rprintf((char*)" %d = [ ", dii);
		    for (djj=0; djj < kLen; djj++) {
			    Rprintf((char*)" %.3f, ", (double) RealXXS[djj * Nlen + dii]);
		    }
		    Rprintf((char*)"\n");
	    }
        }
	    R_FlushConsole();
    }
    if (this->IPLLF > 0) {
	    Rprintf((char*) "RCL: NumTotalRuns = %d, NumEMConv = %d \n", NumTotalRuns, NumEMConv);
	    R_FlushConsole();
    }
    if (LambdaSeqF == 1 || LambdaSeqF == 3) { 
	    this->OnlambdaD = this->StlambdaD;
        this->OnlambdaA = this->StlambdaA;
    } else {
	    this->OnlambdaD = LambdaDK[0];
	    this->OnlambdaA = LambdaAK[0];
    }
    if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL && this->kLen < 200) {
        this->GenGAll();
        this->UpdateGAll(0);
    }
    for (iijjkk = 0; iijjkk < this->kLen; iijjkk++) {
	    this->BBOn1[iijjkk] = this->ppiuse;
    }
	for (tt1 = 0; tt1 < NumTotalRuns; tt1++) {	
		if (LambdaSeqF == 1 || LambdaSeqF == 3) {
		   LambdaDK[tt1] = this->OnlambdaD;
		   LambdaAK[tt1] = this->OnlambdaA;
	    } else {
		    this->OnlambdaD = LambdaDK[tt1];
		    this->OnlambdaA = LambdaAK[tt1];
	    }
		if (SuccessFlag <0 ) { 
			Rprintf((char*) "ConvergentLARS:: RunPILARSAlgorithm, SuccessFlag < 0\n");
			R_FlushConsole(); R_ProcessEvents(); return(-1);
        }		
		this->Ontt1 = tt1+1;
		//this->IPLLF = 2;
		if (this->IPLLF > 1) {
			Rprintf((char*)"CL: RunPILARSAlgorrithm::  within tt1, OnRun: tt1 = %d/%d, OnlambdaA = %.4f\n", 
			    tt1, NumTotalRuns, (double) OnlambdaA);
			R_FlushConsole();
        }
  
		if (tt1 > 2) {
			NumEMCWant = ceiling( pow(MultEMCons, tt1-1) * NumEMConv);
                  NumCDOConvWant = ceiling( pow(MultEMCons, tt1-1) * NumCDOConv);
			if (NumEMCWant < 1) {NumEMCWant = 2;}
        }
        for (tt2 = 0; tt2 < NumEMCWant; tt2++) {
		if (this->IPLLF > 3) {
			Rprintf((char*)"CL: RunPILARSAlgorrithm::  within tt2, OnRun: tt2 = %d/%d,\f\n", 
			    tt2, NumEMCWant, (double) OnlambdaA);
			R_FlushConsole();
        }	        
	       this->Ontt2 = tt2+1;
	       // if (tt2 == 0 && tt1 == 0) {
		   //     this->OnLARS->PrintFlag = 2;
	       // } else {
		   //     this->OnLARS->PrintFlag = 2;
	       // }
	        if (IPLLF > 3) {
	           Rprintf((char*)"CL Run PILARSAlgorithm:: Current OnlambdaD = %.4f, OnlambdaA = %.4f\n", 
	                                    (double) OnlambdaD,(double) OnlambdaA);
	           R_FlushConsole();
            }
	        //this->OnLARS->lambda = 2 * sigmaNoiseSq;
	        /*if (IPLLF > 3) {
		        Rprintf((char*)"Printing BBHatOn1, kLen = %d\n", this->kLen);
		        for (ii = 0; ii < this->kLen; ii++) {
			        Rprintf((char*) "%d=%.3f, ", ii, (double) this->BBOn1[ii]);
		        }
		        Rprintf((char*)"\n");
		        R_FlushConsole();
	        }*/	
	        if (this->LambdaSeqF <= 2 || tt2 > 0 || OrderSeq == NULL || 
	             (this->LambdaSeqF >= 3 && OrderSeq[tt1] >= 0)) {   
		        if (IPLLF > 3) {
		           Rprintf((char*)"Run Algorithm, OrderSeq > 0 ConvertCreateOnMultiCons\n", 
		                                    (double) OnlambdaD,(double) OnlambdaA);
		           R_FlushConsole();
	            }		                       
	          this->ConvertCreateOnMultiCons();
              if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL) {
	                 this->OnLARS->lambda = 2 * sigmaNoiseSq;

	  if (IPLLF > 3 && ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL) { 
        Rprintf((char*)"RunConvergentPIREC: running OnLARS->PrintWeights\n");
		                  R_FlushConsole(); R_ProcessEvents();
		                  this->OnLARS->PrintWeights();   
	  }
	             
	  if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS!= NULL) {
		    this->OnLARS->LARSInit();
	  }
       
	  ErrorLarsFlag = this->OnLARS->LarsLassoAlgorithm2();
	  if (ErrorLarsFlag < 0) {
		  Rprintf((char*)"ErrorLarsFlag =%d\n", ErrorLarsFlag);
		               this->OnBetas[0] = -999;
		   R_FlushConsole();
		   return(ErrorLarsFlag);
	  }

	  this->ConvertOnNusOnBetas();
  }  else if (ChoiceCDOOnLARSFlag == 2 && this->CDO != NULL) {
	      if (IPLLF > 3) { 
		            Rprintf("RunPILARS:  about to SetupCDO()\n");
		            R_FlushConsole(); R_ProcessEvents();
	      }
	      if (tt1 == 0) {
            this->SetupCDO(1);
        } else {
	          this->SetupCDO(0);
        }
        if (CDO->ModifiedXTResid == 1) {
						CDO->MakeXTResid();   
				}                  
     if (IPLLF > 3) {                
       Rprintf("RunPILARS:  about to RunAlgorithm()\n");
		   R_FlushConsole(); R_ProcessEvents();          
	   }           
     CDO->RunAlgorithm( NumCDOConvWant, CDOEpsilon );
     if (IPLLF > 3) {
		    Rprintf("RunPILARS:  about to AfterCDO()\n");
		    R_FlushConsole(); R_ProcessEvents();   
	   }                  
      this->AfterCDO();
     if (IPLLF > 3) {
		    Rprintf("RunPILARS:  Finished with CDO()\n");
		    R_FlushConsole(); R_ProcessEvents();
	   }                     
   } else {
	    Rprintf((char*)"EMPILARS:RunConvergentPILARS, No CDO/OnLARS choice\n");
	    R_FlushConsole(); R_ProcessEvents();
      SuccessFlag = -1;
      return(-1);
   }
	        //ReviseBHatAll(NumEConv, this->kLen, this->BBOn1, this->BBHatsNewStep, 
            //   this->OnBetas,  this->OnlambdaA, this->OnlambdaD, 
            //   this->m1, this->m2);  
      if (IPLLF > 3) {
	       Rprintf((char*) "RunPILARSAlgorithm, about to run BBOn1 \n");
	       R_FlushConsole(); R_ProcessEvents();
      }  
	 this->CalculateBBOn1();
      if (IPLLF > 3) {
	       Rprintf((char*) "RunPILARSAlgorithm, about to run UpdateMaxPi \n");
	       R_FlushConsole(); R_ProcessEvents();
      }  	            
	 this->UpdateMaxPi();
      if (IPLLF > 3) {
	       Rprintf((char*) "RunPILARSAlgorithm, about to run SigmaSq \n");
	       R_FlushConsole(); R_ProcessEvents();
      }  		              
	    if (TDFNu > 0) {
          if (IPLLF > 3) {
					    Rprintf("RunConvergentPiLARS:: TDFNu = %f > 0, Going to Update\n", TDFNu);   
					    R_FlushConsole(); R_ProcessEvents(); 
			    } 	                
          SuccessFlag = RefreshTD2Lasso();
          if (SuccessFlag < 0) {
              Rprintf("RunConvPiLARS: TD2Lasso Refresh Failed\n"); R_FlushConsole();
              R_ProcessEvents();
              return(-1);
          }
      } else {
          if (IPLLF > 3) {
					    Rprintf((char*)"RunConvergentPiLARS:: TDFNu = %f <= 0, didn't do anything\n", TDFNu);   
					    R_FlushConsole(); R_ProcessEvents(); 
			    } 	 
      }		         
              if (SigmaEta >= 0) {
		            this->UpdateSigmaSq(SigmaEta, SigmaBarSq);
		            if (SuccessFlag < 0) {
			            Rprintf((char*) "WB::SuccessFlag went below 0 after UpdateSigmaSq\n");
			                Rprintf((char*) " tt1 = %d, tt2 = %d\n", tt1, tt2);
			                R_FlushConsole(); R_ProcessEvents();
			                return(-1);
		                }
	            }
        
	  }   else if ( tt2 == 0 && OrderSeq != NULL &&
	             (this->LambdaSeqF >= 3 && OrderSeq[tt1] < 0)) { 
               if (IPLLF > 3) {
	                Rprintf((char*) "RunPILARSAlgorithm, about to run BBOn1 \n");
	                R_FlushConsole(); R_ProcessEvents();
                }  		             
            this->CalculateBBOn1();	
            if (IPLLF > 3) {
	                Rprintf((char*) "RunPILARSAlgorithm, about to run CreateOnMultiCons \n");
	                R_FlushConsole(); R_ProcessEvents();
                }
	                                	                         
	        this->ConvertCreateOnMultiCons();
              if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL) {
	                 this->OnLARS->lambda = 2 * sigmaNoiseSq;
	        /*if (IPLLF > 3) {
		        Rprintf((char*)"Printing CurrentMultiCons, kLen = %d\n", this->kLen);
		        for (ii = 0; ii < this->kLen; ii++) {
			        Rprintf((char*) "%d=%.3f, ", ii, (double)(double) this->OnMultiCons[ii]);
		        }
		        Rprintf((char*)"\n");
		        Rprintf((char*) "  MC = c( ");
		        for (ii = 0; ii < this->kLen-1; ii++) {
			        Rprintf((char*) "%.8f, ", (double)(double) this->OnMultiCons[ii]);
		        }
		        Rprintf((char*) " %.8f ) \n", this->OnMultiCons[this->kLen-1]);
		        R_FlushConsole();
	        }*/
	                 if (IPLLF > 3 && this->OnLARS!= NULL) {
		                  this->OnLARS->PrintWeights();   
	                   }
	              //this->OnLARS->PrintFlag = 0;
	                this->OnLARS->LARSInit();
	               /*
	                if (IPLLF > 3) {
		             Rprintf((char*)"Printing WDiag, kLen = %d\n", this->kLen);
		             for (ii = 0; ii < this->kLen; ii++) {
			          Rprintf((char*) "%d=%.3f, ", ii, (double)this->OnLARS->WDiag[ii]);
		             }
		             Rprintf((char*)"\n");
		             Rprintf((char*) "  WDiag = c( ");
		             for (ii = 0; ii < this->kLen-1; ii++) {
			           Rprintf((char*) "%.8f, ", (double)(double) this->OnLARS->WDiag[ii]);
		             }
		             Rprintf((char*) " %.8f ) \n", this->OnLARS->WDiag[this->kLen-1]);		        
		             R_FlushConsole();
	                }*/	        
	                ErrorLarsFlag = this->OnLARS->LarsLassoAlgorithm2();
	                if (ErrorLarsFlag < 0) {
		               Rprintf((char*)"ErrorLarsFlag =%d\n", ErrorLarsFlag);
		               this->OnBetas[0] = -999;
		               R_FlushConsole();
		               return(ErrorLarsFlag);
	                 }

	                  /* if (this->IPLLF > 3) {
		           Rprintf((char*)"FinishLLA, End of tt1=%d, tt2 = %d/%d\n",tt1,  tt2, NumTotalRuns);
		           Rprintf((char*)"lambda = %.4f, OnlambdaD=%.4f, OnlambdaA=%.4f\n", 
		                        (double) this->OnLARS->lambda,
		                        (double) this->OnlambdaD, (double) this->OnlambdaA);
		           Rprintf((char*)"  Now OnNus \n");
	                 for (ii = 0; ii < this->kLen; ii++) {
		              Rprintf((char*)"%d=%.4f,  ", ii, (double) this->OnLARS->OnBetas[ii]);
	                 }
		           Rprintf((char*) "\n  OnNus = c( ");
		           for (ii = 0; ii < this->kLen-1; ii++) {
			          Rprintf((char*) "%.8f, ", (double)(double) this->OnLARS->OnBetas[ii]);
		           }
		           Rprintf((char*) " %.8f ) \n", this->OnLARS->OnBetas[this->kLen-1]);	          
	                 R_FlushConsole();
	                 Rprintf((char*)"\n");
                     }*/
	                 this->ConvertOnNusOnBetas();
                   }  else if (ChoiceCDOOnLARSFlag == 2 && CDO != NULL) {
					   if (IPLLF > 3) {
						   Rprintf((char*)"ConvergentLARS:: Settingup CDO \n");
						    R_FlushConsole(); R_FlushConsole();
					   }   	                   
			          if (tt1 == 1) {                
                         this->SetupCDO(1);
                } else {
	                      this->SetupCDO(0);
                }
					   if (IPLLF > 3) {
						   Rprintf((char*)"ConvergentLARS:: Settingup CDO, Making XTResid \n"); R_FlushConsole();
					   }   
					   if (CDO->ModifiedXTResid == 1) {
						   CDO->MakeXTResid();   
				       }                     
                     CDO->RunAlgorithm( NumCDOConvWant, CDOEpsilon );
                     this->AfterCDO();
                   } else {
	                   Rprintf((char*)"EMPILARS:RunConvergentPILARS, No CDO/OnLARS choice\n");
	                   R_FlushConsole(); R_ProcessEvents();
                       SuccessFlag = -1;
                       return(-1);
                   }
	        //ReviseBHatAll(NumEConv, this->kLen, this->BBOn1, this->BBHatsNewStep, 
            //   this->OnBetas,  this->OnlambdaA, this->OnlambdaD, 
            //   this->m1, this->m2);    
	    this->CalculateBBOn1();
	    if (TDFNu > 0) {
	        if (IPLLF > 3) {
              Rprintf("RunPiALG:  TDFNu = %f, Updateing TD2\n", TDFNu);
              R_FlushConsole();
          }
          RefreshTD2Lasso();
          if (SuccessFlag < 0) {
              Rprintf("RunConvLARS: TD2Lasso Refresh Failed\n"); R_FlushConsole();
              R_ProcessEvents();
              return(-1);
          }                 
      }	            
	    this->UpdateMaxPi();  
	    if (SigmaEta >= 0) {
		      this->UpdateSigmaSq(SigmaEta, SigmaBarSq);
	    }	 
	    if (SuccessFlag < 0) {
		      Rprintf((char*) "SuccessFlag went below 0 after UpdateSigmaSq\n");
		      Rprintf((char*) " tt1 = %d, tt2 = %d\n", tt1, tt2);
		      R_FlushConsole(); R_ProcessEvents();
		      return(-1);
	    }
    }    
      if (this->IPLLF > 2) {
		        Rprintf((char*)"End of tt1=%d, tt2 = %d/%d\n",tt1,  tt2, NumTotalRuns);
		        R_FlushConsole();
		        AllLARSCDODiagnositc();
	    }
  }  // End of tt2 EM Loop
      
  // Update Lambdas    
  if (LambdaSeqF == 1 || LambdaSeqF == 3) {
		  OnlambdaD = OnlambdaD * MultDConst;
			if (this->IPLLF > 1 ) {
				  Rprintf((char*) "StopAInt = %d, tt1 = %d, MultAConst = %.4f, OnlambdaA = %.4f\n",
				      StopAInt, tt1, (double) MultAConst, (double) OnlambdaA );
			}
			if (StopAInt > tt1) {
			    OnlambdaA = OnlambdaA * MultAConst;
			}
	 }  else if (LambdaSeqF == 2|| LambdaSeqF == 4) {
			if (tt1 >= NumTotalRuns -1) { 
			} else { 
			    OnlambdaD = this->LambdaDK[tt1+1];
			    OnlambdaA = this->LambdaAK[tt1+1];
		  }
	 }
   for (ii = 0; ii < this->kLen; ii++) {
	    this->NusKeep[cnt1] = (double) this->OnNus[ii];
	    this->BetasKeep[cnt1] = (double) this->OnBetas[ii];
	    this->BBHatsKeep[cnt1] = (double) this->BBOn1[ii];
	    cnt1++;
   }
   if (this->IPLLF > 2) {
	    Rprintf((char*)"End of tt1 = %d/%d, L1 = %.4f, L2=%.4f\n", 
	        tt1, NumTotalRuns, (double) this->OnlambdaD, (double) this->OnlambdaA);
	    R_FlushConsole();
	    R_ProcessEvents();
   }
   if (this->DConfidenceInt >  0 && this->DConfidenceInt < 1.0 ) {	        
	   this->CalculateCurrentPostPBj();
	   this->SetCopyCurrentPostPBj(tt1);
	   this->CalculateCurrentConfidenceInts(this->DConfidenceInt);
	   this->SetCopyCurrentConfidenceInts(tt1);
   }
     this->PiRecVec[tt1] = this->ppiuse;
        //for (iijjkk = 0; iijjkk < this->kLen; iijjkk++) { this->PiRecVec[tt1] += this->BBOn1[iijjkk]; }
        //this->PiRecVec[tt1] = (this->PiRecVec[tt1] + this->m1) /
        //                         (this->kLen + this->m1 + this->m2);
        //this->ppiuse = this->PiRecVec[tt1];
     this->SigmaRecVec[tt1] = this->sigmaNoiseSq;
     if (this->IPLLF > 3) {
        Rprintf((char*)"End of tt1 = %d/%d\n", tt1, NumTotalRuns);
        Rprintf((char*)"  Now OnBetas \n");
        for (ii = 0; ii < this->kLen; ii++) {
	          Rprintf((char*)"%d=%.4f,  ", ii, (double) this->OnBetas[ii]);
        }
        R_FlushConsole();
        Rprintf((char*)"\n");
        Rprintf((char*)"  Now BBOn1 \n");
        for (ii = 0; ii < this->kLen; ii++) {
	           Rprintf((char*)"%d=%.4f,  ", ii, (double) this->BBOn1[ii]);
        }
        Rprintf((char*)"\n");
        R_FlushConsole();          
        R_FlushConsole();
      }
      if (this->IPLLF > 0) {
	        Rprintf((char*) "End of tt1 = %d/%d, the Current OnBetas =\n", tt1, NumTotalRuns);
	        R_FlushConsole(); R_ProcessEvents();
	        PrintVector(this->OnBetas, this->kLen);
      }
   }
   return(1);	
}

/////////////////////////////////////////////////////////////////////
///  Sub ConvergentLARS functions
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//   Set DiagXTX
//     DiagXTX is a vector containing only the diagonals of 
//     the GRAM matrix XTX.  Since this may be used for weighting
//     Procedures independent of GramXTX, it may be useful to calculate
//     simply the diagonal alone, instead of full matrix, unless of course
//     the matrix itself is already calculated.
int ConvergentLARS::SetDiagXTX() {
       if (this->OnLARS == NULL) {
        Rprintf("ConvergentLARS:: SetDiagXTX, OnLARS NULL, why bother???\n"); 
        SuccessFlag = -1; R_FlushConsole(); R_ProcessEvents(); return(-1);
       }
  if (DiagXTX != NULL) { FFree(&DiagXTX); }
	DiagXTX = (double *) Calloc( this->kLen+1, double);
	if (DiagXTX == NULL) {
		Rprintf((char*) "Error DiagXTX is Null, cannot apportion Space \n");
		R_FlushConsole();
		return(-1);
  }
  int ii;
	if (this->GAll != NULL)	 {
		ii = 0;
		for (ii = 0; ii < this->kLen; ii++) {
		      DiagXTX[ii] = this->GAll[ (ii) * kLen + ii ];	
        }
        return(1);
  }
 // double RunTotal;
  int cn1, kk;
  int One = 1;  double OneD = 1.0; int Zero = 0; double ZeroD = 0.0;
  if (iiWeights == NULL) {
      for (ii = 0; ii < kLen; ii++ ) {
         DiagXTX[ii] = F77_CALL(ddot)(&Nlen, RealXXS + ii * Nlen, &One, 
             RealXXS + ii * Nlen, &One);
      }
  }  else {
      if (kLen < Nlen) {
          for (ii = 0; ii < kLen; ii++) {
              DiagXTX[ii] = 0; cn1 = this->Nlen * ii;
              for (kk = 0; kk < Nlen; kk++) {
                DiagXTX[ii] += RealXXS[cn1] * RealXXS[cn1] * iiWeights[kk];
                cn1++;
              }
          }
      }  else {
         F77_CALL(dscal)(&kLen, &ZeroD, DiagXTX, &One);
         for (ii = 0; ii < Nlen; ii++) {
             F77_CALL(dsbmv)("U", &kLen, &Zero, iiWeights + ii, RealXXS + ii * Nlen, &One,
                       RealXXS + ii * Nlen, &One, &OneD,DiagXTX + ii, &One);          
         }
      }
  }
// for (ii = 0; ii < kLen; ii++) {
//	     cn1 = this->Nlen * ii;
//	     RunTotal = (double) 0;
//	     if (iiWeights != NULL) {
//          
//       }
//	     DiagXTX[ii]
//	     for (kk = 0; kk < Nlen; kk++) {
//		     RunTotal += (double) RealXXS[cn1] * (double) RealXXS[cn1];
//		     cn1++;
//		     cn2++;
//	     }
//	     DiagXTX[ii] = (double) RunTotal;
//    }
    return(1);
}

////////////////////////////////
//   SetupCDO
//
//    Work presumed necessary to start CDO
//
//
int ConvergentLARS::SetupCDO(int SetBeta) {
   if (CDO == NULL) {Rprintf((char*)"ConvergentLARS::SetupCDO, CDO is NULL\n");
                     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1); }
   if (IPLLF > 0 && SetBeta > 0) {
	   Rprintf((char*) "ConvergentLARS::SetupCDO, IPLLF = %d, SetBeta = %d \n",  IPLLF, SetBeta);
	   R_FlushConsole();
   }
   if (IPLLF > 2) {
	   Rprintf((char*)"ConvergentLARS:: Settingup CDO \n"); R_FlushConsole();
   }
   int ii;
   if (FactorGammaDenominator != 1.0) {
	   for (ii = 0; ii < kLen; ii++) {
	     CDO->OnGammas[ii] = (BBOn1[ii] * OnlambdaA + 
	                          (1.0- BBOn1[ii]) *OnlambdaD)  * 
	                          sigmaNoiseSq * 
	        FactorGammaDenominator;
	   }
   } else {
       for (ii = 0; ii < kLen; ii++) {
	     CDO->OnGammas[ii] = (BBOn1[ii] * OnlambdaA + 
	                          (1.0- BBOn1[ii]) *OnlambdaD)  * 
	                          sigmaNoiseSq;
	   }
   }
   if (SetBeta == 1) {
     int One = 1;
     F77_CALL(dcopy)(&kLen, OnBetas, &One, CDO->OnBeta, &One);
	   //for (ii = 0; ii < kLen; ii++) {
	   //  CDO->OnBeta[ii] = OnBetas[ii];
	   //}	   
     CDO->ModifiedXTResid = 1;
   }
   if (IPLLF > 2) {
	   Rprintf((char*) "ConvergentLARS:: SetupCDO The CDO OnGammas Are : \n"); R_FlushConsole();
	   PrintVector(CDO->OnGammas, kLen);
     Rprintf((char*) "ConvergentLARS:: SetupCDO , initial OnBeta Are : \n"); R_FlushConsole();
     PrintVector(CDO->OnBeta, kLen);
     Rprintf((char*) "ConvergentLARS:: Active Beta are : \n"); R_FlushConsole();
     int tti;
     if (CDO->OnKappaS < kLen && CDO->OnKappaS > 0 && CDO->kXFinder != NULL) {
       for (tti = 0; tti < CDO->OnKappaS; tti++) {
          Rprintf(" %d: %d  = %.4f \n", tti, CDO->kXFinder[tti], 
            CDO->OnBeta[CDO->kXFinder[tti]]);
          R_FlushConsole();
       }
     }
     R_FlushConsole(); R_ProcessEvents();
   }

   if (IPLLF > 2) {
	   Rprintf((char*)"ConvergentLARS:: Settingup CDO, finished XTResid \n"); R_FlushConsole();
   }   
   if (IPLLF > 1) {
      Rprintf("ConvergentLARS:: Ending SetupCDO\n"); R_FlushConsole();
   }
   return(1);
}

int ConvergentLARS::AfterCDO() {
   if (CDO == NULL) {Rprintf((char*)"ConvergentLARS::AfterCDO, CDO is NULL\n");
                     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1); }
  //int ii;
   if (IPLLF > 2) {
	   Rprintf((char*)"ConvergentLARS:: After CDO \n");
   } 
     int One = 1;
     F77_CALL(dcopy)(&kLen, CDO->OnBeta, &One, OnNus, &One);     
     F77_CALL(dcopy)(&kLen, CDO->OnBeta, &One, OnBetas, &One);         
  // for (ii = 0; ii < kLen; ii++) {
  //      OnNus[ii] = CDO->OnBeta[ii];
  //      OnBetas[ii] = CDO->OnBeta[ii];
  // }
   if (IPLLF > 2) {
	   Rprintf((char*) "ConvergentLARS:: AfterCDO The CDO OnGammas Are : \n"); R_FlushConsole();
	   PrintVector(CDO->OnGammas, kLen);
       Rprintf((char*) "ConvergentLARS:: AfterCDO , final OnBeta Are : \n"); R_FlushConsole();
       PrintVector(CDO->OnBeta, kLen);
   }   
   return(1);

}


/////////////////////////////////////////////////////////
//  int ConvergentLARS::GenGAll() ;
//  GAll is the Gram Matrix X^TX
//    Much like LARS builtin, one can choose to calculate
//    OrUpdate based upon whim or n * k size.
// GALL is matrix of all Covariate factors, only run once
int ConvergentLARS::GenGAll()   {
   int ii,jj,kk;
   FFree(&GAll);
   this->GAll = NULL;
   this->GAll = (double *) Calloc( this->kLen* this->kLen+1, double);
   if (this->GAll == NULL) {
	   Rprintf((char*)"GenGAll Fail to allocate Data space \n");
	   R_FlushConsole();
	   return(-1);
   }
   double RunTotal;
   int cn1, cn2;
   //SqMat
   SqMat(GAll, Nlen, kLen, RealXXS);
   return(1);
   for (ii = 0; ii < kLen; ii++) {
	  cn1 = this->Nlen * ii;
	    RunTotal = (double) 0;
	    for (kk= 0; kk < Nlen; kk++) {
		     RunTotal += (double) RealXXS[cn1] * (double) RealXXS[cn1];
		     cn1++;
		     cn2++;
	     }
	     GAll[ ii * this->kLen + ii] = (double) RunTotal;
	    if (ii < kLen-1) {
		      for (jj = ii+1; jj < kLen; jj++) {
			     cn1 = this->Nlen * ii;
			     cn2 = this->Nlen * jj;
			     RunTotal = 0;
			     for (kk= 0; kk < Nlen; kk++) {
				     RunTotal += (double) RealXXS[cn1] * (double) RealXXS[cn2];
				     cn1++;
				     cn2++;
			     }
			     GAll[ ii * this->kLen + jj] = RunTotal;
			     GAll[ jj * this->kLen + ii] = RunTotal;
		    }
        }
  }
  return(1);	     
} 

/////////////////////////////////////////////////////////
//  int ConvergentLARS::UpdateGAll(int unweight);
//  GAll is the Gram Matrix X^TX
//    Much like LARS builtin, one can choose to calculate
//    OrUpdate based upon whim or n * k size.
int ConvergentLARS::UpdateGAll(int unweight) {
  if (this->OnLARS == NULL) {
     Rprintf((char*) "ConvergentLARS::UpdateGAll, Why do you want to update GAll CDO?\n");
     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
  }
  if (this->OnLARS->GAll == NULL) {
	  this->OnLARS->GAll = (double *)Calloc( this->kLen * this->kLen, double);
	  if (this->OnLARS->GAll == NULL) {
		  Rprintf((char*)"UpdateGall, No Memory\n");
		  R_FlushConsole();
		  return(-1);
      }
     if (this->OnLARS->XTYAll == NULL) {
	     this->OnLARS->XTYAll = (double *)Calloc( this->kLen, double);
	     if (this->OnLARS->XTYAll == NULL) {
		     Rprintf((char*)"UpdateGall, No memory for XTYAll\n");
		     R_FlushConsole();
		     Free(this->OnLARS->GAll);
		     this->OnLARS->GAll = NULL;
		     return(-1);
	     }
     }
  }
  int cnt1, ii, jj;
  if (unweight == 0) {	  
	  cnt1 = 0;
	  int kLenMult = this->kLen * this->kLen; int One = 1;
	  F77_CALL(dcopy)(&kLenMult, this->GAll, &One, this->OnLARS->GAll, &One);
	  //for (ii = 0; ii < this->kLen; ii++) {
		//  for (jj = 0; jj < this->kLen; jj++) {
		//     this->OnLARS->GAll[cnt1] = (double) this->GAll[cnt1] ;
		//     cnt1++;
	  //    }
	  //}	 
	  cnt1 = 0; 
    tMatTimesVec(this->kLen, this->Nlen,  this->RealXXS,  
                 this->yys, this->OnLARS->XTYAll);  
    //F77_CALL(dsymv)("U", &kLen, OneD,
		//    this->RealXXS, ,
		//const double *x, const int *incx,
		//const double *beta, double *y, const int *incy);
	  //F77_CALL(dcopy)(&kLen, this->GAll, &One this->OnLARS->GAll, &One);	  
	  //for (ii = 0; ii < this->kLen; ii++) { this->OnLARS->XTYAll[ii] = 0; }
	  //for (ii = 0; ii < this->kLen; ii++) {
		//  for (jj = 0; jj < this->Nlen; jj++) {
		//      this->OnLARS->XTYAll[ii] += (double) this->RealXXS[cnt1] * 
		//            (double) this->yys[jj];
		//      cnt1++;
	  //    }
    //  }
      this->OnLARS->GFlag = 1;
  } else {
	      cnt1 = 0;
		  for (ii = 0; ii < this->kLen; ii++) {
			  for (jj = 0; jj < this->kLen; jj++) {
			     this->OnLARS->GAll[cnt1] = (double) this->GAll[cnt1] / (OnMultiCons[ii]  * OnMultiCons[jj] );
			     cnt1++;
		      }
		  }
	  cnt1 = 0;
    tMatTimesVec(this->kLen, this->Nlen,  this->RealXXS,  
                 this->yys, this->OnLARS->XTYAll); 
                 	
    for (ii = 0; ii < kLen; ii++) {
      this->OnLARS->XTYAll[ii] = this->OnLARS->XTYAll[ii] / OnMultiCons[ii];
    }
	  //for (ii = 0; ii < kLen; ii++) { this->OnLARS->XTYAll[ii] = 0; }
	  //for (ii = 0; ii < kLen; ii++) {
		//  for (jj = 0; jj < Nlen; jj++) {
		//      this->OnLARS->XTYAll[ii] += (double) this->RealXXS[cnt1] *  this->yys[jj];
	  //	      cnt1++;
	  //    }
	  //    this->OnLARS->XTYAll[ii] = this->OnLARS->XTYAll[ii] / OnMultiCons[ii];
    //  }
      this->OnLARS->GFlag = 1;		  
  }	   
  return(1);
}

///////////////////////////////////////////////////
//  int ConvergentLARS::UpdateWAll(int weightFlag)
//
//  Updates Weight vector (Diagonal Weight Matrix)
//     Inside the OnLARS object within the ConvergentLARS object.
//     Only if weightFlag == 1 does this operation currently work.
//
int ConvergentLARS::UpdateWAll(int weightFlag) {
   int ii;
  if (this->OnLARS == NULL) {
     Rprintf((char*)"ConvergentLARS::UpdateWAll, why UpdateWAll when doing CDO???\n");
     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
  }
  if (weightFlag == 0) {
	  return(0);
  } else if (weightFlag == 1) {
	  if (this->OnLARS->WDiag == NULL) {
		 this->OnLARS->WDiag = (double *)Calloc( this->kLen +1, double);
			  if (this->OnLARS->WDiag == NULL) {
				  Rprintf((char*)"UpdateWDiag, No Memory\n");
				  R_FlushConsole();
				  return(-1);
		      }
      }
      this->OnLARS->WDiagFlag = 1; this->OnLARS->AvWeights = 0.0;
      for (ii = 0; ii < this->kLen; ii++) {
	       this->OnLARS->WDiag[ii] = (double) 1.0 / this->OnMultiCons[ii];
	       this->OnLARS->AvWeights += this->OnLARS->WDiag[ii];
      }
      this->OnLARS->AvWeights = this->OnLARS->AvWeights / this->kLen;
      return(1);
  } else if (weightFlag == 2) {
	  Rprintf((char*)"WeightFlag = 2 I do not know what to do at this time\n");
	  if (this->OnLARS->WAll == NULL) {
		  this->OnLARS->WAll = (double *)Calloc( this->kLen * this->kLen, double);
		  if (this->OnLARS->WAll == NULL) {
			  Rprintf((char*)"UpdateGall, No Memory\n");
			  R_FlushConsole();
			  return(-1);
	      }
	  }
	  this->OnLARS->WDiagFlag = 2;
	  return(2);
  }
  return(-1);
}

///////////////////////////////////////////////////////////////
//   ConvertCreateOnMultiCons
//
//    If OnMultiCons = BBOn1[jj] * Onlambda_A + ( 1- BBOn1[jj] ) * Onlambda_D
//    Then this creates weight vector, sets up weight vector in the OnLARS
//    And then weights for the running of the algorithm.  
 void ConvergentLARS::ConvertCreateOnMultiCons() {
	  int jj;      int One = 1;
	  //double InvMultiCons;  int cnt1;
    if (this->OnLARS != NULL) {
       for (jj = 0; jj < kLen; jj++) {
			    //cnt1 = jj * Nlen;
			    OnMultiCons[jj] =  (BBOn1[jj] * OnlambdaA + ( 1.0 - BBOn1[jj]) * OnlambdaD );
      }
      if (SigmaMultVec != NULL  && SigmaMultVec[tt1] > 0 && SigmaMultVec[tt1] != 0) {
        F77_CALL(dscal)(&kLen, SigmaMultVec + tt1, OnMultiCons, & One);
      }	
		   for (jj = 0; jj < kLen; jj++) {
			    //InvMultiCons =1.0/ OnMultiCons[jj];
			    OnNus[jj] = OnBetas[jj] * OnMultiCons[jj];
			    OnLARS->OnBetas[jj] = 0;
		  }
		  UpdateWAll(1);
		  this->OnLARS->lambda = 2.0 * this->sigmaNoiseSq; 
   }  else if (CDO != NULL) {
	    if (FactorGammaDenominator != 1) { 
          for (jj = 0; jj < kLen; jj++) {
	              //cnt1 = jj * Nlen;
	              OnMultiCons[jj] = (BBOn1[jj] * OnlambdaA + (1.0-BBOn1[jj]) * OnlambdaD );
	              OnNus[jj] = OnBetas[jj];
	              CDO->OnBeta[jj] = OnBetas[jj];       
	              CDO->OnGammas[jj] =   this->sigmaNoiseSq * OnMultiCons[jj] 
	                                                      * FactorGammaDenominator;
			    }	          
		  } else {
          for (jj = 0; jj < kLen; jj++) {
	          //cnt1 = jj * Nlen;
	          OnMultiCons[jj] = (BBOn1[jj] * OnlambdaA + (1.0-BBOn1[jj]) * OnlambdaD );
	          OnNus[jj] = OnBetas[jj];
	          CDO->OnBeta[jj] = OnBetas[jj];       
	          CDO->OnGammas[jj] =   this->sigmaNoiseSq * OnMultiCons[jj];
			}
      if (SigmaMultVec != NULL  && SigmaMultVec[tt1] > 0 && SigmaMultVec[tt1] != 0) {
        F77_CALL(dscal)(&kLen, SigmaMultVec + tt1, CDO->OnGammas, & One);
      }			    
	 }
  } else {
            Rprintf((char*)"ConvergentLARS::ConvertCreateOnMultiCons, Error, not CDO or OnLARS\n");
            R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1;
  }
	return;
}
     
     
     
///////////////////////////////////////////////////////////////////
//  UpdateMaxPi()
//    For those who need to find maximum pi_A, this asses the result
//   using the current BBOn1 estimates.
//
//
//
int ConvergentLARS::UpdateMaxPi() {
   double SupBHat = 0;
   int ii;
   for (ii = 0; ii < this->kLen; ii++) {
	    SupBHat += BBOn1[ii];
   }
   this->ppiuse = (m1 +1.0+ SupBHat) / ( m1 + m2 + this->kLen-2.0);
   if (this->ppiuse < 0) { 
	   this->ppiuse = (m1) / ( m1 + m2 + this->kLen);
   } else if (this->ppiuse > 1) {
	   this->ppiuse = (m1 + this->kLen) / ( m1 + m2 + this->kLen);
   }
   return(1);
 }	 

///////////////////////////////////////////////////////////////
// SetCopyCurrentPostPBj
//
//  An experimental method for finding confidence intervals for 
//   BBOn1[jj] values.  Does not appear to work, usually as confident
//   for a given Bvalue as the value itself.  
//	 
 int ConvergentLARS::SetCopyCurrentPostPBj(int Run) {
     int jjSt = Run * this->kLen;
     int ii = 0, iimax = this->kLen;
     for (ii = 0; ii < iimax; ii++) {
	      RecordPostPBj[jjSt + ii] = CurrentPostPBj[ii];
     }
     return(1);
 }
 
//////////////////////////////////////////////////////////
//  CurrentPostPBj
//
//   Experimental: Tries to calculate a confidence interval set.
 int ConvergentLARS::CalculateCurrentPostPBj() {
	  double bCurSumXTXBnj;
	  double aCurXTXjj;
	  int jj;
	  double logsqrtpi = .5 * log( M_PI);
	  double LS1A, LS1B, LS2A, LS2B;
	  double CurgammaD, CurgammaA;
          if (this->OnLARS == NULL) {
           Rprintf((char*) "ConvergentLARS::CalculateCurrentPostPBj, can't do CDO yet\n");
           R_FlushConsole(); R_ProcessEvents(); return(-1);
          }
	  for (jj = 0; jj < this->kLen; jj++) {
	    aCurXTXjj = this->DiagXTX[jj] / ( 2 * sigmaNoiseSq);
	    bCurSumXTXBnj = ((double)this->OnLARS->cc[jj] + 
	                     (double) DiagXTX[jj] * (double) this->OnBetas[jj]
	                    )  / ((double)sigmaNoiseSq);
	    LS1B = (bCurSumXTXBnj - OnlambdaD) *   (bCurSumXTXBnj - OnlambdaD) / (4 * aCurXTXjj);    
	    LS1A = Rf_pnorm5(   (bCurSumXTXBnj - OnlambdaD) / sqrt( 2 * aCurXTXjj), 
	                 0, 1, 1, 1);   
	    LS2B = (bCurSumXTXBnj + OnlambdaD) *   (bCurSumXTXBnj + OnlambdaD) / (4 * aCurXTXjj);    
	    LS2A = Rf_pnorm5(   -(bCurSumXTXBnj + OnlambdaD) / sqrt( 2 * aCurXTXjj), 
	                 0, 1, 1, 1); 
	    CurgammaD = exp( LS1A + LS1B + logsqrtpi - .5 * log(aCurXTXjj)) +  
	                exp( LS2A + LS2B + logsqrtpi - .5 * log(aCurXTXjj));
	    LS1B = (bCurSumXTXBnj - OnlambdaA) *   (bCurSumXTXBnj - OnlambdaA) / (4 * aCurXTXjj);    
	    LS1A = Rf_pnorm5(   (bCurSumXTXBnj - OnlambdaA) / sqrt( 2 * aCurXTXjj), 
	                 0, 1, 1, 1);   
	    LS2B = (bCurSumXTXBnj + OnlambdaA) *   (bCurSumXTXBnj + OnlambdaA) / (4 * aCurXTXjj);    
	    LS2A = Rf_pnorm5(   -(bCurSumXTXBnj + OnlambdaA) / sqrt( 2 * aCurXTXjj), 
	                 0, 1, 1, 1); 
	    CurgammaA = exp( LS1A + LS1B + logsqrtpi - .5 * log(aCurXTXjj)) +  
	                exp( LS2A + LS2B + logsqrtpi - .5 * log(aCurXTXjj));	
	    this->CurGammaA[jj] = CurgammaA; this->CurGammaD[jj] = CurgammaD;  
	    CurrentPostPBj[jj] = CurgammaA * OnlambdaA * ppiuse /
	                  (CurgammaA * OnlambdaA * ppiuse + CurgammaD * OnlambdaD * (1-ppiuse));              
	                 		  
      }
      return(1);
 }     
 
///////////////////////////////////////////////////////////////////////////
//  int ConvergentLARS::SetCopyCurrentConfidenceInts(int Run)
//
//    Expeirmental Confidence ints.  Just records this is a vector.
// 
 int ConvergentLARS::SetCopyCurrentConfidenceInts(int Run) {
     int jjSt = Run * 2 * this->kLen;
     int ii = 0, iimax = this->kLen * 2;
     for (ii = 0; ii < iimax; ii++) {
	      RecordConfidenceInts[jjSt + ii] = CurrentConfidenceInts[ii];
     }
     return(1);
 }
 
////////////////////////////////////////////////////////////////////////////
//  int ConvergentLARS::CalculateCurrentConfidenceInts(double DesiredCInt)
//
//  This is a confidence interval for Betas.  But because confidence intervals 
//   for BBOn1 cannot be reliably calculated the result is the same.  
// 
 int ConvergentLARS::CalculateCurrentConfidenceInts(double DesiredCInt) {
	  double bCurSumXTXBnj;
	  double aCurXTXjj;
	  int jj;
	  double logsqrtpi = .5 * log( M_PI);
	  double LS1A, LS1B, LS2A, LS2B, varrho1APlus, varrho1AMinus;
	  double CurBetaHat1;
	  //double PBj0;
    double PBj1;
	  double PBetaG0LBetaHat;
	  double HDesiredCInt = DesiredCInt / 2.0;
	  double MyCalculator = 0;
	  double pnormBetaHatMinusMeanLeftTail, pnormBetaHatPlusMeanRightTail;
	  double pnormPlusTailZero, pnormMinusTailZero;
	   for (jj = 0; jj < this->kLen; jj++) {
	    aCurXTXjj = DiagXTX[jj] / ( 2 * sigmaNoiseSq);
	    bCurSumXTXBnj = ((double)this->OnLARS->cc[jj] + 
	                     (double) DiagXTX[jj] * (double) this->OnBetas[jj]
	                    )  / ((double)sigmaNoiseSq);
	    CurBetaHat1 = (double) OnBetas[jj];                
	    //PBj0 = 1 - this->CurrentPostPBj[jj];
	    PBj1 = this->CurrentPostPBj[jj];
	    
	    LS1B = (bCurSumXTXBnj - OnlambdaA) *   (bCurSumXTXBnj - OnlambdaA) / (4 * aCurXTXjj);  
	    LS1A = Rf_pnorm5(   (bCurSumXTXBnj - OnlambdaA) / sqrt( 2 * aCurXTXjj), 
	                 0, 1, 1, 1);   
	    pnormPlusTailZero = LS1A;  
	    LS2B = (bCurSumXTXBnj + OnlambdaA) *   (bCurSumXTXBnj + OnlambdaA) / (4 * aCurXTXjj);    
	    LS2A = Rf_pnorm5(   -(bCurSumXTXBnj + OnlambdaA) / sqrt( 2 * aCurXTXjj), 
	                 0, 1, 1, 1); 
	    pnormMinusTailZero = LS2A;  	                 
	    varrho1APlus =   exp(LS1A + LS1B + logsqrtpi - .5 * log(aCurXTXjj)  - log(CurGammaA[jj]) );
	    varrho1AMinus =  exp(LS2A + LS2B + logsqrtpi - .5 * log(aCurXTXjj) - log(CurGammaA[jj]) );	    
		   pnormBetaHatMinusMeanLeftTail = Rf_pnorm5( (2 * aCurXTXjj * CurBetaHat1 - 
		                    (bCurSumXTXBnj+ OnlambdaA)) / 
		                    sqrt(2 * aCurXTXjj), 0,1, 0, 1);
		   pnormBetaHatPlusMeanRightTail = Rf_pnorm5 ( ( 
		                    (bCurSumXTXBnj- OnlambdaA- 2.0 * aCurXTXjj * CurBetaHat1)) / 
		                    sqrt(2 * aCurXTXjj), 0,1, 0, 1);		                    
		    if (CurBetaHat1 > 0) {
			   PBetaG0LBetaHat = (pnormPlusTailZero - pnormBetaHatPlusMeanRightTail  ) /  
			                    pnormPlusTailZero  * varrho1APlus * PBj1;
	                MyCalculator = HDesiredCInt / (varrho1APlus * PBj1) * 
				                       		                pnormPlusTailZero ;
			        MyCalculator = MyCalculator + (1.0-pnormBetaHatPlusMeanRightTail); 
			        MyCalculator = Rf_qnorm5(MyCalculator,0,1, 1, 0);
			        MyCalculator = MyCalculator * sqrt(2.0 * aCurXTXjj) + ( bCurSumXTXBnj - OnlambdaA);
			        MyCalculator = MyCalculator / (2 * aCurXTXjj);
				       CurrentConfidenceInts[ 2 * jj+1] = MyCalculator;				                    
			   if (PBetaG0LBetaHat > HDesiredCInt ) {
			        MyCalculator = HDesiredCInt / (varrho1APlus	   * PBj1) *  
			                    pnormPlusTailZero;
			        MyCalculator = MyCalculator + pnormBetaHatPlusMeanRightTail;
			        MyCalculator = Rf_qnorm5(MyCalculator, 0, 1, 1,0);
			        MyCalculator = -MyCalculator * sqrt( 2.0 * aCurXTXjj) + ( bCurSumXTXBnj - OnlambdaA);
			        MyCalculator = MyCalculator / (2.0 * aCurXTXjj);
				    CurrentConfidenceInts[ 2 * jj] = MyCalculator;
				              
		       } else if (PBetaG0LBetaHat + (1-PBj1) > HDesiredCInt) {
			       CurrentConfidenceInts[ 2* jj ] = 0.0;
	           } else {
		           MyCalculator = HDesiredCInt- (PBetaG0LBetaHat + (1.0-PBj1));
		           MyCalculator = MyCalculator / ( varrho1AMinus * PBj1) * pnormMinusTailZero;
		           MyCalculator = MyCalculator + (1.0- pnormMinusTailZero);
		           MyCalculator = sqrt(2.0 * aCurXTXjj) * Rf_qnorm5(MyCalculator,0,1,1,0);
		           MyCalculator = -MyCalculator + ( bCurSumXTXBnj + OnlambdaA);
		           MyCalculator = MyCalculator / (2.0 * aCurXTXjj);
		           CurrentConfidenceInts[2 * jj] = MyCalculator;
	            }
	        } else if (CurBetaHat1 == 0) {
	                MyCalculator = (HDesiredCInt - (1.0-PBj1))/ (varrho1APlus	   * PBj1) *  
			                    pnormPlusTailZero;
			        MyCalculator = MyCalculator + (1.0- pnormPlusTailZero);
			        MyCalculator = Rf_qnorm5(MyCalculator, 0, 1,1,0);
			        MyCalculator = MyCalculator * sqrt( 2.0 * aCurXTXjj) + ( bCurSumXTXBnj - OnlambdaA);
			        MyCalculator = MyCalculator / (2.0 * aCurXTXjj);
				    CurrentConfidenceInts[ 2 * jj+1] = MyCalculator;
				   MyCalculator = HDesiredCInt - (1.0-PBj1);
	               MyCalculator = MyCalculator / ( varrho1AMinus * PBj1) * pnormMinusTailZero;
		           MyCalculator = -MyCalculator + pnormMinusTailZero;
		           MyCalculator = sqrt(2.0 * aCurXTXjj) * Rf_qnorm5(MyCalculator,0,1,1,0);
		           MyCalculator = MyCalculator + ( bCurSumXTXBnj + OnlambdaA);
		           MyCalculator = MyCalculator / (2.0 * aCurXTXjj);
		           CurrentConfidenceInts[2 * jj] = MyCalculator;	
	        } else {		 
		           PBetaG0LBetaHat = (  pnormMinusTailZero - pnormBetaHatMinusMeanLeftTail) /  
			                    pnormMinusTailZero  * varrho1APlus * PBj1;   
		           MyCalculator = HDesiredCInt;
	               MyCalculator = MyCalculator / ( varrho1AMinus * PBj1) * pnormMinusTailZero;
		           MyCalculator = pnormBetaHatMinusMeanLeftTail - MyCalculator;
		           MyCalculator = sqrt(2 * aCurXTXjj) * Rf_qnorm5(MyCalculator,0,1,1,0);
		           MyCalculator = MyCalculator + ( bCurSumXTXBnj + OnlambdaA);
		           MyCalculator = MyCalculator / (2 * aCurXTXjj);
		           CurrentConfidenceInts[2 * jj+1] = MyCalculator;	
		       if (PBetaG0LBetaHat > HDesiredCInt ) {
			        MyCalculator = HDesiredCInt / (varrho1AMinus	   * PBj1) *  
			                    pnormMinusTailZero;
			        MyCalculator = MyCalculator + (pnormBetaHatMinusMeanLeftTail);
			        MyCalculator = Rf_qnorm5(MyCalculator, 0, 1,1,0);
			        MyCalculator = MyCalculator * sqrt( 2 * aCurXTXjj) + ( bCurSumXTXBnj + OnlambdaA);
			        MyCalculator = MyCalculator / (2 * aCurXTXjj);
				    CurrentConfidenceInts[ 2 * jj] = MyCalculator;				              
		       } else if (PBetaG0LBetaHat + (1-PBj1) > HDesiredCInt) {
			       CurrentConfidenceInts[ 2* jj+1 ] = 0.0;
	           } else {
		           MyCalculator = HDesiredCInt- (PBetaG0LBetaHat + (1.0-PBj1));
		           MyCalculator = MyCalculator / ( varrho1APlus * PBj1) * pnormPlusTailZero;
		           MyCalculator = MyCalculator + (1.0-pnormPlusTailZero);
		           MyCalculator = sqrt(2 * aCurXTXjj) * Rf_qnorm5(MyCalculator,0,1,1,0);
		           MyCalculator = MyCalculator + ( bCurSumXTXBnj - OnlambdaA);
		           MyCalculator = MyCalculator / (2 * aCurXTXjj);
		           CurrentConfidenceInts[2 * jj+1] = MyCalculator;
	            }
	          }                
       }
	 
	 return(1);
 }

///////////////////////////////////////////////////////////////////////////////
//  ReviseBHatAll
//
//   Directs calculation of new BHat for all B vectors individually. 
// 
int ReviseBHatAll(int iters, int kLen, double *CurBHat, double *NewBHat, 
               double *CurBetas,  double LambdaA, double LambdaD, 
               double m1, double m2) {
	int nn, ii;
    //Rprintf((char*)"ReviseBHatAll: Printing BBHatAll As First\n");
    //for (ii = 0; ii < kLen; ii++) {
	//    Rprintf((char*)" %.4f",(double) CurBHat[ii]);
	//    if (ii < kLen-1) {Rprintf((char*)",");}
    // }
    //Rprintf((char*)"\n");
    //R_FlushConsole();	
    //Rprintf((char*)"ReviseBHatAll: Printing CurBetas As Given\n");
    //for (ii = 0; ii < kLen; ii++) {
	//    Rprintf((char*)" %.4f", (double) CurBetas[ii]);
	//    if (ii < kLen-1) {Rprintf((char*)",");}
    //}
    //Rprintf((char*)"\n");
    //R_FlushConsole();	    
    for (nn = 0; nn < iters; nn++) {
	    ReviseBHat1( kLen, CurBHat, NewBHat, CurBetas, LambdaA, LambdaD, m1, m2);
	    for (ii = 0; ii < kLen; ii++) {
		     CurBHat[ii] = NewBHat[ii];
	    }	    
    }
    //Rprintf((char*)"ReviseBHatAll: Printing BBHatAll As Estimated\n");
    //for (ii = 0; ii < kLen; ii++) {
	//    Rprintf((char*)" %.4f", (double) CurBHat[ii]);
	//    if (ii < kLen-1) {Rprintf((char*)",");}
    //}
    //Rprintf((char*)"\n");
    //R_FlushConsole();
    return(1);
}

/////////////////////////////////////////////////////////////////////////////
//  ReviseBHat1
//  
//    The BHat vector has to be expected one by one.  
//
int ReviseBHat1(int kLen, double *CurBHat, double *NewBHat, 
               double *CurBetas,  double LambdaA, double LambdaD, 
               double m1, double m2) {
   double CurSBBar = 0;
   double Betasj1, Betasj2;
   //double = lBetasj1, lBetasj2, BetajR;
   int ii;
   for (ii = 0; ii < kLen; ii++) {
	   CurSBBar += CurBHat[ii];
   }
   Betasj2 = 1;
   for (ii = 0; ii < kLen; ii++) {
	   //lBetasj1 = beta( (double) (m1 + CurSBBar - CurBHat[ii] + 1), 
	   //                (double) (m2 + kLen - CurSBBar+CurBHat[ii] - 1));
       //lBetasj2 = beta( (double) (m1 + CurSBBar - CurBHat[ii]), 
       //                (double) (m2 + kLen - CurSBBar+CurBHat[ii]));
       //BetajR = (double) ((int) lBetasj1);
       //Betasj1 = exp( lBetasj1 - BetajR);
       //Betasj2 = exp( lBetasj2 - BetajR);
       Betasj1 = (m2 - CurSBBar + CurBHat[ii] + kLen -1) /
                 (m1 + CurSBBar - CurBHat[ii] );
       //Rprintf((char*)"ReviseBHat1, m1 = %.4f, m2 = %.4f, Betasj1 = %.4f, Betasj2 = %.4f, a1 = %.4f\n", 
       //        (double) m1, (double) m2, (double) Betasj1, (double) Betasj2,
       //        (double) (m1 + CurSBBar - CurBHat[ii] + 1)) ;
       NewBHat[ii] = Betasj1 * LambdaA * exp(- LambdaA * fabs((double) CurBetas[ii]) );
       NewBHat[ii] = NewBHat[ii] / ( NewBHat[ii] + 
                       Betasj2 * LambdaD * exp(-LambdaD * fabs((double) CurBetas[ii]) ) );
   }
   return(1);
}


/////////////////////////////////////////////////////////////////////////////
//  UpdateSigmaSq()
// We're really putting into sigmaNoiseSq 1 / E[ 1 / sigmaNoiseSq | Beta ]
// AKA We don't need sigmaNoiseSq itself, but the expected inverse.
int ConvergentLARS::UpdateSigmaSq(double prioralpha, double priormean) {
	//Rprintf((char*)"UpdateSigmaSq: this->SumCurResids = %.4f \n", 
	//    (double) this->OnLARS->SumCurResids);
	//R_FlushConsole();
  double CurSBBar = 0;
  int ii;
  for (ii = 0; ii < kLen; ii++) {
	  CurSBBar += this->BBOn1[ii]; 
  }
  //CurSBBar -=  - kLen * this->ppiuse * this->OnlambdaA
  //	                       (this->ppiuse * this->OnlambdaA + (1-this->ppiuse) * this->OnlambdaD);
  //this->sigmaNoiseSq = (this->OnLARS->SumCurResids + priormean * prioralpha) / 
  //                   (this->Nlen	 + prioralpha - CurSBBar);
    double TotalResids = 0.0;
    double TempSumYYSq;  double SumRDSquared;
    int One = 1;   int Zero = 0; 
    double OneD = 1; double ZeroD = 0.0;
    double SumWeightsii = 0;
     double OtherSumRDS = 0; 
    //int ii1; 
    //int jj1;
    // double OnYResid =0;  
    // int MarkMe = 0;
     
  //ii1 = 1; 
  //jj1 = 1;
  OtherSumRDS++; 
  //MarkMe = 0; OnYResid = 0;
  if (this->TDFNu > 0 && this->CDO != NULL && kLen > Nlen) {
    TempSumYYSq = 0;
    for (ii = 0; ii < Nlen; ii++) {
      TempSumYYSq += CDO->YY[ii] * CDO->YY[ii] * iiWeights[ii];
    }
    SumRDSquared =   TempSumYYSq;
    if (CDO->XTXFlag == 2) {
      for (ii = 0; ii < CDO->OnKappaS; ii++) {
        SumRDSquared -= CDO->OnBeta[CDO->kXFinder[ii]] * 
          (CDO->XTResid[CDO->kXFinder[ii]] + 
           CDO->XTY[CDO->kXFinder[ii]]);
      }
    } else {
      SumRDSquared -= F77_CALL(ddot)(&kLen, CDO->OnBeta, &One, CDO->XTResid, &One);
      SumRDSquared -= F77_CALL(ddot)(&kLen, CDO->OnBeta, &One, CDO->XTY, &One);
    }
    SumWeightsii = 0; 
    for (ii = 0; ii < Nlen; ii++) {
      SumWeightsii += iiWeights[ii];
    }
     
     /*
      for (ii1 = 0; ii1 < Nlen; ii1++) {
        OnYResid = CDO->YY[ii1];
        if (CDO->kXFinder != NULL) {
          for (jj1 = 0; jj1 < CDO->OnKappaS; jj1++) {
             OnYResid -=  CDO->XX[ Nlen * CDO->kXFinder[jj1] + ii1] *
                          CDO->OnBeta[CDO->kXFinder[jj1]];
          
          }
        } else {
          MarkMe = ii1;
          for (jj1 = 0; jj1 < kLen; jj1++) {
             if (CDO->OnBeta[jj1] != 0) {
                OnYResid -= CDO->XX[ MarkMe] *CDO->OnBeta[jj1];
             }
             MarkMe+=Nlen;
          }
        }
        OtherSumRDS += iiWeights[ii1] * OnYResid * OnYResid;
      }
     
      if (abs(OtherSumRDS - SumRDSquared) > .1) {
        Rprintf("Type 1: SumRDSSquared = %f, OtherSumRDS = %f Error in UpdateSq \n",
           SumRDSquared, OtherSumRDS);
        Rprintf("tt1 = %d, tt2 = %d, Nlen = %d, kLen = %d", tt1,tt2, Nlen, kLen);
        R_FlushConsole();
        Rprintf("\n sigmaNoiseSq was %4f, now will be %.4f\n",
             this->sigmaNoiseSq, (OtherSumRDS + priormean * prioralpha) /
                (SumWeightsii + prioralpha));   
        Rprintf("         SumWeightsii = %.4f, prioralpha = %.4f\n",
                  SumWeightsii, prioralpha);     
      }
      */
        //Rprintf("\n sigmaNoiseSq was %4f, now will be %.4f\n",
        //     this->sigmaNoiseSq, (OtherSumRDS + priormean * prioralpha) /
        //        (SumWeightsii + prioralpha));   
        //Rprintf("         SumWeightsii = %.4f, prioralpha = %.4f\n",
        //          SumWeightsii, prioralpha);
        //R_FlushConsole();      
     this->sigmaNoiseSq = (OtherSumRDS + priormean * prioralpha) / (SumWeightsii + prioralpha);  
     return(1);
  } else if (this->TDFNu > 0 && this->CDO != NULL && kLen <= Nlen) {
     TempSumYYSq = 0;
     if (CDO->rWXVec != NULL)  {
       F77_CALL(dsbmv)("U", &Nlen, &Zero, &OneD, iiWeights, &One,
             CDO->YY, &One, &ZeroD, CDO->rWXVec, &One);
       TempSumYYSq = F77_CALL(ddot)(&Nlen, CDO->rWXVec, &One, CDO->YY, &One);
      } else {
       for (ii = 0; ii < Nlen; ii++) {
         TempSumYYSq += 
          CDO->YY[ii] * CDO->YY[ii] * iiWeights[ii];
       }
     }
     SumRDSquared = TempSumYYSq;
     if (CDO->XTXFlag == 2) {
          for (ii = 0; ii < CDO->OnKappaS; ii++) {
            SumRDSquared -= CDO->OnBeta[CDO->kXFinder[ii]] * 
                (CDO->XTResid[CDO->kXFinder[ii]] + 
                 CDO->XTY[CDO->kXFinder[ii]]);
          }
     } else {
        SumRDSquared -= F77_CALL(ddot)(&Nlen, CDO->XTResid, &One, CDO->OnBeta, & One);
        SumRDSquared -= F77_CALL(ddot)(&Nlen, CDO->XTY, &One, CDO->OnBeta, &One);
     }
     SumWeightsii = 0; for (ii = 0; ii < Nlen; ii++) {
       SumWeightsii += iiWeights[ii];
     }
     /*
     OtherSumRDS = 0;
     for (ii1 = 0; ii1 < Nlen; ii1++) {
        OnYResid = CDO->YY[ii1];
        if (CDO->kXFinder != NULL) {
          for (jj1 = 0; jj1 < CDO->OnKappaS; jj1++) {
             OnYResid -=  CDO->XX[ Nlen * CDO->kXFinder[jj1] + ii1] *
                          CDO->OnBeta[CDO->kXFinder[jj1]];
          
          }
        } else {
          MarkMe = ii1;
          for (jj1 = 0; jj1 < kLen; jj1++) {
             if (CDO->OnBeta[jj1] != 0) {
                OnYResid -= CDO->XX[ MarkMe] *CDO->OnBeta[jj1];
             }
             MarkMe+=Nlen;
          }
        }
        OtherSumRDS += iiWeights[ii1] * OnYResid * OnYResid;
      }
      if (abs(OtherSumRDS - SumRDSquared) > .1) {
        Rprintf("Type 2: SumRDSSquared = %f, OtherSumRDS = %f Error in UpdateSq \n",
           SumRDSquared, OtherSumRDS);
        Rprintf("tt1 = %d, tt2 = %d, Nlen = %d, kLen = %d", tt1,tt2, Nlen, kLen);
        R_FlushConsole();
        Rprintf("\n sigmaNoiseSq was %4f, now will be %.4f\n",
             this->sigmaNoiseSq, (OtherSumRDS + priormean * prioralpha) /
                (SumWeightsii + prioralpha));
        
      }  */
           this->sigmaNoiseSq = (OtherSumRDS + priormean * prioralpha) / 
                     (SumWeightsii + prioralpha);
  
  } else if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL) {
    this->sigmaNoiseSq = (this->OnLARS->SumCurResids + priormean * prioralpha) / 
                     (this->Nlen	 + prioralpha);
  } else if (ChoiceCDOOnLARSFlag == 2 && this->CDO != NULL) {
	  if (IPLLF > 3) {
	  Rprintf((char*) "Trying to updateSigmaSq, with SumYYSq = %.4f\n", (double) SumYYSq);
	  R_FlushConsole(); R_ProcessEvents();
      }
     TotalResids = SumYYSq;
     if (CDO->XTXFlag == 1) {
	  if (IPLLF > 3) {
	  Rprintf((char*) "Trying to updateSigmaSq, XTXFlag == 1\n", (double) SumYYSq);
	  R_FlushConsole(); R_ProcessEvents();	  
      }   
       for (ii = 0; ii < kLen; ii++) {
          if (CDO->OnBeta != NULL && CDO->OnBeta[ii] != -999 && CDO->OnBeta[ii] != 0.0) {
               TotalResids -= CDO->OnBeta[ii] * (CDO->XTResid[ii] + CDO->XTY[ii]);
           }
       }
     } else  if (CDO->XTXFlag == 2) {
	  if (IPLLF > 3) {
	  Rprintf((char*) "Trying to updateSigmaSq, XTXFlag == 2\n", (double) SumYYSq);
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
      Rprintf("ConvergentLARS:: EstimateSigmaSq, CDO Error total resids is less than zero? = %.4f\n",
                   TotalResids);
      Rprintf("\nConvergentLARS:: Current CDO->OnBeta\n");
       PrintVector(CDO->OnBeta, kLen);
       //
      Rprintf("\nConvergentLARS:: Current CDO->OnGamma = %.4f, CDO->OnGammas\n", CDO->OnGamma);
       Rprintf("\n OnGammas = "); PrintVector(CDO->OnGammas, kLen);
      Rprintf("\nConvergentLARS:: Current OnBetas\n");
        Rprintf("\n OnBetas = "); PrintVector(OnBetas, kLen); R_FlushConsole();
      Rprintf("\nConvergentLARS:: CurrentXTResid\n");
        Rprintf("\n XTResid = "); PrintVector(CDO->XTResid, kLen);
      R_FlushConsole(); R_ProcessEvents();
        SetupSumYYSq();
        if (RealXXS != NULL && yys != NULL && CDO->XTY != NULL) {
              tMatTimesVec(kLen, Nlen,  RealXXS, yys, CDO->XTY); 
        }
      Rprintf("\nMaking XTResid for Residuals \n");
        CDO->MakeXTResid();
        Rprintf("ConvergentLARS:: We're going to try this one more time\n"); 
        R_FlushConsole();
        TotalResids = SumYYSq;
     if (CDO->XTXFlag == 1) {
  	  if (IPLLF > 3) {
  	  Rprintf((char*) "Trying to updateSigmaSq, XTXFlag == 1\n", (double) SumYYSq);
  	  R_FlushConsole(); R_ProcessEvents();	  
        }   
         for (ii = 0; ii < kLen; ii++) {
            if (CDO->OnBeta != NULL && CDO->OnBeta[ii] != -999 && CDO->OnBeta[ii] != 0.0) {
                 TotalResids -= CDO->OnBeta[ii] * (CDO->XTResid[ii] + CDO->XTY[ii]);
             }
         }
     } else  if (CDO->XTXFlag == 2) {
	  if (IPLLF > 3) {
	  Rprintf((char*) "Trying to updateSigmaSq, XTXFlag == 2\n", (double) SumYYSq);
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
		      Rprintf("ConvergentLARS:: EstimateSigmaSq, CDO Error total resids is less than zero? = %.4f\n",
		                   TotalResids);
		      R_FlushConsole(); R_ProcessEvents();      
		      Rprintf((char*) "SumYYSq = %.5f, Printing XTY\n", SumYYSq); 
		      PrintVector(CDO->XTY, kLen);
		      //PrintVector(CDO->XTY, kLen);
		      Rprintf((char*) "Printing XTResid\n"); 
		      PrintVector(CDO->XTResid, kLen); R_FlushConsole();       
		     
		      Rprintf((char*) "OnBeta: \n");
		      PrintVector(CDO->OnBeta, kLen);
		      Rprintf((char*) "OnKappaS = %d, the kXFinder: \n", CDO->OnKappaS);
		      Rprintf((char*) "TotalResid Pretend Calc = %.4f\n", TotalResids);
		      PrintVector(CDO->kXFinder, CDO->OnKappaS);
		        int jj; double Resid;
		        TotalResids = 0;
		      if (RealXXS != NULL) {
		      for (ii = 0; ii < Nlen; ii++) {
			      Resid = yys[ii];
			      for (jj = 0; jj < CDO->OnKappaS; jj++) {
				      Resid -= RealXXS[ CDO->kXFinder[jj] * Nlen + ii] * 
				                   CDO->OnBeta[CDO->kXFinder[jj]];
			      }
			      TotalResids += Resid*Resid;
		      }
	          } else{
		          Rprintf("Needed RealXXS to not be NULL, failed \n");
		          R_FlushConsole(); R_ProcessEvents();
		          SuccessFlag = -1; return(-1);
	          }
		      Rprintf((char*) "Calculated it Manually and it is %.4f\n", TotalResids);
		      R_FlushConsole(); R_ProcessEvents();
		       SuccessFlag = -1; return(-1);
	     }
     }
     this->sigmaNoiseSq = (TotalResids + priormean * prioralpha) / (this->Nlen + prioralpha);
     if (this->sigmaNoiseSq < 0 && priormean > 0) { this->sigmaNoiseSq = priormean; }
	  if (IPLLF > 3) {
		  Rprintf((char*) "Finished XTX Calculating SigmaNoiseSq\n", (double) SumYYSq);
		  R_FlushConsole(); R_ProcessEvents();	     
      }     
  } else {
    Rprintf("ConvergentLARS:: Estimate SigmaSq, you've got to be kidding!\n");
    R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
  }
  return(1);	
}	

////////////////////////////////////////////////////////////////////////////////
//   Baseic 2 Norm of yys if necessary
//
//
int ConvergentLARS::SetupSumYYSq() {
   SumYYSq = 0.0;
   int ii; 
   for (ii = 0; ii < Nlen; ii++) {
     SumYYSq += yys[ii] * yys[ii];
   }
   return(1);
}

//////////////////////////////////////////////////////////////////////////////
//   Here we copy sequences rLambdaDK and rLambdaAK to their places.
//
//
int ConvergentLARS::SetupLambdaKSeq( double *rLambdaDK, double *rLambdaAK) {
	   FFree(&LambdaDK); FFree(&LambdaAK);
		 this->LambdaDK = (double*) Calloc(TotalRuns+1, double);
		 this->LambdaAK = (double*) Calloc(TotalRuns+1, double);
		 if (this->LambdaAK == NULL) {
			  Rprintf("SetupLambdaKSeq, no room for LambdaKSeq\n");
			  if (this->LambdaDK != NULL) {
			    Free(this->LambdaDK); this->LambdaDK = NULL;
		      }
		      return(-1);
		}
	   //int ii;
	   int  One = 1;
	   F77_CALL(dcopy)(&TotalRuns, rLambdaDK, &One,
		      this->LambdaDK, &One);
	   F77_CALL(dcopy)(&TotalRuns, rLambdaAK, &One,
		      this->LambdaAK, &One);    
	   // for (ii = 0; ii < this->TotalRuns; ii++) {
		 //  this->LambdaDK[ii] = rLambdaDK[ii];
		 //  this->LambdaAK[ii] = rLambdaAK[ii];
     //  }
       if (this->LambdaSeqF == 1) {
	       this->LambdaSeqF = 2;
       } else if (this->LambdaSeqF == 3) {
	       this->LambdaSeqF = 4;
       }       
       //this->LambdaSeqF = 1;
	   return(1);
}
int ConvergentLARS::SetupOrderSeq( int *rOrderSeq) {
	    FFree(&OrderSeq);
		this->OrderSeq = (int*) Calloc(TotalRuns+1, int);
		 if (this->OrderSeq == NULL) {
			  Rprintf("SetupLambdaKSeq, no room for LambdaKSeq\n");
		      return(-1);
		}
	   int ii;
	   for (ii = 0; ii < this->TotalRuns; ii++) {
		   this->OrderSeq[ii] = rOrderSeq[ii];
       }
       if (this->LambdaSeqF == 1) {
	       this->LambdaSeqF = 3;
       } else if (this->LambdaSeqF == 2) {
	       this->LambdaSeqF = 4;
       }
	   return(1);
}

///////////////////////////////////////////////////////////////////////////////
//   InitFixKa:
//       Useful for ``Fermi-Dirac'' Approximation Two-Lasso
//        This sets up PiRecVec to record changes to Pi necessary
//        to define a fixed KActive set as  algorithm progresses.                                                       
//
int ConvergentLARS::InitFixKa(int rFixKa) {
	     FixKa = rFixKa;
	     ppiuse = ((double)FixKa )/((double) kLen);
	     FFree(&PiRecVec);
         this->PiRecVec = NULL;
	     this->PiRecVec = (double *)Calloc( TotalRuns+2, double);
	     if (this->PiRecVec == NULL) {
		     Rprintf((char*)"InitFixKa:  Cannot Apportion PiRecVec \n");
		     R_FlushConsole();
		     return(-1);
	     }		     	
	     return(1);	  
}


////////////////////////////////////////////////////////////////////////////////
//      AofX:
//          Used in ``Fermi-Dirac'' Approximation Two-Lasso
//          Sums exponentiated Logit probabilities, 
//           necessary for robust computation
double ConvergentLARS::AofX(double Onmu) {
	double logLambdasRat = log( OnlambdaA / OnlambdaD);
	double LambdaDiff = OnlambdaD-OnlambdaA;
	double RetTotal = 0;
	int ii; 
	for (ii = 0; ii < kLen; ii++) {
      if ( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu < - 35) {
		RetTotal += 1;
      } else if ( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu > 35)  {
	    RetTotal += exp(- log( 1.0 + 
	          exp( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu)) );
      }  else {
        RetTotal += exp(- log(1.0 + 
              exp( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu))
                       );		
      }
   }
   return(RetTotal);
}

////////////////////////////////////////////////////////////////////////////////
//      InsertBBOn1ofX:
//          Used in ``Fermi-Dirac'' Approximation Two-Lasso
//          Defines new BBOn1s after declaration of chemical fugacity Onmu. 
int ConvergentLARS::InsertBBOn1ofX(double Onmu) {
	double logLambdasRat = log( OnlambdaA / OnlambdaD);
	double LambdaDiff = OnlambdaD-OnlambdaA;
	int ii; 
	for (ii = 0; ii < kLen; ii++) {
	  if ( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu < - 35) {
		   BBOn1[ii] = 1;
      } else if ( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu > 35)  {
	       BBOn1[ii] = exp(- log( 1 + 
	          exp( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu)) );
      }  else {
         BBOn1[ii] = exp( - log( 1 +
                   exp( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu))
                        );		
      }
   }
   return(1);
}

////////////////////////////////////////////////////////////////////////////////
//      DAofX:
//          Used in ``Fermi-Dirac'' Approximation Two-Lasso
//          Calculates local derivative of AofX as function of fugacity Onmu
double ConvergentLARS::DAofX(double Onmu) {
	double logLambdasRat = log( OnlambdaA / OnlambdaD);
	double LambdaDiff = OnlambdaD-OnlambdaA;
	double RetTotal = 0;
	double LDenom = 0;
	int ii; 
	for (ii = 0; ii < kLen; ii++) {
      if ( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu < - 35) {
		   RetTotal += exp(-logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu);
      } else if ( -logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu > 35)  {
         LDenom = - (log( 1 + exp( -logLambdasRat - LambdaDiff * fabs((double) OnBetas[ii]) - Onmu)));	
         RetTotal += exp(LDenom);	
      } else {
         LDenom = - (log( 1 + exp( -logLambdasRat - LambdaDiff * fabs((double) OnBetas[ii]) - Onmu)));	
         RetTotal += exp(-logLambdasRat - LambdaDiff * fabs((double)OnBetas[ii]) - Onmu
                    + LDenom + LDenom);		
      }
   }
   return(RetTotal);
}

////////////////////////////////////////////////////////////////////////////////
//      SolveFormu:
//          Used in ``Fermi-Dirac'' Approximation Two-Lasso
//          Solves for best fugacity Onmu that fixes Sum OnBBOn1 = kActive
double ConvergentLARS::SolveFormu( double GoalTotal, 
         double SuffEpsilon, int MaxIters, double StartOnmu) {
	double Onmu = StartOnmu;
	//double Prevmu = StartOnmu;
	double OnTotal = AofX(Onmu);
	//double PrevTotal = OnTotal;
	
	double OnDerA;
	double ChangeMe;
	int ii; 
	  if (IPLLF > 5) {
		  Rprintf("SolveFormu: Starting with StartOnmu = %.4e, OnTotal = %.4e, GoalTotal Will be %.4e\n",
		      StartOnmu, OnTotal, GoalTotal);
		  R_FlushConsole(); R_ProcessEvents();
      }
	for (ii = 0; ii < MaxIters; ii++) {
		 if (fabs(OnTotal - GoalTotal) < SuffEpsilon)  {
			 return(Onmu);
         }
		 OnDerA = DAofX(Onmu);
		 ChangeMe = - (OnTotal-GoalTotal) / OnDerA;
		 
	  if (IPLLF > 5) {
		  Rprintf("SolveFormu: ii = %d, Onmu = %.4e, OnTotal = %.4e, OnDerA = %.4e, ChangeMe = %.4e\n",
		      (int) ii, (double) Onmu, (double) OnTotal, (double) OnDerA, (double) ChangeMe);
		  R_FlushConsole(); R_ProcessEvents();
      }
		 if ( ChangeMe < 0 && OnTotal < GoalTotal) {
      if (IPLLF > 5) {
		  Rprintf("SolveFormu: ii = %d, ChangeMe = %.4e < 0, but OnTotal = %.4e < GoalTotal %.4e, will do %.4e\n",
		      (int) ii, (double) ChangeMe , (double) OnTotal, (double) GoalTotal, 
		      -.1 *  (OnTotal-GoalTotal) / fabs(OnDerA));
		  R_FlushConsole(); R_ProcessEvents();
      }
			 ChangeMe = -.1 *  (OnTotal-GoalTotal) / fabs(OnDerA);
         } else if (ChangeMe > 0 && (OnTotal > GoalTotal)) {
      if (IPLLF > 5) {
		  Rprintf("SolveFormu: ii = %d, ChangeMe = %.4e > 0, but OnTotal = %.4e > GoalTotal %.4e, will do %.4e\n",
		      (int) ii, (double) ChangeMe , (double) OnTotal, (double) GoalTotal, 
		      -.1 *  (OnTotal-GoalTotal) / fabs(OnDerA));
		  R_FlushConsole(); R_ProcessEvents();
      }	         
	         ChangeMe = -.1 *  (OnTotal-GoalTotal) / fabs(OnDerA);
         }
         //Prevmu = Onmu;
		 Onmu = Onmu + ChangeMe;
		 //PrevTotal = OnTotal;
		 OnTotal = AofX(Onmu);		
		 if (!R_FINITE(Onmu) || ISNAN(Onmu)) {
			 return(Onmu);
         } 
    }
	return(Onmu);
}

int ConvergentLARS::FixKCalculateBBOn() {
  if (FixKa > kLen || FixKa < 0) {
	  Rprintf("FixKCalculateBBOn:  Cannot make this specious Kcount = %d, kLen=%d \n", 
      FixKa, kLen);
  }
  int ii;
	int SuccFlag = 0;
	double Solvedmu = 0;    
  double CurrentCount = 0; 
	double LogitOnppiuse;
	if (ppiuse > 0.0 && ppiuse < 1.0) {
	  LogitOnppiuse = log(ppiuse / (1-ppiuse));
  } else if (ppiuse == 1.0) {
	   LogitOnppiuse = 20;
  } else if (ppiuse == 0.0) {
	  LogitOnppiuse = -20;
  } else {
	  LogitOnppiuse = 0;
  }
	double OnCurrentppiuseCount = AofX(LogitOnppiuse);
	double ApparentMax = 0;
	double ApparentMin  = 0;
	double CountAtMax = 0;
	double CountAtMin = 0;
	int BreakItYouBoughtIt = 0;
 
	if (IPLLF > 4) {
		Rprintf(" \n \n FixKCalculateBBOn(): tt1 = %d, tt2 = %d, OnlambdaA = %.4e, OnlambdaD = %.4e \n",
			tt1, tt2, (double) OnlambdaA, (double) OnlambdaD);
		  Rprintf(" OnGammas  = "); PrintVector(CDO->OnGammas, kLen);
		  Rprintf("\n OnBetas = "); PrintVector(OnBetas, kLen);
		  Rprintf("\n    Now to Run Calculations \n");
		  R_FlushConsole(); R_ProcessEvents();
      }
	  if (OnlambdaA == OnlambdaD) {
		  ppiuse = ((double)FixKa / ((double)kLen));
		  if (ppiuse <= 0 || ppiuse >= 1) {
			  Rprintf("FixKCalculateBBOn(), we have a ppiuse = %.4e Error, quitting \n", ppiuse);
			  R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return(-1);
          } 
		  LogitOnppiuse = log(ppiuse / (1-ppiuse));
		  Solvedmu = LogitOnppiuse;
		  CurrentCount = AofX(LogitOnppiuse);
		  if (IPLLF > 4) {
			  Rprintf("FixKCalculateBBOn(): tt1 = %d, tt2 = %d, OnlambdaA = %.4e, OnlambdaD = %.4e \n",
			      tt1, tt2, (double) OnlambdaA, (double) OnlambdaD);
			  Rprintf("   Started at Default, ppiuse = %.4e, LogitOnppiuse = %.4e, CurrentCount = %.4e \n",
			      (double) ppiuse, (double) LogitOnppiuse, (double) CurrentCount);
			  R_FlushConsole();
          }
      }  else {       
		  if (IPLLF > 4) {
			    Rprintf("FixKCalculateBBOn():  Starting to Solve: ppiuse = %.4e, OnCurrentppiuseCount = %.4e \n",
			       (double) ppiuse, (double) OnCurrentppiuseCount);
			    R_FlushConsole(); R_ProcessEvents();
		  }
		    Solvedmu = SolveFormu(FixKa, SUFFEPSILON_FD, MAXITERS_FD, LogitOnppiuse);   
	        CurrentCount = AofX(Solvedmu);
		  if (IPLLF > 4) {
			    Rprintf("FixKCalculateBBOn(): Initial Solvedmu = %.4e, CurrentCount = %.4e \n",
			       (double) Solvedmu, (double) CurrentCount);
			    R_FlushConsole(); R_ProcessEvents();
		  }	        
		       if (  ISNAN(Solvedmu) ||  !R_FINITE(Solvedmu) ||
		            (fabs( ((double) FixKa )- CurrentCount ) > SUFFEPSILON_FD * 10 )  ) {			            
	  ApparentMax = 60;
	  ApparentMin  = -60;
	  CountAtMax = AofX(ApparentMax);
	  CountAtMin = AofX(ApparentMin);
	  if (OnCurrentppiuseCount < CountAtMax && OnCurrentppiuseCount > (double) FixKa) {
		  ApparentMax = LogitOnppiuse;
		  CountAtMax = OnCurrentppiuseCount;
      }  
	  if (OnCurrentppiuseCount > CountAtMin && OnCurrentppiuseCount < (double) FixKa) {
		  ApparentMax = LogitOnppiuse;
		  CountAtMax = OnCurrentppiuseCount;
      } 
      if (CountAtMax < FixKa) {
        if (IPLLF > 1) {
        Rprintf("Issue, CountatMax = %d <  FixKa= %d \n", CountAtMax, FixKa);
	      Rprintf("FindBBOn: IPLLF = %f, Apparently FixKa is ", IPLLF);
        Rprintf("already less than Max, probabilities to big to tweak \n");
	      R_FlushConsole(); R_FlushConsole();
	      }
	      LogitOnppiuse = ApparentMax; OnCurrentppiuseCount = CountAtMax;
      } else if (CountAtMin > FixKa) {
        if (IPLLF > 1) {
        Rprintf("Issue, CountAtMin=%d >  FixKa = %d \n", CountAtMin, FixKa);
	      Rprintf("FindBBOn: Apparently FixKa is already less than Max, probabilities to big to tweak \n");
	      R_FlushConsole(); R_FlushConsole();
	      }
	      LogitOnppiuse = ApparentMin; OnCurrentppiuseCount = CountAtMin;	      
      } else {	
	       BreakItYouBoughtIt = 0;	
	               OnCurrentppiuseCount = AofX(LogitOnppiuse);            
			       for (ii = 0; ii < MAXITERS_FD; ii++) {			       
				       if ( fabs( ((double) FixKa )- OnCurrentppiuseCount ) < SUFFEPSILON_FD * 10) {
					       Solvedmu = LogitOnppiuse; CurrentCount = OnCurrentppiuseCount;
					       BreakItYouBoughtIt = 1;
					       break;    
				       }
				       if ( ISNAN(Solvedmu) || !R_FINITE(Solvedmu) ||
				            fabs(((double) FixKa )- CurrentCount ) >
				            fabs(((double) FixKa) - OnCurrentppiuseCount) ) {
					       if (FixKa > OnCurrentppiuseCount) {
						       LogitOnppiuse = (ApparentMax + LogitOnppiuse) / 2.0 ;
						       OnCurrentppiuseCount = AofX(LogitOnppiuse);
						       if (OnCurrentppiuseCount < FixKa) {
							       ApparentMin = LogitOnppiuse;
							       CountAtMin = OnCurrentppiuseCount;
						       } else if (OnCurrentppiuseCount > FixKa) {
							       ApparentMax = LogitOnppiuse;
							       CountAtMax = OnCurrentppiuseCount;							       
						       }
						       Solvedmu = LogitOnppiuse	; 			       
					       } else {
						       LogitOnppiuse = (ApparentMin + LogitOnppiuse) / 2.0 ;
						       OnCurrentppiuseCount = AofX(LogitOnppiuse);
						       if (OnCurrentppiuseCount < FixKa) {
							       ApparentMin = LogitOnppiuse;
							       CountAtMin = OnCurrentppiuseCount;
						       } else if (OnCurrentppiuseCount > FixKa) {
							       ApparentMax = LogitOnppiuse;
							       CountAtMax = OnCurrentppiuseCount;							       
						       }						       
						       Solvedmu = LogitOnppiuse	;
					       }         
					       Solvedmu = SolveFormu(FixKa,SUFFEPSILON_FD , MAXITERS_FD, LogitOnppiuse);   
				       } else {
                               if (CurrentCount < FixKa) {
							       ApparentMin = Solvedmu;
							       CountAtMin = CurrentCount;
						       } else if (CurrentCount > FixKa) {
							       ApparentMax = Solvedmu;
							       CountAtMax = CurrentCount;							       
						       }
					       LogitOnppiuse = Solvedmu; OnCurrentppiuseCount = CurrentCount;
					       Solvedmu = SolveFormu(FixKa,SUFFEPSILON_FD , MAXITERS_FD, LogitOnppiuse);   
			           }
			       	   CurrentCount = AofX(Solvedmu);
			       	   if ( !ISNAN(Solvedmu)  && R_FINITE(Solvedmu) && 
			       	        fabs( ((double) FixKa )- CurrentCount ) < SUFFEPSILON_FD * 10) {
				       	   LogitOnppiuse = Solvedmu;
				       	   OnCurrentppiuseCount = CurrentCount;
				       	   BreakItYouBoughtIt = 2;
				       	   break;
			           }
			           if (LogitOnppiuse > 50) {
				           Solvedmu = LogitOnppiuse;
				           CurrentCount = OnCurrentppiuseCount;
				           BreakItYouBoughtIt = 3;
				           break;   
			           } else if (LogitOnppiuse < -50) {
				           Solvedmu = LogitOnppiuse;
				           CurrentCount = OnCurrentppiuseCount;
				           BreakItYouBoughtIt = 4;
				           break;   				           			           
			           }
		           }
		           if ( ISNAN(LogitOnppiuse) || !R_FINITE(LogitOnppiuse) )  {
		             Rprintf("###############################################\n");
		             R_FlushConsole();
			           Rprintf("#### FixKCalculateBBOn ISNAN(LogitOnppiuse) failed, ");
                 Rprintf("#### we've got FixKa = %d, kLen = %d, ", FixKa, kLen);
                 Rprintf("#### LogitOnppiuse = %.4e, piLogitOnppiuse = %.4e \n",
			              (double) LogitOnppiuse,
                    (double) exp(LogitOnppiuse) / ( 1 + exp(LogitOnppiuse)) );
			           Rprintf(" ii = %d \n", ii);
			           Rprintf(" BreakItYouBoughtIt = %d\n", BreakItYouBoughtIt);			                  
			           R_FlushConsole(); R_ProcessEvents();
			           SuccessFlag = -1; return(-1);
		           }		           
		           if ( fabs( ((double) FixKa )- OnCurrentppiuseCount ) > SUFFEPSILON_FD * 4000
		                  && LogitOnppiuse < 50 && LogitOnppiuse > -50) {
			           Rprintf("####  ERROR on FixKa Failure"); R_FlushConsole();
                 Rprintf("####  FixKCalculateBBOn we failed, we've got ");
                 R_FlushConsole();
                 Rprintf("FixKa = %d, kLen = %d, LogitOnppiuse = %.4e, piLogitOnppiuse = %.4e \n",
			                  FixKa, kLen, (double) LogitOnppiuse,(double) exp(LogitOnppiuse) / ( 1 + exp(LogitOnppiuse)) );
			           Rprintf("#### OnCurrentppiuseCount = %.4e, CurrentCount = %.4e \n",
			              OnCurrentppiuseCount, CurrentCount); R_FlushConsole();
			           Rprintf("#### Solvedmu = %.4e, piSolvedmu = %.4e \n", 
			               Solvedmu, exp(Solvedmu) / ( 1.0 + exp(Solvedmu)));
			           Rprintf("#### ApparentMax = %.4e, CountAtMax = %.4e \n", 
			               ApparentMax, CountAtMax);  R_FlushConsole();	
			           Rprintf("#### ApparentMin = %.4e, CountAtMin = %.4e \n", 
			               ApparentMin, CountAtMin);	R_FlushConsole();		               		               
			           Rprintf("#### ii = %d \n", ii);
			           Rprintf("#### BreakItYouBoughtIt = %d\n", BreakItYouBoughtIt);
			           R_FlushConsole(); R_ProcessEvents();
			           SuccessFlag = -1; return(-1);
		           }
	    }	           
			   } else {
				    LogitOnppiuse = Solvedmu;
				    OnCurrentppiuseCount = CurrentCount;
		       }
		     if (ISNAN(LogitOnppiuse) || !R_FINITE(LogitOnppiuse)) {
		       Rprintf("#### ERROR ANOUNCEMENT \n");  R_FlushConsole();
			     Rprintf("#### 2L:FixKCalculateBBOn: LogitOnppiuse ISNAN!!, LogitOnppiuse = %.4e, FixKa = %d \n",  
			       (double)  LogitOnppiuse, FixKa);  R_FlushConsole();
			     Rprintf("#### 2L:FixKCalculateBBOn: Setting SuccessFlag = -1\n");
			     R_FlushConsole();
			     SuccessFlag =   -1; return(-1);
		     }		       
		     if (LogitOnppiuse > 50) {
			     ppiuse = 1.0;
		     } else if (LogitOnppiuse < -50) {
			     ppiuse = exp(LogitOnppiuse);
		     } else { 
			     ppiuse = exp(LogitOnppiuse) / ( 1.0 + exp(LogitOnppiuse));
		     }

    }
	if (IPLLF > 3) {
	      Rprintf("#### 2L:FixKCalculateBBOn ");
        Rprintf("tt1 = %d, tt2 = %d, FixKCalculateBBOn settles",
          tt1, tt2);
        Rprintf(" upon Solvedmu = %.4e, LogitOnppiuse = %.4e, ppiuse = %.4e, CurrentCount = %.4e\n",
	        (double) Solvedmu, (double) LogitOnppiuse, (double) ppiuse, (double) CurrentCount);
	      R_FlushConsole(); R_ProcessEvents();
    }
    if (PiRecVec != NULL)  {
      PiRecVec[tt1] = ppiuse;
    }
	SuccFlag = InsertBBOn1ofX(LogitOnppiuse);    
	return(SuccFlag);
}



int ConvergentLARS::SetupTD2Lasso(double rTDFNu, double *riiWeights) {
              TDFNu = rTDFNu;
   iiWeights = riiWeights;
   if (CDO->rWXVec != NULL) {
     CDO->RefreshWeights(riiWeights);
   }   else {
     CDO->SetupWeights(riiWeights);  
   }
   return(1);            
}

int ConvergentLARS::RefreshTD2Lasso() {
   int ii;   double ResidOnii;  int One = 1;   double InvSigmaSq = 1.0/ sigmaNoiseSq;
   int jj;
   if (CDO->YY == NULL || CDO->XX == NULL) {
     Rprintf("RefreshTD2Lasso: I can't do that I have NULL values\n");
   }
   if (CDO->XTXFlag == 2 && CDO->OnKappaS > 0 && CDO->OnKappaS < Nlen /2  && CDO->OnKappaS < kLen / 2) {
     F77_CALL(dcopy)(&Nlen, CDO->YY, &One, iiWeights, &One);
     for (jj = 0; jj < CDO->OnKappaS; jj++) {
       if (CDO->OnBeta[CDO->kXFinder[jj]] != 0.0) {
         F77_CALL(daxpy)(&Nlen, &CDO->OnBeta[CDO->kXFinder[jj]], 
           CDO->XX + Nlen * CDO->kXFinder[jj], &One, 
           iiWeights, &One);
       }
     }
     for (ii = 0; ii < Nlen; ii++) {
        iiWeights[ii] = (TDFNu + 1) / (TDFNu + iiWeights[ii] * iiWeights[ii] * InvSigmaSq);
     } 
   } else {
     for (ii = 0; ii < Nlen; ii++) {
        //Rprintf("RTD2Lasso: Updating Coordinate[ii=%d] \n", ii); R_FlushConsole();
        ResidOnii = CDO->YY[ii] - F77_CALL(ddot)(&kLen, CDO->XX + ii, &Nlen, CDO->OnBeta, &One);
        //Rprintf("RTD2Lasso: Updated Coordinate[ii=%d]  to ResidOnii = %.4f\n", ii, ResidOnii); R_FlushConsole();
        iiWeights[ii] = (TDFNu +1)  / ( TDFNu + ResidOnii * ResidOnii * InvSigmaSq); 
        //Rprintf("RTD2Lasso: Updated Coordinate[ii=%d]  to iiWeights[%d] = %f\n",
        //    ii, ii, iiWeights[ii]); R_FlushConsole();      
     }
   }
   if (IPLLF > 3) {
     Rprintf("RefreshTD2Lasso: Got to End of updating all iiWeights\n");
   }  
   //Rprintf("RefreshTD2Lasso: Going to Update Weights now on CDO I Promise!!!\n"); 
   SuccessFlag = CDO->RefreshWeights(iiWeights);
   if (SuccessFlag < 0) {
     Rprintf("RefreshTD2Lasso: Refresh Weights Failed: Abort\n");
     R_FlushConsole();
   }
   return(SuccessFlag);
}
  
   


int ConvergentLARS::RunConvergentAlgorithmInitializeStuff(
   int NumCDOConv, int NumTotalRuns) {
 //int cnt1 = 0;
 //this->IPLLF = 0;
 int ErrorLarsFlag = -1; 
 
 if (ErrorLarsFlag < 0) {ErrorLarsFlag = -1;}//int jj;
 //if (this->OnLARS != NULL) {this->OnLARS->PrintFlag = 4;}
 if (this->IPLLF > 2) {
	   Rprintf((char*) "About to Run ConvergentLARS \n"); 
     R_FlushConsole(); R_ProcessEvents();
	   if (this->OnLARS != NULL) {
		    Rprintf((char*) "this->OnLARS != NULL \n"); 
        R_FlushConsole(); R_ProcessEvents();
	   }
	   if (this->CDO != NULL) {
		    Rprintf((char*) "this->CDO != NULL \n"); 
        R_FlushConsole(); R_ProcessEvents();
		    if (this->CDO->YY != NULL) {
			    Rprintf((char*) "Printing CDO->YY \n"); R_FlushConsole();
		        PrintVector(CDO->YY, Nlen); R_FlushConsole();
	      } 
        if (this->CDO->YY == NULL) {
		        Rprintf((char*) "CDO->YY is NULL\n"); R_FlushConsole();
	      }
	      if (this->CDO->XTY != NULL) {
		        Rprintf((char*) "Printing CDO->XTY \n"); R_FlushConsole();
		        PrintVector(CDO->XTY, kLen); R_FlushConsole();
	      }
	      if (this->CDO->XX != NULL) {
		        Rprintf((char*) "Printing CDO->XX \n"); R_FlushConsole();
		        PrintRMatrix(CDO->XX, Nlen, kLen);
	      } if (this->CDO->YY == NULL) {
		        Rprintf((char*) "CDO->XX is NULL\n"); R_FlushConsole();
	      }
	      if (this->CDO->pXTX != NULL) {
		        Rprintf((char*) "Printing CDO->XTX \n"); R_FlushConsole();
		        PrintRpMatrix(CDO->pXTX, kLen, kLen);
	      }	   
	      R_FlushConsole(); R_ProcessEvents();     
	    }	
  }  
  if (this->IPLLF > 4 && FakeFlag >= 0) {
	    Rprintf((char*) "RunConvergentLARS, RealSet, here's the YY\n");
	    PrintVector(this->yys, Nlen);
	    Rprintf((char*) "\n And here is the XX\n");
	    PrintRMatrix(RealXXS, Nlen, kLen);
        Rprintf((char*) "\n\n");
	    R_FlushConsole();
  }
  if (this->IPLLF > 0) {
	    Rprintf((char*) "RunConvergentLars: NumTotalRuns = %d,", NumTotalRuns);
      Rprintf(" NumEMConv = %d \n",  NumEMConv);
	    R_FlushConsole();
  }
  this->OnlambdaD = this->StlambdaD;
  this->OnlambdaA = this->StlambdaA;
  if (LambdaSeqF == 1 || LambdaSeqF == 3 ) {
	   if (LambdaAK != NULL) {FFree(&LambdaAK); }
	   LambdaAK = (double*) Calloc( NumTotalRuns+2, double);
	   if (LambdaAK == NULL) {Rprintf((char*)"RunConvergentLars: no space for LambdaAK\n");
	     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);}   
	   if (LambdaDK != NULL) {FFree(&LambdaDK); }
	   LambdaDK = (double*) Calloc( NumTotalRuns+2, double);
	   if (LambdaDK == NULL) {Rprintf((char*)"RunConvergentLars: no space for LambdaAK\n");
	     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);}  	     
  }
  if (this->IPLLF > 0) {
	    Rprintf((char*) "LambdaSeqF == %d \n", LambdaSeqF);
	    Rprintf((char*) "StlambdaA = %.4f, MultAConst = %.4f\n", 
	        (double) StlambdaA, (double) MultAConst);
	    Rprintf((char*) "StlambdaD = %.4f, MultDConst = %.4f\n", 
	        (double) StlambdaD, (double) MultDConst);
	    if (LambdaSeqF == 2 || LambdaSeqF == 4) {
		    Rprintf((char*) "Printing LambdaAK \n");
		    PrintVector(LambdaAK, NumTotalRuns);
		    Rprintf((char*) "Printing LambdaDK\n");
		    PrintVector(LambdaDK, NumTotalRuns);   
	    }
	    R_FlushConsole();
   }  
   if (ChoiceCDOOnLARSFlag == 1 && this->OnLARS != NULL && 
       FakeFlag >= 0 && this->kLen < 200 && this->GAll == NULL) {
         this->GenGAll();
        this->UpdateGAll(0);
   }
   if (this->IPLLF > 1) {
	    Rprintf((char*) "RunConvergentLars: Done Initiating GAll \n",
        NumTotalRuns, NumEMConv);
	    R_FlushConsole();
   }
   return(1);

}
int ConvergentLARS::TwoLassoPrintDiagnosticA () {
		 if (CDO != NULL) {
			 Rprintf((char*) "RunConvergentLars:: Just Prior to Run CDO \n");
			 Rprintf((char*) "RunConvergentLars:: tt1 = 0, tt2 = 0; sigmaNoiseSq  = %.8f\n",
			              (double) sigmaNoiseSq);			 
			 Rprintf((char*) "RunConvergentLars:: tt1 = 0, tt2 = 0, CDO->XTX = \n");
			 PrintRpMatrix(CDO->pXTX, kLen, kLen);
			 Rprintf((char*) "RunConvergentLars:: tt1 = 0, tt2 = 0, CDO->XTY = \n");
			 PrintVector(CDO->XTY, kLen);
			 Rprintf((char*) "RunConvergentLars:: tt1 = 0; tt2 = 0, CDO->OnBeta = \n");
			 PrintVector(CDO->OnBeta, kLen);
			 Rprintf((char*) "RunConvergentLars::  CDO->XTResid = \n");
			 PrintVector(CDO->XTResid, kLen);
			 Rprintf((char*) "RunConvergentLars::  CDO->NLen = %d, CDO->OnGamma = %.4f\n",
			               CDO->NLen, CDO->OnGamma);
			 if (CDO->OnGammas != NULL) {
				 Rprintf((char*) "RunConvergentLARS:: OnGammas = \n");
				 PrintVector(CDO->OnGammas, kLen);
       }
      }
   return(1);
}
int ConvergentLARS::RunAlgorithmRunLARS() {
   int ErrorLarsFlag;
	          if (IPLLF >= 4) {
		            Rprintf((char*)"RunConvergentLars:: OnLARS path.\n");
	              R_FlushConsole(); R_ProcessEvents();		        			        
		        } 
			      this->OnLARS->lambda = 2 * sigmaNoiseSq;
			      if (IPLLF > 3) {
				        this->OnLARS->PrintWeights();   
			      }
            	
			      if (IPLLF >= 1) {
			          Rprintf((char*)"About to run LarsINIT\n", 
		                                    (double) OnlambdaD,(double) OnlambdaA);
		            R_FlushConsole();		        
				    }  				        	        
			      this->OnLARS->LARSInit();
			      
            ////////////////////////////////////////////////////////////
            // Runnings Weighted LARS to fit LASSO
            if (IPLLF >= 4) {
			            Rprintf((char*)"About to run LARS Lasso Algorithm\n", 
		                                     (double) OnlambdaD,(double) OnlambdaA);
		              R_FlushConsole();		        
				    }  		        
			      ErrorLarsFlag = this->OnLARS->LarsLassoAlgorithm2();
				    if (ErrorLarsFlag < 0) {
					      Rprintf((char*)"ErrorLarsFlag =%d\n", ErrorLarsFlag);
					      this->OnBetas[0] = -999;
					      R_FlushConsole();
					      return(ErrorLarsFlag);
				    }
            return(1);            
}


int ConvergentLARS::AllLARSCDODiagnositc() {
    int ii;
		          if (this->OnLARS != NULL) {
			            Rprintf((char*)"\n  Now OnLARS->OnNus \n");
				          for (ii = 0; ii < this->kLen; ii++) {
					           Rprintf("%d=%.4f,  ", ii, (double) this->OnLARS->OnBetas[ii]);
				          }
				          R_FlushConsole();
				          Rprintf((char*)"\n");
			              Rprintf((char*) "\n  OnLARS->OnNus = c( ");
					        for (ii = 0; ii < this->kLen-1; ii++) {
						        Rprintf((char*) "%.8f, ", (double)(double) this->OnLARS->OnBetas[ii]);
					        }
					        Rprintf((char*) " %.8f ) \n", this->OnLARS->OnBetas[this->kLen-1]);	          
				            R_FlushConsole();
			        }
			          if (this->CDO != NULL) {
				          Rprintf((char*)"\n  CDO->OnBeta\n");
				          R_FlushConsole(); PrintVector(CDO->OnBeta, kLen);
			          }          
			          Rprintf((char*)"  Now OnBetas \n");
			              R_FlushConsole(); PrintVector(this->OnBetas, kLen);
		             Rprintf((char*)"\n  Now BBOn1 \n");
                          R_FlushConsole(); PrintVector(this->BBOn1, kLen);  
			          R_FlushConsole();	  	          
			          R_FlushConsole();    	                
	  R_FlushConsole();
		return(1);
}

