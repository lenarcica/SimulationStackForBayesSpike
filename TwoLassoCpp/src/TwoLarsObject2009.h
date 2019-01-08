#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <Rdefines.h>
  #include <Rinternals.h>
  #include <stdio.h>
  #define RMATH 0
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


#ifndef CoordinateDescentLassoDD
  #include "CoordinateDescentLasso2014.h"
  #define CoordinateDescentLassoDD 0
#endif

#ifndef MyMatrixOp2009DD
  #include "MyMatrixOp2009.h"
  #define MyMatrixOp2009DD 0
#endif

#define SUFFEPSILON_FD .0000001
#define MAXITERS_FD  140
#define KLENCOPY 200


int LarsConvergencyBITest(double * yys, int NLen, double *xxs, int kLen,   //4
  double *InitBetas,                                   //5
  double *OldBetas, double *OutBetas, double *BetasKeep,  //8
  int NumTotalRuns, int NumEMConv, double MultEMCons,  //11
  double StlambdaD, double StlambdaA,                   //13
  double lambdaDMultC, double lambdaAmultC,           //15
  int lambdaAmultStop, double ppiuse,                 //17
  double sigmaNoiseSq, double *BBHatsKeep,            //19
  double *NusKeep,                                    //20
  double *LambdaDK, double *LambdaAK,                 //22
  double *RecordPostPBj, double DConfidenceInt, double *RecordConfidenceInts,//25
  int LAKSeq, int *OrderSeq, int NumCDOConv, double CDOEpsilon,             //29
  int InitKKs,                                                              //30
  double InverseGammaConstant, int FixKa, double *PiRecVec,                 //33
                      double *iiWeights, double TDFNu, int PrintFlag,  //36
                      double Lambda3, double *Lambda3Seq);           //37,38
                      
int LarsConvergencyBI(double * yys, int NLen, double *xxs, int kLen, 
  double *InitBetas, double *OldBetas, double *OutBetas, double *BetasKeep,
  int NumTotalRuns, int NumEMConv, double MultEMCons,
  double StlambdaD, double StlambdaA, 
  double lambdaDMultC, double lambdaAmultC,
  int lambdaAmultStop, double ppiuse,
  double psigmaNoiseSq, 
  double *BBHatsKeep, 
  double *NusKeep, double *LambdaDK, double *LambdaAK,
  double *RecordPostPBj, double DConfidenceInt, double *RecordConfidenceInts,
  int LAKSeq, int *OrderSeq, int NumCDOConv, double CDOEpsilon,
  int InitKKs, double InverseGammaConstant,
  int FixKa, double *PiRecVec,
  double *iiWeights, double TDFNu, int PrintFlag,
  double *pSigmaVec, double Lambda3, double *Lambda3Seq);

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
  double InverseGammaConstant,
  double *iiWeights, double TDFNu, int PrintFlag,
  double *pSigmaMultVec);

                   
class ConvergentLARS  {
	public: 
	int ChoiceCDOOnLARSFlag;
	LARSObject  *OnLARS;             // Store activity in LARSObject?
  CoordinateDescentObject *CDO;    // Store in a CDO?

	int FakeFlag;
	
	int Nlen; int kLen;              // Nlen: Length of Y, kLen # of covariates 
	
	// Data, in yys, and RealXXS
	//   PassByPointer Flags note whether data is a pointer copy (don't free)
	double *yys;     int yysPassByPointer;
	double *RealXXS; int RealXXSPassByPointer;
	
	//  "GAll"  is t(X) %*% X is sum of squares matrix
	// Depending On what was passed, GAll might be undesired data
	double *GAll;    int GAllPassByPointer;

  //  TotalRuns, NumEMConv
	int TotalRuns; int NumEMConv;
  double *OnBetas;
	double *OnNus;
	double *OnMultiCons;
	double *BBOn1;
	double *BBPrev;
	double *BetasKeep;
	int SuccessFlag;
	double FactorGammaDenominator;

	double *NusKeep;
	double *BBHatsKeep;
	
	double *BBHatsNewStep;	
	double *SigmaRecVec;  int LSigmaRecVec;
	double *PiRecVec;
	
	double m1, m2;
  int OnShr;
  double BVL1, BVL2, sigmaNoiseSq;

  double ppiuse;
  double StlambdaD, StlambdaA;
  double OnlambdaD, OnlambdaA;
    
  double *LambdaDK, *LambdaAK; int LambdaSeqF; int *OrderSeq;
  double MultDConst, MultAConst;

  double *CurrentPostPBj; double *RecordPostPBj;    
  double *CurrentConfidenceInts; double *RecordConfidenceInts;
  
  double *CurGammaA, *CurGammaD;
  
  double *DiagXTX;
  
  //Unworking Confidence Interval Calculations
  int SetCopyCurrentPostPBj(int Run);
  int SetCopyCurrentConfidenceInts(int Run);
  int CalculateCurrentPostPBj();
  int CalculateCurrentConfidenceInts(double DesiredCInt);
 
  double DConfidenceInt;
  int SetupLambdaKSeq( double *rLambdaDK, double *rLambdaAK);
  int SetupOrderSeq( int *rOrderSeq );
     
  int StopAInt;
  
  int Ontt1, Ontt2; int tt1, tt2;
  
  int IPLLF;
  
  int GenGAll();
  int UpdateGAll(int unweight);
  int UpdateWAll(int weightFlag);
  int RunConvergentLARS(int NumTotalRuns, int NumEMConv, double MultEMCons,
    int NumCDOConv, double CDOEpsilon);
  int RunConvergentAlgorithmInitializeStuff(int NumCDOConv, int NumTotalRuns);
  int TwoLassoPrintDiagnosticA();   int AllLARSCDODiagnositc();
  int RunAlgorithmRunLARS();
  
  int RunConvergentPILARS(int NumTotalRuns, int NumEMConv, double MultEMCons, 
    double rm1, double rm2, int NumEConv, double SigmaEta, double SigmaBarSq,
    int NumCDOConv, double CDOEpsilon);  
  int UpdateSigmaSq(double prioralpha, double priormean);  
  int SetDiagXTX();
  int UpdateMaxPi();
  double SumYYSq;
  int SetupSumYYSq();
	int SetupCDO(int SetBeta); int AfterCDO();
	double *SigmaMultVec;

///////////////////////////////////////////////////////////////
// FixKa Versions
//   If FixKa > 0: We are doing Fermi-Dirac version of Lasso
//
	int FixKa;
  int InitFixKa(int rFixKa);
  int InsertBBOn1ofX(double Onmu);
  int FixKCalculateBBOn();
  double AofX(double Onmu);	
  double DAofX(double Onmu);
  double SolveFormu( double GoalTotal, 
         double SuffEpsilon, int MaxIters, double StartOnmu);	
         
  double TDFNu; double *iiWeights; 
  int SetupTD2Lasso(double rTDFNu, double *riiWeights);
  int RefreshTD2Lasso();
  int UpdateT2LWeights();
  
     
  double *Lambda3Seq;
  int SetupLambda3Seq (double *rLambda3Seq) {
    Lambda3Seq = (double*) Calloc(TotalRuns, double);
    if (Lambda3Seq ==NULL) {
      Rprintf("Error: Lambda3Seq could not be specified\n");R_FlushConsole();
      SuccessFlag = -1; return(-1);
    }
    for(int ii = 0; ii < TotalRuns; ii++) {
      Lambda3Seq[ii] = rLambda3Seq[ii];
    }
    return(1);
  }

  ConvergentLARS(){
	  CurrentPostPBj = NULL;
    RecordPostPBj = NULL;    CurrentConfidenceInts = NULL;
    RecordConfidenceInts = NULL; 
    Lambda3Seq = NULL;
    SigmaMultVec = NULL;
    return;}
   
  ConvergentLARS(double * ryys, int rNLen, double *rxxs, int rkLen, 
    double *rInitBetas, double rStlambdaD, double rStlambdaA,
    double rsigmaNoiseSq,
    double rMultDConst, double rMultAConst,
    int rStopAInt, double rppiuse, int rTotalRuns, int InitKKs,
    double InverseGammaConstant, int PrintFlag) {
    CDO = NULL; OnLARS = NULL;  
    ppiuse = rppiuse;
	     
    TotalRuns = rTotalRuns;
	  SuccessFlag = 1;
	  
    LambdaDK = NULL;  RealXXS = NULL; yys = NULL;
	  LambdaAK = NULL; OrderSeq = NULL; OnMultiCons = NULL;
	  OnBetas = NULL; BBOn1 = NULL; OnMultiCons = NULL; BBOn1 = NULL;
	  BBPrev = NULL; BetasKeep = NULL; BBHatsKeep = NULL; BBHatsNewStep = NULL;
	  SigmaRecVec = NULL; PiRecVec = NULL; CurrentPostPBj = NULL;
    RecordPostPBj = NULL;    CurrentConfidenceInts = NULL;
    RecordConfidenceInts = NULL; CurGammaA = NULL; CurGammaD = NULL; DiagXTX=NULL;
    ChoiceCDOOnLARSFlag = 2; LambdaSeqF = 1;  FixKa = -1;  Lambda3Seq = NULL;
    SigmaMultVec = NULL;
           
    iiWeights = NULL; TDFNu = -1;

    NumEMConv = 0;
    SigmaMultVec = NULL;    	                                    

    SuccessFlag = 1; 
   
    ChoiceCDOOnLARSFlag = 1;
    FactorGammaDenominator = InverseGammaConstant; 
	  TotalRuns = rTotalRuns;
	  tt1 = 0; tt2 = 0;
	  
    NusKeep = NULL;
	  OrderSeq = NULL; LambdaSeqF = 1;

    sigmaNoiseSq = rsigmaNoiseSq;
    	     	                    
	  if (PrintFlag >= 0) { IPLLF = PrintFlag;} else {IPLLF = -1;}

	  int ii,jj, cnt1;
	     
	  if (IPLLF > 0) {
		  Rprintf((char*)"Generating ConvergentLARS\n");
		  R_FlushConsole();
	  }
	  int cn1;
	  if (rNLen < 0) { 
		  if (IPLLF >= 0) {
			  Rprintf("ConvergentLARS Load GAll version only \n");
			  R_FlushConsole();
		  }
		  Nlen = - rNLen;
		  FakeFlag = -1;
		  kLen = rkLen;

		  GAll = NULL;
		  GAllPassByPointer = -1;
		  GAll = (double *) Calloc((rkLen+1) * (rkLen) + 1, double);
		  if (GAll == NULL) {
        Rprintf((char*)"LarsCov, Cannot Apportion GAll\b");
		    R_FlushConsole(); SuccessFlag = -1;
		    return;
	    }
	    cn1 = kLen * kLen;
	    for (ii = 0; ii < cn1; ii++) {
		    GAll[ii] = (double) rxxs[ii];
	    }
		  if (IPLLF >= 0) {
			  Rprintf("ConvergentLARS: Loaded GAll\n");
			  R_FlushConsole(); SuccessFlag = -1;
		  }	          
	  }  else {
		  FakeFlag = 1;
		  Nlen = rNLen; kLen = rkLen;
		  yys = NULL; yysPassByPointer = -1;
		  yys = (double *) Calloc( rNLen+2,double);
			if (yys == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion yys \n");  RFP();
				SuccessFlag = -1;
				return;
			}
			RealXXSPassByPointer = -1;
      RealXXS = NULL;		     	      
	    RealXXS = (double *) Calloc( (rNLen+1) * (rkLen+1), double);
		  if (RealXXS == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion RealXXS \n"); RFP();
			  SuccessFlag = -1;
			  return;
		  }	      
	    GAll = NULL;
	    m1 = 0; m2 = 0;
	    for (ii = 0; ii < rNLen; ii++) { yys[ii] = (double) ryys[ii]; }
	      cn1 = rNLen * rkLen;	
	      for (ii = 0; ii < cn1; ii++) { RealXXS[ii] = (double) rxxs[ii]; }
      }
      BetasKeep = NULL;	      
	    BetasKeep = (double *)Calloc((rkLen+1)* TotalRuns, double);
      if (BetasKeep == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion BetasKeep \n"); RFP();
			  SuccessFlag = -1;
			  return;
		  }
      BBHatsKeep = NULL;		     	     
	    BBHatsKeep = (double *)Calloc((rkLen +1)* TotalRuns, double);
      if (BBHatsKeep == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion BBHatsKeep \n"); RFP();
			  SuccessFlag = -1; return;
		  }
		  NusKeep = NULL;
	    NusKeep = (double *)Calloc((rkLen +1)* TotalRuns, double);
      if (NusKeep == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion BBHatsKeep \n"); RFP();
			  SuccessFlag = -1; return;
      }		 
		  OnBetas = NULL;
	    OnBetas = (double *)Calloc(kLen+1, double);
      if (OnBetas == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion OnBetas \n"); RFP();
			  SuccessFlag = -1; return;
		  }
		  if (rInitBetas == NULL) {
			  for (ii = 0; ii < kLen; ii++) {
				  OnBetas[ii] = 0;
	      } 
      } else {
	      for (ii = 0; ii < kLen; ii++) {
		      OnBetas[ii] = rInitBetas[ii];
	      }
      }
	    OnNus= NULL;		     	     
	    OnNus = (double *)Calloc(kLen+1, double);
      if (OnNus == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion OnNus \n");  RFP();
			  SuccessFlag = -1; return;
		  }	     
	    OnMultiCons = (double *)Calloc(kLen+1, double);
	    BBOn1 = NULL;
	    BBOn1 = (double *)Calloc(kLen+1, double);
      if (BBOn1 == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion BBOn1 \n"); RFP();
			  SuccessFlag = -1; return;
		  }	     
	    BBPrev = (double *)Calloc(kLen+1, double);
	    sigmaNoiseSq = (double) rsigmaNoiseSq;
	    StlambdaD = (double) rStlambdaD;
	    StlambdaA = (double) rStlambdaA;
	    OnlambdaD = StlambdaD;
	    OnlambdaA = StlambdaA;
	    MultDConst = (double) rMultDConst;
	    MultAConst = (double) rMultAConst;
		  StopAInt = rStopAInt;
		  ppiuse = (double) rppiuse;
		  LambdaDK = NULL; LambdaAK = NULL; LambdaSeqF = 1;
		  OrderSeq = NULL;
	      	     	     	     

	    OnLARS = new LARSObject(ryys, rNLen, rxxs, rkLen, 
        rInitBetas, (StlambdaD+StlambdaA)/2.0);
	    if (OnLARS == NULL) {
		    Rprintf((char*)"MyConvLars: No room for OnLARS\n");  RFP();
		    SuccessFlag = -1; return;
	    }
	     	
	    BBHatsNewStep = NULL;	SigmaRecVec = NULL; PiRecVec = NULL;    
  
	    CurrentPostPBj = NULL;		     	     
	    CurrentPostPBj = (double *)Calloc((rkLen +1)+1, double);
      if (CurrentPostPBj == NULL) {
	      Rprintf((char*)"LarsConv Cannot Apportion CurrentPostPBj \n"); RFP();
			  return;
		  } 
		  RecordPostPBj = NULL;		     	     
	    RecordPostPBj = (double *)Calloc((rkLen +1)* TotalRuns, double);
      if (RecordPostPBj == NULL) {
	      Rprintf((char*)"LarsConv Cannot Apportion RecordPostPBj \n"); RFP();
			  return;
		  } 
      CurrentConfidenceInts = NULL;
      RecordConfidenceInts = NULL;          	     
	    CurrentConfidenceInts = (double *)Calloc((rkLen +1)*2+1, double);
      if (CurrentConfidenceInts == NULL) {
	      Rprintf((char*)"LarsConv Cannot Apportion CurrentConfidenceInts \n");
			  RFP(); return;
		  } 
		  RecordConfidenceInts = NULL;		     	     
	    RecordConfidenceInts = (double *)Calloc((rkLen +1)*2* TotalRuns, double);
      if (RecordConfidenceInts == NULL) {
	      Rprintf((char*)"LarsConv Cannot Apportion RecordConfidenceInts \n");
			  RFP(); return;
		  } 
		  CurGammaA = NULL;
	    CurGammaA = (double *)Calloc((rkLen +1)+1, double);
      if (CurGammaA == NULL) {
	      Rprintf((char*)"LarsConv Cannot Apportion CurGammaA \n");
			  RFP(); return;
		  }
		  CurGammaD = NULL;
	    CurGammaD = (double *)Calloc((rkLen +1)+1, double);
      if (CurGammaD == NULL) {
	      Rprintf((char*)"LarsConv Cannot Apportion CurGammaD \n");
			  RFP(); SuccessFlag = -1; return;
		  }		
		  if (FakeFlag >= 0) {	 
			  DiagXTX = NULL;
		    this->GenGAll();
		    this->SetDiagXTX();
		    this->UpdateGAll(0);
		    for (ii = 0; ii < kLen; ii++) {
			    this->OnMultiCons[ii] = .5;
		    }
		    if (this->OnLARS->GAllOnB == NULL) {
			    this->OnLARS->GAllOnB = (double *) Calloc(kLen+1, double);
		    }
			  for (ii = 0; ii < this->kLen; ii++) {
				  if (fabs(OnBetas[ii]) > 0) {
					  cnt1 = ii * this->kLen;
				    for (jj= 0; jj < this->kLen; jj++) {
					    this->OnLARS->GAllOnB[jj] += this->OnLARS->GAll[cnt1] * 
                this->OnLARS->OnBetas[ii];  
					    cnt1++;
				    }   
			    }
			  }
	    }
	    this->UpdateWAll((int) 1.0);
	    for (ii = 0; ii < kLen; ii++) { 
		    BBPrev[ii] = .5;       
		    BBOn1[ii] = .5;
	    }  
	    DConfidenceInt = -1.0;   
	    CDO = NULL;    	                      
    }

 ///////////////////////////////////////////////////////////////////
 //   CDO Convergent Lars Input
 //
 //
 ConvergentLARS(int rNLen, double * ryys, double *rxxs, int rkLen, 
   double *rInitBetas, double rStlambdaD, double rStlambdaA,
   double rsigmaNoiseSq,
   double rMultDConst, double rMultAConst,
   int rStopAInt, double rppiuse, int rTotalRuns,
   int rInitKKs, double InverseGammaConstant, 
   int PrintFlag) {

	 OnLARS = NULL; CDO = NULL;
	 LambdaDK = NULL;  RealXXS = NULL; yys = NULL;
	 LambdaAK = NULL; OrderSeq = NULL; OnMultiCons = NULL;
	 OnBetas = NULL; BBOn1 = NULL; OnMultiCons = NULL; BBOn1 = NULL;
	 BBPrev = NULL; BetasKeep = NULL; BBHatsKeep = NULL; BBHatsNewStep = NULL;
	 SigmaRecVec = NULL; PiRecVec = NULL; CurrentPostPBj = NULL;
   RecordPostPBj = NULL;    CurrentConfidenceInts = NULL;
   RecordConfidenceInts = NULL; CurGammaA = NULL; CurGammaD = NULL;                  
	 TDFNu = -1; iiWeights = NULL;
	 SigmaMultVec = NULL;
	 Lambda3Seq=NULL;  DiagXTX=NULL; 
	     
   if (PrintFlag >= 0) { IPLLF = PrintFlag; } else {IPLLF = -1;} // No Printing
   SuccessFlag = 1;    // Flag triggers to -1 when TwoLasso Breaks
  
   ChoiceCDOOnLARSFlag = 2;   // This Will be Coordinate Descent Based TwoLasso
          
   //////////////////  Now Insert Known Quantities
	 TotalRuns = rTotalRuns;
	 tt1 = 0; tt2 = 0;
	
   LambdaSeqF = 1;
	 FactorGammaDenominator = InverseGammaConstant;
	 NumEMConv = 0;
	  
	 int FeedKKs = 0;
	     	                      
   ppiuse = rppiuse;
   NusKeep = NULL; BBPrev = NULL;
        
   ChoiceCDOOnLARSFlag = 2; LambdaSeqF = 1;
   FactorGammaDenominator = InverseGammaConstant;
   NumEMConv = 0;
          
   tt1 = 0; tt2 = 0;
           
   FixKa = -1;

	 int ii; 
	 if (IPLLF > 0) {
		 Rprintf((char*)"Generating ConvergentLARS\n"); RFP();
   }
	 
   int cn1;
	 if (rNLen < 0) { 
		 if (IPLLF >= 0) {
			 Rprintf("ConvergentLARS CDO: Load GAll version only \n"); RFP();
		 }
		 Nlen = - rNLen;
		 FakeFlag = -1;
		 kLen = rkLen;
		 yys = NULL;
		 RealXXS = NULL;
		 GAll = NULL;
		 if (kLen < KLENCOPY) {
		   GAllPassByPointer = -1;
  		 GAll = (double *) Calloc((rkLen+1) * (rkLen) + 1, double);
  		 if (GAll == NULL) {
         Rprintf((char*)"LarsCov, Cannot Apportion GAll\b");  RFP();
  		   SuccessFlag = -1;  return;
  	   }
  	   cn1 = kLen * kLen;
  	   for (ii = 0; ii < cn1; ii++) {
  		   GAll[ii] = (double) rxxs[ii];
  	   }
  	} else {
      GAll = rxxs;
    }
		if (IPLLF >= 0) {
			Rprintf("ConvergentLARS: Loaded GAll\n");  RFP();
		}	          
	}  else {
	  FakeFlag = 1;
		Nlen = rNLen;
		kLen = rkLen;
		yys = NULL;
		if (kLen < KLENCOPY) {
		  yysPassByPointer = -1;
  		yys = (double *) Calloc( rNLen+2,double);
  		if (yys == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion yys \n"); RFP();
        return;
  		}
	    for (ii = 0; ii < rNLen; ii++) { 
        yys[ii] = (double) ryys[ii]; 
      }  			     
      RealXXS = NULL;		
      RealXXSPassByPointer = -1;     	      
  	  RealXXS = (double *) Calloc( (rNLen+1) * (rkLen+1), double);
  		if (RealXXS == NULL) {
        Rprintf((char*)"LarsConv Cannot Apportion RealXXS \n");
  			R_FlushConsole();
  			return;
  		}
	    cn1 = rNLen * rkLen;	
	    for (ii = 0; ii < cn1; ii++) { 
        RealXXS[ii] = (double) rxxs[ii]; 
      }             	 
    }  else {
      yysPassByPointer = 1; RealXXSPassByPointer = 1;
      yys = ryys; RealXXS = rxxs;
    }    
	  GAll = NULL;
	  m1 = 0; m2 = 0;

    SetupSumYYSq();
  }
  BetasKeep = NULL;	      
	BetasKeep = (double *)Calloc((rkLen+1)* TotalRuns, double);
  if (BetasKeep == NULL) {
    Rprintf((char*)"LarsConv Cannot Apportion BetasKeep \n"); RFP();
		return;
  }
      BBHatsKeep = NULL;		     	     
	    BBHatsKeep = (double *)Calloc((rkLen +1)* TotalRuns, double);
      if (BBHatsKeep == NULL) {Rprintf((char*)"LarsConv Cannot Apportion BBHatsKeep \n");
			                       R_FlushConsole(); SuccessFlag = -1;
			                       return;
		  }
		  NusKeep = NULL;
	    NusKeep = (double *)Calloc((rkLen +1)* TotalRuns, double);
      if (NusKeep == NULL) {
         Rprintf((char*)"LarsConv Cannot Apportion BBHatsKeep \n");
			   R_FlushConsole(); SuccessFlag = -1;
			    return;
		  }		 
		  OnBetas = NULL;
	    OnBetas = (double *)Calloc(kLen+1, double);
      if (OnBetas == NULL) {
          Rprintf((char*)"LarsConv Cannot Apportion OnBetas \n");
			    R_FlushConsole(); SuccessFlag = -1;
			    return;
		  }
		  if (rInitBetas == NULL) {
			  for (ii = 0; ii < kLen; ii++) {
				   OnBetas[ii] = 0;
	      } 
      } else {
	      for (ii = 0; ii < kLen; ii++) {
		         OnBetas[ii] = rInitBetas[ii];
	      }
      }
	     
      OnNus= NULL;		     	     
	    OnNus = (double *)Calloc(kLen+1, double);
      if (OnNus == NULL) {
         Rprintf((char*)"LarsConv Cannot Apportion OnNus \n");
			   R_FlushConsole(); SuccessFlag = -1;
			   return;
		  }	     
	    OnMultiCons = (double *)Calloc(kLen+1, double);
	    BBOn1 = NULL;
	    BBOn1 = (double *)Calloc(kLen+1, double);
      if (BBOn1 == NULL) {
         Rprintf((char*)"LarsConv Cannot Apportion BBOn1 \n");
			   R_FlushConsole(); SuccessFlag = -1;
			   return;
		  }	 
		  for (ii = 0; ii < kLen; ii++) {
         BBOn1[ii] = rppiuse;
      }    
	    BBPrev = (double *)Calloc(kLen+1, double);
	    
      sigmaNoiseSq = rsigmaNoiseSq; 

	    StlambdaD = (double) rStlambdaD;
	    StlambdaA = (double) rStlambdaA;
	    OnlambdaD = StlambdaD;
	    OnlambdaA = StlambdaA;
	    MultDConst = (double) rMultDConst;
	    MultAConst = (double) rMultAConst;
		  StopAInt = rStopAInt;
		  ppiuse = (double) rppiuse;
		   LambdaDK = NULL; LambdaAK = NULL; LambdaSeqF =1;
		   OrderSeq = NULL;
	      	     	     	     

	    OnLARS = NULL;
	    if (IPLLF > 1) { 
		     Rprintf((char*)"new ConvergentLars: Opening the CDO\n");
		     R_FlushConsole(); R_ProcessEvents();
	    }
	    if (rNLen < 0) {
		    if (IPLLF > 1) { 
			     Rprintf((char*)"new ConvergentLars: #1 CDO will be XTX and XTY\n");
			     R_FlushConsole(); R_ProcessEvents();
		    }
		       CDO = new CoordinateDescentObject(-Nlen, kLen, rxxs, ryys, 
		          rInitBetas, OnlambdaA / 2.0, (double *)NULL, (double*) NULL, PrintFlag);         		       
		       //CDO = new CoordinateDescentObject(-Nlen, kLen, rxxs, ryys, 
		       //   rInitBetas, OnlambdaA / 2.0, (double *)NULL, (double*) NULL, 
		       //   InverseGammaConstant);
  // CoordinateDescentObject(int rNLen, int rkLen, double *rXXX, double *rYYY,
  //           double *rOnBeta, double rOnGamma, double *rRecordOnBetas, 
  //           double *rOnGammas, int InitKKs, double InverseGammaConstant		          
  // CoordinateDescentObject(int rNLen, int rkLen, double *rXXX, double *rYYY,
  //           double *rOnBeta, double rOnGamma, double *rRecordOnBetas, 
  //            double *rOnGammas)		          
        if (CDO == NULL) {
		       Rprintf((char*)"MyConvLars: No room for CDO \n");
		       R_FlushConsole();
		       R_ProcessEvents();
		       return;
	      }
        CDO->SetUpGammas(BBOn1);		          
      } else if (rInitKKs > 0) {
           if (rInitKKs > kLen) {
             Rprintf("TwoLarsObject2009.h -- Interesting rInitKKs > kLen");
             Rprintf("   rInitKKs = %d, kLen = %d \n", rInitKKs, kLen);
             R_FlushConsole();
             rInitKKs = kLen;
           }
           CDO = new CoordinateDescentObject(Nlen, kLen, rxxs, ryys,
              rInitBetas, OnlambdaA/2.0, (double *) NULL, 
             (double*) NULL,  rInitKKs, InverseGammaConstant, PrintFlag);
           if (CDO == NULL) {
			       Rprintf((char*)"MyConvLars: No room for CDO \n");
			       R_FlushConsole();
			       R_ProcessEvents();
			       return;
	         }
           CDO->SetUpGammas(BBOn1);             
      } else if (Nlen < kLen) {
			     if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Opening the CDO, #2 xx Format\n");
				     R_FlushConsole(); R_ProcessEvents();
			     }	
           CDO = new CoordinateDescentObject(Nlen, kLen, rxxs, ryys, 
                   rInitBetas, OnlambdaA / 2.0, (double*) NULL, (double*) NULL,
                   PrintFlag);                      
              //CDO = new CoordinateDescentObject(Nlen, kLen, rxxs, ryys, 
              //     rInitBetas, OnlambdaA / 2.0, (double*) NULL, (double*) NULL,
              //     InverseGammaConstant); 
           if (CDO == NULL) {
			       Rprintf((char*)"MyConvLars: No room for CDO \n");
			       R_FlushConsole();
			       R_ProcessEvents();
			       return;
	         }
           CDO->SetUpGammas(BBOn1);
      } else if (Nlen < (kLen * 3)) {
           if (Nlen < 0) {
              Rprintf("TwoLarsObject2009.h, Interesting Nlen  = %d < 0??", Nlen);
              R_FlushConsole();
              Nlen = -Nlen;
           }
	         if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Opening the CDO, xx Format\n");
				     R_FlushConsole(); R_ProcessEvents();
			     }
			     FeedKKs = Nlen;
			     if (Nlen > 50) {
              FeedKKs = 50;
           }
           if (kLen < Nlen) {
              FeedKKs = kLen-1;
           }
           //if (Nlen > kLen) {
           //  Rprintf("TwoLarsObject2009.h -- Initeresting Nlen > kLen, but still we feed!");
           //  Rprintf("   Nlen = %d, kLen = %d  \n", Nlen, kLen);
           //  R_FlushConsole();
           //  FeedKKs = kLen;
           //}
    
           CDO = new CoordinateDescentObject(Nlen, kLen, rxxs, ryys,
              rInitBetas, OnlambdaA/2.0, (double *) NULL, 
             (double*) NULL,  FeedKKs, InverseGammaConstant,
             PrintFlag);          
          //   CDO = new CoordinateDescentObject(Nlen, kLen, rxxs, ryys, 
          //         rInitBetas, OnlambdaA / 2.0, (double*) NULL, (double*) NULL,
          //         InverseGammaConstant); 
           if (CDO == NULL) {
			       Rprintf((char*)"MyConvLars: No room for CDO \n");
			       R_FlushConsole();
			       R_ProcessEvents();
			       return;
	        }
          CDO->SetUpGammas(BBOn1);	 
      } else {
         	 if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Opening the CDO, #2 xx Format\n");
				     R_FlushConsole(); R_ProcessEvents();
			     }	
           CDO = new CoordinateDescentObject(Nlen, kLen, rxxs, ryys, 
                   rInitBetas, OnlambdaA / 2.0, (double*) NULL, (double*) NULL,
                   PrintFlag);                      
              //CDO = new CoordinateDescentObject(Nlen, kLen, rxxs, ryys, 
              //     rInitBetas, OnlambdaA / 2.0, (double*) NULL, (double*) NULL,
              //     InverseGammaConstant); 
           if (CDO == NULL) {
			       Rprintf((char*)"MyConvLars: No room for CDO \n");
			       R_FlushConsole();
			       R_ProcessEvents();
			       return;
	         }
           CDO->SetUpGammas(BBOn1);         
      }
      /* else {
			     if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Opening the CDO, #4 creating XTX\n");
				     R_FlushConsole(); R_ProcessEvents();
			     }	 	           
             GAll = (double *) Calloc(kLen * kLen + 4, double);
             if (GAll == NULL) {
                  Rprintf((char*)"MyConvLars: No Room for GAll\n");
                  R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
             }
			     if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Opening the CDO, #4 make GAll\n");
				     R_FlushConsole(); R_ProcessEvents();
			     }	              
             SqMat(GAll, Nlen, kLen, rxxs);
             //CDO = new CoordinateDescentObject(-Nlen, kLen, GAll, ryys,
             //           rInitBetas, OnlambdaA / 2.0, (double*) NULL, (double*) NULL,
             //           InverseGammaConstant);
             CDO = new CoordinateDescentObject(-Nlen, kLen, GAll, ryys,
                        rInitBetas, OnlambdaA / 2.0, (double*) NULL, 
                        (double*) NULL);                        
	         if (CDO == NULL) {
  		       Rprintf((char*)"MyConvLars: No room for CDO \n");
  		       R_FlushConsole();
  		       R_ProcessEvents();
  		       return;
	         }
  			  if (IPLLF > 1) { 
  				     Rprintf((char*)"new ConvergentLars: Opening the CDO, #4 creating XTY\n");
  				     R_FlushConsole(); R_ProcessEvents();
  			  }	 	         
          tMatTimesVec(kLen, Nlen,  rxxs, ryys, CDO->XTY); 
     }	*/ 
     if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Survived CDO, now setup Gammas\n");
				     R_FlushConsole(); R_ProcessEvents();
		 }	    	
     SuccessFlag = CDO->SetUpGammas(BBOn1);
     if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Survived CDO, setup SigmaRecVec/PiRecVec\n");
				     R_FlushConsole(); R_ProcessEvents();
		 }	
                  
	   BBHatsNewStep = NULL;	
	   SigmaRecVec = NULL;
	   PiRecVec = NULL;    
  
	   CurrentPostPBj = NULL;		     	     
	   CurrentPostPBj = (double *)Calloc((rkLen +1)+1, double);
     if (CurrentPostPBj == NULL) {
	        Rprintf((char*)"LarsConv Cannot Apportion CurrentPostPBj \n");
			    R_FlushConsole();
			    return;
		 } 
		 RecordPostPBj = NULL;		     	     
	   RecordPostPBj = (double *)Calloc((rkLen +1)* TotalRuns, double);
     if (RecordPostPBj == NULL) {
	               Rprintf((char*)"LarsConv Cannot Apportion RecordPostPBj \n");
			                       R_FlushConsole();
			                       return;
		 } 
     CurrentConfidenceInts = NULL;
     RecordConfidenceInts = NULL;         
	   //  CurrentConfidenceInts = NULL;		     	     
	   //  CurrentConfidenceInts = (double *)Calloc((rkLen +1)*2+1, double);
     //        if (CurrentConfidenceInts == NULL) {
	   //            Rprintf((char*)"LarsConv Cannot Apportion CurrentConfidenceInts \n");
		 //	                       R_FlushConsole();
		 //	                       return;
		 //    } 
		 //RecordConfidenceInts = NULL;		     	     
	   //  RecordConfidenceInts = (double *)Calloc((rkLen +1)*2* TotalRuns, double);
     //        if (RecordConfidenceInts == NULL) {
	   //            Rprintf((char*)"LarsConv Cannot Apportion RecordConfidenceInts \n");
		 //	                       R_FlushConsole();
		 //	                       return;
		 //    } 
		 CurGammaA = NULL;
	   CurGammaA = (double *)Calloc((rkLen +1)+1, double);
             if (CurGammaA == NULL) {
	               Rprintf((char*)"LarsConv Cannot Apportion CurGammaA \n");
			                       R_FlushConsole();
			                       return;
		 }
		 CurGammaD = NULL;
	   CurGammaD = (double *)Calloc((rkLen +1)+1, double);
     if (CurGammaD == NULL) {
	       Rprintf((char*)"LarsConv Cannot Apportion CurGammaD \n");
			   R_FlushConsole();
			   return;
		 }	
         if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Deciding to make DiagXTX?\n");
				     R_FlushConsole(); R_ProcessEvents();
	     }			     	
		 if (FakeFlag >= 0) {	 
			 DiagXTX = NULL;
		     //this->GenGAll();
		     //this->SetDiagXTX();
		     //this->UpdateGAll(0);
		     if (OnMultiCons != NULL) {
			     for (ii = 0; ii < kLen; ii++) {
				     this->OnMultiCons[ii] = .5;
			     }
	         }
	     }
	     //this->UpdateWAll(1.0);  //No need for weights with CoordinateDescent
	     for (ii = 0; ii < kLen; ii++) { 
		        BBPrev[ii] = .5;       
		        BBOn1[ii] = .5;
	     }  
	     DConfidenceInt = -1.0;     
	     if (IPLLF > 1) { 
				     Rprintf((char*)"new ConvergentLars: Made it to end of loadin\n");
				     R_FlushConsole(); R_ProcessEvents();
	     }	
	     ChoiceCDOOnLARSFlag = 2; OnLARS = NULL;		
	     //IPLLF = 5;
   }

   void CalculateBBOn1() {
	     int ii;
	     double OnlambdaDiff = (OnlambdaD - OnlambdaA) / 2.0;
	     int SuccFlag = 0;
	     if (FixKa > 0) {
		     SuccFlag = FixKCalculateBBOn();
		     if (SuccFlag < 0) {
			     Rprintf(" CalculateBBOn1, we got a return FixKCalculateBBOn is angry \n");
			     R_FlushConsole(); return;
		     }
		     return;
	     }
	     if (OnlambdaD < 100 && OnlambdaA > .001) {
	     for (ii = 0; ii < kLen; ii++) {
		   this->BBOn1[ii] = 
		       ( this->ppiuse)* this->OnlambdaA * 
		                                exp(-fabs(OnBetas[ii]) * this->OnlambdaA);
		   this->BBOn1[ii] = this->BBOn1[ii] / ( this->BBOn1[ii] + 
		    (1.0-this->ppiuse) * this->OnlambdaD * exp(-fabs(OnBetas[ii]) * this->OnlambdaD));	
		    if (ISNAN( BBOn1[ii])) {
			    if (fabs(OnBetas[ii]) > 0) { BBOn1[ii] = 1.0; } else { BBOn1[ii] = 
			         OnlambdaA * this->ppiuse / 
			         (OnlambdaA *this->ppiuse + OnlambdaD * (1.0-this->ppiuse)); }
		    }	     
          }
         } else {
	       
	       for (ii = 0; ii < kLen; ii++) {
		      if (OnlambdaDiff * OnBetas[ii] > 900) {
			      BBOn1[ii] = 1.0;
		      } else {
              this->BBOn1[ii] = 
		       ( this->ppiuse)* this->OnlambdaA * 
		                                exp(fabs(OnBetas[ii]) * (OnlambdaDiff));		      
		      this->BBOn1[ii] = this->BBOn1[ii] / ( this->BBOn1[ii] + 
		      (1.0-this->ppiuse) * this->OnlambdaD * exp(-fabs(OnBetas[ii]) * (OnlambdaDiff)));		        
              if (ISNAN(BBOn1[ii])) {
			    if (fabs(OnBetas[ii]) > 0) { BBOn1[ii] = 1.0; } else { BBOn1[ii] = 
			         OnlambdaA * this->ppiuse / 
			         (OnlambdaA *this->ppiuse + OnlambdaD * (1.0-this->ppiuse)); }
		      }
	         }	   		       
	       } 	         
         }
	     return;
  }
  void ConvertOnNusOnBetas() {
           if (OnLARS == NULL) {
            Rprintf("ConvergentLARS::ConvertOnNusOnBetas, OnLARS NULL\n");
            SuccessFlag = -1; R_FlushConsole();  R_ProcessEvents();
            return;
           }
	    int ii;
	    for (ii = 0; ii < kLen; ii++) {OnNus[ii] = OnLARS->OnBetas[ii];}
	    for (ii = 0; ii < kLen; ii++) {OnBetas[ii] = OnNus[ii] * OnLARS->WDiag[ii];}	       
  }
  
  void ConvertCreateOnMultiCons() ;
	
  ~ConvergentLARS() {
      //if (iiWeights != NULL || TDFNu > 0) {
      //  IPLLF = 8;
      //}
  
       if (yysPassByPointer == 1) {
           yys = NULL;
       }
       if (RealXXSPassByPointer == 1) {
            RealXXS = NULL;
       }
       if (GAllPassByPointer == 1) {
            GAll = NULL;
       }
       //if (XTYPassByPointer == 1) {
       //     XTY = NULL;
       //}
		    if (IPLLF > 0) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Start\n");
			   R_FlushConsole();
	      }  
		    if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete LambdaDK\n");
			   R_FlushConsole();
	      }
	       FFree(&LambdaDK);
	       FFree(&LambdaAK);
		    if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete OnLARS?\n");
			   R_FlushConsole();
	      }
      if (OnLARS != NULL) {
	         delete(OnLARS); OnLARS = NULL;
      }
		  if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete CDO?\n");
			   R_FlushConsole();
	    }         
      if (CDO != NULL) {
		    if (IPLLF > 2) {
			     Rprintf((char*)"2LO->Deleting: Defenitely deleting CDO!\n");
			     R_FlushConsole();
	      }        
              delete (CDO); CDO = NULL;
      }
		  if (IPLLF > 2 && Lambda3Seq != NULL) {
			   Rprintf((char*)"2LO->Deleting: Delete Lambda3Seq\n");
			   R_FlushConsole();
	    } 
      if (Lambda3Seq != NULL ) {
        Free(Lambda3Seq); Lambda3Seq = NULL;
      }
		  if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete BBOn1?\n");
			   R_FlushConsole();
	    }               
        FFree(&BBOn1);
      if (IPLLF > 2 && BBPrev!= NULL) {
        Rprintf("2LO->Free NonNull BBPrev\n"); R_FlushConsole();
      } 
      FFree( &BBPrev);
      if (IPLLF > 2 && OnBetas!=NULL) {
        Rprintf("2LO->Free NonNull OnBetas\n"); R_FlushConsole();
      }
      FFree(&OnBetas);
		  if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete OnMultiCons?\n");
			   R_FlushConsole();
	    } 
        FFree(&OnNus); FFree( &OnMultiCons); FFree( &RealXXS);
        FFree(&yys); FFree(&BBHatsNewStep); 
        
      if (LSigmaRecVec == 0) { 
        FFree(&SigmaRecVec);  
      } else { SigmaRecVec = NULL; }
        FFree(&PiRecVec); FFree(&NusKeep);
		  if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete BBOn1 second?\n");
			   R_FlushConsole();
	    }        
	    if (BBOn1 != NULL) {Free(BBOn1); BBOn1 = NULL;}  //3
		  if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete BBPrev second?\n");
			   R_FlushConsole();
	    } 	    
	  if (BBPrev!=NULL) {Free(BBPrev); BBPrev =NULL;}	 //4
		if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete OnBetas?\n");
			   R_FlushConsole();
	  } 	    
	    if (OnBetas!=NULL) {Free(OnBetas); OnBetas=NULL;}//5
	    if (OnNus!=NULL) {Free(OnNus); OnNus = NULL;}  //6
		  if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete OnMultiCons?\n");
			   R_FlushConsole();
	    } 	    
	    if (OnMultiCons!=NULL){Free(OnMultiCons);OnMultiCons=NULL;}//7
		  if (IPLLF > 2 && BetasKeep != NULL) {
			   Rprintf((char*)"2LO->Free NonNull BetasKeep?\n");
			   R_FlushConsole();
	    } 	    
	    if (BetasKeep!=NULL){Free(BetasKeep); BetasKeep=NULL;}//8
		  if (IPLLF > 2 && BBHatsKeep != NULL) {
			   Rprintf((char*)"2LO0>Free NonNull BBHasKeep?\n");
			   R_FlushConsole();
	    } 	    
	    if (BBHatsKeep!=NULL) {Free(BBHatsKeep);BBHatsKeep=NULL;}//9
	    if (NusKeep!=NULL) {Free(NusKeep); NusKeep = NULL;}
		  if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete Data Fills?\n");
			   R_FlushConsole();
	    } 	    
	    if (kLen < KLENCOPY && RealXXS != NULL) {
		      Free(RealXXS);
		      RealXXS = NULL;//10
	    }
	    if (kLen < KLENCOPY && yys!=NULL){Free(yys); yys=NULL;} //11
	    if (kLen < KLENCOPY && GAll != NULL) {
		       Free(GAll); GAll = NULL;//12
	  	}
		if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, BBHatsNewStep\n");
			   R_FlushConsole();
	  }
        if (BBHatsNewStep != NULL) {Free(BBHatsNewStep); BBHatsNewStep = NULL; } //13	
        if (SigmaRecVec != NULL) {Free(SigmaRecVec); SigmaRecVec = NULL; } //14
        if (PiRecVec != NULL) {Free(PiRecVec); PiRecVec = NULL; } //15
	    if (CurrentConfidenceInts != NULL) {
		        Free(CurrentConfidenceInts); CurrentConfidenceInts = NULL;} // 16
        if (RecordConfidenceInts != NULL) { 
	        Free(RecordConfidenceInts); RecordConfidenceInts = NULL;} // 17;    
		if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Delete DiagXTX?\n");
			   R_FlushConsole();
	  }
	    if (DiagXTX != NULL) {
		    Free(DiagXTX); DiagXTX = NULL; } // 18;  
		if (IPLLF > 2) {
			   Rprintf((char*)"2LO->Deleting: ConvergentLARS, Delete LambdaDK again\n");
			   R_FlushConsole();
	  }		    
		if (LambdaDK != NULL) {
			Free(LambdaDK); LambdaDK = NULL; }
		if (LambdaAK != NULL) {
			Free(LambdaAK); LambdaAK = NULL; }	// 19;   
		if (OrderSeq != NULL) {
			Free(OrderSeq); OrderSeq = NULL; }	// 20;
    if (IPLLF > 2) {
      Rprintf("2LO->Deleting CurrentPostPBj\n"); R_FlushConsole();
    }   	
	    FFree(&CurrentPostPBj);
	  if (IPLLF > 2) {
      Rprintf("2LO->Deleting CurGammaA and CurGammaD\n"); R_FlushConsole();
    }
    if (IPLLF > 2 && CurGammaA != NULL) {
      Rprintf("2LO-> Freeing Non Null CurGammaA\n"); R_FlushConsole();
    }
	    FFree(&CurGammaA); 
    if (IPLLF > 2 && CurGammaD != NULL) {
      Rprintf("2LO-> Freeing Non Null CurGammaD\n"); R_FlushConsole();
    }     FFree(&CurGammaD);
		if (IPLLF > 2 && BetasKeep != NULL) {
			   Rprintf((char*)"2LO->Free NonNull BetasKeep\n");
			   R_FlushConsole();
	  }       	
	  FFree(&BetasKeep);
		if (IPLLF > 2) {
			   Rprintf((char*)"2LO->Deleting: ConvergentLARS, Delete RecordPostPBj\n");
			   R_FlushConsole();
	  }
	  FFree(&RecordPostPBj);    FFree(&CurrentConfidenceInts);
		if (IPLLF > 2 && RecordConfidenceInts != NULL) {
			   Rprintf((char*)"2LO->Deleting: NonNull");
         Rprintf(" ConvergentLARS, Delete RecordConfidenceInts\n");
			   R_FlushConsole();
	  }	  
    FFree(&RecordConfidenceInts); 
    if (IPLLF > 2 && DiagXTX != NULL) {
      Rprintf((char*)"2LO->Free NonNull DiagXTX\n"); R_FlushConsole();
    }
    FFree(&DiagXTX);
    
    if (IPLLF > 2 && OnNus != NULL) {
       Rprintf((char*)"2LO->Free NonNull OnNus \n"); R_FlushConsole();
    }
    FFree(&OnNus); 
    if (IPLLF > 2 && GAll != NULL) {
       Rprintf((char*)"2LO->Free NonNull GAll\n"); R_FlushConsole();
    }
    FFree(&GAll);    
    if (IPLLF > 2 && BBPrev!= NULL ) {
      Rprintf((char*)"2LO-> Free NonNull BBPrev\n"); R_FlushConsole();
    }
    FFree(&BBPrev);
    if (IPLLF > 2 && SigmaMultVec != NULL ) {
      Rprintf((char*)"2LO-> Free NonNull SigmaMultVec\n"); R_FlushConsole();
    }
    if (SigmaMultVec != NULL) {
        FFree(&SigmaMultVec);
    }
	  if (IPLLF > 2 && iiWeights != NULL) {
       Rprintf("Delete: ConvergentLARS: deleting nonNull iiWeights\n");
       R_FlushConsole();
    }
    if (iiWeights != NULL) {
       Free(iiWeights);  iiWeights = NULL;
    }
    if (IPLLF > 2) {
			   Rprintf((char*)"Deleting: ConvergentLARS, Done Deleting\n");
			   R_FlushConsole();
	  }     
    }	
};    



int ReviseBHat1(int kLen, double *CurBHat, double *NewBHat, 
               double *CurBetas,  double LambdaA, double LambdaD, 
               double m1, double m2);
               
                
               

