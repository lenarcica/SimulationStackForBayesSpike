#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <Rdefines.h>
  #include <Rinternals.h>
  #include <stdio.h>
  #define RMATH 0
#endif

#ifndef INCIDENTALOP2009DD
  #include "IncidentalOp2009.h"
  #define INCIDENTALOP2009DD 0 
#endif
#ifndef MyMatrixOp2009DD
  #include "MyMatrixOp2009.h"
  #define MyMatrixOp2009DD 0
#endif


#ifndef DefaultMaximumAllocation
  #define DefaultMaximumAllocation 4096
  #define DefaultStartRemoval   1024
#endif

#ifndef LARGEDATA
  #define LARGEDATA 50000
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

int CoordinateDescentLasso(int NLen, int kLen,  double *XXX, double *YYY,
             double *OnBeta, double OnGamma, double *OnRecordBetas,
             double *FinalBeta, int *pTotalLoops, double MaxEpsilon, 
             double *OnGammas, int InitKKs, double InverseGammaConstant,
             double *iiWeights, int PrintFlag);    
SEXP GLMLasso(SEXP sX, SEXP sy01, SEXP sOnBeta, SEXP sBeta0,
   SEXP sBeta0prev, SEXP sBetaprev, SEXP sOnGamma,
   SEXP sRecordOnBetas,SEXP sBeta0Records,
   SEXP sOnGammas,
   SEXP sInitKKs, SEXP sInverseGammaConstant,
   SEXP snProbWeights, SEXP slogitCurrentProb, 
   SEXP snlogitPrevProb, SEXP rZ, SEXP sPrintFlag,
   SEXP sMaxCDOEpsilon, SEXP sMaxCDORuns,
   SEXP sIter);     
             
class CoordinateDescentObject {
  private: 
  public: 
  
   int ModifiedXTResid;   // Flag is 0 in case XTResid needs to be modified
   double *OnGammas;      // Individual weights for each beta estimate
   int DOnGammas;  // Flag whether pointer points to OnGammas Elsewhere
   int PrintFlag;         // A flag whether to print echo for algorithms
   
   int TriggerReallocFlag;  int SuccessDOP;
   int TriggerInvestigateFlag;
   int DataIntegrity(int DoIn);
   
   int XTXFlag;           // Flag whether XTX is being dynamically scaled
   
   double **pXTX;           // Matrix holds Squared errors between covariates
                          //    XTX_{ij} = sum_t X_{ti} *X_{tj}                         
   double DiagShrinkage;
   
   int XTXPassByPointer;
   double *XTY;           // Vector of correlations with response
                          //    XTY_{j} = sum_t X_ti Y_t
   int XTYPassByPointer;                       
   double *YY, *XX;       // Y is vector of data, XX is matrix of covariates
   int YYPassByPointer; int XXPassByPointer;
   double *XTResid;       // Holds  XTY -  XTX * Beta, current residuals
   double *OnBeta;        // Current BetaEstimate
   int DOnBeta;        // Whether Beta Estimate is a pointer or allocated memory
   double *pBeta0;
   
   int RecordBetasPassByPointer;
   double *OnRecordBetas;  // Record of past Beta Estimates
   double CurrentXTResid;  // Measure of XTResid for current coordinate
   double RecentMove;      // Measures sum of total moves in last iteration
                           //    used to measure convergence
                           
   double OnGamma, OnLambda; // Single Penalties, Gamma = Lambda/2
   //double OnLambda2;  // SecondLevelPenalty
   
   double *InvSXX;        // Vector of the multiplicative inverse of XTX 
   
   int TestCDO();  // Test Elements of CDO for Data integrity. 
   
   int MakeXTResid();     // Updates XTResid matrix from begining
                          //     Could be conceivably sped up by BLAS
   int MemoryFailFlag;
   
   ///////////////////////////////////////////////////////////////                       
   double *iiWeights;    // Useful For Doing WLS Lasso estimates
   int iiWeightsHere;
   double *rWXVec;       // Including doing t-Lasso;
   int SetupWeights(double *riiWeights); 
   int RefreshWeights(double *riiWeights);
   int *iWeightedXtX;
   
   int MaximumAllocation;
   
   double ShowSFunction(int Coord);
   inline double SFunction(int Coord); // Calculates where to Beta for 
                                //   current coordinate
   int UpdateXTY;   // Flag whether XTY needs to be updated
   int SuccessFlag; // Flag whether algorithm is running successfuly
                    //   When SuccessFlag goes negative, algorithm will exit
   int NLen, kLen;  // NLen = # of samples, kLen = # of Covariates
   
   int OnCoord; int OnLoop; // OnCoord, current Beta_i coordinate
                            //   OnLoop,  current Loop through all coordinates                          
   double TotMove;
                            
   ////////////////////////////////////////////////////
   //  RunAlgorithm: Key Algorithm to Running Coordinate Descent Lasso                         
   int RunAlgorithm(int MaxLoops, double MaxEpsilon);
   int RunAlgorithmDiagonostic(int MaxLoops, double MaxEpsilon);
   
   int ReweightCoordinate(int OnCoord);
   
   ////////////////////////////////////////////////////////
   //  UpdateCoord();  
   //    Gets new value for Beta_i for i = OnCoord
   int UpDateCoord(); 
   ////////////////////////////////////////////////////////
   //  MoveCoord();
   //    Selects a new OnCoord in sequence
   inline int MoveCoord(int ii);

   double SFunction3(int OnCoord);    int UpDateCoord3(); 
   double OnLambda3; double *OnLambdas3;
   int SetupOnLambdas3(double *rOnLambdas3) {
      OnLambdas3 = (double *) Calloc(kLen, double);
      if (OnLambdas3 == NULL) {
        Rprintf("SetupOnLambdas3: Failed Setup memory\n"); R_FlushConsole();
        SuccessFlag = -1;
        return(-1);
      }
      for (int ii = 0; ii < kLen; ii++) {
         OnLambdas3[ii] = rOnLambdas3[ii];
      }
      return(1);
   }

   /////////////////////////////////////////////////////////
   //   Dynamic Memory Management
   int *XLC;   // XLC and kXFinder store locations of current X_ij in XTX
   int *kXFinder; 
   int OnKappaS; int OnKappaMem; // Total number of covariates in memory
   int InitKappaMem;
   
   int *RunOrder;
   double OnEp;
   void SetupRunOrder(int *rRunOrder) {
     RunOrder = rRunOrder;
   }
   
   /////////////////////////////////////////////////////////
   //  SetUpGammas: Sets up system for indiviual weighting
   //int SetUpGammas(double *InputGammas);
   int SetUpGammas(double *InputGammas);
   
   //////////////////////////////////////////////////////////
   //  Adds Coordinate NewCoord to dynamic space in memory
   int AddCoord(int NewCoord) ;
   int DoubleMemory(); // Doubles XTX size in memory for more readins
   
   ///////////////////////////////////////////////////////////////
   //  SetUpToPartial:  Opens Coordinate Descent version expanding 
   ///   in dynamic  memory management
   int SetUpToPartial(int InitKKs);
   int AddAllNewNonzeroCoords();
   
   //////////////////////////////////////////////////////////////
   //  FactorGammaDenominator: An inverse of current OnGammas[OnCoord]
   double FactorGammaDenominator;

   //////////////////////////////////////////////////////////////////
   //  Sets up new 
   int SetEpsilonAndLoops(double rMaxEpsilon, double rMaxLoops) {
	   MaxEpsilon = rMaxEpsilon;
	   MaxLoops = (int) rMaxLoops;
	   return(1);
   }
   double GMaxEpsilon() {
	   return(MaxEpsilon);
   }
   int GMaxLoops() {
	   return(MaxLoops);
   }
   int StatusReport();
   public :
   CoordinateDescentObject(){
    return;
   }
   ///////////////////////////////////////////////////////////////
   //  CoordinateDescentObject:: Constructor
   //    If NLen is < 0, constructor assumes one is given pXTX and XTY 
   //     (i.e. the product matrices)
   //     Else, one assumes constructor has been given x matrix and y matrix
   //     Depending on x structure, then appropriate decision is made whether to
   //     Convert x immediately into pXTX or partially into pXTX
   //     Default InitKappas are used.
   CoordinateDescentObject(int rNLen, int rkLen, double *rXXX, double *rYYY,
             double *rOnBeta, double rOnGamma, double *rRecordOnBetas, 
             double *rOnGammas, int rPrintFlag) {
	   int One = 1; double OneD = 1.0; double ZeroD = 0.0; rWXVec = NULL; 
	   OnLambda3=0; OnLambdas3 = NULL;               
	   ModifiedXTResid = 0;   DiagShrinkage  =0.0;   FactorGammaDenominator = 1.0;
	   //FactorGammaDenominator = InverseGammaConstant;
	   
     //OnLambda2 = 0.0;
     TriggerReallocFlag=0;  SuccessDOP = 0;
     TriggerInvestigateFlag=0;
	   RunOrder = NULL;
	   // Currently Echoes are all Off
	   PrintFlag = rPrintFlag;  OnEp =0.0;
	   //PrintFlag = 3;
	   if (PrintFlag > 0) {
       Rprintf("CoordinateDescentObject #1: Starting Constructor\n;");
       R_FlushConsole();
     }
	   /////////////////////////////////////////////////////////
	   //  Set initial values for fields
	   pXTX = NULL; XTY=NULL; YY = NULL; XX = NULL; InvSXX = NULL;
     kLen = rkLen; XTResid = NULL;     TDFNoise = -1.0; TDFSigma = -1.0;
     XLC = NULL; kXFinder = NULL;  OnBeta = NULL;
     kLen = rkLen;
     OnKappaS = rkLen; OnKappaMem = rkLen;
	   SuccessFlag = 1;
	   XTXFlag  = -1;
     OnGammas = NULL;
     pBeta0 =NULL;  MemoryFailFlag = 0;   TDFNoise = -1.0; TDFSigma = -1.0;
     MaximumAllocation = DefaultMaximumAllocation;
     
     YYPassByPointer = 0; XXPassByPointer = 0; RecordBetasPassByPointer = 0;
     XTYPassByPointer = 0; XTXPassByPointer = 0;
       
     iiWeights = NULL;  rWXVec = NULL;  iiWeightsHere = 0;
     iWeightedXtX = NULL;
     
     // When rNLen input is negative, this means that supplied rXXX is acutally
     //   an pXTX symmetric matrix, and rYYY is XTY. 
     if (rNLen < 0) { NLen = -rNLen; } else {NLen = rNLen; }
     
     //int Maxss;
     int ii;        

     ////////////////////////////////////////////////////////////////////////
     //  Allocate OnBeta, If rOnBeta is non null, inputs OnBeta values
     if (rOnBeta != NULL && rOnBeta[0] != -999) {
       DOnBeta = 1;
       OnBeta = rOnBeta;
       //F77_CALL(dcopy)(&kLen, rOnBeta, &One, OnBeta, &One);
		   ModifiedXTResid = 1;        
		   if (PrintFlag > 0) {
         Rprintf("CoordinateDescentObject #1: Loaded DOnBeta = 1\n");
         R_FlushConsole();
       }
     } else {
       DOnBeta = 0;
       OnBeta = (double *) Calloc( kLen+2, double);
       if (OnBeta == NULL) {
		     Rprintf("CoordinateDescentObject #1:: Calloc Fails on OnBeta \n");
		     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
	     }
       int AOne = 1; double ZeroD = 0.0;
       F77_CALL(dscal)(&kLen, &ZeroD, OnBeta, &AOne);	
       DOnBeta = 0;
	   }   
	      
       ///////////////////////////////////////////////////////////////////////
       //  Positive Covariate set must be supplied     
     if (rkLen <= 0) {
	     Rprintf((char*)"new CoordinateDescentObject #1:: Error rkLen < 0 = %d \n", rkLen);
	     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1;
	     return;
     }
       
     ///////////////////////////////////////////////////////////////////
     // Echo Printing
     if (PrintFlag > 1) {
			     Rprintf((char*) "new CoordinateDescentObject: Loading in\n");
			     Rprintf((char*) "rNLen = %d, rkLen = %d, rXXX[0] = %.4f, rYYY[0] = %.4f\n",
			                 rNLen, rkLen, (double) rXXX[0], (double) rYYY[0]);
			     Rprintf((char*) "rOnBeta[0] = %.4f, rOnGamma  = %.4f\n", 
			                      rOnBeta[0], rOnGamma);
			     if (rRecordOnBetas == NULL) {
				       Rprintf((char*) "rRecordOnBetas is NULL\n");
				       R_FlushConsole();
			     } else if (rRecordOnBetas != NULL && rRecordOnBetas[0] < 0 ) {
				       Rprintf((char*)"rRecordOnBetas Negative, rRecordOnBetas[0] = %.4f\n",
				                (double) rRecordOnBetas[0] );
			     } else if (rRecordOnBetas != NULL ) {
				       Rprintf((char*)"rRecordOnBetas Positive, rRecordOnBetas[0] = %.4f\n",
				                (double) rRecordOnBetas[0] );
			     }
			     R_FlushConsole(); R_ProcessEvents();
	     }
	     
	     /////////////////////////////////////////////////
	     //  Load weighted gammas if applicable
       if (rOnGammas != NULL && rOnGammas[0] >= 0) {
         if (PrintFlag > 1) {
           Rprintf((char*) "new CoordinateDescentObject:: rOnGammas != NULL, Setting up Gammas\n");
           R_FlushConsole(); R_ProcessEvents();
         }
         SetUpGammas(rOnGammas);
       } else {
	      if (PrintFlag > 1) {
		      Rprintf((char*)" new CoordinateDescentObject:: rOnGammas == NULL or Less than Zero\n");
		      R_FlushConsole(); R_ProcessEvents();
	      }    
       }
       OnGamma = rOnGamma; OnLambda = OnGamma * 2;
       kLen = rkLen;
       RecentMove = 0;
       	 if (PrintFlag > 1) {
		     Rprintf((char*) "new CoordinateDescentObject #1: Deciding XX, pXTX\n");
		     R_FlushConsole(); R_ProcessEvents();
	     }
     if (rNLen == 0) {
	     Rprintf((char*) "new CoordinateDescentObject #1:: Error, NLen == 0\n");
	     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
     } else if (rNLen < 0) {  // If We're Fed pXTX
	     if (PrintFlag > 1) {
		     Rprintf((char*) "new CoordinateDescentObject #1: NLen < 0, Using Given XTX\n");
		     R_FlushConsole(); R_ProcessEvents();
	     }
	     NLen = - rNLen;  XTXFlag = 1;
	     kLen = rkLen;
	     XTXPassByPointer = 1; 
       pXTX = (double**) Calloc(kLen+1, double *);
	     if (pXTX == NULL) {
		     Rprintf("CoordinateDescentObject #1:: Calloc Fails on pXTX \n");
		     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
	     }
       OnKappaMem = kLen;  OnKappaS = kLen;
       for (int iit = 0; iit < kLen; iit++) {
         pXTX[iit] = rXXX + iit * kLen;
	     }

	     if (PrintFlag > 1) {
		     Rprintf((char*) "new CoordinateDescentObject #1: NLen < 0, Using Given XTY\n");
		     R_FlushConsole(); R_ProcessEvents();
	     }
	     
       XTYPassByPointer = 1; XTY = rYYY;	     
	     //XTY = (double *) Calloc( kLen + 6, double);
	     if (XTY == NULL) {
		     Rprintf("CoordinateDescentObject #1:: Calloc Fails on XTY \n");
		     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
	     }	

	     //Maxss = (kLen * kLen); 
	     //      F77_CALL(dcopy)( &Maxss,rXXX, &One, pXTX, &One );
	     //      F77_CALL(dcopy)( &kLen,rYYY, &One, XTY, &One );
	    //Rprintf("XTY after Fortran \n");
	    //PrintVector(XTY, kLen);
	    // for (ii = 0; ii < kLen; ii++) { pXTY[ii] = (double) rYYY[ii]; }
	    if (PrintFlag > 0) {
        Rprintf("CoordinateDescentObject #1: Memory for XTResid\n"); R_FlushConsole();
      }
      XTResid = NULL; XTResid = (double *) Calloc( kLen, double);
		      	if (XTResid == NULL) {
				     Rprintf("CoordinateDescentObject #1:: Calloc Fails on XTResid \n");
				     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
			}	
			 F77_CALL(dcopy)( &kLen, XTY, &One, XTResid, &One );
	    // for (ii = 0; ii < kLen; ii++) { XTResid[ii] = (double) XTY[ii]; }	
		 InvSXX = (double*) Calloc(kLen+2, double);
		 if (InvSXX == NULL) {
			 R_FlushConsole(); Rprintf("CoordinateDescentObject #1:: Calloc Fails on InvSXX\n");
			 R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
         }
			   for (ii = 0; ii < kLen; ii++) { InvSXX[ii] = ((double)1.0 / 
			                       ( (double) (*pXTX[ii] + ii))); }     
      } else {
	      	 if (PrintFlag > 1) {
			     Rprintf((char*) "new CoordinateDescentObject #1: Loading XX\n");
			     R_FlushConsole(); R_ProcessEvents();
	         }
	       NLen = rNLen; XTXFlag = 0;
	       kLen = rkLen;
	       if (PrintFlag > 1) {
		       Rprintf((char*) "new CoordinateDescentObject #1: NLen > 0, loading XX, N=%d, k=%d\n", 
		         NLen, kLen);
		       R_FlushConsole(); R_ProcessEvents();
	       }

	     //XX = (double *) Calloc( kLen * NLen + 4, double);
       //
	     //if (XX == NULL) {
		   //  Rprintf("CoordinateDescentObject:: Calloc Fails on XX, N = %d, k=%d \n",
		   //      NLen, kLen);
		   //  R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
	     //}
	     if (PrintFlag > 1) {
         Rprintf("new CoordinateDescentObject#1: Passing X by Pointer\n");
         R_FlushConsole();
       }
	     XXPassByPointer = 1;
	     XX = rXXX;
	     if (PrintFlag > 1) {
		     Rprintf((char*) "new CoordinateDescentObject: NLen > 0, loading YY, k=%d\n", kLen);
		     R_FlushConsole(); R_ProcessEvents();
	     }	     
	     //YY = (double *) Calloc( NLen, double);
	     //if (YY == NULL) {
		   //  Rprintf("CoordinateDescentObject:: Calloc Fails on YY \n");
		   //  R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
	     //}
	     
	     YYPassByPointer = 1;
       YY = rYYY;	     
	     //Maxss = (kLen * NLen); 
	    // F77_CALL(dcopy)( &Maxss, rXXX, &One, XX, &One );
	    // F77_CALL(dcopy)( &NLen, rYYY, &One, YY, &One );	     
	     //for (ii = 0; ii < Maxss; ii++) { XX[ii] = rXXX[ii]; }
	     //for (ii = 0; ii < NLen; ii++)  { YY[ii] = rYYY[ii]; }	   
	     if (PrintFlag > 1) {
		     Rprintf((char*) "new CoordinateDescentObject #1: NLen > 0, XX,YY Loaded\n");
		     R_FlushConsole(); R_ProcessEvents();
	     }	   
	     if (PrintFlag > 3) {
		      Rprintf((char*)"Printing XX \n");
		      PrintRMatrix(XX, NLen, kLen);
		      Rprintf((char*)"Printing YY\n");
		      PrintVector(YY, NLen);
	     }			          
	     if (kLen > SMALLESTkLen) {
		     if (PrintFlag > 1) {
			     Rprintf((char*) "new CoordinateDescentObject #1: Setting Up Partial\n");
			     R_FlushConsole(); R_ProcessEvents();
	       }
			   SetUpToPartial(SMALLESTInitKKs);	     		     
	     }  else if (NLen > kLen) {
		     if (PrintFlag > 1) {
			     Rprintf((char*) "new CoordinateDescentObject: NLen > kLen, supplied pXTX version\n");
			     R_FlushConsole(); R_ProcessEvents();
	       }
		     XTXFlag = 0; XTXFlag = 1;
		     
		     XTXPassByPointer = -1;
		     pXTX = (double **) Calloc( kLen + 4, double*);
		     if (pXTX == NULL) {
			     Rprintf("CoordinateDescentObject:: Calloc Fails on supplied XTX  version\n");
			     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
		     }
         for (int jj = 0; jj < kLen; jj++) {
           pXTX[jj] = (double *) Calloc(kLen, double);
         }
		     OnKappaMem = kLen;
		     XTYPassByPointer = -1;
		     XTY = (double *) Calloc( kLen + 4, double);
		     if (XTY == NULL) {
			     Rprintf("CoordinateDescentObject:: Calloc Fails on XTY \n");
			     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
		     }
     //F77_NAME(dsyrk)(const char *uplo, const char *trans,
		 //  const int *n, const int *k,
		 //  const double *alpha, const double *a, const int *lda,
		 // const double *beta, double *c, const int *ldc);		         
		     pSqMat (pXTX, NLen, kLen, XX);
	  if (PrintFlag > 1) {
	       Rprintf((char*) "CoordinateDescentOjbect:PostSquare, XTX is = \n");
	       PrintRpMatrix(this->pXTX, this->kLen, this->kLen);
	       R_FlushConsole();
      }		     
		     tMatTimesVec(kLen, NLen,  (double*) XX, (double*) YY, (double*) XTY);
	  if (PrintFlag > 1) {
	       Rprintf((char*) "CoordinateDescentOjbect:PostSquare, XTY is = \n");
	       PrintVector(this->XTY, this->kLen);
	       R_FlushConsole();
      }
		     XTXFlag = 1;	
	         XTResid = NULL; XTResid = (double *) Calloc( kLen, double);
		      	if (XTResid == NULL) {
				     Rprintf("CoordinateDescentObject:: Calloc Fails on XTResid \n");
				     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
			    }	
			for (ii = 0; ii < kLen; ii++) { XTResid[ii] = (double) XTY[ii]; }	
  			InvSXX = (double*) Calloc(kLen+2, double);
	  	    if (InvSXX == NULL) {
			     Rprintf("CoordinateDescentObject #1 :: Calloc Fails on InvSXX \n");
			     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
		    }
		    for (ii = 0; ii < kLen; ii++) { InvSXX[ii] = (double) 1.0 / (
		                                      (double) *(pXTX[ii] + ii)); }
	     } else {
		     if (PrintFlag > 1) {
			     Rprintf((char*) "new CoordinateDescentObject #1: else-load pXTX\n");
			     R_FlushConsole(); R_ProcessEvents();
	         }
		     XTXFlag = 1;     XTXPassByPointer = 0;
	       pXTX = (double **) Calloc( kLen +1, double*);
		     if (pXTX == NULL) {
			     Rprintf("CoordinateDescentObject #1:: Calloc Fails on pXTX \n");
			     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
		     }
         OnKappaMem = kLen;
         for (int jj = 0; jj < kLen; jj++) {
           pXTX[jj] = (double*) Calloc(kLen, double);
           if (pXTX[jj] == NULL) {
             Rprintf("CoordinateDescentObject #1: Calloc fails inside pXTX[jj=%d]", jj); R_FlushConsole();
             SuccessFlag = -1; return;
           }
         }
		     XTY = (double *) Calloc( kLen + 4, double);  XTYPassByPointer = 0;
		     if (XTY == NULL) {
			     Rprintf("CoordinateDescentObject #1:: Calloc Fails on XTY \n");
			     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
		     }	     
		     pSqMat (pXTX, NLen, kLen, XX);
         F77_CALL(dscal)(&kLen, &ZeroD, XTY, &One);
            F77_CALL(dgemv)("T", &NLen, &kLen,
		        &OneD, XX, &NLen,
		          YY, &One, &ZeroD,
		          XTY, &One);    
         // Rprintf("XTY After Fortran is \n");
         // PrintVector(XTY, kLen);  
        for (ii = 0; ii < kLen; ii++) { XTY[ii] = 0; }
        tMatTimesVec(kLen, NLen,  XX,  YY, XTY);
         //   Rprintf("XTY After Manual is \n");
         //   PrintVector(XTY, kLen);          
   
	    if (PrintFlag > 1) {
	       Rprintf((char*) "CoordinateDescentOjbect #2:PostSquare, XTY is = \n");
	       PrintVector(this->XTY, this->kLen);
	       R_FlushConsole();
      }		     
	    XTResid = (double *) Calloc( kLen, double);
		  if (XTResid == NULL) {
				Rprintf("CoordinateDescentObject #1:: Calloc Fails on XTResid \n");
				R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
			}	
      F77_CALL(dcopy)( &kLen, XTY, &One, XTResid, &One );  			    
		//	for (ii = 0; ii < kLen; ii++) { XTResid[ii] = (double) XTY[ii]; }
  	  InvSXX = (double*) Calloc(kLen+2, double);
	  	if (InvSXX == NULL) {
			     Rprintf("CoordinateDescentObject #1:: Calloc Fails on InvSXX \n");
			     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
		  }
		  for (ii = 0; ii < kLen; ii++) { 
          InvSXX[ii] = ((double) 1.0) /
		                                     ( (double) *(pXTX[ii] + ii)); 
      }       		     
		  	     
	   }  
   }
      
	 if (PrintFlag > 1) {
			Rprintf((char*) "new CoordinateDescentObject #1: Setting up OnBeta\n");
			R_FlushConsole(); R_ProcessEvents();
	 }
      
    OnCoord = 0; OnLoop = 0;
    if (PrintFlag > 1) {
        Rprintf((char*)"CoordinateDescent #1, OnRecordBetas Create?\n");
            R_FlushConsole();
    }        
    OnRecordBetas = NULL;
    RecordBetasPassByPointer = -1;
    if (rRecordOnBetas != NULL) {
	    if (PrintFlag > 1) {
	      Rprintf((char*)"CoordinateDescent #1, rRecordOnBetas is not NULL\n");
	      R_FlushConsole();
	    }
      RecordBetasPassByPointer = 1; 
	    OnRecordBetas = rRecordOnBetas; 	        
		//OnRecordBetas = (double *)Calloc( kLen, double);
		//    	if (OnRecordBetas == NULL) {
		//		     Rprintf("CoordinateDescentObject:: Calloc Fails on OnRecordBetas \n");
	  //		     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
	  //      }	
         }
      
      if (PrintFlag > 1) {
            Rprintf((char*)"CoordinateDescent #1, Ready to go with Gamma = %.4f, NLen = %d, kLen = %d\n",
                      OnGamma, NLen, kLen);
            R_FlushConsole();
      }
      if (PrintFlag > 1) {
	       Rprintf((char*) "CoordinateDescentObject #1:Load, Starting with XTResid = \n");
	       PrintVector(this->XTResid, this->kLen);
	       R_FlushConsole();
      }
        PrintFlag = 0;
        if (ModifiedXTResid == 1) {
	        MakeXTResid();
        }
       if (kLen >= 40000)  {
          int PropMaximum = (int) floor(2147483646.0 / ((2.0)*this->kLen));
          if (PropMaximum < DefaultMaximumAllocation) {
            MaximumAllocation = PropMaximum;
          }
       }
     }	   

////////////////////////////////////////////////////////////////
//  CoordinateDescentObject: Constructor #2
//
//   In this version, input  is only X and Y matrix, this is an adaptable
//    Program where pXTX is created with limiting memory principal.  The inital number
//    Of memory allocated for pXTX is "InitKKs" which should be the anticipated active
//    set size.  Dynamic allocation is to double the memory allocation size every time the
//    barrier is set.
   CoordinateDescentObject(int rNLen, int rkLen, double *rXXX, double *rYYY,
             double *rOnBeta, double rOnGamma, double *rRecordOnBetas, 
             double *rOnGammas, int InitKKs, double InverseGammaConstant, int rPrintFlag) {
     //Rprintf("CDO: XY Input Has anything changed at all?\n\n"); R_FlushConsole();
     //Rprintf("Start create CoordinateDescentObject #2!\n"); R_FlushConsole();
	   PrintFlag = rPrintFlag;
	   if (PrintFlag > 0) {
       Rprintf("CoordinateDescentObject #2: Creating with InitKKs = %d\n", InitKKs);
       R_FlushConsole();
     }
     if (InitKKs <= 0) {
       Rprintf("CoordinateDescentObject #2: InitKKs =%d, that shall not pass goodbye! \n", InitKKs);
       SuccessFlag = -1;
       R_FlushConsole(); return;  FactorGammaDenominator = 1.0;
      
     }
	   XTXPassByPointer = -1; XTYPassByPointer = -1; XXPassByPointer = 0;
	   YYPassByPointer = 0;  RunOrder = NULL;
	   pBeta0 = NULL;  InitKappaMem = InitKKs;
     OnLambda3 = 0;  OnLambdas3 = NULL;  DOnBeta = 0; DOnGammas = 0;
     kLen = rkLen;  NLen = rNLen;  OnEp = 0.0;  DiagShrinkage = 0.0;
     TriggerReallocFlag = 0;  SuccessDOP = 0;  pXTX = NULL;
     TriggerInvestigateFlag = 0;    MemoryFailFlag = 0;
     MaximumAllocation = DefaultMaximumAllocation;   TDFNoise = -1.0; TDFSigma = -1.0;
     //OnLambda2 = 0.0;
   
     if (PrintFlag > 0) {
       Rprintf("CoordinateDescentObject #2: Check InverseGammaconstant \n");
       R_FlushConsole();
     }
	   FactorGammaDenominator = InverseGammaConstant;
	   if (PrintFlag > 0) {
       Rprintf("CoordinateDescentObject #2: We just set FactorGammaDenominator\n");
       R_FlushConsole();
     }
	   if (InitKKs > rkLen) {
       Rprintf("CoordinateDescentObject Start Constructor #2: InitKKs = %d and kLen = %d\n",
          InitKKs, rkLen); R_FlushConsole();
       Rprintf("         Tell Me at Least We Get Here!\n");
          InitKKs = rkLen;    R_FlushConsole();
       SuccessFlag = -1;
       return;
     }
     if (kLen <= 0) {
       Rprintf("CoordinateDescentObject: Constructor #2, no way we take kLen=%d <= 0!\n", kLen);
       R_FlushConsole();  SuccessFlag = -1; return;
     }
     if (InitKKs <= 0) {
       Rprintf("CoordinateDescentObject: StartConstructor #2: InitKKs = %d < 0, not good!\n",
         InitKKs);     R_FlushConsole();
       Rprintf("  kLen is at least %d \n", kLen); R_FlushConsole();
       Rprintf("    We are going to set an error and return!\n"); R_FlushConsole();
       SuccessFlag = -1;
       return;
     }
     if (PrintFlag > 1) {
       Rprintf("CDO #2, InitKappaMem getting set\n"); R_FlushConsole();
     }
	   InitKappaMem = InitKKs;
	   if (PrintFlag > 0) {
		   Rprintf((char*) "Generate CDO #2: pXTXFlag = #2 Forced X,Y Version\n"); R_FlushConsole();
     }
	   pXTX = NULL; XTY=NULL; YY = NULL; XX = NULL; XTResid = NULL;
     XLC = NULL; kXFinder = NULL; OnGammas = NULL; InvSXX = NULL;
     kLen = rkLen; NLen = rNLen;
     iiWeights = NULL;  rWXVec = NULL;  iiWeightsHere = 0;     
     iWeightedXtX = NULL;
     
     if (PrintFlag > 1) {
       Rprintf("CDO #2 for your personal aggrandizement, kLen = %d \n", kLen); R_FlushConsole();
       Rprintf("  But note that InitKKs = %d \n", InitKKs); R_FlushConsole();
     }
     ////////////////////////////////////////////////////////////////////////
     //  Allocate OnBeta, If rOnBeta is non null, inputs OnBeta values
     if (PrintFlag > 1) {
       Rprintf("CDO #2:  Checking OnBeta \n"); R_FlushConsole();
     }
     if (rOnBeta != NULL && rOnBeta[0] != -999) {
       DOnBeta = 1;
       OnBeta = rOnBeta;
       //F77_CALL(dcopy)(&kLen, rOnBeta, &One, OnBeta, &One);
		   ModifiedXTResid = 1; 
       if (PrintFlag > 1) {
         Rprintf("CoordinateDescentObject #2: Using DOnBeta = 1, from supply \n");
         R_FlushConsole();
       }       
     } else {
       OnBeta = (double *) Calloc( kLen+2, double);
       if (OnBeta == NULL) {
		     Rprintf("CoordinateDescentObject #2:: Calloc Fails on OnBeta \n");
		     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
	     }
       int AOne = 1; double ZeroD = 0.0;
       F77_CALL(dscal)(&kLen, &ZeroD, OnBeta, &AOne);	
       DOnBeta = 0;
	   }           
     if (InitKKs <= 0 ) {
       Rprintf((char*) "CoordinateDescentObject #2: InitKKs set to less ");
       Rprintf(" than Zero = %d, false hypothesis\n", InitKKs);
       R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return;
     } else if (InitKKs > kLen) {
        Rprintf((char*)"CoordinateDescentObject Constructor 2: \n"); R_FlushConsole();
        Rprintf("   InitKKs set to more than kLen, no use! \n"); R_FlushConsole();
        Rprintf("  InitKKs = %d \n", InitKKs); R_FlushConsole();
        Rprintf("  kLen = %d \n", kLen); R_FlushConsole();
        Rprintf("  InverseGammaConstant = %d\n", InverseGammaConstant); R_FlushConsole();
        Rprintf((char*)"   InitKKs = %d, kLen = %d, InverseGammaConstant = %f, CDOEpsilon = %f \n",
          InitKKs, kLen, InverseGammaConstant, 0);  R_FlushConsole();
        Rprintf("  Is this okay? \n"); R_FlushConsole();
        Rprintf("    Note, it is rather problematic that this is true \n");
              //R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return; 
      InitKKs = kLen;
    }
      
    if (PrintFlag > 1) {
      Rprintf(" new CoordinateDescentObject #2:: Well setting InitKKs, ");
      Rprintf(" OnKappaS, OnGammas \n"); R_FlushConsole();
    }   
    OnKappaMem = 0;  
    OnKappaS = 0;
    OnGammas = NULL;
    
    if (PrintFlag > 1) {
      Rprintf(" new CoordinateDescentObject #2:: OnGammas Setting\n");
      R_FlushConsole();
    }     
    if (rOnGammas != NULL && rOnGammas[0] >= 0) {
      if (PrintFlag > 1) {
        Rprintf((char*) "new CoordinateDescentObject #2:: rOnGammas != NULL, Setting up Gammas\n");
        R_FlushConsole(); R_ProcessEvents();
      }
      SetUpGammas(rOnGammas);
    } else {
	    if (PrintFlag > 1) {
		    Rprintf((char*)" new CoordinateDescentObject #2:: rOnGammas == NULL or Less than Zero\n");
		      R_FlushConsole(); R_ProcessEvents();
	    }    
    }

	   //int nnii, tt;
	   SuccessFlag = 1;
	   XTXFlag  = 2; 
     OnGamma = rOnGamma; OnLambda = OnGamma * 2;
     kLen = rkLen; NLen = rNLen;
     RecentMove = 0;

	   kLen = kLen;
     if (kLen >= 40000)  {
          int PropMaximum = (int) floor(2147483646.0 / ((2.0)*this->kLen));
          if (PropMaximum < DefaultMaximumAllocation) {
            MaximumAllocation = PropMaximum;
          }
     }      	
	 XXPassByPointer = 1;
   XX = rXXX;     
	 if (PrintFlag > 1) {
	 	 Rprintf((char*) "new CoordinateDescentObject #2: NLen > 0, loading YY\n");
	   R_FlushConsole(); R_ProcessEvents();
	 }
	
   YYPassByPointer = 1;
   YY = rYYY; 
   SetUpToPartial(InitKKs);
                        
   OnCoord = 0; OnLoop = 0;
   if (PrintFlag > 1) {
      Rprintf((char*)"CoordinateDescent #2, OnRecordBetas Create?\n");
      R_FlushConsole();
   }        
   if (PrintFlag > 1) {
     Rprintf((char*)"CoordinateDescent #2, Checking OnRecordBetas\n");
     R_FlushConsole();
   }
   OnRecordBetas = NULL;
   RecordBetasPassByPointer = -1;
    if (rRecordOnBetas != NULL) {
	    if (PrintFlag > 1) {
	      Rprintf((char*)"CoordinateDescent #2, rRecordOnBetas is not NULL\n");
	      R_FlushConsole();
	    } 
	    OnRecordBetas = rRecordOnBetas;
      RecordBetasPassByPointer = 1; 	        
		  //OnRecordBetas = (double *)Calloc( kLen, double);
		  //    	if (OnRecordBetas == NULL) {
		  //		     Rprintf("CoordinateDescentObject:: Calloc Fails on OnRecordBetas \n");
		  //		     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return;
		  //      }	
   }
   if (PrintFlag > 1) {
     Rprintf((char*)"CoordinateDescent #2, Ready to go with Gamma = %.4f, NLen = %d, kLen = %d\n",
       OnGamma, NLen, kLen);
     R_FlushConsole();
   }
   if (ModifiedXTResid == 1) {
	   MakeXTResid();
   }
   if (PrintFlag > 1) {
     Rprintf((char*)"CoordinateDescent #2, Go Home\n");
     R_FlushConsole();
   }
   return;
  }	   
  ~CoordinateDescentObject() {
     PrintFlag = 4;
     if (RecordBetasPassByPointer > 0) {
       OnRecordBetas = NULL; RecordBetasPassByPointer = 1;}
     if (XTXPassByPointer >= 1) {
       Free(pXTX);  pXTX = NULL;
       XTXPassByPointer = 1;
     }
     if (XTYPassByPointer >= 1) {
       XTY = NULL;  XTYPassByPointer = 1;
     }
     if (XXPassByPointer >= 1)  {
       XX = NULL; XXPassByPointer = 1;
     }
     if (YYPassByPointer >= 1) {
       YY = NULL;  YYPassByPointer = 1;
     }
     RunOrder = NULL;
 	   if (DOnBeta > 0) {OnBeta = NULL;} // Do not delete Beta pointed to elsewhere
 	   if (DOnGammas > 0) {OnGammas = NULL;} // Do not Delete OnGammas pointed to elsewhere
     if (iiWeights != NULL) {
       //PrintFlag = 5;    //  iiWeights no longer that important
     }
     if (PrintFlag > 1) {
	     Rprintf((char*) "CDO: SetUpToPartial Freeing pXTX \n"); R_FlushConsole(); R_ProcessEvents();
     }	 	
     FFree(&InvSXX); 	 
	 	 if (XTXPassByPointer <= 0 && pXTX != NULL) { 
      for (int ii = 0; ii < OnKappaMem; ii++) {
        if (pXTX[ii] != NULL) {
          Free(pXTX[ii]); pXTX[ii] = NULL;
        }
      }
      pXTX = NULL;
     }
	 	 FFree(&OnLambdas3);
     if (PrintFlag > 1) {
	     Rprintf((char*) "CDO: SetUpToPartial Freeing XTY \n"); R_FlushConsole(); R_ProcessEvents();
     }	 
  	 if (XTYPassByPointer <= 0 && XTY != NULL) { FFree(&XTY); }
   //  if (PrintFlag > 1) {
	 //    Rprintf((char*) "CDO: SetUpToPartial Freeing XX \n"); R_FlushConsole(); R_ProcessEvents();
   //  }	   	 	 
	 //	FFree(&XX);
   // if (PrintFlag > 1) {
	 //    Rprintf((char*) "CDO: SetUpToPartial Freeing YY \n"); R_FlushConsole(); R_ProcessEvents();
   // }	 
   //	 	 FFree(&YY);
     if (PrintFlag > 1) {
	      Rprintf((char*) "CDO: SetUpToPartial Freeing XTResid \n"); R_FlushConsole(); R_ProcessEvents();
     }	    	 	 
	 	 FFree(&XTResid);
     if (PrintFlag > 1) {
	     Rprintf((char*) "CDO: SetUpToPartial Freeing OnBeta \n"); R_FlushConsole(); R_ProcessEvents();
     }	 	 	
	if (DOnBeta <= 0 && OnBeta != NULL) { FFree(&OnBeta); }
  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing OnRecordBetas \n"); R_FlushConsole(); R_ProcessEvents();
  }	 	    
	if (RecordBetasPassByPointer && OnRecordBetas != NULL) {FFree(&OnRecordBetas);}
  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing InvSXX \n"); R_FlushConsole(); R_ProcessEvents();
  }	 	    
	FFree(&InvSXX);
  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing kXFinder \n"); R_FlushConsole(); R_ProcessEvents();
  }	    
  FFree(&kXFinder);  
  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing XLC \n"); R_FlushConsole(); R_ProcessEvents();
  }	 
  FFree(&XLC); 
  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing OnGammas \n"); R_FlushConsole(); R_ProcessEvents();
  }	                      
  FFree(&OnGammas); 
     
  if (PrintFlag > 1) {
    Rprintf((char*) "CDO: do we free rWXVec?\n"); R_FlushConsole();
  }
  if (rWXVec != NULL) {
    if (PrintFlag > 1) {
      Rprintf("CDO: Yes Free rWXVec. \n"); R_FlushConsole();
    }
    Free(rWXVec);  rWXVec = NULL;
  }
     
  if (PrintFlag > 1) {
    Rprintf("CDO:  Do we Free iiWeights? \n"); R_FlushConsole();
  }
  if (iiWeights != NULL && iiWeightsHere == 1) {
    if (PrintFlag > 1) {
      Rprintf("CDO: Yes Free iiWeights. \n"); R_FlushConsole();
    }
    Free(iiWeights);  iiWeights = NULL;  iiWeightsHere = 0;
  } else if (iiWeights != NULL && iiWeightsHere <= 0) {
    iiWeights = NULL;
  } else {
    if (PrintFlag > 1) {
      Rprintf("CDO: No don't free iiWeights, they are null. \n"); R_FlushConsole();
    }
  }
  if (iWeightedXtX != NULL) {
    if (PrintFlag > 1) {
      Rprintf("CDO: Do we Free iWeightedXtX? \n"); R_FlushConsole();
    }
    Free(iWeightedXtX);
    iWeightedXtX = NULL;
  }
  if( XXPassByPointer <= 0 && XX != NULL) {
    if (PrintFlag > 1) {
      Rprintf("CDO: Yes Free XX, XXPassByPointer = %d. \n", XXPassByPointer); R_FlushConsole();
    } 
    Free(XX); XX = NULL;
  }
   if (PrintFlag > 1) {
      Rprintf("CDO: all we could free is gone. \n", XXPassByPointer); R_FlushConsole();
   }          
 }
 double TDFNoise, TDFSigma;
 int UpdateTDFNoise() {
   // Imagine Epsilon ~ Z / sqrt( sigma^2 * V / TDF) for V sim Inverse Chi-squared.
   //  f(Epsilon|V, sigma) * p(V) propto  exp( - Epsilon_i^2 / (2 sigma^2 V / TDF)) / sqrt( 2 * sigma^2 * V/TDF) * V^(-TDF/2-1) * exp(-1/2V)
   //    Thus V has a posterior given Epsilon  V sim  (TDF + Epsilon^2/sigma^2) * Inverse-ChiSquared(1 + TDF)
   //  The maximum of a inverse chi-squared occurs at d/dV (-TDF/2-1)*log(V) - .5 TDF/V = (-TDF/2 -1)/ V + .5 1/ V^2 = 0
   //   VMax =  1 / ( TDF+2)
   if (TDFNoise < 0 || TDFSigma < 0) {
     Rprintf("UpdateTDFNoise: shouldn't be here, TDFNoise =%f, TDFSigma = %f \n",
       TDFNoise, TDFSigma); R_FlushConsole();
     return(-1);
   }
   if (iiWeights == NULL) {
     Rf_error("UpdateTDFNoise: no, iiWeights is NULL!");
   }
   if (iWeightedXtX == NULL) {
     Rf_error("UpdateTDFNoise: no, iWeightedXtX is NULL!\n");
   }
   if (TDFNoise <= 0  || TDFSigma <= 0) {
     Rf_error("UpdateTDFNoise, no TDFNoise and TDFNoise=%f aren't right!\n", TDFNoise, TDFSigma);
   }
   if (XX == NULL) {
     Rf_error("UpdateTDFNoise, can't do this if XX is NULL!\n");
   }
   if (YY == NULL) {
     Rf_error("UpdateTDFNoise, can't do this if YY is NULL\n");
   }
   for (int ii = 0; ii < NLen; ii++) {
     iiWeights[ii] = 0.0;
   }
   int One = 1;
   if (XTXFlag == 2) {
     for (int ii = 0; ii < OnKappaS; ii++) {
       if (OnBeta[kXFinder[ii]] != 0.0) {
         F77_CALL(daxpy)(&NLen, OnBeta + kXFinder[ii], XX + kXFinder[ii] * NLen, &One, 
           iiWeights, &One);
       }
     }
   } else {
     for (int ii = 0; ii < kLen; ii++) {
       if (OnBeta[ii] != 0.0) {
          F77_CALL(daxpy)(&NLen, OnBeta + ii, XX + ii * NLen, &One, 
           iiWeights, &One);
       }
     }
   }
   for (int ii = 0; ii < NLen; ii++) {
       iiWeights[ii] = YY[ii] - iiWeights[ii];
       iiWeights[ii] = (TDFNoise + 2) / ( TDFNoise + iiWeights[ii] * iiWeights[ii] / TDFSigma);
   }
   return(RefreshWeights(iiWeights));   
 }
 protected :
   double MaxEpsilon;
   int MaxLoops;	
};


class GLMBinomObject: public CoordinateDescentObject  {
  //CoordinateDescentObject CDO;
  private :
  int *y01;
  double *x; double *Z;
  int n, p;
  public :
  double *nlogitCurrentProb;
  double *nProbWeights;
  
  double *pBeta0prev; double *pBetasPrev;
  double *pBeta0Records;
  double *pRecordOnBetas;
  int LengthRecordOnBetas;
  double *nlogitPrevProb;
  int *pIter;  int pIterHere;
  
  public :
  int RunGLMLasso();
  int Sety01Null(){y01 = NULL;return(1);}
  int RefreshCurrentLogitProb();
  int RefreshGLMWeights();
  int PrintFlagGLMB;
  double SumWeights;  double SumYWeights;
  double GLMTotMove;
  
  
  int DataGLMIntegrity(int DoIn);
  
  double CalculateBeta0();
  
  void deletepIter() { 
    if(this->pIter != NULL && pIterHere == 1) {
    Free(this->pIter);  this->pIter=NULL;}}
  GLMBinomObject(int rNLen, int rkLen, double *rX, int *ry01, 
      double *rOnBeta, double *rpBeta0, double *rpBeta0prev, 
      double *rpBetasPrev, double rOnGamma, 
      double *rRecordOnBetas, double *rpBeta0Records,
      double *rOnGammas, int InitKKs, double InverseGammaConstant, 
      double *rnProbWeights, double *rnlogitCurrentProb, double *rnlogitPrevProb,
      double *rZ, int rPrintFlag, double rMaxCDOEpsilon, int rMaxCDORuns,
      int *rpIter, int InMaximumAllocation)
  {
     //  CoordinateDescentObject(rNLen, rkLen, rX, ry01,
     //        rOnBeta, rOnGamma, rRecordOnBetas, 
     //        rOnGammas, InitKKs, InverseGammaConstant)  {
    pIterHere = 1;
	  SuccessFlag = 1;  FactorGammaDenominator = 1.0;
     XTXPassByPointer = -1; XTYPassByPointer = -1; XXPassByPointer = 0;
	   YYPassByPointer = 0;  RunOrder = NULL;
	   pBeta0 = NULL;  InitKappaMem = InitKKs;  //OnLambda2 = 0.0;
     OnLambda3 = 0;  OnLambdas3 = NULL;  DOnBeta = 0; DOnGammas = 0;
     kLen = rkLen;  NLen = rNLen;  OnEp = 0.0;  DiagShrinkage = 0.001;
     TriggerReallocFlag = 0;  SuccessDOP = 0;  pXTX = NULL;
     TriggerInvestigateFlag = 0;    MemoryFailFlag = 0;
     MaximumAllocation = DefaultMaximumAllocation;   TDFNoise = -1.0; TDFSigma = -1.0;
     
    PrintFlagGLMB = rPrintFlag;
    PrintFlag = PrintFlagGLMB - 6; GLMTotMove = 0.0;
    LengthRecordOnBetas = 0;
    if (PrintFlag>0) {
      Rprintf("GLMBinomObject: Loading\n");
      R_FlushConsole();
    }
    if (rpIter == NULL) {
      Rf_error("GLMBinomObject: Error, please supply rpIter!\n");
    }
    if (rpIter != NULL) {
      pIterHere = 0;
      pIter  = rpIter;
    } else {
      pIter = (int*) Calloc(3, int);
    }
    pIter[0] = 0;
    n = rNLen;  NLen = n;
    p = rkLen;  kLen = p;    TDFNoise = -1.0; TDFSigma = -1.0;
    XTXFlag = 2; pXTX=NULL; XTY=NULL; 
    OnKappaMem=InitKKs;InitKappaMem=InitKKs; OnKappaS=0;
    x = rX; XX = rX; y01 = ry01;
    nProbWeights = rnProbWeights; iiWeightsHere   = 0;
    YY =rZ;
    Z = rZ;
    
      XTResid = NULL; InvSXX = NULL; 
      XLC = NULL; kXFinder = NULL; OnGammas = NULL; 
      iiWeights = NULL; iiWeightsHere = 0; rWXVec = NULL;  iWeightedXtX = NULL;
      OnLambdas3 = NULL;  TriggerInvestigateFlag = -1;    
          
    SetEpsilonAndLoops(rMaxCDOEpsilon, rMaxCDORuns); 
    if (rpBeta0Records != NULL && rpBeta0Records[0] >= 0) {
        pBeta0Records = rpBeta0Records;
    } else {
        pBeta0Records = NULL;
    }
    if (rRecordOnBetas != NULL && rRecordOnBetas[0]>= 0) {
        pRecordOnBetas = rRecordOnBetas; OnRecordBetas= NULL;
        LengthRecordOnBetas = (int) fabs(rRecordOnBetas[0]);
    } else {
        pRecordOnBetas = NULL; OnRecordBetas = NULL;
    }
     
    nlogitCurrentProb = rnlogitCurrentProb;   
    pBeta0=rpBeta0;  nlogitPrevProb = rnlogitPrevProb;
    OnBeta = rOnBeta; OnGamma = rOnGamma;
    OnGammas = rOnGammas;  MaximumAllocation = DefaultMaximumAllocation;
    if (InMaximumAllocation >= 100) {
      MaximumAllocation = InMaximumAllocation;
    }
    pBeta0prev = rpBeta0prev;  pBetasPrev =rpBetasPrev;
 
	    XTXPassByPointer = -1; XTYPassByPointer = -1; XXPassByPointer = 1;
	    YYPassByPointer = 1;

	   FactorGammaDenominator = InverseGammaConstant;
	   if (InitKappaMem > rkLen) {
       Rprintf("CoordinateDescentObject Start Constructor #2:\n");
       Rprintf(" InitKKs = %d but kLen = %d\n",
          InitKKs, rkLen); R_FlushConsole();
       Rprintf("         Tell Me at Least We Get Here!\n");
          InitKKs = rkLen;
     }
     
	   if (InitKappaMem > kLen) {InitKappaMem=kLen;}
	   if (InitKappaMem <= 0) {InitKappaMem=2;}
	   OnKappaMem = InitKappaMem;
	   
	   if (PrintFlagGLMB > 0) {
		   Rprintf((char*) "Generate CDO: XTXFlag = #2 Forced Version\n");
     }
      
      
     InvSXX = (double*) Calloc(p+1, double); 
     OnBeta = rOnBeta;    
     ModifiedXTResid=1;     
         
     OnGammas = rOnGammas;
     OnGamma = rOnGamma; OnLambda = OnGamma * 2;
    
     RecentMove = 0;
  	 if (PrintFlagGLMB > 1) {
  		     Rprintf((char*) "new GLMB: Set up Partial Now \n");
  		     R_FlushConsole(); R_ProcessEvents();
  	 }

    if (PrintFlagGLMB > 1) {
       Rprintf(" new GLMB: here are the current ProbWeights: ");
       PrintVector(nProbWeights, kLen);
       Rprintf("\n"); R_FlushConsole();
     }
     if (PrintFlagGLMB >= -4) {
       Rprintf("  new GLMB: At this point we are setting up to Partial. \n");
       if (XTY == NULL) {
         Rprintf("         Note: XTY is null, good. \n"); R_FlushConsole();
       } else {
         Rprintf("         Note: XTY is not null, bad!\n");
       }
       if (pXTX == NULL) {
         Rprintf("         Note: pXTX is null, good. \n"); R_FlushConsole();
       } else {
         Rprintf("         Note: pXTX is not null, bad!\n");
       }
     } 
     SetUpToPartial(InitKappaMem);
     if (PrintFlagGLMB >= -4) {
       Rprintf("  new GLMB: At this point we done setting up Partial. \n");
       if (XTY == NULL) {
         Rf_error("         Note: XTY is null, VERY BAD. \n"); R_FlushConsole();
       } else {
         Rprintf("         Note: XTY is not null, it is: \n");  R_FlushConsole();
         PrintVector(XTY, kLen); Rprintf("\n");   R_FlushConsole();
       }
       if (pXTX == NULL) {
         Rf_error("         Note: pXTX is NULL, VERY BAD. \n"); R_FlushConsole();
       } else {
         Rprintf("         Note: pXTX is not null, its first row is: \n"); R_FlushConsole();
         PrintVector(pXTX[0], kLen); Rprintf("\n");
         Rprintf("         Note OnKappaS = %d, kXFinder is: \n", OnKappaS);
         PrintVector(kXFinder, OnKappaS);
         Rprintf("          Is this good or bad? PrintFlagGLMB = %d \n", PrintFlagGLMB); R_FlushConsole();
         Rprintf("   --     On To More Setup Action. \n"); R_FlushConsole();
       }
     }                         
     OnCoord = 0; OnLoop = 0;
   
 
     if (PrintFlagGLMB > 1) {
        Rprintf((char*)"new GLMB, Ready to MakeXTResid, ");
        Rprintf("  NLen = %d, kLen = %d\n",
                         NLen, kLen);
        R_FlushConsole();
     }
     if (ModifiedXTResid >= 1) {
  	        MakeXTResid();
     }
     
     if (PrintFlagGLMB >1) {
       Rprintf(" new GLMB: going to Setup Weights Now\n");
       R_FlushConsole();
     } 

     SetupWeights(nProbWeights);  iiWeightsHere = 0;
     if (PrintFlagGLMB >= -4) {
       Rprintf("  new GLMB: At this point we done setting up nProbWeights. \n");
       if (XTY == NULL) {
         Rf_error("         Note: XTY is null, VERY BAD. \n"); R_FlushConsole();
       } else {
         Rprintf("         Note: XTY is not null, it is: \n");  R_FlushConsole();
         PrintVector(XTY, kLen); Rprintf("\n");   R_FlushConsole();
       }
       if (pXTX == NULL) {
         Rf_error("         Note: pXTX is NULL, VERY BAD. \n"); R_FlushConsole();
       } else {
         Rprintf("         Note: pXTX is not null, its first row is: \n"); R_FlushConsole();
         PrintVector(pXTX[0], kLen); Rprintf("\n");
         Rprintf("         Note OnKappaS = %d, kXFinder is: \n", OnKappaS);
         PrintVector(kXFinder, OnKappaS);
         Rprintf("          Is this good or bad? \n"); R_FlushConsole();
         
       }
     } 
     if (PrintFlagGLMB > 1) {
       Rprintf(" new GLMB: here are the XTY after the startup: ");
       PrintVector(XTY, kLen);
       Rprintf("\n"); R_FlushConsole();
     }
          
     if (PrintFlagGLMB > 1) {
        Rprintf("GLMBinomObject: FinishedConstructor\n");
        R_FlushConsole();
     }
      
  }
   ~GLMBinomObject() {
   
     PrintFlagGLMB = 3;
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: Calling Destructor\n"); R_FlushConsole();
     }
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: Freeing InvSXX\n"); R_FlushConsole();
     }
     if (InvSXX != NULL) {
        Free(InvSXX); InvSXX = NULL;
     }
     if (PrintFlagGLMB > 2 && pIter != NULL && pIterHere == 1) {
       Rprintf("GLMBinomObject: Freeing nonNull pIter\n"); R_FlushConsole();
     }
     if (pIter != NULL && pIterHere == 1) {
       Free(pIter); pIter = NULL;
     }
     if (PrintFlagGLMB > 2 && XTResid != NULL) {
       Rprintf("GLMBinomObject: Freeing NonNull XTResid\n"); R_FlushConsole();
     }
     if (XTResid != NULL) {
       if (PrintFlagGLMB > 2) {
        Rprintf("GLMBinomObject: Freeing XTResid\n"); R_FlushConsole();
       }
        Free(XTResid); XTResid = NULL;
     }
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: Decide to free pXTX\n"); R_FlushConsole();
     }
     if (pXTX != NULL) {
       if (OnKappaMem >= 1) {
       if (PrintFlagGLMB > 2) {
         Rprintf("GLMBinomObject: Will Have to Free pXTX\n"); R_FlushConsole();
       }
       for (int ii = 0; ii < OnKappaMem; ii++) {
         if (pXTX[ii] != NULL) {
         for (int jj = 0; jj < kLen; jj++) {
           *(pXTX[ii]+jj) = 0.0;
         }
         }
       }
        for (int ii = 0; ii < OnKappaMem; ii++) {
          if (pXTX[ii] != NULL) {
            if (PrintFlagGLMB > 4) {
              if (ii % 10 == 0) { Rprintf("\n"); R_FlushConsole(); }
              Rprintf("Free(%d=%d):", ii, kXFinder[ii]); R_FlushConsole();
            }
            Free(pXTX[ii]);  pXTX[ii] = NULL;
          }
        }
       }
       Free(pXTX); pXTX=NULL;
     }
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: Freeing XTY\n"); R_FlushConsole();
     }     
     if (XTY != NULL) {
        Free(XTY); XTY=NULL;
     }
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: Freeing XLC\n"); R_FlushConsole();
     }
     if (XLC != NULL) {
        Free(XLC); XLC = NULL;
     }
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: Freeing kXFinder\n"); R_FlushConsole();
     }
     FFree(&kXFinder);
     if (PrintFlagGLMB > 2) {
       Rprintf("~GLMBinomObject: Freeing rWXVec\n");  R_FlushConsole();
     }
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: Freeing WXVec\n"); R_FlushConsole();
     }
     FFree(&rWXVec);
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: Setting Things to to destroy to NULL\n"); R_FlushConsole();
     }
     XX = NULL;
     YY = NULL;
     if (iWeightedXtX != NULL) {
       if (PrintFlagGLMB >= -1) {
         Rprintf("GLMBinomObject: Freeing iWeightedXtX \n"); R_FlushConsole();
       }
       Free(iWeightedXtX); iWeightedXtX = NULL;
     }
     iiWeights = NULL; iiWeightsHere = 0; iWeightedXtX = NULL;
     pBeta0=NULL; pBeta0Records=NULL;
     OnBeta=NULL; OnGammas=NULL;
     pBetasPrev=NULL;  pBeta0prev=NULL; OnRecordBetas=NULL; pRecordOnBetas=NULL;    
     if (PrintFlagGLMB > 2) {
       Rprintf("GLMBinomObject: End OF Free\n"); R_FlushConsole();
     }
   }
   //int RefreshCurrentLogitProb();
   //int RefreshGLMWeights();
  
   int SetBeta0(double InBeta0) {
     pBeta0[0] = 0.0;
     return(1);
   }
   double get_GLMTotMove() { return(GLMTotMove); }
   SEXP get_Z() {
     if (Z == NULL) {
       Rf_error("GLMBinomObject: No Z. \n"); R_FlushConsole();
     }
     SEXP sOut = R_NilValue; Rf_protect(sOut = Rf_allocVector(REALSXP, n));
     for (int ii = 0; ii < n; ii++) {REAL(sOut)[ii] = Z[ii]; }
     Rf_unprotect(1); return(sOut);
   }
   SEXP get_y01() {
     if (y01 == NULL) {
       Rf_error("GLMBinomObject: No y01. \n"); R_FlushConsole();
     }
     SEXP sOut = R_NilValue; Rf_protect(sOut = Rf_allocVector(INTSXP, n));
     for (int ii = 0; ii < n; ii++) {INTEGER(sOut)[ii] = y01[ii]; }
     Rf_unprotect(1); return(sOut);
   }
  // CoordinateDescentObject(int rNLen, int rkLen, double *rXXX, double *rYYY,
  //           double *rOnBeta, double rOnGamma, double *rRecordOnBetas, 
  //           double *rOnGammas, int InitKKs, double InverseGammaConstant
};
   





