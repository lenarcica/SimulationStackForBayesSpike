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
#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif

///  DEFCONST.  Currently the LARS algorithm is not allowed
///     to be run more than 2.5 times the kLen, or number of covariates.
# ifndef DEFCONST
   #define DEFCONST 3.5
#endif

int LarsLassoAlgorithmBI( double * yys, int NLen, double *xxs, int kLen, 
                      double *InitBetas, 
                      double *OldBetas, double lambda,double *OutBetas, double *BetasKeep,
                      int *JKeep, double *PenaltyKeep,
                      double *GammasKeep, double *gFKeep,
                      int GFlag, 
                      int WeightLen, double *Weights, int *pFinaltt);
                      
double CalculateGFTAll(double *TFLV2, double *GFT0, double *GFT1, int OnL, int TotK,
          int *ActiveOns, double *OnBetas, double *BetaOlds, int NLen, double *OnResid, 
          double *wa, int *sjA, double *ua, double v1, double v2, double OverallMax);                      
//////////////////////////////////////////////////////////
//   LARSObject
//
//     The containing C++ object for the functions, datastorage
//         and algorithm running our implementation of LARS
//
class LARSObject      {
   public: 	
    int GunFlag;
    int PrintFlag;
    
    int NLen; int kLen; //NLen, length of Y, kLen, # of covariates
        
	  double * yys; //
    double *xxs;  // xxs, double vector of xxs
    	  
    double SumYYSq;
    double SumCurResids;

    
    double *InitBetas;  // Insert Betas
    double *OldBetas;
    double *OutBetas;   // Output Betas
        
    double *BetasKeep;  // matrix of Betas returned by algorithm
    
    double lambda;  // lambda Criterion to stop
    

    int *JKeep;          // Record of coordinates as they enter model
    double *PenaltyKeep; // Matrix of Returned Penalties at each step
    double *GammasKeep;  // Matrix of Returned Gammas
    double *gFKeep;      // Matrix of returned Gamma
    
    double *GAll;             //GAll is XTX matrix
    double *XTYAll, *GAllOnB; // XTYAll is XTYAll for all coeffiicents
    
    int GFlag, FakeFlag;
    
    double *WAll;
    double *WDiag; double * WDiagBeta;
    double *WIV;
    
    double AvWeights;
    int WDiagFlag;
    double *BetaOlds;
    double *wa;
    //double *Wwa;
    double *aa;
    double *OnBetas;
    double *ActiveBetas;
    double *PrevBetas;
    double *gammaList1;
    double *gammaList2;
    double *ua;
    double *cc;
    double *OnMuA;
    double *PrevMuA;
    double *OnResid;
    double *PrevResid;
    double *XXA;
    double *GGAA;
    double *IGGA;
    int MaxIter;
    double CurSumYTXAwa;
    double CurSumOnBetaXTXAwa;
    double CurSumwaXATXAwa;
    int L2Involved;
    int *ActiveOns;
    int *AAInActive;
    double *CurMaxs;
    int *ActiveSigns;
    int *OnActives;
    //int *OldActivesOn, OldOnL;
    double EpsilonBM;
    int MaxIters;
    double GFT0, GFT1, TFLV2;
    double GFTAll, FFeatConst; 
    int maxish;
    int Ontt;
    double CMaxNow;
    double  RunningMin; 
    double CurrentPenalty;
    double PrevPenalty;
    int OnInL, OnL, PrevOnL;
      double AAA;
      int CurrentGoat;
      double CalcAAA();
      int CreateWa();
      int GGAPrint();
      int PrintWeights();
      double LARSgammaLists();
      int CalculateGFTAll(double OverallMax);
      int GenMecc(int jj);
      int AssessCC();
      int ElseCC(int PrevGoal, int GotGoal);
      int MakeGGA();
      int Makeaa();
      int GAllPrint();
      int CalcWithGoat(double FeatConst);
      int FeedFreeBetas(int tt, double FeatConst);
      int LarsLassoAlgorithm2();
          int LARSInit ();
      int CreateWAll(int WeightLen, double *Weights);
      int CreateGAll();
      int CalculateYTXwScores();
      int CalcGAwa();
      
      int SetMeI ( int kSLen, double *GMat ); // Sets up Identity Matrix
      //int SetMeXIn (int kLen, int NLen, int OnL, 
      //        int* ActiveOns, int *ActiveSigns, double *XInputMat , 
      //        double *XOutputMat); // Resigns Matrix according to signs
      int SetMeXIn (int kLen, int NLen, int OnL, int* ActiveOns,  int*ActiveSigns,
              double *XInputMat , 
              double *XOutputMat); // Resigns Matrix according to signs
       int SetMeXIn (int kLen, int NLen, int OnL, 
              int* ActiveOns,  int*ActiveSigns,
              long double *XInputMat , 
              long double *XOutputMat);              
      LARSObject(double * ryys, int rNLen, double *rxxs, int rkLen, 
                      double *rInitBetas, 
                      double rlambda) {
	      FakeFlag = 0;
	      Ontt = 0;
	      GunFlag = 1;
	      //PrintFlag = 4;
	      //PrintFlag = 1;
	      //PrintFlag = 1;
	      PrintFlag = 0;
	      if (PrintFlag > 1) {
		      Rprintf((char*)"LARSObject::  Creating LarsObject NLen = %d, kLen = %d\n", rNLen, rkLen);
		      R_FlushConsole();
		      R_ProcessEvents();
	      }
	      CurrentGoat = -20;
	      L2Involved = 0;
	      GFlag = 0;
	      XTYAll = NULL;
	      GAll = NULL;
	      WAll = NULL;
	      WDiag = NULL; WDiagBeta = NULL;
	      GAllOnB = NULL;
	      WDiagFlag = 0;
	      PrevOnL = 0;
	      EpsilonBM = .0000002;
	      CurrentPenalty = -1.0;
	      PrevPenalty = -1.0;
          int ii,jj,tt;
          WIV = NULL;
          OnResid = NULL; PrevResid = NULL; OnMuA = NULL;PrevMuA = NULL;
	      //int OldOnL = 0;
	      //OldActivesOn = NULL;
	      //int OldActivesOn = (int *) Calloc(rkLen+1, int);
	      lambda = rlambda;	      
	      if (rNLen < 0) {
		     if (PrintFlag >= 1) {
			     Rprintf("LarsObject:: Doing rNLen < 0 Creation\n");
			     R_FlushConsole();
		     }
		     NLen = -rNLen;
		     kLen = rkLen;
		     FakeFlag = -1; GFlag = 1;
		     yys = NULL;
		     xxs = NULL;
		     GAll = NULL;
		     GAll = (double *) Calloc( rkLen *  rkLen + 4,double);
				    	   if (GAll == NULL) {Rprintf((char*)"LarsObject Cannot Apportion GAll \n");
					                       R_FlushConsole();
					                       GunFlag = -1;
					                       return;
				           }	
				jj = kLen * kLen;
				for (ii = 0; ii < jj; ii++) {
					GAll[ii] = (double) rxxs[ii];
		        }
		      XTYAll = (double *) Calloc(kLen + 3, double);
		                  if (XTYAll == NULL) {Rprintf((char*)"LarsObject Cannot Apportion XTYAll\n");
		                       R_FlushConsole(); GunFlag = -1;
		                       return;
	                      }
	           for (ii = 0; ii < kLen; ii++) {
		            XTYAll[ii]  = (double) ryys[ii];
	           }
		       GAllOnB = (double *) Calloc(kLen + 3, double);
		                  if (GAllOnB == NULL) {Rprintf((char*)
		                     "LarsObject Cannot Apportion GAllOnB\n");
		                       R_FlushConsole(); GunFlag = -1;
		                       return;
	                      }	           
		     if (PrintFlag >= 1) {
			     Rprintf("LarsObject:: Finished with xxs specific Creation\n");
			     R_FlushConsole();
		     }	           
	           SumYYSq = GAll[0];
	            OnMuA = NULL;
	            PrevMuA = NULL; PrevResid = NULL;
	            OnResid = NULL; OnMuA = NULL;
	            XXA = NULL;
	      } else {
		     FakeFlag = 1;
	 	      NLen = rNLen;
		      kLen = rkLen; 
		          XTYAll = NULL; GAll = NULL;   
			      //Wwa = NULL;
		          yys = NULL;	      
			      yys = (double *) Calloc( rNLen+2,double);
				    	   if (yys == NULL) {Rprintf((char*)"LarsObject Cannot Apportion yys \n");
					                       R_FlushConsole();
					                       GunFlag = -1;
					                       return;
				           }	
		          xxs = NULL;		                 
			      xxs = (double *) Calloc( (rNLen+1) * (rkLen+1), double);
				    	   if (xxs == NULL) {Rprintf((char*)"LarsObject Cannot Apportion xxs \n");
					                       R_FlushConsole(); GunFlag = -1;
					                       return;	    
		                   }
		      tt = NLen * kLen;
		      for (ii = 0; ii < tt; ii++) { xxs[ii] = (double) rxxs[ii]; }		                   
			  SumYYSq = 0;      
		      for (ii = 0; ii < NLen; ii++) { yys[ii] = (double) ryys[ii]; 
		           SumYYSq += yys[ii] * yys[ii];}	
		      SumCurResids = SumYYSq;	
            OnMuA = NULL;		            
		    OnMuA  = (double *)Calloc(NLen+5, double); // Current MuA
		           if (OnMuA == NULL) {Rprintf((char*)"LarsObject Cannot Apportion OnMuA \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }
		    PrevMuA = NULL; 		    
		    PrevMuA = (double *)Calloc(NLen+5, double);  // Prev MuA
		           if (PrevMuA == NULL) {Rprintf((char*)"LarsObject Cannot Apportion PrevMuA \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           } 		    
		    OnResid = NULL;
		    OnResid = (double *) Calloc(NLen+5, double);
		    	   if (OnResid == NULL) {Rprintf((char*)"LarsObject Cannot Apportion OnResid \n");
			                       R_FlushConsole();
			                       return;
		           }		    
		    PrevResid=NULL;
		    PrevResid = (double *)Calloc(NLen+5, double);
		    	   if (PrevResid == NULL) {Rprintf((char*)"LarsObject Cannot Apportion PrevResid \n");
			                       R_FlushConsole();
			                       return;
		           }
		        XXA = NULL;
			    XXA = (double  *) Calloc( (NLen+1) * (kLen+4)+kLen + NLen, double);
			    if (XXA == NULL) {Rprintf((char*)"LarsObject Cannot Aportion XXA\n");
			           R_FlushConsole(); GunFlag = -1;
			           return;
		        }
	        }



          MaxIters = kLen;
            if(PrintFlag > 2) {
	            Rprintf((char*)"LARSObject: Creating Memory for kLen Items, NLen = %d, kLen=%d\n", 
	                 NLen, kLen);
	            R_FlushConsole();
	            R_ProcessEvents();
            }
            PenaltyKeep = NULL;
            PenaltyKeep = (double *)Calloc((((int) ceil(kLen*DEFCONST))+2) * 4, double);
		    	   if (PenaltyKeep == NULL) {Rprintf((char*)"LarsObject Cannot Apportion PenaltyKeep \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }
		    GammasKeep = NULL;            
            GammasKeep = (double *)Calloc((((int) ceil(kLen*DEFCONST))+3), double);
 		    	   if (GammasKeep == NULL) {Rprintf((char*)"LarsObject Cannot Apportion GammasKeep \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }  
		    gFKeep = NULL;         
            gFKeep = (double*)Calloc((((int) ceil(kLen*DEFCONST))+2), double);
		    	   if (gFKeep == NULL) {Rprintf((char*)"LarsObject Cannot Apportion gFKeep \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           } 
		    BetasKeep = NULL;             
            BetasKeep = (double *)Calloc((((int) ceil(kLen*DEFCONST))+2) * (kLen+1), double);
		    	   if (BetasKeep == NULL) {Rprintf((char*)"LarsObject Cannot Apportion BetasKeep \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }              
            JKeep = NULL;
            JKeep = (int*) Calloc(((int) ceil(kLen*DEFCONST))+3, int);
		    	   if (JKeep == NULL) {Rprintf((char*)"LarsObject Cannot Apportion JKeep \n");
			                       R_FlushConsole();GunFlag = -1;
			                       return;
		           } 
		    OnActives=NULL;                 
		    OnActives = (int *) Calloc( kLen+5, int);
		    	   if (OnActives == NULL) {Rprintf((char*)"LarsObject Cannot Apportion OnActives \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           } 
		    BetaOlds=NULL; 		    
		    BetaOlds = (double *) Calloc(kLen+5, double);
		    	   if (BetaOlds == NULL) {Rprintf((char*)"LarsObject Cannot Apportion BetaOlds \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }  		    
		    InitBetas = NULL;
		    InitBetas = (double *)Calloc(kLen+3, double);
                   if (GammasKeep == NULL) {Rprintf((char*)"LarsObject Cannot Apportion InitBetas \n");
			                       R_FlushConsole();GunFlag = -1;
			                       return;
		           } 
		        OnBetas = NULL;  		        
		        OnBetas = (double *)Calloc(kLen+5, double);
		    	   if (GammasKeep == NULL) {Rprintf((char*)"LarsObject Cannot Apportion GammasKeep \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           } 
		        OutBetas = NULL; 		        
		        OutBetas = (double *) Calloc(kLen + 5, double);
		           if (OutBetas == NULL) {Rprintf((char*)"LarsObject Cannot Apportion OutBetas \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           } 
		        PrevBetas = NULL;
		        PrevBetas = (double *)Calloc(kLen+5, double);
		           if (PrevBetas == NULL) {Rprintf((char*)"LarsObject Cannot Apportion PrevBetas \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }		           
		    if (rInitBetas == NULL) {
			    for (ii = 0; ii < kLen; ii++) {
				    InitBetas[ii] = 0;
				    BetaOlds[ii] = 0;
				    OnBetas[ii] = 0; 
				    OutBetas[ii] = 0; 		
				    PrevBetas[ii] = 0; 				    		    
			    } 
		    } else {
			    for (ii = 0; ii < kLen; ii++) {
				    InitBetas[ii] = (double) rInitBetas[ii];
				    BetaOlds[ii] = (double) rInitBetas[ii];
				    OnBetas[ii] = (double) rInitBetas[ii];
				    OutBetas[ii] = (double) rInitBetas[ii];	
				    PrevBetas[ii] = (double) rInitBetas[ii];			    
			    }
		    }
		     		    
            if(PrintFlag > 2) {
	            Rprintf((char*)"LARSObject: Copying in Initial Betas\n");
	            R_FlushConsole();
	            R_ProcessEvents();
	            Rprintf((char*)"LARSObject: rInitBetas[0] = %.6f, rInitBetas[kLen-1] = %.6f\n",
	                                (double)  rInitBetas[0], double(rInitBetas[kLen-1]));
	            R_FlushConsole();
	            R_ProcessEvents();	            
            }		    
		    for (ii =0;ii < kLen; ii++) {
			      InitBetas[ii] = (double) rInitBetas[ii];
			      BetaOlds[ii] = (double) rInitBetas[ii];
			}
            if(PrintFlag > 2) {
	            Rprintf((char*)"Creating Memory for wa, aa, OnBetas, glists!\n");
	            R_FlushConsole();
	            R_ProcessEvents();
            }		
                wa = NULL;	
		        wa = (double *) Calloc(kLen+5, double);
		    	   if (wa == NULL) {Rprintf((char*)"LarsObject Cannot Apportion wa \n");
			                       R_FlushConsole();GunFlag = -1;
			                       return;
		           }
		        aa = NULL;		        
		        aa = (double *) Calloc(kLen+5, double);
		    	   if (aa == NULL) {Rprintf((char*)"LarsObject Cannot Apportion aa \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }

		        gammaList1 = NULL;		        
		        gammaList1 = (double *)Calloc(kLen+5, double);
		           if (gammaList1 == NULL) {Rprintf((char*)"LarsObject Cannot Apportion gammaList1 \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }
		        gammaList2 = NULL; 		        
		        gammaList2 = (double *) Calloc(kLen+5, double);
		           if (gammaList2 == NULL) {Rprintf((char*)"LarsObject Cannot Apportion gammaList2 \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           } 		        
			    for (ii = 0; ii < kLen; ii++) {
			        OnActives[ii] = 0;
			        wa[ii] = 0;
			        aa[ii] = 0;
			        OnBetas[ii] = InitBetas[ii];
			        PrevBetas[ii] = 0;
			        BetasKeep[ii] = InitBetas[ii];
			        GammasKeep[ii] = 0;
			        gFKeep[ii] = 0;
			    }
           maxish = kLen;
           if (kLen < NLen) {maxish = NLen;}
            if(PrintFlag > 2) {
	            Rprintf((char*)"Creating Memory for ua, cc, PrevMuA\n");
	            R_FlushConsole();
	            R_ProcessEvents();
            } 
            ua = NULL;          
		    ua = (double *) Calloc (maxish+2, double);
		           if (ua == NULL) {Rprintf((char*)"LarsObject Cannot Apportion ua \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           } 
		    cc = NULL;		    
		    cc = (double *) Calloc( maxish+2, double);  // Current Residual Correlations
		    		if (cc == NULL) {Rprintf((char*)"LarsObject Cannot Apportion cc \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }
            MaxIter = (int) ((int) kLen* ((int) DEFCONST) );
		    OldBetas = NULL;
    //Rprintf( (char *) "LarsLasso: Data Loaded in \n");
            CurSumYTXAwa = 0;
            CurSumOnBetaXTXAwa = 0;
            CurSumwaXATXAwa = 0;
            if (FakeFlag >= 0) {
			    if(PrintFlag > 2) {
		            Rprintf((char*)"OnMuA copying now\n");
		            R_FlushConsole();
		            R_ProcessEvents();
		            Rprintf((char*)"OnMuA:  yys[0] = %.5f, xxs[0] = %.6f, OnBetas[0] = %.6f, OnMuA[0] = %.6f\n",
		                (double) ryys[0], (double) rxxs[0],  (double) this->OnBetas[0], (double) OnMuA[0]);
		            R_FlushConsole();
		            R_ProcessEvents();
	            }
			    for (ii = 0; ii < NLen; ii++) { OnMuA[ii] = 0;}
			    int TotF = 0;
			    for (jj = 0; jj < kLen; jj++) {
				    for (ii = 0; ii < NLen; ii++) { 
					      OnMuA[ii] += (double) xxs[TotF] * (double) OnBetas[jj];
					      TotF++;
				    }
			    }
			    for (ii = 0; ii < NLen; ii++) {
			        //OnMuA[ii] = 0;
			        //PrevMuA[ii] = 0;
			        ua[ii] = 0;
			        OnResid[ii] = (double) yys[ii] - (double) OnMuA[ii];
			        PrevResid[ii] = (double) yys[ii] - (double) OnMuA[ii];
			    } 			    
	        }
		    if (PrintFlag > 2) {
			    Rprintf((char*)"Freshening up cc, ua, OnResid, PrevResid\n");
			    R_FlushConsole();
			    R_ProcessEvents();
	     	}
	     	if (FakeFlag >= 0 && PrintFlag > 3) {
		     	Rprintf( (char*) "Printing Initial OnMuA \n");
		     	for (ii = 0; ii < NLen; ii++) {
			     	Rprintf((char*) "      %d=%.6f\n", ii, (double) OnMuA[ii]);
		        }
		        R_FlushConsole();
	        }
		    for (ii = 0; ii < kLen; ii++) {  
			    cc[ii] = 0;
		    }

		    if (FakeFlag >= 0 && PrintFlag > 3) {
		     	Rprintf( (char*) "Printing Initial OnResid -- PrevResid \n");
		     	for (ii = 0; ii < NLen; ii++) {
			     	Rprintf((char*) "      %d=%.6f=%.6f\n", 
			     	      ii, (double) OnResid[ii], (double) PrevResid[ii]);
		        }
		        R_FlushConsole();
	        }

            CurrentPenalty = -1.0;
            PrevPenalty = -1.0; 

            if(PrintFlag > 2) {
	            Rprintf((char*)"Generating XXA, GGAA mem\n");
	            R_FlushConsole();
	            R_ProcessEvents();
            }

		        GGAA = NULL;
			    GGAA = (double *) Calloc( (kLen+2) * (kLen+2) + kLen + 2, double);
			    if (GGAA == NULL) {Rprintf((char*)"LarsObject Cannot Apportion GGAA \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		                          }
		        IGGA = NULL;
			    IGGA = (double *) Calloc((kLen+2) * (kLen+2) + kLen+1, double);   
			    	if (IGGA == NULL) {Rprintf((char*)"LarsObject Cannot Apportion IGGA \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }
		        ActiveOns = NULL;
			    ActiveOns = (int *)Calloc(kLen+5, int);
			    if (ActiveOns == NULL) {Rprintf((char*)"LarsObject Cannot Apportion ActiveOns \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }
		        AAInActive = NULL;
			    AAInActive = (int*)Calloc(kLen+5, int);
			    if (AAInActive == NULL) {Rprintf((char*)"LarsObject Cannot Apportion AAInActive \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }	
		        CurMaxs = NULL;		    
			    CurMaxs = (double *)Calloc(kLen+5, double); 
			    if (CurMaxs == NULL) {Rprintf((char*)"LarsObject Cannot Apportion CurMaxs \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }	
		        ActiveSigns = NULL;		    
			    ActiveSigns = (int *)Calloc(kLen+5, int);
			    if (ActiveSigns == NULL) {Rprintf((char*)"LarsObject Cannot Apportion ActiveSigns \n");
			                       R_FlushConsole(); GunFlag = -1;
			                       return;
		           }			    
			  
              if(PrintFlag > 2) {
	            Rprintf((char*)"Calculating current penalty\n");
	            R_FlushConsole();
	            R_ProcessEvents();
            }
			    PrevPenalty = 0;
			    for (tt = 0; tt < (kLen +2) * 4; tt++) {
				       PenaltyKeep[tt] = 0;
			    }
			    if (FakeFlag >= 0) {
				    for (tt = 0; tt < NLen; tt++) { 
					       PenaltyKeep[0] += (OnResid[tt]) * (OnResid[tt]);
				    }
			    } else {
				    PenaltyKeep[0] = SumYYSq;
			    }
				//       PrevPenalty += (OnResid[tt]) * (OnResid[tt]); 
			      for (tt = 0; tt < kLen; tt++) { 
				      PenaltyKeep[1] += fabs( OnBetas[tt] );
				      PenaltyKeep[2] += fabs( OnBetas[tt] - BetaOlds[tt] );
			      }
			      PenaltyKeep[1] = PenaltyKeep[1] * lambda;
			      PenaltyKeep[2] = PenaltyKeep[2] * 0.0;
			      PenaltyKeep[3] = PenaltyKeep[0] + PenaltyKeep[1] + PenaltyKeep[2];
			      PrevPenalty = PenaltyKeep[3];  
			  if (PrintFlag >= 1) {
				  Rprintf("LarsObject:: Finished Loading \n");
				  R_FlushConsole();
	         }     
    }
    int RefreshLARS();
    ~LARSObject() {
       //PrintFlag = 4;
       //PrintFlag = 4;
       if (PrintFlag > 1) {
	       Rprintf((char*) "LarsObject: Freeing, Cleaning, GunFlag = %d\n", GunFlag);
	       R_FlushConsole();
       }
       
       if (PrintFlag > 2) {
	       Rprintf((char*) "LarsObject: Freeing OldBetas\n", GunFlag);
	       R_FlushConsole();
       }	    	    
       if (OldBetas!= NULL) {Free(OldBetas);OldBetas=NULL;}
       if (PrintFlag > 2) {
	       Rprintf((char*) "LarsObject: Freeing xxs\n", GunFlag);
	       R_FlushConsole();
       } 
       if (GAllOnB != NULL) {Free(GAllOnB); GAllOnB=NULL;}      
       if (xxs != NULL) {Free(xxs);xxs = NULL; }
       if (PrintFlag > 2) {
	       Rprintf((char*) "LarsObject: Freeing yys\n", GunFlag);
	       R_FlushConsole();
       }        
       if (yys!= NULL) {Free(yys);yys= NULL;}
      
       if (PrintFlag > 1) {
	       Rprintf((char*) "Freeing GAll, XTYAll, BetaOlds \n");
	       R_FlushConsole();
       }
      // if (OldActivesOn != NULL) {Free(OldActivesOn); OldActivesOn = NULL;}		
       if (GAll != NULL) {Free(GAll);GAll = NULL;} 
       if (XTYAll != NULL) {Free(XTYAll); XTYAll = NULL;}
       if (WDiag != NULL) {Free(WDiag); WDiag = NULL;}
       if (WDiagBeta != NULL) {Free(WDiagBeta); WDiagBeta = NULL; }
       //if (Wwa != NULL) {Free(Wwa); Wwa = NULL; }
       if (WAll != NULL) {Free(WAll); WAll = NULL; }   
	   if(BetaOlds!= NULL) {Free(BetaOlds ); BetaOlds = NULL; }
		    if(InitBetas!=NULL) {Free(InitBetas); InitBetas=NULL;}
		    if(PenaltyKeep!=NULL){Free(PenaltyKeep); PenaltyKeep=NULL;}
		    if(GammasKeep!=NULL){Free(GammasKeep); GammasKeep=NULL;}
		    if(BetasKeep!=NULL){Free(BetasKeep); BetasKeep=NULL;}
		    if(JKeep!=NULL){Free(JKeep); JKeep=NULL;}
		    if(gFKeep!=NULL){Free(gFKeep); gFKeep=NULL;}
       if (PrintFlag > 1) {
	       Rprintf((char*) "Freeing wa \n");
	       R_FlushConsole();
       }	   
		        if(wa!=NULL) {Free(wa); wa = NULL;}
		        if(aa!=NULL) {Free(aa); aa = NULL;}
		        if(OnBetas!=NULL) {Free(OnBetas); OnBetas=NULL;}
		        if(OutBetas!=NULL) {Free(OutBetas); OutBetas=NULL;}
		        if(PrevBetas!=NULL) {Free(PrevBetas); PrevBetas=NULL;}
		        if(gammaList1!=NULL) {Free(gammaList1); gammaList1=NULL;}
		        if(gammaList2!=NULL) {Free(gammaList2); gammaList2=NULL;}
		    if(ua!=NULL){ Free(ua); ua=NULL;}
		    if(cc!=NULL){Free(cc); cc=NULL;}  // Current Residual Correlations
		    if(OnMuA!=NULL){Free(OnMuA);  OnMuA=NULL;}// Current MuA
		    if(PrevMuA!=NULL){Free(PrevMuA); PrevMuA = NULL;}  // Prev MuA
	         if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing OnResid \n");
		       R_FlushConsole();
	         }		    
		    if(OnResid!=NULL) {Free(OnResid); OnResid=NULL;}
	         if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing PrevResid \n");
		       R_FlushConsole();
	         }			    
		    if(PrevResid!=NULL) {Free(PrevResid);PrevResid=NULL;}
		     if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing XXA \n");
		       R_FlushConsole();
	         }	
		    if(XXA!=NULL) {Free(XXA);XXA = NULL;}
		     if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing GGAA \n");
		       R_FlushConsole();
	         }	
		    if(GGAA!=NULL) {Free(GGAA); GGAA = NULL;}
	         if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing IGGA \n");
		       R_FlushConsole();
	         }			    
		    if(IGGA!=NULL) {Free(IGGA);   IGGA = NULL;}  
		     if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing ActiveOns \n");
		       R_FlushConsole();
	         }	
		    if(ActiveOns!=NULL) {Free(ActiveOns); ActiveOns = NULL;}
	         if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing AAInactive \n");
		       R_FlushConsole();
	         }			    
		    if(AAInActive!=NULL) {Free(AAInActive); AAInActive = NULL;}
	         if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing CurMaxs \n");
		       R_FlushConsole();
	         }			    
		    if(CurMaxs!=NULL) {Free(CurMaxs); CurMaxs=NULL; }

	       if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing ActiveSigns \n");
		       R_FlushConsole();
	       }		    
		    if (ActiveSigns!= NULL) {Free(ActiveSigns); ActiveSigns = NULL;}
	       if (PrintFlag > 1) {
		       Rprintf((char*) "Freeing OnActives \n");
		       R_FlushConsole();
	       }		    
		    if (OnActives != NULL) {Free(OnActives); OnActives = NULL; }
	       if (PrintFlag > 1) {
		       Rprintf((char*) "All Freed \n");
		       R_FlushConsole();
	       }
        if (WIV != NULL) {Free(WIV); WIV = NULL; }			    
		}

}; 




