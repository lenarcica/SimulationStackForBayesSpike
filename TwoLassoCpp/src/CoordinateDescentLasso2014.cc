///////////////////////////////////////////////////////////////////////////////
//  CoordinateDescentLasso.cc
//
//  TwoLassoPackage Implementation of Coordinate Descent Lasso
//
//   These are C functions that implement the Friedman et al. 2006 paper
//    method of lasso regression known as Coordinate Descent.
//
//  Note that Coordinate Descent Lasso is three things:
//    1. a coordinate descent algorthim for lasso regression
//    2. a dynamic memory management scheme that adds key t(X) %*% X features
//        to perform regression 1
//    3. A weighting management that allows for general lizearized regressions
//        to be performed for Regression.
//
//                        Alan Lenarcic 10/30/2013
//

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

#ifndef RMATH
  #include <Rmath.h>
  #include <R.h>
  #include <stdio.h>
  #define RMATH 0
#endif
#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif

///////////////////////////////////////////////////////////////////////////////
//  CoordinateDescentLasso(int NLen, int kLen,  double *XXX, double *YYY,
//             double *OnBeta, double OnGamma, double *OnRecordBetas,
//             double *FinalBeta, int *pTotalLoops, double MaxEpsilon, double *OnGammas,
//             int InitKKs)
//
//  C function that initiates Coordinate Descent Object and runs the algorithm
//    Inputs:
//       NLen = Sample Size
//       kLen = # of Covariates
//       XXX  = Input Matrix of Covariates
//       YYY  = Input Vector of Response
//       OnBeta = Initial Beta Vector, and Output Beta Vector
//       OnGamma = A single weight for all Lasso Penalties
//       OnRecordBetas = An Optional Matrix for recording Output
//       FinalBeta = OutputBeta vector
//       pTotalLoops = Pointer to total max number of loops to run Descent
//       MaxEpsilon = Convergence Criterion to Stop  Loops
//       OnGammas = A vector of weights for weighted Lasso
int CoordinateDescentLasso(int NLen, int kLen,  double *XXX, double *YYY,
  double *OnBeta, double OnGamma, double *OnRecordBetas,
  double *FinalBeta, int *pTotalLoops, double MaxEpsilon, double *OnGammas,
  int InitKKs, double InverseGammaConstant, double *iiWeights,
  int PrintFlag) {
 int PrintF = PrintFlag - 1;            
 CoordinateDescentObject *CDO = NULL;
 double *RecordDecision;
 int ii;
 
 if (PrintF > 1) {  // Print Routines
		Rprintf((char*) "CoordinateDescentLasso function Starting in\n");
		Rprintf((char*) "NLen = %d, kLen = %d, XXX[0] = %.4f, YYY[0] = %.4f\n",
			                 NLen, kLen, (double) XXX[0], (double) YYY[0]);
	  Rprintf((char*) "OnBeta[0] = %.4f, OnGamma  = %.4f\n", 
			                      (double) OnBeta[0], (double) OnGamma);
		if (OnRecordBetas == NULL) {
				       Rprintf((char*) "OnRecordBetas is NULL\n");
				       R_FlushConsole();
		} else if (OnRecordBetas[0] < 0 ) {
				       Rprintf((char*)"OnRecordBetas Negative, OnRecordBetas[0] = %.4f\n",
				                OnRecordBetas[0] );
		} else {
				       Rprintf((char*)"OnRecordBetas Positive, OnRecordBetas[0] = %.4f\n",
				                OnRecordBetas[0] );
		}
		Rprintf((char*)"FinalBeta[0] = %.4f, pTotalLoops[0] = %d, MaxEpsilon = %.4f\n",
			         FinalBeta[0], pTotalLoops[0], MaxEpsilon);
		R_FlushConsole(); R_ProcessEvents();
		if (OnGammas == NULL) {
			     Rprintf("OnGammas is NULL\n");
			     R_FlushConsole();
		} else {
			     Rprintf("OnGammas[0] = %.4f\n", OnGammas[0]);
			     R_FlushConsole(); R_ProcessEvents();
		}
 }
	 
 // Check that we can record path of coordinate descent
 if (OnRecordBetas != NULL) {
	 if (PrintF > 1) {
		 Rprintf((char*)"CoordinateDescentLasso: OnRecordBetas is not NULL\n");
		 R_FlushConsole(); R_ProcessEvents();
   }
 	 if (OnRecordBetas[0] == -999) {
		 if (PrintF > 1) {
			 Rprintf((char*)"CoordinateDescentLasso: OnRecordBetas[0] = -999, RecordDecision will be NULL\n");
			 R_FlushConsole(); R_ProcessEvents();
	   }	 	 
		 RecordDecision = NULL;
     } else {
		   if (PrintF > 1) {
			   Rprintf((char*)"CoordinateDescentLasso: OnRecordBetas[0] != -999, RecordDecision is NOT NULL\n");
			   R_FlushConsole(); R_ProcessEvents();
	     }	 	 
	     RecordDecision = OnRecordBetas;
     }
 } else {
	 RecordDecision = NULL; 
 }
 if (PrintF > 1) {
   Rprintf("CoordinateDescentLasso: InitKKs = %d, but kLen = %d\n", InitKKs, kLen);
   R_FlushConsole();
 }
 if (InitKKs > kLen) {
  Rprintf("CoordinateDescentObject.cc CDLasso, START ISSUE, InitKKs = %d and kLen = %d\n",
          InitKKs, kLen);
  InitKKs = kLen -1;
  R_FlushConsole();
}
 // Create New CoordinateDescentObject
 if (InitKKs > 0) { 
    CDO = new CoordinateDescentObject(NLen, kLen, XXX, YYY,
      OnBeta, OnGamma, RecordDecision, OnGammas, InitKKs, 
      InverseGammaConstant, PrintFlag);
 } else {
    CDO = new CoordinateDescentObject(NLen, kLen, XXX, YYY,
      OnBeta, OnGamma, RecordDecision, OnGammas, PrintFlag);     
 }  
 if (PrintF > 0) {
   Rprintf("CoordinateDescent: Object Loaded \n"); R_FlushConsole();
 }
 if (PrintF > 0) {CDO->PrintFlag = PrintF;}
 if (NLen > 0 && iiWeights != NULL && iiWeights[0] >= 0) {
      CDO->SetupWeights(iiWeights); CDO->iiWeightsHere = 0;
 }         
                              
 if (CDO == NULL) {
	 if (PrintF > 0) {
		 Rprintf((char*) "CoordinateDescentLasso: No Memory for Object \n");
		 R_FlushConsole(); R_ProcessEvents(); return(-1);
     }
  } else if (CDO->SuccessFlag  < 0) {
	  Rprintf((char*) "CoordinateDescentLasso: Not a successful Loadup \n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1);
  }
  CDO->RunAlgorithm(pTotalLoops[0], MaxEpsilon);
  if (FinalBeta != NULL) {
   	  for (ii = 0; ii < kLen; ii++) { FinalBeta[ii] = CDO->OnBeta[ii]; }
  }
  pTotalLoops[0] = CDO->OnLoop;
  //if (RecordDecision != NULL) {
  //	  MaxSS =  TotalLoops * kLen;
  //	  for (ii = 0; ii < MaxSS; ii++) {OnRecordBetas[ii] = CDO->OnRecordBetas[ii]; }
  //}
  return(1);         
}


///////////////////////////////////////////////////////////////////////////
//  RunAlgorithm
//    Runs a setup Coordinate Descent object
//     until convergence (defined by MaxEpsilon)
//     or maximum number of loops (defined by MaxLoops)
int CoordinateDescentObject::RunAlgorithm(int MaxLoops, double MaxEpsilon) {
     /////////////////////////////////////////////////////////////////////
     //  Make sure that correct numbers for XTResid are located given
     //    beta estimate
     if (PrintFlag > 3) {
       Rprintf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");  R_FlushConsole();
       Rprintf("+++++ CoordinateDescentObject():RunAlgorithm \n");  R_FlushConsole();
       Rprintf("+++++      In CoordinateDescentLasso.cc \n");    R_FlushConsole();
       Rprintf("+++++   About to Run Algorithm, ModFlag = %f \n",
         (double) ModifiedXTResid); R_FlushConsole();  
     }
	   if (ModifiedXTResid > 0) {
	      if (PrintFlag > 3) {
          Rprintf("++++ CDO: RunAlgorithm: Must Modify XTResid\n");
          R_FlushConsole();
        }
		    MakeXTResid();
	      ModifiedXTResid = 0;	   
     }
     
    ///////////////////////////////////////////////////////////////////////
    // OnEp records current move, compares to MaxEpsilon after every loop 
	  OnEp = 0.0;
	  TotMove = 0.0;
	  
	  int ii;
	  int tt = 0;
	  if (PrintFlag > 1) {
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, XTXFlag = %d\n", XTXFlag);
	       R_FlushConsole();
    }  
	  if (PrintFlag > 2) { // Echo
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, Starting with XTResid = \n");
	       PrintVector(this->XTResid, this->kLen);
	       R_FlushConsole();
    }
    if (PrintFlag > 1) {
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, PrintFlag = %d\n", (int) PrintFlag);
	       R_FlushConsole();
    }
    if (PrintFlag > 1) {
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, XTXFlag = %d\n", XTXFlag);
	       R_FlushConsole();
	       if (OnGammas == NULL) {
           Rprintf("++++                     We are doing fixed OnGamma version.\n"); R_FlushConsole();
         } else {
           Rprintf("++++                     We are doing OnGammas vector version.\n"); R_FlushConsole();
         }
    }
      
      //////////////////////////////////////////////////////
      // If kLen = 1 then solution should be immediate
    if (kLen == 1) { 
	      if (PrintFlag > 1) {
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, kLen == 1, insta-solution\n", XTXFlag);
	       R_FlushConsole();
        }
	      double DopeGammas = OnGamma;
	      if (OnGammas != NULL) { DopeGammas = OnGammas[0]; }    
	      if (XTY[0] > DopeGammas) { OnBeta[0] = 
	              ( XTY[0] - DopeGammas ) / ((double) *(pXTX[0]+0));                      
          } else if (XTY[0] < - DopeGammas) { OnBeta[0] = (XTY[0] + DopeGammas) / *(pXTX[0]+0); 
          } else { OnBeta[0] = 0.0;
          } 
	      return(1);
   }
      
      //////////////////////////////////////////////////////////////////////////
      //  OnLoop through iterations for multiple coordinates
      for (OnLoop = 0; OnLoop < MaxLoops; OnLoop++) {
   	    if (PrintFlag > 1) {
	   	      Rprintf((char*) "+++++ CDO: RunAlgorithm, Starting Loop %d/%d\n",
	   	                     OnLoop, MaxLoops); 
            R_FlushConsole(); 
            R_ProcessEvents();
   	    }	
   	    if (PrintFlag > 2 && XTXFlag == 2) {
	        Rprintf((char*) "+++++ CDO:RunAlgorithm, XTXFlag = %d, kLen = %d, ",
            XTXFlag, kLen); R_FlushConsole();
          Rprintf((char*) " OnKappaS = %d, OnKappaMem = %d, PrintRpMatrix\n",
	          OnKappaS, OnKappaMem); R_FlushConsole();
	        PrintRpMatrix(pXTX, kLen, OnKappaMem);
	        R_FlushConsole();
        }      
	      OnEp = 0;
	      
        if (OnLoop >= 1 && OnLoop <= 4 && TDFNoise > 0 && TDFSigma > 0) {
          UpdateTDFNoise();
        }
	      ////////////////////////////////////////////////////////////////////////
	      // Loop through coordinates in order, looping OnCoord
	      for (ii = 0; ii < kLen; ii++) {
		      if (PrintFlag > 3) {
		   	      Rprintf((char*) "+++++ CDO: RunAlgorithm, Loop %d/%d, jj = %d, ",
                   OnLoop, MaxLoops, ii);
              Rprintf((char*) " OnBeta[%d] = %.4f\n",
		   	               OnCoord, (double) OnBeta[OnCoord]); 
              R_FlushConsole(); R_ProcessEvents();
	   	    }
          if (TriggerInvestigateFlag >= 1) {
            DataIntegrity(TriggerInvestigateFlag-1);
          }
	   	    
	   	    RecentMove = 0;
	        UpDateCoord();  // Update coordinate OnCoord
          if (MemoryFailFlag == 1) {
            Rprintf("+++++ CDO: RunAlgorithm, MemoryFailFlag %d/%d, ii=%d, OnC=%d, OnKappaS=%d/%d/%d\n",
              OnLoop, MaxLoops, ii, OnCoord, OnKappaS, OnKappaMem, kLen);
            return(-1);
          }
	        OnEp += fabs(RecentMove); TotMove +=fabs(RecentMove);
	   	    if (PrintFlag > 3) {
		   	      Rprintf((char*) "+++++ CDO: RunAlgorithm, Loop %d/%d, jj = %d,",
                      OnLoop, MaxLoops, OnCoord);
              Rprintf((char*) " Move was %.4f\n",
		   	              RecentMove); 
              R_FlushConsole(); R_ProcessEvents();
	   	    }	  	        	      
	   	    MoveCoord(ii);  //Choose next Coordinate in Order
   	      }
   	if (PrintFlag > 1) {
	   	Rprintf((char*) "+++++ CDO: RunAlgorithm, finished Loop %d/%d, ",
        OnLoop, MaxLoops);
      Rprintf((char*) "OnEp = %.6f, TriggerInvestigateFlag = %d\n",
	   	  (double) OnEp); 
      Rprintf("  Here is XTX Matrix now. \n");
      PrintRpMatrix(pXTX, kLen, OnKappaMem); R_FlushConsole();
      Rprintf("  Here is iiWeights: \n");
      if (iiWeights == NULL) {
        Rprintf("   --  No Weights. \n"); R_FlushConsole();
      } else {
        PrintVector(iiWeights, NLen);
      }
      Rprintf("  Is this good? \n"); 
      R_FlushConsole(); R_ProcessEvents();
   	}   	      
   	if (OnRecordBetas != NULL) {
	   	for (ii = 0; ii < kLen; ii++) {
		    OnRecordBetas[tt] = OnBeta[ii]; tt++;
	   	}
   	}    
   	if (MaxEpsilon >= 0 && OnEp <= MaxEpsilon) {
	   	return(1);
   	}
  }
  if (PrintFlag > 0) {
	  Rprintf((char*) "+++++ CDO:RunAlgorithm, Finished\n", XTXFlag);
	  R_FlushConsole();
  }
  return(1);
}

double CoordinateDescentObject::ShowSFunction(int OnCoord) {
  double Out = SFunction(OnCoord);
  return(Out); 
}
inline double CoordinateDescentObject::SFunction(int OnCoord) {
  if (OnCoord >= kLen || OnCoord < 0) {
    Rprintf("ERROR ERROR ERROR ERROR ERROR ERRROR \n");
	  Rprintf((char*) "ERROR CoordinateDescentObject::SFunction Error OnCoord Is False\n");
	  R_FlushConsole();
    OnCoord = 0;
  }
  double NumberOneMove;
	double DopeGamma;
  if (OnGammas != NULL) {
    DopeGamma = (double) OnGammas[OnCoord];
  } else {
    DopeGamma = (double) OnGamma;
  }
  //double DopeGamma = .5 * OnGamma * InvSXX[OnCoord];
  //double PropMove = OnBeta[OnCoord] + XTResid[OnCoord] *InvSXX[OnCoord];
  //double NotMove = OnBeta[OnCoord] *(pXTX[OnCoord] + OnCoord]) + XTResid[OnCoord];
  if (XTXFlag == 1) {
	  if ((double) OnBeta[OnCoord]  == (double) 0.0) {
		  NumberOneMove = (double) XTResid[OnCoord];
	  } else {
      NumberOneMove = ((double) XTResid[OnCoord]) + ((double) *(pXTX[OnCoord] + OnCoord) * 
                                              (double) OnBeta[OnCoord]);
    }
    if (1 == 1|| DiagShrinkage <= 0) {
	    if ( NumberOneMove > DopeGamma ) {
		    return(  (NumberOneMove - DopeGamma ) * ((double)InvSXX[OnCoord]) );
	    } else if (NumberOneMove < -DopeGamma ) {
		    return((NumberOneMove + DopeGamma ) * (double) InvSXX[OnCoord]);
	    } else {
		    return(0.0);
	    }
	  } else {
	    if ( NumberOneMove > DopeGamma ) {
		    return(  (NumberOneMove - DopeGamma ) / 
          ((double) *(pXTX[OnCoord]+OnCoord) + DiagShrinkage) );
	    } else if (NumberOneMove < -DopeGamma ) {
		    return((NumberOneMove + DopeGamma ) * 
          ((double) *(pXTX[OnCoord]+OnCoord) + DiagShrinkage) );
	    } else {
		    return(0.0);
	    }    
    }
  } else if (XTXFlag == 2) {
    if (XLC[OnCoord] >= kLen) {
      Rf_error("Error: XLC[%d=OnCoord]=%d, but kLen=%d!\n", OnCoord, kLen);
    } else if (XLC[OnCoord] >= OnKappaS) {
      Rf_error("Error: XLC[%d=OnCoord]=%d, but OnKappaS=%d, OnKappaMem=%d, kLen=%d!\n",
        OnCoord, OnKappaS, OnKappaMem, kLen);
    } else if (XLC[OnCoord] < (int) 0 || OnBeta[OnCoord] == 0.0) {
      NumberOneMove = ((double) XTResid[OnCoord]);
    } else {
	    NumberOneMove = ((double) XTResid[OnCoord]) 
	      + ((double) *(pXTX[XLC[OnCoord]] + OnCoord)
	      * (double) OnBeta[OnCoord]);
    }
    if (PrintFlag > 5) {
      Rprintf("SFunction(%d), NumberOneMove=%f, for DopeGamma=%f, OnBeta[%d] = %f, OnGammas[%d] = %f, XTResid[%d] = %f\n",
        OnCoord, NumberOneMove, DopeGamma, OnCoord, OnBeta[OnCoord], OnCoord, OnGammas[OnCoord],
        OnCoord, XTResid[OnCoord]); R_FlushConsole();
      Rprintf(" -- Here's XTResid: "); PrintVector(XTResid, kLen); R_FlushConsole();
    }
    if (1==1 || DiagShrinkage <= 0.0) {
	    if ( NumberOneMove > DopeGamma ) {
		    return(  (NumberOneMove - DopeGamma ) * InvSXX[OnCoord] );
	    } else if (NumberOneMove < -DopeGamma ) {
		    return((NumberOneMove + DopeGamma ) * InvSXX[OnCoord]);
	    } else {
		    return(0.0);
	    }
	  } else if (XLC[OnCoord] >= 0  && XLC[OnCoord] < OnKappaS) {
	    if ( NumberOneMove > DopeGamma ) {
		    return(  (NumberOneMove - DopeGamma ) /
          (*(pXTX[XLC[OnCoord]] + OnCoord) + DiagShrinkage));
	    } else if (NumberOneMove < -DopeGamma ) {
		    return((NumberOneMove + DopeGamma ) /
          (*(pXTX[XLC[OnCoord]] + OnCoord) + DiagShrinkage));
	    } else {
		    return(0.0);
	    }    
    } else {
	    if ( NumberOneMove > DopeGamma ) {
		    return(  (NumberOneMove - DopeGamma ) /
          ( 1.0/InvSXX[OnCoord] + DiagShrinkage ));
	    } else if (NumberOneMove < -DopeGamma ) {
		    return((NumberOneMove + DopeGamma ) /
          (1.0/InvSXX[OnCoord] + DiagShrinkage ));
	    } else {
		    return(0.0);
	    }
    
    }
  } 
  return(0);
}


double CoordinateDescentObject::SFunction3(int OnCoord) {
    if (OnCoord > kLen || OnCoord < 0) {
	    Rprintf((char*) "Error OnCoord Is False\n");
	    R_FlushConsole();
    }
    double NumberOneMove;
	double DopeGamma;
    if (OnGammas != NULL) {
         DopeGamma = (double) OnGammas[OnCoord];
    } else {
         DopeGamma = (double) OnGamma;
    }
    if (OnLambdas3 != NULL) {
       DopeGamma += 1.5 * (double) OnLambdas3[OnCoord] *
              fabs(OnBeta[(int)OnCoord]) * 
              fabs(OnBeta[(int)OnCoord]);    
    } else {
       DopeGamma += 1.5 * (double) OnLambda3 *
              fabs(OnBeta[(int)OnCoord]) * 
              fabs(OnBeta[(int)OnCoord]);
    }
    //double DopeGamma = .5 * OnGamma * InvSXX[OnCoord];
    //double PropMove = OnBeta[OnCoord] + XTResid[OnCoord] *InvSXX[OnCoord];
    //double NotMove = OnBeta[OnCoord] * *(pXTX[OnCoord] + OnCoord) + XTResid[OnCoord]; 
    if (XTXFlag == 1) {
	    if ((double) OnBeta[OnCoord]  == (double) 0.0) {
		   NumberOneMove = (double) XTResid[OnCoord];
	    } else {
          NumberOneMove = ((double) XTResid[OnCoord]) + ((double) *(pXTX[OnCoord] + OnCoord) * 
                                              (double) OnBeta[OnCoord]);
      }
	    if ( NumberOneMove > DopeGamma ) {
		    return(  (NumberOneMove - DopeGamma ) * ((double)InvSXX[OnCoord]) );
	    } else if (NumberOneMove < -DopeGamma ) {
		    return((NumberOneMove + DopeGamma ) * (double) InvSXX[OnCoord]);
	    } else {
		    return(0.0);
	    }
    }  else if (XTXFlag == 2) {
         if (XLC[OnCoord] < 0.0 || OnBeta[OnCoord] == 0.0) {
                         NumberOneMove = ((double) XTResid[OnCoord]);
         } else {
	       NumberOneMove = ((double) XTResid[OnCoord]) 
	                       + ((double) *(pXTX[XLC[OnCoord]] + OnCoord) 
	                                           * (double) OnBeta[OnCoord]);
         }
	    if ( NumberOneMove > DopeGamma ) {
		    return(  (NumberOneMove - DopeGamma ) * InvSXX[OnCoord] );
	    } else if (NumberOneMove < -DopeGamma ) {
		    return((NumberOneMove + DopeGamma ) * InvSXX[OnCoord]);
	    } else {
		    return(0.0);
	    }
    } 
    return(0);
}


int CoordinateDescentObject::DoubleMemory() {

   if (PrintFlag >= 2 ||   TriggerReallocFlag >= 1) {
     Rprintf(" CDO: We've got to double memory sorry, OnKappaMem = %d, kLen=%d, OnKappaS=%d, OnLoop=%d\n", OnKappaMem, kLen, OnKappaS, OnLoop);    
     R_FlushConsole();
     R_ProcessEvents();
   }
   if (OnKappaMem < kLen && pXTX[kLen-1] != NULL) {
     Rf_error("DoubleMemory, pXTX[kLen-1=%d] is not null!\n", kLen-1);
   }
   if (OnKappaMem < kLen && pXTX[OnKappaMem] != NULL) {
     Rf_error("DoubleMemory, Error pXTX[OnKappaMem = %d] is not NULL!\n", kLen);
   }
   if (kXFinder == NULL) {
     Rprintf(" CDO: DoubleMemory, uh, oh, kXFinder is NULL!\n"); 
     R_FlushConsole();
     Rf_error(" --  Double Memory, kXFinder NULL error. \n");
   }
   int NewKappaMem = 0;
   int AddKappaMem = InitKappaMem;
   if (XTXFlag <= 1) {
     Rprintf("ERRRORERRORERRORERROR \n"); R_FlushConsole();
     Rprintf("CDO: DoubleMemory: why are we here?  kLen=%d, OnKappaS =%d, OnKappaMem = %d, XTXFlag = %d.", kLen, OnKappaS, OnKappaMem, XTXFlag);
     R_FlushConsole();
     SuccessFlag = -1;  return(-1);
   }
   if (kLen > 1000) {
     if (InitKappaMem <= 1000) {
       AddKappaMem = 1000;
     }
     if (OnKappaMem > 1000) {
       Rprintf("CDO: OnCoord = %d, DoubleMemory OnKappaS=%d, adding %d to OnKappaMem=%d, MaximumAllocation is %d\n",
         OnCoord, OnKappaS, AddKappaMem, OnKappaMem, MaximumAllocation); R_FlushConsole();
       NewKappaMem = OnKappaMem + AddKappaMem;
     } else {
       NewKappaMem = OnKappaMem * 2;
     }
   } else if (OnKappaMem * 2 > kLen) {
        NewKappaMem = kLen+1;   
   } else {
        NewKappaMem = OnKappaMem * 2;
   }
   if (NewKappaMem >= kLen+1) {
      NewKappaMem = kLen+1;
   }
   if (NewKappaMem > MaximumAllocation) {
     Rprintf("CDO::DoubleMemory: No OnKappaS = %d, OnKappaMem=%d, AddKappaMem=%d, MaximumAllocation = %d!\n", OnKappaS, OnKappaMem, AddKappaMem, MaximumAllocation);
     Rprintf("CDO:: NewKappaMem proposed to be %d \n", NewKappaMem);
     int TotalNonBeta = 0;
     for (int ii = 0; ii < kLen; ii++) {
       if (OnBeta[ii] != 0.0) {
         TotalNonBeta++;
       }
     }
     Rprintf("CDO::DoubleMemory, there are %d TotalNon Zero Beta when OnLambda = %f\n", TotalNonBeta, OnLambda); R_FlushConsole();
     Rprintf("CDO::DoubleMemory: We cite a memory failure and quit. \n"); R_FlushConsole();
     MemoryFailFlag  = 1;
     return(-1);
   }
	 if (PrintFlag > 2 || TriggerReallocFlag >= 1) { 
	   Rprintf((char*)"CDO::DoubleMemory, From OnKappaMem = %d, NewKappaMem = %d\n", 
	     OnKappaMem, NewKappaMem);
	   R_FlushConsole();	   R_ProcessEvents();
   } else if (PrintFlag == 1 || PrintFlag == 2) {
     Rprintf("CDO::DoubleMemory, OnKS=%d/%d, NewKappaMem = %d, kLen=%d -- ", OnKappaS, OnKappaMem, NewKappaMem, kLen); R_FlushConsole();
   }
   //int *NewkXFinder = NULL;
   //double *NewInvSXX = NULL;
    //double *NewXTX = NULL;
	if (PrintFlag > 2) { 
	   Rprintf((char*)"CDO::DoubleMemory, Apportioning NewkXFinder, NewKappaMem = %d\n", 
	     OnKappaMem, NewKappaMem);
     Rprintf("Old kxFinder = \n"); R_FlushConsole();
     PrintVector(kXFinder, OnKappaMem);  Rprintf("\n");
	   R_FlushConsole();	   R_ProcessEvents();
   } 
   if (OnKappaS >= NewKappaMem) {
     Rprintf("ERROR ERROR: CDO: DoubleMemory, error, OnKappaS=%d, but OnKappaMem=%d, NewKappaMem=%d!\n",
       OnKappaS, OnKappaMem, NewKappaMem); R_FlushConsole();
       SuccessFlag = -1; return(-1);
   }  
   if (OnKappaMem == NewKappaMem) {                                                                             
     Rprintf("ERROR, ERROR, CDO: DoubleMemory, error, you didn't go beyond memory requirement. OnKM = %d, NKM = %d\n", OnKappaMem, NewKappaMem);
     R_FlushConsole();
     SuccessFlag = -1;  return(-1);
   }
   if (NewKappaMem > kLen) {
     NewKappaMem = kLen;
   }
   //NewkXFinder = (int *) Calloc( NewKappaMem+1, int);
   if (PrintFlag > 2) {
     Rprintf("CDO:: Double Memory, re allocating kXFinder \n"); R_FlushConsole();
     int CTten = 0;
     for (int ttt = 0; ttt < OnKappaMem; ttt++) {
       if (CTten >= 15) {
         Rprintf("\n"); CTten = 0; R_FlushConsole();
       }
       Rprintf("%d, ", kXFinder[ttt]); R_FlushConsole();
       CTten++;
     }
   }
   kXFinder = Realloc(kXFinder, NewKappaMem, int);
   if (PrintFlag > 2) {
     Rprintf("CDO:: Successfuly re-allocated kXFinder. \n"); R_FlushConsole();
   }
   if (kXFinder == NULL) {
     Rprintf((char*) "CoordinateDescentObject::Double Memory, not enough for NewkXFinder\n");
     SuccessFlag = -1; R_FlushConsole(); R_ProcessEvents(); return(-1);
   }
   //NewInvSCC = (double *) Calloc( NewKappaMem, double);
   //if (NewInvSXX == NULL) {
   //  Rprintf((char*) "CoordinateDescentObject::Double Memory, not enough for NewInvSXX \n");
   //  SuccessFlag = -1; R_FlushConsole(); R_ProcessEvents(); return(-1);
   //}
	if (PrintFlag > 2) { 
	   Rprintf((char*)"CDO::DoubleMemory, Now Apportioning NewXTX = %d*%d = %d\n", 
	                     NewKappaMem, kLen, NewKappaMem * kLen);
	   R_FlushConsole();
     R_ProcessEvents();
   }      
   //NewXTX = (double *) Calloc( (NewKappaMem * kLen) + 1, double);
   int OnRightOldKappaMem = OnKappaMem;
   for (int iti = OnRightOldKappaMem; iti  < NewKappaMem; iti++) {
     if (pXTX[iti] != NULL) {
       if (PrintFlag > 0) {
         Rprintf("CDO::Double Memory, OnK=%d, NK=%d, why does iti=%d have nonnull pXTX?, we're going to free it!\n", OnKappaMem, NewKappaMem, iti);
         R_FlushConsole();
       }
       Rprintf("CDO::Double Memory, pXTX[iti=%d] != NULL, but OnKappaMem = %d, NewKappaMem = %d, kLen=%d, this is a disaster end error. \n",
         iti, OnKappaMem, NewKappaMem, kLen); R_FlushConsole();
       Rf_error("CDO:: DoubleMemory Go Seek Error!\n");
       Free(pXTX[iti]);  pXTX[iti] = NULL;
       if (PrintFlag > 0) {
         Rprintf("pXTX[iti=%d] was non null, it's null now. \n", iti); R_FlushConsole();
       }
     }
     if (PrintFlag >= 3) {
       Rprintf(" -- CDO:: Double Allocate iti = %d/%d", iti, NewKappaMem); R_FlushConsole();
     }
     pXTX[iti] = (double*) Calloc(kLen, double); 
     if (pXTX[iti] == NULL) {
       Rprintf((char*) "CoordinateDescentObject::Double Memory, not enough for NewXTX \n");
       SuccessFlag = -1; R_FlushConsole(); R_ProcessEvents(); return(-1);
     }
   }
   this->OnKappaMem = NewKappaMem;
   
   if (PrintFlag > 2) {
     Rprintf("CDO: DoubleMemory, checking Data Integrity.  -- "); R_FlushConsole();
     *(pXTX[NewKappaMem-1] + kLen - 1) = *(pXTX[NewKappaMem-1] + kLen - 1) +0.0;
     *(pXTX[OnRightOldKappaMem-1] + kLen - 1) = *(pXTX[OnRightOldKappaMem-1] + kLen - 1) +0.0; 
     Rprintf(" -- QUICK PASS. \n"); R_FlushConsole();         
   }
   //int ii, jj, nn1;
   //for (ii = 0; ii < OnKappaMem; ii++) {
   //   NewkXFinder[ii] = kXFinder[ii];
   //}
   int ii;
   if (PrintFlag > 2 || TriggerReallocFlag >= 1) {
     Rprintf("CDO: DoubleMemory, now filling holes in kXFinder.\n"); 
     R_FlushConsole();
   }
   for (ii = OnKappaS; ii < NewKappaMem; ii++) {
      kXFinder[ii] = -10;
   }
   
   if (PrintFlag > 2 || TriggerReallocFlag >= 1) {
     Rprintf("CDO: DoubleMemory, now checking integrity!\n"); R_FlushConsole();
     DataIntegrity(1);
     Rprintf("CDO: Memory has Integrity!\n"); R_FlushConsole();   
   }
   if (PrintFlag == 1 || PrintFlag == 2) {
     Rprintf(" -- PASS1 -- "); R_FlushConsole();
     if (OnKappaMem >= 0 && OnKappaMem < NewKappaMem &&
       pXTX[OnKappaMem] == NULL ) {
       Rf_error("DoubleMemory Error, OnKappaMem = %d, NewKappaMem=%d, kLen=%d, no memory there!\n", 
         OnKappaMem, NewKappaMem, kLen);
     }
     if (pXTX[NewKappaMem-1] == NULL) {
       Rf_error("DoubleMemory Error, NewKappaMem = %d, no memory under it!\n", NewKappaMem);
     }
     *(pXTX[OnRightOldKappaMem] + kLen-1) = *(pXTX[OnRightOldKappaMem] + kLen-1);
     *(pXTX[NewKappaMem-1] + kLen-1) = *(pXTX[NewKappaMem-1] + kLen-1);
     Rprintf(" -- ALL PASS . \n"); R_FlushConsole();     
   }
   //for (ii = 0; ii < OnKappaMem; ii++) {
   //   NewInvSXX [ii] = InvSCC[ii];
   //}
   //for (ii = OnKappaMem; ii < NewKappaMem; ii++) {
   //   NewInvSXX [ii] = -1;
   //}
   //nn1 = 0;
   //for (ii = 0; ii < OnKappaMem; ii++) {
   //     for (jj = 0; jj < kLen; jj++) {
   //        NewXTX[nn1] = XTX[nn1]; 
   //        nn1++;
   //     }
   //}
   //int *GoopNewkXFinder = NewkXFinder;
   //NewkXFinder = kXFinder; kXFinder = GoopNewkXFinder;
   //NewkXFinder = kXFinder + NewkXFinder;
   //kXFinder = NewkXFinder - kXFinder;
   //NewkXFinder = NewkXFinder - kXFinder;
   //FFree(&NewkXFinder);
   //NewInvSXX = InvSXX + NewInvSXX ;
   //InvSXX = NewInvSXX - InvSXX ;
   //NewInvSXX = NewInvSXX - InvSXX ;
   //Free(&NewInvSXX);
   //double *GoopNewXTX = NewXTX;
   //NewXTX = XTX; XTX = GoopNewXTX;
   //NewXTX = (double*) ((double*) XTX + (double*) NewXTX);
   //XTX = (double*) ((double*) NewXTX - (double*) XTX);
   //NewXTX = (double*) ((double*) NewXTX - (double*) XTX);
   //FFree(&NewXTX );
   OnKappaMem = NewKappaMem;
   if (PrintFlag >= 2) {
     Rprintf("DoubleMemory: OnKappaMem = %d, NewKappaMem = %d, OnRightOldKappaMem = %d\n",
       OnKappaMem, NewKappaMem, OnRightOldKappaMem); R_FlushConsole();
   }
   int ACheckDouble = DataIntegrity(0);
   if (ACheckDouble < 0) {
     Rprintf("DoubleMemory: We have DataIntegrity broke when OnKappaMem = %d, kLen=%d, OnLoop = %d\n", OnKappaMem, kLen, OnLoop);
     R_FlushConsole();
   }
   return(1);
}
int CoordinateDescentObject::AddAllNewNonzeroCoords() {
  int MaxInOnKappaS = 0, MinInOnKappaS = kLen;
  int ii;
  if (PrintFlag > 1) {
    Rprintf("CDO: AddAllNewNonZeroCoords, kLen = %, OnKappaS = %d", kLen, OnKappaS); R_FlushConsole();
  }
  /*
  if (OnKappaS > 0) {
    for (ii = 0; ii < OnKappaS; ii++) {
      if (kXFinder[ii] < MinInOnKappaS) {
        MinInOnKappaS = kXFinder[ii];
      }
      if (kXFinder[ii] > MaxInOnKappaS) {
        MaxInOnKappaS = kXFinder[ii];
      }
    }
  } else {
    MaxInOnKappaS = -1;  MinInOnKappaS = kLen + 1;
  }
  */
  if (PrintFlag > 1) {
    Rprintf("CDO:AANNZC:  MinInOnKappaS = %d, MaxInOnKappaS = %d\n",
      MinInOnKappaS, MaxInOnKappaS); R_FlushConsole();
  }
  //int MaxAddedCoord = 0;  int MinAddedCoord = kLen;
  int jj = 0;  int DoAdd = 1;
  int NumAdded = 0;
  int CountAd = 0;
  for (ii = 0; ii < kLen; ii++) {
    if (OnBeta[ii] != 0.0 && (XLC[ii] < 0 || XLC[ii] >= kLen)) {
      CountAd++;
    }
  }
  if (jj == -1) {
  }
  if (CountAd + OnKappaS > MaximumAllocation)  {
    Rprintf("CDO:AANNZC: No, CountAd=%d/%d, and OnKappaS=%d, we will not be able to commit MaximumAllocation=%d\n",
      CountAd, OnKappaS, MaximumAllocation); R_FlushConsole(); \
    MemoryFailFlag = 1;
    return(-1);
  }
  if (CountAd <= 0) {
    return(1);
  }
  for (ii = 0; ii < kLen; ii++) {
    DoAdd = 0;
    if (OnBeta[ii] != 0.0  && (XLC[ii] < 0 || XLC[ii] >= kLen)) {
      AddCoord(ii);  NumAdded++;
      if (MemoryFailFlag == 1) {
        Rprintf("CDO:AANNZC fail due to memory \n"); R_FlushConsole();
        return(-1);
      }
      if (PrintFlag >= 3) {
        Rprintf("CDO:AANZC, just finished AddCoord, %d.  OnKappaS=%d/%d\n", ii, OnKappaS, OnKappaMem); R_FlushConsole();
      }
      if ( NumAdded > 0 && NumAdded %  100 == 0) {
        Rprintf("CDO:AANNZC: We've added a lot of coordinates here NumAdded = %d! though kLen=%d\n", NumAdded, kLen); R_FlushConsole();
      }
      if ( NumAdded > 0 && NumAdded %  500 == 0 ) {
        Rprintf("CDO:AANNZC: We've added a lot of coordinates in this step NumAdded = %d! though kLen=%d\n", NumAdded, kLen); R_FlushConsole();
      }
    }
    /*
      if (ii < MinInOnKappaS) {
        if (PrintFlag > 1) {
          Rprintf("CDO:AANNZC: We found coordinate ii = %d to add, MinInOnKappaS = %d\n",
            ii, MinInOnKappaS); R_FlushConsole();
        }
        AddCoord(ii); NumAdded++;
        MaxAddedCoord = ii;  DoAdd = 1;  MinAddedCoord = ii;     NumAdded++;
      } else if (ii > MaxInOnKappaS) {
        if (PrintFlag > 1) {
          Rprintf("CDO:AANNZC: We found coordinate ii = %d to add, MaxInOnKappaS = %d\n",
            ii, MaxInOnKappaS); R_FlushConsole();
        }
        AddCoord(ii);
        MaxAddedCoord = ii; DoAdd = 1;
      } else {
        DoAdd = 1;
        for (jj = 0; jj < OnKappaS; jj++) {
          if (kXFinder[jj] == ii) {
            DoAdd = 0;
            break;
          }
        }
        if (DoAdd == 1) {
          if (PrintFlag > 1) {
            Rprintf("CDO:AANNZC: We found coordinate ii = %d to add, in middle of OnKappa = %d\n",
              ii); R_FlushConsole();
          }
          AddCoord(ii);  NumAdded++;
          MaxAddedCoord = ii;
        }
      }
     */
    if (PrintFlag > 3) {
      if (DoAdd == 1) {
        Rprintf("CDO:AANNZC:  Added %d, %f \n", ii, OnBeta[ii]); R_FlushConsole();
      } else {
        Rprintf("CDO::AANNZO:   Did Not add %d, %f \n", ii, OnBeta[ii]); R_FlushConsole();
      }
    }
  }
  return(1);
}
int CoordinateDescentObject::AddCoord(int NewCoord) {
   if (MemoryFailFlag == 1) { return(-1); }
   int SFlag = 0;
   if (NewCoord < 0 || NewCoord >= kLen) {
     Rf_error("CoordinateDescent: can't add %d if kLen = %d!\n", NewCoord, kLen);
   }
   if (XLC[NewCoord] >= 0) {
     Rprintf((char*)"CDO: Hey, NewCoord = %d is already in the list!\n"); R_FlushConsole();
     return(0);
   }
   int OnPlacement = OnKappaS; 
   if (PrintFlag > 0) { 
	   Rprintf((char*)"CDO:: Loop %d, AddCoord, adding NewCoord = %d/%d to OnPlacement = %d, max=%d\n", 
	     OnLoop, NewCoord, kLen, OnPlacement, OnKappaMem);
	   R_FlushConsole();	   R_ProcessEvents();
   }
   if (OnKappaS >= OnKappaMem-1 && OnKappaMem < kLen) {
     if (PrintFlag >= 0  || kLen >= LARGEDATA) {
       Rprintf((char*)"CDO:: AddCoord: We're going to have to Double Memory, OnKappaS = %d, OnKappaMem = %d!\n",
         OnKappaS, OnKappaMem); R_FlushConsole();
     }
     SFlag = DoubleMemory();
     if (MemoryFailFlag == 1 && SFlag == -1) {
       Rprintf("CDO:AddCoord[NewCoord=%d, OnKappaS=%d]: We concede a MemoryFailFlag.\n", NewCoord, OnKappaS);
       R_FlushConsole();
       return(-1);
     }
     if (SuccessFlag < 0 || SFlag  < 0) {
       Rprintf("CoordinateDescentObject::AddCoord failed to Double Memory\n");
       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
     }
     if (PrintFlag >=2 || kLen >= LARGEDATA){
       Rprintf((char*)"CDO:: AddCoord: We doubled memory, OnKappaMem = %d!\n",
         OnKappaMem); R_FlushConsole();
     }
   }
   if (OnKappaMem -1 < OnPlacement) {
     SuccessFlag = -1;
     Rf_error((char*) "Hey: We can't do this, we have a memory error OnPlacement = %d, OnMem=%d!\n",
       OnPlacement, OnKappaMem); 
   }
   OnKappaS++;
   if (PrintFlag >= 3) {
     Rprintf("CDO: Inserting NewCoord=%d into XLC position XLC: %d\n", NewCoord, OnPlacement); 
     R_FlushConsole();
   }
   if (NewCoord < 0 || NewCoord >= kLen) {
     Rprintf("ERROR CDO: What, after this NewCoord is now %d/%d which is bad!\n", NewCoord, kLen);
     SuccessFlag = -1;
     Rf_error("Hey: we can't do this!\n");
   }
   if (OnPlacement < 0 || OnPlacement >= OnKappaMem) {
     Rf_error("Error CDO: What, after this NewCoord, OnPlacement = %d? Mem = %d bad!\n", OnPlacement,
       OnKappaMem);
   }
   XLC[NewCoord] = OnPlacement;
   kXFinder[OnPlacement] = NewCoord;
   if (PrintFlag >= 2) {
     Rprintf("CDO: AddCoord, just set kXFinder[%d/%d] \n", OnPlacement, OnKappaMem);  R_FlushConsole();
   }
   //Rprintf("CDO: Adding Coord %d as %d = OnPlacement\n", NewCoord, OnPlacement);
   int ii,  nnPlace, nnii;
   //int jj, nnOnK;
   nnPlace =  OnPlacement * kLen;
   int One = 1; int Zero = 0;
   if (PrintFlag >= 3) {
     Rprintf("CDO: We inserted NewCoord = %d, now filling pXTX. \n", NewCoord); R_FlushConsole();
   }
   double OneD = 1.0; double ZeroD = 0.0;  One = 1;
   if (OnKappaMem < 0  || OnKappaS < 0) {
     Rprintf("********************************************************************************************\n");
     Rprintf("** CoordinateDescentLasso.cc::AddCoord(NewCoord=%d):  Max Out memory error !\n");
     Rprintf("** CDO: Oh No, we can't fit all this memory in the computer as allocated \n");
     Rprintf("** kLen = %d, OnKappaMem = %d, kLen*OnKappaMem = %d, this is bad juju, \n",
       kLen, OnKappaMem, kLen * OnKappaMem); R_FlushConsole();
     Rprintf("**  OnKappaS=%d/%d \n", OnKappaS, OnKappaMem);
     Rprintf("**  Bad News for the Double Memory! \n");
     int OnAll = 0;
     for (ii = 0; ii < kLen; ii++) {
       if (this->OnBeta[ii] != 0.0) {
         OnAll++;
       }
     }
     Rprintf("**  We have %d coordinates actually on and non-zero\n");  R_FlushConsole(); 
     int MyFail = 0;
     for (MyFail = 1; MyFail < OnKappaMem; MyFail++) {
       if (MyFail * kLen <= 0) {
          Rprintf("** Doh, if kLen = %d, then %d * %d = %d is first below zero!\n",
            kLen, kLen, MyFail, kLen*MyFail); R_FlushConsole();
          MyFail = OnKappaMem+1;
       }
     }
     Rf_error("**  AddCoord, ends because memory blacks out!\n");
   }
   double AApply = 0.0;
   if (iiWeights != NULL) {
       if (PrintFlag >= 3) {
         Rprintf("CDO: AddCoord: With Weights, Adding a new placement to pXTX + %d * %d!\n",
           OnPlacement, kLen);
       }
       nnPlace = XLC[NewCoord];
       if (nnPlace < 0 || nnPlace >= kLen * OnKappaMem) {
         Rprintf("CDO: AddCoord, nnPlace=%d, but kLen = %d, XLC[NewCoord=%d] = %d, ERROR !\n", nnPlace, kLen, NewCoord, XLC[NewCoord]);
         Rf_error("CDO:  Integers OnKS = %d/%d too big on AddCoord. \n", OnKappaS, OnKappaMem);
       }
       if (nnPlace + kLen < 0 || nnPlace + kLen >= kLen * OnKappaMem) {
         Rprintf("CDO: AddCoord, nnPlace = %d, but kLen = %d, XLC[NewCoord=%d] = %d, ERROR, nnPlace+kLen = %d!\n",
           nnPlace, kLen, NewCoord, XLC[NewCoord], nnPlace+kLen); R_FlushConsole();
         Rf_error("CDO:  Integers OnKS = %d/%d too big on AddCoord. \n", OnKappaS, OnKappaMem);
       }
       if (PrintFlag >= 3) {
         Rprintf("CDO: AddCoord[OnCoord=%d] to place %d/%d, we are Calculating Weighted size now. \n",
           OnCoord, OnPlacement, OnKappaMem); R_FlushConsole();
         Rprintf("  -- Test Weights -- "); R_FlushConsole();
         for (ii = 0; ii < NLen; ii++) {
           iiWeights[ii] = iiWeights[ii] + 0.0;
         }
         Rprintf("  -- PASS \n"); R_FlushConsole();
         Rprintf(" -- Test XX[ii+NLen*NewCoord] -- "); R_FlushConsole();
         for (ii = 0; ii < NLen; ii++) {
           AApply = XX[NLen*NewCoord+ii];
         }
         Rprintf(" -- PASS. \n"); R_FlushConsole();
         Rprintf(" -- Test pXTX[OnPlacement=%d] -- ", OnPlacement); R_FlushConsole();
         
       }
        for (ii = 0; ii < kLen; ii++) { *(pXTX[OnPlacement] + ii) = 0.0; }
       if (PrintFlag >= 3) {
         Rprintf("CDO: AddCoord[OnCoord=%d] to place %d/%d, we have blanked successfully. \n",
           OnCoord, OnPlacement, OnKappaMem); R_FlushConsole();
       }
        One = 1;
        for (ii = 0; ii < NLen; ii++) {
          AApply = iiWeights[ii] * XX[NLen*NewCoord+ii];
          F77_CALL(daxpy)(&kLen, &AApply, XX + ii, &NLen, pXTX[OnPlacement], &One);
        } 
        if (iWeightedXtX != NULL) { iWeightedXtX[OnPlacement] = 1; }
       if (PrintFlag >= 2) {
         Rprintf("CDO: AddCoord[OnCoord=%d] to place %d/%d, All daxpy performed.. \n",
           OnCoord, OnPlacement, OnKappaMem); R_FlushConsole();
       }  
   } else {
     if (PrintFlag >= 3) {
       Rprintf("CDO:AddCoord: we are now adding %d to pXTX no s at %d\n", NewCoord, OnPlacement);
       R_FlushConsole();
     }
     OneD = 1.0; One = 1;
     F77_CALL(dgemv)("T", &NLen, &kLen, &OneD, XX, &NLen, XX + NewCoord * NLen,
            &One, &ZeroD, pXTX[OnPlacement], &One);	
     if (PrintFlag >= 2) {
       Rprintf("CDO:AddCoord: finished adding %d/%d to pXTX at %d/%d\n", NewCoord, kLen, OnPlacement, OnKappaMem);
       R_FlushConsole();
     }
   }
   //for (ii = 0; ii < kLen; ii++) {
   //  if (ii != NewCoord && XLC[ii] >= 0 && XLC[ii] < OnKappaS) {
   //    *(pXTX[OnPlacement] + ii) = *(pXTX[ XLC[ii]] + NewCoord);
   //  } else {
   //    pXTX[nnPlace] = 0;
   //           nnii = ii * NLen; nnOnK = NewCoord * NLen;
   //           for (jj = 0; jj < NLen; jj++) {
   //              XTX[nnPlace] += XX[ nnii] * XX[nnOnK];
   //              nnii++; nnOnK++;
   //           }
   //  } 
   //  nnPlace++;
   //}
   if (PrintFlag >= 3) {
     Rprintf("CDO:AddCoord: Successfully added NewCoord = %d, kLen=%d\n", NewCoord, kLen);
     R_FlushConsole();
   }
   return(1);
}

////////////////////////////////////////////////
//  Moves the Coordinate of Interest in the Loop
inline int CoordinateDescentObject::MoveCoord(int ii) {
  //if (RunOrder != NULL) {
  //  OnCoord = RunOrder[ii];
  //  return(1);
  //}
  if (OnCoord >= kLen-1) {
	  OnCoord = 0;
  } else if (OnCoord < 0) {
    OnCoord = 0;
  } else {
	  OnCoord++;
  }
	return(1);
}
int CoordinateDescentObject::ReweightCoordinate(int iOnCoord) {
  if (iiWeights == NULL) {
    Rprintf("ReweightCoordinate: Can't do if iiWeights is NULL!\n");
    return(-6);
  }   
  if (iWeightedXtX == NULL) {
    Rprintf("ReweightCoordinate: Why do this if iWeightedXtX is NULL?\n");
    return(-6);
  }            
  double Applier = 0.0;
  if (iOnCoord < 0 || iOnCoord >= kLen) {
    Rf_error("ReweightCoordinate: Not of OnCoord = %d!\n", iOnCoord);
  }
  if (XTXFlag >= 2 && XLC[iOnCoord] < 0) {
    Rf_error("ReweightCoordinate: No, OnCoord=%d is not in model yet, XLC[%d]=%d!\n",
      iOnCoord, iOnCoord, XLC[iOnCoord]);
  }
  if (XTXFlag >= 2 && XLC[iOnCoord] >= OnKappaS) {
    Rf_error("ReweightCoordinate:, No iOnCoord = %d, XLC[%d] = %d!\n",
      iOnCoord, iOnCoord, XLC[iOnCoord]);
  }
  int One = 1;
  if (XTXFlag >= 2) {
  for (int ii = 0; ii < kLen; ii++) {
    *(pXTX[XLC[iOnCoord]] + ii) = 0.0;
  }
  for (int ii = 0; ii < NLen; ii++) {
     Applier = XX[ii + iOnCoord * NLen] * iiWeights[ii];
     F77_CALL(daxpy)(&kLen, &Applier, XX + ii, &NLen, pXTX[XLC[iOnCoord]], &One);
  }
    iWeightedXtX[XLC[iOnCoord]] = 1;     
    return(1);             
  }             
  for (int ii = 0; ii < kLen; ii++) {
    *(pXTX[iOnCoord] + ii) = 0.0;
  }
  for (int ii = 0; ii < NLen; ii++) {
     Applier = XX[ii + iOnCoord * NLen] * iiWeights[ii];
     F77_CALL(daxpy)(&kLen, &Applier, XX + ii, &NLen, pXTX[iOnCoord], &One);
  }
  if (iWeightedXtX != NULL) {
    iWeightedXtX[iOnCoord] = 1;
  }     
    return(1);                          
}
////////////////////////////////////////////////
// UpdateCoord:
//   Based upon the current Coordinate to work with
//   Gets a new fit for the beta from the "SFunction"
//   Then Updates Residuals based upon the degree of 
//   the new movement
int CoordinateDescentObject::UpDateCoord() {
  int SFlag = 0;
  if (OnCoord < 0 || OnCoord >= kLen) {
    Rf_error("UpdateCoord, no OnCoord=%d, kLen=%d, cannot Be!\n", OnCoord, kLen);
  }
  if (iiWeights != NULL && XTXFlag >= 2 && OnBeta[OnCoord] != 0.0 && XLC[OnCoord] >= 0 &&
    iWeightedXtX != NULL && iWeightedXtX[XLC[OnCoord]] <= 0) {
    ReweightCoordinate(OnCoord);
  }  else if (iiWeights != NULL && XTXFlag <= 1 && OnBeta[OnCoord] != 0.0 && iWeightedXtX != NULL && iWeightedXtX[OnCoord] <= 0) {
    ReweightCoordinate(OnCoord);
  }
	double NewBeta = SFunction(OnCoord);
  if (NewBeta != 0.0 && iiWeights != NULL && XTXFlag >= 2 && XLC[OnCoord] >= 0 && iWeightedXtX != NULL &&
    iWeightedXtX[XLC[OnCoord]] <= 0) {
    ReweightCoordinate(OnCoord);
  } else if (NewBeta != 0.0 && iiWeights != NULL && XTXFlag <= 1 && iWeightedXtX != NULL && iWeightedXtX[OnCoord] <= 0) {
    ReweightCoordinate(OnCoord);
  }
  if (PrintFlag > 4) {
		   Rprintf((char*) "CDO: UpdateCoord, OldBeta[%d] = %.4f, NewBeta = %.4f, XTResid[%d]=%f, OnGamma[%d]=%f\n",
		   	 OnCoord, (double) OnBeta[OnCoord], NewBeta, OnCoord,
          XTResid[OnCoord], OnCoord, OnGammas != NULL ? OnGammas[OnCoord] : -999); R_FlushConsole(); R_ProcessEvents();
	}
	int ii, jj;
	if ( (double) NewBeta - (double) OnBeta[OnCoord] == (double) 0.0)  {
		RecentMove = 0;
		return(0);
  }
  double DiffBetas = (double) NewBeta - (double) OnBeta[OnCoord];
     if (PrintFlag > 4 && XTXFlag == 2) {
	     Rprintf((char*)"CDOUpdateCoord, kXifnder updating %d, where XLC[%d] = %d\n", 
	                OnCoord, OnCoord, XLC[OnCoord]);
	     R_FlushConsole(); R_ProcessEvents();
     }
  if (XTXFlag == 2 && XLC == NULL) {
	     Rprintf((char*)"CDO UpdateCoord, Error, XLC is NULL\n");
	     R_FlushConsole(); R_FlushConsole();  SuccessFlag = -1; return(-1);
  } else if (XTXFlag ==2 && (DiffBetas != 0.0 || NewBeta != 0.0) &&
    (XLC[OnCoord] < 0 || XLC[OnCoord] >= kLen)) {
    SFlag = AddCoord(OnCoord);
    if (PrintFlag >= 2) {
      Rprintf("CDO UpdateCoord, we added %d to %d/%d, OnKS=%d, and SFlag = %d\n",
        OnCoord, XLC[OnCoord], OnKappaMem, OnKappaS, SFlag); R_FlushConsole();
    }
    if (SFlag == -1 && MemoryFailFlag == 1) {
      Rprintf("CDO UpdateCoord, we concede a memory Fail Flag. \n");  
      R_FlushConsole(); SuccessFlag = -1;
      return(-1);
    }
    if (MemoryFailFlag == 1)  {
      Rprintf("CDO UpdateCoord, AddCoord(%d) returned Memory Fail. \n", OnCoord); R_FlushConsole();  SuccessFlag = -1;
    }
    if (SFlag < 0 || SuccessFlag < 0) {
       Rprintf((char*) "CDO UpDateCoord, Error Adding Coordinate\n");
       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
    }
    if (XLC[OnCoord] < 0 || XLC[OnCoord] >= OnKappaMem) {
       Rprintf((char*) "CDO UpDateCoord, Something bad when we added Coordinate\n");
       R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
    }
  }

    RecentMove = DiffBetas;
    int tt = OnCoord * kLen; int tt1, tt2;
    if (XTXFlag == 1) {
	    tt = OnCoord * kLen;
	    for (ii = 0; ii < kLen; ii++) {
		    XTResid[ii] = ((double) XTResid[ii]) - 
		             ((double) *(pXTX[OnCoord]+ii) * (double) DiffBetas);
		    tt++;
	    }
    } else if (XTXFlag == 2) {
        tt = XLC[OnCoord];                             
        if (XLC[OnCoord] < 0 || XLC[OnCoord] >= OnKappaS) {
          Rf_error("No, way, supposedly OnCoord=%d, but XLC[%d] = %d is not in Data, OnKappaS=%d, OnKappaMem=%d",
            OnCoord, OnCoord, XLC[OnCoord], OnKappaS, OnKappaMem);
        }
        if (iiWeights != NULL && iWeightedXtX != NULL && XLC[OnCoord] >= 0 && iWeightedXtX[XLC[OnCoord]] <= 0) {
          Rprintf("CoordinateDescentLasso2014.cc::UpdateCoord(OnCoord=%d), No Way, XtXFlag=%d, iWeightedXtX[XLC[%d]=%d]=%d!\n",
            OnCoord, XTXFlag, OnCoord, XLC[OnCoord], iWeightedXtX[XLC[OnCoord]]); R_FlushConsole();
          Rf_error("CoordinateDescentLasso2014.14::UpdatecCoord iWeightedXtX Fail. No Really Bad News!");
        }
        for (ii = 0; ii < kLen; ii++) {
  		     XTResid[ii] = ((double)XTResid[ii]) - ((double) *(pXTX[tt]+ii)) *
		                   ((double) DiffBetas);
  	    }
    } else {
	    for (ii = 0; ii < kLen; ii++) {
		    tt1 = ii * NLen;
		    tt2 = OnCoord * NLen;
		    for (jj = 0; jj < NLen; jj++) {
		      XTResid[ii] = XTResid[ii] - DiffBetas * XX[tt1] * XX[tt2];
		      tt1++; tt2++;
	      }
      }
    }
  OnBeta[OnCoord] = NewBeta;
  return(1);
}


////////////////////////////////////////////////
// UpdateCoord3:
//   Based upon the current Coordinate to work with
//   Gets a new fit for the beta from the "SFunction"
//   Then Updates Residuals based upon the degree of 
//   the new movement
int CoordinateDescentObject::UpDateCoord3() {
      int SFlag = 0;
  if (OnBeta[OnCoord] != 0.0 && iiWeights != NULL && XTXFlag >= 2 && XLC[OnCoord] >= 0 && iWeightedXtX != NULL &&
    iWeightedXtX[XLC[OnCoord]] <= 0) {
    ReweightCoordinate(OnCoord);
  } else if (OnBeta[OnCoord] != 0.0 && iiWeights != NULL && XTXFlag <= 1 && iWeightedXtX != NULL && iWeightedXtX[OnCoord] <= 0) {
    ReweightCoordinate(OnCoord);
  }
	double NewBeta = SFunction3(OnCoord);
  if (NewBeta != 0.0 && iiWeights != NULL && XTXFlag >= 2 && XLC[OnCoord] >= 0 && iWeightedXtX != NULL &&
    iWeightedXtX[XLC[OnCoord]] <= 0) {
    ReweightCoordinate(OnCoord);
  } else if (NewBeta != 0.0 && iiWeights != NULL && XTXFlag <= 1 && iWeightedXtX != NULL && iWeightedXtX[OnCoord] <= 0) {
    ReweightCoordinate(OnCoord);
  }
		if (PrintFlag > 4) {
		   Rprintf((char*) "CDO: UpdateCoord, OldBeta[%d] = %.4f, NewBeta = %.4f\n",
		   	 OnCoord, (double) OnBeta[OnCoord], NewBeta); R_FlushConsole(); R_ProcessEvents();
	   	}
	int ii, jj;
	if ( (double) NewBeta - (double) OnBeta[OnCoord] == (double) 0.0)  {
		RecentMove = 0;
		return(0);
      }
     double DiffBetas = (double) NewBeta - (double) OnBeta[OnCoord];
     if (PrintFlag > 4 && XTXFlag == 2) {
	     Rprintf((char*)"CDOUpdateCoord, kXifnder updating %d, where XLC[%d] = %d\n", 
	                OnCoord, OnCoord, XLC[OnCoord]);
	     R_FlushConsole(); R_ProcessEvents();
     }
     if (XTXFlag == 2 && XLC == NULL) {
	     Rprintf((char*)"CDO UpdateCoord, Error, XLC is NULL\n");
	     R_FlushConsole(); R_FlushConsole();  SuccessFlag = -1; return(-1);
     } else if (XTXFlag ==2 && DiffBetas != 0 && XLC[OnCoord] < 0) {
        SFlag = AddCoord(OnCoord);
        if (MemoryFailFlag == 1) {
          Rprintf("CDO UpdateCoord, Memory FailFlag on coordinate OnCoord=%d, OnKappaS=%d/%d\n",
            OnCoord, OnKappaS, OnKappaMem); R_FlushConsole();
          return(-1);
        }         
        if (PrintFlag >= 2) {
          Rprintf("CDO:UpdateCoord3, well, AddCoord(%d) was performed\n", OnCoord); R_FlushConsole();
        }
        if (SFlag < 0 || SuccessFlag < 0) {
           Rprintf((char*) "CDO UpDateCoord, Error Adding Coordinate\n");
           R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
        }
      }

    RecentMove = DiffBetas;
    int tt = OnCoord * NLen; int tt1, tt2;
    if (XTXFlag == 1) {
	    tt = OnCoord * kLen;
	    for (ii = 0; ii < kLen; ii++) {
		    XTResid[ii] = ((double) XTResid[ii]) - 
		             ((double) *(pXTX[OnCoord]+ii) * (double) DiffBetas);
		    tt++;
	    }
    } else if (XTXFlag == 2) {
        tt = XLC[OnCoord];
        if (tt < 0 || tt >= OnKappaS) {
          Rf_error("CDO: UpdateCoord, no, XLC[%d] = %d, but OnKappaS=%d/%d, not good!\n",
            OnCoord, tt, OnKappaS, OnKappaMem);
        }
        if (iiWeights != NULL && iWeightedXtX != NULL && XLC[OnCoord] >= 0 && iWeightedXtX[XLC[OnCoord]] <= 0) {
          Rprintf("CoordinateDescentLasso2014.cc::UpdateCoord(OnCoord=%d), No Way, XtXFlag=%d, iWeightedXtX[XLC[%d]=%d]=%d!\n",
            OnCoord, XTXFlag, OnCoord, XLC[OnCoord], iWeightedXtX[XLC[OnCoord]]); R_FlushConsole();
          Rf_error("CoordinateDescentLasso2014.14::UpdatecCoord iWeightedXtX Fail. No Really Bad News!");
        }
        for (ii = 0; ii < kLen; ii++) {
		    XTResid[ii] = ((double)XTResid[ii]) - ((double)*(pXTX[tt]+ii)) *
		                  ((double) DiffBetas);
		    tt++;
	    }
    } else {
	    for (ii = 0; ii < kLen; ii++) {
		   tt1 = ii * NLen;
		   tt2 = OnCoord * NLen;
		   for (jj = 0; jj < NLen; jj++) {
		       XTResid[ii] = XTResid[ii] - DiffBetas * XX[tt1] * XX[tt2];
		       tt1++; tt2++;
	       }
        }
   }
   OnBeta[OnCoord] = NewBeta;
   return(1);
}

int CoordinateDescentObject::TestCDO() {
  Rprintf("**** CoordinateDescentObject.cc::TestCDO()  InitKappaMem = %d\n", InitKappaMem);
  if (pXTX == NULL) {
    Rprintf("**** TestCDO: No, pXTX is NULL!\n"); R_FlushConsole();
    Rf_error("**** Eliminate TestCDO!\n");
  }
  int ii, jj;  double CBeta = 0.0;
  if (OnKappaMem * kLen <= 0) {
    Rprintf("**** TestCDO: kLen = %d, OnKappaMem = %d, product is %d, that's an error!\n",
      kLen, OnKappaMem, kLen * OnKappaMem); R_FlushConsole();
    Rf_error("**** TestCDO: kLen Error!\n");
  }
  if (NLen * kLen <= 0) {
    Rprintf("**** TestCDO: NLen = %d, kLen = %d, NLen * kLen = %d, too difficult to store this in memory!\n", NLen, kLen, NLen*kLen);
    Rf_error("**** TestCDO: kLen Error !\n");
  }
  if (XTXFlag <= 1) {
    Rprintf("****  Test XTX because XTXFlag = %d\n", XTXFlag);
    Rprintf("****   Test XTX -- "); R_FlushConsole();
    for (ii = 0; ii < kLen; ii++) {
      if (pXTX[ii] == NULL) {
        Rf_error("XTXFlag = %d, but pXTX[ii=%d] was null, no way!\n", XTXFlag, ii);
      }
      for (int jj = 0; jj < kLen; jj++) {
        CBeta += *(pXTX[ii]+jj);
      }
    }
    Rprintf(" PASS, pXTX Total = %f\n", CBeta); R_FlushConsole();
  } else {
    Rprintf("****  Test pXTX[%d,%d<%d] -- ", kLen, OnKappaS, OnKappaMem);
    R_FlushConsole();
    ii = 0;
    if (OnKappaS >= 1) {
      for (ii = 0; ii < OnKappaS; ii++) {
        if (pXTX[ii] == NULL) {
          Rf_error("XTXFlag = %d, but pXTX[ii=%d/%d/%d] is NULL, no way!\n", XTXFlag, ii, OnKappaS, OnKappaMem);
        }
        for (int jj = 0; jj < kLen; jj++) {
          CBeta += *(pXTX[ii]+jj);
          *(pXTX[ii]+jj) = *(pXTX[ii]+jj);
        }
      }
      Rprintf(" OnTotal = %f -- ", CBeta); R_FlushConsole();
    } else {
      Rprintf("  No Active Coordinates "); R_FlushConsole();
    }
    if (OnKappaMem >= 1 && OnKappaMem > OnKappaS) {
      CBeta = 0.0;
      for (ii = OnKappaS; ii < OnKappaMem; ii++) {
        if (pXTX[ii] == NULL) {
          Rf_error("XTXFlag = %d, but pXTX[ii=%d/%d/%d] is NULL, Over OnKappaS, no way!\n",
            XTXFlag, ii, OnKappaS, OnKappaMem);
        }
        for (int jj = 0; jj < kLen; jj++) {
          CBeta += *(pXTX[ii]+jj);
          *(pXTX[ii]+jj) = *(pXTX[ii]+jj);
        }
      }
      Rprintf(" Off Total = %f, FULL PASS\n", CBeta); R_FlushConsole();
    } else {
      Rprintf(" OnKappaMem = %d, still pass? \n", OnKappaMem); R_FlushConsole();
    }
  }
  Rprintf("****  Test XX \n"); R_FlushConsole();
  if (XX == NULL && XTXFlag == 2) {
    Rprintf("  Test CDO Error, XX is NULL!, XTXFlag = %d\n", XTXFlag); R_FlushConsole();
    Rf_error("  XX Is NULL");
  }
  CBeta = 0.0;
  Rprintf("**** Count XX -- "); R_FlushConsole();
  for (ii = 0; ii < NLen * kLen; ii++) {
    CBeta += XX[ii];
  }
  Rprintf("  PASS, Total XX = %f\n", XX); R_FlushConsole();
  Rprintf("**** DiagShrinkage is %f \n", DiagShrinkage); R_FlushConsole();                      
  Rprintf("**** XTXPassByPointer is %d\n", XTXPassByPointer); R_FlushConsole();

  Rprintf("**** Test Y -- "); R_FlushConsole();
  CBeta = 0.0;
  for (ii = 0; ii < NLen; ii++) {
    CBeta += YY[ii];
    YY[ii] = YY[ii];
  }
  Rprintf(" PASS Total = %f\n", CBeta); R_FlushConsole();  
  
  Rprintf("**** Test XTY \n"); R_FlushConsole();
  if (XTY == NULL) {
    Rf_error("*** Doh!, XTY is NULL\n");
  }
  Rprintf("**** Test XTY -- "); R_FlushConsole();
  CBeta = 0.0;
  for (ii = 0; ii < kLen; ii++) {
    CBeta += XTY[ii];
    XTY[ii] = XTY[ii];
  }
  Rprintf(" PASS Total = %f\n", CBeta); R_FlushConsole();
  
  Rprintf("**** Careful Numerical Check of XTY \n"); R_FlushConsole();
  int NumFails = 0; int OnPrint10 = 0; int AGo = 0;
  double MoreGo;
  if (iiWeights == NULL) {
  for (ii = 0; ii < kLen; ii++) {
    CBeta = 0.0;
    for (jj = 0; jj < NLen; jj++) {
      CBeta += XX[jj + NLen*ii] * YY[jj];
    }
    if (fabs(CBeta - XTY[ii]) >= .0001) {
      NumFails++;
    }
  }} else {
   for (ii = 0; ii < kLen; ii++) {
    CBeta = 0.0;
    for (jj = 0; jj < NLen; jj++) {
      CBeta += XX[jj + NLen*ii] * YY[jj] * iiWeights[jj];
    }
    if (fabs(CBeta - XTY[ii]) >= .0001) {
      NumFails++;
    }
   }
  }
  if (NumFails == 0) {
    Rprintf("**** Numerical Check of XTY Passes!\n"); R_FlushConsole();
  } else {
    Rprintf("**** Oh No, we have number of fails NumFails in XTY = %d\n", NumFails);   R_FlushConsole();
    R_FlushConsole();
    OnPrint10 = 0; AGo = 0;
    for (ii = 0; ii > kLen; ii++) {
      CBeta = 0.0;
      for (jj = 0; jj < NLen; jj++) {
        CBeta += XX[jj + NLen*ii] * YY[jj];
      }
      if (fabs(CBeta-XTY[ii]) >= .0001) {
        AGo++;
        if (OnPrint10 >= 10) { Rprintf("\n"); OnPrint10= 0;}
        Rprintf("%d:%d:%f, ", AGo, ii, CBeta-XTY[ii]); 
        R_FlushConsole();
        OnPrint10++;
      }
    }
    Rf_error("**** Not Good if XTY is Broken! \n"); R_FlushConsole();
  }
  
 Rprintf("**** Test OnBeta \n"); R_FlushConsole();
  if (OnBeta == NULL) {
    Rf_error("*** Doh!, OnBeta is NULL\n");
  }
  Rprintf("**** Test OnBeta -- "); R_FlushConsole();
  CBeta = 0.0;
  for (ii = 0; ii < kLen; ii++) {
    CBeta += OnBeta[ii];
    OnBeta[ii] = OnBeta[ii];
  }
  Rprintf(" PASS Total = %f\n", CBeta); R_FlushConsole(); 
    
  Rprintf("**** Test XTResid \n"); R_FlushConsole();
  if (XTResid == NULL) {
    Rf_error("*** Doh!, XTResid is NULL\n");
  }
  Rprintf("**** Test XTResid -- "); R_FlushConsole();
  CBeta = 0.0;
  for (ii = 0; ii < kLen; ii++) {
    CBeta += XTResid[ii];
    XTResid[ii] = XTResid[ii];
  }
  Rprintf(" PASS Total = %f\n", CBeta); R_FlushConsole(); 
  
  Rprintf("*** Careful Numerical Check of pXtX \n"); R_FlushConsole();
  double MyP = 0.0;
  int *FailPerCoordinate = NULL; 
  
  int TotalXtXFails = 0;
  if (XTXFlag >= 2 && OnKappaS >= 1) {
    FailPerCoordinate = (int*) Calloc(OnKappaS, int);
    for (int iti = 0; iti < OnKappaS; iti++) { FailPerCoordinate[iti] = 0; }
    if (iiWeights != NULL) {
      for (int jj = 0; jj < OnKappaS; jj++) {
        if (iWeightedXtX != NULL && iWeightedXtX[jj] <= 0) {
        } else {
        for (int kk = 0; kk < kLen; kk++) {
           MyP = 0.0;
           for (int iti = 0; iti < NLen; iti++) {
             MyP += XX[kXFinder[jj] * NLen + iti] * XX[kk * NLen + iti] * iiWeights[iti];
           }
           if (fabs(*(pXTX[jj] + kk) - MyP) >= .0001) {
             TotalXtXFails++;
             FailPerCoordinate[jj]++;
           }
        }
        }
      }    
    } else {
      for (int jj = 0; jj < OnKappaS; jj++) {
        for (int kk = 0; kk < kLen; kk++) {
           MyP = 0.0;
           for (int iti = 0; iti < NLen; iti++) {
             MyP += XX[kXFinder[jj] * NLen + iti] * XX[kk * NLen + iti];
           }
           if (fabs(*(pXTX[jj] + kk) - MyP) >= .0001) {
             TotalXtXFails++;
             FailPerCoordinate[jj]++;
           }
        }
      }
    }
  } else if (XX != NULL && iiWeights != NULL) {
    FailPerCoordinate = (int*) Calloc(kLen, int);
    for (int iti = 0; iti < kLen; iti++) { FailPerCoordinate[iti] = 0; }
    for (int jj = 0; jj < kLen; jj++) {
      if (iWeightedXtX != NULL && iWeightedXtX[jj] <= 0) {
      } else {
       for (int kk = 0; kk < kLen; kk++) {
         MyP = 0.0;
         for (int iti = 0; iti < NLen; iti++) {
           MyP += XX[jj *NLen + iti] * XX[kk*NLen + iti] * iiWeights[iti];
         }
         if (fabs(*(pXTX[jj] + kk) - MyP) >= .0001) {
             TotalXtXFails++;
             FailPerCoordinate[jj]++;
         }
       }
      }
    }
  }
  if (TotalXtXFails >= 1 && FailPerCoordinate == NULL) {
    Rf_error("TestCDO: Why are we even Here, TotalXtXFails=%d but you never configured FailsPerCoordinate?\n", TotalXtXFails); R_FlushConsole();
  } else if (TotalXtXFails >= 1) {
    Rprintf("*** TestCDO: Error, we see that XtX is in bad error!\n");
    Rprintf("*** There are a total of %d XtX Errors\n", TotalXtXFails);
    if (iWeightedXtX == NULL) {
      Rprintf("*** There is no influence by iWeightedXtX"); R_FlushConsole();
    } else if (iiWeights == NULL) {
      Rprintf("*** There is no weights whatsoever. \n"); R_FlushConsole();
    } else {
      Rprintf("*** iWeightedXtX is Not NULL we are using it. \n");
    }
    if (XTXFlag >= 2) {
      for (int iti = 0; iti < OnKappaS; iti++) {
        if (FailPerCoordinate[iti] > 0) {
          Rprintf("  -- Fail Coordinate[iti=%d/OnKappaS=%d = %d] with %d total Fails. \n",
            iti, OnKappaS, kXFinder[iti], FailPerCoordinate[iti]); R_FlushConsole(); 
        }
      }
    } else {
      for (int iti = 0; iti < kLen; iti++) {
        if (FailPerCoordinate[iti] > 0) {
          Rprintf("  -- Fail Coordinate[iti=%d/p=%d = %d] with %d total Fails. \n",
            iti, kLen, iti, FailPerCoordinate[iti]); R_FlushConsole(); 
        }
      }    
    }
    Free(FailPerCoordinate);
    Rprintf("*** That's a TestCDO error!\n");
    Rf_error("TestCDO: XTX is not running right. \n");
  }

  Rprintf("**** Careful Numerical Check of XTResid \n"); R_FlushConsole();
  NumFails = 0;
  int NumSuccess = 0;
  double CountDiffDouble = 0.0;
  int *iiFails = Calloc(kLen+1, int);
  double *ddFails = Calloc(kLen+1, double);
  Rprintf("**** iiFails and ddFails Allocated.\n"); R_FlushConsole();
  if (XTXFlag <= 1) {
    Rprintf("**** XTXFlag = %d, we are checking \n", XTXFlag); R_FlushConsole();
    for (ii = 0; ii < kLen; ii++) {
      CBeta = 0.0;
      for (jj = 0; jj < kLen; jj++) {
        if (OnBeta[jj] != 0.0) {
          if (jj + kLen * ii >= kLen * kLen) {
            Rf_error("**** XTXFlag Error, oh no jj = %d, ii=%d, kLen=%d why are we here?", jj, ii, kLen);
          }
          CBeta += *(pXTX[ii]+jj) * OnBeta[jj];
        }
      }
      if (fabs(XTY[ii]- CBeta - XTResid[ii]) >= .0001) {
        iiFails[NumFails] = ii;
        ddFails[NumFails] = XTY[ii] - CBeta;
        NumFails++;
      }
    }
  } else {
    int ATotal = 0;
    Rprintf("**** TestCDO: Check for the Non Zero Lengths of CDO\n"); R_FlushConsole();
    for (ii = 0; ii < kLen; ii++) {
      if (OnBeta[ii] != 0.0 && XLC[ii] < 0) {
        ATotal++;
      }
    }
    Rprintf("**** TestCDO, if OnBeta is not zero but XLC[ii] < 0 for a total of ATotal = %d (0 is perfect, no problem)\n", ATotal);
    R_FlushConsole();
    if (ATotal >= 1) {
      Rprintf("***********************************************************************\n");
      Rprintf("***  Test CDO: Note Too Many OnBetas Issue OnLoop=%d/%d              **\n",
        OnLoop, MaxLoops);
      Rprintf("*** Oh no, ATotal = %d, is total non-zero betas not included yet in XLC\n",
        ATotal); 
      Rprintf("*** This is a suggestion of fail for OnKS = %d/%d                    **\n", 
        OnKappaS, OnKappaMem);
      int ACount = 0;
      for (int iti1 = 0; iti1 < kLen; iti1++) {
        if (OnBeta[iti1] != 0.0) {
          ACount++;
        }
      }
      double AGam = 0.0;  double SqGam = 0.0;
      for (int iti2 = 0; iti2 < kLen; iti2++) {
        AGam += OnGammas[iti2];
        SqGam += OnGammas[iti2] * OnGammas[iti2];
      }
      Rprintf("*** Mean of Gammas is %f and Sd of Gammas is %f , SqGam = %f          **\n",
        AGam/kLen, sqrt((SqGam - AGam*AGam/kLen)/(kLen-1)), SqGam);
      Rprintf("*** Choose how to handle this error.                                  **\n");
      Rprintf("************************************************************************\n");
      R_FlushConsole();
    }
    R_FlushConsole();
    if (ATotal == 0) {
      int COn10 = 0;
      Rprintf("**** ATotal TestCDO now performing. COn10 = %d\n", COn10); R_FlushConsole();
      for (ii = 0; ii < kLen; ii++) {
        CBeta = 0.0;
        if ((OnKappaS >= OnKappaMem && OnKappaMem < kLen) ||
          OnKappaS > kLen) {
          Rf_error("**** TestCDO error, OnKappaS=%d, OnKappaMem=%d!\n", OnKappaS, OnKappaMem);
        }
        //if (COn10 >= 20) {
        //  Rprintf("\n"); R_FlushConsole();   COn10 = 0;
        //}
        //Rprintf("%d, ", ii);  COn10++; R_FlushConsole();
        for (jj = 0; jj < OnKappaS; jj++) {
          if (kXFinder[jj] < 0 || kXFinder[jj] >= kLen) {
            Rf_error("**** TestCDO, no kXFinder[jj=%d] = %d, but p = %d, oh no!\n", jj, kXFinder[jj], kLen); 
          }
          CBeta += *(pXTX[jj] + ii) * OnBeta[kXFinder[jj]];
        }
        CountDiffDouble += fabs(XTY[ii] - CBeta - XTResid[ii]);
        if (fabs(XTY[ii] - CBeta - XTResid[ii]) >= .0001) {
          Rprintf("**** TestCDO, Oh no, we get a fail, ii =%d, NumFails=%d\n", ii, NumFails); R_FlushConsole();
          iiFails[NumFails] = ii;
          ddFails[NumFails] = XTY[ii] - CBeta;
          NumFails++;
        } else {
          NumSuccess++;
        }
      }
    } else {
      int kk;  int jjNLen, iiNLen;
      for (ii = 0; ii < kLen; ii++) {
        CBeta = 0.0;
        for (jj = 0; jj < kLen; jj++) {
          if (OnBeta[jj] != 0.0) {
            if (XLC[jj] >= 0 && XLC[jj] <= OnKappaS && (iWeightedXtX == NULL || (iWeightedXtX[XLC[jj]] >= 1))) {
              CBeta += *(pXTX[XLC[jj]]+ii) * OnBeta[jj];
            } else {
              if (iiWeights == NULL) {
              MoreGo = 0.0;
              jjNLen = jj *NLen;  iiNLen = ii * NLen;
              for (kk = 0; kk < NLen; kk++) {
                MoreGo += XX[kk + jjNLen] * XX[kk + iiNLen];
              }
              } else {
              MoreGo = 0.0;
              jjNLen = jj *NLen;  iiNLen = ii * NLen;
              for (kk = 0; kk < NLen; kk++) {
                MoreGo += XX[kk + jjNLen] * XX[kk + iiNLen] * iiWeights[kk];
              }
              }
              CBeta += MoreGo * OnBeta[jj];
            }
          }
        }
        if (fabs(XTY[ii] - CBeta - XTResid[ii]) >= .0001) {
          iiFails[NumFails] = ii;
          ddFails[NumFails] = XTY[ii] - CBeta;
          NumFails++;
        }
      }    
    }
  }
  if (NumFails == 0) {
    Rprintf("**** Numerical Check of XTResid Passes!  NumFails = %d, NumSuccess = %d, CountDiffDouble = %f.\n", 
      NumFails, NumSuccess, CountDiffDouble); R_FlushConsole();
  } else {
    Rprintf("**** Oh No, we have number of fails NumFails in XTResid.  The NumFails is %d\n", NumFails);   R_FlushConsole();
    R_FlushConsole();
    if (NumFails > kLen) {
      Rprintf("**** Oh No, we have NumFails=%d, kLen=%d, this is not good!\n", NumFails, kLen); R_FlushConsole();
    }
    int OnPrint10 = 0;
    for (ii = 0; ii < NumFails; ii++) {
      if (OnPrint10 >= 10) { Rprintf("\n"); OnPrint10=0; R_FlushConsole();}
      Rprintf("%d:%d:%f, ", ii, iiFails[ii], ddFails[ii]-XTResid[iiFails[ii]]); 
      R_FlushConsole();
      OnPrint10++;
    }
    Rprintf("Well we will free iiFails and ddFails.\n"); R_FlushConsole();
    Free(iiFails); Free(ddFails);   iiFails = NULL;  ddFails = NULL;
    Rf_error("**** Not Good if XTResid is Broken! \n"); R_FlushConsole();
  }
  Rprintf("**** TestCDO freeing iiFails, ddFails\n"); R_FlushConsole();
  Free(iiFails); Free(ddFails);
                       
  Rprintf("**** Successful Test of CDO\n"); R_FlushConsole();
  return(1);
}
//////////////////////////////////////////////////////////////////////////
//  If there are multiple Gammas weights
//   this inserts the Gammas weights into the Coordinate Descent Object
//
int CoordinateDescentObject::SetUpGammas(double *InsertGammas) {
	if (PrintFlag > 1 ) {
		Rprintf((char*)"SetUpGammas, kLen = %d, Setting up InsertGammas[0] = %.4f\n",
		       kLen, InsertGammas[0]);
        R_FlushConsole(); R_ProcessEvents();
        Rprintf((char*)"SetUpGammas, Printing InsertGammas\n");
        R_FlushConsole();
        PrintVector(InsertGammas, kLen);
        R_FlushConsole(); R_ProcessEvents();
   }
   if (InsertGammas == NULL) {
     DOnGammas = 0;
     OnGammas = (double*) Calloc(kLen+3, double);
     if (OnGammas == NULL) {
         SuccessFlag = -1;
         Rprintf((char*)"CoordinateDescentObject: Cannot Apportion for Multiple Gammas\n");
         return(-1);
      }
   }  else {
     DOnGammas = 1;
     OnGammas = InsertGammas;
   }
   return(1);
}

/*
int CoordinateDescentObject::SetUpGammas(double *InsertGammas) {
    if (PrintFlag > 1 ) {
		Rprintf((char*)"SetUpGammas, kLen = %d, Setting up InsertGammas[0] = %.4f\n", 
		    kLen, InsertGammas[0]);
        R_FlushConsole(); R_ProcessEvents();
        Rprintf((char*)"SetUpGammas, Printing InsertGammas\n");
        PrintVector(InsertGammas, kLen);
        R_FlushConsole(); R_ProcessEvents();
   }
   if (OnGammas == NULL) {
     OnGammas = (double*) Calloc(kLen+2, double);
     if (OnGammas == NULL) {
         SuccessFlag = -1;
         Rprintf((char*)"CoordinateDescentObject: Cannot Apportion for Multiple Gammas\n");
         return(-1);
      }
   }
   int ii;
   for (ii = 0; ii < kLen; ii++) {
       OnGammas[ii] = (double) InsertGammas[ii];
   }
   return(1);
}
*/

int CoordinateDescentObject::SetupWeights(double *riiWeights) {
      //Rprintf("CDO: Setting Up Weights, XTXFlag = %d\n", XTXFlag);
      //R_FlushConsole();
      if (PrintFlag > 5) {
        Rprintf("SetupWeights, setting up with riiWeights\n"); R_FlushConsole();
        Rprintf("riiWeights[0] = %.4f\n", riiWeights); R_FlushConsole();
      }
      double ZeroD = 0.0; int One = 1;  //double OneD = 1.0;
      //int Zero = 0;
      int jj = 0, ii;
      if (iWeightedXtX == NULL) {
        iWeightedXtX = (int*) Calloc(kLen,int);
        for (int iti = 0; iti < kLen; iti++) { iWeightedXtX[iti] = 0; }
      }
      //int nntt = 0;
      int TotalXTX = 0;
      if (TotalXTX < 0) {
      }
      if (pXTX == NULL || XTY == NULL) {
         if (PrintFlag>1) {
           Rprintf("SetupWeights:  Weird we must apportion pXTX or XTY by setting up to partial again\n");
           R_FlushConsole();
         }
         SetUpToPartial(5);
      }      
     
      //Rprintf("CDO: SetWeights About to equate iiWeights to riiWeights\n", XTXFlag);
      //R_FlushConsole();      
      iiWeights = riiWeights; iiWeightsHere = 0; 
      //Rprintf("CDO: SetWeights Done equateing iiWeights to riiWeights\n", XTXFlag);
      //R_FlushConsole();         
      if (XX == NULL || YY == NULL) {
         Rprintf("You Kidding, XX or YY are NULL, I ain't setting up weights\n");
              SuccessFlag = -1;
              R_FlushConsole();   return(-1);
      }  
      double Applier = 0;
      //Rprintf("CDO: SetWeights About to do XTXFlag = %d related work\n", XTXFlag);
      R_FlushConsole();         
      if (XTXFlag == 2) {
         //Rprintf("CDO: SetWeights XTXFlag was 2 Moving On\n");
         //R_FlushConsole();
         if (PrintFlag>6) {
           Rprintf("  SetupWeights: XTXFlag = %d,NLen = %d, kLen = %d\n",
                 XTXFlag, NLen, kLen);
         }                 
         ZeroD = 0.0;
         F77_CALL(dscal)(&kLen, &ZeroD, InvSXX, &One);
         for (jj = 0; jj < OnKappaS; jj++) {
           F77_CALL(dscal)(&kLen, &ZeroD, pXTX[jj], &One);
         }
         F77_CALL(dscal)(&kLen, &ZeroD, XTY, &One);
         for (ii = 0; ii < NLen; ii++)  {
           Applier = iiWeights[ii] * YY[ii]; 
           F77_CALL(daxpy)(&kLen, &Applier, XX + ii, &NLen,
             XTY, &One);                   
         }
         TotalXTX=0;
         for (jj = 0; jj < kLen; jj++) {
           for (ii = 0; ii < NLen; ii++) {
             InvSXX[jj] += iiWeights[ii] * XX[TotalXTX] * XX[TotalXTX];
             TotalXTX++;
           }
           InvSXX[jj] += DiagShrinkage;
         }
         for (jj = 0; jj <  OnKappaS; jj++) {
           for (ii = 0; ii < NLen; ii++)  {
             Applier = iiWeights[ii] * XX[ kXFinder[jj] * NLen + ii]; 
             F77_CALL(daxpy)(&kLen, &Applier, XX + ii, &NLen, pXTX[jj], &One);   
           }
           if (iWeightedXtX != NULL) {
             iWeightedXtX[jj] = 1;
           }
         }          
      }  else {
        //Rprintf("CDO:SetupWeights, XTXFlag = %d, What to Do\n", XTXFlag);
        //R_FlushConsole();        
        if (PrintFlag >= 4) {
            Rprintf("CDO:SetupWeights, XTXFlag = %d, and NLen < kLen, Starting\n", XTXFlag);
             R_FlushConsole();   R_ProcessEvents();
        }
        F77_CALL(dscal)(&kLen,&ZeroD, InvSXX, &One);
        for (jj = 0; jj < kLen; jj++) {
          F77_CALL(dscal)(&kLen, &ZeroD, pXTX[jj], &One);
        }
        F77_CALL(dscal)(&kLen, &ZeroD, XTY, &One);
        for (ii = 0; ii < NLen; ii++)  {
          F77_NAME(daxpy)(&kLen, &Applier, XX + ii, &NLen,
            XTY, &One);                       
        }
        TotalXTX = 0;
        for (jj = 0; jj < kLen; jj++) {
          for (ii = 0; ii < NLen; ii++) {
            InvSXX[jj] += XX[TotalXTX] * XX[TotalXTX] * iiWeights[ii];
            TotalXTX++;   
          }
          InvSXX[jj] += DiagShrinkage;
        }
        for (jj = 0; jj < kLen-1; jj++) {
          TotalXTX = kLen - jj;
           for (ii = 0; ii < NLen; ii++) {
              Applier = iiWeights[ii] * XX[ii + jj * NLen];
              F77_CALL(daxpy)(&TotalXTX, &Applier, XX+ii, &NLen, pXTX[jj] + jj, &One);
           }
           for (ii = 0; ii < jj; ii++) {        
              *(pXTX[jj] + ii) = *(pXTX[ii] + jj);
           }
           if (iWeightedXtX != NULL) { iWeightedXtX[jj] = 1; }
        }
    }
    if (PrintFlag > 6) {
      Rprintf("  SetupWeights:   About to Setup InvSXX\n");
      R_FlushConsole();
    }
    if (InvSXX == NULL)  {
      Rprintf("InvSXX is NULL, uh oh\n"); R_FlushConsole(); SuccessFlag = -1;
      return(-1);
    }
    for (int iti = 0; iti < kLen; iti++) {
      InvSXX[iti] =  1.0 / InvSXX[iti];
    }
    //Rprintf("CDO: Finished, Setting Up Weights, remaking XTResid\n");
    //  R_FlushConsole();
    if (PrintFlag > 6) {
      Rprintf("  SetupWeights:   Reached End, about to MakeXTResin\n");
      R_FlushConsole();
    }
    MakeXTResid();
    return(1);
}

////////////////////////////////////////////////////////////////////
// Refreshing the Weights
//
//  Given new weights riiWeights for eacy Yi,Xi, we need to reweight entire
//  XTX matrix, XTY vector.  
// 
//  Given that XTXFlag = 2, then XTX matrix is only partial 
//   (contains only previously non-zero ceofficients)
//
int CoordinateDescentObject::RefreshWeights(double *riiWeights) {
  double ZeroD = 0.0; int One = 1;  double OneD = 1.0;
  int Zero = 0;
  int jj=0, ii=0;
  //int nntt = 0;
  int TotalXTX = 0;
  if (TotalXTX < 0) {
  }
      
  iiWeights = riiWeights;
  if (XX == NULL || YY == NULL) {
    Rprintf("You Kidding, XX or YY are NULL, I ain't setting up weights\n");
    SuccessFlag = -1;
    R_FlushConsole();   return(-1);
  }  
  double Applier = 0;
  if (XTXFlag == 2) {
          //Rprintf("CDO: SetWeights XTXFlag was 2 Moving On\n");
         //R_FlushConsole();
         if (PrintFlag>6) {
           Rprintf("  SetupWeights: XTXFlag = %d,NLen = %d, kLen = %d\n",
                 XTXFlag, NLen, kLen);
         }                 
         ZeroD = 0.0;
         F77_CALL(dscal)(&kLen, &ZeroD, InvSXX, &One);
         if (OnKappaS > 0) {
         for (jj = 0; jj < OnKappaS; jj++) {
           F77_CALL(dscal)(&kLen, &ZeroD, pXTX[jj], &One);
         }}
         F77_CALL(dscal)(&kLen, &ZeroD, XTY, &One);
         for (ii = 0; ii < NLen; ii++)  {
           Applier = iiWeights[ii] * YY[ii]; 
           F77_CALL(daxpy)(&kLen, &Applier, XX + ii, &NLen,
             XTY, &One);                   
         }
         TotalXTX=0;
         for (jj = 0; jj < kLen; jj++) {
           for (ii = 0; ii < NLen; ii++) {
             InvSXX[jj] += iiWeights[ii] * XX[TotalXTX] * XX[TotalXTX];
             TotalXTX++;
           }
           InvSXX[jj] += DiagShrinkage;
         }
         if (OnKappaS > 0) {
         for (jj = 0; jj <  OnKappaS; jj++) {
           if (iWeightedXtX != NULL && OnBeta[kXFinder[jj]] == 0.0) {
             iWeightedXtX[jj] = 0;
           } else {
           for (ii = 0; ii < NLen; ii++)  {
             Applier = iiWeights[ii] * XX[ kXFinder[jj] * NLen + ii]; 
             F77_CALL(daxpy)(&kLen, &Applier, XX + ii, &NLen, pXTX[jj], &One);   
           }
             if (iWeightedXtX != NULL) {
               iWeightedXtX[jj] = 1;
             }
           }
         }
         } 
  } else {
        //Rprintf("CDO:SetupWeights, XTXFlag = %d, What to Do\n", XTXFlag);
        //R_FlushConsole();        
        if (PrintFlag >= 4) {
            Rprintf("CDO:SetupWeights, XTXFlag = %d, and NLen < kLen, Starting\n", XTXFlag);
             R_FlushConsole();   R_ProcessEvents();
        }
        F77_CALL(dscal)(&kLen,&ZeroD, InvSXX, &One);
        for (jj = 0; jj < kLen; jj++) {
          F77_CALL(dscal)(&kLen, &ZeroD, pXTX[jj], &One);
        }
        F77_CALL(dscal)(&kLen, &ZeroD, XTY, &One);
        for (ii = 0; ii < NLen; ii++)  {
          F77_NAME(daxpy)(&kLen, &Applier, XX + ii, &NLen,
            XTY, &One);                       
        }
        TotalXTX = 0;
        for (jj = 0; jj < kLen; jj++) {
          for (ii = 0; ii < NLen; ii++) {
            InvSXX[jj] += XX[TotalXTX] * XX[TotalXTX] * iiWeights[ii];
            TotalXTX++;   
          }
          InvSXX[jj] += DiagShrinkage;
        }
        for (jj = 0; jj < kLen-1; jj++) {
          TotalXTX = kLen - jj;
          if (iWeightedXtX != NULL) { iWeightedXtX[jj] = 1; }
           for (ii = 0; ii < NLen; ii++) {
              Applier = iiWeights[ii] * XX[ii + jj * NLen];
              F77_CALL(daxpy)(&TotalXTX, &Applier, XX+ii, &NLen, pXTX[jj] + jj, &One);
           }
           if (jj < kLen-1) {
           for (ii = jj+1; ii < kLen; ii++) {        
              *(pXTX[ii] + jj) = *(pXTX[jj] + ii);
           }
           }
        }
        iWeightedXtX[kLen-1] = 1;
    }
    // Rprintf("Post Refresh Weights: XTX, kLen = %d \n", kLen);
    // PrintRMatrix(XTX, kLen, kLen); R_FlushConsole();
    // Rprintf("And XTY for Good measure\n");
    // PrintVector(XTY, kLen); R_FlushConsole();
    // return(-1); 
         
    for (jj = 0; jj < kLen; jj++) { InvSXX[jj] = 1.0 / InvSXX[jj]; }
    MakeXTResid();
    return(1);  
}

///////////////////////////////////////////////////////////////////////////////
//  If There are memory constraints for Coordinate Descent
//   Then this procedure sets it up such that only partial
//    results of XTX, XTY, XTResid are maintained.
int CoordinateDescentObject::SetUpToPartial(int InitKKs) {
  int One = 1; double OneD = 1.0; double ZeroD = 0.0;
  // First Step is to Clear memory if relevant
  if (XTXPassByPointer > 0) {
    Rprintf("CDO: SetupToPartial: No Way are we doing this XTXPassBy Pointer > 0\n");
    R_FlushConsole();  return(-1);
  }
  if (XTXPassByPointer > 0) {
    Rprintf("CDO: SetupToPartial: No Way are we doing this XTYPassByPointer > 0\n");
    R_FlushConsole();  return(-1);
  }
  if (PrintFlag > 0) {
    Rprintf("CDO: SetupToPartial: InitKKs = %d, kLen = %d\n"); R_FlushConsole();
  }
  if (InitKKs > kLen) {
    Rprintf("CDO: SetupToPartial, ISSUE, InitKKs = %d and kLen = %d\n",
            InitKKs, kLen);
    InitKKs = kLen -1;
    R_FlushConsole();
  }
  if (PrintFlag > 0) {
		   Rprintf((char*) "CDO:  SetUpToPartial  Starting, Clearing Mem\n");
		   R_FlushConsole(); R_ProcessEvents();
  }
  if (pXTX != NULL) {
    if (PrintFlag >= -2) {
      Rprintf("CDO: Hey, pXTX Is Not NULL at SetUpToPartial!\n"); R_FlushConsole();
    }
    for (int ii2 = 0; ii2 < kLen; ii2++) {
      if (pXTX[ii2] != NULL) {
        Free(pXTX[ii2]);  pXTX[ii2] = NULL;
      }
    }
    Free(pXTX); pXTX = NULL;
  }


  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing XTY \n"); R_FlushConsole(); R_ProcessEvents();
  }	   
	FFree(&XTY);  
  
  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing XTResid \n"); R_FlushConsole(); R_ProcessEvents();
  }
  FFree(&XTResid);
     
  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing XLC \n"); R_FlushConsole(); R_ProcessEvents();
  }
  FFree(&XLC);
  
  if (PrintFlag > 1) {
	  Rprintf((char*) "CDO: SetUpToPartial Freeing kXFinder \n"); R_FlushConsole(); R_ProcessEvents();
  }
  FFree(&kXFinder);
              
  if  (PrintFlag > 0) {
		   Rprintf((char*) "CDO:  SetUpToPartial, InitKKs = %d, but p=%d, n=%d\n", InitKKs, kLen, NLen);
		   R_FlushConsole(); R_ProcessEvents();
  }
     
     // InitKKs is a number that is smaller than kLen, the initial size of XTX
  if (InitKKs <= 0 ) {
    Rprintf((char*) "CoordinateDescentObject2014:SetUpToPartial() InitKKs set to less than Zero, false hypothesis\n");
    R_FlushConsole(); R_ProcessEvents();  SuccessFlag = -1; return(-1);
  } else if (InitKKs > kLen) {
    Rprintf((char*)"CoordinateDescentObject2014:SetUpToPartial(): InitKKs set to more than kLen, no use\n");
    Rprintf("   InitKKs = %d, and kLen = %d. \n", InitKKs, kLen);
      // R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1); 
    InitKKs = kLen;
  }
  OnKappaMem = InitKKs;  OnKappaS = 0;
	int Maxss = 0, ii, nnOn;
	   
  if (Maxss < 0) {
  
  }
	int ActiveBetas = 0, jj;
  if (OnBeta == NULL) {
    Rf_error("CoordinateDescentObject2014:SetupToPartial(): can't because OnBeta is invalid. \n");
  }
	for (ii = 0; ii < kLen; ii++) {
		if (OnBeta[ii] != 0) { ActiveBetas++; }
  }
  if (ActiveBetas >= 1 && PrintFlag >= 1) {
    Rprintf("CoordinateDescentObject2014:SetupToPartial(): You know that there were %d active Betas at the start non-zero\n", ActiveBetas);
    R_FlushConsole();
  }
  OnKappaMem+= ActiveBetas;
  if (OnKappaMem > kLen) { OnKappaMem = kLen; }
	     //int nnii, tt;
	SuccessFlag = 1;
	XTXFlag  = 2;
  RecentMove = 0;

	     kLen = kLen;
	     Maxss = (kLen * NLen); 

       XLC = (int *) Calloc( kLen+2, int);
	     if (XLC == NULL) {
		     Rprintf("CoordinateDescentObject:: Calloc Fails on XLC \n");
		     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
	     }	     
       for (ii = 0; ii < kLen; ii++) { XLC[ii] = -10; }

       kXFinder  = (int *) Calloc( OnKappaMem + 2, int);
       if (kXFinder  == NULL) {
		     Rprintf("CoordinateDescentObject:: Calloc Fails on kXFinder  \n");
		     R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
	     }	     
           for (ii = 0; ii < OnKappaMem; ii++) { kXFinder[ii] = -1; } 

        OnCoord = 0; OnLoop = 0;
     if (PrintFlag >= 3) {
       Rprintf((char*)"CoordinateDescent, Apportioning XTX\n");
       R_FlushConsole(); R_ProcessEvents();
     }                 
     pXTX = (double **) Calloc(kLen, double*);
     for (int jj = 0; jj < kLen; jj++) {
       pXTX[ii] = NULL;
     }
		 if (pXTX == NULL) {
			 Rprintf("CoordinateDescentObject:: Calloc Fails on XTX \n");
			 R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
		 }
     for (int jj = 0; jj < OnKappaMem; jj++) {
		   pXTX[jj] = (double *) Calloc( kLen, double);
		   if (pXTX[jj] == NULL) {
		  	 Rprintf("CoordinateDescentObject:: Calloc Fails on pXTX[jj=%d] \n", jj);
		  	 R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
	  	 }
     }

		 if (PrintFlag >= 3) {
       Rprintf((char*)"CoordinateDescent, Apportioning XTY\n");
       R_FlushConsole(); R_ProcessEvents();
     }		     
	  	XTY = (double *) Calloc(kLen+4, double);
	  	if (XTY == NULL) {
   			Rprintf("CoordinateDescentObject:: Calloc Fails on XTY\n");
		  	R_FlushConsole(); SuccessFlag = -1; return(-1);
      }
     if (PrintFlag > 1) {
       Rprintf("CDO: Setup To Partial: XTY Exists, filling with YY[%d]\n",
         NLen); R_FlushConsole();
     }
      F77_CALL(dscal)(&kLen, &ZeroD, XTY, &One);
      F77_CALL(dgemv)("T", &NLen, &kLen,
		     &OneD, XX, &NLen,
		       YY, &One, &ZeroD,
		       XTY, &One);    
       //Rprintf("XTY After Fortran is \n");
     if (PrintFlag > 1) {
       Rprintf("XTY: After Fortran is : "); PrintVector(XTY, kLen);
       Rprintf("\n"); R_FlushConsole();
     }  
     if (PrintFlag > 1) {
       Rprintf("CDO: Setup To Partial: XTYResid Create Length %d\n",
         kLen); R_FlushConsole();
     }
     XTResid = (double *) Calloc( kLen+4, double);
     if (XTResid == NULL) {
		    Rprintf("CoordinateDescentObject:: Calloc Fails on XTResid \n");
		    R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
	   }
	   if (PrintFlag >= 1) {
       Rprintf("CDO: XTYResid was Created Now to fill it!\n");
       R_FlushConsole();
     }
     // Rprintf("XTY After Manual is \n");
     // PrintVector(XTY, kLen);          
      F77_CALL(dcopy)( &kLen, XTY, &One, XTResid, &One );
      //for (ii = 0; ii < kLen; ii++) { XTResid[ii] = (double ) XTY[ii]; }	    
    	jj = 0, nnOn = 0;
		if (ActiveBetas > 0) { 
		  if (PrintFlag >= 1) {
        Rprintf("CDO: SetupToPartial ActiveBetas > 0 \n;"); R_FlushConsole();
      }
      AddAllNewNonzeroCoords();
       if (OnLambda3 > 0) {
        if (PrintFlag >= 2) {
          Rprintf("CDO: SetupToPartial: OnLambda3 supposedly greater than zero, no way\n");
          R_FlushConsole();
        }
		    for (ii = 0; ii < kLen; ii++) {
			    if (OnBeta[ii] != 0) {   
				    OnCoord = ii;
				    UpDateCoord3();
			    }
	      }	    
       }
    }
	   if (ActiveBetas > 0) {
		   if (PrintFlag > 0) {
			   Rprintf((char*) "CDO:: SetupToPartial: ActiveBetas > 0, going to Make XTResid\n");
			   R_FlushConsole();
	     }
		   MakeXTResid();
       } else {		
		     if (PrintFlag > 0) {
			     Rprintf((char*) "CDO:: ActiveBetas == 0, going to Make XTResid = XTY: \n");
			     PrintVector( XTY, kLen); Rprintf("\n");
			     R_FlushConsole();
			     Rprintf((char*) "\nCDO: XX is \n");
			     PrintRMatrix( XX, NLen, kLen);
			     Rprintf((char*)"\nCDO: YY is \n");
			     PrintVector(YY, NLen);Rprintf((char*)"\n"); R_FlushConsole();
	       }	          
		    //for (ii = 0; ii < kLen; ii++) { XTResid[ii] = XTY[ii]; }
       }        

 
  OnCoord = 0; OnLoop = 0;
  if (PrintFlag > 1) {
    Rprintf((char*)"CoordinateDescent, InvSXX Create?\n");
    R_FlushConsole(); R_ProcessEvents();
  }        
  InvSXX = (double*) Calloc(kLen +2, double);
  if (InvSXX == NULL) {
		Rprintf("CoordinateDescentObject:: Calloc Fails on InvSXX \n");
		R_FlushConsole(); R_ProcessEvents(); SuccessFlag = -1; return(-1);
	}
  if (iiWeights == NULL) {
	for (ii = 0; ii < kLen; ii++) { 
    InvSXX[ii] = 0.0;
    nnOn = ii * NLen;  
    for (jj = 0; jj < NLen; jj++) {
      InvSXX[ii] += ((double) XX[nnOn] * XX[nnOn]); nnOn++;
    }
    InvSXX[ii] += DiagShrinkage;
    InvSXX[ii] = ((double) 1.0) / ((double) InvSXX[ii]);
  }} else {
  for (ii = 0; ii < kLen; ii++) { 
    InvSXX[ii] = 0.0;
    nnOn = ii * NLen;  
    for (jj = 0; jj < NLen; jj++) {
      InvSXX[ii] += ((double) XX[nnOn] * XX[nnOn] * iiWeights[jj]); nnOn++;
    }
    InvSXX[ii] += DiagShrinkage;
    InvSXX[ii] = ((double) 1.0) / ((double) InvSXX[ii]);
  }}
  if (PrintFlag > 1) {
    Rprintf((char*)"CoordinateDescent, Ready to go with Gamma = %.4f,", OnGamma);
    Rprintf(" NLen = %d, kLen = %d\n", NLen, kLen);
    R_FlushConsole();
  }
  return(1);	
}

int CoordinateDescentObject::MakeXTResid() {
  int jj;
  int One = 1; double OneD = 1.0; 
  if (OneD == -1.0) {
  
  }
  //int ii;  
  //double ZeroD = 0;
  //double NegOneD = -1.0; 
  //int Zero = 0;
    //Rprintf((char*) "CDO->MakeXTResid() Starting \n");
    //R_FlushConsole(); R_ProcessEvents();
  if (XTResid == NULL)    {
    Rprintf("What the heck, XTResid is null?");
     R_FlushConsole(); return(-1);
  }
  int nntt = 0;
  if (nntt < 0) {
    Rprintf("MakeXTResid \n");  R_FlushConsole();
  }
  if (XTXFlag == 1) {
    //Rprintf((char*) "CDO->MakeXTResid() Starting  XTXFlag = 1\n");
    //R_FlushConsole(); R_ProcessEvents();	  
	    nntt = 0;
	    
//	    BLAS_extern void   /* DCOPY - copy x to y */
//      F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
//     		double *dy, const int *incy);
      F77_CALL(dcopy)( &kLen, XTY, &One, XTResid, &One );
//	   	for (ii = 0; ii < kLen; ii++) {
//		   	XTResid[ii] = ((double) XTY[ii]);
//	    }
    //Rprintf((char*) "CDO->MakeXTResid() Starting nextstep XTX and OnBeta\n");
    //R_FlushConsole(); R_ProcessEvents();		
/* DGEMV - perform one of the matrix-vector operations */
/* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,  */
//BLAS_extern void
//F77_NAME(dgemv)(const char *trans, const int *m, const int *n,
//		const double *alpha, const double *a, const int *lda,
//		const double *x, const int *incx, const double *beta,
//		double *y, const int *incy); 
     double alpha = -1.0;
     for (int jj = 0; jj < kLen; jj++) {
       if (OnBeta[jj] != 0.0) {
         if (iiWeights != NULL && iWeightedXtX != NULL && iWeightedXtX[jj] == 0) {
         
         }
         if (iiWeights == NULL || iWeightedXtX == NULL || (iWeightedXtX != NULL && iWeightedXtX[jj] == 1)) {
           alpha = -1.0 * OnBeta[jj];
           F77_CALL(daxpy)(&kLen, &alpha, pXTX[jj], &One, XTResid, &One);
         } else if (iiWeights != NULL && iWeightedXtX != NULL && (iWeightedXtX != NULL && iWeightedXtX[jj] <= 0)) {
           ReweightCoordinate(jj);
           alpha = -1.0 * OnBeta[jj];
           F77_CALL(daxpy)(&kLen, &alpha, pXTX[jj], &One, XTResid, &One);
         }
       }
       //F77_CALL(dgemv)("T", &kLen, &kLen, &alpha, XTX, 
       //  &kLen, OnBeta, &One, &OneD, XTResid, &One); 
     } 
    //Rprintf("MakeXTResid: Just Did Residual Make vector, noting wrong\n");
    //R_FlushConsole(); R_ProcessEvents();
	   // for (jj = 0; jj < kLen; jj++) {
		 //   if (OnBeta[jj] == 0.0) {
		 //   } else {
		 //	    nntt = kLen * jj;
		//		for (ii = 0; ii < kLen; ii++) {
		//			XTResid[ii] -= ((double) XTX[nntt]) * ((double) OnBeta[jj]); 
		//			nntt++;
		//        }
    //        }
    //    }
  } else if (XTXFlag == 2) {
    //Rprintf((char*) "CDO->MakeXTResid() XTXFlag == 2 \n");
    //R_FlushConsole(); R_ProcessEvents();
    if (XTResid == NULL) {
	    Rprintf((char*) "CDO-> What the heck, XTResid is NULL\n");
	    R_FlushConsole();
    }
    if (XTY== NULL) {
	    Rprintf((char*) "CDO-> What the heck, XTY is NULL\n");
	    R_FlushConsole();
    }
    One = 1;
    F77_CALL(dcopy)( &kLen, XTY, &One, XTResid, &One );
      // for (ii = 0; ii < kLen; ii++) {
	    //   XTResid[ii] = XTY[ii];
      // }
       nntt = 0;
    
    //if (NLen < kLen) {
    //F77_CALL(dscal)(&kLen,&ZeroD, InvSXX, &One);
    double Neggy  = 0.0;   One = 1;
    double ACount = 0.0;
    for (jj = 0; jj < kLen; jj++) {
      if (OnBeta[jj] != 0.0) {
        if (XLC[jj] < 0 || XLC[jj] >= kLen) {
          if (iiWeights == NULL) {
           for (int kk = 0; kk < kLen; kk++) {
             ACount = 0.0;
             for (int it = 0; it < NLen; it++) {
                 ACount +=  XX[it+NLen*jj] * XX[it+NLen * kk];
             }
             XTResid[kk] -= ACount * OnBeta[jj];
           }
          } else {
           for (int kk = 0; kk < kLen; kk++) {
             ACount = 0.0;
             for (int it = 0; it < NLen; it++) {
                 ACount +=  XX[it+NLen*jj] * XX[it+NLen * kk] * iiWeights[it];
             }
             XTResid[kk] -= ACount * OnBeta[jj];
           }          
          }
        } else {
          Neggy = -1.0 * OnBeta[jj];
          if (iiWeights != NULL && iWeightedXtX != NULL && iWeightedXtX[XLC[jj]] <= 0) {
            ReweightCoordinate(jj);
          }
          F77_CALL(daxpy)(&kLen, &Neggy, pXTX[XLC[jj]], &One, XTResid, &One);
        }
      }
    }
    //  for (jj = 0; jj < OnKappaS; jj++)  {
    //    if (kXFinder[jj] < 0 || kXFinder[jj] >= kLen) {
    //      Rf_error("ModifyXTResid: Error: kXFinder[jj=%d] = %d! \n", jj,
    //        kXFinder[jj]);
    //    }
    //    Neggy = -1 * OnBeta[kXFinder[jj]];
    //    F77_CALL(daxpy)(&kLen, &Neggy, XTX + jj * kLen, &One, 
    //          XTResid, &One);
    //  } 
      
    //Rprintf((char*) "CDO->MakeXTResid() Onto XTX and OnBeta with 2\n");
    //R_FlushConsole(); R_ProcessEvents();       
     //  for (jj = 0; jj < OnKappaS; jj++) {
	   //    for (ii = 0; ii < kLen; ii++) {
		 //     XTResid[ii] -= XTX[nntt] * OnBeta[ kXFinder[jj]]; nntt++;
	   //    }
     // }
    //Rprintf((char*) "CDO->MakeXTResid() Finished from it \n");
    //R_FlushConsole(); R_ProcessEvents();       
  }
  ModifiedXTResid = 0;
  return(1);
}

double GetDMax(double *dB, int L) {
  double MyM = R_isnancpp(dB[0]) ? -999.9 : dB[0];
  if (L <= 1) { return(MyM); }
  for (int ii = 1; ii < L; ii++) {
    if (!R_isnancpp(dB[ii]) && MyM < dB[ii]) {
      MyM = dB[ii];
    }
  }
  return(MyM);
}
double GetDMin(double *dB, int L) {
  double MyM = R_isnancpp(dB[0]) ? -999.9 : dB[0];
  if (L <= 1) { return(MyM); }
  for (int ii = 1; ii < L; ii++) {
    if (!R_isnancpp(dB[ii]) && MyM > dB[ii]) {
      MyM = dB[ii];
    }
  }
  return(MyM);
}
int CountNan(double *dB, int L) {
  int CN = 0;
  for (int ii = 0; ii < L; ii++) {
    if (R_isnancpp(dB[ii])) {
      CN++;
    }
  }
  return(CN);
}
int GLMBinomObject::RunGLMLasso() {
  if (PrintFlagGLMB > 0) {
     Rprintf("RunGLMLasso: Ready to start Algorithm, first set pIter.\n"); R_FlushConsole();
  }
  int iiter=0;int jj;
  pIter[0] = iiter+1;
  double deltamove;            int One = 1;

  if (PrintFlagGLMB > 0) {
    Rprintf("Initial y01 vector: "); PrintVector(y01, NLen); 
    Rprintf("\n"); R_FlushConsole();
    double MySelf = 0.0;
    int One = 1;
    MySelf = F77_CALL(dasum)(&p, OnBeta, &One);
    Rprintf("OnBeta has sum %f, Beta0=%f and rest is: ", MySelf, pBeta0[0]); PrintVector(OnBeta, kLen);
    Rprintf("\n"); R_FlushConsole();
  }
  if (PrintFlagGLMB > 0) {
     Rprintf("RunGLMLasso: Starting by getting Weights right\n");
     R_FlushConsole();
  }
  this->RefreshCurrentLogitProb();
  if (PrintFlagGLMB > 1) {
     Rprintf("RunGLMLasso: Refreshing Current Logit Prob\n");
     Rprintf("  Here is the fitted Log Probabilities. \n");
     PrintVector(nlogitCurrentProb, 5);  R_FlushConsole();
     R_FlushConsole();
  } else if (PrintFlagGLMB > 0) {   
    Rprintf("  --  max Log Probabilities = %f, min is %f, Nan = %d\n", GetDMax(nlogitCurrentProb, NLen),
      GetDMin(nlogitCurrentProb, NLen), CountNan(nlogitCurrentProb, NLen));
    R_FlushConsole();
  }
  this->RefreshGLMWeights();
  if (PrintFlagGLMB > 1) {
    Rprintf("Initial nlogitCurrentProb vector, MaxLoops = %d, ", MaxLoops); 
    PrintVector(nlogitCurrentProb, NLen);
    Rprintf("Initial Z vector: "); PrintVector(Z, NLen);
    Rprintf("Beta0 now is %f \n", pBeta0[0]);
    Rprintf("Here are weights. :"); PrintVector(iiWeights, NLen);
    Rprintf("\n"); R_FlushConsole();
    Rprintf(" We are running to test integrity. \n");
    DataGLMIntegrity(1);
    Rprintf(" Guess we passed integrity. \n"); R_FlushConsole();
  }  
  int OldPrint = 0;
  GLMTotMove = 0.0;
  for (iiter = 0; iiter < MaxLoops;iiter++) {
    pIter[0] = iiter+1;
    if (PrintFlagGLMB > 2) {
      Rprintf("RunGLMLasso:  OuterLoop iiter = %d/%d\n", iiter, MaxLoops);
      Rprintf("  Now integrity check: ");
      DataGLMIntegrity(0);
      Rprintf("    Guess we passed Integrity. Check. \n"); R_FlushConsole();
      //if (fabs(OnBeta[0]) > 1000.0) {
      //  Rprintf("RunGLMLasso: Really bad, why is OnBeta so bad? \n"); R_FlushConsole();
      //  Rf_error("Error in RunGLMLasso. \n");
      //}
      R_FlushConsole();
    }
    pBeta0prev[0]=pBeta0[0];
    for (jj =0;jj<p;jj++) {
      pBetasPrev[jj] = OnBeta[jj]; 
    }
    for (jj=0;jj<n;jj++) {
      nlogitPrevProb[jj] = nlogitCurrentProb[jj];
    }
    if (PrintFlagGLMB > 1) {
      Rprintf("RunGLMLasso: About to Run CDO Algorithm %d/%d. \n", iiter, MaxLoops);
      R_FlushConsole();
    }       
    if (PrintFlagGLMB >= 3) {
      OldPrint = this->PrintFlag;
      this->PrintFlag = PrintFlagGLMB;
    }
    //this->PrintFlag = 10;
    //if (this->PrintFlag >= 5) {
    //  Rprintf("RunGLMLasso: I have set the printflag on this. \n"); R_FlushConsole();
    //}
    this->OnLoop = 0;  this->OnCoord = 0;
    //if (this->PrintFlag >= 5 || iiter <= 2) {
    //  Rprintf("RunGLMLasso: OnLoop = %d, OnCoord=%d \n",
    //    this->OnLoop, this->OnCoord); R_FlushConsole();
    //}
    this->RunAlgorithm(MaxLoops, MaxEpsilon);
    if (PrintFlagGLMB >= 3) {
      this->PrintFlag = OldPrint;
      Rprintf("RunGLMLasso: We finished RunAlgorithm. \n"); R_FlushConsole();
    }
    if (GetDMax(OnBeta, NLen) >= 99999999.9 || GetDMin(OnBeta, NLen) <= -999999.9) {
      Rprintf(" -----------------------------------------------------------------------------------------\n"); R_FlushConsole();
      Rprintf("RunGLMLasso: Oh No, we have problematic Beta!\n"); R_FlushConsole();
      for (int ii = 0; ii < kLen; ii++) {
        if (ii % 10 == 0) { Rprintf("\n%d:", ii); } R_FlushConsole();
        Rprintf("%.1f,", OnBeta[ii]); R_FlushConsole();
      }
      Rprintf("\n"); R_FlushConsole();
      Rprintf("  --  Here was On Gammas \n"); R_FlushConsole();
      for (int ii = 0; ii < kLen; ii++) {
        if (ii % 10 == 0) { Rprintf("\n%d:", ii); } R_FlushConsole();
        Rprintf("%.1f,", OnGammas[ii]); R_FlushConsole();
      }
      Rprintf("\n"); R_FlushConsole();
      Rprintf(" -- Here was nlogitCurrentProb \n"); R_FlushConsole();
      for (int ii = 0; ii < NLen; ii++) {
        if (ii % 10 == 0) { Rprintf("\n%d:", ii); } R_FlushConsole();
        Rprintf("%.1f,", nlogitCurrentProb[ii]); R_FlushConsole();
      }      
      Rprintf(" --- \n"); R_FlushConsole();
      SuccessFlag = -1;
      Rprintf(" ---  You're going to have to quit and figure out what went wrong. \n"); R_FlushConsole();
      Rprintf(" -----------------------------------------------------------------------------------------\n");
      R_FlushConsole();
      return(-1);
      
    }
    GLMTotMove += TotMove;
    TotMove = 0.0;
    if (PrintFlagGLMB > 2) {
      Rprintf("RunGLMLasso:  OuterLoop iiter = %d/%d - ", iiter, MaxLoops);
      Rprintf("  Survived CDO Algorithm, Doing Reweighting\n");
      R_FlushConsole();
    }
    this->RefreshCurrentLogitProb();
    if (PrintFlagGLMB > 2) {
      Rprintf("RunGLMLasso: About to Refresh GLMWeights\n");
    }
    this->RefreshGLMWeights();
       deltamove = 0.0;
       deltamove += fabs(pBeta0prev[0] -pBeta0[0]);
       for (jj =0;jj< p; jj++) {
           deltamove += fabs(pBetasPrev[0] -OnBeta[0]);
       }
       for (jj = 0; jj< n; jj++) {
          deltamove+= fabs(nlogitCurrentProb[jj] - nlogitPrevProb[jj]);
       }
      if (PrintFlagGLMB > 2) {
          Rprintf("RunGLMLasso:  OuterLoop iiter = %d/%d - ", iiter, MaxLoops);
          Rprintf("  deltamove = %.7f  but MaxEpsilon = %.7f\n",
             deltamove, MaxEpsilon);
          R_FlushConsole();
       }  
       if (PrintFlagGLMB > 3) {
         Rprintf("  Current OnBeta is \n");
           PrintVector(OnBeta, p);
       }
       if (pRecordOnBetas != NULL) {
         F77_CALL(dcopy)(&p,OnBeta, &One, pRecordOnBetas + p *iiter, &One);
       }
       //if (pBeta0Records != NULL  && (int) ceiling(LengthRecordOnBetas / ((double)kLen)) > iiter) {
       //   pBeta0Records[iiter] = pBeta0[0];
       //}
       if (deltamove <=MaxEpsilon){
          if (PrintFlagGLMB > 0) {
           Rprintf("RunGLMLasso, ending at deltamove=%f \n", deltamove); R_FlushConsole();
          }
          return(1);
       }

  }
  if (PrintFlagGLMB > 0) {
    Rprintf("RunGLMLasso, ending at a maximum load in. \n"); R_FlushConsole();
  }
  return(0);
}

int GLMBinomObject::RefreshGLMWeights() {
  // double ZeroD = 0.0; int One = 1;  double OneD = 1.0;
  // int Zero = 0;
  if (PrintFlagGLMB > 2) {
    Rprintf("GLMBinomObject: RefreshingGLMWeights\n");
    R_FlushConsole();
  }
  int ii;
  //int nntt = 0;
  double MaxIf=0;
  for (ii = 0; ii < n; ii++) {
       if (nlogitCurrentProb[ii] < -30.0) {
         if (nlogitCurrentProb[ii] < -100.0) {
            MaxIf = -100.0;
         }  else {
            MaxIf = nlogitCurrentProb[ii];
         }
         if (y01[ii] == 0) {
            nProbWeights[ii] = exp(MaxIf);
            this->Z[ii] = nlogitCurrentProb[ii] -(1.0 + exp(MaxIf));
         } else {
            nProbWeights[ii] = exp(MaxIf);
            this->Z[ii] = nlogitCurrentProb[ii] + exp(-1.0 * MaxIf);
         }
       } else if (nlogitCurrentProb[ii] > 30.0) {
        if (nlogitCurrentProb[ii] > 100.0) {
            MaxIf = 100.0;
         }  else {
            MaxIf = nlogitCurrentProb[ii];
         }
          if (y01[ii] == 1) {
            nProbWeights[ii] = exp(-MaxIf);
            this->Z[ii] = nlogitCurrentProb[ii] + (1.0 - exp(-MaxIf));
          } else {
            nProbWeights[ii] = exp(-MaxIf);
            this->Z[ii] = nlogitCurrentProb[ii] - exp(MaxIf);
          }    
    } else { 
      MaxIf =  nlogitCurrentProb[ii];
      nProbWeights[ii]  = exp(MaxIf);
      nProbWeights[ii] = nProbWeights[ii] / ( 1.0 + nProbWeights[ii]);
      Z[ii] = nlogitCurrentProb[ii] + (y01[ii] - nProbWeights[ii])/
        (nProbWeights[ii] * (1.0 - nProbWeights[ii]));
      nProbWeights[ii] = (nProbWeights[ii] * (1.0 - nProbWeights[ii]));   
    }
  }
  /*
  double nPW = 0.0;
  for (ii = 0; ii < n; ii++) {
    nPW += nProbWeights[ii];
  }
  int One = 1; double dS = NLen / nPW;
  F77_CALL(dscal)(&n, &dS, nProbWeights, &One);
  */
  if (PrintFlagGLMB > 4) {
    Rprintf(" nProbWeights: MaxIf = %f, Now : ", MaxIf); PrintVector(nProbWeights, n);
    Rprintf(" nProbWeights and Z are inputed, now to refresh weights\n");
    R_FlushConsole();
  }
  this->RefreshWeights(nProbWeights);
  if (pBeta0 != NULL) {
   pBeta0[0] = CalculateBeta0();
  }
  return(1);
}

double GLMBinomObject::CalculateBeta0() {
  SumWeights = 0.0;
  SumYWeights = 0.0;
  for (int ii = 0; ii < n; ii++) {
    SumWeights += nProbWeights[ii];
  }
  for (int ii = 0; ii < n; ii++) {
    SumYWeights += Z[ii] * nProbWeights[ii];
  }
  return(SumYWeights/SumWeights);
}

////////////////////////////////////////////////////////////////////////////////
//  GLMBinomObject():: RefreshCurrentLogitProb
//
//  A binomial object useful for doing logistic regression with Lasso penalties.
//
int GLMBinomObject::RefreshCurrentLogitProb() {
    if (PrintFlagGLMB > 2) {
      Rprintf("GLMBinomObject:RefreshCurrengLogitProb() running\n");
      R_FlushConsole();
    } 
    //double ZeroD = 0.0; 
    int One = 1;  double OneD = 1.0;
    //int Zero = 0;
    int ii;     // int DontJump = 0;
    //int nntt = 0;
    if (PrintFlagGLMB > 4) {
      Rprintf("GLMBinomObject: nlogitCurrentProb going to scale to n=%d\n",n);
      R_FlushConsole();
    }
    double Cp = 0.0;
 
    if (pBeta0 != NULL) {
    for(ii = 0; ii < n; ii++)  {      
      // Rprintf("ii=%d, ", ii); R_FlushConsole();
      // if ( DontJump >= 10 ) {
      //     Rprintf("\n"); R_FlushConsole();
      //     DontJump = 0;
      // }
      // if (ii >= n) {
      //  Rprintf(" Effed Up Shit GoingDown, why are we doing this?, ii= %d", ii);
      //  R_FlushConsole(); SuccessFlag = -1; return(-1);
      // }
       nlogitCurrentProb[ii] = pBeta0[0];
      // DontJump++;
    }}
     
    if (PrintFlagGLMB > 4) {
      Rprintf("GLMBinomObject: XTXFlag about to be looked at %d\n", XTXFlag);
      R_FlushConsole();
    } 
    if (this->XTXFlag == 2) {
      if (PrintFlagGLMB > 4) {
        Rprintf("GLMBinomObject: Doing a Multiplication, OnKappaS = %d\n", OnKappaS);
        Rprintf("kXFinder[0] = %d, or Beta = %f, so this should be fun\n", 
          kXFinder[0], this->OnBeta[kXFinder[0]]);
        Rprintf("pBeta0[0] = %f\n", pBeta0[0]);
        Rprintf("BeforeMultiplication: nlogit: ");
        PrintVector(nlogitCurrentProb, n); Rprintf("\n"); R_FlushConsole();
        R_FlushConsole();
      } 
      if (this->OnKappaS > 0) {
        for (ii = 0; ii < this->OnKappaS;ii++) {
             F77_CALL(daxpy)(&n, this->OnBeta + kXFinder[ii],
		              this->XX + kXFinder[ii] * n, &One,
	                nlogitCurrentProb, &One);
        }        
      }
      if (PrintFlagGLMB > 4) {
        Rprintf("GLMBinomObject: After Multiplication, OnKappaS = %d\n", OnKappaS);
        Rprintf("After Multiplication: nlogit: ");
        PrintVector(nlogitCurrentProb, n); Rprintf("\n"); R_FlushConsole();
        R_FlushConsole();
      } 
    } else {
      if (PrintFlagGLMB > 4) {
        Rprintf("GLMBinomObject: XTXFlag is not 2 HardWay \n");
        R_FlushConsole();
      } 
      /* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y, */
      F77_NAME(dgemv)("N", &n, &p,
		      &OneD, this->XX, &n, this->OnBeta,
          &One, &OneD, nlogitCurrentProb, &One);
    
    } 
    for (int ii = 0; ii < p; ii++) {
      Cp += fabs(nlogitCurrentProb[ii] / ((double) p));
    }
    int UpdateCp = 0;
    if (Cp >= 1000.0) {
      Rprintf("----------------------------------------------------------------------------------\n");
      Rprintf("GLMBinomObject(OnLoop=%d): We will update DiagShrinkage = %f, because Cp is out of line. \n", OnLoop, DiagShrinkage); R_FlushConsole();
      Rprintf("GLMBinomObject(OnLoop=%d), we get a Cp unsustainable point = %f, TestCDOs. \n", OnLoop, Cp);
      //TestCDO();
      //Rprintf(" --- Now Check DataIntegirty \n"); DataIntegrity(0);
      R_FlushConsole();
    }
    while (Cp >= 1000.0) {
      UpdateCp++;
      double AvOnGammas = 0.0;
      if (OnGammas != NULL) {
      for (int ii = 0; ii < p; ii++) {
        AvOnGammas += OnGammas[ii];
      }
       AvOnGammas = AvOnGammas / ((double) p);
      } else {
        AvOnGammas = OnGamma;
      }
      Rprintf("GLMBinomObject(OnLoop=%d): Issue, sum of Betas is %f, unsustainable.  Av(OnGammas) = %f, Increasing Lambda2 to %f, factor update=%d\n", 
        OnLoop, Cp, AvOnGammas, DiagShrinkage, UpdateCp);
      R_FlushConsole();
      double OldDiag = DiagShrinkage;
      if (DiagShrinkage <= 0.0) {
        DiagShrinkage = 0.0001;
      } else {
        DiagShrinkage = DiagShrinkage * 2.0;
      }
      for (int ii = 0; ii < p; ii++) {
        InvSXX[ii] = 1.0/(1.0/InvSXX[ii] + DiagShrinkage - OldDiag);
      }
      double MyShrink = 1.0 / (1.0+DiagShrinkage);
      for (int ii = 0; ii < p; ii++) {
        if (this->OnBeta[ii] != 0.0) {
          this->OnBeta[ii] *= MyShrink;
        }
      }
      for (int ii = 0; ii < n; ii++) {
        nlogitCurrentProb[ii] *= MyShrink;
      }
      Cp = Cp * MyShrink;
    }
    if (UpdateCp >= 1) {
      for (int ii = 0; ii < p; ii++) {
        XTResid[ii] = XTY[ii];
        for (int jj = 0; jj < OnKappaS; jj++) {
          if (OnBeta[kXFinder[jj]] != 0.0) {
            XTResid[ii] -= *(pXTX[jj]+ii) * OnBeta[kXFinder[jj]];
          }
        }
      }
      //Rprintf("GLMBinomObject: XTResid is updated now let's check Data Integrity.");
      DataIntegrity(0);
    }
    if (PrintFlagGLMB > 4) {
      Rprintf("GLMBinomObject:: Printing first 5 members of nlogitCurrentProb\n");
      PrintVector(nlogitCurrentProb, 5);  R_FlushConsole();
      Rprintf("GLMBinomObject:: Last Element = %.6f\n", nlogitCurrentProb);
      R_FlushConsole();
    }
    return(1);       
}

///////////////////////////////////////////////////////////////////////////////
//  GLMLasso(sX, sy01, sOnBeta)
//
//  A function for constructing and performing GLM generalized 
//   linear Logit regression
//   Style Lasso using CoordinateDescentLasso.cc
//
//
SEXP GLMLasso(SEXP sX, SEXP sy01, SEXP sOnBeta, SEXP sBeta0,
   SEXP sBeta0prev, SEXP sBetasPrev, SEXP sOnGamma,
   SEXP sRecordOnBetas,SEXP sBeta0Records,
   SEXP sOnGammas,
   SEXP sInitKKs, SEXP sInverseGammaConstant,
   SEXP snProbWeights, SEXP slogitCurrentProb, 
   SEXP snlogitPrevProb, SEXP sz, SEXP sPrintFlag,
   SEXP sMaxCDOEpsilon, SEXP sMaxCDORuns,
   SEXP sIter) {
  SEXP dimsX = GET_DIM(sX); int n,p;
  int PrintFlag1 =0;
  if (isInteger(sPrintFlag)) { PrintFlag1 = INTEGER(sPrintFlag)[0];} else {
     PrintFlag1 = (int) REAL(sPrintFlag)[0];
  }
  if (PrintFlag1 > 1) {
     Rprintf("Taking Dimension of dimsX\n"); R_FlushConsole();
  }
  if (isInteger(dimsX)){
       n = INTEGER(dimsX)[0];  
       p = INTEGER(dimsX)[1];
  } else{
       n = (int)REAL(dimsX)[0];
       p = (int)REAL(dimsX)[1];
  }
  if (PrintFlag1 > 4) {
     Rprintf("Dimensions are n=%d, p=%d\n", n, p); R_FlushConsole();
  }

  int yLen = (int) length(sy01);  
  if (yLen != n) { 
    Rprintf("yLen = %d, n = %d, not same length, reject!\n", yLen, n); 
    R_FlushConsole();
    return(sy01);
  }
  if (length(sz) < n) {
     Rprintf("Unfortunately sz is not big enough %d, n = %d\n",
        length(sz), n);R_FlushConsole();
  }
  int *py01;  int Ispy01Int = 0;
  if (PrintFlag1 > 4) {
     Rprintf(" About to take length of sy01 \n"); R_FlushConsole();
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // Unfortunately we have to convert Reals to Integers if we are unsure of
  //   SEXP type being delivered.
  //
  // It is often difficult to get user to deliver Integer SEXP
  if (isInteger(sy01)) {  
    py01 = INTEGER(sy01);
  }  else {
    Ispy01Int = 1;
    py01 = (int*) Calloc(length(sy01), int);
    for (int iji = 0; iji < length(sy01); iji++) {
      py01[iji] = (int) REAL(sy01)[iji];
    }
  }
  int BetasLen = length(sOnBeta);
  if (BetasLen < p) {
    Rprintf("BetasLen = %d, p = %d, not same length, reject!\n", BetasLen, p);
    R_FlushConsole(); return(sOnBeta);
  }
  if (length(sBetasPrev) < p) {
    Rprintf("sBetasprev not long enough=%d,not %d\n", 
    length(sBetasPrev), p);  R_FlushConsole();return(sBetasPrev);
  }  
  double *pBetasPrev = REAL(sBetasPrev);  
  double *pBeta0prev =REAL(sBeta0prev);
  if (length(sOnGammas) != p) {
    Rprintf(" sOnGammaslength = %d, too short to p = %d\n", length(sOnGammas),p);
    R_FlushConsole(); return(sOnGammas);
  }
  
  double *pRecordOnBetas =NULL; 
  double *pBeta0Records=NULL;
  int OuterNumLoops;
  if (PrintFlag1 > 4) {
    Rprintf("About to take MaxCDORuns \n"); R_FlushConsole();
  }
  if (isInteger(sMaxCDORuns)) {
     OuterNumLoops = INTEGER(sMaxCDORuns)[0];
  } else {
     OuterNumLoops = (int) REAL(sMaxCDORuns)[0];
  }
  if (PrintFlag1 > 4) {
     Rprintf("  OuterNumLoops = %d\n", OuterNumLoops); R_FlushConsole();
  }

  if (length(sRecordOnBetas) < p * OuterNumLoops) {
    if (length(sRecordOnBetas) <= 1) {
    } else {
       Rprintf("sRecordOnBetas not long enough, turning to NULL \n");
       R_FlushConsole();
    }
    pRecordOnBetas = NULL;
    ///Rprintf(" Fuck You Not Enough Data for Record Betas\n"); R_FlushConsole(); 
    //return(sRecordOnBetas);
  } else {
    pRecordOnBetas=REAL(sRecordOnBetas); 
    pRecordOnBetas[0] = 100.0;  pRecordOnBetas[p + 20] = -999;
    pRecordOnBetas[p * OuterNumLoops - 1] = -10;
    //Rprintf("  At Least I have Record OnBetas \n"); R_FlushConsole();
    //return(sOnBeta);
  }  
  if (Rf_isNull(sBeta0Records)) {
    pBeta0Records = NULL;
  } else if (Rf_length(sBeta0Records) < OuterNumLoops) {
    if (length(sBeta0Records) <= 1) {
    } else {
       //Rprintf("sBeta0Records not long enough, turning to NULL \n");
       // R_FlushConsole();
    }
    pBeta0Records = NULL;
  } else {pBeta0Records=REAL(sBeta0Records);}
  if (PrintFlag1 >4) {
     Rprintf("Looking at snlogitPrevProb \n"); R_FlushConsole();
  }
  if (length(snlogitPrevProb)< n){
    Rprintf("Error, length(snlogitPrevProb)=%d, smaller than n=%d\n",
      length(snlogitPrevProb), n); R_FlushConsole();  return(snlogitPrevProb);
  }         
  int PrintFlag;
  if (isInteger(sPrintFlag)) {
     PrintFlag = INTEGER(sPrintFlag)[0];
  } else {
     PrintFlag = (int) REAL(sPrintFlag)[0];
  }     
  if (PrintFlag > 1) {
    Rprintf("GLMLasso function: Loaded up to PrintFlag\n"); R_FlushConsole();
  }
  if (PrintFlag1 > 4) {
    Rprintf("MaxCDOEpsilon, taking REAL \n"); R_FlushConsole();
  }
  double MaxCDOEpsilon = REAL(sMaxCDOEpsilon)[0]; 
  int MaxCDORuns;     
  if (isInteger(sMaxCDORuns)) {
    MaxCDORuns = INTEGER(sMaxCDORuns)[0];
  } else {
    MaxCDORuns = (int) REAL(sMaxCDORuns)[0];
  }

  if (PrintFlag > 4) {
    Rprintf("InitKKs, about to measure\n"); R_FlushConsole();
  }
  int InitKKs;
  if (isInteger(sInitKKs)) {
    InitKKs = (int) INTEGER(sInitKKs)[0];
  } else {
    InitKKs = (int) REAL(sInitKKs)[0];
  }
  if (PrintFlag > 4) {
    Rprintf("About to create new GLMO \n"); R_FlushConsole();
  }
  int *pIter; int MakepIter;
  if (isInteger(sIter)) {
    pIter = INTEGER(sIter); MakepIter = 0;
  } else {
    pIter = (int*) Calloc(2, int); MakepIter = 1;
  }
  GLMBinomObject *GLMO = new GLMBinomObject(n, p,(double*)REAL(sX), (int*)py01, 
    REAL(sOnBeta),  REAL(sBeta0),  pBeta0prev,
    (double*) pBetasPrev, (double) REAL(sOnGammas)[0],
    pRecordOnBetas,pBeta0Records, (double *) REAL(sOnGammas),
    (int) InitKKs, REAL(sInverseGammaConstant)[0], 
    REAL(snProbWeights), REAL(slogitCurrentProb), REAL(snlogitPrevProb),
    (double*) REAL(sz), (int) PrintFlag, (double) MaxCDOEpsilon, 
    (int) MaxCDORuns, pIter, -1);      
  GLMO->RunGLMLasso();
  if (Ispy01Int == 1) {
    Free(py01); py01 = NULL;  GLMO->Sety01Null();
  }
  if (MakepIter == 1) {
     REAL(sIter)[0] = (double) pIter[0];
     Free(pIter); pIter = NULL;
  }
  delete(GLMO);
  return(sOnBeta);
}

////////////////////////////////////////////////////////////////////////////////
// StatusReport()
//
// Print A Status Report and test the memory of the function 
//    for leaks or overwrites.
int CoordinateDescentObject::StatusReport() {
  Rprintf("CDO-> Status Report\n");
  Rprintf("OnGammas are: "); PrintVector(OnGammas, kLen); Rprintf("\n");
  Rprintf((char*) "Initial OnBeta Are : \n"); R_FlushConsole();
  PrintVector(OnBeta, kLen);  Rprintf("\n");
  Rprintf((char*) "CoordinateDescentObject:: Active Beta are : \n"); R_FlushConsole();
  int tti;
  if (OnKappaS < kLen && OnKappaS > 0 && kXFinder != NULL) {
    for (tti = 0; tti < OnKappaS; tti++) {
      Rprintf(" %d: %d  = %.4f \n", tti, kXFinder[tti], 
        OnBeta[kXFinder[tti]]);
        R_FlushConsole();
      }
  }
  R_FlushConsole(); R_ProcessEvents();
  return(1);
}
int GLMBinomObject::DataGLMIntegrity(int DoIn) {
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("DataGLMIntegrity: Start to check. \n"); R_FlushConsole();
  }
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is pIter good ? "); R_FlushConsole();
  }
  if (pIter != NULL) {
    pIter[0] = pIter[0] + 0.0;
  }
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf(" PASS. \n "); R_FlushConsole();
  }
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is intY good ? "); R_FlushConsole();
  }
  if (y01 != NULL) {
    for (int ii = 0; ii < n; ii++) {
      y01[ii] = y01[ii] + 0;
    }
  } else {
    Rprintf("  Hey, no y01 is NULL!\n"); R_FlushConsole();
  }
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf(" PASS. \n "); R_FlushConsole();
  }
  
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is nlogitCurrentProb good ? "); R_FlushConsole();
  }
  if (nlogitCurrentProb != NULL) {
    for (int ii = 0; ii < n; ii++) {
      nlogitCurrentProb[ii] = nlogitCurrentProb[ii] + 0.0;
    }
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf(" PASS. \n "); R_FlushConsole();
    }
  } else {
    Rprintf("  Hey, no nlogitCurrentProb is NULL!\n"); R_FlushConsole();
  }

  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is nlogitPrevProb good ? "); R_FlushConsole();
  }
  if (nlogitPrevProb != NULL) {
    for (int ii = 0; ii < n; ii++) {
      nlogitPrevProb[ii] = nlogitPrevProb[ii] + 0.0;
    }
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf(" PASS. \n "); R_FlushConsole();
    }
  } else {
    Rprintf("  Hey, no nlogitPrevProb is NULL!\n"); R_FlushConsole();
  }
    
  
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is nProbWeights good ? "); R_FlushConsole();
  }
  if (nProbWeights != NULL) {
    for (int ii = 0; ii < n; ii++) {
      nProbWeights[ii] = nProbWeights[ii] + 0.0;
    }
  } else {
    Rprintf("  Hey, no nProbWeights is NULL!\n"); R_FlushConsole();
  }
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf(" PASS. \n "); R_FlushConsole();
  }
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is Z good ? "); R_FlushConsole();
  }
  if (Z != NULL) {
    for (int ii = 0; ii < n; ii++) {
      Z[ii] = Z[ii] + 0.0;
    }
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf(" PASS. \n "); R_FlushConsole();
    }
  } else {
    Rprintf("  Hey, no Z is NULL!\n"); R_FlushConsole();
  }
  
           
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is pBetasPrev good ? "); R_FlushConsole();
  }
  if (pBetasPrev != NULL) {
    for (int ii = 0; ii < p; ii++) {
      pBetasPrev[ii] = pBetasPrev[ii] + 0.0;
    }
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf(" PASS. \n "); R_FlushConsole();
   }
  } else {
    Rprintf("  Hey, no pBetasPrev is NULL!\n"); R_FlushConsole();
  }
  
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is pBeta0prev good ? "); R_FlushConsole();
  }
  if (pBeta0prev != NULL) {
    pBeta0prev[0] = pBeta0prev[0] + 0.0;
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf(" PASS. \n "); R_FlushConsole();
   }
  } else {
    Rprintf("  Hey, no pBeta0Prev is NULL!\n"); R_FlushConsole();
  }     

    
  int OuterNumLoops = LengthRecordOnBetas;
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is pRecordOnBetas good ? "); R_FlushConsole();
  }
  if (pRecordOnBetas != NULL) {
    for (int iti = 0; iti < p * OuterNumLoops; iti++) {
     pRecordOnBetas[iti] = pRecordOnBetas[iti] + 0.0;
    }
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf(" PASS. \n "); R_FlushConsole();
   }
  } else {
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf("  Hey, no pRecordOnBetas is NULL!\n"); R_FlushConsole();
    }
  }     
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("--  Is pBeta0Records good ? "); R_FlushConsole();
  }
  if (pBeta0Records != NULL) {
    int ALen = LengthRecordOnBetas / kLen;
    if (ALen > 0) {
    //for (int iti = 0; iti < ALen; iti++) {
    // pBeta0Records[iti] = pBeta0Records[iti] + 0.0;
    //}
    }
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf(" PASS. \n "); R_FlushConsole();
   }
  } else {
    if (DoIn >= 1 || PrintFlagGLMB >= 2) {
      Rprintf("  Hey, no pBeta0Records is NULL!\n"); R_FlushConsole();
    }
  }    
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf("-- Check Overall Integrity \n");
  }
  DataIntegrity(DoIn);
  if (DoIn >= 1 || PrintFlagGLMB >= 2) {
    Rprintf(" -- They All Pass. \n"); R_FlushConsole();
  }                
  return(1);
}
int CoordinateDescentObject::DataIntegrity(int DoIn) {
  int DoP = DoIn;
  if (TriggerInvestigateFlag >= 2) {
    DoP = 1;
  }
  if (OnKappaMem < kLen) {
    if (pXTX[OnKappaMem] != NULL) {
      Rprintf("CDO->DataIntegrity, error, pXTX[OnKappaMem=%d] is not NULL!\n", OnKappaMem);
      return(-1);
    }
    if (OnKappaMem < kLen-1 && pXTX[kLen-1] != NULL) {
      Rprintf("CDO->DataIntegrity, error, pXTX[OnKappaMem=%d] is not NULL!\n", OnKappaMem);
      return(-1);   
    }
  }
  if (DoP >= 1) {
    Rprintf("CDO->Data Integrity, starting with OnGammas.\n");
    R_FlushConsole();
  }
  double tG = 0.0;
  if (OnGammas != NULL) {
    for (int ii = 0; ii < kLen; ii++) {
      tG += OnGammas[ii];
    }
  }
 if (DoP >= 1) {
    Rprintf("CDO->Data Integrity, starting with OnBeta.\n");
    R_FlushConsole();
  }
  double tB = 0.0;
  if (tB <= -1.0) {
  }
  if (OnBeta != NULL) {
    for (int ii = 0; ii < kLen; ii++) {
      tG += OnBeta[ii];
    }
  }
  if (DoP >= 1) {
    Rprintf("CDO->Data Integrity, starting XTResid\n"); R_FlushConsole();
  }
  double tXT = 0.0;
  if (XTResid != NULL) {
    for (int ii = 0; ii < kLen; ii++) {
      tXT += XTResid[ii];
    }
  }
  if (DoP >= 1) {
    Rprintf("CDO->Data Integrity checking kXFinder!\n");  R_FlushConsole();
  }
  double tkX = 0.0;
  if (kXFinder != NULL && XTXFlag >= 2) {
    for (int ii = 0; ii < OnKappaS; ii++) {
      if (kXFinder[ii] < 0 || kXFinder[ii] >= kLen) {
        Rprintf("DataIntegerity: What, ii=%d for OnKappaS=%d but flaw in kXFinder!\n",
          ii, OnKappaS); R_FlushConsole();
      }
      tkX += (double) kXFinder[ii];
    }
    if (DoP >= 1) {
      Rprintf("Checking Extended kXFinder!\n");  R_FlushConsole();
    }
    tkX = 0;
    for (int ii = OnKappaS; ii < OnKappaMem; ii++) {
      tkX += (double) kXFinder[ii];
    }
  }
  int RunO = 0;
  if (RunOrder != NULL) {
    if (DoP >= 1) {
      Rprintf("Checking RunOrder\n");   R_FlushConsole();
    }
    for (int ii = 0; ii < kLen; ii++) {
      RunO += RunOrder[ii];
    }
  }
  double tXTX = 0.0;
  if (kXFinder != NULL && XTXFlag >= 2) {
    for (int ii = 0; ii < OnKappaS; ii++) {
      for (int jj = 0; jj < kLen; jj++) {
        tXTX += (double) *(pXTX[ii]+jj);
      }
    }
    if (DoP >= 1) {
      Rprintf("Checking Extended XTX!\n"); R_FlushConsole();
    }
    tXTX = 0;
    for (int ii = OnKappaS; ii < OnKappaMem; ii++) {
      tkX += (double) kXFinder[ii];
    }
  }

  int CountBadInvSXX = 0;
  if (InvSXX == NULL) {
    Rprintf("DataIntegrity(): Hey, InvSXX is NULL!\n"); R_FlushConsole();
  } else {
    double ShouldBeXX;  double ShouldBeInvXX;
    int OnOne = 0; int Onk = 0;
    if (DoP >= 1) {
      Rprintf("DataIntegrity(): Checking InvSXX: XTXFlag = %d\n", XTXFlag); R_FlushConsole();
    }
    if (kXFinder != NULL && XTXFlag >= 2) {
      while(OnOne < kLen) {
        if (XLC[OnOne] < 0 || (iiWeights != NULL && iWeightedXtX != NULL && iWeightedXtX[XLC[OnOne]] <= 0)) {
          ShouldBeXX = 0.0;
          if (iiWeights == NULL) {
             Onk = OnOne * NLen;
             for (int ii = 0; ii < NLen; ii++) {
               ShouldBeXX += XX[Onk +ii] * XX[Onk + ii];
             }
          } else {
            Onk = OnOne * NLen;
            ShouldBeXX = 0.0;
            for (int ii = 0; ii < NLen; ii++) {
              ShouldBeXX += XX[Onk +ii] * XX[Onk + ii] * iiWeights[ii];
            }          
          }
        } else if (XLC[OnOne] >= 0 && XLC[OnOne] < OnKappaS) {
          ShouldBeXX = *(pXTX[XLC[OnOne]] + OnOne);
        }
        ShouldBeInvXX = 1.0 / (ShouldBeXX + DiagShrinkage);
        if (fabs(ShouldBeInvXX - InvSXX[OnOne]) >= .0001) {
          Rprintf("Error Here On OnOne = %d, ShouldBeInvXX=%f, InvSXX[%d] = %f, XLC[%d]=%d, ShouldBeXX=%f, DiagShrink=%f\n", OnOne, ShouldBeInvXX, OnOne, InvSXX[OnOne],
            OnOne, XLC[OnOne], ShouldBeXX, DiagShrinkage);
          R_FlushConsole();
          CountBadInvSXX++;
        }
        OnOne++;
      }
    } else {
      for (int ii = 0; ii < kLen; ii++) {
        ShouldBeXX = *(pXTX[ii] +ii);
        ShouldBeInvXX = 1.0/(ShouldBeXX + DiagShrinkage);
        if (fabs(ShouldBeInvXX - InvSXX[ii]) >= .0001) {
          CountBadInvSXX++;
        }
      }
    }
    if (CountBadInvSXX > 0) {
      Rprintf("ERROR: DataIntegrity: No, CountBadInvSXX = %d, not good\n", CountBadInvSXX);
      Rprintf(" -- You Should probably check out what is wrong with InvSXX. \n");
      CountBadInvSXX = 0;
      if (kXFinder != NULL && XTXFlag >= 2) {
       Rprintf("-- We are Running the Error Test again. \n"); R_FlushConsole();
      while(OnOne < kLen) {
        if (XLC[OnOne] < 0 || (iiWeights != NULL && iWeightedXtX != NULL && iWeightedXtX[XLC[OnOne]] <= 0)) {
          ShouldBeXX = 0.0;
          if (iiWeights == NULL) {
             Onk = OnOne * NLen;
             for (int ii = 0; ii < NLen; ii++) {
               ShouldBeXX += XX[Onk +ii] * XX[Onk + ii];
             }
          } else {
            Onk = OnOne * NLen;
             for (int ii = 0; ii < NLen; ii++) {
               ShouldBeXX += XX[Onk +ii] * XX[Onk + ii] * iiWeights[ii];
             }          
          }
        } else if (XLC[OnOne] >= 0 && XLC[OnOne] < OnKappaS) {
          ShouldBeXX = *(pXTX[XLC[OnOne]] + OnOne);
        }
        ShouldBeInvXX = 1.0 / (ShouldBeXX + DiagShrinkage);
        if (fabs(ShouldBeInvXX - InvSXX[OnOne]) >= .0001) {
         CountBadInvSXX++;
         Rprintf(" -- Error[%d], InvSXX[OnOne=%d] = %f, but we think ShouldBeInvXX=%f, ShouldBeXX=%f (XLC[%d] = %d)\n",
           CountBadInvSXX, OnOne, InvSXX[OnOne], ShouldBeInvXX, ShouldBeXX, OnOne, XLC[OnOne]); R_FlushConsole();
        }
        OnOne++;
      }
    } else {
      for (int ii = 0; ii < kLen; ii++) {
        ShouldBeXX = *(pXTX[ii] +ii);
        ShouldBeInvXX = 1.0/(ShouldBeXX + DiagShrinkage);
        if (fabs(ShouldBeInvXX - InvSXX[ii]) >= .0001) {
          CountBadInvSXX++;
          Rprintf(" -- Error[%d], InvSXX[OnOne=%d] = %f, but we think ShouldBeInvXX=%f, ShouldBeXX=%f (XLC[%d] = %d)\n",
           CountBadInvSXX, OnOne, InvSXX[OnOne], ShouldBeInvXX, ShouldBeXX, OnOne, XLC[OnOne]); R_FlushConsole();
        }
      }
    }
      Rprintf("Did We find any valuable errors? CountBadInvSXX = %d\n", CountBadInvSXX); R_FlushConsole();
      Rf_error("DataIntegrity: Flaw in InvSXX \n");
    }

  }
  SuccessDOP++;
  if (iiWeights != NULL && DoP >= 1) {
    Rprintf("Checking iiWeights. \n"); R_FlushConsole();
  }
  if (iiWeights != NULL && iiWeightsHere == 1) {
    double *NewiiWeights = Calloc(NLen, double);
    for (int ii = 0; ii < NLen; ii++) {
      NewiiWeights[ii] = iiWeights[ii];
    }
    Free(NewiiWeights);
    //iiWeights = NewiiWeights;
  }
  if (iWeightedXtX != NULL && DoP >= 1) {
    Rprintf("Checking iWeightedXtX. \n"); R_FlushConsole();
  }
  if (iWeightedXtX != NULL) {
    int *NewiWeightedXtX = Calloc(kLen, int);
    for (int ii = 0; ii < kLen; ii++) {
      NewiWeightedXtX[ii] = iWeightedXtX[ii];
    }
    Free(iWeightedXtX);
    iWeightedXtX = NewiWeightedXtX;
  }
  if (DoP >= 1) {
    Rprintf("Successful DOP Conclusion!"); R_FlushConsole();
  }
  return(1);
}



///////////////////////////////////////////////////////////////////////////
//  RunAlgorithm
//    Runs a setup Coordinate Descent object
//     until convergence (defined by MaxEpsilon)
//     or maximum number of loops (defined by MaxLoops)
int CoordinateDescentObject::RunAlgorithmDiagonostic(int MaxLoops, double MaxEpsilon) {
     /////////////////////////////////////////////////////////////////////
     //  Make sure that correct numbers for XTResid are located given
     //    beta estimate
     if (PrintFlag > 3) {
       Rprintf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");  R_FlushConsole();
       Rprintf("+++++ CoordinateDescentObject():RunAlgorithm \n");  R_FlushConsole();
       Rprintf("+++++      In CoordinateDescentLasso.cc \n");    R_FlushConsole();
       Rprintf("+++++   About to Run Algorithm, ModFlag = %f \n",
         (double) ModifiedXTResid); R_FlushConsole();  
     }
	   if (ModifiedXTResid > 0) {
	      if (PrintFlag > 3) {
          Rprintf("++++ CDO: RunAlgorithm: Must Modify XTResid\n");
          R_FlushConsole();
        }
		    MakeXTResid();
	      ModifiedXTResid = 0;	   
     }
     
    ///////////////////////////////////////////////////////////////////////
    // OnEp records current move, compares to MaxEpsilon after every loop 
	  OnEp = 0.0;
	  TotMove = 0.0;
    int MoreMove = 0;
    int AfterMoreMove = 0;
	  
	  int ii;
	  int tt = 0;
	  if (PrintFlag > 1) {
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, XTXFlag = %d\n", XTXFlag);
	       R_FlushConsole();
    }  
	  if (PrintFlag > 2) { // Echo
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, Starting with XTResid = \n");
	       PrintVector(this->XTResid, this->kLen);
	       R_FlushConsole();
    }
    if (PrintFlag > 1) {
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, PrintFlag = %d\n", (int) PrintFlag);
	       R_FlushConsole();
    }
    if (PrintFlag > 1) {
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, XTXFlag = %d\n", XTXFlag);
	       R_FlushConsole();
	       if (OnGammas == NULL) {
           Rprintf("++++                     We are doing fixed OnGamma version.\n"); R_FlushConsole();
         } else {
           Rprintf("++++                     We are doing OnGammas vector version.\n"); R_FlushConsole();
         }
    }
      
      //////////////////////////////////////////////////////
      // If kLen = 1 then solution should be immediate
    if (kLen == 1) { 
	      if (PrintFlag > 1) {
	       Rprintf((char*) "+++++ CDO:RunAlgorithm, kLen == 1, insta-solution\n", XTXFlag);
	       R_FlushConsole();
        }
	      double DopeGammas = OnGamma;
	      if (OnGammas != NULL) { DopeGammas = OnGammas[0]; }    
	      if (XTY[0] > DopeGammas) { OnBeta[0] = 
	              ( XTY[0] - DopeGammas ) / ((double) *(pXTX[0]));
          } else if (XTY[0] < - DopeGammas) { OnBeta[0] = (XTY[0] + DopeGammas) / *(pXTX[0]); 
          } else { OnBeta[0] = 0.0;
          } 
	      return(1);
   }
   int WasOnCoord;   
   int IntO = 0;
      //////////////////////////////////////////////////////////////////////////
      //  OnLoop through iterations for multiple coordinates
      for (OnLoop = 0; OnLoop < MaxLoops; OnLoop++) {
   	    if (PrintFlag > 1) {
	   	      Rprintf((char*) "+++++ CDO: RunAlgorithm, Starting Loop %d/%d\n",
	   	                     OnLoop, MaxLoops); 
            R_FlushConsole(); 
            R_ProcessEvents();
   	    }	
   	    if (PrintFlag > 2 && XTXFlag == 2) {
	           Rprintf((char*) "+++++ CDO:RunAlgorithm, XTXFlag = %d, kLen = %d, ",
                           XTXFlag, kLen);
             Rprintf((char*) " OnKappaS = %d, OnKappaMem = %d\n",
	                         OnKappaS, OnKappaMem); R_FlushConsole();
	           PrintRpMatrix(pXTX, kLen, OnKappaMem);
	           R_FlushConsole();
        }      
	      OnEp = 0;
	      //if (((int) OnCoord)/((int) 1000) * 1000 == OnCoord) {
          Rprintf("RunDataInt: OnLoop = %d/%d Beggining Check. -- "); R_FlushConsole();
            DataIntegrity(0);  Rprintf(" --- PASS \n"); R_FlushConsole();
        //}
	      ////////////////////////////////////////////////////////////////////////
	      // Loop through coordinates in order, looping OnCoord
	      for (ii = 0; ii < kLen; ii++) {
          if (ii >= 14000 && ii < 15000) {
            IntO = 0; 
          } else { IntO = 0; }
		      if (PrintFlag > 3) {
		   	      Rprintf((char*) "+++++ CDO: RunAlgorithm, Loop %d/%d, jj = %d, ",
                   OnLoop, MaxLoops, ii); R_FlushConsole();
              if (OnCoord < 0 || OnCoord > kLen) {
                Rprintf("+++++   Ooh OnCoord is bad!\n"); R_FlushConsole();
              }  else {
                Rprintf((char*) " OnBeta[%d] = %.4f\n",
		   	                OnCoord, (double) OnBeta[OnCoord]); 
              }
              R_FlushConsole(); R_ProcessEvents();
	   	    }
          if (TriggerInvestigateFlag >= 1) {
            DataIntegrity(TriggerInvestigateFlag-1);
          }
	   	    
	   	    RecentMove = 0;   WasOnCoord = OnCoord;
          if ( (((int) OnCoord) / ((int) 1000)) * 1000 == OnCoord || IntO == 1) {
            Rprintf(" -- OnLoop = %d/%d, OnCoord = %d/%d, Mem =%d, OnKS=%d --  Try UpDateCoord(), WasOnCoord=%d. \n",
              OnLoop, MaxLoops, OnCoord, kLen, OnKappaMem, OnKappaS, WasOnCoord); R_FlushConsole();
          }
	        UpDateCoord();  // Update coordinate OnCoord
          if (MemoryFailFlag == 1) {
            Rprintf("CoordinateDescentLasso.cc::RunAlgorithmDiagonostic() -- UpdateCoord returned a MemoryFailFlag=%d. quit. \n", MemoryFailFlag); R_FlushConsole();
            return(-1);
          } 
          if (SuccessFlag == -1) {
            Rprintf("CoordinateDescentLasso.cc::RunAlgorithmDiagonostic() -- UpdateCoord returned a SuccessFlag=%d. quit. \n", SuccessFlag); R_FlushConsole();
            return(-1);
          }
          if ( (((int) OnCoord) / ((int) 1000)) * 1000 == OnCoord || IntO == 1) {
            Rprintf(". PASS Update (Move=%f), ", RecentMove); R_FlushConsole();
            Rprintf("OnKappaS=%d/%d, ", OnKappaS, OnKappaMem);  R_FlushConsole();
          }
	        OnEp += fabs(RecentMove); TotMove +=fabs(RecentMove);
	   	    if (PrintFlag > 3) {
		   	      Rprintf((char*) "+++++ CDO: RunAlgorithm, Loop %d/%d, jj = %d,",
                      OnLoop, MaxLoops, OnCoord);
              Rprintf((char*) " Move was %.4f\n",
		   	              RecentMove); 
              R_FlushConsole(); R_ProcessEvents();
	   	    }	  
          if (((int) OnCoord / (int) 1000) * 1000 == OnCoord || IntO == 1 ) {
           Rprintf("Move(%d) = %f, Check Integrity: ", OnCoord, RecentMove); R_FlushConsole();
             DataIntegrity(0);  Rprintf("   PASS ----\n"); R_FlushConsole();  
          } else if (RecentMove > 0 && MoreMove < 200) {
            Rprintf("  --  Move on OnLoop = %d/%d, OnCoord = %d/%d move %f for total %d moved ",
              OnLoop, MaxLoops, OnCoord, kLen, RecentMove, MoreMove); R_FlushConsole();
            Rprintf("OnKappaS=%d/%d, ", OnKappaS, OnKappaMem);
            Rprintf("Check Integrity: "); R_FlushConsole();
              DataIntegrity(0);  Rprintf("   PASS ----\n"); R_FlushConsole();     
            MoreMove++;       
          }	else if (RecentMove > 0) {
            if (AfterMoreMove >= 200) {
              Rprintf("  --  Move on OnLoop = %d/%d, OnCoord = %d/%d move %f for total %d moved ",
                OnLoop, MaxLoops, OnCoord, kLen, RecentMove, MoreMove); R_FlushConsole();
              Rprintf("OnKappaS=%d/%d, ", OnKappaS, OnKappaMem);
              Rprintf("Check Integrity: "); R_FlushConsole();
                DataIntegrity(0);  Rprintf("   PASS ----\n"); R_FlushConsole();     
              MoreMove++;   
              AfterMoreMove=0;
            } else {
              DataIntegrity(0);
            }
            MoreMove++;  AfterMoreMove++;
          }        	      
	   	    MoveCoord(ii);  //Choose next Coordinate in Order
   	    }
      
   	if (OnRecordBetas != NULL) {
	   	for (ii = 0; ii < kLen; ii++) {
		    OnRecordBetas[tt] = OnBeta[ii]; tt++;
	   	}
   	}   
   	if (MaxEpsilon >= 0 && OnEp <= MaxEpsilon) {
      if (PrintFlag > 0) {
        Rprintf("+++++ CDO:RunAlgorithm: MaxEpsilon produces a return, OnEp = %f < %f, OnLoop = %d\n", OnEp, MaxEpsilon, OnLoop);
        R_FlushConsole();
      }
	   	return(1);
   	}
  }
  if (PrintFlag > 0) {
	  Rprintf((char*) "+++++ CDO:RunAlgorithm, Finished\n", XTXFlag);
    R_FlushConsole();
  }
  return(1);
}

