
#ifndef LARSOBJECTDD
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
#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif

////////////////////////////////////////////////////////////////////////
////////////////////////  External Code After Here:


////////////////////////////////////////////////////////////////
//  LarsLassoAlgorithmBI
//
//   This directs the program, constructs a LARSObject, which runs the lars lasso
//    algorithm.  
//
int LarsLassoAlgorithmBI( double * yys, int NLen, double *xxs, int kLen, 
                      double *InitBetas, 
                      double *OldBetas, double lambda,double *OutBetas, double *BetasKeep,
                      int *JKeep, double *PenaltyKeep,
                      double *GammasKeep, double *gFKeep,
                      int GFlag, 
                      int WeightLen, double *Weights, int *pFinaltt) {
	  //int IIPF = 2;
	  int IIPF = -1;
	  int RetFlag = 0;
	if (IIPF > 0) {  
		Rprintf((char*)"LLABI: Creating Lars Object, lambda = %.4f\n", lambda);
		R_FlushConsole();
		R_ProcessEvents(); 
    }                     
	LARSObject *MyLARSObject = new LARSObject(yys, NLen, xxs,  kLen, 
                                    InitBetas,  lambda);
    MyLARSObject->Ontt = 0;
    if (IIPF > 0) {
	    Rprintf((char*)"LLABI: About to create GAll\n");
	    R_FlushConsole();
    }  
    if (GFlag == 1) {
	    MyLARSObject->CreateGAll();
	    MyLARSObject->PrintFlag = 0;
	   
       if (IIPF > 0) { 
          MyLARSObject->GAllPrint();
          Rprintf((char*) "GGFlag = %d: Performing GGAll=1, WFlag = 0, lambda = %.5f\n",
	                                GFlag, (double) MyLARSObject->lambda);
       }
    } else if (GFlag == 0) {
	    MyLARSObject->PrintFlag = 0;
	    if (IIPF > 0) {
	    Rprintf((char*) "GGFlag = %d: Performing GGAll=0, WFlag = 0, lambda = %.5f\n",
	                                GFlag, (double) MyLARSObject->lambda);
        R_FlushConsole();	
        }                                
    }  else if (GFlag == 2) {
	    MyLARSObject->PrintFlag = 0;
	    MyLARSObject->CreateWAll(WeightLen, Weights);
	    if (IIPF > 0) {
		    Rprintf((char*) "GGFlag = %d: Performing GGAll=0, WFlag = 1, lambda = %.5f\n",
		                                GFlag, (double) MyLARSObject->lambda);	
	        R_FlushConsole();	  
        }                                  
    } else if (GFlag == 3) {
	    MyLARSObject->CreateGAll();
	    MyLARSObject->PrintFlag = 0;
       if (IIPF > 0) { 
	        MyLARSObject->GAllPrint();   
		      Rprintf((char*) "GGFlag = %d: Performing GGAll=1, WFlag = 1, lambda = %.5f\n",
		                                GFlag, (double) MyLARSObject->lambda);	  
		      R_FlushConsole(); 
        }  
       MyLARSObject->CreateWAll(WeightLen, Weights);
    }
    R_FlushConsole();
    R_ProcessEvents();
    if (IIPF > 0) {
	    Rprintf((char*)"LLABI: About to run LLA Algorithm\n");
		  R_FlushConsole();
	   	R_ProcessEvents();        
    }
 
    //MyLARSObject->CreateWAll(WeightLen, Weights);
    RetFlag = MyLARSObject->LarsLassoAlgorithm2();   
    if (IIPF > 0) {   
	    Rprintf((char*)"LLABI: Finished LLA, Copying Data, RetFlag = %d\n", RetFlag);
		R_FlushConsole();
		R_ProcessEvents();      
    }                               
	int ii,cn1 = MyLARSObject->kLen * (MyLARSObject->kLen+2);
	pFinaltt[0] = MyLARSObject->Ontt;
	  for (ii = 0; ii < MyLARSObject->kLen; ii++) {
		  OutBetas[ii]=(double)MyLARSObject->OutBetas[ii];                     
		  InitBetas[ii]=(double) MyLARSObject->InitBetas[ii];  
    }
    for (ii = 0; ii < MyLARSObject->Ontt; ii++) {
		  JKeep[ii] = (int) MyLARSObject->JKeep[ii];
		  GammasKeep[ii] = (double) MyLARSObject->GammasKeep[ii];
		  gFKeep[ii] = (double) MyLARSObject->gFKeep[ii];	
    }
    cn1 = (MyLARSObject->Ontt+2) * 4;
    for (ii = 0; ii < cn1; ii++) {
	      PenaltyKeep[ii] = (double) MyLARSObject->PenaltyKeep[ii];
    }
    cn1 = (MyLARSObject->Ontt+2) * (MyLARSObject->kLen);
    for (ii = 0; ii < cn1;ii++) {	  		  
		  BetasKeep[ii]=(double) MyLARSObject->BetasKeep[ii];   		
    }
    if (IIPF > 0) {
         Rprintf((char*)"LLABI: Freeing MyLARSObject\n");
		 R_FlushConsole();
		 R_ProcessEvents();  
    }         
    delete(MyLARSObject);
    if (IIPF > 0) {
         Rprintf((char*)"LLABI: Freeing MyLARSObject deleted\n");
		 R_FlushConsole();
		 R_ProcessEvents();  
    }       
    return(1);  	                      
}


/////////////////////////////////////////////////////////////////
// int LARSObject::LarsLassoAlgorithm2()
//
//   Runs the algorithm over the LARSObject.  
//
//  Steps are described in the algorithm.  
//
int LARSObject::LarsLassoAlgorithm2( 
                     ) {
	int AssessCCFlag = 0;
	/*if (this->PrintFlag > -1) {
      Rprintf( (char *) "Entering LarsLassoAlgorithm: kLen = %d, lambda=%.6f \n", 
               (int) this->kLen, (double) lambda);
      Rprintf( (char *) "Processing mode GFlag = %d, WDiagFlag = %d\n,", 
         (int) GFlag, (int) WDiagFlag);
      R_FlushConsole();	 
    }   */   
    if (ISNAN(this->lambda)) {
	    Rprintf((char*)"Ourlambda is currently NAN, flaw\n");
	    R_FlushConsole();
	    return(-1);
    } 
    /*
    if (this->PrintFlag > 8) {
	    Rprintf((char*) "EnteringLarsLassoAlgorithm, here's the YY\n");
	    int dii, djj;
	    for (dii = 0; dii < NLen; dii++) {
		    Rprintf((char*) "%d=%.4f,", dii, (double) this->yys[dii]);
	    }
	    Rprintf((char*) "\n And here is the XX\n");
	    for (dii = 0; dii < NLen; dii++) {
		    Rprintf((char*)" %d = [ ", dii);
		    for (djj=0; djj < kLen; djj++) {
			    Rprintf((char*)" %.3f, ", (double) xxs[djj * NLen + dii]);
		    }
		    Rprintf((char*)"\n");
	    }
	    R_FlushConsole();
    }*/
    int ii,jj,tt;
    //int kk, OnFun;   
    MaxIters = (int) floor(this->kLen * DEFCONST);
      if (this->kLen > this->NLen-1) {
	        MaxIters =(int) floor(this->NLen * DEFCONST);
      }
    int InvertFlag;
      
    for (tt = 0; tt < this->MaxIters; tt++) {
	   this->Ontt = tt+1;
	    //Rprintf((char*)"On tt = %d / %d \n", tt, this->MaxIters);
	    //R_FlushConsole();
	   /* if (tt == 8 || tt == 9 || tt == 10 || tt == 11 ) {
		    PrintFlag = 5;
	    } else {
		    PrintFlag = 0;
	    }*/
	   //  if ( tt == 7 ) {
		//    PrintFlag = 5;
	   // } else {
		//    PrintFlag = 0;
	   // }	    
	  /*if (this->PrintFlag > 2) {
         Rprintf( (char *) "LLA: tt = %d of %d \n", (int) tt, (int) this->MaxIters);
         R_FlushConsole();	 
      }  
	  if (this->PrintFlag > 2) {
         Rprintf( (char *) "LLA: tt = %d; Running AsessCC \n", (int) tt, (int) this->MaxIters);
         R_FlushConsole();	 
      } */       
	    AssessCCFlag = this->AssessCC();  
	    if (AssessCCFlag < 0) {
		    Rprintf((char*)"AssessCCFlag = %d, Flaw Go Home\n", AssessCCFlag);
		    R_FlushConsole();
		    return(-1);
	    }
	    if (this->PrintFlag > 2) {
         Rprintf( (char *) "LLA: tt = %d; OnInL = %d, OnL = %d\n", (int) tt, (int) this->OnInL, this->OnL);
         R_FlushConsole();	 
      } 
      if ((double) this->OnL > (double) this->PrevOnL + 1.5) {
	         Rprintf( (char*) "OnL > PrevOnL + 1, tt = %d, OnL = %d, PrevOnL = %d, kLen = %d, error \n",
	             tt, this->OnL, this->PrevOnL, this->kLen);
	         R_FlushConsole();
	         R_ProcessEvents();
	         this->ElseCC(this->PrevOnL, this->OnL);
      }
      if ( (double) this->OnL < (double) this->PrevOnL - 1.5) {
	         Rprintf( (char*) "OnL < PrevOnL - 1, tt = %d, OnL = %d, PrevOnL = %d, kLen = %d, error \n",
	             tt, this->OnL, this->PrevOnL, this->kLen);
	         R_FlushConsole();
	         R_ProcessEvents();
	         this->ElseCC(this->PrevOnL, this->OnL);
      }
      this->PrevOnL = this->OnL;
      if (this->OnL > this->kLen || this->OnInL > this->kLen  || this->OnL == 0) {
	        Rprintf( (char *) "OnL = %d, OnInL = %d, kLen = %d, error \n", 
	            this->OnL, this->OnInL, this->kLen);
	        R_FlushConsole();
	        R_ProcessEvents();
	        break;
      }
      /*if (this->PrintFlag > 2) {
         Rprintf( (char *) "LLA: tt = %d: Calculating cc\n", (int) tt);
         R_FlushConsole();	 
      } */ 
      this->CMaxNow = (double) fabs((double) this->CurMaxs[OnL-1]);  
      if (GAll == NULL) { 
            SetMeXIn (this->kLen, this->NLen, this->OnL, this->ActiveOns, 
              this->ActiveSigns, this->xxs , this->XXA);
      }
      this->MakeGGA();

      SetMeI ( this->OnL, this->IGGA );
      /* if (this->PrintFlag > 2) {
         Rprintf( (char *) "LLA: tt = %d: Calculating IGGA\n", (int) tt);
         R_FlushConsole();	 
         }
      */  
        InvertFlag = SolveSymPosMatrix(this->OnL, this->GGAA, this->IGGA);
        if (PrintFlag >= 2) {
	        Rprintf("LarsObject: Post Inversion\n");
	        this->GGAPrint();
        }        
        if (InvertFlag < 0) {
	          Rprintf((char*)"LLA: Invert Failed on tt = %d, OnL = %d\n", tt, this->OnL);
	          Rprintf((char*)"LLA: Suggests flaw in GGAA\n");
	          return(-1);
        }   
        this->CalcAAA();       
        this->CreateWa();
        if (this->PrintFlag > 3) {
	         Rprintf((char*)"LLA: AAA = %.4E \n", (double) this->AAA);
	         R_FlushConsole();
        }
      //if (this->PrintFlag > 2) {
      //   Rprintf( (char *) "LLA: tt = %d: Calculating Matrix Multiplying, ua/aa\n", (int) tt);
      //   R_FlushConsole();	 
      //}     
        //if (this->PrintFlag > 3 && this->GFlag == 1) {
	    //    GAllPrint();
        //} 
        Makeaa();   
        if (this->PrintFlag > 3) {
	        Rprintf((char*) "Printing wa for tt = %d\n", tt);
	          for (ii = 0; ii < this->OnL; ii++) {
		         Rprintf((char*)"ii=%d, wa[%d] = %.5E, ActiveOns[ii] = %d, ActiveSigns[ii] = %d\n", 
		               ii, ii, (double) wa[ii], ActiveOns[ii], ActiveSigns[ii]);  
	          }
	          R_FlushConsole();
	        Rprintf((char*) "\n wa = c(");
	        for (ii = 0; ii < this->OnL-1; ii++) {Rprintf((char*) " %.8f, ",(double)  wa[ii]); }
	        Rprintf((char*) "%.8f)\n", (double) wa[this->OnL-1]);
	        Rprintf((char*) "Printing out AAA = %.4f\n", (double) AAA);
	        if (PrintFlag > 5 && this->kLen < 10) {
		        Rprintf((char*) "Printing aa for tt = %d\n", tt);
		        for (ii = 0; ii < this->kLen; ii++) {
			        Rprintf((char*)"aa[%d] = %.4f\n", ii, (double) aa[ii]);
			        R_FlushConsole();
		        }
	            Rprintf((char*) "\n aa = c(");
		        for (ii = 0; ii < this->kLen-1; ii++) {Rprintf((char*) " %.8f, ", (double) aa[ii]); }
		        Rprintf((char*) "%.8f)\n", (double) aa[this->kLen-1]);	
            }        
        }
        //if (this->PrintFlag > 2) {
        // Rprintf( (char *) "LLA: tt = %d: Calculating gammaLists\n", (int) tt);
        // R_FlushConsole();	 
        //}  
       this->LARSgammaLists();  

       if (GFlag == 1 || GFlag == 3) {
      //if (this->PrintFlag >2) {
      //   Rprintf( (char *) "GFlag = %d, LLA: tt = %d: Calculating YTXwScores\n", 
      //           (int)GFlag, (int) tt);
      //   R_FlushConsole();	 
      //} 	       
          CalculateYTXwScores(); 
       }
      //if (this->PrintFlag > 2) {
      //   Rprintf( (char *) "LLA: tt = %d: Calculating GFTAll\n", (int) tt);
      //   R_FlushConsole();	 
      //}     
    this->CalculateGFTAll(this->RunningMin); 
    if (GFlag == 1) {
       CurSumwaXATXAwa = this->GFTAll * this->GFTAll;
    }
    if(this->PrintFlag > 3 ) {                     
             Rprintf( (char *) "LarsLasso reached OnL = %d, GFTAll = %f, RunningMin = %f, for tt = %d \n", 
                           this->OnL,  (double) this->GFTAll, (double) this->RunningMin, tt);
               R_FlushConsole();  
    } 
              
         this->GammasKeep[tt] = this->GFTAll;
         this->gFKeep[tt] = this->RunningMin; 
         //if (this->PrintFlag > 1) {
         //Rprintf( (char *) "LLA: tt = %d: CWhat to do with GRT All?\n", (int) tt);
         //R_FlushConsole();	 
         //}                   
    if ((double) GFTAll <= (double) 0.0) {
	        //if (PrintFlag > 0) {
	         Rprintf( (char *) "We've Reached a Negative value for GFTAll = %f \n", (double) GFTAll);
          //}
	       //tt = MaxIters-1;
	       if (tt >= MaxIters) {
		       Rprintf((char*)"GFTAll Negative: Too Bad, tt is huge\n"); R_FlushConsole();
	       }
	       R_ProcessEvents();
	       R_FlushConsole();
	       this->FFeatConst= 0.0;
           this->FeedFreeBetas(tt, 0.0);	
	       this->PrevPenalty = (double) this->PenaltyKeep[4 * tt + 3];
	       this->CurrentPenalty = (double) this->PenaltyKeep[4 * (tt+1) + 3];
	        for (ii = 0; ii < kLen; ii++) { this->PrevBetas[ii] = this->OnBetas[ii];}       
	       break;
    }  else if ((double) this->GFTAll < (double) this->RunningMin) {  
      // When the GFTAll exceeds a value, certain issues need to be recalculated
	       CurrentGoat = -20;
	       if (PrintFlag > 1) {
	          Rprintf( (char *) "GFTAll = %f less than Running Min = %f \n", (double) GFTAll, (double)RunningMin);
	          R_FlushConsole();
         }
	       if (tt >= MaxIters) {
		       Rprintf((char *)"GFTAll Negative: Too Bad, tt is huge\n"); R_FlushConsole();
	       }           
	       R_FlushConsole();
	       R_ProcessEvents();
	       this->FFeatConst = this->GFTAll;
	       FeedFreeBetas( tt, (double) this->GFTAll);
	           this->PrevPenalty = (double) this->PenaltyKeep[4 * tt + 3];
	           this->CurrentPenalty = (double) this->PenaltyKeep[4 * (tt+1) + 3];	           
           R_ProcessEvents();
           R_FlushConsole();	        
	       break;
    }  else {
	    this->FFeatConst = (double) this->RunningMin;
	    FeedFreeBetas( tt,  (double) this->RunningMin );
	           this->PrevPenalty = (double) this->PenaltyKeep[4 * tt + 3];
	           this->CurrentPenalty = (double) this->PenaltyKeep[4 * (tt+1) + 3];	           
    } 
      //if (this->PrintFlag > 2) {
      //   Rprintf( (char *) "LLA: tt = %d: Making Penalty Judgement\n", (int) tt);
      //   R_FlushConsole();	 
      //}       
     if ((double) this->PrevPenalty > 0 && (double) this->PrevPenalty < (double) this->CurrentPenalty)  {
	       if (tt == 0) {
		         if (PrintFlag > 1) {
		             Rprintf( (char *) "At tt = 0, Previous Penalty= %.4f ",
                    (double)(this->PrevPenalty));
                 Rprintf((char*) "greater than Zero and Less than Current ");
                 Rprintf((char*) "Penalty %.4f\n", 
		                   (double) this->CurrentPenalty);
	           }
		          // Means that Penalty is so large that nothing is sufficiently correlated
	        } else { 
	             //Rprintf( (char *) "Previous Penalty greater than Zero and less than Current Penalty \n");
		         if (PrintFlag > 1) {
	              Rprintf( (char *) "Previous Penalty= %.4f greater than",
                     (double) this->PrevPenalty);
                Rprintf((char*) "Zero and Less than Current Penalty %.4f\n", 
		                 (double) this->CurrentPenalty);	             
		            Rprintf( (char *) "Previous Penalty parts 0-%.4f, 1- %.4f, ",
		               (double) this->PenaltyKeep[(tt) *4 +0],
		               (double) this->PenaltyKeep[(tt) *4 +1]);
		            Rprintf( (char*) "2- %.4f, 3-%.4f \n",
		               (double) this->PenaltyKeep[(tt) *4 +2],
                   (double) this->PenaltyKeep[(tt) *4 +3]);
		            Rprintf( (char *) "Current Penalty parts 0-%.4f, 1- %.4f, 2- %.4f, 3-%.4f \n",
		               (double) this->PenaltyKeep[(tt+1) *4 +0], (double) this->PenaltyKeep[(tt+1) *4 +1],
		               (double) this->PenaltyKeep[(tt+1) *4 +2], (double) this->PenaltyKeep[(tt+1) *4 +3]);
	           }
	           if (kLen-1 >= MaxIter) {
		           Rprintf("We were going to update an Unalterable Penalty, sorry\n");
		           Rprintf("kLen = %d, MaxIter = %d\n", kLen, MaxIter);
		           R_FlushConsole(); R_FlushConsole(); 
	           } else {    
			       this->PenaltyKeep[(kLen+1) * 4 + 0] = this->PenaltyKeep[(tt) * 4 + 3];
			       this->PenaltyKeep[(kLen+1) * 4 + 1] = this->PrevPenalty;
			       this->PenaltyKeep[(kLen+1) * 4 + 2] = this->CurrentPenalty;
			       this->PenaltyKeep[(kLen+1) * 4 + 3] = - tt;		               
 		       }               
	             R_FlushConsole();
	             R_ProcessEvents();
           }
           Rprintf("Quitting Becaomse PrevPenalty Exceeds Current Penalty");
           Rprintf("Prev = %.4f, Current = %.4f\n",
               this->PrevPenalty, this->CurrentPenalty);
	         break;
       }
     /*   if (this->PrintFlag > 2) {
           Rprintf( (char *) "LarsLasso Finishing Step tt = %d \n", tt);
           R_FlushConsole();
        }
      if (this->PrintFlag > 10) {
         Rprintf( (char *) "LLA: Printing OnBetas:  \n");
         for (int nnii = 0; nnii < this->kLen; nnii++) { 
	         Rprintf((char*) " %d=%.4f", nnii, (double) this->OnBetas[nnii]);
         }
         Rprintf((char*)"\n");         
         R_FlushConsole();	 
      }*/  
      //Rprintf(" Heck, tt = %d, and Ontt = %d\n", tt, Ontt);
    }
    if (this->PrintFlag > 1) {
         Rprintf( (char *) "LLA: LLA all Done, considering stopping time\n");
         R_FlushConsole();	 
    }    
    if (this->PrevPenalty > 0 && this->PrevPenalty < this->CurrentPenalty) {	   
	     Rprintf( (char *) "tt = %d, Previous Penalty parts 0-%.4f, 1- %.4f, 2- %.4f, 3-%.4f \n",
		      tt, (double) this->PenaltyKeep[(tt) *4 +0], (double) this->PenaltyKeep[(tt) *4 +1],
		      (double) this->PenaltyKeep[(tt) *4 +2],  (double) this->PenaltyKeep[(tt) *4 +3]);
		 //this->PenaltyKeep[(this->kLen+1) * 4 + 0] = this->PenaltyKeep[(tt) * 4 + 3];
		 //this->PenaltyKeep[(this->kLen+1) * 4 + 1] = this->PrevPenalty;
		 //this->PenaltyKeep[(this->kLen+1) * 4 + 2] = this->CurrentPenalty;
		 //this->PenaltyKeep[(this->kLen+1) * 4 + 3] = - (double) tt;
		 this->InitBetas[0] = -9999;
	//	 Rprintf( (char *) "tt = %d, Current Penalty parts 0-%.4f, 1- %.4f, 2- %.4f, 3-%.4f \n",
	//	    tt, (double) this->PenaltyKeep[(tt+1) *4 +0], (double) this->PenaltyKeep[(tt+1) *4 +1],
	//	               (double) this->PenaltyKeep[(tt+1) *4 +2], 
	//	               (double) this->PenaltyKeep[(tt+1) *4 +3]); 
         for (jj = 0; jj < this->kLen; jj++) {  
	            this->OutBetas[jj] = (double) this->PrevBetas[jj];     
	     }
    //     Rprintf( (char *) " Freeing with the PrevPenalty \n");
         R_FlushConsole();
         //Free(OnBetas);
         //Free(PrevBetas);
    //   if (this->PrintFlag > 1) {
        // Rprintf( (char *) "LLA: tt = %d=%d: End Case 2, Prev Penalty\n", 
        //          (int) tt, (int) this->MaxIters);
        // R_FlushConsole();	 
    //  }    
         return(2);
    }    else {
	  // if (this->PrintFlag > 1) {
       //  Rprintf( (char *) "LLA: tt = %d/%d: End Case 1, Current Penalty", 
       //     (int) tt, (int) this->MaxIters);
       //  Rprintf(" better than Old Penalty At Least\n");
       //  R_FlushConsole();	 
     //  }    
	     for (jj = 0; jj < this->kLen;jj++) {
           this->OutBetas[jj] = (double) this->OnBetas[jj];  
       }
	     //Rprintf( (char *) " Freeing with the Current Penalty \n");
         R_FlushConsole();
	     //Free(OnBetas);
	     //Free(PrevBetas);
	     return(1);
    }                     
}

/////////////////////////////////////////////////////////////////////////////////
//  CalculateGFTAll
//
//   This calculates the stopping condition for an iteration of the algorithm.
//   It may be that the minimizer of the "gammas", which is the stopping condition before
//   Adding another variable to the dataset,is the true minimizer.  This looks for the first
//   Active variable to turn "off" (return to zero), or the minimizer that stops at Onlambda.
//
//   Note that this comes from the derivative:
//       D/ dgamma \{    .5 sum( (Y-X*beta_0 - gamma * u_a )^2 ) + 
//                                             lambda * sum_j |beta_0_j+ gamma *w_a_j|  \}
//     =    2 u_a %*% ( Y - X* beta - gamma * ua ) + 
//                  lambda * sum_j |w_a_j| * sign(w_a_j * beta_0_j) = 0;
//  Hence one is solving for the correlation of u_a with the residuals as well as the
//    Sum of w_a as it interacts with current beta_0.  
int LARSObject::CalculateGFTAll(double OverallMax) {
	   if(PrintFlag > 3) {
	     Rprintf( (char *) " Calculating GFTAll, NLen = %d, OnL = %d, MaxIters = %d \n", 
	                          NLen, OnL, MaxIters);
	     R_FlushConsole();
         //Rprintf( (char *) "Come On Print Something Anything \n");
	     //   R_FlushConsole();	   
       }
	   int jj;
	   int RunningMinjj1 = -2, RunningMinjj2 = -2;
	   int OnSign = -999;
	   double MaxMove1 = -1, MaxMove2 = -1, ProMove;
       this->TFLV2 = (double) 0.0;
       this->GFT0 = (double) 0.0;
       this->GFT1 = (double) 0.0;
       double LDTFLV2 = 0, LDGFT0 = 0, LDGFT1 = 0;
       if (LDTFLV2 <= -1 || LDGFT0 == -1 || LDGFT1 == -1) {
         Rprintf("**  Calculate LDTFLV2 \n"); R_FlushConsole();
       }
       int goflag = 0, OnGo = 0;
       int MAXGO = 100;
       double MinGammas = 0.0;
       double RetAns = 0.0;
       L2Involved = 0;
       if (PrintFlag > 10) {
         Rprintf( (char *) "Going Through OnBetas\n");
         R_FlushConsole();
         R_ProcessEvents();
         MaxIters = this->kLen; if (MaxIters > 10) { MaxIters = 10; }
         for (jj = 0; jj < MaxIters; jj++) {
	          Rprintf( (char *) "    OnBetas[%d] = %.4f\n", jj, (double) OnBetas[jj]);
         }
         Rprintf( (char *) "Going Through BetaOlds\n");
         R_FlushConsole();
         R_ProcessEvents();
         for (jj = 0; jj < MaxIters; jj++) {
	          Rprintf( (char *) "    BetaOlds[%d] = %.4f\n", jj, (double) BetaOlds[jj]);
         }
       }
       while (goflag == 0 && OnGo < MAXGO) {
	       OnGo++;
	       if (PrintFlag > 3) {
	         Rprintf( (char *) "Updated GFTSearch: Attempting OnGo == %d, RetAns = %.4f \n", OnGo, (double) RetAns);
	         R_FlushConsole();
	         R_ProcessEvents();
           }
	       this->TFLV2 = (double) 0.0; LDTFLV2 = 0.0;
	       this->GFT0 = (double) 0.0;  LDGFT0 = 0.0;
	       this->GFT1 = (double) 0.0;  LDGFT1 = 0.0;  
       if (PrintFlag > 2) {
	       Rprintf((char*) "Updated GFTSearch:  goflag = %d, OnGo = %d, MAXGO = %d\n", goflag, OnGo, MAXGO);
	       R_FlushConsole();
       }        
       for (jj = 0; jj < this->OnL; jj++ ) { 
            if ( this->OnBetas[this->ActiveOns[jj]]== 0) {
		        if(PrintFlag > 3) {
		          Rprintf( (char *) "Beta[%d] = 0\n", ActiveOns[jj]);
	              R_FlushConsole();
	              R_ProcessEvents();
                }
               if (this->wa[jj] < 0 ) { this->GFT0 -= this->wa[jj]; LDGFT0 -= this->wa[jj];
               }  else { this->GFT0 += this->wa[jj]; LDGFT0+= this->wa[jj]; }
	        } else {
		       OnSign = 0;
		       if (this->OnBetas[this->ActiveOns[jj]] > (double) 0.0 && 
		           this->ActiveSigns[jj] > 0) {
			         OnSign = 1;
	           } else if (this->OnBetas[this->ActiveOns[jj]] > (double) 0.0 &&
	                       this->ActiveSigns[jj] < 0 ) {
	                 OnSign = -1;
               } else if (this->OnBetas[this->ActiveOns[jj]] < (double) 0.0 && 
		           this->ActiveSigns[jj] > 0) {
			         OnSign = -1;
	           } else if (this->OnBetas[this->ActiveOns[jj]] < (double) 0.0 &&
	                       this->ActiveSigns[jj] < 1 ) {
	                 OnSign = 1;
               } 
		       if (OnSign == -1) {
			       ProMove = this->OnBetas[this->ActiveOns[jj]] / this->wa[jj];
			       if (ProMove < 0) { ProMove = -1.0 * ProMove; }
			       if (RunningMinjj1 < 0) {
				       MaxMove1 = ProMove;
				       RunningMinjj1 = jj;
			       } else if (MaxMove1  >  ProMove ) {
					   MaxMove1  =  ProMove;
				       RunningMinjj1 = jj;
			       }
                   if ( ProMove > MinGammas) {
	                      if (this->wa[jj] < 0) { 
		                       this->GFT0 += this->wa[jj]; LDGFT0 += this->wa[jj];
                          } else { this->GFT0 -= this->wa[jj]; LDGFT0 -= this->wa[jj]; }

		             } else {
		                if (this->wa[jj] < 0) { 
			               this->GFT0 -= this->wa[jj]; LDGFT0 -= this->wa[jj];
                        }  else { this->GFT0 += this->wa[jj]; LDGFT0 += this->wa[jj];}
	                 }
                } else  {
	                if (this->wa[jj] < 0) { 
		                this->GFT0 -= this->wa[jj]; LDGFT0 -= this->wa[jj]; 
                    }else { this->GFT0 += this->wa[jj]; LDGFT0 += this->wa[jj];}    
                }        		         
		    }	          
       }
         if (PrintFlag > 2) {
           Rprintf((char*)"UpdateGFTAll: LDGFT0 = %.5E \n", (double) LDGFT0);
         }
              this->GFT0 = - (double) this->lambda * (double) GFT0 / (double) 2.0;
              LDGFT0 = - this->lambda * LDGFT0 / ((double) 2.0);
       if (this->GFlag == 0 || this->GFlag == 2) {
           this->GFT1 = 0.0;
           LDGFT1 = 0.0;
	       for (jj = 0; jj <  NLen; jj++) {      
		        this->GFT1 += (double) this->ua[jj] * (double) this->OnResid[jj];   
		        LDGFT1 += this->ua[jj] * this->OnResid[jj]; 
	       }
       } else if (this->GFlag == 1 || this->GFlag == 3) {
	       this->GFT1 = CurSumYTXAwa - CurSumOnBetaXTXAwa;      
	       LDGFT1 = CurSumYTXAwa - CurSumOnBetaXTXAwa;      
       } else {
	       Rprintf((char*)"GFTAll: GFlag error = %d\n", GFlag);
	       R_FlushConsole();
	       return(-1);
       }
         if (PrintFlag > 2) {
           Rprintf("UpdateGFTAll: FixedLDGFT0 = %.5E, LDGFT1 %.5E \n", 
                     (double) LDGFT0, (double) LDGFT1);
         }
        RetAns = (double) this->GFT0 +  (double) this->GFT1;
        RetAns =  LDGFT0 + LDGFT1;
        this->GFT0 = (double) LDGFT0; this->GFT1 = (double) LDGFT1;
        this->GFTAll = RetAns;
        if ((double) OverallMax < (double) RetAns) {
	         goflag = 1; return( (int) RetAns);
        } else if ( (MaxMove1 > 0 && MaxMove1 < (double) RetAns) ) {
	             goflag = 0;
		         MinGammas = MaxMove1;
        } else {
	        goflag = 1; return((int) RetAns);    
        }
      }
      return((int) RetAns);   
       while (goflag == 0 && OnGo < MAXGO) {
	       OnGo++;
	       if (PrintFlag > 3) {
	         Rprintf( (char *) "Attempting OnGo == %d, RetAns = %.4f \n", OnGo, (double) RetAns);
	         R_FlushConsole();
	         R_ProcessEvents();
           }
       this->TFLV2 = (double) 0.0;
       this->GFT0 = (double) 0.0;
       this->GFT1 = (double) 0.0;	       
       for (jj = 0; jj < this->OnL; jj++ ) {
	       if (PrintFlag > 4) {
		       Rprintf( (char *) "Looking jj = %d, OnL = %d, ActiveOns[jj] = %d, TotK = %d \n", jj, OnL,
		                     this->ActiveOns[jj], this->kLen );
		       R_FlushConsole();
		       R_ProcessEvents();
		       Rprintf( (char *) "OnBetas[%d] = %.3f, BetaOlds[%d] = %.3f, wa[%d] = %.4f \n", ActiveOns[jj],
		              (double) OnBetas[ActiveOns[jj]], this->ActiveOns[jj], (double) BetaOlds[ActiveOns[jj]],
		               jj, (double) this->wa[jj]);
		       R_FlushConsole();
		       R_ProcessEvents();
          }
	       if (this->ActiveOns[jj] >= this->kLen) {
		       Rprintf( (char *) "Error ActiveOns[jj] = %d, OnL = %d, TotK = %d \n", 
		            this->ActiveOns[jj], this->OnL, this->kLen);
		       R_FlushConsole();
		       R_ProcessEvents();
		       return(-1);
	       }	
	       if (PrintFlag > 4) {       	 
	         Rprintf( (char *) "TFLV2[0] = %.4f, GFT0[0] = %.4f, GFT1[0] = %.4f\n", 
	                        this->TFLV2, this->GFT0, this->GFT1);
	         R_FlushConsole();
	         R_ProcessEvents(); 
           }     
           if (L2Involved == 1) {
	       if ( (double) this->OnBetas[this->ActiveOns[jj]] - 
	             (double) this->BetaOlds[this->ActiveOns[jj]] == 0) {
		        if (PrintFlag > 4) {
		           Rprintf( (char *) "L2Involved = 1, Beta[%d] - BetaOlds[%d] = 0\n", this->ActiveOns[jj], this->ActiveOns[jj]);
	              R_FlushConsole();
	              R_ProcessEvents();
               }
		       this->TFLV2 += fabs((double) this->wa[jj] );
	       } else {
		        //Rprintf( (char *) "Beta[%d] - BetaOlds[%d] = %f\n", ActiveOns[jj], ActiveOns[jj],
		        //(double) (  this->OnBetas[this->ActiveOns[jj]] - this->BetaOlds[this->ActiveOns[jj]] ));
	            //R_FlushConsole();
	            //R_ProcessEvents();
		       OnSign = GetSign( ((double) this->OnBetas[this->ActiveOns[jj]] - 
		                         (double) this->BetaOlds[this->ActiveOns[jj]]) * 
		                                     this->ActiveSigns[jj] );
		       if (OnSign == -1) {
			       if (RunningMinjj2 < 0) {
				       MaxMove2 = fabs((double) this->OnBetas[this->ActiveOns[jj]] - 
				                       (double) this->BetaOlds[this->ActiveOns[jj]]) / fabs(this->wa[jj]);
				       RunningMinjj2 = jj;
			       } else if (MaxMove2  >  fabs((double) this->OnBetas[this->ActiveOns[jj]] - 
				                       (double) this->BetaOlds[this->ActiveOns[jj]] ) / fabs(this->wa[jj]) ) {
					   MaxMove2  =  fabs((double) this->OnBetas[this->ActiveOns[jj]] - 
				                       (double) this->BetaOlds[this->ActiveOns[jj]]) / fabs(this->wa[jj]);
				       RunningMinjj2 = jj;
			       }
			       if (fabs( (double) this->OnBetas[this->ActiveOns[jj]] - 
		                         (double) this->BetaOlds[this->ActiveOns[jj]] ) / 
		             fabs( (double) this->wa[jj] ) > MinGammas) {
			              this->TFLV2 += fabs((double) this->wa[jj]) * 
			                 GetSign( (double) this->OnBetas[this->ActiveOns[jj]] - 
			                 (double) this->BetaOlds[this->ActiveOns[jj]] );
		             } else {
		                this->TFLV2 += fabs((double) this->wa[jj]);
	                 }
		        } else {
		            this->TFLV2 += fabs((double) this->wa[jj]);
		        }         		        
		   }
	       } 
	       if ( (double) this->OnBetas[this->ActiveOns[jj]]== 0) {
		        if(PrintFlag > 4) {
		          Rprintf( (char *) "Beta[%d] = 0\n", ActiveOns[jj]);
	              R_FlushConsole();
	              R_ProcessEvents();
                }
		       this->GFT0 += fabs((double) this->wa[jj] );
	       } else {
		       OnSign = GetSign( (double) this->OnBetas[this->ActiveOns[jj]]  * 
		                    (double) this->ActiveSigns[jj]);
		       if (OnSign == -1) {
			       if (RunningMinjj1 < 0) {
				       MaxMove1 = fabs((double) this->OnBetas[this->ActiveOns[jj]]) / fabs(this->wa[jj]);
				       RunningMinjj1 = jj;
			       } else if (MaxMove1  >  fabs((double) this->OnBetas[this->ActiveOns[jj]]) /
			                    fabs((double)this->wa[jj]) ) {
					   MaxMove1  =  fabs((double) this->OnBetas[this->ActiveOns[jj]])/
					      fabs((double)this->wa[jj]);
				       RunningMinjj1 = jj;
			       }
                   if (fabs( (double) this->OnBetas[this->ActiveOns[jj]])  / 
		             fabs( (double) this->wa[jj] ) > MinGammas) {
			              this->GFT0 += fabs((double) this->wa[jj]) * 
			                 OnSign;
		             } else {
		                this->GFT0 += fabs((double) this->wa[jj]);
	                 }
	             if (PrintFlag > 4) {
			        Rprintf((char*)"jj = %d, NegSign, OnSign = %d, GFT0 = %.4f, wa[%d] = %.4f, lambda = %.4f\n", 
			                   jj, OnSign, (double) this->GFT0, jj, (double) this->wa[jj], (double) lambda); 
		        }  
		        } else {
		            this->GFT0 += fabs((double) this->wa[jj]);
		        }
     		    /*if (PrintFlag > 4) {
			        Rprintf((char*)"jj = %d, Else/Else, OnSign = %d, GFT0 = %.4f, wa[%d] = %.4f, lambda = %.4f\n", 
			                   jj, OnSign, (double) this->GFT0, jj, (double) this->wa[jj], (double) lambda); 
		        } */ 	        		         
		    }
				/*if (PrintFlag > 4) {
			        Rprintf((char*)"jj = %d, Else, OnSign = %d, GFT0 = %.4f, wa[%d] = %.4f, lambda = %.4f\n", 
			                   jj, OnSign, (double) this->GFT0, jj, (double) this->wa[jj], (double) lambda); 
		        } */     
       }
		              
       this->GFT0 = - (double) this->lambda * (double) GFT0 / (double) 2.0;
       this->TFLV2 = -(double) this->TFLV2 * (double) 0.0 / (double) 2.0;
       if (this->GFlag == 0) {
	       this->GFT1 = 0.0;
	       if (FakeFlag >= 0) {
		       if (ua == NULL || OnResid == NULL) {
		       Rprintf((char*)"CalcGFTAll Error, FakeFlag > 0 but ua, OnResid are NULL\n");
		       R_FlushConsole(); R_ProcessEvents(); return(-1);
	           }
		       for (jj = 0; jj <  NLen; jj++) {      
			        this->GFT1 += (double) this->ua[jj] * (double) this->OnResid[jj];    
		       }
	       }
       } else if (this->GFlag == 1 || this->GFlag == 3) {
	       this->GFT1 = CurSumYTXAwa - CurSumOnBetaXTXAwa;            
       } else {
	       Rprintf((char*)"GFTAll: GFlag error = %d\n", GFlag);
	       R_FlushConsole();
	       return(-1);
       }
        RetAns = (double) this->GFT0 + (double) this->TFLV2 + (double) this->GFT1;
        this->GFTAll = RetAns;
        goflag = 1;
        if ((double) OverallMax < (double) RetAns) {
	         goflag = 1;
        } else if ( (MaxMove1 > 0 && MaxMove1 < (double) RetAns) ||
             (MaxMove2 > 0 && MaxMove2 < (double) RetAns) ) {
	             goflag = 0;
	             if ( (MaxMove1 > 0 ) && ((double) MaxMove1 < (double) MaxMove2) ) {
		               MinGammas = MaxMove1;
	             } else  {
		               MinGammas = MaxMove2;
	             }
        }
      }
         /*if (PrintFlag > 2) {
           Rprintf( (char *) " Returning GFTAll, NLen = %d, OnL = %d \n", NLen, OnL);
           Rprintf(  (char*) "  GFTAll = %.4f, GFT0 = %.4f, GFT1 = %.4f \n", 
                               (double) this->GFTAll, (double) this->GFT0, 
                               (double) this->GFT1);
	       R_FlushConsole();
         }*/
       return((int) RetAns);      
}
 	    

//////////////////////////////////////////////////////////////
//  LARSObject::CalcAAA()
//     
//    A is the inverse square root of the sum:
//                vec(1)^T *W_A^-1 G_A^-1* W_A^-1 * (vec1)
//     This adds up all entries in the Inverse X_A^T X_A matrix G_A
//     Weights are reintroduced as necessary.
double LARSObject::CalcAAA() {
	    //Rprintf( (char *) " Creating AAA, OnL = %d\n ", OnL);
	    //R_FlushConsole();
        this->AAA = 0.0;
        int OnFun = 0,jj, ii;
        double IDiagjj;
        if (WDiagFlag == 1) {
	        for (jj = 0; jj < this->OnL; jj++) {
		         IDiagjj = ((double) 1.0) / WDiag[ActiveOns[jj]];
		         for (ii = 0; ii < this->OnL; ii++) {
			         this->AAA += IDiagjj * this->IGGA[OnFun]  / ( WDiag[ActiveOns[ii]]);
			         OnFun++;
		         }
	        }
        } else {
	        for (jj = 0; jj < this->OnL; jj++) {
		         for (ii = 0; ii < this->OnL; ii++) {
			         this->AAA += this->IGGA[OnFun];
			         OnFun++;
		         }
	        }	        
        }
        //Rprintf( (char *) "AAA = %f \n", (double) 1  / sqrt(AAA));
        this->AAA = 1 / sqrt((double) AAA);
        return( this->AAA );
}

///////////////////////////////////////////////////////////////////////////
//  LARSInit()
//
//    If one chooses to run the LARS algorithm again with same LARSObject,
//  One needs to erase certain values from the system, this operation preps
//  for such a need, returning beta values to zero, and thus returning residuals
//  and such to their maximum values.
//
//
  int LARSObject::LARSInit () {
	      if (PrintFlag > 2) {
		      Rprintf((char*)"LARSObject:: Creating NLen = %d, kLen = %d\n", NLen, kLen);
		      R_FlushConsole();
		      R_ProcessEvents();
	      }
	      EpsilonBM = .00000002;
	      CurrentPenalty = -1.0;
	      PrevPenalty = -1.0;
	      Ontt = 0;
	     
	      int ii, tt; // int jj
	      SumCurResids = SumYYSq;
	      this->PrevOnL = 0;
          MaxIters = kLen;
           if(PrintFlag > 2) {
	            Rprintf((char*)"LARSObject:: Creating Memory for kLen Items, NLen = %d, kLen=%d\n", 
	                 NLen, kLen);
	            R_FlushConsole();
	            R_ProcessEvents();
            }
           
            if(PrintFlag > 2) {
	            Rprintf((char*)"LARSObject: Copying in Initial Betas\n");
	            R_FlushConsole();
	            R_ProcessEvents();
	            Rprintf((char*)"LARSObject: rInitBetas[0] = %.6f, rInitBetas[kLen-1] = %.6f\n",
	                                (double)  InitBetas[0], (double)InitBetas[kLen-1]);
	            R_FlushConsole();
	            R_ProcessEvents();	            
            }		
            if (InitBetas != NULL) {     
			    for (ii =0;ii < kLen; ii++) {
				      OnBetas[ii] = (double) InitBetas[ii];
			          PrevBetas[ii] = 0;
			          BetasKeep[ii] = InitBetas[ii];
				      //BetaOlds[ii] = (double) InitBetas[ii];
				}
            }  else {
			    for (ii =0;ii < kLen; ii++) {
				      OnBetas[ii] = 0;
			          PrevBetas[ii] = 0;
			          BetasKeep[ii] = 0;
				      //BetaOlds[ii] = (double) InitBetas[ii];
				}            
            }
            if(PrintFlag > 2) {
	            Rprintf((char*)"LARSObject:Finished Copying InitialBetas\n");
	            R_FlushConsole();
            }
         /*   if(PrintFlag > 2) {
	            Rprintf((char*)"Creating Memory for wa, aa, OnBetas, glists!\n");
	            R_FlushConsole();
	            R_ProcessEvents();
            }*/			
		   if(PrintFlag > 2) {
	            Rprintf((char*)"LARSObject: Setting OnActives to Zero\n");
	            R_FlushConsole();
            }       
			    for (ii = 0; ii < kLen; ii++) {
			        OnActives[ii] = 0;
			        wa[ii] = 0;
			        aa[ii] = 0;
			        GammasKeep[ii] = 0;
			        gFKeep[ii] = 0;
			    }
			    OnL  = 0;
           maxish = kLen;
           if (kLen < NLen) {maxish = NLen;}
           /* if(PrintFlag > 2) {
	            Rprintf((char*)"Creating Memory for ua, cc, PrevMuA\n");
	            R_FlushConsole();
	            R_ProcessEvents();
            }*/           

    //Rprintf( (char *) "LarsLasso: Data Loaded in \n");

		    /* if(FakeFlag >= 0 && PrintFlag > 2) {
	            Rprintf((char*)"OnMuA copying now\n");
	            R_FlushConsole();
	            R_ProcessEvents();
	            Rprintf((char*)"OnMuA:  yys[0] = %.5f, xxs[0] = %.6f, OnBetas[0] = %.6f, OnMuA[0] = %.6f\n",
	                (double) yys[0], (double) xxs[0],  (double) this->OnBetas[0], 
	                (double) OnMuA[0]);
	            R_FlushConsole();
	            R_ProcessEvents();
            } */
		    if (FakeFlag >= 0) {
			     for (ii = 0; ii < NLen; ii++) { OnMuA[ii] = 0; PrevMuA[ii] = 0;}
		    }
		    //int TotF = 0;
		    if (PrintFlag > 2) {
			    Rprintf("Either Loading OnMuA or loading GAllOnB \n");
			    R_FlushConsole(); R_ProcessEvents();
		    }
		    if (FakeFlag >=0 && (this->GFlag ==0 || this->GFlag == 2)) {
		        if (this->WDiagFlag == 0) {
		         MatTimesVec(NLen, kLen, xxs, OnBetas, OnMuA);
				   // for (jj = 0; jj < kLen; jj++) {
					 //   for (ii = 0; ii < NLen; ii++) { 
					 //	      OnMuA[ii] += (double) xxs[TotF] * (double) OnBetas[jj];
					 //	      TotF++;
					 //    }
				   //  
		          } else {
		          MatTimesVecVec(NLen, kLen, xxs, WDiag, OnBetas, WIV, OnMuA);
  			      //for (jj = 0; jj < kLen; jj++) {
  				    //for (ii = 0; ii < NLen; ii++) { 
  					  //    OnMuA[ii] += (double) xxs[TotF] * 
  					  //          (double) WDiag[jj]* (double) OnBetas[jj];
  					  //    TotF++;
  				    //}
			        // }     
	            }
            }
            if (FakeFlag < 0 && (this->GFlag == 1  || this->GFlag == 3)) {
	            if (this->WDiagFlag == 1) {
	              MatTimesVecVec(kLen, kLen, GAll,
                     WDiag, 
                     OnBetas, WIV, GAllOnB);  
		            //for (ii = 0; ii < kLen; ii++) {
			          //  GAllOnB[ii] = 0.0;
		            //}
		            //for (jj = 0; jj < kLen; jj++) {
			          //   TotF = jj * kLen;
			          //   for (ii = 0; ii < kLen; ii++) {
				        //     GAllOnB[jj] += GAll[TotF + ii] * WDiag[ii] * OnBetas[ii];				             
			          //   }
	              //  }		            
	            } else {
	              MatTimesVec(kLen, kLen,  (double*) GAll, 
                                         (double*) OnBetas, 
                                         (double*) GAllOnB);
		            //for (ii = 0; ii < kLen; ii++) {
			          //  GAllOnB[ii] = 0.0;
		            //}
		            //for (jj = 0; jj < kLen; jj++) {
			          //   TotF = jj * kLen;
			          //   for (ii = 0; ii < kLen; ii++) {
				        //     GAllOnB[jj] += GAll[TotF + ii] * OnBetas[ii];				             
			          //   }
	              //  }			            
		            
	            }
            }
            if (PrintFlag > 2) {
  			    Rprintf("Finished Either Loading OnMuA or loading GAllOnB \n");
  			    R_FlushConsole(); R_ProcessEvents();
		        }
            /*
		    if (PrintFlag > 2) {
			    Rprintf((char*)"Freshening up cc, ua, OnResid, PrevResid\n");
			    R_FlushConsole();
			    R_ProcessEvents();
	     	}
	     	if (PrintFlag > 3) {
		     	Rprintf( (char*) "Printing Initial OnMuA \n");
		     	for (ii = 0; ii < NLen; ii++) {
			     	Rprintf((char*) "      %d=%.6f\n", ii, (double) OnMuA[ii]);
		        }
		        R_FlushConsole();
	        }*/
		    for (ii = 0; ii < kLen; ii++) {  
			    cc[ii] = 0;
		    }
		    if (FakeFlag >= 0 && (this->GFlag == 0 || this->GFlag == 2)) {
			    SumCurResids = 0;
			    int One = 1;   double OneD = 1.0;       double ZeroD = 0.0;
			    F77_CALL(dcopy)(&NLen, yys, &One, OnResid, &One);
			    F77_CALL(daxpy)(&NLen, &OneD,
		           OnMuA, &One,
		           OnResid, &One);
		      F77_CALL(dcopy)(&NLen, OnResid, &One, PrevResid, &One);
		      F77_CALL(dscal)(&NLen, &ZeroD, ua, &One);
			    for (ii = 0; ii < NLen; ii++) {
			        SumCurResids += OnResid[ii]  * OnResid[ii];
			    } 
	        } else {
		        SumCurResids = SumYYSq;
	        }
		    /*if (PrintFlag > 3) {
		     	Rprintf( (char*) "Printing Initial OnResid -- PrevResid \n");
		     	for (ii = 0; ii < NLen; ii++) {
			     	Rprintf((char*) "      %d=%.6f=%.6f\n", 
			     	      ii, (double) OnResid[ii], (double) PrevResid[ii]);
		        }
		        R_FlushConsole();
	        }*/

            CurrentPenalty = -1.0;
            PrevPenalty = -1.0; 

         /*   if(PrintFlag > 2) {
	            Rprintf((char*)"Generating XXA, GGAA mem\n");
	            R_FlushConsole();
	            R_ProcessEvents();
            }
			  
              if(PrintFlag > 2) {
	            Rprintf((char*)"Calculating current penalty\n");
	            R_FlushConsole();
	            R_ProcessEvents();
            }*/
			    PrevPenalty = 0;
			    for (tt = 0; tt < (kLen +2) * 4; tt++) {
				       PenaltyKeep[tt] = 0;
			    }
			    if (FakeFlag >= 0 ) { 
			       PenaltyKeep[0] =  0;
				    for (tt = 0; tt < NLen; tt++) { 
					       PenaltyKeep[0] += (OnResid[tt]) * (OnResid[tt]);
					//       PrevPenalty += (OnResid[tt]) * (OnResid[tt]); 
					}
		        } else {
			        PenaltyKeep[0] = SumCurResids;
		        } 
			      for (tt = 0; tt <kLen; tt++) { 
				      PenaltyKeep[1] += fabs( OnBetas[tt] );
				      PenaltyKeep[2] += fabs( OnBetas[tt] - BetaOlds[tt] );
			      }
			      PenaltyKeep[1] = PenaltyKeep[1] * lambda;
			      PenaltyKeep[2] = PenaltyKeep[2] * 0.0;
			      PenaltyKeep[3] = PenaltyKeep[0] + PenaltyKeep[1] + PenaltyKeep[2];
			      PrevPenalty = PenaltyKeep[3];     
			    return(1);     
    }

////////////////////////////////////////////////////////////////////////
// int LARSObject::CalculateYTXwScores()
//
//    YTXw is particularly useful if we are doing the O(k*k) option
//      In this case we have a previously calculated X^T * Y vector
//    that gives one the ability to calculate a new residual correlation vector
//    to judge the next AssessCC() step.
//
int LARSObject::CalculateYTXwScores() {
	 CurSumYTXAwa = 0.0;
	 CurSumOnBetaXTXAwa = 0.0;
   CurSumwaXATXAwa = 0.0;
     
   //SumCurResids = SumYYSq;
   int jj, jj2, cnt1 = 0;
   if (cnt1 < 0) {
     Rprintf("CalculateYTXwScores, why?\n"); R_FlushConsole();
   }
	 if (WDiagFlag == 0) {
		 if (PrintFlag > 3) {
       Rprintf((char*)"CalculateYTXwScores, to calculate CurSumYTXAwa\n");
       R_FlushConsole();
		 }
		 for (jj = 0; jj < this->OnL; jj++) {
			 CurSumYTXAwa += XTYAll[ActiveOns[jj]] * ActiveSigns[jj] * 
			   this->wa[jj];                 
		 }
		 if (PrintFlag > 3) {
			 Rprintf((char*)"CalculateYTXwScores, to calculate CurSumOnBetastXXAwa\n");
			 R_FlushConsole();
		 }
		 for (jj2 = 0; jj2 < this->OnL; jj2++) {
			 CurSumOnBetaXTXAwa+= GAllOnB[ActiveOns[jj2]] *
		     ActiveSigns[jj2] * wa[jj2];
		 }
		 CurSumwaXATXAwa = this->GFTAll * this->GFTAll;
		 cnt1 = 0;
   } else {
	   if (PrintFlag > 3) {
			 Rprintf((char*)"CalculateYTXwScores, to calculate CurSumYTXAwa\n");
			 R_FlushConsole();
		 }
	 	 for (jj = 0; jj < this->OnL; jj++) {
	     CurSumYTXAwa += 
			   XTYAll[ActiveOns[jj]] * ActiveSigns[jj] * 
			   this->wa[jj] * this->WDiag[ActiveOns[jj]];                 
		 }
		 if (PrintFlag > 3) {
			 Rprintf((char*)"CalculateYTXwScores, to calculate CurSumOnBetastXXAwa\n");
       R_FlushConsole();
		 }
		      
		 for (jj2 = 0; jj2 < this->OnL; jj2++) {
	     CurSumOnBetaXTXAwa+= GAllOnB[ActiveOns[jj2]] *	
       ActiveSigns[jj2] * wa[jj2] *this->WDiag[ActiveOns[jj2]];
		 }
		 cnt1 = 0;
     CurSumwaXATXAwa = this->GFTAll * this->GFTAll;
	 }
   if (PrintFlag >3 ) {
	   Rprintf((char*) "CalculateYTXwScores: CumSumOnBetaXTXAwa = %.4f\n", 
	     (double) CurSumOnBetaXTXAwa);   
   }
   // SumCurResids = SumYYSq;
   return(1);	 
}

///////////////////////////////////////////////////////////////////////////////
//  LARSgammaLists()
//
//  This calculates the gammas for the current step (how far to move our target vector)
//    Before hitting (A-aa[jj])/(C-cc[jj]) or (A+aa[jj])/(C-cc[jj]) maximums.
double LARSObject::LARSgammaLists() {	
	    int jj;
	    int tjj;
	    this->RunningMin = -9;
	   //Rprintf( (char *) "Exploring Gamma Lists for kLen = %d, NLen = %d\n", (int) kLen, (int) NLen);
	   //R_FlushConsole();
	       if (this->kLen < 8) {
           //Rprintf( (char *) "   Writing the Lists \n");
           }
      // for (tjj = 0; tjj < this->OnInL; tjj++) {
	  //     jj = AAInActive[tjj];
	  for (tjj = 0; tjj < this->OnInL; tjj++) {
		  jj = AAInActive[tjj];
	       if ( (double) this->AAA - (double) this->aa[jj] < (double) this->EpsilonBM ) {
	             this->gammaList1[jj] = (double) -1.0;
	             if ((double) CMaxNow + (double) this->cc[jj] < (double) this->EpsilonBM) {
		             this->gammaList2[jj] = (double) -99.0;
	             } else { 
	                this->gammaList2[jj] = ((double) this->CMaxNow + (double) this->cc[jj]) 
	                               / ((double) this->AAA + (double) this->aa[jj]);
                 }
           } else if ((double) this->AAA + (double) this->aa[jj] < (double) this->EpsilonBM ) {
	             if ((double) this->CMaxNow - (double) cc[jj] < (double) this->EpsilonBM) {
		             this->gammaList1[jj] = (double) -99.0;
	             } else {
	                 this->gammaList1[jj] = ((double) this->CMaxNow - (double) this->cc[jj]) / 
	                                 ((double)this->AAA-(double) this->aa[jj]);
                 }
	              this->gammaList2[jj] = (double) -1.0;
           } else if ((double) this->CMaxNow - (double) this->cc[jj] < (double) this->EpsilonBM ) {
	              this->gammaList1[jj] = (double) -99.0;
	              this->gammaList2[jj] = (2.0 * (double) this->CMaxNow ) /
	                                     ((double)this->AAA + this->aa[jj]);
           } else if ((double) this->CMaxNow + (double) this->cc[jj] < (double) this->EpsilonBM) {
	              this->gammaList1[jj] = (2.0 * (double) this->CMaxNow) /
	                                     ((double) this->AAA - this->aa[jj]);
	              this->gammaList2[jj] = (double) -99.0;
           } else {
	              this->gammaList1[jj] = ((double) this->CMaxNow - (double) this->cc[jj]) / 
	                               ((double) this->AAA - (double) this->aa[jj]);
	              this->gammaList2[jj] = ((double) this->CMaxNow + (double) this->cc[jj]) / 
	                               ((double) this->AAA + (double) this->aa[jj]);
           }
           if (this->RunningMin <= (double) 0.0 ) {
	            if (this->gammaList1[jj] > 0.0  && this->gammaList2[jj] > 0.0 && 
	                this->gammaList1[jj] < this->gammaList2[jj]) {
		            this->RunningMin = (double) this->gammaList1[jj];
	            }  else if ((double) this->gammaList1[jj] > (double) 0.0 && 
	                        (double) this->gammaList2[jj] > (double) 0.0  &&
	                        (double) this->gammaList2[jj] < (double) this->gammaList1[jj]) {
		            this->RunningMin = (double) this->gammaList2[jj];
	            } else if (this->gammaList1[jj] > 0.0) {
		            this->RunningMin = this->gammaList1[jj];
	            } else if (this->gammaList2[jj] > 0.0) {
		            this->RunningMin = this->gammaList2[jj];
	            }
           } else {
	          if ((double) this->gammaList1[jj] > (double) 0.0 && 
	               (double) this->gammaList2[jj] > (double) 0.0 && 
	                  (double) this->gammaList1[jj] < (double) this->gammaList2[jj] && 
	                  (double) this->gammaList1[jj] < (double) this->RunningMin) {
		            this->RunningMin = (double) this->gammaList1[jj];
	           } else if ((double) this->gammaList1[jj] > (double) 0.0 && 
	                      (double) this->gammaList2[jj] > (double) 0.0 &&
	                  (double) this->gammaList2[jj] < (double) this->gammaList1[jj] &&
	                  (double) this->gammaList2[jj] < (double) this->RunningMin ) {
		            this->RunningMin = (double) this->gammaList2[jj];
	           } else if ((double) this->gammaList1[jj] > 0.0 && 
	                      (double) this->gammaList1[jj] < (double) this->RunningMin) {
		            this->RunningMin = (double) this->gammaList1[jj];
	           } else if ((double) this->gammaList2[jj] > (double) 0.0 && 
	                      (double) this->gammaList2[jj] < (double) this->RunningMin) {
		            this->RunningMin = (double) this->gammaList2[jj];
	           }
            }
            if ((PrintFlag > 3 && this->kLen < 10) || PrintFlag > 8) {
	           Rprintf((char*) "gammaList1[%d] = %.5f, gammaList2[%d] = %.5f, RunningMin = %.5f\n",
	              jj, (double) gammaList1[jj], jj, (double) gammaList2[jj], (double) RunningMin);
	           R_FlushConsole();	            
            }
           //if (kLen < 8) {
           //   Rprintf( (char *) " %d --   %.4f       %.4f    (curMin = %.4f) \n", 
	       //      jj, (double) gammaList1[jj], (double) gammaList2[jj], (double) RunningMin);
           //}
        }
        for (tjj = 0; tjj < this->OnL; tjj++) {
	       jj = ActiveOns[tjj];
	       if ( (double) this->AAA - (double) this->aa[jj] < (double) this->EpsilonBM ) {
	             this->gammaList1[jj] = (double) -1.0;
	             if ((double) CMaxNow + (double) this->cc[jj] < (double) this->EpsilonBM) {
		             this->gammaList2[jj] = (double) -99.0;
	             } else { 
	                this->gammaList2[jj] = ((double) this->CMaxNow + (double) this->cc[jj]) 
	                               / ((double) this->AAA + (double) this->aa[jj]);
                 }
           } else if ((double) this->AAA + (double) this->aa[jj] < (double) this->EpsilonBM ) {
	             if ((double) this->CMaxNow - (double) cc[jj] < (double) this->EpsilonBM) {
		             this->gammaList1[jj] = (double) -99.0;
	             } else {
	                 this->gammaList1[jj] = ((double) this->CMaxNow - (double) this->cc[jj]) / 
	                                 ((double)this->AAA-(double) this->aa[jj]);
                 }
	              this->gammaList2[jj] = (double) -1.0;
           } else if ((double) this->CMaxNow - (double) this->cc[jj] < (double) this->EpsilonBM ) {
	              this->gammaList1[jj] = (double) -99.0;
	              this->gammaList2[jj] = (2.0 * (double) this->CMaxNow ) /
	                                     ((double)this->AAA + this->aa[jj]);
           } else if ((double) this->CMaxNow + (double) this->cc[jj] < (double) this->EpsilonBM) {
	              this->gammaList1[jj] = (2.0 * (double) this->CMaxNow) /
	                                     ((double) this->AAA - this->aa[jj]);
	              this->gammaList2[jj] = (double) -99.0;
           } else {
	              this->gammaList1[jj] = ((double) this->CMaxNow - (double) this->cc[jj]) / 
	                               ((double) this->AAA - (double) this->aa[jj]);
	              this->gammaList2[jj] = ((double) this->CMaxNow + (double) this->cc[jj]) / 
	                               ((double) this->AAA + (double) this->aa[jj]);
           }
           if (this->RunningMin <= (double) 0.0 ) {
	            if( PrintFlag > 3) {
	               Rprintf((char*)"No Running Min, having to use jj = %d from active Set \n", jj);
                }
	            if (this->gammaList1[jj] > 0.0  && this->gammaList2[jj] > 0.0 && 
	                this->gammaList1[jj] < this->gammaList2[jj]) {
		            this->RunningMin = (double) this->gammaList1[jj];
	            }  else if ((double) this->gammaList1[jj] > (double) 0.0 && 
	                        (double) this->gammaList2[jj] > (double) 0.0  &&
	                        (double) this->gammaList2[jj] < (double) this->gammaList1[jj]) {
		            this->RunningMin = (double) this->gammaList2[jj];
	            } else if (this->gammaList1[jj] > 0.0) {
		            this->RunningMin = this->gammaList1[jj];
	            } else if (this->gammaList2[jj] > 0.0) {
		            this->RunningMin = this->gammaList2[jj];
	            }
           } else {
	          if ((double) this->gammaList1[jj] > (double) 0.0 && 
	               (double) this->gammaList2[jj] > (double) 0.0 && 
	                  (double) this->gammaList1[jj] < (double) this->gammaList2[jj] && 
	                  (double) this->gammaList1[jj] < (double) this->RunningMin) {
		        if( PrintFlag > 3) {
	               Rprintf((char*)"RunningMin, having from jj = %d from active Set \n", jj);
                }
		            this->RunningMin = (double) this->gammaList1[jj];
	           } else if ((double) this->gammaList1[jj] > (double) 0.0 && 
	                      (double) this->gammaList2[jj] > (double) 0.0 &&
	                  (double) this->gammaList2[jj] < (double) this->gammaList1[jj] &&
	                  (double) this->gammaList2[jj] < (double) this->RunningMin ) {
		        if( PrintFlag > 3) {
	               Rprintf((char*)"RunningMin, having from jj = %d from active Set \n", jj);
                }		                  
		            this->RunningMin = (double) this->gammaList2[jj];
	           } else if ((double) this->gammaList1[jj] > 0.0 && 
	                      (double) this->gammaList1[jj] < (double) this->RunningMin) {
		        if( PrintFlag > 3) {
	               Rprintf((char*)"RunningMin, having from jj = %d from active Set \n", jj);
                }		                      
		            this->RunningMin = (double) this->gammaList1[jj];
	           } else if ((double) this->gammaList2[jj] > (double) 0.0 && 
	                      (double) this->gammaList2[jj] < (double) this->RunningMin) {
		        if( PrintFlag > 3) {
	               Rprintf((char*)"RunningMin, having from jj = %d from active Set \n", jj);
                }		                      
		            this->RunningMin = (double) this->gammaList2[jj];
	           }
            }
           //if (kLen < 8) {
           //   Rprintf( (char *) " %d --   %.4f       %.4f    (curMin = %.4f) \n", 
	       //      jj, (double) gammaList1[jj], (double) gammaList2[jj], (double) RunningMin);
           //}
        }
        if (PrintFlag > 3) {
	    	  Rprintf((char*)"Current Running Min = %.5f \n", (double) this->RunningMin); 
	    	  R_FlushConsole();    
        }
        for (jj = 0; jj < this->OnL; jj++) {
	      if (GetSign(this->OnBetas[this->ActiveOns[jj]] * 
	                   this->wa[jj] * this->ActiveSigns[jj]) < 0) {
		      if (fabs(
		            (double)this->OnBetas[this->ActiveOns[jj]] 
		            / (double) this->wa[jj]) < 
		            this->RunningMin) {
			        if (PrintFlag > 3) {
				        Rprintf((char*)"Needed to take a Lasso RunningMin, Onjj = %d\n", ActiveOns[jj]);
			        }
			        this->RunningMin = (double) fabs((double)this->OnBetas[this->ActiveOns[jj]] / 
			                                (double)this->wa[jj]);
			        CurrentGoat = ActiveOns[jj];
		      }                 
	      }  
        } 
      
     if (PrintFlag > 3) {
	     Rprintf((char*)"gammalist1 = c(");
	    for (jj = 0; jj < kLen-1; jj++) {
		    Rprintf((char*) "  %.5f, ", (double) gammaList1[jj]);
	    }
	    Rprintf((char*) " %.5f)\n", (double) gammaList1[kLen-1] );
	    Rprintf((char*)"gammalist2 = c(");
	    for (jj = 0; jj < kLen-1; jj++) {
		    Rprintf((char*) "  %.5f, ", (double) gammaList2[jj]);
	    }
	    Rprintf((char*) " %.5f)\n", (double) gammaList2[kLen-1] );	   
	    
	    Rprintf((char*)"Current Running Min = %.5f \n", (double) this->RunningMin); 
	    R_FlushConsole();
     }
     return(this->RunningMin);        
}

///////////////////////////////////////////////////////////////////////
//  LARSObject::GenMecc
//     "cc" is the current vector of X^T * Residuals
//     This chooses what the current active factors should be.
 int LARSObject::GenMecc(int jj) {
	 double InsertMe = 0;
	 int OnFun, ii;
	 if (CurrentGoat == jj) {
		 CurrentGoat = jj;
		 cc[jj] = 0;
		 return(1);
     }
	 if (FakeFlag < 0 || (this->GFlag == 1 || this->GFlag == 3)) {
		 InsertMe = this->XTYAll[jj];
		 if (WDiagFlag == 0) {
			 OnFun = jj * this->kLen;
			 InsertMe = InsertMe - this->GAllOnB[jj];
         } else if (WDiagFlag == 1) {
	         OnFun = jj * this->kLen;
	         InsertMe = (this->XTYAll[jj]  - this->GAllOnB[jj]) * this->WDiag[jj];
	     } else if (WDiagFlag == 2) {
		     Rprintf((char*)"GenMecc You don't want to deal with NonDiagonal W at this point\n");
		     R_FlushConsole();  
		     return(-1);  
	     }
	     this->cc[jj] = InsertMe;
	     if (ISNAN(this->cc[jj])) {
		     Rprintf((char*)"this->cc[%d] is NAN, Flag in GenMeCC\n", jj);
		     R_FlushConsole();
		     return(-1);
	     }
	     return(1);
     }
     if (FakeFlag >= 0 && (this->GAll == NULL || this->GFlag == 0)) {
        if ( WDiagFlag == 0) {
		    this->cc[jj] = (double) 0;
			    OnFun = this->NLen * jj;
		        for (ii = 0; ii < this->NLen; ii++) { 
			       this->cc[jj] += (double) this->xxs[OnFun] * (double) this->OnResid[ii];
			       OnFun++;
		        }
	    } else if (WDiagFlag == 1) {
		     this->cc[jj] = (double) 0;
		      OnFun = this->NLen * jj;
		        for (ii = 0; ii < this->NLen; ii++) { 
			       this->cc[jj] += (double) this->xxs[OnFun] * (double) this->OnResid[ii];
			       OnFun++;
		        } 
		        this->cc[jj] = this->cc[jj] * (double) this->WDiag[jj];
	    } else if (WDiagFlag == 2) {
		     Rprintf((char*)"GenMecc You don't want to deal with NonDiagonal W at this point\n");
		     R_FlushConsole();
	    }	     
        if (ISNAN(this->cc[jj])) {
		     Rprintf((char*)"this->cc[%d] is NAN, Flag in GenMeCC\n", jj);
		     R_FlushConsole();
		     return(-1);
	     }
	     return(1);
    }
    Rprintf((char*)"GenMecc: Should not be here \n");
    R_FlushConsole();
    return(-1);	 
}

////////////////////////////////////////////////////////////////////////////////////
// int LARSObject::ElseCC(int PrevGoal, int GotGoal)
//
//   In the case that the LARS algorithm defies the usual model expectations and
//     loads +/- 2 new covariates in a step, the ElseCC operation is called to calculate
//     all residual values a little more closely and assess the situation whether 
//     this was a correct operation.  Given that datasets with k >> n are sort of guaranteed to
//     have certain non-linearly independent qualities, this function can expect to be called
//     on datasets of that structure, and does not guarantee a flaw in the run of the algorithm.
//
 int LARSObject::ElseCC(int PrevGoal, int GotGoal) {
	    double EpsilonNM = .00000000001;
	    
	    int cnt1;   double DiagMult;
	    double SumYTXOnB, SumOnBXTXOnB;
	    if (WDiagFlag > 0) { EpsilonNM = EpsilonNM * AvWeights; }
	    int OnFun = 0, jj;
        int ii, ii1, ii2, cno1;
	    //int ii;
	    this->OnL = 0;
	   // Rprintf((char*)"ElseCC: About Make Space for OnMuA XTYAll \n");
	   // R_FlushConsole(); R_ProcessEvents();
	    int NoXMax = this->kLen * this->NLen;
	    if (FakeFlag >= 0) {
		    for (ii2 = 0; ii2 < this->NLen; ii2++) {
			     OnMuA[ii2] = 0; OnResid[ii2] = 0;
		    }
	    }
	    if (FakeFlag >= 0 && XTYAll != NULL) {
		    for (ii1 = 0; ii1 < this->kLen; ii1++) {
			     XTYAll[ii1] = 0.0;
		    }
        }
	   // Rprintf((char*)"ElseCC: about to run Now go fill OnMuA, XTYAll\n");
	   // R_FlushConsole();
	   // R_ProcessEvents();
	   if (FakeFlag >= 0) {
		      if (WDiagFlag == 1) {
		      cno1 = 0;
			      for (ii1 = 0; ii1 < this->kLen; ii1++) {
				      if (this->OnBetas[ii1] != 0) {
					      for (ii2 = 0; ii2 < this->NLen; ii2++) {
			
						          this->OnMuA[ii2] += 
						            (double) (this->WDiag[ii1]) * OnBetas[ii1] * this->xxs[cno1];
					          cno1++;
				          }
		              } else { cno1 += this->NLen; }
		          }
	          } else {
			      cno1 = 0;
			      for (ii1 = 0; ii1 < this->kLen; ii1++) {
				      if (this->OnBetas[ii1] != 0) {
					      for (ii2 = 0; ii2 < this->NLen; ii2++) {					   
						          this->OnMuA[ii2] += (double) (this->OnBetas[ii1]) * this->xxs[cno1];
	
					          cno1++;
				          }
		              } else { cno1 += this->NLen; }	          
	               }
	          }
         }
          if (FakeFlag >= 0 && XTYAll != NULL) {
	   //     Rprintf((char*)"Filling XTYAll \n");
	   //     R_FlushConsole();
	        cno1 = 0;
	        for (ii1 = 0; ii1 < this->kLen; ii1++) {
		        for (ii2 = 0; ii2 < this->NLen; ii2++) {
			        XTYAll[ii1] += (double) this->yys[ii2] * this->xxs[cno1];
			        cno1++;
		        }
	        }    
	    //      Rprintf((char*)"ElseCC: Reweighting XTYAll\n");
	    //      R_FlushConsole();
	          if (WDiagFlag == 1) {
		           for (ii1 = 0; ii1 < this->kLen; ii1++) {
			            XTYAll[ii1] *= WDiag[ii1];
		           }
	          }          
          }

        //  Rprintf((char*)"ElseCC: About to Look for SumCurResids\n");
        //  R_FlushConsole();
        if (FakeFlag >= 0) {
	          SumCurResids = 0.0;
	          for (ii2 = 0; ii2 < this->NLen; ii2++) {
		              this->OnResid[ii2] =  this->yys[ii2] - this->OnMuA[ii2];
		              SumCurResids += this->OnResid[ii2] * this->OnResid[ii2];
	          }
        } else {
			   SumYTXOnB = 0;
			   SumOnBXTXOnB = 0;
			   if (WDiagFlag == 0) {
				   for (ii = 0; ii < kLen; ii++) {GAllOnB[ii] = 0; }
				   for (jj = 0; jj < kLen; jj++) {
					 cnt1 = jj * kLen;
					 if (OnBetas[jj] == 0) {
			         } else {
				         for (ii = 0; ii < kLen; ii++) {
					         GAllOnB[ii] +=GAll[cnt1 + ii] * OnBetas[jj];
					        
				         }
				         SumYTXOnB +=  XTYAll[jj] * OnBetas[jj];
			         }
			       }
			       for (ii = 0; ii < kLen; ii++) { SumOnBXTXOnB += GAllOnB[ii] * OnBetas[ii]; }
		      } else {
				 for (ii = 0; ii < kLen; ii++) {GAllOnB[ii] = 0; }
				     for (jj = 0; jj < kLen; jj++) {
					 cnt1 = jj * kLen;
					 DiagMult = WDiag[jj] * OnBetas[jj];
					 if (OnBetas[jj] == 0) {
			         } else {
				         for (ii = 0; ii < kLen; ii++) {
					         GAllOnB[ii] +=GAll[cnt1 + ii] * DiagMult;
					        
				         }
				         SumYTXOnB +=  XTYAll[jj] * WDiag[jj] * OnBetas[jj];
			         }
			       }
			       //for (ii = 0; ii < kLen; ii++) {GAllOnB[ii] = GAllOnB[ii] * WDiag[ii]; }     
		           for (ii = 0; ii < kLen; ii++) { SumOnBXTXOnB += GAllOnB[ii] * WDiag[ii] * OnBetas[ii]; }
		       }
		       SumCurResids = SumYYSq - 2 * SumYTXOnB + SumOnBXTXOnB;	   	         
        } 
        //  Rprintf((char*)"ElseCC: Calculating cc\n");
        //  R_FlushConsole();   
          if (FakeFlag >= 0 && ((GFlag == 0 || GFlag == 2) || this->NLen < this->kLen)) {
	        cno1=0;
	        for (ii1 = 0; ii1 < this->kLen; ii1++) {
		        cc[ii1] = 0;
		        for (ii2 = 0; ii2 < this->NLen; ii2++) {
			        if (WDiagFlag == 1) {
				       cc[ii1] += ((double) WDiag[ii1] * xxs[cno1] * OnResid[ii2]);  
			        } else {
				       cc[ii1] +=((double)  xxs[cno1] * OnResid[ii2]);    
			        }
			        cno1++;
		        } 	            
            }  
          } else {
	          if (WDiagFlag == 1) {
		           cno1 = 0;
		          for (ii1 = 0; ii1 < this->kLen; ii1++) {
	                  cc[ii1] = XTYAll[ii1];		        
			          for (ii2 = 0; ii2 < this->kLen; ii2++) {
					         cc[ii1] -= ((double) WDiag[ii1] * WDiag[ii2] *
                                        GAll[cno1] * OnBetas[ii2]);
				          cno1++;   
			          }
		          }  
	          } else {
		           cno1 = 0;
		          for (ii1 = 0; ii1 < this->kLen; ii1++) {
	                  cc[ii1] = XTYAll[ii1];		        
			          for (ii2 = 0; ii2 < this->kLen; ii2++) {			      
					      cc[ii1] -= ((double) GAll[cno1] * OnBetas[ii2]);    			       
				        cno1++;   
			          }
		          }  		          
	          }
          }
        //  Rprintf((char*)"ElseCC: About to look through ccjj \n");
        //  R_FlushConsole();
	    //this->PrintFlag = 5;  
        //if (this->PrintFlag > 6 && this->GFlag == 0) {
	    //    Rprintf( (char*) "LLA: printing current OnResid \n");
	    //    for (jj = 0; jj < this->NLen; jj++) {
	    //    Rprintf( (char*)"         %d = %.6f \n", jj, (double) this->OnResid[jj]);
        //}
        //R_FlushConsole();
        //}
	    for (jj = 0; jj < kLen; jj++) {
		/*	  if (this->PrintFlag > 6) {
		         Rprintf( (char *) "LLA: ElseCC, jj = %d /%d\n", jj, this->kLen);
		         R_FlushConsole();	 
		      }*/  		    
		    if (NoXMax <= OnFun) {
			//    Rprintf( (char *) "We're going to steal from too much xxs: NoXMax = %d, OnFun = %d\n",
			//                 NoXMax, OnFun);
			    R_FlushConsole();
		    }
	
	        if (ISNAN(cc[jj])) {
		    //    Rprintf((char*) "LLA: ElseCC, cc[%d] is nan Flaw \n", jj);
		    //    R_FlushConsole();
		        return(-4);
	        }
	        //Rprintf((char*) "LLA: ElseCC, jj = %d/%d, cc = %.4f\n", 
	        //                jj, this->kLen, (double) this->cc[jj]);
	        //R_FlushConsole();
	        if (this->cc[jj] < EpsilonNM && this->cc[jj] > EpsilonNM) {
		        
            } else if (this->OnL == 0) {
		        this->ActiveOns[this->OnL] = jj;
		        this->ActiveSigns[this->OnL] =  (int) GetSign((double) this->cc[jj]);
		        this->CurMaxs[this->OnL] = (double) this->cc[jj];
		        this->OnL++; 
	        } else if (( this->cc[jj] > 0 && this->CurMaxs[0] > 0 && 
	                    this->cc[jj] > this->CurMaxs[0] + EpsilonNM) ||
	                   ( this->cc[jj] > 0 && this->CurMaxs[0] < 0 &&
	                     this->cc[jj] > -this->CurMaxs[0] + EpsilonNM)) {
		          this->OnL = 0;
		          this->ActiveOns[this->OnL] = jj;
		          this->CurMaxs[this->OnL] = (double) this->cc[jj];
		          this->ActiveSigns[this->OnL] = 1;
		          this->OnL++;
	        } else if (( this->cc[jj] > 0 && this->CurMaxs[0] > 0 && 
	                    -this->cc[jj] > this->CurMaxs[0] + EpsilonNM) ||
	                   ( this->cc[jj] < 0 && this->CurMaxs[0] < 0 &&
	                     -this->cc[jj] > -this->CurMaxs[0] + EpsilonNM)) {  
	              this->OnL = 0;
		          this->ActiveOns[this->OnL] = jj;
		          this->CurMaxs[this->OnL] = (double) this->cc[jj];
		          this->ActiveSigns[this->OnL] = 1;
		          this->OnL++;
	        } else if (( this->cc[jj] > 0 && this->CurMaxs[0] > 0 && 
	                    this->cc[jj] > this->CurMaxs[0] - EpsilonNM) ||
	                   ( this->cc[jj] > 0 && this->CurMaxs[0] < 0 &&
	                     this->cc[jj] > -this->CurMaxs[0] - EpsilonNM)) {
		        this->ActiveOns[this->OnL] = jj;
		        this->ActiveSigns[this->OnL] = -1;
		        this->CurMaxs[this->OnL]  = this->cc[jj];
		        this->OnL++;
	        } else if (( this->cc[jj] < 0 && this->CurMaxs[0] > 0 && 
	                    -this->cc[jj] > this->CurMaxs[0] - EpsilonNM) ||
	                   ( this->cc[jj] < 0 && this->CurMaxs[0] < 0 &&
	                     -this->cc[jj] > -this->CurMaxs[0] - EpsilonNM)) {
		        this->ActiveOns[this->OnL] = jj;
		        this->ActiveSigns[this->OnL] = -1;
		        this->CurMaxs[this->OnL]  = this->cc[jj];
		        this->OnL++;
	        }		        		                     
        }
          //Rprintf((char*)"ElseCC: About to look through ccjj \n");
          //R_FlushConsole();
	  /*if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: ElseCC, finished Loop checking OnL = %d/%d\n", 
                      this->OnL, this->kLen);
         R_FlushConsole();	 
      } */  
        if (this->OnL > PrevGoal) {
	        Rprintf((char*)"ElseCC Error, this->OnL = %d, but PrevGoal = %d\n", 
	           this->OnL, PrevGoal);
	        if (this->OnL  == GotGoal) {
		        Rprintf((char*)"GotGoal is still %d \n", GotGoal);
	        } else {
		        Rprintf((char*)"But GotGoal was %d \n", GotGoal);
	        }
	        R_FlushConsole();
	        //Rprintf((char*)"By the way, sending out the Betas, cjj\n");
	        //Rprintf((char*)"First 5 max cjj\n");
	        //ii1 = 5; if (ii1 > this->OnL) { ii1 = this->OnL; }
	        //Rprintf((char*)" Those Indices \n");
	        //for (ii = 0; ii < ii1; ii++) {
		    //    Rprintf((char*)" %.d  ", (double) ActiveOns[ii] );
	        //}
	        //R_FlushConsole(); Rprintf((char*) "The cjj \n");
	        //for (ii = 0; ii < ii1; ii++) {
		    //    Rprintf((char*)" %.7f  ", (double) this->cc[ActiveOns[ii]] );
	        //}
	        //Rprintf((char*)" \n " ); R_FlushConsole(); R_ProcessEvents();
	        //Rprintf((char*)"Corresponding Betas \n");
	        //for (ii = 0; ii < ii1 ; ii++) {
		    //    Rprintf((char*)" %. 5f   ", (double) this->OnBetas[ActiveOns[ii]] );
	        //}   
	        //R_FlushConsole();
	        return(-1);
        }       
        if (this->OnL > this->kLen) {
	        Rprintf( (char *) "Error OnL > kLen\n");
	        return(-1);
        }
      /*if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: ElseCC, filling OnActives, wa \n");
         R_FlushConsole();	 
      } */ 
      // Rprintf((char*)"ElseCC: Configuring wa \n");
      // R_FlushConsole();
        for (jj = 0; jj < this->kLen; jj++) {
	        this->OnActives[jj] = 0;
	        this->wa[jj] = 0; 
        }
      /*if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: ElseCC, filling ActiveSigns \n");
         R_FlushConsole();	 
      } */
        for (jj = 0; jj < this->OnL; jj++) {
	        this->OnActives[ this->ActiveOns[jj] ] = this->ActiveSigns[jj];
        }
      /*if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: ElseCC, InActive \n");
         R_FlushConsole();	 
      } */        
        this->OnInL = 0;
        for (jj = 0; jj < this->kLen; jj++) {
	        if (this->OnActives[jj] == 0) {
		        this->AAInActive[this->OnInL] =  jj;
		        this->OnInL++;
	        }
        }
        /*if (this->OnInL > this->kLen) {
	        Rprintf( (char *) "Problem OnInL[0] = %d > kLen = %d\n", this->OnInL, this->kLen);
	        R_FlushConsole();
	        R_ProcessEvents();
        }*/
        /*if (this->OnL <= 0) {
		    Rprintf( (char *) "ElseCC flaw: OnL == 0 \n");
		    R_FlushConsole();
		    R_ProcessEvents();
	    } */
	    if (PrintFlag > 2) {
		    R_FlushConsole();
		    Rprintf( (char *) "Found OnL = %d, Printing cc \n cc = c(", (int) this->OnL);
		    for (ii = 0; ii < kLen; ii++) {
			    Rprintf( (char *) " %.4f", (double) cc[ii]);
			    if (ii == kLen-1) {Rprintf((char*)" ); \n");}
			    else { Rprintf((char*)", "); }
		    }
	    }
	    if (PrintFlag > 5) {
		    Rprintf( (char *) "\n");
		    Rprintf( (char *) "Printing OnActives \n");
		    for (ii = 0; ii < kLen; ii++) {
			    Rprintf( (char *) " %d = %d, ", ii, this->OnActives[ii]);
		    }
		    R_FlushConsole();
		    Rprintf( (char *) "\n");
		    Rprintf( (char *) "Printing ActiveOns \n");
		    for (ii = 0; ii < this->OnL; ii++) {
			    Rprintf( (char *) " %d = %d, ", ii, this->ActiveOns[ii]);
		    }
		     Rprintf( (char *) "\n");
		    Rprintf( (char *) "Printing ActiveSigns \n");
		    for (ii = 0; ii < this->OnL; ii++) {
			    Rprintf( (char *) " %d = %d, ", ii, this->ActiveSigns[ii]);
		    }
		    Rprintf( (char *) "\n");	  
        } 
       // Rprintf("ElseCC: Finishing Strong\n");
       // R_FlushConsole();  
       // R_ProcessEvents();
        return(1);
     }
  
 ////////////////////////////////////////////////////////////////////////////////
 //   int AssessCC()
 //
 //
 //    At beginning of LARS step, this calculates the "cc", correlated residuals vector
 //   then chooses the active set based upon shared maximum cc values.  Once C and the
 //   active set is picked, this also loads the sign vector for the active set.
 //    
 int LARSObject::AssessCC() {
	    int OnFun = 0, jj;
        int ii;
	    //int ii;
	    this->OnL = 0;
	    int NoXMax = this->kLen * this->NLen;
	    //this->PrintFlag = 5;
	  if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: AssessCC \n");
         R_FlushConsole();	 
      }  
        if (this->PrintFlag > 6 && this->GFlag == 0) {
	        Rprintf( (char*) "LLA: printing current OnResid \n");
	        for (jj = 0; jj < this->NLen; jj++) {
	        Rprintf( (char*)"         %d = %.6f \n", jj, (double) this->OnResid[jj]);
        }
        R_FlushConsole();
        }
        double EpsilonUBM = EpsilonBM;
        if (WDiagFlag == 1) {
	         EpsilonUBM = EpsilonBM * AvWeights;
        }
	    for (jj = 0; jj < kLen; jj++) {
		/*	  if (this->PrintFlag > 6) {
		         Rprintf( (char *) "LLA: AssessCC, jj = %d /%d\n", jj, this->kLen);
		         R_FlushConsole();	 
		      }*/  		    
		    if (NoXMax <= OnFun) {
			    Rprintf( (char *) "We're going to steal from too much xxs: NoXMax = %d, OnFun = %d\n",
			                 NoXMax, OnFun);
			    R_FlushConsole();
		    }
		    GenMecc(jj);
	
	        if (ISNAN(cc[jj])) {
		        Rprintf((char*) "LLA: AssessCC, cc[%d] is nan Flaw \n", jj);
		        R_FlushConsole();
		        return(-4);
	        }
	        //Rprintf((char*) "LLA: AssessCC, jj = %d/%d, cc = %.4f\n", 
	        //                jj, this->kLen, (double) this->cc[jj]);
	        //R_FlushConsole();
	        if (this->OnL == 0) {
		        this->ActiveOns[this->OnL] = jj;
		        this->ActiveSigns[this->OnL] =  (int) GetSign((double) this->cc[jj]);
		        this->CurMaxs[this->OnL] = (double) this->cc[jj];
		        this->OnL++; 
	        } else if ( (double) fabs((double) this->cc[jj]) - (double) fabs( 
	                      (double) this->CurMaxs[0]) > (double) EpsilonUBM ) {
		        this->OnL = 0;
		        this->ActiveOns[this->OnL] = jj;
		        this->CurMaxs[this->OnL] = (double) this->cc[jj];
		        this->ActiveSigns[this->OnL] = GetSign((double) this->cc[jj]);
		        this->OnL++;
	        } else if ((double) fabs((double) this->cc[jj]) - 
	                   (double) fabs((double) this->CurMaxs[0]) > (double) -EpsilonUBM ) {
		        this->ActiveOns[this->OnL] = jj;
		        this->ActiveSigns[this->OnL] = GetSign((double)cc[jj]);
		        this->CurMaxs[this->OnL]  = this->cc[jj];
		        this->OnL++;
	        }
        }
	  /*if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: AssessCC, finished Loop checking OnL = %d/%d\n", 
                      this->OnL, this->kLen);
         R_FlushConsole();	 
      } */         
        if (this->OnL > this->kLen) {
	        Rprintf( (char *) "Error OnL > kLen\n");
	        return(-1);
        }
      /*if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: AssessCC, filling OnActives, wa \n");
         R_FlushConsole();	 
      } */ 
        for (jj = 0; jj < this->kLen; jj++) {
	        this->OnActives[jj] = 0;
	        this->wa[jj] = 0; 
        }
      /*if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: AssessCC, filling ActiveSigns \n");
         R_FlushConsole();	 
      } */
        for (jj = 0; jj < this->OnL; jj++) {
	        this->OnActives[ this->ActiveOns[jj] ] = this->ActiveSigns[jj];
        }
      /*if (this->PrintFlag > 4) {
         Rprintf( (char *) "LLA: AssessCC, InActive \n");
         R_FlushConsole();	 
      } */        
        this->OnInL = 0;
        for (jj = 0; jj < this->kLen; jj++) {
	        if (this->OnActives[jj] == 0) {
		        this->AAInActive[this->OnInL] =  jj;
		        this->OnInL++;
	        }
        }
        /*if (this->OnInL > this->kLen) {
	        Rprintf( (char *) "Problem OnInL[0] = %d > kLen = %d\n", this->OnInL, this->kLen);
	        R_FlushConsole();
	        R_ProcessEvents();
        }*/
        /*if (this->OnL <= 0) {
		    Rprintf( (char *) "AssessCC flaw: OnL == 0 \n");
		    R_FlushConsole();
		    R_ProcessEvents();
	    } */
	    if (PrintFlag > 2) {
		    R_FlushConsole();
		    Rprintf( (char *) "Found OnL = %d, Printing cc \n cc = c(", (int) this->OnL);
		    for (ii = 0; ii < kLen; ii++) {
			    Rprintf( (char *) " %.4f", (double) cc[ii]);
			    if (ii == kLen-1) {Rprintf((char*)" ); \n");}
			    else { Rprintf((char*)", "); }
		    }
	    }
	    if (PrintFlag > 5) {
		    Rprintf( (char *) "\n");
		    Rprintf( (char *) "Printing OnActives \n");
		    for (ii = 0; ii < kLen; ii++) {
			    Rprintf( (char *) " %d = %d, ", ii, this->OnActives[ii]);
		    }
		    R_FlushConsole();
		    Rprintf( (char *) "\n");
		    Rprintf( (char *) "Printing ActiveOns \n");
		    for (ii = 0; ii < this->OnL; ii++) {
			    Rprintf( (char *) " %d = %d, ", ii, this->ActiveOns[ii]);
		    }
		     Rprintf( (char *) "\n");
		    Rprintf( (char *) "Printing ActiveSigns \n");
		    for (ii = 0; ii < this->OnL; ii++) {
			    Rprintf( (char *) " %d = %d, ", ii, this->ActiveSigns[ii]);
		    }
		    Rprintf( (char *) "\n");	  
        } 
        R_FlushConsole();  
        R_ProcessEvents();
        return(1);
     }
     
     
     

////////////////////////////////////////////////////////////////
//  Create Wa
//
//   Wa = A * IGGA * 1
//    this is the unit increase vector for "Beta" during this step.
int LARSObject::CreateWa() {
	    //Rprintf( (char *) "Creating wa \n");
	    //R_FlushConsole();
        int OnFun= 0,jj, ii;
        if (WDiagFlag == 1) {
	        for (jj = 0; jj < this->OnL; jj++) {
		        this->wa[jj] = 0;
		        for (ii = 0; ii < this->OnL; ii++)  {
			        this->wa[jj] +=  (double) this->IGGA[OnFun] /
			                         ((double) WDiag[ActiveOns[ii]]);
			        OnFun++;
		        }
		        this->wa[jj] = (double) this->wa[jj]  / WDiag[ActiveOns[jj]];
		        this->wa[jj] = (double) this->wa[jj] * (double) this->AAA;
	       }
       } else {
            for (jj = 0; jj < this->OnL; jj++) {
		        this->wa[jj] = 0;
		        for (ii = 0; ii < this->OnL; ii++)  {
			        this->wa[jj] +=  (double) this->IGGA[OnFun];
			        OnFun++;
		        }
		        this->wa[jj] = (double) this->wa[jj] * (double) this->AAA;
	       }	       
       }
       //if ( this->OnL < 8) {
       //Rprintf( (char *) "Printing wa \n");
       //for (jj = 0; jj < OnL; jj++) {
	   //     Rprintf( (char *) "  %.4f ", (double) wa[jj]);
       //}
       // Rprintf( (char *) " \n");
       //}
       return(1);
}   

///////////////////////////////////////////////////////////////////////////////////////
//  int LARSObject::Makeaa()
//
//   "aa" is the vector used for the minimum Gamma finding.
//   It relates the perpendicularness of wa to the the full covariate set.
//
int LARSObject::Makeaa() {
	   int ii, jj, cnt1, tii;
	   //if (GFlag == 0) {
	   //  SetMeXIn(this->kLen, this->NLen, this->OnL, this->ActiveOns, this->ActiveSigns, 
	   //                 this->xxs, this->XXA);	   
       //}
       if (FakeFlag < 0 || GFlag == 1 || GFlag == 3) {
	      for (ii = 0; ii < this->kLen; ii++) {
		       aa[ii] = 0;
	      }
	      if (WDiagFlag == 0) {
		      if (PrintFlag > 3) {
	          Rprintf((char*) "aa Calc With WDiagFlag = 0, GFlag = 1 \n");
	          R_FlushConsole();
              }
		    tii = 0;
	        for (tii = 0; tii < this->OnInL; tii++)   {
		        ii = AAInActive[tii];		       
	            cnt1 = ii * this->kLen;
	            if (this->OnL == 1) {
		            this->aa[ii] = ((double) this->GAll[ cnt1 + ActiveOns[0] ]) * 
		                           (double) this->ActiveSigns[0] * 
		                           (double) this->wa[0];
	            } else {
			        for (jj = 0; jj < this->OnL; jj++) {
				        this->aa[ii] =  (double) this->aa[ii] + 
				              ((double) this->GAll[cnt1+ ActiveOns[jj] ])*
	                          (double) this->wa[jj] *
				              (double)  this->ActiveSigns[jj] 
				             ;			                 
			        }
	            } 
	         }
	         if (CurrentGoat > 0) {
		         aa[CurrentGoat] = 0.0;
		         CurrentGoat = -20;
             }	         
	         CalcGAwa();
	        // for (ttii = 0; ttii < this->OnL; ttii++) {
		    //   ii = ActiveOns[ttii];	           
	        // }
	          return(1);
          } else if (WDiagFlag == 1) {
	          if (PrintFlag > 3) {
	          Rprintf((char*) "aa Calc With WDiagFlag = 1, GFlag = 1 \n");
	          R_FlushConsole();
              }
            for (tii = 0; tii < this->OnInL; tii++)   {
	            ii = AAInActive[tii];
	            cnt1 = ii * this->kLen;
	            if (this->OnL == 1) {
		            this->aa[ii] = ((double) this->GAll[ cnt1 + ActiveOns[0]] * 
		                (double) this->WDiag[ActiveOns[0]] *
		               (double) this->ActiveSigns[0] * (double) this->wa[0]);
	            } else {
			        for (jj = 0; jj < this->OnL; jj++) {
				        this->aa[ii] = (double) this->aa[ii] +
				                 ((double) this->GAll[cnt1  + ActiveOns[jj] ]) 
				        			                 * (double) this->WDiag[ActiveOns[jj]] 
				                 * (double) this->wa[jj] 
				                 * (double)  this->ActiveSigns[jj]
	                            ;
			        }
	            }
		        this->aa[ii] = this->aa[ii] * this->WDiag[ii];
	         }
	         if (CurrentGoat > 0) {
		         aa[CurrentGoat] = 0.0;
		         CurrentGoat = -20;
             }
	         CalcGAwa();
	         return(1);          	          
          }  else if (WDiagFlag >= 2) {
	          Rprintf((char*)"Makeaa: Error cannot get W Non Diag \n");
	          return(-1);
          }
      }
    if (WDiagFlag >= 2) {
	    Rprintf((char*)"Makeaa: Cannot do W Non Diag\n");
	    return(-1);
    }
    if (WDiagFlag == 1) {
	    	  if (PrintFlag > 3) {
	          Rprintf((char*) "aa Calc With WDiagFlag = 1, GFlag = 0 \n");
	          R_FlushConsole();
              }
        int nii; int cno1 = 0; int ii;
        for (nii = 0; nii < this->NLen; nii++) { this->ua[nii] = 0; }
        if (XXA == NULL) {
	        Rprintf("LarsObject:: MakeAA Error, FakeFlag >=0, XXA is NULL\n");
	        R_FlushConsole();
	        return(-1);
        }
        for (ii = 0; ii < this->OnL; ii++) {
	       for (nii = 0; nii < this->NLen; nii++) {
		       this->ua[nii] += this->XXA[cno1] * WDiag[ActiveOns[ii]] * this->wa[ii];
		       cno1++;    
	       }    
        }
	    double SSua = 0; 
	    for (nii = 0; nii < this->NLen; nii++) {
		    SSua += this->ua[nii] * this->ua[nii];
	    }
	    if (SSua < 1 - EpsilonBM || SSua > 1 + EpsilonBM ) {
		    Rprintf("Warning Sum Squared ua = %.8E \n");
		    R_FlushConsole();
	    }
    } else {
	          if (PrintFlag > 3) {
	          Rprintf((char*) "aa Calc With WDiagFlag = 0, GFlag = 0 \n");
	          R_FlushConsole();
              }
	    MatTimesVec(this->OnL, this->NLen, this->XXA, this->wa, this->ua);
    }
    tMatTimesVec(this->kLen, this->NLen, this->xxs, this->ua, this->aa);
    if (WDiagFlag == 1) {
	    for (ii = 0; ii < this->kLen; ii++) {
		    aa[ii] = aa[ii] * (double) this->WDiag[ii];
	    }
    }
    if (CurrentGoat > 0) {
	    aa[CurrentGoat] = 0.0;
	    CurrentGoat = -20;
    }
    return(1);
}

////////////////////////////////////////////////////////////////////////////////
//  CalcGAwa()
//
//   Calculates the "aa" vector, useful for assessing the "gamma" in the data.
//   The weighting operation is important because GAll is unweighted matrix.
//   Most of the code is garbage because we have reworded some definitions.
int LARSObject::CalcGAwa() {
  int ii,jj;
  if (PrintFlag > 3) {
	  Rprintf((char*)"Calculating GAwa \n");
      R_FlushConsole();
  }
  if (aa == NULL) {
	  Rprintf("CalcGAwa: aa is NULL cannot follow\n");
	  R_FlushConsole();
	  return(-1);
   }
  for (ii = 0; ii < OnL; ii++) {
	  aa[this->ActiveOns[ii]] = AAA * ActiveSigns[ii];
  }
  return(1);
  /// Garbage Code::
  if (WDiagFlag == 0) {
  for (ii = 0; ii < OnL; ii++) {
	  aa[this->ActiveOns[ii]] = AAA;
  }		
  for (ii = 0; ii < OnL;ii++) {
    if (ActiveSigns[ii] < 0) {
	   for (jj = 0; jj < OnL; jj++) {
	     aa[this->ActiveOns[jj]] += -2* GAll[ ii * kLen + jj] * wa[ii];	   		   
       }    
    }	  
  }
  for (ii = 0; ii < OnL; ii++) {
     if (ActiveSigns[ii] < 0) {
	    aa[this->ActiveOns[ii]] = -1.0 * aa[this->ActiveOns[ii]];    
     }	  
  }
  } else {
	for (ii = 0; ii < OnL; ii++) {
		  aa[this->ActiveOns[ii]] = AAA;
	  }		
	  for (ii = 0; ii < OnL;ii++) {
	    if (ActiveSigns[ii] <0) {
		   for (jj = 0; jj < OnL; jj++) {
		     aa[this->ActiveOns[jj]] += -2* GAll[ ii * kLen + jj] * WDiag[ii] * WDiag[jj] * wa[ii];	   		   
	       }    
	    }	  
	  }
	  for (ii = 0; ii < OnL; ii++) {
	     if (ActiveSigns[ii] < 0) {
		    aa[this->ActiveOns[ii]] = -1.0 * aa[this->ActiveOns[ii]];    
	     }	  
	  }	   
 }
  return(1);
}
   

//////////////////////////////////////////////////////////////
//  The "Goat" is the variable getting sacrificed in a LARS step
//   that happens to involve an eliminated factor.  The Goat is 
//   not allowed to be featured and brought back to life in the next step,
//   though it can be returned to play afterwards.  Because a Goat tends
//   to turn a variable exactly "off", one must do some recalculations of 
//   feature vectors because the desire is to make sure that the active sets
//   in OnBetas, YTXOnB, etc are maintained.
//   
//  This does not allow yet for more than one Goat.  According to the statistical
//  theory, this should not occur, but one hopes that the dataset does not contain
//  unfortunate non-independence that arises with this.
//
//
 int LARSObject::CalcWithGoat( double FeatConst) {
	 /*if (PrintFlag > 3) {
	    Rprintf((char*)"Running the Goat Calculation\n");
	    R_FlushConsole();
     }*/
	 int ii, jj, cnt1;
	 double SumYTXOnB;
	 double SumOnBXTXOnB;
	 double Diagjj;
	 if ((GFlag ==1 || GFlag == 3) && GAllOnB == NULL) {
		 Rprintf((char*)"LARSObject:: CalcWithGoat, Error, GAllOnB is NULL\n");
		 R_FlushConsole(); return(-1);
     }
	 if ((GFlag ==1 || GFlag == 3) && GAll == NULL) {
		 Rprintf((char*)"LARSObject:: CalcWithGoat, Error, GAll is NULL\n");
		 R_FlushConsole(); return(-1);
     }     
	       for (jj = 0; jj < this->OnL; jj++) {
	           this->OnBetas[ this->ActiveOns[jj]] = 
	                 (double) this->OnBetas[ this->ActiveOns[jj]] +(double) FeatConst * 
	                 (double) this->ActiveSigns[ jj ] * (double) this->wa[jj];	
	           if (CurrentGoat == ActiveOns[jj]) {
		           this->OnBetas[this->ActiveOns[jj]] = 0;
	           }            
           }
   if (FakeFlag < 0 || (GFlag == 1 || GFlag == 3)) {
	   SumYTXOnB = 0;
	   SumOnBXTXOnB = 0;
	   if (WDiagFlag == 0) {
		   for (ii = 0; ii < kLen; ii++) {GAllOnB[ii] = 0; }
		   for (jj = 0; jj < kLen; jj++) {
			 cnt1 = jj * kLen;
			 if (OnBetas[jj] == 0) {
	         } else {
		         for (ii = 0; ii < kLen; ii++) {
			         GAllOnB[ii] +=GAll[cnt1 + ii] * OnBetas[jj];
			        
		         }
		         SumYTXOnB +=  XTYAll[jj] * OnBetas[jj];
	         }
	       }
	       for (ii = 0; ii < kLen; ii++) { SumOnBXTXOnB += GAllOnB[ii] * OnBetas[ii]; }
      } else {
		 for (ii = 0; ii < kLen; ii++) {GAllOnB[ii] = 0; }
		     for (jj = 0; jj < kLen; jj++) {
			 cnt1 = jj * kLen;
			 Diagjj = WDiag[jj] * OnBetas[jj];
			 if (OnBetas[jj] == 0) {
	         } else {
		         for (ii = 0; ii < kLen; ii++) {
			         GAllOnB[ii] +=GAll[cnt1 + ii] * Diagjj;
			        
		         }
		         SumYTXOnB +=  XTYAll[jj] * WDiag[jj] * OnBetas[jj];
	         }
	       }
	       //for (ii = 0; ii < kLen; ii++) {GAllOnB[ii] = GAllOnB[ii] * WDiag[ii]; }     
           for (ii = 0; ii < kLen; ii++) { SumOnBXTXOnB += GAllOnB[ii] * WDiag[ii] * OnBetas[ii]; }
       }
       SumCurResids = SumYYSq - 2 * SumYTXOnB + SumOnBXTXOnB;	   
	   return(1);
   }
   if (OnMuA == NULL) {
	   Rprintf("CalcWithGoat: Error FakeFlag >= 0 but OnMuA is NULL\n");
	   R_FlushConsole();  return(-1);
   }
   for (ii = 0; ii < NLen; ii++) { OnMuA[ii] = 0; }
   if (WDiagFlag == 0) { 
	   for  (jj = 0; jj < kLen; jj++) {
		   if (OnBetas[jj] != 0) {
			   cnt1 = jj * NLen;
			   for (ii = 0; ii < NLen; ii++) {
				   OnMuA[ii] += OnBetas[jj] * xxs[cnt1];
				   cnt1++;
		       }
	        }
       }
    } else {
	  for  (jj = 0; jj < kLen; jj++) {
		   if (OnBetas[jj] != 0) {
			   cnt1 = jj * NLen;
			   Diagjj = OnBetas[jj] * WDiag[jj];
			   for (ii = 0; ii < NLen; ii++) {
				   OnMuA[ii] += Diagjj * xxs[cnt1];
				   cnt1++;
		       }
	        }
       }    	    
    }
    SumCurResids = 0.0;
    for (ii = 0; ii < NLen; ii++) {
	    OnResid[ii] = yys[ii] - OnMuA[ii];
	    SumCurResids += OnResid[ii] * OnResid[ii];
    }
   return(1);	 
 }
  

////////////////////////////////////////////////////////////////////////////////
// int MakeGGA() 
//
//  The GGA
//
//
int LARSObject::MakeGGA() {
   int kOnLen = this->OnL;
   int Onnn = 0;
   int ii,jj;
   int WillPostOn;
   int WillPostOff;
   double CurTot;
   int Onn2 = 0;
   if (GGAA == NULL) {
	   Rprintf("LARSObject::MakeGGA Flaw, GGAA is NULL\n");
	   R_FlushConsole(); return(-1);
   }
   if (kOnLen == 1) {
	   GGAA[0] = 0;
	   if (FakeFlag >= 0 && GFlag == 0) {
		   	 Onn2 =  ActiveOns[0] * this->NLen;
	         for (Onnn = 0; Onnn < this->NLen; Onnn++) {
		        this->GGAA[0] += this->xxs[Onn2] * this->xxs[Onn2];
		        Onn2++;
             }
             return(1);
      } else {
	       if (GAll == NULL) {
		       Rprintf("LARSObject:MakeGGA Flaw, GAll is NULL\n");
		       R_FlushConsole(); return(-1);
	       }
	       GGAA[0] = GAll[ActiveOns[0] * this->kLen + ActiveOns[0]];
	       return(1);
      }
   }
   if (FakeFlag < 0 || GFlag == 1 || GFlag == 3) {
	       if (GAll == NULL) {
		       Rprintf("LARSObject:MakeGGA Flaw, GAll is NULL\n");
		       R_FlushConsole(); return(-1);
	       }	   
            for (ii = 0; ii < (kOnLen-1); ii++) {
			     Onnn = ii * kOnLen;
			     Onn2 = ActiveOns[ii] * this->kLen;
			     this->GGAA[Onnn + ii] = 
			         (double) this->GAll[ Onn2+ ActiveOns[ii] ];
			     for (jj = ii+1; jj < kOnLen; jj++) {
				     this->GGAA[Onnn + jj] = 
				        (double) this->GAll[ Onn2 + ActiveOns[jj]] *
				        (double) this->ActiveSigns[ii] * (double) this->ActiveSigns[jj];
				     this->GGAA[jj * kOnLen + ii] = this->GGAA[Onnn + jj];
			     }
		     }
		     this->GGAA[(kOnLen-1) * kOnLen + (kOnLen-1)] =  
		          (double) this->GAll[ ActiveOns[kOnLen-1] *
		                  this->kLen + ActiveOns[kOnLen-1] ];
		     return(1);	   
	   
   }
	       if (xxs == NULL) {
		       Rprintf("LARSObject:MakeGGA Flaw, xxs is NULL\n");
		       R_FlushConsole(); return(-1);
	       }	   
    for (ii = 0; ii < (kOnLen-1); ii++) {
		       Onn2 = this->NLen * ActiveOns[ii];
		       CurTot = 0;
		       for (Onnn = 0; Onnn < this->NLen; Onnn++) {
			       CurTot += this->xxs[Onn2] * this->xxs[Onn2];
			       Onn2++;
		       }
		       this->GGAA[ii * kOnLen + ii] = CurTot;
		       for (jj = ii+1; jj < kOnLen; jj++) {
			       CurTot = 0;
			       WillPostOn = this->ActiveOns[ii] * NLen; 
			       WillPostOff = this->ActiveOns[jj] * NLen;
			       for (Onnn = 0; Onnn < this-> NLen; Onnn++) {
				       CurTot += this->xxs[WillPostOn + Onnn] * this->xxs[WillPostOff + Onnn];
			       }
			       this->GGAA[ ii * kOnLen + jj] = CurTot * 
			                   this->ActiveSigns[ii] * this->ActiveSigns[jj];
			       this->GGAA[ jj * kOnLen + ii] = (double)  this->GGAA[ ii * kOnLen + jj];
		       }
	    }
	    Onn2 = this->NLen * ActiveOns[kOnLen-1];
	    CurTot = 0;
	    for (Onnn= 0; Onnn < this->NLen; Onnn++) {
		    CurTot += this->xxs[Onn2] * this->xxs[Onn2];
		    Onn2++;
	    }
	    GGAA[ (kOnLen-1) * kOnLen + (kOnLen-1) ] = CurTot;
	    return(1);
   if (kOnLen == 1) {
	   GGAA[0] = 0;
	   if (GFlag == 0) {
		  if (WDiagFlag == 0) { 
			 Onn2 =  ActiveOns[0] * this->NLen;
	         for (Onnn = 0; Onnn < this->NLen; Onnn++) {
		        this->GGAA[0] += this->xxs[Onn2] * this->xxs[Onn2];
		        Onn2++;
             }
             return(1);
          } else if (WDiagFlag == 1) {
	         Onn2 = ActiveOns[0] * this->NLen;
		     for (Onnn = 0; Onnn < this->NLen; Onnn++) {
		        this->GGAA[0] += this->xxs[Onn2] * this->xxs[Onn2];
		        Onn2++;
             }          
	         this->GGAA[0] = (double) this->GGAA[0] * (double) WDiag[ActiveOns[0]] * 
	                                         (double) WDiag[ActiveOns[0]]; 
	         return(1);
          } else if (WDiagFlag >= 2) {
	          Rprintf((char*)"GAll Can't deal with WXTXW just yet\n");
	          R_FlushConsole();
	          return(-1);
          }
      } if (GFlag == 1) {
	      if (WDiagFlag == 0) {
		      this->GGAA[0] = (double) this->GAll[ActiveOns[0] * this->kLen + ActiveOns[0]];
		      return(1);
	      } else if (WDiagFlag == 1) {
		      this->GGAA[0] = (double) this->GAll[ActiveOns[0] * this->kLen + ActiveOns[0]] * 
		                  WDiag[ActiveOns[0]] * WDiag[ActiveOns[0]];
		      return(1);
	      } else if (WDiagFlag >= 2) {
		      Rprintf((char*)"GAll Can't deal with WXTXW just yet\n");
		      return(-1);
	      }
      }
    }
	 if (GFlag ==1 ) {
		 if (WDiagFlag == 0) {
		     for (ii = 0; ii < (kOnLen-1); ii++) {
			     Onnn = ii * kOnLen;
			     Onn2 = ActiveOns[ii] * this->kLen;
			     this->GGAA[Onnn + ii] = 
			         (double) this->GAll[ Onn2+ ActiveOns[ii] ];
			     for (jj = ii+1; jj < kOnLen; jj++) {
				     this->GGAA[Onnn + jj] = 
				        (double) this->GAll[ Onn2 + ActiveOns[jj]] *
				        this->ActiveSigns[ii] * this->ActiveSigns[jj];
				     this->GGAA[jj * kOnLen + ii] = this->GGAA[Onnn + jj];
			     }
		     }
		     this->GGAA[(kOnLen-1) * kOnLen + (kOnLen-1)] =  
		          (double) this->GAll[ ActiveOns[kOnLen-1] *
		                  this->kLen + ActiveOns[kOnLen-1] ];
		     return(1);
	     } else if (WDiagFlag == 1) {
           for (ii = 0; ii < (kOnLen-1); ii++) {
	             Onnn = ii * kOnLen;
	             Onn2 = ActiveOns[ii] * this->kLen;
			     this->GGAA[Onnn + ii] = (double) 
			           this->GAll[ Onn2+ ActiveOns[ii]] 
			           * WDiag[ActiveOns[ii]] * WDiag[ActiveOns[ii]];
			     for (jj = ii+1; jj < kOnLen; jj++) {
				     this->GGAA[Onnn  + jj] = 
				       (double) this->GAll[ Onn2 + ActiveOns[jj]] * 
				               WDiag[ActiveOns[ii]] * WDiag[ActiveOns[jj]] *
				       this->ActiveSigns[ii] * this->ActiveSigns[jj];
				     this->GGAA[ jj * kOnLen + ii] = this->GGAA[Onnn + jj];
			     }
		     }
		     this->GGAA[(kOnLen-1) * kOnLen + (kOnLen-1)] =  
		          (double) this->GAll[ ActiveOns[kOnLen-1] *
		                  this->kLen + ActiveOns[kOnLen-1]] *
		          WDiag[ActiveOns[kOnLen-1]] * WDiag[ActiveOns[kOnLen-1]];
		     return(1);		     
	     } else if (WDiagFlag >= 2) {
		     Rprintf((char*)"Gall generate: Can't do Nondiagonal WAll yet\n");
		     R_FlushConsole();
		     return(-1);    
		     
	     }
     }
     if (kOnLen == 1) {
	     Onn2 = this->NLen * ActiveOns[0];
	     CurTot = 0;
	     for (Onnn = 0; Onnn < this->NLen; Onnn++) {
		     CurTot += this->xxs[Onn2] * this->xxs[Onn2];
		     Onn2++;
	     }
	     GGAA[0] = CurTot;
     } else {
	     for (ii = 0; ii < (kOnLen-1); ii++) {
		       Onn2 = this->NLen * ActiveOns[ii];
		       CurTot = 0;
		       for (Onnn = 0; Onnn < this->NLen; Onnn++) {
			       CurTot += this->xxs[Onn2] * this->xxs[Onn2];
			       Onn2++;
		       }
		       this->GGAA[ii * kOnLen + ii] = CurTot;
		       for (jj = ii+1; jj < kOnLen; jj++) {
			       CurTot = 0;
			       WillPostOn = this->ActiveOns[ii] * NLen; 
			       WillPostOff = this->ActiveOns[jj] * NLen;
			       for (Onnn = 0; Onnn < this-> NLen; Onnn++) {
				       CurTot += this->xxs[WillPostOn + Onnn] * this->xxs[WillPostOff + Onnn];
			       }
			       this->GGAA[ ii * kOnLen + jj] = CurTot * 
			                   this->ActiveSigns[ii] * this->ActiveSigns[jj];
			       this->GGAA[ jj * kOnLen + ii] = (double)  this->GGAA[ ii * kOnLen + jj];
		       }
	    }
	    Onn2 = this->NLen * ActiveOns[kOnLen-1];
	    CurTot = 0;
	    for (Onnn= 0; Onnn < this->NLen; Onnn++) {
		    CurTot += this->xxs[Onn2] * this->xxs[Onn2];
		    Onn2++;
	    }
	    GGAA[ (kOnLen-1) * kOnLen + (kOnLen-1) ] = CurTot;
    }
    if (WDiagFlag == 0) {
	    return(1);
    } else if (WDiagFlag == 1) {
	    if (kOnLen == 1) {
		    GGAA[0] = GGAA[0] * WDiag[ActiveOns[0]] * WDiag[ActiveOns[0]];
		    return(1);
	    }  
	    for (ii = 0; ii < kOnLen-1; ii++) {
		    WillPostOn = ii * kOnLen + ii;
		    this->GGAA[ WillPostOn] = this->GGAA[ WillPostOn ] *
		                WDiag[ActiveOns[ii]] * WDiag[ActiveOns[ii]];
		    for (jj = ii+1; jj < kOnLen; jj++) {
			    WillPostOn++;
			    this->GGAA[ WillPostOn] =
			    this->GGAA[ WillPostOn] * WDiag[ActiveOns[ii]] * WDiag[ActiveOns[jj]];
			    this->GGAA[ jj * kOnLen + ii] = this->GGAA[WillPostOn];
		    }
	    }
	    GGAA[ (kOnLen-1) * kOnLen + (kOnLen-1) ] = GGAA[(kOnLen-1) * kOnLen + (kOnLen-1) ] * 
	                     WDiag[ActiveOns[kOnLen-1]] * WDiag[ActiveOns[kOnLen-1]];
	    return(1);
    } else if (WDiagFlag == 2) {
	    Rprintf((char*)"GAll make, still can't deal with WDiagFlag == 2\n");
	    R_FlushConsole();
	    return(-1);
    }
    return(-1);
}


/////////////////////////////////////////////////////
//  LARSObject::PrintWeights()
//
//    Print weights one is currently dealing with.
//
int LARSObject::PrintWeights() {
	if (WDiag== NULL) {return(-1);}
	Rprintf((char*) "Our CurrentWeights, for lambda = %.10f\n", (double) this->lambda);
	Rprintf((char*) "WeightDiag = c( ");
	int ii;
	for (ii = 0; ii < this->kLen-1; ii++) {
		Rprintf((char*) " %.7f, ", (double) WDiag[ii] );
    }
    Rprintf((char*) " %.7f ) \n", (double) WDiag[this->kLen-1]);
    R_FlushConsole();
    R_ProcessEvents();
    return(1);
}





/////////////////////////////////////////////////////////
//  int  SetMeI(int kSLen, double GMat
//      Sets up an identity matrix into GMat
//
int LARSObject::SetMeI ( int kSLen, double *GMat ) { 
  int ii, jj, Onii = 0;
  if (GMat == NULL) {
	  Rprintf((char*)"LARSObject::SetMeI Error, GMat is NULL\n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1); 
  }
  for (	ii = 0; ii < kSLen; ii++) {
	  for (jj = 0; jj < kSLen; jj++) {
		  GMat[Onii] = 0;
		  Onii++;
      }
      GMat[ ii * kSLen + ii] = 1;
  }
  return(1);
}
//////////////////////////////////////////////////////////
//  SetMeXIn
//     Sets up the general X_A matrix with ``double'' inputs
//     Care is made to set it up with usual sign
//     This is used if GFlag = 0;
int LARSObject::SetMeXIn (int kLen, int NLen, int OnL, 
              int* ActiveOns, int *ActiveSigns, double *XInputMat , 
              double *XOutputMat) {
  if (ActiveSigns == NULL) {
	  Rprintf((char*)"LARSObject::SetMeXIn Error, ActiveSigns is NULL\n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1); 
  }
  if (XInputMat == NULL) {
	  Rprintf((char*)"LARSObject::SetMeXIn Error, XInputMat is NULL\n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1); 
  }	 
  if (XOutputMat == NULL) {
	  Rprintf((char*)"LARSObject::SetMeXIn Error, XOutputMat is NULL\n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1); 
  }	     	              
	int ii, jj;
	int WillPostOn;
	int WillPostOff;
    for (ii = 0; ii < OnL; ii++) {
		    WillPostOn  = ActiveOns[ii] * NLen;
		    WillPostOff = ii * NLen;
		    for (jj = 0; jj < NLen; jj++) {
			     XOutputMat[WillPostOff] = ((double)ActiveSigns[ii]) * XInputMat[WillPostOn];
			     WillPostOff++;
			     WillPostOn++;
	        }
    }
	return(1);
}
//////////////////////////////////////////////////////////
//  SetMeXIn
//     Sets up the general X_A matrix with ``long double'' inputs
//     Care is made to set it up with usual sign
//     This is used if GFlag = 0;
int LARSObject::SetMeXIn (int kLen, int NLen, int OnL, int* ActiveOns,  int*ActiveSigns,
              long double *XInputMat , 
              long double *XOutputMat) {
  if (ActiveSigns == NULL) {
	  Rprintf((char*)"LARSObject::SetMeXIn Error, ActiveSigns is NULL\n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1); 
  }
  if (XInputMat == NULL) {
	  Rprintf((char*)"LARSObject::SetMeXIn Error, XInputMat is NULL\n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1); 
  }	 
  if (XOutputMat == NULL) {
	  Rprintf((char*)"LARSObject::SetMeXIn Error, XOutputMat is NULL\n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1); 
  }	
  if (ActiveSigns == NULL) {
	  Rprintf((char*)"LARSObject::SetMeXIn Error, ActiveOns is NULL\n");
	  R_FlushConsole(); R_ProcessEvents(); return(-1); 
  }
	int ii, jj;
	int WillPostOn;
	int WillPostOff;
    for (ii = 0; ii < OnL; ii++) {
		    WillPostOn  = ActiveOns[ii] * NLen;
		    WillPostOff = ii * NLen;
		    for (jj = 0; jj < NLen; jj++) {
			     XOutputMat[WillPostOff] = ActiveSigns[ii] * XInputMat[WillPostOn];
			     WillPostOff++;
			     WillPostOn++;
	        }
    }
	return(1);
}

////////////////////////////////////////////////////////////////////
//  int LARSObject::FeedFreeBetas
//
//    FeedFreeBetas records the results of a lars step in a vector.
//
//
 int LARSObject::FeedFreeBetas(int tt, double FeatConst) {
		       //Rprintf( (char *) "Running FeedFree Beats for FeadConst = %f, tt = %d \n", 
		       //       (double) FeatConst, (int) tt);
		       //R_FlushConsole();
		 if ( (tt+1) *4 +3 >= (floor(this->kLen*DEFCONST)+2) * 4) {
			        Rprintf( (char *) "FeedFreeBetas Not Inputting to Current Penalty, klen = %d \n",
			              this->kLen);
			        return(-1);
		 }
         int jj, ii, cnt1;
         double PumpMe;
         double ResidsAdds1 = 0;
	     double ResidsAdds2 = 0;
	     double ResidsAdds3 = 0;
	       for (jj = 0; jj < this->kLen; jj++) {
	           this->PrevBetas[jj] = (double) this->OnBetas[jj];
           }
           if (GFlag == 0) { for (jj = 0; jj < this->NLen; jj++) {
	                        this->PrevResid[jj] = (double) this->OnResid[jj];
	                        this->PrevMuA[jj] = (double) this->OnMuA[jj];
                       }
           }
	       this->PrevPenalty = (double) this->CurrentPenalty;
	       this->CurrentPenalty = 0;
	       this->PenaltyKeep[ (tt+1) * 4 ] = 0;
	       this->PenaltyKeep[ (tt+1) * 4 + 1] = 0;
	       this->PenaltyKeep[ (tt+1) * 4 + 2] = 0;
	       this->PenaltyKeep[ (tt+1) * 4 +3 ] = 0;
	       if (CurrentGoat >= 0) {
		       this-> CalcWithGoat(FeatConst);   
		       this->CurrentPenalty = this->SumCurResids; 
	       } else {
	           for (jj = 0; jj < this->OnL; jj++) {
	              this->OnBetas[ this->ActiveOns[jj]] = 
	                 (double) this->OnBetas[ this->ActiveOns[jj]] +(double) FeatConst * 
	                 (double) this->ActiveSigns[ jj ] * (double) this->wa[jj];	          
               }		       
		       if (FakeFlag >= 0 && GFlag == 0) {
			       for ( jj = 0; jj < this->NLen; jj++) {
			           this->OnMuA[jj] = (double) this->OnMuA[jj] + 
			                   (double) FeatConst * (double) this->ua[jj];
			           this->OnResid[jj] = (double) this->yys[jj] - (double) this->OnMuA[jj];
			           this->CurrentPenalty += ((double) this->OnResid[jj] ) * 
			                                 ((double) this->OnResid[jj] );
		           }
		           this->SumCurResids = this->CurrentPenalty;
	           } else if (GFlag == 1 || GFlag == 3) {
			          if (this->WDiagFlag == 1) {	            
			             for (jj = 0; jj < this->kLen; jj++) {
				              PumpMe = 0;
				              cnt1 = jj * this->kLen;
				              for (ii = 0; ii < this->OnL; ii++) {
					                 PumpMe += GAll[cnt1+ActiveOns[ii]] * WDiag[ActiveOns[ii]] * wa[ii] * 
					                      ActiveSigns[ii];                        
				              }
				              GAllOnB[jj] += PumpMe * FeatConst;   
		               }
		           } else {
		               for (jj = 0; jj < this->kLen; jj++) {
				              PumpMe = 0;
				              cnt1 = jj * this->kLen;
				              for (ii = 0; ii < this->OnL; ii++) {
					                PumpMe += GAll[cnt1+ActiveOns[ii]] * wa[ii] * 
					                    ActiveSigns[ii];                        
				              }
				              GAllOnB[jj] += PumpMe * FeatConst;   
		              }	              
		         }
	                    // double CurSumYTXAwa;
	                    // double CurSumOnBetaXTXAwa;
	                    // double CurSumwaXATXAwa;
	                    ResidsAdds1 =(double) -2.0 * (double) FeatConst * 
	                                 (double) CurSumYTXAwa;
	                    ResidsAdds2 = (double) 2.0 * (double) FeatConst * 
	                                  (double) CurSumOnBetaXTXAwa;
	                    ResidsAdds3 = FeatConst * FeatConst;
		               this->SumCurResids += ResidsAdds1 + ResidsAdds2 + ResidsAdds3;
			           this->CurrentPenalty = this->SumCurResids;	           
	           } else {
		           Rprintf("LARSObject::FeedFreeBetas, Don't know how we got here\n");
		           R_FlushConsole();
	           }
           }
           this->PenaltyKeep[(tt+1) * 4] = (double) this->CurrentPenalty;
           for (jj = 0; jj < kLen; jj++) {
	           this->CurrentPenalty += this->lambda * fabs( (double) this->OnBetas[jj] );
	           this->PenaltyKeep[(tt+1)*4+1] += this->lambda * fabs( (double) this->OnBetas[jj]);
	           //this->CurrentPenalty[0] += 0.0 * fabs( (double) this->OnBetas[jj] - (double) this->OldBetas[jj] );
	           //this->PenaltyKeep[(tt+1) *4 +2] += 0.0 * fabs( (double) this->OnBetas[jj] - (double) this->OldBetas[jj] );
           }
           for (jj = 0; jj < this->kLen; jj++) {
	            this->BetasKeep[(tt+1) * this->kLen + jj] = (double) this->OnBetas[jj];
           }
           this->JKeep[tt] =this->OnL;
           this->PenaltyKeep[(tt+1) * 4 + 3] = (double) this->CurrentPenalty;
         //  if (PrintFlag > 5) {
	     //      Rprintf((char*)"Penalties: CurrentPenalty = %.4f, PrevPenalty = %.4f\n",
	     //         CurrentPenalty, PrevPenalty);
	     //      R_FlushConsole();
	     //      Rprintf((char*)"SumCurResids = %.4f\n");
	     //      if (GFlag == 1) {		           
	     //        Rprintf((char*) "GFlag = %d, ResidsAdds1 = %.4f, ResidsAdds2 = %.4f, ResidsAdds3 = %.4f\n",
	     //          GFlag, (double) ResidsAdds1, (double) ResidsAdds2, (double) ResidsAdds3);
	     //        Rprintf((char*) "CurSumYTXAwa = %.4f, CurSumOnBetaXTXAwa = %.4f\n", 
	     //           (double) CurSumYTXAwa, (double) CurSumOnBetaXTXAwa);
	     //        R_FlushConsole();
         //      }
         //  }
           return(1);
}


///////////////////////////////////////////////////////////////////
//  int CreateWAll(weightLen, Weights
//
//   Creates a diagonal weight vector (Other weight matrices too complicated to conceive)
//     for the factors.
int LARSObject::CreateWAll(int weightLen, double *Weights) {
   int ii, jj, cnt1;
  if (weightLen != this->kLen) {
	  Rprintf((char*)"CreateWAll, cannot use Weight matrix\n");
	  return(0);
  } 
  
	  if (this->WDiag == NULL) {
		  this->WDiag = (double *)Calloc( this->kLen +1, double);
			  if (this->WDiag == NULL) {
				  Rprintf((char*)"UpdateWDiag, No Memory\n");
				  R_FlushConsole();
				  return(-1);
		      }
    }
    if (this->WDiagBeta == NULL) {
		  this->WDiagBeta = (double *)Calloc( this->kLen +1, double);
			  if (this->WDiagBeta == NULL) {
				  Rprintf((char*)"UpdateWDiag, No Memory for WDiagBeta\n");
				  R_FlushConsole();
				  return(-1);
		      }
    }    
    if (this->WIV == NULL)        {
		  this->WIV = (double *)Calloc( this->kLen +1, double);
			  if (this->WIV == NULL) {
				  Rprintf((char*)"UpdateWDiag, No Memory for WIV\n");
				  R_FlushConsole();
				  return(-1);
		      }    
    }

      //if (this->Wwa == NULL) {
	  //    this->Wwa = (double *)Calloc( this->kLen+1, double);
	  //    if (this->Wwa == NULL) {
	//	      Rprintf((char*)"UpdateWDiag, No Memory for Wwa \n");
	//	      Free(WDiag);
	//	      WDiag = NULL;
	//	      R_FlushConsole();
	//	      return(-1);
	//      }
    //  }
      this->WDiagFlag = 1;
      for (ii = 0; ii < this->kLen; ii++) {
	       this->WDiag[ii] = (double ) Weights[ii];
      }
      AvWeights = 0.0;
      for (ii = 0; ii < this->kLen; ii++) { AvWeights += WDiag[ii]; }
       AvWeights = AvWeights / this->kLen; 

      if (this->GFlag == 1) {
	      int One = 1; double ZeroD = 0.0; double OneD = 0.0;      
	      for (ii = 0; ii < this->kLen; ii++) {
			  WDiagBeta[ii] = WDiag[ii] * OnBetas[ii];
	      }
	      F77_CALL(dsymv)("U", &this->kLen,
		     &OneD, GAll, &this->kLen,
		       this->WDiagBeta, &One, &ZeroD,
		       this->GAllOnB, &One);       
      }
      return(1);
                   
      if (this->GFlag == 1) {
	      for (ii = 0; ii < this->kLen; ii++) {
			  GAllOnB[ii] = 0;
	      }
	      for (ii = 0; ii < this->kLen; ii++) {
		     if (fabs(OnBetas[ii]) > 0) {
			   cnt1 = ii * this->kLen;
		       for (jj= 0; jj < this->kLen; jj++) {
			      GAllOnB[jj] += GAll[cnt1] * WDiag[ii] * OnBetas[ii];
			      cnt1++;
		       }   
	         }
	      } 
      }

}
///////////////////////////////////////////////////////////
//  CreateGAll()
//
//     GAll is the Gram Matrix XTX as stored in Lars Object
//     One may choose not to calculate this unless GFlag = 1;
//     This is calculated best in n > k case.
int LARSObject::CreateGAll() {
  if (this->GAll == NULL) {
	  this->GAll = (double *)Calloc( this->kLen * this->kLen, double);
	  if (this->GAll == NULL) {
		  Rprintf((char*)"UpdateGall, No Memory\n");
		  R_FlushConsole();
		  return(-1);
      }
     this->XTYAll = (double *)Calloc( this->kLen, double);
     if (this->XTYAll == NULL) {
	     Rprintf((char*)"UpdateGall, No memory for XTYAll\n");
	     R_FlushConsole();
	     Free(this->GAll);
	     this->GAll = NULL;
	     return(-1);
     }
  }
  if (GAllOnB == NULL) {
	      this->GAllOnB = (double *)Calloc( this->kLen+1, double);
	      if (this->GAllOnB == NULL)  {
	        Rprintf((char*) "UpdateGAll, No memory for GAllOnB \n");
	        R_FlushConsole();
	        Free(this->GAll); this->GAll = NULL;
	        Free(this->XTYAll); this->XTYAll = NULL;
	        return(-1);
         }
  }
  int cnt1, cntx1, cntx2, ii, jj, kk; 
	  cnt1 = 0;
	  int One = 1; double ZeroD = 0.0; double OneD = 1.0;
    SqMat (GAll, this->NLen, this->kLen, xxs);
      //F77_CALL(dscal)(&this->kLen, &ZeroD, this->GAllOnB, &One);
      //F77_CALL(dgemv)("T", &this->NLen, &this->kLen,
		  //   &OneD, GAll, &this->kLen,
		  //     this->OnBetas, &One, &ZeroD,
		  //     this->GAllOnB, &One);    
      //F77_CALL(dscal)(&this->kLen, &ZeroD, this->GAllOnB, &One);
      F77_CALL(dsymv)("U", &this->kLen,
		     &OneD, GAll, &this->kLen,
		       this->OnBetas, &One, &ZeroD,
		       this->GAllOnB, &One); 
      F77_CALL(dscal)(&this->kLen, &ZeroD,this->XTYAll, &One);
      F77_CALL(dgemv)("T", &this->NLen, &this->kLen,
		     &OneD, xxs, &this->NLen,
		       yys, &One, &ZeroD,
		       XTYAll, &One);                
     // F77_CALL(dscal)(&this->kLen, &ZeroD, this->XTYAll, &One);
     // F77_CALL(dgemv)("T", &this->NLen, &this->kLen,
		 //       &OneD, xxs, &this->NLen,
		 //     this->OnBetas, &One, &ZeroD,
		 //      this->GAllOnB, &One);   		       

      this->SumYYSq  = 0;
      for (ii = 0; ii < this->NLen; ii++) {
	        SumYYSq += (double) this->yys[ii] * (double) this->yys[ii];
      }
      this->GFlag = 1;
      return(1);
      //SumCurResids = SumYYSq;

                 		       
	  for (ii = 0; ii < this->kLen; ii++) {
		  cnt1 = ii * this->kLen + ii;
		  GAll[cnt1] = 0;
		  cntx1 = this->NLen * ii;
		  for (kk = 0; kk < this->NLen; kk++) {
			  GAll[cnt1] += (double) ((double) xxs[cntx1] * (double) xxs[cntx1]);
			  cntx1++;
          }
          GAllOnB[ii] += GAll[cnt1] * OnBetas[ii];
		  if (ii < ((this->kLen)-1)) {
		    for (jj = (ii+1); jj < this->kLen; jj++) {
			  cntx1 = this->NLen * ii; 
			  cntx2 = this->NLen * jj;
			  cnt1 =  ii * this->kLen + jj;
			  this->GAll[cnt1] = 0;
			  for (kk = 0; kk < this->NLen; kk++) {
		          this->GAll[cnt1] += (double) this->xxs[cntx1] * (double) this->xxs[cntx2] ;
		          cntx1++; 
		          cntx2++;
	          }
	          this->GAll[ jj * this->kLen + ii ] = this->GAll[ cnt1];
	           GAllOnB[ii] += GAll[cnt1]  * OnBetas[jj];
	           GAllOnB[jj] += GAll[cnt1]  * OnBetas[ii];
	        }
	       
         }
	}	 
      for (ii = 0; ii < this->kLen; ii++) {
		  GAllOnB[ii] = 0;
      }
      for (ii = 0; ii < this->kLen; ii++) {
	     if (fabs(OnBetas[ii]) > 0) {
		   cnt1 = ii * this->kLen;
	       for (jj= 0; jj < this->kLen; jj++) {
		      GAllOnB[jj] += GAll[cnt1] * OnBetas[ii];  
		      cnt1++;
	       }   
         }
      }
	  cnt1 = 0;
	  for (ii = 0; ii < this->kLen; ii++) { this->XTYAll[ii] = 0; }
	  for (ii = 0; ii < this->kLen; ii++) {
		  for (jj = 0; jj < this->NLen; jj++) {
		      this->XTYAll[ii] += (double) this->xxs[cnt1] * 
		            (double) this->yys[jj];
		      cnt1++;
	      }
      }

  return(1);
}

////////////////////////////////////////////////////////
//  Prints off reduced vector "G" matrix
//
//     See Make GGA
//
//
//
int LARSObject::GGAPrint() {
	   int kOnLen = this->OnL;
	   int ii, jj;
	   if (GGAA == NULL) {
		   Rprintf("LARSOBject::GGAPrint Error, GGAA is NULL\n");
		   R_FlushConsole(); return(-1);
       }
       if (this->PrintFlag > 2) {
	       int PMin = 6;
	       if (kOnLen < 6) {
		       PMin = kOnLen;
	       }
	       Rprintf( (char *) " Printing GGA \n");
	       for (ii =0; ii < PMin; ii++) {
	           Rprintf( (char *) " %d [  ", ii);
	           for (jj = 0; jj < PMin; jj++) {
	   	       Rprintf( (char *) " %.6f ", (double) GGAA[jj * kOnLen + ii]);
	           }
	           Rprintf( (char *) " ]\n");
	      }
	      Rprintf((char*)"\n cbind( c( ");
	       for (ii =0; ii < PMin; ii++) {
	           Rprintf( (char *) " c(  ", ii);
	           for (jj = 0; jj < PMin; jj++) {
	   	       Rprintf( (char *) " %.6f ", (double) GGAA[jj * kOnLen + ii]);
	   	         if (jj < PMin-1) {Rprintf((char*)",");}
	           }
	           Rprintf( (char *) ")");
	           if (ii < PMin-1) {Rprintf((char*)",\n");}
	      }
	      Rprintf((char*)")\n");	      
	      R_FlushConsole();
      }
       if (this->PrintFlag > 3) {
	       int PMin = 6;
	       if (kOnLen < 6) {
		       PMin = kOnLen;
	       }
	       Rprintf( (char *) " Printing IGGA \n");
	       for (ii =0; ii < PMin; ii++) {
	           Rprintf( (char *) " %d [  ", ii);
	           for (jj = 0; jj < PMin; jj++) {
	   	       Rprintf( (char *) " %.3E ", (double) IGGA[jj * kOnLen + ii]);
	           }
	           Rprintf( (char *) " ]\n");
	      }
	      R_FlushConsole();
      }
    return(1);
} 


///////////////////////////////////////////
//  GAllPrint()
//
//  Print the GAll matrix
//
int LARSObject::GAllPrint() {
	   if (GAll == NULL) {
		   Rprintf("LARSObject::GAllPrint Error, GAll is NULL\n");
		   R_FlushConsole(); return(-1);
       }
	   int kOnLen = this->kLen;
	   int ii, jj;
       if (this->PrintFlag > 2) {
	       int PMin = 8;
	       if (kOnLen < 8) {
		       PMin = kLen;
	       }
	       Rprintf( (char *) " Printing GAll \n");
	       for (ii =0; ii < PMin; ii++) {
	           Rprintf( (char *) " %d [  ", ii);
	           for (jj = 0; jj < PMin; jj++) {
	   	       Rprintf( (char *) " %.6f ", (double) this->GAll[jj * kLen + ii]);
	           }
	           Rprintf( (char *) " ]\n");
	      }
	      R_FlushConsole();
      }
    return(1);
} 





///////////////////////////////////////////////////////////////////////////////////////
//  double CalculateGFTAll()
//
//   The GFTAll Algorithm as a stand alone function.
double CalculateGFTAll(double *TFLV2, double *GFT0, double *GFT1, int OnL, int TotK,
          int *ActiveOns, double *OnBetas, double *BetaOlds, int NLen, double *OnResid, 
          double *wa, int *sjA, double *ua, double v1, double v2, double OverallMax) {
	   //Rprintf( (char *) " Calculating GFTAll, NLen = %d, OnL = %d, TotK = %d \n", NLen, OnL, TotK);
	   //R_FlushConsole();
       //Rprintf( (char *) "Come On Print Something Anything \n");
	   //    R_FlushConsole();	   
	   int jj;
	   int RunningMinjj1 = -2, RunningMinjj2 = -2;
	   int OnSign;
	   double MaxMove1 = -1, MaxMove2 = -1;
       TFLV2[0] = (double) 0.0;
       GFT0[0] = (double) 0.0;
       GFT1[0] = (double) 0.0;
       int goflag = 0, OnGo = 0;
       int MAXGO = 100;
       double MinGammas = 0.0;
       double RetAns = 0.0;
       //Rprintf( (char *) "Going Through OnBetas\n");
       //R_FlushConsole();
       //R_ProcessEvents();
       //for (jj = 0; jj < TotK; jj++) {
	   //     Rprintf( (char *) "    OnBetas[%d] = %.4f\n", jj, (double) OnBetas[jj]);
       //}
       //Rprintf( (char *) "Going Through BetaOlds\n");
       //R_FlushConsole();
       //R_ProcessEvents();
       //for (jj = 0; jj < TotK; jj++) {
	   //     Rprintf( (char *) "    BetaOlds[%d] = %.4f\n", jj, (double) BetaOlds[jj]);
       //}       
       while (goflag == 0 && OnGo < MAXGO) {
	       OnGo++;
	       //Rprintf( (char *) "Attempting OnGo == %d, RetAns = %.4f \n", OnGo, (double) RetAns);
	       //R_FlushConsole();
	       //R_ProcessEvents();
       TFLV2[0] = (double) 0.0;
       GFT0[0] = (double) 0.0;
       GFT1[0] = (double) 0.0;	       
       for (jj = 0; jj < OnL; jj++ ) {
	       
	      // Rprintf( (char *) "Looking jj = %d, OnL = %d, ActiveOns[jj] = %d, TotK = %d \n", jj, OnL,
	      //               ActiveOns[jj], TotK );
	       //R_FlushConsole();
	       //R_ProcessEvents();
	       //Rprintf( (char *) "OnBetas[%d] = %.3f, BetaOlds[%d] = %.3f, wa[%d] = %.4f \n", ActiveOns[jj],
	       //       (double) OnBetas[ActiveOns[jj]], ActiveOns[jj], (double) BetaOlds[ActiveOns[jj]],
	       //        jj, (double) wa[jj]);
	       //R_FlushConsole();
	       //R_ProcessEvents();
	       if (ActiveOns[jj] >= TotK) {
		       Rprintf( (char *) "Error ActiveOns[jj] = %d, OnL = %d, TotK = %d \n", ActiveOns[jj], OnL, TotK);
		       R_FlushConsole();
		       R_ProcessEvents();
		       return(-1);
	       }	       	 
	       //Rprintf( (char *) "TFLV2[0] = %.4f, GFT0[0] = %.4f, GFT1[0] = %.4f\n", TFLV2[0], GFT0[0], GFT1[0]);
	       //R_FlushConsole();
	       //R_ProcessEvents();      
	       if ( (double) OnBetas[ActiveOns[jj]] - (double) BetaOlds[ActiveOns[jj]] == 0) {
		        // Rprintf( (char *) "Beta[%d] - BetaOlds[%d] = 0\n", ActiveOns[jj], ActiveOns[jj]);
	            //R_FlushConsole();
	            //R_ProcessEvents();
		       TFLV2[0] += fabs((double) wa[jj] );
	       } else {
		        //Rprintf( (char *) "Beta[%d] - BetaOlds[%d] = %f\n", ActiveOns[jj], ActiveOns[jj],
		        //(double) (  OnBetas[ActiveOns[jj]] - BetaOlds[ActiveOns[jj]] ));
	            //R_FlushConsole();
	            //R_ProcessEvents();
		       OnSign = GetSign( ((double) OnBetas[ActiveOns[jj]] - 
		                         (double) BetaOlds[ActiveOns[jj]]) * sjA[jj] );
		       if (OnSign == -1) {
			       if (RunningMinjj2 < 0) {
				       MaxMove2 = fabs((double) OnBetas[ActiveOns[jj]] - 
				                       (double) BetaOlds[ActiveOns[jj]]) / fabs(wa[jj]);
				       RunningMinjj2 = jj;
			       } else if (MaxMove2  >  fabs((double) OnBetas[ActiveOns[jj]] - 
				                       (double) BetaOlds[ActiveOns[jj]] ) / fabs(wa[jj]) ) {
					   MaxMove2  =  fabs((double) OnBetas[ActiveOns[jj]] - 
				                       (double) BetaOlds[ActiveOns[jj]]) / fabs(wa[jj]);
				       RunningMinjj2 = jj;
			       }
			       if (fabs( (double) OnBetas[ActiveOns[jj]] - 
		                         (double) BetaOlds[ActiveOns[jj]] ) / 
		             fabs( (double) wa[jj] ) > MinGammas) {
			              TFLV2[0] += fabs((double) wa[jj]) * 
			                 GetSign( (double) OnBetas[ActiveOns[jj]] - 
			                 (double) BetaOlds[ActiveOns[jj]] );
		             } else {
		                TFLV2[0] += fabs((double) wa[jj]);
	                 }
		        } else {
		            TFLV2[0] += fabs((double) wa[jj]);
		        }         		        
		   } 
	       if ( (double) OnBetas[ActiveOns[jj]]== 0) {
		        //Rprintf( (char *) "Beta[%d] = 0\n", ActiveOns[jj]);
	            //R_FlushConsole();
	            //R_ProcessEvents();
		       GFT0[0] += fabs((double) wa[jj] );
	       } else {
		       OnSign = GetSign( (double) OnBetas[ActiveOns[jj]]  * (double) sjA[jj]);
		       if (OnSign == -1) {
			       if (RunningMinjj1 < 0) {
				       MaxMove1 = fabs((double) OnBetas[ActiveOns[jj]]) / fabs(wa[jj]);
				       RunningMinjj1 = jj;
			       } else if (MaxMove1  >  fabs((double) OnBetas[ActiveOns[jj]]) / fabs((double)wa[jj]) ) {
					   MaxMove1  =  fabs((double) OnBetas[ActiveOns[jj]])/fabs((double)wa[jj]);
				       RunningMinjj1 = jj;
			       }
                   if (fabs( (double) OnBetas[ActiveOns[jj]])  / 
		             fabs( (double) wa[jj] ) > MinGammas) {
			              GFT0[0] += fabs((double) wa[jj]) * 
			                 OnSign;
		             } else {
		                GFT0[0] += fabs((double) wa[jj]);
	                 }
		        } else {
		            GFT0[0] += fabs((double) wa[jj]);
		        }         		        		         
		    }
		      
       }
		              
       GFT0[0] = - v1 * (double) GFT0[0] / (double) 2.0;
       TFLV2[0] = -(double) TFLV2[0] * (double) v2 / (double) 2.0;
       
       GFT1[0] = 0.0;
     //  if (FakeFlag >= 0 && OnResid != NULL) {
	       if (OnResid == NULL || ua == NULL) {
		       Rprintf((char*)"Lars:CalculateGFTAll error, ua, OnResid are NULL\n");
		       R_FlushConsole(); return(-1);
	       }	       
	       for (jj = 0; jj <  NLen; jj++) {      
		        GFT1[0] += (double) ua[jj] * (double) OnResid[jj];    
	       } 
    //   } else { 
	//       if (XTYAll == NULL || GAllOnB == NULL) {
	//	       Rprintf((char*)"Lars:CalculateGFTAll error, XTYAll, GAllOnB are NULL\n");
	//	       R_FlushConsole(); return(-1);
	//       }
	//       for (jj = 0; jj < OnL; jj++) {
	//	       GFT1[0] += XTYAll[ ActiveOns[jj]] * ActiveSigns[jj] * wa[jj] -
	//	                  GAllOnB[ ActiveOns[jj]] * ActiveSigns[jj] * wa[jj];
	//     } 
    //   }
        RetAns = (double) GFT0[0] + (double) TFLV2[0] + (double) GFT1[0];
        goflag = 1;
        if ((double) OverallMax < (double) RetAns) {
	         goflag = 1;
        } else if ( (MaxMove1 > 0 && MaxMove1 < (double) RetAns) ||
             (MaxMove2 > 0 && MaxMove2 < (double) RetAns) ) {
	             goflag = 0;
	             if ( (MaxMove1 > 0 ) && ((double) MaxMove1 < (double) MaxMove2) ) {
		               MinGammas = MaxMove1;
	             } else  {
		               MinGammas = MaxMove2;
	             }
        }
      }
       //Rprintf( (char *) " Returning GFTAll, NLen = %d, OnL = %d \n", NLen, OnL);
	   //R_FlushConsole();
       return(RetAns);      
}
