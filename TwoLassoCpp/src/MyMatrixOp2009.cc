/******************************************************************************/
/*    MyMatrixOp.cc                                                           */
/*       Code For matrix printing in R                                        */
/*       Uses MyMatrixOp.h                                                    */
/*                                                                            */
/******************************************************************************/
#ifndef MyMatrixOp2009DD
  #include "MyMatrixOp2009.h"
  #define MyMatrixOp2009DD 0
#endif
#ifndef INCIDENTALOP2009DD
  #include "IncidentalOp2009.h"
  #define INCIDENTALOP2009DD 0 
#endif

#ifndef LAPACKDD
  #include <R_ext/Lapack.h>
  #include <R_ext/BLAS.h>
  #define LAPACKDD 0
#endif

//#define  %.3f %.3f
//#define  %.4e %.4e
//#define  %.9f %.8f
//#define  %.9e %.8e

#define PRMATRIX 0

///////////////////////////////////////////////////////////
//   PrintRMatrix
//
//      Prints a reasonably sized, formatted matrix
//
//
//
void PrintVector(long double *Vec, int Len) {
	 const int MaxVec = 10;
	 const int MorePrint = 6;
	 int LLen = Len;
	 if (Len > MaxVec) {
		   LLen = MaxVec;
     }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
	     if (LLen >= MorePrint) {
		     if ( fabs(Vec[ii]) > .001) {
		        Rprintf(" %.3f,", (double) Vec[ii]);
	         } else {
		        Rprintf(" %.4e,", (double) Vec[ii]);
	         }
         } else {
		     if ( fabs(Vec[ii]) > .001) {
		        Rprintf(" %.9f,", (double) Vec[ii]);
	         } else {
		        Rprintf(" %.9e,", (double) Vec[ii]);
	         }	         
         }
     }
     if (Len < MorePrint) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.9f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.9e) \n", (double) Vec[ii]);
         }	 
     } else if (LLen < Len) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f, ... )\n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e, ... )\n", (double) Vec[ii]);
         }
     } else {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e) \n", (double) Vec[ii]);
         }
     }
     R_FlushConsole();
     return;
}

void PrintVectorAll(long double *Vec, int Len) {
	 //const int MaxVec = 10;
	 int LLen = Len;
	 const int MorePrint = 6;
	 //if (Len > MaxVec) {
	//	   LLen = MaxVec;
    // }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
         if (MorePrint > LLen) {
	         Rprintf(" %.9f,", (double) Vec[ii]);
         } else {
	         Rprintf(" %.3f,", (double) Vec[ii]);
         }	         
     }
     if (LLen < Len) {
	     if (MorePrint >= LLen) {
		     if ( fabs(Vec[ii]) > .001) {
		        Rprintf(" %.9f, ... )\n", (double) Vec[ii]);
	         } else {
		        Rprintf(" %.9e, ... )\n", (double) Vec[ii]);
	         }	   
         } else {
			 if ( fabs(Vec[ii]) > .001) {
		        Rprintf(" %.3f, ... )\n", (double) Vec[ii]);
	         } else {
		        Rprintf(" %.4e, ... )\n", (double) Vec[ii]);
	         }	         
	         
         }  
     } else {
	    if (MorePrint >= LLen) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.9f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.9e) \n", (double) Vec[ii]);
         }	 
        } else {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e) \n", (double) Vec[ii]);
         }	 	        
        }    
     }
     R_FlushConsole();
     return;
}
void PrintVector(int *Vec, int Len) {
	 const int MaxVec = 10;
	 int LLen = Len;
	 if (Len > MaxVec) {
		   LLen = MaxVec;
     }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
	     Rprintf(" %d,", (int) Vec[ii]);
     }
     if (LLen < Len) {
	     Rprintf(" %d, ... )\n", (int) Vec[ii]);
     } else {
	     Rprintf(" %d)", (int) Vec[ii]);
     }
     R_FlushConsole();
     return;
}
void PrintVector(short int *Vec, int Len) {
	 const int MaxVec = 10;
	 int LLen = Len;
	 if (Len > MaxVec) {
		   LLen = MaxVec;
     }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
	     Rprintf(" %d,", (int) Vec[ii]);
     }
     if (LLen < Len) {
	     Rprintf(" %d, ... )\n", (int) Vec[ii]);
     } else {
	     Rprintf(" %d)", (int) Vec[ii]);
     }
     R_FlushConsole();
     return;
}

void PrintVectorAll(int *Vec, int Len) {
	 //const int MaxVec = 10;
	 int LLen = Len;
	//if (Len > MaxVec) {
	//	   LLen = MaxVec;
    // }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
	     Rprintf(" %d,", (int) Vec[ii]);
     }
     if (LLen < Len) {
	     Rprintf(" %d, ... )\n", (int) Vec[ii]);
     } else {
	     Rprintf(" %d)", (int) Vec[ii]);
     }
     R_FlushConsole();
     return;
}

void PrintVector(double *Vec, int Len) {
    if (Vec == NULL) {
		Rprintf((char*)"PrintVector:: Error, Vec is NULL\n");
		R_FlushConsole();
		//return(-1);
    }
	 const int MaxVec = 10;
	 const int MorePrint = 6;
	 int LLen = Len;
	 if (Len > MaxVec) {
		   LLen = MaxVec;
     }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
	     if (MorePrint >= LLen) {
		     if ( fabs(Vec[ii]) > .001) {
		        Rprintf(" %.9f,", (double) Vec[ii]);
	         } else {
		        Rprintf(" %.9e,", (double) Vec[ii]);
	         }
         } else {
		     if ( fabs(Vec[ii]) > .001) {
		        Rprintf(" %.3f,", (double) Vec[ii]);
	         } else {
		        Rprintf(" %.9f,", (double) Vec[ii]);
	         }	         
         }
     }
     if (LLen < Len) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f, ... )\n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e, ... )\n", (double) Vec[ii]);
         }	     
     } else {
	    if (MorePrint >= LLen) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.9f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.9e) \n", (double) Vec[ii]);
         }
        } else {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e) \n", (double) Vec[ii]);
         }	        
        }	     
     }
     R_FlushConsole();
     return;
}
 

void PrintVectorTimes(double *Vec, int Len, double Factor) {
    if (Vec == NULL) {
		Rprintf((char*)"PrintVector:: Error, Vec is NULL\n");
		R_FlushConsole();
		//return(-1);
    }
	 const int MaxVec = 10;
	 const int MorePrint = 6;
	 int LLen = Len;
	 if (Len > MaxVec) {
		   LLen = MaxVec;
     }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
	     if (MorePrint >= LLen) {
		     if ( fabs(Vec[ii]) > .001) {
		        Rprintf(" %.9f,", (double) Vec[ii]*Factor);
	         } else {
		        Rprintf(" %.9e,", (double) Vec[ii]*Factor);
	         }
         } else {
		     if ( fabs(Vec[ii]) > .001) {
		        Rprintf(" %.3f,", (double) Vec[ii]*Factor);
	         } else {
		        Rprintf(" %.9f,", (double) Vec[ii]*Factor);
	         }	         
         }
     }
     if (LLen < Len) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f, ... )\n", (double) Vec[ii]*Factor);
         } else {
	        Rprintf(" %.4e, ...) \n", (double) Vec[ii]*Factor);
         }	     
     } else {
	    if (MorePrint >= LLen) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.9f) \n", (double) Vec[ii]*Factor);
         } else {
	        Rprintf(" %.9e) \n", (double) Vec[ii]*Factor);
         }
        } else {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f) \n", (double) Vec[ii]*Factor);
         } else {
	        Rprintf(" %.4e) \n", (double) Vec[ii]*Factor);
         }	        
        }	     
     }
     R_FlushConsole();
     return;
}
 

void PrintVectorElim(double *Vec, int Len, int Elim) {
    if (Vec == NULL) {
		Rprintf((char*)"PrintVector:: Error, Vec is NULL\n");
		R_FlushConsole();
		//return(-1);
    }
	 const int MaxVec = 10;
	 const int MorePrint = 6;
	 int LLen = Len;
	 if (Len > MaxVec) {
		   LLen = MaxVec;
     }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
	     if (ii != Elim) {
		     if (MorePrint >= LLen) {
			     if ( fabs(Vec[ii]) > .001) {
			        Rprintf(" %.9e", (double) Vec[ii]);
		         } else {
			        Rprintf(" %.9e", (double) Vec[ii]);
		         }
	         } else {
			     if ( fabs(Vec[ii]) > .001) {
			        Rprintf(" %.3f", (double) Vec[ii]);
		         } else {
			        Rprintf(" %.9f", (double) Vec[ii]);
		         }	         
	         }
	         if (ii == LLen-2 && LLen < Len && LLen-1 == Elim) {
		         Rprintf(", ... )\n");
	         } else if (ii == LLen-2 && LLen-1 == Elim) {
		         Rprintf(")");
	         } else {
		         Rprintf(",");
	         }
         }
     }
     if (LLen < Len) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f, ... \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e, ... \n", (double) Vec[ii]);
         }	     
     } else {
	    if (MorePrint >= LLen) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.9f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.9e) \n", (double) Vec[ii]);
         }
        } else {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e) \n", (double) Vec[ii]);
         }	        
        }	     
     }
     R_FlushConsole();
     return;
}
 

void PrintVectorAll(double *Vec, int Len) {
    if (Vec == NULL) {
		Rprintf((char*)"PrintVector:: Error, Vec is NULL\n");
		R_FlushConsole();
		//return(-1);
    }
	// const int MaxVec = 10;
    const int Moreprint = 5;
	 int LLen = Len;
	 //if (Len > MaxVec) {
	//	   LLen = MaxVec;
    // }
     Rprintf("c( ");
     int ii; 
     for (ii = 0; ii < LLen-1; ii++) {
	   if (LLen <= Moreprint) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.9f,", (double) Vec[ii]);
         } else {
	        Rprintf(" %.9e,", (double) Vec[ii]);
         }
      } else {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f,", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e,", (double) Vec[ii]);
         }	      
      }
     }
     if (LLen >= Moreprint) {
         if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.9f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.9e) \n", (double) Vec[ii]);
         }		     
     } else if (LLen < Len) {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f, ... \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e, ... \n", (double) Vec[ii]);
         }	     
     } else {
	     if ( fabs(Vec[ii]) > .001) {
	        Rprintf(" %.3f) \n", (double) Vec[ii]);
         } else {
	        Rprintf(" %.4e) \n", (double) Vec[ii]);
         }	     
     }
     R_FlushConsole();
     return;
}
     
 
void PrintRMatrix(long double *Mat, int NR, int NC) {
	if (Mat == NULL) {
		Rprintf((char*)"PrintRMatrix:: Severe Error, Mat is NULL\n");
		R_FlushConsole();
		//return(-1);
    }
	const int maxNR = 8;
	const int maxNC = 8;
	const int MorePrintC = 5;
	int NRR, NCC;
	if (NR > maxNR) {
		 NRR = maxNR;
    } else {
	     NRR = NR;
    }
    if (NC > maxNC) {
		 NCC = maxNC;
    } else {
	     NCC = NC;
    } 
	int ii, jj, ptot;
	if (PRMATRIX == 1) {
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "  %d   ", jj);
    }
    Rprintf( (char *) "\n");
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "-------", jj);
    }
    Rprintf( (char *) "\n");
    R_FlushConsole();    
     
	for (ii = 0; ii < NRR; ii++) {
		ptot = ii;
		Rprintf( (char *) "Row %d:[ ", ii);
		for (jj = 0; jj < NCC;jj++) {
			if ( fabs((double) Mat[ptot]) > .001) {
			  Rprintf( (char *) " %.3f", (double) Mat[ptot]);
		    } else {
			  Rprintf( (char *) " %.4e", (double) Mat[ptot]);
		    }
			R_FlushConsole();
			ptot += NR;
        }
        Rprintf( (char *) "\n");
        R_FlushConsole();
        R_ProcessEvents();
    }
    Rprintf( (char *) "\n");
    }
    R_FlushConsole();
    if ( NC <= NCC && NR <= NRR) {
	   // Rprintf((char*)"MatE <- rbind( \n");
	   Rprintf((char*)" rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
			    if (NCC <= MorePrintC) {
					if ( fabs((double) Mat[jj * NR + ii]) > .001) {
					  Rprintf( (char *) " %.9f, ", (double) Mat[jj * NR + ii]);
				    } else {
					  Rprintf( (char *) " %.9e, ", (double) Mat[jj * NR + ii]);
				    }
		        } else {
					if ( fabs((double) Mat[jj * NR + ii]) > .001) {
					  Rprintf( (char *) " %.3f, ", (double) Mat[jj * NR + ii]);
				    } else {
					  Rprintf( (char *) " %.4e, ", (double) Mat[jj * NR + ii]);
				    }			        
		        }
		    }
		    if (NCC <= MorePrintC) {
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.9f )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.9e )", (double) Mat[(NCC-1)*NR + ii]);
			    }
	        } else {
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.3f )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.4e )", (double) Mat[(NCC-1)*NR + ii]);
			    }		     
	        }
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
		        Rprintf( " ) \n");        
	        } 
	        Rprintf("\n");
	        R_FlushConsole();
       }
   } 
    return;
}
void PrintRMatrix(double *Mat, int NR, int NC) {
	const int maxNR = 8;
	const int maxNC = 8;
	const int MorePrintC = 5;
	int NRR, NCC;
	if (NR > maxNR) {
		 NRR = maxNR;
    } else {
	     NRR = NR;
    }
    if (NC > maxNC) {
		 NCC = maxNC;
    } else {
	     NCC = NC;
    } 
	int ii, jj, ptot;
	if (PRMATRIX == 1) {
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "  %d   ", jj);
    }
    Rprintf( (char *) "\n");
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "-------", jj);
    }
    Rprintf( (char *) "\n");
    R_FlushConsole();    
     
	for (ii = 0; ii < NRR; ii++) {
		ptot = ii;
		Rprintf( (char *) "Row %d:[ ", ii);
		for (jj = 0; jj < NCC;jj++) {
			if ( fabs((double) Mat[ptot]) > .001) {
			  Rprintf( (char *) " %.3f", (double) Mat[ptot]);
		    } else {
			  Rprintf( (char *) " %.4e", (double) Mat[ptot]);
		    }
			R_FlushConsole();
			ptot += NR;
        }
        Rprintf( (char *) "\n");
        R_FlushConsole();
        R_ProcessEvents();
    }
    Rprintf( (char *) "\n");
    }
    //R_FlushConsole();
    //Rprintf( (char *) "\n");
    R_FlushConsole();
    if ( NC <= NCC && NR <= NRR) {
	    //Rprintf((char*)"  <- rbind( \n");
	    Rprintf((char*)"  rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
				if (NCC <= MorePrintC) {
					if ( fabs((double) Mat[jj * NR + ii]) > .001) {
					  Rprintf( (char *) " %.9f,", (double) Mat[jj * NR + ii]);
				    } else {
					  Rprintf( (char *) " %.9e,", (double) Mat[jj * NR + ii]);
				    }			    
	            } else {
					if ( fabs((double) Mat[jj * NR + ii]) > .001) {
					  Rprintf( (char *) " %.3f,", (double) Mat[jj * NR + ii]);
				    } else {
					  Rprintf( (char *) " %.4e,", (double) Mat[jj * NR + ii]);
				    }		            
	            }
		    }
            if (NCC <= MorePrintC) {		    
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.9f )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.9e )", (double) Mat[(NCC-1) * NR + ii]);
			    }	
		    } else {
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.3f )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.4e)", (double) Mat[(NCC-1) * NR + ii]);
			    }				   
		    } 	    
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
		        Rprintf( " ) \n");
	        } 
	        Rprintf("\n");
	        R_FlushConsole();
       }
   }  else {
	    //Rprintf((char*)"  <- rbind( \n");
	    Rprintf((char*)"  rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
				if (NCC <= MorePrintC) {
					if ( fabs((double) Mat[jj * NR + ii]) > .001) {
					  Rprintf( (char *) " %.9f,", (double) Mat[jj * NR + ii]);
				    } else {
					  Rprintf( (char *) " %.9e,", (double) Mat[jj * NR + ii]);
				    }			    
	            } else {
					if ( fabs((double) Mat[jj * NR + ii]) > .001) {
					  Rprintf( (char *) " %.3f,", (double) Mat[jj * NR + ii]);
				    } else {
					  Rprintf( (char *) " %.4e,", (double) Mat[jj * NR + ii]);
				    }		            
	            }
		    }
            if (NCC <= MorePrintC) {		    
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.9f )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.9e )", (double) Mat[(NCC-1) * NR + ii]);
			    }	
		    } else {
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.3f )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.4e)", (double) Mat[(NCC-1) * NR + ii]);
			    }				   
		    } 	    
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
		        Rprintf( " ) \n");
	        } 
	        Rprintf("\n");
	        R_FlushConsole();
       }
   }     
    return;
}          


void PrintRpMatrix(double **pMat, int NR, int NC) {
	const int maxNR = 8;
	const int maxNC = 8;
	const int MorePrintC = 5;
	int NRR, NCC;
	if (NR > maxNR) {
		 NRR = maxNR;
    } else {
	     NRR = NR;
    }
    if (NC > maxNC) {
		 NCC = maxNC;
    } else {
	     NCC = NC;
    } 
	int ii, jj, ptot;
	if (PRMATRIX == 1) {
	Rprintf( (char *) "        ");  R_FlushConsole();
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "  %d   ", jj); R_FlushConsole();
  }
  Rprintf( (char *) "\n");
	Rprintf( (char *) "        ");  R_FlushConsole();
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "-------", jj);
    }
    Rprintf( (char *) "\n");
    R_FlushConsole();    
     
	for (ii = 0; ii < NRR; ii++) {
		ptot = ii;
		Rprintf( (char *) "Row %d:[ ", ii);
		for (jj = 0; jj < NCC;jj++) {
			if ( fabs((double) *(pMat[jj]+ii)) > .001) {
			  Rprintf( (char *) " %.3f", (double) *(pMat[jj]+ii));
		    } else {
			  Rprintf( (char *) " %.4e", (double) *(pMat[jj]+ii));
		    }
			R_FlushConsole();
			ptot += NR;
        }
        Rprintf( (char *) "\n");
        R_FlushConsole();
        R_ProcessEvents();
    }
    Rprintf( (char *) "\n");  R_FlushConsole();
    }
    //R_FlushConsole();
    //Rprintf( (char *) "\n");
    R_FlushConsole();
    if ( NC <= NCC && NR <= NRR) {
	    //Rprintf((char*)"  <- rbind( \n");
	    Rprintf((char*)"  rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
				if (NCC <= MorePrintC) {
					if ( fabs((double) *(pMat[jj] + ii)) > .001) {
					  Rprintf( (char *) " %.9f,", (double) *(pMat[jj]+ii));
				  } else {
					  Rprintf( (char *) " %.9e,", (double) *(pMat[jj]+ii));
				  }			    
	      } else {
					if ( fabs((double) *(pMat[jj]+ ii)) > .001) {
					  Rprintf( (char *) " %.3f,", (double) *(pMat[jj] + ii));
				  } else {
					  Rprintf( (char *) " %.4e,", (double) * (pMat[jj] + ii));
				  }		            
	      }
		    }
        if (NCC <= MorePrintC) {		    
				if ( fabs((double) *(pMat[(NCC-1)] + ii)) > .001) {
				  Rprintf( (char *) " %.9f )", (double) *(pMat[NCC-1] + ii));
			    } else {
				  Rprintf( (char *) " %.9e )", (double) *(pMat[(NCC-1)] + ii));
			    }	
		    } else {
				if ( fabs((double) *(pMat[(NCC-1)] + ii)) > .001) {
				  Rprintf( (char *) " %.3f )", (double) *(pMat[(NCC-1)]+ ii));
			    } else {
				  Rprintf( (char *) " %.4e)", (double) *(pMat[(NCC-1)] + ii));
			    }				   
		    } 	    
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
		        Rprintf( " ) \n");
	        } 
	        Rprintf("\n");
	        R_FlushConsole();
       }
   }  else {
	    //Rprintf((char*)"  <- rbind( \n");
	    Rprintf((char*)"  rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
				if (NCC <= MorePrintC) {
					if ( fabs((double) *(pMat[jj] + ii)) > .001) {
					  Rprintf( (char *) " %.9f,", (double) *(pMat[jj] + ii));
				    } else {
					  Rprintf( (char *) " %.9e,", (double) *(pMat[jj] + ii));
				    }			    
	            } else {
					if ( fabs((double) *(pMat[jj] + ii)) > .001) {
					  Rprintf( (char *) " %.3f,", (double) *(pMat[jj] + ii));
				    } else {
					  Rprintf( (char *) " %.4e,", (double) *(pMat[jj] + ii));
				    }		            
	            }
		    }
            if (NCC <= MorePrintC) {		    
				if ( fabs((double) *(pMat[(NCC-1)] + ii)) > .001) {
				  Rprintf( (char *) " %.9f )", (double) *(pMat[(NCC-1)] + ii));
			    } else {
				  Rprintf( (char *) " %.9e )", (double) *(pMat[(NCC-1)] + ii));
			    }	
		    } else {
				if ( fabs((double) *(pMat[(NCC-1)] + ii)) > .001) {
				  Rprintf( (char *) " %.3f )", (double) *(pMat[(NCC-1)] + ii));
			    } else {
				  Rprintf( (char *) " %.4e)", (double) *(pMat[(NCC-1)]+ ii));
			    }				   
		    } 	    
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
		        Rprintf( " ) \n");
	        } 
	        Rprintf("\n");
	        R_FlushConsole();
       }
   }     
    return;
}


void PrintRMatrixDenseAll(double *Mat, int NR, int NC) {
	const int maxNR = 8;
	const int maxNC = 8;
	const int MorePrintC = 5;
	int NRR, NCC;
	if (NR > maxNR) {
		 NRR = maxNR;
    } else {
	     NRR = NR;
    }
    if (NC > maxNC) {
		 NCC = maxNC;
    } else {
	     NCC = NC;
    } 
	int ii, jj, ptot;
	if (PRMATRIX == 1) {
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "  %d   ", jj);
    }
    Rprintf( (char *) "\n");
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "-------", jj);
    }
    Rprintf( (char *) "\n");
    R_FlushConsole();    
     
	for (ii = 0; ii < NRR; ii++) {
		ptot = ii;
		Rprintf( (char *) "Row %d:[ ", ii);
		for (jj = 0; jj < NCC;jj++) {
			if ( fabs((double) Mat[ptot]) > .001) {
			  Rprintf( (char *) " %.8e", (double) Mat[ptot]);
		    } else {
			  Rprintf( (char *) " %.8e", (double) Mat[ptot]);
		    }
			R_FlushConsole();
			ptot += NR;
        }
        Rprintf( (char *) "\n");
        R_FlushConsole();
        R_ProcessEvents();
    }
    Rprintf( (char *) "\n");
    }
    R_FlushConsole();
    //Rprintf( (char *) "\n");
    R_FlushConsole();
    if ( NC <= NCC && NR <= NRR) {
	   // Rprintf((char*)" <- rbind( \n");
	   Rprintf((char*)" rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
				if (NCC <= MorePrintC) {
					if ( fabs((double) Mat[jj * NR + ii]) > .001) {
					  Rprintf( (char *) " %.9e,", (double) Mat[jj * NR + ii]);
				    } else {
					  Rprintf( (char *) " %.9e,", (double) Mat[jj * NR + ii]);
				    }			    
	            } else {
					if ( fabs((double) Mat[jj * NR + ii]) > .001) {
					  Rprintf( (char *) " %.8e,", (double) Mat[jj * NR + ii]);
				    } else {
					  Rprintf( (char *) " %.8e,", (double) Mat[jj * NR + ii]);
				    }		            
	            }
		    }
            if (NCC <= MorePrintC) {		    
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.9e )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.9e )", (double) Mat[(NCC-1) * NR + ii]);
			    }	
		    } else {
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.8e )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.8e)", (double) Mat[(NCC-1) * NR + ii]);
			    }				   
		    } 	    
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
		        Rprintf( " ) \n");
	        } 
	        Rprintf("\n");
	        R_FlushConsole();
       }
   }     
    return;
}


void PrintRMatrixElim(double *Mat, int NR, int NC, int ElimR, int ElimC) {
	const int maxNR = 8;
	const int maxNC = 8;
	const int MorePrintC = 5;
	int NRR, NCC;
	if (NR > maxNR) {
		 NRR = maxNR;
    } else {
	     NRR = NR;
    }
    if (NC > maxNC) {
		 NCC = maxNC;
    } else {
	     NCC = NC;
    } 
	int ii, jj, ptot;
	if (PRMATRIX==1) {
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		if (jj != ElimC) {
		  Rprintf( (char *) "  %d   ", jj);
	    }
    }
    Rprintf( (char *) "\n");
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "-------", jj);
    }
    Rprintf( (char *) "\n");
    R_FlushConsole();    
     
	for (ii = 0; ii < NRR; ii++) {
		ptot = ii;
		if (ii != ElimR) {
			Rprintf( (char *) "Row %d:[ ", ii);
			for (jj = 0; jj < NCC;jj++) {
				if (jj != ElimC) {
					if ( fabs((double) Mat[ptot]) > .001) {
					  Rprintf( (char *) " %.3f", (double) Mat[ptot]);
				    } else {
					  Rprintf( (char *) " %.4e", (double) Mat[ptot]);
				    }
					R_FlushConsole();
			    }
				ptot += NR;
	        }
	        Rprintf( (char *) "\n");
	        R_FlushConsole();
	        R_ProcessEvents();
        }
    }
    Rprintf( (char *) "\n");
    }
    //R_FlushConsole();
    //Rprintf( (char *) "\n");
    R_FlushConsole();
    if ( NC <= NCC && NR <= NRR) {
	    //Rprintf((char*)"  <- rbind( \n");
	    Rprintf((char*)"  rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    if (ii != ElimR) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
			    if (ElimC != jj) {
					if (NCC <= MorePrintC) {
						if ( fabs((double) Mat[jj * NR + ii]) > .001) {
						  Rprintf( (char *) " %.9f", (double) Mat[jj * NR + ii]);
					    } else {
						  Rprintf( (char *) " %.9e", (double) Mat[jj * NR + ii]);
					    }			    
		            } else {
						if ( fabs((double) Mat[jj * NR + ii]) > .001) {
						  Rprintf( (char *) " %.3f", (double) Mat[jj * NR + ii]);
					    } else {
						  Rprintf( (char *) " %.4e", (double) Mat[jj * NR + ii]);
					    }		            
		            }
		            if (NCC-1 == ElimC && jj == NCC-2) {
			            Rprintf(")");
		            } else {Rprintf((char*) ",");}
	             }
		    }
		    if (ElimC != NCC-1) {
            if (NCC <= MorePrintC) {		    
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.9f )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.9e )", (double) Mat[(NCC-1) * NR + ii]);
			    }	
		    } else {
				if ( fabs((double) Mat[(NCC-1) * NR + ii]) > .001) {
				  Rprintf( (char *) " %.3f )", (double) Mat[(NCC-1) * NR + ii]);
			    } else {
				  Rprintf( (char *) " %.4e)", (double) Mat[(NCC-1) * NR + ii]);
			    }				   
		    } 	
	        }    
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
		        Rprintf( " ) \n");
	        } 
	        Rprintf("\n");
            }
	        R_FlushConsole();
       }
   }     
    return;
}

void PrintRMatrix(int *Mat, int NR, int NC) {
	const int maxNR = 8;
	const int maxNC = 8;
	int NRR, NCC;
	if (NR > maxNR) {
		 NRR = maxNR;
    } else {
	     NRR = NR;
    }
    if (NC > maxNC) {
		 NCC = maxNC;
    } else {
	     NCC = NC;
    } 
	int ii, jj, ptot;
	if (PRMATRIX == 1) {
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "  %d   ", jj);
    }
    Rprintf( (char *) "\n");
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "-------", jj);
    }
    Rprintf( (char *) "\n");
    R_FlushConsole();    
     
	for (ii = 0; ii < NRR; ii++) {
		ptot = ii;
		Rprintf( (char *) "Row %d:[ ", ii);
		for (jj = 0; jj < NCC;jj++) {
			Rprintf( (char *) " %d", (int) Mat[ptot]);
			R_FlushConsole();
			ptot += NR;
        }
        Rprintf( (char *) "\n");
        R_FlushConsole();
        R_ProcessEvents();
    }
    Rprintf( (char *) "\n");
    }
    //R_FlushConsole();
    //Rprintf( (char *) "\n");
    R_FlushConsole();
    if ( NC <= NCC && NR <= NRR) {
	    //Rprintf((char*)" <- rbind( \n");
	    Rprintf((char*)"  rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
			    Rprintf(" %d, ", (int) Mat[jj * NR + ii ]);
		    }
		    Rprintf(" %d )", (int) Mat[(NCC-1) * NR + ii]);
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
               Rprintf("\n");
	        } 
	        Rprintf("\n");
	        R_FlushConsole();
       }
   }  
   Rprintf(") \n");   
    return;
}



void PrintRMatrix(short int *Mat, int NR, int NC) {
	const int maxNR = 8;
	const int maxNC = 8;
	int NRR, NCC;
	if (NR > maxNR) {
		 NRR = maxNR;
    } else {
	     NRR = NR;
    }
    if (NC > maxNC) {
		 NCC = maxNC;
    } else {
	     NCC = NC;
    } 
	int ii, jj, ptot;
	if (PRMATRIX == 1) {
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "  %d   ", jj);
    }
    Rprintf( (char *) "\n");
	Rprintf( (char *) "        ");
	for (jj = 0; jj < NCC; jj++) {
		Rprintf( (char *) "-------", jj);
    }
    Rprintf( (char *) "\n");
    R_FlushConsole();    
     
	for (ii = 0; ii < NRR; ii++) {
		ptot = ii;
		Rprintf( (char *) "Row %d:[ ", ii);
		for (jj = 0; jj < NCC;jj++) {
			Rprintf( (char *) " %d", (int) Mat[ptot]);
			R_FlushConsole();
			ptot += NR;
        }
        Rprintf( (char *) "\n");
        R_FlushConsole();
        R_ProcessEvents();
    }
    Rprintf( (char *) "\n");
    }
    //R_FlushConsole();
    //Rprintf( (char *) "\n");
    R_FlushConsole();
    if ( NC <= NCC && NR <= NRR) {
	    //Rprintf((char*)"  rbind( \n");
	    Rprintf((char*)"  rbind( \n");
	    for (ii = 0; ii < NRR; ii++) {
		    Rprintf("c( ");
		    for (jj = 0; jj < NCC-1; jj++) {
			    Rprintf(" %d, ", (int) Mat[jj * NR + ii ]);
		    }
		    Rprintf(" %d )", (int) Mat[(NCC-1) * NR + ii]);
			if (ii < NRR-1) {
				Rprintf(",");
	        } else {
               Rprintf("\n");
	        } 
	        Rprintf("\n");
	        R_FlushConsole();
       }
   }    
   Rprintf(") \n"); 
    return;
}


   int SignZero(double Try, double CZ) {
	if (Try <= CZ && Try >= -CZ ) {
		 return(0);
    } else if (Try > CZ) {
	    return(1);
    } else {
	    return(-1);
    }
}
////////////////////////////////////////////////////////////////////////////////
//  Self Implemented Matrix Multiplication
//
//    You probably don't want to use any of these functions but they are here
//     for backwards compatibility, testing things out in case Lapack Access
//     fails.s
//
////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
//  InvertMatrix:  Matrix Inversion Algorithm
//
//     Quick solution to 1 and 2 dimensional matrix cases
//     Row reduction to do all inversions.  
//
//
 int InvertMatrix(int MatLeng, long double *InputMatrix, long double *OutputMatrix) {
	 if (InputMatrix == NULL) {
		Rprintf((char*)"InvertMatrix:: Severe Error, InputMatrix is NULL\n");
		R_FlushConsole();
		return(-1);
    }
	 if (OutputMatrix == NULL) {
		Rprintf((char*)"InvertMatrix:: Severe Error, OutputMatrix is NULL\n");
		R_FlushConsole();
		return(-1);
    }    
  int ii,jj,kk, cno = 0;
  //int i1,j1;
  int KiiS, Kiijj;
  double OnDivisor;
  //double COT;
  if (MatLeng == 1) {
	     if (InputMatrix[0] == 0) {
		     Rprintf("InvertMatrix Error: Trying to invert 1/ 0\n");
		     R_FlushConsole(); R_ProcessEvents(); return(-1);
	     }
	     OutputMatrix[0] = ((long double)1.00) / ((long double)InputMatrix[0]); return(1);
  }
  if (MatLeng == 2) {
	  OnDivisor = InputMatrix[0] * InputMatrix[3] - InputMatrix[1] * InputMatrix[2];
	  if (OnDivisor == 0) { 
		  Rprintf("InvertMatrix Error: OnDivisor == 0 \n");
		  R_FlushConsole(); R_ProcessEvents(); return(-1);
      }
      OnDivisor = ((long double)1.0) / ((long double)OnDivisor);
      OutputMatrix[0] = InputMatrix[3] * OnDivisor; OutputMatrix[3] = InputMatrix[0] * OnDivisor;
      OutputMatrix[1] = -InputMatrix[1] * OnDivisor;
      OutputMatrix[2] = -InputMatrix[2] * OnDivisor;
      return(1);
  }
  for (ii = 0; ii < MatLeng; ii++) {
    for (jj = 0; jj <MatLeng; jj++) {
       OutputMatrix[cno] = 0;
       cno++;
    }
  }
  cno = 0;
  OutputMatrix[0] = 1;
  for (ii = 0; ii < MatLeng; ii++) {
     OutputMatrix[cno] =1;
     cno += MatLeng + 1;
  }
  /*  if (MatLeng < 8) { 
     Rprintf( (char *) " Start Inverse Matrix Solution:  \n");
     for (ii = 0; ii < MatLeng; ii++) {
	     Rprintf( (char *) " [  ");
	     for (jj = 0; jj < MatLeng; jj++) {
		      Rprintf( (char *) "  %.3f ", (double) InputMatrix[ MatLeng * jj + ii ] );
	     }
	     Rprintf( (char *) " ] \n");
     }
     R_FlushConsole();
  } */
  KiiS = -MatLeng;
  for (ii = 0; ii < MatLeng; ii++) {
     KiiS += MatLeng;
     if ( InputMatrix[ KiiS + ii ] == 0.0000) {
       Rprintf( (char *) "Invert Matrix Flaw in Invert Matrix, Zero Determinant, Dim = %d \n", 
              (int) MatLeng);
       R_FlushConsole();
       R_ProcessEvents();
       return(-1);
     } else {
         OnDivisor = InputMatrix[ KiiS + ii ];
     }
     Kiijj = ii;
     for ( jj = 0; jj < MatLeng; jj++) {
          if (jj >= ii) {
	           InputMatrix[ Kiijj ]  =  InputMatrix[Kiijj] / OnDivisor;
	      }
          OutputMatrix[ Kiijj]  = OutputMatrix[Kiijj] / OnDivisor;
          Kiijj += MatLeng;
     }
        for (jj  = 0; jj < MatLeng; jj++) {
             if ( InputMatrix[ KiiS + jj] == 0) {
             } else if (jj == ii) {
             } else {
                OnDivisor = InputMatrix[ KiiS+jj] / InputMatrix [KiiS+ ii];
                Kiijj = 0;
                for (kk = 0; kk < MatLeng; kk++) {
	                if (kk >= ii) {
                       InputMatrix[ Kiijj+jj] = InputMatrix[Kiijj+jj] - InputMatrix[Kiijj+ii] * OnDivisor;
                    }
                    OutputMatrix[ Kiijj+jj] = OutputMatrix[Kiijj+jj] - OutputMatrix[Kiijj + ii] * OnDivisor;
                    Kiijj += MatLeng;                                    
                }
             }     
        
        }
        /*
  if (MatLeng < 8) { 
     Rprintf( (char *) " Iteration %d Inverse Matrix Process:  \n", ii);
     for (i1 = 0; i1 < MatLeng; i1++) {
	     Rprintf( (char *) " [  ");
	     for (j1 = 0; j1 < MatLeng; j1++) {
		      Rprintf( (char *) "  %.4f ", (double) OutputMatrix[ MatLeng * j1 + i1 ] );
	     }
	     Rprintf( (char *) " ] \n");
     }
     R_FlushConsole();
  }
   if (MatLeng < 8) { 
     Rprintf( (char *) " Iteration %d Source Matrix :  \n", ii);
     for (i1 = 0; i1 < MatLeng; i1++) {
	     Rprintf( (char *) " [  ");
	     for (j1 = 0; j1 < MatLeng; j1++) {
		      Rprintf( (char *) "  %.3f ", (double) InputMatrix[ MatLeng * j1 + i1 ] );
	     }
	     Rprintf( (char *) " ] \n");
     }
     R_FlushConsole();
  } */
  }
  //if (MatLeng < 8) { 
  //   Rprintf( (char *) " Inverse Matrix Solution:  \n");
  //   for (ii = 0; ii < MatLeng; ii++) {
  //     Rprintf( (char *) " [  ");
  //	     for (jj = 0; jj < MatLeng; jj++) {
  //		      Rprintf( (char *) "  %.6f ", (double) OutputMatrix[ MatLeng * jj + ii ] );
  //	     }
  //	     Rprintf( (char *) " ] \n");
  //   }
  //   R_FlushConsole();
  //} 
  return(1);
}

 int InvertMatrix(int MatLeng, double *InputMatrix, double *OutputMatrix) {
	 if (InputMatrix == NULL) {
		Rprintf((char*)"InvertMatrix:: Severe Error, InputMatrix is NULL\n");
		R_FlushConsole();
		return(-1);
    }
	 if (OutputMatrix == NULL) {
		Rprintf((char*)"InvertMatrix:: Severe Error, OutputMatrix is NULL\n");
		R_FlushConsole();
		return(-1);
    } 	 
  int ii,jj,kk, cno = 0;
  //int i1,j1;
  int KiiS, Kiijj;
  double OnDivisor;
  //double COT;
  if (MatLeng == 1) {
	     if (InputMatrix[0] == 0) {
		     Rprintf("InvertMatrix Error: Trying to invert 1/ 0\n");
		     R_FlushConsole(); R_ProcessEvents(); return(-1);
	     }
	     OutputMatrix[0] = ((long double)1.00) / ((long double)InputMatrix[0]); return(1);
  }
  if (MatLeng == 2) {
	  OnDivisor = InputMatrix[0] * InputMatrix[3] - InputMatrix[1] * InputMatrix[2];
	  if (OnDivisor == 0) { 
		  Rprintf("InvertMatrix Error: OnDivisor == 0 \n");
		  R_FlushConsole(); R_ProcessEvents(); return(-1);
      }
      OnDivisor = ((long double)1.0) / ((long double)OnDivisor);
      OutputMatrix[0] = InputMatrix[3] * OnDivisor; OutputMatrix[3] = InputMatrix[0] * OnDivisor;
      OutputMatrix[1] = -InputMatrix[1] * OnDivisor;
      OutputMatrix[2] = -InputMatrix[2] * OnDivisor;
      return(1);
  }
  for (ii = 0; ii < MatLeng; ii++) {
    for (jj = 0; jj <MatLeng; jj++) {
       OutputMatrix[cno] = 0;
       cno++;
    }
  }
  cno = 0;
  OutputMatrix[0] = 1;
  for (ii = 0; ii < MatLeng; ii++) {
     OutputMatrix[cno] =1;
     cno += MatLeng + 1;
  }
  /*  if (MatLeng < 8) { 
     Rprintf( (char *) " Start Inverse Matrix Solution:  \n");
     for (ii = 0; ii < MatLeng; ii++) {
	     Rprintf( (char *) " [  ");
	     for (jj = 0; jj < MatLeng; jj++) {
		      Rprintf( (char *) "  %.3f ", (double) InputMatrix[ MatLeng * jj + ii ] );
	     }
	     Rprintf( (char *) " ] \n");
     }
     R_FlushConsole();
  } */
  KiiS = -MatLeng;
  for (ii = 0; ii < MatLeng; ii++) {
     KiiS += MatLeng;
     if ( InputMatrix[ KiiS + ii ] == 0.0000) {
       Rprintf( (char *) "Invert Matrix Flaw in Invert Matrix, Zero Determinant, Dim = %d \n", 
              (int) MatLeng);
       R_FlushConsole();
       R_ProcessEvents();
       return(-1);
     } else {
         OnDivisor = InputMatrix[ KiiS + ii ];
     }
     Kiijj = ii;
     for ( jj = 0; jj < MatLeng; jj++) {
          if (jj >= ii) {
	           InputMatrix[ Kiijj ]  =  InputMatrix[Kiijj] / OnDivisor;
	      }
          OutputMatrix[ Kiijj]  = OutputMatrix[Kiijj] / OnDivisor;
          Kiijj += MatLeng;
     }
        for (jj  = 0; jj < MatLeng; jj++) {
             if ( InputMatrix[ KiiS + jj] == 0) {
             } else if (jj == ii) {
             } else {
                OnDivisor = InputMatrix[ KiiS+jj] / InputMatrix [KiiS+ ii];
                Kiijj = 0;
                for (kk = 0; kk < MatLeng; kk++) {
	                if (kk >= ii) {
                       InputMatrix[ Kiijj+jj] = InputMatrix[Kiijj+jj] - InputMatrix[Kiijj+ii] * OnDivisor;
                    }
                    OutputMatrix[ Kiijj+jj] = OutputMatrix[Kiijj+jj] - OutputMatrix[Kiijj + ii] * OnDivisor;
                    Kiijj += MatLeng;                                    
                }
             }     
        
        }
        /*
  if (MatLeng < 8) { 
     Rprintf( (char *) " Iteration %d Inverse Matrix Process:  \n", ii);
     for (i1 = 0; i1 < MatLeng; i1++) {
	     Rprintf( (char *) " [  ");
	     for (j1 = 0; j1 < MatLeng; j1++) {
		      Rprintf( (char *) "  %.4f ", (double) OutputMatrix[ MatLeng * j1 + i1 ] );
	     }
	     Rprintf( (char *) " ] \n");
     }
     R_FlushConsole();
  }
   if (MatLeng < 8) { 
     Rprintf( (char *) " Iteration %d Source Matrix :  \n", ii);
     for (i1 = 0; i1 < MatLeng; i1++) {
	     Rprintf( (char *) " [  ");
	     for (j1 = 0; j1 < MatLeng; j1++) {
		      Rprintf( (char *) "  %.3f ", (double) InputMatrix[ MatLeng * j1 + i1 ] );
	     }
	     Rprintf( (char *) " ] \n");
     }
     R_FlushConsole();
  } */
  }
  //if (MatLeng < 8) { 
  //   Rprintf( (char *) " Inverse Matrix Solution:  \n");
  //   for (ii = 0; ii < MatLeng; ii++) {
  //     Rprintf( (char *) " [  ");
  //	     for (jj = 0; jj < MatLeng; jj++) {
  //		      Rprintf( (char *) "  %.6f ", (double) OutputMatrix[ MatLeng * jj + ii ] );
  //	     }
  //	     Rprintf( (char *) " ] \n");
  //   }
  //   R_FlushConsole();
  //} 
  return(1);
}


//////////////////////////////////////////////////////
//  MatTimesVec Functions
//
//    Multiply a Matrix times a vector, Use Lapack
//
int MatTimesVec(int kOnLen, int NLen, double *XXA, double *InputV, 
       double *OutPutV) {
	//Rprintf( (char *) "MatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	if (XXA == NULL) {
		Rprintf("MatTimesVec Error XXA is NULL\n");
		R_FlushConsole();
		return(-1);
    }
	if (InputV == NULL) {
		Rprintf("MatTimesVec Error InputV is NULL\n");
		R_FlushConsole();
		return(-1);
    }   
	if (OutPutV == NULL) {
		Rprintf("MatTimesVec Error OutPutV is NULL\n");
		R_FlushConsole();
		return(-1);
    }  
     int One = 1; double ZeroD = 0.0; double OneD = 1.0;
    //	F77_CALL(dscal)(&kOnLen, &ZeroD, OutPutV, &One);
      F77_CALL(dgemv)("N", &NLen, &kOnLen,
		       &OneD, XXA, &NLen,
		       InputV, &One, &ZeroD,
		       OutPutV, &One); 
      return(1);    
	//R_FlushConsole();
	int Onnn = 0, jj;
	int OnFun = 0;
	for (Onnn = 0; Onnn < NLen; Onnn++) {
		OutPutV[Onnn] = (long double) 0;
    }
    OnFun = 0;
    for (jj = 0; jj < kOnLen; jj++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Onnn] += (long double) XXA[OnFun]  * (long double) InputV[jj];
		    OnFun++;
        }
    }
    return(1);
}
//////////////////////////////////////////////////////
//  MatTimesVec Functions
//
//    Multiply a Matrix times a vector, hardly hot stuff.
//
int MatTimesVec(int kOnLen, int NLen, long double *XXA, long double *InputV, long double *OutPutV) {
	//Rprintf( (char *) "MatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	if (XXA == NULL) {
		Rprintf("MatTimesVec Error XXA is NULL\n");
		R_FlushConsole();
		return(-1);
    }
	if (InputV == NULL) {
		Rprintf("MatTimesVec Error InputV is NULL\n");
		R_FlushConsole();
		return(-1);
    }   
	if (OutPutV == NULL) {
		Rprintf("MatTimesVec Error OutPutV is NULL\n");
		R_FlushConsole();
		return(-1);
    }  
     //int One = 1; double ZeroD = 0.0;  
    //	F77_CALL(dscal)(&kOnLen, &ZeroD, OutPutV, &One);
     // F77_CALL(dgemv)("N", &NLen, &kOnLen,
		 //    &OneD, XX, &NLen,
		 //      InputV, &One, &ZeroD,
		 //      OutPutV, &One); 
     // return(1);    
	//R_FlushConsole();
	int Onnn = 0, jj;
	int OnFun = 0;
	for (Onnn = 0; Onnn < NLen; Onnn++) {
		OutPutV[Onnn] = (long double) 0;
    }
    OnFun = 0;
    for (jj = 0; jj < kOnLen; jj++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Onnn] += (long double) XXA[OnFun]  * (long double) InputV[jj];
		    OnFun++;
        }
    }
    return(1);
}
int MatTimesVecVec(int kOnLen, int NLen, long double *XXA, long double *WeightV, 
                  long double *InputV, long double *OutPutV) {
	//Rprintf( (char *) "MatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	if (XXA == NULL) {
		Rprintf("MatTimesVecVec Error XXA is NULL\n");
		R_FlushConsole();
		return(-1);
    }
	if (WeightV == NULL) {
		Rprintf("MatTimesVecVec Error WeightV is NULL\n");
		R_FlushConsole();
		return(-1);
    }       
	if (InputV == NULL) {
		Rprintf("MatTimesVecVec Error InputV is NULL\n");
		R_FlushConsole();
		return(-1);
    }   
	if (OutPutV == NULL) {
		Rprintf("MatTimesVecVec Error OutPutV is NULL\n");
		R_FlushConsole();
		return(-1);
    }
     	
	//R_FlushConsole();
	int Onnn = 0, jj;
	int OnFun = 0;
	for (Onnn = 0; Onnn < NLen; Onnn++) {
		OutPutV[Onnn] = (long double) 0;
    }
    OnFun = 0;
    for (jj = 0; jj < kOnLen; jj++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Onnn] += (long double) XXA[OnFun]  * (long double) InputV[jj] * 
		                     (long double) WeightV[jj];
		    OnFun++;
        }
    }
    return(1);
}

int MatTimesVecVec(int kOnLen, int NLen, double *XXA, double *WeightV, 
                  double *InputV, double *WIV, double *OutPutV) {
	//Rprintf( (char *) "MatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	if (XXA == NULL) {
		Rprintf("MatTimesVecVec Error XXA is NULL\n");
		R_FlushConsole();
		return(-1);
    }
	if (WeightV == NULL) {
		Rprintf("MatTimesVecVec Error WeightV is NULL\n");
		R_FlushConsole();
		return(-1);
    }       
	if (InputV == NULL) {
		Rprintf("MatTimesVecVec Error InputV is NULL\n");
		R_FlushConsole();
		return(-1);
    }   
	if (OutPutV == NULL) {
		Rprintf("MatTimesVecVec Error OutPutV is NULL\n");
		R_FlushConsole();
		return(-1);
    }
     int Zero = 0; int One = 1; 
    double OneD = 1;    double ZeroD = 0.0;
F77_CALL(dsbmv)("U", &kOnLen, &Zero,
		&OneD, WeightV, &One,
		InputV, &One,
		&ZeroD, WIV, &One );
F77_NAME(dgemv)("N", &NLen, &kOnLen,
		&OneD, XXA, &NLen,
		WIV, &One, &ZeroD,
		OutPutV, &One);      
     	
	//R_FlushConsole();
	int Onnn = 0, jj;
	int OnFun = 0;
	for (Onnn = 0; Onnn < NLen; Onnn++) {
		OutPutV[Onnn] = (long double) 0;
    }
    OnFun = 0;
    for (jj = 0; jj < kOnLen; jj++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Onnn] += (long double) XXA[OnFun]  * (long double) InputV[jj] * 
		                     (long double) WeightV[jj];
		    OnFun++;
        }
    }
    return(1);
}

int SMatTimesVecVec(int kOnLen, int NLen, double *XXA, double *WeightV, 
                  double *InputV, double *WIV, double *OutPutV) {
   int One = 1;    int Zero = 0;
   double OneD = 1.0; double ZeroD = 0.0;
	//Rprintf( (char *) "MatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	if (XXA == NULL) {
		Rprintf("MatTimesVecVec Error XXA is NULL\n");
		R_FlushConsole();
		return(-1);
    }
	if (WeightV == NULL) {
		Rprintf("MatTimesVecVec Error WeightV is NULL\n");
		R_FlushConsole();
		return(-1);
    }       
	if (InputV == NULL) {
		Rprintf("MatTimesVecVec Error InputV is NULL\n");
		R_FlushConsole();
		return(-1);
    }   
	if (OutPutV == NULL) {
		Rprintf("MatTimesVecVec Error OutPutV is NULL\n");
		R_FlushConsole();
		return(-1);
    }
F77_CALL(dsbmv)("U", &kOnLen, &Zero,
		&OneD, WeightV, &One,
		InputV, &One,
		&ZeroD, WIV, &One );
F77_CALL(dsymv)("U", &kOnLen, &OneD,
		XXA, &kOnLen,
		WIV, &One,
		&ZeroD, OutPutV, &One);

     	
	//R_FlushConsole();
	int Onnn = 0, jj;
	int OnFun = 0;
	for (Onnn = 0; Onnn < NLen; Onnn++) {
		OutPutV[Onnn] = (long double) 0;
    }
    OnFun = 0;
    for (jj = 0; jj < kOnLen; jj++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Onnn] += (long double) XXA[OnFun]  * (long double) InputV[jj] * 
		                     (long double) WeightV[jj];
		    OnFun++;
        }
    }
    return(1);
}

int tMatTimesVec(int kOnLen, int NLen,  double *XX,long  double *InputV, long double *OutPutV) {
	//Rprintf( (char *) "tMatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	//R_FlushConsole();	
	int Okkk = 0, Onnn;
	int OnFun = 0;
	for (Okkk = 0; Okkk < kOnLen; Okkk++) {
		OutPutV[Okkk] = (long double) 0;
    }
    OnFun = 0;
    for (Okkk = 0; Okkk < kOnLen; Okkk++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Okkk] += (long double) XX[OnFun]  * (long double) InputV[Onnn];
		    OnFun++;
        }
    }
    return(1);
}


int tMatTimesVec(int kOnLen, int NLen,  long double *XX,  long double *InputV, long double *OutPutV) {
	//Rprintf( (char *) "tMatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	//R_FlushConsole();	
	int Okkk = 0, Onnn;
	int OnFun = 0;
	for (Okkk = 0; Okkk < kOnLen; Okkk++) {
		OutPutV[Okkk] = (long double) 0;
    }
    OnFun = 0;
    for (Okkk = 0; Okkk < kOnLen; Okkk++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Okkk] += (long double) XX[OnFun]  * (long double) InputV[Onnn];
		    OnFun++;
        }
    }
    return(1);
}
int tMatTimesVec(int kOnLen, int NLen,  double *XX,  
                 double *InputV, double *OutPutV) {
	//Rprintf( (char *) "tMatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	//R_FlushConsole();	
	 int One = 1; double OneD = 1; double ZeroD = 0.0;
	 
	    F77_CALL(dscal)(&kOnLen, &ZeroD, OutPutV, &One);
      F77_CALL(dgemv)("T", &NLen, &kOnLen,
		     &OneD, XX, &NLen,
		       InputV, &One, &ZeroD,
		       OutPutV, &One); 
      return(1);   
	int Okkk = 0, Onnn;
	int OnFun = 0;
	for (Okkk = 0; Okkk < kOnLen; Okkk++) {
		OutPutV[Okkk] = (long double) 0;
    }
    OnFun = 0;
    for (Okkk = 0; Okkk < kOnLen; Okkk++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Okkk] += (long double) XX[OnFun]  * (long double) InputV[Onnn];
		    OnFun++;
        }
    }
    return(1);
}
int tMatTimesVec(int kOnLen, int NLen,  long double *XX,  long double *InputV, double *OutPutV) {
	//Rprintf( (char *) "tMatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
	//R_FlushConsole();	
	int Okkk = 0, Onnn;
	int OnFun = 0;
	for (Okkk = 0; Okkk < kOnLen; Okkk++) {
		OutPutV[Okkk] = (double) 0;
    }
    OnFun = 0;
    for (Okkk = 0; Okkk < kOnLen; Okkk++) {
	    for (Onnn = 0; Onnn < NLen; Onnn++) {
		    OutPutV[Okkk] += (double) ((long double) XX[OnFun]  * (long double) InputV[Onnn]);
		    OnFun++;
        }
    }
    return(1);
}

//int tMatTimesVec(int kOnLen, int NLen,  double *XX,  double *InputV, double *OutPutV) {
//	//Rprintf( (char *) "tMatTimesVec, kOnLen = %d, NLen = %d\n", kOnLen, NLen);
//	//R_FlushConsole();	
//	int Okkk = 0, Onnn;
//	int OnFun = 0;
//	for (Okkk = 0; Okkk < kOnLen; Okkk++) {
//		OutPutV[Okkk] = ( double) 0;
//  }
//    OnFun = 0;
//    for (Okkk = 0; Okkk < kOnLen; Okkk++) {
//	    for (Onnn = 0; Onnn < NLen; Onnn++) {
//		    OutPutV[Okkk] += (double) XX[OnFun]  * (double) InputV[Onnn];
//		    OnFun++;
//        }
//    }
//    return(1);
//}



///////////////////////////////////////////////////////////////////////////////
//  C Coded Cholesky decomposition
//     I'm not sure that this is used in any code at present, though useful
//
int Cholesky(int Nlen, long double *InputMatrix, long double *OutputMatrix) {
	//long double Holdaii;
	long double Remainingtotal;
	long double OnTotal;
	int lociiii = 0;
	int ii,jj,kk;
	int locjjjj;
	int lociikk, locjjkk;
	int lociijj;
	for (ii = 0; ii < Nlen; ii++)  {
		jj = 0;
		Remainingtotal = InputMatrix[lociiii];
		locjjjj = 0;
		lociijj = ii;
		//Rprintf( (char *) "  %d   [", ii );
		if ( ii > 0 ) { 
		   for (jj = 0; jj < ii; jj++) {
			if (jj == ii) {
				 break;
	        }

	        lociikk = ii;
	        locjjkk = jj;
	        OnTotal = InputMatrix[lociijj];
	        if (jj == 0) {
		        OutputMatrix[lociijj] = OnTotal/(OutputMatrix[locjjjj]);
	        } else {
                for (kk = 0; kk  < jj; kk++) {
	                OnTotal = OnTotal - (OutputMatrix[lociikk] * OutputMatrix[locjjkk]);
	                lociikk += Nlen;
	                locjjkk += Nlen;
                }
                OutputMatrix[lociijj] = OnTotal / OutputMatrix[locjjjj];
            }
            //Rprintf( (char *) "  %d = %.4f   ", (int) lociijj, (double) OutputMatrix[lociijj]);
            //R_FlushConsole();
            Remainingtotal = Remainingtotal - (OutputMatrix[lociijj] * OutputMatrix[lociijj]);
            lociijj += Nlen;
            locjjjj += Nlen + 1;
          }
        }
        if (Remainingtotal < 0) {
	         Rprintf( (char *) "\nCholesky Flaw at ii = %d, Nlen = %d, Remainintotal = %.4f, lociiii = %d, lociijj = %d\n",
	                       ii, Nlen, (double) Remainingtotal, lociiii, lociijj);
	         Rprintf( (char *) "Start total was %.4f", InputMatrix[lociiii]);
	         R_FlushConsole();
	         R_ProcessEvents();
	         //Rprintf( (char *) "Working MatrixL:\n");
	         //PrintRMatrix(InputMatrix, Nlen,Nlen);
	         //Rprintf( (char *) "Chol Finished:\n");
	         //PrintRMatrix(OutputMatrix, Nlen, Nlen);
	         return(-1);
        }
        OutputMatrix[lociiii] = sqrt(Remainingtotal);
        //Rprintf( (char *) "   %d = %.4f   \n", lociiii, (double) OutputMatrix[lociiii]);
        //R_FlushConsole();
        lociiii+=Nlen+1;
	
    }
    return(1);
}



///////////////////////////////////////////
//  SetToZero
//    Sets Square matrix all to zero.   
void SetToZero(int Nlen, long double *Matrix) {
	int Nleno = Nlen * Nlen,ii;
	for (ii = 0; ii < Nleno; ii++) {
		Matrix[ii] = 0.0;
    }
return;
}


////////////////////////////////////////////////////////////////////////////
//   int MatTMat (.... )
//   Multiples Matrix to it's transpose creates corresponding Square Gram Matrix.
//
//    X^TX = G

int MatTMat (long double *RetMat, int n1, int n2, double *NMat, int m1, int m2, double * MMat) {
  if (n2 != m1) {
	  Rprintf( (char *) "MatTMat Error n2 = %d, m1 = %d \n", n2, m1);
	  R_FlushConsole();
	  return(-1);
  }
  int OnNMat, OnMMat, OnNMatT, OnMMatT;
  double WorkAns;
  int ii,jj,kk;
  for (ii = 0; ii < n1; ii++) {
	  OnNMat = ii;
	  for (jj = 0; jj < m2; jj++) {
		  WorkAns = 0;
		  OnMMat = jj * m1;
		  OnNMatT = OnNMat;
		  OnMMatT = OnMMat;
		  for (kk = 0; kk < n2; kk++) {
			  WorkAns+= NMat[OnNMatT] * MMat[OnMMatT];
			  OnNMatT += n1;
			  OnMMatT++;
	      }
	      RetMat[ n1 * jj + ii ] = (long double) WorkAns;
      }
   }
   return(1);
}
///////////////////////////////////////////////////////////////
// int SqMat1 (long double *RetMat, int n1, int n2, double *InpMat)
//
//    Squares a Matrix with one theory method
//
//
int SqMat1 (long double *RetMat, int n1, int n2, double *InpMat) {
    int ii,jj;
    double WorkAns;
    int Onijk1, Onijk2;
    int kk;
    for (ii = 0; ii < n1; ii++) {
	    WorkAns = 0;
	    Onijk1 = ii;
	    for (kk = 0; kk < n2; kk++)  {
		    WorkAns+=(double) InpMat[Onijk1] * (double) InpMat[Onijk1];
		    Onijk1 += n1;
        }
        RetMat[ ii * n1 + ii ] = WorkAns;
        if (ii < n1-1) {
	        for (jj = ii+1; jj < n1; jj++) {
		        Onijk1 = ii;
		        Onijk2 = jj;
		        WorkAns = 0;
		        for (kk = 0; kk < 2; kk++) {
			        WorkAns+=(double) InpMat[Onijk1] * (double) InpMat[Onijk2];
			        Onijk1 += n1;
			        Onijk2 += n1;
		        }
		        RetMat[ ii * n1 + jj ] = WorkAns;
		        RetMat[ jj * n1 + ii ] = WorkAns;
	        }
        }
    }
    return(1);		        
}
/////////////////////////////////////////////////////////////
//  int SqMat () 
//
//   Squares the RetMat putting it into InpMat
//
//
//
//
int SqMat (long double *RetMat, int n1, int n2, double *InpMat) {
    int ii,jj;
    double WorkAns;
    int Onijk1, Onijk2;
    int kk;
    for (ii = 0; ii < n2; ii++) {
	    WorkAns = 0;
	    Onijk1 = ii * n1;
	    for (kk = 0; kk < n1; kk++)  {
		    WorkAns+=(double) InpMat[Onijk1] * (double) InpMat[Onijk1];
		    Onijk1 ++;
        }
        RetMat[ ii * n2 + ii ] = WorkAns;
        if (ii < n2-1) {
	        Onijk2 = (ii+1) * n1;
	        for (jj = ii+1; jj < n2; jj++) {
		        Onijk1 = ii * n1;
		        WorkAns = 0;
		        for (kk = 0; kk < n1; kk++) {
			        WorkAns+=(double) InpMat[Onijk1] * (double) InpMat[Onijk2];
			        Onijk1 ++;
			        Onijk2 ++;
		        }
		        RetMat[ ii * n2 + jj ] = (long double) WorkAns;
		        RetMat[ jj * n2 + ii ] = (long double) WorkAns;
	        }
        }
    }
    return(1);		        
}
int ReSquareLU(double *RetMat, int kLen) {
   int nniijj = 0; int nnjjii = 0;
   int jj; int ii;
   for (jj = 1; jj < kLen; jj++) {
      nniijj = kLen * jj; nnjjii = jj;
      for (ii = 0; ii <jj; ii++) {
          RetMat[nniijj] = RetMat[nnjjii];
          nniijj++; nnjjii+=kLen;
      }
   }
   return(1);
}
int PackSymMatrixLP(double *InMat, double *OutMat, int kLen) {
   int nniijj = 0; int nnjjii = 0;
   int jj; int ii; int nnO = 0;
   for (jj = 1; jj < kLen; jj++) {
      nniijj = kLen * jj; nnjjii = jj;
      for (ii = 0; ii <jj; ii++) {
          OutMat[nnO] = InMat[nnjjii];
          nniijj++; nnjjii+=kLen;
          nnO++;
      }
   }
   return(1);
}
int PackSymMatrixUP(double *InMat, double *OutMat, int kLen) {
   int nniijj = 0; int nnjjii = 0;
   int jj; int ii; int nnO = 0;
   for (jj = 1; jj < kLen; jj++) {
      nniijj = kLen * jj; nnjjii = jj;
      for (ii = 0; ii <jj; ii++) {
          OutMat[nnO] = InMat[nniijj];
          nniijj++; nnjjii+=kLen;
          nnO++;
      }
   }
   return(1);
}
int ReSquareUL(double *RetMat, int kLen) {
   int nniijj = 0; int nnjjii = 0;
   int jj; int ii;
   for (jj = 1; jj < kLen; jj++) {
      nniijj = kLen * jj; nnjjii = jj;
      for (ii = 0; ii <jj; ii++) {
          RetMat[nnjjii] = RetMat[nniijj];
          nniijj++; nnjjii+=kLen;
      }
   }
   return(1);
}


// What's Going Down?
int pSqMat (double **RetMat, int n1, int n2, double *InpMat) {
    ////Rprintf("SqMat: BLAS Version \n");
    double OneD = 1.0;  double ZeroD = 0.0;
/* DSYRK - perform one of the symmetric rank k operations */
/* C := alpha*A*A' + beta*C or C := alpha*A'*A + beta*C */
   char MyTrans = 'T';
   double MyTotal = 0.0;
   int One = 1;
   int NLeft = 0;
   for (int ii = 0; ii < n2-1; ii++) {
     for (int jj = 0; jj < n2; jj++) {
       *(RetMat[ii] + jj) = 0.0;
     }
     if (ii > 0) {
       for (int jj = 0; jj < ii; jj++) {
         *(RetMat[ii] + jj) = *(RetMat[jj] + ii);
       }
     }  
     MyTotal = 0.0;
     for (int jj = 0; jj < n1; jj++) {
       MyTotal += InpMat[ii*n1 + jj] * InpMat[ii*n1+jj];
     }
     NLeft = n2 - ii - 1;
     if (NLeft > 0) {
       F77_CALL(dgemv)(&MyTrans, &n1, &n2, &OneD, InpMat + n1*(ii+1), &n1, InpMat+ii*n1, &One, &ZeroD, RetMat[ii]+ii+1, &One);
     }
   }
 // Rprintf("SqMat After Blas \n");
 //  PrintRMatrix( RetMat, n2, n2);
 //  Rprintf("AfterFixing \n");
 //  PrintRMatrix(RetMat, n2,n2);
 // R_FlushConsole(); R_ProcessEvents();
  return(1);		
}
// What's Going Down?
int SqMat (double *RetMat, int n1, int n2, double *InpMat) {
    ////Rprintf("SqMat: BLAS Version \n");
    double OneD = 1.0;  double ZeroD = 0.0;
/* DSYRK - perform one of the symmetric rank k operations */
/* C := alpha*A*A' + beta*C or C := alpha*A'*A + beta*C */
F77_CALL(dsyrk)("U", "T", &n2,&n1, &OneD, InpMat,
     &n1, &ZeroD, RetMat, &n2 );
 // Rprintf("SqMat After Blas \n");
 //  PrintRMatrix( RetMat, n2, n2);
 //  Rprintf("AfterFixing \n");
   ReSquareUL(RetMat, n2);
 //  PrintRMatrix(RetMat, n2,n2);
 // R_FlushConsole(); R_ProcessEvents();
  return(1);		
    int ii,jj;
    double WorkAns;
    int Onijk1, Onijk2;
    int kk;
    for (ii = 0; ii < n2; ii++) {
	    WorkAns = 0;
	    Onijk1 = ii * n1;
	    for (kk = 0; kk < n1; kk++)  {
		    WorkAns+=(double) InpMat[Onijk1] * (double) InpMat[Onijk1];
		    Onijk1 ++;
        }
        RetMat[ ii * n2 + ii ] = WorkAns;
        if (ii < n2-1) {
	        Onijk2 = (ii+1) * n1;
	        for (jj = ii+1; jj < n2; jj++) {
		        Onijk1 = ii * n1;
		        WorkAns = 0;
		        for (kk = 0; kk < n1; kk++) {
			        WorkAns+=(double) InpMat[Onijk1] * (double) InpMat[Onijk2];
			        Onijk1 ++;
			        Onijk2 ++;
		        }
		        RetMat[ ii * n2 + jj ] = WorkAns;
		        RetMat[ jj * n2 + ii ] = WorkAns;
	        }
        }
    }
  Rprintf("SqMat After NonBlas \n");
  PrintRMatrix( RetMat, n2, n2);
  R_FlushConsole(); R_ProcessEvents();    
    return(1);		        
}



int CopyMatrix(int n1, int n2, long double *InputMatrix, long double *OutputMatrix) {
   int MM2 = n1 * n2;
   int nn;
   if ((OutputMatrix!= NULL) && InputMatrix != NULL) {
		   for (nn = 0; nn < MM2; nn++) {
			   OutputMatrix[nn] = InputMatrix[nn];	
		   }
		   return(1);	
   }
   return(-1);
}        


int ReduceFormatToOne(int KappaT, int Onk2, long double *OnDiag, long double *OffDiag, 
                    long double *OutW11) {
   int onii = 0, onjj = onii+1, inii = 0, injj = inii+1;
   int nnii = 0; int nnjj = 0; int nn;
   int STO = KappaT * (KappaT-1)/2;
   nnii = injj * (KappaT-1) + inii;
   nnjj = inii * (KappaT -1) + injj;
   for (nn = 0; nn < STO; nn++) {
	   if (onii == Onk2) {
		   nn += (KappaT-1) - onii-1;
		   onii++; onjj = onii+1;
       } else if (onjj == Onk2) {
	       onjj++;
       } else if (onjj >= KappaT-2) {
	       OutW11[nnii] = OffDiag[nn];
	       OutW11[nnjj] = OffDiag[nn];
	       onii = onii +1; onjj = onii+1;
	       inii = inii +1; injj = inii +1;
	       nnii = injj * (KappaT-1) + inii;
	       nnjj = inii * (KappaT-1) + injj;
       } else {
	       OutW11[nnii] = OffDiag[nn];    
	       OutW11[nnjj] = OffDiag[nn];
	       onjj++; injj++; nnjj++;
	       nnii += (KappaT-1);
       }
   }
   nnii = 0;
   for (nn = 0; nn < KappaT-1; nn++) {
	   if (nn != Onk2) {
	      OutW11[nnii] = OnDiag[nn];	      
       }
           nnii +=KappaT;
   }
   return(1);
}

int OnToOff(int KappaT, int Onk2, double *InptSS,
                    long double *OnDiag, long double *OffDiag) {
   int onii = 0, onjj = onii+1;
   int nnjj = 0; int nn;
   int STO = KappaT * (KappaT-1)/2;
   nnjj = onii * (KappaT) + nnjj;
   for (nn = 0; nn < STO; nn++) {
        if (onjj >= KappaT-1) {
	       OffDiag[nn] = InptSS[nnjj]; onii++;
	       onjj = onii+1;
	       nnjj = onii * (KappaT) + onjj;
       } else {
	       OffDiag[nn] = InptSS[nnjj];    
	       onjj++; nnjj++;
       }
   }
   nnjj = 0;
   for (nn = 0; nn < KappaT-1; nn++) {
		  OnDiag[nn] = InptSS[nnjj];
	      nnjj +=KappaT+1;
   }
   return(1);
}



////////////////////////////////////////////////////////////////////////////
// Solve a Symmetric Positive Matrix
//   dimk, is dimension of the matrix
//   Input, is pointer to input matrix (will be unchanged)
//   Output, is pointer to output matrix (will be changed.);
int SolveSymPosMatrix(int dimk, double *Input, double *Output) {
    int Pointdimk[1];  Pointdimk[0] = dimk;
    int ii, jj;  int FlagP[1]; FlagP[0] = 0;
    //int CurCount  = 0;
    //  Because there is no guarantee that Outmatrix is zeroed, which we need,
    int dimkTdimk = dimk * dimk; int One = 1;
    //    We loop through an zero OutMatrix at our discretion.
    //F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
	  //	double *dy, const int *incy);
    F77_CALL(dcopy)(&dimkTdimk, Input,  &One, Output, &One);
    //for (jj = 0; jj < dimk; jj++) {  for (ii = 0; ii < dimk; ii++) {
    //      Output[CurCount] = 0;
    //      if ( ii >= jj) { Output[CurCount] = Input[CurCount]; }
    //      CurCount++;
    //    }
    //}
    ////////////////////////////////////////////////////////////////////
    //  Inversion by Cholesky transform and then inverse program
    // 
    F77_CALL(dpotrf)((char*) "L", (const int *) Pointdimk, 
                      (double*) Output, 
                      (const int *) Pointdimk, 
                      (int *) FlagP);   
    F77_CALL(dpotri)("L", Pointdimk, Output, Pointdimk, FlagP);  
    int Lociijj = dimk+1;
    int Locjjii = 1;
    for (jj = 1; jj < dimk; jj++) {
       Locjjii = jj; Lociijj = dimk * jj;
       for (ii = 0; ii <= jj; ii++) {
          Output[Lociijj] = Output[Locjjii];
          Lociijj++; Locjjii+=dimk;
       }
    }
    return(FlagP[0]);
}

	
