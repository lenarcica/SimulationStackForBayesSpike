/* ========================================================================== */
/*                                                                            */
/*   IncidentalOp2009.cc                                                       */
/*   (c) 2010 Alan Lenarcic                                                   */
/*                                                                            */
/*   C code not used often, mostly to free memory and set matrix positions.   */
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



#ifndef INCIDENTALOP2009DD
  #include "IncidentalOp2009.h"
  #define INCIDENTALOP2009DD 0 
#endif


////////////////////////////////////
//   int GetSign(double value)
//
//    Returns vanilla integer sign
//     (no running minimum value);
int GetSign(double Value) {
	if (Value < 0) {
		return(-1);
    } else if (Value == 0) {
	    return(0);
    } else {
	    return(1);
    }
} 
  
////////////////////////////////////
//  ceiling
//
//     Lowest integer greater than Inpt;
int ceiling (double Inpt) {
	if ((int) Inpt == Inpt ) {
		 return((int)Inpt);
    } else {
	     return( ( (int) Inpt) + 1 ) ;
    }
}

/////////////////////////////////////////////////
//  SetMeI functions are turned off.
//
//
//
int SetMeI ( int kSLen, long double *GMat ) { 
  int ii, jj, Onii = 0;
  for (	ii = 0; ii < kSLen; ii++) {
	  for (jj = 0; jj < kSLen; jj++) {
		  GMat[Onii] = 0;
		  Onii++;
      }
      GMat[ ii * kSLen + ii] = 1;
  }
  return(1);
}
int SetMeXIn (int kLen, int NLen, int OnL, 
              int* ActiveOns, int *ActiveSigns, double *XInputMat , 
              long double *XOutputMat) {
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
int SetMeXIn (int kLen, int NLen, int OnL, int* ActiveOns,  int*ActiveSigns,
              long double *XInputMat , 
              long double *XOutputMat) {
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

void FFree(void **MyPointer) {
  if (MyPointer[0] != NULL) {
	  Free(MyPointer[0]);
	  MyPointer[0] = NULL;
   }
   return;  		
}

void FFree(double **MyPointer) {
  if (MyPointer[0] != NULL) {
	  Free(MyPointer[0]);
	  MyPointer[0] = NULL;
   }
   return;  		
}

void FFree(long double **MyPointer) {
  if (MyPointer[0] != NULL) {
	  Free(MyPointer[0]);
	  MyPointer[0] = NULL;
   }
   return;  		
}

void FFree(int **MyPointer) {
  if (MyPointer[0] != NULL) {
	  Free(MyPointer[0]);
	  MyPointer[0] = NULL;
   }
   return;  		
}

void FFree(short int **MyPointer) {
  if (MyPointer[0] != NULL) {
	  Free(MyPointer[0]);
	  MyPointer[0] = NULL;
   }
   return;  		
}

int EqualToZero(long double A) {
	const long double MinMin = .00000000000000000001;
	if (A < MinMin && A > -MinMin) {
		return(1);
    } else {
        return(0);
    }
	
}
int EqualToZero(double A) {
	const double MinMin = .00000000000000000001;
	if (A < MinMin && A > -MinMin) {
		return(1);
    } else  {
        return(0);
    }
	
}

int IntMin(int A1, int A2) {
   if (A1 < A2) { return(A1); }
   return(A2);	
}

