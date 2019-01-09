/* ========================================================================== */
/*                                                                            */
/*   HorseShoeC.cc                                                            */
/*   (c) 2010 Alan Lenarcic                                                   */
/*                                                                            */
/*   Used to Integrate Lambdaj Densities                                      */
/*     This attempt to implement Horseshoe is not very good and it would  be  */
/*     recommend users look into users of TwoSimR5 for more professional      */
/*     versions of this estimator.                                            */
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
#ifndef HorseShoeCDef
  #include "HorseShoeC.h"
  #define HorseShoeCDef 1
#endif 

int GetAllLambdaj(SEXP SlZVals, SEXP SdlStandard, SEXP SBetaj, SEXP StauSq,
  SEXP SFillerDens, SEXP SUnifDraw, SEXP SOutLambdaj, SEXP Sexpn2lz2) {
 int Needs = length(SBetaj);
 int ii;
 for (ii = 0; ii < Needs; ii++) {
   REAL(SOutLambdaj)[ii]= GetLambdaj(length(SlZVals), REAL(SdlStandard),
       REAL(SBetaj)[ii], REAL(StauSq)[0], REAL(SlZVals),
       REAL(SFillerDens), REAL(SUnifDraw)[0], REAL(Sexpn2lz2));
 }
 return(1);
}

int AllGetManyLambdaj(SEXP SlZVals, SEXP SdlStandard, SEXP SOneBetaj, SEXP StauSq,
  SEXP SFillerDens, SEXP SortedUnifDraw, SEXP SOutLambdaj, SEXP Sexpn2lz2) {
 
  //Rprintf("To Run AllGetManyLambdaj  Fun \n"); R_FlushConsole();
  GetManyLambdaj(length(SlZVals), REAL(SdlStandard),
       REAL(SOneBetaj)[0], REAL(StauSq)[0], REAL(SlZVals),
       REAL(SFillerDens), length(SortedUnifDraw),
       REAL(SortedUnifDraw), REAL(SOutLambdaj), REAL(Sexpn2lz2));
 return(1);
}

double GetLambdaj(int LengthZ2, double *dlStandard, double Betaj, 
     double tauSq, double *lZ2, double *FillerDens, double GetDraw, double *expn2lz2) {
  if (GetDraw == 0) {
    return(0.0);
  }
  if (GetDraw >= 1.0) {
     return( exp(lZ2[LengthZ2-1]));
  }
  int ii;
  int One = 1;
  F77_CALL(dcopy)(&LengthZ2, dlStandard, &One,
  		  FillerDens, &One);
  double Multiplier = - Betaj * Betaj / ( 2 * tauSq);
  for (ii = 0; ii <LengthZ2; ii++) {
     FillerDens[ii] *= exp( Multiplier * expn2lz2[ii]);
  }                               
  double Integrator = F77_CALL(dasum)(&LengthZ2, FillerDens, &One);
  double OnVal = 0.0;
  double UnIntegrator = 1.0 / Integrator;
  
  for(ii = 0; ii < LengthZ2; ii++) {
     OnVal += UnIntegrator * FillerDens[ii];
     if (OnVal > GetDraw) {
        return(exp(lZ2[ii]));
     }
  }
  return(100000);	
}   

int GetManyLambdaj(int LengthZ2, double *dlStandard, double Betaj, 
     double tauSq, double *lZ2, double *FillerDens, 
     int CountDraws, double *SortedDraw, double *Output, double *expn2lz2) {

  //Rprintf("Running GetManyLambdaj \n"); R_FlushConsole();
  int ii,jj=0;            if (CountDraws < 0) {return(-1);}
  int One = 1;
  F77_CALL(dcopy)(&LengthZ2, dlStandard, &One,
  		  FillerDens, &One);
  double Multiplier = - Betaj * Betaj / ( 2 * tauSq);
  for (ii = 0; ii <LengthZ2; ii++) {
     FillerDens[ii] *= exp( Multiplier * expn2lz2[ii]);
  } 
  //Rprintf("Done With FillerDens \n"); R_FlushConsole();                              
  double Integrator = F77_CALL(dasum)(&LengthZ2, FillerDens, &One);
  //Rprintf("Integrator Sum is %.4f\n", Integrator); R_FlushConsole();
  double UnIntegrator = 1.0 / Integrator;
  F77_CALL(dscal)(&LengthZ2, &UnIntegrator, FillerDens, &One);
  
  Integrator = F77_CALL(dasum)(&LengthZ2, FillerDens, &One);
  //Rprintf("AfterSum, Integrator is %.4f\n", Integrator); R_FlushConsole();
  
  double GetDraw;
  int Onii = -1;
  double OnIntegrate = 0.0;
  for (jj = 0; jj < CountDraws; jj++) {
     GetDraw = SortedDraw[jj];
     if (GetDraw == 0.0) {
         Output[jj] = 0;
     } else if (GetDraw == 1.0) {
         Output[jj] = 1.0;
     }else {
       while(OnIntegrate < GetDraw && Onii < LengthZ2-1) {
          Onii++;
          OnIntegrate += FillerDens[Onii];    
       }
       Output[jj] = exp(lZ2[Onii]);
     } 
  }  
  return(100000);	
}
         
