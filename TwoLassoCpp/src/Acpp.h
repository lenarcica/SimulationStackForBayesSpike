/* ========================================================================== */
/*                                                                            */
/*   Acpp.h                                                                   */
/*   (c) 2014 Alan Lenarcic                                                   */
/*                                                                            */
/*   Description:                                                             */
/*     Code attempting to duplicate old Rcpp functionality                    */
/* ========================================================================== */

//SEXP stats = PROTECT(
//R_FindNamespace(
//mkString("stats")));

#ifndef ACPPCONFIGUREME
  #include "AcppConfigureMe.h"
  #define ACPPCONFIGUREME 0
#endif 


#define UpLo 'L'
#define TransWant 'N'
#ifndef RNeeds 
 #define RNeeds 0
 #include <R.h>
 #include <Rmath.h>
 #include <Rinternals.h>
 #include <R_ext/BLAS.h>
 #include <R_ext/Lapack.h>
 #include <Rdefines.h>
#endif



//#define StartListLength = 100

#ifndef AOBJECTH
 #define AOBJECTH 0
class AObject {
  SEXP MyP;
  int InVectorLoc;
  SEXP EnvironS;
  public :
  AObject(SEXP aMyP) {
    int tProtect= 0;
    Rprintf(" -- AObject: New Allocate Start. \n"); R_FlushConsole();
    if (Rf_isNull(aMyP)) {
      MyP = R_NilValue;  InVectorLoc = -1;
      EnvironS = R_NilValue;
      return;
    } else {
       SEXP EnvAcpp = R_NilValue;
       Rprintf(" -- AObject: Now Starting to Find Package NameSpace"); R_FlushConsole();
       Rf_protect(EnvAcpp = R_FindNamespace(Rf_mkString(PCKGNAME)));  tProtect++;
       //Rf_protect(EnvAcpp = R_FindNamespace(Rf_mkString(PCKGNAME)));  tProtect++;
       if (Rf_isNull(EnvAcpp)) {
         Rprintf(" -- Oh no, doesn't look like we got the namespace right!\n"); R_FlushConsole();
       }
       if (!Rf_isEnvironment(EnvAcpp)) {
         Rf_error("Error: We had sought a namespace for PCKGNAME %s, not good!\n", PCKGNAME);
       }
       Rprintf(" --  Got NameSpace for Acpp = %s, about to InsertIntoAList.\n", PCKGNAME); R_FlushConsole();
       SEXP sInsertIntoAList =  R_NilValue;
         Rf_protect(sInsertIntoAList = Rf_findVarInFrame( EnvAcpp, Rf_install(".InsertIntoAList")));
         tProtect++;
       EnvironS = R_NilValue;
       EnvironS  = Rf_findVarInFrame( EnvAcpp, Rf_install("ACPPSAVESPACE"));
       if (Rf_isNull(EnvironS)) {
         Rprintf("Acpp: Tried to fine ACPPSAVESPACE in ACppTest package space but it is NULL. \n");
         R_FlushConsole();
       } else {
         //Rprintf("Acpp: Note, we found EnvironS using findVarInFrame on ACppTest.\n"); R_FlushConsole();
       }
       if (!Rf_isEnvironment(EnvironS)) {
         Rprintf("Acpp::Acpp.h:: Not good, ACPPSAVESPACE turned out to be not an Environment!\n");
         SEXP FindACPPSAVE = R_NilValue;
         Rf_protect(FindACPPSAVE = Rf_findVarInFrame( EnvAcpp, Rf_install(".GetACPPSAVESPACE")));
         tProtect++;
         SEXP GetEnvironS = Rf_protect(Rf_lcons( FindACPPSAVE, R_NilValue)); tProtect++;
         Rprintf(" -- AObject: new, about to call InsertIntoAList, wish me luck!\n"); R_FlushConsole();
         EnvironS = R_NilValue;
         EnvironS = Rf_eval(GetEnvironS, EnvAcpp); 
         Rf_unprotect(2); tProtect--;  tProtect--;
       }
       if (!Rf_isEnvironment(EnvironS)) {
         Rprintf("Acpp::Acpp.h:: Environment Error, EnvironS is ACPPSAVESPACE of EnvAcpp, but it is not Enviornment!\n");
         Rf_error("EnvironS Error!\n");
       }
       if (Rf_isNull(sInsertIntoAList)) {
         Rf_error("-- AObject: Error, did not find .InsertIntoAList!\n");
       }
       if (Rf_isNull(sInsertIntoAList)) {
         Rprintf(" -- AObject: new, sInsertIntoAList came back NULL, there is no function!\n"); R_FlushConsole();
       }
       Rprintf(" -- Presumably sInsertIntoAList is Not NULL!\n"); R_FlushConsole();
       //if (Rf_isNull(sTBSRooSetVerbose)) {
       //  Rprintf("BayesSpikeCL: SetVerbose, Error: no sTBSRooSetVerbose");
       //}
       Rprintf("  Evidently we found sTBSRooSetVerbose \n"); R_FlushConsole();
       SEXP call = Rf_protect(Rf_lcons( sInsertIntoAList, Rf_cons(aMyP, R_NilValue))); tProtect++;
       Rprintf(" -- AObject: new, about to call InsertIntoAList, wish me luck!\n"); R_FlushConsole();
       SEXP RTM = R_NilValue;
         Rf_protect(RTM = Rf_eval(call, EnvAcpp)); tProtect++;
       Rprintf(" -- AObject: finished calling InsertIntoAList");
       if (Rf_isNull(RTM) || Rf_length(RTM) <= 0) {
         Rprintf("Acpp: Creation of new Aob failed in .InsertIntoAList call. \n");
       }
       int iRTM = -1;
       if (Rf_isInteger(RTM) ) { iRTM = INTEGER(RTM)[0]; 
       } else if (Rf_isReal(RTM)) { iRTM = (int)REAL(RTM)[0]; }
       if (iRTM >= 0) { InVectorLoc = iRTM; }
       //Rprintf("Acpp:: As far as we know, the object was inserted into location %d\n", 
       //  InVectorLoc);  R_FlushConsole();
       Rf_unprotect(tProtect);    tProtect = 0;
       //Rprintf(" -- Trying to get ActualAList \n"); R_FlushConsole();
       SEXP ActualAList = 
         Rf_protect(Rf_findVarInFrame( EnvironS, Rf_install(".AcppList")));
       tProtect++;
       SEXP sGetAcppList = R_NilValue;
       SEXP Newcall = R_NilValue;
       if (Rf_isNull(ActualAList) || Rf_isReal(ActualAList) || Rf_isInteger(ActualAList) ||
         Rf_isEnvironment(ActualAList) || Rf_isSymbol(ActualAList) || Rf_length(ActualAList) == 1) {
         Rprintf("Acpp:: Oops, Ef_install from EnvironS didn't work!. ActualAList try again!\n");
         R_FlushConsole();
         ActualAList = R_NilValue;
         Rf_unprotect(1);  tProtect--;
         sGetAcppList =  R_NilValue;
           Rf_protect(sGetAcppList = Rf_findVarInFrame( EnvAcpp, Rf_install(".GetAcppList")));
           tProtect++;
          Newcall = Rf_protect(Rf_lcons( sGetAcppList, R_NilValue)); tProtect++;
          //Rprintf(" -- AObject: new, about to call InsertIntoAList, wish me luck!\n"); R_FlushConsole();
          ActualAList = R_NilValue;
          Rf_protect(ActualAList = Rf_eval(Newcall, EnvAcpp)); tProtect++;
       } 
       if (Rf_isNull(ActualAList)) {
         Rf_error("Acpp::new, sorry we don't get an ActualAList! Came back NULL\n");
       }
       if (Rf_length(ActualAList) <= 0) {
         Rf_error("Acpp::new, no, Actual A list has null length!\n");
       }
       if (Rf_isReal(ActualAList)) {
         Rf_error("Acpp::new, why is ActualAList REAL? First member %f \n",
           REAL(ActualAList)[0]);
       }
       if (Rf_isInteger(ActualAList)) {
         Rf_error("Acpp::Acpp.h:: new, why is ActualAList INTEGER? First member %f \n",
           INTEGER(ActualAList)[0]);
       }
       if (Rf_isEnvironment(ActualAList)) {
         Rf_error("Acpp::Acpp.h:: new Error, hey, ActualAList is an environment!\n");
       }
       if (Rf_isExpression(ActualAList)) {
         Rf_error("Acpp::Acpp.h:: new Big Error, how is ActualAList an expression!\n");
       }
       if (Rf_isSymbol(ActualAList)) {
         Rprintf("Acpp::Acpp.h:: new What an Error, why is ActualAList a symbol?\n"); 
         R_FlushConsole();
         if (Rf_isNull(ActualAList)) {
           Rprintf("Acpp:: ActualAList is a NULL \n"); R_FlushConsole();
         }
         Rf_error("Acpp::Acpp.h:: What is Wrong with AcppList as a symbol?\n");
       }
       if (Rf_length(ActualAList) <= InVectorLoc) {
         Rprintf("Acpp:: Error, sorry, InVectorLoc = %d but AList has length %d!\n",
           InVectorLoc, Rf_length(ActualAList));
         if (Rf_isReal(ActualAList) && Rf_length(ActualAList) >= 1) {
          Rprintf("  Real ActualAList[0] is %f \n", REAL(ActualAList)[0]);
         } else if (Rf_isInteger(ActualAList) && Rf_length(ActualAList) >= 1) {
           Rprintf("  INTEGER ActualList[0] is %d \n", INTEGER(ActualAList)[0]);
         } else {
          Rprintf("  Aparently Actual A List is neither Integer or Real but length %d \n",
            Rf_length(ActualAList)); 
         }
         R_FlushConsole();
         Rf_unprotect(tProtect);  Rf_error("Acpp: new error from .AcppList!\n");
       }
       Rprintf("Supposedly ActualAList is length %d\n", Rf_length(ActualAList));
       R_FlushConsole();
       MyP = VECTOR_ELT(ActualAList, InVectorLoc); 
       if (Rf_isReal(MyP) && Rf_isInteger(aMyP)) {
         Rprintf("Acpp::new Sorry, MyP is Real but you gave us aMyP as Int!\n"); R_FlushConsole();
       } else if (Rf_isInteger(MyP) && Rf_isReal(aMyP)) {
         Rprintf("Acpp::new Sorry, MyP is Integer but you gave us aMyP as Real!\n"); R_FlushConsole();
       } else if ( ((Rf_isInteger(MyP) && Rf_isInteger(aMyP)) ||
         (Rf_isReal(MyP) && Rf_isReal(aMyP)))  && Rf_length(aMyP) != Rf_length(MyP)) {
         Rprintf("Acpp::new, bad result, MyP is length %d but given aMyP has length %d!\n",
           Rf_length(MyP), Rf_length(aMyP)); R_FlushConsole();  
       }
       /*
       Rprintf("Acpp::new, we feel we have inserted a valid MyP"); R_FlushConsole();
       if (Rf_isReal(MyP)) {
         Rprintf("-- ActualAList[%d] is an actual real length %d let's print it. \n",
           InVectorLoc, Rf_length(MyP)); R_FlushConsole();
         Rprintf("MyP = c("); R_FlushConsole();
         if (Rf_length(MyP) >= 1) {
           for (int ii = 0; ii < Rf_length(MyP); ii++) {
             Rprintf("%.6f", REAL(MyP)[ii]);
             if (ii == Rf_length(MyP) -1) { Rprintf(")\n");
             } else if ((ii+1) % 8 == 0) { Rprintf(",\n");
             } else { Rprintf(", "); }
             R_FlushConsole();
           }
         }
       }
       if (Rf_isInteger(MyP)) {
         Rprintf("-- ActualAList[%d] is an actual integer length %d let's print it. \n",
           InVectorLoc, Rf_length(MyP)); R_FlushConsole();
         Rprintf("MyP = c("); R_FlushConsole();
         if (Rf_length(MyP) >= 1) {
           for (int ii = 0; ii < Rf_length(MyP); ii++) {
             Rprintf("%.6f", INTEGER(MyP)[ii]);
             if (ii == Rf_length(MyP) -1) { Rprintf(")\n");
             } else if ((ii+1) % 8 == 0) { Rprintf(",\n");
             } else { Rprintf(", "); }
             R_FlushConsole();
           }  
         }
       }  */
       Rf_unprotect(tProtect); tProtect = 0;
    }
  }
  ~AObject() {
       SEXP EnvAcpp = R_NilValue;
       int tProtect = 0;
       /*
       Rprintf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
       Rprintf("~~  ~AObject(): We are actually deleting the object. \n");
       Rprintf("~~  ~ Let's look. \n"); R_FlushConsole();
       */
       Rf_protect(EnvAcpp = R_FindNamespace(Rf_mkString(PCKGNAME)));  tProtect++;
       //if (EnvAcpp != EnvironS) {
       //  Rprintf("Acpp:: Uh oh, Acpp does no returncorrect EnviornS = %d, EnvAcpp = %d!\n",
       //    EnvironS, EnvAcpp);
       //}
       //Rprintf("~~  ! We have DeleteFromAList function. \n"); R_FlushConsole();
       SEXP sDeleteFromAList = 
         Rf_protect(Rf_findVarInFrame( EnvAcpp, Rf_install(".DeleteFromAList")));
       tProtect++;
       SEXP sInVectorLoc = R_NilValue;
       Rf_protect(sInVectorLoc = Rf_allocVector(INTSXP, 1)); tProtect++;
       INTEGER(sInVectorLoc)[0] = InVectorLoc;
 
       EnvironS  = Rf_findVarInFrame( EnvAcpp, Rf_install("ACPPSAVESPACE"));
       if (!Rf_isEnvironment(EnvironS)) {
         //Rprintf("Acpp::Acpp.h:: Not good, On Delete ACPPSAVESPACE turned out to be not an Environment!\n");
         SEXP FindACPPSAVE = R_NilValue;
         Rf_protect(FindACPPSAVE = Rf_findVarInFrame( EnvAcpp, Rf_install(".GetACPPSAVESPACE")));
         tProtect++;
         SEXP GetEnvironS = Rf_protect(Rf_lcons( FindACPPSAVE, R_NilValue)); tProtect++;
         //Rprintf(" -- AObject: new, about to call InsertIntoAList, wish me luck!\n"); R_FlushConsole();
         EnvironS = R_NilValue;
         EnvironS = Rf_eval(GetEnvironS, EnvAcpp); 
         Rf_unprotect(2); tProtect--;  tProtect--;
       }
       SEXP sAcppList = 
         Rf_protect(Rf_findVarInFrame( EnvironS, Rf_install(".AcppList")));   
       tProtect++;
       if (Rf_isNull(sAcppList)) {
         Rprintf(" ~~ AObject Delete, we found Acpp List is NULL. \n"); R_FlushConsole();
       }   
       if (InVectorLoc < 0 && Rf_isNull(MyP)) {
         // Nulls aren't actually preserved in list but can be deleted.
         Rf_unprotect(tProtect); tProtect = 0;
         return;
       }
       if (InVectorLoc < 0 || InVectorLoc >= Rf_length(sAcppList)) {
         Rprintf(" ~~ AObject Delete: Note that InVectorLoc is broken = %d, length AcppList is %d!\n",
           InVectorLoc, Rf_length(sAcppList)); R_FlushConsole();
         Rf_unprotect(tProtect); tProtect = 0;
         Rf_error(" ~~ No AObject Delete: we are not going to delete this member.\n");
       } 
       /*else {
         SEXP IWant = VECTOR_ELT(sAcppList, InVectorLoc);
         if (IWant != MyP) {
           Rprintf(" ~~ AObject Delete: note IWant is AcppList[%d] is not the same as MyP!\n",
             InVectorLoc);  R_FlushConsole();
           Rprintf(" ~~ IWant-real=%d, MyP-real = %d\n",
             Rf_isReal(IWant) ? 1 : 0, Rf_isReal(MyP) ? 1 : 0); R_FlushConsole();
         }
         SET_VECTOR_ELT(sAcppList, InVectorLoc, R_NilValue);
         Rprintf(" ~~ AObject, we manually deleted from InVectorLoc");
         Rf_unprotect(tProtect); MyP = R_NilValue; InVectorLoc = 0;  EnvironS = R_NilValue;
         Rprintf(" ~~ AObject: Note we are quiting delete now. \n"); R_FlushConsole();
         Rf_unprotect(tProtect); tProtect = 0;
         Rprintf(" ~~ This is when we set new sAcppList \n");   R_FlushConsole();
         Rf_setVar(Rf_install(".AcppList"), sAcppList, EnvironS);
         Rprintf(" ~~ Okay, successful delete!\n"); R_FlushConsole();
         return;
       }             */
       
       //if (Rf_isNull(sTBSRooSetVerbose)) {
       //  Rprintf("BayesSpikeCL: SetVerbose, Error: no sTBSRooSetVerbose");
       //}
       //Rprintf("  Evidently we found sTBSRooSetVerbose \n"); R_FlushConsole();
       //Rprintf(" ~~ AObject, attempting to delete using DeleteFromAList\n"); R_FlushConsole();
       SEXP call = Rf_protect(Rf_lcons( sDeleteFromAList, Rf_cons(sInVectorLoc, R_NilValue)));
       tProtect++;
       SEXP RTM = R_NilValue;
       //Rprintf("~~ Now Calling DeleteFromAList function. \n"); R_FlushConsole();
         Rf_protect(RTM = Rf_eval(call, EnvAcpp)); tProtect++;
       if (Rf_isNull(RTM) || Rf_length(RTM) <= 0) {
         Rprintf("Acpp: Deletion of AMyP failed in .DeleteFromAList call %d. \n", 
           InVectorLoc);
       }
       int iRTM = -1;
       MyP = R_NilValue;
       if (Rf_isInteger(RTM) ) { iRTM = INTEGER(RTM)[0]; 
       } else if (Rf_isReal(RTM)) { iRTM = (int)REAL(RTM)[0]; }
       if (iRTM >= 0) { InVectorLoc = iRTM; }
       //Rprintf("Acpp:: As far as we know, the object deleted at location %d\n", 
       //  InVectorLoc);  R_FlushConsole();
       Rf_unprotect(tProtect); tProtect = 0;
  }
  SEXP asSexp() {
    return(MyP);
  }
  int getInVectorLoc() {
     return(InVectorLoc);
  }
  SEXP OtherAsSexp() {
     SEXP EnvAcpp = R_NilValue;   int tProtect = 0;
     Rf_protect(EnvAcpp = R_FindNamespace(Rf_mkString(PCKGNAME)));  tProtect++;
     EnvironS = Rf_findVarInFrame(EnvAcpp, Rf_install("ACPPSAVESPACE"));
     if (Rf_isNull(EnvironS)) {
       Rprintf(" -- OtherAsSexp(): Error, EnvironS is NULL!\n");
     }
     if (Rf_isReal(EnvironS) || Rf_isInteger(EnvironS)) {
       Rprintf(" -- OtherAsSexp(): This is not good it's real or integer!\n"); R_FlushConsole();
     }
     SEXP ActualAList = 
       Rf_protect(Rf_findVarInFrame( EnvironS, Rf_install(".AcppList")));
     tProtect++;
     if (Rf_isNull(ActualAList) || Rf_length(ActualAList)) {
       Rf_error("Acpp::OtherAsSexp, we pulled null ActulAList!\n");
     }
     if (Rf_length(ActualAList) <= InVectorLoc) {
       Rf_error("Acpp::OtherAsSexp, no InVectorLoc = %d but length AList = %d!\n",
         InVectorLoc, Rf_length(ActualAList));
     }
     return(VECTOR_ELT(ActualAList, InVectorLoc));
  }

};

#endif