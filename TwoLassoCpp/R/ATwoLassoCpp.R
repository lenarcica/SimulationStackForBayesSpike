################################################################################
##
##   ATwoLassoCpp.R
##
##  (c) Alan Lenarcic 2009-2019:
##      Work for Lab of Edoardo Airoldi, Harvard Statistics
##      Work completed in lab of William Valdar, UNC Genetics and used
##      to help make comparison code in Lenarcic and Valdar paper against
##      methods such as those implemented around Arg-max estimators like 2Lasso
##    
##      It was found that although GroupBayes->"BayesSpike" was competitive
##      and similar to leading Argmax penalties, that 2Lasso is similarly very
##      competitive and thus Gibbs sampler approaches do not necessary
##      reject the value of arg-max methods.
##
##   R wrapper functions for TwoLassoCpp Package
##
##   This contains TwoLassoCpp() R function which is the main function
##    for processing calling TwoLasso functions in 2011+ version of package.
##   This shows Rcpp Modules interface to access object fields from R prompt.
##  
##   .onLoad() is the on package "library(TwoLassoCpp)" declarations.
##   .onAttach() is for re-attach purposes.
##
##   There are many default parameters that need to be set 
##    (or can be taken from current Environment)
##
##   As a warning, my  GetG0Text(), SetGText() functions are functions
##    designed to write and lock functions and variables both to globalenv()
##    and TWOLASSONAMESPACE environments.  This uses an "eval(parse(text=))"
##    interface to write code and lock it inside a function.
##
##
##      Alan Lenarcic 10/30/2013

## LICENSE INFO: R CODE
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/
#
#  Note, TwoSimR5 code is predominantly built around running existing 
#  selector impleentations which exist as R packages, most of which have
#  a GNU license.                                    

TWOLASSONAMESPACE <- environment()

.onLoad <- function(libname, pkgname) {
## load the module and store it in our namespace
  print("onLoad: initializing"); flush.console();
  if (!exists("flush.console")) {
    flush.console <- function() { return(1); }
  }
  print("Initallizing package TwoLassoCpp"); flush.console();
  if (!exists("TWOLASSONAMESPACE")) {
    try(TWOLASSONAMESPACE <- environment(), silent=TRUE);
    try(eval(parse(text=SetGText("TWOLASSONAMESPACE", S=1))));
    try(eval(parse(text=SetGText("TWOLASSONAMESPACE", envir="TWOLASSONAMESPACE", S=1))));
  } else {
    try(eval(parse(text=SetGText("TWOLASSONAMESPACE", S=1))));
    try(eval(parse(text=SetGText("TWOLASSONAMESPACE", envir="TWOLASSONAMESPACE", S=1))));  
  }
  try(library(Rcpp, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
  print("Declaring Globals"); flush.console();
   DeclareGlobals();
  
  print("Now Running the loading of onLoadAcpp"); flush.console();
  .onLoadAcpp(libname, pkgname);
  print("Running Second OnLoadFunction!");
   SecondLoadfunction();
  TryO = NULL;
  try(eval(parse(text=GetG0Text("modTwoLassoSexpCL", "globalenv()", S=1)))); 
  try(eval(parse(text=GetG0Text("modTwoLassoSexpCL", "TWOLASSONAMESPACE", S=1)))); 
  try( TryO <- bindingIsActive("modTwoLassoSexpCL", TWOLASSONAMESPACE), silent = TRUE);
  if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
    try(unlockBinding( "modTwoLassoSexpCL" , TWOLASSONAMESPACE ), silent=TRUE);
  }
  try(modTwoLassoSexpCL <- Module( "modTwoLassoSexpCL" ) );
  ##try(assign( "modTwoLassoSexpCL", Module( "modTwoLassoSexpCL" ), 
  ##  TWOLASSONAMESPACE ), silent=TRUE);
  try(eval(parse(text=GetG0Text("TwoLassoSexp", "globalenv()", S=1))), silent=TRUE);
  try(TwoLassoSexp <- modTwoLassoSexpCL$TwoLassoSexp);
  try(eval(parse(text=LockGText("TwoLassoSexp", "globalenv()", S=1))), silent=TRUE);
  try(eval(parse(text=LockGText("TwoLassoSexp", "TWOLASSONAMESPACE", S=1))), silent=TRUE);
  try(lockBinding( "modTwoLassoSexpCL", TWOLASSONAMESPACE ), silent=TRUE);
  try(eval(parse(text=LockGText("modTwoLassoSexpCL", "globalenv()", S=1))), silent=TRUE);
  try(eval(parse(text=LockGText("modTwoLassoSexpCL", "TWOLASSONAMESPACE", S=1))), silent=TRUE);


  print("Completed TwoLassoCpp OnLoad"); flush.console();
  ## FreeTBSRoo deletes TBSRoo from locked background environment when
  ##  the RCpp object for the regressions is deleted.
  ##  TBSRoo will still not be deleted if user keeps pointers to TBSRoo in R
  ##  environment. 
  if (!exists(".Two2LassoCppMiniValidate")) {
    try(eval(parse(text=GetG0Text(".Two2LassoCppMiniValidate", "TWOLASSONAMESPACE", S=1))));
  }
  try(eval(parse(text=SetGText(".Two2LassoCppMiniValidate", "globalenv", S=1))));

  try(eval(parse(text=GetG0Text("TwoLassoStayFitFunction", "globalenv()", S=1))), silent=TRUE);
  TwoLassoStayFitFunction <- function(X,Y,DoLogit=FALSE, SampleInput,StartBeta=-1,
    Verbose = 0, RFastestCFoldValidate2Lasso = .FastestCFoldValidate2Lasso, StandardFlag=FALSE,
    MaximumAllocation = 4096, TestX = NULL, ...) {
    if (!exists("RFastestCFoldValidate2Lasso") || 
      is.null(RFastestCFoldValidate2Lasso) ) {
      RFastestCFoldValidate2Lasso = Verbose;
      eval(parse(text=SetGText("RFastestCFoldValidate2Lasso", S=1)));    
    }  else if (!exists(".FastestCFoldValidate2Lasso") || 
      is.null(.FastestCFoldValidate2Lasso) ) {
      .RFastestCFoldValidate2Lasso = RFastestCFoldValidate2Lasso;    
      eval(parse(text=SetGText(".FastestCFoldValidate2Lasso", S=1)));
    }
    if (!is.null(RFastestCFoldValidate2Lasso) && RFastestCFoldValidate2Lasso > 2) {
      print(paste("RFastestCFoldValidate2Lasso: TwoLassoStayFitRun SI: (",
        paste(SampleInput, collapse=", "),")", sep="")); flush.console();
    }
    eval(parse(text=GetG0Text(".M2Stuff", "globalenv()", S=1)));
    LambdaAKInputs = .M2Stuff$LambdaAKInputs;
    LambdaDKInputs = .M2Stuff$LambdaDKInputs;   
    LambdaAK = .M2Stuff$LambdaASeq;
    LambdaDK = .M2Stuff$LambdaDSeq;
    
    XMean <- colMeans(X);
    X <- t(t(X)-XMean); 
    if (!exists("DoLogit")) {
        if (all(Y %in% c(0,1))) {  DoLogit = TRUE;
        } else {
          DoLogit = FALSE;
        }
    } else if (any(!(Y %in% c(0,1)))) {
      if (DoLogit == TRUE) {
        print("ATwoLassoCpp.R:TwoLassoStayFitFunction() Why is DoLogit = TRUE when there are non 0,1's in Y?");
        flush.console(); DoLogit=FALSE;
      } else {
        DoLogit = FALSE;
      }
    } else {
      if (DoLogit == FALSE) {
        print("ATwoLassoCpp.R:TwoLassoStayFitFunction(): if DoLogit = FALSE why are all Y's in 0,1?"); flush.console();
      } 
    }   
   if (is.null(LambdaDKInputs) || length(LambdaDKInputs) <= 0) {
     print("TwoLassoStayFitFunction;  Uh oh, LambdaDKInputs is Null.");
     flush.console();
   }
   if (is.null(LambdaAKInputs) || length(LambdaAKInputs) <= 0) {
     print("TwoLassoStayFitFunction;  Uh oh, LambdaAKInputs is Null.");
     flush.console();
   }
   if (Verbose > 2) {
     print("Launching Two2LassoMiniValidate")
   }
   try(eval(parse(text=GetG0Text(".HiddenTwoLassoOb", "globalenv()"))), silent=TRUE);  
   .HiddenTwoLassoOb <- TwoLassoCpp:::.Two2LassoCppMiniValidate(X = X, Y = Y, 
      DoLogit=DoLogit,
      PiAVectorInputs = SampleInputs[,1], 
      SigmaVectorInputs = SampleInputs[,2], 
      LambdaDK = LambdaDK, LambdaAK = LambdaAK, 
      LambdaDKInputs = LambdaDKInputs, LambdaAKInputs = LambdaAKInputs,
      OrderSeq = OrderSeq,  MaxCauchy = MaxCauchy * .5, 
      CauchyEpsilon = CauchyEpsilon * sqrt(p), 
      StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
      XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
      InitKKs = -5, WLSWeights = -1, TDFNu = -1,
      Verbose = RFastestCFoldValidate2Lasso-3, RecordFlag = 0,  Groupers = NULL,
      L2ShrinkagePrior = L2ShrinkagePrior, 
      RecordL2Shrinkage = RecordL2Shrinkage,StartBetaMatrix=StartBetaMatrix,
      RunFlag = TRUE,MaximumAllocation = MaximumAllocation, ...)
    try(eval(parse(text=LockGText(".HiddenTwoLassoOb", "globalenv()"))), silent=TRUE);
    try(eval(parse(text=SetGText(".HiddenTwoLassoOb", "TWOLASSONAMESPACE",S=1))), silent=TRUE);
    if (!is.null(TestX)) {
      TestX <- t(t(TestX) - XMean)
      FitY <- TestX %*% .HiddenTwoLassoOb$RecordBetaCVFinish
      if (DoLogit == TRUE) {
        FitY <- FitY + matrix(rep(.HiddenTwoLassoOb$TLS$RecordBeta0CV, each=NROW(TestX)), NROW(TestX), NCOL(FitY));
      }   
      return(list(FitBeta = .HiddenTwoLassoOb$RecordBetaCVFinish, FitY = FitY));
    }
    return(.HiddenTwoLassoOb$RecordBetaCVFinish);
  }
  try(eval(parse(text=LockGText("TwoLassoStayFitFunction", "globalenv()"))), silent=TRUE);  
  try(eval(parse(text=SetGText("TwoLassoStayFitFunction", "TWOLASSONAMESPACE", S=1))), silent=TRUE); 
}
.onAttach <- function(.Library, TwoLasso){
  ##print(".onAttach: Now Declare Globals."); flush.console();
  ## DeclareGlobals();
  ## SecondLoadfunction();
   print("After OnAttach Attaching package TwoLassoCpp, declare Globals"); flush.console();
   try(DeclareGlobals());
   print("OnAttach: Attempting to Second Load Function. "); flush.console();
   try(SecondLoadfunction());
  TWOLASSONAMESPACE <- environment()
  TryO = NULL;
  try( TryO <- bindingIsActive("modTwoLassoSexpCL", TWOLASSONAMESPACE), silent = TRUE);
  if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
   unlockBinding( "modTwoLassoSexpCL" , TWOLASSONAMESPACE )
  }
  try(assign( "modTwoLassoSexpCL", Module( "modTwoLassoSexpCL" ), 
    TWOLASSONAMESPACE ), silent=TRUE);
  try(lockBinding( "modTwoLassoSexpCL", TWOLASSONAMESPACE ), silent=TRUE);
  try(eval(parse(text=LockGText("modTwoLassoSexpCL", "globalenv()"))), silent=TRUE);
  try(eval(parse(text=LockGText("modTwoLassoSexpCL", "TWOLASSONAMESPACE"))), silent=TRUE);
  print("Completed TwoLassoCpp OnAttach"); flush.console();
  ## FreeTBSRoo deletes TBSRoo from locked background environment when
  ##  the RCpp object for the regressions is deleted.
  ##  TBSRoo will still not be deleted if user keeps pointers to TBSRoo in R
  ##  environment. 
    try(eval(parse(text=GetG0Text("TwoLassoStayFitFunction", "globalenv()", S=1))), silent=TRUE);
  TwoLassoStayFitFunction <- function(X,Y,DoLogit = FALSE, SampleInput,StartBeta=-1,
    Verbose = 0, RFastestCFoldValidate2Lasso = .FastestCFoldValidate2Lasso,
    StandardFlag=0,TestX = NULL, ...) {
    if (!exists("RFastestCFoldValidate2Lasso") || 
      is.null(RFastestCFoldValidate2Lasso) ) {
      RFastestCFoldValidate2Lasso = Verbose;
      eval(parse(text=SetGText("RFastestCFoldValidate2Lasso", S=1)));    
    }  else if (!exists(".FastestCFoldValidate2Lasso") || 
      is.null(.FastestCFoldValidate2Lasso) ) {
      .RFastestCFoldValidate2Lasso = RFastestCFoldValidate2Lasso;    
      eval(parse(text=SetGText(".FastestCFoldValidate2Lasso", S=1)));
    }
    if (!is.null(RFastestCFoldValidate2Lasso) && RFastestCFoldValidate2Lasso > 2) {
      print(paste("RFastestCFoldValidate2Lasso: TwoLassoStayFitRun SI: (",
        paste(SampleInput, collapse=", "),")", sep="")); flush.console();
    }
    if (!exists("DoLogit")) {
        if (all(Y %in% c(0,1))) {  DoLogit = TRUE;
        } else {
          DoLogit = FALSE;
        }
    } else if (any(!(Y %in% c(0,1)))) {
      if (DoLogit == TRUE) {
        print("ATwoLassoCpp.R:TwoLassoStayFitFunction() Why is DoLogit = TRUE when there are non 0,1's in Y?");
        flush.console(); DoLogit=FALSE;
      } else {
        DoLogit = FALSE;
      }
    } else {
      if (DoLogit == FALSE) {
        print("ATwoLassoCpp.R:TwoLassoStayFitFunction(): if DoLogit = FALSE why are all Y's in 0,1?"); flush.console();
      } 
    }
    XMean <- colMeans(X);
    X <- t(t(X) - XMean);
    eval(parse(text=GetG0Text(".M2Stuff", "globalenv()", S=1)));
    LambdaAKInputs = .M2Stuff$LambdaAKInputs;
    LambdaDKInputs = .M2Stuff$LambdaDKInputs;   
    LambdaAK = .M2Stuff$LambdaASeq;
    LambdaDK = .M2Stuff$LambdaDSeq; 
   
   if (is.null(LambdaDKInputs) || length(LambdaDKInputs) <= 0) {
     print("TwoLassoStayFitFunction;  Uh oh, LambdaDKInputs is Null.");
     flush.console();
   }
   if (is.null(LambdaAKInputs) || length(LambdaAKInputs) <= 0) {
     print("TwoLassoStayFitFunction;  Uh oh, LambdaAKInputs is Null.");
     flush.console();
   }
   if (Verbose > 2) {
     print("Launching Two2LassoMiniValidate")
   }
   try(eval(parse(text=GetG0Text(".HiddenTwoLassoOb", "globalenv()", S=1))), silent=TRUE);  
   if (!is.null(.HiddenTwoLassoOb)  && (!is.numeric(.HiddenTwoLassoOb) &&
     .HiddenTwoLassoOb[1] == 0)) {
     try(.HiddenTwoLassoOb$TCpp <- NULL);  
     gc();
     .HiddenTwoLassoOb <- NULL;
     eval(parse(text=SetGText(".HiddenTwoLassoOb", "globalenv()", S=1)));
     gc();
   }
   .HiddenTwoLassoOb <- TwoLassoCpp:::.Two2LassoCppMiniValidate(X = X, Y = Y, 
      DoLogit = DoLogit,
      PiAVectorInputs = SampleInputs[,1], 
      SigmaVectorInputs = SampleInputs[,2], 
      LambdaDK = LambdaDK, LambdaAK = LambdaAK, 
      LambdaDKInputs = LambdaDKInputs, LambdaAKInputs = LambdaAKInputs,
      OrderSeq = OrderSeq,  MaxCauchy = MaxCauchy * .5, 
      CauchyEpsilon = CauchyEpsilon * sqrt(p), 
      StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
      XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
      InitKKs = -5, WLSWeights = -1, TDFNu = -1,
      Verbose = RFastestCFoldValidate2Lasso-3, RecordFlag = 0,  Groupers = NULL,
      L2ShrinkagePrior = L2ShrinkagePrior, 
      RecordL2Shrinkage = RecordL2Shrinkage,StartBetaMatrix=StartBetaMatrix,
      RunFlag = TRUE,...)
    try(eval(parse(text=LockGText(".HiddenTwoLassoOb", "globalenv()"))), silent=TRUE);
    if (!is.null(TestX)) {
      TestX <- t(t(TestX) - XMean)
      FitY <- TestX %*% .HiddenTwoLassoOb$RecordBetaCVFinish
      if (DoLogit == TRUE) {
        FitY <- FitY + matrix(rep(.HiddenTwoLassoOb$TLS$RecordBeta0CV, each=NROW(TestX)), NROW(TestX), NCOL(FitY));
      }
      return(list(FitBeta = .HiddenTwoLassoOb$RecordBetaCVFinish, FitY = FitY));
    }    
    return(.HiddenTwoLassoOb$RecordBetaCVFinish);
  }
  try(eval(parse(text=LockGText("TwoLassoStayFitFunction", "globalenv()", S=1))), silent=TRUE);
  eval(parse(text=SetGText("TwoLassoStayFitFunction", "TWOLASSONAMESPACE", S=1)));
}



################################################################################
## SetGText() - Alan Lenarcic Package Helper function
##
##  Inside any function in R, SetGText constructs a set of R code that
##   tries to set an object with text name "AText" to environment 
##   "envir" (which is default Global) from the value in the current function.  
##    Thus the usage is 
##     eval(parse(text=SetGText("MyWantObject", "globanenv()", S=1)));
##
##  S = 1 is for silent.  AText is a character string for what is wanted.
SetGText <- function(AText, envir = "globalenv()", S=0) {
  if (S==0) {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    try(assign( \"", AText, "\", ", AText, ", ", envir, " ), silent=TRUE);
    ", sep=""));
  } else {
   return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    try(assign( \"", AText, "\", ", AText, ", ", envir, " ), silent=TRUE);
    ", sep="")); 
  }
}

################################################################################
## GetG0Text() - Alan Lenarcic Package Helper function
##
##  Inside any function in R, GetG0Text constructs a set of R code that
##   tries to pull and unlock an object with text name "AText" from environment 
##   "envir" (which is default Global) and then give it to the environment of
##   this function.  Thus the usage is 
##     eval(parse(text=GetG0Text("MyWantObject", "globanenv()", S=1)));
##
##  S = 1 is for silent.  AText is a character string for what is wanted.
GetG0Text <- function(AText, envir = "globalenv()", S=1) {
  if (S== 0) {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));
  } else {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));   
  }
}

################################################################################
##  LockGTest() - Alan Lenarcic Package Helper function
##
##  Inside any function in R, LockGTest constructs a set  of R code that
##   tries to set  and lock an object with text name "AText" to environment 
##   "envir" (which is default Global) from the value in the current function.  
##    Thus the usage is 
##     eval(parse(text=LockGTest("MyWantObject", "globanenv()", S=1)));
##
##  S = 1 is for silent.  AText is a character string for what is wanted.
LockGText <- function(AText, envir="globalenv()", S=0)  {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=TRUE); 
    try(assign( \"", AText, "\", ", AText, ", ", envir, " ), silent=TRUE);
    try(lockBinding( \"", AText, "\", ", envir, " ), silent=TRUE);
    ", sep="");
}

TwoLassoSetupDefaults <- function() {
  return("
  if (!exists(\"XTX\")) { XTX = NULL;}
  if (!exists(\"pCoefs\")) { pCoefs = -1;}
  if (!exists(\"WLSWeights\")) { WLSWeights = NULL; }
  if (!exists(\"nSample\")) { nSample = -1.0; }
  if (!exists(\"Verbose\")) { Verbose = 0; }
  if (!exists(\"XtY\")) { XtY = NULL; }
  if (!exists(\"InverseGammaConstant\")) { InverseGammaConstant = 1.0; }
  if (!exists(\"StartBetas\")) { StartBetas = NULL; }
  if (!exists(\"FixKa\")) { FixKa = -1; }
  if (!exists(\"tt2\")) { tt2 = 0;}
  if (!exists(\"ReturnBetas\")) { ReturnBetas = NULL; }
  if (!exists(\"dfTNoise\")) { dfTNoise = 0.0; }
  if (!exists(\"StandardFlag\")) { StandardFlag = 0; }
  if (!exists(\"tt1\")) { tt1 = 0; }
  if (!exists(\"TargetMinimumBeta\")) { TargetMinimumBeta <- 1.0; }
  
  if (!exists(\"DoLogit\") || is.null(DoLogit)) {
    if (all(Y %in% c(0,1))) { DoLogit = TRUE; } else {
      DoLogit = FALSE;
    }
  } else if (all((Y -min(Y,na.rm=TRUE)) %in% c(0,1))) {
    if (DoLogit == TRUE) {
       Y <- Y - min(Y, na.rm=TRUE);  
    } else {
      print(paste(\"ATwoLassoCpp.R:TwoLassoSetupDefaults: Why is Y very Logit like but DoLogit = ?\", DoLogit, sep=\"\")); flush.console();
    }
  } else if (any(!(Y[!is.na(Y)] %in% c(0,1)))) {
     if (DoLogit == TRUE) {
        print(\"ATwoLassoCpp.R::TwoLassoSetpDefaults: Hey, there are varied Y but DoLogit = TRUE? \"); flush.console();
        print(\"We Paste Bad Y to OnY \"); flush.console();
        OnY <- Y;
        eval(parse(text=SetGText(\"OnY\", \"globalenv()\", S=1)));
        print(\"POSSIBLE ERROR POSSIBLE ERROR POSSIBLE ERROR POSSIBLE ERROR POSSIBLE ERROR POSSIBLE ERROR POSSIBLE ERROR \");
        print(\"POSSIBLE ERROR: Setting DoLogit FALSE \"); flush.console();
        1 <- 2;
        DoLogit = FALSE;
     } 
  }  else if (all(Y[!is.na(Y)] %in% c(0,1))) {
    DoLogit <- TRUE;
  }
  if (!exists(\"LambdaDK\") || is.null(LambdaDK)) {
    LambdaDK <- c(1,2,3);
  }
  if (!exists(\"LambdaAK\") || is.null(LambdaAK)) {
    LambdaAK <- c(1,.5,1/3);
  }
  if (!exists(\"OrderSeq\") || is.null(OrderSeq)) {
    OrderSeq <- c(4,4,4);
  }
  if (!exists(\"OnLambda2\")) { OnLambda2 <- 0.0; }
  if (!exists(\"Lambda2\")) { Lambda2 = -1; }
  if ((!is.numeric(OnLambda2) || OnLambda2 == 0.0) && (is.numeric(Lambda2) && Lambda2[1] > 0)) {
    OnLambda2 = Lambda2[1];
  } else if (!is.numeric(OnLambda2)) {
    print(\"TwoLassoSetupDefaults: Defective OnLambda2 was provided. Set to 0.0; \"); flush.console();
    OnLambda2 = 0.0;
  }
  
  if (!exists(\"SigmaSq\") || is.null(SigmaSq) || is.na(SigmaSq) || is.nan(SigmaSq)) { 
    print(\"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\"); flush.console();
    print(\"Hey you supplied non existant SigmaSq or NULL SigmaSq!\"); flush.console();
    SigmaSq = 1.0; }
  if (!exists(\"DoLogit\")) { DoLogit = FALSE; }
  if (!exists(\"PiA\")) { PiA = .5; }
  if (!exists(\"IndexFirstRandomEffect\")) { IndexFirstRandomEffect = NULL; }
  if (!exists(\"tauEndList\")) { tauEndList = NULL; }
  if (!exists(\"CauchyEpsilon\")) { CauchyEpislon = .00001; }
  if (!exists(\"MaxCauchy\")) { MaxCauchy = 200; }
  if (!exists(\"HoldOn\")) { HoldOn = TRUE; }
  if (!exists(\"SigmaPrior\")) { SigmaPrior = NULL; }
  if (!exists(\"PiAPrior\")) { PiAPrior = NULL; }
  if (!exists(\"mVec\")) { mVec = NULL; }
  if (!exists(\"PiAVector\")) { PiAVector = NULL; }
  if (!exists(\"SigmaVector\")) { SigmaVector = NULL; }
  if (!exists(\"TDFNu\")) { TDFNu = 0; }
  if (!exists(\"RecordSigma\")) { RecordSigma = FALSE; }
  if (!exists(\"SigmaVectorInputs\")) { SigmaVectorInputs = NULL; }
  if (!exists(\"PiAVectorInputs\")) { PiAVectorInputs = NULL; }
  if (!exists(\"RunFlag\")) { RunFlag = -1; }
  if (!exists(\"RecordBetaCVFinish\")) { RecordBetaCVFinish = NULL; }
  if (!exists(\"OnGammas\")) { OnGammas = NULL; }
  if (!exists(\"L2ShrinkageRecords\")) { L2ShrinkageRecords = NULL; }
  if (!exists(\"L2ShrinkagePrior\")) { L2ShrinkagePrior = NULL; } 
  if (!exists(\"ForceXTXFlag\")) { ForceXTXFlag = -1; } 
  if (!exists(\"RecordFlag\")) { RecordFlag = 0; }    
  if (!exists(\"LogitNoise\")) { LogitNoise = 1; }
  if (!exists(\"LambdaAKVectorInputs\")) {
    if (exists(\"LambdaAKInputs\")) { LambdaAKVectorInputs = LambdaAKInputs; }
    LambdaAKVectorInputs = NULL;
  } 
  if (!exists(\"MaximumAllocation\")) { MaximumAllocation = 4096; }
  if (!exists(\"CauchyEpsilon\")) { CauchyEpsilon = .000001; } 
  if (!exists(\"MaxCauchy\")) { MaxCauchy = 200; } 
  if (!exists(\"CVQuitTime\")) { CVQuitTime = -1; }
       
  if (!exists(\"LambdaDKVectorInputs\")) {
    if (exists(\"LambdaDKInputs\")) { LambdaDKVectorInputs = LambdaDKInputs; }
    LambdaDKVectorInputs = NULL;
  }  
  ##eval(parse(text=TwoLassoCpp:::MyExistsF(\"L2ShrinkageRecords\", \"NULL\")));
  ##eval(parse(text=TwoLassoCpp:::MyExistsF(\"L2ShrinkagePrior\", \"NULL\")));
  
  if  (is.numeric(MaxCauchy) && is.numeric(CauchyEpsilon) &&
    (MaxCauchy > 0 || CauchyEpsilon > 0)) {
    if (MaxCauchy <= 0) {
      MaxCauchy = 50;
    }
    if (CauchyEpsilon < 0) {
      CauchyEpsilon = .001;
    }
  } else {
    MaxCauchy = 50; CauchyEpsilon = .0001;   
  }
  if (is.null(PiAPrior)) { PiAPrior = mVec; }
  GroupSexp = NULL; TCSPointer = NULL;  
  if (!exists(\"InitKKs\") || InitKKs <= 0) { InitKKs = 5; }
  if (!exists(\"InitKKs\") || is.null(InitKKs) ||  InitKKs <= 0) {
    print(paste(\"TwoLassoCpp: Invalid InitKKs = \", InitKKs, \" supplied, use default InitKKs = 5;\",
      sep=\"\"));  InitKKs = 5;
  }
  
  if (!exists(\"tt1\") || is.null(tt1) || !is.numeric(tt1) || tt1 <= 0) { tt1 = 1; }
  if (!exists(\"tt2\") || is.null(tt2) || !is.numeric(tt2) || tt2 <= 0) { tt2 = 1; }
   if ((!exists(\"PiA\") || is.null(PiA)) && exists(\"PiAVectorInputs\") && !is.null(PiAVectorInputs) &&
     length(PiAVectorInputs) >= 1) {
    PiA = PiAVectorInputs[1];      
  }

  if (!exists(\"OnLambda2\")) { OnLambda2 <- 0.0; }
  if (!exists(\"Lambda2\")) { Lambda2 = -1; }
  if ((!is.numeric(OnLambda2) || OnLambda2 == 0.0) && (is.numeric(Lambda2) && Lambda2[1] > 0)) {
    OnLambda2 = Lambda2[1];
  } else if (!is.numeric(OnLambda2)) {
    print(\"TwoLassoSetupDefaults: Defective OnLambda2 was provided. Set to 0.0; \"); flush.console();
    OnLambda2 = 0.0;
  }
  if (exists(\"DoLogit\") && DoLogit == TRUE) {
    if (OnLambda2 == 0.0) {
      OnLambda2 <- .001;
    } else if (OnLambda2) {
      print(paste(\"TwoLassoSetupDefaults: Warning, DoLogit is true but you set OnLambda2 = \", OnLambda2, sep=\"\")); flush.console();
    }
  }
  if (OnLambda2 < 0.0) {
    OnLambda2 <- 0.0;
  }
  ");
}
TwoLassoFormatX <- function() {
  return("
  if (!is.null(X) && !is.null(dim(X))) {
    pCoefs = dim(X)[2]  
  } else if (!is.null(XtX) && !is.null(dim(XtX))) {
    pCoefs = dim(XtX)[2]  
  } else {
    print(\"UnSatisfactory X or XtX given\");
    return(-1);
  }

  if(!is.loaded(\"LarsConvergencyShell\") ) {
    LoadTwoCC();
  } 
  if (Verbose >= 1) {
    print(paste(\"LambdaDK is \", paste(LambdaDK, collapse=\", \"), sep=\"\"));
    flush.console();
  }
  if (is.null(LambdaAK) || LambdaAK[1] == -999) {
    if (StLambdaA < 0 & PiA > 0 && PiA < 1.0) {
       StLambdaA =  .6 * log((1-PiA)/PiA);      
    } else if (StLambdaA < 0 && FixKa > 0) {
      pDtry <- FixKa / pCoefs; 
      StLambdaA <-  .6 * log((1-pDtry)/pDtry);
    } else if (StLambdaA < 0) {
      print(\"Why did you not give us StLambdaA\");
    }
  } 
  if (LambdaDK[1] == -999) {
    if (StLambdaD < 0 & PiA > 0 && PiA < 1.0) {
      StLambdaD =  .6 * log((1-PiA)/PiA);      
    } else if (StLambdaD < 0 && FixKa > 0) {
      pDtry <- FixKa / pCoefs; 
      StLambdaD =  .6 * log((1-pDtry)/pDtry);
    } else if (StLambdaD < 0) {
      print(\"Why did you not give us StlambdaD\");
    }
  }
  if (Verbose >= 2) {
    print(paste(\"Well we have LambdaAK = (\", paste(LambdaAK, collapse=\", \"),
      \")\", sep=\"\"));
    print(paste(\"         and LambdaDK = (\", paste(LambdaDK, collapse=\", \"), 
      \")\", sep=\"\"));   
    flush.console(); 
  }
  
  ##############################################################################
  ### Flag whether to attempt to Logit analysis
  ##   We will have to analyze Y vector to confirm that Logit is 
  ##     an appropriate analysis
  ##
  sDoLogit = 0;
  if (is.logical(DoLogit) && DoLogit == FALSE) {
    sDoLogit = 0;
  }  else if (is.logical(DoLogit) && DoLogit == TRUE) { 
    sDoLogit = 1; 
  } else if (is.numeric(DoLogit) && DoLogit[1] >= 1) {
    sDoLogit = 1;
  } else if (length(DoLogit) == 1) {
      sDoLogit = as.integer(DoLogit);
  } else {
        sDoLogit = 0;
  }
  if (sDoLogit == 1 && all(Y %in% c(0,1))) {
    ##if (!all(X[,1] == 1.0)) {
    ##  X <- cbind(rep(1, NROW(X)), X);
    ##  NoShrinkColumns = 1;
    ##} else { NoShrinkColumns = NULL;
    ##}
    NoShrinkColumns = NULL;
  }  else {
    NoShrinkColumns = NULL;
  }
  if (!exists(\"SigmaSq\") || is.null(SigmaSq) || 
    !is.numeric(SigmaSq) || is.na(SigmaSq[1]) || is.nan(SigmaSq[1]) || SigmaSq[1] < 0) { 
    if (!is.null(SigmaVector) && is.numeric(SigmaVector) && SigmaVector[1] > 0) {
      if (Verbose >= 3) {
        print(paste(\"TwoLassoCpp: Getting SigmaSq from SigmaVector = \", SigmaVector[1], sep=\"\"));
        flush.console();
      }
      SigmaSq = SigmaVector[1];
    } else {
      if (Verbose >= 3) {
        print(paste(\"TwoLassoCpp: Setting SigmaSq to 1 because of nothing better.\")); flush.console();
      }
      SigmaSq = 1.0; 
    }
  }  else {
    if (Verbose >= 3) {
      print(paste(\"TwoLassoCpp: You supplied a manual SigmaSq = \", SigmaSq, sep=\"\"));
      flush.console();
    }
  }  
  if (is.numeric(dfTNoise) && dfTNoise > 0) { 
     TDFNu = dfTNoise; 
  } else if (is.numeric(TDFNu) && TDFNu > 0) {
    dfTnoise <- TDFNu
  } else {
    TDFNu <- -1;  dfTNoise = -1;
  }
  ");
}

TwoLassoSetupLambda <- function() {
  return("
    if (!is.numeric(PiA) || PiA[1] == -1) {
    if (is.null(PiAVector)) { PiA = PiAVector[1]; } else {
      PiA = .5;
    }
  }
  if (any(PiA >= 1)) { PiA[PiA >= 1] <- .99; }
  if (any(PiA <= 0.0)) { PiA[PiA <= 0.0] <- 0.01; }
  if (any(is.na(PiA))) { PiA[is.na(PiA)] <- .5; }

  if (is.vector(LambdaDK) == FALSE || length(LambdaDK) == 1 
    || LambdaDK == -999) {
    LambdaDK = StLambdaDK * LambdaDMultiplier^(0:(TotalRuns-1));
    LambdaAK = StLambdaDK * LambdaAMultiplier^(0:(TotalRuns-1));
    LambdaAStop = round(LambdaAStop);
    LambdaAK[(1:TotalRuns) > LambdaAStop] = LambdaAK[LambdaAStop]
  } else {
    TotalRuns = min(length(LambdaDK), length(LambdaAK));
    LambdaDK= LambdaDK[1:TotalRuns];  LambdaAK = LambdaAK[1:TotalRuns];
  }
  if (is.vector(OrderSeq) == TRUE && length(OrderSeq) >= TotalRuns) {
    OrderSeq = OrderSeq[1:TotalRuns];
  } else if (length(OrderSeq) == 1 && OrderSeq == -999) {
    OrderSeq <- rep(8, TotalRuns);
  }
  
  LambdaAList <- list(LambdaAK=LambdaAK, LambdaDK = LambdaDK, OrderSeq=OrderSeq);
  ##try(TLS$SetupLambda(LambdaAK, LambdaDK, OrderSeq)); 
  ##if (Verbose >= 1) {
  ##  print(paste(\"TLS$LambdaDK is \", paste(TLS$LambdaDK, collapse=\", \"), sep=\"\"));
  ##  flush.console();
  ##}
  OnGamma1 = .5*(PiA * LambdaAK[1] + (1.0 - PiA) * LambdaDK[1]);
  OnGamma1 = OnGamma1[1];
  if (is.na(OnGamma1) || is.nan(OnGamma1)) {
    print(paste(\"TwoLassoCpp: Error nan OnGamma1, PiA = (\", paste(PiA, collapse=\",\"),
      \") and LambdaAK = (\", paste(round(LambdaAK,4), collapse = \", \"),
      \") and LambdaDK = (\", paste(round(LambdaDK, 4), collapse=\",\"), \")\", sep=\"\"));
    flush.console();
  }
  ");
}
TwoLassoFormatXBeta <- function() {
  return("
    if (!exists(\"StartBeta\") || is.null(StartBeta)) {
    OnBeta = rep(0, pCoefs);  StartBeta = rep(0, pCoefs);
  } else {
    OnBeta = StartBeta + 0.0;
  }

  ");
}
TwoLassoSetupOne <- function() {
   return("
    if (is.matrix(XtX) == FALSE || length(XtX) == 1 || XtX == -1) {
      print(\"TwoLassoCpp, Error no proper X given, XtX False\"); return;  
    }
    if (dfTNoise > 0) {
      print(\"TwoLassoCpp: Error, Without X, we can't use dfTNoise\"); return;
    }
    if ( (is.matrix(XtY) == FALSE && is.vector(XtY) == FALSE) ||
      length(XtY) == 1 || XtY == -1) {
      print(\"TwoLasso, Error no proper X given, XtY False\"); return;
    }
    if (length(nSample) > 1 || nSample == -1) {
      print(\"TwoLasso, Error, no proper NLen given\");
      return(-1);
    }

    XTXFlag = 1;
    XtY = as.double(XtY);   p = length(XtY);
    TnSample = -nSample; pCoefs = length(XtX[1,]);
    WLSWeights = -1;  dfTNoise = -1;
    sdX = sqrt(diag(XtX)); sdY = NULL;
    ListFlags = as.numeric(c(Verbose[1], TnSample[1], PiA[1], FixKa[1], OnGamma1[1], SigmaSq[1], InitKKs[1], RecordFlag,
      CauchyEpsilon,  MaxCauchy, 0,0));
    try(names(ListFlags) <- c(\"Verbose\", \"TnSample\", \"PiA\", \"FixKa\", \"OnGamma1\", \"SigmaSq\", \"InitKKs\",
        \"RecordFlag\", \"CauchyEpsilon\", \"MaxCauchy\", \"sDoLogit\", \"ForceXTXFlag\")); 
    eval(parse(text=SetGText(\"ListFlags\", \"globalenv()\", S=1)));
    if (Verbose >= 1) {
      print(\"TwoLassoCpp: Now Create TLS: TwoLassoSexp from XtY, XtX. \"); flush.console();
    }
    ");
}
TwoLassoSetupCppTwo <- function() {
    return("
    XTXFlag = 0;
      XTXFlag = 1;
      XtX = as.double( t(X)%*% X );
      XtY = as.double( t(X) %*% Y );  p = length(XtY);
      TnSample = -length(Y); pCoefs = length(X[1,]); 
      sdX = apply(X,2,sd); sdY = sd(as.vector(Y));  


      SigmaSq = as.numeric(SigmaSq);
      if (is.na(SigmaSq) || is.nan(SigmaSq)) {
        SigmaSq = 1.0
      }
      ListFlags = as.numeric(c(Verbose=Verbose, TnSample=TnSampple, 
        PiA=PiA[1], FixKa=FixKa, OnGamma1=OnGamma1, SigmaSq=SigmaSq, 
        InitKKs=InitKKs, 
        RecordFlag=RecordFlag, CauchyEpsilon=CauchyEpsilon, 
        MaxCauchy=MaxCauchy, sDoLogit=sDoLogit,0));
      eval(parse(text=SetGText(\"ListFlags\", \"globalenv()\", S=1)));
     ## eval(parse(text=SetGText(\"ListFlags\", S=1)));
     ## print(paste(\"ListFlags is length: \", length(ListFlags), \" and (\",
     ##   paste(ListFlags, collapse=\", \"), sep=\"\")); flush.console();
     ##   print(paste(\"RecordFlag = \", RecordFlag, sep=\"\"));
     ##   print(paste(\"MaxCauchy = \", MaxCauchy, sep=\"\"));
     ##   print(paste(\" sDoLogit = \", sDoLogit, sep=\"\"));
      if (is.null(XtX)) {
        tryCatch(\"TwoLassoCpp: oops, well XtX is NULL in Y,X to XtX-XtY version\");
      }
      if (is.null(XtY)) {
        tryCatch(\"TwoLassoCpp: oops, well XtY is NULL in Y,X XtX to XtY version\");
      }
    ");
}
TwoLassoSetupCppThree <- function() {
    return("
      TnSample = length(Y); pCoefs = length(X[1,]); p = pCoefs;
      if (length(WLSWeights) != length(Y)) {
        WLSWeights = -1;
      }
      if (pCoefs > 10000) {
        InitKKs = max(InitKKs, round(pCoefs *mean(PiA)));
        if (InitKKs > pCoefs) {InitKKs = pCoefs;}

      }
      if (InitKKs < 0) {
        InitKKs = min(50, length(X[1,]));
      }
      if (InitKKs > 1000) {
        InitKKs <- 1000;
      }
      if (length(dfTNoise) == 1 && dfTNoise > 0) {
        WLSWeights = rep(1, length(Y))
      } else {
        dfTNoise = -1;
      }
      sdX = apply(X,2,sd); sdY = sd(as.vector(Y));
      if (InitKKs <= 0) { InitKKs = 5;}
      ##if (InitKKs < pCoefs) { InitKKs = pCoefs; }
      print(paste(\"SETTING UP ListFlags, SigmaSq[1] = \", SigmaSq[1], sep=\"\")); flush.console();
      ListFlags = as.numeric(c(Verbose[1], TnSample[1], PiA[1], FixKa[1], OnGamma1[1], SigmaSq[1], InitKKs, 
        RecordFlag, CauchyEpsilon, MaxCauchy, sDoLogit, ForceXTXFlag));
      try(names(ListFlags) <- c(\"Verbose\", \"TnSample\", \"PiA\", \"FixKa\", \"OnGamma1\", \"SigmaSq\", \"InitKKs\",
        \"RecordFlag\", \"CauchyEpsilon\", \"MaxCauchy\", \"sDoLogit\", \"ForceXTXFlag\"));
      eval(parse(text=SetGText(\"ListFlags\", \"globalenv()\", S=1)));
      if (is.null(Y)) {
        tryCatch(\"TwoLassoCpp: oops, well Y is NULL in Y,X version\");
      }
      if (is.null(X)) {
        tryCatch(\"TwoLassoCpp: oops, well X is NULL in Y,X version\");
      }
     if (Verbose >= 1) {
       print(\"TwoLassoCpp: CreateTwoLassoSexp from Y, X full version\"); flush.console();
     }
    ");
}
TLSSetupSigmaGammas <- function() {
 return("
  try(TLS$OnSigma <- SigmaSq)
  try(TLS$MaximumAllocation <- MaximumAllocation)

  if (length(PiA) >= 2) {
    try(TLS$OnRPiA <- PiA[2]);
  }
  if (length(PiAPrior) >= 4) {
    try(TLS$RPiAPrior <- PiAPrior[3:4])
  }
  if (exists(\"OnLambda2\") && is.numeric(OnLambda2) && OnLambda2[1] > 0.0) {
    try(TLS$OnLambda2 <- OnLambda2);
  }
  if (Verbose > 1) {
    print(paste(\"Setting Up Lambda with vectors length A:\",
      length(LambdaAK), \", D: \", length(LambdaDK), \", \", 
      \"O: \", length(OrderSeq), sep=\"\"));                                                                       
    flush.console();
  }
  ARText <- \"
  if (DoLogit == TRUE) {
    try(TLS$PrintFlagGLMB <- TLS$Verbose - 3);
  }
  \";
  try(eval(parse(text=ARText)));
  if (InverseGammaConstant != 1.0) {
    try(TLS$InverseGammaConstant <- InverseGammaConstant);
  }
  if (Verbose >= 1){
    print(paste(\"Before Setup LambdaDK is \", paste(LambdaDK, collapse=\", \"), sep=\"\"));
    flush.console();
  }

  if (!exists(\"tt1\")) { tt1 = 1; }  
  if (!exists(\"tt2\")) { tt2 = 1; }
  try(TLS$Setuptt(as.integer(tt1),as.integer(tt2)));
  TLS$tt1 <- 0;  TLS$tt2 <- 0;
  
  if (exists(\"NoShrinkColumns\") && !is.null(NoShrinkColumns) &&
    length(NoShrinkColumns) >= 1) {
    try(TLS$NoShrinkColumns <- NoShrinkColumns);      
  }
  
  if (TDFNu > 0) {
    try(TLS$SetupTDFNu(TDFNu, rep(1, length(Y))) );
  }
  if (!is.null(tauEndList) && length(tauEndList) >= 1 && IndexFirstRandomEffect >= 1) {
    ## subtract one for zero index
    try(TLS$SetupGroupTwoLasso(c(IndexFirstRandomEffect, tauEndList)-1));
  }

  


  if (!is.null(RecordSigma) && ((is.numeric(RecordSigma) && RecordSigma > 0) || RecordSigma != FALSE)) {
    if (Verbose > 1) {
      print(\"RecordSigma > 0, setting up RecordSigma in TLS\"); flush.console();
    }
    if (length(LambdaAK)<= 0) {
      print(paste(\"ERROR Setting up TLS$Sigma, length(LambdaAK) is \", length(LambdaAK), sep=\"\")); flush.console();  
    }
    if (!is.null(TLS$OnSigma) && length(TLS$OnSigma) >= 1  && !is.na(TLS$OnSigma) && is.numeric(TLS$OnSigma) && TLS$OnSigma >= 0 ) {
      TLS$Sigma = rep(TLS$OnSigma, length(LambdaAK)+1);
    } else {
      TLS$Sigma = rep(1.0, length(LambdaAK)+1);
    }
    if (DoLogit == 1 && LogitNoise != 1.0)  {
      try(TLS$Sigma <- rep(LogitNoise, length(LambdaAK)+1));
      try(TLS$LogitNoise <- LogitNoise);
    }
  }
  ");
}
TLSSetupVectorInputs <- function() {
  return("
  if (!exists(\"SigmaVector\")) {
    SigmaVector = NULL;
  } 
  if (!exists(\"CVQuitTime\") || CVQuitTime == -1) {
    try(eval(parse(text=GetG0Text(\"CVQuitTime\", \"globalenv()\", S=1))));
    if (is.null(CVQuitTime) || CVQuitTime == 0) {
      CVQuitTime = -1;
    } 
  }
  try(TLS$CVQuitTime <- CVQuitTime-1);
  if (!is.null(SigmaVector)) {
    if (length(SigmaVector) == length(LambdaAK)) {
      try(TLS$SigmaVector <- SigmaVector);
    } else {
      printf(paste(\"Error: Let SigmaVector have length \", length(LambdaAK), \",\",
        \"  not \", length(SigmaVector), \", continueing...\", sep=\"\"));
      flush.console();
    }
    try(TLS$Sigma <- as.numeric(SigmaVector[1]));
  } else {
    if (DoLogit == 1 && LogitNoise != 1.0) {
      try(TLS$Sigma <- as.numeric(LogitNoise));    
    } else {
      try(TLS$Sigma <- as.numeric(1.0));
    }
  }
  if (!exists(\"PiAVector\")) { PiAVector = NULL; }
  if (!is.null(PiAVector)) {
    if (length(PiAVector) == length(PiAVector)) {
      try(TLS$PiAVector <- PiAVector);
    } else {
      printf(paste(\"Error: Let PiAVector have length \", length(LambdaAK), \",\",
        \"  not \", length(PiAVector), \", continueing...\", sep=\"\"));
      flush.console();
    }
  }

  if (exists(\"TargetMinimumBeta\") && length(TargetMinimumBeta) == 1 &&
    TargetMinimumBeta != 1.0 && TargetMinimumBeta > 0.0) {
   try(TLS$TargetMinimumBeta <- TargetMinimumBeta);   
 }
  if ( (!is.null(SigmaVectorInputs)) || !is.null(PiAVectorInputs)) {
    if (is.null(SigmaVectorInputs) || length(SigmaVectorInputs) < PiAVectorInputs) {
      print(paste(\"  Please Let SigmaVector[\", length(SigmaVectorInputs), 
        \"] Inputs have same length as PiAVectorInputs[\", 
        length(PiAVectorInputs), \"]\", sep=\"\")); flush.console();
      SigmaVectorInputs = rep(SigmaSq[1], length(PiAVectorInputs));
    }
    if(is.null(PiAVectorInputs)  || length(PiAVectorInputs) < SigmaVectorInputs) {
      print(paste(\"  PiAShort!: Please Let SigmaVector[\", length(SigmaVectorInputs), 
        \"] Inputs have same length as PiAVectorInputs[\", 
        length(PiAVectorInputs), \"]\", sep=\"\")); flush.console();
        PiAVectorInputs = rep(PiA[1], length(SigmaVectorInputs));
    }
    if (length(PiAVectorInputs) != length(SigmaVectorInputs)) {
      print(paste(\"PiAVector/SigmaVector: Doh! after that, still short! \",
        length(PiAVectorInputs), \", \", length(SigmaVectorInputs), \".\", sep=\"\"));
      flush.console();
    }
    if (is.null(RecordBetaCVFinish) || length(RecordBetaCVFinish) !=
      p * length(SigmaVectorInputs) ) {
      RecordBetaCVFinish = matrix(0, p, length(SigmaVectorInputs));
    }
    if (length(LambdaAKVectorInputs) != length(LambdaAK) * length(SigmaVectorInputs)) {
      LambdaAKVectorInputs = matrix(rep(LambdaAK, length(SigmaVectorInputs)),
        length(LambdaAK), length(SigmaVectorInputs));
    }
    if (length(LambdaDKVectorInputs) != length(LambdaDK) * length(SigmaVectorInputs)) {
      LambdaDKVectorInputs = matrix(rep(LambdaDK, length(SigmaVectorInputs)),
        length(LambdaDK), length(SigmaVectorInputs));
    }
    try(TLS$SetupInputs(PiAVectorInputs,  SigmaVectorInputs,
      LambdaAKVectorInputs,  LambdaDKVectorInputs,
      StartBetaMatrix, RecordBetaCVFinish))
  } else {
    RecordBetaCVFinish = NULL;
  }
  ");

}

###############################################################################
## TLSSetupRecordings()
##
##   Setup Recordining Principles in TLS
TLSSetupRecordings <- function() {
  return("
      RecBBOn1 = TLS$sRecBBOn1;  
  RecOnBeta = TLS$sRecOnBeta;
     
  if ((is.logical(RecordFlag) && RecordFlag[1] == TRUE) || RecordFlag  == 1) {
    RecBBOn1 = matrix(0, p, length(LambdaAK));
    RecOnBeta = matrix(0, p, length(LambdaAK));
    try(TLS$SetupRecords(RecBBOn1, RecOnBeta));   
  } else if ( RecordFlag >=2 ) {
    RecBBOn1 = array(0, dim=c(p, length(LambdaAK), max(OrderSeq)));
    RecOnBeta = array(0, dim=c(p, length(LambdaAK), max(OrderSeq)));
    try(TLS$SetupRecords(RecBBOn1, RecOnBeta));   
  } else { RecordFlag = 0; RecBBOn = NULL;  RecOnBeta = NULL; }
  if (!is.null(PiAPrior)) {
    try(TLS$PiAPrior <- PiAPrior[1:2]);
  }
  if (DoLogit == FALSE && !is.null(SigmaPrior)) {
    try(TLS$SigmaPrior <- SigmaPrior);
  }
  try(TLS$tt1 <- 0);  try(TLS$tt2 <- 0);
  ");
}

TwoLassoTwoLassoSexpNew <- function() {
  return("
    if ((is.null(X) || is.matrix(X) == FALSE || length(X) == 1 || X == -1) && DoLogit == 0) {
    eval(parse(text=TwoLassoSetupOne()));
    TLS <- new(TwoLassoSexp, ListFlags, XtY, XtX, OnBeta, LambdaAList); 
    if (Verbose > 0) {
      print(paste(\"Two2Lasso: We're going to supply XtX, XtY with TnSample = \", 
         TnSample, sep=\"\"));
      flush.console();
    }
  } else if (( is.vector(Y) == FALSE && is.matrix(Y) == FALSE) || 
    length(Y) == 1 || Y == -1) {
    print(\"TwoLasso, Error no proper Y given\"); 
    print(\" Y is : \"); print( Y ); return(-1);
  } else if (length(Y) > 1000 && NCOL(X) < length(Y) &&
      (length(WLSWeights) != length(Y) ||
        length(dfTNoise) == 1 && dfTNoise > 0  ) && DoLogit == 0) {
      eval(parse(text=TwoLassoSetupCppTwo()));
      TLS <- new(TwoLassoSexp, ListFlags, XtY, XtX, OnBeta, LambdaAList); 
      if (Verbose > 0) {
        print(paste(\"Two2Lasso: Y too long, We're going to supply XtX, XtY with TnSample = \", TnSample, sep=\"\"));
        flush.console();
      }
  }  else {
      eval(parse(text=TwoLassoSetupCppThree()));
      TLS <- new(TwoLassoSexp, ListFlags, Y, X, OnBeta, LambdaAList); 
      if (Verbose > 0) {
        print(paste(\"Two2Lasso: We're going to supply X, Y with TnSample = \",   
          TnSample, sep=\"\"));
        flush.console();
      }
  }      

  if (DoLogit == 1) {
    try(TLS$LogitNoise <- LogitNoise);
  }
  
  ")
}

################################################################################
##  TwoLassoCpp(X, Y, ...)
##
##  TwoLassoCpp is the key R wrapper function in 2011+ version of TwoLasso
##
##  This function creates a TwoLassoSexp object TLS (TwoLassoObject2011.h)
##   using Rcpp modules.  This module object has C functions to operate
##   on it which can be called from within R terminal.
##
##
TwoLassoCpp <- function(
     X = NULL, Y= NULL, XtX = NULL, XtY = NULL,
     DoLogit = FALSE,
     ReturnBetas = NULL, StartBeta = NULL,
     InverseGammaConstant= 1.0,
     dfTNoise = -1, StandardFlag = 0,
     CauchyEpsilon= .00001, MaxCauchy= 200,
     pCoefs = -1, nSample = -1, FixKa = -1,
     tt1 = NULL, tt2 = NULL, WLSWeights= NULL, 
     Verbose= 0, LambdaAK = c(1,.5,1/3), 
     LambdaDK= c(1,2,3), OrderSeq= c(4,4,4), TargetMinimumBeta = 1.0, PiA = .5, SigmaSq = 1,
     PiAVector= NULL, SigmaVector = NULL, PiAPrior = NULL,
     mVec= NULL, SigmaPrior= NULL,  InitKKs = 20, GroupSexp= NULL, 
     OnGammas = NULL, RecordFlag = 0, RecordBetaCVFinish = NULL,
     SigmaVectorInputs = NULL, PiAVectorInputs = NULL,
     LambdaAKVectorInputs= NULL, LambdaDKVectorInputs = NULL,
     L2ShrinkagePrior = NULL, L2ShrinkageRecords = NULL,  StartBetaMatrix = NULL,
     RunFlag = -1, RecordSigma = FALSE, TDFNu = 0, HoldOn = TRUE, 
     IndexFirstRandomEffect = NULL,
     tauEndList = NULL, ForceXTXFlag = 0, QuitBeforeRunReturnTLS = FALSE, CVQuitTime = -1, 
     MaximumAllocation = 4096, LogitNoise = 1, OnLambda2 = 0.0, Lambda2 = -1, ...) {
        
  eval(parse(text=TwoLassoSetupDefaults()));
  eval(parse(text=TwoLassoFormatX()));
  ## eval(parse(text=TwoLassoCpp:::MyExistsF("RunFlag", -1)));
                                   
  if (Verbose > 0) {
    print("TwoLassoCpp:: Start"); flush.console();
  }    

  if (Verbose >= 3) {
    try(print(paste("TwoLassoCpp: we have that Sigmasq on input is ", SigmaSq, sep="")));
    flush.console();
  }
  eval(parse(text=TwoLassoSetupLambda()));

  
  eval(parse(text=GetG0Text("TWOLASSONAMESPACE", envir="globalenv()", S=1)));
  if (is.environment(TWOLASSONAMESPACE)) {
    suppressWarnings(try(eval(parse(text=GetG0Text("TwoLassoSexp", envir="TWOLASSONAMESPACE", S=1)))));
  } else {
    suppressWarnings(eval(parse(text=GetG0Text("TwoLassoSexp", envir="globalenv()", S=0))));
  }
  ##TwoLassoSexp <- modTwoLassoSexpCL$TwoLassoSexp

  if (Verbose > 0) {
    print("TwoLassoCpp: Creating Rcpp C++ class"); flush.console();
  }  
  
  eval(parse(text=TwoLassoFormatXBeta()));
  
  ##Actual Allocation of TwoLassoSexp C++ Class.
  eval(parse(text=TwoLassoTwoLassoSexpNew()));
  
  eval(parse(text=TLSSetupSigmaGammas()));
  
  TLS$RefreshBBOn1();
  TLS$ReUpdateOnGammas();
  

  
  ##TLS$UpdateBBOn1();
  ##TLS$UpdateBeta();
  eval(parse(text=TLSSetupVectorInputs()));
  eval(parse(text=TLSSetupRecordings()));

  if (QuitBeforeRunReturnTLS == TRUE) {
    print("QuitBeforReunReturnTLS: Here you go");
    return(TLS);
  }
  ABZ <- -1;
  if (exists("RecordBetaCVFinish") && !is.null(RecordBetaCVFinish) && length(RecordBetaCVFinish) > 1) {
    if (Verbose >= 1 || TLS$p >= 50000) {
      print("ATwoLassoCpp.R::TwoLassoCpp():  Now Run TLS$RunCrossValidate Option"); flush.console();
      try(print(paste("  PiA vary from (", min(PiAVectorInputs), ", ", max(PiAVectorInputs), ")",
        " and Sig vary from (", min(SigmaVectorInputs), ", ", max(SigmaVectorInputs), ")", sep="")));
      flush.console();
    }
    CVTLS <- TLS;
    eval(parse(text=SetGText("CVTLS", "globalenv()", S=1)));
    try(ABZ <- TLS$RunCrossValidate());
    if (exists("NukeCrossValidate") && NukeCrossValidate == 1) {
      print("ATwoLassoCpp.R::TwoLassoCpp(): NukeCrossValidate is set, we are saving TLS$CV after RunCrossValidate to OnCVTLS");
      eval(parse(text=SetGText("CVTLS", "globalenv()", S=1)));
      return(CVTLS);
    }
    if (Verbose >= 1 || TLS$p >= 50000) {
      print("ATwoLassoCpp.R::TwoLassoCpp(): TLS$RunCrossValidate completed. "); flush.console();
    }
    if (ABZ == -1) {
      try(print(paste("ERROR: ATwoLassoCpp.R::TwoLassoCpp(), We ran CrossValidate and got ABZ failure: Memory=",
        TLS$MemoryFailFlag, sep=""))); flush.console();
      eval(parse(text=SetGText("CVTLS", "globalenv()", S=1)));
      return(-1);
    }
  } else if (!is.null(L2ShrinkagePrior)) {                                                 
    if (is.null(L2ShrinkageRecords) || length(L2ShrinkageRecords) < length(LambdaAK)) {
      L2ShrinkageRecords = rep(0, length(LambdaAK));
    }
    try(TLS$SetupL2Shrinkage(L2ShrinkagePrior, L2ShrinkageRecords));
    try(ABZ <- TLS$RunTwoLassoRegression(RecordFlag));   ## Run Algorithm
  } else {
    try(ABZ <- TLS$RunTwoLassoRegression(RecordFlag));   ## Run Algorithm
  }
  if (ABZ < 0) {
    print("TwoLassoCpp: TLS RuntwoLassoRegression ends in error, returning TLS! ");
    ErrorTLS <- TLS;
    eval(parse(text=SetGText("ErrorTLS", "globalenv()", S=1)));
    flush.console();
    try(TLS$FailureFlag <- 1);
    return(TLS);
  }

  if (StandardFlag == 1) {
    USOnBeta = TLS$Beta;
    Beta = TLS$Beta * sdX / sdY;
    if (!is.null(TLS$RecOnBeta)) {
      try(RecOnBeta <- TLS$RecOnBeta * sdX / sdY);
      try(USRecOnBeta <- TLS$RecOnBeta);
    } else {
      RecOnBeta = NULL; USRecOnBeta = NULL; RecBBOn1 = NULL;
    }
  } else {
    USOnBeta = NULL;
    Beta = TLS$Beta;
    if (!is.null(TLS$RecOnBeta)) {
      try(USRecOnBeta <- NULL);
      try(RecOnBeta <- TLS$RecOnBeta);
    } else {
      USRecOnBeta = NULL;
      RecOnBeta = NULL;
    }
  }
  BBOn1 <- NULL;
  try(BBOn1 <- TLS$BBOn1);
  if (exists("RecordBetaCVFinish") && !is.null(RecordBetaCVFinish) && length(RecordBetaCVFinish) > 1) {
  return(
    list(FinalBeta = Beta, USOnBeta = USOnBeta, RecOnBeta = RecOnBeta,
      RecBBOn1 = RecBBOn1,
      USRecOnBeta = USRecOnBeta, TLS = TLS, DidIFail = 0,BBOn1 = BBOn1,
      FailureFlag = 0,RecordBetaCVFinish = RecordBetaCVFinish,
      DoLogit=DoLogit
    )
  );  
  }
  return(
    list(FinalBeta = Beta, USOnBeta = USOnBeta, RecOnBeta = RecOnBeta,
      USRecOnBeta = USRecOnBeta, TLS = TLS, DidIFail = 0,BBOn1 = BBOn1,
      FailureFlag = 0, DoLogit=DoLogit
    )
  );
}

##############################################################################
##  ProptotypeRunTwoLassoRegression
##
##   A prototype for using a TLS (TwoLassoObject2011.cc) object
##     Alan Lenarcic 10/30/2013
##
ProptotypeRunTwoLassoRegression <- function(TLS, tt1 = 1) {
  TLS$tt1 <- tt1;
  TLS$OnSigma
  
  ##TLS$UpdateSigmaSq();
  TLS$ReUpdateOnGammas()
  
  OldBeta <- TLS$Beta;
  for (ii in 1:TLS$p) {  
    TLS$OnCoord = ii-1;
    TLS$UpdateCoord();
  }
  MyMove <- sum(abs(OldBeta - TLS$Beta));

  if (TLS$OrderSeq[TLS$tt1+1] >= 0) {
    TLS$UpdateBeta();

    TLS$UpdateBBOn1();
    if (TLS$SigmaBar > 0) {
      TLS$UpdateSigmaSq();
    }
    if (TLS$PiAPrior[1] > 0) {
      UpdateHatPiA();
    }
  }
  for (tt2 in 1:10) {
    TLS$UpdateBeta();

    TLS$UpdateBBOn1();
    if (TLS$SigmaBar > 0) {
      TLS$UpdateSigmaSq();
    }
    if (TLS$PiAPrior[1] > 0) {
      UpdateHatPiA();
    }        
  }
}

StartCFoldValidateText <- function() {
 return("
      eval(parse(text=MyExistsF(\"nPi\", 20)));
  eval(parse(text=MyExistsF(\"nSigma\", 10)));
  eval(parse(text=MyExistsF(\"cFolds\", 5)));
  eval(parse(text=MyExistsF(\"OrderSeq\", \"c(20,1)\")));
  eval(parse(text=MyExistsF(\"SABS\", .05)));
  eval(parse(text=MyExistsF(\"MaxCauchy\", 100)));
  eval(parse(text=MyExistsF(\"CauchyEpsilon\", .00001)));
  eval(parse(text=MyExistsF(\"L2ShrinkagePrior\", \"NULL\")));
  eval(parse(text=MyExistsF(\"RandomStartBetas\", 2)));
  eval(parse(text=MyExistsF(\"RecordL2Shrinkage\", \"FALSE\")));
  eval(parse(text=MyExistsF(\"KillBeforeFinal\", \"FALSE\")));
  eval(parse(text=MyExistsF(\"CVQuitTime\", -1)));
  eval(parse(text=MyExistsF(\"MaximumAllocation\", 4096)));
  eval(parse(text=MyExistsF(\"OnLambda2\", .001)));  
  eval(parse(text=MyExistsF(\"Lambda2\", -1)));  
  if (!exists(\"DoLogit\")) {
    if (all(Y %in% c(0,1))) {
      print(paste(\"Hey, DoLogit Does Not exist but Y is 0,1 setting to TRUE\")); 
      flush.console();
      DoLogit <- TRUE;
    } else {
      DoLogit <- FALSE;
    }
  }
    
  if (!exists(\"OnLambda2\")) { OnLambda2 = 0.0; }
  if (!exists(\"Lambda2\")) { Lambda2 = -1; }
  if (Lambda2 >= 0 && OnLambda2 == 0.0) {
    OnLambda2 <- Lambda2;
  }
  p = length(X[1,]);  n = length(Y);
  SampleInputs = TwoLassoSampleInputs(X,Y,nPi, nSigma, DoLogit=DoLogit);

    if (Verbose > 2) {
    print(\"-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): Creating M2Stuff\"); flush.console();
  }
  try(eval(parse(text=GetG0Text(\".M2Stuff\", \"globalenv()\",S=1))),silent=TRUE);
  .M2Stuff = M2Items(X, sigma=SampleInputs[,2], 
    OrderSeq = OrderSeq, SABS = SABS, Verbose = Verbose);
  eval(parse(text=SetGText(\".M2Stuff\", \"globalenv()\", S=1)));
  if (Verbose > 1) {
    print(\"-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): - gotten, paste to global\"); flush.console();
  }
  eval(parse(text=SetGText(\".M2Stuff\", \"globalenv()\")));
  LambdaAKInputs = .M2Stuff$LambdaAKInputs;
  LambdaDKInputs = .M2Stuff$LambdaDKInputs;   
  LambdaAK = .M2Stuff$LambdaASeq;
  LambdaDK = .M2Stuff$LambdaDSeq; 
  StartBetaMatrix <- NULL;
  if (is.null(CVQuitTime)) { CVQuitTime = -1; }
  if ( CVQuitTime >= -1) {
     eval(parse(text=SetGText(\"CVQuitTime\", \"globalenv()\", S=1)));
  }
 ");
 }
CFoldAllocateStartBetaMatrix <- function() {
  return("
    if (is.numeric(RandomStartBetas) && RandomStartBetas > 0) { 
     if (Verbose > 1) {
        print(\"-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): - Using RandomStartBetas\"); flush.console();
     }
     MyXeigen = NULL;
     if (p >= 10000) {
        print(\"-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): We cannot use Eigen values to choose a start. \"); flush.console();
        MyXeigen = NULL;
     } else {
       try(MyXeigen <- eigen(t(X) %*% X));
     }
     if (p > 10000) {
       try(MS <- (t(X) %*% Y ) / colSums(X*X));
       try(StartBetaMatrix <- matrix(0, p, ceiling(RandomStartBetas)));
       for (kk in 1:ceiling(RandomStartBetas)) {
         try(IID <- sort(sample(1:p, prob=abs(MS), size=100, replace=FALSE)));
         IDO <- colSums(X[,IID]*X[,IID]);
         IS <- colSums((Y-X[,IID] %*% MS[IID])^2) / NROW(X);
         try(StartBetaMatrix[IID,kk] <- MS[IID] + sqrt(IS)* rnorm(length(IID),0,1) / sqrt(IDO));
       }
     } else if (is.null(MyXeigen)) {
       try(StartBetaMatrix <- matrix(rnorm(p * ceiling(RandomStartBetas)),
         p, ceiling(RandomStartBetas)));
     } else {
       Di = MyXeigen$values;
       Di[abs(Di) >= .0001] = 1/ Di[abs(Di) >= .0001];
       try(Di[Di < 0] <- 0);
       try(Di[is.na(Di)] <- 0);
       PointMean <- t(MyXeigen$vectors) %*% diag(Di) %*% 
         MyXeigen$vectors  %*% t(X) %*% Y;
       SigmaEst =  sum( (Y - X %*% PointMean)^2) / (length(Y)-1);
       try(StartBetaMatrix <- matrix(rnorm(p * ceiling(RandomStartBetas)),
         p, ceiling(RandomStartBetas)));
       try(StartBetaMatrix <- matrix(rep(PointMean, RandomStartBetas),
         length(PointMean), RandomStartBetas) +
         sqrt(SigmaEst) * (t(MyXeigen$vectors) %*% 
           sqrt(diag(Di))) %*% MyXeigen$vectors %*% StartBetaMatrix);
       eval(parse(text=SetGText(\"StartBetaMatrix\", \"globalenv()\", S=1)));
     } 
     if (is.null(StartBetaMatrix)  || length(StartBetaMatrix) < p * ceiling(RandomStartBetas)) {
      print(\"-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp: Failure to not allocate the StartBetaMatrix. \"); flush.console();
      return(-1);
     }
  } else {
    StartBetaMatrix = NULL;
    eval(parse(text=SetGText(\"StartBetaMatrix\", \"globalenv()\", S=1)));
    if (Verbose > 1) {
      print(\"-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): Note, StartBetaMatrix set intentionally to NULL. \"); flush.console();
    }
  }
  
  if (Verbose > 1) {
    print(\"-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): Successful allocation StartBetaMatrix\"); flush.console();
  }
  ");  
    
}
TwoLassoCppFitStartHitFunction <- function() {
 return("
    eval(parse(text=MyExistsF(\"OrderSeq\", \"c(20,1)\")));
    if (!exists(\"OrderSeq\")) {
      eval(parse(text=GetG0Text(\".OrderSeq\", \"globanenv()\", S=1)));
      if (is.null(.OrderSeq) || (is.numeric(.OrderSeq) && .OrderSeq[1] == 0)) {
        .OrderSeq = c(20,1);  eval(parse(text=SetGText(\".OrderSeq\", \"globalenv()\", S=1)));
      }
      OrderSeq = .OrderSeq;
    }
    eval(parse(text=GetG0Text(\".M2Stuff\", \"globalenv()\", S=1)));
    LambdaAKInputs = .M2Stuff$LambdaAKInputs;
    LambdaDKInputs = .M2Stuff$LambdaDKInputs;   
    LambdaAK = .M2Stuff$LambdaASeq;
    LambdaDK = .M2Stuff$LambdaDSeq; 
    PiAVectorInputs = SampleInputs[,1];  SigmaVectorInputs = SampleInputs[,2];
   
   if (is.null(LambdaDKInputs) || length(LambdaDKInputs) <= 0) {
     print(\"ATwoLassoCpp.R:::CFoldValidateTwoLassoCpp():: TwoLassoStayFitFunction;  Uh oh, LambdaDKInputs is Null.\");
     flush.console();
   }
   if (is.null(LambdaAKInputs) || length(LambdaAKInputs) <= 0) {
     print(\"ATwoLassoCpp.R:::CFoldValidateTwoLassoCpp():: TwoLassoStayFitFunction;  Uh oh, LambdaAKInputs is Null.\");
     flush.console();
   }
   if (Verbose > 2) {
     print(\"ATwoLassoCpp.R:::CFoldValidateTwoLassoCpp():: Launching Two2LassoMiniValidate\"); flush.console();
   }
   try(eval(parse(text=GetG0Text(\".HiddenTwoLassoOb\", \"globalenv()\"))), silent=TRUE);  
   .HiddenTwoLassoOb = NULL;
    try(eval(parse(text=SetGText(\".HiddenTwoLassoOb\", \"globalenv()\"))), silent=TRUE); 
 ");
 }
CFoldInstallTwoLassoFit <- function() {
  return("
  eval(parse(text=GetG0Text(\".OrderSeq\", \"globalenv()\", S=1)));
  .OrderSeq = .M2Stuff$OrderSeq;  OrderSeq = .M2Stuff$OrderSeq;
  eval(parse(text=SetGText(\".OrderSeq\", \"globalenv()\", S=1)));
  #############################################################################
  ##  Set TwoLassoCppFitFunction
  ##     The \"CFoldValidate\" accepts two functions
  ##      a \"FitFunction\" and a \"RunOnAChange\" function.
  ##
  ##     With these functions, CFoldValidate automatically fits \"cFold\" number of
  ##  folds of the data.  The FitFunction does the heaviest lifting, generating
  ##  a fit.  Then the RunOnAChange function is designed, if possible to change
  ##  global parameters if that is faster than resetting data.  At this point
  ##  we will deal with a NULL action for RunOnAChange, because we do not consider
  ##  TwoLassoCpp objects stable enough to try to change X data within, and thus
  ##  Allocating another one is more efficient. 
  ##
  ##  Both TwoLassoCppFitFunction and RunOnAChange are assigned to globalenv()
  ##   here, and are then used on CFoldValidate.
  ##
  if (Verbose >= 1) {
    print(\"-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): Placing TwoLassoCppFitFunction. \"); flush.console();
  }
  try(eval(parse(text=GetG0Text(\"TwoLassoCppFitFunction\", \"globalenv()\", S=1))), silent=TRUE);
  TwoLassoCppFitFunction <- function(X,Y,SampleInput,StartBeta=-1,
    Verbose = 0,  StandardFlag=FALSE, TestX = NULL, ...) {
    eval(parse(text=TwoLassoCppFitStartHitFunction()));
   eval(parse(text=GetG0Text(\"CVQuitTime\", \"globalenv()\", S=1)));
   XMean <- colMeans(X);
   X <- t(t(X) - XMean);
   eval(parse(text=GetG0Text(\".HiddenTwoLassoOb\", \"globalenv()\")));
   .HiddenTwoLassoOb <- TwoLassoCpp:::TwoLassoCpp(X = X, Y = Y, 
      DoLogit = DoLogit,
      PiAVectorInputs = PiAVectorInputs, 
      SigmaVectorInputs = SigmaVectorInputs, 
      LambdaDK = LambdaDK, LambdaAK = LambdaAK, 
      LambdaDKInputs = LambdaDKInputs, LambdaAKInputs = LambdaAKInputs,
      OrderSeq = OrderSeq,  MaxCauchy = MaxCauchy * .5, 
      CauchyEpsilon = CauchyEpsilon * sqrt(p), 
      StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
      XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
      InitKKs = 5, WLSWeights = -1, TDFNu = -1,
      Verbose = Verbose, RecordFlag = 0,  Groupers = NULL,
      L2ShrinkagePrior = L2ShrinkagePrior, 
      RecordL2Shrinkage = RecordL2Shrinkage,StartBetaMatrix=StartBetaMatrix,
      RunFlag = TRUE,  CVQuitTime = CVQuitTime, MaximumAllocation = MaximumAllocation, OnLambda2 = OnLambda2,...);
    print(\"Finished .HiddenTwoLassoOb:TwoLassoCppFitFunction, Check CV Integrity \"); flush.console();
    .HiddenTwoLassoOb$TLS$CVDataIntegrity(3);
    print(\"  --- CVDataIntegrity Check Passed. \"); flush.console();
    if (is.numeric(.HiddenTwoLassoOb) && .HiddenTwoLassoOb[1] == -1) {
      print(\"ERROR: CFoldValidateTwoLassoCpp::Error.HiddenTwoLassoOb returned a NILL, check CVTLS\"); flush.console();
      MyT <- \"1 = 2\";
      eval(parse(text=MyT));
    }
    try(eval(parse(text=LockGText(\".HiddenTwoLassoOb\", \"globalenv()\"))), silent=TRUE);
    if (!is.null(TestX)) {
      TestX <- t(t(TestX) - XMean)
      FitY <- TestX %*% .HiddenTwoLassoOb$RecordBetaCVFinish
      if (DoLogit == TRUE) {
        FitY <- FitY + matrix(rep(.HiddenTwoLassoOb$TLS$RecordBeta0CV, each=NROW(TestX)), NROW(TestX), NCOL(FitY));
      }
      return(list(FitBeta = .HiddenTwoLassoOb$RecordBetaCVFinish, FitY = FitY));
    }    
    return(.HiddenTwoLassoOb$RecordBetaCVFinish);
  }   
  ");       
}

CFoldValidateInstallRunOnAChange <- function() {
  return("
  ## Assign \"RunOnAChange\" function to be NULL in this case:
  try(eval(parse(text=LockGText(\"TwoLassoCppFitFunction\", \"globalenv()\"))), silent=TRUE);
  try(parse(eval(text=GetG0Text(\"RunOnAChange\", \"globalenv()\",S=1))), silent=TRUE);
  RunOnAChange <- function(XTrain, YTrain,...) {  
  }
  try(parse(eval(text=SetGText(\"RunOnAChange\", \"gloablenv()\"))), silent=TRUE);
  if (Verbose >= 2) {
    print(\"ATwoLassoCpp.R:::CFoldValidateTwoLassoCpp():: Run CFoldValidate\"); flush.console();
  }
  try(eval(parse(text=SetGText(\"SampleInputs\", \"globalenv()\", S=1))));
  try(eval(parse(text=GetG0Text(\"MyCFold\", \"globalenv()\", S=1))), silent=TRUE);  
  ");
}
    
CFoldValidateTwoLassoCpp <- function(X,Y,DoLogit=FALSE, nPi=20, nSigma=10, cFolds =5,
  OrderSeq = c(20,1), SABS=.05, MaxCauchy = 100, CauchyEpsilon = .00001,
  Verbose = -1, RFastestCFoldValidate2Lasso = NULL,
  L2ShrinkagePrior = NULL,  RecordL2Shrinkage = FALSE,
  RandomStartBetas = 2, KillBeforeFinal = FALSE, CVQuitTime = -1,
  MaximumAllocation = 4096, OnLambda2 = .0001, Lambda2 = -1) {
  eval(parse(text=StartCFoldValidateText()));
  eval(parse(text=CFoldAllocateStartBetaMatrix()));
  if (DoLogit == FALSE && all(Y %in% c(0,1))) {
    print("ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): Hey, DoLogit is FALSE but Y is all 0,1");
    flush.console();
  }
  
  eval(parse(text=CFoldInstallTwoLassoFit()));
  eval(parse(text=CFoldValidateInstallRunOnAChange()));
  
  if (Verbose >= 1) {
    print("-- ATwoLassoCpp.R::CFoldValidateTwoLassoCpp(): About to Run CFoldValidate. "); flush.console();
  }
  MyCFold = CFoldValidate(X,Y,DoLogit=DoLogit, TwoLassoCppFitFunction,SampleInputs,
    cFolds = cFolds, RunOnAChange = RunOnAChange,
    CFoldVerbose=Verbose - 3,
    FitEmAll = TRUE);
  if (NCOL(X) >= 40000) {
    print(paste("Pasting MyCFold to globalenv()", sep="")); flush.console();
    eval(parse(text=SetGText("MyCFold", envir="globalenv()", S=1)));
  }

  if (Verbose >= 2) {
    print("ATwoLassoCpp.R:::CFoldValidateTwoLassoCpp():: Ran CFoldValidate.  Now Run M2Stuff"); flush.console();
  }
  try(parse(eval(text=SetGText("MyCFold", "gloablenv()"))), silent=TRUE); 
   M2Stuff <- M2Items(X, sigma=MyCFold$BestSampleInputs[2], 
      OrderSeq = OrderSeq, SABS = SABS);
  NewLambdaAK = M2Stuff$LambdaASeq;  NewLambdaDK = M2Stuff$LambdaDSeq;
    NewOrderSeq = M2Stuff$OrderSeq;  
  if (Verbose >= 0  || NCOL(X) >= 20000) {
    print(paste("ATwoLassoCpp.R:::CFoldValidateTwoLassoCpp():: Set now run FinalM2Out.  With Verbose = ", 
      Verbose-2, sep="")); flush.console();
    print(paste("      Successful CFoldValidate, with BestSampleInputs PiA=",
      MyCFold$BestSampleInputs[1], " and SigmaSq = ", MyCFold$BestSampleInputs[2], sep=""));
    flush.console();
    AMin <- min(5, length(NewLambdaAK));
    print(paste("      LambdaAK will be (", paste(NewLambdaAK[1:AMin], collapse=", "), ")", sep=""));
    try(print(paste("      LambdaDK will be (", paste(NewLambdaDK[1:AMin], collapse=", "), ")", sep="")));
    flush.console();
  } 
  if (!is.logical(KillBeforeFinal) || KillBeforeFinal == TRUE) {
    print("ATwoLassoCpp.R:::CFoldValidateTwoLassoCpp(): We are Killing before Final, this is not a normal exit.");
    flush.console();
    eval(parse(text=SetGText("MyCFold", "globalenv()", S=1)));
    eval(parse(text=SetGText("M2Stuff", "globalenv()", S=1)));
    return(MyCFold);
  }
  PiA <- MyCFold$BestSampleInputs[1];
  SigmaSq <- MyCFold$BestSampleInputs[2];
  
  if (Verbose >= 1) {
    print("----------------------------------------------------------------------------");
    print("---  ATwoLassoCpp.r::CFoldValidateTWoLassoCpp(): We have best inputs. "); flush.console();
    print(paste("---  Bets PiA = ", PiA, ",  Best SigmaSq = ", SigmaSq, sep="")); flush.console();
    print("---  About to run FinalM2Out. "); flush.console();
    print("----------------------------------------------------------------------------"); flush.console();
  }
  try(eval(parse(text=GetG0Text("FinalM2Out", "globalenv()", S=1))), silent=TRUE);
  FinalM2Out  <- TwoLassoCpp(X = X, Y = Y, 
    PiA = PiA, SigmaSq = SigmaSq, 
    LambdaDK = NewLambdaDK, LambdaAK = NewLambdaAK, OrderSeq = NewOrderSeq,
    MaxCauchy = MaxCauchy, CauchyEpsilon = CauchyEpsilon, 
    StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
    XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
    InitKKs = 5, WLSWeights = -1, TDFNu = -1,
    Verbose = Verbose-2, RecordFlag = 1, HoldOn = TRUE, 
    Groupers = NULL,  DoLogit=DoLogit,
    L2ShrinkagePrior = L2ShrinkagePrior, RecordL2Shrinkage = RecordL2Shrinkage,     
    RunFlag = 1, MaximumAllocation = MaximumAllocation);
  FinalM2Out$MyCFold = MyCFold;
  try(eval(parse(text=SetGText("FinalM2Out", "globalenv()", S=1))), silent=TRUE);
  if (Verbose >= 1) {
    print("ATwoLassoCpp.R:::CFoldValidateTwoLassoCpp():: All done, FinalM2Out"); flush.console();
  }
  return(FinalM2Out);  
}



ZeroOutTwoCpp <- function(TLS) {

print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
print("$$ ZeroOutTwoCpp"); flush.console();
print("$$ Hey: TLS$ is large before a sample we are going to blank up and sample small");
flush.console();

TLS$Beta <- rep(0,TLS$p);
MyTryProbs <- NULL;
MyTryProbTau <- NULL;  MyTryProbFixed <- NULL;
if (TLS$OnChainIter >= 2) {
  try(MyTryProbs <- MBS$ProbVector);
  if (!is.null(MyTryProbs)) {
    if (is.null(TLS$tauEndList)) {
      MyTryProbFixed <- log(MyTryProbs/(1.0-MyTryProbs));  
      MyTryProbFixed[MyTryProbFixed > 0 & !is.finite(MyTryProbFixed)] <- 99999;
      MyTryProbFixed[MyTryProbFixed < 0 & !is.finite(MyTryProbFixed)] <- -99999;
    } else if (MBS$FirstRandom > 1) {
      MyTryProbFixed <- log(MyTryProbs[1:MBS$iFirstRandom]/(1.0-MyTryProbs[1:MBS$iFirstRandom]))
      MyTryProbFixed[MyTryProbFixed > 0 & !is.finite(MyTryProbFixed)] <- 99999;
      MyTryProbFixed[MyTryProbFixed < 0 & !is.finite(MyTryProbFixed)] <- -99999;
    }
    if (!is.null(MBS$tauEndList)) {
      MyTryProbTau <- log(MyTryProbs[MBS$FirstRandom:length(MyTryProbs)] /
        (1.0 - MyTryProbs[TLS$FirstRandom:length(MyTryProbs)]));
      MyTryProbTau[MyTryProbTau > 0 & !is.finite(MyTryProbTau)] <- 99999;
      MyTryProbTau[MyTryProbTau < 0 & !is.finite(MyTryProbTau)] <- -99999;
    }
  }  
}
try(MBS$BlankAllNewCoords());
  if (length(MBS$tauEndList) >= 1) {
    if (!is.null(MyTryProbTau)) {
      MyProb <- MyTryProbTau;
    } else {
      MBS$SampleTausOnly <- 1;
      MBS$SampleNewTaus();
      MyProb <- MBS$ProbTau;
    }
    GoTau  <- MBS$OnTau;
    print(paste("$$ ZeroOutFunction: We sample Taus and the number of nonzero are ", length(GoTau[GoTau>0]),
     "/", length(GoTau), sep="")); flush.console();
    Fac = 1;
    if (length(MBS$tauEndList) >= 2)   {
      TauLens <- MBS$tauEndList[2:length(MBS$tauEndList)] -
        MBS$tauEndList[1:(length(MBS$tauEndList)-1)];
      TauLens <- c(TauLens, MBS$tauEndList- MBS$FirstRandom+1);
      Fac <- floor(mean(TauLens));
    } else if (length(MBS$tauEndList) == 1) {
      Fac <- MBS$tauEndList[1] - MBS$FirstRandom+1;
    }
    OnPiAT = MBS$OnPiA[1];
    if (length(MBS$OnPiA) == 2) { OnPiAT <- MBS$OnPiA[2]; }
    KeepCount <- round(5*length(GoTau)  * OnPiAT);
    if (KeepCount >= .5 * length(GoTau)) { KeepCount = round(.5* length(GoTau)); }
    if (KeepCount >= length(GoTau[GoTau>0])) { KeepCount <- length(GoTau[GoTau>0]); }
    if (KeepCount <= 0) { KeepCount = 1; }     
    if (KeepCount < length(GoTau[GoTau>0])) {
      Expt <- GiveMeNewBalance(lX=MyProb/2.0, WantTot = KeepCount, StartMove=1.0, MaxTTs = 300);
      if (length(Expt[Expt > 0]) < KeepCount) {
        KeepCount <- length(Expt[Expt > 0]);
      }
      print(paste("$$ after calculating Expt we have its sum is ", round(sum(Expt),4),
        "/", length(Expt), sep="")); flush.console();
      Keepers <- sample(1:length(GoTau), prob=Expt, replace=FALSE, size=KeepCount);
      print(paste("$$ Now we reduce active set to length(Keepers) = ", length(Keepers), sep=""));
      flush.console();
      DoTau <- rep(0, length(GoTau));
      DoTau[Keepers] <- GoTau[Keepers];
      try(MBS$DoAddCoordsOnSetBeta <- 0);
      try(MBS$BlankAllNewCoords());
      try(MBS$AssignTau(DoTau));
      print(paste("$$ Reduced DoTau length(Keepers) = ", length(Keepers), " is added", sep=""));
      print(paste("$$ After AssignTau, AllNewCoords = ", MBS$AllNewCoords, sep="")); flush.console();
      flush.console();
    } else {
      print(paste("$$ KeepCount = ", KeepCount, " but length(GoTau) = ",
        length(GoTau[GoTau>0]), "/", length(GoTau), " so no change.", sep=""));
      flush.console(); 
      try(MBS$BlankAllNewCoords());
      MBS$AssignTau(GoTau);
      print(paste("$$ After AssigningTau to GoTau length ", length(GoTau), " we ",
        " have AllNewCoords = ", MBS$AllNewCoords, sep=""));
    }
  }
  AllRandomCoords <- MBS$AllNewCoords;
  if(MBS$FirstRandom > 1 || length(MBS$tauEndList)==0) {
      if (length(MBS$tauEndList) == 0) {   pF = MBS$p;  
      } else { pF <- MBS$FirstRandom-1;}
      if (!is.null(MyTryProbFixed)) {
        BProbs <- MyTryProbFixed;
        if (Verbose >= 1) {
          print(paste("$$ ZeroOutFunction(): Probability of Historical fixed variables ",
            "has Historical MIPS summing to ", round(sum(dLogit(BProbs)),4), "/", pF, sep=""));
         flush.console();
        }
      } else {
        MBS$SampleTausOnly <- 1;
        print(paste("$$  ZeroOutFunction() : We're refreshing fixed effects. ")); flush.console();
        MBS$SampleFixedB();
        MBS$SampleTausOnly <- 0;
        BProbs <- MBS$ProbFixed;
        if (Verbose >= 1) {
          print(paste("$$ ZeroOutFunction() : MIPs taken from MBS SampleFixedB()", 
            round(sum(dLogit(BProbs)),4), "/", pF, sep=""));
          flush.console();
        }
      }     
      flush.console();
      OnPiAF <- MBS$OnPiA[1];
      KeepCount <- round(2.5*pF * OnPiAF);
      if (KeepCount <= 0) { KeepCount = 1; }
      if (KeepCount > pF) { KeepCount <- pF; }
      
      eBProbs <- GiveMeNewBalance(lX=BProbs/2.0, WantTot = KeepCount, StartMove=1.0, MaxTTs = 300);
      if (length(eBProbs[eBProbs > 0])  < KeepCount) {
        KeepCount <- length(eBProbs[eBProbs > 0]); 
      }
      if (Verbose >= 1) {
        print("$$ ZeroOutFunction() Sampling Keepers.");
      }
      Keepers <- sample(1:pF, prob=eBProbs, replace=FALSE, size=KeepCount);
      DKeepers <- (1:pF)[!(Keepers %in% (1:pF))]; 
      print(paste("$$  We have chosen to eliminate a count of Keepers ", length(DKeepers)));
      flush.console();
      if (MBS$FirstRandom <= 0 || is.null(MBS$tauEndList) || length(MBS$tauEndList) <= 0) {
        ABeta <- MBS$Beta;
        ABeta[MBS$Beta == 0.0 && MBS$BFixed >= 1] <- 1;
      } else if (MBS$FirstRandom > 1) {
        ABeta <- c(rep(1, MBS$FirstRandom-1), rep(0, MBS$p- MBS$FirstRandom+1)) * MBS$Beta;
        ABeta[(1:length(MBS$BFixed))[MBS$Beta == 0.0 && MBS$BFixed >= 1]] <- 1;
      }
      if (length(DKeepers) > 0) {
        ABeta[DKeepers] <- 0.0;
      }
      try(MBS$DoAddCoordsOnSetBeta <- 0);
     
      try(MBS$Beta <- ABeta);
      try(MBS$SetAllNewCoords(AllRandomCoords+MBS$NewFixedCoords, MBS$NewFixedCoords));
    }
    print(paste("$$  Adding initially ", MBS$AllNewCoords, "/", MBS$p, " new coords, maximum allotment is ",
      MBS$MaximumAllocation, sep=""));
    print(paste("$$  Resizing now to run one regression and hopefully eliminate noise variables. ", sep=""));
    try(MBS$AddAllNewCoords());
    try(MBS$RefreshOrderedActive(1));
    try(MBS$PrepareForRegression()); 
    print(paste("$$  Now sampling PropBeta which will have length ", MBS$NumActive, sep=""));
    try(MBS$SamplePropBeta());
    CNBeta <- MBS$PropBeta[abs(MBS$PropBeta) > 999999]
    if (any(is.nan(MBS$PropBeta)) || length(CNBeta) >= 1) {
      print(paste("$$ CNBeta has errors in PropBeta, please new work.  Triggering an error", sep="")); flush.console();
      -1 <- 2;
    }
    try(MBS$FillsBetaFromPropBetaAndCompute()); 
    try(MBS$UpdateSigma()); 
    try(print(paste("$$ ZeroOut Sampler, predictability after this is var(Y) = ",
      round(var(MBS$Y),6), " and YResid Squared is ", round(MBS$YResidSq,6), 
      " so Sigma = ", MBS$OnSigma, sep=""))); flush.console();
    print(paste("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$", sep=""));
    flush.console();
    return(1);
}

TwdLogit <- function(lX) {
  NI <- exp(lX) / (1.0+exp(lX));
  NI[lX < -10] = exp(-abs(lX[lX < -10]));
  NI[lX > 10] = 1.0 - exp(-abs(lX[lX > 10]));
  return(NI);
}
##  GiveMeNewBalance: Helper function
##
##  Let lX be a real vector
##    We are looking for a vector LB = exp(lX + alpha) / (1.0 + exp(lX+alpha))
##     for some alpha such that  sum(LB) = WantTot
##
TwGiveMeNewBalance <- function(lX, WantTot = 1, StartMove = 1, MaxTTs = 300, epsilon = .00001) {
  if (WantTot <= 0) {
    print(paste("Error:GiveMeNewBalance, does not work if WantTot=", WantTot, sep=""));
    return(-1);
  }
  if (length(lX) < WantTot) {
    print(paste("Error:GiveMeNeBalance, there are only ", length(lX),
      "members of lX not good enough to make ", WantTot, sep=""));
  }
  alphaOld = 0;
  AG <- TwdLogit(lX);
  if (length(AG[AG >= 1.0]) >= WantTot) {
    print(paste("GiveMeNewBalance, cannot be corrected, there are already ",
      " length(AG[AG >= 1.0]) = ",
      length(AG[AG >= 1.0]), " non zero, sample a subset then already!",
      sep="")); 
      flush.console();
    return(TwdLogit(lX));
  }
  if (length(AG[AG <= 0.0]) >= length(lX) - WantTot) {
    print(paste("GiveMeNewBalance, cannot be corrected, there are already ",
      "length(AG[AG <= 0.0]) = ",
      length(AG[AG <= 0.0]), " are zero, be happy its small!",
      sep="")); 
      flush.console();
    return(TwdLogit(lX));
  }
  LBOld <- sum(TwdLogit(lX));
  if (LBOld > WantTot) { StartMove = -1.0 *abs(StartMove) }
  alphaNew = alphaOld + StartMove;
  LBNew <- sum(TwdLogit(lX+alphaNew));
  tt = 0;
  if (LBNew > WantTot && LBOld > WantTot) {
    while(LBNew > WantTot && tt < MaxTTs) {
      alphaNew = alphaNew - abs(StartMove);
      LBNew <- sum(TwdLogit(lX+alphaNew)); tt <- tt+1;
      if (tt >= MaxTTs) {
        print(paste("GiveMeNewBalance, cannot succeed, tt = ", tt, 
          " and WantTot = ", WantTot, sep=""));  flush.console();
        return(-1);
      }
    }
    if (LBNew > WantTot) {
      print(paste("Error: GiveMeNewBalance in 300 moves we could not make WantTot = ",
        WantTot, " for LBNew = ", LBNew, " for alphaNew = ", alphaNew, sep=""));
      flush.console();  return(-1);
    }
  } else if (LBNew < WantTot && LBOld < WantTot) {
    while(LBNew < WantTot && tt < MaxTTs) {
      alphaNew = alphaNew + abs(StartMove);
      LBNew <- sum(TwdLogit(lX+alphaNew)); tt <- tt+1;
      if (tt >= MaxTTs) {
        print(paste("GiveMeNewBalance, cannot succeed, tt = ", tt, 
          " and WantTot = ", WantTot, sep=""));  flush.console();
        return(-1);
      }
    }
    if (LBNew < WantTot) {
      print(paste("Error: GiveMeNewBalance in 300 moves we could not make large WantTot = ",
        WantTot, " for LBNew = ", LBNew, " for alphaNew = ", alphaNew, sep=""));
      print(paste(" We will return a smaller form. ")); flush.console();
      return(TwdLogit(lX+1.0));
      flush.console();  return(-1);
    }
  }
  tt = 0;
  if (alphaNew < alphaOld) {
    LBMid <- LBNew;  alphaMid <- alphaNew;
    alphaNew <- alphaOld;  LBNew <- LBOld;
    alphaOld <- alphaMid;  LBOld <- LBMid;
  }
  while(abs(LBNew- WantTot) > epsilon && tt < MaxTTs) {
    alphaMid <- (alphaNew+alphaOld)  / 2.0;
    LBMid <- sum(TwdLogit(lX+alphaMid));  tt <- tt+1;
    if (abs(LBMid-WantTot) <= epsilon) {
      return(TwdLogit(alphaMid+lX));
    } else if (LBMid > WantTot) {
      alphaNew <- alphaMid;  LBNew <- LBMid;
    } else {
      alphaOld <- alphaMid;  LBOld <- LBMid;
    } 
  }
  return(TwdLogit(alphaMid+lX));  
}

print("ATwoLassoCpp.R: complete Loading")