################################################################################
##
##   TwoLarsDirected.R
##
##  (c) Alan Lenarcic 2009-2010:
##
##  Work primarily from Alan Lenarcic's 2009 thesis research, and this is the
##  package code that was used at that time.  Some of this is to run simulation
##  code and attempt to compare estimators, both implemented in this package
##  and elsewhere.
##
##
##     This uses files  
##         TwoLarsObject.cc,  LARSObject.cc, CoordinateDescentObject.cc,
##         MyMatrixOp.cc, IncidentalOp.cc, CovMatrixEstimation.cc
##           And their .h files.
##
##      Though the relevant algorithm in TwoLarsObject.cc  is
##      "LarsConvergencyBI" which calls  "RunConvergentLars",
##      And they are both called by EM2Lasso
##
##     Simulation code is in function "SimForMe2", much of the code contained
##      is for plotting, running, and saving many simulations of regression
## 
##     ## Other code is avi
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



DeclareGlobals <- function() {
  Verbose = 3;
  if (Verbose >= 1) {
    print("DeclareGlobals, start: InitiateS. "); flush.console();
  }
  eval(parse(text=GetG0Text("InitiateS", S=1)));
  InitiateS = "BFOUR";
  eval(parse(text=SetGText("InitiateS", S=1)));
  try(eval(parse(text=SetGText("InitiateS", envir="TWOLASSONAMESPACE", S=1))));
  eval(parse(text=GetG0Text("InitiateS2", S=1)));
  InitiateS2 = "BSIXB8";
  eval(parse(text=SetGText("InitiateS2", S=1)));
  try(eval(parse(text=SetGText("InitiateS2", envir="TWOLASSONAMESPACE", S=1))));

  if (Verbose >= 1) {
    print("DeclareGlobals, Load required libraries lars, corpcor."); flush.console();
  }
  try(library(lars, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE)
  try(library(corpcor, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE)

  eval(parse(text=GetG0Text("ALT", S=1)));
  ALT = 7 * 6 / 2;
  eval(parse(text=SetGText("ALT", S=1)));
  try(eval(parse(text=SetGText("ALT", envir="TWOLASSONAMESPACE", S=1))));

  eval(parse(text=GetG0Text("ALT2", S=1)));
  ALT2 = 10 * 11 / 2;
  eval(parse(text=SetGText("ALT2", S=1)));
  try(eval(parse(text=SetGText("ALT2", envir="TWOLASSONAMESPACE", S=1))));
    
  eval(parse(text=GetG0Text("ALT3", S=1)));
  ALT3 = 10;
  eval(parse(text=SetGText("ALT3", S=1)));
  try(eval(parse(text=SetGText("ALT3", envir="TWOLASSONAMESPACE", S=1))));
  
  if (Verbose >= 1) {
    print("DeclareGlobals, start: RunSetSaveHomes. "); flush.console();
  }
  MyAT <- -1;
  MyTT <- "
    try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
    try(AlanDirectories:::SetSaveHomes());
    MyAT == 1;
  "
  try(eval(parse(text=MyTT)));
  if (MyAT == -1) {
    SaveHome <- "";
  }
  
  if (Verbose >= 1) {
    print("DeclareGlobals, start: Set PathMe. "); flush.console();
  }
   eval(parse(text=GetG0Text("PathMe2", S=1)));
   try(PathMe2 <- paste(SaveHome, "//2008Summer", sep=""));
   try(dir.create(PathMe2, showWarnings  = FALSE));
   try(PathMe2 <- paste(PathMe2, "//LarsProject", sep=""));
   try(dir.create(PathMe2, showWarnings=FALSE));
   try(PathMe2 <- paste(PathMe2, "//code", sep=""));
   try(dir.create(PathMe2, showWarnings=FALSE));
   try(eval(parse(text=SetGText("PathMe2", S=1))));
   try(eval(parse(text=SetGText("PathMe2", envir="TWOLASSONAMESPACE", S=1))));
     
   eval(parse(text=GetG0Text("PathMeFiles", S=1)));
   if (!exists("PathMe")) {
     try(PathMe <- MakePathMe());
   }
   if (!exists("PathMe")) {
     PathMe = "";
   }
   if (!exists("PathMe")) {
     DefaultTwoLassoCppDir = ""; PathMe = "";
     try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE), 
       silent=TRUE); try(AlanDirectories:::SetSaveHomes());
     try(eval(parse(text=GetG0Text("PathMe"))));
     try(PathMe <- DefaultTwoLassoCppDir);
     try(eval(parse(text=SetG0Text("PathMe"))));                              
     
   }
   if (exists("PathMe")) {
     try(PathMeFiles <- paste(PathMe, "//SavedOutput", sep=""));
   } else {
     try(PathMeFiles <- "SavedOutput", sep="");
   }
   try(dir.create(PathMeFiles, showWarnings=FALSE));
   try(eval(parse(text=SetGText("PathMeFiles", S=1))), silent=TRUE);
   try(eval(parse(text=SetGText("PathMeFiles", envir="TWOLASSONAMESPACE", S=1))),
     silent=TRUE);
   
  eval(parse(text=GetG0Text("n1", S=1)));
  n1 = 100;
  eval(parse(text=SetGText("n1", S=1)));
  try(eval(parse(text=SetGText("n1", envir="TWOLASSONAMESPACE", S=1))));

  if (Verbose >= 1) {
    print("DeclareGlobals, start: BetaBarMin. "); flush.console();
  }
  eval(parse(text=GetG0Text("BetaBarMin", S=1)));
  BetaBarMin = 4;
  eval(parse(text=SetGText("BetaBarMin", S=1)));
  try(eval(parse(text=SetGText("BetaBarMin", envir="TWOLASSONAMESPACE", S=1))));

  eval(parse(text=GetG0Text("MinMin", S=1)));
  MinMin = .005;
  eval(parse(text=SetGText("MinMin", S=1)));
  try(eval(parse(text=SetGText("MinMin", envir="TWOLASSONAMESPACE", S=1))));
      
  eval(parse(text=GetG0Text("PathMe", S=1)));
  if (substr(.Library,1, nchar("/nas")) == "/nas") {
    try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
    PathMe = "";
    try(AlanDirectories:::SetSaveHomes(), silent=TRUE);
    try(PathMe <- DefaultTwoLassoCppDir, silent=TRUE);
  } else {
    try(PathMe <- paste( .Library, "/TwoLasso/", sep="") );
  }
  eval(parse(text=SetGText("PathMe", S=1)));
  try(eval(parse(text=SetGText("PathMe", envir="TWOLASSONAMESPACE", S=1))));

  if (Verbose >= 1) {
    print("DeclareGlobals, start: SavedOutputDirectory. "); flush.console();
  }  
  eval(parse(text=GetG0Text("SavedOutPutDirectory", S=1)));
  try(SavedOutPutDirectory <- PathMe );
  eval(parse(text=SetGText("SavedOutPutDirectory", S=1)));
  try(eval(parse(text=SetGText("SavedOutPutDirectory", envir="TWOLASSONAMESPACE", S=1))));
  
   LoadTwoCC();
}

SecondLoadfunction <- function() {
  ################################################################
  ##  When Errors are created in Simulations we add to these lists
  ##     for debugging purposes
  ##
  ##
  eval(parse(text=GetG0Text("FailToFitSim", S=1)));
  try(FailToFitSim <- list() );
  eval(parse(text=SetGText("FailToFitSim", S=1)));
  try(eval(parse(text=SetGText("FailToFitSim", envir="TWOLASSONAMESPACE", S=1))));

  eval(parse(text=GetG0Text("FailToFitEMOb", S=1)));
  try(FailToFitEMOb <- list() );
  eval(parse(text=SetGText("FailToFitEMOb", S=1)));
  try(eval(parse(text=SetGText("FailToFitEMOb", envir="TWOLASSONAMESPACE", S=1))));

  eval(parse(text=GetG0Text("OnFails", S=1)));
  try(OnFails <- 0 );
  eval(parse(text=SetGText("OnFails", S=1)));
  try(eval(parse(text=SetGText("OnFails", envir="TWOLASSONAMESPACE", S=1))));
  
}
##
MakePathMe <- function() {
  eval(parse(text=GetG0Text("PathMe", S=1)));
  try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
  try(AlanDirectories:::SetSaveHomes());
  PathMe = "";
  try(PathMe <- DefaultTwoLassoCppDir, silent=TRUE);
  if (!exists("DefaultTwoLassoCppDir") || DefaultTwoLassoCppDir == "" ||
    PathMe == "") {
    try(PathMe <- paste( .Library, "/TwoLassoCpp/", sep="") );
  }
  try(dir.create(PathMe, showWarnings=FALSE));
  eval(parse(text=SetGText("PathMe", S=1)));
  try(eval(parse(text=SetGText("PathMe", envir="TWOLASSONAMESPACE", S=1))));
  
  return(PathMe);
}
LoadTwoCC <- function() {
  ##dyn.unload(paste(PathMe,"TwoLarsOperatorProgram.dll", sep="")) 
  if (!exists(".Library")) {
    .Library <- "";
  }
  FilesInDIR <- unlist(list.files(.Library));
  if (any(FilesInDIR == "TwoLasso")) {
    try(eval(parse(text=TwoLassoCpp:::GetG0Text("PathMe", S=1))));
    PathMe <- MakePathMe();
    eval(parse(text=SetGText("PathMe", S=1)));
  } else {
    eval(parse(text=GetG0Text("PathMe", S=1)));
    try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
    try(AlanDirectories:::SetSaveHomes());
    if (!exists("SAVEHOME")) {
      try(SAVEHOME <- paste(.Library, "/TwoLassoCpp/", sep=""));
      try(PathMe <- SAVEHOME);
    } else {
      try(PathMe <- DefaultTwoLassoCppDir);
      try(dir.create(PathMe, showWarnings=FALSE));
    }
    try(dir.create(PathMe, showWarnings=FALSE));
    eval(parse(text=SetGText("PathMe", S=1)));
    try(eval(parse(text=SetGText("PathMe", envir="TWOLASSONAMESPACE", S=1))));
  }
  
  #################
  ### Load Library Lars  
  if (any(FilesInDIR == "lars")) {
    try(library(lars, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
  } else {
    print("TwoLasso Requires Installation of Package \"lars\"");
    flush.console();
  }
}

###########################################################
### PenaltyFunction: 
###
###   Standard Lasso Penalty function, adding a second  L1 lag penalty
PenaltyFunction <- function(yys, Betas, xxs, v1, v2, BetOlds) {
  sum((as.vector(yys) - as.vector(xxs %*% Betas))^2) + 
    v1 * sum(abs(Betas)) + v2* sum(abs(Betas-BetOlds));
}



##########################################3
###  meanC
###   typical Mean of columns function
meanC <- function(MatX) {
  Ret <- vector("numeric", length(MatX[1,]));
  for (ii in 1:length(Ret)) { Ret[ii] <- mean(MatX[,ii]) }
  return(Ret);
}
   
################################################
####  StandardizeXX();
####  Standardizes matrix vector XXX if applicable
####
StandardizeXX <- function(XX) {
  if ( length(dim(XX)) == 2) {
        XX <- t( ( t(XX) - meanC(XX)) / apply(XX,2,sd) );
  }  else {
        XX <-  (XX - mean(XX)) / sd(as.vector(XX));
  }
  return(XX);
}   

##########################################################
## EM-RIDGE Method: Coded in R, with more EM iterations
##   based upon finding the "pi_A" value.  This is semi
##   archaic, "Maxmizes the partial model", then maximizes
##   p after a few iterations, 
##   and a better coded method would include
##   a solution coded within the algorithm.  However nothing
##   seems perhaps wrong with the code.
EMRIDGESolveP <- function(yys, xxs, n1 = 128, n2 = 40, EMiters = 10, 
     pstart = -999, p1p = 1, p2p = 2, sigmaalpha = .5, sigmatau = 1,
     tauDsqStart=1, tauAsqStart=4, tauDMult=-999, tauAMult = -999,
     tauAMultStop=0, BStartBetas=-1, SigmaSqNoiseSt=.5,
     ZeroCriterion =4) {
     if(!is.loaded("EMRIDGEShell2") ) {
        LoadTwoCC();
     }
     LoadTwoCC();
     p1p = abs(p1p);
     p2p = abs(p2p);
     if (pstart < 0 || pstart > 1) {
         pstart = p1p / (p1p + p2p);
     }
     puse = pstart;
     Returner = list( ReturnBetasPI = matrix(0, EMiters, length(XX[1,])),
                      BBHatsKeepPI = matrix(0, EMiters, length(XX[1,])),
          sigmaNoiseFit = vector("numeric", EMiters), 
          pfit = vector("numeric", EMiters), 
          Elpfit = vector("numeric", EMiters), El1mpfit = vector("numeric", EMiters),
          sigmaInvNoiseFit = vector("numeric", EMiters), 
          ReturnBetas = vector("numeric", length(xxs[1,])) );
    ## );
      XX = xxs; YY = yys;
     SigmaNoiseForInput = SigmaSqNoiseSt;
     E1OverSigma = 1 / SigmaNoiseForInput;
     Elp = log(puse);
     El1mp = log( 1- puse);
     for ( ii in 1:EMiters) {
     MyResult <- EMRIDGE(yys=yys, xxs=xxs, TotalRuns = 128, NumEMConv = 40, ppiuse = puse,
         tauDsqStart = tauDsqStart, tauAsqStart = tauAsqStart, tauDMult = tauDMult, 
         tauAMult = tauAMult, tauAMultStop = tauAMultStop, BStartBetas = -1, SigmaSqNoise =  SigmaNoiseForInput, 
         EFlag = TRUE, 
         logp = Elp, log1mp = El1mp, SigmaNoiseInv = E1OverSigma );   
         NewNoise <-  sum( (YY - XX %*% MyResult$BEndBetas)^2 );
         E1OverSigma <-  (sigmaalpha + length(YY)-1) / (sigmatau + NewNoise );
         SigmaNoiseForInput <- 1 / E1OverSigma;
         ETrueNoise <- (sigmatau + NewNoise) / (sigmaalpha + length(YY) -3);
         NumNo = length(MyResult$BEndBetas[round(MyResult$BEndBetas, ZeroCriterion) == 0]);
         pp <- 1:99999 / 100000
         Elp <- sum (log(pp) * dbeta( pp, p1p + (length(XX[1,]) - NumNo), p2p + NumNo)  * 
                       (pp[2] - pp[1]));
         El1mp <- sum(log(1-pp) * dbeta( pp, p1p + (length(XX[1,]) - NumNo), p2p + NumNo) * 
                        (pp[2] - pp[1]) );
         puse <- sum ( pp * dbeta( pp, p1p + (length(XX[1,]) - NumNo), p2p + NumNo)  * 
                       (pp[2] - pp[1]));
         Returner$ReturnBetasPI[ii,] <- MyResult$BEndBetas;
         Returner$sigmaNoiseFit[ii] <- ETrueNoise;
         Returner$pfit[ii] <-puse;
         Returner$Elpfit[ii] <- Elp;
         Returner$El1mpfit[ii] <- El1mp;
         Returner$sigmaInvNoiseFit <- E1OverSigma;
         ##Returner$BBHatsKeepPI[ii,] <- 
         xxx = 1:length(MyResult$ReturnBetasAll[,1,n2]);
         yylim = c(min(MyResult$ReturnBetasAll[,,n2]), 
                   max(MyResult$ReturnBetasAll[,,n2]));
         kLen = length(XX[1,]);
         cols <- rainbow(kLen);
         plot(x=xxx,
             y=MyResult$ReturnBetasAll[1,,n2], main=paste("On EMiter = ", ii, sep=""), type="n",
             ylim = yylim);
           for (jjk in 1:kLen) { 
              lines(x=xxx, y=MyResult$ReturnBetasAll[jjk,,n2],  col = cols[jjk],
              lty = jjk, lwd=2);
           }
     }
     Returner$ReturnBetas <- MyResult;
     RetCEV$Type = "EMPIRIDGE";
     return(Returner);
}

########################################################################3
##  EMRIDGE  Algorithm
##
##    Uses converging Gaussian distriubtions to pick apart active variables.
##
##  There are optional function versions with "NOPRINT" and "EFlag" that simplify output.
##   
EMRIDGE <- function(xxs, yys, TotalRuns=128, NumEMConv=4, ppiuse = .5,
     tauDsqStart = -999, tauAsqStart = -999, tauDMult = -999, tauAMult = -999, 
     tauAMultStop = 0, BStartBetas = -1, sigmaNoiseSq = .5, SigmaSqNoise= -9, EFlag = FALSE,
     logp = -1, log1mp = -1, SigmaNoiseInv = -1, NOPRINT = TRUE, n1 = -9, n2 = -9,
     NumCDOConv = 40, CDOEpsilon = .00001, TauSqASeq = -1, TauSqDSeq = -1, OrderSeq = -1) { 
TOTALMAXK <- 50;
     XX <- xxs; YY <- yys;  puse <- ppiuse; n1 <- TotalRuns;
     if (n2 < 0) {n2 <- NumEMConv;} else {NumEMConv = n2;}
     if (n1 < 0) { n1 <- TotalRuns; } else {TotalRuns = n1;}
     if (SigmaSqNoise  == -9) {SigmaSqNoise <- sigmaNoiseSq; }
     if(!is.loaded("CalculateGFTAllShell") ) {
        LoadTwoCC();
     }
     if (puse >= 1) {puse = .99; logp = log(puse);  log1mp = log(1-puse);}
     if (puse <= 0) {puse = .01; logp = log(puse);  log1mp = log(1-puse);}
     if (tauDsqStart == -999) {
        XEffect <- sd(as.vector(YY)) / mean(apply(XX,2,sd));
        tauDsqStart = (.5 * XEffect)^2;
        tauAsqStart = (2 * XEffect)^2;
        SigmaSqNoise = (.25 * sd(as.vector(YY)))^2;
        SigmaNoiseInv = 1 / SigmaSqNoise;
     }
     if (tauDMult == -999) {
        tauDMult = 2^(-1/8);
     }
     if (tauAMult == -999) {
        tauAMult = 2^(1/8);
     } 
     if (length(OrderSeq) == 1) {
         OrderSeq = rep(NumEMConv, TotalRuns);
     }
     if (length(TauSqASeq) == 1) {
         TauSqASeq = tauAsqStart * tauAMult^(0:(TotalRuns));
         if (tauAMultStop > 0 & tauAMultStop < TotalRuns) {
             TauSqASeq[tauAMultStop:TotalRuns] <- tauAsqStart * tauAMult^tauAMultStop;
         }
         TauSqDSeq = tauDsqStart * tauDMult^(0:(TotalRuns));         
     }
     if (length(TauSqDSeq) == 1) {
         TauSqASeq = tauAsqStart * tauAMult^(0:(TotalRuns));
         if (tauAMultStop > 0 & tauAMultStop < TotalRuns) {
             TauSqASeq[tauAMultStop:TotalRuns] <- tauAsqStart * tauAMult^tauAMultStop;
         }
         TauSqDSeq = tauDsqStart * tauDMult^(0:(TotalRuns));         
     }     
     
     if (logp == -1) { logp = log(ppiuse); }
     if (log1mp == -1) { log1mp = log(1-ppiuse); }
     if (SigmaNoiseInv == -1) { SigmaNoiseInv = 1 / SigmaSqNoise; }
     BStartBBs <- puse + 0 * (1:length(XX[1,]));
     BStartBetas <- 0 * (1:length(XX[1,]));
  if (det(t(XX) %*% XX) > .000001) {
       BStartBetas <- solve(t(XX) %*% XX) %*% t(XX) %*% YY;
       BStartBBs <- 0 * (1:length(XX[1,])) + ppiuse;
  } else  if (tauDsqStart != tauAsqStart) {
xxSt <- sqrt( abs( ( .5 * log(2 * pi * tauDsqStart) - .5 * log(2*pi * tauAsqStart)
    + log(puse) - log(1-puse) ) /
        ( 1 / (2 * tauDsqStart) - 1 / (2 * tauAsqStart) ) ) );
    BStartBetas = sign( t(XX) %*% YY) * xxSt;
    BStartBetas <-1/diag(t(XX) %*% XX) * t(XX) %*% YY
    BStartBBs <- puse + 0 * (1:length(XX[1,]));
  }    else {
    BStartBBs <- puse + 0 * (1:length(XX[1,]));
    BStartBetas <-  0 * (1:length(XX[1,]));
    BStartBetas <-1/diag(t(XX) %*% XX) * t(XX) %*% YY
  }

##kLen = length(BStartBetas);
   kLen = length(XX[1,]);
   NLen = length(YY);
##n1 = 128; 
##n2 = 40;
if (NOPRINT == TRUE) {
  ReturnBetasAllEM <- c(-999,-999);
  BBHatsAllEM <- c(-999,-999);
  ##print("NOPRINT IS TRUE we Won't have any ReturnBetasAll or BBHatsAll");
} else if (kLen > TOTALMAXK && NumCDOConv >= 0) {
  ReturnBetasAllEM <- as.double(vector("numeric", ((kLen+1) * (n1+2))+2));
  BBHatsAllEM <- as.double(vector("numeric", ((kLen+1) * (n1+2)) + 2));
  ##print(paste("ReturnBetasAll length = ", length(ReturnBetasAll), 
  ##   " and length(BBHatsAll) = ", length(BBHatsAll), sep=""));
} else if (EFlag == FALSE) {
  ReturnBetasAllEM <- as.double(vector("numeric", kLen * n1 * (n2+1)));
  BBHatsAllEM <- as.double(vector("numeric", kLen * n1 * (n2+1)));
  ##print("EMRIDGE: Performing EMDYMANICRIDGEShell2");
} else if (kLen > TOTALMAXK && NumCDOConv < 0) {
  ReturnBetasAllEM <- as.double(vector("numeric", kLen * (n1+1) * (n2+1)));
  BBHatsAllEM <- as.double(vector("numeric", kLen * (n1+1) * (n2+1)));
} else {
  ReturnBetasAllEM <- as.double(vector("numeric", kLen * (n1+1) * (n2+1)));
  BBHatsAllEM <- as.double(vector("numeric", kLen * (n1+1) * (n2+1)));
  ##print("EMRIDGE: Performing EMDYMANICRIDGEShell2");
}
flush.console();
     if (is.na(puse) || is.nan(puse) || !is.finite(puse)) {
      print(paste("Heck, puse is ", puse, " warning!"));
       puse = .5;
     }
     if (is.na(logp) || is.nan(logp) || !is.finite(logp)) {
      print(paste("Heck, logp is ", logp, " warning!"));
       logp = log(puse);
     }
    if (is.na(log1mp) || is.nan(log1mp) || !is.finite(log1mp)) {
      print(paste("Heck, log1mp is ", log1mp, " warning!"));
       log1mp = log(1-puse);
     }     
     if (puse >= 1) {puse = .99; logp = log(puse);  log1mp = log(1-puse);}
     if (puse <= 0) {puse = .01; logp = log(puse);  log1mp = log(1-puse);}
     if (is.na(tauDsqStart) || is.nan(tauDsqStart)) {
       print("tauDsqStart is NA sorry");
       tauDsqStart = 1;
     }
     if (!is.finite(tauDsqStart) ) {
       print("tauDsqStart is Infinite sorry\n");
       tauDsqStart = 1;
     }
     if (is.na(tauAsqStart) || is.nan(tauAsqStart) || !is.finite(tauAsqStart)) {
       print(paste("tauAsqStart is", tauAsqStart, "sorry\n", sep=""));
       tauDqStart = 1;
     }
     if (is.na(tauAMult) || is.nan(tauAMult) || !is.finite(tauAMult)) {
       print(paste("tauAMult is", tauAMult, "sorry", sep=""));
       print("tauAMult is NA sorry");
       tauAMult = .99;
     }
     if (is.na(tauDMult) || is.nan(tauDMult)) {
       print(paste("tauDMult is", tauDMult, "sorry", sep=""));
       print("tauDMult is NA sorry");
       tauDMult = 1.01;
     }     
     if (is.na(SigmaSqNoise) || is.nan(SigmaSqNoise) || !is.finite(SigmaSqNoise)) {
        print(paste("SigmaSqNoise is", SigmaSqNoise, "sorry", sep=""));
        SigmaSqNoise = 1;
     }
##print("EMRIDGE: Performing EMDYMANICRIDGEShell2 #2");
RetCEV <- .C("EMDYMANICRIDGEShell2", NLen = as.integer(length(YY)), kLen=as.integer(length(XX[1,])), 
         YY = as.double(YY), XX = as.double(XX), 
         ChangeNSteps = as.integer(n1), ConvergenceNSteps = as.integer(n2), 
         SigmaSqNoise = as.double(SigmaSqNoise),
         ppiuse = as.double(puse), 
         tauDsqStart =as.double(tauDsqStart), 
         tauAsqStart = as.double(tauAsqStart), 
         tauDMult = as.double(tauDMult), tauAMult = as.double(tauAMult), 
         tauAMultStop = as.integer(tauAMultStop), 
         TauSqASeq = as.double(TauSqASeq), TauSqDSeq = as.double(TauSqDSeq),
         BStartBBs = as.double(BStartBBs),
         BStartBetas = as.double (BStartBetas),
         ReturnBetasAllEM = as.double(ReturnBetasAllEM),
         ReturnBetas = as.double(vector("numeric", kLen+5)),
         BBHatsAllEM = as.double(BBHatsAllEM),
         BBHatsFinal = as.double(vector("numeric", kLen +5)),
         logp = as.double(logp), log1mp = as.double(log1mp), 
         SigmaNoiseInv = as.double(SigmaNoiseInv),
         NumCDOConv = as.integer(NumCDOConv), CDOEpsilon = as.double(CDOEpsilon),
         OrderSeq = as.integer(OrderSeq),
         CountDescentOperations = as.integer(vector("numeric", TotalRuns+2))
      ); 
      kLen = RetCEV$kLen;
      RetCEV$ReturnBetas = RetCEV$ReturnBetas[1:kLen];
      RetCEV$BBHatsFinal = RetCEV$BBHatsFinal[1:kLen];
      RetCEV$CountDescentOperations = RetCEV$CountDescentOperations[1:RetCEV$ChangeNSteps];
      RetCEV$puse = RetCEV$ppiuse;
  n1 = RetCEV$ChangeNSteps;
  RetCEV$tauAK <- RetCEV$tauAsqStart[1] * RetCEV$tauAMult^(0:(n1-1));
  if (RetCEV$tauAMultStop[1] < n1 && RetCEV$tauAMultStop > 0) {
      RetCEV$tauAK[(RetCEV$tauAMultStop+1):n1] <- RetCEV$tauAK[RetCEV$tauAMultStop];
  }
  RetCEV$tauDK <- RetCEV$tauDsqStart * RetCEV$tauDMult^(0:(n1-1));
  MP <- round( length(RetCEV$tauDK) /2);
  try(MP <- sort( abs(RetCEV$SigmaSqNoise[1] - 1/RetCEV$tauDK+1/RetCEV$tauAK), index=TRUE)$ix[1]);
  try(RetCEV$MP <- MP);
                                       
if (NOPRINT == TRUE) {
  RetCEV$ReturnBetasAllEM <- NULL;  RetCEV$BBHatsAllEM <- NULL;
  RetCEV$MIP <- NULL;
} else if (kLen > TOTALMAXK && NumCDOConv >= 0) {
  RetCEV$ReturnBetasAllEM <- t(matrix(RetCEV$ReturnBetasAllEM[1:(kLen * n1)], n1, kLen));
  RetCEV$BBHatsAllEM <- t(matrix(RetCEV$BBHatsAllEM[1:(kLen * n1)], n1, kLen));  
  try(RetCEV$MIP <- RetCEV$BBHatsAllEM[,MP]);
} else if (EFlag == FALSE) {
  RetCEV$ReturnBetasAllEM <- array(RetCEV$ReturnBetasAllEM[1:(kLen * n1 * n2)], dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") );
  RetCEV$BBHatsAllEM <- array(RetCEV$BBHatsAllEM[1:(kLen*n1*n2)], dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") ); 
  try(RetCEV$MIP <- RetCEV$BBHatsAllEM[,n2,MP]);
} else if (kLen > TOTALMAXK && NumCDOConv < 0) {
  RetCEV$ReturnBetasAllEM <- array(RetCEV$ReturnBetasAllEM[1:(kLen * n1 * n2)], dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") );
  RetCEV$BBHatsAllEM <- array(RetCEV$BBHatsAllEM[1:(kLen*n1*n2)], dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") ); 
  try(RetCEV$MIP <- RetCEV$BBHatsAllEM[,n2,MP]);
} else {
  RetCEV$ReturnBetasAllEM <- array(RetCEV$ReturnBetasAllEM[1:(kLen * n1 * n2)], 
                                      dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") );
  RetCEV$BBHatsAllEM <- array(RetCEV$BBHatsAllEM[1:(kLen*n1*n2)], 
                              dim = c(kLen, n2, n1), dimnames = 
                                    c("k ", "ConvIter", "ChangeIter") ); 
  try(RetCEV$MIP <- RetCEV$BBHatsAllEM[,n2,MP]);
   
}      
 

  RetCEV$Type = "EMRIDGE";
  if (RetCEV$BStartBetas[1] == -999) {
     print("RetCEV  Signals RIDGE Failure \n");
     RetCEV$FailureFlag = 1;
  } else {
     RetCEV$FailureFlag = 0;
  }
  return(RetCEV);    
}

  
########################################################################3
##  EMOLDRIDGE  Algorithm
##     
##    Uses converging Gaussian distriubtions to pick apart active variables.
##
##    There are optional function versions with "NOPRINT" and "EFlag" that simplify output.
##    This one calls a version of EMTWORIDGE not using Coordinate Descent
EMOLDRIDGE <- function(xxs, yys, TotalRuns=128, NumEMConv=4, ppiuse = .5,
     tauDsqStart = -999, tauAsqStart = -999, tauDMult = -999, tauAMult = -999, 
     tauAMultStop = 0, BStartBetas = -1, sigmaNoiseSq = .5, SigmaSqNoise= -9, EFlag = FALSE,
     logp = -1, log1mp = -1, SigmaNoiseInv = -1, NOPRINT = TRUE, n1 = -9, n2 = -9) { 
     XX <- xxs; YY <- yys;  puse <- ppiuse; n1 <- TotalRuns;
     if (n2 < 0) {n2 <- NumEMConv;}
     if (n1 < 0) { n1 <- TotalRuns; }
     if (SigmaSqNoise  == -9) {SigmaSqNoise <- sigmaNoiseSq; }
     if(!is.loaded("CalculateGFTAllShell") ) {
        LoadTwoCC();
     }
     if (tauDsqStart == -999) {
        XEffect <- sd(as.vector(YY)) / mean(apply(XX,2,sd));
        tauDsqStart = (.5 * XEffect)^2;
        tauAsqStart = (2 * XEffect)^2;
        SigmaSqNoise = (.25 * sd(as.vector(YY)))^2;
        SigmaNoiseInv = 1 / SigmaSqNoise;
     }
     if (tauDMult == -999) {
        tauDMult = 2^(-1/8);
     }
     if (tauAMult == -999) {
        tauAMult = 2^(1/8);
     } 
     
     if (logp == -1) { logp = log(ppiuse); }
     if (log1mp == -1) { log1mp = log(1-ppiuse); }
     if (SigmaNoiseInv == -1) { SigmaNoiseInv = 1 / SigmaSqNoise; }
     BStartBBs <- puse + 0 * (1:length(XX[1,]));
     BStartBetas <- 0 * (1:length(XX[1,]));
  if (det(t(XX) %*% XX) > .000001) {
       BStartBetas <- solve(t(XX) %*% XX) %*% t(XX) %*% YY;
       BStartBBs <- 0 * (1:length(XX[1,])) + ppiuse;
  } else  if (tauDsqStart != tauAsqStart) {
xxSt <- sqrt( abs( ( .5 * log(2 * pi * tauDsqStart) - .5 * log(2*pi * tauAsqStart)
    + log(puse) - log(1-puse) ) /
        ( 1 / (2 * tauDsqStart) - 1 / (2 * tauAsqStart) ) ) );
    BStartBetas = sign( t(XX) %*% YY) * xxSt;
    BStartBetas <-1/diag(t(XX) %*% XX) * t(XX) %*% YY
    BStartBBs <- puse + 0 * (1:length(XX[1,]));
  }    else {
    BStartBBs <- puse + 0 * (1:length(XX[1,]));
    BStartBetas <-  0 * (1:length(XX[1,]));
    BStartBetas <-1/diag(t(XX) %*% XX) * t(XX) %*% YY
  }

##kLen = length(BStartBetas);
   kLen = length(XX[1,]);
   NLen = length(YY);
##n1 = 128; 
##n2 = 40;
if (NOPRINT == TRUE) {
  ##print("EMRIDGE: Performing EMRIDGENoRecShell");
RetCEV <- .C("EMRIDGENoRecShell", NLen = as.integer(NLen), kLen=as.integer(kLen),
         YY = as.double(YY), XX = as.double(XX), 
         ChangeNSteps = as.integer(n1), ConvergenceNSteps = as.integer(n2), 
         ConvergeCloseEnough = as.double(.000001), 
         SigmaSqNoise = as.double(SigmaSqNoise),
         ppiuse = as.double(ppiuse), 
         tauDsqStart = as.double(tauDsqStart), 
         tauAsqStart = as.double(tauAsqStart), 
         tauDMult = as.double(tauDMult), tauAMult = as.double(tauAMult), 
         tauAMultStop = as.integer(tauAMultStop), 
         BStartBBs = as.double(BStartBBs),
         BStartBetas = as.double (BStartBetas),
         ReturnBetas = as.double(vector("numeric", kLen *2)),
         BBHatsFinal = as.double(vector("numeric", kLen * 2)),
         logp = as.double(logp), log1mp = as.double(log1mp), 
         SigmaNoiseInv = as.double(SigmaNoiseInv)
      ); 
      RetCEV$ReturnBetas = RetCEV$ReturnBetas[1:kLen];
      RetCEV$BBHatsFinal = RetCEV$BBHatsFinal[1:kLen];
     ##print("EMRIDGE: Finished with .C out to RetCEV");
} else if (EFlag == FALSE) {
##print("EMRIDGE: Performing EMRIDGEShell2 #1");
RetCEV <- .C("EMRIDGEShell2", NLen =as.integer(NLen), kLen=as.integer(kLen), 
         YY = as.double(YY), XX = as.double(XX), 
         ChangeNSteps = as.integer(n1), ConvergenceNSteps = as.integer(n2), 
         SigmaSqNoise = as.double(SigmaSqNoise),
         ppiuse = as.double(ppiuse), 
         tauDsqStart =as.double(tauDsqStart), 
         tausqStart = as.double(tauAsqStart), 
         tauDMult = as.double(tauDMult), sigma2Mult = as.double(tauAMult), 
         tauAMultStop = as.integer(tauAMultStop), 
         BStartBBs = as.double(BStartBBs),
         BStartBetas = as.double (BStartBetas),
         ReturnBetasAllEM = as.double(vector("numeric", kLen * n1 *n2 )),
         ReturnBetas = as.double(vector("numeric", kLen)),
         BBHatsAllEM = as.double(vector("numeric", kLen * n1 * n2)),
         BBHatsFinal = as.double(vector("numeric", kLen)),
         logp = as.double(logp), log1mp = as.double(log1mp),
         SigmaNoiseInv = as.double(SigmaNoiseInv)
      ); 
      RetCEV$ReturnBetasAllEM <- array(RetCEV$ReturnBetasAllEM, dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") );
      RetCEV$BBHatsAllEM <- array(RetCEV$BBHatsAllEM, dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") ); 
} else {
##print("EMRIDGE: Performing EMRIDGEShell2 #2");
RetCEV <- .C("EMRIDGEShell2", NLen = as.integer(length(YY)), kLen=as.integer(length(XX[1,])), 
         YY = as.double(YY), XX = as.double(XX), 
         ChangeNSteps = as.integer(n1), ConvergenceNSteps = as.integer(n2), 
         SigmaSqNoise = as.double(SigmaSqNoise),
         ppiuse = as.double(ppiuse), 
         tauDsqStart =as.double(tauDsqStart), 
         tauAsqStart = as.double(tauAsqStart), 
         tauDMult = as.double(tauDMult), tauAMult = as.double(tauAMult), 
         tauAMultStop = as.integer(tauAMultStop), 
         BStartBBs = as.double(BStartBBs),
         BStartBetas = as.double (BStartBetas),
         ReturnBetasAllEM = as.double(vector("numeric", kLen * (n1+2) *(n2+2) )),
         ReturnBetas = as.double(vector("numeric", kLen+5)),
         BBHatsAllEM = as.double(vector("numeric", kLen * (n1+2) * (n2+2))),
         BBHatOutVec = as.double(vector("numeric", kLen +5)),
         logp = as.double(logp), log1mp = as.double(log1mp), 
         SigmaNoiseInv = as.double(SigmaNoiseInv)
      ); 
      RetCEV$ReturnBetas = RetCEV$ReturnBetas[1:kLen];
      RetCEV$BBHatsFinal = RetCEV$BBHatsFinal[1:kLen];
      RetCEV$ReturnBetasAllEM <- array(RetCEV$ReturnBetasAllEM[1:(kLen * n2 * n1)], dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") );
      RetCEV$BBHatsAllEM <- array(RetCEV$BBHatsAllEM[1:(kLen * n2 * n1)], dim = c(kLen, n2, n1), dimnames = 
                                      c("k ", "ConvIter", "ChangeIter") ); 
} 
  RetCEV$tauAK <- tauAsqStart * tauAMult^(0:(n1-1));
  if (tauAMultStop < n1 && tauAMultStop > 0) {
      RetCEV$tauAK[(tauAMultStop+1):n1] <- RetCEV$tauAK[tauAMultStop];
  }
  RetCEV$tauDK <- tauDsqStart * tauDMult^(0:(n1-1));
  RetCEV$Type = "EMRIDGE";
  if (RetCEV$BStartBetas[1] == -999) {
     print("RetCEV  Signals RIDGE Failure \n");
     RetCEV$FailureFlag = 1;
  } else {
     RetCEV$FailureFlag = 0;
  }
  return(RetCEV);    
}



####################################################################
####  Lars Prototype
####
####     This is a R-Coded prototype for the LARS maximization algoirthm.
####     We can choose what our BetStart should be based upon whether
####        v1(BetStart[jj]) "< or >"  v2(BetStart[jj])
####     Keep the BetStart that is a smaller initial penalty. 
LarsPrototype <- function(xxs,yys, BetStart = -1, lambda = 0, BetOlds = -1)  {
     if(!is.loaded("CalculateGFTAllShell") ) {
        LoadTwoCC();
     }
   if ( length(BetStart) ==1 && BetStart == -1) {
       OnBetas <- xxs[1,] * 0;
   } else {
       OnBetas <- BetStart;
   }
   BetaOlds = BetOlds;
   if ( length(BetOlds) == 1 && BetOlds == -1) {
       BetOlds <- xxs[1,] * 0;
       BetaOlds <- xxs[1,] * 0;
   }
   OnVecAW <- vector("numeric", length(OnBetas));

     LarsKeep <- list(USReturnBetasAll = matrix(0, length(xxs[1,])+1, length(xxs[1,])), 
                      USReturnBetas = vector("numeric", length(xxs[1,])),                    
                      CCMaxsKeep = vector("numeric", length(xxs[1,])), 
                      KeepJs = vector("numeric", length(xxs[1,])),
                      Penalty = matrix(0, length(xxs[1,]) + 1, 4) ,
                      AlterVec = matrix(0, length(xxs[1,]), length(xxs[1,])),
                      GAllKeep = vector("numeric", length(xxs[1,])),
                      GammasKeep = vector("numeric", length(xxs[1,])),
                      OnBetas = OnBetas,
                      OnVecAW = OnVecAW,
                      meanXXs = meanC(xxs), meanYYs = meanC(yys),
                      sdXXs = apply(xxs,2,sd), sdYYs = sd(as.vector(yys)),
                      ReturnBetas = 1,
                      ReturnBetasAll = matrix(0, length(xxs[1,])+1, length(xxs[1,]))
                      );
   ##yys <- StandardizeXX(yys);
   ##xxs <- StandardizeXX(xxs);
   xxs = matrix(as.vector(xxs), length(xxs[,1]), length(xxs[1,]));   
   xxs <- t( t(xxs) - meanC(xxs));
   xxs = matrix(as.vector(xxs), length(xxs[,1]), length(xxs[1,]));
   yys <- matrix(as.vector(yys - mean(yys)), length(yys),1);
   if (length(xxs[,1]) != length(yys[,1])) {
      print(paste("dim(xxs)  =(", dim(xxs)[1],", ", dim(xxs)[2],") and dim(yys)",
         " = (", dim(yys)[1],", ", dim(yys)[2], ")", sep="") );
         print("xxs is");
      print(xxs);
         print("yys is");
      print(yys);
      return(1);
   }
   LarsKeep$USReturnBetasAll[1,] <- OnBetas;
   onmuA <-  yys * 0;
   OnJ <- 1;
   epsilonBM <- .0002;
   v1 = lambda;
   v2=0;
   
   LarsKeep$Penalty[1,] = c(PenaltyFunction(yys, OnBetas, xxs, v1, v2, BetOlds),
     sum( (as.vector(yys) - as.vector(xxs %*% OnBetas))^2),
     sum( v1 * abs(OnBetas)),
     sum( v1 * abs( BetOlds - OnBetas ))
     );
          Marker <- NULL;
                               
   for (tt in 1:length(xxs[1,])) {
   OnResid = as.vector(yys) - as.vector(xxs %*% OnBetas);
    if ( dim(t(xxs))[2] != dim( 
      as.matrix( matrix((as.vector(yys)- as.vector(onmuA)), length(yys),1))
        )[1]
         ){
        print("The WORLD Has Gone insane Here is yys!");
        print(yys);
        print(" And Here is xxs!");
        print(xxs);
   }
   cc <- as.matrix(matrix(t(xxs), length(xxs[1,]), length(xxs[,1]))) %*% 
        as.matrix( matrix((as.vector(yys)- as.vector(onmuA)), length(yys),1));
   CCxxix <- sort( abs(cc), decreasing = TRUE, index=TRUE);
   OnJ <- length(cc[abs(abs(cc) - abs(CCxxix$x[1])) < epsilonBM]);
   AAix <- sort( CCxxix$ix[1:OnJ]);
     if (length(Marker) > 0) {
       for (tt in 1:length(Marker)) {
           AAix <- AAix[ AAix != Marker[tt] ];
       }
     }
     OnJ <- length(AAix);
     CCxx <- cc[AAix];
   LarsKeep$KeepJs[tt] <- OnJ;
   LarsKeep$CCMaxsKeep[tt] <- abs(max(CCxx));
   AAixC <- sort(CCxxix$ix[(OnJ+1):length(CCxxix$x)]);
   AAixC <- c(AAixC, Marker); AAixC <- sort(unique(AAixC));
   sjA <- sign(cc[AAix]);
   OnVecA <- (1:length(OnBetas)) * 0;
   OnVecA[ AAix ] = sjA
   
   XXA <-t( t( xxs[, AAix]) * sjA);
   GGA <- t(XXA) %*% XXA;
   GGAs <- solve(GGA);
   AAA <- 1/ sqrt( sum(GGAs));
   wa <- AAA * GGAs %*% matrix(1,length(AAix),1);
   ua <- XXA %*% wa;
   aaa <- t(xxs) %*% ua;
   OnVecAW <- ( 1:length(OnBetas)) * 0;
   OnVecAW[AAix] <- wa * sjA;
 
   gammaList1 <- (max(abs(CCxx)) - cc) / (AAA - aaa);
   gammaList2 <- (max(abs(CCxx))+cc) / (AAA+aaa);
   ZeroEpsilon <- .000001
   gammaList1[ (max(abs(CCxx))-cc) < ZeroEpsilon ] <- -9;
   gammaList2[ (max(abs(CCxx))+cc) < ZeroEpsilon ] <- -9; 
   gammaList1[ (max(AAA) - aaa) < ZeroEpsilon] <- -8;
   gammaList2[ (max(AAA)+aaa) < ZeroEpsilon] <- -8;   
   gF <- min( gammaList1[ is.na(gammaList1) == FALSE & gammaList1 > ZeroEpsilon], 
              gammaList2[ is.na(gammaList2) == FALSE & gammaList2 > ZeroEpsilon] );
   OnResid = as.vector(yys) - as.vector(xxs %*% OnBetas);    
   GLD <- .C("CalculateGFTAllShell", MyAns = as.double(0.0),
          TFLV2 = as.double(0.0), GFT0 = as.double(0.0), GFT1 = as.double(0.0),
          OnL = as.integer(length(wa)), TotK = as.integer(length(OnBetas)), 
          ActiveOns = as.integer((AAix)-1), 
          OnBetas = as.double(OnBetas), BetaOlds = as.double(BetaOlds), 
          NLen = as.integer(length(OnResid)), OnResid = as.double( OnResid), 
          wa = as.double(wa), sjA = as.integer(sjA), ua = as.double(ua), lambda = as.double(v1), 
          v2 = as.double(v2), OverallMax = as.double(0.0) );
   GFT2 = GLD$MyAns;  
   TFLV2 = GLD$TFLV2;
   GFT0 = GLD$GFT0;
   GFT1 = GLD$GFT1;        
   LarsKeep$GAllKeep[tt] <- GFT2;
   LarsKeep$USReturnBetasAll[tt+1,] <- OnBetas
   LarsKeep$AlterVec[tt,] <- OnVecAW;
   LarsKeep$GammasKeep[tt] = gF;
   if (GFT2 < 0) {
       LarsKeep$GammasKeep[tt] = 0;
   } else  if (GFT2 != 0 && (abs(GFT2)  < gF)) {         
      LarsKeep$GammasKeep[tt]  =  min(gF, GFT2);
   } else{
       LarsKeep$GammasKeep[tt] = gF;
   }
   if (GFT2 < 0) {
        print(paste("Taking negative GFT2 = ", GFT2, sep="") );
        gF = 0;
        OnVecAW <- sign(LarsKeep$BetasKeep[tt,]);
        OnBetas <- OnBetas;
        for (ttL in tt:length(xxs[1,])) {
             LarsKeep$BetasKeep[ttL+1,] <- OnBetas;
        }
        LarsKeep$GammasKeep[tt] = 0;
        print("End");
        break;
   } else if ((GFT2 != 0) && (abs(GFT2) < gF)) {
        print(paste("Taking GFT2 = ", GFT2, sep="") );
        gF <- GFT2;
        OnVecAW <- (1:length(OnBetas)) * 0;
        OnVecAW[ AAix ] <- sjA * wa;
        Marker <- (1:length(OnBetas))[OnBetas != 0 & 
                              abs(OnBetas + gF * OnVecAW) < epsilonBM];
        OnBetas <- OnBetas + gF * OnVecAW;           
        onmuA <- onmuA + gF * ua;
        for (ttL in tt:length(xxs[1,])) {
           LarsKeep$BetasKeep[ttL+1,] <- OnBetas;
        }
        LarsKeep$GammasKeep[tt] = gF;
        print("End");
        break;
   } else {
       LarsKeep$GammasKeep[tt] = gF;
        OnVecAW <- (1:length(OnBetas)) * 0;
        OnVecAW[ AAix ] <- sjA * wa;
        OnBetas <- OnBetas + gF * OnVecAW;  
        LarsKeep$USReturnBetasAll[tt+1,] = OnBetas;         
        onmuA <- onmuA + gF * ua;
   }

   LarsKeep$Penalty[tt+1,] = c(PenaltyFunction(yys, OnBetas, xxs, v1, v2, BetOlds),
     sum( (as.vector(yys) - as.vector(xxs %*% OnBetas))^2),
     sum( v1 * abs(OnBetas)),
     sum( v1 * abs( BetOlds - OnBetas ))
     );
   if (LarsKeep$Penalty[tt,1] < LarsKeep$Penalty[tt+1,1] ) {
        LarsKeep$GammasKeep[tt] = 0;
        LarsKeep$OnVecAW <- LarsKeep$AlterVec[tt-1,]
        LarsKeep$OnBetas <- LarsKeep$BetasKeep[tt,];
        OnVecAW <- LarsKeep$OnVecAW;
        OnBetas <- LarsKeep$OnBetas;
        print("End");
   }

  }                           
     LarsKeep$ReturnBetas = LarsKeep$OnBetas * LarsKeep$sdYYs / LarsKeep$sdXXs;
     LarsKeep$ReturnBetasAll = t( t(LarsKeep$USReturnBetasAll)  * LarsKeep$sdYYs / LarsKeep$sdXXs );
       LarsKeep$OnVecAW <- OnVecAW;
     LarsKeep$USReturnBetas <- OnBetas;
  LarsKeep$Type <- "Lars"
  return(LarsKeep);
} 


#########################################################################
####   LarsCC2
####    A function using CC code to run Lars Lasso Algorithm
####     Makes use of the multiple version (weighted/unweighted, O(nk) or O(k^2)
####     coded version in LarsObject.cc
####
#####   The user is free to choose BetStart so that decision is not made
#####    by this algorithm
 LarsCC2 <- function(xxs,yys, BetStart = -1, lambda = 0, BetOlds = -1,
     StandardFlag = 0, Weights = NULL, GFlag = -1, PrintGFlag = FALSE,
     WLSWeights = -1)  {
     
     DEFCONST = 3.5
     if(!is.loaded("LarsLassoAlgorithmShell2") ) {
        LoadTwoCC();
     }

   if ( length(BetStart) ==1 && BetStart == -1) {
       OnBetas <- xxs[1,] * 0;
   } else if (length(BetStart) < length(xxs[1,])) {
       print("BetStart is flawed can't use to trigger OnBetas");
       OnBetas <- xxs[1,] * 0;
   ##} else if (length(BetStart) == length(xxs[1,])) {
   ##  OnBetas <- BetStart;
   ##  OnBetas[ v1 * abs(BetStart) < v2 * abs(BetStart) ] <- 0;
   ##  #### We probably have to think more about this condition
   } else {   
       OnBetas <- BetStart;
   }
   if ( length(BetOlds) == 1 && BetOlds == -1) {
       BetOlds <- xxs[1,] * 0;
   }
   if (GFlag <= -1) {
	   if (length(xxs[1,]) < length(xxs[,1])) {
	        if (PrintGFlag == TRUE) {
	          print("LarsCC2, k Length Advantageous:: Setting GFlag = 1");
	        }
	        GFlag = 1;
	   } else {
	        if (PrintGFlag == TRUE) {
	          print("LarsCC2, Nlength Advantageous:: SettingGFlag = 0");
	        }
	        GFlag = 0;
	   }  
   }
   kLen = length(xxs[1,]);
   if (length(Weights) !=  kLen) {
        Weights = (1:kLen) * 0 + 1.0;
        WFlag = 0;
   } else {
        WFlag = 1;
   }
   GFFlag = GFlag  + 2* WFlag;
   ##print(paste("GFlag = ", GFlag, ",  WFlag = ", WFlag, ", GFFlag = ", GFFlag, sep=""));
   KFlag = 0;
   sdXXs <- apply(xxs,2,sd);
   sdYYs <- sd(as.vector(yys));
   RunLength = floor(min(c(length(xxs[1,]), length(yys))) * DEFCONST); 
     LarsKeep <- list(BetasKeep = matrix(0, RunLength+1, length(xxs[1,])),
                      KeepJs = vector("numeric", RunLength+4),
                      Penalty = matrix(0, RunLength+4, 4) ,
                      AlterVec = matrix(0, RunLength+4, length(xxs[1,])),
                      GammasKeep = vector("numeric", RunLength+4),
                      OnBetas = OnBetas, 
                      meanXXs = meanC(xxs), meanYYs = mean(yys),
                      sdXXs = apply(xxs,2,sd), sdYYs = sd(as.vector(yys)),
                      BetasFull = 1, 
                      lambda = lambda       
                      );
   if (StandardFlag == 1) {
      if (PrintGFlag  == TRUE) { print("LarsLasso2: StandardFlag is active"); flush.console(); }
      yys <- StandardizeXX(yys);
      xxs <- StandardizeXX(xxs); 
      BetaOldsM <- BetOlds   * LarsKeep$sdXXs / LarsKeep$sdYYs;  
      OnBetas <- OnBetas * LarsKeep$sdXXs / LarsKeep$sdYYs;
   } else {
      BetaOldsM <- BetOlds
      OnBetas <- OnBetas
   }
   epsilonZero = .00000001
         DEFCONST <- 3.5
   if (length(yys) == length(xxs[1,]) && abs(det(xxs)) < epsilonZero) {
          ##print("We need to do Determinant Correction");
           Cors <-  t(yys) %*% xxs   
          KFlag = 1;
          AASucks <- sort( abs(Cors), index=TRUE, decreasing = FALSE);
          AAElim <- AASucks$ix[1];
          AAKeep <- sort(AASucks$ix[2:length(AASucks$ix)]);
          ##OnBetas <- OnBetas[AAKeep];
          ##xxs <- xxs[,AAKeep];
          kLenK <- length(xxs[1,]) - 1;
          HalfLars <-  LarsCC2(xxs[,AAKeep],yys, 
             BetStart = -1, lambda = 0, BetOlds = -1,
             StandardFlag = 0, Weights = NULL, GFlag = -1, PrintGFlag = FALSE);
          LarsKeep = HalfLars;
          if (!is.null(HalfLars$USWReturnBetas)) {
             LarsKeep$USWReturnBetas = vector("numeric", kLen);
              LarsKeep$USWReturnBetas[AAKeep] = HalfLars$USWReturnBetas;
          }
          if (!is.null(HalfLars$USWReturnBetasAll)) {
             LarsKeep$USWReturnBetasAll = matrix(0,
                 length(HalfLars$USWReturnBetasAll[,1]), kLen);
             LarsKeep$USWReturnBetasAll[,AAKeep] = HalfLars$USWReturnBetasAll;
          }          
          if (!is.null(HalfLars$ReturnBetas)) {
             LarsKeep$ReturnBetas = vector("numeric", kLen);
             LarsKeep$ReturnBetas[AAKeep] = HalfLars$ReturnBetas;
          }
          if (!is.null(HalfLars$ReturnBetasAll)) {
             LarsKeep$ReturnBetasAll = matrix(0,
                 length(HalfLars$ReturnBetasAll[,1]), kLen);
             LarsKeep$ReturnBetasAll[,AAKeep] = HalfLars$ReturnBetasAll;
          } 
          if (!is.null(HalfLars$USReturnBetas)) {
             LarsKeep$USReturnBetas = vector("numeric", kLen);
             LarsKeep$USReturnBetas[AAKeep] = HalfLars$USReturnBetas;
          }
          if (!is.null(HalfLars$USReturnBetasAll)) {
             LarsKeep$USReturnBetasAll = matrix(0,
                 length(HalfLars$USReturnBetasAll[,1]), kLen);
             LarsKeep$USReturnBetasAll[,AAKeep] = HalfLars$ReturnBetasAll;
          }                         
   } else {               
    kLen <- length(xxs[1,])    
      kLenK = kLen;   
          LLS <- .C("LarsLassoAlgorithmShell2",
              yys = as.double(yys), NLen = as.integer(length(yys)), 
              xxs = as.double(as.vector(xxs)), kLen = as.integer(kLenK), 
              InitBetas = as.double(OnBetas), OldBetas = as.double(BetaOldsM),
              lambda = as.double(lambda), 
              USWReturnBetas = as.double(vector("numeric", kLenK)), 
              USWReturnBetasAll = as.double( vector("numeric", (RunLength+2) * (kLenK))
                                   ),
              JKeep = as.integer(vector("numeric", RunLength+2)),
              PenaltyKeep = as.double(vector("numeric", (RunLength + 3)*4)),
              GammasKeep = as.double(vector("numeric", RunLength +3)),
              gFKeep = as.double(vector("numeric", RunLength +3)), 
              GFFlag = as.integer(GFFlag),
              WeightLen = as.integer(length(Weights) ),
              Weights = as.double(Weights),
              Finaltt = as.integer(1), 
              WLSWeights = as.double(WLSWeights)                   
            ); 
     Finaltt = LLS$Finaltt;
     if (length(WLSWeights) > 1) {
        LarsKeep$WLSWeights = LLS$WLSWeights[1:length(yys)];
     }
     LarsKeep$USWReturnBetasAll = t(matrix(LLS$USWReturnBetasAll[1:(kLenK*(Finaltt+1))], kLenK, Finaltt+1 ) );
     LarsKeep$USWReturnBetas = LLS$USWReturnBetas[1:kLenK];
       if (WFlag == 1) {
          LarsKeep$USReturnBetasAll = t( t(LarsKeep$USWReturnBetasAll ) * Weights);
          LarsKeep$USReturnBetas = LarsKeep$USWReturnBetas * Weights;
       } else {
           LarsKeep$USReturnBetas = LarsKeep$USWReturnBetas;
           LarsKeep$USReturnBetasAll = LarsKeep$USWReturnBetasAll;                   
       }
     LarsKeep$ReturnBetasAll = LarsKeep$USReturnBetasAll;
     LarsKeep$ReturnBetas = LarsKeep$USReturnBetas;
     LarsKeep$OnBetas  = LLS$USWReturnBetas[1:kLenK];   
     LarsKeep$KeepJs = LLS$JKeep[1:Finaltt];
     LarsKeep$Penalty = t(matrix(LLS$PenaltyKeep[1:(4*(Finaltt+1))], 4, Finaltt+1));
     if (StandardFlag == 1) {
        if (PrintGFlag == TRUE) { print("DeStandardizing"); flush.console(); }
        LarsKeep$ReturnBetas = LarsKeep$USReturnBetas * LarsKeep$sdYYs / LarsKeep$sdXXs;
        LarsKeep$ReturnBetasAll = t( t(LarsKeep$USReturnBetasAll)  * LarsKeep$sdYYs / LarsKeep$sdXXs );
     } else {
        LarsKeep$ReturnBetas= LarsKeep$USReturnBetas
        LarsKeep$ReturnBetasAll = LarsKeep$USReturnBetasAll;
     }
     LarsKeep$InitBetas = LLS$InitBetas;
     LarsKeep$GammasKeep[1:Finaltt] = LLS$GammasKeep[Finaltt];
     LarsKeep$gFKeep = LLS$gFKeep;
         if (length(Weights) != kLen) {
            OnBetas= LLS$USWReturnBetas[1:kLenK];
          } else {
          OnBetas = LLS$USWReturnBetas[1:kLenK] * Weights[1:kLenK];
         }       
     LarsKeep$Type <- "Lars"
  }
  if (StandardFlag == 1) {
     yys = yys * sdYYs
     xxs = t(sdXXs * t(xxs)  );
  }
  if (is.matrix(LarsKeep$USWReturnBetasAll)) {
    LarsKeep$Lambdas = IdentifyLambdas(LarsKeep, xxs, yys, Weights);
  }
  LarsKeep$Type <- "Lars"
  return(LarsKeep);
}

IdentifyLambdas <- function(LarsKeep, xxs, yys, Weights = NULL) {
    LambdaList = vector("numeric", length(LarsKeep$USWReturnBetasAll[,1]));
    for (ii in 1:length(LarsKeep$USWReturnBetasAll[,1])) {
	     OnJJ <-  (1:length(LarsKeep$USWReturnBetasAll[ii,])
                                   )[ LarsKeep$USWReturnBetasAll[ii,] != 0 ];
	     WhatHaveYou = 
          as.vector(2*t(xxs[,OnJJ])  %*% (yys-xxs[,OnJJ] %*%
             matrix( as.numeric(LarsKeep$USWReturnBetasAll[ii,OnJJ]),
                length(OnJJ),1 )            
              ))/
             as.vector(sign(LarsKeep$USWReturnBetasAll[ii,OnJJ]));
       if (!is.null(Weights)) {
          WhatHaveYou = WhatHaveYou / Weights[OnJJ];
       }
        MDD <-  WhatHaveYou[is.na(WhatHaveYou) == FALSE]
        if (length(MDD) <= 0) {
           LambdaList[ii] = max(LambdaList);
        }  else {
          LambdaList[ii] = median(MDD);
        }
    }

    return(LambdaList);
}

###################################################################
##  M2Lasso
##      Creates the "M"-type Max-LambdaD, Min-LambdaA fit
##
OldM2Lasso <- function (xxs=-1, yys=-1, ppiuse = .5, noise = -1,
  FixK1 = -100, sigmaNoiseSq = -1, NLen = -1, XTX = -1, XTY =-1,
  R = 0, StartBeta = -1) {
  if (length(xxs) == 1 && length(XTX) == 1) {
    print("M2Lasso -- Bad Entry, XTX and xxs are unincluded");
    return(-1);
  }
  if (length(xxs) == 1 && length(XTX) > 1) {
    if (NLen == -1) {
      print("M2Lasso -- must supply me with Nlen for XTX");
    }  
    if (sigmaNoiseSq < 0) {
      print("M2Lasso -- must supply me with sigma for XTX");
    }
    kLen = length(XTX[1,])
    USENOISE <- sigmaNoiseSq;

    SDD <- sqrt(min(apply(xxs,2,sd))^2 * (NLen-1) );
    LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD, 
      sqrt(3.7 / USENOISE) * SDD *  NLen);
    SABS <- min( abs(SMS$BetasReal[SMS$BetasReal != 0])  );
    SDA <- SDD / sqrt(NLen-1);
    LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) *
      sqrt( 2 * SDA / USENOISE)
    LambdaASeq = c( 1.1 * LambdaAS1, .05 * LambdaAS1);
    if (LambdaASeq[1] >= LambdaDSeq[2]) {
      LambdaASeq[1] = LambdaDSeq[1] + LambdaASeq[1];
      LambdaDSeq[1] = LambdaASeq[1] - LambdaDSeq[1];
      LambdaASeq[1] = LambdaASeq[1] - LambdaDSeq[1];
    }
    LambdaASeq = c(LambdaASeq[1], LambdaASeq);
    LambdaDSeq = c(LambdaDSeq[1], LambdaDSeq);        
    OrderSeq = c(20, 4,1);
    if (R > 0 && R <= .25) {
      SigmaVec = c(USENOISE * NLen^R, USENOISE, USENOISE);
      LambdaASeq = c(LambdaASeq[1] * N^(.5 * R), LambdaASeq[1], 
        LambdaASeq[length(LambdaASeq)]);
      LambdaDSeq = c(LambdaDSeq[1] * N^(-.5 * R), LambdaDSeq[1],
        LambdaDSeq[length(LambdaDSeq)] );
    } else {
      SigmaVec = -1;
    }   
    HitBatWingFast <- EM2Lasso(XTX=XTX, XTY=XTY, NLen = NLen, ppiuse = ppiuse, 
      FixKa = FixKa, sigmaNoiseSq = USENOISE,
      RatWant = .2, StlambdaD=LambdaDSeq[1], 
      StlambdaA = LambdaASeq[1], LambdaDK = LambdaDSeq, LambdaAK = LambdaASeq, 
      OrderSeq = OrderSeq, TotalRuns = 5, 
      NumEMConv = 4, MultEMCons = .99, 
      NumCDOConv = 40, CDOEpsilon = .000001,
      SigmaVec = SigmaVec, StartBeta = StartBeta); 
   return(HitBatWingFast); 
  }
  if (sigmaNoiseSq < 0) {
    if (length(yys) > 1) {
      sigmaNoiseSq <- .5 * var(yys);
    } else {
      sigmaNoiseSq = .5;
    } 
  }
  if (length(dim(xxs)) == 2) {
    kLen = length(xxs[1,])
  }
  NLen = length(yys);
  USENOISE <- sigmaNoiseSq; kLen <- max(length(XTY), length(xxs[1,]) );

	SDD <- sqrt(min(apply(xxs,2,sd))^2 * (length(xxs[,1])-1) );
	LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD, 
    sqrt(3.7 / USENOISE) * SDD *  NLen);
	SABS <- min( abs(SMS$BetasReal[SMS$BetasReal != 0])  );
	SDA <- SDD / sqrt(NLen-1);
	LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) * 
    sqrt( 2 * SDA / USENOISE)
  LambdaASeq = c( 1.1 * LambdaAS1, .05 * LambdaAS1);
  if (LambdaASeq[1] >= LambdaDSeq[2]) {
    LambdaASeq[1] = LambdaDSeq[1] + LambdaASeq[1];
    LambdaDSeq[1] = LambdaASeq[1] - LambdaDSeq[1];
    LambdaASeq[1] = LambdaASeq[1] - LambdaDSeq[1];
  }
  LambdaASeq = c(LambdaASeq[1], LambdaASeq);
  LambdaDSeq = c(LambdaDSeq[1], LambdaDSeq);        
	OrderSeq = c(20, 4,1);
  if (R > 0 && R <= .25) {
    SigmaVec = c(USENOISE * NLen^R, USENOISE, USENOISE);
    LambdaASeq = c(LambdaASeq[1] * N^(.5 * R), LambdaASeq[1], 
      LambdaASeq[length(LambdaASeq)]);
    LambdaDSeq = c(LambdaDSeq[1] * N^(-.5 * R), LambdaDSeq[1],
      LambdaDSeq[length(LambdaDSeq)] );
  }  else {
    SigmaVec = -1;
  }   	      
  HitBatWingFast <- EM2Lasso(xxs=xxs, yys=yys, ppiuse = ppiuse, 
	  FixKa = FixKa, sigmaNoiseSq = USENOISE,
	  RatWant = .2, StlambdaD=LambdaDSeq[1], 
	  StlambdaA = LambdaASeq[1], LambdaDK = LambdaDSeq,
    LambdaAK = LambdaASeq, OrderSeq = OrderSeq, 
    TotalRuns = 5, 
	  NumEMConv = 4, MultEMCons = .99, 
    NumCDOConv = 40, CDOEpsilon = .000001,
    SigmaVec = SigmaVec, StartBeta = StartBeta);
  return(HitBatWingFast);	                         
}

##################################################################
## EM2Lasso: Basic Coding of the EM2Lasso Algorithm
##   EM2Lasso is an algorithm solving the Two-Lasso problem for 
##   a number of specified lambda values
##    (specified by multiplying initial values by a multiplicative constant)
##
##   Does not estimate pi_A or sigma
##
##   By Default uses CoordinateDescent Code from CoordinateDescent.cc
##   However, It can use LARS algorithm in stead from LarsC.cc
EM2Lasso <- function(xxs = -1, yys = -1, NLen = -1, ppiuse = 3/7, sigmaNoiseSq = -999, 
  StlambdaD = -999, StlambdaA = -999, lambdaDMultC = 2^(2/8), 
  lambdaAMultC = 2^(-2/8), lambdaAmultStop = 20, TotalRuns = 20, 
  NumEMConv = 4, MultEMCons = .99, StartBeta = -999, StandardFlag = 0,
  RatWant = .2, DConfidenceInt = -1.0, 
  LambdaDK = -999, LambdaAK = -999,
  OrderSeq = -999, NumCDOConv = 50, CDOEpsilon = .0000001, XTX = -1, XTY = -1,
  InverseGammaConstant = 1, FixKa = -100, InitKKs = -5, WLSWeights = -1, TDFNu = -1,
  PrintFlag = -1, Lambda3 = 0, Lambda3Seq = -999,
  SigmaVec = -1) {    
  if(!is.loaded("LarsConvergencyShell") ) {
    LoadTwoCC();
  } 
  if (ppiuse   >= 1) {ppiuse = .99;}
  if (ppiuse <=0) {ppise = .01;}
  if (is.na(ppiuse) || is.nan(ppiuse)) { ppiuse = .5;}

  if (is.matrix(xxs) == FALSE || length(xxs) == 1 || xxs == -1) {
    if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == -1) {
      print("EM2Lasso, Error no proper X given, XTX False"); return;  
    }
    if ( (is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) ||
      length(XTY) == 1 || XTY == -1) {
      print("EM2Lasso, Error no proper X given, XTY False"); return;
    }
    if (length(NLen) > 1 || NLen == -1) {
      print("EM2Lasso, Error, no proper NLen given");
      return(-1);
    }
    XTXFlag = 1;
    XXX = as.double(XTX); YYY = as.double(XTY);
    TNLen = - NLen; kLen = length(XTX[1,]);
    WLSWeights = -1;  TDFNu = -1;
  } else if (( is.vector(yys) == FALSE && is.matrix(yys) == FALSE) || 
    length(yys) == 1 || yys == -1) {
    print("EM2Lasso, Error no proper Y given"); 
    print("   yys is : ");
    print( yys );
    return(-1);
  } else {
    XTXFlag = 0;
    if (length(yys) > 1000 && length(xxs[1,]) < length(yys) &&
        ( length(WLSWeights) != length(yys) ||
          length(TDFNu) == 1 && TDFNu > 0  )
         ) {
      XTXFlag = 1;
      XXX = as.double( t(xxs)%*% xxs);
      YYY = as.double( t(xxs) %*% yys);
      TNLen = -length(yys); kLen = length(xxs[1,]);     
    }  else {
      XXX = as.double(xxs); YYY = as.double(yys);
      TNLen = length(yys); kLen = length(xxs[1,]);
      if (length(WLSWeights) != length(yys)) {
        WLSWeights = -1;
      }
      if (kLen > 10000) {
        InitKKs = round(kLen *ppiuse);
        if (InitKKs > kLen) {InitKKs = kLen;}
      }
      if (length(TDFNu) == 1 && TDFNu > 0) {
        WLSWeights = rep(1, length(yys))
      } else {
        TDFNu = -1;
      }
    }      
  }   
   if (InitKKs < 0 && kLen > 1000) {
         InitKKs = min(round(kLen * 1.5 * ppiuse), kLen-1);
   } 
          
   if (sigmaNoiseSq < 0) {
      sigmaNoiseSq = (.4 * sd (yys) / mean(apply(xxs,2,sd)))^2;
   }
   if (LambdaAK[1] == -999) {
      if (StlambdaA < 0 & ppiuse > 0 && ppiuse < 1.0) {
        StlambdaA =  .6 * log((1-ppiuse)/ppiuse);      
      } else if (StlambdaA < 0 && FixKa > 0) {
          pDtry <- FixKa / kLen; 
          StlambdaA =  .6 * log((1-pDtry)/pDtry);
      } else if (StlambdaA < 0) {
         print("Why did you not give us StlambdaA");
      }
   } 
   if (LambdaDK[1] == -999) {
      if (StlambdaD < 0 & ppiuse > 0 && ppiuse < 1.0) {
        StlambdaD =  .6 * log((1-ppiuse)/ppiuse);      
      } else if (StlambdaD < 0 && FixKa > 0) {
          pDtry <- FixKa / kLen; 
          StlambdaA =  .6 * log((1-pDtry)/pDtry);
      } else if (StlambdaD < 0) {
         print("Why did you not give us StlambdaD");
      }
   }    

   if (InverseGammaConstant < 0) {
      print(paste("InverseGammaConstant is less than zero is : ", InverseGammaConstant, sep=""));
      return(-1);
   }
   if (length(StartBeta) < kLen) {
     StartBeta = (1:kLen) * 0;
   }  

   sdYYs <- sd(as.vector(yys));
   sdXXs <- apply(xxs,2,sd);   
   if (StandardFlag == 1) {
        yys = StandardizeXX(yys);
        xxs = StandardizeXX(xxs);
        sigmaNoiseSq = sigmaNoiseSq  / sdYYs^2;
        StlambdaD = StlambdaD * sdYYs;
        StlambdaA = StlambdaA * sdYYs;
        ##print(paste("We Standardized Everything, sdYYs =", 
        ##  round(sd(as.vector(yys))), ",  sdXXs = ", round(mean(sd(xxs))),
        ##       ",  sigmaNoiseSq = ", sigmaNoiseSq, sep=""));
   }
    if (InitKKs > kLen) {
      InitKKs = kLen -1;
   } 

      OnBetas = StartBeta;
      LAKSeq = 2;
      if (is.vector(LambdaDK) == FALSE || length(LambdaDK) == 1 || LambdaDK == -999) {
        LambdaDK = vector("numeric", TotalRuns);
        LAKSeq = 1;
        LambdaAK = vector("numeric", TotalRuns);
      } else {
        TotalRuns = length(LambdaDK);
        StlambdaA = LambdaAK[1]; StlambdaD = LambdaDK[1];
      }
      if (is.vector(OrderSeq) == TRUE && length(OrderSeq) >= TotalRuns) {
         LAKSeq = LAKSeq + 2;
      } else if (length(OrderSeq) == 1 && OrderSeq == -999) {
         OrderSeq <- vector("numeric", TotalRuns);
      } else if (is.vector(OrderSeq) == FALSE) {
         OrderSeq <- vector("numeric", TotalRuns);
      } else if (is.vector(OrderSeq) == FALSE || length(OrderSeq) == 1 || OrderSeq != -999) {
         OrderSeq <- vector("numeric", TotalRuns);
      } else {       
         LAKSeq = LAKSeq +2;
      }
      if (length(Lambda3Seq) < TotalRuns) {
         Lambda3Seq = -999;
      }
    if (InitKKs > kLen) {
      InitKKs = kLen -1;
   } 
   if(length(SigmaVec) < TotalRuns) {
       SigmaVec = -1;
   }
   ##print(paste("LambdaAK = ", paste(LambdaAK, collapse=", "), " and SigmaVec = ", 
   ##  paste(SigmaVec, collapse=", "), sep=""));  flush.console();
     LCS =.C("LarsConvergencyShell", YYY = as.double(YYY), NLen = as.integer(TNLen),
              XXX= as.double(XXX), kLen = as.integer(kLen), 
              InitBetas = as.double(OnBetas[1:kLen]), OldBetas=as.double(StartBeta),
              USReturnBetas=as.double(vector("numeric", kLen+3)), 
              USReturnBetasAll =as.double(vector("numeric", kLen * TotalRuns)),
              NumTotalRuns=as.integer(TotalRuns), NumEMConv=as.integer(NumEMConv),
              MultEMCons=as.double(MultEMCons), StlambdaD=as.double(StlambdaD),
              StlambdaA=as.double(StlambdaA), lambdaDMultC=as.double(lambdaDMultC),
              lambdaAMultC=as.double(lambdaAMultC),
              lambdaAmultStop=as.integer(lambdaAmultStop), ppiuse = as.double(ppiuse),
              sigmaNoiseSq = as.double(sigmaNoiseSq),
              BBHatsAll=as.double(vector("numeric", kLen * TotalRuns)),
              NusKeep = as.double(vector("numeric", kLen * (TotalRuns+1))),
              LambdaDK = as.double(LambdaDK),
              LambdaAK = as.double(LambdaAK),
              RecordPostPBj = as.double(vector("numeric", (TotalRuns+1) * kLen) ), 
              DConfidenceInt = as.double(DConfidenceInt),             
              RecordConfidenceInts = as.double(vector("numeric", (TotalRuns+1) * kLen * 2) ), 
              LAKSeq = as.integer(LAKSeq), OrderSeq = as.integer(OrderSeq),
              NumCDOConv = as.integer(NumCDOConv), CDOEpsilon = as.double(CDOEpsilon),
              InitKKs = as.integer(round(InitKKs)),
              InverseGammaConstant = as.double(InverseGammaConstant),
              FixKa = as.integer(FixKa), 
              PiRecVec = as.double( rep(FixKa / kLen, TotalRuns * 40)),
              WLSWeights = as.double(WLSWeights), TDFNu = as.double(TDFNu),
              PrintFlag = as.integer(PrintFlag),
              SigmaVec = as.double(SigmaVec),
              Lambda3 = as.double(Lambda3), Lambda3Seq = as.double(Lambda3Seq)
            );
        LCS$BBHatsAll=t(matrix(LCS$BBHatsAll[1:(kLen*TotalRuns)],kLen,TotalRuns));
        LCS$USReturnBetasAll=t(matrix(LCS$USReturnBetasAll[1:(kLen*TotalRuns)],kLen,TotalRuns));  
        LCS$NusKeep=t(matrix(LCS$NusKeep[1:(kLen*TotalRuns)],kLen,TotalRuns)); 
        LCS$RecordPostPBj=t(matrix(LCS$RecordPostPBj[1:(kLen*TotalRuns)], kLen, TotalRuns) );
        RecordConfidenceInts =  array( dim=c(kLen,2, TotalRuns), dimnames=c("Bj", "LjRj", "Runs"));
        for (ii in 1:TotalRuns) {
           RecordConfidenceInts[,,ii] = 
                t(matrix(LCS$RecordConfidenceInts[ ((ii-1) * TotalRuns) + 1:(kLen*2)],2,kLen));
        }

        EMLarsKeep <- list( NLen = LCS$NLen, kLen = LCS$kLen, 
           USReturnBetas = LCS$USReturnBetas[1:kLen],  USReturnBetasAll = LCS$USReturnBetasAll, 
           NusKeep = LCS$NusKeep,
           NumTotalRuns = LCS$NumTotalRuns, NumEMConv = LCS$NumEMConv, 
           MultEMCons = LCS$MultEMCons, LambdaDK = LCS$LambdaDK, LambdaAK = LCS$LambdaAK,
           lambdaDMultC = LCS$lambdaDMultC, lambdaAMultC = LCS$lambdaAMultC, 
           lambdaAmultStop = LCS$lambdaAmultStop, ppiuse = LCS$ppiuse, 
           USsigmaNoiseSq = LCS$sigmaNoiseSq, 
           sigmaNoiseSq = LCS$sigmaNoiseSq,
           sdYYs = sdYYs, sdXXs = sdXXs,
           BBHatsAll = LCS$BBHatsAll,
           BBHatsFinal = LCS$BBHatsAll[TotalRuns,],
           ReturnBetasAll = 2, ReturnBetas = 2,
           RecordPostPBj=LCS$RecordPostPBj, RecordConfidenceInts = RecordConfidenceInts,
           ppiuse = LCS$ppiuse, LAKSeq = LCS$LAKSeq, OrderSeq=LCS$OrderSeq,
           NumCDOConv = as.integer(NumCDOConv), CDOEpsilon = as.double(CDOEpsilon),
           InverseGammaConstant = as.double(InverseGammaConstant),
           XTX = XTX, XTY = XTY, xxs = xxs, yys = yys,
           FailureFlag = 0, FixKa = FixKa,
           PiRecVec = LCS$PiRecVec[1:TotalRuns], TDFNu = as.double (TDFNu),
           SigmaVec = SigmaVec );
        if (length(WLSWeights) > 1) {
          EMLarsKeep$WLSWeights = LCS$WLSWeights[1:length(yys)];
        }
        if (LCS$USReturnBetas[1] == -999) {
            Ontt1 = round(LCS$USReturnBetas[2]);
            Ontt2 = round(LCS$USReturnBetas[3]);
        } else {
            Ontt1 = TotalRuns;
            Ontt2 = NumEMConv;
        }
        
        if (StandardFlag == 1) {
          EMLarsKeep$ReturnBetasAll = LCS$USReturnBetasAll;
          EMLarsKeep$ReturnBetasAll[1:Ontt1, 1:kLen] = t( t(
                      LCS$USReturnBetasAll[1:Ontt1, 1:kLen]  * sdYYs / sdXXs));
          if (LCS$USReturnBetas[1] != -999) {
                EMLarsKeep$ReturnBetas = LCS$USReturnBetas[1:kLen] * sdYYs / sdXXs;
          } else {
                EMLarsKeep$ReturnBetas = LCS$USReturnBetas[1:kLen] * 0 - 999;
          }
            EMLarsKeep$sigmaNoiseSq = sigmaNoiseSq * (sdYYs)^2; 
        } else {
           EMLarsKeep$ReturnBetasAll = LCS$USReturnBetasAll;
	          if (LCS$USReturnBetas[1] != -999) {
	                EMLarsKeep$ReturnBetas = LCS$USReturnBetas[1:kLen];
	          } else {
	                EMLarsKeep$ReturnBetas = LCS$USReturnBetas[1:kLen] * 0 - 999;
	          }       
        }
     if (LCS$InitBetas[1] == -999) {
        EMLarsKeep$FailureFlag = 1;
     }
     EMLarsKeep$NusKeep = NULL;
     EMLarsKeep$Type <- "EM2Lasso";
     return(EMLarsKeep);              
} 


##################################################################################
##  EMPI2Lasso
##
##   A convergent solution to the Two-Lasso problem for a sequence
##   of growing lambda1 (lambdaD) and shrinking lambda2 (lambdaA) vlaues
##   pi_A is estimated along with sigma (at the user's request, by setting
##    SigmaEta and SigmaBarSq to something other than -1, -999.
##
##   m1 and m2 are parameters for the Beta prior on piA
##   SigmaEta and SigmaBarSq are Gamma parameters for the sigmaSq noise
EMPI2Lasso <- function(xxs=-1, yys=-1, NLen = -1,
          StlambdaD = -999, StlambdaA = -999, lambdaDMultC = 2^(2/8), 
          lambdaAMultC = 2^(-2/8), lambdaAmultStop = 20, TotalRuns = 20, 
          NumEMConv = 4, MultEMCons = .99, StartBeta = -999, StandardFlag = 0,
          RatWant = .2,
          pPpiuseSt = -999, sigmaNoiseSqSt = -999, 
          m1 = 1, m2 = 1, NumEConv = 6, SigmaEta = -1, SigmaBarSq = - 999, 
          DConfidenceInt = -1.0, NumCDOConv = 50, CDOEpsilon = .000001, 
          LambdaAK = -999, LambdaDK = -999, OrderSeq = -999, XTX = -1, XTY =-1,
          InitKKs = -5,
          InverseGammaConstant = 1, WLSWeights = -1, TDFNu = -1, PrintFlag = -1,
          SigmaMultVec = -1) {
   if (pPpiuseSt > 0 && pPpiuseSt < 1) {
      ppiuse <- pPpiuseSt;
   } else {
      ppiuse <- m1 / ( m1 + m2);
   }
     if(!is.loaded("LarsConvergencyPIEstShell") ) {
        LoadTwoCC();
     } 
 if (is.matrix(xxs) == FALSE || length(xxs) == 1 || xxs == -1) {
       if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == -1) {
          print("EMPI2LassoS, Error no proper X given, XTX False"); return;  
       }
       if ( (is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) || length(XTY) == 1 || XTY == -1) {
          print("EMPI2Lasso, Error no proper X given, XTY False"); return;
       }
       if (length(NLen) > 1 || NLen == -1) {
           print("EMPI2Lasso, Error, no proper NLen given");
           return;
       }
       XTXFlag = 1;
       XXX = as.double(XTX); YYY = as.double(XTY);
       WLSWeights = -1;  TDFNu = -1;
       TNLen = - NLen; kLen = length(XTX[1,]);
   } else if ((is.vector(yys) == FALSE && is.matrix(yys) == FALSE) || length(yys) == 1 || yys == -1) {
        print("EMPI2Lasso, Error no proper Y given"); 
        print( " The yys is ");
        print(yys);
        return(-1);
   } else {
       XTXFlag = 0;
       XXX = as.double(xxs); YYY = as.double(yys);
       TNLen = length(yys); kLen = length(xxs[1,]);
       if (length(WLSWeights) != length(yys)) {
          WLSWeights = -1;
       }
       if (length(TDFNu) == 1 && TDFNu > 0) {
             WLSWeights = rep(1, length(yys));
       }  else {
            TDFNu = -1;
       }
       ##if (InitKKs < 0) {
       ##   print("Cannot do xx version if xx is non matrix\n");
       ##   return;
       ##}       
   }  
   if (InitKKs < 0 && kLen > 1000) {
         InitKKs = round( ppiuse * kLen * 1.5);
         if (InitKKs > kLen) {InitKKs = kLen-1;}
   }                    
   if (sigmaNoiseSqSt < 0) {
      sigmaNoiseSqSt = (.4 * sd (yys) / mean(apply(xxs,2,sd)))^2;
   }
   if (SigmaBarSq == -999) {
         SigmaBarSq = sigmaNoiseSqSt;
   }
   if (InverseGammaConstant < 0) {
      print(paste("InverseGammaConstant is less than zero is : ", InverseGammaConstant, sep=""));
      return(-1);
   }
   if (pPpiuseSt == -999)  {
      pPpiuseSt = 4 / length(xxs[1,]);
   }
   if (SigmaEta < 0) {
      SigmaEta <- 3;
   }
   if (StlambdaD < 0 || StlambdaA < 0) {
       StlambdaD =  .6 * log((1-pPpiuseSt)/pPpiuseSt);
       StlambdaA =  .5 * log((1-pPpiuseSt)/pPpiuseSt);      
   }  
   if (length(StartBeta) < kLen) {
     StartBeta = (1:kLen) * 0;
   }  
   sdYYs <- sd(as.vector(yys));
   sdXXs <- apply(xxs,2,sd);   
   if (StandardFlag == 1) {
        yys = StandardizeXX(yys);
        xxs = StandardizeXX(xxs);
        sigmaNoiseSqSt = sigmaNoiseSqSt  / sdYYs^2;
        SigmaBarSq = SigmaBarSq / sdYYs^2;
        StlambdaD = StlambdaD  * sdYYs;
        StlambdaA = StlambdaA * sdYYs;        
        ##yys = yys * sdYYs;
        ##sigmaNoiseSqSt = sigmaNoiseSqSt * sdYYs^2;
        ##StlambdaD = StlambdaD  / sdYYs;
        ##StlambdaA = StlambdaA / sdYYs;
        ##SigmaBarSq = SigmaBarSq * sdYYs^2;
        ##print(paste("We Standardized Everything, sdYYs =", round(sd(as.vector(yys))), ",  sdXXs = ", round(mean(sd(xxs))),
        ##       ",  sigmaNoiseSq = ", sigmaNoiseSq, sep=""));
   }
      LAKSeq = 2;
      if (is.vector(LambdaDK) == FALSE || length(LambdaDK) == 1 || LambdaDK == -999) {
        LambdaDK = vector("numeric", TotalRuns);
        LAKSeq = 1;
        LambdaAK = vector("numeric", TotalRuns);
      } else {
        TotalRuns = length(LambdaDK);
        StlambdaA = LambdaAK[1]; StlambdaD = LambdaDK[1];
      }
      if (is.vector(OrderSeq) == TRUE && length(OrderSeq) >= TotalRuns) {
         LAKSeq = LAKSeq + 2;
      } else if (length(OrderSeq) == 1 && OrderSeq == -999) {
         OrderSeq <- vector("numeric", TotalRuns);
      } else if (is.vector(OrderSeq) == FALSE) {
         OrderSeq <- vector("numeric", TotalRuns);
      } else if (is.vector(OrderSeq) == FALSE || length(OrderSeq) == 1 || OrderSeq != -999) {
         OrderSeq <- vector("numeric", TotalRuns);
      } else {       
         LAKSeq = LAKSeq +2;
      }   
      OnBetas = (1:kLen) * 0;
   if (InitKKs > kLen) {
      InitKKs = kLen -1;
   } 
   if (length(SigmaMultVec) < TotalRuns) {
     SigmaMultVec = -1;   
   }
     LCS =.C("LarsConvergencyPIEstShell", YYY = as.double(YYY), NLen = as.integer(TNLen),
              XXX= as.double(XXX), kLen = as.integer(kLen), 
              InitBetas = as.double(OnBetas[1:kLen]), OldBetas=as.double(StartBeta),
              USReturnBetas=as.double(vector("numeric", kLen+3)), 
              USReturnBetasAll =as.double(vector("numeric", kLen * TotalRuns)),
              NumTotalRuns=as.integer(TotalRuns), NumEMConv=as.integer(NumEMConv),
              MultEMCons=as.double(MultEMCons), StlambdaD=as.double(StlambdaD),
              StlambdaA=as.double(StlambdaA), lambdaDMultC=as.double(lambdaDMultC),
              lambdaAMultC=as.double(lambdaAMultC),
              lambdaAmultStop=as.integer(lambdaAmultStop),
              BBHatsAll=as.double(vector("numeric", kLen * TotalRuns)),
              NusKeep = as.double(vector("numeric", kLen * (TotalRuns+1))),
              LambdaDK = as.double(LambdaDK),
              LambdaAK = as.double(LambdaAK), LAKSeq = as.integer(LAKSeq),
              OrderSeq = as.integer(OrderSeq),
              RecordPostPBj = as.double(vector("numeric", (TotalRuns+1) * kLen) ), 
              DConfidenceInt = as.double( DConfidenceInt ),             
              RecordConfidenceInts = as.double(vector("numeric", (TotalRuns+1) * kLen * 2) ),
              m1 = as.double(m1), m2 = as.double(m2), NumEConv = as.integer(NumEConv),
              SigmaEta = as.double(SigmaEta), SigmaBarSq = as.double(SigmaBarSq),
              PiRecVec = as.double(vector("numeric", TotalRuns)), 
              SigmaRecVec = as.double(vector("numeric", TotalRuns)),
              NumCDOConv = as.integer(NumCDOConv), CDOEpsilon = as.double(CDOEpsilon),
              InitKKs = as.integer(InitKKs),
              InverseGammaConstant = as.double(InverseGammaConstant),
              WLSWeights = as.double(WLSWeights), TDFNu = as.double(TDFNu),
              PrintFlag = as.integer(PrintFlag),
              SigmaMultVec = as.double(SigmaMultVec)
        );
        LCS$BBHatsAll=t(matrix(LCS$BBHatsAll[1:(kLen*TotalRuns)],kLen,TotalRuns));
        LCS$USReturnBetasAll=t(matrix(LCS$USReturnBetasAll[1:(kLen*TotalRuns)],kLen,TotalRuns));  
        LCS$NusKeep=t(matrix(LCS$NusKeep[1:(kLen*TotalRuns)],kLen,TotalRuns));
        LCS$RecordPostPBj=t(matrix(LCS$RecordPostPBj[1:(kLen*TotalRuns)], kLen, TotalRuns) );
        RecordConfidenceInts =  array( dim=c(kLen,2, TotalRuns), dimnames=c("Bj", "LjRj", "Runs"));
        for (ii in 1:TotalRuns) {
           RecordConfidenceInts[,,ii] = 
                t(matrix(LCS$RecordConfidenceInts[ ((ii-1) * TotalRuns) + 1:(kLen*2)],2,kLen));
        }  
        EMLarsKeep <- list( NLen = LCS$NLen, kLen = LCS$kLen, 
           USReturnBetas = LCS$USReturnBetas[1:kLen],  USReturnBetasAll = LCS$USReturnBetasAll,
           NusKeep = LCS$NusKeep,
           NumTotalRuns = LCS$NumTotalRuns, NumEMConv = LCS$NumEMConv, 
           MultEMCons = LCS$MultEMCons,
           StlambdaA = LCS$StlambdaA, StlambdaD = LCS$StlambdaD,
           lambdaDMultC = LCS$lambdaDMultC, lambdaAMultC = LCS$lambdaAMultC, 
           lambdaAmultStop = LCS$lambdaAmultStop,
           SigmaBarSq = LCS$SigmaBarSq,
           SigmaEta = LCS$SigmaEta, 
           sdYYs = sdYYs, sdXXs = sdXXs,
           NusKeep = LCS$NusKeep[1:LCS$NumTotalRuns],
           BBHatsAll = LCS$BBHatsAll, BBHatsFinal = LCS$BBHatsAll[TotalRuns,], 
           ReturnBetasAll = 2, ReturnBetas = 2, 
           PiRecVec = LCS$PiRecVec, SigmaRecVec = LCS$SigmaRecVec,
           m1 = m1, m2 = m2, SigmaEta = SigmaEta, SigmaBarSq = SigmaBarSq,
           pPpiuseSt = pPpiuseSt,
           RecordPostPBj = LCS$RecordPostPBj,
           RecordConfidenceInts = RecordConfidenceInts,
           FailureFlag = 1, LambdaAK = LCS$LambdaAK[1:LCS$NumTotalRuns],
           LambdaDK = LCS$LambdaDK[1:LCS$NumTotalRuns],
           OrderSeq = LCS$OrderSeq[1:LCS$NumTotalRuns],
           NumCDOConv = LCS$NumCDOConv, CDOEpsilon = LCS$CDOEpsilon,
           InverseGammaConstant = as.double(InverseGammaConstant), TDFNu = as.double(TDFNu) );
        if (length(WLSWeights) > 1) {
          EMLarsKeep$WLSWeights = LCS$WLSWeights[1:length(yys)];
        }           
        if ( LCS$USReturnBetas[1] == -999) {
            Ontt1 = round(LCS$USReturnBetas[2]);
            Ontt2 = round(LCS$USReturnBetas[3]);
        } else {
            Ontt1 = TotalRuns;
            Ontt2 = NumEMConv;
        }
        
        if (StandardFlag == 1) {
          EMLarsKeep$USReturnBetasAll = LCS$USReturnBetasAll;
          EMLarsKeep$UsReturnBetas = LCS$USReturnBetas;
          EMLarsKeep$ReturnBetasAll = LCS$USReturnBetasAll;
          EMLarsKeep$ReturnBetasAll[1:Ontt1, 1:kLen] = t( t(
                      LCS$USReturnBetasAll[1:Ontt1, 1:kLen]  * sdYYs / sdXXs));
          if (LCS$USReturnBetas[1] != -999) {
                EMLarsKeep$ReturnBetas = LCS$USReturnBetas[1:kLen] * sdYYs / sdXXs;
          } else {
                EMLarsKeep$ReturnBetas = LCS$USReturnBetas[1:kLen] * 0 - 999;
          }
            EMLarsKeep$sigmaNoiseSq = sigmaNoiseSqSt * (sdYYs)^2; 
        } else {
           EMLarsKeep$ReturnBetasAll = LCS$USReturnBetasAll;
	          if (LCS$USReturnBetas[1] != -999) {
	                EMLarsKeep$ReturnBetas = LCS$USReturnBetas[1:kLen];
	          } else {
	                EMLarsKeep$ReturnBetas = LCS$USReturnBetas[1:kLen] * 0 - 999;
	          }       
        }
     if (LCS$InitBetas[1] == -999) {
       EMLarsKeep$FailureFlag = 1;
     } else {
        EMLarsKeep$FailureFlag = 0;
     }
     EMLarsKeep$Type <- "EMPI2Lasso";
     return(EMLarsKeep);              
}                           
                                   

###########################################################
## XLOperation ()
##
##   This is the Abramovich method of posterior Medians
##   This uses the pseudoinverse to plan parameter distributions
## 
XLOperationOld <- function(XX, YY, probpi, sigmasqNoise, TauOther = 60) {
    try(library(corpcor, warn.conficts=FALSE, quietly=TRUE), silent=TRUE);
    XTX <- t(XX) %*% XX;
    psXTX <- pseudoinverse(XTX);
    diag <- cbind( 1:length(XX[1,]), 1:length(XX[1,]));
    ReturnBetas <- psXTX %*% t(XX) %*% YY;
    ErrorEst <- sigmasqNoise * psXTX[diag];
    ErrorEst[ErrorEst < 0] = 0;
    tausqA = TauOther;
    ## UnOddsprobNew <- .5 *
    ##    ( 1 / (sigmasqNoise / length(YY) ) -
    ##      1 / ( TauOther + sigmasqNoise/length(YY)) ) * ReturnBetas^2
    ##OddsprobNew <- UnOddsprobNew;
    ##OddsprobNew[UnOddsprobNew < 100 &
    ##            UnOddsprobNew > -100 ] <- exp( UnOddsprobNew[UnOddsprobNew < 100 &
    ##                         UnOddsprobNew > - 100] )  *
    ##  sqrt(  (sigmasqNoise / length(YY)) /
    ##         ( TauOther + sigmasqNoise/length(YY)) ) * probpi / ( 1 - probpi) ;
    ##OddsprobNew[ UnOddsprobNew >= 100 | UnOddsprobNew <= - 100] <- 1;
    ##ProbNew <- OddsprobNew / ( 1 + OddsprobNew);    
    ##ProbNew[ UnOddsprobNew >= 100] <- 1.0;
    ##ProbNew[ UnOddsprobNew <= -100] <- 0.0;   
   BestGuestA <- psXTX %*% t(XX) %*% YY;
   VarOnA = ( diag(XTX)/sigmasqNoise + 1 / TauOther )^(-1)
      ##nDF = NLen - length(BetaProposed[BetaPropsed != 0] );
   logOddsPB1OverB0 = log(probpi / ( 1- probpi))  +  
               .5 * - log( diag(XTX) / sigmasqNoise )  +
               .5 * - log( tausqA + sigmasqNoise/ diag(XTX) )+
              .5 *  (diag(XTX)/ sigmasqNoise -  ( tausqA + sigmasqNoise/diag(XTX))^(-1)) *
              BestGuestA^2; 
      PosteriorProbB1 = rep(0, length(BestGuestA));
      if (length(PosteriorProbB1[abs(logOddsPB1OverB0) <= 60])  > 0 ) {
      PosteriorProbB1[abs(logOddsPB1OverB0) <= 60] = 
             exp( logOddsPB1OverB0[abs(logOddsPB1OverB0) <= 60]) / 
             ( 1 + exp( logOddsPB1OverB0[abs(logOddsPB1OverB0) <= 60]));
      }
      if (length( PosteriorProbB1[logOddsPB1OverB0 < -60]) > 0) {
        PosteriorProbB1[logOddsPB1OverB0 < -60] = 
               exp( logOddsPB1OverB0[logOddsPB1OverB0 < -60]);   
      }    
      if (length(PosteriorProbB1[logOddsPB1OverB0 > 60] )  > 0) {
          PosteriorProbB1[logOddsPB1OverB0 > 60] = 1;       
      }           
      PosteriorProbB0 = 1 - PosteriorProbB1;   
      ProbNew = PosteriorProbB1;   
    IntEstsX <- (ProbNew) * pnorm( -ReturnBetas / sqrt(ErrorEst),  0, 1 );
    IntEstsX[ErrorEst == 0] = (ProbNew[ErrorEst == 0]) *( 1 - sign(ReturnBetas[ErrorEst ==0])  ) / 2;
    BBHatsFinal <- IntEstsX * 0 + 1;
    BBHatsFinal[ IntEstsX <= .5 & IntEstsX + (1-ProbNew) >= .5] = 0;
    XLTReturn <- list(ReturnBetas = ReturnBetas, BBHatsFinal = BBHatsFinal, Type="XL");
    return(XLTReturn);
    
}
  XLPostAnalysisOld <- function( BetaProposed, xxs=-1, yys=-1, NLen= -1, ppiuse, sigmaNoiseSq,
         tausqA, XTX = -1, XTY = -1, YtY = 0) { 
 if (is.matrix(xxs) == FALSE || length(xxs) == 1 || xxs == -1) {
       if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == -1) {
          print("XLPostAnalysis, Error no proper X given, XTX False"); return;  
       }
       if ( (is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) || length(XTY) == 1 || XTY == -1) {
          print("XLPostAnalysis, Error no proper X given, XTY False"); return;
       }
       if (length(NLen) > 1 || NLen == -1) {
           print("XLPostAnalysis, Error, no proper NLen given");
           return(-1);
       }
       XTXFlag = 1;
       XXX = as.double(XTX); YYY = as.double(XTY);
       TNLen = - NLen; kLen = length(XTX[1,]);
   } else if (( is.vector(yys) == FALSE && is.matrix(yys) == FALSE) || length(yys) == 1 || yys == -1) {
        print("XLPostAnalysis, Error no proper Y given"); 
        print("   yys is : ");
        print( yys );
        return(-1);
   } else {
       XTXFlag = 0; NLen = length(yys);
       XXX = as.double(xxs); YYY = as.double(yys);
       XTX = t(xxs) %*% xxs;
       XTY <- t(xxs) %*% yys;
       YtY <- sum(yys*yys);
       TNLen = length(yys); kLen = length(xxs[1,]);
       ##if (InitKKs < 0) {
       ##   print("Cannot do xx version if xx is non matrix\n");
       ##   return;
       ##}       
   }   
      AllSumResidSq =  YtY - 2 * t(BetaProposed) %*% XTY + t(BetaProposed) %*% XTX %*% BetaProposed;          
      TurnOffSumResidSq = AllSumResidSq + 2 * as.vector(BetaProposed) * as.vector(XTY) - BetaProposed^2 * diag(XTX);
      XTYResidpD = as.vector(XTY) - as.vector(XTX %*% BetaProposed) + as.vector(diag(XTX) * BetaProposed);
      BestGuestA = ( diag(XTX)/sigmaNoiseSq + 1 / tausqA )^(-1) * XTYResidpD/ sigmaNoiseSq;
      VarOnA = ( diag(XTX)/sigmaNoiseSq + 1 / tausqA )^(-1)
      ##nDF = NLen - length(BetaProposed[BetaPropsed != 0] );
   logOddsPB1OverB0 = log(ppiuse / ( 1- ppiuse))  +  
               .5 * - log( diag(XTX) / sigmaNoiseSq )  +
               .5 * - log( tausqA + sigmaNoiseSq/ diag(XTX) )+
              .5 *  (diag(XTX)/ sigmaNoiseSq -  ( tausqA + sigmaNoiseSq/diag(XTX))^(-1)) *
              BestGuestA^2; 
      PosteriorProbB1 = rep(0, length(BestGuestA));
      if (length(PosteriorProbB1[ abs(logOddsPB1OverB0) < 60 ] )   > 0) {
      PosteriorProbB1[ abs(logOddsPB1OverB0) < 60 ] <-
          exp( logOddsPB1OverB0[ abs(logOddsPB1OverB0) < 60 ]) / 
           ( 1 + exp( logOddsPB1OverB0[ abs(logOddsPB1OverB0) < 60 ]));
      }
      if (length(PosteriorProbB1[ logOddsPB1OverB0 < - 60 ]) > 0) {
        PosteriorProbB1[ logOddsPB1OverB0 < - 60 ] <- 
             exp( logOddsPB1OverB0[ logOddsPB1OverB0 < -60 ]) 
      }
      if (length(PosteriorProbB1[ logOddsPB1OverB0 > 60 ])       > 0) {
            PosteriorProbB1[ logOddsPB1OverB0 > 60 ] <- 1;    
      }
      PosteriorProbB0 = 1 - PosteriorProbB1;     
      ##ReturnBetas[PosteriorProbB1 < .5] = 0;
      PLessThanZero = PosteriorProbB1* pnorm(0, BestGuestA, sqrt( ( 1/tausqA + diag(XTX)/ sigmaNoiseSq )^(-1) ) );
      ReturnBetas = BestGuestA;
      ReturnBetas[ PLessThanZero > .5] = BestGuestA[ PLessThanZero > .5 ];
      ReturnBetas[ PLessThanZero + PosteriorProbB0 < .5] = 
                    BestGuestA[ PLessThanZero + PosteriorProbB0 < .5];
      ReturnBetas[ PLessThanZero < .5 & PLessThanZero + PosteriorProbB0 > .5 ] <- 0;
      ReturnBBs = ReturnBetas; ReturnBBs[abs(ReturnBetas) > 0.0] <- 1;
      return(list(ReturnBetas = ReturnBetas, BBHatsFinal = ReturnBBs, BetaProposed = BetaProposed,
          BestGuestA = BestGuestA, logOddsPB1OverB0 = logOddsPB1OverB0,
          PosteriorProbB1 = PosteriorProbB1));
  }

###########################################################
## XLOperation ()
##
##   This is the Abramovich method of posterior Medians
##   This uses the pseudoinverse to plan parameter distributions
## 
XLOperation <- function(XX, YY, probpi, sigmasqNoise = -1, TauOther = 60, Interval = -1, ConfidenceQuantiles = NULL) {
    sigmaNoiseSq <- sigmasqNoise;
    try(library(corpcor, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
    Time1 <- proc.time();
    XTX <- t(XX) %*% XX;
    if (!exists("sigmasqNoise") || is.null(sigmasqNoise) || (is.numeric(sigmasqNoise) && sigmasqNoise[1] <= 0)) {
      if (length(YY) <= NCOL(XX)) {
        print("**********************************************************");
        print("** XLOperation: No, cannot go without sigmasqNoise if length(YY) >= NCOL(XX)");
        print("** Thus  we really need you to supply sigmasqNoise in this case. ");
        print("** I'm defaulting to sigmasqNoise = 1 as a warning. ")
        sigmasqNoise = 1;
      } else {
        EESq <- (sum(YY^2) - sum(YY * as.vector(XX %*% solve(XTX) %*% t(XX) %*% YY)));
        sigmasqNoise = EESq / ( length(Y) - 1 - NCOL(XX));
      }
    }

    psXTX <- pseudoinverse(XTX);
    diag <- cbind( 1:length(XX[1,]), 1:length(XX[1,]));
    ReturnBetas <- psXTX %*% t(XX) %*% YY;
    ErrorEst <- sigmasqNoise * psXTX[diag];
    ##print("Test ErrorEst "); flush.console();
    if (length(ErrorEst[ErrorEst < 0]) > 0) {
        ErrorEst[ErrorEst < 0] = 0;
    }
    tausqA = TauOther;
    ## UnOddsprobNew <- .5 *
    ##    ( 1 / (sigmasqNoise / length(YY) ) -
    ##      1 / ( TauOther + sigmasqNoise/length(YY)) ) * ReturnBetas^2
    ##OddsprobNew <- UnOddsprobNew;
    ##OddsprobNew[UnOddsprobNew < 100 &
    ##            UnOddsprobNew > -100 ] <- exp( UnOddsprobNew[UnOddsprobNew < 100 &
    ##                         UnOddsprobNew > - 100] )  *
    ##  sqrt(  (sigmasqNoise / length(YY)) /
    ##         ( TauOther + sigmasqNoise/length(YY)) ) * probpi / ( 1 - probpi) ;
    ##OddsprobNew[ UnOddsprobNew >= 100 | UnOddsprobNew <= - 100] <- 1;
    ##ProbNew <- OddsprobNew / ( 1 + OddsprobNew);    
    ##ProbNew[ UnOddsprobNew >= 100] <- 1.0;
    ##ProbNew[ UnOddsprobNew <= -100] <- 0.0;   
   BestGuestA <- psXTX %*% t(XX) %*% YY;
   VarOnA = ( diag(XTX)/sigmasqNoise + 1 / TauOther )^(-1)
      ##nDF = NLen - length(BetaProposed[BetaPropsed != 0] );
   logOddsPB1OverB0 = log(probpi / ( 1- probpi))  +  
               .5 * - log( diag(XTX) / sigmasqNoise )  +
               .5 * - log( tausqA + sigmasqNoise/ diag(XTX) )+
              .5 *  (diag(XTX)/ sigmasqNoise -  ( tausqA + sigmasqNoise/diag(XTX))^(-1)) *
              BestGuestA^2; 
      PosteriorProbB1 = rep(0, length(BestGuestA));
      if (length(PosteriorProbB1[abs(logOddsPB1OverB0) <= 60])  > 0 ) {
        PosteriorProbB1[abs(logOddsPB1OverB0) <= 60] = 
               exp( logOddsPB1OverB0[abs(logOddsPB1OverB0) <= 60]) / 
               ( 1 + exp( logOddsPB1OverB0[abs(logOddsPB1OverB0) <= 60]));
      }
      if (length( PosteriorProbB1[logOddsPB1OverB0 < -60]) > 0) {
        PosteriorProbB1[logOddsPB1OverB0 < -60] = 
               exp( logOddsPB1OverB0[logOddsPB1OverB0 < -60]);   
      }    
      if (length(PosteriorProbB1[logOddsPB1OverB0 > 60] )  > 0) {
           PosteriorProbB1[logOddsPB1OverB0 > 60] = 1;       
      }           
      PosteriorProbB0 = 1 - PosteriorProbB1;   
    ProbNew = PosteriorProbB1;
    ##print("Test IntEstsX"); flush.console();   
    IntEstsX <- (ProbNew) * pnorm( -ReturnBetas / sqrt(ErrorEst),  0, 1 );
    if (length(ErrorEst[ErrorEst == 0]) > 0) {
      IntEstsX[ErrorEst == 0] = (ProbNew[ErrorEst == 0]) *
            ( 1 - sign(ReturnBetas[ErrorEst ==0])  ) / 2;

    }
    BBHatsFinal <- IntEstsX * 0 + 1;
    Time2 <- proc.time();
    PointTime <- Time2[1:3] - Time1[1:3];
    ##print("Test BBHatsFinal");
    if (length(BBHatsFinal[ IntEstsX <= .5 & IntEstsX + (1-ProbNew) >= .5]) > 0) {
       BBHatsFinal[ IntEstsX <= .5 & IntEstsX + (1-ProbNew) >= .5] = 0;
    }
    if (exists("sigmasqNoise") && !exists("sigmaNoiseSq")) {
      sigmaNoiseSq  = sigmasqNoise;
    }
    CITime = NULL;
    if (exists("ConfidenceQuantiles") && !is.null(ConfidenceQuantiles) &&
      is.numeric(ConfidenceQuantiles)) {
      CITime1 <- proc.time();
      dSq <- diag(XTX)
      XTR <- t(XX) %*% (YY - XX %*% ReturnBetas) + dSq * ReturnBetas;
      ConfidenceMatrix <- matrix(0, NCOL(XX), length(ConfidenceQuantiles));
      UnshrunkConfidenceMatrix <- matrix(0, NCOL(XX), length(ConfidenceQuantiles));
      for (ii in 1:length(ConfidenceQuantiles)) {
         UnshrunkConfidenceMatrix[,ii] <-  (XTR
           ) / dSq + sqrt(sigmaNoiseSq/dSq) * qnorm(ConfidenceQuantiles[ii]);
      }
      colnames(UnshrunkConfidenceMatrix) <- ConfidenceQuantiles;
      rownames(UnshrunkConfidenceMatrix) <- 1:NROW(UnshrunkConfidenceMatrix);
      DL <- pnorm(-(XTR/dSq) / sqrt(sigmaNoiseSq/dSq));
      for (ii in 1:length(ConfidenceQuantiles)) {
        AT <- ProbNew*DL > ConfidenceQuantiles[ii]
        ConfidenceMatrix[AT,ii] <-
          XTR[AT]/dSq[AT] + sqrt(sigmaNoiseSq/dSq[AT]) * qnorm(ConfidenceQuantiles[ii] / ProbNew[AT])
        DT <- ProbNew *DL + (1-ProbNew) <= ConfidenceQuantiles[ii];
        ConfidenceMatrix[!AT & !DT,ii] <- 0.0;
        ConfidenceMatrix[DT, ii] <-  XTR[DT]/dSq[DT] - sqrt(sigmaNoiseSq/dSq[DT]) * qnorm((1-ConfidenceQuantiles[ii])/
          ProbNew[DT])
      }
      colnames(ConfidenceMatrix) <- ConfidenceQuantiles;
      rownames(ConfidenceMatrix) <- 1:NROW(UnshrunkConfidenceMatrix);
      CIEst <- list();
      CIEst[[1]] <- ConfidenceMatrix;
      CIEst[[2]] <- UnshrunkConfidenceMatrix;
      try(names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix"));
      CITime2 <- proc.time();
      CITime <- CITime2[1:3] - CITime1[1:3];
    } else {
      CITime <- NULL;  CIEst <- NULL; 
    }

    XLTReturn <- list(ReturnBetas = ReturnBetas, BBHatsFinal = BBHatsFinal, Type="XL",
      PointTime = PointTime, ProbNew = ProbNew, PosteriorProbB1 = PosteriorProbB1,
      CITime = CITime, CIEst = CIEst, ConfidenceQuantiles=ConfidenceQuantiles);
    ##print("Gone to Test Interval");
      if (!exists("Interval") || !is.numeric(Interval) || is.null(Interval) ||
           length(Interval) < 0 || Interval < 0) {
      } else  if (Interval < 1 && Interval > 0) {
         XLTReturn$CI = XLCreditability( BetaProposed= ReturnBetas, AllSumResidSq, 
         TurnOffSumResidSq, BestGuestA, VarOnA, logOddsPB1OverB0, PosteriorProbB1,
         Interval = Interval);
         XLTReturn$Interval = Interval;    
      }
    return(XLTReturn);
    
}

XLPostAnalysis <- function( BetaProposed, xxs=-1, yys=-1, NLen= -1, 
  ppiuse, sigmaNoiseSq, tausqA, XTX = -1, XTY = -1, YtY = 0, Interval = -1,
  ConfidenceQuantiles = NULL) { 
  Time1 <- proc.time();
  if (is.matrix(xxs) == FALSE || length(xxs) == 1 || xxs == -1) {
    if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == -1) {
      print("XLPostAnalysis, Error no proper X given, XTX False"); return;  
    }
    if ((is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) ||
      length(XTY) == 1 || XTY == -1) {
      print("XLPostAnalysis, Error no proper X given, XTY False"); return;
    }
    if (length(NLen) > 1 || NLen == -1) {
      print("XLPostAnalysis, Error, no proper NLen given");
      return(-1);
    }
    XTXFlag = 1;
    XXX = as.double(XTX); YYY = as.double(XTY);
    TNLen = - NLen; kLen = length(XTX[1,]);
  } else if (( is.vector(yys) == FALSE && is.matrix(yys) == FALSE) 
    || length(yys) == 1 || yys == -1) {
    print("XLPostAnalysis, Error no proper Y given"); 
    print("   yys is : ");
    print( yys );
    return(-1);
  } else {
    XTXFlag = 0; NLen = length(yys);
    Y = as.double(yys);
    if (!exists("XTX") || length(XTX) == 1) {
      XTX <- t(xxs) %*% xxs;
    }
    if (!exists("XTY") || length(XTY) == 1) {
      XTY <- t(xxs) %*% yys
    }
    if (length(dim(xxs)) != 2) {
      print("XLPostAnalysis: Error, defective xxs given no columns!");
      flush.console();
      return(-1);
    }
    if (NROW(xxs) != length(Y)) {
      print(paste("XLPostAnalysis: Error, NROW(xxs) = ", NROW(xxs),
        " but length(Y) = ", length(Y), sep="")); flush.console();
      return(-1);
    }
    if (NCOL(xxs) != length(BetaProposed)) {
      print(paste("XLPostAnalysis: ERROR, xxs not right dimension to multiply by BetaProposed."));
      print(paste("NCOL(xxs) = ", NCOL(xxs), " but length(BetaProposed) = ",
        length(BetaProposed), sep="")); flush.console();
      print("Go inspect error!"); return(-1);
    }
    YSResid = sum( (Y - xxs %*% BetaProposed)^2 );
    diagXTX = colSums(xxs^2);
    ATryText <- "Fun = -1;
      XTYResidPD = as.vector(t(xxs) %*% (as.vector(Y) - as.vector(xxs %*% BetaProposed))) + as.vector(diagXTX * BetaProposed) ;
      Fun = 1;
    "
    try(eval(parse(text=ATryText)));
    if (Fun == -1) {
      print("XLPostAnalysis Error, we still failed to generate XTYResidPD!");
      flush.console();
      print(paste("length(BetaPropsed) = ", length(BetaProposed), sep=""));
      flush.console();
      return(-1);
    }
    YSResidj =  YSResid  + 2 * BetaProposed * XTYResidPD - BetaProposed^2 * diagXTX;
    BSqOverA = (XTYResidPD/ (2*sigmaNoiseSq) )^2 / ( diagXTX/ (2*sigmaNoiseSq) + 1/(2*tausqA) );
    C = YSResidj /  (2*sigmaNoiseSq)
    if (min(YSResidj) < 0) {
      print("Error: min YSResidj < 0!!"); flush.console();
    }  
    Qj = diagXTX / sigma^2 + 1/tausqA;
    piDj = log(1-ppiuse) - .5 * length(Y) * log( 2*sigmaNoiseSq )   -C^2 / 2;
    piAj = log(ppiuse) - .5 * length(Y) * log( 2* sigmaNoiseSq) - C^2 /2 + 
      BSqOverA  +  .5 * log( 2 * pi / Qj);
    BNj = piAj / ( piAj + piDj);
    
  }   
  AllSumResidSq =  sum(yys^2) - 2 * t( BetaProposed) %*% XTY + 
    t(BetaProposed) %*% XTX %*% BetaProposed;          
  TurnOffSumResidSq = AllSumResidSq + 2 * as.vector(BetaProposed) * 
    as.vector(XTY) - BetaProposed^2 * diag(XTX);
  XTYResidpD = XTY - XTX %*% BetaProposed + diag(XTX) * BetaProposed;
  BestGuestA = ( diag(XTX)/sigmaNoiseSq + 1 / tausqA )^(-1) * XTYResidpD/ sigmaNoiseSq;
  VarOnA = ( diag(XTX)/sigmaNoiseSq + 1 / tausqA )^(-1)
  ##nDF = NLen - length(BetaProposed[BetaPropsed != 0] );
  logOddsPB1OverB0 = log(ppiuse / ( 1- ppiuse))  +  
    .5 * - log( diag(XTX) / sigmaNoiseSq )  +
    .5 * - log( tausqA + sigmaNoiseSq/ diag(XTX) )+
    .5 *  (diag(XTX)/ sigmaNoiseSq -  ( tausqA + sigmaNoiseSq/diag(XTX))^(-1)) *
    BestGuestA^2; 
    
  PosteriorProbB1 = rep(0, length(BestGuestA));
  if (length(PosteriorProbB1[ abs(logOddsPB1OverB0) < 60 ] )   > 0) {
    PosteriorProbB1[ abs(logOddsPB1OverB0) < 60 ] <-
      exp( logOddsPB1OverB0[ abs(logOddsPB1OverB0) < 60 ]) / 
        ( 1 + exp( logOddsPB1OverB0[ abs(logOddsPB1OverB0) < 60 ]));
  }
  if (length(PosteriorProbB1[ logOddsPB1OverB0 < - 60 ]) > 0) {
    PosteriorProbB1[ logOddsPB1OverB0 < - 60 ] <- 
      exp( logOddsPB1OverB0[ logOddsPB1OverB0 < -60 ]) 
  }
      if (length(PosteriorProbB1[ logOddsPB1OverB0 > 60 ])       > 0) {
            PosteriorProbB1[ logOddsPB1OverB0 > 60 ] <- 1;    
      }
      PosteriorProbB0 = 1 - PosteriorProbB1;     
      ##ReturnBetas[PosteriorProbB1 < .5] = 0;
      PLessThanZero = PosteriorProbB1* pnorm(0, BestGuestA, sqrt( ( 1/tausqA + diag(XTX)/ sigmaNoiseSq )^(-1) ) );
      ReturnBetas = BestGuestA;
      ReturnBetas[ PLessThanZero > .5] = BestGuestA[ PLessThanZero > .5 ];
      ReturnBetas[ PLessThanZero + PosteriorProbB0 < .5] = 
                    BestGuestA[ PLessThanZero + PosteriorProbB0 < .5];
      ReturnBetas[ PLessThanZero < .5 & PLessThanZero + PosteriorProbB0 > .5 ] <- 0;
      ReturnBBs = ReturnBetas; ReturnBBs[abs(ReturnBetas) > 0.0] <- 1;
      Time2 <- proc.time();
      PointTime <- Time2[1:3] - Time1[1:3];

    if (exists("ConfidenceQuantiles") && !is.null(ConfidenceQuantiles) &&
      is.numeric(ConfidenceQuantiles) && length(ConfidenceQuantiles) >= 1) {
      CITime1 <- proc.time();
      dSq <- diag(XTX)
      XTR <- t(xxs) %*% (yys - xxs %*% ReturnBetas) + dSq * ReturnBetas;
      ConfidenceMatrix <- matrix(0, NCOL(xxs), length(ConfidenceQuantiles));
      UnshrunkConfidenceMatrix <- matrix(0, NCOL(xxs), length(ConfidenceQuantiles));
      for (ii in 1:length(ConfidenceQuantiles)) {
         UnshrunkConfidenceMatrix[,ii] <-  (XTR
           ) / dSq + sqrt(sigmaNoiseSq/dSq) * qnorm(ConfidenceQuantiles[ii]);
      }
      colnames(UnshrunkConfidenceMatrix) <- ConfidenceQuantiles;
      rownames(UnshrunkConfidenceMatrix) <- 1:NROW(UnshrunkConfidenceMatrix);
      DL <- pnorm(-(XTR/dSq) / sqrt(sigmaNoiseSq/dSq));
      for (ii in 1:length(ConfidenceQuantiles)) {
        AT <- PosteriorProbB1*DL > ConfidenceQuantiles[ii]
        ConfidenceMatrix[AT,ii] <-
          XTR[AT]/dSq[AT] + sqrt(sigmaNoiseSq/dSq[AT]) * qnorm(ConfidenceQuantiles[ii] / PosteriorProbB1[AT])
        DT <- PosteriorProbB1 *DL + (1-PosteriorProbB1) <= ConfidenceQuantiles[ii];
        ConfidenceMatrix[!AT & !DT,ii] <- 0.0;
        if (sum(DT) >= 1) {
        ConfidenceMatrix[DT, ii] <-  XTR[DT]/dSq[DT] - sqrt(sigmaNoiseSq/dSq[DT]) * qnorm((1-ConfidenceQuantiles[ii])/
          PosteriorProbB1[DT])
        }
      }
      colnames(ConfidenceMatrix) <- ConfidenceQuantiles;
      rownames(ConfidenceMatrix) <- 1:NROW(UnshrunkConfidenceMatrix);
      CIEst <- list();
      CIEst[[1]] <- ConfidenceMatrix;
      CIEst[[2]] <- UnshrunkConfidenceMatrix;
      try(names(CIEst) <- c("CredibilityMatrix", "UnshrunkCredibilityMatrix"));
      CITime2 <- proc.time();
      CITime <- CITime2[1:3] - CITime1[1:3];
    } else {
      CITime <- NULL;  CIEst <- NULL; 
    }
    
      returner = list(ReturnBetas = ReturnBetas, BBHatsFinal = ReturnBBs, BetaProposed = BetaProposed,
          BestGuestA = BestGuestA, logOddsPB1OverB0 = logOddsPB1OverB0,
          PosteriorProbB1 = PosteriorProbB1, ConfidenceQuantiles=ConfidenceQuantiles,
          CIEst=CIEst, CITime=CITime, PointTime=PointTime);
      ##if (Interval < 1 && Interval > 0) {
      ##   returner$CI = XLCreditability( BetaProposed = BetaProposed, AllSumResidSq, 
      ##   TurnOffSumResidSq, BestGuestA, 
      ##   VarOnA, logOddsPB1OverB0, PosteriorProbB1,  Interval = Interval);
      ##   returner$Interval = Interval;    
      ##}
      return(returner);
  }
  Creditability <- function( muj, sigj, Bj, alphO2) {
    LB = rep(0, length(muj)); UB = rep(0, length(muj));
    qNL = qnorm( alphO2 / Bj, muj, sigj);
    LB[qNL < 0] = qNL;
    qNL = qnorm( 1.0 - alphO2 / (1-Bj), muj, sigj);
    UB[qNL > 0] = qNL;
    
    qTL = pnorm(0, muj, sigj) * Bj + (1-Bj);    

  
  }
  XLCreditability <- function( BetaProposed, AllSumResidSq, 
         TurnOffSumResidSq, BestGuestA, VarOnA, 
         logOddsPB1OverB0, PosteriorProbB1, Interval = .95) {
     ## AllSumResidSq =  YtY - 2 * t(BetaProposed) %*% XTY + t(BetaProposed) %*% XTX %*% BetaProposed;          
     ## TurnOffSumResidSq = AllSumResidSq + 2 * as.vector(BetaProposed) * as.vector(XTY) - BetaProposed^2 * diag(XTX);
     ## XTYResidpD = XTY - XTX %*% BetaProposed + diag(XTX) * BetaProposed;
     ## BestGuestA = ( diag(XTX)/sigmaNoiseSq + 1 / tausqA )^(-1) * XTYResidpD/ sigmaNoiseSq;
     ##  VarOnA = ( diag(XTX)/sigmaNoiseSq + 1 / tausqA )^(-1)
      ##nDF = NLen - length(BetaProposed[BetaPropsed != 0] );
      
      PosteriorProbB0 = 1 - PosteriorProbB1;     
      ##ReturnBetas[PosteriorProbB1 < .5] = 0;
      PLessThanZero = PosteriorProbB1* pnorm(0, BestGuestA, sqrt( VarOnA )^(-1)  );
      PMoreThanZero = PosteriorProbB1 - PLessThanZero;
      LB = rep(0, length(BetaProposed));
      UB = rep(0, length(BetaProposed));
      ITALL <- 1:length(BetaProposed);
      Workit <- ITALL[ (PLessThanZero <= (1.0-Interval)/2 & PLessThanZero + PosteriorProbB0 > (1.0-Interval)/2) ];
      if (length(Workit) > 0) 
      {   LB[Workit] = 0;     
         Workit <- ITALL[ PLessThanZero + PosteriorProbB0 < (1.0-Interval)/2 ];
           #### PosteriorProbB0 + PosteriorProbB1* pnorm( X, BestGuestA, sqrt( VarOnA)  ) = 
           ####    (1.0-Interval)/2;
         LB[Workit] = qnorm(  (  (1.0-Interval)/ 2 
                - PosteriorProbB0[Workit]) / PosteriorProbB1[Workit],
                BestGuestA[Workit], sqrt(  VarOnA[Workit] ) ); 
      }
      Workit <- ITALL[ PLessThanZero >= (1.0-Interval) /2 ];
      if (length(Workit) > 0) {
        LB[Workit] = qnorm( (1.0-Interval)/(2*PosteriorProbB1[Workit]),
            BestGuestA[Workit], sqrt( VarOnA[Workit] ) );
           ## PosteriorProbB1 * pnorm( X, BestGuestA, sqrt( ( 1/tausqA + diag(XTX)/ sigmaNoiseSq )^(-1) ) ) = (1.0-Interval)/2;
      }
      Workit = ITALL[ (PMoreThanZero + PosteriorProbB0 < (1.0-Interval)/2) ] ;
        ## PosteriorProbB1 * pnorm(X,BestGuestA, sqrt( ( 1/tausqA + diag(XTX)/ sigmaNoiseSq )^(-1) ) ) =
        ##    (1- (1.0-Interval)/2);
      if (length(Workit) > 0) {
       UB[Workit] = qnorm(   ( 1- (1.0-Interval)/2) / PosteriorProbB1[Workit],
                    BestGuestA[Workit], sqrt( VarOnA[Workit] ) );
      }
      Workit = ITALL[PMoreThanZero < (1.0-Interval)/2 & PMoreThanZero + PosteriorProbB0 >= (1.0-Interval)/2 ];
         UB[Workit] = 0;
      Workit = ITALL[PMoreThanZero >= (1.0-Interval)/2];
          ## PosteriorProbB0 + PosteriorProbPosteriorProbB1
          ##    * pnorm( X, BestGuestA, sqrt( ( 1/tausqA + diag(XTX)/ sigmaNoiseSq )^(-1) ) ) = 1.0 -(1.0-Interval)/2;
      if (length(Workit) > 0) {
          UB[Workit] = qnorm( (1.0 -(1.0-Interval)/2 
                               - PosteriorProbB0[Workit]) / PosteriorProbB1[Workit] ,
          BestGuestA[Workit], sqrt( VarOnA[Workit] ) )
      }  
      return(cbind(LB,UB)) 
  }
##############################################################
##   SimMeStuff()
##    This function simulates a sparse linear regression dataset
##    The covariates will be correlated.
##    Noise is Sigma noise level.  BRealVecs is vector of first set of non-zero values
##    ppCov,sigmaFactCov are parameters for correlating covariates.
SimMeStuff <- function(NN = -999, NR = -999, 
     Noise = -999, BRealVecs = -999, ppCov = -999, 
     sigmaFactCov = -999) {
 if (NR == -999) { NR = 20; }
 if (NN == -999) { NN = round(NR * .8); }
 
 if (Noise == -999) {Noise <- .5; }

 if (length(BRealVecs) == 1 && BRealVecs == -999) { 
    BRealVecs <- c(4, 3, -2.5, 1.5);   
 }    
 
  BReal <- 1:NR * 0;
  if (length(BRealVecs) > NR) { BRealVecs = BRealVecs[1:NR] }
  BReal[1:length(BRealVecs)] <- BRealVecs;
  puse <- min(length(BRealVecs) / length(BReal), 1);
  if (ppCov == -999) {
    ppCov  <- .8;
  }
  if (sigmaFactCov == -999) {  sigmaFactCov = .4; }
     
if (ppCov < 1 && ppCov >= 0) { 
    WCOR <- matrix( rnorm(NR * NR,sigmaFactCov,sigmaFactCov) * ( 2 * rbinom(NR*NR,1,.5) - 1) * 
                       rbinom( NR * NR, 1,ppCov), NR, NR);
     for (ii in 1:(NR)) {
	     WCOR[ii, ii:NR] <- 0;
	     if (sum(WCOR[ii,]^2) > 1) {
	              ##print("Hit it ");
	              WCOR[ii,ii] <-  2 *( sum(WCOR[ii,(1:(ii-1))])^2 - 1);
	              WCOR[ii,] <- WCOR[ii,] / sqrt(sum(WCOR[ii,]^2));
	     } else {
	                WCOR[ii,ii] <- sqrt( 1 - sum(WCOR[ii,]^2) )
	     }
     }  
	LCovX <- WCOR;
	CovX <- WCOR %*% t(WCOR);
	XX <- LCovX  %*% matrix( rnorm(NN * length(CovX[1,]),0,1), length(CovX[1,]), NN);
    XX <- t(XX);
} else if (ppCov == -20) {
   WCor <- matrix( sigmaFactCov, NR, NR);
   XX <- matrix( rnorm(NN * NR, 0, 1), NN, NR );
   XX <- XX + rnorm(NN, 0, 1);
}  else {
    WCOR <- matrix( sigmaFactCov, NR, NR);
    WCOR[cbind(1:NR, 1:NR)] <- 1;
	LCovX <- t(chol(WCOR));
	CovX <- WCOR;  
    XX <- LCovX  %*% matrix( rnorm(NN * length(CovX[1,]),0,1), length(CovX[1,]), NN);
    XX <- t(XX);	
}
       SigmaSqNoise = Noise^2;

		YY <-    XX %*% BReal + Noise * rnorm(NN, 0,1 );
	
     MyRet = list(XX=XX, YY=YY, puse = puse, 
          RealVec = BReal, 
          SigmaSqNoise = Noise^2, 
          NoiseLevelSigmaNotSq = Noise,
          BReal = BReal);
     return(MyRet);
}

#############################################
############# SignZero()
#############   Gives the sign vector (-1,0,or 1) of real numbers)
#############  Just the Sign of the value;
SignZero <- function(aa, CZ) {
 raa <- aa * 0;
 raa[aa > abs(CZ)] <- 1;
 raa[aa < -abs(CZ)] <- -1;
 return(raa);
}

 
##################################################################
##  InspectAllFunction
##
##     For small k values, inspects all possible calculations for L0 minimizer
InspectAllFunction <- function(YYtake, XXtake, PenKeep = 1, NKeep = 4) {  
   if(length(XXtake[1,]) > 12) {
   return(-1);
   }       
   oddif <- function(num) {
     retns <- num;
     num[floor(num) %%2 == 0] <- 0;
     num[floor(num) %%2 == 1] <- 1;
     return(num);
   }
   KLen = length(XXtake[1,]);
   KLenTo2 <- 2^(KLen);
   Outers <- 1:KLenTo2;
   MassiveCombos <- matrix(0, 2^(KLen), KLen);
   for (ii in 1:KLen) {
     MassiveCombos[ ,ii] <- oddif( (Outers-1) / 2^(ii-1));
    }
   
    PenMatrix <- matrix(0, KLenTo2, 3);
    BetMatrix <- matrix(0, KLenTo2, KLen);
    for (kk in 1:KLenTo2) {
        xxUse <- (1:KLen)[MassiveCombos[kk,] == 1];
        XXtt <- XXtake[,xxUse];
        XXttXXtt <- t(XXtt) %*% XXtt;
        if (length(xxUse) == 0) {
            BetEsts <- 1:KLen * 0;
            SumS2R <- sum ( YYtake^2 );
            BetMatrix[kk,] <- BetEsts; 
        } else if (SignZero(det(XXttXXtt), .000001) == 0) {
            BetMatrix[kk,]<- - 99999;
        } else if (length(xxUse) > 0) { 
            BetEsts <- solve( t(XXtt)%*% XXtt) %*% t(XXtt) %*%YYtake;
            SumS2R <- sum( (YYtake - XXtt%*%BetEsts)^2 );
            BetMatrix[kk,xxUse] <- BetEsts;            
        } else {
            BetEsts <- 1:KLen * 0;
            SumS2R <- sum ( YYtake^2 );
            BetMatrix[kk,] <- BetEsts;            
        }
        Penalty <- length(xxUse);
        PenMatrix[kk,] <- c(SumS2R + Penalty, SumS2R, Penalty);

    }
    Pens <- PenMatrix[,2] + PenKeep * PenMatrix[,3]
    StIndex <- sort(Pens, index=TRUE)$ix;
    Rlist <- list( UseId = MassiveCombos[StIndex[1:NKeep],],
                   BetRet = BetMatrix[StIndex[1:NKeep],],
                   PenKeeps = c( Pens[StIndex[1:NKeep]], PenMatrix[StIndex[1:NKeep],2:3] ) );
    return(Rlist);
}


 


#############################################################
###  The GiveGiveMe functions 
###   generally assess differences between multiple vectors.
###
###
###################################################
##  GiveGiveMe 
##    First assesses beta members active in first seq and not in second
##    Then assess betamembers active in second seq and not in first
GiveGiveMe <- function(Arrs) {
     
    if ((is.vector(Arrs) && !is.matrix(Arrs)) ||length(dim(Arrs)) != 2) {
       print(paste("dim(Arrs) = ", dim(Arrs), sep=""));
       return(-1);
    }
    Arrs <- matrix(as.numeric(Arrs), length(Arrs[,1]), length(Arrs[1,]));
    ALTER = -20;
      ALTOO = -20;
    if (ALTOO <= 0) {
         ALTOO = length(Arrs[,1]);
    }
    print(paste( "ALTOO =", ALTOO, sep=""));
    RetAD <- matrix(0, 2, ALTOO * (ALTOO-1) / 2);
    ndn <- 1;
    ##if (length(dim(Arrs)) != 2) {
    ##   print(paste("dim(Arrs) = ", dim(Arrs), sep=""));
    ##}
    iiList <- 1:(ALTOO-1);
    for (ii in 1:(ALTOO-1)) {
      jjList <- (ii+1):ALTOO;
     for( jj in jjList) {
        AD1 <- Arrs[ii,];
        AD2 <- Arrs[jj,];
        ADD <- AD1 - AD2;
        ##print(paste("ii = ", ii, " and jj = ", jj, " and jjList is ", sep=""));
        ##print(jjList);
        if (length(AD1) < 1 || length(AD2) < 1) {
           RetAD[1,ndn]  <- -10;
           RetAD[2,ndn] <- -10;
        } else if (is.nan(AD1[1])  || is.na(AD1[1]) ) {
           RetAD[1,ndn] <- - 10;
           RetAD[2,ndn] <- -10;
        } else if (is.nan(AD2[1])  || is.na(AD2[1]) ) {
           RetAD[1,ndn] <- - 10;
           RetAD[2,ndn] <- -10;       
        } else if (AD1[1] < -1) {
           RetAD[1,ndn] <- - 10;
           RetAD[2,ndn] <- -10; 
        } else if (AD2[1] < -1) {
           RetAD[1,ndn] <- - 10;
           RetAD[2,ndn] <- -10;              
       ## } else if ((is.nan(AD1[1])) || is.nan(AD2[1]) || AD1[1] < -1 || AD2[1] < -1) { 
       ##    RetAD[1,ndn] <- - 10;
       ##    RetAD[2,ndn] <- -10;
        } else {
           RetAD[1,ndn] <- length(ADD[ADD==1]);
           RetAD[2,ndn] <- length(ADD[ADD==-1]);
        }
        ndn <- ndn +1;
     }
    }
    return(RetAD);
}

###################################################
##  GiveGiveMe 
##    
##    Calculates total sum squared Beta difference differences.
##
GiveGiveMeDiffD <- function(ArrsCoefs) {
    ALTER = -20;
    if (ALTER <= 0) {
        ALTER <<- length(ArrsCoefs[,1]);
    }
      ALTOO = -20;
    if (ALTOO <= 0) {
         ALTOO = length(ArrsCoefs[,1]);
    }    
    RetAD <- vector("numeric", ALTOO * (ALTOO-1) / 2);
    ndn <- 1;
    for (ii in 1:(ALTOO-1)) {
     for( jj in (ii+1):ALTOO) {
        AD1 <- ArrsCoefs[ii,];
        AD2 <- ArrsCoefs[jj,];
        ADD <- AD1 - AD2;
        if (length(AD1) < 1 || length(AD2) < 1) {
           RetAD[ndn]  <- -10;
           RetAD[ndn] <- -10;
        } else if (is.nan(AD1[1])  || is.na(AD1[1]) ) {
           RetAD[ndn] <- - 10;
           RetAD[ndn] <- -10;
        } else if (is.nan(AD2[1])  || is.na(AD2[1]) ) {
           RetAD[ndn] <- - 10;
           RetAD[ndn] <- -10;       
        } else if (AD1[1] < -1) {
           RetAD[ndn] <- - 10;
           RetAD[ndn] <- -10;         
        ##} else if ((is.nan(AD1[1])) || is.nan(AD2[1]) || AD1[1] < -999 || AD2[1] < -999) { 
        ##   RetAD[ndn] <- - 10;
        } else {
           RetAD[ndn] <- sum(ADD^2);
        }
        ndn <- ndn +1;
     }
    }
    return(RetAD);
}

###############################################
##  TwoDigits is a Printing Algorithm
##      Just prints out the 10s place of data
##      useful representation for printing.
TwoDigits <- function(Nums, NDigits) {
   sabbyN <- abs(Nums);
   logss <- floor(log(sabbyN,10));
   SSAS <- round(Nums / 10^(logss), NDigits);
   return( SSAS * 10^logss);
}


#################################################
##   LarsGiveMeClose
##
##   Reduces a vector from LarsOb set to compile active and inactive vector.
##   Looks for first active vector set with IntClose Betas already activated.
LarsGiveMeClose <- function(LarsOb, IntClose, MinMin = .00001) {
ii = 1;
for (ii in 1:length(LarsOb$beta[,1])) {
  VVT <- LarsOb$beta[ii,];
  VVT[abs(VVT) > MinMin] <- 1;
  VVT[abs(VVT) <= MinMin] <- 0;
  TOT <- sum(VVT);
  if (sum(TOT) >= IntClose) {
     return(ii);
  }
}
return(-1);
}

##   Looks for first active vector set with IntClose Betas already activated.
MyLarsGiveMeClose <- function(MyLarsOb, IntClose, MinMin = .00001) {
ii = 1;
  if (is.null(MyLarsOb$ReturnBetasAll) && !is.null(MyLarsOb$beta)) {
      LenLen <- length(MyLarsOb$beta[,1]);
      MyLarsOb$Lambdas = MyLarsOb$lambda;
  } else if (!is.null(MyLarsOb$ReturnBetasAll)) {
      LenLen <- length(MyLarsOb$ReturnBetasAll[,1]);
  } else {
     print("MyLarsOb is not complete\n");
     return(1);
  }
  if (length(MyLarsOb$Lambdas) <= 0) {
    print("MyLarsGiveMeClose:  No Lambdas what so ever");
    flush.console();  return(1);
  }  else if (length(MyLarsOb$Lambdas) < LenLen) {
         LenLen =   length(MyLarsOb$Lambdas);
  }
  if (LenLen <= 0) {
     print("MyLarsGiveMeClose:  No Lambdas what so ever");
    flush.console();  return(1);
  }
for (ii in 1:LenLen) {
  if (is.null(MyLarsOb$ReturnBetasAll) && !is.null(MyLarsOb$beta)) {
      VVT <- MyLarsOb$beta[ii,];
  } else if (!is.null(MyLarsOb$ReturnBetasAll)) {
      VVT <- MyLarsOb$ReturnBetasAll[ii,];
  } else {
     print("MyLarsOb is not complete\n");
     return(1);
  }
  VVT[abs(VVT) > MinMin] <- 1;
  VVT[abs(VVT) <= MinMin] <- 0;
  TOT <- sum(VVT);
  if (ii > length(MyLarsOb$Lambdas)) {
    print(paste("Come On, ii = ", ii, " but you idiot you selected length(Lambdas) = ", 
       length(MyLarsOb$Lambdas) )    );
    flush.console();
    if (ii -1 > 0) { return(ii-1); } else { return(1); }
  }
  if (is.na(MyLarsOb$Lambdas[ii])) {
    print(paste("MyLarsOb$Lambdas[", ii, "] is NA", sep=""));
    print(MyLarsOb$Lambdas);
  } else if (is.nan(MyLarsOb$Lambdas[ii])) {
    print(paste("MyLarsOb$Lambdas[", ii, "] is NAN", sep=""));
    print(MyLarsOb$Lambdas);
  }  else if (is.null(MyLarsOb$Lambdas[ii])) {
    print(paste("MyLarsOb$Lambdas[", ii, "] is NAN", sep=""));
    print(MyLarsOb$Lambdas);
  }  else if (!is.finite(MyLarsOb$Lambdas[ii])) {
    print(paste("MyLarsOb$Lambdas[", ii, "] is NAN", sep=""));
    print(MyLarsOb$Lambdas);
  }  else if (sum(TOT) >= IntClose && !is.null(MyLarsOb$Lambdas)) {
      iion = ii;
      if (MyLarsOb$Lambdas[ii] < 0) {
            for (ii in iion:1) {
               if (MyLarsOb$Lambdas[ii] > 0) {
                  return(ii);
               }
            }
      }
  } else if (sum(TOT) >= IntClose) {
      iion = ii;
      if (MyLarsOb$lambda[ii] < 0) {
            for (ii in iion:1) {
               if (MyLarsOb$lambda[ii] > 0) {
                  return(ii);
               }
            }
      }  
  }
  if (sum(TOT) >= IntClose) {
     return(ii);
  }
}
return(length(MyLarsOb$ReturnBetasAll[,1]));
return(-1);
}

##############################################################
##############################################################

################################################################################3
##  Defines current SavedOutPutDirectory in the computer file.
##
LoadSavedOutPutDirectory <- function() {
   FilesInDIR <- unlist(list.files(.Library));
   if (!exists("PathMe")) {
     PathMe <- MakePathMe();
     eval(parse(text=SetGText("PathMe", S=1)));
   }
     if (any(FilesInDIR == "TwoLasso")) {
         eval(parse(text=SetGText(PathMe, S=1)));  
         FilesInDIR <- unlist(list.files(PathMe));
         if (any(FilesInDIR == "SavedOutput")) {
           PathMeFiles <<- paste(PathMe, "//SavedOutput", sep="");
         } else {
           PathMeFiles <<- paste(PathMe, "//SavedOutput", sep="");
           dir.create(PathMe, showWarnings = FALSE, recursive = FALSE);
         }
     } else {
         PathMeFiles <<- "c://Stat//2008Summer//LarsProject//code//SavedOutput//"
     }
     PathMeFiles <<- "c://Stat//2008Summer//LarsProject//code//SavedOutput//"
     SavedOutPutDirectory <<- PathMeFiles;
     return(SavedOutPutDirectory);
}
##SavedOutPutDirectory <- "c://Stat//2008Summer//LarsProject//Code//SavedOutputSigExp//"


####################################################################
##  tSeq  Printing out rounded numbers that doesn't blow out titles
##   These are involved in the names inevitably used for saved files.
####################################################################
##  tSeq  Printing out rounded numbers that doesn't blow out titles
##   These are involved in the names inevitably used for saved files.
tSeq <- function(Number) {

CN <- sub("e", "e", as.character(Number));
CN <- sub("-", "n", as.character(CN));
CN <- sub("\\+", "f", as.character(CN));
CN <- sub("\\.", "p", as.character(CN));
CN <- sub("0p", "p", as.character(CN));
   return(CN);
}


#########################################################################################
##  NamingClenture()
##
##   Naming standard for files.
##
NamingClenture <- function(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps) {
  NC <- paste(InitiateS, "A", tSeq(ppiuse), "B", tSeq(sigmaFactCov), "C", tSeq(ppCov), 
        "D", tSeq(NNOn), 
        "E", tSeq(NROn), "F", tSeq(NoiseOn), "G", tSeq(NumReps), ".csv", sep="");
  return(NC);
}

#########################################################################################
## NamingClenture2
##
##   Naming standard for files 2.
##
##
NamingClenture2 <- function(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps,
    BRealInput = c(1,-1,1,-1,1,-1), LARSSeekFlag = 0) {
  if (!is.vector(BRealInput)) {
      UseStart <- "NonV";
  } else if (length(BRealInput) == 6 && all(BRealInput == c(1,-1,1,-1,1,-1))) {
      UseStart <- InitiateS2;
  } else if (length(BRealInput) == 4 && all(BRealInput == c(4,3,-2.5,1))) {
      UseStart <- InitiateS;
  } else {
      UseStart <- "Alt";
  }
  if (LARSSeekFlag == 0) {
	  NC <- paste(UseStart, "A", tSeq(ppiuse), "B", tSeq(sigmaFactCov), "C", tSeq(ppCov), 
	        "D", tSeq(NNOn), 
	        "E", tSeq(NROn), "F", tSeq(NoiseOn), "G", tSeq(NumReps), ".csv", sep="");
  } else {
	  NC <- paste(UseStart, "LS", tSeq(LARSSeekFlag), 
	        "A", tSeq(ppiuse), "B", tSeq(sigmaFactCov), "C", tSeq(ppCov), 
	        "D", tSeq(NNOn), 
	        "E", tSeq(NROn), "F", tSeq(NoiseOn), "G", tSeq(NumReps), ".csv", sep="");
  
  }
  return(NC);
}


################################################################################
##  GetForMeAll
##
##    Under the standardied "SimMeStuff1" version, runs the simulations and saves in files
##     for simulation counts of all the combinations of ppiuse, sigmaFact, ppCov, NNOn, NROn, etc.
##
GetForMeAll <- function(ppiuseV, sigmaFactV, ppCovV, NNOnV, NROnV, NoiseOnV, NumReps,
     PrintFlag = 0, musize=10) {
  StuffMeMatrix <- matrix(0, NumReps * length(ppiuseV) * length(sigmaFactV) * 
                                length(ppCovV) * length(NNOnV) * 
                                length(NROnV) * length(NoiseOnV), 7 + ALT * 3  +1);
  SavedOutPutDirectory <- LoadSavedOutPutDirectory();
  FilesInDIR <- unlist(list.files(SavedOutPutDirectory));
  countt <- 0;
  for (ii1 in 1:length(ppiuseV)) {
     ppiuse <- ppiuseV[ii1];
     for (ii2 in 1:length(sigmaFactV)) {
         sigmaFactCov <- sigmaFactV[ii2];
         for (ii3 in 1:length(ppCovV)) {
            ppCov <- ppCovV[ii3];
            for (ii4 in 1:length(NNOnV)) {
               NNOn <- NNOnV[ii4];
               for (ii5 in 1:length(NROnV)) {
                 NROn = NROnV[ii5];
                 for (ii6 in 1:length(NoiseOnV)) {
                    NoiseOn = NoiseOnV[ii6];
                    GMF <- GetForMeFunction(ppiuse, sigmaFactCov, ppCov, NNOn, 
                       NROn, NoiseOn, NumReps, 
                       FilesInDIR, musize);
                    StuffMeMatrix[ (countt+1):(countt + NumReps),  ] <- as.matrix(GMF);
                    countt <- countt + NumReps;
                    if (PrintFlag > 0 && countt %% PrintFlag == 0) {
                       plot(GMF[, 8]~GMF[,9], main=paste(
                          "ppiuse = ", ppiuse, ", sigFC = ", sigmaFactCov, 
                          ", ppC = ", ppCov, "\n NNOn = ", NNOn,
                          ", NROn = ", NROn, ", NoiseOn = ", NoiseOn, 
                          "countt = ", countt, sep="") );
                    }
                 }
               }
             }
          }
     }
  }
  return(StuffMeMatrix);
}

################################################################################
##  GetForMeAll2
##
##    Under the standardied "SimMeStuff2" version, runs the simulations and saves in files
##     for simulation counts of all the combinations of ppiuse, sigmaFact, ppCov, NNOn, NROn, etc.
##
GetForMeAll2 <- function(ppiuseV, sigmaFactV, ppCovV, NNOnV, NROnV, NoiseOnV, NumReps,
     PrintFlag = 0, BRealInput = c(1,-1,1,-1,1,-1), LARSSeekFlag = 0, musize=10) {
                    ALT2 <<- 10 * 11 / 2;  ALT3 <<- 10 ;
  StuffMeMatrix <- matrix(0, NumReps * length(ppiuseV) * length(sigmaFactV) * 
                                length(ppCovV) * length(NNOnV) * 
                                length(NROnV) * length(NoiseOnV), 7 + ALT2 * 3 + ALT3 * 3 );
  SavedOutPutDirectory <- LoadSavedOutPutDirectory();                                
  FilesInDIR <- unlist(list.files(SavedOutPutDirectory));
  countt <- 0;
  for (ii1 in 1:length(ppiuseV)) {
     ppiuse <- ppiuseV[ii1];
     for (ii2 in 1:length(sigmaFactV)) {
         sigmaFactCov <- sigmaFactV[ii2];
         for (ii3 in 1:length(ppCovV)) {
            ppCov <- ppCovV[ii3];
            for (ii4 in 1:length(NNOnV)) {
               NNOn <- NNOnV[ii4];
               for (ii5 in 1:length(NROnV)) {
                 NROn = NROnV[ii5];
                 for (ii6 in 1:length(NoiseOnV)) {
                    NoiseOn = NoiseOnV[ii6];
                    GMF <- GetForMeFunction2(ppiuse = ppiuse, 
                       sigmaFactCov = sigmaFactCov, ppCov = ppCov, NNOn = NNOn, 
                       NROn = NROn, NoiseOn = NoiseOn, NumReps = NumReps, 
                       FilesInDIR = FilesInDIR, BRealInput = BRealInput,
                       LARSSeekFlag = LARSSeekFlag, musize);
                    StuffMeMatrix[ (countt+1):(countt + NumReps),  ] <- as.matrix(GMF);
                    countt <- countt + NumReps;
                    if (PrintFlag > 0 && countt %% PrintFlag == 0) {
                       plot(GMF[, 8]~GMF[,9], main=paste(
                          "ppiuse = ", ppiuse, ", sigFC = ", sigmaFactCov, 
                          ", ppC = ", ppCov, "\n NNOn = ", NNOn,
                          ", NROn = ", NROn, ", NoiseOn = ", NoiseOn, 
                          "countt = ", countt, "\n",
                          "LARSSeekFlag = ", LARSSeekFlag, sep="") );
                    }
                 }
               }
             }
          }
     }
  }
  return(StuffMeMatrix);
}


##################################################################
##  GetForMeAAA  
##    Given a stack of vectors OneVV which should define ppiuse -> SigmaNoise values
##   this then simulates enough samples using SimForMe (the 4 vector version).
##
GetForMeAAA <- function(  OneVV, NumReps, PrintFlag = 0, musize = 10) {
 ##ppiuseV, sigmaFactV, ppCovV, NNOnV, NROnV, NoiseOnV, NumReps,
 
   ##if(is.null(ALT2)) {
    DeclareGlobals();
  ##}
  
  StuffMeMatrix <- matrix(0, length(OneVV[,1]) * NumReps, 7 + ALT * 3 +1 );
  SavedOutPutDirectory <- LoadSavedOutPutDirectory();
  FilesInDIR <- unlist(list.files(SavedOutPutDirectory));
  countt <- 0;
  for (ii in 1:length(OneVV[,1])) {
     ppiuse <- OneVV[ii,1]; sigmaFactCov <- OneVV[ii,2];
     ppCov <- OneVV[ii,3]; NNOn <- OneVV[ii,4]; NROn <- OneVV[ii,5];
     NoiseOn <- OneVV[ii,6];
     GMF <- GetForMeFunction(ppiuse, sigmaFactCov, ppCov, NNOn, 
                       NROn, NoiseOn, NumReps, 
                       FilesInDIR );
     StuffMeMatrix[ (countt+1):(countt + NumReps),  ] <- as.matrix(GMF);
     countt <- countt + NumReps;
      if (PrintFlag > 0 && countt %% PrintFlag == 0) {
                       plot(GMF[, 8]~GMF[,9], main=paste(
                          "ppiuse = ", ppiuse, ", sigFC = ", sigmaFactCov, 
                          ", ppC = ", ppCov, "\n NNOn = ", NNOn,
                          ", NROn = ", NROn, ", NoiseOn = ", NoiseOn, 
                          "countt = ", countt, sep="") );
      }
      
  }
  return(StuffMeMatrix);
}

##################################################################
##  GetForMeAAA2 
##    Given a stack of vectors OneVV which should define ppiuse -> SigmaNoise values
##   this then simulates enough samples using SimForMe2 (the 6 vector version).
##
GetForMeAAA2 <- function(  OneVV, NumReps, PrintFlag = 0, BRealInput = c(1,-1,1,-1,1,-1),
  LARSSeekFlag = 0, musize =10) {
  ##if(is.null(ALT2)) {
    DeclareGlobals();
  ##}
  ##  if(is.null(ALT2)) {
    DeclareGlobals();
  ## }
 ##ppiuseV, sigmaFactV, ppCovV, NNOnV, NROnV, NoiseOnV, NumReps,
             ## ALT2 <<- 10 * 11 / 2; 
             ## ALT3 <<- 10 ;
              if (PrintFlag < 0) {
                 PrintOutFlags <- 0;
              } else if (PrintFlag == 0) {
                 PrintOutFlags <- 1
              } else {
                 PrintOutFlags <- PrintFlag;
              }
  StuffMeMatrix <- matrix(0, length(OneVV[,1]) * NumReps, 7 + ALT2 * 3 + ALT3 * 3);
  SavedOutPutDirectory <- LoadSavedOutPutDirectory();
  FilesInDIR <- unlist(list.files(SavedOutPutDirectory));
  countt <- 0;
  for (ii in 1:length(OneVV[,1])) {
     ppiuse <- OneVV[ii,1]; sigmaFactCov <- OneVV[ii,2];
     ppCov <- OneVV[ii,3]; NNOn <- OneVV[ii,4]; NROn <- OneVV[ii,5];
     NoiseOn <- OneVV[ii,6];
     GMF <- GetForMeFunction2(ppiuse = ppiuse, sigmaFactCov = sigmaFactCov,
                       ppCov = ppCov, NNOn = NNOn, 
                       NROn = NROn, NoiseOn = NoiseOn, NumReps = NumReps, 
                       FilesInDIR = FilesInDIR, PrintOutFlags = PrintOutFlags,
                       BRealInput = BRealInput, LARSSeekFlag = LARSSeekFlag);
         if( dim(as.matrix(GMF))[2] != length(StuffMeMatrix[1,])) {
          GMF = as.matrix(GMF);
          print(paste("GetForMeAAA2: Stuff Error at OneVV[", ii, "], we have : ", sep=""));
          print(paste("Def = 10", "  ALT2 = ", ALT2, "  ALT3 = ", ALT3, sep=""));
          print(paste(" StuffMeMatrixDims = c( ", length(StuffMeMatrix[,1]), ", ", 
                     length(StuffMeMatrix[1,]), ")", sep=""));
          print(paste(" GMF = c( ", length(GMF[,1]), ", ", 
                     length(GMF[1,]), ")", sep=""));  
          print(paste("  NumReps = ", NumReps, sep=""));
          print( paste("  countt = ", countt, " and countt+NumReps=", countt+NumReps, sep=""));
          return(StuffMeMatrix);
     }                 
     StuffMeMatrix[ (countt+1):(countt + NumReps),  ] <- as.matrix(GMF);
     countt <- countt + NumReps;
      if (PrintFlag > 0 && countt %% PrintFlag == 0) {
                       plot(GMF[, 8]~GMF[,9], main=paste(
                          "ppiuse = ", ppiuse, ", sigFC = ", sigmaFactCov, 
                          ", ppC = ", ppCov, "\n NNOn = ", NNOn,
                          ", NROn = ", NROn, ", NoiseOn = ", NoiseOn, 
                          "countt = ", countt, "\n",
                          "LARSSeekFlag = ", LARSSeekFlag, sep="") );
      }
      
  }
  return(StuffMeMatrix);
}



##################################################################
##  GetForMeFunction
##  
##   Given a single active values of ppiuse through NoiseOn, this simulates
##   NumReps and saves them temporarily to a file name, and finally
##   (conditional on a full run), saves all to a final file.  If the program is 
##   stopped temporarily, the SimForMe function temporarily saves files to a "-a" file
##
GetForMeFunction <- function(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, 
                     FilesInDIR = -999, musize = 10) {
   SavedOutPutDirectory <- LoadSavedOutPutDirectory();
   NCN <- NamingClenture(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps);
   if (length(FilesInDIR) == 1 && FilesInDIR == -999) {
       FilesInDIR <-unlist(list.files(SavedOutPutDirectory));
   }
   if (any( NCN == FilesInDIR) ) {
      MyTable <- read.table( paste(SavedOutPutDirectory, NCN, sep=""), sep=",", header=TRUE);
      return(MyTable);
   } else if (any (NCN == paste(FilesInDIR, "a", sep=""))) {
      MyTableI <- as.matrix(read.table( paste(SavedOutPutDirectory, NCN, "a", sep=""), sep=",", 
                                header = TRUE));
      MTL <- length(MyTableI[MyTableI[,1] != 0,1]);
         FileNamea = paste(SavedOutPutDirectory, NCN, "a", sep="");
      MyTable <- SimForMe2(ppiuse = ppiuse, 
      sigmaFactCov = sigmaFactCov, ppCov=ppCov, NNOn = NNOn, NROn = NROn, NoiseOn=NoiseOn, 
         NumReps=NumReps, SigmaEta = 2, SigmaBar = -999, musize =4, 
               InputTable = MyTableI, FileName = FileNamea); 
      write.table(MyTable, paste(SavedOutPutDirectory, NCN, sep=""), sep=", ", eol="\n", na = "NA", quote=FALSE, 
                 row.names=FALSE, col.names=colnames(MyTable));
      return(MyTable);                                       
   } else {
      FileNamea = paste(SavedOutPutDirectory, NCN, "a", sep="");
      MyTable <- SimForMe2(ppiuse = ppiuse, 
      sigmaFactCov = sigmaFactCov, ppCov=ppCov, NNOn = NNOn, NROn = NROn, NoiseOn=NoiseOn, 
         NumReps=NumReps, SigmaEta = 2, SigmaBar = -999, musize =4, 
             InputTable=NULL, FileName = FileNamea);
      write.table(MyTable, paste(SavedOutPutDirectory, NCN, sep=""), sep=", ", eol="\n", na = "NA", quote=FALSE, 
                 row.names=FALSE, col.names=colnames(MyTable));
      return(MyTable);
   }
}
##GetForMeFunction <- function(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, 
##                     FilesInDIR = -999) {
##   NCN <- NamingClenture(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps);
##   if (length(FilesInDIR) == 1 && FilesInDIR == -999) {
##       FilesInDIR <-unlist(list.files(SavedOutPutDirectory));
##   }
##   if (any( NCN == FilesInDIR) ) {
##      MyTable <- read.table( paste(SavedOutPutDirectory, NCN, sep=""), sep=",", header=FALSE);
##      return(MyTable);
##   } else {
##      MyTable <- SimForMe(ppiuse = ppiuse, 
##      sigmaFactCov = sigmaFactCov, ppCov=ppCov, NNOn = NNOn, NROn = NROn, NoiseOn=NoiseOn, 
##         NumReps=NumReps);
##      write.table(MyTable, paste(SavedOutPutDirectory, NCN, sep=""), sep=", ", eol="\n", na = "NA", quote=FALSE, 
##                 row.names=FALSE, col.names=FALSE);
##      return(MyTable);
##   }
##}

##################################################################
##  GetForMeFunction2
##  
##   Given a single active values of ppiuse through NoiseOn, this simulates
##   NumReps and saves them temporarily to a file name, and finally
##   (conditional on a full run), saves all to a final file.  If the program is 
##   stopped temporarily, the SimForMe2 function temporarily saves files to a "-a" file
##
GetForMeFunction2 <- function(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, 
                     FilesInDIR = -999, PrintFlag = 0, PrintOutFlags = 1,
                     BRealInput = c(1,-1,1,-1,1,-1), LARSSeekFlag = 0, musize = 10) {
              if (PrintOutFlags >0) {
              } else if (PrintFlag < 0) {
                 PrintOutFlags <- 0;
              } else if (PrintFlag == 0) {
                 PrintOutFlags <- 1
              } else {
                 PrintOutFlags <- PrintFlag;
              }                     
   SavedOutPutDirectory <- LoadSavedOutPutDirectory();
   NCN <- NamingClenture2(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps,
       BRealInput = BRealInput, LARSSeekFlag = LARSSeekFlag);
   if (length(FilesInDIR) == 1 && FilesInDIR == -999) {
       FilesInDIR <-unlist(list.files(SavedOutPutDirectory));
   }
   if (any( NCN == FilesInDIR) ) {
      MyTable <- read.table( paste(SavedOutPutDirectory, NCN, sep=""), sep=",", header=TRUE);
      return(MyTable);
   } else if (any (paste(NCN, "a", sep="") == FilesInDIR)) {
      MyTableI <- as.matrix(read.table( paste(SavedOutPutDirectory, NCN, "a", sep=""), sep=",", 
                                header = TRUE));
      MTL <- length(MyTableI[MyTableI[,1] != 0,1]);
         FileNamea = paste(SavedOutPutDirectory, NCN, "a", sep="");
     print(paste("GetForMe2: InputTable for has nonzero length = ", MTL, " of ", NumReps, sep=""));
     print(paste("GetForMe2: FileNamea = ", FileNamea, sep=""));
     flush.console();       
      MyTable <- SimForMe2(ppiuse = ppiuse, 
      sigmaFactCov = sigmaFactCov, ppCov=ppCov, NNOn = NNOn, NROn = NROn, NoiseOn=NoiseOn, 
         NumReps=NumReps, SigmaEta = 2, SigmaBar = -999, musize =4, 
               InputTable = MyTableI, FileName = FileNamea, PrintOutFlags = PrintOutFlags, 
               LARSSeekFlag = LARSSeekFlag, BRealInput); 
      write.table(MyTable, paste(SavedOutPutDirectory, NCN, sep=""), sep=", ", eol="\n", na = "NA", quote=FALSE, 
                 row.names=FALSE, col.names=colnames(MyTable));
      return(MyTable);                                       
   } else {
      FileNamea = paste(SavedOutPutDirectory, NCN, "a", sep="");
      print(paste("GetForMe2: InputTable Not Found, will start new", sep=""));
      print(paste("GetForMe2: FileNamea = ", FileNamea, sep=""));       
      flush.console();  
      MyTable <- SimForMe2(ppiuse = ppiuse, 
      sigmaFactCov = sigmaFactCov, ppCov=ppCov, NNOn = NNOn, NROn = NROn, NoiseOn=NoiseOn, 
         NumReps=NumReps, SigmaEta = 2, SigmaBar = -999, musize = musize, 
             InputTable=NULL, FileName = FileNamea, PrintOutFlags = PrintOutFlags,
             LARSSeekFlag = LARSSeekFlag, BRealInput);
      write.table(MyTable, paste(SavedOutPutDirectory, NCN, sep=""), sep=", ", eol="\n", na = "NA", quote=FALSE, 
                 row.names=FALSE, col.names=colnames(MyTable));
      return(MyTable);
   }
}

##############################################################################
##  SimForMe2()
##
##   An integrated assessment of many/most model selection mechansims for many
##   possible choises of n, k, noise, and fixed active vector (1,-1,1,-1,1,-1,1) for first
##    seven covariates.
##
##  Simulates LARS Fixed, CP Lars solution, Ming and Yuan LARS solution, Limit-Ridge,
##   a Fast EM2Lasso, a Limit-Lasso, and Limit-Pi-Lasso (semi-estimates Pi), FD-Lasso,
##   Marginal Median, and Marginal Median Lasso 
##
SimForMe2 <- function(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, 
                            musize = 10, SigmaEta = 2, SigmaBar =-999, 
                            InputTable = NULL, FileName = NULL, PrintOutFlags = 0,
                            BRealInput = c(1,-1,1,-1,1,-1),
                    LARSSeekFlag = 0) {
       try(library(lars, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
       try(library(corpcor, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
        ## PrintOutFlags = 2;
  ##if(is.null(ALT2)) {
    DeclareGlobals();
  ##}
             ##ALT2 <- 10 * 11 / 2;  ALT3 <- 10 ;

      ReturnMatrix <- matrix(0, NumReps,  3 * ALT2 + 3 * ALT3 + 7);
      ii4Start <- 1;
      print(paste("SimForMe2: LARSSeekFlag = ", LARSSeekFlag, sep=""));
      if (is.null(InputTable)) {  
      } else if ( length(InputTable[1,]) != length(ReturnMatrix[1,])) {
         print(paste("SimForMe2: InputTable Not shaped like Return Matrix, dim(InputTable) = (",
                      dim(InputTable)[1], ", ", dim(InputTable)[2], ")", sep=""));
         print(paste("FileName: ", FileName, sep=""));      
      } else {
         LIT <- length(InputTable[InputTable[,1] != 0,1]);
         ReturnMatrix[1:LIT,] <- 
             InputTable[1:LIT,];
         ii4Start <- LIT +1;
         print(paste("SimForMe2: InputTable for has nonzero length = ", ii4Start-1, " of ", NumReps, sep=""));
         print(paste("FileName: ", FileName, sep=""));
      }
      if (length(BRealInput) > 1) { 
         BRealVecs = BRealInput;
      } else { 
         BRealVecs <- c(1,-1,1,-1,1,-1,1);
      }
      PPract <- length(BRealVecs) / NROn;
      CLLIST <- c("Real", "LS2", "LSCP", "LSMY", "LimRidge", "F2Lasso", "Lim2Lasso", 
                      "Lim2PILasso", "XLMedian", "LimFDLasso", "XLMLasso");
        AltList1 <- rep(0, ALT2); AltList2 <- rep(0, ALT2); AltList3 <- rep(0,ALT2); nii = 1;
        for (kii in 1:(length(CLLIST)-1)) {
            aS <- CLLIST[kii];
            for (jjj in (kii+1):length(CLLIST) ) {
               bS <- CLLIST[jjj];
               AltList1[nii] <- paste("In", aS, "Not", bS, sep="");
               AltList2[nii] <- paste("In", bS, "Not", aS, sep="");    
               AltList3[nii] <- paste("SumDiff", aS, "and", bS, sep="");   
               nii <- nii +1;                     
            }
        }
        colnames(ReturnMatrix) <- c( "ppiuse", "sigmaFactCov", "ppCov", "NNOn", "NROn", 
                  "NoiseOn", "NumReps", AltList1, AltList2, AltList3, 
                   paste("tUser", CLLIST[2:length(CLLIST)], sep=""),
                   paste("tSystem", CLLIST[2:length(CLLIST)], sep=""),
                   paste("tElapsed", CLLIST[2:length(CLLIST)], sep="")  );
     ## Just a TestSim        
      SMS <- SimMeStuff(NN = NNOn, NR = NROn, 
			     Noise = NoiseOn, BRealVecs = BRealVecs, ppCov = ppCov, 
			     sigmaFactCov = sigmaFactCov);
      if (SigmaEta > 0) {
         SigmaBar <- NoiseOn^2;
      }                   
     for (ii4 in ii4Start:NumReps) {
        if (PrintOutFlags >= 1) {
          print(paste("ii4 = ", ii4, ", starting beginning"));
          flush.console();
        }
	     SMS <- SimMeStuff(NN = NNOn, NR = NROn, 
			     Noise = NoiseOn, BRealVecs = BRealVecs, ppCov = ppCov, 
			     sigmaFactCov = sigmaFactCov);
	     SMSReal = as.vector(SMS$BReal);
	     SMSReal[abs(SMS$BReal) <= MinMin] = 0;
	     SMSReal[abs(SMS$BReal) > MinMin] = 1;
	     if (SigmaBar == -999) {
	       SigmaBarS <- (.5 * sd(as.vector(SMS$YY)))^2
	     } else {
	       SigmaBarS <- SigmaBar;
	     }
	     ##print(paste("ii4 = ", ii4, ", starting with Lars2Fit"));
	    T1 <- proc.time();
	    ##HitBatWing <- Batwing(RatWant=.4, SMS$XX, SMS$YY, 
	    ##                      BetaBarMin = BetaBarMin, puse = SMS$puse);
	    try(library(lars, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
	    HitLars2Ob <- lars(SMS$XX, SMS$YY);
	    ITO <- LarsGiveMeClose(HitLars2Ob, IntClose = round(SMS$puse * NROn) , MinMin = .00001);    
	    Lars2Fit <- HitLars2Ob$beta[ITO,];
	    Lars2Fit[abs(Lars2Fit) >= MinMin] <- 1;
	    Lars2Fit[abs(Lars2Fit) < MinMin] <- 0;  
	    Lars2Fit2 <- Lars2Fit;
        
	      if (sum(Lars2Fit) > 0 || length(SMS$XX[1,Lars2Fit == 1]) > 0) { 
	         XXSM <-  SMS$XX[,Lars2Fit == 1];
	         XXSXX <- t(XXSM) %*% XXSM;
	         if (any(is.null(XXSXX)) || length(dim(XXSXX)) != 2 || dim(XXSXX)[1] != dim(XXSXX)[2] || 
	                dim(XXSXX)[1] < 1 || any(is.nan(XXSXX))) { 
	             Lars2Fit2 <- Lars2Fit; BetaSM <- 0;
	         } else if (length(XXSXX) == 1) {
	              BetaSM <- 1 / XXSXX * sum(XXSM * SMS$YY);
	         } else if ( !is.null(dim(XXSXX)) && dim(XXSXX)[1] > 1 && dim(XXSXX)[1] == dim(XXSXX)[2] && det(XXSXX) == 0 ) {
	          try(library(corpcor, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
	          BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;
	        } else if (!is.null(dim(XXSXX)) ){
	          BetaSM <- solve(t(XXSM) %*% XXSM) %*% t(XXSM) %*% SMS$YY;
	        }
	          Lars2Fit2[Lars2Fit != 0 ] <- BetaSM;
	      } else {
	         Lars2Fit2 <- Lars2Fit;
	      }

	    T2 <- proc.time();
	        Lars2FitTime <- T2 - T1;
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", starting with HitLars"));
	          flush.console();
	      }
	 	    HitLars <- LarsCC2(SMS$XX,SMS$YY, BetStart = -1, lambda = 1/ BetaBarMin * 
	                log((1-SMS$puse)/SMS$puse), BetOlds = -1,
	     StandardFlag = 1);  
	     
	    T1 <- proc.time();    
	      HitLars2Ob <- lars(SMS$XX, SMS$YY); 
	      ##Cp <-  vector("numeric", length(HitLars2Ob$beta[,1]));
	      ##for (jj in 1:length(Cp)) {
	      ##  Cp[jj] <- ( sum ( (SMS$YY - SMS$XX %*% as.vector(HitLars2Ob$beta[jj,])  )^2 ) / SMS$SigmaSqNoise 
	      ##          - length(SMS$YY) + 2 * length(HitLars2Ob$beta[jj, HitLars2Ob$beta[jj,] != 0]));
	      ##}
	      Cp <- HitLars2Ob$Cp;
	      if (any(is.nan(Cp)) || any(is.na(Cp))) {
		       Cp <-  vector("numeric", length(HitLars2Ob$beta[,1]));
	          for (jj in 1:length(Cp)) {
	              Cp[jj] <- ( sum ( (SMS$YY - SMS$XX %*% as.vector(HitLars2Ob$beta[jj,])  )^2 ) / SMS$SigmaSqNoise 
	                   - length(SMS$YY) + 2 * length(HitLars2Ob$beta[jj, HitLars2Ob$beta[jj,] != 0]));
	          }      
	      }
	      if (any(is.nan(Cp)) || any(is.na(Cp))) {
            LarsCpFitAll <- HitLars2Ob$beta[1,] * 0 - 999;
	          LarsCpFitOnes <- LarsCpFitAll; 
	          LarsCpFitMP <- LarsCpFitAll;	      
	      } else if (length(Cp) != 1 && Cp[1] != -999) {
	          STL <- sort(Cp, index=TRUE);
		      LarsCpFitAll <- HitLars2Ob$beta[STL$ix[1],];
		      LarsCpFitOnes <- LarsCpFitAll;
		      LarsCpFitOnes[abs(LarsCpFitAll) >= MinMin] <- 1;
		      LarsCpFitOnes[abs(LarsCpFitAll) < MinMin] <- 0;
		      LarsCpFitMP <- LarsCpFitOnes;
		         
		         if (sum(LarsCpFitOnes) > 0 || length(SMS$XX[1, LarsCpFitOnes == 1]) > 0)  { 
		            XXSM <- SMS$XX[, LarsCpFitOnes == 1];
		             XXSXX <- t(XXSM) %*% XXSM;
		             if (any(is.null(XXSXX)) || length(dim(XXSXX)) != 2 
		                || dim(XXSXX)[1] != dim(XXSXX)[2] || 
		                dim(XXSXX)[1] < 1 || any(is.nan(XXSXX))) { 
		                LarsCpFitMP <- LarsCpFitAll;
		                BetaSM <- 0;
		             } else if (length(XXSXX) == 1) {
		              BetaSM <- 1 / XXSXX * sum(XXSM * SMS$YY);
		             } else if (!is.null(dim(XXSXX)) && dim(XXSXX)[1] > 1 && 
		                        dim(XXSXX)[1] == dim(XXSXX)[2] && det(XXSXX) == 0 ) {
		               try(library(corpcor, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
		               BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;
		             } else if (!is.null(dim(XXSXX))) {
		                BetaSM <- solve(t(XXSM) %*% XXSM) %*% t(XXSM) %*% SMS$YY;
		              }
		            LarsCpFitMP[LarsCpFitOnes != 0] <- BetaSM;
		         } else {
		            LarsCpFitMP <- LarsCpFitOnes;
		         }	        
	      } else   if (any(is.nan(HitLars2Ob$Cp)) || any(is.na(HitLars2Ob$Cp))) {
	          LarsCpFitAll <- HitLars2Ob$beta[1,] * 0;
	          LarsCpFitOnes <- LarsCpFitAll; 
	          LarsCpFitMP <- LarsCpFitAll;
	      } else {
		   STL <- sort(HitLars2Ob$Cp, index=TRUE);
		      LarsCpFitAll <- HitLars2Ob$beta[STL$ix[1],];
		      LarsCpFitOnes <- LarsCpFitAll;
		      LarsCpFitOnes[abs(LarsCpFitAll) >= MinMin] <- 1;
		      LarsCpFitOnes[abs(LarsCpFitAll) < MinMin] <- 0;
		      LarsCpFitMP <- LarsCpFitOnes;
		         
		         if (sum(LarsCpFitOnes) > 0 || length(SMS$XX[1, LarsCpFitOnes == 1]) > 0)  { 
		            XXSM <- SMS$XX[, LarsCpFitOnes == 1];
		             XXSXX <- t(XXSM) %*% XXSM;
		             if (any(is.null(XXSXX)) || length(dim(XXSXX)) != 2 
		                || dim(XXSXX)[1] != dim(XXSXX)[2] || 
		                dim(XXSXX)[1] < 1 || any(is.nan(XXSXX))) { 
		                LarsCpFitMP <- LarsCpFitAll;
		                BetaSM <- 0;
		             } else if (length(XXSXX) == 1) {
		              BetaSM <- 1 / XXSXX * sum(XXSM * SMS$YY);
		             } else if (!is.null(dim(XXSXX)) && dim(XXSXX)[1] > 1 && 
		                        dim(XXSXX)[1] == dim(XXSXX)[2] && det(XXSXX) == 0 ) {
		               try(library(corpcor, warn.conflicts=TRUE, quietly=TRUE), silent=TRUE);
		               BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;
		             } else if (!is.null(dim(XXSXX))) {
		                BetaSM <- solve(t(XXSM) %*% XXSM) %*% t(XXSM) %*% SMS$YY;
		              }
		            LarsCpFitMP[LarsCpFitOnes != 0] <- BetaSM;
		         } else {
		            LarsCpFitMP <- LarsCpFitOnes;
		         }
	      }
	      T2 <- proc.time();
	      LarsCpT <- T2 - T1;
	      
	    T1 <- proc.time();  
             HitLars2Ob <- lars(SMS$XX, SMS$YY);
	    ITO2 <- LarsGiveMeClose(HitLars2Ob, IntClose = min(max(round(SMS$puse * NROn), 
	         min(NNOn-1, round(.5 * NROn))), length(HitLars2Ob$beta[,1]) )
	                   , MinMin = .00001);  
	    if ( length(SMS$XX[1,]) < ITO2) {
	        ITO2 <- length(SMS$XX[1,]) - 1;
	    } else  if (length(HitLars2Ob$beta) == 1 || dim(HitLars2Ob$beta)[2] <= 1) {
	       ITO2 <- -999;
	    } else if ( length(HitLars2Ob$beta[,1]) < ITO2) {
	        ITO2 <- -999;
	    }
	    if (ITO2  > 0) {
	      Hitb2 <- HitLars2Ob$beta[ITO2,];
	      MtMtM <- matrix(as.numeric(Hitb2), length(Hitb2), 1);
	      if (dim(SMS$XX)[2] != length(MtMtM)) {
	        print( paste(" dim(SMS$XX)[2] = ", dim(SMS$XX)[2], " and ", 
	             "length(MtMtM) = ", length(MtMtM), sep="") );
	      }
	         TDOALLX <- t(SMS$XX)%*%(SMS$YY - SMS$XX %*% MtMtM)
	         TDONX <- sum(TDOALLX[ Hitb2 != 0 ] );
	         BTTry <- abs(TDONX / length(Hitb2[Hitb2 != 0] ));
	    } else {
	         BTTry <- 1 / (sd(as.vector(SMS$YY)) / max(apply(SMS$XX,2,sd))	 );
	    }  
	    T2 <- proc.time();
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", starting with lars, HitLars2Ob"));
	          flush.console();
	      }	      
	    ##print(paste("ii4 = ", ii4, ", starting with lars"));
	     T1 <- proc.time();
		       HitLars2Ob <- lars(SMS$XX, SMS$YY);
		       q <- PPract;
		       lambdaLook <- 1 * (1-q) / q  * 2 * sqrt( 2 * SMS$SigmaSqNoise);
		       LTS <- max(HitLars2Ob$lambda);
		       if (LTS < lambdaLook) { 
		            Lars4Fit <- HitLars2Ob$beta[1,] * 0;
		            Lars4Fit4A <- HitLars2Ob$beta[1,] * 0;
		       } else { 
		           ITO22 <- (1:length(HitLars2Ob$lambda))[HitLars2Ob$lambda > lambdaLook];
		           ITOOn <- (1:length(HitLars2Ob$lambda))[HitLars2Ob$lambda == lambdaLook];
			       ITO33 <- (1:length(HitLars2Ob$lambda))[HitLars2Ob$lambda < lambdaLook];
			       minLambda <- min(HitLars2Ob$lambda);
			       maxLambda <- max(HitLars2Ob$lambda);
			       if (minLambda > lambdaLook) {
			          Lars4Fit <- rep(1, length(HitLars2Ob$beta[1,]));
			          Lars4Fit4A <- HitLars2Ob$beta[ length(HitLars2Ob$beta[,1]),];
			       } else if (maxLambda < lambdaLook) {
			          Lars4Fit <- rep(0, length(HitLars2Ob$beta[1,]));
			          Lars4Fit4A <- rep(0, length(HitLars2Ob$beta[1,]));
			       } else if ( length(HitLars2Ob$lambda[HitLars2Ob$lambda > 0]) <= 0) {
			          Lars4Fit <- HitLars2Ob$beta[1,] * 0;
			          Lars4Fit4A <- rep(0, length(HitLars2Ob$beta[1,]));
			       } else if (length(ITOOn) > 0) {
			          Lars4Fit <- HitLars2Ob$beta[ ITOOn[1] ];
			          Lars4Fit[abs(Lars4Fit) >= MinMin] <- 1;
			          Lars4Fit[abs(Lars4Fit) < MinMin] <- 0;
			          Lars4Fit4A <- HitLars2Ob$beta[ITOOn[1],];
			       } else if (length(ITO33) >= 1) {
			         Lars4Fit <- HitLars2Ob$beta[ min(ITO33),];
			         Lars4Fit[abs(Lars4Fit) >= MinMin] <- 1;
			         Lars4Fit[abs(Lars4Fit) < MinMin] <- 0;	
			         Lars4Fit4A <- HitLars2Ob$beta[ min(ITO33),];		       			          
			       } else if (max(ITO22)+1 <= length(HitLars2Ob$beta[,1]) ) {
			         Lars4Fit <- HitLars2Ob$beta[ max(ITO22)+1,];
			         Lars4Fit[abs(Lars4Fit) >= MinMin] <- 1;
			         Lars4Fit[abs(Lars4Fit) < MinMin] <- 0;
			         Lars4Fit4A <- HitLars2Ob$beta[max(ITO22),];
			       } else {
			          print("I Can't figure out what's going");
			          Lars4Fit <- rep(-999, length(HitLars2Ob$beta[1,]));
			       }
			   }
		       Lars4Fit4 <- Lars4Fit;
		           if (sum(Lars4Fit) > 0) { 
			             XXSM <- SMS$XX[, Lars4Fit == 1];
			             XXSXX <- t(XXSM) %*% XXSM;
			             if (any(is.null(XXSXX)) || length(dim(XXSXX)) != 2 || dim(XXSXX)[1] != dim(XXSXX)[2] || 
	                           dim(XXSXX)[1] < 1 || any(is.nan(XXSXX))) { 
	                           Lars4Fit4 <- Lars4Fit;
	                           BetaSM <- 0;
	                    } else if (length(XXSXX) == 1) {
	                           BetaSM <- 1 / XXSXX * sum(XXSM * SMS$YY);
		                } else if (!is.null(dim(XXSXX)) && dim(XXSXX)[1] > 1 && dim(XXSXX)[1] == dim(XXSXX)[2] && 
		                            det(XXSXX) == 0 ) {
		                 try(library(corpcor, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
		                 BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;
		                } else if (!is.null(dim(XXSXX))){
		                  BetaSM <- solve(t(XXSM) %*% XXSM) %*% t(XXSM) %*% SMS$YY;
		                }
			             Lars4Fit4[Lars4Fit != 0 ] <- BetaSM; 	       
		            } else {
		               Lars4Fit4 <- Lars4Fit;
		            }
           T2 <- proc.time();
           TLars4 <- T2 - T1;	   
	                        
	       ##USENOISE <- sqrt(SMS$SigmaSqNoise);  ##USENOISE <- SMS$SigmaSqNoise
	       USENOISE <- SMS$SigmaSqNoise;
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", starting with EMRIDGE"));
	          flush.console();
	      }	 
        T1 <- proc.time();
     
        ##print(paste("ii4 = ", ii4, ", starting with EMRIDGE"));
        
	    EMCFit = EMRIDGE(yys=SMS$YY, xxs=SMS$XX, n1 = 100, n2 = 4, 
	          ppiuse = SMS$puse,
	     tauDsqStart = .4, tauAsqStart = .8, 
	     tauDMult = 2^(-1/8), tauAMult = 2^(1/8), 
	     tauAMultStop = 50, BStartBetas = -1, 
	     SigmaSqNoise = USENOISE, EFlag = FALSE,
	     logp = log(SMS$puse), log1mp = log(1-SMS$puse), 
	     SigmaNoiseInv = 1/USENOISE, NOPRINT = TRUE);

	  if (EMCFit$FailureFlag ==1 ) {
	     print("EMCFit Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
         FailToFitEMOb[[OnFails]] <<- EMCFit;	     
	  }                
	   
	  EMCFitBBs <- as.vector(EMCFit$BBHatsFinal);
	  EMCFitBBs[abs(EMCFit$BBHatsFinal) <= MinMin] = 0;
	  EMCFitBBs[abs(EMCFit$BBHatsFinal) > MinMin] = 1;
	   T2 <- proc.time();
	     EMFitT <- T2 - T1;
	   ##print(paste("ii4 = ", ii4, ", starting with HitBatWingFast"));
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", starting with HitBatWingFast"));
	          flush.console();
	      }	  
	  T1 <- proc.time();
	  if (LARSSeekFlag == 0) {
	      StlambdaA = .5/ (2*USENOISE);
	      StlambdaD = .5 / (2*USENOISE);
      
	  } else if (LARSSeekFlag >= 1) {
	      MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)
	      rnd <- round(SMS$puse * NROn * LARSSeekFlag); rnd = min( NROn, rnd);
	      ITO <- MyLarsGiveMeClose(MyFit, 
	            IntClose = rnd, MinMin = 0.000000001);   
	     StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	     StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	     ##StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE);
	     ##StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE);
	  }
	  HitBatWingFast <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, 
	    FixKa = -100, sigmaNoiseSq = USENOISE,
	    RatWant = .2, StlambdaD=StlambdaD, StlambdaA=StlambdaA, 
	    lambdaDMultC = 3^(1/2), lambdaAMultC = 3^(-1/4),
	    lambdaAmultStop = 5, TotalRuns = 5, 
	    NumEMConv = 4, MultEMCons = .99, NumCDOConv = 40, CDOEpsilon = .000001);

	  if (HitBatWingFast$FailureFlag ==1 ) {
	     print("HitBatWingFast Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
	     FailToFitEMOb[[OnFails]] <<-HitBatWingFast;
	  }  	                         
	                                                  
	  #if (length(HitBatWingFast) == 1) {
	  #  BatWingFitBBsFast = SMSReal * 0 - 10;
	  #} else if ( length(HitBatWingFast$ReturnBetas) > 1) {
	  #BatWingFitBBsFast <- as.vector(HitBatWingFast$ReturnBetas);
	  #BatWingFitBBsFast[abs(HitBatWingFast$ReturnBetas) <= 1/3] = 0;
	  #BatWingFitBBsFast[abs(HitBatWingFast$ReturnBetas) > 1/3] = 1;     
	  #} else {
	  #  BatWingFitBBsFast = SMSReal * 0 - 999;
	  #}
	  #if ((length(BatWingFitBBsFast)==1 && is.nan(BatWingFitBBsFast)) || is.nan(BatWingFitBBsFast[1])) {
	  #    BatWingFitBBsFast = (1:length(BatWingFitBBsFast)) * 0 - 999;
	  #}
	  if (length(HitBatWingFast$ReturnBetas) > 1) {
	  LarsFitBBsFast <- as.vector(HitBatWingFast$ReturnBetas);
	  LarsFitBBsFast[abs(HitBatWingFast$ReturnBetas) <= MinMin] = 0;
	  LarsFitBBsFast[abs(HitBatWingFast$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBsFast = SMSReal * 0 - 10;
	  }
	  T2 <- proc.time();
	     EM2LassoFastT <- T2 - T1;
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", starting with EM2Lasso 1"));
	          flush.console();
	      }	  	  
	     ##print(paste("ii4 = ", ii4, ", starting with longer EM2Lasso"));
	  T1 <- proc.time();
	      HitBatWing1 <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, 
	                         sigmaNoiseSq = USENOISE,
	                         RatWant = .2, StlambdaD=StlambdaD, 
	                         StlambdaA=StlambdaA, 
	                         lambdaDMultC = 2^(1/8), lambdaAMultC = 2^(-1/8),
	                         lambdaAmultStop = 50, TotalRuns = 100, 
	                         NumEMConv = 4, MultEMCons = .99, NumCDOConv = 40, CDOEpsilon = .000001,
	                         FixKa = -100);
	  if (HitBatWing1$FailureFlag ==1 ) {
	     print(paste("HitBatWing: LimitLasso: normal version, SMS$puse = ",
	                SMS$puse, " Failed ", sep="")); 
	     OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
         FailToFitEMOb[[OnFails]] <<- HitBatWing1;	     
	  }  	                         
	  #if (length(HitBatWing1) == 1) {
	  #  BatWingFitBBs1 = SMSReal * 0 - 10;
	  #} else if ( length(HitBatWing1$ReturnBetas) > 1) {
	  #BatWingFitBBs1 <- as.vector(HitBatWing1$ReturnBetas);
	  #BatWingFitBBs1[abs(HitBatWing1$ReturnBetas) <= MinMin] = 0;
	  #BatWingFitBBs1[abs(HitBatWing1$ReturnBetas) > MinMin] = 1;     
	  #} else {
	  #  BatWingFitBBs1 = SMSReal * 0 - 999;
	  #}
	  #if ((length(BatWingFitBBs1)==1 && is.nan(BatWingFitBBs1)) || is.nan(BatWingFitBBs1[1])) {
	  #    BatWingFitBBs1 = (1:length(BatWingFitBBs1)) * 0 - 999;
	  #}
	  if (length(HitBatWing1$ReturnBetas) > 1) {
		  LarsFitBBs1 <- as.vector(HitBatWing1$ReturnBetas);
		  LarsFitBBs1[abs(HitBatWing1$ReturnBetas) <= MinMin] = 0;
		  LarsFitBBs1[abs(HitBatWing1$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBs1 = SMSReal * 0 - 10;
	  }
	   T2 <- proc.time();
	   LarsFit1T <- T2 - T1;
	  T1 <- proc.time();
	      FixKa = length( SMS$BReal[ abs(SMS$BReal) != 0 ] );
	      HitBatWingFD <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = -SMS$puse, 
	                         FixKa = FixKa,
	                         sigmaNoiseSq = USENOISE,
	                         RatWant = .2, StlambdaD=StlambdaD/10, 
	                         StlambdaA=StlambdaA/10, 
	                         lambdaDMultC = 2^(1/8), lambdaAMultC = 2^(-1/8),
	                         lambdaAmultStop = 50, TotalRuns = 100, 
	                         NumEMConv = 4, MultEMCons = .99, NumCDOConv = 40, CDOEpsilon = .000001);
	  if (HitBatWingFD$FailureFlag ==1 ) {
	     print(paste("HitBatWingFD Failed, we tried FixKa = ", HitBatWingFD$FixKa, sep="")); 
	     OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
         FailToFitEMOb[[OnFails]] <<- HitBatWing1;	     
	  }  	                         
	  #if (length(HitBatWingFD) == 1) {
	  #  BatWingFitBBsFD = SMSReal * 0 - 10;
	  #} else if ( length(HitBatWing1$ReturnBetas) > 1) {
	  #BatWingFitBBsFD <- as.vector(HitBatWingFD$ReturnBetas);
	  #BatWingFitBBsFD[abs(HitBatWingFD$ReturnBetas) <= MinMin] = 0;
	  #BatWingFitBBsFD[abs(HitBatWingFD$ReturnBetas) > MinMin] = 1;     
	  #} else {
	  #  BatWingFitBBs1 = SMSReal * 0 - 999;
	  #}
	  #if ((length(BatWingFitBBs1)==1 && is.nan(BatWingFitBBs1)) || is.nan(BatWingFitBBs1[1])) {
	  #    BatWingFitBBs1 = (1:length(BatWingFitBBs1)) * 0 - 999;
	  #}
	  if (length(HitBatWingFD$ReturnBetas) > 1) {
		  LarsFitBBsFD <- as.vector(HitBatWingFD$ReturnBetas);
		  LarsFitBBsFD[abs(HitBatWingFD$ReturnBetas) <= MinMin] = 0;
		  LarsFitBBsFD[abs(HitBatWingFD$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBsFD = SMSReal * 0 - 10;
	  }
	   T2 <- proc.time();
	   LarsFitFDT <- T2 - T1;	   
	  ##print(paste("ii4 = ", ii4, ", starting with HitBatWingPI"));
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", starting with EMPI2Lasso"));
	          flush.console();
	      }	 	  
	   T1 <- proc.time();	                         
	     HitBatWingPI <- EMPI2Lasso(xxs=SMS$XX, yys=SMS$YY, pPpiuseSt = SMS$puse, 
	                         sigmaNoiseSq = USENOISE,
	                         RatWant = .2, StlambdaD= StlambdaD, 
	                         StlambdaA=StlambdaA, 
	                         lambdaDMultC = 2^(1/8), lambdaAMultC = 2^(-1/8),
	                         lambdaAmultStop = 50, TotalRuns = 100, 
	                         NumEMConv = 4, MultEMCons = .99,
	                         m1 = musize * SMS$puse, m2 = musize * (1-SMS$puse),
	                         SigmaEta = SigmaEta, SigmaBarSq = SigmaBarS,
	                         StandardFlag = 0, NumCDOConv = 50, CDOEpsilon = .000001);	
	  if (HitBatWingPI$FailureFlag ==1 ) {
	     print("HitBatWingPI Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
         FailToFitEMOb[[OnFails]] <<- HitBatWingPI;	    	     
	  }  
	  if (length(HitBatWingPI) == 1) {
	    BatWingFitBBsPI = SMSReal * 0 - 10;
	  } else if ( length(HitBatWingPI$ReturnBetas) > 1) {
		  BatWingFitBBsPI <- as.vector(HitBatWingPI$ReturnBetas);
		  BatWingFitBBsPI[abs(HitBatWingPI$ReturnBetas) <= MinMin] = 0;
		  BatWingFitBBsPI[abs(HitBatWingPI$ReturnBetas) > MinMin] = 1;     
	  } else {
	    BatWingFitBBsPI = SMSReal * 0 - 999;
	  }
	  if ((length(BatWingFitBBsPI)==1 && is.nan(BatWingFitBBsPI)) || is.nan(BatWingFitBBsPI[1])) {
	      BatWingFitBBsPI = (1:length(BatWingFitBBsPI)) * 0 - 999;
	  }
	  if (length(HitBatWingPI$ReturnBetas) > 1) {
		  LarsFitBBsPI <- as.vector(HitBatWingPI$ReturnBetas);
		  LarsFitBBsPI[abs(HitBatWingPI$ReturnBetas) <= MinMin] = 0;
		  LarsFitBBsPI[abs(HitBatWingPI$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBsPI = SMSReal * 0 - 10;
	  }
	    T2 <- proc.time();
	    LarsFitPIT <- T2 - T1;
	    ##print(paste("ii4 = ", ii4, ", starting with GiveGiveMes"));
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", starting with GiveGiveMes"));
	          flush.console();
	      }	 
	     T1 <- proc.time();
	     XLFit <- XLOperation(SMS$XX, SMS$YY, SMS$puse, SMS$SigmaSqNoise, 
	         TauOther = 60)
	     XLFitBBs <- XLFit$BBHatsFinal;
	     XLFitBBs[abs(XLFit$BBHatsFinal) <= MinMin] = 0;
	     XLFitBBs[abs(XLFit$BBHatsFinal) > MinMin] = 1;   
	     XLFitBBs <- as.vector(XLFitBBs);
	       XLFitBetas <- XLFit$BBHatsFinal * XLFit$ReturnBetas;
	       ##if (sum(XLFitBBs) == 0) {
	       ##    XLFitBetas =XLFitBBs;
	       ##} else if (sum(XLFitBBs) < length(SMS$YY) && (sum(XLFitBBs > 0) && det(t(SMS$XX[, XLFitBBs]) %*% SMS$XX[,XLFitBBs]) > 0)) {
	       ##    XLFitBetas <- XLFitBBs;
	       ##    XLFitBetas[XLFitBBs == 1] <- pseudoinverse(t(SMS$XX[, XLFitBBs]) %*% SMS$XX[,XLFitBBs]) %*% 
	       ##                          t(SMS$XX[,XLFitBBs]) %*%  SMS$YY;
	       ##} else {
	       ##    XLFitBetas <- XLFitBBs;
	       ##    XLFitBetas[XLFitBBs == 1] <-  (t(SMS$XX[, XLFitBBs]) %*% SMS$YY ) / 
	       ##        ( diag(t(SMS$XX[, XLFitBBs]) %*% SMS$XX[,XLFitBBs]) );
	       ##}
	            
	     T2 <- proc.time();
	     XLFitPIT <- T2 - T1;
	     T1 <- proc.time();
	     XLPostFit <- XLPostAnalysis(BetaProposed = HitBatWing1$ReturnBetas, 
	          xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, sigmaNoiseSq = SMS$SigmaSqNoise, 
	         tausqA = 60)
	     XLPostFitBBs <- XLPostFit$BBHatsFinal;
	     XLPostFitBBs[abs(XLPostFit$BBHatsFinal) <= MinMin] = 0;
	     XLPostFitBBs[abs(XLPostFit$BBHatsFinal) > MinMin] = 1;  
	     XLPostFitBBs <- as.vector(XLPostFitBBs);	     
	     XLPostFitBetas <- XLPostFit$ReturnBetas; 	            
	     T2 <- proc.time();
	     XLPostFitPIT <- T2 - T1;	     
	     
	               	    
	     GG1 <- GiveGiveMe( 
	        Arrs = rbind(SMSReal, Lars2Fit, LarsCpFitOnes, Lars4Fit,
                 EMCFitBBs, LarsFitBBsFast, LarsFitBBs1, LarsFitBBsPI, XLFitBBs,
                  LarsFitBBsFD, XLPostFitBBs) 
	     );
	     GG2 <- GiveGiveMeDiffD( rbind(SMS$BReal,
	               Lars2Fit2, LarsCpFitMP, Lars4Fit4, 
	               EMCFit$ReturnBetas,
	               HitBatWingFast$ReturnBetas, 
	               HitBatWing1$ReturnBetas, HitBatWingPI$ReturnBetas, as.vector(XLFitBetas),
	               HitBatWingFD$ReturnBetas, as.vector(XLPostFitBetas)));
	    TimeVec <- rbind(Lars2FitTime, LarsCpT, TLars4, EMFitT, EM2LassoFastT, LarsFit1T, 
	                                 LarsFitPIT, XLFitPIT, LarsFitFDT, XLPostFitPIT);
	    if (PrintOutFlags >= 1) {
	      plot(TimeVec, main=paste("Iteration = ", ii4, sep=""), 
	         xlab=paste("NN = ", length(SMS$YY), ", NR =", length(SMS$XX[1,]), 
	            " NoiseSq = ", round(SMS$SigmaSqNoise,2)),
	         ylab=paste("BReal = c(", round(SMS$BReal[1],1), ", ", round(SMS$BReal[2],1), ", ", 
	            round(SMS$BReal[3],1), ", ", round(SMS$BReal[4],1), ")", sep=""));
	    }
	    ##NN = NNOn, NR = NROn, 
		##	     Noise = NoiseOn, BRealVecs = BRealVecs, ppCov = ppCov, 
		##	     sigmaFactCov = sigmaFactCov
	    TimeVec <- TimeVec[, 1:3];
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", putting in ReturnMatrix"));
	          flush.console();
	      }	 	    
	   ReturnMatrix[ii4,] <- c(ppiuse, sigmaFactCov, ppCov, NNOn, NROn, NoiseOn, NumReps, 
	                       t(GG1), GG2, TimeVec);
	   if (!is.null(FileName )) {
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", Writing a table"));
	          flush.console();
	      }	 	   
	     ## write.table(ReturnMatrix, FileName, sep=", ", eol="\n", na = "NA", quote=FALSE, 
         ##        row.names=FALSE, col.names=colnames(ReturnMatrix));
          maxp <- max(2,ii4);
	      write.table(ReturnMatrix[1:maxp,], FileName, sep=", ", eol="\n", na = "NA", quote=FALSE, 
                 row.names=FALSE, col.names=colnames(ReturnMatrix[1:maxp,]));         
       }
	      if (PrintOutFlags >= 2) {  
	        print(paste("ii4 = ", ii4, ", Moving to next ii4"));
	          flush.console();
	      }	        
    }
	      if (PrintOutFlags >= 1) {  
	        print(paste("ii4 = ", ii4, ", returning final ReturnMatrix"));
	          flush.console();
	      }	     
    return(ReturnMatrix);
}

######################################
##  meanall
##
##   Classic mean of each vector.
 meanall <- function(MMT) {
    vecR <- vector("numeric", length(MMT[1,]));
     for (ii in 1:length(vecR)) {
       vecR[ii] <- mean(MMT[,ii]);
     }
     return(vecR);
 }


##############################################################################
##  PrintMyPlotter
##
##   Designed to print out in 2d space depictions of accuracy in hitting right active set
##   (No wrongly included and no wrongly eliminated.
##   The first version is targeted towards a 4 active vector set
##   The second version is targeted toward a 7 active beta set.
 PrintMyPlotter <- function(GFMAA,  OneVV, NumReps = 100, PrMeVec, TMM = FALSE,
     colList, AlphList, xlims = -999, ylims = -999) {
   ##rbind(SMSReal, 
   ##         L0Fit2, EMCFitBBs, BatWingFitBBs, LarsFitBBs, Lars2Fit, XLFitBBs );
   cto = 0;
   if ( length(OneVV[,1]) == 12) {
      iin <- 3;
      jjn <- 4;
      par(mfrow=c(3,4));
   } else if ( length(OneVV[,1]) == 8) {
      iin <- 2;
      jjn <- 4;
      par(mfrow=c(2,4));
   } else if ( length(OneVV[,1]) == 4) {
      iin <- 2;
      jjn <- 2;
      par(mfrow=c(2,2));
   } else if ( length(OneVV[,1]) == 9) {
      iin <- 3;
      jjn <- 3;
      par(mfrow=c(3,3));
   } else {
      return;
   }   
       for (ii in 1:iin) {
          for (jj in 1:jjn) {
                      cto <- cto +1;
            if (length(xlims) == 1 && xlims == -999) {
               xlimss = c(0,4);
            } else if (length(xlims) == 1) { 
               xlimss = c(0, xlims);
            } else if (length(xlims) == 2) {
               xlimss = c(xlims[1], xlims[2]);
            } else if (length(xlims) < length(OneVV[,1])) {
               xlimss = c(0,4);
            } else if (length(xlims) == length(OneVV[,1])) {
               xlimss = c(0, xlims[  cto]);
            } else if (length(xlims) < 2 * length(OneVV[,1])) {
               xlimss = c(0,4);
            } else if (length(xlims) == 2 * length(OneVV[,1])) {
               xlimss = xlims[  cto ,];
            } else { 
               xlimss = c(0,4);
            }
            if (length(ylims) ==1 && ylims == -999) {
               ylimss =  c(0, min(c(6, max(OneVV[cto,5] - 4,1))))
            } else if (length(ylims) == 1) { 
               ylimss = c(0, ylims);
            } else if (length(ylims) == 2) {
               ylimss = c(ylims[1], ylims[2]);
            } else if (length(ylims) < length(OneVV[,1])) {
               ylimss =  c(0, min(c(6, max(OneVV[cto,5] - 4,1))))
            } else if (length(ylims) == length(OneVV[,1])) {
               ylimss = c(0, ylims[  cto + jj]);
            } else if (length(ylims) < 2 * length(OneVV[,1])) {
               ylimss =  c(0, min(c(6, max(OneVV[cto,5] - 4,1))))
            } else if (length(ylims) == 2 * length(OneVV[,1])) {
               ylimss = ylims[  cto ,];
            } else { 
               ylimss =  c(0, min(c(10, max(OneVV[cto,5] - 4,1))))
            }            
            
            

            SubSPlot <- GFMAA[ GFMAA[,1] == OneVV[cto,1] & GFMAA[,2] == OneVV[cto,2] & 
                   GFMAA[,3] == OneVV[cto,3] & GFMAA[,4] == OneVV[cto,4] & 
                   GFMAA[,5] == OneVV[cto,5] & GFMAA[,6] == OneVV[cto,6], ];
            if (max(meanall(SubSPlot[, 7 + ALT + PrMeVec])) > max(ylimss)) {
                ylimss[2] = max(meanall(SubSPlot[, 7 + ALT + PrMeVec]));
            }
            plot(x=c(0,1), y=c(0,1), type="n", 
                   main=paste("n = ", OneVV[cto,4], ", kppa = ", OneVV[cto,5], 
                              ", sig^2 = ", OneVV[cto,6], sep=""), ylab="", xlab="",
                   xlim=xlimss, ylim = ylimss
                 );
            for (tt in 1:length(PrMeVec)) {
                if (TMM == FALSE) {
                   PrV1 <-SubSPlot[, 7 + PrMeVec[tt]];
                   PrV1 <- PrV1[!is.na(PrV1) & PrV1 >= 0 ];
                   PrV2 <-SubSPlot[, 7 + ALT + PrMeVec[tt]];
                   PrV2 <- PrV2[!is.na(PrV2) & PrV2 >= 0 ];                   
                } else {
                   PrV1 <-SubSPlot[, 7 + PrMeVec[tt]];
                   PrV1[is.na(PrV1) | PrV1 < 0 ] <- max(PrV1[!is.na(PrV1) & PrV1 >= 0 ]);
                   PrV2 <-SubSPlot[, 7 + ALT + PrMeVec[tt]];
                   PrV2[is.na(PrV2) | PrV2 < 0]  <- max(PrV1[!is.na(PrV2) & PrV2 >= 0 ]);                                
                }
                points( x=mean(PrV1), y= mean(PrV2), col= colList[tt],
                         pch=16, lwd=8 ); 
                points( x=mean(PrV1)+.01, y= mean(PrV2)-.01, col= "black",
                         pch=AlphList[tt] );                                       
            }
        }
     }
  return(1);
 } 
 

#################################################################################3
##   PrintHists2
##
##   Prints in Histograms counts of missing and ill assessed data.
##
 Pnames <- c("Sim Real", "L0 Minimizer", "Limit-Ridge", "Limit-Lasso", "Lasso-Wide", 
             "Lasso-Fixed", "Posterior Median")
 PrintHists2 <- function(GFMA1, TCX = c(3,4,6), TCY = c(1,2,3), SixTeller = FALSE) {
    par(mfrow=c(length(TCY), length(TCX)) );
    for (iix in 1:length(TCX)) {
         for (jjy in 1:length(TCY)) {
        if (TCY[jjy] == 1) {
           DatALL <- as.vector(GFMA1[, 7 + TCX - 1]);
        } else  if (TCY[jjy] == 2) {
           DatALL <- as.vector(GFMA1[, 7 + ALT + TCX - 1]);       
        } else {
           DatALL <- as.vector(GFMA1[, 7 + ALT*2 + TCX - 1]);              
        }
            NameT <- Pnames[TCX[iix]];
            if (TCY[jjy] == 1) { 
                DatT <- GFMA1[ ,7 + TCX[iix]-1];
                yaxname <- "Type II missing"
                bks <- -.5 + -20:100
            } else if (TCY[jjy] == 2) {
                DatT <- GFMA1[ , 7 + ALT + TCX[iix]-1];
                yaxname <- "Type I False +"
                bks <- -.5 + -20:100
            } else if (TCY[jjy] == 3) {
                if (SixTeller == TRUE && TCX[iix] == 6) {
                   DatT <- GFMA1[, 7 + 3 * ALT + 1 ];
                } else { 
                   DatT <- GFMA1[, 7 + 2 * ALT + TCX[iix] - 1];
                }
                yaxname <- "SumSqDiff";
                bks <- -5 + (0:30)*5;
            }
            if (is.null(bks)) {
               hist(DatT, main="", xlab=NameT, ylab= yaxname, xlim=c(-.5, max(DatALL)));
               TMM <- mean(DatT);
               lines(x=c(TMM,TMM), c(0,10000), lwd=3, col="red", lty=3);
            } else {
               TMM <- mean(DatT);
               hist(DatT, main="", xlab=NameT, ylab= yaxname, breaks = bks, 
                  xlim=c(-.5, max(DatALL)+.5));
               lines(x=c(TMM,TMM), c(0,10000), lwd=3, col="red", lty=3);                  
            } 
         }
     }
    return;
}


###################################################################################
##  PrintHistsV2
##
##   Assess histograms representing counts of errors in each dataset against Real set
##
PnamesV2 <- c("Sim Real", "Lasso-Fixed 7", "Lasso-Cp", "Lasso Min-Yuan", 
              "Limit-Ridge", "Two-Lasso - Fast", "Limit-2-Lasso", "Limit-2PI-Lasso",
              "XL Median", "Fermi Limit-Lasso", "XLM Lasso");
 PrintHistsV2 <- function(GFMA1, TCX = c(3,4,6), TCY = c(1,2,3), SixTeller = FALSE) {
  ## if(is.null(ALT2)) {
    DeclareGlobals();
  ## }
    ##ALT2 <<- 10 * 11 / 2;
    par(mfrow=c(length(TCY), length(TCX)) );
    for (iix in 1:length(TCX)) {
         for (jjy in 1:length(TCY)) {
        if (TCY[jjy] == 1) {
           DatALL <- as.vector(GFMA1[, 7 + TCX - 1]);
        } else  if (TCY[jjy] == 2) {
           DatALL <- as.vector(GFMA1[, 7 + ALT2 + TCX - 1]);       
        } else {
           DatALL <- as.vector(GFMA1[, 7 + ALT2*2 + TCX - 1]);              
        }
            NameT <- Pnames[TCX[iix]];
            if (TCY[jjy] == 1) { 
                DatT <- GFMA1[ ,7 + TCX[iix]-1];
                yaxname <- "Type II missing"
                bks <- -.5 + -20:100
            } else if (TCY[jjy] == 2) {
                DatT <- GFMA1[ , 7 + ALT2 + TCX[iix]-1];
                yaxname <- "Type I False +"
                bks <- -.5 + -20:100
            } else if (TCY[jjy] == 3) {
                if (SixTeller == TRUE && TCX[iix] == 6) {
                   DatT <- GFMA1[, 7 + 3 * ALT2 + 1 ];
                } else { 
                   DatT <- GFMA1[, 7 + 2 * ALT2 + TCX[iix] - 1];
                }
                yaxname <- "SumSqDiff";
                bks <- -5 + (0:30)*5;
            }
            if (is.null(bks)) {
               hist(DatT, main="", xlab=NameT, ylab= yaxname, xlim=c(-.5, max(DatALL)));
               TMM <- mean(DatT);
               lines(x=c(TMM,TMM), c(0,10000), lwd=3, col="red", lty=3);
            } else {
               TMM <- mean(DatT);
               hist(DatT, main="", xlab=NameT, ylab= yaxname, breaks = bks, 
                  xlim=c(-.5, max(DatALL)+.5));
               lines(x=c(TMM,TMM), c(0,10000), lwd=3, col="red", lty=3);                  
            } 
         }
     }
    return;
}   

############################################################################
##  GiveMeIDOnOutput
##     Helperfunction locates specific indexes for PrintMyPlotter functions
##
##   
 GiveMeIDOnOutput <- function(GFMAA, OneVV, VNum, PrMeNum) {
  ##if(is.null(ALT2)) {
    DeclareGlobals();
  ##}
  ##ALT2 <- 10 * 11 / 2;
  cto = VNum
  SubSPlot <- GFMAA[ GFMAA[,1] == OneVV[cto,1] & GFMAA[,2] == OneVV[cto,2] & 
                   GFMAA[,3] == OneVV[cto,3] & GFMAA[,4] == OneVV[cto,4] & 
                   GFMAA[,5] == OneVV[cto,5] & GFMAA[,6] == OneVV[cto,6], ];
           
                   PrV1 <-SubSPlot[, 7 + VNum];
                   ##PrV1 <- PrV1[!is.na(PrV1) & PrV1 >= 0 ];
                   PrV2 <-SubSPlot[, 7 + ALT2 + VNum];
                   ##PrV2 <- PrV2[!is.na(PrV2) & PrV2 >= 0 ];  
                   PrV3 <-SubSPlot[, 7 + 2*ALT2 + VNum];                                    
         return(cbind(PrV1, PrV2, PrV3))        
 }

######################################################################
##   PrintMyPlotterV2  
##
##   Useful in printing out data in Points distance format.
##
##
##
## 
 PrintMyPlotterV2 <- function(GFMAA,  OneVV, NumReps = 100, PrMeVec, TMM = FALSE,
     colList, AlphList, xlims = -999, ylims = -999) {
##PnamesV2 <- c("Sim Real", "Lasso-Fixed 7", "Lasso-Cp", "Lasso Min-Yuan", 
##              "Limit-Ridge", "2-Lasso-Fast", "Limit-2-Lasso", "Limit-2PI-Lasso", "Posterior Median",
##              "Limit-Lasso FD);
  PnV2 <- c("Lasso-Fixed", "Lasso-Cp", "Lasso w=1 LY", 
              "Limit-Ridge", "2-Lasso-9X", "Limit-2-Lasso", "Limit-2PI-Lasso", "Marg Median",
              "Limit-Lasso FD", "MargM Lim-Lasso");
  PnV3 <- c("Lasso-\nFixed", "Lasso-\nCp", "Lasso \nw=1 LY", 
              "Limit-\nRidge", "2-Lasso-\n9X", "Limit-2\n-Lasso", "Limit-2PI\n-Lasso", "Marg\n Median",
              "Limit-\nLasso FD", "MargM Lim\n-Lasso");              
   cto = 0;
   ALT2 <- 10 * 11 / 2;
   if ( is.vector(OneVV)) {
        iin <- 1; jjn <- 1;
        par(mfrow=c(1,1));
          OneVV <- t( matrix( c(OneVV, OneVV), length(OneVV), 2));
   } else if ( length(OneVV[,1]) == 12) {
      iin <- 3;
      jjn <- 4;
       matlayout <- rbind(
     c( 24,15,25 ,16,26 ,17,27),
     c(1,2,3,4,5,6,7),
     c(28 ,18,29 ,19,30 ,20,31),
     c(8,9,10,11,12,13,14),
     c(32 ,21,33 ,22,34 ,23,35)   );
nf <- layout(matlayout, c(1,.02, 1,.02, 1,.02, 1), c(1,.02, 1, .02,1), TRUE)
       matlayout <- rbind(
     c(1, 25,16,26 ,17,27 ,18,28),
     c(1, 2,3,4,5,6,7,8),
     c(1, 29 ,19,30 ,20,31 ,21,32),
     c(1, 9,10,11,12,13,14,15),
     c(1, 33 ,22,34 ,23,35 ,24,36)   );
nf <- layout(matlayout, c(.5, 1,.02, 1,.02, 1,.02, 1), c(1,.02, 1, .02,1), TRUE)
  par(plt=c(0,1,0,1));
  plot(1,1, xlab="", ylab="", axes=FALSE, type="n", ylim=c(0,1), xlim=c(0,1));
      nNames <- PnV3[PrMeVec];
      ylen <- length(AlphList) + 4;
      yspots <- (ylen - .5  - 1:length(AlphList))/ylen;
    for (ii in 1:length(AlphList)) {
      points(x=.09, y=yspots[ii], pch=19, lwd=17-ii, col=colList[ii]);
      points(x=.092, y=yspots[ii]+.000, pch=AlphList[ii], col="black");
      text( x=.55, y=yspots[ii]-.001, labels=paste(AlphList[ii], "-", nNames[ii]) );
    }
    text(x=.5, y=yspots[length(AlphList)]-1/ylen, labels=("X-axis = Type II\nTrue Missed"));
    text(x=.5, y=yspots[length(AlphList)]-2.5/ylen, labels=("Y-axis = Type I\nFalse Pos"));    
##layout.show(nf)
    par(plt=c(0,1,0,1));
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8), type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
     par(plt=c(0,1,0,1));
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8), type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       par(plt=c(0,1,0,1)); 
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");      
      par(plt=c(.13,.975,.155,.83))
   } else if ( length(OneVV[,1]) == 8) {
      iin <- 2;
      jjn <- 4;
      par(mfrow=c(2,4));
   } else if ( length(OneVV[,1]) == 4) {
      iin <- 2;
      jjn <- 2;
      par(mfrow=c(2,2));
   } else if ( length(OneVV[,1]) == 9) {
      iin <- 3;
      jjn <- 3;
      par(mfrow=c(3,3));
             matlayout <- rbind(
     c( 24,15,25 ,16,26 ,17,27),
     c(1,2,3,4,5,6,7),
     c(28 ,18,29 ,19,30 ,20,31),
     c(8,9,10,11,12,13,14),
     c(32 ,21,33 ,22,34 ,23,35)   );
nf <- layout(matlayout, c(1,.02, 1,.02, 1,.02, 1), c(1,.02, 1, .02,1), TRUE)
       matlayout <- rbind(
     c(1, 18,12,19 ,13,20),
     c(1, 2,3,4,5,6),
     c(1, 21 ,14,22 ,15,23),
     c(1, 7,8,9,10,11),
     c(1, 24 ,16,25 ,17,26)   );
nf <- layout(matlayout, c(.5, 1,.02, 1,.02, 1), c(1,.02, 1, .02,1), TRUE)
  par(plt=c(0,1,0,1));
  plot(1,1, xlab="", ylab="", axes=FALSE, type="n", ylim=c(0,1), xlim=c(0,1));
      nNames <- PnV3[PrMeVec];
      ylen <- length(AlphList) + 4;
      yspots <- (ylen - .5  - 1:length(AlphList))/ylen;
    for (ii in 1:length(AlphList)) {
      points(x=.09, y=yspots[ii], pch=19, lwd=17-ii, col=colList[ii]);
      points(x=.092, y=yspots[ii]+.000, pch=AlphList[ii], col="black");
      text( x=.55, y=yspots[ii]-.001, labels=paste(AlphList[ii], "-", nNames[ii]) );
    }
    text(x=.5, y=yspots[length(AlphList)]-1/ylen, labels=("X-axis = Type II\nTrue Missed"));
    text(x=.5, y=yspots[length(AlphList)]-2.5/ylen, labels=("Y-axis = Type I\nFalse Pos"));    
##layout.show(nf)
    par(plt=c(0,1,0,1));
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8), type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    ##plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    ##plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
     par(plt=c(0,1,0,1));
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8), type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    ##plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    ##plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       par(plt=c(0,1,0,1)); 
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      ##plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
      ##    type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      ##plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
      ##    type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      ##plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
      ##    type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");      
      par(plt=c(.13,.975,.155,.83))
   } else if ( length(OneVV[,1]) == 16) {
      iin <- 4; jjn <- 4;
      par(mfrow=c(4,4));
   } else {
      return;
   }   
       for (ii in 1:iin) {
          for (jj in 1:jjn) {
                      cto <- cto +1;
            if (length(xlims) == 1 && xlims == -999) {
               xlimss = c(0,7);
            } else if (length(xlims) == 1) { 
               xlimss = c(0, xlims);
            } else if (length(xlims) == 2) {
               xlimss = c(xlims[1], xlims[2]);
            } else if (length(xlims) < length(OneVV[,1])) {
               xlimss = c(0,7);
            } else if (length(xlims) == length(OneVV[,1])) {
               xlimss = c(0, xlims[  cto]);
            } else if (length(xlims) < 2 * length(OneVV[,1])) {
               xlimss = c(0,7);
            } else if (length(xlims) == 2 * length(OneVV[,1])) {
               xlimss = xlims[  cto ,];
            } else { 
               xlimss = c(0,7);
            }
            if (length(ylims) ==1 && ylims == -999) {
               ylimss =  c(0, min(c(6, max(OneVV[cto,5] - 4,1))))
            } else if (length(ylims) == 1) { 
               ylimss = c(0, ylims);
            } else if (length(ylims) == 2) {
               ylimss = c(ylims[1], ylims[2]);
            } else if (length(ylims) < length(OneVV[,1])) {
               ylimss =  c(0, min(c(6, max(OneVV[cto,5] - 4,1))))
            } else if (length(ylims) == length(OneVV[,1])) {
               ylimss = c(0, ylims[  cto + jj]);
            } else if (length(ylims) < 2 * length(OneVV[,1])) {
               ylimss =  c(0, min(c(6, max(OneVV[cto,5] - 4,1))))
            } else if (length(ylims) == 2 * length(OneVV[,1])) {
               ylimss = ylims[  cto ,];
            } else { 
               ylimss =  c(0, min(c(10, max(OneVV[cto,5] - 4,1))))
            }            
            
            

            SubSPlot <- GFMAA[ GFMAA[,1] == OneVV[cto,1] & GFMAA[,2] == OneVV[cto,2] & 
                   GFMAA[,3] == OneVV[cto,3] & GFMAA[,4] == OneVV[cto,4] & 
                   GFMAA[,5] == OneVV[cto,5] & GFMAA[,6] == OneVV[cto,6], ];
            if (max(meanall(SubSPlot[, 7 + ALT2 + PrMeVec])) > max(ylimss)) {
                ylimss[2] = max(meanall(SubSPlot[, 7 + ALT2 + PrMeVec]));
            }
            ##plot(x=c(0,1), y=c(0,1), type="n", 
            ##       main=paste("n = ", OneVV[cto,4], ", kppa = ", OneVV[cto,5], 
            ##                  ", sig^2 = ", OneVV[cto,6], sep=""), ylab="", xlab="",
            ##       xlim=xlimss, ylim = ylimss
            ##     );
            plot(x=c(0,1), y=c(0,1), type="n",, ylab="", xlab="",
                   xlim=xlimss, ylim = ylimss
                 ); 
            title(main=substitute(list(n,kappa,sigma)== group("(",list(a, x,y),")"),              
                           list(a=OneVV[cto,4], x=OneVV[cto,5], y=OneVV[cto,6])))                                
            for (tt in 1:length(PrMeVec)) {
                if (TMM == FALSE) {
                   PrV1 <-SubSPlot[, 7 + PrMeVec[tt]];
                   PrV1 <- PrV1[!is.na(PrV1) & PrV1 >= 0 ];
                   PrV2 <-SubSPlot[, 7 + ALT2 + PrMeVec[tt]];
                   PrV2 <- PrV2[!is.na(PrV2) & PrV2 >= 0 ];                   
                } else {
                   PrV1 <-SubSPlot[, 7 + PrMeVec[tt]];
                   PrV1[is.na(PrV1) | PrV1 < 0 ] <- max(PrV1[!is.na(PrV1) & PrV1 >= 0 ]);
                   PrV2 <-SubSPlot[, 7 + ALT2 + PrMeVec[tt]];
                   PrV2[is.na(PrV2) | PrV2 < 0]  <- max(PrV1[!is.na(PrV2) & PrV2 >= 0 ]);                                
                }
                points( x=mean(PrV1), y= mean(PrV2), col= colList[tt],
                         pch=19, lwd=(17-tt) ); 
                points( x=mean(PrV1)+.01, y= mean(PrV2)-.01, col= "black",
                         pch=AlphList[tt] );                                       
            }
        }
     }
  return(1);
 }
 
 

##########################################################################
##  CoordinateDescent
##    Performs CoordinateDescent Lasso algorithm on Linear Vector Set
##    Uses CoordinateDescent.cc compiled C++ code
CoordinateDescent <- function(xx = -1, yy = -1, XTX = -1, XTY = -1, NLen = -1,
   OnBeta = -1, OnGamma = -1, OnLambda = -1,
   RecordBetasFlag = FALSE, OnGammas = -999, InitKKs = -5,
   NumCDOConv = -999, CDOEpsilon = -999, WLSWeights = -1,
   TotalLoops = -999, MaxEpsilon = -999,
   Verbose = 0) {
   
   TNLen = abs(NLen);
   ## if(!is.loaded("CoordinateDescentLassoShell2009") ) {
   ##     LoadTwoCC();
   ##  }
   XTXFlag = -1;
   if (Verbose > 0) {
     print("CoordinateDescent: Starting"); flush.console();
   }
   if (TotalLoops == -999 && MaxEpsilon == -999 && NumCDOConv == -999 && CDOEpsilon == -999) {
      print("CoordinateDescentError, No Number of Loops or Epsilon Value given\n");
      print("We'll go 200 times");
      TotalLoops = 200; MaxEpsilon = .0000001; NumCDOConv = 200;
      CDOEpsilon = .00000001;
      ##return(-1);
   } else if (TotalLoops == -999 || MaxEpsilon == -999) {
       TotalLoops = NumCDOConv; MaxEpsilon = CDOEpsilon;
   } else if (NumCDOConv == -999 && CDOEpsilon == -999) {
       NumCDOConv = TotalLoops; CDOEpsilon = MaxEpsilon;
   }
   if (is.matrix(xx) == FALSE || length(xx) == 1 || xx == -1) {
       if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == -1) {
          print("CoordinateDescent, Error no proper X given, XTX False"); return;  
       }
       if ( (is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) || length(XTY) == 1 || XTY == -1) {
          print("CoordinateDescent, Error no proper X given, XTY False"); return;
       }
       if (length(NLen) > 1 || NLen == -1) {
           print("CoordinateDescent, Error, no proper NLen given");
           return;
       }
       if (length(dim(XTX)) != 2) {
           print("CoordinateDescentError, no good dimensions to XTX");
           flush.console();
           return(-20);
       }
       XTXFlag = 1;
       XXX = as.double(XTX); YYY = as.double(XTY);
       WLSWeights = -1;
       TNLen = - NLen; kLen = length(XTX[1,]);
   } else if ((is.matrix(yy) == FALSE  && is.vector(yy) == FALSE) ||
                        length(yy) == 1 || yy == -1) {
        print("CoordinateDescent, Error no proper Y given"); return;
   } else {
       XTXFlag = 0;
       XXX = as.double(xx); YYY = as.double(yy);
       TNLen = length(yy); kLen = length(xx[1,]);
       if (length(WLSWeights) != length(yy) || min(WLSWeights) < 0) {
              WLSWeights = -1;
       }
      ## if (InitKKs > 0) {
      ##    print("Cannot do xx version if xx is non matrix\n");
      ##    return;
      ## }       
   }  
   if (OnGamma < 0 && OnLambda > 0) {
     OnGamma <- OnLambda / 2
   } else if (OnGamma < 0 && OnLambda < 0 && length(OnGammas) != kLen) {
      print("CoordinateDescentError, OnGamma not given\n"); return();
   }                                      
   if (length(OnBeta) != kLen) {
     OnBeta = rep(0, kLen);
   }  else {
     OnBeta =   as.vector(matrix(as.numeric(OnBeta), 1, length(OnBeta)));
   }
   if (length(OnGammas) == 1 && OnGammas == -1) {
      OnGammas = rep(1, kLen);  OnGammas[1] = -999;
   }
   if (RecordBetasFlag == TRUE) {
      OnRecordBetas = rep(0, (kLen+1) * (TotalLoops +3));
   } else {
      OnRecordBetas = -999;
   }
 
   InitKKs <- min(c(kLen-1, InitKKs));
   ##if (InitKKs > kLen) {InitKKs = kLen-1;}
   ###int CoordinateDescentLasso(int NLen, int kLen,  double *XXX, double *YYY,
   ###          double *OnBeta, double OnGamma, double *OnRecordBetas,
   ###          double *FinalBeta, int TotalLoops, double MaxEpsilon)
   if (InitKKs < 0) {
     if (Verbose > 0) {
       print(paste("CoordinateDescent: InitKKs : ", InitKKs, " < 0 version!"));
       flush.console();
     }
   CDO <- .C("CoordinateDescentLassoShell2009", NLen = as.integer(TNLen),
       kLen = as.integer(kLen), XXX = as.double(XXX), YYY = as.double(YYY),
       OnBeta = as.double(OnBeta), OnGamma = as.double(OnGamma), 
       OnRecordBetas = as.double(OnRecordBetas), 
       FinalBeta = as.double(OnBeta), TotalLoops = as.integer(TotalLoops),
       MaxEpsilon = as.double(MaxEpsilon), OnGammas = as.double(OnGammas),
       InitKKs = as.integer(0),
       InverseGammaConstant = as.double(-1),
       WLSWeights = as.double(WLSWeights),
       PrintFlag = as.integer(Verbose -1));
   } else {
     if (Verbose > 0) {
       print(paste("CoordinateDescent: InitKKs : ", InitKKs, " > 0 version!"));
       flush.console();
     }
   CDO <- .C("CoordinateDescentLassoShell2009", NLen = as.integer(TNLen),
       kLen = as.integer(kLen), XXX = as.double(XXX), YYY = as.double(YYY),
       OnBeta = as.double(OnBeta), OnGamma = as.double(OnGamma), 
       OnRecordBetas = as.double(OnRecordBetas), 
       FinalBeta = as.double(OnBeta), TotalLoops = as.integer(TotalLoops),
       MaxEpsilon = as.double(MaxEpsilon), OnGammas = as.double(OnGammas),
       InitKKs = as.integer(InitKKs),
       InverseGammaConstant = as.double(-1),
       WLSWeights = as.double(WLSWeights),
       PrintFlag = as.integer(Verbose -1)
      );   
   }
   if (Verbose > 0) {
     print("CoordinateDescent: Finished .C Algorithm"); flush.console();
   }
   if (length(OnRecordBetas) > 1) {
    CDO$OnRecordBetas = t( matrix(CDO$OnRecordBetas[1: ( kLen * TotalLoops)],
                                      kLen, TotalLoops) );
   }
   if (XTXFlag == 0) {
     CDO$xx = xx; CDO$yy = yy;
     CDO$XXX = XXX; CDO$YYY = YYY;
     CDO$XTY = XTY; CDO$XTX = XTX;
   } else {
     CDO$xx = xx; CDO$yy = yy;
     CDO$XXX = XXX; CDO$YYY = YYY;
     CDO$XTY = XTY; CDO$XTX = XTX;
   }
   CDO$UsedLoops = CDO$TotalLoops; CDO$TotalLoops = TotalLoops;
   CDO$ReturnBetas = CDO$FinalBeta;
   CDO$OnBeta = NULL;
   return(CDO);
}



##########################################################################
##  CoordinateDescentPrototype
##    Performs CoordinateDescent Lasso algorithm on Linear Vector Set
##    Uses only Native R code
CoordinateDescentPrototype <- function(xx = -1, yy = -1, XTX = -1, XTY = -1, NLen = -1,
   TotalLoops = 20, MaxEpsilon = .0001, OnBeta = -1, OnGamma = 1,
   RecordBetasFlag = FALSE, OnGammas = -999) {
   XTXFlag = -1;
   if (is.matrix(xx) == FALSE || length(xx) == 1 || xx == -1) {
       if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == -1) {
          print("CoordinateDescent, Error no proper X given, XTX False"); return;  
       }
       if ( (is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) || length(XTY) == 1 || XTY == -1) {
          print("CoordinateDescent, Error no proper X given, XTY False"); return;
       }
       if (length(NLen) > 1 || NLen == -1) {
           print("CoordinateDescent, Error, no proper NLen given");
           return;
       }
       XTXFlag = 1;
       XXX = XTX; YYY = XTY;
       TNLen = - NLen; kLen = length(XTX[1,]);
   } else if ( (is.matrix(yy) == FALSE && is.vector(yy) == FALSE) || length(yy) == 1 || yy == -1) {
        print("CoordinateDescent, Error no proper Y given"); return;
   } else {
       XTXFlag = 0;
       XXX = xx; YYY = yy;  XTX <- t(xx) %*% xx; XTY <- t(xx) %*% yy;
       TNLen = NLen; kLen = length(xx[1,]);
   }  
   if (length(OnBeta) != kLen) {
      OnBeta = rep(0, kLen);
   }
   if (length(OnGammas) == 1 || OnGammas == -1) {
      OnGammas = rep(1, kLen);  OnGammas[1] = -999;
      OnGammas = rep(OnGamma,kLen);
   }
   if (RecordBetasFlag == TRUE) {
      OnRecordBetas = matrix(0, TotalLoops, kLen);
   } else {
      OnRecordBetas = -999;
   }
   if ( !is.matrix(XTX) ||  XTX == -1)  {
      XTX <- t(xx) %*% xx;  XTY <- t(xx) %*% yy;
   }
   XTYResid <- XTY;
   OnBetas =  rep(0, kLen);
   PrevBetas = OnBetas;
   ##print(paste("kLen = ", kLen, " and begining OnBetas = \n", sep=""));
   ##print(OnBetas);
   for (ii in 1:TotalLoops) {
        PrevBetas = OnBetas;
        for (jj in 1:kLen) {
            PRD = XTYResid[jj] + XTX[jj,jj] *OnBetas[jj];
            PropBeta = 0;
            if ( PRD > OnGammas[jj] ) {
               PropBeta = (PRD- OnGammas[jj]) / XTX[jj,jj];
            } else if (PRD < - OnGammas[jj] ) {
               PropBeta = (PRD + OnGammas[jj]) / XTX[jj,jj];
            } else { 
               PropBeta = 0;
            }
            XTYResid = XTYResid - XTX[,jj] * (PropBeta - OnBetas[jj]);
            OnBetas[jj] = PropBeta;
        }
        if (is.matrix(OnRecordBetas)) {
          OnRecordBetas[ii,] <- OnBetas;
        }
        if (sum(abs(OnBetas - PrevBetas)) <= MaxEpsilon) {
             break;
        }
   }

   CDO <- list( NLen = NLen, kLen = kLen, xx = xx, yy = yy, XTX = XTX, XTY = XTY,
          OnBeta = OnBetas, OnGamma = OnGamma, OnRecordBetas = OnRecordBetas,
          FinalBeta = OnBetas, TotalLoops = 20, MaxEpsilon = .0003,
          OnGammas = OnGammas, ReturnBetas = OnBetas, NumLoops = ii);
   return(CDO);
}


 
 
 
 
 
 ###   PrintCode
###  Designed functions in R for printing
PrintSpreadFunction <- function (EMObject, ITFlag = TRUE, weightFlag = FALSE, 
    USFlag = FALSE, logFlag = FALSE, mfrows=c(2,1), symmet=TRUE,
    lDAxis = FALSE, headerWANT = FALSE) {
  if (EMObject$Type=="Lars") {
     print("No Functionality for Lars Object, try LARS package");
     return;
  }
 
   logScale = FALSE; revScale = FALSE;
  if (EMObject$Type == "EMLARS" || EMObject$Type == "EMLars" ||
        EMObject$Type == "EM2Lasso" ) {

     pphave = EMObject$ppiuse;
     StartD = EMObject$StlambdaD;
     StartA = EMObject$StlambdaA;
     SeqD = EMObject$LambdaDK;
     SeqA = EMObject$LambdaAK;
     SDSeqD = sqrt(2) / SeqD;
     SDSeqA = sqrt(2) / SeqA;
     sigSeq = EMObject$sigmaNoiseSq + SDSeqD * 0 ;
     if (weightFlag == TRUE) {
        BetaSeq = EMObject$USWReturnBetas;
     } else if (USFlag == TRUE) {
        BetaSeq = EMObject$USReturnBetas;
     } else {
        BetaSeq = EMObject$ReturnBetasAll;
     }
     BBSeq = EMObject$BBHatsAll;
     xx <- 1:length(BBSeq[,1]);
     piSeq = EMObject$ppiuse + SDSeqD * 0;
     MyText <- "EM2Lasso";
     ##round(sqrt(sigsq), 2)
     ##S.E = substitute(
     ##  expression(EMLARS, sigma = a, pi = b, lambda_A = c, lambda_D = d, 
     ##      MTP = e), list( 
     ##      a = round(sqrt(sigSeq[length(sigSeq)]),2),
     ##      b=round(piSeq[length(piSeq)],3), 
     ##      c = round(piSeq[length(piSeq)],3),
     ##      d =  round(1/SDSeqD[1],3),
     ##      e = round( 1/SDSeqA[1],3),
     ##      f= round((SDSeqD[1]/SDSeqD[2]),3)
     ##      ))
     Mymain=paste(MyText, "\n", ", sig=", round(sqrt(sigSeq[length(sigSeq)]), 2), 
        ", pi=", round(piSeq[length(piSeq)],3), 
        ", lbdaD[1]= ", 
        round( SeqD[1],3), ", lbdaA[1]=", round( SeqA[1],3), ", MTP=", 
        round( (SDSeqD[1]/SDSeqD[2]),3), sep="");
     if (ITFlag == FALSE && logFlag == TRUE) {
        xx <- EMObject$LambdaDK;
     } else if (ITFlag == FALSE && logFlag == FALSE) {
        xx <- EMObject$LambdaDK;  logScale = TRUE;
     }
  }
  if (EMObject$Type =="EMPILARS" || EMObject$Type == "EMPiLars" ||
     EMObject$Type=="EMPILars" ||
     EMObject$Type=="EM2PILasso" || EMObject$Type == "EM2PiLasso") {
    sigsq = EMObject$sigmaNoiseSq;
     pphave = EMObject$ppiuse;
     StartD = EMObject$StlambdaD;
     StartA = EMObject$StlambdaA;
     SeqD = EMObject$LambdaDK;
     SeqA = EMObject$LambdaAK;
     SDSeqD = sqrt(2) / SeqD;
     SDSeqA = sqrt(2)/ SeqA;
     if (weightFlag == TRUE) {
        BetaSeq = EMObject$USWReturnBetasAll;
     } else if (USFlag == TRUE) {
        BetaSeq = EMObject$USReturnBetasAll;
     } else {
        BetaSeq = EMObject$ReturnBetasAll;
     }
     BBSeq = EMObject$BBHatsAll;
       xx <- 1:length(BBSeq[,1]);
     piSeq = EMObject$PiRecVec; 
     sigSeq = EMObject$SigRecVec;
     MyText <- "EMPILARS";
     Mymain=paste(MyText, "\n", ", sigEst=", round(sqrt(sigSeq[length(sigSeq)]), 2), 
        ", piEst=", round(piSeq[length(piSeq)],3), 
        ", lbdaD[1]= ", 
        round(SeqD[1],3), ", lbdaA[1]=", round(SeqA[1],3), ", MTP=", 
        round( (SDSeqD[1]/SDSeqD[2])^2,3), sep="");
     if (ITFlag == FALSE && logFlag == TRUE) {
        xx <- EMObject$LambdaDK;
     } else if (ITFlag == FALSE && logFlag == FALSE) {
        xx <- EMObject$LambdaDK;  logScale = TRUE;
     }       
  }
  if (EMObject$Type == "EMRIDGE") {
     sigsq = EMObject$SigmaSqNoise;
     pphave = EMObject$ppiuse;
     SDSeqD = sqrt(EMObject$tauDK);
     SDSeqA = sqrt(EMObject$tauAK);
     ddBetaSeq = dim(EMObject$ReturnBetasAllEM);
       if (length(ddBetaSeq) == 3) {
          BetaSeq = t(EMObject$ReturnBetasAllEM[,ddBetaSeq[2],]);
       } else {
          BetaSeq = t(EMObject$ReturnBetasAllEM);
       }
     ddBBHatsAll = dim(EMObject$BBHatsAllEM);
        if (length(ddBBHatsAll) == 3) {
          BBSeq = t(EMObject$BBHatsAll[,ddBBHatsAll[2],]);
        } else {
          BBSeq = t(EMObject$BBHatsAll);
        }
          xx <- 1:length(BBSeq[,1]);
     piSeq = SDSeqD * 0 + pphave;
     sigSeq = SDSeqD * 0 + sigsq;
     MyText <- "EMRIDGE"
     
     Mymain=paste(MyText, "\n", ", sig=", round(sqrt(sigsq), 2), ", pi=", round(pphave,3), 
        ", tauD[1]= ", 
        round(SDSeqD[1],3), ", tauA[1]=", round(SDSeqA[1],3), ", MTP=", 
        round(SDSeqD[2]/SDSeqD[1],3), sep="");
     if (ITFlag == FALSE && logFlag == TRUE) {
        xx <- EMObject$tauDK;  revScale = TRUE;
     } else if (ITFlag == FALSE && logFlag == FALSE) {
        xx <- EMObject$tauDK;  logScale = TRUE; revScale = TRUE;
     }        
   }
    ABCD = c("A","B","C","D","E","F","G","H","I", "J", "K", "L", "M", "N","O","P",
     "Q", "R", "S","T", "U", "V", "W", "X", "Y", "Z");
    ##par(mfrow=c(2,1))

MyCols <- rainbow(length(BBSeq[1,]));
   if (length(mfrows)== 1 && mfrows==FALSE)  {
    print("FALSE mfrows")
   } else if (length(mfrows) ==1 && mfrows < 2)  {
        par(mfrow=c(2,1));
   } else if (length(mfrows) == 1 ) {
        par(mfrow=c(mfrows,1));
   } else if (mfrows[1] * mfrows[2] < 2) {
       par(mfrow=c(2,1))
   } else {
      par(mfrow=mfrows); 
   }
   if (length(headerWANT) > 1 || headerWANT != FALSE) {
     if (length(headerWANT) >=2)  {
       MyMainA = headerWANT[1];
     } else {
       MyMainA = headerWANT; 
     }
   }  else {
       MyMainA = "";
   }
   if (revScale == FALSE && logScale == TRUE) { 
     plot( x=c(min(xx), max(xx)), y=c(min(BBSeq), max(BBSeq)),
           main=MyMainA, xlab="", ylab="", type="n", log="x");
   } else if (revScale == FALSE && logScale == FALSE) {
     plot( x=c(min(xx), max(xx)), y=c(min(BBSeq), max(BBSeq)),
           main=MyMainA, xlab="", ylab="", type="n");  
   } else if (revScale == TRUE && logScale == FALSE) {   
     plot( x=c(max(xx), min(xx)), y=c(min(BBSeq), max(BBSeq)),
           main=MyMainA, xlab="", ylab="", type="n");    
   } else if (revScale == FALSE && logScale == TRUE) { 
     plot( x=c(max(xx), min(xx)), y=c(min(BBSeq), max(BBSeq)),
           main=MyMainA, xlab="", ylab="", type="n", log="x");
   }
   if (MyMainA == "") {
     padjA = -1;
     padjB = -.8;
    sigsq = EMObject$sigmaNoiseSq;
     mtext(substitute(sigma == a, 
       list(a = round(sqrt(sigsq),2)
           ) ), side=3, cex=1.25, adj = .05, padj = padjA
           )
    piA = pphave ;
     mtext(substitute(pi == a,
       list(a = round(piA,3) ) ),
        side =3, cex=1.25, adj = .26, padj = padjA );
    ##mtext("Start", cex=1.25, side=3, adj=.5, padj = padjB);
    if (EMObject$Type=="EM2Lasso" || EMObject$Type=="EM2PiLasso"
       || EMObject$Type == "EMLARS" || EMObject$Type=="EMLars" ||
       EMObject$Type == "EMPiLars")  {
    mtext(substitute(lambda[0] == a ,
       list(a = round(EMObject$LambdaDK[1],3) ) ),
        side =3, cex=1.25, adj = .55, padj = padjB );
    mtext(paste("Mlt = ", round(
         EMObject$LambdaDK[2] /
         EMObject$LambdaDK[1],3), sep=""), cex=1.25, 
         adj = .86, padj = padjA )
    } else if (EMObject$Type=="EMRIDGE") {
      mtext(substitute(tau == a ,
       list(a = round(EMObject$tauD[1],3) ) ),
        side =3, cex=1.25, adj = .55, padj = padjB );
      mtext(paste("Mlt = ", round(
         EMObject$tauDK[2] /
         EMObject$tauAK[1],3), sep=""), cex=1.25, 
         adj = .86,padj = padjA )        
    }
   }
  mtext(expression(hat(B)), side=2,adj=3.5, cex=1.5,  padj=0,
       las=1)
	for (ii in 1:length(MyCols)) {
	    lines(x=xx, y= BBSeq[,ii], lty=ii, lwd=3, col=MyCols[ii]);
	}
  
     idd <-  floor(length(xx) * .5 / length(MyCols));
     idx <- floor((0:(length(MyCols)-1) * idd)) + floor( .4 * length(xx));
     ##if (length(xx) * .45 > length(MyCols) * 2) {
     ##  idx <- 2 *floor(0:(length(MyCols)-1)) * idd + floor( .4 * length(xx));
     ##}
     ##if (length(xx) * .45 > length(MyCols) * 3) {
     ##  idx <- 3 *floor(0:(length(MyCols)-1)) * idd + floor( .4 * length(xx));
     ##}     
     xxpts <- xx[idx];
       pps <- round(length(xx)/2);
     idxii = sort( BBSeq[pps,], index=TRUE, decreasin=TRUE)$ix;
     for (ii in 1:length(MyCols)) {
          yypts = BBSeq[idx[ii], idxii[ii]];
          xpt = xx[idx[ii]];
          points( x= xpt, y= yypts, pch=19, lwd=12, col="orange");
          ##points( x= (xpt+.001), y= yypts+.001, pch=ABCD[idxii[ii]], ps=8, col="green" );
          text( x= (xpt+.001), y= yypts+.001, 
            labels=substitute(bold(xx), list(xx=ABCD[idxii[ii]])), cex=1.1, col="black" );
     }
   if (symmet == TRUE) {
      yLims <- c( -(max(abs(BetaSeq))), max(abs(BetaSeq))  );
   } else {
      yLims <-  c(min(BetaSeq), max(BetaSeq))
   } 
   if (length(headerWANT) > 1 ||headerWANT != FALSE) {
     if (length(headerWANT) >=2)  {
       MyMainB = headerWANT[2];
     } else {
       MyMainB = headerWANT; 
     }
   }  else {
       MyMainB = "";
   }
   if (revScale == FALSE && logScale == TRUE) { 
     plot( x=c(min(xx), max(xx)), y=yLims,
           main=MyMainB, xlab="", ylab=, type="n", log="x");
   } else if (revScale == FALSE && logScale == FALSE) {
     plot( x=c(min(xx), max(xx)), y=yLims,
           main=MyMainB, xlab="", ylab="", type="n");    
   } else if (revScale == TRUE && logScale == FALSE) {
     plot( x=c(max(xx), min(xx)), y=yLims,
           main=MyMainB, xlab="", ylab="", type="n");    
   } else if (revScale == FALSE && logScale == TRUE) { 
     plot( x=c(max(xx), min(xx)), y=yLims,
           main=MyMainB, xlab="", ylab="",
            type="n", log="x");
   }
   mtext(expression(hat(beta)), side=2, adj=3.5, padj = 0, cex=1.5,
       las=1)
	##plot( x=c(min(xx), max(xx)), y=c(min(BetaSeq), max(BetaSeq)),
	##    main="", xlab="", ylab="", type="n", log="x");
	   lines(x=xx, y=2*SDSeqD, lwd=1, col="black", lty=2);
	   lines(x=xx, y=-2*SDSeqD, lwd=1, col="black", lty=2);
	   lines(x=xx, y=2*SDSeqA, lwd=1, col="black", lty=2);
	   lines(x=xx, y=-2*SDSeqA, lwd=1, col="black", lty=2);	   
	for (ii in 1:length(MyCols)) {
	    lines(x=xx, y= BetaSeq[,ii], lty=ii, lwd=3, col=MyCols[ii]);
	}
     idd <-  floor(length(xx) * .5 / length(MyCols));
     idx <- floor((0:(length(MyCols)-1) * idd)) + floor( .4 * length(xx));
     ##if (length(xx) * .45 > length(MyCols) * 2) {
     ##  idx <- 2 *floor(0:(length(MyCols)-1)) * idd + floor( .4 * length(xx));
     ##}
     ##if (length(xx) * .45 > length(MyCols) * 3) {
     ##  idx <- 3 *floor(0:(length(MyCols)-1)) * idd + floor( .4 * length(xx));
     ##}  
     if (length(MyCols) <= 26) {    
     xxpts <- xx[idx];
     idxii = sort( BetaSeq[pps,], index=TRUE, decreasin=TRUE)$ix;
     for (ii in 1:length(MyCols)) {
          yypts = BetaSeq[idx[ii], idxii[ii]];
          xpt = xx[idx[ii]];
          points( x= xpt, y= yypts, pch=20, lwd=12, col="orange");
          ##points( x= (xpt+.001), y= yypts+.001, pch=ABCD[idxii[ii]], ps=8, col="green" );
          text( x= (xpt+.001), y= yypts+.001, 
            labels=substitute(bold(xx), 
            list(xx=ABCD[idxii[ii]])), cex=1.1, col="black" );          
     }
     }

   return();
     
  
  
  }
 
 
 #######################################################################
 ##  SimCor : A quick, C++ based way to simulate vector random normals
 ##   with local correlation structure
 ##    Cor(X_i,X_j) = rho^(abs|i-j|); 
 SimCor <- function(lengthV = 2, rho = .3) {
     RTVec <- rnorm(lengthV, 0, 1);
     LTO <- .C("SimCor", lengthV = as.integer(lengthV), 
        rho= as.double(rho), RTVec = as.double(RTVec));
     return(LTO$RTVec);
 } 
 
 
 ######################################################################
##  GLMLogitLasso
##
##    Fits Logistic Regression using Selective LASSO penalty
##       CoordinateDescent is the algorithm of choice
##
##   y01 is a zero or 1 response vector
##   X is the matrix of covariates
##     Usually, a intercept is created for the X matrix
GLMLogitLasso <- function(X, y01, OnBeta = -999, Beta0 = 0,  OnGamma = 1,
   OnGammas = -1, RecOnBetas = FALSE,
   PrintFlag = 0, MaxCDOEpsilon = .00001,MaxCDOLoops = 100,
   InitKKs = 10, InverseGammaConstant =1)              {
   if (length(dim(X)) != 2) {
     print("GLMLogitLasso: X is not right dimension"); return(X);
   }
   if (length(OnBeta) != length(X[1,])) {
      OnBeta =  rep(0, length(X[1,]));
   }
   Beta0prev = 0;
   BetasPrev  = rep(0, length(OnBeta));  BetasPrev[1:length(OnBeta)]= OnBeta;
   if (length(OnGammas) <length(OnBeta)) {
       OnGammas = rep(OnGamma[1], length(OnBeta));
   }
   if (RecOnBetas == TRUE) {
      RecordOnBetas = matrix(0,p,MaxCDOLoops);
      Beta0Records = vector("numeric", MaxCDOLoops);
      RecordAllBetas = rbind(Beta0Records, RecordOnBetas);
   }  else {
      RecordOnBetas = -999; Beta0Records= -999;
      RecordAllBetas= -999;
   }
   ProbWeights = vector("numeric",n);
   ProbWeights = rep(1,n);
   logitCurrentProb = vector("numeric", n);
   logitCurrentProb = rep(0,n);
   logitPrevProb = rep(0,n);
   z = y01 *1.1;   Iter = 0;
   
   sendBetas = c(Beta0, OnBeta);
   prevBetas = c(Beta0prev, BetasPrev);
   AllGammas=c(0,OnGammas);
   XUse = cbind(rep(1, length(X[,1])), X);
   
   GLMLassoOb <- .Call("GLMLassoShell", XUse=XUse, y01=y01, 
     sendBetas=sendBetas, Beta0=Beta0, Beta0prev=Beta0prev, 
     prevBetas=prevBetas, OnGamma=OnGamma,
     RecordAllBetas=RecordAllBetas, Beta0Records=Beta0Records,
     AllGammas=AllGammas, InitKKs=InitKKs, 
     InverseGammaConstant=InverseGammaConstant,
     ProbWeights=ProbWeights, logitCurrentProb=logitCurrentProb, 
     logitPrevProb=logitPrevProb, z=z, 
     PrintFlag=PrintFlag, MaxCDOEpsilon=MaxCDOEpsilon, 
     MaxCDOLoops=MaxCDOLoops,Iter = Iter);
  
   RTOb <- list(ReturnBetas = sendBetas, XUse =XUse, y01=y01,
     OnGamma = OnGamma,
     RecordAllBetas=RecordAllBetas, 
     OnGammas = AllGammas,
     InitKKs = InitKKs,
     InverseGammaConstant = InverseGammaConstant,
     ProbWeights = ProbWeights,
     logitCurrentProb=logitCurrentProb, 
     logitPrevProb=logitPrevProb, z=z, 
     PrintFlag=PrintFlag, MaxCDOEpsilon=MaxCDOEpsilon, 
     MaxCDOLoops=MaxCDOLoops,Iter = Iter,
     prevBetas=prevBetas, Beta0=Beta0,
     Beta0Records = Beta0Records);
    
   return(RTOb);
}


######################################################################
## GLMG2Lasso:  Logit Regression using Two-Lasso
##
##     y01 is vector of zero-one responses, X is covariates
##    
##    piA is supposed pre activation, sigmaSq is a blurring factor
GLMG2Lasso <- function(X, y01, StartBeta=-999, StartBeta0=0, piA =.5, sigmaSq=2,
  LambdaAK = -999,LambdaDK = -999,OrderSeq =-999,
  StLambdaA = 1, StLambdaD = 1, MultLambdaA = .98,MultLambdaD=.98,
  TotalRuns = 100, MaxLambdaASteps = 50,
  RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
  MaxCDOEpsilon=.00001, 
  MaxCDOLoops=100, FixKa = -1, PriorPi = -1, InMaximumAllocation = -1)     {
  if (all(X[,1] == 1)) {
    X = X[,2:length(X[1,])];
  }
  p = length(X[1,]);
    
  if (length(StartBeta) != dim(X)[2]) {
    OnBeta = rep(0, dim(X)[2]);
  }
  BBHats= rep(piA,length(OnBeta));
  
  ProbWeights = rep(1,n);
  logitCurrentProb = vector("numeric", n);
  logitCurrentProb = rep(0,n);
  logitPrevProb = rep(0,n);
  z = y01 *1.1;   Iter = 0;
  Beta0 =StartBeta0;
  if (length(StartBeta) < p) {
    OnBeta = rep(0,p);
  }   else {OnBeta= StartBeta;}
   
  Beta0prev = rep(0, length(Beta0)); BetasPrev = rep(0, length(OnBeta));
  sendBetas = c(Beta0, OnBeta);
  prevBetas = c(Beta0prev, BetasPrev);
  BBHatsAll= rep(piA,p+1);
  
  if (LambdaAK[1] < 0) {
    LambdaAK = StLambdaA* MultLambdaA^(0:TotalRuns);
  } 
  if (LambdaDK[1] < 0) {
    LambdaDK = StLambdaD * MultLambdaD^(0:TotalRuns);
  }
  if (length(OrderSeq) < length(LambdaDK)) {
    OrderSeq = c(OrderSeq, 
      rep(OrderSeq[1], length(LambdaDK)- length(OrderSeq)));
  }
  if (length(LambdaAK) < length(LambdaDK)) {
    LambdaAK = c(LambdaAK, rep(LambdaAK[length(LambdaAK)], 
      length(LambdaDK)- length(LambdaAK)));
  }
  OnGammas = BBHats * LambdaAK[1] + (1-BBHats) * LambdaDK[1];
  AllGammas=c(0,OnGammas);
  XUse = cbind(rep(1, length(X[,1])), X);
  prevBetas = c(Beta0prev,BetasPrev);
  SigmaRecVec = rep(sigmaSq, length(LambdaDK));
  InverseGammaConstant = 1;
  
  PiRecVec = rep(piA, length(LambdaDK)); 
  if (RecordBetaFlag == TRUE) {
    RecordOnBetas = matrix(0,p, length(LambdaDK));
    RecordBBHats = matrix(0,p, length(LambdaDK));
    RecordBeta0 = vector("numeric", length(LambdaDK));
    RecordBBHatsAll = rbind(rep(1,length(LambdaDK)), RecordBBHats);
    RecordAllBetas = rbind(rep(1,length(LambdaDK)), RecordOnBetas);
  } else {
    RecordOnBetas = -999;  RecordBBHats=-999; RecordBeta0=-999;
    RecordAllBetas = -999;  RecordBBHatsAll= -999;
  }
  Iter = 0;  
  
  G2GLMO <- .Call("RunG2GLMLasso", X=XUse, y01=y01, 
    sendBetas=sendBetas, Beta0=Beta0,
    AllGammas=AllGammas, BBHatsAll=BBHatsAll, 
    RecordAllBetas = RecordAllBetas,
    Beta0prev=Beta0prev, prevBetas = prevBetas,
    RecordBeta0=RecordBeta0, RecordBBHatsAll=RecordBBHatsAll,
    LambdaDK=LambdaDK, LambdaAK=LambdaAK,
    piA=piA, PiRecVec = PiRecVec, sigmaSq=sigmaSq,
    SigmaRecVec=SigmaRecVec,
    OrderSeq=OrderSeq, ProbWeights=ProbWeights,
    logitCurrentProb=logitCurrentProb, 
    logitPrevProb=logitPrevProb, z=z, 
    InverseGammaConstant = InverseGammaConstant, InitKKs = InitKKs,  
    PrintFlag=PrintFlag, MaxCDOEpsilon=MaxCDOEpsilon, 
     MaxCDOLoops=MaxCDOLoops,Iter = Iter, FixKa = FixKa, sPriorS = PriorPi,
     InMaximumAllocation=InMaximumAllocation);
  ReturnBetasAll = sendBetas;  ReturnBetas = sendBetas;
  if (all(XUse[,1] == 1)) {
    Beta0 = ReturnBetas[1];
    ReturnBetas = ReturnBetas[2:length(ReturnBetas)];
    BBHatsAll = BBHatsAll[2:length(BBHatsAll)];
  }
  RetOb <- list(
    ReturnBetas=ReturnBetas,
    X=XUse, y01=y01, 
    sendBetas=sendBetas, Beta0=Beta0,
    AllGammas=AllGammas, BBHatsAll=BBHatsAll, 
    RecordAllBetas = RecordAllBetas,
    Beta0prev=Beta0prev, prevBetas = prevBetas,
    RecordBeta0=RecordBeta0, RecordBBHatsAll=RecordBBHatsAll,
    LambdaDK=LambdaDK, LambdaAK=LambdaAK,
    piA=piA, PiRecVec = PiRecVec, sigmaSq=sigmaSq,
    SigmaRecVec=SigmaRecVec,
    OrderSeq=OrderSeq, ProbWeights=ProbWeights,
    logitCurrentProb=logitCurrentProb, 
    logitPrevProb=logitPrevProb, z=z, 
    InverseGammaConstant = InverseGammaConstant, InitKKs = InitKKs,  
    PrintFlag=PrintFlag, MaxCDOEpsilon=MaxCDOEpsilon, 
    MaxCDOLoops=MaxCDOLoops,Iter = Iter, PriorPi = PriorPi,
    FailureFlag = FALSE);   
  if (FixKa >0) {RetOb$FixKa = FixKa;}
  return(RetOb);
}  
     
     
     
     
logZeroLambdaZ2 = (-20000:10000) / 1000;
minlZ2 = min(logZeroLambdaZ2)
DensLambdaZ2 = exp( -exp(2 * logZeroLambdaZ2) / 2);
SampleLambdaj0 <- function(NNeed) {
   lZ2 = sample(logZeroLambdaZ2, size=NNeed, prob = DensLambdaZ2, replace=TRUE);
   Z2 = exp(lZ2)
   Z2[lZ2== minlZ2]  =  0;
   Z1 =  sqrt(rchisq(NNeed,1));
   return(Z2/ Z1);
}

###################################################################3
### HorseShoe Implementation
###  An Implementation of Carvalho et al's horseshoe Gibbs Sampler
###     Issue is that it is slower and not a penalty
dlStandard = 1 / ( 1 +exp(2 * logZeroLambdaZ2));
FillerDens = dlStandard * 0;
nZeroLambdaZ2Sq = exp(-2 * logZeroLambdaZ2);
DrawAllLambdaj<-function(Betaj, tauSq, printFlag = 0) {
  Lambdaj = Betaj * 0;
  NZBeta = Betaj[Betaj != 0];
  if (length(NZBeta) > 0) {
    DrawLambdaNZj = NZBeta * 0;
    UnifDraw = runif(length(NZBeta));
    DrawLambdaNZj = .Call("OutGetAllLambdaj", logZeroLambdaZ2, dlStandard, NZBeta, tauSq,
       FillerDens, UnifDraw, DrawLambdaNZj, nZeroLambdaZ2Sq);
    Lambdaj[Betaj != 0] = DrawLambdaNZj;
  }
  if (length(Betaj[Betaj ==0]) > 0) {
     DrawLambdajZ = SampleLambdaj0(length(Betaj[Betaj==0]));
     Lambdaj[Betaj == 0] = DrawLambdajZ; 
  }
  return(Lambdaj);
}
#####################################################################
##  Draw Many Lambdaj is simply a test function designed
## To see if we are sampling from the posterior of Lamdaj | Betaj
##
##
DrawManyLambdaj <- function(Betaj, tauSq, NCount) {
   Lambdaj = rep(0, NCount);
   UnifDraws = runif(NCount);
   SUnifDraws = sort(UnifDraws, index=TRUE);
   Lambdaj = .Call("OutGetManyLambdaj", logZeroLambdaZ2, dlStandard, Betaj, tauSq,
       FillerDens, SUnifDraws$x, Lambdaj, nZeroLambdaZ2Sq);
   Lambdaj = Lambdaj[SUnifDraws$ix];
   return(Lambdaj);
}


SampleGibbsHorseShoe <- function (X, Y, tauSq, LengthSamples=250, SigmaSq, printFlag = 0){
   if (printFlag >= 1) {
     print("SampleGibbsHorseShoe: we are starting. "); flush.console();
   }
   RecordBetajs = matrix(0, length(X[1,]), LengthSamples);
   RecordLambdajs = matrix(0, length(X[1,]), LengthSamples);
   OnLambdaj = rep(1, length(X[1,]));
   XtX = t(X) %*% X;
   Q = XtX + diag(rep(SigmaSq/tauSq, length(X[1,])));
   LS = NULL;
   try(LS <- svd(Q));
   if (is.null(LS) || max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
      SQ <- pseudoinverse( Q);  
   }   else {
      SQ <- solve(Q);
   }   
   tXY= t(X) %*% Y;
   StartBetaj = SQ %*% tXY;
   OnBetaj = StartBetaj;
   if (printFlag >= 1) {
     print(paste("SampleGibbsHorseShoe: Ready to study samples for ",
       LengthSamples, sep="")); flush.console();
   }
   for (tjjin in 1:LengthSamples) {
     OnLambdaj = DrawAllLambdaj(OnBetaj, tauSq);
     OnBetaj  = SampleBetaj(X,Y, OnLambdaj, tauSq, tXY, XtX, SigmaSq);
     RecordBetajs[,tjjin] = OnBetaj;
     RecordLambdajs[,tjjin] = OnLambdaj;
     if (printFlag > 0 && tjjin %% ceiling(50/printFlag) == 0) {
       print(paste("SampleGibbsHorseShoe: and now tjjin=",
         tjjin, "/", LengthSamples, ".", sep=""));
       flush.console();
       par(mfrow=c(1,2));
          Ints <- 1:tjjin;
         plot(RecordBetajs[1,1:tjjin]~1:Ints, type="l", main="Beta 1");
         plot(RecordLambdajs[1,1:tjjin]~1:Ints, type="l", main="Lambda 1");         
     }
   }
   return(list(RecordBetajs = RecordBetajs, RecordLambdajs = RecordLambdajs));
}		                

SampleBetaj <- function(X,Y, Lambdaj, tauSq, tXY, XtX, sigmaSq) {
   if (is.null(Lambdaj)) {
     print("Error, Lambdaj is NULL");
   }
   Betaj = Lambdaj * 0;
   Betaj[Lambdaj == 0] = 0;
   if (sum(Lambdaj == 0) == length(Lambdaj)) {
     return(Betaj);
   }
   if (length(Lambdaj[Lambdaj != 0]) > 0) {
     nQ = XtX[Lambdaj != 0, Lambdaj != 0] + 
           diag( sigmaSq / (tauSq * Lambdaj[Lambdaj != 0]));
     if (length(nQ) == 1) {
       SnQ = 1/ nQ;
     }  else {
       LS = NULL;
       SnQ = NULL;
     
       try(SnQ <- solve(nQ), silent = TRUE);
       if (is.null(SnQ)) {
         try(SnQ <- pseudoinverse(nQ), silent=TRUE);
       }
       if (is.null(SnQ)) {
         dimnQ = dim(nQ);
         if (length(dim(nQ)) == 2 && (dim(nQ))[1] > 0 ) {
           SnQ = NULL;
           try(SnQ <- matrix(0, (dim(nQ))[1], (dim(nQ))[2]));
           if (is.null(SnQ)) {
             print("We can't help It, SnQ Sucks sucks Sucks!")
           }
           diag(SnQ) = 1 / diag(nQ);
         } else if (length(nQ) == 1)  {
           SnQ = matrix(1/nQ, 1,1);
         } else {
           print("Dimension of nQ is wrong.")
           print(paste("length of nQ == ", length(nQ), sep=""));
           SnQ = nQ * 0;
           return(Betaj);
         }
       }
     }
     ##try(LS <- svd(nQ));
     ##if (is.null(LS) || max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
     ##   SnQ =NULL;
     ##   try(SnQ <- pseudoinverse( nQ));
     ##   if (is.null(SnQ)) {
     ##     SnQ =matrix(0, length(nQ[1,]), length(nQ[,1]));
     ##     diag(SnQ) = 1 / diag(nQ);
     ##   }  
     ##}   else {
     ##   SnQ = solve(nQ);
     ##}
     MeanNBetaj = SnQ %*% tXY[Lambdaj != 0];
     sdNBetaj = NULL;
     try(sdNBetaj <-  t(chol(SnQ)));
     if (is.null(sdNBetaj)) {
        sdNBetaj = SnQ * 0;
        diag(sdNBetaj) = sqrt(diag(SnQ));
     }
     Betaj[Lambdaj != 0] = MeanNBetaj + sqrt(sigmaSq) * sdNBetaj  %*% rnorm( length(MeanNBetaj))
   }
   return(Betaj);
}


## Spike and Slab Gibbs Sampler
SampleGibbsSpikeSlab <- function (X, Y, tauASq, piA, LengthSamples, SigmaSq, printFlag = 0){
   RecordBetajs = matrix(0, length(X[1,]), LengthSamples);
   RecordBjs = matrix(0, length(X[1,]), LengthSamples);
   OnLambdaj = rep(1, length(X[1,]));
   XtX = t(X) %*% X;
   n = length(Y); p = length(X[1,]);
   Q = XtX + diag(rep(SigmaSq/tauASq, length(X[1,])));
   LS = NULL;
   try(LS <- svd(Q));
   SQ = NULL;
   if (is.null(LS) || max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
      try(SQ <- pseudoinverse( Q));  
   }   else {
      try(SQ <- solve(Q));
   }   
   tXY= t(X) %*% Y;
   StartBetaj = SQ %*% tXY;
   OnBetaj = StartBetaj;
   for (tjjin in 1:LengthSamples) {
     if (n < p ) {
        XtYResid = t(X) %*% (as.vector(Y) - as.vector(X %*% OnBetaj));
     } else {
        XtYResid = tXY - XtX %*% OnBetaj;
     }
     OnBj = SpikeAndSlabSampleBj(piA, tauAsq = tauASq, SigmaSq, XtYResid, OnBetaj, XtX)
     OnLambdaj = OnBj * tauASq;
     OnBetaj  = SampleBetaj(X,Y, OnLambdaj, 1, tXY, XtX, SigmaSq);
     RecordBetajs[,tjjin] = OnBetaj;
     RecordBjs[,tjjin] = OnBj;
     if (printFlag > 0) {
       par(mfrow=c(1,3));
          Ints <- 1:tjjin;
         plot(RecordBetajs[1,1:tjjin]~1:Ints, type="l", main="Beta 1");
         plot(RecordBjs[1,1:tjjin]~1:Ints, type="l", main="Bj 1");   
         if (tjjin > 1) {
           plot(colMeans(RecordBjs[,1:tjjin]), type="l", main=" mean Bj ");      
         }
     }
   }
   return(list(RecordBetajs = RecordBetajs, RecordBjs = RecordBjs));
}		          
                                               
SpikeAndSlabSampleBj <- function(piA, tauAsq, sigmaSq, XtYResid, Beta, XtX) {
    probNonZero = exp( (XtYResid + diag(XtX) * Beta)^2 / (diag(XtX) + sigmaSq/tauAsq ))*
                piA * sqrt( 2 * pi * sigmaSq) / sqrt( diag(XtX) + sigmaSq / tauAsq )    
    probZero = (1-piA)
    ##TotProb = probNonZero + probZero;
    probNonZero = probNonZero / (probNonZero + probZero)
    probNonZero[!is.finite(probNonZero)] = 1;
    return(rbinom(length(probNonZero), 1, probNonZero));
}
