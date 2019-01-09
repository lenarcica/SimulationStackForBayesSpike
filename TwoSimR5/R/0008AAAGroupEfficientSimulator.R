###############################################################################
##  0008AAAGroupEfficientSimulator.R
##    R scripts for Group methods used in TwoSimR5
##
##   (c) Alan Lenarcic 2009-2019
##     Written for simulation work related to methods paper for Valdar Lab, UNC Genetcis
##
##  These are the commands in TwoSimR52 that take the information in a "SMS"
##   object and then run the following algorithms to produce an estimator
##   based upon a group selection hypothesis.  This includes grpreg, grplasso, etc.
##
##   This is for "Grouped Effects" methods in the way "EfficientSimulator" is for Fixed effects.
##
##   There are many estimators in R for sparse group model selection, and
##   this code is intended to try and unify the input formats across a number
##   of estimators.  All of these estimators were current as of the last
##   time I had access to UNC Killdevil in Summer 2017.  However,
##   not all estimator results were used in Lenarcic and Valdar paper, mostly
##   because better or more canonical examples produced similar results.
##   Clearly, this covers only a tiny subset of available options, but
##   a purpose of this was to focus on Bayesian and Penalized estimators
##   which had a sparse setting, and to look specifically in Type 1 Type 2 
##   performance when the Beta versus noise setting was .25 up to 2.5.
##   All of these would be R packages with implementations sufficient for single
##   thread analyses, and no need for CUDA or multi-thread to compete. At the
##   time of the TwoSimR5 development (2010-2016), most methods fell into this category.
##
##   While it is hopeful that code like these simulations could help other
##    projects doing simulation performance tests in some other region of
##    the p,n,k,sigma,Beta,Cov(X), etc. choices every group likely has
##    a different region of interest and specialty, and TwoSimR5 seems
##    to demonstrate a role for all of these estimators, depending on the
##    location of noise space one belongs to, as well as the amount of prior
##    information one can come into the experiment holding.
##
##    Helper File for 221 project
##
##    This file corrects issues in TwoLasso package "SimForMe2" and represents
##      a fair testing situation for producing Stodden plots comparing TwoLasso
##      estimators to LARS estimates
##
##    For the most part you will only be using 
##      "GenerateLarsFixed" and "GenerateLimitLasso"
##      estimates and plotting to graphs.  However the function
##                       "SimulateCompare()"
##      gives a method for generating a table rating performance of
##      all of the LimitLasso estimators.
##    
##   SimMeData is a function at the end that generates data to rate LimitLasso
##     and TwoLasso against each other
##
##   Key inputs are a function "CorrelationXmatrix" which generates a 
##     Correlation matrix.  Choose one of CorrelationXmatrix1 or
##     CorrleationXmatrix2, and consider reasonable inputs.
##     DefaultCorrelationXmatrix is set yere to CorrelationXmatrix1
##     with an input of PosCorr = .3
##
##   Generating the Betas vector comes from a function declared
##     "GenerateBetaVec".  Choose one of the GenerateBetaVec functions
##     from GenerateBetaVec1, GenerateBetaVec2, GenerateBetaVec3
##     default is GenerateBetaVec2
##
##   An issue is that I've discovered some programming issue in underlying 
##     C code.  Though the main algorithms run as I intended them to,
##     I wrote my own helper version of the LARS algorithm in the function
##     LarsCC2() which is in desperate need of Lapack optimizing.  I'll be
##     working on an update of the underlying package C++ code for the next
##     few weeks.  So depending on how the preliminary results come out on this
##     assignment, the best runs will get into the next package update,
##     possibly getting in a publication.
##
##   


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

###############################################################################

##
##


## d/dx log(x) e^(-x) = - log(x) e^{-x}  + e^(-x)/x
## dv =  dx/x, u = e^(-x)
## v = log(x),  du = e^(-x) dx
##   log(x) e^(-x) - int log(x) e^(-x) dx
if (FALSE)  {
  install.packages("grpreg");
  install.packages("grplasso");
  install.packages("SGL");
  install.packages("standGL");
}

GenerateGroupTwoLassoSpread <- function(SMS, DoCI = DefaultDoCI) {
  library(TwoLassoCpp);
    USENOISE <- SMS$SigmaSqNoise; kLen <- SMS$kLen;
    ##AFilePrint("0008AAAGroupEfficientSimulator:  Running 2 Lasso Spread.")
	  T1 <- proc.time();
	      SDD <- sqrt(min(apply(SMS$XX,2,sd))^2 * (length(SMS$XX[,1])-1) );
	      LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD
                         , 
                       sqrt(3.7 / USENOISE) * SDD *  max( (2*length(SMS$XX[,1])), 200)
                           );
       ## LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD
       ##                  , 
       ##                sqrt(3.7 / USENOISE) * SDD *  200
       ##                    );                           
	      SABS <- min( abs(SMS$BetasReal[SMS$BetasReal != 0])  );
	      SDA <- SDD / sqrt(length(SMS$XX[,1])-1);
	      LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) * sqrt( 2 * SDA^2 / USENOISE)
        LambdaASeq = c( 1.1 * LambdaAS1, .05 * LambdaAS1);
        if (LambdaASeq[1] >= LambdaDSeq[2]) {
           LambdaASeq[1] = LambdaDSeq[1] + LambdaASeq[1];
           LambdaDSeq[1] = LambdaASeq[1] - LambdaDSeq[1];
           LambdaASeq[1] = LambdaASeq[1] - LambdaDSeq[1];
        }
    if (SMS$LogitRegression == TRUE) { DoLogit=TRUE; } else {
      DoLogit=FALSE;
    }
    LambdaASeq = c(LambdaASeq[1], LambdaASeq);
    LambdaDSeq = c(LambdaDSeq[1], LambdaDSeq);        
	OrderSeq = c(20, 4,1);
	puseGroups <- SMS$puseGroups;
        ## Do this to debug TwoLassoCpp to look for write error!
    puseGroups <- SMS$puseGroups;
	if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == FALSE) { 
     HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, DoLogit=FALSE, PiA = puseGroups, 
       SigmaSq = USENOISE * log(NCOL(SMS$XX)), LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = 1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = FALSE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior = c(USENOISE, 2),
       PiAPrior = c(puseGroups * sqrt(NCOL(SMS$XX)), (1-puseGroups) * sqrt(NCOL(SMS$XX))),
       IndexFirstRandomEffect = SMS$FirstGroupIndex, 
       tauEndList = SMS$EndGroupIndices);   
	   ## HitBatWingFast <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, 
	   ##   FixKa = -100, sigmaNoiseSq = USENOISE,
	   ##   RatWant = .2, StlambdaD=LambdaDSeq[1], 
	   ##   StlambdaA = LambdaASeq[1], LambdaDK = LambdaDSeq,
     ##   LambdaAK = LambdaASeq, OrderSeq = OrderSeq, 
     ##   TotalRuns = 5, NumEMConv = 4, MultEMCons = .99, NumCDOConv = 40, 
     ##   CDOEpsilon = .000001, TDFNu = SMS$tNoiseDF);
      if (HitBatWingFast$DidIFail >= 1) {
        AFilePrint("GenerateGroupTwoLasso Failure, we return TLS object. ");
        return(HitBatWingFast);
      }
    } else if (SMS$LogitRegression == TRUE) {
      ##LambdaAK = StlambdaA * lambdaAMultC^(0:4);
      ##LambdaDK = StlambdaD * lambdaDMultC^(0:4);
      ##OrderSeq = rep(10,5);
    if (FALSE) {
      HitBatWingFast = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
        StartBeta=-999, StartBeta0=0, piA = puseGroups, sigmaSq=2,
        LambdaAK = LambdaASeq,LambdaDK = LambdaDSeq,OrderSeq = OrderSeq,
        StLambdaA = LambdaASeq[1], StLambdaD = LambdaDSeq[1], 
        MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
        TotalRuns = 5, MaxLambdaASteps = 50,
        RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
        MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = -1,
        IndexFirstRandomEffect = SMS$FirstGroupIndex, 
        tauEndList = SMS$EndGroupIndices,verbose=3)  
    }
    if(SMS$n < SMS$p) {
      LogitNoise <- SMS$n / SMS$p;
    } else {
      LogitNoise <- 1.0;
       }
       HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, DoLogit=TRUE, PiA = SMS$puseGroups, 
        SigmaSq = USENOISE, LambdaAK = LambdaASeq,
        LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
        InitKKs = max(c(5,SMS$puse*length(SMS$YY))), RecordFlag = 1, HoldOn = FALSE, 
        CauchyEpsilon = .00001, MaxCauchy = 150,
        TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = puseGroups, 
        IndexFirstRandomEffect = SMS$FirstGroupIndex, 
        tauEndList = SMS$EndGroupIndices, DoLogit = TRUE, LogitNoise=LogitNoise);
    }
    AL <- length(HitBatWingFast$TLS$LambdaAK)-1;
    MIPRandom <-  (
      HitBatWingFast$TLS$LambdaDK[AL] - HitBatWingFast$TLS$GroupLambdaRecord[,AL])/
      ( HitBatWingFast$TLS$LambdaDK[AL] - HitBatWingFast$TLS$LambdaAK[AL]  );
    MIPRandom[MIPRandom >= 1.0] <- 1.0;
    MIPRandom[MIPRandom <= 0.0] <- HitBatWingFast$TLS$BackGroupsBBOn1[MIPRandom <= 0.0];

	  if (!is.null(HitBatWingFast$FailureFlag) && 
      length(HitBatWingFast$FailureFlag) >= 1 &&
      HitBatWingFast$FailureFlag == 1 ) {
	     AFilePrint("HitBatWingFast Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
	     FailToFitEMOb[[OnFails]] <<-HitBatWingFast;
	  }  	
    if (length(HitBatWingFast$FinalBeta) > 1) {
	  LarsFitBBsFast <- as.vector(HitBatWingFast$FinalBeta);
	  LarsFitBBsFast[abs(HitBatWingFast$FinalBeta) <= MinMin] = 0;
	  LarsFitBBsFast[abs(HitBatWingFast$FinalBeta) > MinMin] = 1;  
	  } else {
	    LarsFitBBsFast = SMS$BetasReal * 0 - 10;
	  }
	  T2 <- proc.time();
	     EMLARSFastT <- T2 - T1; 
    if (DoCI == TRUE) {
      try(rm(HitBatWingFast));
      CITime1 <- proc.time();
      if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) { 
       HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, DoLogit=DoLogit, PiA = SMS$puseGroups, 
       SigmaSq = USENOISE * log(NCOL(SMS$XX)), LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF,  SigmaPrior = c(USENOISE, 2),
       PiAPrior = c(puseGroups * sqrt(NCOL(SMS$XX)), (1-puseGroups) * sqrt(NCOL(SMS$XX))),
       IndexFirstRandomEffect = SMS$FirstGroupIndex, 
       tauEndList = SMS$EndGroupIndices); 
      }  
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      try(HitBatWingFast$TLS$ConfidenceQuantiles <- GetConfidenceIntervals);
      try(HitBatWingFast$TLS$LambdaIndexForConfidenceIntervals <- length(LambdaDSeq)-2);
      
      CIEst <- list();
      MIPMIP <-  HitBatWingFast$TLS$RecBBOn1[, length(LambdaDSeq)-1];
      try(MIPRandom <- HitBatWingFast$TLS$BackGroupsBBOn1);
      try(CIEst[[1]]<- HitBatWingFast$TLS$ConfidenceMatrix + 0.0);
      if (HitBatWingFast$DidIFail == 1) {
        AFilePrint("GenerateGroupTwoLasso: GetConfidenceIntervals fail!  Inspect!"); flush.console();
        return(HitBatWingFast);
      }
      try(CIEst[[2]]<- HitBatWingFast$TLS$UnshrunkConfidenceMatrix + 0.0);
      if (HitBatWingFast$DidIFail == 1) {
        AFilePrint("GenerateGroupTwoLasso: GetConfidenceIntervals fail!  Inspect!"); flush.console();
        return(HitBatWingFast);
      }
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix");
      CIEst <- CleanCIEst(CIEst);
      CIQuantiles <- HitBatWingFast$TLS$ConfidenceQuantiles + rep(0.0,
        length(HitBatWingFast$TLS$ConfidenceQuantiles) );
      CITime2 <- proc.time();  CITime <- CITime2-CITime1;
    } else {
      CIEst <- NULL;   CIQuantiles <- NULL;  CITime = NULL;  MIPMIP=NULL;
    }                                                                                                             
  return(list( type="GroupTwoLassoSpread", BetaFit = HitBatWingFast$FinalBeta,
    BBHatsFit = LarsFitBBsFast, FitTime =EMLARSFastT,
	OtherBetas = HitBatWingFast$TLS$Beta,
    EMLARSObject = HitBatWingFast,CIEst=CIEst, 
    CIQuantiles = CIQuantiles, CITime = CITime,
    MIP = MIPMIP, MIPReport = MIPMIP, MIPRandom=MIPRandom) );
}


###############################################################################
###  GenerateGroupBayesSpike()
###
###    A function for using BayesSpike Group inference
###
###    This is Group BayesSpike
###     
###
###    
###     
###    
###
###   
### 
GenerateGroupBayesSpike <- function(SMS, tauSqA = 1.0, CountSamp = 1000, 
  NumChains = 3, PriorStrength = 0, SigmaPriorStrength = 0,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo = "Auto",
  DoLogitPostPreProb = 1) {
  library(BayesSpike);
  kLen = SMS$kLen;
  ##AFilePrint("ENetFixed Run"); flush.console();
  T1 <- proc.time();
  if (!exists("NumChains")) { NumChains = 3; }
  if (!exists("CountSamp")) { CountSamp = 1000; }
  if (!exists("PriorStrength")) { PriorStrength = 0; }
  if (!exists("tauSqA")) { tauSqA = .5; }
  if (!exists("CutOff")) { CutOff = .1; }
  if (!exists("DoCI")) { DoCI = FALSE; }
  if (!exists("SigmaPriorStrength")) { SigmaPriorStrength <- 0; }
  if (!exists("DoLogitPostPreProb")) { DoLogitPostPreProb = 1; }
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE)  {
    DoLogitPreProb = FALSE; DoLogitPostProb = FALSE;
     dfRobit = -1; PreRunMaxGibbsIters = 0; 
      AlterWeightFlag = FALSE;
  } else {
  if (DoLogitPostPreProb == 2) {
    DoLogitPreProb = TRUE;  DoLogitPostProb = FALSE;
     dfRobit = -1; PreRunMaxGibbsIters = 25; 
     AlterWeightFlag = TRUE;
  } else {
    DoLogitPreProb = FALSE;  DoLogitPostProb = TRUE;
    dfRobit = -1; PreRunMaxGibbsIters = 25; 
    AlterWeightFlag = TRUE;
  }
  }
  
  if (tauSqA != 1.0) {
    tauPriordf = 1; tauPriorMean = tauSqA;
  } else if (DoMedian == FALSE) {
    tauPriordf = 2; tauPriorMean = 2; 
    ##PriorStrength = SMS$p;  SigmaPriorStrength = SMS$n;
  } else {
    tauPriordf = 1; tauPriorMean = 1;         
  }
  if (NCOL(SMS$X) >= 1100) {
    SigmaSq <- .5 * var(SMS$Y);
  } else {
    SigmaSq <- SMS$SigmaSqNoise;
  }
  
  if (is.logical(AutoInfo) && AutoInfo == TRUE) {
    AutoInfo = "Auto"
  } else if (is.logical(AutoInfo) && AutoInfo == FALSE) {
    AutoInfo = "Info"
  } else if (is.character(AutoInfo) && !(AutoInfo %in% c("Auto", "Info"))) {
    print(paste("Error AUTOINFO NOT SET!", sep=""));
  } else if (is.character(AutoInfo) && AutoInfo %in% c("Auto", "Info")) {
  } else {
    AutoInfo = "Auto";
  }

  nEndGroup = round(length(SMS$EndGroupIndices));
  if (AutoInfo == "Info") {
  if (SigmaPriorStrength <= 0) { SigmaPriorStrength = SMS$n; }  
  if (NCOL(SMS$X) >= 1100 && PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puseGroups + 2, PriorStrength*(1.0-SMS$puseGroups) +
      round(length(SMS$EndGroupIndices)));
    if (SigmaPriorStrength <= 0) { SigmaPriorStrength = SMS$n; }
    SigmaPrior <- c(SigmaPriorStrength, SMS$SigmaSqNoise);
  } else if (NCOL(SMS$X) >= 1000 && PriorStrength <= 0) {
    PiAPrior = c( max(1,SMS$puseGroups* SMS$n),
      max(1,SMS$n * (1-SMS$puseGroups))); 
    SigmaPrior <- c(SMS$n,SMS$SigmaSqNoise);   
  } else if (PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puseGroups, PriorStrength * (1.0-SMS$puseGroups));
    SigmaPrior = c(SigmaPriorStrength, SMS$SigmaSqNoise);
  } else {
    SigmaPrior = c(1,1); PiAPrior = c(-1,-1);
  }
  } else {
  if (SigmaPriorStrength <= 0) { SigmaPriorStrength = SMS$n; }
  if (NCOL(SMS$X) >= 1100 && PriorStrength > 0) {
    PiAPrior = c(PriorStrength * 1/SMS$kLen + 1, PriorStrength * round(length(SMS$EndGroupIndices)) + 1
      );
    SigmaPrior <- c(SigmaPriorStrength, SMS$SigmaSqNoise);
  } else if (NCOL(SMS$X) >= 1000 && PriorStrength <= 0) {
    PiAPrior = c(1/nEndGroup + 1,
        nEndGroup+1); 
    SigmaPrior <- c(SMS$n,SMS$SigmaSqNoise);   
  } else if (PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puseGroups, PriorStrength * (1.0-SMS$puseGroups));
    SigmaPrior = c(SigmaPriorStrength, SMS$SigmaSqNoise);
  } else {
    SigmaPrior = c(1,1); PiAPrior = c(-1,-1);
  }    
    
  }
  eval(parse(text=GetG0Text("BSSaveDir")));
  if (is.character(BSSaveDir) && (nchar(BSSaveDir) < nchar(UniqueProcessIdentifier) ||
    substr(BSSaveDir, nchar(BSSaveDir) - nchar(UniqueProcessIdentifier)+1,
    nchar(BSSaveDir)) != UniqueProcessIdentifier)) {
    BSSaveDir <- paste(BSSaveDir, "//", UniqueProcessIdentifier, sep="");      
  }

  if ((kLen <= 250 || is.null(BSSaveDir) || is.numeric(BSSaveDir) || !is.character(BSSaveDir)) && DoLogitPostProb == FALSE ) {
    GibbsSS <- BayesSpike:::BayesSpikeRegression(X=SMS$X, Y=SMS$Y, 
      tauASq = tauSqA, PiAStart= PiAPrior[1] / sum(PiAPrior), MaxGibbsIters = CountSamp,
      NumChains=NumChains,     tauPriordf = tauPriordf, tauPriorMean = tauPriorMean, 
      SigmaSq = SMS$SigmaSqNoise, Verbose = 0, PiAPrior = PiAPrior,
      SigmaPrior = SigmaPrior, DoSave=FALSE,
      DoRecord=c(1,0,0,0,0,0,1), IndexFirstRandomEffect = SMS$FirstGroupIndex,
      tauEndList = SMS$EndGroupIndices,
      DoLogitPreProb = DoLogitPreProb, DoLogitPostProb = DoLogitPostProb,
      dfRobit = dfRobit, PreRunMaxGibbsIters = PreRunMaxGibbsIters, 
      AlterWeightFlag = AlterWeightFlag);
    

    ## MIPTAU = GibbsSS$RunProbVector / GibbsSS$TotEveryProbVector;
     if (DoMedian == TRUE) {
     	HL <- GibbsSS$CodaList
     	NHL <- list();
	  for (ii in 1:length(HL)) {
        NHL[[ii]] <- as.mcmc(HL[[ii]][round(NROW(HL[[ii]])/4):NROW(HL[[ii]]),]);
      }
     try(NHL <- as.mcmc.list(NHL));
     suppressWarnings(MySum <- summary(NHL));
     MakeAns <- MySum[[2]][1:GibbsSS$p,3];
     KeepIters <- (1:NCOL(SMS$X))[MySum[[2]][1:GibbsSS$p,3] != 0.0];
     HSAns <- MySum[[2]][1:GibbsSS$p,3]; 
     MakeAns <- MySum[[1]][1:GibbsSS$p,1];
     MIP <- GibbsSS$MIP;
     MIPRandom <- MIP;

     } else {
       MIP <- GibbsSS$ProbVector;
       MIPRandom <- MIP;
       MakeAns <- GibbsSS$QMedianBeta;
       HSAns <- MakeAns;
       KeepIters <- (1:GibbsSS$p)[HSAns != 0.0];
     }
     T2 <- proc.time();

     if (DoMedian == FALSE) {
       HL <- GibbsSS$CodaList;
       for (ii in 1:length(HL)) {
        NHL[[ii]] <- as.mcmc(HL[[ii]][round(NROW(HL[[ii]])/4):NROW(HL[[ii]]),]);
      }
      try(NHL <- as.mcmc.list(NHL));
      suppressWarnings(MySum <- summary(NHL));
      MakeAns <- MySum[[2]][1:GibbsSS$p,3];  
     }     
     Keep <- rep(0, NCOL(SMS$X));
     if (GibbsSS$FirstRandom > 1) {
       Keep[1:(GibbsSS$FirstRandom-1)] <- MIP[1:(GibbsSS$FirstRandom-1)]
     }
     AK <- GibbsSS$FirstRandom-1; AST = GibbsSS$FirstRandom-1;
     if (GibbsSS$FirstRandom <= 0) { AK = 0;  AST = 0; }
     for (ii in 1:length(GibbsSS$OnTau)) {
       Keep[AK:(GibbsSS$tauEndList[ii]+1)] = MIP[AST+ii];
       AK <- GibbsSS$tauEndList[ii]+1;
     }
     MIPMIP <- Keep;
  } else {
    if (FALSE) {
      GibbsSS <- BayesSpike:::BayesSpikeRegression(X=SMS$X, Y=SMS$Y, 
      tauASq = tauSqA, PiAStart=SMS$puse, MaxGibbsIters = 200,
      NumChains=NumChains,     tauPriordf = tauPriordf, tauPriorMean = tauPriorMean, 
      SigmaSq = as.numeric(SigmaSq), Verbose = 2, PiAPrior = PiAPrior,
      SigmaPrior = SigmaPrior, DoSave=TRUE, DoRecord = c(0,0,0,0,0,0,0),
      SaveDir = BSSaveDir, IndexFirstRandomEffect = SMS$FirstGroupIndex,
      tauEndList = SMS$EndGroupIndices, ZeroOutBeforeSample = TRUE,
      StartRunProbVector = 10,
      DoLogitPreProb = DoLogitPreProb, DoLogitPostProb = DoLogitPostProb,
      dfRobit = dfRobit, PreRunMaxGibbsIters = PreRunMaxGibbsIters, 
      AlterWeightFlag = AlterWeightFlag); 
    }
    if (NCOL(SMS$X) >= 10000) {
      MaxGibbsIters <- 200;
    } else {
      MaxGibbsIters <- 500;
    }
    GibbsSS <- try(BayesSpike:::BayesSpikeRegression(X=SMS$X, Y=SMS$Y, 
      tauASq = tauSqA, PiAStart=PiAPrior[1]/sum(PiAPrior), MaxGibbsIters = MaxGibbsIters,
      NumChains=NumChains,
      SigmaSq = as.numeric(SigmaSq), Verbose = 0, PiAPrior = PiAPrior,
      tauPriordf = tauPriordf, tauPriorMean = tauPriorMean,
      SigmaPrior = SigmaPrior, DoSave=TRUE, DoRecord = c(0,0,0,0,0,0,0),
      SaveDir = BSSaveDir, IndexFirstRandomEffect = SMS$FirstGroupIndex,
      tauEndList = SMS$EndGroupIndices, ZeroOutBeforeSample = TRUE,
      StartRunProbVector = 10,
      DoLogitPreProb = DoLogitPreProb, DoLogitPostProb = DoLogitPostProb,
      dfRobit = dfRobit, PreRunMaxGibbsIters = PreRunMaxGibbsIters, 
      AlterWeightFlag = AlterWeightFlag, burnin = 25)); 
     
     MIP <- GibbsSS$ProbVector;  HL <- NULL;
     MIPRandom <- MIP;  
     Keep <- rep(0, NCOL(SMS$X));
     if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE && GibbsSS$FirstRandom > 1) {
       HL <- NULL;
        MIPRandom <- MIP[2:length(MIP)];
         Keep <- rep(0, NCOL(SMS$X));
        if (GibbsSS$FirstRandom > 1) {
          Keep[1:(GibbsSS$FirstRandom-1)] <- MIP[1:(GibbsSS$FirstRandom-1)]
        }
        AK <- GibbsSS$FirstRandom-1; AST = GibbsSS$FirstRandom-1;
          for (ii in 1:length(GibbsSS$OnTau)) {
              if (MIPRandom[ii] > .5) {
                ##print(paste("On ii = ", ii, " we have tauEndList[ii]+1=",
                ##  GibbsSS$tauEndList[ii]+1)); flush.console();
                ##break;
              }
              Keep[(AK+1):(GibbsSS$tauEndList[ii]+1)] <- MIPRandom[ii];
              AK <-GibbsSS$tauEndList[ii]+1;
          }
     }  else if (GibbsSS$FirstRandom > 1) {
       HL <- NULL;
       try(GibbsSS$burnin <- 25);      
       Keep[1:(GibbsSS$FirstRandom-1)] <- MIP[1:(GibbsSS$FirstRandom-1)]
            AK <- GibbsSS$FirstRandom-1; AST = GibbsSS$FirstRandom-1;
       if (GibbsSS$FirstRandom <= 1) { AK = 0;  AST = 0; }  else {
         AST <- GibbsSS$FirstRandom-1
       }
       for (ii in 1:length(GibbsSS$OnTau)) {
         Keep[(AK + 1):(GibbsSS$tauEndList[ii]+1)] = MIPRandom[AST+ii];
         AK <- GibbsSS$tauEndList[ii]+1;
       }
     } else {
      HL <- NULL;
      try(GibbsSS$burnin <- 25);
      AK <- GibbsSS$FirstRandom-1; AST <- 0;
      for (ii in 1:length(GibbsSS$OnTau)) {
         Keep[(AK + 1):(GibbsSS$tauEndList[ii]+1)] = MIPRandom[AST+ii];
         AK <- GibbsSS$tauEndList[ii]+1;
      }
     }


     

     MIPMIP <- Keep;
     if (TRUE || DoMedian == TRUE) {
     SubSetCoords <- (1:length(Keep))[Keep >= .4]; 
     if(length(SubSetCoords) <= 0) {
       SubSetCoords <- sort(sort(Keep, decreasing=TRUE, 
         index=TRUE)$ix[1:round(SMS$kActiveLen* (SMS$EndGroupIndices[1]-SMS$FirstGroupIndex))] )
     }
     GibbsSS$SubSetCoords <- SubSetCoords;

     if (length(SubSetCoords) >= 1 && DoLogitPostProb == TRUE) {
       eval(parse(text=SetGText("GibbsSS", "globalenv()", S=1)));
       LSCL <- GibbsSS$TBSR5$LoadSubCodaList(StartIter = 1, EndIter = GibbsSS$TBSR5$MaxGibbsIters);
       ALongBetaMatrix <- AStackCoda(GibbsSS$SubCodaList, GibbsSS$burnin);
       eval(parse(text=SetGText("ALongBetaMatrix", "globalenv()", S=1)));
       print(paste("Dim ALongBetaMatrix is (", paste(dim(ALongBetaMatrix), collapse=", ", sep=""), ")", sep=""));
       flush.console();
        try(WeightCodas <- GibbsSS$AlterWeightCodaList)
        try(StackWeightCodas <- AStackCoda(WeightCodas, GibbsSS$burnin));
        if (sd(StackWeightCodas) >= 1.5) {
          StackWeightCodas[StackWeightCodas > quantile(StackWeightCodas,.99)] <- quantile(StackWeightCodas, .99);
         StackWeightCodas <- StackWeightCodas / sd(StackWeightCodas) * 1.5;
        }

        ReWeight <- exp(StackWeightCodas-max(StackWeightCodas));
        ReWeight[ReWeight > (quantile(ReWeight,.75) + 
          1.5 * (quantile(ReWeight,.75) - quantile(ReWeight,.25)))]  <- 
          (quantile(ReWeight,.75) + 
          1.5 * (quantile(ReWeight,.75) - quantile(ReWeight,.25)))
        ReWeight <- ReWeight/ sum(ReWeight);
      if (is.null(ReWeight)) {
        print("Note ReWeight is null, not good on alternate weight system.");
        flush.console();
      }  else {
        print(paste("Note: ReWeight has a dimension ", length(ReWeight), " in case that's good for you. ", sep=""));
        eval(parse(text=SetGText("ReWeight", S=1)));
        
      }
        MyQuant <- QuantileMe(ReWeight, ALongBetaMatrix, QuantileLocs = .5, Verbose = 0);
       HSAns <- rep(0, GibbsSS$p);
       HSAns[SubSetCoords] <- MyQuant;
       MakeAns <- rep(0, GibbsSS$p);
       if (length(SubSetCoords) == 1) {
         MakeAns[SubSetCoords] <- sum(as.vector(ReWeight) * ALongBetaMatrix);
       } else {
         MakeAns[SubSetCoords] <- colSums(as.vector(ReWeight) * ALongBetaMatrix); 
       }
     }  else {
     ACoda2 <- GibbsSS$SubCodaList;
     MySumP <- summary(ACoda2);
     KeepIters <- SubSetCoords[MySumP[[2]][,3] != 0.0];
     HSAns <- rep(0, GibbsSS$p);
     HSAns[SubSetCoords] <- MySumP[[2]][,3];
     MakeAns <- rep(0, GibbsSS$p);
     MakeAns[SubSetCoords] <- MySumP[[1]][,1];
     }
     } else {
       HSAns <- GibbsSS$QMedianBeta;
       MakeAns <- HSAns;
     }
     T2 <- proc.time();
     if (FALSE && DoMedian == FALSE) {
      SubSetCoords <- (1:length(Keep))[Keep >= .5]; 
      if(length(SubSetCoords) <= 0) {
        SubSetCoords <- sort(Keep, decreasing=TRUE, 
          index=TRUE)$ix[1:round(SMS$kActiveLen* (SMS$EndGroupIndices[1]-SMS$FirstGroupIndex))]
      }
      GibbsSS$SubSetCoords <- SubSetCoords;
      ACoda2 <- GibbsSS$SubCodaList;
      MySumP <- summary(ACoda2);
      KeepIters <- SubSetCoords[MySumP[[2]][,3] != 0.0];           
     }
  }
   
    GibbsSS$TBSR5$unsave(TRUE);

  ActiveHS = HSAns;
  ActiveHS[ActiveHS != 0] = 1;
  if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
    if ((SMS$p +1) == length(ActiveHS)) {
      try(ActiveHS <- ActiveHS[2:length(ActiveHS)]);
      try(HSAns <- HSAns[2:length(HSAns)]);
      try(MakeAns <- MakeAns[2:length(MakeAns)]); 
      if (length(MIPMIP) == (SMS$p +1)) {
        try(MIPMIP <- MIPMIP[2:length(MIPMIP)]);
      }
    }
  }
	HSFitTime <- T2 - T1;
  if (DoCI == TRUE) {
      if (NCOL(SMS$X) >= 10000) {
        CountSamp <- 1000;
      } else {
        CountSamp <- 1000;
      }
	  eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
         SubAns <- c(0, MakeAns);
      } else { SubAns <- MakeAns; }
	  CITime1 <- proc.time();
	  GibbsSS <- BayesSpike:::BayesSpikeRegression(X=SMS$X, Y=SMS$Y, 
      BetaStart = SubAns, NoNoiseBetaStart=TRUE,
      tauASq = tauSqA, PiAStart=PiAPrior[1]/sum(PiAPrior), MaxGibbsIters = CountSamp,
      NumChains=NumChains,
      SigmaSq = as.numeric(SigmaSq), Verbose = 2, PiAPrior = PiAPrior,
      tauPriordf = tauPriordf, tauPriorMean=tauPriorMean,
      SigmaPrior = SigmaPrior, DoSave=TRUE, DoRecord = c(0,0,0,0,0,0,0), 
      ZeroOutBeforeSample = TRUE,
      StartRunProbVector = 10,
      SaveDir = BSSaveDir, IndexFirstRandomEffect = SMS$FirstGroupIndex,
      tauEndList = SMS$EndGroupIndices, DoLongCI = TRUE,
      DoLogitPreProb = DoLogitPreProb, DoLogitPostProb = DoLogitPostProb,
      dfRobit = dfRobit, PreRunMaxGibbsIters = PreRunMaxGibbsIters, 
      AlterWeightFlag = AlterWeightFlag, burnin = 25); 
      
      try(GibbsSS$burnin <- 25);
      GibbsSS$TBSR5$KeepPosteriorQuantiles <- GetConfidenceIntervals;
	  CIEst <- list();
	   try(CIEst[[1]] <- GibbsSS$TBSR5$BetaSymmetricQuantiles);
	   try(CIEst[[2]] <- GibbsSS$TBSR5$BetaSymmetricUnshrinkQuantiles);
	   try(CIEst[[3]] <- GibbsSS$TBSR5$BetaHPDQuantiles);
	   try(CIEst[[4]] <- GibbsSS$TBSR5$BetaHPDUnshrinkQuantiles);   
 	   try(CIEst[[5]] <- GibbsSS$TBSR5$BetaSymmetricLongQuantiles);
	   try(CIEst[[6]] <- GibbsSS$TBSR5$BetaHPDLongQuantiles);   	   
    CITime2 <- proc.time();
    CIEst <- CleanCIEst(CIEst);
      try(names(CIEst) <- c("BetaSymmetricQuantiles", "BetaSymmetricUnshrinkQuantiles",
        "BetaHPDQuantiles", "BetaHPDUnshrinkQuantiles", 
        "BetaSymmetricLongQuantiles", "BetaHPDLongQuantiles"));
      CIQuantiles <- GibbsSS$TBSR5$KeepPosteriorQuantiles;
    CITime <- CITime2-CITime1;
    if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
      for (ii in 1:length(CIEst)) {
        try(CIEst[[ii]] <- CIEst[[ii]][2:NROW(CIEst[[ii]]),]);
        try(rownames(CIEst[[ii]]) <- paste("Beta:", 1:NROW(CIEst[[ii]]), sep=""));
      }
    }
    GibbsSS$TBSR5$unsave(TRUE);
   } else {
     CIEst <- NULL;   CIQuantiles <- NULL;  CITime <- NULL;
   }  
	##AFilePrint("Finished RunningENetFixed"); flush.console();
	return(list( type="GroupBayesSpike", BetaFit = HSAns, BBHatsFit = ActiveHS,
	  FitTime =HSFitTime, OtherBetas = MakeAns, HitOb = HL,
    CIEst=CIEst, CIQuantiles=CIQuantiles,CITime=CITime, MIP = MIPMIP,
    MIPRandom = MIPRandom, MIPReport = MIPMIP, DoMedian = DoMedian) );                  
}

###############################################################################
###  GenerateGroupGrpReg()
###
###    A function for using Lasso GroupGrg inference
###
###    This is Group Lasso
###     
###
###    
###     
###    
###
###   
### 
GenerateGroupLassoReg <- function(SMS, Interval = Interval,
  Criterion = "BIC") {
  MyT <- "ATT <- FALSE; 
    library(grplasso, warn.conflicts=FALSE, quietly=TRUE);  
    ATT <- TRUE;"
  try(eval(parse(text=MyT)));
  if (ATT == FALSE) {
    AFilePrint("GenerateGroupLassoReg: FAIL no grplasso package! installed.");
     1 = 2;
  }
	T1 <- proc.time();
  if (is.character(Criterion)) {
    if (Criterion[1] %in% c("1", "BIC", "bic", "Bic")) {
      Criterion = "BIC"
    } else if (Criterion[1] %in% c("2", "AIC", "Aic", "aIC")) {
      Criterion = "AIC"
    } else if (Criterion[1] %in% c("GCV", "3", "gcv", "GCv", "Gcv")) {
      Criterion = "GCV";
    } else {
      AFilePrint(paste("GenerateGroupLassoReg: Hey You GAVE A VERY BAD character Criterion = ",
        Criterion, sep="")); flush.console();
      One = Two;
    }   
  } else if (is.numeric(Criterion)) {
    if (Criterion[1] == 1) {
      Criterion = "BIC"
    } else if (Criterion[1] == 2) {
      Criterion = "AIC"
    } else if (Criterion[1] == 3) {
      Criterion = "GCV";
    } else {
      AFilePrint(paste("GenerateGroupLassoReg: Hey you GAVE A VERY BAD, numeric Criterion = ",
        Criterion, ". ", sep="")); flush.console();
      One = Two;
    }
  } 
	

	## grpreg         

  if (!exists("LogitRegression")) {
    LogitRegression <- SMS$LogitRegression;
  }
  if (LogitRegression == FALSE) {
	  lambda <- lambdamax(SMS$XX, y = SMS$YY, index = SMS$index, penscale = sqrt,
                    model = LinReg(), center=FALSE) * 0.5^(0:5)
    fit <- grplasso(x=SMS$XX, y=SMS$YY, index=SMS$index, 
         weights = rep(1, length(SMS$YY)), offset = rep(0,
         length(SMS$YY)), lambda, coef.init = rep(0, ncol(SMS$XX)),
         penscale = sqrt, model = LinReg(), center = FALSE,
         standardize = TRUE, control = grpl.control());
    SS <- colSums((fit$fitted-
      matrix(rep(SMS$YY, NCOL(fit$fitted)), 
      length(SMS$YY), NCOL(fit$fitted)))^2)
  } else {
    lambda <- lambdamax(SMS$XX, y = SMS$YY, index = SMS$index, penscale = sqrt,
                    model = LinReg(), center=FALSE) * 0.5^(0:5)
    fit <- grplasso(x=SMS$XX, y=SMS$YY, index=SMS$index, 
         weights = rep(1, length(SMS$YY)), offset = rep(0,
         length(SMS$YY)), lambda, coef.init = rep(0, ncol(SMS$XX)),
         penscale = sqrt, model = LogReg(), center = FALSE,
         standardize = TRUE, control = grpl.control()); 
    SS <- colSums(log(fit$fitted[SMS$YY>= .5]))+ colSums(log(1-fit$fitted[SMS$YY <= .5]));
  }
  AFitT <-fit$coefficients;
  AFitT[fit$coefficients != 0.0] <- 1;
  kks <- colSums(AFitT);
  BICC <- SS+ kks * log(length(SMS$YY));
  Winner <- sort.int(BICC, index.return=TRUE)$ix[1]
 
  T2 <- proc.time();
  GPLFit <- fit$coefficients[,Winner];
  GPLFitBBs <- rep(0, length(GPLFit));
  GPLFitBBs[GPLFit != 0.0] = 1.0;
  
  MIPRandom <- sort(unique(SMS$index));
  for (kk in 1:length(MIPRandom)) {
    MIPRandom[kk] <- mean(GPLFitBBs[SMS$index==kk])
  }
  
  if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE &&
    NCOL(SMS$XX)+1 == length(GPLFit)) {
    GPLFit <- GPLFit[2:length(GPLFit)];    
    GPLFitBBs <- GPLFitBBs[2:length(GPLFitBBs)];
  }


	GPLFitTime <- T2 - T1;
	return(list( type="GroupLassoReg", BetaFit = GPLFit,
    BBHatsFit = GPLFitBBs, FitTime = GPLFitTime,
	  OtherBetas = GPLFit, GroupRegFitObject = fit, CIEst = NULL, CITime = NULL,
    CI = NULL, MIPReport = GPLFitBBs, MIPRandom=MIPRandom) );        		     
}

###############################################################################
###  GenerateGroupGrpReg()
###
###    A function for using GroupGrp  GroupLasso inference
###
###    Function "grpreg" has numerous stopping criterions, need to run them all.
###     
###
###    
###     
###    
###
GenerateGroupGrpReg <- function(SMS, Interval = Interval,
  Criterion = "BIC") {
    library(grpreg);
	T1 <- proc.time();
  if (is.character(Criterion)) {
    if (Criterion[1] %in% c("1", "BIC", "bic", "Bic")) {
      Criterion = "BIC"
    } else if (Criterion[1] %in% c("2", "AIC", "Aic", "aIC")) {
      Criterion = "AIC"
    } else if (Criterion[1] %in% c("GCV", "3", "gcv", "GCv", "Gcv")) {
      Criterion = "GCV";
    } else {
      AFilePrint(paste("GenerateGroupGrpReg: Hey You GAVE A VERY BAD character Criterion = ",
        Criterion, sep="")); flush.console();
      One = Two;
    }   
  } else if (is.numeric(Criterion)) {
    if (Criterion[1] == 1) {
      Criterion = "BIC"
    } else if (Criterion[1] == 2) {
      Criterion = "AIC"
    } else if (Criterion[1] == 3) {
      Criterion = "GCV";
    } else {
      AFilePrint(paste("GenerateGroupGrpReg: Hey you GAVE A VERY BAD, numeric Criterion = ",
        Criterion, ". ", sep="")); flush.console();
      One = Two;
    }
  }
	
	## grpreg         
  MyT <- "ATT <- FALSE; 
    library(grpreg, warn.conflicts=FALSE, quietly=TRUE);  
    ATT <- TRUE;"
  try(eval(parse(text=MyT)));
  if (ATT == FALSE) {
    AFilePrint("GenerateGroupGrpReg: FAIL no grpreg package! installed.");
     1 = 2;
  }
  if (!exists("LogitRegression")) {
    LogitRegression <- SMS$LogitRegression;
  }
  if (LogitRegression == FALSE) {
    fit <- grpreg(SMS$XX,SMS$YY,index=SMS$index,penalty="grLasso");
  } else {
    fit <- grpreg(SMS$XX,SMS$YY,index=SMS$index,penalty="grLasso"); 
  }
  
  AS <- select(fit, Criterion);
  T2 <- proc.time();
  GPLFit <- AS$beta;
  GPLFitBBs <- rep(0, length(GPLFit));
  GPLFitBBs[GPLFit != 0.0] = 1.0;

  MIPRandom <- sort(unique(SMS$index));
  for (kk in 1:length(MIPRandom)) {
    MIPRandom[kk] <- mean(GPLFitBBs[SMS$index==kk])
  }

  if (!is.null(SMS$LogitRegression) && SMS$LogitRegression==TRUE) {
     GPLFit <- GPLFit[2:length(GPLFit)];
     GPLFitBBs <- GPLFitBBs[2:length(GPLFitBBs)];
  }
	GPLFitTime <- T2 - T1;
	return(list( type="GroupRegLasso", BetaFit = GPLFit,
    BBHatsFit = GPLFitBBs, FitTime = GPLFitTime,
	  OtherBetas = GPLFit, GroupRegFitObject = fit, CIEst = NULL, CITime = NULL,
    MIPReport = GPLFitBBs, MIPRandom=MIPRandom) );        		     
}

###############################################################################
###  GenerateStandGLReg()
###
###    A function for using StandGL  GroupLasso inference
###
###    Function "standGL" bases performance off of crossvalidation routine.
###     
###
###    
###     
###    
###
GenerateStandGLReg <- function(SMS, Interval = Interval) {
	T1 <- proc.time();
	
	## grpreg         
  MyT <- "ATT <- FALSE; 
    library(standGL, warn.conflicts=FALSE, quietly=TRUE);  
    ATT <- TRUE;"
  try(eval(parse(text=MyT)));
  if (ATT == FALSE) {
    AFilePrint("GenerateGroupGrpReg: FAIL no grpreg package! installed.");
     1 = 2;
  }
  
  if (!exists("LogitRegression")) {
    LogitRegression <- SMS$LogitRegression;
  }
  if (LogitRegression == FALSE) {
    cv.fit <- cv.standGL(X=SMS$XX,y=SMS$YY,index=SMS$index,family="linear");
  } else {
    cv.fit <- cv.standGL(X=SMS$XX,y=SMS$YY,index=SMS$index,family="logit"); 
  }
  
  lambdaFit <- cv.fit$lambda.min;
  BetaWant <- cv.fit$fit$beta[, cv.fit$fit$lam.path == lambdaFit];
  T2 <- proc.time();
  StandGLFit <- BetaWant;
  StandGLFitBBs <- rep(0, length(BetaWant));
  StandGLFitBBs[BetaWant != 0.0] = 1.0;

  MIPRandom <- sort(unique(SMS$index));
  for (kk in 1:length(MIPRandom)) {
    MIPRandom[kk] <- mean(StandGLFitBBs[SMS$index==kk])
  }

	StandGLFitTime <- T2 - T1;
	return(list( type="StandGLCVLasso", BetaFit = StandGLFit,
    BBHatsFit = StandGLFitBBs, FitTime = StandGLFitTime,
	  OtherBetas = StandGLFit, StandGLFitObject = cv.fit, CIEst = NULL,
    CITime = NULL, MIPReport = StandGLFitBBs, MIPRandom = MIPRandom) );        		     
}


#######################################################################3
##  GenerateSGLReg:
##
##    Tests SGL-package Regularized Group Lasso
##
##    This algorithm uses cross validation.
##
GenerateSGLReg <- function(SMS, Interval = Interval,
  Criterion = "BIC") {
  MyData <- list(y=SMS$YY, x=SMS$XX, n=length(SMS$YY));
	T1 <- proc.time();
  library(SGL)	
	## grpreg         
  MyT <- "ATT <- FALSE; 
    library(SGL, warn.conflicts=FALSE, quietly=TRUE);  
    ATT <- TRUE;"
  try(eval(parse(text=MyT)));
  if (ATT == FALSE) {
    AFilePrint("GenerateGroupGrpReg: FAIL no grpreg package! installed.");
     1 = 2;
  }
  
  
  cv.fit = cvSGL(data=MyData, index=SMS$index, type = "linear")
  AFitT <- cv.fit$fit$beta;
  AFitT[cv.fit$fit$beta != 0.0] <- 1.0;
  SS <- (colSums(AFitT) - SMS$kActiveLen)^2
  Minerizer <- sort.int(SS, index.return=TRUE)$ix[1];
  ##lambdaFit <- cv.fit$lambda.min;
  BetaWant <- cv.fit$fit$beta[, Minerizer];
  T2 <- proc.time();
  StandGLFit <- BetaWant;
  StandGLFitBBs <- rep(0, length(BetaWant));
  StandGLFitBBs[BetaWant != 0.0] = 1.0;

  MIPRandom <- sort(unique(SMS$index));
  for (kk in 1:length(MIPRandom)) {
    MIPRandom[kk] <- mean(StandGLFitBBs[SMS$index==kk])
  }

  if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE &&
    length(StandGLFit) == NCOL(SMS$XX)+1) {
    StandGLFit <- StandGLFit[2:length(StandGLFit)];
    StandGLFitBBs <- StandGLFitBBs[2:length(StandGLFitBBs)];    
  }

	SGLFitTime <- T2 - T1;
	return(list( type="SGLCVLasso", BetaFit = StandGLFit,
    BBHatsFit = StandGLFitBBs, FitTime = SGLFitTime,
	  OtherBetas = StandGLFit, StandGLFitObject = cv.fit, CI = NULL,
    CIEst = NULL, CITime = NULL, MIPReport = StandGLFitBBs, MIPRandom=MIPRandom) );        		     
}


WriteFileInfoT <- function(FileDirectory, FileName = "Category.R", NameFile="Category.R", ...) {
  eval(parse(text=GetG0Text("NN", "globalenv", S=1)));
  eval(parse(text=GetG0Text("n", "globalenv", S=1)));
  eval(parse(text=GetG0Text("sigma", "globalenv", S=1)));
  eval(parse(text=GetG0Text("NR", "globalenv", S=1)));
  eval(parse(text=GetG0Text("p", "globalenv", S=1)));
  eval(parse(text=GetG0Text("k", "globalenv", S=1)));
  if (!is.null(NameFile) && is.character(NameFile) && NameFile != "Category.R" && 
    (is.character(FileName) && FileName == "Category.R")) {
    FileName = FileName;  
  }

  MyF <- file(description = paste(FileDirectory, "//", FileName, sep=""), 
    open = "w", blocking = TRUE,
    encoding = getOption("encoding"));
  writeLines(con=MyF, tex=paste("NN = ", NN, ";", sep=""));
  writeLines(con=MyF, tex=paste("n = ", n, ";", sep=""));
  writeLines(con=MyF, tex=paste("sigma = ", sigma, ";", sep=""));  
  writeLines(con=MyF, tex=paste("NR = ", NR, ";", sep="")); 
  writeLines(con=MyF, tex=paste("p = ", p, ";", sep="")); 
  writeLines(con=MyF, tex=paste("k = ", k, ";", sep=""));   
  close(MyF);
}
