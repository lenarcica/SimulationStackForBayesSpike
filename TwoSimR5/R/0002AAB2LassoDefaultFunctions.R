###############################################################################
##  0002AAB2LassoDefaultFunctions  -- Alan Lenarcic
##
##    (c) 2009-2019
##
##    This R code installs versions of simulation code defined in 0001AAAEfficientSimulator.r
##   into the R prompt with certain default parameter options defined.  It also
##   gives names to the output of these functions that can be used in retrieval
##   including Latex format names.  Each "function * parameter setting" is given
##   it's own "Generate_XXXX()" function name which curries in the implied
##   estimator approach.  For instance, a researcher might want to 
##   input "oracle" k into Elastic Net, that is tell the Net exactly the model
##   size and it will pick the best fit with this value.  Instead the
##   researcher might give a "noisy k", or some number close to what they believe
##   is probably k.  Both would use the same code from 0001AAAEfficientSimulator
##   however, the two approaches to supplying input parameters are simulated 
##   differently.  The same is true for hyper-prior approaches for GroupBayes=BayesSpike
##
##   There are many estimators in R for sparse fixed model selection, and
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

Interval = "";  piAPriorData = 0;

GenerateLarsFixedR <- function(SMS) {
  GenerateLarsFixed(SMS); 
}
GenerateLarsFixedNoiseR <- function(SMS) {
  GenerateLarsFixedR(SMS$SMSNoise); 
}
GenerateENetFixedR <- function(SMS) {
  GenerateLarsFixed(SMS); 
}
GenerateENetFixedNoiseR <- function(SMS) {
  GenerateLarsFixed(SMS$SMSNoise); 
}
GenerateENetCVMinR <- function(SMS) {
  GenerateENetCVMin(SMS); 
}
GenerateLarsCpR <- function(SMS=NULL) {
  GenerateLarsCp(SMS); 
}
GenerateLarsCpNoiseR <- function(SMS=NULL) {
  if (is.null(SMS$SMSNoise)) {
    return(GenerateLarsCp(SMS)); 
  } else {
    return(GenerateLarsCp(SMS$SMSNoise)); 
  }
}
GenerateLassoStoddenR <- function(SMS) {
  GenerateLassoStodden(SMS, LARSSeekFlag);
}
GenerateLassoStoddenNoiseR <- function(SMS) {
  GenerateLassoStodden(SMS$SMSNoise, LARSSeekFlag);
}  
GenerateLimitLassoR <- function(SMS) {
  GenerateLimitLasso(SMS, LARSSeekFlag);
}

GenerateBayesVarSelR <- function(SMS) {
  GenerateBayesVarSel(SMS);
}

GenerateLimitLassoNoiseR <- function(SMS) {
  GenerateLimitLasso(SMS$SMSNoise, LARSSeekFlag);
}

GenerateBayesLassoR <- function(SMS) {
 GenerateBayesLasso(SMS);
}

GeneratepiALimitLassoR <- function(SMS) {
  GeneratepiALimitLasso(SMS, LARSSeekFlag, 
    musize = piAPriorData * n, SigmaEta = sigmaPriorData * n, 
    SigmaBarS = SMS$SigmaSqNoise);
}

GeneratepiALimitLassoNoiseR <- function(SMS) {
  GeneratepiALimitLasso(SMS$SMSNoise, LARSSeekFlag, 
    musize = piAPriorData * n, SigmaEta = sigmaPriorData * n, 
    SigmaBarS = SMS$SigmaSqNoise);
}


GeneratepiAM2LassoR <- function(SMS) {
  GeneratepiAM2Lasso(SMS, LARSSeekFlag, musize = piAPriorData * n, 
    SigmaEta = sigmaPriorData * n, SigmaBarS = SMS$SigmaSqNoise);
}
GeneratepiAM2LassoNoiseR <- function(SMS) {
  GeneratepiAM2Lasso(SMS$SMSNoise, LARSSeekFlag, musize = piAPriorData * n, 
    SigmaEta = sigmaPriorData * n, SigmaBarS = SMS$SigmaSqNoise);
}
GeneratepiAR2LassoR <- function(SMS) {
  GeneratepiAM2Lasso(SMS, LARSSeekFlag, 
    musize = piAPriorData * n, SigmaEta = sigmaPriorData * n,
    SigmaBarS = SMS$SigmaSqNoise, R = .5);
}
GeneratepiAR2LassoNoiseR <- function(SMS) {
  GeneratepiAM2Lasso(SMS$SMSNoise, LARSSeekFlag, 
    musize = piAPriorData * n, SigmaEta = sigmaPriorData * n,
    SigmaBarS = SMS$SigmaSqNoise, R = .5);
}
GenerateFDLimitLassoR <- function(SMS) {
  GenerateFDLimitLasso(SMS, LARSSeekFlag);
}
GenerateFDLimitLassoNoiseR <- function(SMS) {
  GenerateFDLimitLasso(SMS$SMSNoise, LARSSeekFlag);
}
GenerateFDM2LassoR <- function(SMS) {
  GenerateFDM2Lasso(SMS, LARSSeekFlag);
}
GenerateFDM2LassoNoiseR <- function(SMS) {
  GenerateFDM2Lasso(SMS$SMSNoise, LARSSeekFlag);
}
GenerateFDR2LassoR <- function(SMS){
  GenerateFDM2Lasso(SMS, LARSSeekFlag,R = .5);
}
GenerateHorseShoeR <- function(SMS)  {
  eval(parse(text=GetG0Text("DefaultHorseShoeSamples", "globalenv()", S=1)));
  GenerateHorseShoe(SMS, tauSq = 1, CountSamp = DefaultHorseShoeSamples) 
}
GenerateHorseShoeOldR <- function(SMS)  {
  eval(parse(text=GetG0Text("DefaultHorseShoeSamples", "globalenv()", S=1)));
  GenerateHorseShoeOld(SMS, tauSq = 1, CountSamp = DefaultHorseShoeSamples) 
}
GenerateSpikeAndSlabR <- function(SMS) {
  eval(parse(text=GetG0Text("DefaultHorseShoeSamples", "globalenv()", S=1)));
  GenerateSpikeAndSlab(SMS, tauSqA = 1, CountSamp = DefaultHorseShoeSamples);
}
GenerateSpikeAndSlabNoiseR <- function(SMS) {
  eval(parse(text=GetG0Text("DefaultHorseShoeSamples", "globalenv()", S=1)));
  GenerateSpikeAndSlab(SMS$SMSNoise, tauSqA = 1, CountSamp = DefaultHorseShoeSamples);
}
GenerateXLLimitLassoR <- function(SMS) {
  GenerateXLLimitLasso(SMS,GenerateLimitLassoR(SMS), Interval=Interval); 
}
GenerateXLLimitLassoNoiseR <- function(SMS) {
  GenerateXLLimitLasso(SMS,GenerateLimitLassoR(SMS$SMSNoise), Interval=Interval); 
}
GenerateXLM2LassoR <- function(SMS) {
  GenerateXLLimitLasso(SMS,GenerateTwoLassoSpreadR(SMS), Interval=Interval); 
}
GenerateXLM2LassoNoiseR <- function(SMS) {
  GenerateXLLimitLasso(SMS$SMSNoise,GenerateTwoLassoSpreadR(SMS$SMSNoise), Interval=Interval); 
}
GenerateXLR2LassoR <- function(SMS) {
  GenerateXLLimitLasso(SMS, GenerateR2Lasso(SMS, R = .5),
    Interval = Interval);
}
GenerateR2LassoR <- function(SMS) {
  GenerateR2Lasso(SMS, R = .5);
}
GenerateIgnorantLimitLassoR <- function(SMS) {
  GenerateLimitLasso(SMS, LARSSeekFlag, Ignorant =TRUE);
}
GenerateB2LassoR <- function(SMS) {
  GenerateB2Lasso(SMS);
}
GenerateTwoLasso9XR <- function(SMS) {
  GenerateTwoLasso9X(SMS, LARSSeekFlag);
}
GenerateTwoLasso9XNoiseR <- function(SMS) {
  GenerateTwoLasso9X(SMS$SMSNoise, LARSSeekFlag);
}
GenerateEMRidgeR <- function(SMS) {
  GenerateEMRidge(SMS, LARSSeekFlag);
}
GenerateTwoLassoSpreadR <- function(SMS)  {
  GenerateTwoLassoSpread(SMS); 
}

GenerateBayesSpikeInfoR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Info", TestAtGoal = FALSE); 
}
GenerateBayesSpikeAutoR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", TestAtGoal = FALSE); 
}

GenerateBayesSpikeAutoSQRTR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
    NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
    CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="AutoSqrt", TestAtGoal = FALSE);  
}

GenerateBayesSpikeAutoSQRTMeanR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
    NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
    CutOff = .1, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="AutoSqrt", TestAtGoal = FALSE);  
}

GenerateBayesSpikeAutoLogitSlowR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto",
  DoLogitPostPreProb = 2); 
}
GenerateBayesSpikeInfoLogitSlowR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Info",
  DoLogitPostPreProb = 2); 
}

GenerateBayesSpikeAutoMeanR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="Auto", TestAtGoal = FALSE); 
}

GenerateBayesSpikeAutoMeanR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="Auto", TestAtGoal = FALSE); 
}

GenerateBayesSpikeInfoMeanR <- function(SMS) {
  GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="Info", TestAtGoal = FALSE); 
}

GenerateBayesSpikeInfoMeanNoiseR <- function(SMS) {
  GenerateBayesSpike(SMS$SMSNoise, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength,  SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="Info", TestAtGoal = FALSE); 
}

GenerateBayesSpikeInfoNoiseR <- function(SMS) {
  GenerateBayesSpikeInfoR(SMS$SMSNoise); 
}

  
GenerateEMRidgeNoiseR <- function(SMS) {
  GenerateEMRidge(SMS$SMSNoise, LARSSeekFlag);
}

GenerateGroupBayesSpikeAutoR <- function(SMS) {
   GenerateGroupBayesSpike(SMS, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto");      
}

GenerateGroupBayesSpikeInfoR <- function(SMS) {
   GenerateGroupBayesSpike(SMS, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Info");      
}
GenerateGroupBayesSpikeAutoLogitSlowR <- function(SMS) {
   if (SMS$p > 5000) {
      print("GenerateGroupBayesSpikeInfoLogitSlowR: Too Much data to do this!"); flush.console();
      return(NULL);
   }
   return(GenerateGroupBayesSpike(SMS, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto",
   DoLogitPostPreProb = 1));      
}
GenerateGroupBayesSpikeInfoLogitSlowR <- function(SMS) {
   if (SMS$p > 5000) {
      print("GenerateGroupBayesSpikeInfoLogitSlowR: Too Much data to do this!"); flush.console();
      return(NULL);
   }
   return(GenerateGroupBayesSpike(SMS, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Info",
   DoLogitPostPreProb = 1));      
}
GenerateGroupBayesSpikeInfoNoiseR <- function(SMS) {
   GenerateGroupBayesSpike(SMS$SMSNoise, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Info");      
}
GenerateGroupBayesSpikeAutoMeanR <- function(SMS) {
   GenerateGroupBayesSpike(SMS, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="Auto");      
}
GenerateGroupBayesSpikeAutoSQRTMeanR <- function(SMS) {
   GenerateGroupBayesSpike(SMS, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="AutoSqrt");      
}
GenerateGroupBayesSpikeInfoMeanR <- function(SMS) {
   GenerateGroupBayesSpike(SMS, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="Info");      
}
GenerateGroupBayesSpikeInfoMeanNoiseR <- function(SMS) {
   GenerateGroupBayesSpike(SMS$SMSNoise, DoCI = DefaultDoCI, DoMedian = FALSE, AutoInfo="Info");      
}
GenerateGroupTwoLassoSpreadNI <- function(SMS) {
  GenerateGroupTwoLassoSpread(SMS$SMSNoise);
}

GenerateBayesSpikeAutoTemp1R <- function(SMS) {
  ## 1 50000,10000
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 6200, burnin=200, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", TemperatureList=1, NumTemp = -1, TestAtGoal=TRUE)); 

  ##GenerateGroupBayesSpike(SMS, DoCI = TRUE, DoMedian = FALSE, AutoInfo="Auto", NumTemp=1,TestAtGoal=TRUE);
  ##GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  ##CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  ##NumTemp = 1, TestAtGoal = FALSE,  ...)
}

GenerateBayesSpikeAutoTemp2R <- function(SMS) {
  #1 50000 10000
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 6200, burnin =200, 
  NumChains = 1, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", TemperatureList=c(1.25,1),TestAtGoal=TRUE,
  EEMergeEvery = 10, RevertTemperatureEvery = 50, EEProbSortWidth = .25)); 
  ##GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  ##CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  ##NumTemp = 1, TestAtGoal = FALSE,  ...)
}

GenerateBayesSpikeAutoTemp2c10R <- function(SMS) {
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", NumTemp=2,TestAtGoal=TRUE,
  EEMergeEvery = 10)); 
  ##GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  ##CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  ##NumTemp = 1, TestAtGoal = FALSE,  ...)
}

GenerateBayesSpikeAutoTemp2c05R <- function(SMS) {
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", NumTemp=2,TestAtGoal=TRUE,
  EEMergeEvery = 5)); 
  ##GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  ##CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  ##NumTemp = 1, TestAtGoal = FALSE,  ...)
}

GenerateBayesSpikeAutoTemp2c20R <- function(SMS) {
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", NumTemp=2,TestAtGoal=TRUE,
  EEMergeEvery = 20)); 
  ##GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  ##CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  ##NumTemp = 1, TestAtGoal = FALSE,  ...)
}


GenerateBayesSpikeAutoTemp3c10R <- function(SMS) {
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", NumTemp=3,TestAtGoal=TRUE,
  EEMergeEvery = 10)); 
  ##GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  ##CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  ##NumTemp = 1, TestAtGoal = FALSE,  ...)
}

GenerateBayesSpikeAutoTemp3c05R <- function(SMS) {
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", NumTemp=3,TestAtGoal=TRUE,
  EEMergeEvery = 5)); 
  ##GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  ##CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  ##NumTemp = 1, TestAtGoal = FALSE,  ...)
}

GenerateBayesSpikeAutoTemp3c20R <- function(SMS) {
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", NumTemp=3,TestAtGoal=TRUE,
  EEMergeEvery = 20)); 
  ##GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  ##CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  ##NumTemp = 1, TestAtGoal = FALSE,  ...)
}

GenerateBayesSpikeAutoTemp3R <- function(SMS) {
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 6200, burnin=200, 
  NumChains = 1, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", TemperatureList=c(1.5,1), NumTemp = -1, TestAtGoal=TRUE,
  EEMergeEvery = 10, RevertTemperatureEvery = 50, EEProbSortWidth = .25)); 

}
GenerateBayesSpikeAutoTemp4R <- function(SMS) {
  return(GenerateBayesSpike(SMS, tauSqA = 1, CountSamp = 6200,burnin = 200, 
  NumChains = 1, PriorStrength = SMS$PriorStrength, SigmaPriorStrength = SMS$SigmaPriorStrength,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, AutoInfo="Auto", TemperatureList=c(1.5,1.25,1), NumTemp = -1, TestAtGoal=TRUE,
  EEMergeEvery = 10, RevertTemperatureEvery = 50, EEProbSortWidth = .25)); 
}
GenerateSGLRegNI <- function(SMS) {
  GenerateSGLReg(SMS$SMSNoise);   
}
GenerateTwoLassoSpreadNoiseR <- function(SMS)  {
  GenerateTwoLassoSpread(SMS$SMSNoise); 
}
GenerateSCADCVMR <- function(SMS)  {
  return( GenerateSCAD(SMS, LARSSeekFlag, DoConvexMin = 1)); 
}
GenerateSCADFourVR <- function(SMS)  {
  return( GenerateSCAD(SMS, LARSSeekFlag, DoConvexMin = 0)); 
}
GenerateSCADFourVNoiseR <- function(SMS)  {
  return( GenerateSCAD(SMS$SMSNoise, LARSSeekFlag, DoConvexMin = 0)); 
}
GenerateSCADMinR <- function(SMS)  {
  ### WLT method
  return( GenerateSCAD(SMS, LARSSeekFlag, DoConvexMin =-1)); 
}

GenerateMCPMinR <- function(SMS)  {
  return( GenerateMCP(SMS, LARSSeekFlag, DoConvexMin =-1)); 
}

GenerateMCPCVMR <- function(SMS)  {
  return( GenerateMCP(SMS, LARSSeekFlag, DoConvexMin = 1)); 
}

GT <- function() {
  GenerateLarsFixed(stuff);
}
GenerateXLMMedianR <- function(SMS) {
  GenerateXLMMedian(SMS,Interval)
}
GenerateXLMMedianNoiseR <- function(SMS) {
  GenerateXLMMedian(SMS$SMSNoise,Interval)
}
GenerateLassoQLYR <- function(SMS) {
  GenerateLassoQLY(SMS);
}
GenerateLassoQLYNoiseR <- function(SMS) {
  GenerateLassoQLY(SMS$SMSNoise);
}


WriteFileInfoT <- function(FileDirectory, FileName = "Category.R", NameFile = "Category.R", ...) {
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


Generate2LassoCVR <- function(SMS) {
  Generate2LassoCV(SMS);
}




OlderUseFunctions <- function() {
UseFunctions <- c(GenerateLarsFixed, GenerateLarsCp, GenerateLassoStoddenR,
  GenerateLimitLassoR, GeneratepiALimitLassoR, 
  GenerateFDLimitLassoR,  GenerateXLMMedianR,                 
  GenerateXLLimitLassoR, GenerateIgnorantLimitLassoR,
  GenerateLassoQLYR, GenerateTwoLasso9XR, GenerateEMRidgeR,
  GenerateTwoLassoSpreadR, GenerateSCADCVMR,
  GeneratepiAM2LassoR,GenerateFDM2LassoR,
  GenerateXLM2LassoR,GenerateSCADFourVR, GenerateSCADMinR,
  GenerateENetFixedR, GenerateENetCVMinR,
  GenerateR2LassoR, GeneratepiAR2LassoR,
  GenerateFDR2LassoR, GenerateXLR2LassoR,
  GenerateHorseShoeR, GenerateSpikeAndSlabR, GenerateB2LassoR,
  GenerateBayesLassoR, Generate2LassoCVR, GenerateBayesVarSelR,
  GenerateHorseShoeOldR);                 
NameFunctions <- c("LarsFixed", "LarsCp", "LassoStodden",
    "LimitLasso", "piALimitLasso", "FDLimitLasso",
    "XLMMedian", "XLLimitLasso", "IgnorantLimitLasso", "LassoLY", "TwoLasso9X", 
    "LimitRidge", "2LassoSpread", "SCADminConvex", "piAM2Lasso", "FDM2Lasso", "XLM2Lasso",
    "SCAD4xStodden", "SCADMinLRY", "ENetFixed", "ENetCV",
    "R2Lasso", "R2PiLasso", "R2FDLasso", "R2XLLasso", "HorseShoe",
    "SpikeAndSlab", "B2Lasso", "BayesLasso", "2LassoCV",
    "BayesVarSel", "HorseShoeAL");
}
DefineDefaultUseFunction <- function () {
eval(parse(text=GetG0Text("UseFunctions", S=1)));
eval(parse(text=GetG0Text("DefaultUseFunctions", S=1)));
  DefaultAllUseFunctions <- c(GenerateLarsFixed, GenerateLarsCp, GenerateLassoStoddenR,
  GenerateLimitLassoR, GeneratepiALimitLassoR, 
  GenerateFDLimitLassoR,  GenerateXLMMedianR,                 
  GenerateXLLimitLassoR, GenerateIgnorantLimitLassoR,
  GenerateLassoQLYR, GenerateTwoLasso9XR, GenerateEMRidgeR,
  GenerateTwoLassoSpreadR, GenerateSCADCVMR,
  GeneratepiAM2LassoR,GenerateFDM2LassoR,
  GenerateXLM2LassoR,GenerateSCADFourVR, GenerateSCADMinR,
  GenerateENetFixedR, GenerateENetCVMinR,
  GenerateHorseShoeR, GenerateSpikeAndSlabR, GenerateB2LassoR, GenerateBayesLassoR,
  Generate2LassoCVR, GenerateISpike,
  GenerateBayesSpikeInfoR, GenerateBayesSpikeAutoR,
  GenerateLarsFixedNoiseR, GenerateLarsCpNoiseR, GenerateLassoStoddenNoiseR,
  GenerateLimitLassoNoiseR, GeneratepiALimitLassoNoiseR, 
  GenerateFDLimitLassoNoiseR,  GenerateXLMMedianNoiseR,                 
  GenerateXLLimitLassoNoiseR,
  GenerateLassoQLYNoiseR, GenerateTwoLasso9XNoiseR, GenerateEMRidgeNoiseR,
  GenerateTwoLassoSpreadNoiseR,
  GeneratepiAM2LassoNoiseR,GenerateFDM2LassoNoiseR,
  GenerateXLM2LassoNoiseR,  GenerateSCADFourVR,
  GenerateENetFixedNoiseR, 
  GenerateBayesSpikeInfoNoiseR, GenerateBayesSpikeAutoMeanR,
  GenerateBayesSpikeInfoMeanR, GenerateBayesSpikeInfoMeanNoiseR,
  GenerateMCPCVMR, GenerateMCPMinR,
  GenerateBayesSpikeAutoLogitSlowR, GenerateBayesSpikeInfoLogitSlowR,
  GenerateBayesVarSelR, GenerateHorseShoeOldR,
  GenerateBayesSpikeAutoSQRTR,
  GenerateBayesSpikeAutoSQRTMeanR,
  GenerateBayesSpikeAutoTemp1R,
  GenerateBayesSpikeAutoTemp2R,
  GenerateBayesSpikeAutoTemp3R,
  GenerateBayesSpikeAutoTemp4R,
  GenerateBayesSpikeAutoTemp2c10R,
  GenerateBayesSpikeAutoTemp2c05R,
  GenerateBayesSpikeAutoTemp2c20R,
  GenerateBayesSpikeAutoTemp3c05R,
  GenerateBayesSpikeAutoTemp3c10R,
  GenerateBayesSpikeAutoTemp3c20R
  );       
  DefaultUseFunctions <- DefaultAllUseFunctions              
UseFunctions <- DefaultUseFunctions;       
eval(parse(text=SetGText("UseFunctions", S=1)));
eval(parse(text=SetGText("DefaultUseFunctions", S=1)));
eval(parse(text=GetG0Text("TWOSIMENVIRONMENT", S=1)));
if (!is.environment(TWOSIMENVIRONMENT)) {
  TWOSIMENVIRONMENT <- evironment();
} 
eval(parse(text=SetGText("UseFunctions", 
  envir="TWOSIMENVIRONMENT", S=1))); 
eval(parse(text=SetGText("DefaultUseFunctions", 
  envir="TWOSIMENVIRONMENT", S=1))); 
  
eval(parse(text=GetG0Text("DefaultNameFunctions", S=1)));      
eval(parse(text=GetG0Text("DefaultAllGroupFunctions", S=1)));
eval(parse(text=GetG0Text("NameFunctions", S=1)));  
eval(parse(text=GetG0Text("OriginalOldwd", S=1)));
OriginalOldwd <- getwd();
eval(parse(text=SetGText("OriginalOldwd", S=1)));           
DefaultNameFunctions <- c("LarsFixed", "LarsCp", "LassoStodden",
    "LimitLasso", "piALimitLasso", "FDLimitLasso",
    "XLMMedian", "XLLimitLasso", "IgnorantLimitLasso", "LassoLY", "TwoLasso9X", 
    "LimitRidge", "2LassoSpread", "SCADminConvex", "piAM2Lasso", "FDM2Lasso", "XLM2Lasso",
    "SCAD4xStodden", "SCADMinLRY", "ENetFixed", "ENetCV",
    "HorseShoe",
    "OurFirstSpike", "B2Lasso", "BayesLasso", "2LassoCV", "IshwaranSpike",
    "GenerateBayesSpikeInfoPrior", "GenerateBayesSpikeAutoPrior",
    "LarsFixedNI", "LarsCpNI", "LassoStoddenNI",
    "LimitLassoNI", "piALimitLassoNI", "FDLimitLassoNI",
    "XLMMedianNI", "XLLimitLassoNI",
    "LassoLYNI", "TwoLasso9XNI", 
    "LimitRidgeNI", "2LassoSpreadNI",
    "piAM2LassoNI", "FDM2LassoNI", "XLM2LassoNI", "SCAD4xStoddenNI", 
    "ENetFixedNI",
    "GenerateBayesSpikeInfoPriorNI", "GenerateBayesSpikeAutoMean",
    "GenerateBayesSpikeInfoMean",     "GenerateBayesSpikeInfoMeanNI",
    "MCPCVMR", "MCPMinR",
    "GenerateBayesSpikeAutoLogitSlow", "GenerateBayesSpikeInfoLogitSlow",
    "BayesVarSel", "HorseShoeAL",
    "GenerateBayesSpikeAutoSQRTR",
    "GenerateBayesSpikeAutoSQRTMeanR",
  "GenerateBayesSpikeAutoTemp1R",
  "GenerateBayesSpikeAutoTemp2R",
  "GenerateBayesSpikeAutoTemp3R",
  "GenerateBayesSpikeAutoTemp4R",
  "GenerateBayesSpikeAutoTemp2c10R",
  "GenerateBayesSpikeAutoTemp2c20R",
  "GenerateBayesSpikeAutoTemp2c05R",
  "GenerateBayesSpikeAutoTemp3c10R",
  "GenerateBayesSpikeAutoTemp3c20R",
  "GenerateBayesSpikeAutoTemp3c05R"     
  ); 
NameFunctions <- DefaultNameFunctions;
eval(parse(text=SetGText("NameFunctions", S=1)));
eval(parse(text=SetGText("NameFunctions", 
  envir="TWOSIMENVIRONMENT", S=1))); 
eval(parse(text=SetGText("DefaultNameFunctions", S=1)));
eval(parse(text=SetGText("DefaultNameFunctions", 
  envir="TWOSIMENVIRONMENT", S=1))); 
 
eval(parse(text=GetG0Text("DefaultAllGroupFunctions", S=1))); 
DefaultAllGroupFunctions <- c(
  GenerateGroupBayesSpikeAutoR, 
  GenerateGroupBayesSpikeInfoR, GenerateGroupBayesSpikeInfoNoiseR,
  GenerateGroupBayesSpikeAutoMeanR,GenerateGroupBayesSpikeInfoMeanR,
  GenerateGroupBayesSpikeInfoMeanNoiseR,
  GenerateGroupTwoLassoSpread,
  GenerateGroupLassoReg,
  GenerateGroupGrpReg, GenerateStandGLReg, GenerateSGLReg,
  GenerateGroupTwoLassoSpreadNI,
  GenerateSGLRegNI,
  GenerateGroupBayesSpikeAutoLogitSlowR,
  GenerateGroupBayesSpikeInfoLogitSlowR
);

eval(parse(text=SetGText("DefaultAllGroupFunctions", S=1)));
eval(parse(text=SetGText("DefaultAllGroupFunctions", 
  envir="TWOSIMENVIRONMENT", S=1))); 
DefaultAllGroupNameFunctions <- c(
  "GroupBayesSpikeAuto", "GroupBayesSpikeInfo",
  "GroupBayesSpikeSpikeInfoNoise", "GroupBayesSpikeAutoMean",
  "GroupBayesSpikeInfoMean", "GroupBayesSpikeInfoMeanNoise",
  "GenerateGroupTwoLassoSpread", 
  "GenerateGroupLassoReg", "GenerateGroupGrpReg",
  "GenerateStandGLReg", "GenerateSGLReg",
  "GenerateGroupTwoLassoSpreadNI",
  "GenerateSGLRegNI",
  "GenerateGroupBayesSpikeAutoLogitSlowR",
  "GenerateGroupBayesSpikeInfoLogitSlowR"
);
eval(parse(text=SetGText("DefaultAllGroupNameFunctions", S=1)));
eval(parse(text=SetGText("DefaultAllGroupNameFunctions", 
  envir="TWOSIMENVIRONMENT", S=1))); 
eval(parse(text=GetG0Text("OnSimType")));
  if (!is.null(OnSimType) && is.character(OnSimType) &&
  OnSimType %in% c("Group", "GroupRobit", "GroupLogit")) {
     NameFunctions <- DefaultAllGroupNameFunctions;
     eval(parse(text=SetGText("NameFunctions", S=1)));
     eval(parse(text=SetGText("NameFunctions", 
       envir="TWOSIMENVIRONMENT", S=1))); 
     UseFunctions  <- DefaultAllGroupFunctions;
     eval(parse(text=SetGText("UseFunctions", S=1)));
     eval(parse(text=SetGText("UseFunctions", 
       envir="TWOSIMENVIRONMENT", S=1))); 
  }
}

