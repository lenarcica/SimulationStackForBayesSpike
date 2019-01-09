################################################################################
## RunAProcessExperiment in ProcessRunFile.r 
##  (c) Alan Lenarcic 2009-2019
##
##  RunAProcessExperiment() run a estimator(by "FunctionName") against a simulation contained in "SMSName"
##   The estimator is fit against the data and then assessments of performance are conducted and the results
##   saved to disk in a separate file and directory so that the next time RunAProcessExperiment is run
##   it knows that the result does not need to be run.
##
## This will try to run an estimator on a given dataset and conclude selection properties.
##   The TwoSimR5 approach to conducint simulations is
##     1. Decide upon which N_M estimators to use in a problem, install functions
##      through EfficientSimulator and declare parameter settings with DeclareAllTstFunctions
##     2. Define a simulation (n,p,k, sigma, covarianceX, size of Beta...)
##     2. Simulate N_S 500-1000 indepenent simulations from t  (AssessFirstTime() which runs "SimMeData()" as needed)
##     3. RunAProcessExperiment()" Open N_M * (N_S/N_divisor) threads on Killdevil to attempt to solve each simulation with each estimator
##        (Killdevil requires 1 minute or more processes with 1 hour kill time, so some threads will solve many problems,
##           so each thread gets to solve multiple randomly picked estimator sand simulation problems)
##        Each solution must:
##          a. fit the estimator
##          b. Assess the L1, L2, etc error of the estimator
##          c. Assess Credibility or Model Inclusion performance if applciable
##          d. Save results on this estimator somewhere to remind program to not rerun.
##     4. Rerun threads when an estimator may have failed due to a timeout or memory error thus every
##       estimator gets at least two chances to complete a run.  Most time they don't need any, but sometimes
##       a problem is begun right before a 1 hour timeout.  As problems get larger certain R packages hit
##       memory limits, and this gives us a chance to double check.
##     5.Collect the results of all of these simulations by estimator
##     6.Produce Latex usable tables of results or other interesting RData.   
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
RunAProcessExperiment =  function(SMSName, FunctionName, verbose = -1, DoCI = DefaultDoCI, 
    ForceAFail=FALSE, ...) {
  eval(parse(text=InitiateEvaluations()));
  eval(parse(text=TestForVerbose()));
  eval(parse(text=LoadSMSNow()));  
  eval(parse(text=FunctionUnderstandProcesses()));
  eval(parse(text=EvaluateProgress())); 
  eval(parse(text=EstablishFunctionAndIdentifiers()));
  
  RTI = unlist(UseFunctions[IFunction]);
  eval(parse(text=SetGText("SMS", S=1)));
    ##AFilePrint(paste(" Doing Analysis NameFunctions[", ii, "] = ",
    ##   NameFunctions[ii], sep="")); 
  if (p > 1000) {
    AFilePrint(paste("Doing Analysis NameFunctions[", IFunction, "] = ",
        OnNameFunctions[IFunction], sep=""));  flush.console();
  }
  limlas = NULL;
     
  try(MyOldwd <- getwd());
  try(setwd(CrashDir));
  try(limlas <- (RTI[[1]])(SMS));   ## Now Run the Function!
  try(setwd(MyOldwd));
  ## AL <- GenerateGroupGrpReg(SMS)
  eval(parse(text=SetGText("limlas", S=1)));
  if (verbose > 0) {
    AFilePrint(paste("We just ran function: ", NameFunctions[IFunction], sep=""));
    flush.console();
  }
  if (exists("ForceAFail") && is.logical(ForceAFail) && ForceAFail == TRUE) {
    AFilePrint(paste("We ran function ", NameFunctions[IFunction], 
      " but deliberately call it a failure. "));
    limlas <- NULL;
  } else {
    ForceAFail <- FALSE;
  }
  if (is.null(limlas)) {
    eval(parse(text=GetG0Text("NameFunctions")));
    AFilePrint(paste("We've got a limlas error, we are on ii = ", IFunction, " and NameFunctions = ",
       NameFunctions[IFunction], sep="")); flush.console();
    beta.hat = SMS$BetasReal * 0 - 666;
    l1error = -2;
    l2error = -2;
    Type1 = -2; Type2 = -2; 
    L.zero <- -2 
    time = c(NA,NA,NA);
    Hit = "D"
    FittedOb = NULL;
  } else if (!is.null(limlas)) {
    FittedOb = limlas;
    beta = SMS$BetasReal
  	beta.hat = limlas$BetaFit
  	if (exists("IFunction") && IFunction == 1) {
      BetaStart = beta.hat;
      eval(parse(text=SetGText("BetaStart", S=1)));
    }
  	
  	# l1 and l2 errors
  	l1error <- lnorm(beta.hat - beta, 1) / lnorm(beta, 1)
  	l2error <- lnorm(beta.hat - beta, 2) / lnorm(beta, 2)
    if (beta.hat[1] == -999) {
      l2error = -999;
      l1error = -999;
    }
    # Type1 and Type2 errors

    if (!is.null(limlas$BBHatsFit)) {
      if (!is.null(SMS$TrueOnGroups) && length(SMS$TrueOnGroups) == SMS$pGroups) {
         TOn <- rep(0,SMS$pGroups);
         ish <- SMS$FirstGroupIndex;
         for (jti in 1:SMS$pGroups) {
            if (any(limlas$BBHatsFit[ish:SMS$EndGroupIndices[jti]] != 0)) {
              TOn[jti] = 1;
            }
            ish <- SMS$EndGroupIndices[jti]+1;
         }
         Type1 <- sum(((SMS$TrueOnGroups - TOn)==-1)*1)
         Type2 <- sum(((-1*SMS$TrueOnGroups + TOn)==-1)*1)
      } else {
        Type1 <- sum(((SMS$BBsReal - limlas$BBHatsFit)==-1)*1)
        Type2 <- sum(((-1*SMS$BBsReal + limlas$BBHatsFit)==-1)*1)
      }
    } else {Type2 <- NA; Type1 <- NA; }
    # L_0 Error
    L.zero <- 0.5*Type1 + 0.5*Type2  	
    if (SMS$BBsReal[1] <= -1) {
      Type1 = -2;
      Type2 = -2;
      L.zero = -2;
    }
  

  	if (length(limlas$FitTime) >= 3) {
     	time <- as.numeric(limlas$FitTime[1:3])
    } else {
      time <- c(NA, NA, NA);
    }
    eval(parse(text=GetG0Text("CurrentSmallContainDir")));
    if (p > 5000) {
      try(paste("    Onjobii = ", jobii, 
        " Finished ", IFunction, " = ", NameFunctions[IFunction], 
        " in time ", sum(time[1:2]), sep=""));
      try(dir.create(CurrentSmallContainDir, showWarnings = FALSE));
      try(FilesInDir <- unlist(list.files(CurrentSmallContainDir)));
      MyTT <- "
      if (any(FilesInDir == \"Progress.txt\")) {
        try(fp <- file(paste(CurrentSmallContainDir, \"//Progress.txt\", sep=\"\"), \"a\"))
      }   else {
        try(fp <- file(paste(CurrentSmallContainDir, \"//Progress.txt\", sep=\"\"), \"w\"))
      }
      ";
      try(eval(parse(text=MyTT)));
      try(writeLines(text=paste("Onjobii = ", jobii, 
                         " Finished ", IFunction, " = ", NameFunctions[IFunction], 
                         " in time ", sum(time[1:2]), sep=""),
                   con = fp));
      try(close(fp));
    }
   
   Hit = "+"
   }
    if (!is.numeric(l1error) || is.na(l1error) || !is.finite(l1error) || length(l1error) != 1) { l1error = NA; }
    if (!is.numeric(l2error) || is.na(l2error) || !is.finite(l2error) || length(l2error) != 1) { l2error = NA; }
    if (!is.numeric(Type1) || is.na(Type1) || !is.finite(Type1) || length(Type1) != 1) { Type1 = NA; }
    if (!is.numeric(Type2) || is.na(Type2) || !is.finite(Type2) || length(Type2) != 1) { Type2 = NA; }
    if (!is.numeric(L.zero) || is.na(L.zero) || !is.finite(L.zero) || length(L.zero) != 1) { L.zero = NA; }        
   ## SimVec = c(l1error, l2error, Type1, Type2, L.zero, as.numeric(time[1:3]));
   MyV = c(l1error, l2error, Type1, Type2, L.zero, time);
   if (!is.numeric(IFunction) || IFunction > length(NameFunctions) ||
     IFunction <= 0) {
     AFilePrint(paste("Oh No, IFunction = ", IFunction, " this won't work.!", sep=""));  
     flush.console();
   }
   
   NoNoiseBetaStart <- NULL;
   TemperatureList <- NULL;
   if (!is.null(limlas$NoNoiseBetaStart)) {
     NoNoiseBetaStart <- limlas$NoNoiseBetaStart;
   }
   if (!is.null(limlas$TemperatureList)) {
     TemperatureList <- limlas$TemperatureList;
   } 
   eval(parse(text=GetG0Text("CurrentSmallContainDir", "globalenv()", S=1)));
   eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));

   setwd(CurrentSmallContainDir);
   try(dir.create(paste("Processes/", ZZMyID, sep=""), showWarnings=FALSE, recursive=TRUE)); 
    NowToSaveTime <- TwoSimR5:::GetRealTime();
    FinishesTime <- NowToSaveTime;
   SimulationIdentifier=IDUS; 
   RunProcessIdentifier <- UniqueProcessIdentifier; 
   RunProcessIdentifer <- RunProcessIdentifier; 
   NameFunction <- FunctionName;
   eval(parse(text=SetGText("MyV", "globalenv()", S=1)));
   eval(parse(text=SetGText("RunProcessIdentifer", "globalenv()", S=1)));
   eval(parse(text=SetGText("IDUS", "globalenv()", S=1)));
   FileSavedMyVName <- paste("Processes//", ZZMyID, 
     "//S", FunctionName, "S", IDUS, ".RData", sep="");
   eval(parse(text=SetGText("FileSavedMyVName", "globalenv()", S=1)));
   eval(parse(text=GetG0Text("sigma", "globalenv()", S=1)));
   SigmaNoise <- sigma;
   eval(parse(text=GetG0Text("CovRho", "globalenv()", S=1)));
   UseRho <- CovRho;
   try(save(MyV=MyV, RunProcessIdentifier = RunProcessIdentifier, 
     RunProcessIdentifer = RunProcessIdentifer, NowToSaveTime = NowToSaveTime,
     MyID = MyID,  SimulationIdentifier = SimulationIdentifier,
     UniqueSimulationIdentifier=IDUS, NameFunction=NameFunction,SigmaNoise = SigmaNoise,
     sigma = sigma, CovRho = CovRho, UseRho=UseRho,
     FunctionName=FunctionName, ZZMyID=ZZMyID,  FinishesTime=FinishesTime, StartsTime=StartsTime,
     TemperatureList=TemperatureList,
     NoNoiseBetaStart=NoNoiseBetaStart,
     file=paste("Processes//", ZZMyID, 
     "//S", FunctionName, "S", IDUS, ".RData", sep="")));
   if (FALSE) {
     BackMyV <- MyV;  MyV <- NULL;
     try(load(paste("Processes//", ZZMyID, 
     "//S", FunctionName, "S", IDUS, ".RData", sep="")))
   }  
     
     



   eval(parse(text=GetG0Text("TotalThralls", S=1)));
   if (TotalThralls <= 0) { TotalThralls <= MyID; }
   MyL <- unlist(list.files(paste("Processes//", ZZMyID, sep="")));
   if (any(MyL == "MyProcesses.RData")) {
     load(paste("Processes//", ZZMyID, "//MyProcesses.RData", sep=""));
     ##ALD <- length(ProgressData);
     if (!exists("PrintOnOne") ||  PrintOnOne > length(ProgressNames) ||
       !(ProgressNames[[PrintOnOne]][1] %in%  c(SMS$UniqueSimulationIdentifier,
         paste("SMS", SMS$UniqueSimulationIdentifier, sep=""))) ||
       !(ProgressNames[[PrintOnOne]][2] == FunctionName)) {
       ATK <- -1;
       for (KK in 1:length(ProgressNames)) {
         if (ProgressNames[[KK]][1] %in%  c(SMS$UniqueSimulationIdentifier,
         paste("SMS", SMS$UniqueSimulationIdentifier, sep="")) &&
           ProgressNames[[KK]][2] == FunctionName
         ) {
            ATK <- KK;  PrintOnOne <- KK;
            break; 
         }
       }
       if (ATK <= 0) {
          PrintOnOne <- length(ProgressNames)+1;
          ProgressNames[[PrintOnOne]] <- c(SMS$UniqueSimulationIdentifier, FunctionName);
          OurTime <- GetRealTime();
          ProgressData[[PrintOnOne]] <-  c( 
            StartDays=OurTime$JustDays, StartSeconds=OurTime$JustSeconds,
            AttemptSucceedFail = -1, ProcessTime1=0, ProcessTime2=0, ProcessTime3=0,
            L2Accuracy=0.0, FinishDays=0, FinishSeconds=0)
       }
     }
     SimulationIdentifier <- SMS$UniqueSimulationIdentifier;
     if (as.character(ProgressNames[[PrintOnOne]][1]) != as.character(SimulationIdentifier)) {
        AFilePrint("Weird, we loaded Progress but didn't get SMS US back!");
     }
     if (as.character(ProgressNames[[PrintOnOne]][2]) != as.character(FunctionName)) {
        AFilePrint("Weird, we loaded Progress but didn't get SMS US back!");
     }
     FinishTime = GetRealTime(); 
     if (is.na(l2error)) {
       ProgressData[[PrintOnOne]] <- c(ProgressData[[PrintOnOne]][1:2], -1,-1,-1,-1,-1, FinishTime$JustDays,
         FinishTime$JustSeconds); 
     }  else {
       ProgressData[[PrintOnOne]] <- c(ProgressData[[PrintOnOne]][1:2], 1,time[1:3],l2error, FinishTime$JustDays,
         FinishTime$JustSeconds);         
     }
     save(ProgressNames=ProgressNames, ProgressData=ProgressData, 
      RunProcessIdentifier=RunProcessIdentifier,
      file=paste("Processes//", ZZMyID, "//MyProcesses.RData", sep=""));
   } else {
      AFilePrint("Hey Weird Error there was no processes file in Table!");
   }
   
   NowToSaveTime <- TwoSimR5:::GetRealTime();
   if (!exists("ProgressNames") || is.null(ProgressNames)) {
    ProgressNames <- list();
    SimulationIdentifier <- SMS$UniqueSimulationIdentifier;
    ProgressNames[[1]] <- c(SimulationIdentifier=SimulationIdentifier, FunctionName=FunctionName);
    ProgressData <- list();
    ProgressData[[1]] <-c( 
      StartDays=OurTime$JustDays, StartSeconds=OurTime$JustSeconds,
      AttemptSucceedFail = -1, ProcessTime1=0, ProcessTime2=0, ProcessTime3=0,
      L2Accuracy=MyV[2], FinishDays=0, FinishSeconds=0)
     save(ProgressNames=ProgressNames, ProgressData=ProgressData, 
      RunProcessIdentifier= RunProcessIdentifier,
      file=paste("Processes//", ZZMyID, "//MyProcesses.RData", sep=""));
   };
      
   MyProcessCIFit <- NULL;
   if (is.logical(DoCI) && DoCI == TRUE  && !is.null(FittedOb$CIEst) ) {
     RunProcessIdentifier <- UniqueProcessIdentifier;
     SimulationIdentifier <- SMS$UniqueSimulationIdentifier;
     try(MyProcessCIFit <- ProcessAssessCIFit(SMS = SMS, FittedOb = FittedOb,
     RunProcessIdentifier = RunProcessIdentifier,
     UniqueSimulationIdentifier = SimulationIdentifier, 
     NameFunction = FunctionName, 
     Hit = Hit,
     quoteMore = paste(quoteMore, " - ", sep=""), 
     verbose = max(verbose -1, verbose), AlreadyLocked=FALSE)); 
   }
   MyProcessBetaCIFit <- NULL;
   if (is.logical(DoCI) && DoCI == TRUE  && !is.null(FittedOb$CIEst) ) {
     RunProcessIdentifier <- UniqueProcessIdentifier;
     SimulationIdentifier <- SMS$UniqueSimulationIdentifier;
     try(MyProcessBetaCIFit <- ProcessAssessCIFitAllNonZeroBeta(SMS = SMS, FittedOb = FittedOb,
     RunProcessIdentifier = RunProcessIdentifier,
     UniqueSimulationIdentifier = SimulationIdentifier, 
     NameFunction = FunctionName, 
     Hit = Hit,
     quoteMore = paste(quoteMore, " - ", sep=""), 
     verbose = max(verbose -1, verbose), AlreadyLocked=FALSE)); 
   }
   
   AssesMIP <- NULL;
   if (!is.null(FittedOb$MIPReport) ) {
     SimulationIdentifier <- SMS$UniqueSimulationIdentifier
     try(AssesMIP <- ProcessAssessMIPFit(SMS = SMS, FittedOb = FittedOb,
     RunProcessIdentifier = RunProcessIdentifier,
     UniqueSimulationIdentifier = SimulationIdentifier, 
     NameFunction = FunctionName, 
     Hit = Hit,
     quoteMore = paste(quoteMore, " - ", sep=""), 
     verbose = max(verbose -1, verbose), AlreadyLocked=FALSE)); 
   }
   
   ## Nothing to unlock here
   ##UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MDir);
   if (verbose >= 3) {
     AFilePrint("RunOneSimExperiment:  All Finished, return(1)!"); flush.console();
   } 
   try(setwd(CurrentSmallContainDir)); 
   
   if (any(names(limlas) == "GibbsSS")) {
        DSaveDir <- limlas$GibbsSS$TBSR5$SaveDir;
     try(print(paste("Going to Delete GibbsSS BS directory: ",
       limlas$GibbsSS$TBSR5$SaveDir, sep=""))); flush.console();
     try(limlas$GibbsSS$TBSR5$unsave(TRUE), silent=TRUE);
     try(unlink(DSaveDir));
     try(limlas$GibbsSS <- NULL);
   }
   if (any(names(limlas) == "BackGibbsSS")) {
     DSaveDir <- limlas$BackGibbsSS$TBSR5$SaveDir;
     try(print(paste("Going to Delete BackGibbsSS BS directory: ",
       limlas$BackGibbsSS$TBSR5$SaveDir, sep=""))); flush.console();
     try(limlas$BackGibbsSS$TBSR5$unsave(TRUE), silent=TRUE);
     try(unlink(DSaveDir));
     try(limlas$BackGibbsSS <- NULL);
   }
   if (any(names(limlas) == "NGibbsSS")) {
     DSaveDir <- limlas$NGibbsSS$TBSR5$SaveDir;
     try(print(paste("Going to Delete NGibbsSS BS directory: ",
       limlas$NGibbsSS$TBSR5$SaveDir, sep=""))); flush.console();
     try(limlas$NGibbsSS$TBSR5$unsave(TRUE), silent=TRUE);
     try(unlink(DSaveDir));
     try(limlas$NGibbsSS <- NULL);
   }
   try(setwd(OriginalOldwd));
   return(1);
};

SetBackToMandatory <- function() {
  eval(parse(text=GetG0Text("MandatoryOnFunction", "globalenv()", S=1)));  
  eval(parse(text=GetG0Text("SubMitList", "globalenv()", S=1)));
  print(paste("SetBackToMandatory: Running we have MandatoryOnFunction is ", MandatoryOnFunction, sep=""));
  flush.console();
  MyGoods <- (1:NROW(SubMitList))[(SubMitList[,1] == MandatoryOnFunction |
     SubMitList[,2] == MandatoryOnFunction |
     SubMitList[,1] == paste("Generate", MandatoryOnFunction, sep="") |
     SubMitList[,2] == paste("Generate", MandatoryOnFunction, sep="") |
     SubMitList[,1] == paste("Generate", MandatoryOnFunction, "NI", sep="") |
     SubMitList[,2] == paste("Generate", MandatoryOnFunction, "NI", sep="") |
     SubMitList[,2] == paste(MandatoryOnFunction, "NI", sep="")
     ) &
     as.character(SubMitList[,3]) %in% c("0", 0, " 0", "0 ", " 0 ", "1", "1 ", " 1")]; 
  if (length(MyGoods) == 0) {
    print(paste("Error SetBackToMandatory we can't find in SubMitList an Mandatory = ", MandatoryOnFunction, sep=""));
    print(paste("SubMit List is: "));
    print(SubMitList);
    flush.console();
    print("Error In Run Process!"); flush.console();
    quit();
  } 
  ABB <- MyGoods[1];
  eval(parse(text=SetGText("ABB", "globalenv()", S=1)));
  SMSName <- SubMitList[ABB,1];
  FunctionName  <- SubMitList[ABB,2];
  eval(parse(text=SetGText("SMSName", "globalenv()", S=1)));
  eval(parse(text=SetGText("FunctionName", "globalenv()", S=1)));
}
InitiateEvaluations <- function() {
  return("
  eval(parse(text=GetG0Text(\"OriginalOldwd\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"UniqueProcessIdentifier\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"MandatoryOnFunction\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"NameFunctions\", \"globalenv()\", S=1)));
  if (any(SubMitList[,2] == \"GroupBayesSpike\")) {
    SubMitList[SubMitList[,2] == \"GroupBayesSpike\",2] <- \"GenerateGroupBayesSpike\";
  }
  if (any(SubMitList[,2] == \"GroupTwoLassoSpread\")) {
    SubMitList[SubMitList[,2] == \"GroupTwoLassoSpread\",2] <- \"GenerateGroupTwoLassoSpread\";
  }               
  eval(parse(text=GetG0Text(\"CrashDir\", \"globalenv()\", S=1)));
  if (!is.null(MandatoryOnFunction) && is.character(MandatoryOnFunction) &&
    MandatoryOnFunction != \"\" && MandatoryOnFunction %in% NameFunctions) {
    SetBackToMandatory();
    eval(parse(text=GetG0Text(\"FunctionName\", \"globalenv()\", S=1)));    
    eval(parse(text=GetG0Text(\"ABB\", \"globalenv()\", S=1)));  
    eval(parse(text=GetG0Text(\"SMSName\", \"globalenv()\", S=1)));  
    eval(parse(text=GetG0Text(\"SubMitList\", \"globalenv()\", S=1)));  
  }
  DefineDefaultUseFunction();
  RunProcessIdentifier <- UniqueProcessIdentifier;
  ##eval(parse(text=GetG0Text(\"ListOfParams\")));  
  if (substr(SMSName, nchar(SMSName)-nchar(\".RData\")+1, nchar(SMSName)) == \".RData\") {
    SMSName <- substr(SMSName, 1,nchar(SMSName)-nchar(\".RData\"));
  }
  if (substr(SMSName, 1, nchar(\"SMS\")) == \"SMS\") {
    NoSMSSMSName <- substr(SMSName, nchar(\"SMS\")+1, nchar(SMSName));
    SMSName <- SMSName;
  } else {
    NoSMSSMSName <- SMSName;
    SMSName <- paste(\"SMS\", SMSName, sep=\"\");    
  }
  try(setwd(OriginalOldwd));
  try(eval(parse(text=AOverallPaste(\"ListOfParams\"))), silent=TRUE);
  ");
}
TestForVerbose <- function() {
  return("
    #Use custom version of Alan's SimMeData
  if (!exists(\"FunctionName\")) { 
     print(\"No: You've got to supply a function!\");
     return(-1);
  }
  if (exists(\"verbose\") && !is.null(verbose) && is.numeric(verbose) &&
	  length(verbose) == 1 && verbose == -1) {
	  MyT <- \"verbose <- verbose;\"
	  try(eval(parse(text=MyT)), silent=TRUE);
  } else if (!exists(\"verbose\") || is.null(verbose)  || !is.numeric(verbose)) {
    eval(parse(text=GetG0Text(\"verbose\")));
    if (is.null(verbose) || !is.numeric(verbose) || verbose[1] == 0) {
      MyT <- \"verbose <- verbose\";
      try(eval(parse(text=MyT)), silent=TRUE);
    } else if (is.logical(verbose)) {
      if (verbose == TRUE) {
        if (verbose >= 1) {
        } else {
          try(verbose <- as.integer(1));
        }
      }  else {
        try(verbose <- as.integer(0));
      }
    } else if (is.numeric(verbose) && verbose[1] == -1) {
      MyT <- \"verbose <- verbose\";
      try(eval(parse(text=MyT)), silent=TRUE);
    } else {
      try(verbose <- verbose);
    }
  } else {
    MyT <- \"verbose <- verbose\";
    try(eval(parse(text=MyT)), silent=TRUE);
  }
  ");
}
LoadSMSNow <- function() {
  return("
   eval(parse(text=GetG0Text(\"CurrentLargeContainDir\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"CurrentSmallContainDir\", \"globalenv()\", S=1)));
  setwd(CurrentLargeContainDir);
  MyLLS <- unlist(list.files(\"SMS\"));
  if (!(SMSName %in% MyLLS) && !(NoSMSSMSName %in% MyLLS)) {
    AErrorPrint(\"Error: No we can't find the SMS you are looking for!\"); flush.console();
    return(-1);
  }

  SMS <- NULL;
  setwd(CurrentSmallContainDir);
  try(dir.create(\"SMS\", showWarnings=FALSE, recursive=TRUE));
  try(dir.create(paste(\"SMS/\", \"SMS\", NoSMSSMSName, sep=\"\"), showWarnings=FALSE, recursive=TRUE));
  TooManyFail <- 10;  CountFails <- 0;
  while( 
  TwoSimR5:::LockMeIn(verbose=as.numeric(verbose), 
    quoteMore=\"Reading in an SMS\", LFDir = paste(\"SMS/SMS\", NoSMSSMSName, sep=\"\")) ==FALSE 
    && TooManyFail > CountFails
    ) {
    Sys.sleep(runif(1,0,4));  
    CountFails <- CountFails+1;    
  }
  setwd(CurrentLargeContainDir);
  try(load(paste(\"SMS/SMS\", NoSMSSMSName, \"//\", NoSMSSMSName, \".RData\", sep=\"\")));
  setwd(CurrentSmallContainDir)
  try(TwoSimR5:::UnLockMe(LFDir = paste(\"SMS/\", SMSName, sep=\"\"), 
       verbose=verbose, quoteMore=\"Loaded an SMS\"));
  if (is.null(SMS)) {
    AFilePrint(\"Error We tried to load SMS but got NULL SMS!\");
    AFilePrint(paste(\"NoSMSSMSName = \", NoSMSSMSName, sep=\"\"));
    return(-1);
  }
  ");
}

FunctionUnderstandProcesses <- function() {
  return("
   eval(parse(text=GetG0Text(\"MyID\", \"globalenv()\", S=1)));
   eval(parse(text=GetG0Text(\"TotalThralls\", S=1)));
   if (TotalThralls <= 0) { TotalThralls <= MyID; }
   ZZMyID <- AZeroOut(as.numeric(MyID),10^ceiling(log(TotalThralls+1,10)))
   MyL <- unlist(list.files(paste(\"Processes//\", ZZMyID, \"\", sep=\"\")));
       
  dir.create(paste(\"Processes//\", ZZMyID, sep=\"\"), showWarnings=FALSE, recursive=TRUE);
  MyL <- unlist(list.files(paste(\"Processes//\", ZZMyID, sep=\"\")));
  if (any(MyL == \"MyProcesses.RData\")) {
    ProgressNames <- NULL;  ProgressData <- NULL;
    BackNameFunction <- FunctionName;  BackSMSName <- SMSName;
    try(load(paste(\"Processes//\", ZZMyID, \"//MyProcesses.RData\", sep=\"\")));
    if (is.null(ProgressNames)) {
      try(eval(parse(text=RecoverProgress(MyL))));
    }
    ##ProgressNames=ProgressNames, ProgressData=ProgressData, 
  }
  OurTime = GetRealTime();  StartsTime <- OurTime;
   try(dir.create(\"Processes\", showWarnings=FALSE, recursive=TRUE));
   try(dir.create(paste(\"Processes//\", ZZMyID, sep=\"\"), showWarnings=FALSE, recursive=TRUE));

   IDUS <- SMS$UniqueSimulationIdentifier;
  ");
}

EvaluateProgress <- function() {
  return("
    if (!exists(\"ProgressNames\") || is.null(ProgressNames)) {
    ProgressNames <- list();
    SimulationIdentifier   <- SMS$UniqueSimulationIdentifier
    ProgressNames[[1]] <- c(SimulationIdentifier=SimulationIdentifier, FunctionName=FunctionName);
    ProgressData <- list();
    ProgressData[[1]] <-c( 
      StartDays=OurTime$JustDays, StartSeconds=OurTime$JustSeconds,
      AttemptSucceedFail = -1, ProcessTime1=0, ProcessTime2=0, ProcessTime3=0,
      L2Accuracy=0.0, FinishDays=0, FinishSeconds=0)
    PrintOnOne <- 1;
  } else {
     ATK <- -1
     SimulationIdentifier  <- SMS$UniqueSimulationIdentifier
     for (KK in 1:length(ProgressNames)) {
       if (ProgressNames[[KK]][1] == SimulationIdentifier  &&
         ProgressNames[[KK]][1] == FunctionName) {
         ATK = KK; break;   
       }
     }
     if (ATK >= 1) {
       AFilePrint(\"Hey: We are writing over a previous attempt\");  
       PrintOnOne <- ATK;
     } else {     
       PrintOnOne <- length(ProgressNames)+1;
     }
       ProgressNames[[PrintOnOne]] <- 
         c(SimulationIdentifier=SimulationIdentifier,
           FunctionName=FunctionName);
       ProgressData[[PrintOnOne]] <- 
         c(StartDays=OurTime$JustDays, StartSeconds=OurTime$JustSeconds,
         AttemptSucceedFail = -1, ProcessTime1=0, ProcessTime2=0, ProcessTime3=0,
         L2Accuracy=0.0, FinishDays=0, FinishSeconds=0);

  }
  setwd(CurrentSmallContainDir);
  save(ProgressNames=ProgressNames, ProgressData=ProgressData, 
    file=paste(\"Processes//\", ZZMyID, \"//MyProcesses.RData\", sep=\"\"));
  ");
}
EstablishFunctionAndIdentifiers <- function() {
  return("
  eval(parse(text=GetG0Text(\"UseFunctions\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"NameFunctions\", \"globalenv()\", S=1)));
  MyNameFunction = FunctionName;
  NameFunction <- FunctionName;
  IFunction = (1:length(NameFunctions))[NameFunctions == MyNameFunction];
  eval(parse(text=SetGText(\"IFunction\", S=1)));
  eval(parse(text=SetGText(\"MyNameFunction\", S=1)));
  eval(parse(text=SetGText(\"SMS\", S=1)));  
  SimulationIdentifier <- SMS$UniqueSimulationIdentifier
  if (is.null(SimulationIdentifier)) {
    AFilePrint(\"RunExperimentRun:  Error: SMS$UniqueSimulationIdentifier is NULL! \"); flush.console();
    eval(parse(text=GetG0Text(\"OriginalOldwd\", \"globalenv()\", S=1)));
    try(setwd(OriginalOldwd));
    return(-666);
  }
  if (verbose > 0) {
    AFilePrint(paste(\"We're going to run function: \", NameFunctions[IFunction], sep=\"\"));
    flush.console();
  }
  ");
}


ProcessAssessMIPFit = function(SMS, FittedOb,
     RunProcessIdentifier,
     UniqueSimulationIdentifier, 
     NameFunction, 
     Hit = "+",
     quoteMore = " - ", 
     verbose = 0, DontRecord=FALSE, AlreadyLocked=FALSE) {
     if (is.null(FittedOb$MIPReport)) {
       AFilePrint(paste("AssessMIPFit: Hey, we tried to do MIP fit for ", NameFunction, " but got NULL!"));
     } 
     MIPReport <- FittedOb$MIPReport;
     MyTestT <- " 
     if (is.logical(verbose) && verbose == FALSE) {
       verbose = FALSE;
     } else if (is.logical(verbose) && verbose == TRUE) { verbose = 2; 
     } else if (!is.numeric(verbose)) {
       printf(\"AssessMIPFit: Hey, this isn't a good verbose;\"); verbose = 2;
     } else {
       verbose = verbose-2;
     } ";
     try(eval(parse(text=MyTestT)));
     if (verbose >= 1) {
       AFilePrint(paste("AssessMIPFit; Starting for ", NameFunction, ", Hit = ", Hit, sep=""));
       flush.console();
    }
    if (is.null(FittedOb$MIPReport) || length(FittedOb$MIPReport) <= 0) {
      AFilePrint(paste("AssessMIPFit: No MIPReport Given for ", NameFunction, sep="")); flush.console();
      return(-1);
    }
    RTIF <- NULL;
    Type1 = NULL; Type2=NULL; SphereLoss = NULL; AUC=NULL; MarginalHellinger=NULL;
    try(Type1 <- Type1Function(FittedOb$MIPReport, SMS$BetasReal));
    try(Type2 <- Type2Function(FittedOb$MIPReport, SMS$BetasReal)); 
    try(RTIF <- RankTotalIntegrateFunction(FittedOb$MIPReport, SMS$BetasReal));
    try(SphereLoss <- SphereLossFunction(FittedOb$MIPReport, SMS$BetasReal));
    try(AUC <- AUCFunction(FittedOb$MIPReport, SMS$BetasReal));
    try(MarginalHellinger <- MarginalHellingerFunction(FittedOb$MIPReport, SMS$BetasReal));
    
    uniqueBetas <- sort(unique(SMS$BetasReal));
    iMIPReport <- rep(0, length(uniqueBetas));
    for (ii in 1:length(uniqueBetas)) {
      idx <- (1:length(MIPReport))[SMS$BetasReal == uniqueBetas[ii]];
      if (length(idx)==1) {
        iMIPReport[ii] <- MIPReport[idx];
      } else if (length(idx) > 1) {
        iMIPReport[ii] <- mean(MIPReport[idx]);
      }
    }
    names(iMIPReport) <- uniqueBetas;
    RTIFRandom = NULL;
    Type1Random=NULL; Type2Random=NULL; SphereLossRandom=NULL;
    AUCRandom=NULL; MarginalHellingerRandom=NULL; 
    if (!is.null(SMS$TrueOnGroups) && length(SMS$TrueOnGroups) >= 1) {
      if (!is.null(FittedOb$MIPRandom) && length(FittedOb$MIPRandom) == length(SMS$TrueOnGroups)) {
      try(Type1Random <- Type1Function(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(Type2Random <- Type2Function(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(SphereLossRandom <- SphereLossFunction(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(AUCRandom <- AUCFunction(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(RTIFRandom <- RankTotalIntegrateFunction(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(MarginalHellingerRandom <- MarginalHellingerFunction(FittedOb$MIPRandom, SMS$TrueOnGroups));
      }
    }
  
  try(IndexOnBetas <- (1:length(SMS$BetasReal))[SMS$BetasReal!=0.0]);
  try(IndexOffBetas <- (1:length(SMS$BetasReal))[SMS$BetasReal==0.0]);
  try(BetasOn <- SMS$BetasReal[IndexOnBetas]);
  try(FitBetaOn <- FittedOb$BetaFit[IndexOnBetas]);
  try(MIPBetaOn <- FittedOb$MIPReport[IndexOnBetas]);
  try(MeanMIPBetaOff <- mean(FittedOb$MIPReport[IndexOffBetas])); 
  MaxMIPBetaOff <- NULL;
  try(MaxMIPBetaOff <- max(FittedOb$MIPReport[IndexOffBetas])); 
  NoNoiseBetaStart <- FALSE;
  try(NoNoiseBetaStart <- FittedOb$NoNoiseBetaStart);
  TemperatureList <- NULL; 
  if (is.null(FittedOb$BetaStart)) {
     BetaStartOn <- NULL; IndexOnBetaStart <- NULL;
  } else {
    IndexOnBetaStart <- (1:length(FittedOb$BetaStart))[FittedOb$BetaStart != 0];
    BetaStartOn <- FittedOb$BetaStart[IndexOnBetaStart];
  } 
  if (!is.null(FittedOb$TemperatureList)) {
    TemperatureList <- FittedOb$TemperatureList;
  }
  AllMIPReturn <- list(Type1=Type1, Type2=Type2, SphereLoss=SphereLoss,
  AUC=AUC, MarginalHellinger=MarginalHellinger, Type1Random=Type1Random,
  Type2Random=Type2Random, SphereLossRandom=SphereLossRandom, RTIF=RTIF, RTIFRandom=RTIFRandom,
  AUCRandom=AUCRandom,MarginalHellingerRandom = MarginalHellingerRandom,
  iMIPReport=iMIPReport,
  IndexOnBetas = IndexOnBetas, 
  BetasOn=BetasOn, MIPBetaOn = MIPBetaOn, MeanMIPBetaOff = MeanMIPBetaOff, 
  MaxMIPBetaOff = MaxMIPBetaOff,
  TemperatureList=TemperatureList, IndexOnBetaStart=IndexOnBetaStart, BetaStartOn=BetaStartOn,
  NoNoiseBetaStart=NoNoiseBetaStart, TemperatureList=TemperatureList
  )
  if (exists("DontRecord") && is.logical(DontRecord) && DontRecord == TRUE) {
    return(AllMIPReturn);
  }
  eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
  if (verbose >= 2) {
      AFilePrint(paste("----  AssessMIPFit[", UniqueSimulationIdentifier, 
        "] --- CurrentSmallContainDir = ", CurrentSmallContainDir, " for NameFunction = ", NameFunction, sep=""));
      flush.console();
  }
 BADUniqueProcessIdentifier <- 1;
  ATryGo <- "
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  }
  if (length(UniqueProcessIdentifier) > 1) {
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1]);
  }
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  } else if (is.null(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NULL\");
    try(UniqueProcessIdentifier <- \"\"); 
  } else if (length(UniqueProcessIdentifier) > 1) {
    UniqueProcessIdentifier <- UniqueProcessIdentifier[1]; 
  } else if (is.na(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NA\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.numeric(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NUMERIC\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.character(UniqueProcessIdentifier)) {
    if (UniqueProcessIdentifier[1] == \"\")  {
      AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is Blank!\");
    } else {
      BADUniqueProcessIdentifier <- 0;
    }
  } else {
    try(UniqueProcessIdentifier <- \"\");
  }
  ";
  try(eval(parse(text=ATryGo)))
  if (BADUniqueProcessIdentifier == 1) {
    AFilePrint("AssessMIPFit: UnqiueProcessIdentifier is NULL!"); flush.console();  
  }
  if (verbose >= 2) {
    AFilePrint(paste("  --- AssedMIPFit, ", NameFunction, " looking for HMS. ", sep="")); flush.console();
  }
  eval(parse(text=GetG0Text("CurrentSmallContainDir", "globalenv()", S=1)));
  try(setwd(CurrentSmallContainDir));
  eval(parse(text=GetG0Text("ZZMyID")));
  SimulationIdentifier <- SMS$UniqueSimulationIdentifier;
  FunctionName <- NameFunction;
  TemperatureList<- NULL; BetaStartOn <- NULL;  IndexOnBetaStart = NULL;  NoNoiseBetaStart <- NULL;
  if (!is.null(FittedOb$BetaStart)) {
    IndexOnBetaStart <- (1:length(FittedOb$BetaStart))[FittedOb$BetaStart != 0];
    BetaStartOn <-   FittedOb$BetaStart[IndexOnBetaStart];
    NoNoiseBetaStart <- TRUE;
  }
  if (!is.null(FittedOb$TemperatureList)) {
    TemperatureList <- FittedOb$TemperatureList;
  }
  save(AllMIPReturn=AllMIPReturn, NameFunction=NameFunction,
     FunctionName = FunctionName, RunProcessIdentifier=RunProcessIdentifier,
     Hit=Hit, SimulationIdentifier = SimulationIdentifier,
     TemperatureList=TemperatureList, 
     file=paste("Processes//", ZZMyID, 
     "//MIP", NameFunction, "S", SimulationIdentifier, ".RData", 
     sep=""));
  
  return(AllMIPReturn);
};

ProcessAssessCIFit = function(SMS, FittedOb,
     RunProcessIdentifier,
     UniqueSimulationIdentifier, 
     NameFunction, 
     Hit = "+",
     quoteMore = " - ", 
     verbose = 0, DontRecord=FALSE, AlreadyLocked=FALSE) {
     if (is.null(FittedOb$CIEst)) {
       AFilePrint(paste("AssessCIFit: Hey, we tried to do CI fit for ", NameFunction, " but got NULL!"));
     } 
     MyTestT <- " 
     if (is.logical(verbose) && verbose == FALSE) {
       verbose = FALSE;
     } else if (is.logical(verbose) && verbose == TRUE) { verbose = 2; 
     } else if (!is.numeric(verbose)) {
       printf(\"AssessCIFit: Hey, this isn't a good verbose;\"); verbose = 2;
     } else {
       verbose = verbose-2;
     } ";
     try(eval(parse(text=MyTestT)));
     if (verbose >= 1) {
       AFilePrint(paste("AssessCIFit; Starting for ", NameFunction, ", Hit = ", Hit, sep=""));
       flush.console();
     }
     Betas <- SMS$BetasReal;
     IndexOnBetas <- (1:length(Betas))[Betas!=0.0];
     IndexOffBetas <- (1:length(Betas))[Betas==0.0];
     IndexSmallBetas <- (1:length(Betas))[Betas!=0.0 &
       Betas^2 <= SMS$SigmaNoiseReal / length(SMS$NLen)];
     CIQuantiles <- FittedOb$CIQuantiles; 
     if (CIQuantiles[1] == .5) {
       AK <- 1; 
     } else { AK <- 0;}
     QuantileCov <- rep(0, length(CIQuantiles)/2);
     ArAt <- AK;
     for (jj in 1:length(QuantileCov)) {
       QuantileCov[jj] <- CIQuantiles[ArAt + 2] - CIQuantiles[ArAt+1];
       ArAt <- ArAt +2;
     }
     CITime = FittedOb$CITime;
     
     AllCIReturn <- list();
     for (tjt in 1:length(FittedOb$CIEst)) {
       if (any(CIQuantiles == .5)) {
         MedianAbsDiff <-sum(abs(Betas - FittedOb$CIEst[[tjt]][,CIQuantiles == .5]))
         MedianOnAbsDiff <- sum(abs(Betas[IndexOnBetas] - FittedOb$CIEst[[tjt]][IndexOnBetas,CIQuantiles == .5]))
         MedianOffAbsDiff <- sum(abs(Betas[IndexOffBetas] - FittedOb$CIEst[[tjt]][IndexOffBetas,CIQuantiles == .5]))
         MedianSqDiff <- sum((Betas - FittedOb$CIEst[[tjt]][,CIQuantiles == .5])^2)
         MedianOnSqDiff <- sum((Betas[IndexOnBetas] - FittedOb$CIEst[[tjt]][IndexOnBetas,CIQuantiles == .5])^2)
         MedianOffSqDiff <- sum((Betas[IndexOffBetas] - FittedOb$CIEst[[tjt]][IndexOffBetas,CIQuantiles == .5])^2)
         if (length(IndexSmallBetas) >= 1) {
           MedianSmallAbsDiff <- sum(abs(Betas[IndexSmallBetas] - FittedOb$CIEst[[tjt]][IndexSmallBetas,CIQuantiles == .5]))
           MedianSmallSqDiff <- sum((Betas[IndexSmallBetas] - FittedOb$CIEst[[tjt]][IndexSmallBetas,CIQuantiles == .5])^2)            
         } else {
           MeanSmallAbsDiff <- NULL; MeanSmallSqDiff = NULL;  
           MedianSmallAbsDiff <- NULL; MedianSmallSqDiff = NULL;         
         }          
       } else {
         MedianAbsDiff = NULL;  MedianOnAbsDiff = NULL;
         MedianOffAbsDiff = NULL;  MedianSmallAbsDiff = NULL;
         MedianSqDiff = NULL;  MedianOnSqDiff = NULL;
         MedianOffSqDiff = NULL;  MedianSmallSqDiff = NULL;
       }
       ArAt <- AK;
       CoverageQuantiles <- rep(0, length(QuantileCov));
       CoverageOnQuantiles <- rep(0, length(QuantileCov));
       CoverageOffQuantiles <- rep(0, length(QuantileCov));
       CoverageSmallQuantiles <- rep(0, length(QuantileCov));              
       MeanWidthQuantiles <- rep(0, length(QuantileCov));
       MeanOnWidthQuantiles <- rep(0, length(QuantileCov));
       MeanOffWidthQuantiles <- rep(0, length(QuantileCov));
       MeanSmallWidthQuantiles <- rep(0, length(QuantileCov));
                     
       MeanSqWidthQuantiles <- rep(0, length(QuantileCov));
       MeanSqOnWidthQuantiles <- rep(0, length(QuantileCov));
       MeanSqOffWidthQuantiles <- rep(0, length(QuantileCov));
       MeanSqSmallWidthQuantiles <- rep(0, length(QuantileCov)); 
       
       CINormalLoss <- rep(0, length(QuantileCov));
       CIUnifLoss <- rep(0, length(QuantileCov));
       CIOnNormalLoss <- rep(0, length(QuantileCov));
       CIOnUnifLoss <- rep(0, length(QuantileCov));
       CIOffNormalLoss <- rep(0, length(QuantileCov));
       CIOffUnifLoss <- rep(0, length(QuantileCov));
       CISmallNormalLoss <- rep(0, length(QuantileCov));
       CISmallUnifLoss <- rep(0, length(QuantileCov));
                          
       for (jj in 1:length(QuantileCov)) {
         CoverageQuantiles[jj] <- sum( FittedOb$CIEst[[tjt]][,ArAt+1] <=
            Betas & FittedOb$CIEst[[tjt]][,ArAt+2] >= Betas) /length(Betas);
         CoverageOnQuantiles[jj] <- sum( FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+1] <=
            Betas[IndexOnBetas] & FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+2] >= Betas[IndexOnBetas])/
            length(IndexOnBetas);
         CoverageOffQuantiles[jj] <- sum( FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+1] <=
            Betas[IndexOffBetas] & FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+2] >= Betas[IndexOffBetas])/
            length(IndexOffBetas);
         CINormal <- CINewLoss(FittedOb$CIEst[[tjt]][, ArAt+1], 
           FittedOb$CIEst[[tjt]][, ArAt+2], Betas, QuantileCov[jj]);
         CINormalLoss[jj] <- sum(CINormal);
         CIOnNormalLoss[jj] <- sum(CINormal[IndexOnBetas]);
         CIOffNormalLoss[jj] <- sum(CINormal[IndexOffBetas]);
         try(CIUnif <- CIUnifLoss(FittedOb$CIEst[[tjt]][, ArAt+1], 
           FittedOb$CIEst[[tjt]][, ArAt+2], Betas, QuantileCov[jj]));
         CIUnifLoss[jj] <- sum(CIUnif);
         CIOnUnifLoss[jj] <- sum(CIUnif[IndexOnBetas]);
         CIOffUnifLoss[jj] <- sum(CIUnif[IndexOffBetas]);
        
         MeanWidthQuantiles[jj] <- sum(mean(FittedOb$CIEst[[tjt]][,ArAt+2]-
           FittedOb$CIEst[[tjt]][,ArAt+1]));
         MeanOnWidthQuantiles[jj] <- sum(mean(FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+1]));           
         MeanOffWidthQuantiles[jj] <- sum(mean(FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+1]));
         MeanSqWidthQuantiles[jj] <- sum(mean((FittedOb$CIEst[[tjt]][,ArAt+2]-
           FittedOb$CIEst[[tjt]][,ArAt+1])^2));
         MeanSqOnWidthQuantiles[jj] <- sum(mean((FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+1])^2));           
         MeanSqOffWidthQuantiles[jj] <- sum(mean((FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+1])^2));           


         if (length(IndexSmallBetas) >= 1) {
           CoverageSmallQuantiles[jj] <- sum( FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+1] <=
            Betas[IndexSmallBetas] & FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+2] >= Betas[IndexSmallBetas])/
            length(IndexSmallBetas);
           MeanSmallWidthQuantiles[jj] <- sum(mean(FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+1]));                 
           MeanSqSmallWidthQuantiles[jj] <- sum(mean((FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+1])^2)); 
           try(CISmallNormalLoss[jj] <- sum(CINormal[IndexSmallBetas]));
           try(CISmallUnifLoss[jj] <- sum(CIUnif[IndexSmallBetas]));

         } else {
           CoverageSmallQuantiles <- NULL; MeanSqSmallWidthQuantiles = NULL;
           MeanSmallWidthQuantiles = NULL;    
           CISmallNormalLoss <- NULL;  CISmallUnifLoss <- NULL;      
         } 
         ArAt <- ArAt+2;           
       }
       APiece <- list(
        NameOfInterval = names(FittedOb$CIEst)[tjt],
        NameFunction = NameFunction, Hit=Hit,
        CIQuantiles = FittedOb$CIQuantiles,
        MedianAbsDiff = MedianAbsDiff,
        MedianOnAbsDiff = MedianOnAbsDiff,
        MedianOffAbsDiff = MedianOffAbsDiff, MedianSqDiff = MedianSqDiff,
        MedianOnSqDiff = MedianOnSqDiff, MedianOffSqDiff = MedianOffSqDiff,
        MedianSmallAbsDiff = MedianSmallAbsDiff,  MedianSmallSqDiff = MedianSmallSqDiff,   
         CoverageQuantiles = CoverageQuantiles,
         MeanWidthQuantiles = MeanWidthQuantiles,
         MeanSqWidthQuantiles = MeanSqWidthQuantiles,
         CoverageOnQuantiles =CoverageOnQuantiles,
         MeanOnWidthQuantiles = MeanOnWidthQuantiles,
         MeanSqOnWidthQuantiles = MeanSqOnWidthQuantiles,
         CoverageOffQuantiles = CoverageOffQuantiles,
         MeanOffWidthQuantiles = MeanOffWidthQuantiles,
         MeanSqOffWidthQuantiles = MeanSqOffWidthQuantiles,
         CoverageSmallQuantiles = CoverageSmallQuantiles,
         MeanSmallWidthQuantiles = MeanSmallWidthQuantiles,
         MeanSqSmallWidthQuantiles = MeanSqSmallWidthQuantiles,
         CINormalLoss =  CINormalLoss,
         CIUnifLoss = CIUnifLoss,
         CIOnNormalLoss =  CIOnNormalLoss,
         CIOnUnifLoss = CIOnUnifLoss,
         CIOffNormalLoss = CIOffNormalLoss, 
         CIOffUnifLoss = CIOffUnifLoss,
         CISmallNormalLoss = CISmallNormalLoss,
         CISmallUnifLoss = CISmallUnifLoss);
       AllCIReturn[[tjt]] <- APiece;
   }
   
  if (is.logical(DontRecord) && DontRecord == TRUE) {
    return(AllCIReturn);
  }
  eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
  if (verbose >= 2) {
      AFilePrint(paste("----  AssessCIFit[", UniqueSimulationIdentifier, 
        "] --- CurrentSmallContainDir = ", CurrentSmallContainDir, " for NameFunction = ", NameFunction, sep=""));
      flush.console();
  }
   BADUniqueProcessIdentifier <- 1;
  ATryGo <- "
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  }
  if (length(UniqueProcessIdentifier) > 1) {
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1]);
  }
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  } else if (is.null(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NULL\");
    try(UniqueProcessIdentifier <- \"\"); 
  } else if (length(UniqueProcessIdentifier) > 1) {
    UniqueProcessIdentifier <- UniqueProcessIdentifier[1]; 
  } else if (is.na(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NA\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.numeric(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NUMERIC\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.character(UniqueProcessIdentifier)) {
    if (UniqueProcessIdentifier[1] == \"\")  {
      AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is Blank!\");
    } else {
      BADUniqueProcessIdentifier <- 0;
    }
  } else {
    try(UniqueProcessIdentifier <- \"\");
  }
  ";
  try(eval(parse(text=ATryGo)))
  if (BADUniqueProcessIdentifier == 1) {
    AFilePrint("AssessCIFit: UnqiueProcessIdentifier is NULL!"); flush.console();  
  }
  if (verbose >= 2) {
    AFilePrint(paste("  --- AssedCIFit, ", NameFunction, " looking for HMS. ", sep="")); flush.console();
  }
  CIEst = NULL;
  if (!exists("FittedOb") || is.null(FittedOb)) {
    AFilePrint("Why did FittedOb AssesCI? that not happen what?"); flush.console(); 
  } else if (is.null(FittedOb$CIEst)) {
    AFilePrint("What: Why is FittedOb$CIEst null?"); flush.console();
  } else {
    CIEst <- FittedOb$CIEst;
  }
  eval(parse(text=GetG0Text("CurrentSmallContainDir", "globalenv()", S=1)));
  try(setwd(CurrentSmallContainDir));
  eval(parse(text=GetG0Text("ZZMyID")));
  SimulationIdentifier   <- SMS$UniqueSimulationIdentifier;
  eval(parse(text=GetG0Text("UniqueProcessIdentifier")));
  RunProcessIdentifier <- UniqueProcessIdentifier;
  FunctionName <- NameFunction;
  print(paste("About to save AllCIReturn to ",
    paste("Processes//", ZZMyID, 
     "//CI", FunctionName, "S", SimulationIdentifier, ".RData", sep=""), sep=""));
  flush.console();
  try(save(AllCIReturn=AllCIReturn, NameFunction=NameFunction, 
    FunctionName = FunctionName, Betas=Betas, CIQuantiles=CIQuantiles,
    CIEst = CIEst, Hit=Hit, CITime=CITime, RunProcessIdentifier = RunProcessIdentifier,
    SimulationIdentifier = SimulationIdentifier, file=paste("Processes//", ZZMyID, 
     "//CI", FunctionName, "S", SimulationIdentifier, ".RData", sep="")));
     
  eval(parse(text=GetG0Text("OriginalOldwd", "globalenv()", S=1)));
  try(setwd(OriginalOldwd));
  if (verbose >= 3) {
    AFilePrint(paste("  --- AssessedCI: Finish ", NameFunction, " length = ", length(AllCIReturn), sep="")); flush.console();
  }
  return(AllCIReturn);
}



FTCIEst <- function(FittedCIs, CIQuantiles, RealBetas, NameOfFitted,
  NameFunction, Hit, MIPReport) {
  TableBetas <- table(RealBetas);
  uBetas <- sort(unique(RealBetas));
  AvFittedCIs <- matrix(NA, length(TableBetas), length(CIQuantiles),
    dimnames=list(Betas=names(TableBetas), Quantiles=CIQuantiles));
  for (ii in 1:length(TableBetas)) {
    ASet <-  FittedCIs[RealBetas==uBetas[ii],]
    if (length(ASet) == length(CIQuantiles)) {
       AvFittedCIs[ii,] <- as.vector(ASet);
    }  else {
        AvFittedCIs[ii,] <- colMeans(ASet)
    }
    
  }  
  UniqueBetas <- sort(unique(RealBetas));
  NBeta <- length(TableBetas);
  MedianAbsDiffs <- rep(0, NBeta);
  MedianSqDiffs <- rep(0,NBeta);
  QuantileCov <- rep(0, floor(length(CIQuantiles) / 2));
  AK <- 0;  OS <- 1;
  if (CIQuantiles[1] == .5) { AK <- 1;}
  while(AK < length(CIQuantiles)) {
    QuantileCov[OS] <- CIQuantiles[AK+2]-CIQuantiles[AK+1]; AK <- AK+2; OS <- OS + 1;
  }
  iCoverageQuantiles <- matrix(0, length(TableBetas), length(QuantileCov),
    dimnames=list(Betas=names(TableBetas), Quantiles=QuantileCov));
  iMeanWidthQuantiles <-  matrix(0, length(TableBetas), length(QuantileCov),
    dimnames=list(Betas=names(TableBetas), Quantiles=QuantileCov));
  iMeanSqWidthQuantiles <- matrix(0, length(TableBetas), length(QuantileCov),
    dimnames=list(Betas=names(TableBetas), Quantiles=QuantileCov));
  iCINormalLoss <- matrix(0, length(TableBetas), length(QuantileCov),
    dimnames=list(Betas=names(TableBetas), Quantiles=QuantileCov));
  iCIUnifLoss <- matrix(0, length(TableBetas), length(QuantileCov),
    dimnames=list(Betas=names(TableBetas), Quantiles=QuantileCov));
  iMIPReport <- rep(0, length(UniqueBetas));
  try(names(iMIPReport) <- names(UniqueBetas));
  for (ii in 1:length(UniqueBetas)) {
    IndexMyBeta <- (1:length(RealBetas))[abs(RealBetas -UniqueBetas[ii]) <= .001]; 
    if (length(IndexMyBeta) > 1) {
      iMIPReport[ii] <- mean(MIPReport[IndexMyBeta])
    } else if (length(IndexMyBeta) == 1) {
      iMIPReport[ii] <- MIPReport[IndexMyBeta];
    }
  }
  Betas <- RealBetas;
  for (ii in 1:length(TableBetas)) {
    IndexMyBeta <- (1:length(RealBetas))[RealBetas == UniqueBetas[ii]];
    if (any(CIQuantiles == .5)) {
      MedianAbsDiffs[ii] <-sum(abs(Betas[IndexMyBeta] - FittedCIs[IndexMyBeta,CIQuantiles == .5]))
      MedianSqDiffs[ii] <- sum((Betas[IndexMyBeta] - FittedCIs[IndexMyBeta,CIQuantiles == .5])^2)
    } else {
      MedianAbsDiffs[ii] <- 0;  MedianSqDiffs[ii] <- 0;

    }
    OnAK <- 0;  OnOS <- 1; if (CIQuantiles[1] == .5) { OnAK <- 1; }
    for (OnOS in 1:length(QuantileCov)) {
      iCoverageQuantiles[ii, OnOS] <- sum( FittedCIs[IndexMyBeta,OnAK+1] <=
        Betas[IndexMyBeta] & 
        FittedCIs[IndexMyBeta,OnAK+2] >= Betas[IndexMyBeta]) / 
        length(Betas[IndexMyBeta]);
      iCINormal <- CINewLoss(FittedCIs[IndexMyBeta, OnAK+1], 
        FittedCIs[IndexMyBeta, OnAK+2], Betas[IndexMyBeta], QuantileCov[OnOS]);
      iCINormalLoss[ii, OnOS] <- sum(iCINormal);
      try(iCIUnif <- CIUnifLoss(FittedCIs[IndexMyBeta, OnAK+1], 
        FittedCIs[IndexMyBeta, OnAK+2], Betas[IndexMyBeta], QuantileCov[OnOS]));
      iCIUnifLoss[ii, OnOS] <- sum(iCIUnif);
        
      iMeanWidthQuantiles[ii, OnOS] <- mean(FittedCIs[IndexMyBeta,OnAK+2]-
        FittedCIs[IndexMyBeta,OnAK+1]);
      iMeanSqWidthQuantiles[ii, OnOS] <- mean((FittedCIs[IndexMyBeta,OnAK+2]-
           FittedCIs[IndexMyBeta,OnAK+1])^2);
      OnAK <- OnAK+2;        
    }
  }
  APiece <- list(
      NameOfInterval = NameOfFitted, NameFunction = NameFunction, Hit=Hit,
      CIQuantiles = CIQuantiles,
      MedianAbsDiffs = MedianAbsDiffs,  MedianSqDiffs = MedianSqDiffs,   
      iCoverageQuantiles = iCoverageQuantiles,
      iMeanWidthQuantiles = iMeanWidthQuantiles,
      iMeanSqWidthQuantiles = iMeanSqWidthQuantiles,
      iCINormalLoss =  iCINormalLoss,
      iCIUnifLoss = iCIUnifLoss, iMIPReport = iMIPReport, FittedCIs = AvFittedCIs, 
      TableBetas = TableBetas, UniqueBetas = UniqueBetas);
  return(APiece);
}

ProcessAssessCIFitAllNonZeroBeta = function(SMS, FittedOb,
     RunProcessIdentifier,
     UniqueSimulationIdentifier, 
     NameFunction, 
     Hit = "+",
     quoteMore = " - ", 
     verbose = 0, DontRecord=FALSE, AlreadyLocked=FALSE) {
     TableNonZeroBeta <- table(SMS$BetasReal);
     if (is.null(FittedOb$CIEst)) {
       AFilePrint(paste("AssessCIFit: Hey, we tried to do CI fit for ", NameFunction, " but got NULL!"));
     } 
     MyTestT <- " 
     if (is.logical(verbose) && verbose == FALSE) {
       verbose = FALSE;
     } else if (is.logical(verbose) && verbose == TRUE) { verbose = 2; 
     } else if (!is.numeric(verbose)) {
       printf(\"AssessCIFit: Hey, this isn't a good verbose;\"); verbose = 2;
     } else {
       verbose = verbose-2;
     } ";
     try(eval(parse(text=MyTestT)));
     if (verbose >= 1) {
       AFilePrint(paste("AssessCIFit; Starting for ", NameFunction, ", Hit = ", Hit, sep=""));
       flush.console();
     }
     
     Betas <- SMS$BetasReal;
     IndexOnBetas <- (1:length(Betas))[Betas!=0.0];
     if (SMS$p <= 1200) {
       FitBeta <- FittedOb$BetaFit;
     } else {
       FitBeta <- NULL;
     }
     
     
     IndexOffBetas <- (1:length(Betas))[Betas==0.0];
     IndexSmallBetas <- (1:length(Betas))[Betas!=0.0 &
       Betas^2 <= SMS$SigmaNoiseReal / length(SMS$NLen)];
     CIQuantiles <- FittedOb$CIQuantiles; 
     if (CIQuantiles[1] == .5) {
       AK <- 1; 
     } else { AK <- 0;}
     QuantileCov <- rep(0, length(CIQuantiles)/2);
     ArAt <- AK;
     for (jj in 1:length(QuantileCov)) {
       QuantileCov[jj] <- CIQuantiles[ArAt + 2] - CIQuantiles[ArAt+1];
       ArAt <- ArAt +2;
     }
     CITime = FittedOb$CITime;
     
     AllCIReturn <- list();
     AllNames <- rep("", length(FittedOb$CIEst));
     ##for (ii in 1:length(TableNonZeroBeta)) {
     for (tjt in 1:length(FittedOb$CIEst)) {
       APiece <- FTCIEst(FittedCIs = FittedOb$CIEst[[tjt]], CIQuantiles = FittedOb$CIQuantiles, 
       SMS$BetasReal, NameOfFitted = names(FittedOb$CIEst)[tjt],
       NameFunction = NameFunction, Hit = Hit, FittedOb$MIPReport );
       AllCIReturn[[tjt]] <- APiece;
       try(names(AllCIReturn)[tjt] <- names(FittedOb$CIEst)[tjt]);
       AllNames[tjt] <- names(FittedOb$CIEst)[tjt];
     }
    try(names(AllCIReturn) <- AllNames);
    
   
  if (exists("DontRecord") && is.logical(DontRecord) && DontRecord == TRUE) {
    return(AllCIReturn);
  }
  eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
  if (verbose >= 2) {
      AFilePrint(paste("----  AssessCIFit[", UniqueSimulationIdentifier, 
        "] --- CurrentSmallContainDir = ", CurrentSmallContainDir, " for NameFunction = ", NameFunction, sep=""));
      flush.console();
  }
  BADUniqueProcessIdentifier <- 1;
  ATryGo <- "
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  }
  if (length(UniqueProcessIdentifier) > 1) {
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1]);
  }
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  } else if (is.null(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NULL\");
    try(UniqueProcessIdentifier <- \"\"); 
  } else if (length(UniqueProcessIdentifier) > 1) {
    UniqueProcessIdentifier <- UniqueProcessIdentifier[1]; 
  } else if (is.na(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NA\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.numeric(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NUMERIC\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.character(UniqueProcessIdentifier)) {
    if (UniqueProcessIdentifier[1] == \"\")  {
      AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is Blank!\");
    } else {
      BADUniqueProcessIdentifier <- 0;
    }
  } else {
    try(UniqueProcessIdentifier <- \"\");
  }
  ";
  try(eval(parse(text=ATryGo)))
  if (BADUniqueProcessIdentifier == 1) {
    AFilePrint("AssessCIFit: UnqiueProcessIdentifier is NULL!"); flush.console();  
  }
  if (verbose >= 2) {
    AFilePrint(paste("  --- AssedCIFit, ", NameFunction, " looking for HMS. ", sep="")); flush.console();
  }
  CIEst = NULL;
  if (!exists("FittedOb") || is.null(FittedOb)) {
    AFilePrint("Why did FittedOb AssesCI? that not happen what?"); flush.console(); 
  } else if (is.null(FittedOb$CIEst)) {
    AFilePrint("What: Why is FittedOb$CIEst null?"); flush.console();
  } else {
    CIEst <- FittedOb$CIEst;
  }
  TemperatureList <- NULL;  BetaStartOn <- NULL;   IndexOnBetaStart <- NULL;
  NoNoiseBetaStart <- NULL;
  if (!is.null(FittedOb$BetaStart)) { 
    IndexOnBetaStart <- (1:length(FittedOb$BetaStart))[FittedOb$BetaStart != 0]
    BetaStartOn <- FittedOb$BetaStart[IndexOnBetaStart];
    NoNoiseBetaStart <- TRUE;
  }
  if (!is.null(FittedOb$TemperatureList)) {
    TemperatureList <- FittedOb$TemperatureList;
  }
  eval(parse(text=GetG0Text("CurrentSmallContainDir", "globalenv()", S=1)));
  try(setwd(CurrentSmallContainDir));
  eval(parse(text=GetG0Text("ZZMyID")));
  SimulationIdentifier   <- SMS$UniqueSimulationIdentifier;
  eval(parse(text=GetG0Text("UniqueProcessIdentifier")));
  RunProcessIdentifier <- UniqueProcessIdentifier;
  FunctionName <- NameFunction;
  print(paste("About to save AllCIReturn to ",
    paste("Processes//", ZZMyID, 
     "//BetaCI", FunctionName, "S", SimulationIdentifier, ".RData", sep=""), sep=""));
  flush.console();
  
  ##Betas <- SMS$BetasReal;
  ##IndexOnBetas <- (1:length(Betas))[Betas!=0.0];
    
  try(BetasOn <- Betas[IndexOnBetas]);
  try(FitBetaOn <- FittedOb$BetaFit[IndexOnBetas]);
  try(MIPBetaOn <- FittedOb$MIPReport[IndexOnBetas]);
  try(MeanMIPBetaOff <- mean(FittedOb$MIPReport[IndexOffBetas])); 
  if (!exists("IndexOnBetas")) {
    print("BetaCI: ERROR IndexOnBetas doesn't exist not generated? "); flush.console();
  }    
  if (!exists("FitBetaOn")) {
    print("BetaCI: ERROR FitBetaOn doesn't exist not generated? "); flush.console();
  }   
  if (!exists("MIPBetaOn")) {
    print("BetaCI: ERROR MIPBetaOn doesn't exist not generated? "); flush.console();
  }   
  if (!exists("MeanMIPBetaOff")) {
    print("BetaCI: ERROR MeanMIPBetaOff doesn't exist not generated? "); flush.console();
  }   
  try(save(AllCIReturn=AllCIReturn, NameFunction=NameFunction, 
    FunctionName = FunctionName, IndexOnBetas = IndexOnBetas, 
    BetasOn=BetasOn, MIPBetaOn = MIPBetaOn, MeanMIPBetaOff = MeanMIPBetaOff, 
    TemperatureList=TemperatureList, IndexOnBetaStart=IndexOnBetaStart, BetaStartOn=BetaStartOn,
    NoNoiseBetaStart=NoNoiseBetaStart, 
    CIQuantiles=CIQuantiles, QuantileCov = QuantileCov, 
    CIEst = CIEst, Hit=Hit, CITime=CITime, RunProcessIdentifier = RunProcessIdentifier, 
    SimulationIdentifier = SimulationIdentifier, file=paste("Processes//", ZZMyID, 
     "//BetaCI", FunctionName, "S", SimulationIdentifier, ".RData", sep="")));
   
  print("Successfully Saved Final BetaCI. "); flush.console();  
  eval(parse(text=GetG0Text("OriginalOldwd", "globalenv()", S=1)));
  try(setwd(OriginalOldwd));
  if (verbose >= 3) {
    AFilePrint(paste("  --- AssessedCI: Finish ", NameFunction, " length = ", length(AllCIReturn), sep="")); flush.console();
  }
  return(AllCIReturn);
}



RecoverProgress <- function(MyL) {
    eval(parse(text=GetG0Text("ZZMyID", "globalenv()", S=1)));
    ASS <- MyL[substr(MyL,1,1)=="S"];  
      eval(parse(text=GetG0Text("SubMitList", "globalenv()", S=1)));
      ProgressNames <- list()
      ProgressData <- list();
      AGo <- 0;
      for (ktk in 1:length(ASS)) {
         MyV = NULL;
         try(load(paste("Processes//", ZZMyID, "//", ASS[ktk], sep="")));
         ##ASLoc <- (1:NROW(SubMitList))[
         ##  (SubMitList[,1] == SimulationIdentifier  | 
         ##   SubMitList[,1] == paste("SMS", SimulationIdentifier, sep="") |
         ##   SubMitList[,1] == paste("SMS", SimulationIdentifier, ".RData", sep="")
         ##   ) &
         ##  SubMitList[,2] == NameFunction];
         ##if (length(ASLoc) == 1) {
         ##   ProgressNames[ASLoc,] <- c(SimulationIdentifier, NameFunction);
         ##   ProgressData[ASLoc,] <- MyV;
         ##}
         if (!is.null(MyV)) {
            AGo <- AGo+1;
            ProgressNames[[AGo]] <- c(SimulationIdentifier, NameFunction); 
            if (is.na(MyV[1])) {
              AttemptSucceedFail <- -1;  
            } else { AttemptSucceedFail <- 1;}
            if (!exists("StartsTime")) { StartsTime <- GetRealTime(); }
            if (!exists("FinishesTime")) { FinishesTime <- GetRealTime(); }
            ProgressData[[AGo]] <- c( 
             StartDays=StartsTime$JustDays, StartSeconds=StartsTime$JustSeconds,
             AttemptSucceedFail = AttemptSucceedFail, ProcessTime1=MyV[6], 
             ProcessTime2=MyV[7], ProcessTime3=MyV[8],
             L2Accuracy=MyV[2], FinishDays=FinishesTime$JustDays,
             FinishSeconds=FinishesTime$JustSeconds)
         }
      }  
      try(save(ProgressNames=ProgressNames, ProgressData=ProgressData,
        file=paste("Processes//", ZZMyID, "//MyProcesses.RData", sep=""))
      );
      eval(parse(text=SetGText("ProgressNames", "globalenv()", S=1)));
      eval(parse(text=SetGText("ProgressData", "globalenv()", S=1)));  
      RTText <- "
      eval(parse(text=GetG0Text(\"ProgressNames\", \"globalenv()\", S=1)));
      eval(parse(text=GetG0Text(\"ProgressData\", \"globalenv()\", S=1)));        
      "  
      return(RTText);
}