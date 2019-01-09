################################################################################
##  LookToGetNewSimulation()
##    (c) 2009-2018 Alan Lenarcic
##
##   As commanded by "AssessFirstTime()" this will look into which Estimators and Simulatiors
##  haven't been conducted and uses SimMeData() to get new simulations.
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
LookToGetNewSimulation <- function(verbose = NULL, quoteMore = "LookForNewSimulation : ",
  ForceFunction = NULL, TCS = NULL, ...) {
  eval(parse(text=GetG0Text("CurrentLargeContainDir", S=1)));
  eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
  if (is.numeric(CurrentLargeContainDir) && CurrentLargeContainDir[1] == 0) {
    eval(parse(text=GetG0Text("CurrentLargeContainDir", S=1)));
    eval(parse(text=GetG0Text("n", S=1)));
    eval(parse(text=GetG0Text("p", S=1)));
    eval(parse(text=GetG0Text("k", S=1)));
    eval(parse(text=GetG0Text("sigma", S=1)));
    eval(parse(text=GetG0Text("NumReps", S=1)));
    eval(parse(text=GetG0Text("LARSSeekFlag", S=1)));
    eval(parse(text=GetG0Text("DefaultGenerateBetaVec", S=1)));
    eval(parse(text=GetG0Text("DefaultCorrelationXMatrix", S=1)));
    eval(parse(text=GetG0Text("AllDir", S=1)));
    CurrentLargeContainDir <- paste(AllDir, "//",  DirNamingClenture( 
      DefaultGenerateBetaVec, DefaultCorrelationXmatrix,
      NLen = n, kLen = p, kActiveLen = k, SigmaNoise = sigma, 
      NumReps = NumReps, LARSSeekFlag  = LARSSeekFlag), sep="");
      dir.create(CurrentLargeContainDir, showWarnings = FALSE);   
    eval(parse(text=SetGText("CurrentLargeContainDir", S=1))); 
  }
  if (is.numeric(CurrentSmallContainDir) && CurrentSmallContainDir[1] == 0) {
    eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
    eval(parse(text=GetG0Text("n", S=1)));
    eval(parse(text=GetG0Text("p", S=1)));
    eval(parse(text=GetG0Text("k", S=1)));
    eval(parse(text=GetG0Text("sigma", S=1)));
    eval(parse(text=GetG0Text("NumReps", S=1)));
    eval(parse(text=GetG0Text("LARSSeekFlag", S=1)));
    eval(parse(text=GetG0Text("DefaultGenerateBetaVec", S=1)));
    eval(parse(text=GetG0Text("DefaultCorrelationXMatrix", S=1)));
    eval(parse(text=GetG0Text("AllDir", S=1)));
    CurrentSmallContainDir <- paste(AllDir, "//",  DirNamingClenture( 
      DefaultGenerateBetaVec, DefaultCorrelationXmatrix,
      NLen = n, kLen = p, kActiveLen = k, SigmaNoise = sigma, 
      NumReps = NumReps, LARSSeekFlag  = LARSSeekFlag), sep="");
      dir.create(CurrentSmallContainDir, showWarnings = FALSE);   
    eval(parse(text=SetGText("CurrentSmallContainDir", S=1))); 
  }
  CurrentTotalCompleted <- 0;
  
  if (is.null(verbose)) { 
    eval(parse(text=GetG0Text("verbose")));
    if (verbose == 0) {
      verbose = 0;
    } 
  } else { verbose = as.numeric(verbose); }
  MySummaryPerformanceDirectory = paste(
    CurrentLargeContainDir, "//", "SumPerformance", sep="");
  dir.create(MySummaryPerformanceDirectory, showWarnings=FALSE, recursive=TRUE);
  if (verbose > 0) {
    AFilePrint("Looking for New Simulation Start with LockMeIn"); 
    flush.console(); 
  }
  try(Oldwd <- getwd());
  try(setwd(CurrentSmallContainDir));
  while(TwoSimR5:::LockMeIn(verbose=as.numeric(verbose), quoteMore=quoteMore, 
    LFDir = "SumPerformance", LockIn="SumLock")==FALSE){
    try(Sys.sleep(runif(1,0,4))); }
  try(setwd(Oldwd));
  try(MySummarySimulationsList <- NULL);
  try(eval(parse(text=SetGText("MySummarySimulationsList", 
    "globalenv()", S=1))))
  MyListfiles = list.files(MySummaryPerformanceDirectory);
  try(MySummarySimulationsList <- NULL);

  try(setwd(CurrentSmallContainDir));
  try(load(paste(
    "SumPerformance//SummaryPerformanceRData.RData", sep="")));
  try(eval(parse(text=TwoSimR5:::RecoveryTSSLText(CurrentLargeContainDir))));
  
  if (is.null(MySummarySimulationsList)) {
    paste(AFilePrint("**GetSims: we did not load in MySummarySimulationsList. "));
  }
  GetSims <- NULL;
  
  try(GetSims <- MySummarySimulationsList$GetMeASimulation(
    AlreadyLocked=TRUE, ForceFunction = ForceFunction, TCS = TCS));
  try(eval(parse(text=SetGText("GetSims", "globalenv()", S=1))))
  if (!exists("GetSims") || is.null(GetSims) || is.numeric(GetSims) || 
    !is.list(GetSims) || is.null(GetSims$SMS)) {
    AFilePrint("*************************************************************");
    AFilePrint(paste("** Error inspection: LocateASim: GetSims was returned ",
      "as a NULL, we fail", sep=""));
    AFilePrint("** SaveMySummarySimulationList to terminal: ");
    eval(parse(text=SetGText("MySummarySimulationsList", "globalenv()", S=1)));
    return(-1);
    tryCatch("We are in it bad!");
  }

    Oldwd <- getwd();
    try(setwd(CurrentSmallContainDir));
    try(secure.save.MSS(AObject=MySummarySimulationsList, 
      ObjectName="MySummarySimulationsList",
      file=paste("SumPerformance//SummaryPerformanceRData.RData", sep="")));
    try(setwd(Oldwd)); 
  try(setwd(CurrentSmallContainDir));
  UnLockMe(verbose=verbose, quoteMore=quoteMore, 
    LFDir = "SumPerformance", LockIn="SumLock");
  try(setwd(Oldwd)); 
    
  try(eval(parse(text=TryFillCurrentTotalCompleted())));

 
  if (is.character(GetSims) && GetSims[1] == "AllComplete") {
    if (verbose >= 2) {
      AFilePrint("We received a Getsims is AllComplete message from GetMeASimulation");
      AFilePrint("End Get Sims. ");
    }
    return(111);
  }
  if ((is.character(GetSims) && Getsims[1] == "AFail") ||
    (is.numeric(GetSims) && length(GetSims) == 1 && GetSims <= -1)) {
    AFilePrint("*****************************************************************");
    AFilePrint("*****************************************************************");
    AFilePrint("LookToGetNewSimulation: Fail From GetSims!");
    AFilePrint("LookToGetNewSimulation: We were delivered -1 from Get Sims");
    AFilePrint(paste("Our MySummaryPerformanceDirectory is: ",
       MySummaryPerformanceDirectory, sep=""));
    AErrorPrint(paste("LookToGetNewSimulation: MySummaryPerformanceDirectory",
      " results in Error.", sep="")); flush.console();
    return(-1);
    tryCatch("Fleeing LookToGetNewSimulation (In GetANewSimulation.r)")
  }
    
  ##UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory);

  
  #LockedSMSFile <- LockInSMS(SMS, CurrentContainDir = CurrentContainDir, OnFunction=NameFunctions[ISample]);    

  return(list(SMS=GetSims$SMS, UniqueSimulationIdentifier=GetSims$UniqueSimulationIdentifier,
             InSimulationsList=GetSims$InSimulationsList, 
             OnFunction = GetSims$OnFunction,
             OnNumFunction=GetSims$OnNumFunction));
}

##############################################################################
## GenerateBetaLongQuantiles
SummarySimulationList$methods(
  GetMeASimulation = function(AlreadyLocked=TRUE,ForceFunction=NULL, 
    TCS = NULL,...) {
     
     eval(parse(text=GetG0Text("TargetTotalSimulations")));
     if (!is.null(ForceFunction)) {
       if (is.character(ForceFunction)) {
         AT <- (1:length(.self$NameFunctions))[
           .self$NameFunctions == ForceFunction];
         if (length(AT) == 1) {
           OnNumForceFunction = AT;
         } else {
           OnNumForceFunction = NULL;  ForceFunction = NULL;
         }
       } else if (is.numeric(ForceFunction)) {
         if (round(ForceFunction) >= 1   && 
           round(ForceFunction) <= length(.self$NameFunctions)) {
           OnNumForceFunction <- round(ForceFunction);
           ForceFunction <- .self$NameFunctions[OnNumForceFunction]
         }  else {
           OnNumForceFunction = NULL;  Forcefunction = NULL;
         }
       }
     } else {OnNumForceFunction <- NULL; }
     if (is.null(.self$TheSummarySimulationList) ||
        length(.self$TheSummarySimulationList) <= 0) {
        SMS <- GenerateANewSimulation(AlreadyLocked=AlreadyLocked)
        UniqueSimulationIdentifier <- SMS$UniqueSimulationIdentifier;
        .self$TheSummarySimulationList[[1]] <- SummarySimulationRecord$new(
          UniqueSimulationIdentifier=SMS$UniqueSimulationIdentifier,
          SavedUniqueSimulationRDataFile = paste(.self$LargeContainingDir, 
          "//SMS", SMS$UniqueSimulationIdentifier, ".RData", sep=""));
        try(names(.self$TheSummarySimulationList)[1] <- 
          SMS$UniqueSimulationIdentifier);
        names(.self$TheSummarySimulationList)[1] <-  
          SMS$UniqueSimulationIdentifier;
        VLD <- .self$ValidSimulations;
        SampleAFunction <- sample(VLD, prob=rep(1, length(VLD)),
          size=1, replace=FALSE);
        if (!is.null(OnNumForceFunction)) {
          SampleAFunction <- OnNumForceFunction;
        }
        if (SampleAFunction <= 1) {
           AFilePrint(paste("Error In GetASimulation!  ",
             "SampleAFunction = ", SampleAFunction, sep=""));
        }
        TryAddText <- "
        SuccAdd <- FALSE;   AOn <- -1;
        AOn <- .self$TheSummarySimulationList[[1]]$AddAFunction(
          as.integer(SampleAFunction), SSL=.self);
        SuccAdd <- TRUE;
        ";
        try(eval(parse(text=TryAddText)));
        if (is.character(AOn) && AOn == "AFail") {
          AFilePrint("GetANewSimulation: We get A FAIL!");
        }
        AFilePrint(paste("***  AfterAddAFunction: Well \n Well ",
          "\n Well: We get AOn = ", AOn, sep=""));
        if (SuccAdd == FALSE || (is.numeric(AOn) && 
          AOn[1] != SampleAFunction) || !is.numeric(AOn) ||
          (is.numeric(AOn) && AOn < 0) || 
          (is.character(AOn) && AOn == "AFail")) {
          AFilePrint(paste("**********************************",
            "************************", sep=""));
          AFilePrint(paste("GetANewSimulation.r On[", 1, 
            "] Try add function ", SampleAFunction, " fails!", sep=""));
          AText <- "
          TheSummarySimulationList <- .self$TheSummarySimulationList;
          MySummarySimulationStore <- .self;
          ";
          try(eval(parse(text=AText)));
          eval(parse(text=SetGText("TheSummarySimulationList", S=1)));
          eval(parse(text=SetGText("MySummarySimulationStore", S=1)));  
          eval(parse(text=SetGText("SampleAFunction", S=1)));
          return("AFail")        
          tryCatch("GetANewSimulation.r Fail!")
        }
        return(list(SMS=SMS, UniqueSimulationIdentifier=
          SMS$UniqueSimulationIdentifier,
          InSimulationsList=1, OnFunction = 
          .self$NameFunctions[SampleAFunction],
          OnNumFunction=SampleAFunction,
          CurrentTotalCompleted = CurrentTotalCompleted));
     } 
     VLD <- NULL; TC = NULL; TK = NULL;
     try(VLD <- .self$ValidSimulations);
     try(TC <- .self$GetTotalCompletedFunction());
     try(TK <- .self$MyTotalKilledFunctions);
     try(eval(parse(text=TryFillCurrentTotalCompletedTC())));
     if (all(TC[VLD] >=  TargetTotalSimulations)) {
       return("AllComplete");
     }
    if (exists("OnNumForceFunction") &&
    !is.null(OnNumForceFunction) && OnNumForceFunction >= 1) {
     if (TC[OnNumForceFunction] + TK[OnNumForceFunction] >= 
       length(.self$TheSummarySimulationList)) {
        SMS <- GenerateANewSimulation(AlreadyLocked=AlreadyLocked)
        UniqueSimulationIdentifier <- SMS$UniqueSimulationIdentifier;
        NewN <- length(.self$TheSummarySimulationList)
        .self$TheSummarySimulationList[[NewN]] <- SummarySimulationRecord$new(
          UniqueSimulationIdentifier=SMS$UniqueSimulationIdentifier,
          SavedUniqueSimulationRDataFile = paste(.self$LargeContainingDir, 
          "//SMS", SMS$UniqueSimulationIdentifier, ".RData", sep=""));
        try(names(.self$TheSummarySimulationList)[NewN] <-
          SMS$UniqueSimulationIdentifier);
        names(.self$TheSummarySimulationList)[NewN] <-  
          SMS$UniqueSimulationIdentifier;
        SampleAFunction <- OnNumForceFunction;
        TryAddText <- "
        SuccAdd <- FALSE;   AOn <- -1;
        AOn <- .self$TheSummarySimulationList[[NewN]]$AddAFunction(
          as.integer(SampleAFunction), SSL = .self);
        SuccAdd <- TRUE;
        ";
        try(eval(parse(text=TryAddText)));  
        return(list(SMS=SMS, UniqueSimulationIdentifier=
          SMS$UniqueSimulationIdentifier,
          InSimulationsList=NewN, OnFunction = 
          .self$NameFunctions[SampleAFunction],
          OnNumFunction=SampleAFunction));      
     }
  }
  if (all(TC[VLD] + TK[VLD] >= length(.self$TheSummarySimulationList))) {
        SMS <- GenerateANewSimulation(AlreadyLocked=AlreadyLocked)
        try(NewSim <- length(.self$TheSummarySimulationList)+1);
        try(.self$TheSummarySimulationList[[NewSim]] <- 
          SummarySimulationRecord$new(SMS$UniqueSimulationIdentifier,
          SavedUniqueSimulationRDataFile = paste(.self$LargeContainingDir, 
          "//SMS", SMS$UniqueSimulationIdentifier, ".RData", sep="")));
        try(names(.self$TheSummarySimulationList)[NewSim] <- 
          SMS$UniqueSimulationIdentifier);
        SampleAFunction <- sample(VLD, 
          prob=.self$SampleFunctionProbabilities[VLD],
          size=1);
        if (!is.null(OnNumForceFunction)) { 
          SampleAFunction <- OnNumForceFunction; }
        TryAddText <- "
        SuccAdd <- FALSE;
        AOn <- -1;
        AOn <- .self$TheSummarySimulationList[[NewSim]]$AddAFunction(
          as.integer(SampleAFunction), SSL = .self);
        SuccAdd <- TRUE;
        ";
        try(eval(parse(text=TryAddText)));
        AFilePrint(paste("***  AfterAddAFunction: Well \n Well ",
          "\n Well: We get AOn = ", AOn, sep=""));
        if (is.character(AOn) && AOn == "AFail") {
          AFilePrint("GetANewSimulation: We get A FAIL!");
        }
        if (SuccAdd == FALSE || (is.numeric(AOn) && 
           AOn[1] != SampleAFunction) || !is.numeric(AOn)
           ||
          (is.numeric(AOn) && AOn < 0) || 
          (is.character(AOn) && AOn == "AFail")) {
          AFilePrint("ERROR ERROR, TC > SSL ********************************");
          AFilePrint(paste("TC > TSSL: GetANewSimulation.r On[", NewSim, 
            "] Try add function ", SampleAFunction, " fails!", sep=""));
          ATText <- "
          TheSummarySimulationList <- .self$TheSummarySimulationList;
          MySummarySimulationStore <- .self;
          ";
          try(eval(parse(text=ATText)));
          eval(parse(text=SetGText("TheSummarySimulationList", S=1)));
          eval(parse(text=SetGText("MySummarySimulationStore", S=1)));  
          eval(parse(text=SetGText("SampleAFunction", S=1)));        
          eval(parse(text=SetGText("NewSim", S=1)));  
          return("AFail")
          tryCatch("GetANewSimulation.r Fail!")
        }
        return(list(SMS=SMS, 
          UniqueSimulationIdentifier=SMS$UniqueSimulationIdentifier,
          InSimulationsList=NewSim, 
          OnFunction = .self$NameFunctions[SampleAFunction],
          OnNumFunction=SampleAFunction));
     }
     MyT <- "
     try(PerformanceMatrix <- NULL);
     ";
     try(eval(parse(text=MyT)));
     try(PerformanceMatrix <- .self$GenPerformanceMatrix());
     if (is.null(PerformanceMatrix)) {
       AFilePrint(paste("GenPerformanceMatrix Failed, save ",
         "TheSummarySimulationList to terminal", sep=""));
       ATText <- "
       TheSummarySimulationList <- .self$TheSummarySimulationList
       MySummarySimulationListOnSave <- .self;
       ";
       try(eval(parse(text=ATText)));
       eval(parse(text=SetGText("TheSummarySimulationList", S=1)))
       eval(parse(text=SetGText("MySummarySimulationListOnSave", S=1)))       
     }
     if (!is.null(PerformanceMatrix) && !is.null(OnNumForceFunction) && 
       OnNumForceFunction >= 0) {
       SampleAFunction <- OnNumForceFunction; 
       if (any(PerformanceMatrix[,SampleAFunction] >= 0)) {
          AVT <- rep(0, NROW(PerformanceMatrix));
          if (any(PerformanceMatrix[,SampleAFunction] == 0)) {
            AVT[PerformanceMatrix[,SampleAFunction] == 0] <- 1;
          } else {
            AID <- PerformanceMatrix[,SampleAFunction] >= 0; 
            AVT[AID] <- 
              max(PerformanceMatrix[AID,SampleAFunction]) - 
              PerformanceMatrix[AID,SampleAFunction]+1;
          }
          if (!any(AVT > 0)) {
             AFilePrint(paste("Error in Search PerformanceMatrix for ",
               " force function: ", ForceFunction, sep=""));
             AFilePrint(paste("Performance was : [",
               paste(PerformanceMatrix[,SampleAFunction], collapse = ", "),
               "]", sep=""));
             AFilePrint(paste("Something is wrong in the ",
               "GetANewsimulation and you shoule come and fight it.", sep=""));
             AErrorPrint("Error Trying to get ForceFunction Simulation!");
             return(-1);
          } else {
             SampleASimulation <- sample(
               1:NROW(PerformanceMatrix), prob=AVT,size=1);
          }
       } else {
         SMS <- GenerateANewSimulation(AlreadyLocked=AlreadyLocked)
         try(NewSim <- length(.self$TheSummarySimulationList)+1);
         try(.self$TheSummarySimulationList[[NewSim]] <- 
           SummarySimulationRecord$new(SMS$UniqueSimulationIdentifier,
           SavedUniqueSimulationRDataFile = paste(.self$LargeContainingDir, 
           "//SMS", SMS$UniqueSimulationIdentifier, ".RData", sep="")));
          try(names(.self$TheSummarySimulationList)[NewSim] <- 
          SMS$UniqueSimulationIdentifier);
          SampleASimulation <- NewSim;
       } 
       TryAddText <- "
         SuccAdd <- FALSE;
         AOn <- -1;
         AOn <- .self$TheSummarySimulationList[[
           SampleASimulation]]$AddAFunction(as.integer(SampleAFunction), 
             SSL=.self);
         SuccAdd <- TRUE;
         ";
         try(eval(parse(text=TryAddText)));        
         try(UniqueSimulationIdentifier <- 
           .self$TheSummarySimulationList[[
             SampleASimulation]]$UniqueSimulationIdentifier);
         load(paste(.self$MyBaseSimulationStorageDirectory, "//SMS", 
           .self$TheSummarySimulationList[[
           SampleASimulation]]$UniqueSimulationIdentifier,
           ".RData", sep=""));
         return(list(SMS=SMS, 
           UniqueSimulationIdentifier = SMS$UniqueSimulationIdentifier,
           InSimulationsList=SampleASimulation, OnFunction = 
           .self$NameFunctions[SampleAFunction],
           OnNumFunction = SampleAFunction));
     }
     if (!is.null(PerformanceMatrix)) {
       if (any(PerformanceMatrix[,VLD] == 0)) {
         ATOs <- apply(PerformanceMatrix,2, function(x) { length(x[x==0]) });
         SampleAFunction <- sample(VLD, prob=ATOs[VLD], size=1);
         ARD <- (1:NROW(PerformanceMatrix));
         ATTs <- rep(0, length(ARD));
         ATTs[PerformanceMatrix[,SampleAFunction] == 0] <- 1;
         SampleASimulation <- sample(ARD, prob=ATTs, size=1);
         TryAddText <- "
         SuccAdd <- FALSE;
         AOn <- -1;
         AOn <- .self$TheSummarySimulationList[[
           SampleASimulation]]$AddAFunction(as.integer(SampleAFunction),
             SSL = .self);
         SuccAdd <- TRUE;
         ";
         try(eval(parse(text=TryAddText)));
         AFilePrint(paste("***  AfterAddAFunction: Well \n Well \n ",
           "Well: We get AOn = ", AOn, sep=""));
         if (is.character(AOn) && AOn == "AFail") {
           AFilePrint("GetANewSimulation: We get A FAIL!");
         }
         if (SuccAdd == FALSE || (is.numeric(AOn) 
           && AOn[1] != SampleAFunction) || !is.numeric(AOn)  ||
          (is.numeric(AOn) && AOn < 0) || 
          (is.character(AOn) && AOn == "AFail")) {
           AFilePrint(paste("SampleASimulation = ", SampleASimulation, 
             ": GetANewSimulation.r On[", SampleASimulation, 
             "] Try add function ", SampleAFunction, " fails!", sep=""));
           ATText <- "
           TheSummarySimulationList <- .self$TheSummarySimulationList;
           MySummarySimulationStore <- .self;
           ";
           try(eval(parse(text=ATText)));
           try(eval(parse(text=SetGText("ATOs", S=1))));
           try(eval(parse(text=SetGText("PerformanceMatrix", S=1))));
           eval(parse(text=SetGText("TheSummarySimulationList", S=1)));
           try(TSSL <- .self$TheSummarySimulationList, silent=TRUE);
           eval(parse(text=SetGText("TSSL", S=1)));
           eval(parse(text=SetGText("MySummarySimulationStore", S=1)));  
           eval(parse(text=SetGText("SampleAFunction", S=1)));  
           eval(parse(text=SetGText("SampleASimulation", S=1)));      
           return("AFail")
           tryCatch("GetANewSimulation.r Fail!")
         }
         try(UniqueSimulationIdentifier <- 
           .self$TheSummarySimulationList[[
             SampleASimulation]]$UniqueSimulationIdentifier);
         load(paste(.self$MyBaseSimulationStorageDirectory, "//SMS", 
           .self$TheSummarySimulationList[[
            SampleASimulation]]$UniqueSimulationIdentifier,".RData", sep=""));
         return(list(SMS=SMS, 
           UniqueSimulationIdentifier = SMS$UniqueSimulationIdentifier,
           InSimulationsList=SampleASimulation, OnFunction = 
             .self$NameFunctions[SampleAFunction],
           OnNumFunction = SampleAFunction));
       } else if (NROW(PerformanceMatrix) < TargetTotalSimulations) {
         SMS <- GenerateANewSimulation(AlreadyLocked=AlreadyLocked)
         UniqueSimulationIdentifier <- SMS$UniqueSimulationIdentifier;
         try(NewSim <- length(.self$TheSummarySimulationList)+1);
         try(.self$TheSummarySimulationList[[NewSim]] <- 
           SummarySimulationRecord$new(SMS$UniqueSimulationIdentifier));
         try(names(.self$TheSummarySimulationList)[NewSim] <- 
           .self$TheSummarySimulationList[[NewSim]]$UniqueSimulationIdentifier);
         SampleAFunction <- sample(VLD, 
           prob=.self$SampleFunctionProbabilities[VLD],
           size=1);
         if (!is.null(OnNumForceFunction)) { 
           SampleAFunction <- OnNumForceFunction; 
           ForcingFunction<-TRUE;           
         }  else { ForcingFunction <- FALSE;   }
           TryAddText <- "
             SuccAdd <- FALSE;   AOn <- -1;
             AOn <- .self$TheSummarySimulationList[[NewSim]]$AddAFunction(
               as.integer(SampleAFunction), ForcingFunction=ForcingFunction,
               SSL = .self);
            SuccAdd <- TRUE;
         ";
         try(eval(parse(text=TryAddText)));
         ##  AFilePrint(paste("***  AfterAddAFunction: 
         ## Well \n Well \n Well: We get AOn = ", AOn, sep=""));
         if (is.character(AOn) && AOn == "AFail") {
           AFilePrint("GetANewSimulation: We get A FAIL!");
         }
         if (SuccAdd == FALSE || (is.numeric(AOn) && 
          AOn[1] != SampleAFunction) || !is.numeric(AOn)  ||
          (is.numeric(AOn) && AOn < 0) || 
          (is.character(AOn) && AOn == "AFail")) {
           AFilePrint(paste("NewSim AddOn = ", NewSim, 
             ": GetANewSimulation.r On[", SampleASimulation, 
             "] Try add function ", SampleAFunction, " fails!", sep=""));
           ATText <- "
           TheSummarySimulationList <- .self$TheSummarySimulationList;
           MySummarySimulationStore <- .self;
           ";
           try(eval(parse(text=ATText)));
           try(eval(parse(text=SetGText("ATOs", S=1))));
           try(eval(parse(text=SetGText("PerformanceMatrix", S=1))));
           eval(parse(text=SetGText("TheSummarySimulationList", S=1)));
           eval(parse(text=SetGText("MySummarySimulationStore", S=1)));  
           eval(parse(text=SetGText("SampleAFunction", S=1)));  
           ##eval(parse(text=SetGText("SampleASimulation", S=1)));      
           eval(parse(text=SetGText("NewSim", S=1))); 
           return("AFail") 
           tryCatch("GetANewSimulation.r Fail!")
         }
         return(list(SMS=SMS, 
           UniqueSimulationIdentifier=SMS$UniqueSimulationIdentifier,
           InSimulationsList=NewSim, 
           OnFunction = .self$NameFunctions[SampleAFunction],
           OnNumFunction=SampleAFunction));
       } else {
         if (FALSE) {
         AFilePrint("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ");
         AFilePrint(paste("*** GetANewSimulation Hey we are at a ",
           "MaxDuplicateTrials call but it doesn't make sense", sep=""));
           AFilePrint("***We are not testing this out yet.  ")
           AFilePrint("***Save PerformanceMatrix to ask question!");
           eval(parse(text=SetGText("PerformanceMatrix", "globalenv()", S=1)));
           AFilePrint("*** Figure out source of error! ");
           return("AFail")
           tryCatch("GetANewSimulation.r Fail")
         }
         eval(parse(text=GetG0Text("MaxDuplicateTrials", "globalenv()", S=1)));
         if (MaxDuplicateTrials <= 1) { MaxDuplicateTrials <- 2; }
         LLens <- apply(PerformanceMatrix, 2, function(x) { length(x[x>=0]) });
         Acts <-  apply(PerformanceMatrix, 2, function(x) { 
           sum(x[x>=0], na.rm=TRUE) });
         Maxs <-  apply(PerformanceMatrix, 2, function(x) { 
           AT <- x[x>=0];
           if (length(AT) >= 1) { return(max(AT)); }
           return(0); }); 
         AllMaxs <- max(Maxs[VLD]);   
         DOs <- LLens * AllMaxs - Acts+1;  
         if (!any(DOs[VLD]>0)) {
           ATT <- paste("Hey, Maximum filled table can't find a ",
             "valid DO greater than zero
             VLD = (", paste(VLD, collapse=", "), ")
             DOs = (", paste(DOs, collapse=", "), ")
             LLens = (", paste(LLens, collapse=", "), ")
             Acts =  (", paste(LLens, collapse=", "), ")
             Maxs= (", paste(Maxs, collapse=", "), ")", sep="");
             
           AFilePrint(ATT);
           AErrorPrint(ATT);
         }   
         SampleAFunction <- sample(VLD, prob=DOs[VLD], size=1);
         ART <- PerformanceMatrix[,SampleAFunction];
         maxART <- max(ART);
         TrySims <- 1:length(ART);
         SProb <- rep(0, length(ART));
         if (maxART < 0) {
            AFilePrint("ERROR ERROR ");
            AFilePrint(paste("Hey, We sampled ", SampleAFunction, 
              " but ART is (", paste(ART, collapse=", "), ")", sep=""));
            AErrorPrint(paste("Hey, We sampled ", SampleAFunction, 
              " but ART is (", paste(ART, collapse=", "), ")", sep=""));
         } else if (maxART == 0) {
           SProb[ART == 0] <- 1;
         } else if (any(ART) == 0) {
           SProb[ART == 0] <- 1;
         } else {
           SProb[ART >= 0] = maxART - ART[ART >= 0] +1;
         }
         if (!any(SProb > 0)) {
            AFilePrint("ERROR ERROR");
            AFilePrint(paste("Hey, We sampled ", SampleAFunction, 
              " but ART is (", paste(ART, collapse=", "), ")", sep=""));
            AFilePrint(paste("  but SProb Sucks! (",
              paste(SProb, collapse=", "), ")", sep=""));
            ErrorText <- paste("Hey, We sampled ", SampleAFunction, 
              " but ART is (", paste(ART, collapse=", "), ")
              but SProb Sucks! (",
              paste(SProb, collapse=", "), ")", sep=""); 
            AErrorPrint(ErrorText);
            return(-1);             
         }
        SampleASimulation <- sample(TrySims, prob=SProb, size=1);
        try(.self$TheSummarySimulationList[[
          SampleASimulation]]$AddAFunction(as.integer(SampleAFunction)),
            SSL = .self);
           load(paste(.self$MyBaseSimulationStorageDirectory, "//SMS", 
           try(.self$TheSummarySimulationList[[
             SampleASimulation]]$UniqueSimulationIdentifier,".RData", sep="")));
           return(list(SMS=SMS, UniqueSimulationIdentifier=
             SMS$UniqueSimulationIdentifier,
             InSimulationsList=length(.self$TheSummarySimulationList), 
             OnFunction = .self$NameFunctions[SampleAFunction],
             OnNumFunction=SampleAFunction));
         }
           AFilePrint(paste("GetASimulation: How the Heck did we get ",   
             "Here to return -1!!", sep=""));
           AErrorPrint(paste("GetASimulation: How the Heck did we get ",   
             "Here to return -1!!", sep=""));
           return(-1);
       }
       AFilePrint(paste("GetASimulation: Null PerformanceMatrix ",
         ", How the Heck did we get ",   
         "Here to return Nothing!!", sep=""));
       AErrorPrint(paste("GetASimulation: Null PerformanceMatrix ",
         "How the Heck did we get ",   
             "Here to return Nothing!!", sep=""));
       return(-2);
  }     
);


GenerateANewSimulation <- function(AlreadyLocked=TRUE) {                                                                                                    
  eval(parse(text=GetG0Text("n", S=1)));
  eval(parse(text=GetG0Text("p", S=1)));
  eval(parse(text=GetG0Text("k", S=1)));

  eval(parse(text=GetG0Text("sigma", S=1)));
  eval(parse(text=GetG0Text("tNoiseDF", S=1)));
  eval(parse(text=GetG0Text("GenerateBetaVec", S=1)));
  eval(parse(text=GetG0Text("CorrelationXmatrix", S=1)));
  eval(parse(text=GetG0Text("tNoiseDF", S=1)));
  eval(parse(text=GetG0Text("LogitRegression", S=1)));
  eval(parse(text=GetG0Text("Beta0",S=1)));
  eval(parse(text=GetG0Text("ExperimentName", S=1)));
  eval(parse(text=GetG0Text("jobii", S=1)));
  eval(parse(text=GetG0Text("WorkingRow",S=1)));
  eval(parse(text=GetG0Text("TargetTotalSimulations",S=1)));
  eval(parse(text=GetG0Text("NameFunctions",S=1)));
  eval(parse(text=GetG0Text("OnSimType", S=1)));
  INeed = 1:length(NameFunctions);
  ISample = sample(INeed, size=1);
  SMS = NULL;
  if (verbose > 1) {
    AFilePrint("Look for a new simulation, had to configure a new SMS"); flush.console(); 
  }
  eval(parse(text=GetG0Text("ALargeContainDir")));
  eval(parse(text=GetG0Text("ASmallContainDir")));
  MySummaryPerformanceDirectory <-  paste(ALargeContainDir, 
  "//SumPerformance", sep="");

  if (OnSimType == "Group") {
    eval(parse(text=GetG0Text("GroupSize", S=1)));
    eval(parse(text=GetG0Text("pGroups", S=1)));
    eval(parse(text=GetG0Text("kGroups", S=1)));
    eval(parse(text=GetG0Text("tNoiseDF", S=1)));
    eval(parse(text=GetG0Text("sigma", S=1)));
    eval(parse(text=GetG0Text("SigmaNoise", S=1)));
    eval(parse(text=GetG0Text("Beta0", S=1)));
    eval(parse(text=GetG0Text("GenerateBetaGroupVec", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", S=1)));
    eval(parse(text=GetG0Text("LogitRegression", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", S=1)));
    eval(parse(text=GetG0Text("MeanCenterXX", S=1)));
    eval(parse(text=GetG0Text("jobii", S=1)));
    SMS <- SimGroupData(n = n, pGroups = pGroups, GroupSize=GroupSize,
     kGroups = kGroups, sigma = sigma,
     SigmaNoise = SigmaNoise, tNoiseDF = tNoiseDF, GenerateBetaGroupVec = GenerateBetaGroupVec, 
     CorrelationXmatrix = CorrelationXmatrix,
     LogitRegression = FALSE, Beta0 = Beta0, ExperimentName="", jobii= jobii, 
     WorkingRow = WorkingRow, AlreadyLocked=AlreadyLocked, MeanCenterXX = MeanCenterXX,
     UniqueProcessIdentifier = UniqueProcessIdentifier, ISample = ISample,
     DontSave=TRUE)
    
  } else if (OnSimType == "Robit") {
     try(SMS <- SimMeData(n = n, p = p, k = k, sigma = sigma,
     GenerateBetaVec = GenerateBetaVec, 
     CorrelationXmatrix = CorrelationXmatrix, tNoiseDF = tNoiseDF,
     LogitRegression = TRUE, Beta0 = Beta0, 
     ExperimentName=ExperimentName, jobii= jobii, WorkingRow = WorkingRow,
     AlreadyLocked=AlreadyLocked,
     UniqueProcessIdentifier=UniqueProcessIdentifier, ISample=ISample,
     DontSave = TRUE));
  } else if (OnSimType == "GroupRobit") {
    eval(parse(text=GetG0Text("GroupSize", S=1)));
    eval(parse(text=GetG0Text("pGroups", S=1)));
    eval(parse(text=GetG0Text("kGroups", S=1)));
    eval(parse(text=GetG0Text("tNoiseDF", S=1)));
    eval(parse(text=GetG0Text("sigma", S=1)));
    eval(parse(text=GetG0Text("SigmaNoise", S=1)));
    eval(parse(text=GetG0Text("Beta0", S=1)));
    eval(parse(text=GetG0Text("GenerateBetaGroupVec", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", S=1)));
    eval(parse(text=GetG0Text("LogitRegression", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", S=1)));
    eval(parse(text=GetG0Text("MeanCenterXX", S=1)));
    eval(parse(text=GetG0Text("jobii", S=1)));
    SMS <- SimGroupData(n = n, pGroups = pGroups, GroupSize=GroupSize,
     kGroups = kGroups, sigma = sigma,
     SigmaNoise = SigmaNoise, tNoiseDF = tNoiseDF, GenerateBetaGroupVec = GenerateBetaGroupVec, 
     CorrelationXmatrix = CorrelationXmatrix,
     LogitRegression = TRUE, Beta0 = Beta0, ExperimentName="", jobii= jobii, 
     WorkingRow = WorkingRow, AlreadyLocked=AlreadyLocked, MeanCenterXX = MeanCenterXX,
     UniqueProcessIdentifier = UniqueProcessIdentifier, ISample = ISample,
     DontSave=TRUE)
  } else {
  try(SMS <- SimMeData(n = n, p = p, k = k, sigma = sigma,
     GenerateBetaVec = GenerateBetaVec, 
     CorrelationXmatrix = CorrelationXmatrix, tNoiseDF = tNoiseDF,
     LogitRegression = LogitRegression, Beta0 = Beta0, 
     ExperimentName=ExperimentName, jobii= jobii, WorkingRow = WorkingRow,
     AlreadyLocked=AlreadyLocked,
     UniqueProcessIdentifier=UniqueProcessIdentifier, ISample=ISample,
     DontSave = TRUE));
   }
   if (AlreadyLocked == FALSE) {
     try(MOld <- getwd());
     try(setwd(CurrentSmallContainDir));
     while(LockMeIn(verbose=as.numeric(verbose), quoteMore=quoteMore, LFDir = "SumPerformance", 
       LockIn="SumLock")==FALSE){
       try(Sys.sleep(runif(1,0,4)));
     }
     try(setwd(MOld));
   }
   try(MOld <- getwd());
   try(setwd(CurrentLargeContainDir));  
   save(SMS = SMS, file=paste("SumPerformance", "//SMS", SMS$UniqueSimulationIdentifier, ".RData", sep=""));
   if (AlreadyLocked == FALSE) {
     try(setwd(CurrentSmallContainDir));
     try(UnLockMe(verbose=as.numeric(verbose), quoteMore=quoteMore, LFDir = "SumPerformance", 
       LockIn="SumLock"));
     try(setwd(MOld));
   }
   try(setwd(MOld));
   return(SMS);
}

##############################################################################
## HitASimulation
##
##   After a successful SMS simulation run on function NameFunction, this finds the directory
##    to save this UniqueSimulationIdentiifer setup and record the success.
##
HitAFunctionForSimulation <- function(UniqueProcessIdentifier, UniqueSimulationIdentifier, NameFunction, Hit = "+",
  quoteMore = "HitASimulation", verbose=0, TCS = NULL) {
  if (is.logical(verbose) && verbose == TRUE) {verbose = 1;}
  if (is.logical(verbose) && verbose == FALSE) {verbose =0;}
  if (verbose >= 2) {
    AFilePrint(paste("HitAFunctionForSimulation: Running for Name = ",
      NameFunction, " and UniqueProcess = ",
      UniqueProcessIdentifier, sep="")); flush.console();
    AFilePrint(paste(" --- with UniqueSimulationIdentifier = ", 
      UniqueSimulationIdentifier, sep=""));
    flush.console();
  }
  if (Hit!="+") { Hit = "D" 
    if (verbose >= 2) {
      AFilePrint(paste(
        "----  HitAFunctionForSimulation[", UniqueSimulationIdentifier, 
        "] --- OOF, Hit = ", Hit, " for NameFunction = ", NameFunction, 
        sep=""));
      flush.console();
    }
  }
  CurrentTotalCompleted <- 0;
  eval(parse(text=GetG0Text("CurrentTotalCompleted", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("CurrentLargeContainDir", S=1)));
  eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
  if (verbose >= 2) {
      AFilePrint(paste(
        "----  HitAFunctionForSimulation[", UniqueSimulationIdentifier, 
        "] --- CurrentLargeContainDir = ", CurrentLargeContainDir, 
        " for NameFunction = ", NameFunction, sep=""));
      flush.console();
  }
  MySummaryPerformanceDirectory = paste(
    CurrentLargeContainDir, "//", "SumPerformance", sep="");
  dir.create(MySummaryPerformanceDirectory, showWarnings=FALSE);
  try(quoteMore <- paste("Hit A Simulation, NameFunction = ", NameFunction,
    " and Hit = ", Hit, sep=""));
  if (verbose >= 1) {
      AFilePrint(paste(
        "----  HitAFunctionForSimulation[", UniqueSimulationIdentifier, 
        "] --- Perform LockMeIn", sep=""));
      flush.console();
  }
  MyTestT <- FALSE; 
  TryOutT <- "
    Oldwd <- getwd();
    try(setwd(CurrentSmallContainDir));
    while(LockMeIn(verbose=as.numeric(verbose), quoteMore=quoteMore, 
      LFDir = \"SumPerformance\", NoSpaceLFDir=TRUE, LockIn=\"SumLock\")==FALSE){
      MyTestT <- FALSE;
    }
    try(setwd(Oldwd));
    MyTestT <- TRUE;
  ";
  try(eval(parse(text=TryOutT)));
  if (MyTestT == FALSE) {
    AFilePrint(paste(
      "----  HitAFunctionForSimulation[", UniqueSimulationIdentifier, 
      "] --- Hey we tried to lock in but failed!  ",
      "I will return NULL now!", sep=""), LockIn="SumLock");
    return(NULL);
    flush.console();
  }
  if (verbose >= 1) {
    AFilePrint(paste("---- HitAFunctionForSimulation[", 
      UniqueSimulationIdentifier, 
       "] --- We have finished LockMeIn", sep=""));
  }
  if (verbose >= 1) {
      AFilePrint(paste("----  HitAFunctionForSimulation[", 
        UniqueSimulationIdentifier, 
        "] --- Completed LockMeIn. Get listfiles:", sep=""));
      flush.console();
  }
  eval(parse(text=SetGText("MySummaryPerformanceDirectory", 
    "globalenv()", S=1)));
  try(MyListFiles <- unlist(list.files(MySummaryPerformanceDirectory)));
  if (verbose >= 1) {
      AFilePrint(paste("----  HitAFunctionForSimulation[",
         UniqueSimulationIdentifier, 
         "] --- MyListFiles receieved length ", length(MyListFiles), sep=""));
      flush.console();
  }
  if (verbose >= 1) {
      AFilePrint(paste("----  HitAFunctionForSimulation[", 
        UniqueSimulationIdentifier, 
        "] --- Run GetSummaryPerformanceFile", sep=""));
      flush.console();
  }
  if (is.null(MyListFiles)) {
    AFilePrint(paste("HitAFunctionForSimulation: NOTE WARNING: ",
      "MyListfiles is NULL!", sep=""));
  }
  if (length(MyListFiles) <= 0) {
    AFilePrint(paste("HitAFunctionForSimulation: NOTE WARNING: ",
      "MyListfiles has no length!", sep=""));
  }
  if (is.null(MyListFiles) || length(MyListFiles)<= 0 || 
    !("SummaryPerformanceRData.RData" %in% MyListFiles)) {
    AFilePrint(paste("HitAFunctionForSimulation: Error Summary ",
      "PerformanceRData not in directory: how could you do that?", sep=""));
    AFilePrint("  We have that summary performance was: ");
    AFilePrint(MySummaryPerformanceDirectory);
    AFilePrint(paste("  and MyListFiles is: (", 
      paste(MyListFiles, collapse=", "), ")", sep=""));
    eval(parse(text=SetGText("MyListFiles", "globalenv()", S=1)));
    eval(parse(text=SetGText("MySummaryPerformanceDirectory", 
      "globalenv()", S=1)));
    AFilePrint(paste("  We cannot be running HitAFunctionForSimulation ",
      "without files.", sep=""));
    AErrorPrint(paste("HitAFunctionForSimulation: Error Summary ", 
      "PerformanceRData not in directory: how could you do that?", sep=""));
    tryCatch("HitAFunctionForSimulation");
  }
  if (verbose >= 1) {
    AFilePrint(paste(" Hit AFunction For Simulation: right now ", 
      "getting MySummarySimulationsList", sep=""));
  }
  try(MySummarySimulationsList <- NULL);
  eval(parse(text=SetGText("MySummarySimulationsList", S=1)));
  try(MyOldOldwd <- getwd());
  try(setwd(CurrentSmallContainDir));
  try(load(paste(
    "SumPerformance//SummaryPerformanceRData.RData", sep="")));
  try(eval(parse(text=TwoSimR5:::RecoveryTSSLText(CurrentLargeContainDir))));
  try(ISimID <- (1:length(MySummarySimulationsList$TheSummarySimulationList))[
    names(MySummarySimulationsList$TheSummarySimulationList) == 
    UniqueSimulationIdentifier]);
  try(setwd(MyOldOldwd));
  if (is.null(ISimID) || length(ISimID) != 1) {
    AFilePrint(paste("HitAFunctionForSimulation: after load, we don't ",
      "calculate ISimID for ", UniqueSimulationIdentifier, sep=""));
    AFilePrint("Will Paste Names of MySummarySimulationsList to global!");
    ATText <- "
    TheSummarySimulationList <- 
      MySummarySimulationsList$TheSummarySimulationList;
    namesTheSummarySimulationList <- names(TheSummarySimulationList);
    ";
    try(eval(parse(text=ATText)));
    eval(parse(SetGText("ISimID", "globalenv()", S=1)));
    eval(parse(text=SetGText("namesTheSummarySimulationList", 
      "globalenv()", S=1)));
    eval(parse(text=SetGText("TheSummarySimulationList", "globalenv()", S=1)));
    eval(parse(text=SetGText("UniqueSimulationIdentifier", 
      "globalenv()", S=1)));  
    eval(parse(text=SetGText("MySummarySimulationsList", "globalenv()", S=1)));  
    tryCatch("HitAFunctionForSimulation: We failed!");  
  }
  if (Hit == "+") {
    if (verbose >= 1) {
      AFilePrint(paste(" Hit AFunction about to run successful Complete ", 
        sep=""));
    }
    AR <- MySummarySimulationsList$TheSummarySimulationList[[
      ISimID]]$CompleteAFunction(ANameFunction=NameFunction, 
      SSL =MySummarySimulationsList);
    if (AR %in% c(10, 11) || AR != 1) {
       ReturnAHit = "L"; 
    } else {
      ReturnAHit = "+";
    }
  } else {
    if (verbose >= 1) {
      AFilePrint(paste(" Hit AFunction about to run a fail", 
        sep=""));
    }
    AR <- MySummarySimulationsList$TheSummarySimulationList[[
     ISimID]]$FailAFunction(ANameFunction=NameFunction,
       SSL = MySummarySimulationsList);
    ReturnAHit = "F";
  }
  if (verbose >= 3) {
    AFilePrint(paste("  --- HitASimulation finish, now we will",
      " Unlock.", sep="")); 
    flush.console();
  }
  try(eval(parse(text=TryFillCurrentTotalCompleted())));
  Oldwd <- getwd();
    eval(parse(text=GetG0Text("CurrentSmallContainDir", "globalenv()", S=1)));
    try(setwd(CurrentSmallContainDir));
    try(setwd("SumPerformance"));
    try(secure.save.MSS(AObject=MySummarySimulationsList, 
      ObjectName="MySummarySimulationsList",
      file=paste("SummaryPerformanceRData.RData", sep="")));
  try(setwd(CurrentSmallContainDir));
  UnLockMe(verbose=verbose, quoteMore=paste(
    quoteMore, " - HitASimulation", sep=""), 
    LFDir = "SumPerformance", LockIn="SumUnLock");
  try(setwd(Oldwd));
  if (verbose >= 3) {
    AFilePrint(paste("  --- HitASimulation: AllFinish.", sep="")); 
  }
  return(ReturnAHit);
}


TryFillCurrentTotalCompleted <- function() {  
  MyText <- "
  eval(parse(text=GetG0Text(\"CurrentTotalCompleted\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"ValidSimulations\", \"globalenv()\", S=1)));
  try(CurrentTotalCompleted <- MySummarySimulationsList$
    GetTotalCompletedFunction());
  try(ValidSimulations <- MySummarySimulationsList$ValidSimulations);
  eval(parse(text=SetGText(\"CurrentTotalCompleted\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"ValidSimulations\", \"globalenv()\", S=1)));
 
  eval(parse(text=SetGText(\"CurrentTotalKilled\", \"globalenv()\", S=1))); 
  try(CurrentTotalKilled <- MySummarySimulationsList$MyTotalKilledFunctions);
  eval(parse(text=SetGText(\"CurrentTotalKilled\", \"globalenv()\", S=1)));
  if (!exists(\"TCS\") || is.null(TCS)) {
    eval(parse(text=GetG0Text(\"TCS\", \"globalenv()\", S=1)));
  }
  if (exists(\"TCS\") && !is.null(TCS) && !is.numeric(TCS)) {
    try(TCS$CurrentTotalCompleted <- as.integer(CurrentTotalCompleted));
    try(TCS$ValidSimulations <- as.integer(ValidSimulations));
  }
  ";
}

TryFillCurrentTotalCompletedSELF <- function() {  
  MyText <- "
  eval(parse(text=GetG0Text(\"CurrentTotalCompleted\", \"globalenv()\", S=1)));
  try(CurrentTotalCompleted <- 
    MySummarySimulationsList$GetTotalCompletedFunction());
  eval(parse(text=GetG0Text(\"CurrentTotalKilled\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"CurrentTotalCompleted\", \"globalenv()\", S=1)));
    try(.self$CurrentTotalCompleted <- as.integer(CurrentTotalCompleted));
  eval(parse(text=GetG0Text(\"ValidSimulations\", \"globalenv()\", S=1)));
  try(ValidSimulations <-
    MySummarySimulationsList$ValidSimulations);
  eval(parse(text=SetGText(\"ValidSimulations\", \"globalenv()\",S=1)));
  try(CurrentTotalKilled <- MySummarySimulationsList$MyTotalKilledFunctions);
  eval(parse(text=SetGText(\"CurrentTotalKilled\", \"globalenv()\", S=1)));
  try(.self$ValidSimulations <- as.integer(ValidSimulations));
  ";
}

TryFillCurrentTotalCompletedTC <- function() {  
  MyText <- "
  eval(parse(text=GetG0Text(\"CurrentTotalCompleted\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"ValidSimulations\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"CurrentTotalKilled\", \"globalenv()\", S=1)));
  try(CurrentTotalCompleted <- TC);
  try(CurrentTotalKilled <- TK);
  try(ValidSimulations <- VLD);
  eval(parse(text=SetGText(\"CurrentTotalCompleted\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"CurrentTotalKilled\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"ValidSimulations\", \"globalenv()\", S=1)));
  if (!exists(\"TCS\") || is.null(TCS)) {
    eval(parse(text=GetG0Text(\"TCS\", \"globalenv()\", S=1)));
  }
  if (exists(\"TCS\") && !is.null(TCS) && !is.numeric(TCS)) {
    try(TCS$CurrentTotalCompleted <- as.integer(CurrentTotalCompleted));
    try(TCS$ValidSimulations <- as.integer(ValidSimulations));
  }
  ";
}