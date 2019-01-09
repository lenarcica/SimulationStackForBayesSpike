################################################################################
## 2LassoMaster.r
##
##
##  Master file for TwoSimR5.
##
##     (c) 2009-2019 Alan Lenarcic
##     The code was written to support simulations for Lenarcic and Valdar methods
##      paper.
##
##     This work begain in lab of Edoardo Airoldi to attempt simulations to rate
##     "2Lasso" penalty, but became a larger project for Bayesian and ArgMax penalties
##
##    2LassoMaster gives "AssessFirstTime" which is an algorithm run to get Estimators to fit.
##
##   The TwoSimR5 approach to conducint simulations is
##     1. Decide upon which N_M estimators to use in a problem, install functions
##      through EfficientSimulator and declare parameter settings with DeclareAllTstFunctions
##     2. Define a simulation (n,p,k, sigma, covarianceX, size of Beta...)
##     2. Simulate N_S 500-1000 indepenent simulations from t
##     3. Open N_M * (N_S/N_divisor) threads on Killdevil to attempt to solve each simulation with each estimator
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


##  AssessFirstTime is run by a RunSeekScriptMake.r file
##
##  This assess how many simulations are necessary, how many have yet
##   to be fit by algorithms, simulating more simulations or setting up
##   Processes to run algorithms as necessary.
##
AssessFirstTime <- function(ALargeContainDir, ASmallContainDir, TargetTotal = -1, TotalThralls=100, verbose = NULL) {
  if (is.null(verbose)) {
    eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
  }  
  
  eval(parse(text=GetG0Text("RenameFunctions", "globalenv()",S=1)));
  if (!is.null(RenameFunctions) && !is.numeric(RenameFunctions) && length(RenameFunctions) >= 1) {
    UpdateNameFunctions();
    eval(parse(text=GetG0Text("UseFunctions", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("NameFunctions", "globalenv()", S=1)));
  }
  if (verbose >= 1) {
    AFilePrint("AssessFirstTime:   StartAssessFirstTime  "); flush.console();
  }  
  try(eval(parse(text=StartAssessFirstTime())));

  if (verbose >= 1) {
    print(paste("AssessFirstTime: About to SetupMySimulations")); flush.console();
  }
  setwd(ASmallContainDir);
  eval(parse(text=SetupMySimulations()));
  
  if (verbose >= 1) {
    print(paste("AssessFirstTime: About to Update Length Sims")); flush.console();
  }
  try(eval(parse(text=UpdateLengthSims())));
  
  if (verbose >= 1) {
    print(paste("AssessFirstTime: About to SetupMySimulations 2nd time")); flush.console();
  }
  eval(parse(text=SetupMySimulations()));
      
  try(eval(parse(text=SetupSuccessKillDie())));
  try(eval(parse(text=CutMySimulations())));
 
  MyProcesses <- list.files("Processes");
  
  eval(parse(text=GetG0Text("NowDelete", "globalenv()", S=1)));
  if (is.null(NowDelete) || (is.numeric(NowDelete) && NowDelete[1] == 0)) {
    NowDelete = FALSE;
  }
  
  if (length(MyProcesses) >= 1) {
    if (verbose >= 1) {
      AFilePrint(paste("AssessFirst:  MyProcesses is length ", length(MyProcesses), 
        " perform end save process.", sep="")); flush.console();
    }
  for (jj in 1:length(MyProcesses)) {
    ZZMyID <- MyProcesses[jj];
    ListInProcess <- unlist(list.files(
      paste("Processes//", MyProcesses[jj], sep="") ));
    if (verbose >= 2) {
     print(paste("  We have jj = ", jj, ", and has length(ListInProcess) = ",
       length(ListInProcess), sep="")); flush.console();
    }
    if (length(ListInProcess) >= 1) {
      if ("MyProcesses.RData" %in% ListInProcess) {
        try(eval(parse(text=InTheMiddleProcessText())));
      } else {
        try(unlink(paste("Processes//", MyProcesses[jj], sep=""), recursive=TRUE));
      }
    }
    try(eval(parse(text=EndSaveProcess())));
  }
  }
 
  if (verbose >= 2) {
    AFilePrint("Setting Up Kill List"); flush.console();
  }
  try(eval(parse(text=SetupKillList())));

 
  if (verbose >= 0) {
    AFilePrint(paste("Setup MyCompleteProcessTable for ", length(NameFunctions), " functions to dir: ",
      ASmallContainDir, sep="")); flush.console();
  }
  MyCompletesProcessTable <- matrix(NA, TargetTotal, length(NameFunctions)); 
  try(rownames(MyCompletesProcessTable) <- MySimNoSMS[1:TargetTotal]);
  ##MyV = c(l1error, l2error, Type1, Type2, L.zero, time);
  MySuccessArray <-  NULL; 
  if (!exists("MySuccessTimesTable")) {
    MySuccessTimesTable <- matrix("", TargetTotal, length(NameFunctions)); 
  }
  try(rownames(MySuccessTimesTable) <- MySimNoSMS[1:TargetTotal]);
  try(colnames(MySuccessTimesTable) <- NameFunctions);
  for (ii in 1:length(NameFunctions)) {
    setwd(ASmallContainDir)
    try(eval(parse(text=SetupCompleteProcessTable())));
  }

  ## This will copy processes to Directories
  if (verbose >= 1) {
    AFilePrint("  --  SaveProcessWorkToGlobal");
  }
  try(eval(parse(text=SaveProcessWorkToGlobal())));

  if (verbose >= 1) {
    AFilePrint("  --  Ready Now for Thralls. ");
  }  
  try(eval(parse(text=ReadyForThralls())));
  if (length(ABB) <= 0) {
    AFilePrint("We Have completed, no ABB left to construct!");
    return(1);
  }
  AFill <- 0;  OGo <- LengthGo;
  if (verbose >= 2) {
    AFilePrint(paste("  About to setup SampleSort for ", TotalThralls, 
      " Total Thralls - and length(SampleSort) = ", length(SampleSort),
      sep="")); 
    AFilePrint(paste("  Starting with LengthGo = ", LengthGo, sep=""));
  }
  for (ii in 1:TotalThralls) {
     if (AFill+LengthGo > length(SampleSort)) {
        OGo <- length(SampleSort)-AFill;
     }
     MyKeepSampleSort <-  SampleSort[AFill + 1:OGo];
     DoSMS <- MyKeepSampleSort %% TargetTotal+1;
     DoNameFunc <-   (MyKeepSampleSort-1) %/% TargetTotal+1;
     SubMitList <- cbind(MySimNoSMS[DoSMS], NameFunctions[DoNameFunc],
       rep(0, length(DoSMS)));
     MyRequest <- paste("Request",
        AZeroOut(as.numeric(ii),10^ceiling(log(TotalThralls+1,10))), ".RData", sep="");
     eval(parse(text=GetG0Text("ASmallContainDir")));
     MyTryText <- "
       GoodGo = 0;
       OldOnWd <- getwd();
       setwd(ASmallContainDir);
       GoodGo = 1;
     ";
     try(eval(parse(text=MyTryText)));
     if (GoodGo == 0) {
        print("***********************************************************************  Error");
        print("**  Error ")
        print("** We have an Error trying to set directory to ASmallContainDir, how will we create Submission!");
        flush.console();
        retun(-999);
     }
     CreateDirText = "
     GoodGo = 0;
     dir.create(\"Submissions\", showWarnings=FALSE, recursive=TRUE);
     GoodGo = 1;
     ";
     try(eval(parse(text=CreateDirText)));
     if (GoodGo == 0) {
       print("One Fail on create Submissions directory, taking a pause."); flush.console();
       Sys.sleep(1);
       try(eval(parse(text=CreateDirText)));
     }
     if (GoodGo == 0) {
        print("*******************************************************************************");
        print("** Error ");
        print("** We were in the right directory but couldn't create a Submissions directory!"); flush.console();
        return(-999);
     }
     
     SaveDirectoryText = "
     GoodGo = 0;
     CurrentLargeContainDir <- ALargeContainDir;
     CurrentSmallContainDir <- ASmallContainDir;   
     save(SubMitList=SubMitList, ASmallContainDir=ASmallContainDir, 
       ALargeContainDir=ALargeContainDir, 
       CurrentLargeContainDir=CurrentLargeContainDir, 
       CurrentSmallContainDir=CurrentSmallContainDir, file=paste(\"Submissions//\", MyRequest, sep=\"\"));
     GoodGo = 1;
     ";
     try(eval(parse(text=SaveDirectoryText)));
     if (GoodGo == 0) {
        print("**************************************************************************************");
        print("  Error ");
        print("** We tried to save SubMitList to a directory and failed."); flush.console();
        return(-999);
     }
     AFill <- AFill + OGo;
     MyLS <-rep(paste("R --vanilla --args ", argeese, " ", ii, " 1 1 1 1 1 1 1  ", " < ", RunDirFile, sep=""),
       OGo);
     AFF <- NULL;
     try(AFF <- file(paste(AllDirToSave, "/RBayesAlgorithm", MyNameOfInputFile, ii, ".txt", sep="")));
     AFilePrint(paste("We saved ", OGo, " functions to ", AllDirToSave, 
       "/RBayesAlgorithm", MyNameOfInputFile, ii, ".txt", sep=""));
     if (!is.null(AFF)) {
        writeLines(text=MyLS, con=AFF);
        try(close(AFF));
     }

  }
  print("At the End, please use bsub line : ");
  print(paste("bash ", AllDirToSave, "/RBayesAlgorithm", MyNameOfInputFile, "${LSB_JOBINDEX}.txt", sep=""));
  eval(parse(text=GetG0Text("OriginalOldwd")));
  try(setwd(OriginalOldwd));
}



StartAssessFirstTime <-  function() {
  ART <- " 
  if (!exists(\"TotalThralls\") || is.null(TotalThralls) || TotalThralls <= 0) {
    TotalThralls <- 5;
  }
  if (!exists(\"TargetTotal\") || is.null(TargetTotal) || TargetTotal <= 0) {
     eval(parse(text=GetG0Text(\"TargetTotal\", \"globalenv()\", S=1)))
  }
  if (is.null(TargetTotal) || TargetTotal <= 0) {
    TargetTotal <- 20;
  }
  MOld <- getwd();
  
  eval(parse(text=GetG0Text(\"DuplicationEETest\", \"globalenv()\", S=1)));
  if (is.logical(DuplicationEETest) && DuplicationEETest == TRUE) {
    OnSimType <- \"DuplicationEE\";
  }  else {
    DuplicationEETest <- FALSE;
    eval(parse(text=SetGText(\"DuplicationEETest\", \"globalenv()\", S=1)));
  }
  AllDirToSave <- MakeDirSaves();
  KillList <- NULL; 
  
  if (exists(\"DuplicationEETest\") && is.numeric(DuplicationEETest) && DuplicationEETest[1] == 0) {
    DuplicationEETest <- FALSE;
    eval(parse(text=SetGText(\"DuplicationEETest\", \"globalenv()\", S=1)));
  } 
  if (is.logical(DuplicationEETest) && DuplicationEETest == TRUE) {
    OnSimType <- \"DuplicationEE\";
    eval(parse(text=SetGText(\"OnSimType\", \"globalenv()\",S=1)));
    eval(parse(text=GetG0Text(\"RenameFunctions\", \"globalenv()\",S=1)));
    #NameFunctions <- c(DefaultAllNameFunctions);
    #, 
    #    c(\"GenerateGroupBayesSpikeAutoTemp1R\",
    #    \"GenerateGroupBayesSpikeAutoTemp2R\",
    #    \"GenerateGroupBayesSpikeAutoTemp3R\",
    #    \"GenerateGroupBayesSpikeAutoTemp4R\"
    #    )
    #);
    #UseFunctions <- c(DeclareAllDefaultFunctions,
    #        c(GenerateGroupBayesSpikeAutoTemp1R,
    #    GenerateGroupBayesSpikeAutoTemp2R,
    #    GenerateGroupBayesSpikeAutoTemp3R,
    #    GenerateGroupBayesSpikeAutoTemp4R
    #    )
    #);
    ##eval(parse(text=SetGText(\"NameFunctions\", \"globalenv()\", S=1)));
    ##eval(parse(text=SetGText(\"UseFunctions\", \"globalenv()\", S=1)));    
    if (!(\"GenerateBayesSpikeAutoTemp1R\" %in% RenameFunctions)) {
      RenameFunctions <- c(\"GenerateBayesSpikeAutoTemp1R\",
        \"GenerateBayesSpikeAutoTemp2R\",
        \"GenerateBayesSpikeAutoTemp3R\",
        \"GenerateBayesSpikeAutoTemp4R\"
        );
      eval(parse(text=SetGText(\"RenameFunctions\", \"globalenv()\", S=1)));
    }
      
  }

  eval(parse(text=GetG0Text(\"UniqueProcessIdentifier\", \"globalenv()\", S=1)));
  BackProcessIdentifier <- UniqueProcessIdentifier;
  eval(parse(text=SetGText(\"BackProcessIdentifier\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"TooManyFails\", \"globalenv()\", S=1)));
  
  try(dir.create(paste(ALargeContainDir, \"//SMS\", sep=\"\"), showWarnings=FALSE, recursive=TRUE));
  try(dir.create(paste(ASmallContainDir, \"//SmallSMS\", sep=\"\"), showWarnings=FALSE, recursive=TRUE));
  setwd(ASmallContainDir);
  
  eval(parse(text=GetG0Text(\"RenameFunctions\", \"globalenv()\", S=1)));
  if (length(RenameFunctions) >= 1 && (is.character(RenameFunctions) ||
    (is.numeric(RenameFunctions[1]) && RenameFunctions[1] != 0.0))) {
    print(paste(\"2LassoMaster.r: We have a RenameFunctions so we are updating NameFunctions to length \", 
      length(RenameFunctions), sep=\"\")); flush.console();
    UpdateNameFunctions(); 
    print(paste(\"2LassoMaster.r: UpdateNameFunctions() has returned. \", sep=\"\")); flush.console();            
  }
  ";
  return (ART);
}

UpdateLengthSims <- function() {
  ART <- "

  MySimNoSMS <- MySimNoSMS[MySimNoSMS != \"\"];
  if (length(MySimNoSMS) < TargetTotal) {
    if (verbose >= 1) {
      AFilePrint(paste(\" We are going to have to create \", 
        TargetTotal - length(MySimNoSMS), \", because Target Total = \",
        TargetTotal, \"  and length(MySimNoSMS) = \", length(MySimNoSMS), \"\"));
    }
    for (ii in (length(MySimNoSMS)+1):TargetTotal) {
       SMS <- MasterGenerateANewSimulation();  
       try(setwd(ALargeContainDir));  
       try(SMS$UniqueSimulationIdentifier <- paste(SMS$UniqueSimulationIdentifier, \"N\",
         AZeroOut(as.numeric(ii),10^ceiling(log(TargetTotal+1,10))), sep=\"\"));
       try(SMS$SimulationIdentifier <- SMS$UniqueSimulationIdentifier);
       try(dir.create(paste(\"SMS//SMS\", SMS$UniqueSimulationIdentifier, sep=\"\"), showWarnings=FALSE, recursive=TRUE));
       ##AFilePrint(paste(\"Create and save to directory: \",
       ##  SMS$UniqueSimulationIdentifier, sep=\"\"))
       setwd(paste(\"SMS//SMS\", SMS$UniqueSimulationIdentifier, sep=\"\"));
       try(save(SMS = SMS, file=paste(\"\", SMS$UniqueSimulationIdentifier, \".RData\", sep=\"\")));
       try(setwd(ASmallContainDir));
       try(dir.create(\"SmallSMS\", recursive=TRUE, showWarnings=FALSE));
       SimulationIdentifier <- -1;
       try(SimulationIdentifier <- SMS$UniqueSimulationIdentifier);
       ListSimsParams <- c(n=length(SMS$Y), p=NCOL(SMS$X), puse = SMS$puse);
       try(save(SimulationIdentifier=SimulationIdentifier, 
         ListSimsParams=ListSimsParams,
         file=paste(\"SmallSMS//SMS\", SimulationIdentifier, \".RData\", sep=\"\")));
    }    
  } else {
    if (verbose >= 1) {
      AFilePrint(paste(\" We are not going to create \", 
        TargetTotal - length(MySimNoSMS), \", because Target Total = \",
        TargetTotal, \"  and length(MySimNoSMS) = \", length(MySimNoSMS), \"\"));
    }
  }
  
  setwd(ASmallContainDir);
  dir.create(\"Processes\", showWarnings=FALSE, recursive=TRUE);
  dir.create(\"SubmissionsProcesses\", showWarnings=FALSE, recursive=TRUE);
  eval(parse(text=GetG0Text(\"NameFunctions\")));
";
return(ART);
}

SetupSuccessKillDie <- function() {
  ART <- "
  try(setwd(ASmallContainDir));
  for (ii in 1:length(NameFunctions)) {
    try(dir.create(NameFunctions[ii], showWarnings=FALSE, recursive=TRUE));
  }
  KillList <- NULL;
  dir.create(\"SumPerformance\", showWarnings=FALSE, recursive=TRUE);
  if (any(unlist(list.files(\"SumPerformance\")) == \"MyFails.RData\")) {
    if (verbose >= 1) {
      AFilePrint(\"  --  We are Loading from MyFails.RData\");
    }
    try(load(paste(\"SumPerformance//MyFails.RData\", sep=\"\")));
    if (!exists(\"MyFails\") || !exists(\"MyFailList\")) {
    if (!exists(\"MyFails\")) {
        print(paste(\" We attempted to load MyFails, but it did not exist, regenerate. \"));
        flush.console();
        MyFails <- rep(0, length(NameFunctions));    
        names(MyFails) <- NameFunctions;
    }
    if (!exists(\"MyFailList\")) {
        MyFailList <- NULL;
    }
    save(MyFails=MyFails, MyFailList=MyFailList, file=paste(\"SumPerformance//MyFails.RData\", sep=\"\"));
    }
  } else {
    if (verbose >= 1) {
      AFilePrint(\"  --  We are Creating MyFails.RData\");
    }
    MyFails <- rep(0, length(NameFunctions));
    names(MyFails) <- NameFunctions;
    MyFailList <- NULL;
    save(MyFails=MyFails, MyFailList=MyFailList, file=paste(\"SumPerformance//MyFails.RData\", sep=\"\"));
  }
  if (any(unlist(list.files(\"SumPerformance\")) == \"MySuccesses.RData\")) {
    try(load(paste(\"SumPerformance//MySuccesses.RData\", sep=\"\")));
    if (!exists(\"MySuccessProcessTable\") || is.null(MySuccessProcessTable)) {
      if (verbose >= 2) {
         AFilePrint(\"  --  We are Creating MySuccessProcess Table. \");
      }
      MySuccessProcessTable <- matrix(\"\", TargetTotal, length(NameFunctions)); 
    }
    if (!exists(\"MySuccessTimesTable\")  | !exists(\"MySuccesses\")) {
    if (!exists(\"MySuccessTimesTable\") || is.null(MySuccessTimesTable)) {
      MySuccessTimesTable <- matrix(\"\", TargetTotal, length(NameFunctions)); 
    } else if (exists(\"MySuccessTimesTable\") && NROW(MySuccessTimesTable) < TargetTotal) {
      MySuccessTimesTable <- rbind(MySuccessTimesTable,   
        matrix(0,TargetTotal-NROW(MySuccessTimesTable), NCOL(MySuccessTimesTable)))
    }
    if (!exists(\"MySuccesses\") || is.null(MySuccesses)) {
      MySuccesses <- rep(0, length(NameFunctions));
      names(MySuccesses) <- NameFunctions;
    }
    save(MySuccesses=MySuccesses, MySuccessProcessTable=MySuccessProcessTable, 
      MySuccessTimesTable = MySuccessTimesTable,
      file=paste(\"SumPerformance//MySuccesses.RData\", sep=\"\"));
    }
  } else {
    MySuccesses <- rep(0, length(NameFunctions));
    names(MySuccesses) <- NameFunctions;
    MySuccessProcessTable <- matrix(\"\", TargetTotal, length(NameFunctions));
    MySuccessTimesTable <- matrix(0, TargetTotal, length(NameFunctions));
    MySimNoSMS <- MySimNoSMS[MySimNoSMS != \"\"];
    try(rownames(MySuccessProcessTable) <- MySimNoSMS[1:TargetTotal]);
    try(colnames(MySuccessProcessTable) <- NameFunctions); 
    save(MySuccesses=MySuccesses, MySuccessProcessTable=MySuccessProcessTable, 
      MySuccessTimesTable = MySuccessTimesTable,
      file=paste(\"SumPerformance//MySuccesses.RData\", sep=\"\"));
  }
  if (any(unlist(list.files(\"SumPerformance\")) == \"MyKills.RData\")) {
    try(load(paste(\"SumPerformance//MyKills.RData\", sep=\"\")));
    if (!exists(\"MyKills\") || !exists(\"MyKillList\")) {
    if (!exists(\"MyKills\")) {
    MyKills <- rep(0, length(NameFunctions));
    names(MyKills) <- NameFunctions        
    }
    if (!exists(\"MyKillList\")) { MyKillList <- NULL; }
      try(save(MyKills=MyKills, MyKillList=MyKillList, KillList=KillList, 
        file=paste(\"SumPerformance//MyKills.RData\", sep=\"\")));
    }
  } else {
    MyKills <- rep(0, length(NameFunctions));
    names(MyKills) <- NameFunctions;
    MyKillList <- NULL;
    save(MyKills=MyKills, MyKillList=MyKillList, KillList=KillList, file=paste(\"SumPerformance//MyKills.RData\", sep=\"\"));
  }
  if (any(unlist(list.files(\"SumPerformance\")) == \"MyDies.RData\")) {
    try(load(paste(\"SumPerformance//MyDies.RData\", sep=\"\")));
    if (!exists(\"MyDies\") || !exists(\"MyDieList\")) {
      MyDies <- rep(0, length(NameFunctions));
      names(MyDies) <- NameFunctions;
      MyDieList <- NULL;
      save(MyDies=MyDies, MyDieList=MyDieList, file=paste(\"SumPerformance//MyDies.RData\", sep=\"\"));
    }
  } else {
    MyDies <- rep(0, length(NameFunctions));
    names(MyDies) <- NameFunctions;
    MyDieList <- NULL;
    save(MyDies=MyDies, MyDieList=MyDieList, file=paste(\"SumPerformance//MyDies.RData\", sep=\"\"));
  }
  try(MySims <- substr(MySimulations, nchar(\"SMS\")+1, nchar(MySimulations)));
  if (verbose >= 1) {
      AFilePrint(\"  --  We are done with SetupSuccessKillDie. \");
  }
  ";
  return(ART);
    
}

InTheMiddleProcessText <- function() {
    ART <- "
        if (verbose >= 3) {
          AFilePrint(paste(\"   We start InTheMiddleProcessText about to load \",
            \"MyProcesses.RData for Process[\", jj, \"] named \",
             MyProcesses[jj], sep=\"\"));
        }
        try(load(paste(\"Processes//\", MyProcesses[jj], \"//MyProcesses.RData\", sep=\"\")));
        TableProgressNames <- StackProgresses(ProgressNames)
        TableProgressData <- StackProgresses(ProgressData);
        MySFiles <- ListInProcess[substr(ListInProcess,1,1)==\"S\"];
        MySFiles <- MySFiles[!is.na(MySFiles)];
        OtherTableProgressNames <- matrix(\"N\", NROW(TableProgressNames),
          NCOL(TableProgressNames));
        OtherTableProgressData <- matrix(0, NROW(TableProgressNames),
          NCOL(TableProgressData));
        if (length(MySFiles) >= 1) {
        for (iti in 1:length(MySFiles)) {
          MyV <- NULL;
          NameFunction <- NULL;  FunctionName <- NULL;
          try(rm(NameFunction));  try(rm(FunctionName));
          UniqueSimulationIdentifier <- NULL;
          try(rm(UniqueSimulationIdentifier));
          SimulationIdentifier <- NULL;
          try(rm(SimulationIdentifier));
          RunProcessIdentifier <- NULL;
          try(rm(RunProcessIdentifier));
          try(load(paste(\"Processes//\", MyProcesses[jj], \"//\", MySFiles[iti], sep=\"\")));
          if (!exists(\"NameFunction\") || is.null(NameFunction)) {
            if (!exists(\"FunctionName\")) {
              print(paste(\"Whoa: Before we were goint to work on \", 
                \"file MySFiles[iti=\", iti, \"] = \", MySFiles[iti], 
                \" we don't have Namefunction or FunctionName. \", sep=\"\"));
              flush.console();
            }
            try(NameFunction <- FunctionName);
          }
          if (!exists(\"UniqueSimulationIdentifier\") || is.null(UniqueSimulationIdentifier)) {
            if (!exists(\"SimulationIdentifier\")) {
              print(paste(\"Whoa: Before we were goint to work on \", 
                \"file MySFiles[iti=\", iti, \"] = \", MySFiles[iti], 
                \" we don't have SimulationIdentifier or UniqueSimulationIdentifier. \", sep=\"\"));
              flush.console();
            }
            try(UniqueSimulationIdentifier <- SimulationIdentifier);
          }
          if (!exists(\"MyV\") || is.null(MyV)) {
            AFilePrint(paste(\"Hey Bad load of MyProcesses[\", jj, \"]: SFile = \",
              MySFiles[ii], sep=\"\"));
          }  else {
            ABY <- (1:NROW(TableProgressNames))[TableProgressNames[,1] == SimulationIdentifier & 
              TableProgressNames[,2] == FunctionName];
            if (is.na(MyV[1])) { AttemptSucceedFail = -1;} else {
               AttemptSucceedFail = 1; 
            }
            if (length(ABY) == 1) {
              OtherTableProgressNames[ABY, ] <- c(UniqueSimulationIdentifier, FunctionName);
              OtherTableProgressData[ABY, ] <- as.numeric(c( 
              TableProgressData[ABY,1:2],
              AttemptSucceedFail, MyV[6:8],
              L2Accuracy=MyV[2], TableProgressData[ABY,8:9]))
            }
            eval(parse(text=GetG0Text(\"UniqueProcessIdentifier\", \"globalenv()\", S=1)));
            if (!exists(\"RunProcessIdentifier\") && exists(\"RunProcessIdentifer\")) {
              RunProcessIdentifier <- RunProcessIdentifer;
            } else if (!exists(\"RunProcessIdentifier\") && exists(\"UniqueProcessIdentifier\")) {
              RunProcessIdentifier <- UniqueProcessIdentifier;
            }
            if (exists(\"MyV\") && !is.na(MyV[1])) {
              NameFunction <- FunctionName;
              if (!exists(\"StartsTime\")) { StartsTime <- GetRealTime(); }
              if (!exists(\"FinishesTime\")) { FinishesTime <- GetRealTime(); }            
              save(MyV=MyV, RunProcessIdentifier = RunProcessIdentifier, 
               UniqueSimulationIdentifier=SimulationIdentifier, NameFunction=NameFunction,
               FunctionName=FunctionName,
               NowToSaveTime=NowToSaveTime, FinishesTime=FinishesTime, StartsTime=StartsTime,
               file=paste(NameFunction, \"//\",
               \"//S\", FunctionName, \"S\", SimulationIdentifier, \".RData\", sep=\"\"));
            }
            try(NowToSaveTime <- NULL);
            try(rm(NowToSaveTime));
             
            if (!exists(\"SimulationIdentifier\")) {
              print(\"Woah, No way are we going to achieve here with simulation identifier not exists!\");
              flush.console();
            }
            if (paste(\"CI\", NameFunction, \"S\", SimulationIdentifier, \".RData\", sep=\"\") %in% 
               ListInProcess) {
              AllCIReturn <- NULL;
              try(load(paste(\"Processes//\", MyProcesses[jj], \"//\",
                \"CI\", NameFunction, \"S\", SimulationIdentifier, \".RData\", sep=\"\")));
              if (!is.null(AllCIReturn)) {
               try(save(AllCIReturn=AllCIReturn, NameFunction=NameFunction, 
                 FunctionName = FunctionName,
                 Betas=Betas, CIQuantiles=CIQuantiles, SimulationIdentifier=SimulationIdentifier,
                 CIEst = CIEst, Hit=Hit, CITime=CITime, 
                 RunProcessIdentifier=RunProcessIdentifier,
                 file=paste(FunctionName, \"//\",
                 \"//CI\", FunctionName, \"S\", SimulationIdentifier, 
                 \".RData\", sep=\"\")));
              } 
              rm(AllCIReturn); rm(Betas);  rm(CIQuantiles); rm(CIEst); rm(Hit); rm(CITime);
            }
            if (paste(\"BetaCI\", NameFunction, \"S\", SimulationIdentifier, \".RData\", sep=\"\") %in% 
               ListInProcess) {
              AllCIReturn <- NULL;
              try(load(paste(\"Processes//\", MyProcesses[jj], \"//\",
                \"BetaCI\", NameFunction, \"S\", SimulationIdentifier, \".RData\", sep=\"\")));
              if (!is.null(AllCIReturn)) {
               if (!exists(\"IndexOnBetas\")) {IndexOnBetas <- NULL; }
               if (!exists(\"BetasOn\")) { BetasOn <- NULL; }
               if (!exists(\"MIPBetaOn\")) { MIPBetaOn <- NULL; }
               if (!exists(\"MeanMIPBetaOff\")) { MeanMIPBetaOff <- NULL; }
               if (!exists(\"TemperatureList\")) { TemperatureList <- NULL; }
               if (!exists(\"IndexOnBetaStart\")) { IndexOnBetaStart <- NULL; }
               if (!exists(\"BetaStartOn\")) { BetaStartOn <- NULL; }   
               if (!exists(\"NoNoiseBetaStart\")) { NoNoiseBetaStart <- NULL; }    
               if (!exists(\"CIQuantiles\")) { CIQuantiles <- NULL; }   
               if (!exists(\"CIEst\")) { CIEst <- NULL; }            
               try(save(AllCIReturn=AllCIReturn, NameFunction=NameFunction, 
                FunctionName = FunctionName, 
                IndexOnBetas = IndexOnBetas, 
                BetasOn=BetasOn, MIPBetaOn = MIPBetaOn, MeanMIPBetaOff = MeanMIPBetaOff, 
                TemperatureList=TemperatureList, IndexOnBetaStart=IndexOnBetaStart, BetaStartOn=BetaStartOn,
                NoNoiseBetaStart=NoNoiseBetaStart, 
                CIQuantiles=CIQuantiles,
                CIEst = CIEst, Hit=Hit, CITime=CITime, RunProcessIdentifier = RunProcessIdentifier,
                SimulationIdentifier = SimulationIdentifier,
                 file=paste(FunctionName, \"//\",
                 \"//BetaCI\", FunctionName, \"S\", SimulationIdentifier, 
                 \".RData\", sep=\"\")));
              } 
              rm(AllCIReturn); rm(Betas);  rm(CIQuantiles); rm(CIEst); rm(Hit); rm(CITime);
            }
            if (paste(\"MIP\", NameFunction, \"S\", SimulationIdentifier, \".RData\", sep=\"\") %in% 
               ListInProcess) {
              AllMIPReturn <- NULL;
              try(load(paste(\"Processes//\", MyProcesses[jj], \"//\",
                \"MIP\", NameFunction, \"S\", SimulationIdentifier, \".RData\", sep=\"\")));
              if (!is.null(AllMIPReturn)) {
                try(save(AllMIPReturn=AllMIPReturn, NameFunction=NameFunction, 
                  FunctionName=FunctionName, SimulationIdentifier=SimulationIdentifier, 
                  RunProcessIdentifier=RunProcessIdentifier, 
                  file=paste(FunctionName, \"//MIP\", NameFunction, \"S\", SimulationIdentifier, \".RData\", sep=\"\")));
              } 
              rm(AllMIPReturn);  rm(Hit);
            }
            MyV <- NULL;  try(rm(MyV));
          }
        } 
        }  
        if (exists(\"NowDelete\") && NowDelete == TRUE) {
          for (ii in 1:length(MySFiles)) {
            try(unlink(paste(\"Processes//\", MyProcesses[jj], \"//\", MySFiles[ii], sep=\"\"), recursive=TRUE));
          }
          try(unlink(paste(\"Processes//\", MyProcesses[jj], sep=\"\"), recursive=TRUE));
        }  
        for (itTPN in 1:NROW(TableProgressNames)) {
          if (TableProgressData[itTPN,3] >= 1) {
            try(APR <- (1:NROW(MySuccessProcessTable))[TableProgressNames[itTPN,1] == rownames(MySuccessProcessTable)
              | paste(\"SMS\",TableProgressNames[itTPN,1], sep=\"\") == rownames(MySuccessProcessTable) ]);
            try(APC <- (1:NCOL(MySuccessProcessTable))[TableProgressNames[itTPN,2] == colnames(MySuccessProcessTable)]);
            if (length(APR) >= 1) {
              try(MySuccessProcessTable[APR[1],APC[1]] <- RunProcessIdentifier);
              try(MySuccessTimesTable[APR[1],APC[1]] <- TableProgressData[itTPN,5]);
              try(MySuccesses[APC[1]] <- MySuccesses[APC[1]]+1)
            }
          } else if (TableProgressData[itTPN,3] == 0) {
            ABad <- (1:length(NameFunctions))[NameFunctions==TableProgressNames[itTPN,2]]
            try(MyDies[ABad] <- MyDies[ABad]+1);
            if (is.null(MyDieList) || length(MyDieList) == 0) {
              MyDieList <- list();
            }
            LT <- length(MyDieList)+1;
            MyDieList[[LT]] <- c(TableProgressNames[itTPN,], RunProcessIdentifier, 
              TableProgressData[itTPN,1:3]);
          } else if (TableProgressData[itTPN,3] < 0) {
            ABad <- (1:length(NameFunctions))[NameFunctions==TableProgressNames[itTPN,2]]
            try(MyFails[ABad] <- MyFails[ABad]+1);
            if (is.null(MyFailList) || length(MyFailList) == 0) {
              MyFailList <- list();
            }
            LT <- length(MyFailList)+1;
            MyFailList[[LT]] <- c(TableProgressNames[itTPN,], RunProcessIdentifier, 
              TableProgressData[itTPN,1:3]);
          }   
        }      
    ";
    return(ART);
}

EndSaveProcess <- function() {
  ART <- "
   setwd(ASmallContainDir);
    if (!exists(\"MyFails\")) {
        print(\"Really, no MyFails at this point, why?\"); flush.console(); 
    } else {
      save(MyFails=MyFails, MyFailList=MyFailList, file=paste(\"SumPerformance//MyFails.RData\", sep=\"\"));
    }
    save(MySuccesses=MySuccesses, MySuccessProcessTable=MySuccessProcessTable, 
      MySuccessTimesTable = MySuccessTimesTable,
      file=paste(\"SumPerformance//MyFails.RData\", sep=\"\"));
    if (!exists(\"MyKillList\")) { MyKillList <- NULL; }
    save(MyKills=MyKills, MyKillList=MyKillList, KillList=KillList, file=paste(\"SumPerformance//MyKills.RData\", sep=\"\"));
    save(MyDies=MyDies, MyDieList=MyDieList, file=paste(\"SumPerformance//MyDies.RData\", sep=\"\"));     
  ";
  return(ART);  
    
}

SetupKillList <- function() {
  ART <- "
   if (any(MyKills+MyDies >= TooManyFails)) {
    AFilePrint(\"We have many MyKills, MyDies, the rejectable are: \");
    AFilePrint(paste(\" Rejectable: (\",
      paste(NameFunctions[MyKills+MyDies >= TooManyFails], collapse=\", \"),
      \")\", sep=\"\"));
    KillList <- NameFunctions[MyKills+MyDies >= TooManyFails]
  } else {
    if (!exists(\"KillList\")) { KillList <- NULL; }
    KillList=KillList;
  }
  ";
  return(ART);  
}

SetupCompleteProcessTable <- function() {
  ART <- "
  
    setwd(ASmallContainDir)
    MyS <-  list.files(NameFunctions[ii]);
    MyS <- MyS[substr(MyS,1,1) == \"S\"];
    MyS <- MyS[substr(MyS,nchar(MyS) - nchar(\".RData\")+1, nchar(MyS)) == \".RData\"]
    AFilePrint(paste(\"On NameFunction: \", NameFunctions[ii], \" with \", 
      length(MyS), \" files. \", sep=\"\"));
      CAdd <- 0;
    
    MySMSSolved <- substr(MyS, nchar(NameFunctions[ii])+3,nchar(MyS));
    AFF <- substr(MySMSSolved, nchar(MySMSSolved)-nchar(\".RData\")+1,
      nchar(MySMSSolved)) == \".RData\";
    MySMSSolved[AFF] <- substr(MySMSSolved[AFF], 1, nchar(MySMSSolved[AFF])-nchar(\".RData\"));
    
    if (length(MySMSSolved) >= 1) {
    MyCompletesProcessTable[rownames(MyCompletesProcessTable) %in% MySMSSolved,ii] <- 1;
    for (jj in 1:length(MySMSSolved)) {
      RunProcessIdentifier <- NULL;  
      SimulationIdentifier <- NULL;
      UniqueSimulationIdentifier <- NULL;
      setwd(paste(ASmallContainDir, \"//\", NameFunctions[ii], sep=\"\"));
      try(MyV <- NULL);
      eval(parse(text=SetGText(\"MyV\", \"globalenv()\", S=1)));
      try(rm(MyV));
      try(load(paste( \"S\", NameFunctions[ii], \"S\", MySMSSolved[jj], \".RData\", sep=\"\")));   
      if (!exists(\"MyV\") || is.null(MyV) || !is.numeric(MyV)) {
        AFilePrint(paste(\"Error: we loaded MyS[jj] = \", MyS[jj], \" for jj = \", jj, sep=\"\"));
        AFilePrint(paste(\"  But nothing hit, ii=\", ii, \" and NameFunction = \",
          NameFunctions[ii], sep=\"\"));
      } 
      SimLocation <- (1:NROW(MyCompletesProcessTable))[SimulationIdentifier ==
        rownames(MyCompletesProcessTable)];
      if (is.null(MySuccessArray)) {
        print(paste(\"About to create MySuccess array, note TargetTotal = \", 
          TargetTotal, \" name functions are len \", length(NameFunctions),sep=\"\")); flush.console();
        MySuccessArray <- array(-1, dim=c(TargetTotal, length(NameFunctions),8),
          dimnames = list(Sim=1:TargetTotal, Function=NameFunctions, Result=1:8)); 
      }
      eval(parse(text=SetGText(\"MyV\", \"globalenv()\", S=1)));
      if (exists(\"MyV\") && !is.null(MyV)) {
        MySuccessArray[SimLocation,ii,] <- MyV;
        MySuccessTimesTable[SimLocation,ii] <- MyV[length(MyV)];
        if (!is.null(MyV) && !is.na(MyV[1]) && MyV[1] >= 0) {
          CAdd <- CAdd+1;
        }
      } else if (exists(\"MyV\") && is.null(MyV)) {
        AFilePrint(paste(\" Hey NULL MyV found for MyS = \", MyS, sep=\"\"));
      } else {
        AFilePrint(paste(\" Hey no MyV found for MyS = \", MyS, sep=\"\"));
      }
    }  
    }
   if (is.null(MySuccessArray)) {
        print(paste(\"Second NULL Success array, note TargetTotal = \", 
          TargetTotal, \" name functions are len \", length(NameFunctions),sep=\"\")); flush.console();
        MySuccessArray <- array(-1, dim=c(TargetTotal, length(NameFunctions),8),
          dimnames = list(Sim=(1:TargetTotal), Function=NameFunctions, Result=1:8)); 
   }
  ";
  return(ART);  
    
}

SaveProcessWorkToGlobal <- function() {
  ART <- "
  try(setwd(ASmallContainDir));
  try(dir.create(\"SumPerformance\", recursive=TRUE, showWarnings=FALSE));
  try(save(MySuccessArray=MySuccessArray, MyCompletesProcessTable=MyCompletesProcessTable,
    TargetTotal=TargetTotal, NameFunctions=NameFunctions, 
    MySuccessTimesTable=MySuccessTimesTable, MySuccesses=MySuccesses,
    file=paste(\"SumPerformance/SuccessArray.RData\", sep=\"\")));
  eval(parse(text=SetGText(\"MySuccessArray\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"MyCompletesProcessTable\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"MySuccessTimesTable\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"MySuccesses\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"MyKills\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"MyKillList\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"MyDies\", \"globalenv()\", S=1)));
  eval(parse(text=SetGText(\"MyDieList\", \"globalenv()\", S=1)));
  try(eval(parse(text=SetGText(\"KillList\", \"globalenv()\", S=1))));
  if (exists(\"KillList\") && length(KillList) >= 1) {
    for (ii in 1:length(KillList)) {
      AKK <- (1:length(NameFunctions))[NameFunctions==KillList[ii]];
      AFF <- (1:length(MySuccessArray[,1,1]))[MySuccessArray[,AKK,1] < 0];
      MySuccessArray[AFF,AKK,1] <- -666;
    }
  }
  eval(parse(text=GetG0Text(\"OriginalOldwd\")));
  try(setwd(OriginalOldwd));
  ";
  return(ART);  
    
}

ReadyForThralls <- function() {
  ART <- "
    ABB <- (1:(TargetTotal*length(NameFunctions)))[MySuccessArray[,,1] < 0
    & MySuccessArray[,,1]  != -666]; 
  SampleSort <- sample(ABB, size=length(ABB), replace=FALSE);
  LengthGo <- ceiling(length(SampleSort)/TotalThralls);
  dir.create(\"Submissions\", showWarnings=FALSE, recursive=TRUE);
  MySubs <- unlist(list.files(\"Submissions\"));
  if (length(MySubs) >= 1) {
  for (ii in 1:length(MySubs)) {
    try(unlink(paste(\"Submissions//\", MySubs, sep=\"\"), recursive=TRUE));
  }
  }

  OGo <- LengthGo;
  eval(parse(text=GetG0Text(\"argeese\", \"globalenv()\", S=1)));
  eval(parse(text=GetG0Text(\"MyNameOfInputFile\", \"globalenv()\", S=1)));
  if (is.null(MyNameOfInputFile) || is.numeric(MyNameOfInputFile)) {
    MyNameOfInputFile <- (unlist(strsplit(argeese, \"\\\\.\")))[1]
  }
  ";
  return(ART);  
    
}

SetupMySimulations <- function() {
  ATT <- "
  try(setwd(ASmallContainDir));
  MySimulations <- list.files(paste(\"SmallSMS\", sep=\"\"));
  MySimulations <- MySimulations[MySimulations != \"SMS.RData\"]
  MySimulations <- MySimulations[substr(MySimulations,1,nchar(\"SMS\")) == \"SMS\"];
  MySimNoSMS <-  substr(MySimulations, nchar(\"SMS\")+1, nchar(MySimulations));
  MySimulations <- MySimulations[MySimNoSMS != \"\"];
  MySimNoSMS <- MySimNoSMS[MySimNoSMS != \"\"];
  MySims <- MySimulations;
  AFFT <- substr(MySimNoSMS, nchar(MySimNoSMS)-nchar(\".RData\")+1,
    nchar(MySimNoSMS)) == \".RData\";
  if (length(AFFT) >= 1) {
    MySimNoSMS[AFFT] <- substr(MySimNoSMS[AFFT], 1, 
      nchar(MySimNoSMS[AFFT])-nchar(\".RData\"))
  }
  if (verbose >= 1) {
    if (length(MySimNoSMS) == 0) {
      AFilePrint(\"  -- We have now that MySimNoSMS has no data. \");
    } else {
      AFilePrint(paste(\"  -- We have now that MySimNoSMS has length \", 
        length(MySimNoSMS), \". \", sep=\"\"));        
    }
  }
  eval(parse(text=GetG0Text(\"OriginalOldwd\")));
  try(setwd(OriginalOldwd));
  ";
  return(ATT);   
}
CutMySimulations <- function() {
  ATT <- "
    MySimulations <- MySimulations[1:TargetTotal];
    MySims <- MySims[1:TargetTotal];  
    MySimNoSMS <- MySimNoSMS[1:TargetTotal];  
  "
  return(ATT);
}