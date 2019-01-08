MyFunctionRecord = setRefClass("MyFunctionRecord")

if (FALSE) {
ATestClass= setRefClass("ATestClass",
  fields = list(AValue="integer", BValue="numeric", CValue="character",
  DValue = function(AA=NULL) { if (is.null(AA)) { return(CValue); 
    } else {CValue <- AA; }}),
  methods=list(
    initialize=function(AValue=0,BValue=1, CValue="Oh Boy") {
    try(.self$AValue <- as.integer(AValue));
    try(.self$BValue <- as.real(BValue));
    try(.self$CValue <- as.character(CValue));
    return(.self);
  })
  );
  MyNew <- ATestClass$new(2,3.5, "acb");
  try(save(MyNew=MyNew, file=paste("c:/Stat/ASaveTest.RData", sep="")));
  rm(MyNew);
  load(paste("c:/Stat/ASaveTest.RData", sep=""))
  MyNew$AValue;  MyNew$DValue;
}
SummarySimulationList = setRefClass("SummarySimulationList",
 fields = list(
   TheSummarySimulationList = "list",
   NameFunctions = "character",
   .MyTotalBeatenFunctions = "integer", 
   .MyTotalFailFunctions = "integer",
   .MyTotalCompleteFunctions = "integer",
   .MyTotalKilledFunctions = "integer",
   .LastCheckTotalKilled = "integer",
   MyTotalKilledFunctions=function() {
     if (length(.self$.LastCheckTotalKilled) <= 0 || .self$.LastCheckTotalKilled[1] <= 0) {
        try(.self$.LastCheckTotalKilled <- as.integer(0));
     }
     if (is.null(.self$.MyTotalKilledFunctions)  || length(.self$.MyTotalKilledFunctions) <= 0 ||
       .self$.LastCheckTotalKilled[1] >= 40) {
       try(.self$.MyTotalKilledFunctions <- as.integer(.self$GetTotalKilledFunction()));
       try(names(.self$.MyTotalKilledFunctions) <- .self$NameFunctions);
       try(.self$.LastCheckTotalKilled[1] <- as.integer(0));
     }
     try(.self$.LastCheckTotalKilled <- as.integer(.self$.LastCheckTotalKilled+1));
     return(.self$.MyTotalKilledFunctions);
   },
   SampleFunctionProbabilities = "numeric", 
   LargeContainingDir = "character",  
   SmallContainingDir = "character",  
   CurrentLargeContainDir = function(XX = NULL) {
     if (is.null(XX)) {
        return(.self$LargeContainingDir);
     }
     .self$LargeContainingDir <- XX;
   },
   CurrentSmallContainDir = function(XX = NULL) {
     if (is.null(XX)) {
        return(.self$SmallContainingDir);
     }
     .self$SmallContainingDir <- XX;
   },
   MyBaseSimulationStorageDirectory = "character",
   CountCheckValidSimulations = "integer",
   CurrentTotalCompleted = function() {
     AR <- NULL;
     try(AR <- .self$GetTotalCompletedFunction());
     return(AR);
   },
  .ValidSimulations = "integer",
   ValidSimulations = function() {
     if (is.null(.self$CountCheckValidSimulations) ||
       is.na(.self$CountCheckValidSimulations) || length(.self$CountCheckValidSimulations) <= 0 || 
       .self$CountCheckValidSimulations <= 0) {
         try(.self$CountCheckValidSimulations <- as.integer(1))
     }
     if (is.null(.self$.ValidSimulations) || length(.self$.ValidSimulations) <= 0) {
         .self$.ValidSimulations <- 1:length(NameFunctions);
     }
     if (.self$CountCheckValidSimulations >= 20) {
       try(.self$CountCheckValidSimulations <- as.integer(1)); 
       .self$.ValidSimulations <- 1:length(NameFunctions);
       FailOuts <- .self$MyTotalKilledFunctions;
       eval(parse(text=GetG0Text("TooManyFailsNowKillFunction", "globalenv()", S=1)));
       if (is.null(TooManyFailsNowKillFunction) || TooManyFailsNowKillFunction <= 0) {
         TooManyFailsNowKillFunction <- 25;
       }
       .self$.ValidSimulations <- (1:length(FailOuts))[FailOuts <= TooManyFailsNowKillFunction]
     } else {
       try(.self$CountCheckValidSimulations <- as.integer(.self$CountCheckValidSimulations +1));
     }
      return(.self$.ValidSimulations);
    },
    GoodRecover = "logical"
   ),
   methods = list(
   initialize=function(MySmallBaseDirectory, MyLargeBaseDirectory) {
     ##AFilePrint("Start Creation of SummaryRecover");
     eval(parse(text=GetG0Text("CurrentSmallContainDir")));
     try(.self$SmallContainingDir <- as.character(CurrentSmallContainDir));
     eval(parse(text=GetG0Text("CurrentLargeContainDir")));
     try(.self$LargeContainingDir <- as.character(CurrentLargeContainDir));    
     ## AFilePrint("SummaryRecover: recover BaseStorage");
     try(.self$MyBaseSimulationStorageDirectory <- 
       as.character(paste(CurrentLargeContainDir, "//SumPerformance", sep="")));
     try(MyOld <- getwd());
     ##AFilePrint(paste("SummaryRecover: Try to Setwd to SmallContainDir: ", CurrentSmallContainDir, sep=""));
     try(setwd(CurrentSmallContainDir));
     try(dir.create("SumPerformance", showWarnings=FALSE, recursive=TRUE));
     try(setwd(MyOld));
     ##AFilePrint(paste("SummaryRecover: Try to set GoodRecover.", sep=""));
     try(.self$GoodRecover <- as.logical(TRUE));
     ##AFilePrint(paste("SummaryRecover: Try to get NameFunctions.", sep=""));
     eval(parse(text=GetG0Text("NameFunctions")));
     try(.self$NameFunctions <- NameFunctions);
     ##AFilePrint(paste("SummaryRecover: Try to set TheSummarySimulationList to a list.",  sep=""));
     try(.self$TheSummarySimulationList <- list());
     ##AFilePrint(paste("SummaryRecover: Try to set MyTotalBeatenFunctions.", sep=""));
     try(.self$.MyTotalBeatenFunctions <- as.integer(0))
     try(.self$.MyTotalFailFunctions <- as.integer(0))
     try(.self$.MyTotalCompleteFunctions <- as.integer(0))  
     try(.self$CountCheckValidSimulations <- as.integer(0)); 
     ##AFilePrint("SummaryRecover: Set SampleFunctionProbabilities "); flush.console();
     try(.self$SampleFunctionProbabilities <- rep(1, length(.self$NameFunctions)));   
     try(.self$.LastCheckTotalKilled  <- as.integer(0));
     ##AFilePrint("SummaryRecover: Set MyTotalKilledFunctions "); flush.console();
     try(.self$.MyTotalKilledFunctions <- vector("integer",0))
     ##AFilePrint("SummaryRecover: Return .self "); 
     return(.self);
   },
   GetTotalFailFunction=function() {
     if (is.null(.self$TheSummarySimulationList) || length(.self$TheSummarySimulationList) <= 0) {
        MyRep <- rep(0, length(.self$NameFunctions));
        names(MyRep) <- .self$NameFunctions;
        return(MyRep);
     }
     MyRep <- rep(0, length(.self$NameFunctions));
     names(MyRep) <- .self$NameFunctions;
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        if (ii <= length(.self$TheSummarySimulationList) &&
          !is.null(.self$TheSummarySimulationList[[ii]]) &&
          length(.self$TheSummarySimulationList[[ii]]$MyFails) ==
          length(MyRep)) {
          try(MyRep <- MyRep + .self$TheSummarySimulationList[[ii]]$MyFails);
        }
     }
     try(names(MyRep) <- .self$NameFunctions);
     ATryText <- "
     SampleFunctionProbabilities <- rep(1.0, length(.self$NameFunctions));
     try(names(SampleFunctionProbabilities) <- .self$NameFunctions)
     try(.self$SampleFunctionProbabilities <- SampleFunctionProbabilities);
     ";
     try(eval(parse(text=ATryText)));
     return(MyRep);
   },
   GetSimFailFunction=function() {
     if (is.null(.self$TheSummarySimulationList) || length(.self$TheSummarySimulationList) <= 0) {
        MyRep <- rep(0, length(.self$NameFunctions));
        names(MyRep) <- .self$NameFunctions;
        return(MyRep);
     }
     MyRep <- rep(0, length(.self$NameFunctions));
     names(MyRep) <- .self$NameFunctions;
     ABIt <- 0;
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        AText <- "
          if (ii <= length(.self$TheSummarySimulationList) &&
            !is.null(.self$TheSummarySimulationList[[ii]]) &&
            length(.self$TheSummarySimulationList[[ii]]$MyFails) ==
            length(MyRep)) {
            Aad <- .self$TheSummarySimulationList[[ii]]$MyFails;
            Aad[Aad >= 1] <- 1;
            try(MyRep <- MyRep + .self$TheSummarySimulationList[[ii]]$MyFails);
          }
          if (is.null(.self$TheSummarySimulationsList[[ii]])) {
            ABIt <- 1;
          }
        ";
        try(eval(parse(text=AText)));
     }
     try(names(MyRep) <- .self$NameFunctions);
     if (ABIt == 1) {
        try(.self$ElimNullSSL());
     }
     return(MyRep);
   },
   GetTotalKilledFunction=function() {
     if (is.null(.self$TheSummarySimulationList) || length(.self$TheSummarySimulationList) <= 0) {
        MyRep <- rep(0, length(.self$NameFunctions));
        names(MyRep) <- .self$NameFunctions;
        return(MyRep);
     }
     MyRep <- rep(0, length(.self$NameFunctions));
     names(MyRep) <- .self$NameFunctions;
     ABIt <- 0;
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        AText <- "
          Aad <- rep(0, length(.self$NameFunctions));
          if (ii <= length(.self$TheSummarySimulationList) &&
            !is.null(.self$TheSummarySimulationList[[ii]]) &&
            length(.self$TheSummarySimulationList[[ii]]$AreKilled) ==
            length(Aad)) {
            Aad <- .self$TheSummarySimulationList[[ii]]$AreKilled;
            Aad[Aad == TRUE] <- 1;
            try(MyRep <- MyRep + Aad);
          }
          if (is.null(.self$TheSummarySimulationList[[ii]])) {
            ABIt <- 1;
          }
        ";
        try(eval(parse(text=AText)));
        ##if (length(MyRep) <= 0) {
        ##  print(paste("Error On ii = ", ii, sep=""));
        ##  break;
        ##}
     }
     try(names(MyRep) <- .self$NameFunctions);
     if (ABIt == 1) {
        try(.self$ElimNullSSL());
     }
     return(MyRep);
   },
   GetTotalCompletedFunction=function() {
     if (is.null(.self$TheSummarySimulationList) || length(.self$TheSummarySimulationList) <= 0) {
        MyRep <- rep(0, length(.self$NameFunctions));
        names(MyRep) <- .self$NameFunctions;
        return(MyRep);
     }
     MyRep <- rep(0, length(.self$NameFunctions));
     names(MyRep) <- .self$NameFunctions;
     ABIt <- 0;
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        if (!is.null(.self$TheSummarySimulationList[[ii]]) &&
          length(.self$TheSummarySimulationList[[ii]]$AreCompleted) == 
          length(MyRep)) {
          try(MyRep <- MyRep + .self$TheSummarySimulationList[[ii]]$AreCompleted);     
        }   
        if (is.null(.self$TheSummarySimulationList[[ii]])) {
          ABIt <- 1;
        }
     }
     try(names(MyRep) <- .self$NameFunctions)
     if (ABIt == 1) {
        try(.self$ElimNullSSL());
     }
     return(MyRep);
   },
   GenTotalCompleted=function() {return(.self$GetTotalCompletedFunction()); },
   GetTotalCompleted=function() {return(.self$GetTotalCompletedFunction()); },
   GenCurrentTotalCompleted=function() {return(.self$GetTotalCompletedFunction()); },
   GetCurrentTotalCompleted=function() {return(.self$GetTotalCompletedFunction()); },
   GetTotalAttemptedFunction=function() {
     if (is.null(.self$TheSummarySimulationList) || length(.self$TheSummarySimulationList) <= 0) {
        MyRep <- rep(0, length(.self$NameFunctions));
        names(MyRep) <- .self$NameFunctions;
        return(MyRep);
     }
     MyRep <- rep(0, length(.self$NameFunctions));
     names(MyRep) <- .self$NameFunctions;
     ABIt <- 0;
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        if (!is.null(.self$TheSummarySimulationList[[ii]]) &&
          length(.self$TheSummarySimulationList[[ii]]$AttemptsOnAFunction) ==
          length(MyRep)) {
          try(MyRep <- MyRep + 
            .self$TheSummarySimulationList[[ii]]$AttemptsOnAFunction);
        }
        if (is.null(.self$TheSummarySimulationList)) {
          ABIt <- 1;
        }
     }
     try(names(MyRep) <- .self$NameFunctions);
     if (ABIt == 1) {
        try(.self$ElimNullSSL());
     }
     return(MyRep);
   },
   GetTotalBeatFunction=function() {
     if (is.null(.self$TheSummarySimulationList) || length(.self$TheSummarySimulationList) <= 0) {
        MyRep <- rep(0, length(.self$NameFunctions));
        names(MyRep) <- .self$NameFunctions;
        return(MyRep);
     }
     MyRep <- rep(0, length(.self$NameFunctions));
     names(MyRep) <- .self$NameFunctions;
     ABIt <- 0;
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        if (!is.null(.self$TheSummarySimulationList[[ii]]) && 
          length(.self$TheSummarySimulationList[[ii]]$MyBeatOuts) ==
          length(MyRep)) {
          try(MyRep <- MyRep + .self$TheSummarySimulationList[[ii]]$MyBeatOuts);
        }
        if (is.null(.self$TheSummarySimulationList[[ii]])) {
          ABIt <- 1;
        }
     }
     try(names(MyRep) <- .self$NameFunctions)
     if (ABIt == 1) {
        try(.self$ElimNullSSL());
     }
     return(MyRep);
   },
   GenPerformanceMatrix = function() {
     if (is.null(.self$TheSummarySimulationList) || length(.self$TheSummarySimulationList) <= 0) {
        return(NULL);
     }
     PerformanceMatrix <- matrix(0, length(.self$TheSummarySimulationList),
       length(.self$NameFunctions));
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        CC <- rep(0, length(.self$NameFunctions));
        CD <- rep(FALSE, length(.self$NameFunctions));
        if (!is.null(.self$TheSummarySimulationList[[ii]]) &&
          length(.self$TheSummarySimulationList[[ii]]$AreCompleted) == 
          length(CC)) {
          CC <- .self$TheSummarySimulationList[[ii]]$AreCompleted;    
        }
        if (!is.null(.self$TheSummarySimulationList[[ii]]) &&
          length(.self$TheSummarySimulationList[[ii]]$AreKilled) == 
          length(CD)) {
          CD <- .self$TheSummarySimulationList[[ii]]$AreKilled;
        }
        if (!is.null(.self$TheSummarySimulationList[[ii]]) &&
          length(.self$TheSummarySimulationList[[ii]]$AreAttempting) == 
          length(CD)) {
          PerformanceMatrix[ii,] <- .self$TheSummarySimulationList[[ii]]$AreAttempting;
        }
        PerformanceMatrix[ii,CC >= 1] <- - 5;
        PerformanceMatrix[ii,CD == TRUE] <- -666;
     }
     return(PerformanceMatrix);
   },
   ElimNullSSL = function() {
     ACC <- 0;
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        if (is.null(.self$TheSummarySimulationList[[ii]])) {
          ACC <- ACC+1;
        }
     }
     print(paste("ElimNullSSL: there were ", ACC, " NULL Simulations! ", sep=""));
     NewTSL <- list();
     OnIn <- 0;
     for (ii in 1:length(.self$TheSummarySimulationList)) {
        if (!is.null(.self$TheSummarySimulationList[[ii]])) {
          OnIn <- OnIn +1;
          NewTSL[[OnIn]] <- .self$TheSummarySimulationList[[ii]]  
        }
     }
     try(.self$TheSummarySimulationList <- NewTSL);
   }
   )
);

SummarySimulationRecord = setRefClass("SummarySimulationRecord",
  fields = list(UniqueSimulationIdentifier="character",
  SavedUniqueSimulationRDataFile="character",
  TriedSimulations="vector",
  TotalAttempted="integer",
  TotalCompleted="integer",
  TotalWorking="integer",
  MyAttempts="vector",
  FirstSuccessfulAttempt = "ANY",
  FailedAttempts = "vector",
  AttemptsOnAFunction="integer",
  MyFails="vector",
  MyBeatOuts="integer",
  AreCompleted="logical",
  Completed="logical",
  AreKilled="logical",
  WinningUniqueProcess="character",

  AreAttempting = function() {
     if (is.null(.self$MyAttempts)) {
        return(NULL);
     }
     AttemptVec <- rep(0, length(.self$MyAttempts));
     for (ii in 1:length(.self$MyAttempts)) {
        if (.self$AreCompleted[ii] == TRUE) {
          AttemptVec[ii] <- -1;
        } else if (.self$AreKilled[ii] == TRUE) {
          AttemptVec[ii] <- -666;
        } else { 
          AttemptVec[ii] <- 0;
          if (!is.null(.self$MyAttempts[[ii]]) && length(.self$MyAttempts[[ii]]) >= 1  && 
            !(is.numeric(.self$MyAttempts[[ii]]) &&
              (.self$MyAttempts[[ii]] == 1 || .self$MyAttempts[[ii]] <= 0))) {
              for (jj in 1:length(.self$MyAttempts[[ii]])) {
               if (!(is.numeric(.self$MyAttempts[[ii]][[jj]]) && 
                 .self$MyAttempts[[ii]][[jj]][1] < 0) &&
                 .self$MyAttempts[[ii]][[jj]]$IFailed == TRUE) {
                
               } else {
                AttemptVec[ii] <- AttemptVec[ii]+1;
              }
            }
          }
        }
    }
    return(AttemptVec);
  }
  ),
  methods=list(
    initialize=function(UniqueSimulationIdentifier = "DefaultSimulation",
     SavedUniqueSimulationRDataFile = "SavedUniqueSimulationRDataFile") { 
     .self$UniqueSimulationIdentifier <- as.character(UniqueSimulationIdentifier);
     .self$SavedUniqueSimulationRDataFile <- as.character(SavedUniqueSimulationRDataFile);      
     eval(parse(text=GetG0Text("NameFunctions")));   
     AS <- vector("list", length(NameFunctions));  
     try(names(AS) <- NameFunctions)
     try(.self$MyAttempts  <- AS);
     try(.self$FailedAttempts <- AS);
     try(.self$TotalAttempted <- as.integer(0));
     try(.self$TotalCompleted <- as.integer(0));
     try(.self$TotalWorking <- as.integer(0)); 
     try(.self$Completed <- as.logical(FALSE));
     try(.self$FirstSuccessfulAttempt <- NULL);
     try(.self$WinningUniqueProcess <- rep("", length(NameFunctions))); 
     try(.self$AreCompleted <- rep(FALSE, length(NameFunctions)));
     try(.self$AreKilled <- rep(FALSE, length(NameFunctions)));
     ATryFailText <- "
     MyFails <- rep(0, length(.self$MyAttempts));
     names(MyFails) <- NameFunctions;
     try(.self$MyFails <- MyFails); 
     MyBeatOuts <- rep(0, length(.self$MyAttempts));
     try(names(MyBeatOuts)   <- NameFunctions);
     try(.self$MyBeatOuts <- as.integer(MyBeatOuts));
     try(names(.self$MyBeatOuts) <- NameFunctions);
     ";
     try(eval(parse(text=ATryFailText)));
     AttemptsText <- "
     AttemptsOnAFunction = rep(0, length(.self$MyAttempts));
     try(names(AttemptsOnAFunction) <- NameFunctions)
     try(.self$AttemptsOnAFunction <- as.integer(AttemptsOnAFunction));
     try(names(.self$AttemptsOnAFunction) <- NameFunctions);
     ";
     try(eval(parse(text=AttemptsText)))
     return(.self);     
    },
    AddAFunction = function(ANameFunction = "DefaultANameFunction", ForcingFunction = FALSE,
      SSL = NULL, UniqueProcessIdentifier=NULL, NoWrite=FALSE) {
      if (is.null(SSL)) {
        AFileError("AddAFunction: Error you must non NULL SSL!");
        return("AFail");
      }
      if (is.integer(ANameFunction) || is.numeric(ANameFunction)) {
        AOn <- ANameFunction;
        ANameFunction <- names(.self$MyAttempts)[AOn];
      } else {
        if (!(ANameFunction %in% names(.self$MyAttempts))) {
          AFilePrint(paste("AddAFunction:  No function named ", ANameFunction,
            " in a file. ", sep="")); flush.console();
          AErrorPrint(paste("AddAFunction:  No function named ", ANameFunction,
            " in a file. ", sep="")); flush.console();
          tryCatch("AddAFunction though ANameFunction");
        }
        AOn <-  (1:length(.self$MyAttempts))[
          names(.self$MyAttempts) == ANameFunction];
      }
      ANameFunction <- paste(unlist(strsplit(ANameFunction, " ")), collapse="");
      if (.self$AreCompleted[AOn] == TRUE) {
        AFilePrint(paste("AddAFunction:\n       for OnSimulation: ", 
          .self$UniqueSimulationIdentifier, "\n    function ", ANameFunction,
          "\n was already long completed!  ", sep="")); flush.console();
        AErrorPrint(paste("AddAFunction:\n       for OnSimulation: ", 
          .self$UniqueSimulationIdentifier, "\n    function ", ANameFunction,
          "\n was already long completed!  ", sep="")); flush.console();
        return("AFail");
        tryCatch("AddAFunction Fails to add this function!");
      }
      if (NoWrite == FALSE && ForcingFunction == FALSE  && .self$AttemptsOnAFunction[AOn] - .self$MyFails[AOn] >= 1) {
        AFilePrint("*******************************************************************************");
        AFilePrint("*** AddAFunction: Wait a minute: we are in testing phase and this shouldn't happen");
        AFilePrint(paste("***  You are seeking to add ", .self$UniqueSimulationIdentifier, sep=""));
        AFilePrint(paste("***  To run function ", ANameFunction, " with AOn = ", AOn, " again ", sep=""));
        AFilePrint("***  Not Sure what to do at this point!  ");
        AFilePrint(paste("*** AttemptsOnAFunction[AOn=", AOn, "] = ", .self$AttemptsOnAFunction[AOn], sep=""));
        AFilePrint(paste("*** MyFails[AOn=", AOn, "] = ", .self$MyFails[AOn], sep=""));
        eval(parse(text=SetGText("AOn", S=1)))
        try(AFilePrint(paste("**  AttemptsOnAFunction are: (", 
          paste(.self$AttemptsOnAFunction, collapse=", "), ")", sep="")));
        AFilePrint("**************************************************************");
       ## return("AFail");
       ## tryCatch("*** AddAFunction Fails for this duplicate Call!");
      }
      try(.self$AttemptsOnAFunction[AOn] <- as.integer(.self$AttemptsOnAFunction[AOn]+1));
      if (is.null(UniqueProcessIdentifier)) {
        eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()")));
      }
      if (is.null(.self$MyAttempts[[AOn]]) || length(.self$MyAttempts[[AOn]]) <= 0) {
        .self$MyAttempts[[AOn]] <- list();
         .self$MyAttempts[[AOn]][[1]] <- MySummarySimulationAttempt$new(
          UniqueProcessIdentifier=UniqueProcessIdentifier, NameFunction = ANameFunction)
      } else {
        .self$MyAttempts[[AOn]][[length(.self$MyAttempts[[AOn]])+1]] <- 
          MySummarySimulationAttempt$new(
          UniqueProcessIdentifier=UniqueProcessIdentifier, NameFunction = ANameFunction);               
      }
      try(.self$TotalWorking <- as.integer(.self$TotalWorking+1));
      try(.self$TotalAttempted <- as.integer(.self$TotalAttempted +1));
      
    Oldwd <- getwd();
    if (NoWrite==FALSE) {
      MyOnDir <- paste(SSL$SmallContainingDir, sep="");
      try(setwd(MyOnDir));
      try(dir.create(paste(ANameFunction, "//", "MyAttempts", sep=""), showWarnings=FALSE,
        recursive=TRUE));
      while(TwoSimR5:::LockMeIn(verbose=0, 
        quoteMore = paste("LockIn New Function", ANameFunction, sep=""), 
        LFDir = ANameFunction, 
        NoSpaceLFDir = TRUE,
        MaxLockTime = DefaultMaxLockTime, 
        MoreLock = "", UniqueProcessIdentifier = UniqueProcessIdentifier) == FALSE) {
        try(Sys.sleep(4));
      }
      RDTry <- NULL;
      if (any(unlist(list.files(paste(ANameFunction, "//MyAttempts", sep="")))
         == paste(UniqueProcessIdentifier, ".txt", sep=""))) {
        try(RDTry <- secure.read.table(paste(ANameFunction, 
          "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), nCols=7, header=TRUE));
      }
      RealTime <- GetRealTime();
      if (is.null(RDTry) || NROW(RDTry) <= 0) {
        try(NewTable <- matrix(0, 1, 7));
        rownames(NewTable) <- .self$UniqueSimulationIdentifier;
        NewTable[1,1:2] <- c(RealTime$JustSeconds, RealTime$JustDays)
        try(NewTable[1,3] <- AOn);
        NewTable[1,4] <- .self$AttemptsOnAFunction[AOn];
        NewTable[1,5] <- as.integer(.self$AreCompleted[AOn]);
        try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, 
          "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
          sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", "AttemptsOnAFunction", 
          "AreCompleted", "EndSeconds", "EndDays"),
          row.names=.self$UniqueSimulationIdentifier));
      } else {
        if (any(rownames(RDTry) == .self$UniqueSimulationIdentifier &
          RDTry[,3] == AOn)) {
          try(NewTable <- RDTry)
          try(ATT <- (1:length(rownames(RDTry)))[rownames(RDTry) == .self$UniqueSimulationIdentifier &
            RDTry[,3] == AOn]);
          try(RDTry[ATT,4] <-RDTry[ATT,4]+1);
          try(RDTry[ATT,1:2]  <- c(RealTime$JustSeconds, RealTime$JustDays)); 
          try(TwoSimR5:::secure.write.table(RDTry, file=paste(ANameFunction, "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
            sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", "AttemptsOnAFunction", "AreCompleted",
            "EndSeconds", "EndDays"),
            row.names=rownames(RDTry) ));
        } else {
          try(NewTable <- matrix(0, NROW(RDTry)+1,7));
          if (NROW(RDTry) == 1) {
            try(NewTable[1,1:7] <- as.numeric(unlist(RDTry[1, 1:7])));
          } else {
            try(NewTable[1:NROW(RDTry),1:7] <- as.numeric(unlist(RDTry[1:NROW(RDTry), 1:7])));            
          }
          try(NewTable[NROW(NewTable),1:2] <- c(RealTime$JustSeconds, RealTime$JustDays));
          try(NewTable[NROW(NewTable),3] <- AOn);
          try(NewTable[NROW(NewTable),4] <- .self$AttemptsOnAFunction[AOn]);
          try(NewTable[NROW(NewTable),5] <- as.integer(.self$AreCompleted[AOn]));
          try(rownames(NewTable) <- c(rownames(RDTry), .self$UniqueSimulationIdentifier)); 
          try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
            sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", "AttemptsOnAFunction", "AreCompleted",
            "EndSeconds", "EndDays"),
            row.names=rownames(NewTable) ));
        } 
      }
        
        try(TwoSimR5:::UnLockMe(LFDir = ANameFunction, 
          verbose=verbose, quoteMore="Final Saw no problem Sizeup",
          NoSpaceLFDir=TRUE));
        try(setwd(Oldwd));
      }
      return(AOn);
    },
    CompleteAFunction = function(ANameFunction="DefaultNameFunction", SSL = NULL,
      UniqueProcessIdentifier=NULL, NoWrite=FALSE)  {
      eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
      if (is.null(SSL)) {
         AErrorPrint("CompleteAFunction: must supply SSL now")
      }
      if (is.integer(ANameFunction) || is.numeric(ANameFunction)) {
        AOn <- round(ANameFunction);
        ANameFunction <- names(.self$MyAttempts)[AOn];
      } else {
        if (!(ANameFunction %in% names(.self$MyAttempts))) {
          AFilePrint(paste("CompleteAFunction:  No function named ", ANameFunction,
            " in a file. ", sep="")); flush.console();
          AErrorPrint(paste("CompleteAFunction:  No function named ", ANameFunction,
            " in a file. ", sep="")); flush.console();
          tryCatch("AddAFunction though ANameFunction");
        }
        try(.self$TotalCompleted <- as.integer(.self$TotalCompleted+1));
        AOn <-  (1:length(.self$MyAttempts))[
        names(.self$MyAttempts) == ANameFunction];
      }
      if (is.null(UniqueProcessIdentifier)) {
        try(eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1))));
      }
      Oldwd <- getwd();
      if (.self$AreCompleted[AOn] == TRUE) {
        
        if (NoWrite==FALSE) {
        MyOnDir <-  paste(SSL$SmallContainingDir, sep="");
        try(setwd(MyOnDir)); 
        try(dir.create(paste(ANameFunction, "//MyAttempts", sep=""), showWarnings=FALSE, recursive=TRUE));
        while(TwoSimR5:::LockMeIn(verbose=0, quoteMore = "On Recover", LFDir = ANameFunction,
          MaxLockTime = DefaultMaxLockTime, MoreLock = "", UniqueProcessIdentifier = UniqueProcessIdentifier,
          NoSpaceLFDir=TRUE) == FALSE) {
          for (ii in 1:100) {}
        } 
        RDTry <- NULL;
        if (any(unlist(list.files(paste(ANameFunction, "//MyAttempts", sep="")))
           == paste(UniqueProcessIdentifier, ".txt", sep=""))) {
          try(RDTry <- secure.read.table(paste(ANameFunction, 
            "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), nCols = 7, header=TRUE));
        }
        RealTime <- GetRealTime();
        if (is.null(RDTry) || NROW(RDTry) <= 0) {
          NewTable <- matrix(0, 1, 7);
          rownames(NewTable) <- .self$UniqueSimulationIdentifier;
          NewTable[1,1:2] <- c(RealTime$JustSeconds, RealTime$JustDays);
          try(NewTable[1,3] <- AOn);
          NewTable[1,4] <- .self$AttemptsOnAFunction[AOn];
          NewTable[1,5] <- as.integer(666);
          try(NewTable[6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));
          try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, "//MyAttempts//", 
            UniqueProcessIdentifier, ".txt", sep=""), 
            sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", "AttemptsOnAFunction", "AreCompleted",
              "EndSeconds", "EndDays"),
            row.names=.self$UniqueSimulationIdentifier));
        } else {
          if (any(rownames(RDTry) == .self$UniqueSimulationIdentifier &
            RDTry[,3] == AOn)) {
              try(ATT <- (1:length(rownames(RDTry)))[rownames(RDTry) == .self$UniqueSimulationIdentifier &
                RDTry[,3] == AOn]);
              try(NewTable <- RDTry);
              if (NewTable[ATT,5] <= 0) {
                try(NewTable[ATT,5] <- as.integer(666));
              }  else {
                try(MyFFDir <- paste(ANameFunction, sep=""));
                try(LTT <- unlist(list.files(MyFFDir)));              
                if ("OutTable.RData" %in% LTT) {
                  try(load(paste(ANameFunction, "//OutTable.RData", sep="")));
                  RowOfSim=0;
                  try(RowOfSim <- (1:NROW(RDD))[rownames(RDD) == .self$UniqueSimulationIdentifier]);
                  try(RDDRow <- RDD[RowOfSim,]);
                } else {
                  RowOfSim = 0;  RDD = NULL;
                  Line = 0; RDDRow = NULL;
                }
                MyEP <- paste("***************************************************
                  ****************************************************************
                  Hey, we wanted to save complete to NewTable Function AOn = ", AOn, "
                  but we found that it was already solved with ATT = ", ATT, "
                  and NewTable[ATT,] <- ", paste(NewTable[ATT,], collapse=", "),"
                  UniqueProcessIdentifier = ", UniqueProcessIdentifier, "
                  UniqueSimulationIdentifier = ", .self$UniqueSimulationIdentifier, "
                  list.files(HOMEDIR) = (", paste("\"", unlist(list.files()), "\"", collapse=",
                      ", sep=""), ")
                  And RDDRow = [", paste(names(RDDRow), collapse=", "), "]
                               (", paste(RDDRow, collapse=", "), ")
                  What is going on?
                  *****************************************************************
                \n", sep="");
                AErrorPrint(MyEP);
                AFilePrint(MyEP);
                return(-1);
              } 
              try(NewTable[ATT,6:7] <- c(RealTime$JustSeconds,RealTime$JustDays));
          } else {
            try(NewTable <- matrix(0, NROW(RDTry)+1,7));
            if (NROW(RDTry) == 1) {
              try(NewTable[1,1:7] <- as.numeric(unlist(RDTry[1, 1:7])));
            } else {
              try(NewTable[1:NROW(RDTry),1:7] <- as.numeric(unlist(RDTry[1:NROW(RDTry), 1:7])));            
            }
            try(NewTable[NROW(NewTable),1:2] <- c(RealTime$JustSeconds, RealTime$JustDays));  
            try(NewTable[NROW(NewTable),6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));  
            try(NewTable[NROW(NewTable),3] <- AOn);
            try(NewTable[NROW(NewTable),4] <- .self$AttemptsOnAFunction[AOn]);
            try(NewTable[NROW(NewTable),5] <- as.integer(666));
            try(rownames(NewTable) <- c(rownames(RDTry), .self$UniqueSimulationIdentifier)); 
          }

           try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, 
            "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
            sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", "AttemptsOnAFunction", 
              "AreCompleted", "EndSeconds", "EndDays"),
            row.names=rownames(NewTable) )); 
         }
         try(TwoSimR5:::UnLockMe(LFDir = paste(ANameFunction, sep=""), 
            verbose=verbose, quoteMore="Finished Complete On"));
         try(setwd(Oldwd));
         }
         return(10);
      }
      if ((is.integer(.self$MyAttempts[[AOn]]) || is.numeric(.self$MyAttempts[[AOn]])) &&
        .self$MyAttempts[[AOn]][1] == 1) {
        AFilePrint(paste("CompleteAFunction:  Hey AOn = ", AOn, 
          " and we have integer object but AreCompleted isn't true!",
          sep=""));
        return(11);
      }
      try(.self$AreCompleted[AOn] <- TRUE);
      if (all(.self$AreCompleted == TRUE)) {
        try(.self$Completed <- TRUE);
      }
      if (is.null(UniqueProcessIdentifier)) {
        eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
      }
      try(.self$MyBeatOuts[AOn] <- as.integer(0));
      if (length(.self$MyAttempts[[AOn]]) >= 1) {
      for (ii in 1:length(.self$MyAttempts[[AOn]])) {
        ATT <- "
        if (.self$MyAttempts[[AOn]][[ii]]$UniqueProcessIdentifier == UniqueProcessIdentifier) {
           .self$WinningUniqueProcess[AOn] <- UniqueProcessIdentifier; 
           .self$MyBeatOuts[AOn] <- as.integer(ii-1);
           try(.self$FirstSuccessfulAttempt <- .self$MyAttempts[[AOn]][[ii]]);
           .self$MyAttempts[[AOn]] <- as.integer(1);
           break;
        } 
        ";
        try(eval(parse(text=ATT)));
      }
      }
      try(.self$TotalWorking <- as.integer(.self$TotalWorking-1));
      try(.self$TotalCompleted <- as.integer(.self$TotalCompleted+1));

      MyOnDir <-  paste(SSL$SmallContainingDir, sep="");
      Oldwd <- getwd(); 
      if (NoWrite==FALSE) {
      try(setwd(MyOnDir));
      while(TwoSimR5:::LockMeIn(verbose=0, quoteMore = "CompleteAFunction", LFDir = ANameFunction,
        MaxLockTime = DefaultMaxLockTime, MoreLock = "", UniqueProcessIdentifier = UniqueProcessIdentifier,
        NoSpaceLFDir = TRUE) == FALSE) {
         for (ii in 1:100) {}
      } 
      dir.create(paste(ANameFunction,"//MyAttempts",sep=""),
         showWarnings=FALSE, recursive=TRUE);
      RDTry <- NULL;
      if (any(unlist(list.files(paste(ANameFunction, "//MyAttempts", sep=""))) 
        == paste(UniqueProcessIdentifier, ".txt", sep=""))) {
        try(RDTry <- secure.read.table(paste(ANameFunction, "//MyAttempts//", 
          UniqueProcessIdentifier, ".txt", sep=""), nCols=7, header=TRUE));
      }
      RealTime <- GetRealTime();
      if (is.null(RDTry) || NROW(RDTry) <= 0) {
        NewTable <- matrix(0, 1, 7);
        rownames(NewTable) <- .self$UniqueSimulationIdentifier;
        NewTable[1,1:2] <- c(RealTime$JustSeconds, RealTime$JustDays);
        try(NewTable[1,3] <- AOn);
        NewTable[1,4] <- .self$AttemptsOnAFunction[AOn];
        NewTable[1,5] <- as.integer(1);
        try(NewTable[6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));
        try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, 
          "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
          sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", "AttemptsOnAFunction", "AreCompleted",
            "EndSeconds", "EndDays"),
          row.names=.self$UniqueSimulationIdentifier));
      } else {
        if (any(rownames(RDTry) == .self$UniqueSimulationIdentifier &
          RDTry[,3] == AOn)) {
            try(ATT <- (1:length(rownames(RDTry)))[rownames(RDTry) == .self$UniqueSimulationIdentifier &
              RDTry[,3] == AOn]);
            try(NewTable <- RDTry);
            try(NewTable[ATT,5] <- as.integer(1));
            try(NewTable[ATT,6:7] <- c(RealTime$JustSeconds,RealTime$JustDays));
        } else {
            try(NewTable <- matrix(0, NROW(RDTry)+1,7));
            if (NROW(RDTry) == 1) {
              try(NewTable[1,1:7] <- as.numeric(unlist(RDTry[1, 1:7])));
            } else {
              try(NewTable[1:NROW(RDTry),1:7] <- as.numeric(unlist(RDTry[1:NROW(RDTry), 1:7])));            
            }
          try(NewTable[NROW(NewTable),1:2] <- c(RealTime$JustSeconds, RealTime$JustDays));  
          try(NewTable[NROW(NewTable),6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));  
          try(NewTable[NROW(NewTable),3] <- AOn);
          try(NewTable[NROW(NewTable),4] <- .self$AttemptsOnAFunction[AOn]);
          try(NewTable[NROW(NewTable),5] <- as.integer(1));
          try(rownames(NewTable) <- c(rownames(RDTry), .self$UniqueSimulationIdentifier)); 
        }
        try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, "//MyAttempts//",
          UniqueProcessIdentifier, ".txt", sep=""), 
          sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", "AttemptsOnAFunction", 
            "AreCompleted", "EndSeconds", "EndDays"),
          row.names=rownames(NewTable) )); 
      }
      try(TwoSimR5:::UnLockMe(LFDir = ANameFunction, 
          verbose=verbose, quoteMore="Done Complete Function", NoSpaceLFDir=TRUE));
      try(setwd(Oldwd));
      }
      return(1);
    },
    KillAFunction= function(ANameFunction="DefaultNameFunction", SSL = NULL,
      UniqueProcessIdentifier=NULL, NoWrite=FALSE, ...)  {
      if (is.null(SSL)) {
         AErrorPrint("KillAFunction: Must supply SSL!");
      }
      if (!(ANameFunction %in% names(.self$MyAttempts))) {
        AFilePrint(paste("KillAFunction:  No function named ", ANameFunction,
          " in a file. ", sep="")); flush.console();
        AErrorPrint(paste("KillAFunction:  No function named ", ANameFunction,
          " in a file. ", sep="")); flush.console();
        tryCatch("KillAFunction though ANameFunction");
      }
      if (is.numeric(ANameFunction) && ANameFunction[1] > 0 && 
        ANameFunction[1] <= length(.self$MyAttempts)) {
        AOn <- round(ANameFunction[1]);
        ANameFunction <- names(.self$MyAttempts)[AOn];  
      } else {
        AOn <-  (1:length(.self$MyAttempts))[
          names(.self$MyAttempts) == ANameFunction];
      }
      if (.self$AreCompleted[AOn] == TRUE) {
        AFilePrint("*************************************************************");
        AFilePrint("**  Weird Weird Weird");
        AFilePrint(paste("**  We were told to kill ", ANameFunction, " but it actually completed! ", sep=""));
        AFilePrint("** We won't kill this one! ");
        return(-1);
      } else  {
        try(.self$AreKilled[AOn] <- TRUE);
      } 
      return(1);
    },
    FailAFunction = function(ANameFunction="DefaultNameFunction", SSL = NULL,
      UniqueProcessIdentifier=NULL, NoWrite=FALSE, ...)  {
      eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
      if (is.null(SSL)) {
         AErrorPrint("FailAFunction: Must supply SSL!");
      }
      if (verbose >= 3) {
        AFilePrint(paste("FailAFunction: start for ANameFunction = ",
          ANameFunction, sep=""));
      }
      if (!(ANameFunction %in% names(.self$MyAttempts))) {
        AFilePrint(paste("FailAFunction:  No function named ", ANameFunction,
          " in a file. ", sep="")); flush.console();
        AErrorPrint(paste("FailAFunction:  No function named ", ANameFunction,
          " in a file. ", sep="")); flush.console();
        tryCatch("FailAFunction though ANameFunction");
      }
      if (is.numeric(ANameFunction) && ANameFunction[1] > 0 && 
        ANameFunction[1] <= length(.self$MyAttempts)) {
        AOn <- round(ANameFunction[1]);
        ANameFunction <- names(.self$MyAttempts)[AOn];  
      } else {
        AOn <-  (1:length(.self$MyAttempts))[
          names(.self$MyAttempts) == ANameFunction];
      }
      if (is.null(UniqueProcessIdentifier)) {
        try(eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1))));
      }
      if (verbose >= 1) {
         AFilePrint(paste("FailAFunction: failing function ", ANameFunction, sep=""));
      }
      IDK <- 0;
      if (.self$AreCompleted[AOn] == TRUE) {
            IDK <- -333;
      } else if (.self$AreKilled[AOn] == TRUE) {
            IDK <- -666;
      } 
      OldOldwd <- getwd();
      if (.self$AreCompleted[AOn] == TRUE || .self$AreKilled[AOn] == TRUE) {  
        if (NoWrite==FALSE) {
        MyOnDir <-  paste(SSL$SmallContainingDir, sep="");
        if (verbose >= 2) {
            AFilePrint("FailAFunction: We are about to setwd to small contain dir ");
        } 
        try(setwd(SSL$SmallContainDir));
        if (verbose >= 2) {
            AFilePrint("FailAFunction: we just setwd to small contain dir, now to run LockMeIn ");
        } 
        while(TwoSimR5:::LockMeIn(verbose=0, quoteMore = "FailFunction", LFDir = 
           ANameFunction,
           MaxLockTime = DefaultMaxLockTime, MoreLock = "", 
           UniqueProcessIdentifier = UniqueProcessIdentifier,
           NoSpaceLFDir=TRUE) == FALSE) {
          for (ii in 1:100) {}
        } 
        try(dir.create(paste(ANameFunction, "//MyAttempts", sep=""), 
          showWarnings=FALSE, recursive=TRUE));
        RDTry <- NULL;
        if (any(unlist(list.files(paste(ANameFunction, "//MyAttempts", sep=""))) 
          == paste(UniqueProcessIdentifier, ".txt", sep=""))) {
          try(RDTry <- secure.read.table(paste(ANameFunction, 
            "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), nCols=7, header=TRUE));
        }
        RealTime <- GetRealTime();
        if (is.null(RDTry) || NROW(RDTry) <= 0) {
          try(NewTable <- matrix(0, 1, 7));
          try(rownames(NewTable) <- .self$UniqueSimulationIdentifier);
          try(NewTable[1,1:2] <- c(RealTime$JustSeconds, RealTime$JustDays))
          try(NewTable[1,3] <- AOn);
          try(NewTable[1,4] <- .self$AttemptsOnAFunction[AOn]);
          try(NewTable[1,5] <- as.integer(IDK));
          try(NewTable[1,6:7] <- c(RealTime$JustSeconds, RealTime$JustDays))
          try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, 
            "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
            sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", "AttemptsOnAFunction", "AreCompleted",
            "EndSeconds", "EndDays"),
            row.names=.self$UniqueSimulationIdentifier));
         } else {
          if (NROW(RDTry) >= 1 && any(rownames(RDTry) == .self$UniqueSimulationIdentifier &
            RDTry[,3] == AOn)) {
            try(ATT <- (1:length(rownames(RDTry)))[rownames(RDTry) == .self$UniqueSimulationIdentifier &
              AOn == RDTry[,3]]);
            try(NewTable <- RDTry);
            if (NewTable[ATT,5] <= 0) {
              try(NewTable[ATT,5] <- as.integer(IDK));
            } else {
              MyEP <- paste("***************************************************
                Hey, we wanted to save fail (though Complete) to NewTable Function AOn = ", AOn, "
                but we found that it was already solved with ATT = ", ATT, "
                and NewTable[ATT,] <- ", paste(NewTable[ATT,], collapse=", "),"
                What is going on?
              ", sep="");
              try(setwd(OldOldwd));
              AErrorPrint(MyEP);
              AFilePrint(MyEP);
              return(-1);
            }
            try(NewTable[ATT,6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));
          } else {
            try(NewTable <- matrix(0, NROW(RDTry)+1,7));
            if (NROW(RDTry) == 1) {
              try(NewTable[1,1:7] <- as.numeric(unlist(RDTry[1, 1:7])));
            } else {
              try(NewTable[1:NROW(RDTry),1:7] <- as.numeric(unlist(RDTry[1:NROW(RDTry), 1:7])));            
            }
            try(NewTable[NROW(NewTable),1:2] <- c(RealTime$JustSeconds, RealTime$JustDays));
            try(NewTable[NROW(NewTable),6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));
            try(NewTable[NROW(NewTable),3] <- AOn);
            try(NewTable[NROW(NewTable),4] <- .self$AttemptsOnAFunction[AOn]);
            try(NewTable[NROW(NewTable),5] <- as.integer(IDK));
            try(rownames(NewTable) <- c(rownames(RDTry), .self$UniqueSimulationIdentifier)); 
          }
          try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, 
            "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
            sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", 
            "AttemptsOnAFunction", "AreCompleted", "EndSeconds", "EndDays"),
            row.names=rownames(NewTable) )); 
         }
         try(TwoSimR5:::UnLockMe(LFDir =  ANameFunction, 
            verbose=verbose, quoteMore="Final Saw no problem Sizeup", NoSpaceLFDir=TRUE));
         if (verbose >= 2) {
            AFilePrint("FailAFunction: about to setwd to OldOldwd.");
         } 
         try(setwd(OldOldwd));
         }
         return(10);
      }
      if (is.null(UniqueProcessIdentifier)) {
        eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
      }
      try(.self$MyFails[AOn] <- as.integer(.self$MyFails[AOn]+1));
      MyMoveOut <- 0;
      if (length(.self$MyAttempts[[AOn]]) >= 1) {
      for (ii in 1:length(.self$MyAttempts[[AOn]])) {
        ATT <- "
        if (.self$MyAttempts[[AOn]][[ii]]$UniqueProcessIdentifier == UniqueProcessIdentifier) {
           .self$MyAttempts[[AOn]][[ii]]$IFailed<-TRUE;
           MyMoveOut <- ii;
           break;
        } 
        ";
        try(eval(parse(text=ATT)));
      }
      }
      if (TooManyFailsPerOneSimFunction > 0 && .self$MyFails[AOn] >= TooManyFailsPerOneSimFunction) {
          AFilePrint("****************************************************************************");
          AFilePrint(paste("**  Too many fails for function:  ", ANameFunction, sep=""));
          AFilePrint(paste("**  We will only tolerate ", TooManyFailsPerOneSimFunction, sep=""));
          try(AFilePrint(paste("**  But this is ", .self$MyFails[AOn], sep="")));
          AFilePrint(paste("***  Critical too many fails on AOn = ", AOn, " named ", ANameFunction, sep=""));
          AFilePrint("*** Time to kill this Method ");  flush.console();
          try(.self$AreKilled[AOn] <- as.logical(TRUE))
          try(.self$AreKilled[AOn] <- as.logical(TRUE));
          try(AFilePrint(paste("**  After our set we have AreKilled[", AOn, "] ",
            " is ", .self$AreKilled[AOn], sep="")));
        
          if (!is.null(SSL)) {
            ATryText <- "
             if (length(SSL$.MyTotalKilledFunctions) >= 1) {
               try(MyD <- SSL$.MyTotalKilledFunctions);
               try(MyD[AOn] <- as.integer(MyD[AOn]+1));
               try(SSL$.MyTotalKilledFunctions <- as.integer(MyD));
               try(names(SSL$.MyTotalKilledFunctions) <- SSL$NameFunctions);
             }
            ";
            try(eval(parse(text=ATryText)));
          }
      }
      if (MyMoveOut <= 0) {
        AFilePrint("**********************************************************");
        AFilePrint(paste("*** Hey, FailAFunction could not find process identifier: ", 
          UniqueProcessIdentifier, sep=""));
        try(AErrorPrint(paste("*** Hey, FailAFunction could not find process identifier: ", 
          UniqueProcessIdentifier, sep="")));
      } else {
      if (length(.self$MyAttempts[[AOn]]) == 1) {
        if (is.null(.self$FailedAttempts[[AOn]])) {
          .self$FailedAttempts[[AOn]] <- list();
          .self$FailedAttempts[[AOn]][[1]] <- 
            .self$MyAttempts[[AOn]][[1]];            
        } else {
          .self$FailedAttempts[[AOn]][[length(.self$FailedAttempts[[AOn]])+1]] <- 
            .self$MyAttempts[[AOn]][[1]];
        }
        .self$MyAttempts[[AOn]] <- list();
        names(.self$MyAttempts)[AOn] <- ANameFunction;
        names(.self$FailedAttempts)[AOn] <- ANameFunction;
      } else {
        NewMyAttempts <- vector("list", length(.self$MyAttempts[[AOn]])-1);
        if (MyMoveOut > 1) {
          for (jj in 1:(MyMoveOut-1)) {
            NewMyAttempts[[jj]] <- .self$MyAttempts[[AOn]][[jj]];
          }
        }
        if (MyMoveOut < length(.self$MyAttempts[[AOn]])) {
          for (jj in (MyMoveOut+1):length(.self$MyAttempts[[AOn]])) {
            NewMyAttempts[[jj-1]] <-  .self$MyAttempts[[AOn]][[jj]];
          }
        }
       if (is.null(.self$FailedAttempts[[AOn]])) {
          .self$FailedAttempts[[AOn]] <- list();
          .self$FailedAttempts[[AOn]][[1]] <- 
            .self$MyAttempts[[AOn]][[1]];            
        } else {
          .self$FailedAttempts[[AOn]][[length(.self$FailedAttempts[[AOn]])+1]] <- 
            .self$MyAttempts[[AOn]][[1]];
        }
        .self$MyAttempts[[AOn]] <- NewMyAttempts;
        names(.self$MyAttempts)[AOn] <- ANameFunction;
        names(.self$FailedAttempts)[AOn] <- ANameFunction;
        eval(parse(text=GetG0Text("TooManyFailsPerOneSimFunction", "globalenv()", S=1)));
      }
    }
      try(.self$TotalWorking <- as.integer(.self$TotalWorking-1));  
      
      MyOnDir <-  paste(SSL$SmallContainingDir, sep="");
      if (verbose >= 3) {
        AFilePrint(paste("About to put functions in ", ANameFunction, ".", sep=""));
      }
      OldOldwd <- getwd(); 
      if (verbose >= 2) {
            AFilePrint("FailAFunction: End of Fail, about to set to Small ContainingDir ");
      } 
      eval(parse(text=GetG0Text("CurrentSmallContainDir", "globalenv()", S=1)));
      if (CurrentSmallContainDir != SSL$SmallContainingDir) {
        AFilePrint(paste("FailAFunction: uh oh, we have currentSmall = ",
          CurrentSmallContainDir, sep=""));
        AFilePrint(paste("FailAFunction: but SSL$SmallContain = ",
          SSL$SmallContainingDir, sep=""));
      }
      if (NoWrite==FALSE) {
      try(setwd(SSL$SmallContainingDir));
      if (verbose >= 2) {
        AFilePrint("FailAFunction: we just set SmallContainDir");
      }
      while(TwoSimR5:::LockMeIn(verbose=0, quoteMore = "Fail End Function", 
        LFDir = ANameFunction,
        MaxLockTime = DefaultMaxLockTime, MoreLock = "", 
        UniqueProcessIdentifier = UniqueProcessIdentifier, 
        NoSpaceLFDir=TRUE) == FALSE) {
        for (ii in 1:100) {}
      } 
      if (verbose >= 2) {
        AFilePrint("FailAFunction: we just set locked in SmallContainDir");
      }
      try(dir.create(paste(ANameFunction, "//MyAttempts", sep=""), 
        showWarnings=FALSE, recursive=TRUE));
      RDTry <- NULL;
      if (any(unlist(list.files(paste(ANameFunction, "//MyAttempts", sep=""))) 
        == paste(UniqueProcessIdentifier, ".txt", sep=""))) {
        try(RDTry <- secure.read.table(paste(ANameFunction, "//MyAttempts//",
          UniqueProcessIdentifier, ".txt", sep=""), nCols=7, header=TRUE));
      }
      RealTime <- GetRealTime();
      if (is.null(RDTry) || NROW(RDTry) <= 0) {
        try(NewTable <- matrix(0, 1, 7));
        try(rownames(NewTable) <- .self$UniqueSimulationIdentifier);
        try(NewTable[1,1:2] <- c(RealTime$JustSeconds, RealTime$JustDays))
        try(NewTable[1,3] <- AOn);
        try(NewTable[1,4] <- .self$AttemptsOnAFunction[AOn]);
        if (.self$AreKilled[AOn] == TRUE) {
          try(NewTable[1,5] <- as.integer(-666));
        } else {
          try(NewTable[1,5] <- as.integer(-1));            
        }
        try(NewTable[1,6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));
        try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction, 
          "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
          sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", 
            "AttemptsOnAFunction", "AreCompleted",
            "EndSeconds", "EndDays"),
          row.names=.self$UniqueSimulationIdentifier));
      } else {
        if (any(rownames(RDTry) == .self$UniqueSimulationIdentifier &
          RDTry[,3] == AOn)) {
            try(ATT <- (1:length(rownames(RDTry)))[rownames(RDTry) == .self$UniqueSimulationIdentifier
              & RDTry[,3] == AOn]);
            try(NewTable <- RDTry);
            if (NewTable[ATT,5] <= 0) {
              try(NewTable[ATT,5] <- NewTable[NROW(NewTable),5]-1);
            } else {
              MyEP <- paste("***************************************************
                Hey, we wanted to save fail to NewTable Function (Never Complete) AOn = ", AOn, "
                but we found that it was already solved with ATT = ", ATT, "
                and NewTable[ATT,] <- ", paste(NewTable[ATT,], collapse=", "),"
                What is going on?
              ", sep="");
              AErrorPrint(MyEP);
              AFilePrint(MyEP);
              return(-1);
            }
            try(NewTable[ATT,6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));
        } else {
            try(NewTable <- matrix(0, NROW(RDTry)+1,7));
            if (NROW(RDTry) == 1) {
              try(NewTable[1,1:7] <- as.numeric(unlist(RDTry[1, 1:7])));
            } else {
              try(NewTable[1:NROW(RDTry),1:7] <- as.numeric(unlist(RDTry[1:NROW(RDTry), 1:7])));            
            }
          try(NewTable[NROW(NewTable),1:2] <- c(RealTime$JustSeconds, RealTime$JustDays));
          try(NewTable[NROW(NewTable),6:7] <- c(RealTime$JustSeconds, RealTime$JustDays));
          try(NewTable[NROW(NewTable),3] <- AOn);
          try(NewTable[NROW(NewTable), 4] <- .self$AttemptsOnAFunction[AOn]);
          if (.self$AreKilled[AOn] == TRUE) {
            try(NewTable[NROW(NewTable),5] <- -666);
          } else {
            try(NewTable[NROW(NewTable),5] <- NewTable[NROW(NewTable),5]-1);
          }
          try(rownames(NewTable) <- c(rownames(RDTry), .self$UniqueSimulationIdentifier)); 
        }
        try(TwoSimR5:::secure.write.table(NewTable, file=paste(ANameFunction,
          "//MyAttempts//", UniqueProcessIdentifier, ".txt", sep=""), 
          sep=", ", col.names=c("StartSeconds", "StartDays", "FunctionNum", 
          "AttemptsOnAFunction", "AreCompleted",
          "EndSeconds", "EndDays"),
          row.names=rownames(NewTable) )); 
      } 
      try(TwoSimR5:::UnLockMe(LFDir = ANameFunction, 
        verbose=verbose, quoteMore="Final Saw no problem Sizeup",
        NoSpaceLFDir=TRUE));
      try(setwd(OldOldwd));
      }
      ##try(.self$Completed <- as.logical(FALSE)); 
      return(8)
      }
  )
);
  
MySummarySimulationAttempt = setRefClass("MySummarySimulationAttempt",
  fields = list(UniqueProcessIdentifier="character",
  NameFunction="character",
  StartTime="numeric",
  IFailed="logical",
  Finished="logical"),
  methods=list(
    initialize= function(UniqueProcessIdentifier="DefaultProcess", NameFunction="DefaultNameFunction"){
      try(.self$UniqueProcessIdentifier <- as.character(UniqueProcessIdentifier));
      ARealTime <- GetRealTime();
      try(.self$NameFunction <- NameFunction);
      try(.self$StartTime <- c(ARealTime$AllSeconds, 
        ARealTime$JustSeconds, ARealTime$JustDays));
      try(.self$Finished <- FALSE);
      try(.self$IFailed <- FALSE);
      return(.self);
    } 
  )
);
    