################################################################################
##AccesOperatingSimulation.r
##
##  (c) 2009-2019 by Alan Lenarcic
##    Written in lab of William Valdar, UNC Genetics
##    used to conduct simulations on UNC Killdevil LSF server to perform
##    simulations in many parameter settings and save material to
##    disk setup on server.
##
## The following does not seem to be used and was "FALSE"d out before
##   the writing of this comment and was not active in simulations.  Alan
## 
##
if (FALSE) {
LockSMSFile <- function(SMS, CurrentContainDir, OnFunction, verbose=0) {
  if (is.logical(verbose) && verbose == TRUE) {
    verbose = 1;
  } else if (is.logical(verbose) && verbose == FALSE) {
    verbose = 0;
  }
  MySaveFile <- paste("SaveOn", SMS$UniqueSimulationIdentifier, ".RData");
  AllFiles <- unlist(list.files(paste(CurrentContainDir, "//", OnFunction, sep="")));
  if (MySaveFile %in% AllFiles) {
    try(load(paste(CurrentContainDir, "//", OnFunction, "//", MySaveFile, sep="")));
    if (!exists("MyFileAttemptsList")) {
      AFilePrint(paste("Hey Error in LockSMSFile OnFunction = ", OnFunction, "
         and SMS = ", SMS$UniqueSimulationIdentifier, sep=""));
      AFilePrint("Could Not Find OnGoes");
      AErrorPrint("LockSMSFile: Not happy with a load attempt error!");
      tryCatch("LockSMSFile: Treating this as an error;")
    } 
    OnGoes <- length(MyFileAttemptsList)+1;
    OnGoes <- OnGoes+1;
  }  else {
    MyFileAttemptsList <- list();  OnGoes = 1;
  }
  ARealTime <- GetRealTime();
  MyTime <- c(ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays);
  MyFileAttemptsList[[OnGoes]] <- list()
}   
}
if (FALSE) {

LookForNewSimulation <- function(verbose = NULL, quoteMore = "LookForNewSimulation : ",
  ForceFunction = NULL) {
  eval(parse(text=GetG0Text("CurrentContainDir", S=1)));
  if (is.numeric(CurrentContainDir) && CurrentContainDir[1] == 0) {
    eval(parse(text=GetG0Text("CurrentContainDir", S=1)));
    eval(parse(text=GetG0Text("n", S=1)));
    eval(parse(text=GetG0Text("p", S=1)));
    eval(parse(text=GetG0Text("k", S=1)));
    eval(parse(text=GetG0Text("sigma", S=1)));
    eval(parse(text=GetG0Text("NumReps", S=1)));
    eval(parse(text=GetG0Text("LARSSeekFlag", S=1)));
    eval(parse(text=GetG0Text("DefaultGenerateBetaVec", S=1)));
    eval(parse(text=GetG0Text("DefaultCorrelationXMatrix", S=1)));
    eval(parse(text=GetG0Text("AllDir", S=1)));
    CurrentContainDir <- paste(AllDir, "//",  DirNamingClenture( 
      DefaultGenerateBetaVec, DefaultCorrelationXmatrix,
      NLen = n, kLen = p, kActiveLen = k, SigmaNoise = sigma, 
      NumReps = NumReps, LARSSeekFlag  = LARSSeekFlag), sep="");
      dir.create(CurrentContainDir, showWarnings = FALSE);   
    eval(parse(text=SetGText("CurrentContainDir", S=1))); 
  }
  if (is.null(verbose)) { 
    eval(parse(text=GetG0Text("verbose")));
    if (verbose == 0) {
      verbose = 0;
    } 
  } else { verbose = as.numeric(verbose); }
  MySummaryPerformanceDirectory = paste(CurrentContainDir, "//", "SumPerformance", sep="");
  dir.create(MySummaryPerformanceDirectory, showWarnings=FALSE, recursive=TRUE);
  if (verbose > 0) {
    AFilePrint("Looking for New Simulation Start with LockMeIn"); flush.console(); 
  }
  while(LockMeIn(verbose=as.numeric(verbose), quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory)==FALSE){}
  MyListfiles = list.files(MySummaryPerformanceDirectory);
  MyGG = NULL;
  try(MyGG <- GetSummaryPerformanceFile(CurrentContainDir, verbose = max(verbose -1,0), 
    quoteMore = quoteMore, AlreadyLocked=TRUE));
  if (verbose > 0) {
    AFilePrint("Look for a New Simulation: Got Summary Performance File "); flush.console(); 
  }
  if (is.null(MyGG) || (is.numeric(MyGG) && MyGG[1] == -1)) {
     AFilePrint(paste("LookForNewSimulation: SummaryPerformanceFile returned an error ", sep=""));
     if (is.null(MyGG)) {
       AFilePrint(paste("LookForNewSimulation: MyGG is NULL!")); flush.console();
     }
     flush.console();
     return(-1);
  } 
  eval(parse(text=GetG0Text("TargetTotalSimulations", "globalenv()", S=1)));

  DDR = NROW(MyGG$isDira);
  MyListfiles = MyListfiles[substr(MyListfiles, 1, nchar("SMS")) == "SMS", drop=FALSE]
  MyListSMS = MyListfiles[MyListfiles!="SummaryPerformance.txt" & MyListfiles != "LockFile.txt", drop=FALSE];
  MyListSMS = substr(MyListSMS, 4, nchar(MyListSMS)-6);
  CountOccupied <-  NULL;
  if (DDR > 0) {
    CountOccupied <- rep(0, NCOL(MyGG$isDira));
    CountFailed <- rep(0, NCOL(MyGG$isDira));
    ### Check to see that some algorithm hasn't generated too many fails that it is untrustworth.
    for (ii in 1:NCOL(MyGG$isDira)) {
      RT = length(MyGG$isDira[substr(MyGG$isDira[,ii],1,1) %in% c("-"),ii]);
      CountOccupied[ii] <- RT;
      CountFailed[ii] = length(MyGG$isDira[substr(MyGG$isDira[,ii],1,1) %in% c("D"),ii]);
    }
    MMRT <- sort(CountOccupied, decreasing=TRUE)[2];
    for (ii in 1:NCOL(MyGG$isDira)) {
     RT = CountOccupied[ii];
     if (RT-MMRT + CountFailed[ii] > TooManyFails) {
       MyGG$isDira[substr(MyGG$isDira[,ii],1,1) == "0",ii] = 
         paste("D", UniqueProcessIdentifier, sep="");  
       NewTable = cbind(MyGG$iDList, MyGG$isDira);
         colnames(NewTable) = c("iD", colnames(MyGG$isDira));
       filename = paste(CurrentContainDir, "//SumPerformance//", "SummaryPerformance.txt", sep="");
       asep = ",";
       Goal = NULL;
       ##try(eval(parse(text=
       ##  "try(write.table(NewTable, file = filename,
       ##  append=FALSE, quote=FALSE, sep=asep, col.names=colnames(NewTable), row.names=FALSE);
       ##  Goal = 1;")));
       MySummaryPerformanceTable <- NewTable;
       try(save(MySummaryPerformanceTable = MySummaryPerformanceTable, 
        file=paste(CurrentSmallContainDir, "//SumPerformance//", 
        "SummaryPerformanceRData.RData", sep="")));
       if (is.null(Goal)) {
          AFilePrint(paste("Write to NewTable fails for ii = ", ii, sep=""));
          AFilePrint(paste("filename was ", filename, sep=""));  flush.console();
       }
     }
  }
  }
  if (DDR > 0) {
    ## Go through active "SMS" files, assess which can be deleted from directory
    if (length(MyListSMS) > 0) {
      WentWithLines = MyGG$iDList[MyGG$iDList %in% MyListSMS, drop=FALSE];
      if (length(WentWithLines) > 0) {
        for (ii in 1:length(WentWithLines)) {
          OnLine = (1:length(MyGG$iDList))[ MyGG$iDList == WentWithLines[ii]];
          if (all(substr(MyGG$isDira[OnLine,],1,1) %in% c("+", "D"))) {
            unlink(paste(CurrentContainDir, "//SumPerformance//", "SMS", WentWithLines[ii], ".RData", sep="")); 
          } else if (any(substr(MyGG$isDira[OnLine,],1,1) == "0")) {
            ## This loads the relevant SMS file
            ## Then it changes Summary table to have "-UniqueProcessIdentifier" for the chosen
            ##   sample function;  
            load(paste(CurrentContainDir, "//SumPerformance//", "SMS", WentWithLines[ii], ".RData", sep=""));
            eval(parse(text=SetGText("SMS", S=1)));
            eval(parse(text=GetG0Text("NameFunctions",S=1)));
            INeed = (1:length(MyGG$isDira[OnLine,]))[ substr(MyGG$isDira[OnLine,],1,1) == "0"];
            ISample = (INeed)[sample(1:length(INeed), size=1)];
            MyGG$isDira[OnLine,ISample] = paste("-", UniqueProcessIdentifier, sep="");
            NewTable = cbind(MyGG$iDList, MyGG$isDira);
            colnames(NewTable) = c("iD", colnames(MyGG$isDira));
            Goal = NULL; asep =","; 
            filename = paste(MySummaryPerformanceDirectory, "//SummaryPerformance.txt", sep="")
            try(eval(parse(text="write.table(NewTable, filename,
              sep=asep, col.names=colnames(NewTable), row.names=FALSE, quote=FALSE, append=FALSE);
              Goal = 1")));
            if (is.null(Goal)) {
              AFilePrint(paste(" We tried with Went WithLines ii =  ", ii, " but failed", sep=""));
              AFilePrint(paste(" We're trying to write to zero")); flush.console();
            }
        
            LockedSMSFile <- LockInSMS(SMS, CurrentContainDir = CurrentContainDir, OnFunction=NameFunctions[ISample]);              
            UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory);
            eval(parse(text=SetGText("ISample",S=1)));
            return(list(SMS=SMS,MyNameFunction = NameFunctions[ISample],ISample=ISample, LockedSMSFile=LockedSMSFile));
          }
        }
      }
    }    
  }
  if (length(MyGG$isDira) >= 1) {
    MyGGN = names(MyGG$isDira);
    for (tii in 1:length(MyGG$isDira)) {
       MyGG$isDira[tii] = paste(unlist(strsplit(as.character(MyGG$isDira[tii]), " ")), collapse=""); 
    }
    names(MyGG$isDira) = MyGGN;
  }
  if (DDR > 0 && any(MyGG$iDList == "0") && any(MyGG$iDList != "0")) {
    RTN = (1:(length(MyGG$iDList)))[MyGG$iDList != "0"];
    NewTable = cbind(MyGG$iDList[RTN], MyGG$isDira[RTN,]);
    colnames(NewTable) = c("iD", colnames(MyGG$isDira));
    SumPerformanceFile = paste(CurrentContainDir, "//SumPerformance//", "SummaryPerformance.txt", sep="");
    asep = ",";
    Goal = NULL;
    try(eval(parse(text="write.table(NewTable, file = SumPerformanceFile,
      append=FALSE, quote=FALSE, sep=asep, col.names=colnames(NewTable), 
      row.names=FALSE); Goal = 1")));
    if (is.null(Goal)) {
       AFilePrint("DDR > 0, we tried to write to SumPerformanceFile but failed ");
       AFilePrint("  New Table was :");
       AFilePrint(NewTable);
       flush.console();
       return(-1);
    }
    WorkingRow = length(MyGG$iDList[RTN]) + 1;
  } else {
    WorkingRow = length(MyGG$iDList) +1; 
  }
  eval(parse(text=SetGText("WorkingRow", S=1)));
  
  eval(parse(text=GetG0Text("n", S=1)));
  eval(parse(text=GetG0Text("p", S=1)));
  eval(parse(text=GetG0Text("k", S=1)));
  eval(parse(text=GetG0Text("sigma", S=1)));
  eval(parse(text=GetG0Text("GenerateBetaVec", S=1)));
  eval(parse(text=GetG0Text("CorrelationXmatrix", S=1)));
  eval(parse(text=GetG0Text("tNoiseDF", S=1)));
  eval(parse(text=GetG0Text("LogitRegression", S=1)));
  eval(parse(text=GetG0Text("Beta0",S=1)));
  eval(parse(text=GetG0Text("ExperimentName", S=1)));
  eval(parse(text=GetG0Text("jobii", S=1)));
  eval(parse(text=GetG0Text("WorkingRow",S=1)));
  eval(parse(text=GetG0Text("TargetTotalSimulations",S=1)));
  if (WorkingRow > TargetTotalSimulations) {
    UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory);
    return(NULL); 
  }
  eval(parse(text=GetG0Text("NameFunctions",S=1)));
  INeed = 1:length(NameFunctions);
  ISample = sample(INeed, size=1);
  SMS = NULL;
  if (verbose > 1) {
    AFilePrint("Look for a new simulation, had to configure a new SMS"); flush.console(); 
  }
  if (exists("DropHome") && substr(DropHome, 1, nchar("/Users/")) == "/Users/") {
    try(system("purge"));
  }
  try(SMS <- SimMeData(n = n, p = p, k = k, sigma = sigma,
     GenerateBetaVec = GenerateBetaVec, 
     CorrelationXmatrix = CorrelationXmatrix, tNoiseDF = tNoiseDF,
     LogitRegression = LogitRegression, Beta0 = Beta0, 
     ExperimentName=ExperimentName, jobii= jobii, WorkingRow = WorkingRow,AlreadyLocked=TRUE,
     UniqueProcessIdentifier=UniqueProcessIdentifier, ISample=ISample));
  if (is.null(SMS)) {
    AFilePrint("Crazy Dude, SMS SimMeData failed with a NULL "); flush.console();
    UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory);
    return(-1);
  }
  LockedSMSFile <- LockInSMS(SMS, CurrentContainDir = CurrentContainDir, OnFunction=NameFunctions[ISample]);    
  UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory); 
  return(list(SMS=SMS, MyNameFunction=NameFunctions[ISample], ISample=ISample, LockedSMSFile=LockedSMSFile));
}
}