
        
UpdateNameFunctions <- function() {
  if (OnSimType %in% c("Group", "GroupRobit", "GroupLogit") || grepl("Group", OnSimType)[1] == TRUE) {
    if (exists("DefaultAllGroupFunctions", env=globalenv())) {
       eval(parse(text=GetG0Text("DefaultAllGroupFunctions", "globalenv()", S=1)));
       eval(parse(text=GetG0Text("DefaultAllGroupNameFunctions", "globalenv()", S=1)));
       UseFunctions <-  DefaultAllGroupFunctions;
       NameFunctions <- DefaultAllGroupNameFunctions;
    } else if (exists("DefaultAllGroupFunctions", env=TWOSIMENVIRONMENT)) {
      eval(parse(text=GetG0Text("DefaultAllGroupFunctions", "TWOSIMENVIRONMENT)", S=1)));
      eval(parse(text=GetG0Text("DefaultAllGroupNameFunctions", "TWOSIMENVIRONMENT", S=1)));
      UseFunctions <- DefaultAllGroupFunctions;
      NameFunctions <- DefaultAllGroupNameFunctions;   
    }

  } else  {
        if (exists("DeclareAllDefaultFunctions", env=globalenv())) {
       eval(parse(text=GetG0Text("DefaultAllUseFunctions", "globalenv()", S=1)));
       eval(parse(text=GetG0Text("DefaultAllNameFunctions", "globalenv()", S=1)));
       UseFunctions <-  DefaultAllUseFunctions;
       NameFunctions <- DefaultAllNameFunctions;
    } else if (exists("DeclareAllDefaultFunctions", env=TWOSIMENVIRONMENT)) {
      eval(parse(text=GetG0Text("DefaultAllUseFunctions", "TWOSIMENVIRONMENT)", S=1)));
      eval(parse(text=GetG0Text("DefaultAllNameFunctions", "TWOSIMENVIRONMENT", S=1)));
      UseFunctions <- DefaultAllUseFunctions;
      NameFunctions <- DefaultAllNameFunctions;   
    }
  }
  eval(parse(text=GetG0Text("RenameFunctions", "globalenv()", S=1)));
  ##eval(parse(text=GetG0Text("NameFunctions", "globalenv()", S=1)));
  ##eval(parse(text=GetG0Text("UseFunctions", "globalenv()", S=1)));
  if ( (length(RenameFunctions) == 1 && is.numeric(RenameFuncions) && RenameFunctions[1] == 0) ||
     is.null(RenameFunctions)) {
    eval(parse(text=SetGText("UseFunctions", "globalenv()", S=1)));
    eval(parse(text=SetGText("NameFunctions", "globalenv()", S=1))); 
    return(1);       
  }
  if (length(RenameFunctions) >= 1  && !is.numeric(RenameFunctions) &&
    RenameFunctions[1] != 0) {
    MyMatch <- match(RenameFunctions, NameFunctions);    
    MyKeep <- sort(MyMatch[!is.na(MyMatch)]);
    if (length(MyKeep) >= 1) {
        print(paste("Reducing NameFunctions in size from length ", length(NameFunctions),
          " to length ", length(MyKeep), sep="")); flush.console();
      try(NameFunctions <- NameFunctions[MyKeep]);
      try(UseFunctions <- UseFunctions[MyKeep]);
      eval(parse(text=SetGText("NameFunctions", "globalenv()", S=1)));
      eval(parse(text=SetGText("UseFunctions", "globalenv()", S=1)));
      print(paste("-- We have successfully limited UseFunctions to length ", length(UseFunctions), sep=""));
      flush.console(); return;
    } else {
      print(paste("UpdateNameFunctions: Oops, well, Although RenameFunctions Exists, we don't have any matches!", sep=""));
      return;
    }
  }
}
MasterGenerateANewSimulation <- function() {                                                                                                    
  eval(parse(text=GetG0Text("n", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("p", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("k", "globalenv()", S=1)));
  AlreadyLocked <- FALSE;

  eval(parse(text=GetG0Text("sigma", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("tNoiseDF", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("GenerateBetaVec", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("CorrelationXmatrix", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("tNoiseDF", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("LogitRegression","globalenv()",  S=1)));
  eval(parse(text=GetG0Text("Beta0", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("ExperimentName", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("jobii", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("WorkingRow","globalenv()", S=1)));
  eval(parse(text=GetG0Text("TargetTotalSimulations","globalenv()", S=1)));
  eval(parse(text=GetG0Text("NameFunctions","globalenv()", S=1)));
  eval(parse(text=GetG0Text("OnSimType", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("DuplicationEETest", "globalenv()", S=1)));
  INeed = 1:length(NameFunctions);
  ISample = sample(INeed, size=1);
  SMS = NULL;
  if (verbose > 1) {
    AFilePrint("Look for a new simulation, had to configure a new SMS"); flush.console(); 
  }
  eval(parse(text=GetG0Text("ALargeContainDir", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("ASmallContainDir", "globalenv()", S=1)));
  MySummaryPerformanceDirectory <-  paste(ALargeContainDir, 
  "//SumPerformance", sep="");
  if (n >= 200 || p >= 1000) {
    try(gc());
     if (exists("DropHome") && substr(DropHome, 1, nchar("/Users/")) == "/Users/") {
       try(system("purge"));
     } 
  }
  eval(parse(text=GetG0Text("ElseGetX", "globalenv()", S=1)));
  if (exists("DuplicationEETest") && is.numeric(DuplicationEETest) && DuplicationEETest[1] == 0) {
    DuplicationEETest <- FALSE;
    eval(parse(text=SetGText("DuplicationEETest", "globalenv()", S=1)));
  } 
  if (is.logical(DuplicationEETest) && DuplicationEETest == TRUE) {
    OnSimType == "DuplicationEE";
  }
  if (OnSimType == "Group") {
    eval(parse(text=GetG0Text("GroupSize", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("pGroups", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("kGroups", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("tNoiseDF", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("sigma", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("SigmaNoise", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("Beta0", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("GenerateBetaGroupVec", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("LogitRegression", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("MeanCenterXX", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("jobii", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("ElseGetX", "globalenv()", S=1)));
    SMS <- SimGroupData(n = n, pGroups = pGroups, GroupSize=GroupSize,
     kGroups = kGroups, sigma = sigma,  ElseGetX = ElseGetX,
     SigmaNoise = SigmaNoise, tNoiseDF = tNoiseDF, GenerateBetaGroupVec = GenerateBetaGroupVec, 
     CorrelationXmatrix = CorrelationXmatrix,
     LogitRegression = FALSE, Beta0 = Beta0, ExperimentName="", jobii= jobii, 
     WorkingRow = WorkingRow, AlreadyLocked=TRUE, MeanCenterXX = MeanCenterXX,
     UniqueProcessIdentifier = UniqueProcessIdentifier, ISample = ISample,
     DontSave=TRUE)
    
  } else if (OnSimType == "Robit") {
     eval(parse(text=GetG0Text("dfRobit")));
     try(SMS <- SimMeData(n = n, p = p, k = k, sigma = sigma,
     GenerateBetaVec = GenerateBetaVec, ElseGetX = ElseGetX,
     CorrelationXmatrix = CorrelationXmatrix, tNoiseDF = tNoiseDF,
     LogitRegression = FALSE, Beta0 = Beta0, 
     ExperimentName=ExperimentName, jobii= jobii, WorkingRow = WorkingRow,
     AlreadyLocked=TRUE, dfRobit=dfRobit,
     UniqueProcessIdentifier=UniqueProcessIdentifier, ISample=ISample,
     DontSave = TRUE));
  } else if (OnSimType == "GroupRobit") {
    eval(parse(text=GetG0Text("GroupSize", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("pGroups", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("kGroups", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("tNoiseDF", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("sigma", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("dfRobit")));
    eval(parse(text=GetG0Text("SigmaNoise", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("Beta0", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("GenerateBetaGroupVec", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("LogitRegression", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("MeanCenterXX", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("jobii", "globalenv()", S=1)));
    SMS <- SimGroupData(n = n, pGroups = pGroups, GroupSize=GroupSize,
     kGroups = kGroups, sigma = sigma,  ElseGetX = ElseGetX,
     SigmaNoise = SigmaNoise, tNoiseDF = tNoiseDF, GenerateBetaGroupVec = GenerateBetaGroupVec, 
     CorrelationXmatrix = CorrelationXmatrix,  dfRobit=dfRobit,
     LogitRegression = FALSE, Beta0 = Beta0, ExperimentName="", jobii= jobii, 
     WorkingRow = WorkingRow, AlreadyLocked=TRUE, MeanCenterXX = MeanCenterXX,
     UniqueProcessIdentifier = UniqueProcessIdentifier, ISample = ISample,
     DontSave=TRUE)
  } else if (OnSimType == "Logit") {
     try(SMS <- SimMeData(n = n, p = p, k = k, sigma = sigma,
     GenerateBetaVec = GenerateBetaVec,  ElseGetX = ElseGetX,
     CorrelationXmatrix = CorrelationXmatrix, tNoiseDF = tNoiseDF,
     LogitRegression = TRUE, Beta0 = Beta0, 
     ExperimentName=ExperimentName, jobii= jobii, WorkingRow = WorkingRow,
     AlreadyLocked=TRUE,
     UniqueProcessIdentifier=UniqueProcessIdentifier, ISample=ISample,
     DontSave = TRUE));
  } else if (OnSimType == "GroupLogit") {
    eval(parse(text=GetG0Text("GroupSize", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("pGroups", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("kGroups", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("tNoiseDF", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("sigma", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("SigmaNoise", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("Beta0", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("GenerateBetaGroupVec", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("LogitRegression", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("CorrelationXmatrix", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("MeanCenterXX", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("jobii", "globalenv()", S=1)));
    SMS <- SimGroupData(n = n, pGroups = pGroups, GroupSize=GroupSize,
     kGroups = kGroups, sigma = sigma,
     SigmaNoise = SigmaNoise, tNoiseDF = tNoiseDF, GenerateBetaGroupVec = GenerateBetaGroupVec, 
     CorrelationXmatrix = CorrelationXmatrix,  ElseGetX = ElseGetX,
     LogitRegression = TRUE, Beta0 = Beta0, ExperimentName="", jobii= jobii, 
     WorkingRow = WorkingRow, AlreadyLocked=TRUE, MeanCenterXX = MeanCenterXX,
     UniqueProcessIdentifier = UniqueProcessIdentifier, ISample = ISample,
     DontSave=TRUE)
  } else if (OnSimType == "DuplicationEE") {
    eval(parse(text=SetGText("OnSimType", "globalenv()", S=1)));
    try(SMS <- SimMeData(n = n, p = p/2, k = k, sigma = sigma,
     GenerateBetaVec = GenerateBetaVec,    ElseGetX = ElseGetX,
     CorrelationXmatrix = CorrelationXmatrix, tNoiseDF = tNoiseDF,
     LogitRegression = LogitRegression, Beta0 = Beta0, 
     ExperimentName=ExperimentName, jobii= jobii, WorkingRow = WorkingRow,
     AlreadyLocked=TRUE,
     UniqueProcessIdentifier=UniqueProcessIdentifier, ISample=ISample,
     DontSave = TRUE));
    SMS <- DoubleSimulation(SMS, AlreadyLocked=TRUE, DontSave=TRUE);
  } else {
  try(SMS <- SimMeData(n = n, p = p, k = k, sigma = sigma,
     GenerateBetaVec = GenerateBetaVec,    ElseGetX = ElseGetX,
     CorrelationXmatrix = CorrelationXmatrix, tNoiseDF = tNoiseDF,
     LogitRegression = LogitRegression, Beta0 = Beta0, 
     ExperimentName=ExperimentName, jobii= jobii, WorkingRow = WorkingRow,
     AlreadyLocked=TRUE,
     UniqueProcessIdentifier=UniqueProcessIdentifier, ISample=ISample,
     DontSave = TRUE));
   }


   eval(parse(text=GetG0Text("OriginalOldwd", "globalenv()", S=1)));
   try(setwd(OriginalOldwd));
   return(SMS);
}


StackProgresses <- function(ProgressNames) {
    MyOut <- matrix(0, length(ProgressNames), length(ProgressNames[[1]]));
    for (ii in 1:length(ProgressNames)) {
      MyOut[ii,] <- ProgressNames[[ii]];
    }
   return(MyOut);
}

MakeDirSaves <- function() {
  eval(parse(text=GetG0Text("OnSimType", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("n", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("p", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("k", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("pGroups", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("kGroups", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("sigma", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("piSV", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("sigSV", "globalenv()", S=1)));  
  eval(parse(text=GetG0Text("ADDDIRCONT", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("argeese", "globalenv()", S=1)));  
  
  eval(parse(text=GetG0Text("TotalSlaves", "globalenv()", S=1)));  
  eval(parse(text=GetG0Text("TotalThralls", "globalenv()", S=1)));  
  if (!is.null(TotalSlaves) && length(TotalSlaves) == 1 && TotalSlaves != 0
    && TotalSlaves > TotalThralls) {
    AFilePrint(paste("  We need to replace TotalThralls with TotalSlaves = ", TotalSlaves, sep=""));
    TotalThralls <- TotalSlaves;
    eval(parse(text=SetGText("TotalThralls", "globalenv()", S=1)));      
  }
  
  if (length(list.files("c:/Users/AlanX230")) >= 1) {
     DirToSave <- "c:/Stat/TwoSimR52";
     ASmallContainDir <- "c:/Stat/TwoSimR52/Small";
     ALargeContainDir <- "c:/Stat/TwoSimR52/Large";
     RunDirFile <- paste(DropHome, "//TwoLassoPackage//RunSeekScriptFile.r", sep="");
     APrintDir <- "c:/Stat/TwoSimR52/Print";
     AErrorDir <- "c:/Stat/TwoSimR52/Error";
     CrashDir <- "c:/Stat/Crash";
  } else if (length(list.files("~/SpikeTest")) >= 1) {
     DirToSave <- "~/TR5Script";
     RunDirFile <- paste("~/SpikeTest//RunSeekScriptFile.r", sep="");
     ASmallContainDir <- "/netscr/alenarc/TwoSimR52";
     if (length(list.files("/lustre/scr/a/l/alenarc")) <= 0) {
       ALargeContainDir <- "/netscr/alenarc/TwoSimR52/Large";        
     } else {
       ALargeContainDir <- "/lustre/scr/a/l/alenarc/TwoSimR52";
     }
     APrintDir <- "/netscr/alenarc/TwoSimR5Print";
     AErrorDir <- "/netscr/alenarc/TwoSimR5Error"; 
     CrashDir <- "/netscr/alenarc/TwoSimR52/Crash"; 
  } else if (length(unlist(list.files("/Users/lenarcic/Documents/"))) > 0) {
    DirToSave <-  "/Users/lenarcic/Documents/TwoSimR52";
    dir.create(DirToSave, showWarnings=FALSE, recursive=TRUE);
    RunDirFile <- paste(DropHome, "//TwoLassoPackage//RunSeekScriptFile.r", sep="");
    APrintDir <- "/Users/lenarcic/Documents/TwoSimR52//Print"
    AErrorDir <- "/Users/lenarcic/Documents/TwoSimR52//Error"
    ALargeContainDir <- "/Users/lenarcic/Documents/TwoSimR52/Large";
    ASmallContainDir <- "/Users/lenarcic/Documents/TwoSimR52/Small";   
    CrashDir <- "/Users/lenarcic/Documents/TwoSimR52/Crash" 
  }
  eval(parse(text=GetG0Text("argeese", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("ADDDIRCONT", "globalenv()", S=1)));
  if (is.numeric(ADDDIRCONT)) {
    ADDDIRCONT <- "DefaultDir";
  }
  if (is.character("argeese")) {
    AV <- unlist(strsplit(argeese, "\\."))[1];
    if (nchar(AV) <= 0) {
      AV <- ADDDIRCONT
    }
  } else {
    ADDDIRCONT
  }
  if (OnSimType %in% c("Group", "GroupRobit", "GroupLogit")) {
    MyDirToGo <- paste(ADDDIRCONT, "//N", tSeq(n), "pG", tSeq(pGroups),
      "kG", tSeq(kGroups), "sig", tSeq(sigma), "piSV", tSeq(piSV),
      "sigSV", tSeq(sigSV), "Rho", tSeq(CovRho), sep="");
  } else {
    MyDirToGo <- paste(ADDDIRCONT, "//N", tSeq(n), "p", tSeq(p),
      "k", tSeq(k), "sig", tSeq(sigma), "piSV", tSeq(piSV),
      "sigSV", tSeq(sigSV), "Rho", tSeq(CovRho), sep="");    
  }
  if (OnSimType %in% c("Robit", "GroupRobit")) {
    MyDirToGo <- paste(MyDirToGo, "R", sep="");
  } else if (OnSimType %in% c("Logit", "GroupLogit")) {
    MyDirToGo <- paste(MyDirToGo, "L", sep=""); 
  }
  if (OnSimType == "DuplicationEE") {
    MyDirToGo <- unlist(strsplit(MyDirToGo, "//"));
    MyDirToGo <- paste(paste(MyDirToGo[1:(length(MyDirToGo)-1)], sep="", collapse="//"), 
    "//DD", MyDirToGo[length(MyDirToGo)], sep="");
  }
  ALargeContainDir <- paste(ALargeContainDir, "//", MyDirToGo, sep="");
  ASmallContainDir <- paste(ASmallContainDir, "//", MyDirToGo, sep="");

  APrintDir <- paste(APrintDir, "//", MyDirToGo, sep="");
  AErrorDir <- paste(AErrorDir, "//", MyDirToGo, sep="");
  
  dir.create(CrashDir, showWarnings=FALSE, recursive=TRUE);
  dir.create(ALargeContainDir, showWarnings=FALSE, recursive=TRUE);
  dir.create(ASmallContainDir, showWarnings=FALSE, recursive=TRUE);
  dir.create(DirToSave, showWarnings=FALSE, recursive=TRUE);
  dir.create(APrintDir, showWarnings=FALSE, recursive=TRUE);
  dir.create(AErrorDir, showWarnings=FALSE, recursive=TRUE); 
   
  eval(parse(text=SetGText("CrashDir", "globalenv()", S=1))); 


  LLACrash <- unlist(list.files(CrashDir));
  if (any(substr(LLACrash, 1, nchar("core")) == "core")) {
     AK <- LLACrash[substr(LLACrash,1, nchar("core")) == "core"];
     Oldwd <- getwd();
     setwd(CrashDir);
     for (ii in 1:length(AK)) {
        try(unlink(AK[ii]));
     }
     setwd(Oldwd);
  }

  LLASim <- unlist(list.files(DirToSave));
  if (any(substr(LLASim, 1, nchar("core")) == "core")) {
     AK <- LLASim[substr(LLASim,1, nchar("core")) == "core"];
     Oldwd <- getwd();
     setwd(DirToSave);
     for (ii in 1:length(AK)) {
        try(unlink(AK[ii]));
     }
     setwd(Oldwd);
  }

  LLALarge <- unlist(list.files(ALargeContainDir));
  if (any(substr(LLALarge, 1, nchar("core")) == "core")) {
     AK <- LLALarge[substr(LLALarge,1, nchar("core")) == "core"];
     Oldwd <- getwd();
     setwd(ALargeContainDir);
     for (ii in 1:length(AK)) {
        try(unlink(AK[ii]));
     }
     setwd(Oldwd);
  }
  LLASmall <- unlist(list.files(ASmallContainDir));
  if (any(substr(LLASmall, 1, nchar("core")) == "core")) {
     AK <- LLASmall[substr(LLASmall,1, nchar("core")) == "core"];
     Oldwd <- getwd();
     setwd(ASmallContainDir);
     for (ii in 1:length(AK)) {
        try(unlink(AK[ii]));
     }
     setwd(Oldwd);
  }

  LLAPrint <- unlist(list.files(APrintDir));
  if (any(substr(LLAPrint, 1, nchar("core")) == "core")) {
     AK <- LLAPrint[substr(LLAPrint,1, nchar("core")) == "core"];
     Oldwd <- getwd();
     setwd(APrintDir);
     for (ii in 1:length(AK)) {
        try(unlink(AK[ii]));
     }
     setwd(Oldwd);
  }
  LLAError <- unlist(list.files(AErrorDir));
  if (any(substr(LLAError, 1, nchar("core")) == "core")) {
     AK <- LLAError[substr(LLAError,1, nchar("core")) == "core"];
     Oldwd <- getwd();
     setwd(AErrorDir);
     for (ii in 1:length(AK)) {
        try(unlink(AK[ii]));
     }
     setwd(Oldwd);
  }
 
  eval(parse(text=GetG0Text("JustPrintFile")));
  if (is.character(JustPrintFile) && JustPrintFile %in% 
    c("DoNotPrint", "Do Not Print", "donotprint")) {
    JustPrintFile <- "DoNotPrint";      
  } else {
    JustPrintFile <- paste("P", UniqueProcessIdentifier, ".txt", sep="");
  }
  eval(parse(text=GetG0Text("NameFunctions")));
  if (is.null(NameFunctions) || (is.numeric(NameFunctions) && NameFunctions[1] == 0)) {
    DefineDefaultUseFunction();
     eval(parse(text=GetG0Text("RenameFunctions", "globalenv()",S=1)));
     if (!is.null(RenameFunctions) && !is.numeric(RenameFunctions) && length(RenameFunctions) >= 1) {
      UpdateNameFunctions();
      eval(parse(text=GetG0Text("UseFunctions", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("NameFunctions", "globalenv()", S=1)));
    }
  }

  eval(parse(text=SetGText("JustPrintFile", "globalenv()", S=1)));
  eval(parse(text=SetGText("ALargeContainDir", "globalenv()", S=1)));
  eval(parse(text=SetGText("ASmallContainDir", "globalenv()", S=1)));
  eval(parse(text=SetGText("APrintDir", "globalenv()", S=1)));
  eval(parse(text=SetGText("AErrorDir", "globalenv()", S=1)));
  Oldwd <- getwd();
  setwd(APrintDir);
  if (JustPrintFile %in% unlist(list.files(APrintDir))) {
     try(ACCOn <- file(JustPrintFile, "at"));
  } else {
     try(ACCOn <- file(JustPrintFile, "at"));    
  }
  try(writeLines(text="  Start a run  ", con=ACCOn));
  try(close(ACCOn));
  setwd(Oldwd);
  
  
  dir.create(paste(DirToSave, "//", MyDirToGo, sep=""), 
    showWarnings=FALSE, recursive=TRUE);
  AllDirToSave <- paste(DirToSave, "//", MyDirToGo, sep="");
  eval(parse(text=SetGText("AllDirToSave", "globalenv()", S=1)));
  eval(parse(text=SetGText("DirToSave", "globalenv()", S=1)));
  eval(parse(text=SetGText("MyDirToGo", "globalenv()", S=1)));
  eval(parse(text=SetGText("RunDirFile", "globalenv()", S=1)));
  

    setwd(ASmallContainDir);
    eval(parse(text=GetG0Text("TotalThralls", "globalenv()", S=1)));
    if (TotalThralls <= 0) {
      TotalThralls <- length(unlist(list.files("Submissions")));
      eval(parse(text=SetGText("TotalThralls", "globalenv()", S=1)));
    }
    ZZMyID <- TwoSimR5:::AZeroOut(as.numeric(MyID),10^ceiling(log(TotalThralls+1,10)))
    eval(parse(text=SetGText("ZZMyID", "globalenv()", S=1)));
    setwd(Oldwd);
  MyNameOfInputFile <- unlist(strsplit(argeese, "\\."))[1];
  eval(parse(text=SetGText("MyNameOfInputFile", "globalenv()", S=1)));
  return(AllDirToSave);
}

ScriptGetArgs <- function ()
{
    eval(parse(text = GetG0Text("Myargs", S = 1)))
    eval(parse(text = GetG0Text("argeese", S = 1)))
    eval(parse(text = GetG0Text("MyDefaultArgs", S = 1)))
    if (is.null(MyDefaultArgs) || length(MyDefaultArgs) == 1 &&
        MyDefaultArgs[1] == 0) {
        Myargs = commandArgs(TRUE)
    } else {
        Myargs <- MyDefaultArgs
    }
    if (length(Myargs) >= 1 && is.character(Myargs[1]) && !is.na(Myargs[1])) {
        argeese = Myargs[1]

        if (length(Myargs) <= 1) {
            Myargs <- c(argeese, 1, 1, 1, 1, 1, 1)
        }
        DeclareAllDefaultFunctions();               
        eval(parse(text=GetG0Text("DefaultAllGroupNameFunctions", "globalenv()", S=1)));
        eval(parse(text=GetG0Text("DefaultAllNameFunctions", "globalenv()", S=1)));
        AnyNameFunction <- c(DefaultAllGroupNameFunctions, DefaultAllNameFunctions);
        eval(parse(text=SetGText("AnyNameFunction", "globalenv()", S=1)));

        if (is.character(Myargs[2]) && Myargs[2] %in% paste("Generate", AnyNameFunction, sep="")) {
          MandatoryOnFunction <- AnyNameFunction[paste("Generate", AnyNameFunction, sep="") == Myargs[2]];
          ForceFunction <- MandatoryOnFunction;
          print(paste("Setting MandatoryOnFunction to ", Myargs[2], sep=""));
          flush.console();
          if (length(Myargs) == 2) {
            Myargs <- c(argeese, 1,1,1,1,1,1,1);
          } else if (length(Myargs) >= 3) {
            Myargs <- c(Myargs[1], Myargs[3:length(Myargs)],1,1,1,1,1);
          }  
        } else if  (is.character(Myargs[2]) && paste("Generate", Myargs[2], sep="") %in% AnyNameFunction) {
          MandatoryOnFunction <- AnyNameFunction[AnyNameFunction == paste("Generate", Myargs[2], sep="")];
          ForceFunction <- MandatoryOnFunction;
          print(paste("Setting MandatoryOnFunction to ", Myargs[2], sep=""));
          flush.console();
          if (length(Myargs) == 2) {
            Myargs <- c(argeese, 1,1,1,1,1,1,1);
          } else if (length(Myargs) >= 3) {
            Myargs <- c(Myargs[1], Myargs[3:length(Myargs)],1,1,1,1,1);
          }  
        } else if (is.character(Myargs[2]) && Myargs[2] %in%
          AnyNameFunction) {
          MandatoryOnFunction <- Myargs[2];
          ForceFunction <- MandatoryOnFunction;
          MandatoryOnFunction <- Myargs[2];
          print(paste("Setting MandatoryOnFunction to ", Myargs[2], sep=""));
          flush.console();
          if (length(Myargs) == 2) {
            Myargs <- c(argeese, 1,1,1,1,1,1,1);
          } else if (length(Myargs) >= 3) {
            Myargs <- c(Myargs[1], Myargs[3:length(Myargs)],1,1,1,1,1);
          }
        } else {
          MandatoryOnFunction <- NULL;  
        }
        eval(parse(text=SetGText("MandatoryOnFunction", "globalenv()", S=1)));
        try(print(paste("We will be setting argeese for paper file is ",
            argeese, sep = "")))
    } else if (length(list.files("c:/Stat")) > 0) {
        argeese <- paste(DropHome, "//TwoLassoPackage//ParTabPaper.R",
            sep = "")
        Myargs = c(argeese, 1, 1, 1, 1, 1)
    } else if (length(list.files("~/DownloadPackages")) > 0) {
        argeese <- "~/DownloadPackages/ParTabPaper.R"
        Myargs = c(argeese, 1, 1, 1, 1, 1)
    } else if (length(list.files("~/MyPackages")) > 0) {
        argeese <- "~/SpikeTest/LongTestn100p1000.R"
        Myargs = c(argeese, 1, 1, 1, 1, 1)
    } else {
        print("Please Set RunGetArgs Better for this Server!")
        flush.console()
    }
    if (is.na(argeese)) {
        print("Please Set RunGetArgs Better for this Server, argeese is NA!")
        flush.console()
    }
    eval(parse(text = SetGText("Myargs", S = 1)))
    eval(parse(text = SetGText("argeese", S = 1)))
    MyID = suppressWarnings(try(as.numeric(Myargs[2])));
    
    
    
    if (is.na(MyID)) { MyID <- 1;}
    eval(parse(text=SetGText("MyID", S=1)));
}


ClearOutLargeAndSmall <- function() {
  setwd(ASmallContainDir);
  ASS <- unlist(list.files());
  for (kk in 1:length(ASS)) {
    ABB <- list.files(ASS[kk])
    if (length(ABB) >= 1) {
      for (kkt in 1:length(ABB)) {
        try(unlink(paste(ASS[kk], "//", ABB[kkt], sep=""), recursive = TRUE));
      }
    }
    try(unlink(ASS[kk], recursive=TRUE));
  }

  setwd(ALargeContainDir);
  ASS <- unlist(list.files());
  if (any(ASS == "SMS")) {
     setwd("SMS");
     ASS <- unlist(list.files());
     if (length(ASS) >= 1) {
     for (kk in 1:length(ASS)) {
        try(setwd(ALargeContainDir));
        try(setwd("SMS"));
        try(setwd(ASS[kk]));
        BB <- unlist(list.files());
        if (length(BB) >= 1) {
          for (tt in 1:length(BB)) {
            try(unlink(BB, recursive=TRUE));
          }
        }
     }
     try(setwd(ALargeContainDir))
     try(setwd("SMS"));
     try(unlink(ASS[kk], recursive=TRUE));
    }
    try(setwd(ALargeContainDir));
    try(unlink("SMS", recursive=TRUE));
  }
  setwd(OriginalOldwd)
}


RunSeekScriptPreamble <- function() {
return("
    TwoSimR5:::ScriptGetArgs();
CheckUniqueProcessIdentifier();    ## In 2LassoSaveRoutines.$

try(library(TwoLassoCpp, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);

##TwoSimR5:::ScriptGetArgs();

verbose = 1;
Practice = 1;
source(paste(as.character(argeese), sep=\"\"));

UpdateNameFunctions();

##MyID <- as.numeric(Myargs[1]);
TwoSimR5:::MakeDirSaves();

ABStart <- proc.time();
AStart <- proc.time();
ALTFiles <- NULL;
try(ALTFiles <- unlist(list.files(paste(ASmallContainDir, \"//Submissions\", sep=\"\"))));
ALT2Files <- NULL;
try(ALT2Files <- unlist(list.files(paste(ALargeContainDir, \"//Submissions\", sep=\"\"))));
TotalSlaves <- length(ALTFiles);
TotalSlaves2 <- length(ALT2Files);
WantFile <- paste(\"Request\", ZZMyID, \".RData\", sep=\"\");
if (WantFile %in% ALTFiles) {
  setwd(ASmallContainDir);
  try(load(paste(\"Submissions//\", WantFile, sep=\"\")));
  setwd(OriginalOldwd);
} else if (WantFile %in% ALT2Files) {
  setwd(ALargeContainDir);
  try(load(paste(\"Submissions//\", WantFile, sep=\"\")));
  setwd(OriginalOldwd);
} else {
  print(paste(\"Oh No No RData file for MyID = \", MyID, sep=\"\"));
  q();
}


if (length(Myargs) >= 2 && is.character(Myargs)) {
  AUT <- Myargs[2];
  if (AUT %in% SubMitList[,2]) {
    ForceABB <-  min((1:NROW(SubMitList))[SubMitList[,2] == AUT]);
    ForceFunction <- AUT;
  } else if (paste(\"Generate\", AUT, sep=\"\") %in% SubMitList[,2]) {
    ForceABB <-  min((1:NROW(SubMitList))[SubMitList[,2] == paste(\"Generate\", AUT, sep=\"\")]);    
    ForceFunction <- SubMitList[ForceABB,2];
  } else if (AUT %in% paste(\"Generate\", SubMitList[,2], sep=\"\")) {
    ForceABB <-  min((1:NROW(SubMitList))[AUT == paste(\"Generate\", SubMitList[,2], sep=\"\")]);    
    ForceFunction <- SubMitList[ForceABB,2];
  } else {
    ForceFunction <- NULL;  ForceABB <- -1;
  }
} else {
  ForceFunction <- NULL;
  ForceABB <- -1;
}

if (any(SubMitList[,2] == \"GroupBayesSpike\")) {
  SubMitList[SubMitList[,2] == \"GroupBayesSpike\",2] <- \"GenerateGroupBayesSpike\";
}
if (any(SubMitList[,2] == \"GroupTwoLassoSpread\")) {
  SubMitList[SubMitList[,2] == \"GroupTwoLassoSpread\",2] <- \"GenerateGroupTwoLassoSpread\";
}
if (!any(SubMitList[,3] %in% c(\"0\", 0, \" 0\"))) {
    AFilePrint(\"Hey, we're all done! \");
    flush.console();
    return(-1);
}
if (!exists(\"MandatoryOnFunction\")) {
  print(paste(\"Woah: MandatoryOnFunction does not exist!\")); flush.console();
} else if (is.null(MandatoryOnFunction)) {
  print(paste(\"MandatoryOnFunction Is Null\", sep=\"\")); flush.console();
} else {
  print(paste(\"MandatoryOnFunction Is: \", MandatoryOnFunction, sep=\"\"));
}
ABB <- min((1:NROW(SubMitList))[SubMitList[,3] %in% c(\"0\", 0, \" 0\")])

if (exists(\"ForceABB\") && exists(\"ForceFunction\") 
  && !is.null(ForceABB) && is.numeric(ForceABB) && ForceABB > 0) {
  ABB <- ForceABB;  MandatoryOnFunction <- ForceFunction;
}



SMSName <- SubMitList[ABB,1];
FunctionName  <- SubMitList[ABB,2];
eval(parse(text=GetG0Text(\"UniqueProcessIdentifier\"))); 
  if (is.character(BSSaveDir) && (nchar(BSSaveDir) < nchar(UniqueProcessIdentifier) ||
    substr(BSSaveDir, nchar(BSSaveDir) - nchar(UniqueProcessIdentifier)+1,
    nchar(BSSaveDir)) != UniqueProcessIdentifier)) {
    BSSaveDir <- paste(BSSaveDir, \"//\", UniqueProcessIdentifier, sep=\"\");      
  }
eval(parse(text=SetGText(\"UniqueProcessIdentifier\")));
");
}


StartBayesSpikeErrorPrior <- function() {
  return("
    library(BayesSpike);
  kLen = SMS$kLen;
  ##AFilePrint(\"ENetFixed Run\"); flush.console();
  T1 <- proc.time();
  if (!exists(\"NumChains\")) { NumChains = 3; }
  if (!exists(\"CountSamp\")) { CountSamp = 1000; }
  if (!exists(\"PriorStrength\")) { PriorStrength = 0; }
  if (!exists(\"tauSqA\")) { tauSqA = 1; }
  if (!exists(\"CutOff\")) { CutOff = .1; }
  
  
  if (kLen >= 1100 && PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puse + 2, PriorStrength*(1.0-SMS$puse) +
      round(kLen^.75));
    SigmaPrior = c(PriorStrength, SMS$SigmaSqNoise);
  } else if (SMS$kLen >= 1000 && PriorStrength <= 0) {
    PiAPrior = c(2, SMS$kLen);
    SigmaPrior = c(SMS$n, SMS$SigmaSqNoise);
  } else if (PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puse, PriorStrength * (1.0-SMS$puse));
    SigmaPrior = c(PriorStrength, SMS$SigmaSqNoise);
  } else {
    SigmaPrior = c(SMS$n,SMS$SigmaSqNoise); PiAPrior = c(2,SMS$kLen);
  }
  eval(parse(text=GetG0Text(\"BSSaveDir\")));
  eval(parse(text=GetG0Text(\"UniqueProcessIdentifier\")));
  if (is.character(\"BSSaveDir\") && substr(BSSaveDir, nchar(BSSaveDir)-nchar(UniqueProcessIdentifier)+1,
    nchar(BSSaveDir)) != UniqueProcessIdentifier) {
    BSSaveDir <- paste(BSSaveDir, \"//\", UniqueProcessIdentifier, sep=\"\");
    dir.create(BSSaveDir, showWarnings=FALSE, recursive=TRUE);
  }
  ")
}
StartTestToTryNuke <- function() {
  return("
    verbose = verbose; DoCI = TRUE;

  eval(parse(text=InitiateEvaluations()));
  eval(parse(text=TestForVerbose()));
  eval(parse(text=LoadSMSNow()));  
  eval(parse(text=FunctionUnderstandProcesses()));
  eval(parse(text=EvaluateProgress())); 
  eval(parse(text=EstablishFunctionAndIdentifiers()));
  
  RTI = unlist(UseFunctions[IFunction]);
  eval(parse(text=SetGText(\"SMS\", S=1)));
    ##AFilePrint(paste(\" Doing Analysis NameFunctions[\", ii, \"] = \",
    ##   NameFunctions[ii], sep=\"\")); 
  if (p > 1000) {
    AFilePrint(paste(\"Doing Analysis NameFunctions[\", IFunction, \"] = \",
        OnNameFunctions[IFunction], sep=\"\"));  flush.console();
  }
  limlas = NULL;
  
  ");
}

BayesSpikeBeginnerTestHelpFunction <- function() {
    return("
    
   library(BayesSpike);
  kLen = SMS$kLen;
  ##AFilePrint(\"ENetFixed Run\"); flush.console();
  T1 <- proc.time();
  if (!exists(\"NumChains\")) { NumChains = 3; }
  if (!exists(\"CountSamp\")) { CountSamp = 1000; }
  if (!exists(\"PriorStrength\")) { PriorStrength = 0; }
  if (!exists(\"tauSqA\")) { tauSqA = 1; }
  if (!exists(\"CutOff\")) { CutOff = .1; }
  
  if (NCOL(SMS$X) >= 1100) {
    SigmaSq <- .5 * var(SMS$Y);
  } else {
    SigmaSq <- SMS$SigmaSqNoise;
  }
  if (NCOL(SMS$X) >= 1100 && PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puse + 2, PriorStrength*(1.0-SMS$puse) +
      round(length(SMS$EndGroupIndices)));
    SigmaPrior <- c(PriorStrength, SMS$SigmaSqNoise);
  } else if (NCOL(SMS$X) >= 1000 && PriorStrength <= 0) {
    PiAPrior = c(2,
      round(length(SMS$EndGroupIndices)+1)); 
    SigmaPrior <- c(NROW(SMS$X),SMS$SigmaSqNoise);   
  } else if (PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puseGroups, PriorStrength * (1.0-SMS$puseGroups));
    SigmaPrior = c(PriorStrength, SMS$SigmaSqNoise);
  } else {
    SigmaPrior = c(1,1); PiAPrior = c(-1,-1);
  }

   if (NCOL(SMS$X) >= 10000) {
      MaxGibbsIters <- 200;
    } else {
      MaxGibbsIters <- 500;
    }

    X=SMS$X; Y=SMS$Y; 
      tauASq = tauSqA; PiAStart=SMS$puse; MaxGibbsIters = MaxGibbsIters;
      NumChains<-NumChains;
      eval(parse(text=GetG0Text(\"Verbose\")));
      eval(parse(text=GetG0Text(\"Verbose\", env=\"globalenv()\", S=1)))
      eval(parse(text=GetG0Text(\"Verbose\", env=\"TWOSIMENVIRONMENT\", S=1)));
      SigmaSq = as.numeric(SigmaSq); 
      try(Verbose <- 0); 
      PiAPrior = PiAPrior;
      SigmaPrior = SigmaPrior; DoSave=TRUE; DoRecord = c(0,0,0,0,0,0,0);
      SaveDir = BSSaveDir; IndexFirstRandomEffect = SMS$FirstGroupIndex;
      tauEndList = SMS$EndGroupIndices; ZeroOutBeforeSample = TRUE;
      StartRunProbVector = 10
      
    ")   
}

SecondStartBayesSpike <- function() {
return("
    eval(parse(text=EvalReDeclarers()));
eval(parse(text=EvalXandY()));
eval(parse(text=EvalTemperatureEE()));
eval(parse(text=CheckXandY()));
eval(parse(text=CheckSaveDir()));
eval(parse(text=CreateTBSR5()));
eval(parse(text=SetupToCreateMBS()));

eval(parse(text=CreateAndLockInMBS()));  ## Real Generation of MBS
eval(parse(text=SetupOnPiAAndHowSample()));
eval(parse(text=SetupDoTimePartsList()));
eval(parse(text=SetupTauEndList()));
eval(parse(text=SetupAlternateNoise()));
eval(parse(text=SetupInsertPointers()));
eval(parse(text=BSSetupEigenList()));

eval(parse(text=BayesSpikeSetupMBSOnSigmaAndStuff()));
eval(parse(text=BayesSpikeMergeAndPriors()));
eval(parse(text=BayesSpikeSetupTaus()));
eval(parse(text=BayesSpikeChecksForMemory()));

eval(parse(text=BayesSpikeSetupProbs()));
eval(parse(text=BayesSpikeUpToEarlyEndStep()));
eval(parse(text=BayesSpikeSetupDependencies()));
eval(parse(text=BayesSpikeSetupItersAndSave()));
eval(parse(text=BayesSpikeSetupSaveFileToBeDestroyed()));
eval(parse(text=BayesSpikeCodaListSetter()));
eval(parse(text=BayesSpikeSetupSaveToWrite()));

eval(parse(text=BayesSpikeSetupPriorProbFixedTau()));
eval(parse(text=BayesSpikeSetupTemperatureLists()));
TTii = 1
   eval(parse(text=BayesSpikeStartSetupOnIterACoda()));
   eval(parse(text=BayesSpikeCodaListSetAndCheck()));
   eval(parse(text=BayesSpikeSetupEarlyEndttNow()));
iti = 1;
   eval(parse(text=BayesSpikeOnChainIterInitiate()));
   try(MBS$tt <- 0);
");

    
}

SecondZeroOutHelpTest <- function() {
  return("
    if (length(MBS$tauEndList) >= 1) {
    if (!is.null(MyTryProbTau)) {
      MyProb <- MyTryProbTau;
    } else {
      MBS$SampleTausOnly <- 1;
      MBS$SampleNewTaus();
      MyProb <- MBS$ProbTau;
    }
    if (MBS$FirstRandom >= MBS$p) {
      print(paste(\"$$ ZeroOutFunction: ERROR ERROR ERROR, FirstRandom = \", MBS$FirstRandom, sep=\"\"));
      flush.console();
      print(\"ERRORERRORERRORERRORERRORERRORERROR\"); flush.console();
      MyT = \"1 = 2\";
      eval(parse(text=MyT));
      return(-1);
    }
    GoTau  <- MBS$OnTau;
    print(paste(\"$$ ZeroOutFunction: We sample Taus and the number of nonzero are \", length(GoTau[GoTau>0]),
     \"/\", length(GoTau), sep=\"\")); flush.console();
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
      print(paste(\"$$ after calculating Expt we have its sum is \", round(sum(Expt),4),
        \"/\", length(Expt), sep=\"\")); flush.console();
      Keepers <- sample(1:length(GoTau), prob=Expt, replace=FALSE, size=KeepCount);
      print(paste(\"$$ Now we reduce active set to length(Keepers) = \", length(Keepers), sep=\"\"));
      flush.console();
      DoTau <- rep(0, length(GoTau));
      DoTau[Keepers] <- GoTau[Keepers];
      try(MBS$DoAddCoordsOnSetBeta <- 0);
      try(MBS$BlankAllNewCoords());
      print(paste(\"AssignTau: about to Assign \")); flush.console();
      try(MBS$AssignTau(DoTau));
      print(paste(\"$$ Reduced DoTau length(Keepers) = \", length(Keepers), \" is added\", sep=\"\"));
      print(paste(\"$$ After AssignTau, AllNewCoords = \", MBS$AllNewCoords, sep=\"\")); flush.console();
      flush.console();
    } else {
      print(paste(\"$$ KeepCount = \", KeepCount, \" but length(GoTau) = \",
        length(GoTau[GoTau>0]), \"/\", length(GoTau), \" so no change.\", sep=\"\"));
      flush.console(); 
      try(MBS$BlankAllNewCoords());
      print(paste(\"$$ Now Assigning Tau to GoTau\")); flush.console();
      MBS$AssignTau(GoTau);
      print(paste(\"$$ After AssigningTau to GoTau length \", length(GoTau), \" we \",
        \" have AllNewCoords = \", MBS$AllNewCoords, sep=\"\"));
      flush.console();
    }
  }
  print(\"$$ Well Looking into the length of MBS$tauEndList\"); flush.console();
  if (length(MBS$tauEndList) >= 1) {
    AllRandomCoords <- MBS$AllNewCoords;
  } else {
    AllRandomCoords <- NULL;
  }
  if((MBS$FirstRandom > 1 && MBS$FirstRandom <= MBS$p) || length(MBS$tauEndList) <= 0) {
      if (length(MBS$tauEndList) == 0) {   pF = MBS$p;  
      } else { pF <- MBS$FirstRandom-1;}
      if (!is.null(MyTryProbFixed)) {
        BProbs <- MyTryProbFixed;
        if (Verbose >= 1) {
          print(paste(\"$$ ZeroOutFunction(): Probability of Historical fixed variables \",
            \"has Historical MIPS summing to \", round(sum(dLogit(BProbs)),4), \"/\", pF, sep=\"\"));
         flush.console();
        }
      } else {
        MBS$SampleTausOnly <- 1;
        print(paste(\"$$  ZeroOutFunction() : We're refreshing fixed effects. \")); flush.console();
        MBS$SampleFixedB();
        print(paste(\"$$ Now to SampleTausOnly\")); flush.console();
        MBS$SampleTausOnly <- 0;
        BProbs <- MBS$ProbFixed;
        if (Verbose >= 1) {
          print(paste(\"$$ ZeroOutFunction() : MIPs taken from MBS SampleFixedB()\", 
            round(sum(dLogit(BProbs)),4), \"/\", pF, sep=\"\"));
          flush.console();
        }
      }     
      
      OnPiAF <- MBS$OnPiA[1];
      KeepCount <- round(25*pF * OnPiAF);
      if (KeepCount <= 0) { KeepCount = 1; }
      if (KeepCount > pF) { KeepCount <- pF; }
      
      eBProbs <- GiveMeNewBalance(lX=BProbs/2.0, WantTot = KeepCount, StartMove=1.0, MaxTTs = 300);
      if (length(eBProbs[eBProbs > 0])  < KeepCount) {
        KeepCount <- length(eBProbs[eBProbs > 0]); 
      }                                                   
      if (Verbose >= 1) {
        print(\"$$ ZeroOutFunction() Sampling Keepers.\");
      }
      Keepers <- sample(1:pF, prob=eBProbs, replace=FALSE, size=KeepCount);
      DKeepers <- (1:pF)[!(( 1:pF)  %in% Keepers )]; 
      print(paste(\"$$  We have chosen to eliminate a count of Keepers \", length(DKeepers)));
      flush.console();
      if (MBS$FirstRandom <= 0 || is.null(MBS$tauEndList) || length(MBS$tauEndList) <= 0) {
        ABeta <- MBS$Beta;
        ABeta[MBS$Beta == 0.0 && MBS$BFixed >= 1] <- 1;
      } else if (MBS$FirstRandom > 1 && MBS$FirstRandom <= MBS$p) {
        ABeta <- c(rep(1, MBS$FirstRandom-1), rep(0, MBS$p- MBS$FirstRandom+1)) * MBS$Beta;
        ABeta[(1:length(MBS$BFixed))[MBS$Beta == 0.0 && MBS$BFixed >= 1]] <- 1;
      }
      
      if (length(DKeepers) > 0) {
        ABeta[DKeepers] <- 0.0;
      }
      try(MBS$DoAddCoordsOnSetBeta <- 0);
     
      try(MBS$RefreshBeta(ABeta));
      ##try(MBS$SetAllNewCoords(length(AllRandomCoords)+length(Keepers), length(Keepers)));
    }
    print(\"About to Count All New Coords() \"); flush.console();
    if (length(AllRandomCoords) + length(Keepers) >= 0) {
      try(MBS$CountAllNewCoords());
    }
  ")
}

FirstTwoLassoRegression <- function() {
  return("
   library(TwoLassoCpp);
    USENOISE <- SMS$SigmaSqNoise; kLen <- SMS$kLen;
    ##
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
	      LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) * sqrt( 2 * SDA / USENOISE)
        LambdaASeq = c( 1.1 * LambdaAS1, .05 * LambdaAS1);
        if (LambdaASeq[1] >= LambdaDSeq[2]) {
           LambdaASeq[1] = LambdaDSeq[1] + LambdaASeq[1];
           LambdaDSeq[1] = LambdaASeq[1] - LambdaDSeq[1];
           LambdaASeq[1] = LambdaASeq[1] - LambdaDSeq[1];
        }
        LambdaASeq = c(LambdaASeq[1], LambdaASeq);
        LambdaDSeq = c(LambdaDSeq[1], LambdaDSeq);        
	      OrderSeq = c(20, 4,1);
	    puseGroups <- SMS$puseGroups;
        ## Do this to debug TwoLassoCpp to look for write error!
        puseGroups <- SMS$puseGroups;
        
        eval(parse(text=GetG0Text(\"Verbose\", \"globalenv()\", S=1)));
        X=SMS$XX; Y=SMS$YY; PiA = puseGroups; 
       SigmaSq = USENOISE * log(NCOL(SMS$XX)); LambdaAK = LambdaASeq;
       LambdaDK = LambdaDSeq; OrderSeq = OrderSeq; 
       try(Verbose <- 2); 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5); RecordFlag = FALSE; HoldOn = FALSE; 
       TDFNu = SMS$tNoiseDF; SigmaPrior = c(USENOISE, 2);
       PiAPrior = c(puseGroups * sqrt(NCOL(SMS$XX)), (1-puseGroups) * sqrt(NCOL(SMS$XX)));
       IndexFirstRandomEffect = SMS$FirstGroupIndex; 
       tauEndList = SMS$EndGroupIndices
  ")
}

TextStartpiALasso <- function() {
return("

  if (!exists(\"SigmaEta\")) { SigmaEta = sigmaPriorData * SMS$n;}
  if (!exists(\"SigmaBarS\")) { SigmaBarS = SMS$SigmaSqNoise;}
  if (!exists(\"musize\")) { musize = piAPriorData * SMS$n;}
  if (!exists(\"R\")) { R = 0; }
  if (!exists(\"DoCI\")) { DoCI = DefaultDoCI; }
  library(TwoLassoCpp);
  USENOISE <- SMS$SigmaSqNoise;	  kLen <- SMS$kLen;  
	T1 <- proc.time();
	SDD <- sqrt(min(apply(SMS$XX,2,sd))^2 * (length(SMS$XX[,1])-1) );
	LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD, 
    sqrt(3.7 / USENOISE) * SDD *  max((2*length(SMS$XX[,1])),200));                          
	SABS <- min( abs(SMS$BetasReal[SMS$BetasReal != 0])  );
	SDA <- SDD / sqrt(length(SMS$XX[,1])-1);
	LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) * 
    sqrt( 2 * SDA^2 / USENOISE)
  LambdaASeq = c( 1.1 * LambdaAS1, .05 * LambdaAS1);
  if (LambdaASeq[1] >= LambdaDSeq[2]) {
    LambdaASeq[1] = LambdaDSeq[1] + LambdaASeq[1];
    LambdaDSeq[1] = LambdaASeq[1] - LambdaDSeq[1];
    LambdaASeq[1] = LambdaASeq[1] - LambdaDSeq[1];
  }
  LambdaASeq = c(LambdaASeq[1], LambdaASeq);
  LambdaDSeq = c(LambdaDSeq[1], LambdaDSeq);        
	OrderSeq = c(20, 4,1);
	if (R > 0 && R <= .5) {
	  NN = length(SMS$YY);
    SigmaMultVec = c( NN^R, 1, 1);                         
    LambdaASeq = c(LambdaASeq[1] * NN^(.5 * R), LambdaASeq[1], 
      LambdaASeq[length(LambdaASeq)]);
    LambdaDSeq = c(LambdaDSeq[1] * NN^(-.5 * R), LambdaDSeq[1],
      LambdaDSeq[length(LambdaDSeq)] );
  } else {
    SigmaMultVec = -1;
  }
  SigmaPrior <- c(SigmaBarS, SigmaEta);
  PiAPrior <- musize * c(SMS$puse, 1- SMS$puse);
 ");    
}
TextBayesSpikeMaster <- function() {
return("

X=SMS$X; Y=SMS$Y; 
      tauASq = tauSqA; PiAStart=SMS$puse;InitKKs=6; MaxGibbsIters = CountSamp;
      NumChains=NumChains;
      SigmaSq = SMS$SigmaSqNoise; Verbose = 0; tauEndList = NULL; PiAPrior = PiAPrior;
      SigmaPrior = SigmaPrior; DoSave=TRUE; DoSaveTBSR5=FALSE; DoRecord = c(0,0,0,0,0,0,0);
      SaveDir = BSSaveDir; IndexFirstRandomEffect = -1;
      StartRunProbVector = 5;
      DoLogitPreProb = DoLogitPreProb; DoLogitPostProb = DoLogitPostProb;
      dfRobit = dfRobit; PreRunMaxGibbsIters = PreRunMaxGibbsIters; 
      AlterWeightFlag = AlterWeightFlag
      
   
eval(parse(text=EvalReDeclarers()));
eval(parse(text=EvalXandY()));
eval(parse(text=EvalTemperatureEE()));
eval(parse(text=CheckXandY()));
eval(parse(text=CheckSaveDir()));
eval(parse(text=CreateTBSR5()));
eval(parse(text=SetupToCreateMBS()));

eval(parse(text=CreateAndLockInMBS()));  ## Real Generation of MBS
eval(parse(text=SetupOnPiAAndHowSample()));
eval(parse(text=SetupDoTimePartsList()));
eval(parse(text=SetupTauEndList()));
eval(parse(text=SetupAlternateNoise()));
eval(parse(text=SetupInsertPointers()));
eval(parse(text=BSSetupEigenList()));

eval(parse(text=BayesSpikeSetupMBSOnSigmaAndStuff()));
eval(parse(text=BayesSpikeMergeAndPriors()));
eval(parse(text=BayesSpikeSetupTaus()));
eval(parse(text=BayesSpikeChecksForMemory()));

eval(parse(text=BayesSpikeSetupProbs()));
eval(parse(text=BayesSpikeUpToEarlyEndStep()));
eval(parse(text=BayesSpikeSetupDependencies()));
eval(parse(text=BayesSpikeSetupItersAndSave()));
eval(parse(text=BayesSpikeSetupSaveFileToBeDestroyed()));
TempRunFlag <- 1;
if (Run == FALSE) {
  TTii = 1; if (!exists(\"ItTemps\")) { ItTemps = c(1); }
  eval(parse(text=BayesSpikeCodaListSetAndCheck())); 
  eval(parse(text=BayesSpikeCodaListSetter())); 
  ##return(MBS);
  if (is.null(TempRunFlag) || (is.numeric(TempRunFlag) && TempRunFlag < 0)) {
    print(\"Run Fail return MBS. \"); flush.console(); return(MBS);
  }
} else {
  eval(parse(text=BayesSpikeCodaListSetter()));
  if (is.null(TempRunFlag) || (is.numeric(TempRunFlag) && TempRunFlag < 0)) {
    print(\"Run Fail return GibbsSS. \"); flush.console(); return(MBS);
  }
  ##return(MBS);
}
eval(parse(text=BayesSpikeSetupSaveToWrite()));

eval(parse(text=BayesSpikeSetupPriorProbFixedTau()));
eval(parse(text=BayesSpikeSetupTemperatureLists()));
");    
}
TextStartBayesSpike <- function() {
  
return("
  eval(parse(text=SetGText(\"SMS\", S=1)));
    ##AFilePrint(paste(\" Doing Analysis NameFunctions[\", ii, \"] = \",
    ##   NameFunctions[ii], sep=\"\")); 
  if (p > 1000) {
    AFilePrint(paste(\"Doing Analysis NameFunctions[\", IFunction, \"] = \",
        OnNameFunctions[IFunction], sep=\"\"));  flush.console();
  }
  limlas = NULL;
     
  try(MyOldwd <- getwd());
  try(setwd(CrashDir));
  
  ##SMS; tauSqA = 1; CountSamp = 1000; NumChains = 3;
  PriorStrength = SMS$PriorStrength; SigmaPriorStrength = SMS$SigmaPriorStrength;
  CutOff = 0.1; DoCI = DefaultDoCI; DoMedian = TRUE; AutoInfo = \"Auto\";
 
  library(BayesSpike);
  kLen = SMS$kLen;
  ##AFilePrint(\"ENetFixed Run\"); flush.console();
  T1 <- proc.time();
  if (!exists(\"NumChains\")) { NumChains = 3; }
  if (!exists(\"CountSamp\")) { CountSamp = 1000; }
  if (!exists(\"PriorStrength\")) { PriorStrength = 0; }
  if (!exists(\"tauSqA\")) { tauSqA = 1; }
  if (!exists(\"CutOff\")) { CutOff = .1; }
  if (!exists(\"AutoInfo\")) { AutoInfo = \"Auto\"; }
  if (!exists(\"DoLogitPostPreProb\")) { DoLogitPostPreProb = 1; }
  if (SMS$LogitRegression == FALSE)  {
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
  
  if (is.logical(AutoInfo) && AutoInfo == TRUE) {
    AutoInfo = \"Auto\"
  } else if (is.logical(AutoInfo) && AutoInfo == FALSE) {
    AutoInfo = \"Info\"
  } else if (is.character(AutoInfo) && !(AutoInfo %in% c(\"Auto\", \"Info\"))) {
    print(paste(\"Error AUTOINFO NOT SET!\", sep=\"\"));
  } else if (is.character(AutoInfo) && AutoInfo %in% c(\"Auto\", \"Info\")) {
  } else {
    AutoInfo = \"Auto\";
  }

  if (AutoInfo == \"Info\") {  
  if (kLen >= 1100 && PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puse + 2, PriorStrength*(1.0-SMS$puse) +
      round(kLen^.75));
    SigmaPrior = c(SigmaPriorStrength, SMS$SigmaSqNoise);
  } else if (SMS$kLen >= 1000 && PriorStrength <= 0) {
    PiAPrior = c(2, SMS$kLen);
    SigmaPrior = c(SMS$n, SMS$SigmaSqNoise);
  } else if (PriorStrength > 0) {
    PiAPrior = c(PriorStrength * SMS$puse, PriorStrength * (1.0-SMS$puse));
    SigmaPrior = c(SigmaPriorStrength, SMS$SigmaSqNoise);
  } else {
    SigmaPrior = c(SMS$n,SMS$SigmaSqNoise); PiAPrior = c(2,SMS$kLen);
  }

  } else {
    if (kLen >= 1100 && PriorStrength > 0) {
    PiAPrior = c(PriorStrength * 1/SMS$kLen + 1, PriorStrength +
      1);
    SigmaPrior = c(SigmaPriorStrength, .2*var(SMS$Y));
  } else if (SMS$kLen >= 1000 && PriorStrength <= 0) {
    PiAPrior = c(2, SMS$kLen+1);
    SigmaPrior = c(SMS$n, .2*var(SMS$Y));
  } else if (PriorStrength > 0) {
    PiAPrior = c(PriorStrength * 1 / SMS$kLen+1, PriorStrength +1);
    SigmaPrior = c(SigmaPriorStrength, SMS$SigmaSqNoise);
  } else {
    SigmaPrior = c(-1,-1); PiAPrior = c(-1,-1);
  }
  }
  eval(parse(text=GetG0Text(\"BSSaveDir\")));
  eval(parse(text=GetG0Text(\"UniqueProcessIdentifier\")));
  if (is.character(BSSaveDir) && (nchar(BSSaveDir) < nchar(UniqueProcessIdentifier) ||
    substr(BSSaveDir, nchar(BSSaveDir) - nchar(UniqueProcessIdentifier)+1,
    nchar(BSSaveDir)) != UniqueProcessIdentifier)) {
    BSSaveDir <- paste(BSSaveDir, \"//\", UniqueProcessIdentifier, sep=\"\");      
  }
");
    
}