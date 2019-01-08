DefaultNumDetails = 9;  ANum = 0.0;
.onLoad <- function(libname, pkgname) { 
ANum <- 0.0;
try(eval(parse(text=GetG0Text("ANum", "globalenv()", S=1))));
try(eval(parse(text=SetGText("ANum", "globalenv()", S=1))));
try(library(methods, warn.conflicts = FALSE, quietly = TRUE), silent=TRUE);
try(library(R.methodsS3, warn.conflicts = FALSE, quietly = TRUE), silent=TRUE);
print("Loading R.oo"); flush.console();
  try(library(R.oo, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
print("Loading TwoLasso"); flush.console();
  try(library(TwoLasso, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
print("Loading TwoLassoCpp"); flush.console();
  try(library(TwoLassoCpp, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
print("Loading BayesSpike"); flush.console();
  try(library(Rcpp, warn.conflicts=FALSE, quietly=TRUE));
  ##try(library(BayesSpike, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
  print("Welcome to TwoSim OnLoad "); flush.console();
CovRho = 0;  eval(parse(text=SetGText("CovRho", envir = "globalenv()",S=1)));
BM = 0;  eval(parse(text=SetGText("BM", S=1)));

eval(parse(text=GetG0Text("TWOSIMENVIRONMENT", S=1)));
TWOSIMENVIRONMENT = environment();
eval(parse(text=SetGText("TWOSIMENVIRONMENT", 
  envir = "TWOSIMENVIRONMENT", S=1)));
eval(parse(text=SetGText("TWOSIMENVIRONMENT", 
  envir = "globalenv()", S=1)));
eval(parse(text=GetG0Text("RSIMNAMESPACE", S=1)));
if (is.null(RSIMNAMESPACE) || is.numeric(RSIMNAMESPACE)  ||
  !is.environment(RSIMNAMESPACE) ) {
  RSIMNAMESPACE = environment(); 
}
try(eval(parse(text=SetGText("TWOSIMENVIRONMENT", 
  envir = "RSIMNAMESPACE", S=1))), silent=TRUE);
try(eval(parse(text=SetGText("RSIMNAMESPACE", 
  envir = "RSIMNAMESPACE", S=1))), silent=TRUE); 
try(eval(parse(text=SetGText("RSIMNAMESPACE", 
  envir = "TWOSIMENVIRONMENT", S=1))), silent=TRUE);   
eval(parse(text=SetGText("RSIMNAMESPACE", 
  envir = "globalenv()", S=1)));  
    
  
##argies = 0;

eval(parse(text=GetG0Text("OnSimType", envir="globalenv()", S=1)));
if (is.null(OnSimType) || !is.character("OnSimType")) {
  OnSimType <- "Default";
}
eval(parse(text=SetGText("OnSimType", "globalenv()", S=1)));
ExperimentName = "DefaultE";
eval(parse(text=SetGText("ExperimentName", envir = "globalenv()",S=1)));

UniqueProcessIdentifier = NULL;
eval(parse(text=SetGText("UniqueProcessIdentifier", envir = "globalenv()",S=1)));

DefaultGenerateBetaVec  = GenerateBetaVec2;  BM = 2;

eval(parse(text=SetGText("DefaultGenerateBetaVec",envir = "globalenv()", S=1)));

DefaultCorrelationXmatrix <- function(kLen) {
  CorrelationXmatrix1(kLen, CovRho);
}
eval(parse(text=SetGText("DefaultCorrelationXmatrix",S=1)));
eval(parse(text=SetGText("DefaultGenerateBetaVec",S=1)));
ADDDIRCONT = "DefaultADDDIRCONT"
eval(parse(text=SetGText("ADDDIRCONT",envir = "globalenv()",S=1)));

eval(parse(text=GetG0Text("piSV", envir="globalenv()", S=1)));
piSV = 0; eval(parse(text=SetGText("piSV", envir = "globalenv()",S=1)));
eval(parse(text=SetGText("piSV", S=1)));

eval(parse(text=GetG0Text("piAPriorData", envir="globalenv()", S=1)));
piAPriorData = 0; eval(parse(text=SetGText("piAPriorData", envir = "globalenv()",S=1)));
eval(parse(text=SetGText("piAPriorData", S=1)));
eval(parse(text=GetG0Text("jobii", envir="globalenv()", S=1)));
jobii = 1;
eval(parse(text=SetGText("jobii", S=1)));

print("IdenfityAllDirFunctions");flush.console();
IdentifyAllDirFunction();


print("DeclaringAllTestFunctions");flush.console();
DeclareAllTestFunctions();


print("SetAllComparisonRCodeFunctions");flush.console();
SetAllComparisonRCodeFunctions();

print("SetDefaultOtherParameters"); flush.console();
SetDefaultOtherParameters();
 
print("SetDefaultExperimentName"); flush.console(); 
ExperimentName = "DefaultE";
eval(parse(text=SetGText("ExperimentName",envir = "globalenv()",S=1)));

quoteMore = "DefaultQuoteMore"
eval(parse(text=SetGText("quoteMore", envir = "globalenv()",S=1)));


print("SetEfficientSumulatorDefaults"); flush.console();
SetEfficientSimulatorDefaults();      ## In 0007AListOfTwoLassoDefaultParameters.R
print("SetDefaultSimulationFunctions"); flush.console();
GenerateDefaultSimulationFunctions(); ## in 0006GenerateDefaultSimulationFunctions.R

TargetTotalSimulations = 1000;
eval(parse(text=SetGText("TargetTotalSimulations",envir = "globalenv()",S=1)));

DefaultForceFunction=NULL;
eval(parse(text=SetGText("DefaultForceFunction",envir = "globalenv()",S=1)));

TargetTotalSimulations = 400;
eval(parse(text=SetGText("TargetTotalSimulations",envir = "globalenv()",S=1)));
eval(parse(text=GetG0Text("TooManyFails", envir="globalenv()", S=1)));
if (is.null(TooManyFails) || TooManyFails <= 0)  {
TooManyFails = TargetTotalSimulations+1;
}
eval(parse(text=SetGText("TooManyFails", envir = "globalenv()",S=1)));

eval(parse(text=GetG0Text("TooManyFailsNowKillFunction", envir="globalenv()", S=1)));
if (is.null(TooManyFailsNowKillFunction) || TooManyFailsNowKillFunction <= 0)  {
TooManyFailsNowKillFunction = 20;
}
eval(parse(text=SetGText("TooManyFailsNowKillFunction", envir = "globalenv()",S=1)));

TooManyFailsPerOneSimFunction = 4;
eval(parse(text=SetGText("TooManyFailsPerOneSimFunction", envir = "globalenv()",S=1)));

piSV = 0; eval(parse(text=SetGText("piSV", envir = "globalenv()",S=1)));
sigSV = 0; eval(parse(text=SetGText("sigSV", envir = "globalenv()",S=1)));
LARSSEEKMETHOD=1; eval(parse(text=SetGText("LARSSEEKMETHOD", S=1)));
LARSSeekFlag = 1; eval(parse(text=SetGText("LARSSeekFlag", S=1))); 
piAPriorData=0; eval(parse(text=SetGText("piAPriorData", S=1)));
sigmaPriorData=0; eval(parse(text=SetGText("sigmaPriorData", S=1)));
DefaultNoiseDF=0; eval(parse(text=SetGText("DefaultNoiseDF", S=1)));

TableOutName="OutTable.csv"; eval(parse(text=SetGText("TableOutName", S=1)));

DefaultNumDetails = 9; eval(parse(text=SetGText("DefaultNumDetails", S=1)));

TwoSimLoadIn();

eval(parse(text=GetG0Text("MyStartProcTime", "globalenv()", S=1)));
MyStartProcTime <- proc.time();
eval(parse(text=SetGText("MyStartProcTime", "globalenv()", S=1)));


eval(parse(text=GetG0Text("MaxRunProcTime", "globalenv()", S=1)));
if (MaxRunProcTime == 0) {
MaxRunProcTime <- 60 * 60 * .75
}
eval(parse(text=SetGText("MaxRunProcTime", "globalenv()", S=1)));

verbose = FALSE; eval(parse(text=SetGText("verbose", S=1)));

eval(parse(text=GetG0Text("UpdateAllFunctions", envir="globalenv()", S=1)));
UpdateAllFunctions <- function() {
  if (!exists("CovRho", globalenv())) {
    print("UpdateAllFunctions: Warning CovRho not in globalenv!"); flush.console();
    TwoSim:::.OnLoad("TwoSim", "TwoSim"); 
  }
  eval(parse(text=GetG0Text("NameFunctions", S=1)));
  SubSetNameFunctions = NameFunctions;
  DeclareAllTestFunctions(SubSetNameFunctions=SubSetNameFunctions);
  eval(parse(text=GetG0Text("piSV", envir="globalenv()", S=1)));
  eval(parse(text=GetG0Text("piSV", S=1)));
  eval(parse(text=GetG0Text("sigSV", S=1)));
  eval(parse(text=SetGText("piSV", S=1)));
  eval(parse(text=SetGText("sigSV", S=1)));  
  eval(parse(text=GetG0Text("LARSSEEKMETHOD", S=1)));
  eval(parse(text=GetG0Text("LARSSeekFlag", S=1)));
  eval(parse(text=GetG0Text("piAPriorData", S=1)));
  eval(parse(text=GetG0Text("sigmaPriorData", S=1)));
  eval(parse(text=GetG0Text("DefaultNoiseDF",S=1)));
  eval(parse(text=SetGText("piAPriorData", S=1)));
  eval(parse(text=SetGText("sigmaPriorData", S=1)));
  eval(parse(text=SetGText("DefaultNoiseDF",S=1)));
  
  ATryText <- "try(library(AlanDirectories, 
    warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);  
  AlanDirectories:::SetSaveHomes();
  try(dir.create(SAVEHOME, showWarnings=FALSE));
  "
  try(eval(parse(text=ATryText)));
  eval(parse(text=GetG0Text("TwoSimSaveDir", envir="globalenv()", S=1)));
    TwoSimSaveDir <- paste(SAVEHOME, "//TwoSimSaveHome", sep="");
    try(dir.create(TwoSimSaveDir, showWarnings=FALSE));
    eval(parse(text=SetGText("TwoSimSaveDir", envir="globalenv()", S=1)));
    eval(parse(
      text=GetG0Text(
      "BayesSpikeTwoSimSaveDir", envir="globalenv()", S=1)));
    BayesSpikeTwoSimSaveDir <- 
      paste(TwoSimSaveDir, "//BayesSpikeSaves", sep="")
    try(dir.create(BayesSpikeTwoSimSaveDir, showWarnings=FALSE));
    eval(parse(
      text=SetGText("BayesSpikeTwoSimSaveDir", envir="globalenv()", S=1)));

if (!exists("Myargs", globalenv())) {  
  eval(parse(text=GetG0Text("Myargs", envir="globalenv()", S=1)));  
  Myargs <- commandArgs(TRUE)
  eval(parse(text=SetGText("Myargs", "globalenv()", S=1)))
} else {
  eval(parse(text=GetG0Text("Myargs", "globalenv()", S=1)))
}
argies <- Myargs;
if (length(argies)>= 8) {
  if (as.numeric(argies[8])  == 1) {
    DefaultGenerateBetaVec  = GenerateBetaVec1;  BM = 1;
  } else if (as.numeric(argies[8])  == 2) {
    DefaultGenerateBetaVec  = GenerateBetaVec2;  BM = 2;
  }  else if (as.numeric(argies[9]) == 3) {
    DefaultGenerateBetaVec  = GenerateBetaVec3;  BM = 3;
  }   else {
    DefaultGenerateBetaVec  = GenerateBetaVec2;  BM = 2;
  }
} else {
  DefaultGenerateBetaVec  = GenerateBetaVec2;  BM = 2;
}

eval(parse(text=GetG0Text("DefaultMaxLockTime", "globalenv()", S=1)));
if (is.null(DefaultMaxLockTime) || DefaultMaxLockTime <= 0) {
  if (length(list.files("c:/Stat")) >= 1) {
    DefaultMaxLockTime <- 5;    
  } else {
    DefaultMaxLockTime <- 100;    
  }
}
eval(parse(text=SetGText("DefaultMaxLockTime", "globalenv()", S=1)));

eval(parse(text=GetG0Text("OriginalOldwd", S=1)));
OriginalOldwd <- getwd();
eval(parse(text=SetGText("OriginalOldwd", S=1)));   


eval(parse(text=GetG0Text("MaxDuplicateTrials", "globalenv()", S=1)));
if (is.null(MaxDuplicateTrials) || MaxDuplicateTrials <= 0) {
  MaxDuplicateTrials <- 2;
}
eval(parse(text=SetGText("MaxDuplicateTrials", "globalenv()", S=1)));

eval(parse(text=GetG0Text("InitiateDuplicate", "globalenv()", S=1)));
if (is.null(InitiateDuplicate) || InitiateDuplicate <= 0) {
  InitiateDuplicate <- 24 * 60 * 60 /2;
}
eval(parse(text=SetGText("InitiateDuplicate", "globalenv()", S=1)));

eval(parse(text=GetG0Text("CallFailure", "globalenv()", S=1)));
if (is.null(CallFailure) || CallFailure <= 0) {
  CallFailure <- 24 * 60 * 60 * 2;
}
eval(parse(text=SetGText("CallFailure", "globalenv()", S=1)));

eval(parse(text=SetGText("DefaultGenerateBetaVec", "globalenv()", S=1)))
if (length(Myargs) >= 2) {
  try(eval(parse(text=GetG0Text("NInstance", "globalenv()", S=1))));
  try(NInstance <- as.numeric(Myargs[2]));
  if (is.null(NInstance) || is.na(NInstance)) { NInstance <- 666; }
} else {
  try(eval(parse(text=GetG0Text("NInstance"))));
}
try(eval(parse(text=SetGText("NInstance", "globalenv()", S=1))));
eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
if (!exists("UniqueProcessIdentifier") || is.null(UniqueProcessIdentifier) ||
  length(UniqueProcessIdentifier) != 1 || is.na(UniqueProcessIdentifier)) {
  UniqueProcessIdentifier <- DeclareUniqueProcessIdentifier(NInstance, TimeRun =-1);
  eval(parse(text=SetGText("UniqueProcessIdentifier", "globalenv()", S=1)));
}
if (length(Myargs) == 0) {
  if (length(list.files("c:/Stat")) > 0) {
    Myargs = c("c:/Stat/TwoLassoPackage/ParTabPaper.R", 1); 
  } else if (length(list.files("/n/scratch06/airoldi_scratch/alenarcic/")) > 0) {
    Myargs = c("/n/home13/alenarcic/TwoLassoTest/Palak/ParTabPaper.R", 1);
  } else if (length(list.files("/Users/lenarcic/Dropbox/")) > 0) { 
    Myargs = c("/Users/lenarcic/DropBox/TwoLassoPackage/ParTabPaper.R", 1);
  } else if (length(list.files("/Users/alenarc/Dropbox/")) > 0) {
    Myargs = c("/Users/alenarc/DropBox/TwoLassoPackage/ParTabPaper.R", 1);
  }
  eval(parse(text=SetGText("Myargs", envir="globalenv()", S=1)));
  eval(parse(text=GetG0Text("argeese", envir="globalenv()", S=1)));
  argeese = Myargs[1]; 
  eval(parse(text=SetGText("argeese", envir="globalenv()", S=1)));
  eval(parse(text=GetG0Text("NInstance", envir="globalenv()", S=1)));
  try(NInstance <- Myargs[2]);
  eval(parse(text=SetGText("NInstance", envir="globalenv()", S=1)));
  } else {
  eval(parse(text=GetG0Text("Myargs", envir="globalenv()", S=1)));
  eval(parse(text=GetG0Text("argeese", envir="globalenv()", S=1)));
  argeese = Myargs[1]; 
  eval(parse(text=SetGText("argeese", envir="globalenv()", S=1)));
  try(eval(parse(text=GetG0Text("NInstance", envir="globalenv()", S=1))));
  try(NInstance <- Myargs[2]);
  try(eval(parse(text=SetGText("NInstance", envir="globalenv()", S=1))));
  }
 
  
  ATryText <- "try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE),
    silent=TRUE);  
    AlanDirectories:::SetSaveHomes();
    try(dir.create(SAVEHOME, showWarnings=FALSE));
  "
  try(eval(parse(text=ATryText)));
  eval(parse(text=GetG0Text("TwoSimSaveDir", envir="globalenv()", S=1)));
    TwoSimSaveDir <- paste(SAVEHOME, "//TwoSimSaveHome", sep="");
    try(dir.create(TwoSimSaveDir, showWarnings=FALSE));
    eval(parse(text=SetGText("TwoSimSaveDir", envir="globalenv()", S=1)));
    eval(parse(
      text=GetG0Text(
      "BayesSpikeTwoSimSaveDir", envir="globalenv()", S=1)));
    BayesSpikeTwoSimSaveDir <- 
      paste(TwoSimSaveDir, "//BayesSpikeSaves", sep="")
    try(dir.create(BayesSpikeTwoSimSaveDir, showWarnings=FALSE));
    eval(parse(
      text=SetGText("BayesSpikeTwoSimSaveDir", envir="globalenv()", S=1)));
      
    eval(parse(text=GetG0Text("OurPrintFile", envir="globalenv()", S=1)));
    eval(parse(text=GetG0Text("OurErrorFile", envir="globalenv()", S=1)));
    eval(parse(text=GetG0Text("JustPrintFile", envir="globalenv()", S=1)));
    eval(parse(text=GetG0Text("DefaultTwoSimPrintDir", envir="globalenv()", S=1)));
    eval(parse(text=GetG0Text("DefaultTwoSimErrorDir", envir="globalenv()", S=1)));
    if (!is.null(DefaultTwoSimPrintDir) && is.character(DefaultTwoSimPrintDir) && DefaultTwoSimPrintDir != "") {
      if (!is.na(argeese) && is.character(argeese) && argeese[1] != "") {
        AT <- unlist(strsplit(argeese[1], "\\."))[1];
        if (substr(AT, 1, nchar("c:")) %in% c("c:", "C:") ||
          substr(AT, 1, nchar("/netscr")) %in% c("/netscr", "/lustre")) {
          APrintDir <- paste( AT, "//Print", sep="");
          AErrorDir <- paste( AT, "//Error", sep="");  
          APrintLockDir <- paste(APrintDir, "//Lock", sep="");      
        }  else {
          APrintDir <- paste(DefaultTwoSimPrintDir, "//", AT, sep="");
          AErrorDir <- paste(DefaultTwoSimErrorDir, "//", AT, sep="");
          APrintLockDir <- paste(APrintDir, "//Lock", sep=""); 
        }
        dir.create(APrintDir, showWarnings=FALSE, recursive=TRUE);
        dir.create(APrintLockDir, showWarnings=FALSE, recursive=TRUE);
        dir.create(AErrorDir, showWarnings=FALSE, recursive=TRUE);
        eval(parse(text=SetGText("APrintDir", "globalenv()", S=1)));
        eval(parse(text=SetGText("AErrorDir", "globalenv()", S=1)));
        eval(parse(text=SetGText("APrintLockDir", "globalenv()", S=1)));
      } else {
        APrintDir <- DefaultTwoSimPrintDir; dir.create(DefaultTwoSimPrintDir, showWarnings=FALSE, recursive=TRUE)
        APrintLockDir <- paste(APrintDir, "//Lock", sep="");  dir.create(APrintLockDir, showWarnings=FALSE, recursive=TRUE);
        AErrorDir <- DefaultTwoSimErrorDir; dir.create(DefaultTwoSimErrorDir, showWarnings=FALSE, recursive=TRUE);
        eval(parse(text=SetGText("APrintDir", "globalenv()", S=1)));
        eval(parse(text=SetGText("AErrorDir", "globalenv()", S=1)));
        eval(parse(text=SetGText("APrintLockDir", "globalenv()", S=1)));
      }
      ATime <- date();  ATime <- paste(unlist(strsplit(ATime, " ")), collapse="");
      ATime <- paste(unlist(strsplit(ATime, ":")), collapse="");
      PTime <- proc.time(); PTime <- paste(PTime[1:3], collapse="T");
      PTime <- paste(unlist(strsplit(PTime, "\\.")), collapse="p");
      R <- sample(1:100, size= 1);
      AFileName <- paste(ATime, "D", PTime, "R", R, ".txt", sep="")
      OurPrintFile = paste(APrintDir, "//", AFileName, sep="");
      OurErrorFile = paste(AErrorDir, "//", AFileName, sep="");
      JustPrintFile <- paste(UniqueProcessIdentifier, AFileName, sep="");
      eval(parse(text=SetGText("AFileName", "globalenv()", S=1)));
      eval(parse(text=SetGText("JustPrintFile", envir="globalenv()", S=1)));
      eval(parse(text=SetGText("OurPrintFile", "globalenv()", S=1)));
      eval(parse(text=SetGText("OurErrorFile", "globalenv()", S=1)));
      print(paste("Note: TwoSimAPrint we will save to ", OurPrintFile, sep=""));
      print(paste("Note APrintDir is ", APrintDir, sep="")); flush.console();
      print(paste("Note AErrorDir is ", AErrorDir, sep="")); flush.console();
      print(paste("Note DefaultTwoSimPrintDir is ", DefaultTwoSimPrintDir, sep="")); flush.console();
      print(paste("Note: Also DefaultTwoSimPrintDir: ", DefaultTwoSimPrintDir, sep=""));
      flush.console();
      flush.console();
      AF <- NULL;
      try(Oldwd <-getwd());
      try(setwd(APrintDir));
      try(AF <- file(JustPrintFile, "wt"));
      try(setwd(Oldwd));
      if (!is.null(AF)) {
        try(writeLines(text=paste("Hello: Starting TwoSimR5 Analysis stream: ", AFileName, sep=""),
          con=AF));
        try(writeLines(text=paste("JustPrintFile: ", JustPrintFile, sep="")));
        try(writeLines(text=paste("APrintDir: ", APrintDir, sep="")));
        try(close(AF));
      }
    } else {
      OurPrintFile = "";
      print("  All Printing will be done to the screen.  "); flush.console();
      eval(parse(text=SetGText("OurPrintFile", "globalenv()", S=1)));
    }
    MyPrintEventLocation <- 0;
    eval(parse(text=SetGText("MyPrintEventLocation", "globalenv()", S=1)));
    AFilePrint <- function(AText, LockIn="No") {
      eval(parse(text=GetG0Text("APrintDir", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("JustPrintFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("OurPrintFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("IWillStillNotPrint", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("MyPrintEventLocation", "globalenv()", S=1)));
      APTime <- proc.time()[3]
      if (is.null(MyPrintEventLocation) || !is.numeric(MyPrintEventLocation)) {
        try(MyPrintEventLocation <- 10000001);
      }
      try(MyPrintEventLocation <- MyPrintEventLocation +1);
      if (is.numeric(JustPrintFile) || is.null(JustPrintFile) || 
        (is.character(JustPrintFile) && JustPrintFile == "")) {
        try(print(paste(MyPrintEventLocation, ":", APTime, " : ", AText, sep=""))); 
        flush.console(); 
        return(1);
      }
      AF <- NULL;
      Oldwd <- getwd();
      ATT <- setwd(APrintDir)
      try(AF <- file(JustPrintFile, "at"));
      if (!is.null(AF)) {
        try(writeLines(text=
          paste(MyPrintEventLocation, ":", APTime, " : ", AText, sep=""), con=AF), silent=TRUE);
        try(close(AF));
      }
      if (is.character(LockIn) && LockIn %in% c("Lock", "UnLock", "L", "U", "SumLock", "SumUnlock",
        "SumUnLock")) {
        try(ATT <- setwd(APrintLockDir));
        try(AF <- NULL);
        try(AF <- file(JustPrintFile, "at"));
        if (!is.null(AF)) {
          try(writeLines(text=
            paste(MyPrintEventLocation, ":", APTime, " : ", LockIn, ": ", AText, sep=""), con=AF), silent=TRUE);
          try(close(AF));
        }
      }
      if (is.null(IWillStillNotPrint) ||  !exists("IWillStillNotPrint") || 
        (is.numeric(IWillStillNotPrint)  && IWillStillNotPrint[1] != 1) ||
        (is.logical(IWillStillNotPrint) && IWillStillNotPrint == FALSE)) {
         try(print(paste(MyPrintEventLocation, ":", APTime, " : ", AText, sep="")), silent=TRUE); flush.console();  
      }
      eval(parse(text=SetGText("MyPrintEventLocation", "globalenv()", S=1)));
      try(setwd(Oldwd));
      return(1);
    }
    eval(parse(text=GetG0Text("IWillStillNotPrint", "globalenv()", S=1)));
    if (is.null(IWillStillNotPrint)) { IWillStillNotPrint = 0; }
    eval(parse(text=SetGText("IWillStillNotPrint", "globalenv()", S=1)));
    try(eval(parse(text=SetGText("AFilePrint", "globalenv()", S=1))));
    try(eval(parse(text=SetGText("AFilePrint", "TWOSIMENVIRONMENT", S=1))));
    eval(parse(text=GetG0Text("AErrorPrint", "globalenv()", S=1)));
    AErrorPrint <- function(AText) {
      eval(parse(text=GetG0Text("OurErrorFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("OurPrintFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("JustPrintFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("APrintDir", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("AErrorDir", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("IWillStillNotPrint", "globalenv()", S=1)));
      eval(parse(text=SetGText("MyPrintEventLocation", "globalenv()", S=1)));
      STTText <- unlist(strsplit(AText, "\n"));
      if (is.numeric(JustPrintFile) || is.null(JustPrintFile) || 
        (is.character(JustPrintFile) && JustPrintFile == "")) {
         print("*** ERROR ERROR ERROR ERROR ERROR ERROR ERROR")
         print(paste("*** ERROR: OurPrintFile is ", OurPrintFile, sep=""));
         for (ii in 1:length(STTText)) {
           print(paste("*** ERROR : ", MyPrintEventLocation, " : ", STTText[ii], sep="")); 
           flush.console(); 
         }
         print("***"); flush.console();
         print("***************************************************************");
         flush.console();
         return(1);
      }
      AF <- NULL;
      Oldwd <- getwd();
      try(setwd(AErrorDir));
      try(AF <- file(OurErrorFile, "at"));
      try(setwd(Oldwd));
      try(writeLines(text="*** ERROR ERROR ERROR ERROR ERROR ERROR ERROR", con=AF));
      try(writeLines(text=paste("*** ERROR at MyPrintEventLocation = ", MyPrintEventLocation, sep="")));
      try(writeLines(text=paste("*** ERROR: the Our PrintFile is ", OurPrintFile, sep="")));
      for (ii in 1:length(STTText)) {
         try(writeLines(text=paste("*** ERROR:", STTText[ii], sep=""), con=AF));
      }
    
      try(close(AF));
      if (is.null(IWillStillNotPrint) ||  !exists("IWillStillNotPrint") || 
        (is.numeric(IWillStillNotPrint)  && IWillStillNotPrint[1] != 1) ||
        (is.logical(IWillStillNotPrint) && IWillStillNotPrint == FALSE)) {
        print("*** ERROR ERROR ERROR ERROR ERROR ERROR ERROR")
        print(paste("*** ERROR: OurPrintFile is ", OurPrintFile, sep=""));
        ##print(paste("*** ERROR : ", MyPrintEventLocation, AText, sep=""));         
        for (ii in 1:length(STTText)) {
           print(paste("*** ERROR : ", MyPrintEventLocation, " : ", STTText[ii], sep="")); 
           flush.console(); 
        }
      }
      return(1);
    }
    eval(parse(text=SetGText("AErrorPrint", "globalenv()", S=1)));

}
eval(parse(text=SetGText("UpdateAllFunctions", envir="globalenv()", S=1)));
if (!exists("Myargs", globalenv())) {  
  eval(parse(text=GetG0Text("Myargs", envir="globalenv()", S=1)));  
  Myargs <- commandArgs(TRUE)
  eval(parse(text=SetGText("Myargs", "globalenv()", S=1)))
} else {
  eval(parse(text=GetG0Text("Myargs", "globalenv()", S=1)))
}
if (length(Myargs) == 0) {
if (length(list.files("c:/Stat")) > 0) {
  Myargs = c("c:/Stat/TwoLassoPackage/ParTabPaper.R", 1); 
} else if (length(list.files("/n/scratch06/airoldi_scratch/alenarcic/")) > 0) {
  Myargs = c("/n/home13/alenarcic/TwoLassoTest/Palak/ParTabPaper.R", 1);
} else if (length(list.files("/Users/lenarcic/Dropbox/")) > 0) {
  Myargs = c("/Users/lenarcic/DropBox/TwoLassoPackage/ParTabPaper.R", 1);
} else if (length(list.files("/Users/alenarc/Dropbox/")) > 0) {
  Myargs = c("/Users/alenarc/DropBox/TwoLassoPackage/ParTabPaper.R", 1);
}
}
eval(parse(text=SetGText("Myargs", S=1)));
eval(parse(text=GetG0Text("argeese", S=1)));
argeese = Myargs[1]; 
eval(parse(text=SetGText("argeese", S=1)));
try(eval(parse(text=GetG0Text("NInstance", S=1))));
try(NInstance <- Myargs[2]);
try(eval(parse(text=SetGText("NInstance", "globalenv()", S=1))));

SetDefaultScore("NumDetails", DefaultNumDetails);
SetDefaultScore("NumPrint", 40);
SetDefaultScore("n_job", 40);





#################################################
## Now to Load ParTabPapLarge
## argeese= "c:/Stat/TwoLassoPackage/ParTabPaper.R"

eval(parse(text=GetG0Text("DefaultCorrelationXmatrix", S=1)));
eval(parse(text=GetG0Text("JBList", S=1)));
eval(parse(text=GetG0Text("MyPrintRoutines", S=1)));
eval(parse(text=GetG0Text("ADDDIRCONT", S=1)));
eval(parse(text=GetG0Text("NNList", S=1)));
eval(parse(text=GetG0Text("NRList", S=1)));
eval(parse(text=GetG0Text("NumReps", S=1)));
eval(parse(text=GetG0Text("CovRho", S=1)));
eval(parse(text=GetG0Text("DirectorForImage", S=1)));
eval(parse(text=GetG0Text("NNList", S=1)));
eval(parse(text=GetG0Text("LARSSEEKMETHOD", S=1)));
eval(parse(text=GetG0Text("LARSSeekFlag", S=1)));
eval(parse(text=GetG0Text("piAPriorData", S=1)));
eval(parse(text=GetG0Text("sigmaPriorData", S=1)));
eval(parse(text=GetG0Text("Interval", S=1)));
eval(parse(text=GetG0Text("GimmeList", S=1)));
eval(parse(text=GetG0Text("SigmaList", S=1)));
eval(parse(text=GetG0Text("OneVV2", S=1)));
NNList = c(25,50,100);
NRList = c(100,100,100);
Apt = unlist(strsplit(argeese, "/"));
if (length(Apt) > 1) {
  NewDir = paste(Apt[1:(length(Apt)-1)], collapse="/");
} else {
  NewDir = ""; 
}
if (any(unlist(list.files(NewDir)) == argeese)) {
try(source(paste(as.character(argeese), sep="")));
}
eval(parse(text=SetGText("MyPrintRoutines", S=1)));
eval(parse(text=SetGText("ADDDIRCONT", S=1)));
eval(parse(text=SetGText("NNList", S=1)));
eval(parse(text=SetGText("NRList", S=1)));
eval(parse(text=SetGText("GimmeList", S=1)));
eval(parse(text=SetGText("SigmaList", S=1)));
eval(parse(text=SetGText("OneVV2", S=1)));
eval(parse(text=SetGText("LARSSEEKMETHOD", S=1)));
eval(parse(text=SetGText("LARSSeekFlag", S=1)));
eval(parse(text=SetGText("piAPriorData",S=1)));
eval(parse(text=SetGText("sigmaPriorData",S=1)));
eval(parse(text=SetGText("Interval",S=1)));
eval(parse(text=SetGText("GimmeList",S=1)));
eval(parse(text=SetGText("SigmaList",S=1)));
eval(parse(text=SetGText("DefaultCorrelationXmatrix",S=1)));
eval(parse(text=SetGText("JBList",S=1)));

eval(parse(text=GetG0Text("tSeq",S=1)));
if (is.null(tSeq) || (is.numeric(tSeq))) {
  tSeq <- function (Number) 
{
    CN <- sub("e", "e", as.character(Number))
    CN <- sub("-", "n", as.character(CN))
    CN <- sub("\\+", "f", as.character(CN))
    CN <- sub("\\.", "p", as.character(CN))
    CN <- sub("0p", "p", as.character(CN))
    return(CN)
}
eval(parse(text=SetGText("tSeq",S=1)));

  DefineDefaultUseFunction();
  SetDefaultOtherParameters();
}

UpdateMethodDir(piSV = piSV, sigSV = sigSV, 
  LARSSEEKMETHOD=LARSSEEKMETHOD, LARSSeekFlag = LARSSeekFlag, 
  piAPriorData=piAPriorData, sigmaPriorData=sigmaPriorData,
  DefaultNoiseDF=DefaultNoiseDF);
  SetAllComparisonRCodeFunctions();
UpdateAllFunctions();

##print("Args 1 is ");
## print(args[1]); 
## argeese = args[1];
## args[1] = 0;
## args = as.numeric(args);  
## SkipLoad = 0;
## if (length(args) >= 3) {
##   SkipLoad = as.numeric(args[3]);
## } 
## source(paste(as.character(argeese), sep=""));
## print(paste("length(args) = ", length(args), sep=""));
print("You should now get an argeese"); flush.console();


}

OurPrintFile <- "";
OurErrorFile <- "";
IWillStillNotPrint <- 0;

AFilePrint <- function(AText, LockIn="No") {
      eval(parse(text=GetG0Text("OurPrintFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("JustPrintFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("APrintDir", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("IWillStillNotPrint", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("MyPrintEventLocation", "globalenv()", S=1)));
      if (is.null(MyPrintEventLocation) || !is.numeric(MyPrintEventLocation)) {
        try(MyPrintEventLocation <- 1000001);
      }
      APTime <- proc.time()[3];
      try(MyPrintEventLocation <- MyPrintEventLocation +1);
      if (is.character(JustPrintFile) && JustPrintFile %in% c("Do Not Print",
        "DoNotPrint", "donotprint", "DoNotPrint.txt", "donotprint.txt", "Do Not Print.txt")) {
        return(1);  
      }
      if (is.numeric(JustPrintFile) || is.null(JustPrintFile) || 
        (is.character(JustPrintFile) && JustPrintFile == "" ||
         JustPrintFile %in% c("Do Not Print", "DoNotPrint", "donotprint",
         "DoNotPrint.txt", "donotprint.xt", "Do Not Print.txt"))) {
        ##try(print(paste(MyPrintEventLocation, ":", APTime, " : ", AText, sep="")), silent=TRUE); 
        ##flush.console(); 
        return(1);
      }
      AF <- NULL;
      try(Oldwd <- getwd());
      try(setwd(APrintDir));
      ALLFiles <- unlist(list.files());
      if (any(ALLFiles == JustPrintFile)) {
        try(AF <- file(JustPrintFile, "at"));
      } else {
        try(AErrorPrint(paste("Hey: we tried to APrintFile. \n Dir is ", APrintDir,
          "\n and File is ", JustPrintFile, "\n But no print file  in there!", sep="")));
        try(AF <- file(JustPrintFile, "wt"));
      }  
      try(setwd(Oldwd));
      try(writeLines(text=
        paste(MyPrintEventLocation, ":", APTime, " : ", AText, sep=""), con=AF), silent=TRUE);
      try(close(AF));
      if (is.character(LockIn) && LockIn %in% c("Lock", "UnLock", "L", "U",
        "SumLock", "SumUnlock", "SumUnLock")) {
        try(ATT <- setwd(APrintLockDir));
        if (JustPrintFile %in% unlist(list.files())) {
          try(AF <- file(JustPrintFile, "wt"));
        } else {
          try(AF <- file(JustPrintFile, "at"));            
        }
        if (!is.null(AF)) {
          try(writeLines(text=
            paste(MyPrintEventLocation, ":", APTime, " : ", LockIn, ": ", AText, sep=""), con=AF), silent=TRUE);
          try(close(AF));
        }
        try(setwd(Oldwd));
      }
      if (is.null(IWillStillNotPrint) ||  !exists("IWillStillNotPrint") || 
        (is.numeric(IWillStillNotPrint)  && IWillStillNotPrint[1] != 1) ||
        (is.logical(IWillStillNotPrint) && IWillStillNotPrint == FALSE)) {
         try(print(paste(MyPrintEventLocation, ":", APTime, " : ", AText, sep="")), silent=TRUE); flush.console();  
      }
      eval(parse(text=SetGText("MyPrintEventLocation", "globalenv()", S=1)));
      return(1);
    }
    AErrorPrint <- function(AText) {
      eval(parse(text=GetG0Text("OurErrorFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("OurPrintFile", "globalenv()", S=1)));
      eval(parse(text=SetGText("JustPrintFile", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("AErrorDir", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("APrintDir", "globalenv()", S=1)));
      eval(parse(text=GetG0Text("IWillStillNotPrint", "globalenv()", S=1)));
      eval(parse(text=SetGText("MyPrintEventLocation", "globalenv()", S=1)));
      STTText <- unlist(strsplit(AText, "\n"));
      if (is.numeric(JustPrintFile) || is.null(JustPrintFile) || 
        (is.character(JustPrintFile) && JustPrintFile == "")) {
         print("*** ERROR ERROR ERROR ERROR ERROR ERROR ERROR")
         print(paste("*** ERROR: JustPrintFile is ", JustPrintFile, sep=""));
         for (ii in 1:length(STTText)) {
           print(paste("*** ERROR : ", MyPrintEventLocation, STTText[ii], sep="")); 
           flush.console(); 
         }
         print("***"); flush.console();
         print("***************************************************************");
         flush.console();
         return(1);
      }
      AF <- NULL;
      Oldwd <- getwd();
      try(setwd(AErrorDir));
      try(AF <- file(JustPrintFile, "at"));
      try(setwd(Oldwd));
      try(writeLines(text="*** ERROR ERROR ERROR ERROR ERROR ERROR ERROR", con=AF));
      try(writeLines(text=paste("*** ERROR at MyPrintEventLocation = ", MyPrintEventLocation, sep="")));
      try(writeLines(text=paste("*** ERROR: the Our PrintFile is ", OurPrintFile, sep="")));
      for (ii in 1:length(STTText)) {
         try(writeLines(text=STTText[ii], con=AF));
      }
    
      try(close(AF));
      if (is.null(IWillStillNotPrint) ||  !exists("IWillStillNotPrint") || 
        (is.numeric(IWillStillNotPrint)  && IWillStillNotPrint[1] != 1) ||
        (is.logical(IWillStillNotPrint) && IWillStillNotPrint == FALSE)) {
        print("*** ERROR ERROR ERROR ERROR ERROR ERROR ERROR")
        print(paste("*** ERROR: OurPrintFile is ", OurPrintFile, sep=""));
        ##print(paste("*** ERROR : ", MyPrintEventLocation, AText, sep=""));         
        for (ii in 1:length(STTText)) {
           print(paste("*** ERROR : ", MyPrintEventLocation, " : ", STTText[ii], sep="")); 
           flush.console(); 
        }
      }
      return(1);
    }
    try(eval(parse(text=SetGText("AErrorPrint", "globalenv()", S=1))));