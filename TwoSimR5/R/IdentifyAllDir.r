##IdentifyAllDir.R
IdentifyAllDirFunction <- function() {
  try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
eval(parse(text=GetG0Text("AllLargeDir", S=1)));
eval(parse(text=GetG0Text("AllSmallDir", S=1)));
if (length(unlist(
  list.files("/n/scratch06/airoldi_scratch/alenarcic/"))) > 0) {
  AllLargeDir = "/n/scratch06/airoldi_scratch/alenarcic/TwoSimR5//Large"
  AllSmallDir = "/n/scratch06/airoldi_scratch/alenarcic/TwoSimR5//Small"
  AllDir = "/n/scratch06/airoldi_scratch/alenarcic/TwoSimR5"
  dir.create(AllDir, showWarnings=FALSE, recursive=TRUE)  
  dir.create(AllLargeDir, showWarnings=FALSE, recursive=TRUE) 
  dir.create(AllSmallDir, showWarnings=FALSE, recursive=TRUE) 
} else if (length(unlist(list.files("C:/Stat"))) > 0) {
  AllDir = "c:/Stat/TwoSimR5/Saves";
  try(dir.create(AllDir, showWarnings=FALSE, recursive=TRUE));
  AllLargeDir = "c:/Stat/TwoSimR5/Saves/Large";
  try(dir.create(AllLargeDir, showWarnings=FALSE, recursive=TRUE));
  AllSmallDir = "c:/Stat/TwoSimR5/Saves/Small";
  try(dir.create(AllSmallDir, showWarnings=FALSE, recursive=TRUE));  
} else if (length(unlist(list.files("/Users/lenarcic/Documents/"))) > 0) {
  AllDir = "/Users/lenarcic/Documents/TwoSimR5"
  try(dir.create(AllDir, showWarnings=FALSE, recursive=TRUE));
  AllLargeDir = "/Users/lenarcic/Documents/TwoSimR5/Large"
  try(dir.create(AllLargeDir, showWarnings=FALSE, recursive=TRUE));
  AllSmallDir = "/Users/lenarcic/Documents/TwoSimR5/Small"
  try(dir.create(AllSmallDir, showWarnings=FALSE, recursive=TRUE));
} else if (length(unlist(list.files("~/CudaSDK"))) > 0) {
  AllDir = "~/Temp/TwoSimR5";
  try(dir.create(AllDir, recursive=TRUE, showWarnings=FALSE));
  AllLargeDir = "~/Temp/TwoSimR5/Large";
  try(dir.create(AllLargeDir, recursive=TRUE, showWarnings=FALSE));
  AllSmallDir = "~/Temp/TwoSimR5/Small";
  try(dir.create(AllSmallDir, recursive=TRUE, showWarnings=FALSE));
} else if (length(unlist(list.files("/Users/alenarc/Documents/"))) > 0) {
  AllDir = "/Users/alenarc/Documents/TwoSimR5"
  dir.create(AllDir, showWarnings=FALSE)
  AllLargeDir = "/Users/alenarc/Documents/TwoSimR5/Large"
  dir.create(AllLargeDir, showWarnings=FALSE)
  AllSmallDir = "/Users/alenarc/Documents/TwoSimR5/Small"
  dir.create(AllSmallDir, showWarnings=FALSE)
} else if (length(unlist(list.files("C:/Users/AlanX220/Documents/"))) > 0) {
  AllDir = "C:/Users/AlanX220/Documents/TwoSimR5"
  dir.create(AllDir, showWarnings=FALSE, recursive=TRUE)
  AllLargeDir = "C:/Users/AlanX220/Documents/TwoSimR5/Large"
  dir.create(AllDir, showWarnings=FALSE, recursive=TRUE)
  AllSmallDir = "C:/Users/AlanX220/Documents/TwoSimR5/Small"
  dir.create(AllDir, showWarnings=FALSE, recursive=TRUE);
} else if (length(unlist(list.files("C:/Users/AlanX230/Documents/"))) > 0) {
  AllDir = "C:/Users/AlanX230/Documents/TwoSimR5"
  dir.create(AllDir, showWarnings=FALSE, recursive=TRUE)
  AllLargeDir = "C:/Users/AlanX230/Documents/TwoSimR5/Large"
  dir.create(AllDir, showWarnings=FALSE, recursive=TRUE)
  AllSmallDir = "C:/Users/AlanX230/Documents/TwoSimR5/Small"
  dir.create(AllDir, showWarnings=FALSE, recursive=TRUE);
} else if (length(unlist(list.files("~/SpikeTest"))) > 0) {
  AllDir = "/netscr/alenarc/TwoSimR5";
  dir.create(AllDir, showWarnings=FALSE, recursive=TRUE);
  AllLargeDir = "/lustre/scr/a/l/alenarc/TwoSimR5";
  dir.create(AllLargeDir, showWarnings=FALSE, recursive=TRUE);
  AllSmallDir = "/netscr/alenarc/TwoSimR5/Small";
  dir.create(AllSmallDir, showWarnings=FALSE, recursive=TRUE);
} else {
  AllDir = "/Temp/";  AllLargeDir = "/Temp/Large";  AllSmallDir = "/Temp/Small";
  try(dir.create(AllLargeDir, showWarnings=FALSE, recursive=TRUE));
  try(dir.create(AllSmallDir, showWarnings=FALSE, recursive=TRUE));
  try(dir.create(AllDir, showWarnings=FALSE, recursive=TRUE));
}
  eval(parse(text=GetG0Text("DefaultTwoSimR5SaveDir", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("DefaultTwoSimR5SmallDir", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("DefaultTwoSimR5LArgeDir", "globalenv()", S=1)));
  if (!is.null(DefaultTwoSimR5SmallDir) && DefaultTwoSimR5SmallDir != 0) {
    AllDir = DefaultTwoSimR5SmallDir;
    AllSmallDir = DefaultTwoSimR5SmallDir;
    if (!is.null(DefaultTwoSimR5LargeDir)  && DefaultTwoSimR5LargeDir != 0) {
      AllLargeDir = DefaultTwoSimR5LargeDir;
    }
  } else if (!is.null(DefaultTwoSimR5LargeDir) && DefaultTwoSimR5LargeDir != 0) {
    AllDir = DefaultTwoSimR5LargeDir;
    AllLargeDir = DefaultTwoSimR5LargeDir;
    if (!is.null(DefaultTwoSimR5SmallDir)  && DefaultTwoSimR5SmallDir != 0) {
      AllSmallDir = DefaultTwoSimR5SmallDir;
    }
  } 
  try(eval(parse(text=SetGText("AllDir", envir="globalenv()", S=1))), silent=TRUE);
  try(eval(parse(text=SetGText("AllLargeDir", envir="globalenv()", S=1))), silent=TRUE);
  try(eval(parse(text=SetGText("AllSmallDir", envir="globalenv()", S=1))), silent=TRUE);
  try(eval(parse(text=GetG0Text("RSIMNAMESPACE", S=1))),  silent=TRUE);
  TRYRSIMText = "
  if (!is.numeric(RSIMNAMESPACE) &&  !is.null(RSIMNAMESPACE) &&
    !(is.numeric(RSIMNAMESPACE) &&  RSIMNAMESPACE == 0) &&
    is.environment(RSIMNAMESPACE)) {
    try(eval(parse(text-SetGText(\"AllDir\", 
      envir=\"RSIMNAMESPACE\", S=1))), silent=TRUE);
    try(eval(parse(text-SetGText(\"AllSmallDir\", 
      envir=\"RSIMNAMESPACE\", S=1))), silent=TRUE);
    try(eval(parse(text-SetGText(\"AllLargeDir\", 
      envir=\"RSIMNAMESPACE\", S=1))), silent=TRUE);
  }
  ";
  try(eval(parse(text=TRYRSIMText)), silent=TRUE);
  
TableOutName = "TableOut.csv";
 ADDDIRCONT = "";

 FTM = -1;
 KillerFun = -1;
 
eval(parse(text=GetG0Text("DefaultAContainDir", "globalenv()", S=1)));
if (!exists("DefaultAContainDir") || is.null(DefaultAContainDir) ||
  is.numeric(DefaultAContainDir)) {
  DefaultAContainDir <- AllDir;  
  eval(parse(text=SetGText("DefaultAContainDir", "globalenv()", S=1)));
}
OriginalDir = AllDir;
eval(parse(text=SetGText("OriginalDir", envir="globalenv()", S=1)));

eval(parse(text=GetG0Text("DefaultASmallContainDir", "globalenv()", S=1)));
if (!exists("DefaultASmallContainDir") || is.null(DefaultASmallContainDir) ||
  is.numeric(DefaultASmallContainDir)) {
  DefaultASmallContainDir <- AllSmallDir;  
  eval(parse(text=SetGText("DefaultASmallContainDir", "globalenv()", S=1)));
}
OriginalSmallDir = AllSmallDir;
eval(parse(text=SetGText("OriginalSmallDir", envir="globalenv()", S=1)));

eval(parse(text=GetG0Text("DefaultALargeContainDir", "globalenv()", S=1)));
if (!exists("DefaultALargeContainDir") || is.null(DefaultALargeContainDir) ||
  is.numeric(DefaultALargeContainDir)) {
  DefaultALargeContainDir <- AllLargeDir;  
  eval(parse(text=SetGText("DefaultALargeContainDir", "globalenv()", S=1)));
}
OriginalLargeDir = AllLargeDir;
eval(parse(text=SetGText("OriginalLargeDir", envir="globalenv()", S=1)));


eval(parse(text=GetG0Text("LARSSEEKMETHOD", envir="globalenv()", S=1)));
eval(parse(text=GetG0Text("LARSSeekFlag", envir="globalenv()", S=1)));
eval(parse(text=GetG0Text("piAPriorData", envir="globalenv()", S=1)));
eval(parse(text=GetG0Text("sigmaPriorData", envir="globalenv()", S=1)));
eval(parse(text=GetG0Text("TargetTablesDir", envir="globalenv()", S=1)));
eval(parse(text=GetG0Text("piSV", envir="globalenv()", S=1)));
eval(parse(text=GetG0Text("sigSV", envir="globalenv()", S=1)));


MethodDir = paste( "LSM", tSeq(LARSSEEKMETHOD),
  "LSF", tSeq(LARSSeekFlag),
  "piAP", tSeq(piAPriorData),
  "sPD", tSeq(sigmaPriorData), sep="");
TargetTablesDir = paste(OriginalSmallDir, "//", "Tables", ADDDIRCONT, MethodDir, sep="");
dir.create(TargetTablesDir, showWarnings=FALSE);
if (piSV != 0 || sigSV != 0) {
  MethodDir = paste(MethodDir, "piSV", tSeq(piSV), "sigSV", tSeq(sigSV), sep="");
}
if (exists("DefaultNoiseDF") && !is.null(DefaultNoiseDF)  && 
  DefaultNoiseDF > 0) {
  MethodDir = paste(MethodDir, "TDF", tSeq(DefaultNoiseDF), sep="");  
}
AllDir = paste(OriginalDir, "//", ADDDIRCONT, "CR", tSeq(CovRho), MethodDir, sep="");
AllLargeDir = paste(OriginalLargeDir, "//", ADDDIRCONT, "CR", tSeq(CovRho), MethodDir, sep="");
AllSmallDir = paste(OriginalSmallDir, "//", ADDDIRCONT, "CR", tSeq(CovRho), MethodDir, sep="");

  dir.create(AllDir, showWarnings = FALSE, recursive=TRUE);
  dir.create(AllSmallDir, showWarnings = FALSE, recursive=TRUE);
  dir.create(AllLargeDir, showWarnings = FALSE, recursive=TRUE);
  eval(parse(text=SetGText("TargetTablesDir", S=1)));
  eval(parse(text=SetGText("AllDir",S=1)));
  eval(parse(text=SetGText("AllLargeDir",S=1)));
  eval(parse(text=SetGText("AllSmallDir",S=1)));
  eval(parse(text=SetGText("MethodDir",S=1)));
  eval(parse(text=SetGText("OriginalDir",S=1)));
  eval(parse(text=SetGText("OriginalLargeDir",S=1)));
  eval(parse(text=SetGText("OriginalSmallDir",S=1)));
  eval(parse(text=SetGText("FTM",S=1)));
  eval(parse(text=SetGText("KillerFun",S=1)));
}

########################################################################
## If "OriginalDir" is target directory, UpdateMethodDir
##   creates a directory specific to the current flag sequences.
##
##
UpdateMethodDir <- function(piSV = piSV, sigSV = sigSV, 
  LARSSEEKMETHOD=LARSSEEKMETHOD, LARSSeekFlag = LARSSeekFlag, 
  piAPriorData=piAPriorData, sigmaPriorData=sigmaPriorData,
  DefaultNoiseDF=DefaultNoiseDF) {
  if (!exists("OriginalDir", globalenv()) ||
    !exists("OriginalSmallDir", globalenv())) {
    IdentifyAllDirFunction();
  }
  AFilePrint("UpdateMethodDir:  On Start "); flush.console();
  eval(parse(text=GetG0Text("ADDDIRCONT", "globalenv()", S=1)));
  if (is.null(ADDDIRCONT) || (is.numeric(ADDDIRCONT) && ADDDIRCONT == 0)) {
    ADDDIRCONT = "DefaultADDDIRCONT";
  }
  MethodDir = paste( "LSM", tSeq(LARSSEEKMETHOD),
    "LSF", tSeq(LARSSeekFlag),
    "piAP", tSeq(piAPriorData),
    "sPD", tSeq(sigmaPriorData), sep="");
  ##print(paste("UpdateMethodDir: Initial MethodDir = ", MethodDir, sep="")); flush.console();
  if (piSV != 0 || sigSV != 0) {
    print(paste("UpdateMethodDir: piSV is ", piSV, " and sigSV is ", sigSV, sep="")); flush.console();
    MethodDir = paste(MethodDir, "piSV", tSeq(piSV), "sigSV", tSeq(sigSV), sep="");
  }
  ## print(paste("UpdateMethodDiri: after piSV, MethodDir = ", MethodDir, sep="")); flush.console();
  if (exists("DefaultNoiseDF") && !is.null(DefaultNoiseDF)  && is.numeric(DefaultNoiseDF) &&
    DefaultNoiseDF > 0) {
    MethodDir = paste(MethodDir, "TDF", tSeq(DefaultNoiseDF), sep="");  
  }
  dir.create(paste(OriginalSmallDir, "//", ADDDIRCONT, sep=""), showWarnings=FALSE, recursive=TRUE);
  dir.create(paste(OriginalSmallDir, "//", ADDDIRCONT, "//", "CR", tSeq(CovRho), sep=""), showWarnings=FALSE, recursive=TRUE);

  dir.create(paste(OriginalLargeDir, "//", ADDDIRCONT, sep=""), showWarnings=FALSE, recursive=TRUE);
  dir.create(paste(OriginalLargeDir, "//", ADDDIRCONT, "//", "CR", tSeq(CovRho), sep=""), showWarnings=FALSE, recursive=TRUE);
  
  ## print(paste("UpdateMethodDiri: after tDF noise, MethodDir = ", MethodDir, sep="")); flush.console();
  AllDir = paste(OriginalSmallDir, "//", ADDDIRCONT, "//", "CR", tSeq(CovRho), "//", MethodDir, sep="");
  AllSmallDir = paste(OriginalSmallDir, "//", ADDDIRCONT, "//", "CR", tSeq(CovRho), "//", MethodDir, sep="");
  AllLargeDir = paste(OriginalLargeDir, "//", ADDDIRCONT, "//", "CR", tSeq(CovRho), "//", MethodDir, sep="");
  dir.create(AllDir, showWarnings = FALSE);
  dir.create(AllSmallDir, showWarnings = FALSE);
  dir.create(AllLargeDir, showWarnings = FALSE);
  eval(parse(text=SetGText("MethodDir")));
  eval(parse(text=SetGText("AllDir")));
  eval(parse(text=SetGText("AllLargeDir")));
  eval(parse(text=SetGText("AllSmallDir")));  
  
}  