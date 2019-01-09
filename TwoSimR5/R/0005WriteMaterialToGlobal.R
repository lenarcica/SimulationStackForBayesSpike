###############################################################################
##  0005WriteMaterialToGlobal.r  -- Alan Lenarcic
##
##    01-09-2019
##
##   Saves imporant variable settings to Global Environment
##
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


SetRSIText <- function(AText) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", RSIMNAMESPACE), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", RSIMNAMESPACE ) 
    }
    assign( \"", AText, "\", ", AText, ", RSIMNAMESPACE);
    ", sep="");
}

SetGText <- function(AText, envir = "globalenv()", S = 0) {
  if (S == 0) {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    ", sep=""));
  } else {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    ", sep=""));    
  }
}


SetRSIText <- function(AText) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", RSIMNAMESPACE ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", RSIMNAMESPACE ) 
    }
    assign( \"", AText, "\", ", AText, ", RSIMNAMESPACE) ;
    lockBinding( \"", AText, "\", RSIMNAMESPACE );
    ", sep="");
}

LockRSIText <- function(AText) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", RSIMNAMESPACE), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", RSIMNAMESPACE ) 
    }
    try(unlockBinding( \"", AText, "\", RSIMNAMESPACE ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", RSIMNAMESPACE);
    lockBinding( \"", AText, "\", RSIMNAMESPACE);
    ", sep="");
}


LockGText <- function(AText, envir="globalenv()", S=1)  {
  if (S==1) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep="");
  } else {
    paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=FALSE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=FALSE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep="");    
  }
}

GetRSIText <- function(AText) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", RSIMNAMESPACE), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", RSIMNAMESPACE ) 
    }
    ", AText, " <- get(\"", AText, "\", RSIMNAMESPACE);
    ", sep="");
}
GetGText <- function(AText, envir = "globalenv()",S=1) {
  if (S==0) {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    ", sep=""));
  } else {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = NULL;
    }
    ", sep=""));    
    
  }
}

GetG0Text <- function(AText, envir = "globalenv()", S=0 ) {
  if (S== 0) {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));
  } else {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));    
  }
}
SetSimKS <- function() {
  eval(parse(text=GetG0Text("SimSigmaNoise", S=1)));
  eval(parse(text=GetG0Text("SimKActiveLen", S=1)));
  SimKActiveLen <- function (kActiveLen, NLen, kLen)  {
   OneNum = round( kActiveLen + sqrt(kLen) * piSV * rnorm(1,0,1) )
   if (OneNum >= kLen) { OneNum = kLen -1; }
   if (OneNum >= NLen) { OneNum = NLen -1; }
   if (OneNum <= 0) { OneNum = 1; }
   return(OneNum);
  }
  SimSigmaNoise <- function (SigmaNoise) {
    if (sigSV > 0) {
      ##print("SuckSuckSuckSuck") ;
      return(sqrt(SigmaNoise^2 * rchisq(1, sigSV) /sigSV));
    } else{
      return(SigmaNoise);
    }
  }
  eval(parse(text=SetGText("SimSigmaNoise")));
  eval(parse(text=SetGText("SimKActiveLen")));
}

SetCurrentContainDir <- function() {
  eval(parse(text=GetG0Text("n", S=1)));
  eval(parse(text=GetG0Text("p", S=1)));
  eval(parse(text=GetG0Text("k", S=1)));
  eval(parse(text=GetG0Text("AllDir", S=1)));
  eval(parse(text=GetG0Text("AllSmallDir", S=1)));
  eval(parse(text=GetG0Text("AllLargeDir", S=1)));
  eval(parse(text=GetG0Text("sigma", S=1)));
  eval(parse(text=GetG0Text("NumReps", S=1)));
  eval(parse(text=GetG0Text("LARSSeekFlag", S=1)));
  eval(parse(text=GetG0Text("DefaultGenerateBetaVec", S=1)));
  eval(parse(text=GetG0Text("DefaultCorrelationXmatrix", S=1)));
  eval(parse(text=GetG0Text("AContainDir", S=1)));
  eval(parse(text=GetG0Text("CurrentContainDir", S=1)));
  eval(parse(text=GetG0Text("ALargeContainDir", S=1)));
  eval(parse(text=GetG0Text("CurrentLargeContainDir", S=1)));
  eval(parse(text=GetG0Text("ALargeContainDir", S=1)));
  eval(parse(text=GetG0Text("CurrentLargeContainDir", S=1)));
  CurrentContainDir <- paste(AllSmallDir, "//",  DirNamingClenture( 
    DefaultGenerateBetaVec, DefaultCorrelationXmatrix,
    NLen = n, kLen = p, kActiveLen = k, SigmaNoise = sigma, 
    NumReps = NumReps, LARSSeekFlag  = LARSSeekFlag), sep="");
    dir.create(CurrentContainDir, showWarnings = FALSE);  
  CurrentSmallContainDir <- paste(AllSmallDir, "//",  DirNamingClenture( 
    DefaultGenerateBetaVec, DefaultCorrelationXmatrix,
    NLen = n, kLen = p, kActiveLen = k, SigmaNoise = sigma, 
    NumReps = NumReps, LARSSeekFlag  = LARSSeekFlag), sep="");
    dir.create(CurrentContainDir, showWarnings = FALSE); 
  CurrentLargeContainDir <- paste(AllLargeDir, "//",  DirNamingClenture( 
    DefaultGenerateBetaVec, DefaultCorrelationXmatrix,
    NLen = n, kLen = p, kActiveLen = k, SigmaNoise = sigma, 
    NumReps = NumReps, LARSSeekFlag  = LARSSeekFlag), sep="");
    dir.create(CurrentContainDir, showWarnings = FALSE); 
    dir.create(CurrentLargeContainDir, showWarnings = FALSE); 
    dir.create(CurrentSmallContainDir, showWarnings = FALSE); 
    eval(parse(text=SetGText("CurrentContainDir", S=1)));
    eval(parse(text=SetGText("CurrentLargeContainDir", S=1)));
    eval(parse(text=SetGText("CurrentSmallContainDir", S=1)));
  AContainDir <- CurrentContainDir;             
  ASmallContainDir <- CurrentSmallContainDir;   
  ALargeContainDir <- CurrentLargeContainDir;   
  eval(parse(text=SetGText("CurrentContainDir", S=1)));
  eval(parse(text=SetGText("CurrentLargeContainDir", S=1)));
  eval(parse(text=SetGText("CurrentSmallContainDir", S=1)));
  eval(parse(text=SetGText("AContainDir", S=1)));
  eval(parse(text=SetGText("ALargeContainDir", S=1)));
  eval(parse(text=SetGText("ASmallContainDir", S=1)));
  ##source("c://Stat//2008Summer//LarsProject//code//EfficientSimulator.R")
  ##   Customized SimMeData


  ###########################################################################
  ## LockFile Mechanics
  ##     ## We rename some defaults here, though it's not necessary
  eval(parse(text=GetG0Text("DefaultAContainDir", S=1)));
  eval(parse(text=GetG0Text("DefaultASmallContainDir", S=1)));
  eval(parse(text=GetG0Text("DefaultLFDir", S=1)));
  eval(parse(text=GetG0Text("DefaultLockFile", S=1)));
  DefaultAContainDir <- CurrentSmallContainDir;
  DefaultASmallContainDir <- CurrentSmallContainDir;
  DefaultLFDir <- DefaultASmallContainDir;
  DefaultLockFile <- paste(CurrentSmallContainDir, "//", "LockFile.txt", sep="");  
  eval(parse(text=SetGText("DefaultAContainDir", S=1)));
  eval(parse(text=SetGText("DefaultASmallContainDir", S=1)));
  eval(parse(text=SetGText("DefaultLFDir", S=1)));
  eval(parse(text=SetGText("DefaultLockFile", S=1)));
}

SetDefaultScore <- function(NameVar, Val) {
if (is.numeric(Val)) {
MyCatT <-  paste(
"eval(parse(text=GetG0Text(\"", NameVar, "\", S=1)));
if (is.numeric(", NameVar, ") && ", NameVar, " == 0) {
  ", NameVar, "<- ", Val, ";
  eval(parse(text=SetGText(\"", NameVar, "\", S=1)));
} else if (is.numeric(", NameVar, ") && ", NameVar, " != ", Val, ") {
  eval(parse(text=SetGText(\"", NameVar, "\", S=1)));
} else {
  ", NameVar, " <- ", Val, ";
  eval(parse(text=SetGText(\"", NameVar, "\", S=1)));
}", sep="");
eval(parse(text=MyCatT));
} else if (is.character(Val)) {
MyCatT <-  paste(
"eval(parse(text=GetG0Text(\"", NameVar, "\", S=1)));
if (is.numeric(", NameVar, ") && ", NameVar, " == 0) {
  ", NameVar, "<- \"", Val, "\";
  eval(parse(text=SetGText(\"", NameVar, "\", S=1)));
} else if (is.character(", NameVar, ") && ", NameVar, " != \"", Val, "\") {
  eval(parse(text=SetGText(\"", NameVar, "\", S=1)));
} else {
  ", NameVar, " <- \"", Val, "\";
  eval(parse(text=SetGText(\"", NameVar, "\", S=1)));
}", sep="");
eval(parse(text=MyCatT));
}
}

SetDefaultOtherParameters <- function() {

SetDefaultScore("SkipLoad", 0);
SetDefaultScore("Inverval", .95);
SetDefaultScore("MyPrintRoutines", "");
SetDefaultScore("LARSSEEKMETHOD", 1);
SetDefaultScore("LARSSeekFlag", -.8);
SetDefaultScore("piAPriorData", .1);
SetDefaultScore("sigmaPriorData", .1);
SetDefaultScore("SigmaPriorData", .1);
SetDefaultScore("TableOutName", "TableOut.csv");
SetDefaultScore("Inverval", .95);
SetDefaultScore("DefaultHorseShoeSamples", 250);
SetDefaultScore("NumDetails", DefaultNumDetails);
SetDefaultScore("Inteval", .95);
SetDefaultScore("Inverval", .95);
DefaultIntAssessOutName = paste("IntAssessOut", tSeq(Interval), ".csv", sep="");
SetDefaultScore("IntAssessOutName", DefaultIntAssessOutName);
SetDefaultScore("NumIntervalAssess", 12);
SetDefaultScore("NumChains", 3);
SetDefaultScore("MandatoryOnFunction", NULL);


}

RunGetArgs <- function() {

eval(parse(text=GetG0Text("Myargs", S=1)));
eval(parse(text=GetG0Text("argeese", S=1)));
eval(parse(text=GetG0Text("MyDefaultArgs", S=1)));
if (is.null(MyDefaultArgs) || length(MyDefaultArgs) == 1 && MyDefaultArgs[1] == 0) {
Myargs = commandArgs(TRUE);
} else {
  Myargs <- MyDefaultArgs;
}
if (length(Myargs) >= 1 && is.character(Myargs[1])  && !is.na(Myargs[1])){ 
  argeese = Myargs[1];
  if (length(Myargs) >= 2 && is.character(Myargs)) {
    ForceFunction <- Myargs[2];
  }
  if (length(Myargs) >= 5) {
    Myargs <- c(argeese,1,1,1,1);
  }
  try(print(paste("We will be setting argeese for paper file is ", argeese, sep="")))
} else if (length(list.files("c:/Stat")) > 0) {
  argeese <- paste(DropHome, "//TwoLassoPackage//ParTabPaper.R", sep="");
  Myargs = c(argeese, 1,1,1,1);
} else if (length(list.files("~/DownloadPackages")) > 0) {
  argeese <- "~/DownloadPackages/ParTabPaper.R";
  Myargs = c(argeese, 1,1,1,1);
} else if (length(list.files("~/MyPackages")) > 0) {
  argeese <- "~/KillTwoLasso/ParTabPaper.R";
  Myargs = c(argeese, 1,1,1,1);
} else {
  print("Please Set RunGetArgs Better for this Server!"); flush.console();
}
if (is.na(argeese)) {
  print("Please Set RunGetArgs Better for this Server, argeese is NA!");
  flush.console();
}
eval(parse(text=SetGText("Myargs", S=1)));
eval(parse(text=SetGText("argeese", S=1)));

}

  MyAttacher <- function(AStringList) {
  AOverallPaste <- paste(",
      for (ii in 1:length(", AStringList, ")) {
        MyName <- names(", AStringList, ")[ii]
          MyPaster <- paste(\"
          eval(parse(text=GetG0Text(\\\"\", MyName, \"\\\", S=1)));
            \", MyName, \" <- ", AStringList, "[[ii]];
          eval(parse(text=SetGText(\\\"\", MyName, \"\\\", S=1))); 
          \", sep=\"\");
          eval(parse(text=MyPaster));
      }
      ", sep="");        
     return(AOverallPaste);
  }
  
SetListOfParams <- function() {
  eval(parse(text=GetG0Text("ListOfParams", S=1)));
  eval(parse(text=GetG0Text("n", S=1)));
  eval(parse(text=GetG0Text("k", S=1))); 
  eval(parse(text=GetG0Text("sigma", S=1)));
  eval(parse(text=GetG0Text("p", S=1)));
  eval(parse(text=GetG0Text("piSV", S=1)));
  eval(parse(text=GetG0Text("sigSV", S=1))); 
  ListOfParams <- list(n=n, k=k, p = p, sigma=sigma, NN = n, 
    piSV=piSV, sigSV=sigSV, kLen = p,
    kActiveLen=k, SigmaNoise = sigma);       
  try(eval(parse(text=SetGText("ListOfParams", S=1))), silent=TRUE);
}

SetContainDir <- function() {
  eval(parse(text=GetG0Text("AllDir", S=1)));
  eval(parse(text=GetG0Text("DirectorForImage", S=1)));
  eval(parse(text=GetG0Text("LARSSEEKMETHOD", S=1))); 
  eval(parse(text=GetG0Text("LARSSeekFlag", S=1)));
  eval(parse(text=GetG0Text("piAPriorData", S=1)));
  eval(parse(text=GetG0Text("piSV", S=1))); 
  eval(parse(text=GetG0Text("sigSV", S=1)));
  eval(parse(text=GetG0Text("DefaultNoiseDF", S=1)));
  eval(parse(text=GetG0Text("ContainDir", S=1)));
  ContainDir <-  paste(AllDir, "//", DirectorForImage, "LSM", tSeq(LARSSEEKMETHOD), 
  "LSF", tSeq(LARSSeekFlag), "pPr", tSeq(piAPriorData), 
  "piSV", tSeq(piSV), "sigSV", tSeq(sigSV), sep="");
  eval(parse(text=GetG0Text("DefaultNoiseDF", S=1)));
  if (exists("DefaultNoiseDF") && !is.null(DefaultNoiseDF)  && 
    is.numeric(DefaultNoiseDF) &&
    DefaultNoiseDF > 0) {
    ContainDir = paste(ContainDir, "TDF", tSeq(DefaultNoiseDF), sep="");  
  }
  eval(parse(text=SetGText("ContainDir", S=1)));
  dir.create(ContainDir, showWarnings=FALSE);
}
  
GetNumPrint <- function() {
  eval(parse(text=GetG0Text("UseFunctions", S=1)));
  eval(parse(text=GetG0Text("n", S=1)));
  eval(parse(text=GetG0Text("p", S=1)));
  eval(parse(text=GetG0Text("TryPrint", S=1)));
  eval(parse(text=GetG0Text("NumPrint", S=1)));        
  if (length(UseFunctions) * n * p > 300) {
    TryPrint = 2 + round(40 / ( n * p * length(UseFunctions)^(1/4) ) );
    if (TryPrint > NumPrint) {
    
    } else {
      NumPrint = TryPrint;
    }
  }
  eval(parse(text=SetGText("TryPrint", S=1)));
  eval(parse(text=SetGText("NumPrint", S=1)));
  
}