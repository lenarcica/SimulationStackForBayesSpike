


SetRSIText <- function(AText) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", RSIMNAMESPACE), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", RSIMNAMESPACE ), silent=TRUE); 
    }
    try(assign( \"", AText, "\", ", AText, ", RSIMNAMESPACE), silent=TRUE);
    ", sep="");
}

SetGText <- function(AText, envir = "globalenv()", S = 0) {
  if (!is.character("AText")) {
    print(paste("Hey, SetGText, you passed something bad to me. ")); flush.console();
  }
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
    try(assign( \"", AText, "\", ", AText, ", ", envir, " ), silent=TRUE);
    ", sep=""));    
  }
}


SetRSIText <- function(AText) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", RSIMNAMESPACE ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", RSIMNAMESPACE ), silent=TRUE); 
    }
    try(assign( \"", AText, "\", ", AText, ", RSIMNAMESPACE), silent=TRUE); ;
    try(lockBinding( \"", AText, "\", RSIMNAMESPACE ), silent=TRUE);
    ", sep="");
}

LockRSIText <- function(AText) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", RSIMNAMESPACE), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", RSIMNAMESPACE ), silent=TRUE); 
    }
    try(unlockBinding( \"", AText, "\", RSIMNAMESPACE ), silent=TRUE); 
    try(assign( \"", AText, "\", ", AText, ", RSIMNAMESPACE), silent=TRUE);
    try(lockBinding( \"", AText, "\", RSIMNAMESPACE), silent=TRUE);
    ", sep="");
}



LockGText <- function(AText, envir="globalenv()", S=1)  {
  if (!is.character("AText")) {
    print(paste("Hey, LockGText, you passed something bad to me. ")); flush.console();
  }
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
      try(unlockBinding( \"", AText, "\", RSIMNAMESPACE ), silent=TRUE); 
    }
    ", AText, " <- get(\"", AText, "\", RSIMNAMESPACE);
    ", sep="");
}
GetGText <- function(AText, envir = "globalenv()",S=1) {
  if (!is.character(AText)) {
    print("GetG0Text: some error, AText is not a character"); flush.console();
  }
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    ", sep="");
}

GetG0Text <- function(AText, envir = "globalenv()", S=0 ) {
  if (!is.character(AText)) {
    print("GetG0Text: some error, AText is not a character"); flush.console();
  }
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
  eval(parse(text=GetG0Text("sigma", S=1)));
  eval(parse(text=GetG0Text("NumReps", S=1)));
  eval(parse(text=GetG0Text("LARSSeekFlag", S=1)));
  eval(parse(text=GetG0Text("DefaultGenerateBetaVec", S=1)));
  eval(parse(text=GetG0Text("DefaultCorrelationXmatrix", S=1)));
  eval(parse(text=GetG0Text("AContainDir", S=1)));
  eval(parse(text=GetG0Text("CurrentContainDir", S=1)));
  CurrentContainDir <- paste(AllDir, "//",  DirNamingClenture( 
    DefaultGenerateBetaVec, DefaultCorrelationXmatrix,
    NLen = n, kLen = p, kActiveLen = k, SigmaNoise = sigma, 
    NumReps = NumReps, LARSSeekFlag  = LARSSeekFlag), sep="");
    dir.create(CurrentContainDir, showWarnings = FALSE);  
  eval(parse(text=SetGText("CurrentContainDir", S=1)));
  AContainDir <- CurrentContainDir;             
  eval(parse(text=SetGText("CurrentContainDir", S=1)));
  eval(parse(text=SetGText("AContainDir", S=1)));
  ##source("c://Stat//2008Summer//LarsProject//code//EfficientSimulator.R")
  ##   Customized SimMeData


  ###########################################################################
  ## LockFile Mechanics
  ##     ## We rename some defaults here, though it's not necessary
  eval(parse(text=GetG0Text("DefaultAContainDir", S=1)));
  eval(parse(text=GetG0Text("DefaultLFDir", S=1)));
  eval(parse(text=GetG0Text("DefaultLockFile", S=1)));
  DefaultAContainDir <- CurrentContainDir;
  DefaultLFDir <- DefaultAContainDir;
  DefaultLockFile <- paste(CurrentContainDir, "//", "LockFile.txt", sep="");  
  eval(parse(text=SetGText("DefaultAContainDir", S=1)));
  eval(parse(text=SetGText("DefaultAContainDir", S=1)));
  eval(parse(text=SetGText("DefaultLockFile", S=1)));
}