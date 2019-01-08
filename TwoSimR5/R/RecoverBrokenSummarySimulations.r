if (FALSE) {
  AFile <- "c:/Stat/TryTable.csv";
  AMatrix <- rbind(c(1,2,4,5),c(6,7,8,9));
  col.names=c("One", "Two", "Three", "Four");
  row.names=c("OneRow", "TwoRow");
  sep=", ";
  secure.write.table(file=AFile, AMatrix, col.names=col.names, row.names=row.names,
    sep=", ");
 
    try(ACon <- file(AFile, "wt")); 
    try(writeLines(text=paste("\"", col.names, "\"", collapse=", ", sep=""),
           con=ACon));
        for (tt in 1:NROW(AMatrix)) {
          writeLines(text="Garbage", con=ACon);
          try(writeLines(text=paste("\"", row.names[tt], "\", ",
            paste(AMatrix[tt,], collapse=", "), sep=""), con=ACon));
        }
    writeLines(text="Garbage", con=ACon);
        try(close(ACon));  ART <- NULL;
        
    try(ASS <- secure.read.table(file=AFile, nCols=4,header=TRUE))
    
    AFile <- c("c:/Stat/Long");
    AMatrix <- matrix(rnorm(105*5),105,5);
    colnames(AMatrix) <- c("a", "b", "c", "d", "e");
    rownames(AMatrix) <- 1:NROW(AMatrix);
    write.table(AMatrix, file=AFile, sep=", ", col.names=colnames(AMatrix),
      row.names=rownames(AMatrix), quote=FALSE, eol="\n");
      
    AFF <- file(AFile, "rt");
    A1 <- readLines(AFF,1);
    Done = 0;
    RT <- 1;
    OneOn <- readLines(AFF, RT);
    AOn <- 1;
    ReadLen <- length(OneOn);
    OnRead <- length(OneOn);
    while(Done == 0) {
      print(paste("AOn = ", AOn, " On read = ", OnRead, " and ReadLen = ", ReadLen, sep=""));   flush.console();
      OneOn <- readLines(AFF,RT);
      if (!is.null(OneOn) && length(OneOn) >= 1) {
        OnRead <- length(OneOn);
        ReadLen <- ReadLen + OnRead;
      } else {
        Done = 1;
      }
      AOn <- AOn +1;      
    }
    close(AFF);
      
  
}
secure.write.table <- function(file = NULL, AMatrix, col.names=NULL, row.names=NULL,
  sep=", ", AFile = NULL) {
  AGiveTries <- 100;
  if (is.null(AFile)  && !is.null(file)) { 
    AFile <- file;
  }
  if (is.null(col.names)) { ncolnames = FALSE; } else {
    ncolnames = col.names;
  }
  if (is.null(row.names)) { nrownames = FALSE; } else {
    nrownames = row.names;
  }
  try(write.table(AMatrix, file=AFile, col.names=ncolnames,
    row.names=nrownames, sep=sep, quote=FALSE, eol="\n"), silent=TRUE);
  ART <- NULL;
  if (!is.null(col.names) && (length(ncolnames) >= 2 ||
      !(is.logical(ncolnames) && ncolnames == FALSE))) {
      try(ART <- read.table(file=AFile, header=TRUE,
      sep=","), silent=TRUE);
    } else {
      try(ART <- read.table(file=AFile, header=FALSE,
      sep=","), silent=TRUE);        
  }
  Good <- FALSE;
  ATText <- "
  if (!is.null(ART) && NROW(ART) == NROW(AMatrix) &&
    NCOL(AMatrix) == NCOL(ART)) {
    Good <- TRUE;
    return(TRUE);      
  }
  if (!is.null(ART) && (is.logical(row.names) && row.names[1] == TRUE) &&
    NCOL(AMatrix) == NCOL(ART)-1) {
    Good <- TRUE;
    return(TRUE);  
  }
  if (!is.null(ART) && (length(row.names) == NROW(ART)) &&
    NCOL(AMatrix) == NCOL(ART)-1) {
    Good <- TRUE;
    return(TRUE);  
  }
  ";
  try(eval(parse(text=ATText)));
  if (Good) { return(TRUE); }
  if (is.null(ART) || NROW(ART) != NROW(AMatrix) ||
    NCOL(AMatrix) != NCOL(ART)) {
    AFilePrint(paste("secure.write.table: hey, bad save for ", AFile, sep=""));      
    GoodWrite <- 0;  AGiveTries<- 100;
    while(GoodWrite == 0 && AGiveTries >= 0) {
      ACon <- NULL;
      try(ACon <- file(AFile, "wt"), silent=TRUE);
      if (is.null(ACon)) {
        AFilePrint(paste("   secure.write: fail again (CON) tries = ", AGiveTries, sep=""));
      } else {
        if (is.null(col.names)  || (is.logical(col.names) && col.names==FALSE)) {  
        } else if (is.logical(col.names) && col.names==TRUE) {
          try(writeLines(text=paste("\"", colnames(AMatrix), "\"", collapse=", ", sep=""),
             con=ACon), silent=TRUE);
        } else {
          try(writeLines(text=paste("\"", col.names, "\"", collapse=", ", sep=""),
             con=ACon), silent=TRUE);            
        }
        if (length(row.names) == NROW(AMatrix) && (NROW(AMatrix) >= 2  ||
          (!is.logical(row.names) && !is.null(row.names))) ) {
          for (tt in 1:NROW(AMatrix)) {
            try(writeLines(text=paste("\"", row.names[tt], "\", ",
              paste(AMatrix[tt,], collapse=", "), sep=""), con=ACon));
          }
        } else if (is.logical(row.names) && row.names==TRUE) {
           for (tt in 1:NROW(AMatrix)) {
            try(writeLines(text=paste("\"", rownames(AMatrix)[tt], "\", ",
              paste(AMatrix[tt,], collapse=", "), sep=""), con=ACon));
          }           
        } else if ((is.logical(row.names) && row.names==FALSE) || is.null(row.names)) {
           for (tt in 1:NROW(AMatrix)) {
            try(writeLines(text=paste("",
              paste(AMatrix[tt,], collapse=", "), sep=""), con=ACon));
           }           
        } else {
          AFilePrint(paste("Hey: I don't secure.write.table know what ",
            "to do about row.names!!! = (",
            paste(row.names, collapse=", ", sep=""), ");", sep=""));
          return(FALSE);
        } 
        try(close(ACon), silent=TRUE);  ART <- NULL;
        if (!is.null(col.names) && (length(ncolnames) >= 2 ||
            !(is.logical(ncolnames) && ncolnames == FALSE))) {
          try(ART <- read.table(file=AFile, header=TRUE,
            sep=","), silent=TRUE);
        } else {
          try(ART <- read.table(file=AFile, header=FALSE,
            sep=","), silent=TRUE);            
        }
        if (!is.null(ART) && NCOL(ART) == NCOL(AMatrix) && 
          NCOL(ART) == NCOL(AMatrix)) {
          AFilePrint(paste("    secure.write: good hit on AGiveTries = ", AGiveTries, sep=""));
          GoodWrite <- 1;
          return(TRUE);
        }
        if (!is.null(ART) && (is.logical(row.names) && row.names[1] == TRUE) &&
         NCOL(AMatrix) == NCOL(ART)-1) {
         Good <- TRUE;
         return(TRUE);  
        }
        if (!is.null(ART) && (length(row.names) == NROW(ART)) &&
          NCOL(AMatrix) == NCOL(ART)-1) {
          Good <- TRUE;
          return(TRUE);  
        }
      } 
      AGiveTries <- AGiveTries-1;
    }
  }
  return(FALSE)
}

secure.read.text.lines <- function(file, lines=20, FurtherText="") {
  try(Atries <- 5);  try(AFile <- file);
  eval(parse(text=SetGText("AFile", "globalenv()", S=1)));
  while(Atries >= 0) {
    ACon <- NULL; try(ACon <- file(AFile, "rt"));
    if (!is.null(ACon)) {
      try(ART <- readLines(ACon,lines), silent=TRUE);
      try(close(ACon));
      if (length(ART) >= 1) {
         return(ART);
      } else {
        try(AFilePrint(paste("secure.read.text.lines: ", AFile, " no lines read. try = ", Atries, " : ", FurtherText, sep="")));
      }
    } else {
       try(AFilePrint(paste("secure.read.text.lines: ", AFile, " No connection established. try = ", Atries, " : ", FurtherText, sep="")));
    }
    try(Atries <- Atries - 1);
  }  
  return(NULL);
}
secure.read.table <- function(file=NULL, nCols=0, header=TRUE,
  sep=", ", FurtherText="", AFile=NULL) {
  AGiveTries <- 8;
  if (is.null(AFile) && !is.null(file)) {
    try(AFile <- file);
  }
  try(eval(parse(text=SetGText("AFile", "globalenv()", S=1))));
  ART <- NULL;
  try(ART <- read.table(file=AFile, header=header,
    sep=","), silent=TRUE);
  Good <- FALSE;
  ATText <- "
  if (!is.null(ART) && (nCols <= 0 || NCOL(ART) == nCols)) {
    Good <- TRUE;
    return(ART);      
  }
  if (header==FALSE && (nCols >= 1 && nCols == NCOL(ART)-1)) {
     Good <- TRUE;
     return(ART);
  }
  ";
  try(eval(parse(text=ATText)));
  if (Good) { return(ART); }
  AFilePrint(paste("secure.read.table: hey, bad read for ", AFile, " : ", FurtherText, sep=""));      
  AGiveTries <- 5;  ReturnMatrix <- NULL;
  while(is.null(ReturnMatrix) && AGiveTries >= 0) {
    ACon <- NULL; ReturnMatrix <- NULL;
    try(ReturnMatrix <- secure.read.table.try(AFile = AFile, nCols=nCols, header=header, AGiveTries=AGiveTries));
    if (is.null(ReturnMatrix)) {
      try(AFilePrint(paste("secure.read.table:  Fail again, AGiveTries = ",
        AGiveTries, " : ", FurtherText, sep="")));
    }
    AGiveTries <- AGiveTries -1;
  }
  if (!is.null(ReturnMatrix)) {
    AFilePrint(paste("secure.read.table: we will now secure.write because ",
      "we had a successful hit on: ", AFile, ".", sep=""));
    try(AFilePrint(paste("secure.read.table: Line 1 = ",
       rownames(ReturnMatrix)[1,1], ": ", 
       paste(ReturnMatrix[1,], collapse=", "), sep="")));
    try(AFilePrint(paste("secure.read.table: Line 2 = ",
       rownames(ReturnMatrix)[1,1], ": ", 
       paste(ReturnMatrix[1,], collapse=", "), sep="")));       
    ATT <- " WeGot <- FALSE;
      WeGot <- secure.write.table(AMatrix=ReturnMatrix, file=AFile, col.names=colnames(ReturnMatrix),
      row.names=rownames(ReturnMatrix), sep=\", \"));
    ";
    try(eval(parse(text=ATT)));
    if (is.null(WeGot) || WeGot == FALSE) {
      try(AFilePrint(paste("secure.read.table: we were not able to secure.write a revision for ", AFile, sep="")));
    }
  }
  return(ReturnMatrix);
}

secure.read.table.header.1 <- function() {
  RetText <- "
      if (header==TRUE) {
        try(ColNames <- readLines(con=ACon,1));  OnLine <- OnLine+1;
        if (is.null(ColNames)) {
            AFilePrint(\"Hey secure read header, no header to read!\");
            try(close(ACon));
            return(NULL);
        } else {
          try(AQ <- paste(unlist(strsplit(ColNames, \"\\\"\")), collapse=\"\"));
          try(AQ <- paste(unlist(strsplit(AQ, \" \")), collapse=\"\"));          
          try(new.col.names <- unlist(strsplit(AQ, \",\")));
          if (nCols > 0 && nCols != length(new.col.names)) {
            AFilePrint(paste(\" Hey Error, secure.read.table.header.1: \",
              \"nCols = \", nCols, \"  but length(new.col.names) = \",
              length(new.col.names), sep=\"\"));
            try(close(ACon));  return(NULL); 
          } else if (nCols <= 0) {
            nCols <- length(new.col.names);
          }
        }
      } else {
        new.col.names = NULL;
      }  
  ";
}
secure.read.table.header.2 <- function() {
  RetText <- "
      if (header==TRUE) {
        try(ColNames <- readLines(con=ACon,1));  OnLine <- OnLine+1;
        if (is.null(ColNames)) {FailThisTime = 1;} else{
          try(AQ <- paste(unlist(strsplit(ColNames, \"\\\"\")), collapse=\"\"));
          try(AQ <- paste(unlist(strsplit(AQ, \" \")), collapse=\"\"));          
          try(new.col.names <- unlist(strsplit(AQ, \",\")));
          try(colnames(ReturnMatrix) <- new.col.names);
          if (nCols > 0 && nCols != length(new.col.names)) {
            AFilePrint(paste(\" Hey Error, secure.read.table.header.2: \",
              \"nCols = \", nCols, \"  but length(new.col.names) = \",
              length(new.col.names), \" 
              This is bad because we already succeeded! \", sep=\"\"));
            try(close(ACon));
            return(NULL); 
          }
        }
      } else {
        new.col.names = NULL;
      }  
  ";
}
secure.read.table.first.line <- function() {
  Ret <- "
    try(RLL <- readLines(ACon,1));   OnLine <- OnLine+1;
    if (is.null(RLL) || length(RLL) <= 0) {
      AFilePrint(paste(\"Hey, we tried to secure.read \", AFile, 
        \" But there is no first line\", sep=\"\"));
      NoLine <- TRUE;
      GoodLines <- 0;
    } else if (nchar(RLL) == 0) {
      NoLine <- FALSE;  GoodLines <- 0;
    } else {
        NoLine <- FALSE;
        AT <- paste(unlist(strsplit(RLL[1], \" \")), collapse=\"\");
        AT <- paste(unlist(strsplit(AT, \"\\\"\")), collapse=\"\");
        AQ <- unlist(strsplit(AT, \",\"));
        if (is.null(new.col.names)  || nCols <= 0) {
          if ( !is.na(suppressWarnings(as.numeric(AQ[1]))) ) {
             nCols <- length(AQ);  nColRows <- nCols;
             GoodLines <- 1;
          } else if (length(AQ) == 1) {
            GoodLines <- 0;
            AFilePrint(paste(\"  Hey line :\", OnLine, 
              \" is not good for first row. \", sep=\"\"));
          } else {
            nCols <- length(AQ)-1;  GoodLines <- 1;
            RowNames <- rep(\"\", length(RLL)); nColRows <- nCols+1;
          }
        } else {
           NoLine <- FALSE;
           if (nCols <= 0) { nCols = length(new.col.names); }
           if (length(AQ) == nCols+1) {  nColRows <- nCols+1;
             GoodLines <- 1;
           } else if (length(AQ) != nCols) {
             AFilePrint(paste(\"Error, AQ = (\", 
               paste(AQ, collapse=\", \"), \") but nCols=\", nCols, sep=\"\"));
           } else { nColRows <- nCols;  GoodLines <- 1;
        }
      }   
    } 
    ";
}
secure.read.table.a.line <- function() {
  Ret <- "
    try(RLL <- readLines(ACon,1)); OnLine <- OnLine+1;
    if (is.null(RLL) || length(RLL) <= 0) {
      NoLine <- TRUE;
    } else if (nchar(RLL) == 0) {
      NoLine <- FALSE;
    } else {
        NoLine <- FALSE;
        AT <- paste(unlist(strsplit(RLL[1], \" \")), collapse=\"\");
        AT <- paste(unlist(strsplit(AT, \"\\\"\")), collapse=\"\");
        AQ <- unlist(strsplit(AT, \",\"));
        if (length(AQ) == nColRows) {
          GoodLines <- GoodLines+1;
        } else {
          AFilePrint(paste(\"  We found OnLine = \", OnLine, 
            \" to be faulty, read \", GoodLines, \": GoodLines.\", sep=\"\"));  
        }
    }    
    ";
}
secure.read.table.a.line.fill <- function() {
  Ret <- "
    try(RLL <- readLines(ACon,1)); OnLine <- OnLine+1;
    if (is.null(RLL) || length(RLL) <= 0) {
      NoLine <- TRUE;
    } else if (nchar(RLL) == 0) {
      NoLine <- FALSE;
    } else {
        NoLine <- FALSE;
        AT <- paste(unlist(strsplit(RLL[1], \" \")), collapse=\"\");
        AT <- paste(unlist(strsplit(AT, \"\\\"\")), collapse=\"\");
        AQ <- unlist(strsplit(AT, \",\"));
        if (length(AQ) == nColRows) {
          GoodLines <- GoodLines+1;
          if (nCols == nColRows) {
            try(ReturnMatrix[GoodLines, ] <- suppressWarnings(as.numeric(AQ)));
          } else {
            try(rownames(ReturnMatrix)[GoodLines] <- suppressWarnings(as.character(AQ[1])));
            try(ReturnMatrix[GoodLines,] <- suppressWarnings(as.numeric(AQ[2:length(AQ)])));
          }
        } else {
          AFilePrint(paste(\"  We found OnLine = \", OnLine, 
            \" to be faulty, read \", GoodLines, \": GoodLines.\", sep=\"\"));  
        }
    }    
    ";
}

secure.read.table.try <- function(AFile = NULL, nCols=0, header=FALSE, AGiveTries=0) {
    if (is.null(AFile)) {
      AFilePrint("secure.read.table.try: you supplied NULL AFile ");
    }
    ACon <- NULL;
    try(ACon <- file(AFile, "rt"), silent=TRUE);
    if (is.null(ACon)) {
      try(AFilePrint(paste("   secure.read: fail again (CON) tries = ", 
        AGiveTries, " file = ", AFile, sep="")), silent=TRUE);
      FailThisTime <- 1;
    } else {
      OnLine <- 0;
      try(eval(parse(text=secure.read.table.header.1())));
      GoodLines <- 0; NoLine <- FALSE;
      while(GoodLines == 0) {
        try(eval(parse(text=secure.read.table.first.line())));
        if (NoLine == TRUE) {
          AFilePrint(paste("secure.read: sorry, ", AFile, 
            " simply does not have a good line anywhere!", sep=""));
          try(close(ACon));
          return(NULL);
        }
      }
      while(NoLine == FALSE) {
        try(eval(parse(text=secure.read.table.a.line())));
      }
      try(close(ACon), silent=TRUE);
      OldGoodLines <- GoodLines;
      ReturnMatrix <- matrix(0,GoodLines, nCols);
      rownames(ReturnMatrix) <- c("", NROW(ReturnMatrix));
      ACon <- NULL; try(ACon <- file(AFile, "rt"), silent=TRUE);
      if (is.null(ACon)) {
        AFilePrint("Error secure.read.table: we couldn't open connection on second go around!");
        return(NULL);
      }
      OnLine <- 0;
      try(eval(parse(text=secure.read.table.header.2())));
      GoodLines <- 0; NoLine <- FALSE;
      while(NoLine == FALSE) {
        try(eval(parse(text=secure.read.table.a.line.fill())));
      }
      try(close(ACon));
      return(ReturnMatrix);
    }
  return(NULL);
}
TryToRecoverSummarySimulationList <- function(AContainDir, TCS = NULL,...) {
  eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("DefaultMaxLockTime")));
  eval(parse(text=GetG0Text("UniqueProcessIdentifier")));
  eval(parse(text=GetG0Text("CurrentLargeContainDir")));
  eval(parse(text=GetG0Text("CurrentSmallContainDir")));
  MySSDir <- paste(ALargeContainDir, "//SumPerformance", sep="");
  Oldwd <- getwd();
  try(setwd(CurrentSmallContainDir));
  while(TwoSimR5:::LockMeIn(verbose=0, quoteMore = "Set SumPerformance On Recover", LFDir = "SumPerformance",
    MaxLockTime = DefaultMaxLockTime, MoreLock = "", UniqueProcessIdentifier = UniqueProcessIdentifier) == FALSE) {
      for (ii in 1:100) {}
  } 
  AFilePrint("TryToRecoverSummarySimulation: We are locked in.");
  try(setwd(Oldwd));
  ISucceed <- 0;
  try(setwd(CurrentSmallContainDir));
  if (any(unlist(list.files("SumPerformance")) == "SummaryPerformanceRData.RData")) {
    try(ALLF <- unlist(list.files("SumPerformance")));
    try(load("SumPerformance/SummaryPerformanceRData.RData"));
    ISucceed <- 0;
    ARText <- paste("
    ARLen <- length(ALLF[substr(ALLF, 1, nchar(\"SMS\")) == \"SMS\"]);
    if (is.null(MySummarySimulationsList)) { ISucceed <- 0; 
    } else if (ARLen ==  length(MySummarySimulationsList$TheSummarySimulationList)) {
      ISucceed <- 1;
    }
    ", sep="");
    try(eval(parse(text=ARText)));
    try(setwd(Oldwd));
  }
  if (ISucceed == 1) {
    
    try(TwoSimR5:::UnLockMe(LFDir = paste(CurrentSmallContainDir, sep=""), 
        verbose=verbose, quoteMore="Final Saw no problem Sizeup"));
    try(TwoSimR5:::UnLockMe(LFDir = paste(CurrentSmallContainDir, "//", "SumPerformance", sep=""), 
       verbose=verbose, quoteMore="Final Saw no problem"));    
    AFilePrint("TryToRecover: Saw no problem.")
    return(TRUE);   
  }
  if (verbose >= 1) { 
    AFilePrint("TryToRecoverSummarySimulationList: Oh no recovering!")
  }
  MyNewSummarySimulationList <- SummarySimulationList$new(MySmallBaseDirectory = CurrentSmallContainDir,
    MyLargeBaseDirectory=CurrentLargeContainDir);
  if (verbose >= 1) {
    AFilePrint("TryToRecoverSummarySimulationList: Run Recover")
  }
  try(MyNewSummarySimulationList$Recover());
  if (verbose >= 1) {
    AFilePrint("TryToRecoverSummarySimulationList: Finished")
  }
  try(AFilePrint(paste("TryToRecover, Functon Recover Completed with GoodRecover",
    MySummarySimulationsList$GoodRecover, sep="")));
  if (MyNewSummarySimulationList$GoodRecover == FALSE) {
    try(AFilePrint(paste("TryToRecover: End where the Recover() failed!", sep="")));
  }
  setwd(CurrentSmallContainDir)
  MySummarySimulationsList = MyNewSummarySimulationList;
  
  try(secure.save.MSS(AObject=MySummarySimulationsList, 
    ObjectName="MySummarySimulationsList",
    file=paste("SumPerformance//SummaryPerformanceRData.RData", sep="")));
  try(TwoSimR5:::UnLockMe(LFDir = paste(CurrentSmallContainDir, sep=""), 
        verbose=verbose, quoteMore="Recover Successful"));
  try(TwoSimR5:::UnLockMe(LFDir = paste(CurrentSmallContainDir, "//", "SumPerformance", sep=""), 
       verbose=verbose, quoteMore="Recover Successful"));
  TC <- NULL; TK <- NULL;  VLD <- NULL;
  try(TC <- MyNewSummarySimulationList$GetTotalCompletedFunction());
  try(TK <- MyNewSummarySimulationList$MyTotalKilledFunctions);
  try(VLD <- MyNewSummarySimulationList$ValidSimulations)
  try(eval(parse(text=TryFillCurrentTotalCompletedTC())));
  try(AFilePrint(paste("TryToRecover, at end we have that GoodRecover = ",
    MySummarySimulationsList$GoodRecover, sep="")))
  return(TRUE);
}

RecoveryTSSLText <- function(MyContainDir = "CurrentSmallContainDir") {
  MyR <- paste("
  IFail <- 0;
  eval(parse(text=GetG0Text(\"CurrentLargeContainDir\")));
  eval(parse(text=GetG0Text(\"CurrentSmallContainDir\")));
  if (is.null(MySummarySimulationsList)) { IFail = 1; }
  MyListfiles <- unlist(list.files(paste(\"", CurrentLargeContainDir, "\",
    \"//SumPerformance\", sep=\"\")))
  if (length(MySummarySimulationsList$TheSummarySimulationList) == 0 &&
    any(substr(MyListfiles,1,nchar(\"SMS\")) == \"SMS\")) { 
    IFail = 1;
  }
  if (IFail == 1) {
    TryToRecoverSummarySimulationList(\"", CurrentSmallContainDir, "\")
    try(setwd(CurrentSmallContainDir));
    while(TwoSimR5:::LockMeIn(verbose=as.numeric(verbose), 
      quoteMore=quoteMore, 
      LFDir = \"SumPerformance\")==FALSE){
      for (jj in 1:100) {}}
    try(MySummarySimulationsList <- NULL);
    try(setwd(CurrentSmallContainDir));
    try(load(paste( 
      \"SumPerformance//SummaryPerformanceRData.RData\", sep=\"\")));
    if (is.null(MySummarySimulationsList)) {
      eval(parse(text=GetG0Text(\"UniqueProcessIdentifier\", 
      \"globalenv()\", S=1)))
      ErrorText <- paste(\"Hey we are process \",
        UniqueProcessIdentifier, \"
        and we tried to load in directory: \",
        MySummaryPerformanceDirectory, \"
        But didn't find anything at all!\", sep=\"\");
      try(setwd(CurrentSmallContainDir));
      try(TwoSimR5:::UnLockMe(LFDir = \"SumPerformance\", 
        verbose=verbose, quoteMore=\"Initial Sizeup\"));
      AFilePrint(ErrorText);
      AErrorPrint(ErrorText);
      return(-1);
    }
  }
  ", sep=""); 
}

SummarySimulationList$methods(
  Recover = function() {
    eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("DefaultMaxLockTime", "globalenv()", S=1)));
    try(.self$GoodRecover <- as.logical(FALSE));
    MyFF <- unlist(list.files(.self$MyBaseSimulationStorageDirectory))
    MyNames <- .self$NameFunctions;
    BeginOldwd <- getwd();
    MySMS <- MyFF[substr(MyFF,1,nchar("SMS")) == "SMS"];
    MyNamesSMS <- rep("", length(MySMS));
    if (verbose >= 1) {
      AFilePrint("RECOVERYRECOVERYRECOVERYRECOVERY")
      AFilePrint("We are Running Recover, hope we have good luck.")
    }

    OnP <- 0;
    for (ii in 1:length(MySMS)) {
      try(setwd(.self$SmallContainingDir));
      while(TwoSimR5:::LockMeIn(verbose=0, quoteMore = "On Recover", LFDir = "SumPerformance",
        MaxLockTime = DefaultMaxLockTime, MoreLock = "", UniqueProcessIdentifier = NULL) == FALSE) {
        Sys.sleep(1);
      } 
      try(setwd(.self$LargeContainingDir));
      load(paste("SumPerformance", "//", MySMS[ii], sep=""));  
      try(ABR <- SummarySimulationRecord$new(
          UniqueSimulationIdentifier=SMS$UniqueSimulationIdentifier,
          SavedUniqueSimulationRDataFile = paste("SumPerformance", 
          "//SMS", SMS$UniqueSimulationIdentifier, ".RData", sep="")));
      if (!is.null(ABR)) {
        OnP <- OnP+1;
        try(.self$TheSummarySimulationList[[OnP]] <- ABR);
        try(names(.self$TheSummarySimulationList)[OnP] <- SMS$UniqueSimulationIdentifier);
        MyNamesSMS[ii] <- .self$TheSummarySimulationList[[OnP]]$UniqueSimulationIdentifier;
      }
    }
    try(setwd(BeginOldwd));
    if (verbose >= 1) {
      AFilePrint(paste("We are recover: ", length(.self$TheSummarySimulationList), 
        " sims.", sep=""));
    }
    for (jj in 1:length(MyNames)) {
      try(setwd(.self$SmallContainingDir));
      while(TwoSimR5:::LockMeIn(verbose=0, quoteMore = "On Recover", LFDir = "SumPerformance",
        MaxLockTime = DefaultMaxLockTime, MoreLock = "", UniqueProcessIdentifier = NULL,
        NoSpaceLFDir=TRUE) == FALSE) {
          try(Sys.sleep(runif(1,0,4)));
      }
      try(setwd(.self$SmallContainingDir));  
      try(dir.create(MyNames[jj], showWarnings=FALSE, recursive=TRUE));
      while(TwoSimR5:::LockMeIn(verbose=0, quoteMore = "On Recover", LFDir = MyNames[jj],
        MaxLockTime = DefaultMaxLockTime, MoreLock = "", UniqueProcessIdentifier = NULL,
        NoSpaceLFDir=TRUE) == FALSE) {
          try(Sys.sleep(runif(1,0,4)));
      }  
      if (verbose >= 2) {
        AFilePrint(paste("Recover Function: ", MyNames[jj], 
          "", sep=""));
      }
      Got666 = 0;
      OnName <- MyNames[jj];
      MySSDir <-  paste(OnName, "//MyAttempts", sep="");
      try(dir.create(MySSDir, showWarnings=FALSE, recursive=TRUE));
      setwd(MySSDir);
      ISFiles <- unlist(list.files());
      if (length(ISFiles) >= 1) {
        ISFiles <- ISFiles[substr(ISFiles, nchar(ISFiles)-nchar(".txt")+1, nchar(ISFiles)) == ".txt"]
      }
      if (length(ISFiles) >= 1) {
        Got666 = 0;
        for (kk in 1:length(ISFiles)) {
           RDTry <- NULL;
           try(RDTry <- read.table(paste(ISFiles[kk], sep=""), sep=",", header=TRUE));
           if (is.null(RDTry) || NROW(RDTry) <= 0) {
            AT <- paste(" Uh Oh, Recover doesn't go well for jj = ", jj, ", kk = ", 
              kk,"
              I wonder why if ISFiles = (", paste(ISFiles, collapse=", "), ")
              We are on OnName = ", OnName, sep="");
             AFilePrint(AT);
             AErrorPrint(AT);
           }
           if (!is.null(RDTry)  && NROW(RDTry)>= 1) {
             UNP <- substr(ISFiles[kk], 1,nchar(ISFiles[kk]) - nchar(".txt"));
             for (tjj in 1:NROW(RDTry)) {
                ATS <- rownames(RDTry)[tjj];
                ATi <- (1:length(MyNamesSMS))[MyNamesSMS == ATS]
                ATi <- ATi[1];
                if (RDTry[tjj,5] != 666 && RDTry[tjj,5] > 0 && ATi[1] >= 1) {
                  try(.self$TheSummarySimulationList[[ATi]]$AddAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE)); 
                  try(.self$TheSummarySimulationList[[ATi]]$CompleteAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE)); 
                  Got666 = -1;
                } else if (RDTry[tjj,5] == 666 && RDTry[tjj,5] > 0 && ATi[1] >= 1) {
                  try(.self$TheSummarySimulationList[[ATi]]$AddAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE)); 
                  if (Got666== 0) { Got666 = 1;  ATiWin = ATi;  UNPWin=UNP; }
                } else if (!RDTry[tjj,5] %in% c(-666,-333) && RDTry[tjj,5] < 0 && ATi[1] >= 1) {
                  try(.self$TheSummarySimulationList[[ATi]]$AddAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE)); 
                  try(.self$TheSummarySimulationList[[ATi]]$FailAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE)); 
                } else if (RDTry[tjj,5] == -666 && RDTry[tjj,5] < 0 && ATi[1] >= 1) {
                  try(.self$TheSummarySimulationList[[ATi]]$AddAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE));
                  try(.self$TheSummarySimulationList[[ATi]]$FailAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE)); 
                  try(.self$TheSummarySimulationList[[ATi]]$KillAFunction(
                     ANameFunction=OnName,
                     SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE
                    ));
               } else if (RDTry[tjj,5] == -333 && RDTry[tjj,5] < 0 && ATi[1] >= 1) {
                  try(.self$TheSummarySimulationList[[ATi]]$AddAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE));
                  try(.self$TheSummarySimulationList[[ATi]]$FailAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE)); 
                } else if (RDTry[tjj,5] != 666 && RDTry[tjj,5] == 0 && ATi[1] >= 1) {
                  try(.self$TheSummarySimulationList[[ATi]]$AddAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNP, NoWrite=TRUE)); 
                }
             }
           }

        }
      }
      if (Got666 >= 1) {
              try(.self$TheSummarySimulationList[[ATiWin]]$CompleteAFunction(
                    ANameFunction = OnName, 
                    SSL = .self, UniqueProcessIdentifier=UNPWin, NoWrite=TRUE)); 
      }
      try(setwd(.self$SmallContainingDir));
      try(TwoSimR5:::UnLockMe(LFDir =  MyNames[jj], 
        verbose=verbose, quoteMore="Initial Sizeup", NoSpaceLFDir=TRUE));
      setwd(BeginOldwd);
    }
    try(.self$GoodRecover <- as.logical(TRUE));
    try(setwd(BeginOldwd));
    ##try(save(MySummarySimulationsList=.self, 
    ##  file=paste("SumPerformance", 
    ##  "//SummaryPerformanceRData.RData", sep="")));
    ##try(TwoSimR5:::UnLockMe(LFDir = paste(.self$ContainingDir, sep=""), 
    ##    verbose=verbose, quoteMore="Initial Sizeup"));
    ##try(TwoSimR5:::UnLockMe(LFDir = paste(.self$ContainingDir, "//", "SumPerformance", sep=""), 
    ##    verbose=verbose, quoteMore="Initial Sizeup"));
  }
)

secure.save.MSS <- function(AObject,ObjectName,file) {
  MaxTries <- 50;
  BackObject <- AObject;
  AGoodSave  <- FALSE;
  AFile <- file;
  while(MaxTries >= 0 && AGoodSave == FALSE)  {
    TryGet <- paste("AllGood <- FALSE;
      ", ObjectName, " <-  BackObject  
      AllGood <- TRUE;", sep="");
    try(eval(parse(text=TryGet)));
    if (AllGood == TRUE) {
      TryGet <- paste("AllGood <- FALSE;
      save(", ObjectName, "=", ObjectName, ",
        file=\"", AFile, "\");  
      AllGood <- TRUE;", sep=""); 
      try(eval(parse(text=TryGet)));
      if (AllGood == TRUE) {
        TryGet <- paste("AllGood <- FALSE; 
          ", ObjectName, " <- NULL;
         load(file=\"", AFile, "\"); 
         AllGood <- TRUE;
        ", sep="");
        try(eval(parse(text=TryGet)));
        if (AllGood == TRUE) {
           TryGet <- paste("AllGood <- FALSE;  AGoodSave <- FALSE;
             if (!is.null(", ObjectName, ")) {
                if (length(", ObjectName, "$TheSummarySimulationList) ==
                  length(BackObject$TheSummarySimulationList)) {
                  AGoodSave <- TRUE;  
                }
             }
             AllGood <- TRUE;
           ", sep="");
           try(eval(parse(text=TryGet)));
        }
      }
    }
    MaxTries <- MaxTries -1;
    if (AGoodSave == FALSE) {
      AFilePrint(paste("Hey: save.secure.MSS, we have a fail on MaxTries = ", 
        MaxTries, sep=""));
    }
  }  
  if (MaxTries >= 1) {
    return(1);
  } else {
    return(0);
  } 
}