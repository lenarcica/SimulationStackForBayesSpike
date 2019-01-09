################################################################################
##  2LassoSaveRoutines
##
##  (c) 2009-2019 by Alan Lenarcic
##    Written in lab of William Valdar, UNC Genetics
##    used to conduct simulations on UNC Killdevil LSF server to perform
##    simulations in many parameter settings and save material to
##    disk setup on server.
##
##
##  A "TableCompareSaver" object attempts to read information from the
##  save files directoires, which will be saved as RData (version 0.01 on)
##
##
##
##
##
##
##
##
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

try(library(methods, warn.conflicts = FALSE, quietly=TRUE), silent=TRUE);
try(library(R.methodsS3, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE); 
try(library(R.oo, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
DefaultNumPrint = 25;
##DefaultNumDetails = 9;
DefaultTableOutName = "TableOut.csv";
options(stringsAsFactors = FALSE)

###############################################################################
##  DefaultWriteFileInfoT:
##
##    Apparently writes a default.R file somewhere in an R Directory
DefaultWriteFileInfoT <- function(FileDirectory = "", 
  NameFile = "DefaultCategory.R", FileName ="",...) {
  if (NameFile=="DefaultCategory.R" && FileName != "") {
    NameFile <- FileName;
  }
  MyF = NULL;
  try(MyF <- file(description = paste(FileDirectory, "//", NameFile, sep=""), 
   open = "w", blocking = TRUE,
   encoding = getOption("encoding")));
  if (is.null(MyF)) {
    AFilePrint(paste("DefaultWriteFileInfoT did not open for NameFile=", NameFile,
      sep="")); 
    AFilePrint(paste("And FileDirectory = ", FileDirectory, sep=""));
    flush.console();
    return(-1);
  }
  writeLines(con=MyF, text=paste("DefaultNN = ", 1, ";", sep=""));  
  close(MyF);
}

###############################################################################
## HitASimulation
##
##   After a successful SMS simulation run on function NameFunction, this finds the directory
##    to save this UniqueSimulationIdentiifer setup and record the success.
##
HitASimulation <- function(UniqueProcessIdentifier, UniqueSimulationIdentifier, NameFunction, Hit = "+",
  quoteMore = "HitASimulation", verbose=0) {
  if (is.logical(verbose) && verbose == TRUE) {verbose = 1;}
  if (is.logical(verbose) && verbose == FALSE) {verbose =0;}
  if (verbose >= 2) {
    AFilePrint(paste("HitASimulation: Running for Name = ", NameFunction, " and UniqueProcess = ",
      UniqueProcessIdentifier, sep="")); flush.console();
    AFilePrint(paste(" --- with UniqueSimulationIdentifier = ", UniqueSimulationIdentifier, sep=""));
    flush.console();
  }
  if (Hit!="+") { Hit = "D" 
    if (verbose >= 2) {
      AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- OOF, Hit = ", Hit, " for NameFunction = ", NameFunction, sep=""));
      flush.console();
    }
  }
  eval(parse(text=GetG0Text("CurrentLargeContainDir", S=1)));
eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
  if (verbose >= 2) {
      AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- CurrentLargeContainDir = ", CurrentLargeContainDir, " for NameFunction = ", NameFunction, sep=""));
      AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- CurrentSmallContainDir = ", CurrentSmallContainDir, " for NameFunction = ", NameFunction, sep=""));
      flush.console();
  }
  MySummaryPerformanceDirectory = paste(CurrentLargeContainDir, "//", "SumPerformance", sep="");
  dir.create(MySummaryPerformanceDirectory, showWarnings=FALSE);
  if (verbose >= 3) {
      AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- Perform LockMeIn", sep=""));
      flush.console();
  }
  Oldwd <- getwd();
  setwd(CurrentSmallContainDir)
  while(LockMeIn(verbose=as.numeric(verbose), quoteMore=quoteMore, LFDir = "SumPerformance", LockIn="SumLock")==FALSE){}
  if (verbose >= 3) {
      AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- Completed LockMeIn. Get listfiles:", sep=""));
      flush.console();
  }
  setwd(Oldwd);
  MyListFiles = list.files(MySummaryPerformanceDirectory);
  if (verbose >= 3) {
      AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- MyListFiles receieved length ", length(MyListFiles), sep=""));
      flush.console();
  }
  if (verbose >= 3) {
      AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- Run GetSummaryPerformanceFile", sep=""));
      flush.console();
  }
  setwd(CurrentSmallContainDir);
  MyGG = GetSummaryPerformanceFile("SumPerformance", 
    verbose = max(verbose -1,0), quoteMore = quoteMore, AlreadyLocked=TRUE);
  if (verbose >= 3) {
      AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- Ran GetSummaryPerformanceFile, received MyGG", sep=""));
      flush.console();
  }
  DDR = length(MyGG$isDira);
  if (verbose >= 3) {
      try(AFilePrint(paste("----  HitASimulation[", UniqueSimulationIdentifier, 
        "] --- Got DDR = ",  paste(as.character(DDR), collapse=", "), sep="")), silent=TRUE);
      flush.console();
  }
  setwd(CurrentSmallContainDir)
  GetSmall <- unlist(list.files("SumPerformance"));
  if (is.numeric(CurrentLargeContainDir) && CurrentLargeContainDir[1] == 0) {
     CurrentLargeContainDir = "";
  }
  try(setwd(CurrentLargeContainDir));
  GetLarge <- unlist(list.files("SumPerformance"));
  MyListSMS = GetLarge[GetLarge!="SummaryPerformanceRData.RData" & 
    MyListFiles!="SummaryPerformance.txt" & MyListFiles != "LockFile.txt"];
  ##MyListSMS = MyListfiles[MyListfiles!="LockFile.txt"];
  MyListSMS = substr(MyListSMS, 4, nchar(MyListSMS)-6);
  if (verbose >= 3) {
    AFilePrint(paste("  --- we get MyListSMS = ", 
      paste(MyListSMS, collapse=", "), sep=""));  flush.console();
  }
  if (exists("UniqueProcessIdentifier") && is.list(UniqueProcessIdentifier)) {
    try(UniqueProcessIdentifier <- unlist(UniqueProcessIdentifier));
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1])
  }
  if (!exists("UniqueProcessIdentifier")) {
    AFilePrint("UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!");
    try(UniqueProcessIdentifier <- "");
  } else if (is.null(UniqueProcessIdentifier)) {
    AFilePrint("UhOh:HitASimulation: UniqueProcessIdentifier is NULL");
    try(UniqueProcessIdentifier <- "");  
  } else if (is.na(UniqueProcessIdentifier)) {
    AFilePrint("UhOh:HitASimulation: UniqueProcessIdentifier is NA");
    try(UniqueProcessIdentifier <- "");  
  } else if (is.numeric(UniqueProcessIdentifier)) {
    AFilePrint("UhOh:HitASimulation: UniqueProcessIdentifier is NUMERIC");
    try(UniqueProcessIdentifier <- "");  
  } else if (is.character(UniqueProcessIdentifier)) {
    if (UniqueProcessIdentifier[1] == "")  {
      AFilePrint("UhOh:HitASimulation: UniqueProcessIdentifier is Blank!");
    }
  } else {
    try(UniqueProcessIdentifier <- "");
  }

  if (verbose >= 2) {
    AFilePrint("  --- HitASimulation, looking for HMS. "); flush.console();
  }
  if (any(MyListSMS == UniqueSimulationIdentifier)) {
    ID = (1:length(MyGG$iDList))[MyGG$iDList == UniqueSimulationIdentifier];
    
    TTN = MyGG$isDira[ID,]
    TTN = as.vector(TTN);
    if (length(TTN) <= 0) {
      AFilePrint("HitASimulation: We have ERROR! No TTN length < 0"); flush.console();
    }
    MTA = rep("", length(TTN));
    for (ii in 1:length(TTN)) {
      AD = as.character(unlist(TTN[ii]));
      if (length(AD) > 1) { AD = AD[length(AD)] }
      AD = strsplit(AD, " ");  AD = paste(unlist(AD), collapse="");
      MTA[ii] = AD
    }
    if (any(substr(MTA,2,nchar(MTA)) == UniqueProcessIdentifier)) {
       eval(parse(text=GetG0Text("NameFunctions")));
       RTT = (1:length(NameFunctions))[NameFunctions == NameFunction];
       MyGG$isDira[ID,RTT] = paste(Hit, UniqueProcessIdentifier, sep=""); 
    }
    load(paste(CurrentLargeContainDir, "//SumPerformance//", "SMS", UniqueSimulationIdentifier, ".RData", sep=""));
    CurrentSolvesList = c(CurrentSolvesList, UniqueProcessIdentifier);
    save(SMS = SMS, CurrentSolvesList=CurrentSolvesList, 
      file=paste(MySummaryPerformanceDirectory, "//SMS", UniqueSimulationIdentifier,".RData", sep=""));
    ##load(paste(CurrentContainDir, "//SumPerformance//", 
    ##  "SummaryPerformanceRData.RData", sep=""));
    NewTable = cbind(MyGG$iDList, MyGG$isDira);
    colnames(NewTable) = c("iD", colnames(MyGG$isDira));
    MySummaryPerformanceTable = as.data.frame(NewTable, stringsAsFactors = FALSE);

    AText <- "
      try(setwd(CurrentSmallContainDir));
      try(setwd(\"SumPerformance\"));
      try(secure.save.MSS(AObject=MySummarySimulationsList, 
        ObjectName=\"MySummarySimulationsList\",
        file=paste(\"SummaryPerformanceRData.RData\", sep=\"\")));
      try(setwd(Oldwd));
        MySuccSave <- 1;
    ";
    try(eval(parse(text=AText)));
    ##try(write.table(NewTable, file = paste(CurrentContainDir, "//SumPerformance//", "SummaryPerformance.txt", sep=""),
    ##  append=FALSE, quote=FALSE, sep=", ", col.names=colnames(NewTable), row.names=FALSE)); 
    ##MySummaryPerformanceTable <- NewTable;
    ##try( save(MySummaryPerformanceTable = MySummaryPerformanceTable,
    ##  file = paste(CurrentContainDir, "//SumPerformance//", 
    ##  "SummaryPerformanceRData.RData", sep="")) );
  }
  if (verbose >= 3) {
    AFilePrint(paste("  --- HitASimulation finish, Unlock.", sep="")); 
    flush.console();
  }

  eval(parse(text=GetG0Text("CurrentLargeContainDir", "globalenv()", S=1)));
  try(setwd(CurrentSmallContainDir));
  UnLockMe(verbose=verbose, quoteMore=paste(quoteMore, " - HitASimulation", sep=""), LFDir = "SumPerformance",
    LockIn="SumUnLock");
  try(setwd(Oldwd));
  if (verbose >= 3) {
    AFilePrint(paste("  --- HitASimulation: AllFinish.", sep="")); flush.console();
  }
  return(1);
}

SaveCurrentNewSimulation <- function(SMS, verbose = FALSE, quoteMore = "SaveCurrentNewSim", 
  UniqueSimulationIdentifier,
  AlreadyLocked = TRUE, UniqueProcessIdentifier = UniqueProcessIdentifier, 
  ISample = 1) {
  if(is.null(verbose)) {
  print("Hey: Never Run Me, I'm old and don't work!"); flush.console();
  2 <- 1;
    verbose = 0;
  } else if (is.logical(verbose) && verbose ==FALSE) {verbose = 0;
  } else if (is.logical(verbose) && verbose == TRUE) {verbose = 1;
  } else { verbose = as.numeric(verbose); }
  eval(parse(text=GetG0Text("CurrentLargeContainDir", S=1)));
  eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
    
  Oldwd <- getwd();
  
  MySummaryPerformanceDirectory = paste(CurrentLargeContainDir, "//", "SumPerformance", sep="");
  dir.create(MySummaryPerformanceDirectory, showWarnings=FALSE);
  if (AlreadyLocked == FALSE) {
     setwd(CurrentSmallContainDir);
     while(LockMeIn(verbose=as.numeric(verbose), quoteMore=quoteMore, LFDir = "SumPerformance", LockIn="SumLock")==FALSE){
     }
  }
  CurrentSolvesList = NULL;
  if (!is.na(CurrentLargeContainDir) && !is.numeric(CurrentLargeContainDir) && CurrentLargeContainDir == "") {
    try(setwd(CurrentLargeContainDir));
    save(SMS = SMS, CurrentSolvesList=CurrentSolvesList, file=paste("SumPerformance", 
    "//SMS", UniqueSimulationIdentifier,".RData", sep=""));
  }
  
  MyGG = GetSummaryPerformanceFile("SumPerformance", 
    verbose = max(verbose -1,0), quoteMore = quoteMore, AlreadyLocked=TRUE);
  DDR <- NULL;
  try(DDR <- length(MyGG$isDira));
  if (length(DDR) > 0) {
  if (any(MyGG$iDList == UniqueSimulationIdentifier)) {
    AFilePrint(paste("Woah, I'm an idiot, looking through Summary Performance File, already 1 simulation described: ", UniqueSimulationIdentifier,
      sep="")); flush.console();
    AFilePrint("I can't continue this! ");
    AFileError(paste("Woah, I'm an idiot, looking through Summary Performance File, already 1 simulation described: ", UniqueSimulationIdentifier,
      sep="")); flush.console();
    tryCatch("Error Duplicate with this Simulation Identifier");
  }
  }
  eval(parse(text=GetG0Text("NameFunctions")));
  DDR2 = length(NameFunctions);
  ##if (is.matrix(MyGG$isDira)) { DDR = length(isDira[1,]); }

  MyGG2 = list( iDList = c(MyGG$iDList, UniqueSimulationIdentifier), isDira = rbind(MyGG$isDira, rep(0, DDR2)) );
  MySummaryPerformanceTable <- NULL;
  try(setwd(CurrentSmallContainDir));
  try(load(paste("SumPerformance",
    "//SummaryPerformanceRData.RData", sep="")));
  try(eval(parse(text=TwoSimR5:::RecoveryTSSLText(CurrentSmallContainDir))));
  
  NewP = matrix(c(UniqueSimulationIdentifier, rep(0, DDR2)),1,1+DDR2);
  NewP[1,ISample+1] = paste("-", UniqueProcessIdentifier, sep="");
  ANewTable <- rbind(MySummaryPerformanceTable, NewP);
  colnames(ANewTable) <- colnames(MySummaryPerformanceTable);
  try(setwd(CurrentSmallContainDir));
  try(secure.save.MSS(MySummaryPerformanceTable = MySummaryPerformanceTable, 
    file=paste("SumPerformance",
    "//SummaryPerformanceRData.RData", sep="")));
    
  ##write.table(NewP, file = paste(MySummaryPerformanceDirectory, "//SummaryPerformance.txt", sep=""),
  ##  append = TRUE, quote=FALSE, sep=",", col.names=FALSE, row.names=FALSE);
  if (!exists("MySummaryPerformanceTable") || is.null(MySummaryPerformanceTable)) {
    load(paste(MySummaryPerformanceDirectory,
    "//SummaryPerformanceRData.RData", sep=""));
    try(eval(parse(text=TwoSimR5:::RecoveryTSSLText(CurrentLargeContainDir))));
    if (FALSE) {
    if (!exists("MySummaryPerformanceTable") || is.null(MySummaryPerformanceTable)) {
      MySummaryPerformanceTable <- NewP;
      One = 0;
      ATryT <- "
      save(MySummaryPerformanceTable = MySummaryPerformanceTable,
        file=paste(MySummaryPerformanceDirectory,
        \"//SummaryPerformanceRData.RData\", sep=\"\"));
      One = 1;
      ";
      try(eval(parse(text=ATryT)));s
      if (One == 0) {
      AFilePrint("ERRORERRORERRORERRORERRORERRORERRORERRORERROR"); flush.console();
      AFilePrint("SaveNewSimulation FAIL On Empty MySummaryPerformanceTable Try Save "); flush.console();
      AFilePrint(paste("  Directory name: ", MySummaryPerformanceDirectory, sep=""));
      }  
    } else {
      One = 0;
      ATryT <- "
      MySummaryPerformanceTable <- rbind(MySummaryPerformanceTable, NewP);  
      save(MySummaryPerformanceTable = MySummaryPerformanceTable,
        file=paste(MySummaryPerformanceDirectory,
        \"//SummaryPerformanceRData.RData\", sep=\"\"));  
      One = 1; 
      ";
      try(eval(parse(text=ATryT)));
      if (One == 0) {
      AFilePrint("ERRORERRORERRORERRORERRORERRORERRORERRORERROR"); flush.console();
      AFilePrint("SaveNewSimulation FAIL On Second MySummaryPerformanceDirectory TrySave "); flush.console();
      AFilePrint(paste("  Directory name: ", MySummaryPerformanceDirectory, sep=""));
      }     
    }
    }
  } else { 
    One = 1;
    ATryT <- "
    try(MySummaryPerformanceTable <- rbind(MySummaryPerformanceTable, NewP));
    save(MySummaryPerformanceTable = MySummaryPerformanceTable, 
      file = paste(MySummaryPerformanceDirectory,
      \"//SummaryPerformanceRData.RData\", sep=\"\")
      );
    One = 1;
    ";
    try(eval(parse(text=ATryT)));
    if (One == 0) {
      AFilePrint("ERRORERRORERRORERRORERRORERRORERRORERRORERROR"); flush.console();
      AFilePrint("SaveNewSimulation FAIL On Regular First MySummaryPerformanceDirectory TrySave "); flush.console();
      AFilePrint(paste("  Directory name: ", MySummaryPerformanceDirectory, sep=""));
    }     
  }
  if (AlreadyLocked == FALSE) {
     UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory);
  }
}

GetSummaryPerformanceFile <- function(ALargeContainDir = "", ASmallContainDir="", verbose = TRUE, 
  quoteMore = "GetSumPerformanceFile", AlreadyLocked=FALSE, TCS = NULL, ...) {
  if (is.logical(verbose) && verbose == TRUE) {
    verbose = 1;
  } else if (is.logical(verbose)) {
    verbose = 0;
  } else {
    verbose = as.integer(verbose);
  }
  Oldwd <- getwd();
  if (verbose >= 3) {
    AFilePrint(paste("GetSummaryPerformanceFile: Load in ALargeContainDir = ",
      ALargeContainDir, sep="")); flush.console();
  }
  if (!exists("ALargeContainDir") || ALargeContainDir == "") {
    eval(parse(text=GetG0Text("CurrentLargeContainDir",S=1)));
    if (!(is.numeric(CurrentLargeContainDir) && CurrentLargeContainDir == 0)) {
      ALargeContainDir = CurrentLargeContainDir; 
    }
  }
  if (!exists("ASmallContainDir") || ASmallContainDir == "") {
    eval(parse(text=GetG0Text("CurrentSmallContainDir",S=1)));
    if (!(is.numeric(CurrentSmallContainDir) && CurrentSmallContainDir == 0)) {
      ASmallContainDir = CurrentSmallContainDir; 
    }
  }
  MySummaryPerformanceDirectory = paste(ALargeContainDir, "//", 
    "SumPerformance//", sep="");
  dir.create(MySummaryPerformanceDirectory, recursive=TRUE, showWarnings=FALSE);
  if (AlreadyLocked==FALSE) {while(LockMeIn(verbose=as.numeric(verbose), 
    quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory)==FALSE){}}
  try(setwd(ALargeContainDir));
  ASSubLarge <- unlist(list.files("SumPerformance"));
  try(setwd(ASmallContainDir));
  ASSubSmall <- unlist(list.files("SumPerformance"));
  if (any(substr(ASSubLarge, 1, nchar("SMS")) == "SMS") &&
    !(any(ASSubSmall == "SummaryPerformanceRData.RData"))) {
    TryToRecoverSummarySimulationList(ALargeContainDir, TCS = TCS)
    while(TwoSimR5:::LockMeIn(verbose=as.numeric(verbose), 
      quoteMore=quoteMore, 
      LFDir = MySummaryPerformanceDirectory)==FALSE){
    for (jj in 1:100) {}}
    try(MySummarySimulationsList <- NULL);
    try(setwd(ASmallContainDir));
    try(load(paste( 
      "SumPerformance/SummaryPerformanceRData.RData", sep="")));
    try(TwoSimR5:::UnLockMe(LFDir = MySummaryPerformanceDirectory, 
       verbose=verbose, quoteMore="Initial Sizeup"));
  } else if (any(ASSubSmall=="SummaryPerformanceRData.RData")) {
    try(setwd(ASmallContainDir));
    load(
      paste("SumPerformance", "//", "SummaryPerformanceRData.RData", sep=""));
    try(eval(parse(text=TwoSimR5:::RecoveryTSSLText(CurrentLargeContainDir))));
    if (!is.null(MySummaryPerformanceTable)) {
      MyRD <- MySummaryPerformanceTable;
    } else { MyRD <- NULL; }
  } else{ 
    MyRD <- NULL;
  }
  try(setwd(Oldwd));
  if (exists("MyRD") && !is.null(MyRD)) {
    IDRow = MyRD[,1, drop=FALSE];
    MyWorkRows = as.data.frame(MySummaryPerformanceTable[, 2:length(MySummaryPerformanceTable[1,])], stringsAsFactors=FALSE);
    IDReturnRows = (1:length(IDRow))[as.character(as.vector(IDRow)) != "-1"];
    if (length(IDReturnRows) == 0) {
      return(list(iDList = NULL, isDira = NULL, IDReturnRows = NULL)); 
    }
    if (AlreadyLocked==FALSE) {UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory)};
    isDira = MyWorkRows[IDReturnRows,, drop=FALSE];
    if (length(IDReturnRows) == 1) {
      MyN = names(isDira);
      for (iti in 1:length(isDira)) {
        isDira[iti] = paste(unlist(strsplit(as.character(isDira[iti]), " ")), collapse=""); 
      }
      names(isDira) = MyN;
    } else {
      MyN = colnames(isDira);
      for (tt in 1:length(IDReturnRows)) {
        for (iti in 1:length(MyN)) {
           isDira[tt,iti] = paste(unlist(strsplit(as.character(isDira[tt,iti]), " ")), collapse=""); 
        }
      }
      colnames(isDira) = MyN;
    }
    if (AlreadyLocked==FALSE) {UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory)};
    return(list(iDList = IDRow[IDReturnRows], isDira = isDira, IDReturnRows = IDReturnRows));    
  }
  MySummaryPerformanceTable = matrix(-1, 1, length(NameFunctions) + 1);
  colnames(MySummaryPerformanceTable) = c("iD", NameFunctions);
  MySuccSave <- 0;
  AText <- "
    Oldwd <- getwd();
    try(setwd(ASmallContainDir));
    try(secure.save.MSS(AObject=MySummarySimulationsList, 
      ObjectName=\"MySummarySimulationsList\",
      file=paste(\"SumPerformance//SummaryPerformanceRData.RData\", sep=\"\")));
    try(setwd(Oldwd));
  MySuccSave <- 1;
  ";
  try(eval(parse(text=AText)));
  if (MySuccSave == 0) {
    AFilePrint("GetSummaryPerformanceFile fales to save a new blank MySummaryPerformanceTable");
    AErrorPrint("GetSummaryPerformanceFile fales to save a new blank MySummaryPerformanceTable");
    tryCatch("Error in GetSummaryPerformanceFile");
  }
  
  if (AlreadyLocked==FALSE) {UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory)};
  return(list(iDList = NULL, isDira = NULL));
        
  if ((!exists("MySummaryPerformanceTable") || is.null(MySummaryPerformanceTable) ||
    !is.matrix(MySummaryPerformanceTable)) && any(unlist(list.files(MySummaryPerformanceDirectory)) 
    == "SummaryPerformance.txt")) {
    MyRD = NULL;
    try(MyRD <- read.table(paste(MySummaryPerformanceDirectory, "//",
      "SummaryPerformance.txt",sep=""), 
      header=TRUE,stringsAsFactors = FALSE, sep=","));
    try(save)
    if (is.null(MyRD)) {
      AFilePrint("ERRORERRORERRORERRORERRORERRORERRORERRORERROR"); flush.console();
      AFilePrint("GetSummaryPerformanceFile FAIL On Read "); flush.console();
      AFilePrint(paste("  Directory name: ", MySummaryPerformanceDirectory, sep=""));
      AFilePrint(paste("  FileName = SummaryPerformance.txt", sep=""));
      flush.console();
      MyRD = file(paste(MySummaryPerformanceDirectory, "//",
        "SummaryPerformance.txt",sep=""), "rt");
      if (is.null(MyRD)) {
        AFilePrint("Couldn't even read file using file()!"); flush.console();
        return(-1000);
      }
      ARTL = "-1";
      AFilePrint("FILE IS ---------------------------------"); flush.console();
      ARTL = readLines(MyRD, n=1);
      ALTLine = 1;
      while(!is.null(ARTL) &&  is.character(ARTL) && nchar(ARTL) > 0 && ARTL != "") {
        AFilePrint(paste("Line ", ALTLine, " : ", ARTL, sep="")); flush.console();
        ALTLine = ALTLine + 1;
        ARTL = readLines(MyRD, n=1);
      }
      close(MyRD);
      AFilePrint("-----------------------------------------"); flush.console();
      AFilePrint("  ---   See where the error mistake is? "); flush.console();
      return(-1000);
    } 
    IDRow = MyRD[,1, drop=FALSE];
    MyWorkRows = as.data.frame(MyRD[, 2:length(MyRD[1,])], stringsAsFactors=FALSE);
    IDReturnRows = (1:length(IDRow))[as.character(IDRow) != "-1"];
    if (length(IDReturnRows) == 0) {
      return(list(iDList = NULL, isDira = NULL, IDReturnRows = NULL)); 
    }
    if (AlreadyLocked==FALSE) {UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory)};
    isDira = MyWorkRows[IDReturnRows,, drop=FALSE];
    if (length(IDReturnRows) == 1) {
      MyN = names(isDira);
      for (iti in 1:length(isDira)) {
        isDira[iti] = paste(unlist(strsplit(as.character(isDira[iti]), " ")), collapse=""); 
      }
      names(isDira) = MyN;
    } else {
      MyN = colnames(isDira);
      for (tt in 1:length(IDReturnRows)) {
        for (iti in 1:length(MyN)) {
           isDira[tt,iti] = paste(unlist(strsplit(as.character(isDira[tt,iti]), " ")), collapse=""); 
        }
      }
      colnames(isDira) = MyN;
    }
   return(list(iDList = IDRow[IDReturnRows], isDira = isDira, IDReturnRows = IDReturnRows));
  } else {
    eval(parse(text=GetG0Text("NameFunctions")));
    if (length(NameFunctions) == 1 && is.numeric(NameFunctions) && NameFunctions == 0) {
      eval(parse(text=GetG0Text("NameFunctions", "TwoSim"))); 
    }
    MySummaryPerformanceTable = matrix(-1, 1, length(NameFunctions) + 1);
    colnames(MySummaryPerformanceTable) = c("iD", NameFunctions);
    write.table(MySummaryPerformanceTable, paste(MySummaryPerformanceDirectory, "//", "SummaryPerformance.txt", sep=""), col.names=TRUE,
      quote=FALSE, sep=",", row.names=FALSE, append=FALSE);
    try(setwd(ASmallContainDir));
    try(save(MySummaryPerformanceTable = MySummaryPerformanceTable, 
      file=paste("SumPerformance", "//", "SummaryPerformanceRData.RData", sep=""), 
      col.names=TRUE,
      quote=FALSE, sep=",", row.names=FALSE, append=FALSE));
    try(setwd(Oldwd));
    if (AlreadyLocked==FALSE) {UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MySummaryPerformanceDirectory)};  
    return(list(iDList = NULL, isDira = NULL));
  }
}


TableCompareSaverR5 = setRefClass("TableCompareSaverR5", 
  fields = list(
    ALargeContainDir = "character",
    ASmallContainDir = "character",
    LargeContainDir = function(XX = NULL) {
      if (is.null(XX)) { return(.self$ALargeContainDir); }
      try(.self$ALargeContainDir <- as.character(XX));
      return(.self$ALargeContainDir);
    },
    SmallContainDir = function(XX = NULL) {
      if (is.null(XX)) { return(.self$ASmallContainDir); }
      try(.self$ASmallContainDir <- as.character(XX));
      return(.self$ASmallContainDir);
    },
    ContainDir = function(XX = NULL) {
      if (is.null(XX)) { return(.self$ASmallContainDir); }
      try(.self$ASmallContainDir <- as.character(XX));
      return(.self$ASmallContainDir);
    },
    OnNameFunctions = "ANY", NameFunctions = function(XX = NULL) {
      if (is.null(XX)) { return(.self$OnNameFunctions); }
      try(.self$OnNameFunctions <- as.character(XX));
      return(.self$OnNameFunctions);
    },
    OnNameFunction = function(XX = NULL) {
       if (is.null(XX)) { return(.self$OnNameFunctions); }
      try(.self$OnNameFunctions <- as.character(XX));
      return(.self$OnNameFunctions);   
    },
    OnSimType = "character",
    FileCons = "ANY", FileLenS = "integer",
    FileMyDirectories="ANY",
    CurrentTotalCompleted="integer",
    ValidSimulations = "integer",
    ListTotalCompleted="list",
    FileUniqueSimulationIds = "ANY", FileUniqueProcessIds = "ANY",
    FileList = function() {
      return(list(FileCons = .self$FileCons,
        FileLenS = .self$FileLenS,
        FileMyDirectories = .self$FileMyDirectories,
        FileUniqueSimulationIds = .self$FileUniqueSimulationIds,
        FileUniqueProcessIds = .self$FileUniqueProcessIds       
        ));
    },
    experimentFunction = "ANY", verbose="integer",
    TargetLength = "integer", NumPrint="integer", NumDetails = "integer",
    jobii = "integer", TableOutName = "character", 
    MySummarySimulationList = "ANY",
    .RDataTableOutName = "character",
    RDataTableOutName = function() {
      OutNow <- "";
      MyTextToTry <- "
      if (nchar(.self$.RDataTableOutName) <= 0 || is.na(.self$.RDataTableOutName) ||  
        .self$.RDataTableOutName == \"NONAME\") {
        AS <- paste(.self$TableOutName, \".RData\", sep=\"\");
        if (is.na(.self$TableOutName) || is.na(nchar(.self$TableOutName)) ||
          .self$TableOutName == \"\" || nchar(.self$TableOutName) <= 3) {
          print(\"RDataTableOutName: False no TableOutName!\"); flush.console();
        } else if (substr(.self$TableOutName, nchar(.self$TableOutName)-3,
          nchar(.self$TableOutName)) %in% c(\".csv\", \".txt\")) {
          AS <-  paste(substr(.self$TableOutName,1,nchar(.self$TableOutName)-4), \".RData\", sep=\"\");
        } 
        try(.self$.RDataTableOutName <- AS);
      }
      ";
      try(eval(parse(text=MyTextToTry)));
      return(.self$.RDataTableOutName);
    },
    NotBS = "logical", ListOfParams = "ANY",
    UniqueProcessIdentifier = "character",
    OnUseFunctions = "ANY",
    WriteFileInfoT = "function",
    UseFunctions =  function(XX = NULL) {
      if (is.null(XX)) { return(.self$OnUseFunctions); }
      try(.self$OnUseFunctions <- XX);
      return(.self$UseFunctions);
    },
    OriginalOldwd = "character"
  ),
  methods = list(
    initialize = function(ALargeContainDir = "",
    ASmallContainDir = "",
  OnNameFunctions = NULL, 
  verbose = 0, 
  TargetLength = 800, 
  NumPrint = DefaultNumPrint, NumDetails = DefaultNumDetails,
  jobii = jobii, TableOutName = DefaultTableOutName, 
  WriteFileInfoT = DefaultWriteFileInfoT, NotBS = FALSE,
  ListOfParams = NULL, NInstance = 0, ...) {
  try(.self$OriginalOldwd <- as.character(getwd()));
  try(Sys.sleep(runif(1,0,4)));
  EvalPractice <- "
  if (!exists(\"OnNameFunctions\")) { OnNameFunctions = NULL; }
  if (!exists(\"verbose\")) { verbose = 0; }
  if (!exists(\"TargetLength\")) { 
    eval(parse(text=GetG0Text(\"TargetTotalSimulations\", \"globalenv()\", S=1)));
    if (TargetTotalSimulations >= 1) {
      TargetLength = TargetTotalSimulations;
    } 
  }
  if (!exists(\"NumPrint\")) { NumPrint = DefaultNumPrint; }
  if (!exists(\"NumDetails\")) { NumDetails = DefaultNumDetails; }
  if (!exists(\"TableOutName\")) { TableOutName = DefaultTableOutName; }
  if (!exists(\"WriteFileInfoT\")) { WriteFileInfoT = DefaultWriteFileInfoT; }
  if (!exists(\"NotBS\")) { NotBS = FALSE; }
  if (!exists(\"ListOfParams\")) { ListOfParams <- NULL; }
  ";
  try(eval(parse(text=EvalPractice)));
  BADUniqueProcessIdentifier <- 1;
  try(.self$TableOutName <- "NONAME");
  if (is.character(TableOutName) && nchar(TableOutName) >= 1) {
    try(.self$TableOutName <- TableOutName);
    try(.self$.RDataTableOutName <- "")
  }
  try(eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1))));
  try(.self$TargetLength <- as.integer(400));
  ATryGo <- "
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  }
  if (exists(\"UniqueProcessIdentifier\") && is.list(UniqueProcessIdentifier)) {
    try(UniqueProcessIdentifier <- unlist(UniqueProcessIdentifier));
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1])
  }
  if (length(UniqueProcessIdentifier) > 1) {
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1]);
  }
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  } else if (is.null(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NULL\");
    try(UniqueProcessIdentifier <- \"\"); 
  } else if (length(UniqueProcessIdentifier) > 1) {
    UniqueProcessIdentifier <- UniqueProcessIdentifier[1]; 
  } else if (is.na(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NA\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.numeric(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is NUMERIC\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.character(UniqueProcessIdentifier)) {
    if (UniqueProcessIdentifier[1] == \"\")  {
      AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier is Blank!\");
    } else {
      BADUniqueProcessIdentifier <- 0;
    }
  } else {
    try(UniqueProcessIdentifier <- \"\");
  }
  ";
  try(eval(parse(text=ATryGo)))
  if (BADUniqueProcessIdentifier == 1) {
     try(eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)))); 
     if (is.numeric(UniqueProcessIdentifier) && 
       UniqueProcessIdentifier == 0) {
       MyEvalText = "UniqueProcessIdentifier = 
         \"DefaultUniqueProcessIdentifier\";";
       try(eval(parse(text=MyEvalText)), silent=TRUE);
       try(.self$UniqueProcessIdentifier <- UniqueProcessIdentifier, silent=TRUE);
     }
  }
  if (!exists("NInstance") || NInstance == 0) {
    try(eval(parse(text=GetG0Text("NInstance", "globalenv()", S=1)))); 
  } 
  try(.self$OnSimType <- "DefaultSim");
  try(eval(parse(text=GetG0Text("LogitRegression", S=1))));
  try(.self$LogitRegression <- FALSE);
  eval(parse(text=GetG0Text("OnSimType", "globalenv()", S=1))); 
  if (is.null(OnSimType) || OnSimType == 0 || OnSimType %in% c("Default", "DefaultSim")) {
    try(.self$OnSimType <- "Default");
    .self$LogitRegression <- FALSE
  } else if (OnSimType == "Group")  {
    try(.self$LogitRegression <- FALSE);
    try(.self$OnSimType <- "Group");
  } else if (OnSimType == "Robit") {
    try(.self$LogitRegression <- TRUE);
    try(.self$OnSimType <- "Robit");
  } else if (OnSimType == "GroupRobit") {
    try(.self$LogitRegression <- TRUE);
    try(.self$OnSimType <- "GroupRobit");
  } else if (OnSimType == "GroupLogit") {
    try(.self$LogitRegression <- TRUE);
    try(.self$OnSimType <- "GroupLogit");
  }
  MyS <- " OnSimType <- .self$OnSimType; "
  try(eval(parse(text=MyS)));
  MyS <- " LogitRegression <- .self$LogitRegression; "
  try(eval(parse(text=MyS)));

  try(.self$CurrentTotalCompleted <- as.integer(0));

  try(.self$ListTotalCompleted <- list(a=as.integer(0)));
  
  try(.self$.RDataTableOutName <- "");
  try(.self$UniqueProcessIdentifier <- 
      as.character(DeclareUniqueProcessIdentifier(NInstance)));  
  AUPI <- "
    UniqueProcessIdentifier <- .self$UniqueProcessIdentifier;
    try(eval(parse(text=SetGText(\"UniqueProcessIdentifier\", \"globalenv()\",
      S=1))));
  "
  try(eval(parse(text=AUPI)));
  AT <- -1;
  MyTryText <- ".self$jobii <- as.integer(jobii);  AT = 1;"
  try(eval(parse(text=MyTryText)));
  if (AT == -1) {
    eval(parse(text=GetG0Text("jobii")));
    try(.self$jobii <- as.integer(jobii));
  }

  if (!exists("verbose")) { 
    try(.self$verbose <- as.integer(1), silent=TRUE); 
    MyT <- "verbose <- 1";
    try(eval(parse(text=MyT)), silent=TRUE);
  } else if (is.logical(verbose)) {
    if (verbose == TRUE) { try(.self$verbose <- as.integer(1), silent=TRUE); }
    if (verbose == FALSE){ try(.self$verbose <- as.integer(0), silent=TRUE); }
  }  else if (is.null(verbose)) { 
    try(.self$verbose <- as.integer(0), silent=TRUE); 
  } else {
    try(.self$verbose <- as.integer(verbose[1]), silent=TRUE);
  }
  if (length(.self$verbose) <= 0) {
    try(.self$verbose <- as.integer(0), silent=TRUE);
  }
  
  if (is.null(ALargeContainDir)  || !is.character(ALargeContainDir) ||
    nchar(ALargeContainDir) <= 0) {
    eval(parse(text=GetG0Text("ALargeContainDir", "globalenv()", S=1)));
  }
  try(.self$ALargeContainDir <- as.character(ALargeContainDir));
  if (is.null(ASmallContainDir)  || !is.character(ASmallContainDir) ||
    nchar(ASmallContainDir) <= 0) {
    eval(parse(text=GetG0Text("ASmallContainDir", "globalenv()", S=1)));
  }
  try(.self$ASmallContainDir <- as.character(ASmallContainDir));
  if (.self$verbose >= 1) {
    AFilePrint("TableCompareSaver: Initializing");
    AFilePrint(paste("  Note TableOutName = ", TableOutName, sep=""));
  }
  try(MyOld <- getwd());
  try(setwd(.self$ASmallContainDir));
  try(dir.create("SumPerformance", showWarnings=FALSE, recursive=TRUE));
  try(setwd(MyOld));
  if (.self$verbose >= 1 && NotBS == FALSE) {
   AFilePrint(paste("NotBS=FALSE; ALargeContainDir = ", .self$ALargeContainDir, sep=""));
   AFilePrint(paste("NotBS=FALSE; ASmallContainDir = ", .self$ASmallContainDir, sep=""));
   AFilePrint(paste("NotBS=FALSE; NameFunctions = ", .self$OnNameFunctions, sep=""));
   AFilePrint(paste("NotBS=FALSE;",  sep=""));
   AFilePrint(paste("NotBS=FALSE; TargetLength = ", TargetLength, sep=""));
   ##AFilePrint(paste("NotBS=FALSE; FileCons = ", aFileList$FileCons, sep=""));
   flush.console();
  }
  if (NotBS == FALSE) {
   if (.self$verbose >= 1) {
     AFilePrint(paste("Trying to return when NotBS = FALSE"));
     flush.console();
   }
   try(.self$NotBS <- as.logical(FALSE));
   try(setwd(.self$ASmallContainDir)); 
   try(UnLockMyProcess(verbose=0, quoteMore="UnLock SumPerformance on NoBS", LFDir = "SumPerformance",
      NoSpaceLFDir =FALSE, SumLock="No", LockIn="No"));
   try(setwd(.self$OriginalOldwd));
   return(.self);
  } else {
    try(.self$NotBS <- as.logical(TRUE));
  }
  ##if (is.null(experimentFunction) && NotBS == TRUE) {
  ##  AFilePrint("experimentFunction is Blank! quitting");
  ##  quit();
  ##}
   if (.self$verbose >= 3) {
     AFilePrint(paste("TwoSimR5 Create: Placing OnNameFunctions"));
     flush.console();
   }
  if (.self$OnSimType %in% c("Group", "GroupRobit")) {
    eval(parse(text=GetG0Text("UseFunctions")));
    eval(parse(text=GetG0Text("DefaultAllGroupFunctions")));
    eval(parse(text=GetG0Text("DefaultAllGroupNameFunctions")));
    eval(parse(text=GetG0Text("NameFunctions")));
    if (!any(NameFunctions %in% DefaultAllGroupNameFunctions)) {
      MyTT <- "
      try(NameFunctions <- DefaultAllGroupNameFunctions, silent=TRUE);
      try(UseFunctions <-  DefaultAllGroupFunctions, silent=TRUE);
      try(.self$NameFunctions <- NameFunctions);
      try(.self$UseFunctions <- DefaultAllGroupFunctions)  
      eval(parse(text=SetGText(\"NameFunctions\", \"globalenv()\", S=1)));
      eval(parse(text=SetGText(\"UseFunctions\", \"globalenv()\", S=1)));
      ";  
      try(eval(parse(text=MyTT)));
    } else {
      MyT <- "
       UseFunctions <- UseFunctions[NameFunctions %in% DefaultAllGroupNameFunctions];
       NameFunctions <- NameFunctions[NameFunctions %in% DefaultAllGroupNameFunctions];
       try(.self$NameFunctions <- NameFunctions);
       try(.self$UseFunctions <- UseFunctions);
      ";
      try(eval(parse(text=MyT)));
       eval(parse(text=SetGText("UseFunctions", "globalenv()", S=1)));
       eval(parse(text=SetGText("NameFunctions", "globalenv()", S=1)));
    }
  } else if (!exists("OnNameFunctions") || is.null(OnNameFunctions)) {
    eval(parse(text=GetG0Text("NameFunctions", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("UseFunctions", "globalenv()", S=1)));
    try(.self$OnUseFunctions <- UseFunctions);
    if (is.character(NameFunctions)) {
      try(.self$OnNameFunctions <- as.character(NameFunctions))
    } else {
      AFilePrint("Tried to find Name Functions but couldn't! no Char!");
      try(.self$OnNameFunctions <- as.character(""));
    }
    try(.self$OnUseFunctions <- UseFunctions);
  }  else {
    eval(parse(text=GetG0Text("NameFunctions", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("UseFunctions", "globalenv()", S=1)));
    eval(parse(text=SetGText("OnNameFunctions", "globalenv()", S=1)));
    MON <- NULL;
    try(MON <- match(OnNameFunctions, NameFunctions));
    if (is.null(MON)) {
      AFilePrint("TCS Initialize: we were unable to match OnNameFunctions and get MON!");
      tryCatch("We're going to not run past this and let you find the problem.");
      flush.console();
    }
    idM = (1:length(OnNameFunctions))[ !is.na(MON) ];
    try(.self$OnNameFunctions <- as.character(OnNameFunctions[idM]));
    try(.self$OnUseFunctions <- UseFunctions[MON[!is.na(MON)]]);
  }
  try(.self$ValidSimulations <- as.integer(1:length(.self$OnUseFunctions)));
  if (is.character(.self$ASmallContainDir)  && .self$ASmallContainDir == "") {
    AFilePrint(paste("Please Supply SmallContainDir!")); flush.console();
  }  
  if (is.character(.self$ALargeContainDir)  && .self$ALargeContainDir == "") {
    AFilePrint(paste("Please Supply LargeContainDir!")); flush.console();
  }
  CurrentSmallContainDir <- .self$ASmallContainDir;
  CurrentLargeContainDir <- .self$ALargeContainDir;
   if (.self$verbose >= 1) {
     AFilePrint(paste("TwoSimR5 Create: Check out FileCons"));
     flush.console();
   }
  if (.self$verbose >= 1) {
    AFilePrint(paste("TwoSimR5 create: To Try to Get SummaryPerformancList", sep=""));
    flush.console();
  }
  SummaryPerformanceDirectory <- paste(.self$ALargeContainDir, "//", "SumPerformance", sep="");
  dir.create(SummaryPerformanceDirectory, recursive=TRUE, showWarnings=FALSE);
  dir.create(paste(.self$ASmallContainDir, "//", "SumPerformance", sep=""),
    recursive=TRUE, showWarnings=FALSE);
  eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
  if (verbose >= 1) {
    AFilePrint(paste("New TableCompareSaver: about to conduct Lock in of LFDir = ", SummaryPerformanceDirectory, sep="")); 
  }
  try(Oldwd <- getwd());
  try(setwd(.self$ASmallContainDir));
  LockInText <- " NotGood = TRUE;  OutLock <- FALSE;  CountLocks <- 0;
  while(OutLock == FALSE && CountLocks <= 50) {
    OutLock <- LockMeIn(verbose=.self$verbose, quoteMore = \"Start TCS TableCompareSaver\", 
      LFDir = \"SumPerformance\",
      UniqueProcessIdentifier = UniqueProcessIdentifier, LockIn=\"SumLock\");
    if (OutLock == FALSE) { try(Sys.sleep(runif(1,0, 4))); }
    CountLocks <- CountLocks+1;
    if (CountLocks %% 10 == 0) {
      AFilePrint(paste(\" Hey: Lock Me in first SumPerformance \",
        \"not working.  on count = \", CountLocks, sep=\"\"));
    }
    if (CountLocks >= 50) {
      AFilePrint(paste(\"  Hey: Not working, too long to get SumPerformance directory\")); 
    }
  } 
    NotGood <- FALSE;
  "
  try(eval(parse(text=LockInText)));
  if (NotGood == TRUE || OutLock == FALSE) {
    AFilePrint("TCS Start ***********************************************************");
    AFilePrint(paste("TCS we went for CountLocks = ", CountLocks, sep=""));
    AFilePrint(paste("TCS: But we failed in our attempt to lock in, we should not continue."));
    AFilePrint(paste("TCS: We are now executing a deliberate error to stop this program", sep=""));
    ABD <- " BadOutLock <- NULL;  rm(BadOutLock);";
    eval(parse(text=ABD));
    BadOutLock <- BadOutLock+1;
  }
  
  try(setwd(Oldwd));
  try(setwd(.self$ALargeContainDir));
  ASSubLarge <- unlist(list.files("SumPerformance"));
  try(setwd(.self$ASmallContainDir));
  ASSubSmall <- unlist(list.files("SumPerformance"));
  try(setwd(Oldwd));
  if (any(substr(ASSubLarge, 1, nchar("SMS")) == "SMS") &&
    !(any(ASSubSmall == "SummaryPerformanceRData.RData"))) {
    TryToRecoverSummarySimulationList(CurrentLargeContainDir, TCS = .self)
    try(MySummarySimulationsList <- NULL);
    try(setwd(.self$ASmallContainDir));
    try(load(paste( 
      "SumPerformance//SummaryPerformanceRData.RData", sep="")));
    try(setwd(.self$ASmallContainDir) )
    try(TwoSimR5:::UnLockMe(LFDir = "SumPerformance", 
       verbose=verbose, quoteMore="Initial Sizeup after a recovery", LockIn="SumUnLock"));
    try(setwd(Oldwd));
  } else if (any(ASSubSmall=="SummaryPerformanceRData.RData")) {
    try(setwd(.self$ASmallContainDir), silent=TRUE);
    try(load(paste("SumPerformance", "//SummaryPerformanceRData.RData", sep="")));
    try(eval(parse(text=TwoSimR5:::RecoveryTSSLText(CurrentLargeContainDir))));
    if (!exists("MySummarySimulationsList")) {
      AFilePrint("***********************************************************");
      AFilePrint(paste("TCS Load in: Error: On load from directory ", SummaryPerformanceDirectory, sep=""));
      AFilePrint("TCS Load in: We tried to read SummaryPerformanceRData but no MySummaryPerformanceList");
      AErrorPrint("TCS Load in: Treat fail to load SummaryPerformanceRData as Early return error!");      
      try(setwd(.self$ASmallContainDir)); 
      try(UnLockMyProcess(verbose=0, quoteMore="UnLock SumPerformance on NoBS", LFDir = "SumPerformance",
        NoSpaceLFDir =FALSE, SumLock="No", LockIn="No"));
      try(setwd(.self$OriginalOldwd));
      return(.self);
    }
    try(.self$MySummarySimulationsList <- MySummarySimulationsList);
    try(eval(parse(text=TryFillCurrentTotalCompletedSELF())));
    MyTryUnlock <- " NotGood <- TRUE; 
      try(setwd(.self$ASmallContainDir));
      UnLockMe(verbose=.self$verbose, quoteMore=\"TCS Start\", LFDir = \"SumPerformance\",
        MoreLock = \"\",
      LockIn=\"SumUnLock\");
      try(setwd(Oldwd));
      NotGood <- FALSE;";
    try(eval(parse(text=MyTryUnlock)));
    if (NotGood) {
      try(AFilePrint("ASSub we tried to load in and then failed on the Unlock."));
    }
    try(setwd(Oldwd)); 
    if (.self$verbose >= 1) {
      AFilePrint("TCS new, we succesffuly loaded in old SummaryPerformanceRData!");
    }
  } else {
    if (.self$verbose >= 1) {
      AFilePrint("TCS new, we will need to create new Summary Simulation List");
    }
    EvalTry <-
      "MySummarySimulationsList <- SummarySimulationList$new(MySmallBaseDirectory=.self$ASmallContainDir,
        MyLargeBaseDirectory=.self$ALargeContainDir);"
    try(eval(parse(text=EvalTry)));
    try(.self$MySummarySimulationsList <- MySummarySimulationsList);
    try(setwd(.self$ASmallContainDir));
    try(secure.save.MSS(AObject=MySummarySimulationsList, 
      ObjectName="MySummarySimulationsList",
      file=paste("SumPerformance//SummaryPerformanceRData.RData", sep="")));
    try(setwd(.self$ASmallContainDir))
    try(UnLockMe(verbose=.self$verbose, quoteMore="TCS Start", LFDir = "SumPerformance",  
      LockIn="SumUnLock")); 
    try(setwd(Oldwd));
    if (.self$verbose >= 1) {
      AFilePrint("Successfuly created the new SummarySimulationsList!");
    } 
    try(eval(parse(text=TryFillCurrentTotalCompletedSELF())));
  }
  try(.self$ValidSimulations <- as.integer(MySummarySimulationsList$ValidSimulations))
  try(.self$CurrentTotalCompleted <- as.integer(MySummarySimulationsList$GetTotalCompletedFunction()));
  AText <- "
  try(CurrentTotalCompleted <- .self$CurrentTotalCompleted);  ";
  try(eval(parse(text=AText)));
  try(eval(parse(text=SetGText("CurrentTotalCompleted"))));
  ATText <- "
  try(ValidSimulations <- .self$ValidSimulations);
  ";
  try(eval(parse(text=ATText)));
  try(eval(parse(text=SetGText("ValidSimulations"))));
  if (.self$NotBS == TRUE) {
    if (FALSE) {
    aFileList = GiveFileCons(OnNameFunctions = .self$OnNameFunctions,
      verbose=verbose, ASmallContainDir = .self$ASmallContainDir,
      ALargeContainDir = .self$ALargeContainDir,
      RDataTableOutName = .self$RDataTableOutName,
      TableOutName = .self$TableOutName);
    }
    aFileList <- list(FileCons=NULL, FileLenS = NULL);
  } else {
    if (.self$verbose >= 1) {
      AFilePrint(paste("NotBS = ", NotBS, " we're giving a list")); flush.console();
    }
    aFileList = list(FileCons = NULL, FileLenS = NULL);
    
  }

  try(.self$FileCons <- aFileList$FileCons, silent=TRUE);
  try(.self$FileLenS <- as.integer(aFileList$FileLenS), silent=TRUE);
  try(.self$FileMyDirectories <- aFileList$FileMyDirectories, silent=TRUE);
  try(.self$FileUniqueSimulationIds <- aFileList$FileUniqueSimulationIds, silent=TRUE);
  try(.self$FileUniqueProcessIds <- aFileList$FileUniqueProcessIds, silent=TRUE);

  if (verbose == TRUE)  {
    AFilePrint("TableCompareSaver: Finished getting FileList"); flush.console();
  }
  if (is.null(OnNameFunctions)) {
    AFilePrint("I Cannot Use Null OnNameFunctions, Sorry ");
  }
  ##if (is.null(experimentFunction) && NotBS == TRUE) {
  ##  AFilePrint("You're Going to have to supply experimentFunction");
  ##  flush.console();
  ##} else {
  ##}
   if (.self$verbose >= 4) {
     AFilePrint(paste("TwoSimR5 Create: Check out AContainDir"));
     flush.console();
   }
  if (!is.character(ASmallContainDir)) {
    AFilePrint("Error: TableCompareSaverR5: AContainDir is not a string!");
      flush.console();
    try(.self$ASmallContainDir <- as.character(""), silent=TRUE);
    try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
    eval(parse(text=GetG0Text("DefaultTwoSimR5SmallDir", S=1)));
    if (is.numeric(DefaultTwoSimR5SmallDir) && DefaultTwoSimR5SmallDir == 0) {
      DefaultTwoSimR5SmallDir = "";
    }
    try(.self$ASmallContainDir <- as.character(DefaultTwoSimR5SmallDir));
  }  else {
    try(.self$ASmallContainDir <- as.character(ASmallContainDir));
  }
  if (.self$ASmallContainDir == "") {
    AFilePrint(paste("Hey, you supplied ASmallContainDir = ", ASmallContainDir, sep=""));
    AFilePrint(paste(" --  and we recorded it as ", .self$ASmallContainDir, sep=""));
    flush.console();
  }
  try(dir.create(.self$ASmallContainDir, showWarnings = FALSE, recursive = TRUE, mode = "0777"))

  if (!is.character(ALargeContainDir)) {
    AFilePrint("Error: TableCompareSaverR5: ALargeContainDir is not a string!");
      flush.console();
    try(.self$ALargeContainDir <- as.character(""), silent=TRUE);
    try(library(AlanDirectories, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
    eval(parse(text=GetG0Text("DefaultTwoSimR5LargeDir", S=1)));
    if (is.numeric(DefaultTwoSimR5LargeDir) && DefaultTwoSimR5LargeDir == 0) {
      ABText <- "
        ALargeContainDir <- \"\";
        .self$ALargeContainDir <- as.character(ALargeContainDir);
      "
      try(eval(parse(text=ABText)));
    }
    try(.self$ALargeContainDir <- as.character(DefaultTwoSimR5LargeDir));
  }  else {
    try(.self$ALargeContainDir <- as.character(ALargeContainDir));
  }
  if (.self$ALargeContainDir == "") {
    AFilePrint(paste("Hey, you supplied ALargeContainDir = ", ALargeContainDir, sep=""));
    AFilePrint(paste(" --  and we recorded it as ", .self$ALargeContainDir, sep=""));
    flush.console();
  }
  try(dir.create(.self$ALargeContainDir, showWarnings = FALSE, recursive = TRUE, mode = "0777"))
  
  if (is.function(WriteFileInfoT)) {
    try(.self$WriteFileInfoT <- as.function(WriteFileInfoT), silent=TRUE);
  } else {
    try(.self$WriteFileInfoT <- as.function(DefaultWriteFileInfoT), silent=TRUE);
  }

  if (.self$verbose >= 3) {
     AFilePrint(paste("TwoSimR5 Create: Do WritefileInfoT"));
     flush.console();
  }
  if (.self$NotBS == TRUE) {
    try(.self$WriteFileInfoT(
      FileDirectory=.self$AContainDir, NameFile="Category.R"));
  }

  if (.self$verbose >= 1) {
    AFilePrint(paste("MYE: About to Extend with NotBS = ", NotBS, sep=""));
    flush.console();
  }  
  if (.self$verbose >= 3) {
     AFilePrint(paste("TwoSimR5 Create: Setting TargetLength"));
     flush.console();
  } 
  TargetLengthText <- "
  if (is.numeric(TargetLength) && length(TargetLength) == 1 &&
    TargetLength[1] == 800) {
    eval(parse(text=GetG0Text(\"TargetTotalSimulations\")));
    if (is.numeric(TargetTotalSimulations) && length(TargetTotalSimulations) >= 1 &&
      TargetTotalSimulations >= 1) {
      try(TargetLength <-TargetTotalSimulations);
    }
    if (is.null(TargetLength) || !is.numeric(TargetLength) ||
      length(TargetLength) > 1 || TargetLength <= 0) {
      try(.self$TargetLength <- as.integer(800));  
    } else {
      try(.self$TargetLength <- as.integer(TargetLength));
    }  
  } else {
    try(.self$TargetLength <- as.integer(TargetLength));
  }
  ";
  try(eval(parse(text=TargetLengthText)));
  if (!exists("DefaultNumPrint")) {
    AFilePrint("Error: DefaultNumPrint Does not exists please fix defaults!");
  }
  ABD = -1;
  if (.self$verbose >= 3) {
     AFilePrint(paste("TwoSimR5 Create: Set NumPrint"));
     flush.console();
  }
  ElTryO <- "
  if (!exists(\"NumPrint\") || !is.numeric(NumPrint) || length(NumPrint) != 1 ||
    NumPrint <= 0) {
    eval(parse(text=GetG0Text(\"NumPrint\")));
    if (is.null(NumPrint) || !is.numeric(NumPrint) || length(NumPrint) != 1 ||
      NumPrint == 0) {
      AFilePrint(\"No Good NumPrint, couldn't get it!\"); flush.console();
      .self$NumPrint <- as.integer(5);
      ABD = 1;
    } else {
      .self$NumPrint <- as.integer(NumPrint);
      ABD = 1;
    } 
  } else {
     .self$NumPrint <- as.integer(NumPrint);
     ABD = 1;
  } ";
  try(eval(parse(text=ElTryO)));
  if (ABD == -1) {
    AFilePrint("Error: Tried to set NumPrint but a failure happened!!"); flush.console();
    try(.self$NumPrint <- as.integer(5))
  }
  ABD = -1;
  if (.self$verbose >= 1) {
     AFilePrint(paste("TwoSimR5: Set NumDetails"));
     flush.console();
  }
  SetDefaultScore("NumDetails", DefaultNumDetails);
  try(.self$NumDetails <- as.integer(NumDetails), silent=TRUE);  


  if (.self$verbose >= 1) {
     AFilePrint(paste("TwoSimR5:  Set NameFunctions.", sep=""));
     flush.console();
  }
  if (!is.character(OnNameFunctions)) {
    AFilePrint("Error: TableCompareSaverR5: onNameFunctions not char!"); 
      flush.console();
    try(.self$OnNameFunctions <- as.character(OnNameFunctions));
  } else {
    try(.self$OnNameFunctions <- as.character(OnNameFunctions));
  }
  if (!is.integer(aFileList$FileLenS)) {
     ##AFilePrint("Error: TableCompareSaverR5: FileLenS not integer!"); 
     ## flush.console();
     try(.self$FileLenS <- as.integer(aFileList$FileLenS));
  } else {
    try(.self$FileLenS <- as.integer(aFileList$FileLenS));
  }
  try(.self$experimentFunction <- NULL);
  if (exists("aFileList") && !is.null(aFileList) &&
    !is.null(aFileList) && is.character(aFileList$TableOutName)) {
    try(.self$TableOutName <- as.character(aFileList$TableOutName));  
  } else if (exists("TableOutName") && is.character(TableOutName)) {
    try(.self$TableOutName <- as.character(TableOutName));
  } else if (!exists("aFileList") || is.null(aFileList) || 
    is.null(aFileList$TableOutName) || !is.character(aFileList$TableOutName)) {
    AFilePrint("Error: TableCompareSaverR5: No good TableOutName!");
    flush.console();
    try(.self$TableOutName <- as.character("TableOutName"));
  } 
  try(.self$ListOfParams <- ListOfParams); 
  if (.self$verbose >= 3) {
    AFilePrint("TwoSimR5 Create: Setting in FileCons");
    flush.console();
  }
  if (FALSE) {
  if (!exists("aFileList") || is.null(aFileList)) {
    AFilePrint("Error: TableCompareSaverR5: aFileList was not generated!");
    flush.console();
    .self$FileCons <- NULL;
  } else if (is.null(aFileList$FileCons)) {
    AFilePrint("Error: TableCompareSaverR5: aFileList FileCons was NULL!");
    flush.console();
    .self$FileCons <- NULL;
  } else {
    try(.self$FileCons <- aFileList$FileCons);
  }
  }

  if (.self$verbose >= 1) {
    AFilePrint(paste("TwoSimR5 create: Made it to install end when NotBS = ", .self$NotBS, sep=""));
    flush.console();
  }
  try(setwd(.self$ASmallContainDir)); 
  try(UnLockMyProcess(verbose=0, quoteMore="UnLock SumPerformance on NoBS", LFDir = "SumPerformance",
        NoSpaceLFDir =FALSE, SumLock="No", LockIn="No"));
  try(setwd(.self$OriginalOldwd));  
  return(.self);
}
));

####################################################################
##  CheckUniqueProcessIdentifier();
##
## Attempts to get UniqueProcessIdentifier from the  globalenv for use
##
CheckUniqueProcessIdentifier <- function() {
  eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
  if (exists("UniqueProcessIdentifier") && is.list(UniqueProcessIdentifier)) {
    try(UniqueProcessIdentifier <- unlist(UniqueProcessIdentifier));
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1])
  }
  if (!exists("UniqueProcessIdentifier")) {
    AFilePrint("CheckUniqueProcessIdentifier, sorry, UniqueProcessIdentifier is not existant!");
  } else if (is.null(UniqueProcessIdentifier)) {
    AFilePrint("CheckUniqueProcessIdentifier, sorry, UniqueProcessIdentifier is NULL! ");
  } else if (length(UniqueProcessIdentifier) != 1) {
    AFilePrint(paste("CheckUniqueProcessIdentifier, sorry, UniqueProcessIdentifier has length: ",  
      length(UnqiueProcessIdentifier), sep=""));
  } else if (is.na(UniqueProcessIdentifier)) {
    AFilePrint("CheckUniqueProcessIdentifier, sorry, UniqueProcessIdentifier is NA! ");
  } else if (is.numeric(UniqueProcessIdentifier)) {
    AFilePrint("CheckUniqueProcessIdentifier, sorry, UniqueProcessIdentifier is numeric! ");
  } else if (is.character(UniqueProcessIdentifier) && UniqueProcessIdentifier == "") {
    AFilePrint("CheckUniqueProcessIdentifier, sorry, UniqueProcessIdentifier is Empty String! ");
  } else if (is.character(UniqueProcessIdentifier) && UniqueProcessIdentifier != "") {
    return(1);
  } else {
    try(AFilePrint(paste("CheckUniqueProcessIdentifier: I don't know what is going on with: ",
      UniqueProcessIdentifier, sep="")));
  }
  return(-1);
}


IntroduceToLockList <- function(quoteMore="", LFDir="", MoreLock="",
  UniqueProcessIdentifier = "WhatNone", NoSpaceLFDir = "",
  SumLock="No", Success="") {
  if ((is.character(UniqueProcessIdentifier) && UniqueProcessIdentifier == "WhatNone")  ||
   is.null(UniqueProcessIdentifier) || nchar(UniqueProcessIdentifier)==0){
    eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
  }
  eval(parse(text=GetG0Text("LockList", "globalenv()", S=1)));
  if (!is.list(LockList) || (is.numeric(LockList) && LockList[1] == 0)) {
    LockList <- list();
  } 
  eval(parse(text=GetG0Text("LastLockTime", "globalenv()",S=1)));
  PrevLockTime <- LastLockTime;
  LastLockTime <- proc.time();
  eval(parse(text=SetGText("LastLockTime", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("APrintDir", "globalenv()", S=1)));
  dir.create(paste(APrintDir, "//LockRData", sep=""), recursive=TRUE, showWarnings=FALSE);
  Oldwd <- getwd();
  if (length(LockList) >= 1) {
    for (ii in 1:length(LockList)) {
      AL <- LockList[[ii]]
      if (!is.null(AL$LFDir) && AL$LFDir == LFDir &&
          !is.null(AL$UniqueProcessIdentifier) && 
          AL$UniqueProcessIdentifier == UniqueProcessIdentifier) {
          try(LockList[[ii]]$Success<- Success);
          try(LockList[[ii]]$SumLock <- SumLock);
          eval(parse(text=SetGText("LockList", "globalenv()", S=1)));
          try(setwd(APrintDir));
          try(save(LockList=LockList, PrevLockTime=PrevLockTime,
            LastLockTime=LastLockTime,  UniqueProcessIdentifier=UniqueProcessIdentifier,
            file=paste("LockRData//", UniqueProcessIdentifier, ".RData",sep="")));
          try(setwd(Oldwd));
          return(2);
      }   
    }
  }
  New <- length(LockList)+1;
  LockList[[New]] <- list(Time=proc.time(), LFDir=LFDir, quoteMore=quoteMore,
    MoreLock=MoreLock,UniqueProcessIdentifier=UniqueProcessIdentifier, 
    NoSpaceLFDir=NoSpaceLFDir, SumLock=SumLock, Success=Success);
   eval(parse(text=SetGText("LockList", "globalenv()", S=1))); 
    try(setwd(APrintDir));
    try(save(LockList=LockList, PrevLockTime=PrevLockTime,
    LastLockTime=LastLockTime, UniqueProcessIdentifier=UniqueProcessIdentifier,
    file=paste("LockRData//", UniqueProcessIdentifier, ".RData",sep="")));
  try(setwd(Oldwd));  
  return(1);
}

RemoveFromLockList <- function(quoteMore="", LFDir="", MoreLock="",
  UniqueProcessIdentifier = "WhatNone") {
  if ((is.character(UniqueProcessIdentifier) && UniqueProcessIdentifier == "WhatNone")  ||
   is.null(UniqueProcessIdentifier) || nchar(UniqueProcessIdentifier)==0){
    eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
  }
  eval(parse(text=GetG0Text("LockList", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("LastLockTime", "globalenv()", S=1)));
  PrevLockTime <- LastLockTime;
   LastLockTime <- proc.time();
  eval(parse(text=SetGText("LastLockTime", "globalenv()", S=1)));
  Oldwd <- getwd();
  if (!is.list(LockList) || (is.numeric(LockList) && LockList[1] == 0)) {
    LockList <- list();
    eval(parse(text=SetGText("LockList")));
    return(0);
  }   
  if (is.list(LockList) && length(LockList) >= 1) {
     ID <- -1;
     for (ii in 1:length(LockList)) {
        AL <- LockList[[ii]]
        if (
          !is.null(AL$LFDir) && AL$LFDir == LFDir &&
          !is.null(AL$UniqueProcessIdentifier) && 
          AL$UniqueProcessIdentifier == UniqueProcessIdentifier) {
          ID <- ii;    
          break;
        }
     }
     if (ID <= 0) {
        return(0);
     }
     NewL <- list();
     KK <- 0;
     for (ii in 1:length(LockList)) {
        if (ii != ID) {
          KK <- KK+1;
          NewL[[KK]] <- LockList[[ii]];
        }
     }
     LockList <- NewL;
     if (length(NewL) == 0) {
        LockList <- list();
     }
     eval(parse(text=SetGText("LockList", "globalenv()", S=1)));
     eval(parse(text=GetG0Text("APrintDir", "globalenv()", S=1)));
     try(setwd(APrintDir));
     try(dir.create("LockRData", recursive=TRUE, showWarnings=FALSE));
     if (length(LockList) > 0) {
       try(save(LockList=LockList, PrevLockTime=PrevLockTime,
        LastLockTime=LastLockTime, UniqueProcessIdentifier=UniqueProcessIdentifier,
        file=paste("LockRData//", UniqueProcessIdentifier, ".RData",sep="")));
     } else {
       try(unlink(paste("LockRData//", UniqueProcessIdentifier, ".RData", 
         sep=""))); 
     }
     try(setwd(Oldwd));  
     
     return(1);
  }
  return(-1);
}

LockMeIn <- function(verbose=0, quoteMore = "", LFDir = DefaultAContainDir,
  MaxLockTime = DefaultMaxLockTime, MoreLock = "", UniqueProcessIdentifier = NULL,
  NoSpaceLFDir=FALSE, LockIn="No", SumLock="No")  {
  Oldwd <- getwd();
  if (!exists("MaxLockTime")) { MaxLockTime <- DefaultMaxLockTime}
  if (!exists("MoreLock")) { MoreLock <- ""; }
  if (!exists("MoreLockTime")) { MoreLockTime <- DefaultMaxLockTime; }
  if (!exists("quoteMore")) { quoteMore="";}
  if (!exists("LFDir")) { LFDir <- DefaultAContainDir; }
  if (!exists("NoSpaceLFDir")) { NoSpaceLFDir <- FALSE; }
  if (is.logical(verbose) && verbose == TRUE) {
    verbose = 1;
  } else if (is.logical(verbose) && verbose == FALSE) { verbose = 1; 
  } else if (is.null(verbose)) {  verbose = 0;
  } else { verbose = as.numeric(verbose); }
  
  if (!is.null(LockIn) && is.character(LockIn) && !(LockIn %in% c("No", "NO"))) {
    SumLock <- LockIn;
  }
  if (!is.null(SumLock) && is.character(SumLock) && !(SumLock %in% c("No", "NO"))) {
    LockIn <- SumLock;
  }
  if (SumLock %in% c("Summary", "SumPerformance") || LFDir == "SumPerformance") {
    SumLock = "SumLock";
  }
  if (is.logical(NoSpaceLFDir) && NoSpaceLFDir[1] == TRUE) {
    try(LFDir <- paste(unlist(strsplit(LFDir, " ")), collapse=""));
  }

  BADUniqueProcessIdentifier <- 1;
  ATryGo <- "
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh: LockMeIn: FirstLook: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  }
  if (length(UniqueProcessIdentifier) > 1) {
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1]);
  }
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  } else if (is.null(UniqueProcessIdentifier)) {
    ##AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier is NULL\");
    ##try(UniqueProcessIdentifier <- \"\"); 
    try(eval(parse(text=GetG0Text(\"UniqueProcessIdentifier\", \"globalenv()\", S=1))));
    if (!is.null(UniqueProcessIdentifier) && length(UniqueProcessIdentifier) == 1 &&
      !is.na(UniqueProcessIdentifier) && is.character(UniqueProcessIdentifier) &&
      UniqueProcessIdentifier != \"\") {
      BADUniqueProcessIdentifier <- 0;
    }
  } else if (length(UniqueProcessIdentifier) > 1) {
    UniqueProcessIdentifier <- UniqueProcessIdentifier[1]; 
  } else if (is.na(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier is NA\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.numeric(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier is NUMERIC\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.character(UniqueProcessIdentifier)) {
    if (UniqueProcessIdentifier[1] == \"\")  {
      AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier is Blank!\");
    } else {
      BADUniqueProcessIdentifier <- 0;
    }
  } else {
    ##try(UniqueProcessIdentifier <- \"\");
  }
  ";
  try(eval(parse(text=ATryGo)))
  if (is.null(BADUniqueProcessIdentifier) ||
    BADUniqueProcessIdentifier != 0) {
    eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
  }

  ATryGo <- "
    BADUniqueProcessIdentifier <- 1;
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh: LockMeIn: FirstLook: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
  }
  if (length(UniqueProcessIdentifier) > 1) {
    try(UniqueProcessIdentifier <- UniqueProcessIdentifier[1]);
  }
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier, second load, still does not exist!!\");
    UniqueProcessIdentifier <- 1;
  } else if (is.null(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier, afer GetG0Load, is NULL\");
    try(UniqueProcessIdentifier <- \"\"); 
  } else if (length(UniqueProcessIdentifier) > 1) {
    UniqueProcessIdentifier <- UniqueProcessIdentifier[1]; 
  } else if (is.na(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier, after GetG0Load, is NA\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.numeric(UniqueProcessIdentifier)) {
    AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier, after GetG0Load, is NUMERIC\");
    try(UniqueProcessIdentifier <- \"\");  
  } else if (is.character(UniqueProcessIdentifier)) {
    if (UniqueProcessIdentifier[1] == \"\")  {
      AFilePrint(\"UhOh: LockMeIn: UniqueProcessIdentifier, after GetG0Load, is Blank!\");
    } else {
      BADUniqueProcessIdentifier <- 0;
    }
  } else {
    ##try(UniqueProcessIdentifier <- \"\");
  }
  ";
  if (is.null(BADUniqueProcessIdentifier) ||
    BADUniqueProcessIdentifier != 0)  {
    try(eval(parse(text=ATryGo)));
    if (is.null(BADUniqueProcessIdentifier) ||
      BADUniqueProcessIdentifier != 0)  {
      AFilePrint(paste("UhOh: LockMeIn: On Second try UniqueProcessIdentifier is ",
        UniqueProcessIdentifier, ", bad goings on. ", sep=""));  flush.console();
      UniqueProcessIdentifier <- "DefaultBadLockFileProcess!";
      eval(parse(text=SetGGtext("UniqueProcessIdentifier", "globalenv()", S=1)));
    }
  }
  try(LockFile <- paste(LFDir, "//LockFile", MoreLock, ".txt", sep=""));
  if (verbose >= 1) {
    AFilePrint(paste("LockMeIn, to lock to lockfile ", LockFile, 
      " and MaxLockTime = ", MaxLockTime, ", ", quoteMore, sep=""), LockIn=SumLock);
  }
  try(Oldwd <- getwd());
  MyTryText <- paste("
    Tryer=TRUE;
    if (LFDir != \"\") { setwd(LFDir); }
    Tryer<-FALSE;", sep="");
  try(eval(parse(text=MyTryText)));
  if (Tryer) {
    eval(parse(text=SetGText("LFDir", "globalenv()", S=1)));
    AFilePrint("********************************************************", LockIn=SumLock);
    AFilePrint("** LockMeIn: Error on attempt to set to LFDir  ", LockIn=SumLock);
    AFilePrint(paste("** LockMeIn: Error, LFDir is ", LFDir, sep=""), LockIn=SumLock);
    try(AFilePrint(paste("** LockMeIn: Error, and OnDir is ", getwd(), sep=""), LockIn=SumLock));
    AFilePrint(paste("** LockMeIn: and the quoteMore is ", quoteMore, sep=""), LockIn=SumLock);
    AErrorPrint(paste("LockFileHey at beginning for list.files we tried to setwd to 
      ", LFDir, " but it failed miserably!", sep=""))
    try(IntroduceToLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
      UniqueProcessIdentifier = UniqueProcessIdentifier, NoSpaceLFDir = NoSpaceLFDir,
      SumLock=SumLock, Success="Horrible Set LFDir Error"));
    return(FALSE);
  }
  try(FilesInDir <- unlist(list.files()));
  try(setwd(Oldwd));
  try(ARealTime <- GetRealTime());
  
  MaybeOne <- NULL;
  TrySetLFDirText <- "
    AGood <- TRUE; 
    if (LFDir != \"\") {
      setwd(LFDir);
    }
    AGood <- FALSE;
  ";
  try(eval(parse(text=TrySetLFDirText)));
  if (AGood) {
     AFilePrint("****************************************************************", LockIn=SumLock);
     AFilePrint(paste("** LockMeIn - Error, tried to setwd to LFDir = ", LFDir, sep=""), LockIn=SumLock);
     try(eval(parse(text=SetGText("LFDir", "globalenv()", S=1))));
     try(eval(parse(text=SetGText("Oldwd", "globalenv()", S=1))));
     AFilePrint(paste("** LockMeIn - Error, Oldwd =  ", Oldwd, sep=""), LockIn=SumLock);
     AFilePrint(paste("** LockMeIn - Error, quoteMore = ", quoteMore, sep=""), LockIn=SumLock);
     AFilePrint(paste("** LockMeIn - Error, What are we doing now, lets commit an error!", sep=""), LockIn=SumLock);
     AErrorPrint(paste("LockMeIn-error, setting to LFDir for quoteMore = ", quoteMore, sep=""), LockIn=SumLock);
     ATT <- "1=2";  try(setwd(Oldwd));
     eval(parse(text=ATT));
  }
  TryToWriteToEmptyFiles <- "
  if (length(FilesInDir) == 0 || all(FilesInDir != paste(\"LockFile\", 
    MoreLock, \".txt\", sep=\"\"))) {
    AFile <- paste(\"LockFile\", MoreLock, \".txt\", sep=\"\")
    AMatrix <- rbind(c(ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays),
      c(ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays));
    try(TwoSimR5:::secure.write.table(file=AFile, AMatrix, col.names=NULL, 
      row.names=c(UniqueProcessIdentifier, UniqueProcessIdentifier),
      sep=\", \"));
    try(setwd(Oldwd));
    MaybeOne <- 2;
    try(IntroduceToLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
      UniqueProcessIdentifier = UniqueProcessIdentifier, NoSpaceLFDir = NoSpaceLFDir,
      SumLock=SumLock, Success=\"Write To Empty Files\"));
    return(TRUE);
  } else {
    MaybeOne <- 1;
  }
  ";
  try(eval(parse(text=TryToWriteToEmptyFiles)));
  try(setwd(Oldwd));
  if (!is.null(MaybeOne) && MaybeOne == 2) { return(TRUE); }
  if (is.null(MaybeOne)) {
    AFilePrint("LockMeIn:  Error on LockMe In!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", LockIn=SumLock)
    AFilePrint("LockMeIn: Hey some attempt to LockFile More lock and load in. The all match failed!", LockIn=SumLock);
    eval(parse(text=SetGText("Oldwd", "globalenv()", S=1)));
    try(AFilePrint(paste("LockMeIn: and Oldwd was ", Oldwd, sep=""), LockIn=SumLock));
    AFilePrint(paste("LockMeIn: Error: The LFDir is ", LFDir, sep=""), LockIn=SumLock);
    try(eval(parse(text=SetGText("LFDir", "globalenv()", S=1))));
    try(IntroduceToLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
      UniqueProcessIdentifier = UniqueProcessIdentifier, NoSpaceLFDir = NoSpaceLFDir,
      SumLock=SumLock, Success="LockFile Lock and Load In Fail"));
    AErrorPrint("LockMeIn: Hey some attempt to LockFile More lock and load in. The all match failed!", LockIn=SumLock);
    MyT <- "2=1";
    eval(parse(text=MyT));
    tryCatch("Try to fail out on Oldwd! for LockFile bad save.")
    return(FALSE);
  }

  if (length(FilesInDir) > 0 && any(FilesInDir == paste("LockFile", MoreLock, ".txt", sep=""))) {
    ARealTime <- GetRealTime();
    RDDefault = rbind(c(UniqueProcessIdentifier, ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays),
     c(UniqueProcessIdentifier, ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays));
    RD = NULL;
    ATrySetLFDir <- "
      AGood <- TRUE; 
      if (LFDir != \"\") { setwd(LFDir); }
      AGood <- FALSE;
    ";
    try(eval(parse(text=ATrySetLFDir)));
    if (AGood) {
        AFilePrint("********************************************************************", LockIn=SumLock)
        try(AFilePrint("** LockMeIn: Error Hey, Files in Dir, we tried to set directory for secure.read.table "), LockIn=SumLock);
        try(AFilePrint(paste("** LockMeIn: Error, LFDir, : ", LFDir, sep=""), LockIn=SumLock));
        try(AFilePrint(paste("** LockMeIn: Error, and oldwd is ", Oldwd, sep=""), LockIn=SumLock));
        try(AFilePrint(paste("** LockMeIn: I don't know how this will work", sep=""), LockIn=SumLock));
    }
    try(RD <- secure.read.table(file=paste("LockFile", MoreLock, ".txt", sep=""), 
      header=FALSE, nCols = 3, FurtherText = paste("LFDir=", LFDir, sep="")), silent=TRUE);
    try(setwd(Oldwd));
    if (is.null(RD)) {
      AFilePrint(paste("LockMeIn: Error, on ", quoteMore, pase=""), LockIn=SumLock);
      AFilePrint(paste("LockMeIn: E      Was Looking at : ", LFDir, sep=""), LockIn=SumLock);
      AFilePrint(paste("LockMeIn: Try to set to globalenv: "), LockIn=SumLock);
      try(eval(parse(text=SetGText("LFDir", "globalenv()", S=1))));
      try(eval(parse(text=SetGText("Oldwd", "globalenv()", S=1))));      
      AFilePrint(paste("LockMeIn: E     length FilesInDir = ", length(FilesInDir), sep=""), LockIn=SumLock);
      AFilePrint(paste("LockMeIn: E     Files in Dir are (", paste(FilesInDir, collapse=", "), ")", sep=""), LockIn=SumLock);
      AFilePrint(paste("LockMeIn: E     NumFile = ", (1:length(FilesInDir))[ FilesInDir == "LockFile.txt"],
        sep=""), LockIn=SumLock);
      try(setwd(LFDir));
      SecureText <- "";  SecureRead <- "";
      try(SecureText <- secure.read.text.lines(file=paste("LockFile", MoreLock, ".txt", sep=""), lines = 2,
        FurtherText = paste("LFDir=", LFDir, sep="") ))
      try(SecureRead <- secure.read.table(file=paste("LockFile", MoreLock, ".txt", sep=""), nCols = 3, header=FALSE,
        FurtherText = paste("LFDir=", LFDir, sep="") ));
      try(setwd(Oldwd));
      AFilePrint("LockMeIn: E     The Text of the file is ", LockIn=SumLock);
      try(AFilePrint(SecureText));
      eval(parse(text=SetGText("SecureText", "globalenv()", S=1)));
      eval(parse(text=SetGText("SecureRead", "globalenv()", S=1)));
      try(AFilePrint(paste("LockMeIn: E     The Table is 1: ", SecureRead[1,], sep="")));
      try(AFilePrint(paste("LockMeIn: E     The Table is 2: ", SecureRead[1,], sep="")));
      AErrorPrint(paste("LockMeIn: Read error, file exists but only read NULL!  We are printing over", quoteMore, pase=""));
      ####AFilePrint("For Temporary, we quit");
      ##AFilePrint("File Reads: ");
      ##AFilePrint(cat(file = LockFile));
      ##quit();
      ARealTime <- GetRealTime();
      AMatrix <- rbind(c(ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays),
        c(ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays));
      try(unlink(LockFile), silent=TRUE); 
      AFile <-  AFile <- paste("LockFile", MoreLock, ".txt", sep="");
      if (LFDir != "") {
        try(setwd(LFDir));
      }
      TwoSimR5:::secure.write.table(file=AFile, AMatrix=AMatrix, col.names=NULL, 
        row.names=c(UniqueProcessIdentifier, UniqueProcessIdentifier),
        sep=", "); 
      try(IntroduceToLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
      UniqueProcessIdentifier = UniqueProcessIdentifier, NoSpaceLFDir = NoSpaceLFDir,
      SumLock=SumLock, Success="RD Read Fail"));
      try(setwd(Oldwd))     
      return(TRUE);
    }
    if (as.character(RD[1,1]) == as.character(UniqueProcessIdentifier)) {
      ARealTime <- GetRealTime();
      AMatrix <- matrix(0,2,3);
      try(AMatrix[1,] <- suppressWarnings(as.numeric(RD[1,2:4])));
      try(AMatrix[2,] <- c(ARealTime$AllSeconds, ARealTime$JustSeconds, 
        ARealTime$JustDays));
      AFile <- paste("LockFile", MoreLock, ".txt", sep="")
      if (LFDir != "") { try(setwd(LFDir)); }
      rownames(AMatrix) <- c(UniqueProcessIdentifier, UniqueProcessIdentifier);
      try(TwoSimR5:::secure.write.table(AMatrix = AMatrix, file=AFile,
        col.names=NULL, rownames(AMatrix), 
        sep=", "), silent=TRUE);
      try(setwd(Oldwd));
      try(IntroduceToLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
       UniqueProcessIdentifier = UniqueProcessIdentifier, NoSpaceLFDir = NoSpaceLFDir,
       SumLock=SumLock, Success="Copy over previous lock in"));
      return(TRUE);      
    }
    ARealTime <- GetRealTime();
    DiffDiff <- (ARealTime$JustDays-as.numeric(RD[2,4])) * 24 * 60 * 60  +
      (ARealTime$JustSeconds - as.numeric(RD[2,3]));
    if (abs( DiffDiff ) < MaxLockTime) {
      if (verbose >= 2) {
        AFilePrint(paste("LockMeIn: File Exists apparently, undertime.", 
          quoteMore, sep=""), LockIn=SumLock);
      }
      AT <- 0;
      Sys.sleep(runif(1,0, abs(MaxLockTime-abs(DiffDiff)) * .25));
      return(FALSE);
    }
    AMatrix <- rbind(
      c(ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays), 
      c(ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays));
    try(row.names(AMatrix) <- c(UniqueProcessIdentifier, UniqueProcessIdentifier));
    AFilePrint(paste("LockMeIn: DiffDiff = ", DiffDiff, " but MaxLockTime = ", MaxLockTime, sep=""), LockIn=SumLock);
    AFilePrint(paste("LockMeIn: We are Overriding Lock!  :", LockFile, sep=""), LockIn=SumLock);
    try(AFilePrint(paste("LockMeIn: Current Dir is ", getwd(), sep=""), LockIn=SumLock));
    try(AFilePrint(paste("LockMeIn: Oldwd is ", Oldwd, sep=""), LockIn=SumLock));
    try(AFilePrint(paste("LockMeIn: Now try to setwd to LFDir = ", LFDir, sep=""), LockIn=SumLock));
    try(AFilePrint(paste("LockMeIn and MoreTxt is ", MoreLock, sep=""), LockIn=SumLock));
    TryToLockInText <- "
      NotGood <- TRUE;
      if (LFDir != \"\") { setwd(LFDir); 
      }
      NotGood <- FALSE";
    try(eval(parse(text=TryToLockInText)));
    if (NotGood) {
      AFilePrint("LockMeIn: Error, no we tried to override but did not setwd!", LockIn=SumLock);
      AFilePrint("LockMeIn Error: is doomed to fail.  ", LockIn=SumLock);
    }
    AFile <- paste("LockFile", MoreLock, ".txt", sep="")
    TryOverrideText <- "
      NotGood <- FALSE
      TwoSimR5:::secure.write.table(AMatrix = AMatrix, file=AFile,
        col.names=NULL, row.names=c(UniqueProcessIdentifier, UniqueProcessIdentifier), 
        sep=\", \");
      NotGood <- TRUE;
    ";
    try(eval(parse(text=TryOverrideText)));
    if (NotGood == FALSE) {
        AFilePrint("LockMeIn: Override, we did not successfully override lock in.", LockIn=SumLock);
    } else {
        try(AFilePrint("LockMeIn: Overrride, we just successfully override locked in.", LockIn=SumLock));
    }
    try(setwd(Oldwd));
    try(IntroduceToLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
      UniqueProcessIdentifier = UniqueProcessIdentifier, NoSpaceLFDir = NoSpaceLFDir,
      SumLock=SumLock, Success="Override go over time"));
    return(TRUE);
  }
  ARealTime <- GetRealTime();  
  MyTime <-  c(ARealTime$AllSeconds, ARealTime$JustSeconds, ARealTime$JustDays);
  AMatrix <- rbind(c( as.numeric(MyTime)), 
   c( as.numeric(MyTime)));
  if (verbose >= 3) {
    AFilePrint(paste("Wrote a new Lock File to: ", LockFile, sep=""),LockIn=SumLock);
    ##flush.console();
  }
  ATryGood <- "
    AGood <- TRUE;
    if (LFDir != \"\") { try(setwd(LFDir)); }
    AGood <- FALSE;
  ";
  try(eval(parse(text=ATryGood)));
  if (AGood) {
    AFilePrint("*********************************************************", LockIn=SumLock);
    try(AFilePrint(paste("** LockMeIn a setwd(LFDir) at end fails.  ", LFDir, sep=""),LockIn=SumLock));
    try(AFilePrint(paste("** LockMeIn our Oldwd was ", Oldwd, sep=""), LockIn=SumLock));
    try(eval(parse(text=SetGText("LFDir", "globalenv()", S=1))));
    try(eval(parse(text=SetGText("Oldwd", "globalenv()", S=1))));
  }
  AFile <- paste("LockFile", MoreLock, ".txt", sep="")
  try(TwoSimR5:::secure.write.table(AMatrix = AMatrix, file=AFile,
        col.names=NULL, row.names=c(UniqueProcessIdentifier, UniqueProcessIdentifier), 
        sep=", "), silent=TRUE);
  try(IntroduceToLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
      UniqueProcessIdentifier = UniqueProcessIdentifier, NoSpaceLFDir = NoSpaceLFDir,
      SumLock=SumLock, Success="LockIn New"));
  try(setwd(Oldwd));
  return(TRUE);
}
UnLockMe <- function(verbose=0, quoteMore="", LFDir = DefaultLFDir,
  NoSpaceLFDir =FALSE, MoreLock="", SumLock="No", LockIn="No") {
  if (!exists("NoSpaceLFDir")) { NoSpaceLFDir <- FALSE; }
  if (is.logical(verbose) && verbose == TRUE) {
    verbose = 1;
  } else if (is.logical(verbose) && verbose == FALSE) { verbose = 1; 
  } else if (is.null(verbose)) {  verbose = 0;
  } else { verbose = as.numeric(verbose); }
  eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
  if (is.logical(NoSpaceLFDir) && NoSpaceLFDir[1] == TRUE) {
    try(LFDir <- paste(unlist(strsplit(LFDir, " ")), collapse=""));
  }
  if (!is.null(LockIn) && is.character(LockIn) && !(LockIn %in% c("No", "NO"))) {
    SumLock <- LockIn;
  }
  if (!is.null(SumLock) && is.character(SumLock) && !(SumLock %in% c("No", "NO"))) {
    LockIn <- SumLock;
  }
  if (SumLock %in% c("SumLock", "SumUnLock", "UnLock") || LFDir == "SumPerformance") {
    SumLock <- "SumUnLock";
  }

  Oldwd <- getwd();
  TryText <- "
     NotGood <- TRUE;
     setwd(LFDir);
     NotGood <- FALSE;
  ";
  try(eval(parse(text=TryText)));
  if (NotGood == TRUE) {
    AFilePrint(paste(" LockMeIn: Error: Hey we tried to set(LFDir) = ", LFDir, sep=""), LockIn=SumLock);
    AFilePrint(paste(" LockMeIn: Error: Where we are is in ", getwd(), sep=""), LockIn=SumLock);
    AFilePrint(paste(" LockMeIn: Error: quoteMore = ", quoteMore, sep=""), LockIn=SumLock);
  }
  TryText <-"
     NotGood <- TRUE;
     try(LockFile <- paste(\"LockFile\", MoreLock, \".txt\", sep=\"\"));
     FilesInDir = unlist(list.files());
     if (length(FilesInDir) > 0 && any(FilesInDir == paste(\"LockFile\", MoreLock, \".txt\", sep=\"\"))) {
     if (verbose>= 2) {
       AFilePrint(paste(\"UnLockMe: Unlocking \", LockFile, sep=\"\"), LockIn=SumLock);
       AFilePrint(paste(\" ---   And : \", quoteMore, sep=\"\"), LockIn=SumLock);
     }
     try(unlink(LockFile));
     }
     setwd(Oldwd);
     NotGood <- FALSE";
  try(eval(parse(text=TryText)));
  if (NotGood) {
    AFilePrint(paste("Hey we did not succeed in Unlocking from setwd to \"", LFDir, "\"", sep=""),LockIn=SumLock);
    eval(parse(text=SetGText("LFDir", "globalenv()", S=1)));
    eval(parse(text=SetGText("Oldwd", "globalenv()", S=1)));
    AFilePrint(" Look for issue. ");
  } else {
    try(RemoveFromLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
      UniqueProcessIdentifier = UniqueProcessIdentifier));
    return(TRUE);
  }

  try(LockFile <- paste(LFDir, "//LockFile", MoreLock, ".txt", sep=""));
  FilesInDir = unlist(list.files(LFDir));
  if (length(FilesInDir) > 0 && any(FilesInDir == paste("LockFile", MoreLock, ".txt", sep=""))) {
    if (verbose>= 2) {
      AFilePrint(paste("UnLockMe: Unlocking ", LockFile, sep=""), LockIn=SumLock);
      AFilePrint(paste(" ---   And : ", quoteMore, sep=""), LockIn=SumLock);
    }
    try(unlink(LockFile));
  }
  try(RemoveFromLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
    UniqueProcessIdentifier = UniqueProcessIdentifier));
  return(TRUE);
}



UnLockMyProcess <- function(verbose=0, quoteMore="", LFDir = DefaultLFDir,
  NoSpaceLFDir =FALSE, MoreLock="", SumLock="No", LockIn="No") {
  AFilePrint(paste("UnLockMyProcess: Start fo LFDir = ", LFDir, sep=""));
  eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
  if (!exists("NoSpaceLFDir")) { NoSpaceLFDir <- FALSE; }
  if (is.logical(verbose) && verbose == TRUE) {
    verbose = 1;
  } else if (is.logical(verbose) && verbose == FALSE) { verbose = 1; 
  } else if (is.null(verbose)) {  verbose = 0;
  } else { verbose = as.numeric(verbose); }

  if (is.logical(NoSpaceLFDir) && NoSpaceLFDir[1] == TRUE) {
    try(LFDir <- paste(unlist(strsplit(LFDir, " ")), collapse=""));
  }
  if (!is.null(LockIn) && is.character(LockIn) && !(LockIn %in% c("No", "NO"))) {
    SumLock <- LockIn;
  }
  if (!is.null(SumLock) && is.character(SumLock) && !(SumLock %in% c("No", "NO"))) {
    LockIn <- SumLock;
  }
  if (SumLock %in% c("SumLock", "SumUnLock", "UnLock") || LFDir == "SumPerformance") {
    SumLock <- "SumUnLock";
  }

  Oldwd <- getwd();
  TryText <- "
     NotGood <- TRUE;
     setwd(LFDir);
     NotGood <- FALSE;
  ";
  try(eval(parse(text=TryText)));
  if (NotGood == TRUE) {
    AFilePrint(paste(" UnLockMyProcess: Error: Hey we tried to set(LFDir) = ", LFDir, sep=""), LockIn=SumLock);
    AFilePrint(paste(" UnLockMyProcess: Error: Where we are is in ", getwd(), sep=""), LockIn=SumLock);
    AFilePrint(paste(" UnLockMyProcess: Error: quoteMore = ", quoteMore, sep=""), LockIn=SumLock);
  }
  TryText <-"
     NotGood <- TRUE;
     try(LockFile <- paste(\"LockFile\", MoreLock, \".txt\", sep=\"\"));
     FilesInDir = unlist(list.files());
     if (length(FilesInDir) > 0 && any(FilesInDir == paste(\"LockFile\", MoreLock, \".txt\", sep=\"\"))) {
     if (verbose>= 2) {
       AFilePrint(paste(\"UnLockMyProcess: Looking To Unlock \", LockFile, sep=\"\"), LockIn=SumLock);
       AFilePrint(paste(\" ---   And : \", quoteMore, sep=\"\"), LockIn=SumLock);
     }
     try(ASS <- secure.read.table(file=paste(\"LockFile\", MoreLock, \".txt\", sep=\"\"), 
      header=FALSE, nCols = 3, FurtherText = paste(\"UnLockMyProcess: LFDir=\", LFDir, sep=\"\")), silent=TRUE);
        if (is.null(ASS)) {
        } else if (NCOL(ASS) == 4 && as.character(ASS[1,1]) == as.character(UniqueProcessIdentifier)) {
          try(unlink(LockFile));  
        } else if (NCOL(ASS)==3 && rownames(ASS)[1] !=  UniqueProcessIdentifier) {     
        } else {
          try(unlink(LockFile));  
        }
     }
     setwd(Oldwd);
     NotGood <- FALSE";
  try(eval(parse(text=TryText)));
  if (NotGood) {
    AFilePrint(paste("Hey we did not succeed in Unlock MyProcess from setwd to \"", LFDir, "\"", sep=""),LockIn=SumLock);
    eval(parse(text=SetGText("LFDir", "globalenv()", S=1)));
    eval(parse(text=SetGText("Oldwd", "globalenv()", S=1)));
    AFilePrint(" Look for issue. ");
  } else {
    try(RemoveFromLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
      UniqueProcessIdentifier = UniqueProcessIdentifier));
    try(setwd(Oldwd));
    return(TRUE);
  }
  try(setwd(Oldwd));
  try(RemoveFromLockList(quoteMore=quoteMore, LFDir=LFDir, MoreLock=MoreLock,
    UniqueProcessIdentifier = UniqueProcessIdentifier));
  return(TRUE);
}


PRString <- function() {
  paste("Onjobii = ", jobii, "; n = ", n, "; k = ", k, "; p = ", p,
    "; sigma = ", sigma, sep="");
}

ReadOneFile <- function(OnNameFunction = NULL, AContainDir = NULL, TableOutName = NULL, 
  RDataTableOutName = NULL, verbose = FALSE,
  AlreadyLocked = FALSE) {
  try(eval(parse(text=GetG0Text("MyColNames", "globalenv()", S=1))));

if (is.logical(verbose) && verbose == TRUE) {
  verbose = 1;
} else if (is.logical(verbose) && verbose == FALSE) {
  verbose = 0;
} else {
  verbose = as.integer(verbose);
}
if (verbose >= 2) {
 AFilePrint(paste("ReadOneFile: About to Start", sep="")); flush.console();
}
if (verbose >= 2) {
  AFilePrint(paste("ReadOneFile: TableOutName = ", RDataTableOutName, sep=""));
  flush.console();
}
AD = AContainDir;
try(Oldwd <- getwd());
if (AlreadyLocked == FALSE) {
  if (verbose >= 2) {
    AFilePrint(paste("ReadOneFile: about to LockIn LFDir = ",  AContainDir, sep=""));
    flush.console();
  }
  ATryLock <- "
  DidLock <- FALSE; 
  setwd(AContainDir);
  while(LockMeIn(verbose=verbose, 
    quoteMore=\"ReadOneFile\", LFDir = OnNameFunction) ==FALSE ) { 
    DidUnlock <- FALSE; Sys.sleep(runif(1,0,4));
  }
  DidLock <- TRUE;
  setwd(Oldwd);
  ";
  try(eval(parse(text=ATryLock)));
  if (DidLock == FALSE) {
    AFilePrint("ReadOneFile: we tried to unlock but we got a fail!");
    AErrorPrint("ReadOneFile: Fail on try to lock directory on read one file!");
  }
} 
  ii = (1:length(NameFunctions))[NameFunctions == OnNameFunction];
  if (verbose == TRUE) {
    AFilePrint(paste("GiveAFileCons: Starting and looking into OnNameFunctions ii = ", 
      ii, " is ", OnNameFunction, sep="")); flush.console();
  }
    
  OnDir <- paste(AContainDir, "//", OnNameFunction, sep="");
  if (verbose >= 3) {
    AFilePrint(paste("ReadOneFile: Create OnDir =  ",  OnDir, sep=""));
    flush.console();
  }
  dir.create(OnDir, showWarnings = FALSE)
  FilesInDir = unlist(list.files(OnDir));
  if (length(FilesInDir) >= 1 && any(FilesInDir == RDataTableOutName)) {
    if (verbose >= 4) {
      AFilePrint(paste("ReadOneFile: Got it in FilesInDir, About to Read TableOutName = ", TableOutName, sep=""));
      flush.console();
    }
    RDD = NULL;
    AB <- warnings();
    ##try(RDD <- read.table( paste(OnDir, "//", TableOutName, sep=""), sep=",",
    ##  header=TRUE, na.strings="NA", fill = TRUE, flush=TRUE,
    ##  stringsAsFactors = FALSE, col.names=MyColNames));
    try(setwd(AContainDir));
    try(load(paste(OnNameFunction, "//", RDataTableOutName, sep="")));
    if (AlreadyLocked == FALSE) { 
      try(setwd(AContainDir));
      UnLockMe(LFDir = OnNameFunction, verbose=verbose, quoteMore="Initial Sizeup");
      try(setwd(Oldwd));
      AlreadyLocked <- TRUE;
    }
    try(setwd(Oldwd));

    
    AW <- warnings();
    if (length(AW)-length(AB) >= 1 || length(MyColNames) != NCOL(RDD)) {
      AFilePrint("ISSUE -----------------------------------------------------");
      AFilePrint(" ----- ");
      AFilePrint(paste("length(MyColNames) = ", length(MyColNames), sep=""));
      flush.console();
      AFilePrint(paste("NCOL(RDD) = ", NCOL(RDD), sep=""));
      AFilePrint(" ----- "); flush.console();
      AFilePrint(paste("  MyColNames = (", paste("\"", MyColNames, "\"", collapse=", "),
        ")", sep="")); flush.console();
      AFilePrint(" ---- ");
      AFilePrint(" RDD  = "); flush.console();
      if (NROW(RDD) < 5) {
        AFilePrint(RDD[1:NROW(RDD),]); flush.console();
      } else {
        AFilePrint(RDD); flush.console();
      }
      AFilePrint("The AW warnings are: "); flush.console();
      AFilePrint(AW);
      AFilePrint("  We'll figure out RDD problem next. "); flush.console();
    }
    if (verbose >= 4) {
      AFilePrint(paste("ReadOneFile: We Read RDD, checking dimension.", sep=""));
      flush.console();
    }
    if (!is.null(RDD) && NCOL(RDD) == length(MyColNames)) {  
      FileSLen = NROW(RDD);
      FileCon = RDD;  
      try(FileMyDirectory <- MDir, silent=TRUE);
      try(FileMySimulationIds <- NULL);
      try(FileMySimulationIds <- MyPiDs, silent=TRUE);
      try(FileUniqueProcessIds <- NULL);
      try(FileUniqueProcessIds <- MyUniqueProcesses, silent=TRUE);    
    } else if (is.null(RDD) || dim(RDD)[2] != length(MyColNames)) {
      if (AlreadyLocked == FALSE) {
        try(setwd(AContainDir));
        try(UnLockMe(LFDir = OnNameFunction, 
        verbose=verbose, quoteMore="ReadOneFile got a null dim"));
        try(setwd(Oldwd));
        AlreadyLocked=TRUE;
      }
      if (is.null(RDD)) {
        AFilePrint("ReadOneFile: We Have a Multi-line Problem we must correct");
        AErrorPrint(paste("ReadOneFile: Hey NULL, RDD.", sep=""));
        tryCatch("ReadOneFile: I don't want to continue with defective RDD");
      } else if (is.null(dim(RDD))) {
        AFilePrint("ReadOneFile: dim(RDD) is NULL.");
        AErrorPrint(paste("ReadOneFile: Hey NULL, RDD.", sep=""));        
        tryCatch("ReadOneFile: I don't want to continue with defective RDD");
      } else if (dim(RDD)[2] != length(MyColNames)) {
        AFilePrint(paste("ReadOneFile: dim(RDD) is bad (", 
          paste(dim(RDD), collapse=", "), ")", sep=""));
        AErrorPrint(paste("ReadOneFile: Hey NULL, RDD. (",
          paste(dim(RDD), collapse=", "), ")", sep=""));  
        tryCatch("ReadOneFile: I don't want to continue with defective RDD");        
      } else {
        AFilePrint(paste("ReadOneFile: dim(RDD) is bad (", 
          paste(dim(RDD), collapse=", "), ")", sep=""));
        AErrorPrint(paste("ReadOneFile: Hey NULL, RDD. (",
          paste(dim(RDD), collapse=", "), ")", sep=""));  
        tryCatch("ReadOneFile: I don't want to continue with defective RDD");            
      }
      AFilePrint(paste("Frankly ReadOneFile: I'm not willing to continue with defective RDD for ", OnNameFunction, sep=""));
      tryCatch("ReadOneFile: Can't contine from here!");
      RDD2 = matrix(as.numeric(RDD[1:length(RDD[,1]), 1:(length(MyColNames))]),
        length(RDD[,1]), NumDetails);
      ##write.table(RDD2, paste(OnDir, "//", TableOutName, sep=""), 
      ##  sep=",", eol="\n", na="NA", append = FALSE, 
      ##  row.names=FALSE, col.names=MyColNames);
      RDD = RDD2;
      FileSLen = dim(RDD)[1];
      FileCon = RDD;
      try(FileMyDirectory <- MDir, silent=TRUE);
      try(FileMySimulationIds <- NULL);
      try(FileMySimulationIds <- MyPiDs, silent=TRUE);
      try(FileUniqueProcessIds <- NULL);
      try(FileUniqueProcessIds <- MyUniqueProcesses, silent=TRUE);
    } else {
      FileSLen = 0; FileCon = NULL;
      try(FileMyDirectory <- MDir, silent=TRUE);
      try(FileMySimulationIds <- NULL);
      try(FileMySimulationIds <- MyPiDs, silent=TRUE);
      try(FileUniqueProcessIds <- NULL);
      try(FileUniqueProcessIds <- MyUniqueProcesses, silent=TRUE);
    }
    if (verbose >= 4) {
      AFilePrint(paste("ReadOneFile: After Read In FileSLen = ", 
        paste(FileSLen, collapse=", "), ".", sep=""));
      flush.console();
    }
  } else {
    FileSLen = 0;  FileCon = NULL;
    try(FileMyDirectory <- "");
    try(FileMyDirectory <- MDir, silent=TRUE);
    try(FileMySimulationIds <- NULL);
    try(FileMySimulationIds <- MyPiDs, silent=TRUE);
    try(FileUniqueProcessIds <- NULL);
    try(FileUniqueProcessIds <- MyUniqueProcesses, silent=TRUE);
  }         
if (verbose >= 2) {
  AFilePrint("FileCons: All Finished, About to UnLockMe"); flush.console();
}
if (AlreadyLocked == FALSE) {
  try(setwd(AContainDir));
  UnLockMe(LFDir = OnNameFunction, verbose=verbose, quoteMore="Initial Sizeup");
  try(setwd(Oldwd));   AlreadyLocked=TRUE;
}
if (verbose >= 3) {
  AFilePrint(paste("ReadOneFile: About to Create FileList ", sep=""));
  flush.console();
}
FileList = list(FileCon = FileCon, FileSLen = FileSLen,
  FileMyDirectory=FileMyDirectory, FileMySimulationIds = FileMySimulationIds, 
  FileUniqueProcessIds = FileUniqueProcessIds,OnNameFunction =OnNameFunction);
if (verbose >= 3) {
  AFilePrint(paste("ReadOneFile: FileList created Quit. ", sep=""));
  flush.console();
}
return(FileList);
}


GiveFileCons <- function(OnNameFunctions = NULL, ASmallContainDir = NULL, 
  ALargeContainDir = NULL,
  RDataTableOutName = NULL, TableOutName=NULL, verbose=FALSE) {

if (((!exists("RDataTableOutName") || is.null(RDataTableOutName)) && 
     (!exists("TableOutName") || is.null(TableOutName))) || 
     !exists("OnNameFunctions") ||
    !exists("ASmallContainDir") || !exists("ALargeContainDir") ||
  is.null(OnNameFunctions) || is.null(ASmallContainDir) ||
  is.null(ALargeContainDir)) {
  AFilePrint("GiveFileCons: We got bad Input !");
  if (!exists("RDataTableOutName") || is.null(RDataTableOutName) ||
    (is.character(RDataTableOutName) && RDataTableOutName[1] == "")) {
    AFilePrint("GiveFileCons RDataTableOutName is not good!")
  }
  if (!exists("TableOutName") || is.null(TableOutName) ||
    (is.character(TableOutName) && TableOutName[1] == "")) {
    AFilePrint("GiveFileCons TableOutName is not good!")
  }
  if (!exists("OnNameFunctions") || is.null(OnNameFunctions)) {
    AFilePrint("GiveFileCons OnNameFunctions is not good!")
  }
  if (!exists("ALargeContainDir") || is.null(ALargeContainDir)) {
    AFilePrint("GiveFileCons ALargeContainDir is not good!")
  }
  if (!exists("ASmallContainDir") || is.null(ASmallContainDir)) {
    AFilePrint("GiveFileCons ASmallContainDir is not good!")
  }
  return(list(FileCons = NULL, FileLenS = NULL, 
        FileMyDirectories = NULL,
        FileUniqueSimulationIds = NULL,
        FileUniqueProcessIds = NULL  ));
}
FileCons = list();
FileLenS = rep(0, length(OnNameFunctions));
FileMyDirectories <- rep("", length(OnNameFunctions));
FileUniqueSimulationIds <- list();
FileUniqueProcessIds <- list();
      
if (is.numeric(verbose)) { if (verbose >= 2) { 
  verbose <- TRUE; 
} else {
  verbose <- FALSE; 
}}
if (verbose == TRUE) {
 AFilePrint(paste("GiveFileCons: About to Start", sep="")); flush.console();
}
if (verbose == TRUE) {
  AFilePrint(paste("GiveFileCons: TableOutName = ", TableOutName, sep=""));
}
AD = ASmallContainDir;
ALD = ALargeContainDir;

if ( (!exists("RDataTableOutName") || is.null(RDataTableOutName)) &&
  (exists("TableOutName") && !is.null(TableOutName))  ) {
  RDataTableOutName <- paste(TableOutName, ".RData", sep="");
  if (substr(TableOutName, nchar(TableOutName)-3, nchar(TableOutName)) %in% 
    c(".csv", ".txt", ".tst", ".Txt", ".Csv")) {
    RDataTableOutName <- paste(
      substr(TableOutName, 1, nchar(TableOutName)-4), ".RData", sep="");
  }
}

for (ii in 1:length(OnNameFunctions)) {
  if (verbose == TRUE) {
    AFilePrint(paste("GiveFileCons: Starting and looking into OnNameFunctions ii = ", 
      ii, " is ", OnNameFunctions[ii], sep="")); flush.console();
  }
    
  if (is.na(OnNameFunctions[ii]) || as.character(OnNameFunctions[ii]) == "NA") {
    AFilePrint(paste("GivefileCons for ii = ", ii," you had OnNamefunctions = ",
      OnNameFunctions[ii], sep=""));
    tryCatch("Cannot continue into this directory with FileCons!");
  }
  OnDir <- paste(ASmallContainDir, "//", OnNameFunctions[ii], sep="");
  MDir <- OnDir;
  dir.create(OnDir, showWarnings = FALSE, recursive=TRUE);
  TwoAddDir <- OnNameFunctions[ii];
  Oldwd <- getwd();
  try(setwd(ASmallContainDir));
  while(
    LockMeIn(verbose=verbose, 
    quoteMore=
    paste("GiveFileCons for Name Function ", OnNameFunctions[ii], sep=""), 
    LFDir = OnNameFunctions[ii], NoSpaceLFDir = TRUE) ==FALSE ) {}
  try(FilesInDir <- unlist(list.files(OnNameFunctions[ii]))); 
  try(setwd(Oldwd));
  if (any(FilesInDir == RDataTableOutName)) {
    RDD = NULL;
    ##assign("last.warning", NULL, envir = baseenv());
    AB = warnings();
    ##try(RDD <- read.table( paste(OnDir, "//", TableOutName, sep=""), sep=",",
    ##  header=TRUE, na.strings="NA", fill = TRUE, flush=TRUE,
    ##  stringsAsFactors = FALSE, col.names=MyColNames));
    try(Oldwd <- getwd());  try(setwd(ASmallContainDir));
    try(RDataRDD <- load(paste(OnNameFunctions[ii], "//", RDataTableOutName, sep="")));
    try(setwd(Oldwd));
    AW = warnings();
    if (length(AW)-length(AB) >= 1 || length(MyColNames) != NCOL(RDD)) {
      AFilePrint("**********************************************************");
      AFilePrint("**** ReadTable Issue Give FileCons "); flush.console();
      AFilePrint(paste("*** Trying to read table: ", TableOutName, sep=""));
      flush.console();
      AFilePrint(paste("***    in directory: ", OnDir, sep=""));
      AFilePrint(paste("*** ii = ", ii, " for OnNameFunction = ", 
        OnNameFunctions[ii], sep=""));
      flush.console();
      AFilePrint(paste("*** NCOL(RDD) = ", NCOL(RDD), " but length(MyColNames) = ",
        length(MyColNames), sep="")); flush.console();
      AFilePrint(paste("*** MyColNames = (", paste(MyColNames, collapse=", "),
        ")", sep="")); flush.console();
      AFilePrint(paste("*** RDD = ", sep=""));   flush.console();
      if (NROW(RDD) > 5) { AFilePrint(RDD[1:5,]); } else {AFilePrint(RDD); }
      AFilePrint("*** Now you have to go diagnose that error. "); flush.console();
      AText <- paste("**********************************************************\n",
         "**** ReadTable Issue Give FileCons   \n",
         "**** Trying to read table: ", TableOutName, " \n",
         "****    in directory: ", OnDir, " \n",
         "****    ii = ", ii, " for OnNameFunction = ", OnNameFunctions[ii],  "\n",
         "****  NCOL(RDD) = ", NCOL(RDD), " but length(MyColNames) = ", MyColNames, "\n",
         "****  MyColNames = (", paste(MyColNames, collapse=", "), ") \n",
         "****  Now You have to go diagnose that error. "); flush.console();
      AErrorPrint(AText);   
      tryCatch("GiveFileCons: I cannot continue here. ");
    }
    if (!is.null(RDD) && dim(RDD)[2] == NumDetails && dim(RDD)[1] > 0) {  
      try(FileLenS[ii] <- dim(RDD)[1]); try(FileCons[[ii]] <- RDD);
      try(FileMyDirectories[ii] <- MDir);
      try(FileUniqueSimulationIds[[ii]] <- MyPiDs);
      try(FileUniqueProcessIds[[ii]] <- MyUniqueProcesses);
    } else if (is.null(RDD) || dim(RDD)[2] > NumDetails) {
      if (is.null(RDD)) {
        AFilePrint("GiveFileCons: We Have a Multi-line Problem, NULL RDD");
        AErrorPrint("GiveFileCons: NULL RDD!");
        tryCatch("GiveFileCons: I don't want to continue!");
      } else if (is.null(dim(RDD)))  {
        AFilePrint("GiveFileCons: We Have a Multi-line Problem, dim NULL RDD");
        AErrorPrint("GiveFileCons: dim NULL RDD!");
        tryCatch("GiveFileCons: I don't want to continue!");
      } else if (dim(RDD)[2] > NumDetails) {
        AFilePrint(paste("GiveFileCons: We Have a Multi-line Problem, NumDetails = ", NumDetails, " but dim(RDD)[2] = ", 
          dim(RDD)[2], sep=""));
        AErrorPrint(paste("GiveFileCons: We Have a Multi-line Problem, NumDetails = ", NumDetails, " but dim(RDD)[2] = ", 
          dim(RDD)[2], sep=""));
        tryCatch("GiveFileCons: I don't want to continue!");  
      }
      RDD2 = matrix(as.numeric(RDD[1:length(RDD[,1]), 1:(NumDetails+1)]),
        length(RDD[,1]), NumDetails+1);
      if (is.null(MyColNames)) {
        AFilePrint("Oh, no we do not have MyColNames set up!");
      } else if (length(MyColNames) != NCOL(RDD2)) {
        AFilePrint(paste("Oh no, MyColNames length is ", length(MyColNames),
          "  and NCOL(RDD2) = ", NCOL(RDD2), sep="")); flush.console();
      } else if (!is.character(MyColNames)) {
        AFilePrint("My ColNames is not character! "); flush.console();
      }
      ABText <-  "
        ATFail = TRUE;
        write.table(RDD2, paste(OnDir, \"//\", TableOutName, sep=\"\"), 
          sep=\", \", eol=\"\\n\", na=\"NA\", append = FALSE, 
          row.names=FALSE, col.names=MyColNames);
        save(RDD=RDD2, RDD2=RDD2, MyColNames=MyColNames, 
          file=paste(OnDir, \"//\", RDataTableOutName, sep=\"\")); 
        ATFail = FALSE;
        ";
        try(eval(parse(text=ABText)));
        if (ATFail == TRUE) {
          AFilePrint("GiveFileCons: Oh no, we have a fail with col.names=FALSE, row.names=MyColNames");
          flush.console();
          AFilePrint(paste("Length MyColNames is ", length(MyColNames), sep=""));
          flush.console();
          AFilePrint(paste("They are: ", AFilePrint(MyColNames, collapse=", "), sep=""));
          flush.console();
          AFilePrint(paste("dim(RDD2) = (", paste(dim(RDD2), collapse=", "), ")", sep=""));
          AFilePrint(paste("We were trying to write to ", .self$TableOutName, sep=""));
          AFilePrint("We will now cause deliberate failure."); flush.console();
          PleaseDontLookForMe <- NULL;
          rm(PleaseDontLookForMe);
          OhNoWeFail <- PleaseDontLookForMe;
        }
      RDD = RDD2;
      try(FileLenS[ii] <- dim(RDD)[1]); try(FileCons[[ii]] <- RDD);
      try(FileMyDirectories[ii] <- MDir);
      try(FileUniqueSimulationIds[[ii]] <- MyPiDs);
      try(FileUniqueProcessIds[[ii]] <- MyUniqueProcesses);
      
      ##FileLenS[ii] = dim(RDD)[1]; FileCons[[ii]] = RDD;
    } else {
      try(FileLenS[ii] <- 0); try(FileCons[[ii]] <- -1);
      try(FileMyDirectories[ii] <- MDir);
      try(FileUniqueSimulationIds[[ii]] <- -1);
      try(FileUniqueProcessIds[[ii]] <- -1);
    }
    try(FileLenS[ii] <- dim(RDD)[1]); 
    TryText <- "
    FileCons[[ii]] <- -1;
    if (!is.null(RDD) && NROW(RDD) >= 1) {
      try(FileCons[[ii]] <- RDD);
    } else {
      try(FileCons[[ii]] <- -1);
    }
    ";
    try(eval(parse(text=TryText)));
    try(FileMyDirectories[ii] <- MDir);
    try(FileUniqueSimulationIds[[ii]] <- MyPiDs);
    try(FileUniqueProcessIds[[ii]] <- MyUniqueProcesses);
  } else {
    try(FileLenS[ii] <- 0); try(FileCons[[ii]] <- -1);
    try(FileMyDirectories[ii] <- MDir);
    try(FileUniqueSimulationIds[[ii]] <- -1);
    try(FileUniqueProcessIds[[ii]] <- -1);
  }  
  try(Oldwd <- getwd()); try(setwd(ASmallContainDir));      
  try(UnLockMe(LFDir = TwoAddDir, verbose=verbose, quoteMore="Initial Sizeup"));
  try(setwd(Oldwd));
}
if (verbose == TRUE) {
  AFilePrint("FileCons: All Finished, About to UnLockMe"); flush.console();
}

FileList = list(FileCons = FileCons, FileLenS = FileLenS,FileMyDirectories=FileMyDirectories,
  FileUniqueSimulationIds = FileUniqueSimulationIds, FileUniqueProcessIds = FileUniqueProcessIds);
return(FileList);
}


TableCompareSaverR5$methods(
  RunExperimentsAndSave = function(verbose=-1,...) {
  if (exists("verbose") && !is.null(verbose) && verbose == TRUE) {
    try(.self$verbose <- as.integer(1), silent=TRUE);
  } else if (is.numeric(verbose) &&  verbose > 0) {
    try(.self$verbose <- as.integer(verbose), silent=TRUE);
  } else if (is.logical(verbose) && verbose == FALSE) {
    try(.self$verbose <- as.integer(0), silent=TRUE);
  } else if (is.numeric(verbose) && length(verbose) == 1 &&
    verbose != -1  && verbose <= 0) {
    try(.self$verbose <- as.integer(verbose), silent=TRUE);  
  }
  try(setwd(.self$OriginalOldwd));
  if (.self$verbose >= 1) { 
    if (FALSE) {
    AFilePrint(paste(PRString(), "; Starting AnotherExperiment, 
      NumReps = ", .self$TargetLength,
    " with ", min(.self$FileLenS), " already done", sep="")); 
    }
    flush.console();
  } 
 
  IntervalAssess = -1;
  if (.self$verbose>=1) {
    AFilePrint("Running Experiments And Save Start!!! ");
    flush.console();
  }
  if (.self$verbose >= 1 || .self$verbose >0) {
    AFilePrint(paste(" Before we run, TableOutName = ", .self$TableOutName, sep=""));  
    flush.console();
  }
  if (FALSE) {
  if (min(.self$FileLenS) >= .self$TargetLength) {
    if (.self$verbose >= 3) {
      AFilePrint(paste("we get that .self$FileLens = ", .self$FileLenS,
        " and .self$TargetLength = ", .self$Targetlength, sep=""));
      flush.console();
    }
    return;
  }
  }
  AD <- 100;
  ##while(min(.self$FileLenS) < .self$TargetLength) {
  try(ABTrue <- TRUE);
  while(ABTrue)  { 
    if (.self$verbose >= 3) {
      AFilePrint(paste("Looping: we get that .self$FileLens = ", .self$FileLenS,
        " and .self$TargetLength = ", .self$Targetlength, sep=""));
    }
    .self$RunOneSetAndPrint();
    ##AFilePrint(paste(" AFter RunOneSetAndPrint, min(FileLenS) = ", 
    ##  min(this$FileLenS), sep="")); flush.console();
    
     if (AD != 100) { ABTrue <- FALSE }  
     TryValid <- " NotGood <- TRUE;
     if (!is.null(.self$ValidSimulations) && 
       length(.self$ValidSimulations) >= 1 &&
       !all(.self$CurrentTotalCompleted[.self$ValidSimulations] < 
         .self$TargetLength[1])) {
        ABTrue <- FALSE;  
      } 
      NotGood <- FALSE";
      try(eval(parse(text=TryValid)));
      if (NotGood) {
        AFilePrint("ERROR ERROR on Run Experiment");
        AFilePrint("Run Experiment, we have an issue with TCS Run Experiments");
        try(AFilePrint("The Valid Sims are: "));
        try(AFilePrint(paste("Run Experiment: ValidSims are (",
          paste(.self$ValidSimulations), collapse=", "), ")", sep=""))
        AErrorPrint("  Error: we have a Valid Exeperiment fail in RunExperiment");
      }
      if (MaxRunProcTime <= (proc.time()[3] -MyStartProcTime[3])) {
        ABTrue <- FALSE;
      }
  }
  return;
});

TableCompareSaverR5$methods(
  RunExperiment = function(
  verbose = -1,...) {
  if (is.logical(verbose) && verbose == TRUE) { 
    if (.self$verbose >= 1) {
    } else {
      try(.self$verbose <- as.integer(verbose), silent=TRUE); }
  } else if (is.logical(verbose) && verbose == FALSE) {
    try(.self$verbose <- as.integer(0), silent=TRUE);
  } else if (is.numeric(verbose) && verbose >= 1 && .self$verbose < 1) {
    try(.self$verbose <- as.integer(verbose), silent=TRUE);
  } else if (is.numeric(verbose) && length(verbose) == 1 && verbose == -1) {
  
  } else {
    try(.self$verbose <- as.integer(verbose), silent=TRUE);
  } 
  try(setwd(.self$OriginalOldwd)); 
  eval(parse(text=GetG0Text("Practice")));
  AD = 100;
  if (is.null(.self$TargetLength) || !is.numeric(.self$TargetLength)) {
    AErrorPrint("Hey An error before check lengths, TargetLength is NULL or not numeric!");
  }
  if (length(.self$TargetLength) <= 0) {
    AErrorPrint("Hey an Error, TargetLength in TCS has zero length!"); 
  }
  if (.self$TargetLength[1] <= 0) {
    AErrorPrint("Hey an Error, Target Length <= 0");
  }

  eval(parse(text=GetG0Text("MaxRunProcTime", "globalenv()", S=1)));
  eval(parse(text=SetGText("MyStartProcTime", "globalenv()", S=1)));
  eval(parse(text=GetGText("p", "globalenv()", S=1)));
  if (p >= 2000) {
     MaxRunProcTime <- 23/24 * 24 * 60 * 60;
  }
  if (p >= 200000) {
     MaxRunProcTime <- 6.8 * 24 * 60 * 60;
  }
  ABTrue <- TRUE;
  while(ABTrue)  { 
    RT = NULL;
    if (.self$verbose >= 3) {
      AFilePrint(paste("Looping: Trying to run with FileLenS = ",
        .self$FileLenS, " and TargetLength = ", .self$TargetLength, sep=""));
      flush.console();
    }
    try(RT <- .self$RunOneSimExperiment());
    if (is.null(RT) || (is.numeric(RT) && RT < 0)) {
      AFilePrint("RunExperiment: RunOneSimExperiment Returned Null, inspect the error!"); flush.console();
      return(-1);
    } else if (Practice == 1 && RT == -1) {
      AFilePrint("RunExperiment Encounters Try Error !"); flush.console();
      return(-1);
    }
    if (FALSE && AD != 100) { 
      ABTrue <- FALSE; 
      try(AFilePrint(paste(
        "RunExperiment: Kill ABTrue because AD == ", 
        AD, sep="")));
    }  
     TryValid <- " NotGood <- TRUE;
     if (!is.null(.self$ValidSimulations) && 
       length(.self$ValidSimulations) >= 1 &&
       !all(.self$CurrentTotalCompleted[.self$ValidSimulations] < 
         .self$TargetLength[1])) {
        try(AFilePrint(paste(\"RunExperiment: Kill ABTrue because \",
          \"we hit Valid Simulations\", sep=\"\")));
        ABTrue <- FALSE;  
      } 
      NotGood <- FALSE";
      try(eval(parse(text=TryValid)));
      if (NotGood) {
        AFilePrint("ERROR ERROR on Run Experiment");
        AFilePrint("Run Experiment, we have an issue with TCS Run Experiments");
        try(AFilePrint("The Valid Sims are: "));
        try(AFilePrint(paste("Run Experiment: ValidSims are (",
          paste(.self$ValidSimulations), collapse=", "), ")", sep=""))
        AErrorPrint("  Error: we have a Valid Exeperiment fail in RunExperiment");
      }
      try(AFilePrint(paste("At End of sim MaxProcTime = ", MaxRunProcTime, 
        ", and Run Time is ", proc.time()[3]-MyStartProcTime[3], sep="")));
      if (MaxRunProcTime <= (proc.time()[3] -MyStartProcTime[3])) {
        AFilePrint("We have Gone Over Time");
        ABTrue <- FALSE;
      }
    if (MaxRunProcTime <= (proc.time()[3] -MyStartProcTime[3])) {
      AFilePrint("RunExperiment: Quit RunOneSimExperiment because We've run out of Time!");
      AFilePrint(paste("MaxRunProcTime = ", MaxRunProcTime, 
        ",  But proc.time()[3] = ", proc.time()[3], "  and MyStartProcTime[3] = ",
        MyStartProcTime[3], sep=""));
    }
  }
  if (MaxRunProcTime <= (proc.time()[3] -MyStartProcTime[3])) {
      AFilePrint("RunExperiment: All at last Quit RunOneSimExperiment because We've run out of Time!");
      AFilePrint(paste("MaxRunProcTime = ", MaxRunProcTime, 
        ",  But proc.time()[3] = ", proc.time()[3], "  and MyStartProcTime[3] = ",
        MyStartProcTime[3], sep=""));
  }  
},
  RunExperimen = function(verbose = -1,...) {
    return(.self$RunExperiment(verbose, ...));
  }
);

TableCompareSaverR5$methods(
  RunOneSetAndPrint = function(UniqueIdentifier, verbose = -1,...) {
  if (is.numeric(verbose) && length(verbose) == 1 && verbose == -1) {
  
  } else if (is.logical(verbose)) {
    if (verbose == TRUE) { 
      if (.self$verbose >= 1) {
      } else {
        .self$verbose <- as.integer(1); 
      } 
    } else {
      .self$verbose <- as.integer(0);
    }
  } else if (is.numeric(verbose)) {
    try(.self$verbose <- as.integer(verbose[1]));
  }
  if (.self$verbose >= 3) {
    AFilePrint("RunOneSetAndPrint: Starting."); flush.console();
  }
  

  MatrixBDDL = matrix(0, .self$NumPrint, 
    .self$NumDetails * length(.self$OnNameFunctions));  

  if (is.null(.self$FileCons)) {
    AFilePrint("RunOneSetAndPrint Wants FileCons!! can't do it");
    flush.console();
  }  
  if (is.null(.self$FileLenS)) {
    AFilePrint("RunOneSetAndPrint: Need FileLenS!! "); flush.console();
  }
  if (is.null(.self$OnNameFunctions) || length(.self$OnNameFunctions) <= 0 ||
    (length(.self$OnNameFunctions) == 1 && .self$OnNameFunctions == "")) {
    AFilePrint("RunOneSetAndPrint: Need OnNameFunctions")
  }
  ##if (is.null(this$experimentFunction)) {
  ##  AFilePrint("RunOneSetAndPrint: this$experimentFunction is NULL!"); flush.console(); 
  ##}
 
  ##for (njj in 1:this$NumPrint) {
  ##  PullOff = this$experimentFunction(njj, this$ListOfParams);
  ##  MatrixBDDL[njj,] = PullOff$LS
  ##  if (FALSE) {
  ##    try(UpdateIntervalAssess());
  ##  }    
  ##}
    
##############################################################################
##  Wait until LockFile gives opening, then write to table  
##                
##AFilePrint(paste(PRString(), "; min(FileLenS) == ", min(this$FileLenS), 
##    "  Finished ", this$NumPrint, " Experiments", sep=""));  flush.console();   
##AFilePrint(paste(PRString(),
##  "; min(FileLenS) == ", min(this$FileLenS), 
##    ";   Reading Tables To Check ", sep=""));  
  while(LockMeIn(
    LFDir = .self$ASmallContainDir, 
    verbose=.self$verbose, quoteMore="ReadAndRunOnePrint") ==FALSE ) {} 
  if (.self$verbose >= 2) {
    AFilePrint("Succeeded bypassing LockMeIn for ReadAndRunOnePrint"); flush.console();
  }             

################################################################################
##  Now begin saving to directories
  for (ii in 1:length(.self$OnNameFunctions)) {
    OnDir <- paste(.self$ASmallContainDir, "//", .self$OnNameFunctions[ii], sep="");
    dir.create(OnDir, showWarnings = FALSE)
    FilesInDir = unlist(list.files(OnDir));

#############Debugging PrintCode    
if (.self$verbose >= 1) {
  AFilePrint(paste(" # of Files is ", length(FilesInDir), " in ", OnDir, sep="" ));
  AFilePrint(paste(" The File Names are: ", paste(FilesInDir, collapse = ", "), sep=""));
  AFilePrint(paste("And self$TableOutName = ", .self$TableOutName, sep=""));
  AFilePrint(paste("any Trick: ", any(FilesInDir == .self$TableOutName), sep=""));
}
    if (any(FilesInDir == .self$TableOutName)) {
      RDD = NULL;
      ##AFilePrint("Trying to read in RDD, which can be read");  flush.console();
      
      #########################################################################
      ## Moment of truth, try to read the table
      try(RDD <- read.table( paste(OnDir, "//", .self$TableOutName, sep=""),stringsAsFactors = FALSE, 
        sep=",", header=FALSE, fill = TRUE, flush=TRUE,
        col.names=1:.self$NumDetails));
      if (is.null(RDD) && .selfverbose >= 1) {
        AFilePrint(paste("ReadOneAndPrint: We Read Null RDD for : ", .self$TableOutName, sep=""));
      }
      if (is.null(RDD)) {
        AFilePrint("Sorry RDD was NULL");
        flush.console();
      } else {
      ##  AFilePrint(paste("dimRDD = c(", paste(dim(RDD), collapse=", "), ")", sep="")); 
      ##  flush.console();
      }
      if (!is.null(RDD) && dim(RDD)[2] == .self$NumDetails) { 
        if (.self$verbose >= 1) {
          AFilePrint(paste("ReadOneAndPrint: Appending FileLenS, already",
            dim(RDD)[1], sep=""));
        }
        try(.self$FileLenS[ii] <- dim(RDD)[1]);
        ##AFilePrint("This Seems to have been a good Table!!"); flush.console();
        
        ########################################################################
        ## Moment of truth, try to append the table
        try(write.table(
          matrix( round( MatrixBDDL[1:.self$NumPrint,
            1:.self$NumDetails + (ii-1) * .self$NumDetails],8),
            .self$NumPrint, .self$NumDetails), 
            file = paste(OnDir, "//", .self$TableOutName, sep=""), 
            append = TRUE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double")));
          .self$FileLenS[ii] <- .self$FileLenS[ii] + .self$NumPrint;       
      } else if (is.null(RDD) || dim(RDD)[2] > .self$NumDetails) {
        ##AFilePrint("We Have a Multi-line Problem we must correct");
        RDD2 = matrix(as.numeric(RDD[1:length(RDD[,1]), 1:.self$NumDetails]),
        length(RDD[,1]), .self$NumDetails);
        if (.self$verbose >= 1) {
          AFilePrint(paste("ReadOneAndPrint: MultiLine, has wrong dim (",
            paste(dim(RDD2), collapse = ", "), ")",
            "yet NumPrint = ", .self$NumPrint, sep=""));
        }
        ABText <-  "
        ATFail = TRUE;
        write.table(RDD2, paste(OnDir, \"//\", .self$TableOutName, sep=\"\"), 
          sep=\", \", eol=\"\\n\", na=\"NA\", append = FALSE, 
          row.names=FALSE, col.names=FALSE); 
        ATFail = FALSE;
        ";
        try(eval(parse(text=ABText)));
        if (ATFail == TRUE) {
          AFilePrint("Oh no, we have a fail with col.names=FALSE, row.names=FALSE");
          flush.console();
          AFilePrint(paste("dim(RDD2) = (", paste(dim(RDD2), collapse=", "), ")", sep=""));
          AFilePrint(paste("We were trying to write to ", .self$TableOutName, sep=""));
          AFilePrint("We will now cause deliberate failure."); flush.console();
          PleaseDontLookForMe <- NULL;
          rm(PleaseDontLookForMe);
          OhNoWeFail <- PleaseDontLookForMe;
        }
        RDD = RDD2
        .self$FileLenS[ii] <- dim(RDD)[1];
        try(write.table( matrix( round( MatrixBDDL[1:.self$NumPrint,
            1:.self$NumDetails + (ii-1) * .self$NumDetails],8),
            .self$NumPrint, .self$NumDetails) , 
          file = paste(OnDir, "//", .self$TableOutName, sep=""), 
          append = TRUE, quote = FALSE, sep = ", ",
          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
          col.names = FALSE, qmethod = c("escape", "double")));
        .self$FileLenS[ii] <- .self$FileLenS[ii] + .self$NumPrint;
      } else {
        ##AFilePrint("The Length was short!");
        try(.self$FileLenS[ii] <- 0);
        AFilePrint(paste(PRString(), "; min(.self$FileLenS) == ", 
          min(.self$FileLenS), ";   File for OnNameFunctions[", ii, "] = ",
          .self$OnNameFunctions[ii], " is bad dimensions.", sep=""));   
          AFilePrint(paste("dim(RDD) = c(", paste( dim(RDD), collapse=", "),
          ")", sep=""));
        AFilePrint(paste("ReadOneAndPrint: Unfortunately can't work with wrong dim (",
            paste(dim(RDD2), collapse = ", "), ")",
            "yet NumPrint = ", .self$NumPrint, sep=""));
        if (.self$verbose >= 1) {
          AFilePrint(paste("ReadOneAndPrint: Unfortunately can't work with wrong dim (",
            paste(dim(RDD2), collapse = ", "), ")",
            "yet NumPrint = ", .self$NumPrint, sep=""));
        }
        quit();
        try(write.table( 
          round(MatrixBDDL[ , (ii-1) * .self$NumDetails + 1:.self$NumDetails], 8),
          paste(OnDir, "//", TableOutName, sep=""), sep=",",
          row.names=FALSE, col.names=FALSE,
          quote=FALSE, na="NA", append=FALSE));
        .self$FileLenS[ii] = length(MatrixBDDL[,1]);                          
      }
      ##AFilePrint("End of this strange Else");
    } else {
      .self$FileCons[[ii]] <- as.integer(-1);
      try(.self$FileLenS[ii] <- as.integer(0), silent=TRUE);
      AFilePrint(paste(PRString(), "; min(FileLenS) == ", min(.self$FileLenS), 
        ";   File for OnNameFunctions[", ii, "] = ",
        .self$OnNameFunctions[ii], " doesn't exist", sep=""));     
      ##AFilePrint(paste("ReadOneAndPrint: What the flying F!! We had fail to load file", sep=""));   
      ##quit();
      if (.self$verbose >= 1) {
        AFilePrint(paste("ReadOneAndPrint: Failed To Open ", sep=""));
      }       
      try(write.table( 
        matrix(round( 
          MatrixBDDL[ 1:.self$NumPrint, (ii-1) * .self$NumDetails + 
            1:.self$NumDetails], 8),
          .self$NumPrint, .self$NumDetails),
          paste(OnDir, "//", TableOutName, sep=""), sep=",",
          row.names=FALSE, col.names=FALSE,
          quote=FALSE, na="NA", append=FALSE));
      try(.self$FileLenS[ii] <- length(MatrixBDDL[,1]));  
    }
     
}
 
UnLockMe(LFDir = .self$ASmallContainDir, verbose=.self$verbose, quoteMore="  After Writing");
});


###############################################################################
##  Functions to Assess Confidence intervals, don't work very well
##
##
UpdateInvervalAssess <- function(verbose=FALSE) {
  if (verbose==TRUE) {AFilePrint("Updating Invertval Assess");}
  if (length(IntervalAssess) == 1 && is.numeric(IntervalAssess)) {
    IntervalAssess = PullOff$IntAssess;    
  } else {
    for (ii in 1:length(PullOff$IntAssess)) {
      Gotcha <<- 0;
      for (jj in 1:length(IntervalAssess)) {
        if (IntervalAssess[[jj]]$nameF == PullOff$IntAssess[[ii]]$nameF) {
          IntervalAssess[[jj]]$AssessedCI <<- 
            rbind(IntervalAssess[[jj]]$AssessedCI,
              PullOff$IntAssess[[ii]]$AssessedCI);
          Gotcha = jj;
          break; 
          NumAssessed = NumAssessed + 1; 
        }     
      }
    }         
  }
  return;
}

################################################################################
##  GiveRDTL tries to read the tables for the given functions and
##   piece them together in a table.
##
##  However, R tends to load numeric data into data-structures that do not
##  convert to the correct numeric columns.  As of 12/11/2010 the solution
##  seems to be by deliberately specifying for all read.table columns
##  to use numeric typing.  A lot of code was experimented with, trying
##  to unlist data and re-convert it to numeric data, however, it appears
##  to fail and generate false numbers
##
##
TableCompareSaverR5$methods(
  GiveRDTL = function(verbose=-1) {
  if (is.logical(verbose) && verbose == TRUE) {
    if (.self$verbose < 1) { .self$verbose <- 1; }
  } else if (is.numeric(verbose) && length(verbose) ==1 && verbose== -1) {
  
  } else if (is.logical(verbose)) {
    if (verbose == TRUE) {
      if (.self$verbose >= 1) {
      } else {
        try(.self$verbose <- as.integer(1));
      }
    } else {
      try(.self$verbose <- as.integer(0));
    }
  } else {
    try(.self$verbose <- as.integer(verbose[1]));
  }

   if (.self$verbose >= 1) {
     AFilePrint("GiveRDTL, Starting"); flush.console();
   }
   jobs = NULL;
   RDTL = list();
   MyTryT <- "FileLenS = rep(0, length(.self$OnNameFunctions));
     .self$FileLenS <- FileLenS; ";
   Oldwd <- getwd();
   try(eval(parse(text=MyTryT)), silent=TRUE);

   jobs = matrix(0, .self$TargetLength, 
       .self$NumDetails * length(.self$OnNameFunctions))
   for (ii in 1:length(.self$OnNameFunctions)) {
        try(setwd(.self$ASmallContainDir));
        while(LockMeIn(LFDir=.self$OnNameFunctions[ii], verbose=.self$verbose, 
          paste(quoteMore = "GiverRDTL: jobii = ", .self$jobii, sep="") ) == FALSE ) {}
       if (.self$verbose >= 2) {
         AFilePrint(paste("GiveRDTL, on ii = ", ii, sep="")); flush.console();
       }
       try(dir.create(.self$OnNameFunctions[ii], recursive=TRUE, showWarnings = FALSE));
       try(RDTL[[ii]] <- 0);
       try(MyL <- unlist(list.files(.self$OnNameFunctions[ii])));
       if (.self$TableOutName %in% MyL) {
         try(RDTL[[ii]] <- read.table(paste(.self$OnNameFunctions[ii], "//", 
           .self$TableOutName, sep=""), stringsAsFactors = FALSE, 
           header=TRUE, sep=",", na.strings=c("NA", "NaN", "na"), 
           fill=TRUE, flush=TRUE,
           ));
      }
      try(UnLockMe(LFDir = .self$OnNameFunctions[ii], verbose=.self$verbose,
       quoteMore = paste("GetRDTL, jobii = ", jobii, sep=""))); 
      try(setwd(Oldwd)); 
     ##if (.self$OnNameFunctions[ii] == "2LassoSpread") {
     ##  AFilePrint("Well at first loop: Where are We Now?  GiveRDTL on 2Lasso Spread");
   ##  AFilePrint("RDTL[[ii]] First Line is");
   ##  AFilePrint(RDTL[[ii]][1,]) ;
   ##  AFilePrint("Entire jobs now is ");
   ##  AFilePrint(jobs[1,])
   ##  AFilePrint(""); AFilePrint("");
   ##}
     if (length(RDTL) < ii || is.null(RDTL[[ii]]) || is.null(dim(RDTL)) ||
       dim(RDTL[[ii]])[2] != .self$NumDetails) {
       if (.self$verbose == TRUE) {
       AFilePrint(paste("GiveRDTL Error ii = ", ii, " for OnDir = ", OnDir, sep=""));
       AFilePrint(paste("    Our TableOutName is ", .self$TableOutName, sep="")); flush.console();
       if (length(RDTL) < ii) {
         AFilePrint(paste("RDTL has length ", length(RDTL), " which is less than ii = ", ii, sep=""));
       } else if (is.null(RDTL[[ii]])) {
         AFilePrint("RDTL[[ii]] is null!");
       } else if (is.null(dim(RDTL))) {
         AFilePrint("RDTL[[ii]] has no dimension!"); flush.console();
       } else {
         AFilePrint(paste(".self$NumDetails = ", .self$NumDetails, 
           ", but RDTL[[ii]] has dim (", 
           paste(dim(RDTL[[ii]]), collapse=", "), ")", sep="")); flush.console();
       }
       AFilePrint(paste(PRString(), "; min(FileLenS) == ", 
         min(.self$FileLenS), sep=""));         
       AFilePrint(paste("ReadingRDTL Namefunctions[", ii, "] = ",
         .self$OnNameFunctions[ii],"ReadLengthProblem!", sep=""));
       AFilePrint(paste(" dim(RDTL[[", ii, "]]) = c(", 
         paste( dim(RDTL[[ii]]), collapse=", "), sep=""));
       quit();
     }
     if (dim(RDTL[[ii]])[1] < .self$TargetLength) {
       AFilePrint(paste(PRString(),  "; min(FileLenS) == ", 
         min(.self$FileLenS), sep=""));         
       AFilePrint(paste("ReadingRDTL Namefunctions[", ii, "] = ",
         .self$OnNameFunctions[ii], "ReadLengthProblem!", sep=""));
           AFilePrint(paste(" dim(RDTL[[", ii, "]]) = c(", 
             paste( dim(RDTL[[ii]]), collapse=", "), sep=""));
           quit();
     }
   }
   try(JTB <- matrix(0, length(RDTL[[ii]][,1]), .self$NumDetails ));
   for (jj in 1:.self$NumDetails)  {
     try(JTB[,jj] <- as.numeric(unlist(RDTL[[ii]][,jj])));
   }
   ##AFilePrint("JTB[1,] is ");
   ##AFilePrint(JTB[1,]);
   ##AFilePrint("RDTL[[ii]] is ");
   ##AFilePrint(RDTL[[ii]][1,]);
   try(RDTL[[ii]] <- JTB);
   ##RDTL[[ii]] = data.matrix(RDTL[[ii]]);
   ##if (.self$OnNameFunctions[ii] == "2LassoSpread") {
   ##  AFilePrint("Well after DataMatrix Change GiveRDTL on 2Lasso Spread");
   ##  AFilePrint("RDTL[[ii]] First Line is");
   ##  AFilePrint(RDTL[[ii]][1,]) ;
   ##  AFilePrint("Entire jobs now is ");
   ##  AFilePrint(jobs[1,])
   ##  AFilePrint(""); AFilePrint("");
   ##}
   ##AFilePrint(RDTL[[ii]][1,])
   if ( ii == 1 ) {
     ##AFilePrint(paste("ii is 1 lets look at RDTL[[ii]]"));
     ##AFilePrint(RDTL[[ii]]);
     ##try(jobs[1:.self$TargetLength, 1:.self$NumDetails] <- 
     ##  (RDTL[[ii]])[1:.self$TargetLength, 1:.self$NumDetails]);
     

     ##jobs = matrix(as.numeric(unlist(jobs)), length(jobs[,1]), length(jobs[1,]));
   } else {       

     ##try(jobs[1:.self$TargetLength, (ii-1) * .self$NumDetails + 1:.self$NumDetails] <- 
     ##  (RDTL[[ii]])[1:.self$TargetLength, 1:.self$NumDetails]); 
     ##AFilePrint("Talking about jobs:");
     ##AFilePrint("Printing unlist(jobs): ");
     ##AFilePrint(unlist(jobs)[1:20]);
     ##AFilePrint("But jobs[1:20,1:NumDetails*2] was");
     ##AFilePrint(jobs[1:20,1:(.self$NumDetails*2)]);

     ##jobs = matrix(as.numeric(unlist(jobs)), length(jobs[,1]), length(jobs[1,]));
     ##AFilePrint("But jobs[1:20,1:NumDetails*2] after sprint was");
     ##AFilePrint(jobs[1:20,]);    
     ##quit();
   }

   ##AFilePrint(paste(" As we expand, dim(jobs) = ", 
   ##  paste(dim(jobs), collapse=", "), ")", sep=""));
  } 
    
  ##jobs = matrix(0, .self$TargetLength, 
  ##  .self$NumDetails * length(this$OnNameFunctions));
  for (ii in 1:length(.self$OnNameFunctions)) {
    ##AFilePrint(paste("ii = ", ii)); flush.console();
    ##AFilePrint(paste("dim(jobs) == ", paste(dim(jobs), collapse=", "), sep=""));
    ##AFilePrint(paste("dim(RDTL[[ii]]) = ", paste(dim(RDTL[[ii]]), collapse=", "),
    ##  sep="")); flush.console();
    ##AFilePrint(paste(" And Length RDTL[[ii]][1,] is ", length(RDTL[[ii]][1,]), sep=""));
    ##flush.console();
    try(jobs[1:(.self$TargetLength), (ii-1) * .self$NumDetails + 
      1:(.self$NumDetails)] <- 
      (RDTL[[ii]])[1:(.self$TargetLength), 1:(.self$NumDetails)]);
   ##if (this$OnNameFunctions[ii] == "2LassoSpread") {
   ##  AFilePrint("GiveRDTL on 2Lasso Spread");
   ##  AFilePrint("RDTL[[ii]] First Line is");
   ##  AFilePrint(RDTL[[ii]][1,]) ;
   ##  AFilePrint("Entire jobs now is ");
   ##  AFilePrint(jobs[1,])
   ##}
  }
  return(jobs);
});

SaveIntAsses <- function() {
 if (FALSE && any(FilesInDir ==  IntAssessOutName)) {
        ##try(RDDIA = read.table( paste(OnDir, "//", IntAssessOutName, sep=""), 
        ##         sep=",",
        ##         header=FALSE, fill = TRUE, flush=TRUE,
        ##         col.names=1:NumIntervalAssess));
        if (dim(RDDIA)[2] == NumIntervalAssess) {  
          ##FileLenS[ii] = dim(RDD)[1]; 
        } else if (dim(RDDIA)[2] > NumIntervalAssess) {
          AFilePrint("We Have a IA Multi-line Problem we must correct");
            RDD2IA = matrix(as.numeric(RDDIA[1:length(RDDIA[,1]), 
            1:NumIntervalAssess]), length(RDDIA[,1]), NumIntervalAssess);
          ##try(write.table(RDD2IA, paste(OnDir, "//", IntAssessOutName, sep=""), 
          ##   sep=", ", eol="\n", na="NA", append = FALSE, 
          ##   row.names=FALSE, col.names=FALSE));
          RDDIA = RDD2IA
          ##FileLenS[ii] = dim(RDD)[1];
        } else {
          ##FileLenS[ii] = 0;
          AFilePrint(paste(PRString(), 
            "; min(FileLenS) == ", min(FileLenS), 
            ";   IntAsses File for NameFunctions[", ii, "] = ",
            this$OnNameFunctions[ii], " is bad dimensions.", sep=""));   
          AFilePrint(paste("dim(RDDI) = c(", 
            paste( dim(RDDIA), collapse=", "), ")", sep=""));                      
        }
  }
  return;
}

GetRealTime <- function() {    
  AllMonths <- c(31,28,31,30, 31, 30, 31, 31, 30, 31, 30, 31)
  op <- options(digits.secs=6)
  OurTime <- as.character(Sys.time())
  options(op)
  MyN <- unlist(strsplit(as.character(OurTime), " "));
  A1 <- as.numeric(unlist(strsplit(MyN[1], "-")));
  A2 <- as.numeric(unlist(strsplit(MyN[2], ":")));
  LastYear <- A1[1]-1;
  Num400s <- LastYear %/% 400;
  Num100s <- LastYear %/% 100;
  NumLeaps <- LastYear %/% 4;
  isLeap <- A1[1] %% 4 == 0 & !(A1[1] %% 100 == 0 & !(A1[1] %% 400 == 0));
  FormerLeaps <- NumLeaps - Num100s + Num400s;
  
  LastDays <- 365 * LastYear +  FormerLeaps; 
  Yesterday <- A1[3] -1;
  LastMonth <- A1[2] -1;
  if (LastMonth == 0) {
    NumDays <- Yesterday;
  } else {
    NumDays <- Yesterday + sum(AllMonths[1:LastMonth]);
  }
  if (isLeap == TRUE) {
    if (LastMonth >= 2) {
      NumDays <- NumDays+1;
    }
  }
  
  JustSeconds <- 60 * 60 * A2[1] + 60 * A2[2] +A2[3] 
  JustDays <- (LastDays + NumDays)
  AllSeconds <- JustDays * 60 * 60 * 24 + JustSeconds;
  return(
    list(AllSeconds=AllSeconds, JustSeconds=JustSeconds, JustDays=JustDays, 
    OurTime=OurTime, FormerLeaps=FormerLeaps, isLeap=isLeap,NumDays=NumDays+LastDays));
}
      