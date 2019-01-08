CollectIntoTable <- function(OrderedRhos, OrderedSigmas, 
  NeedColumns, WantLines = LocateLines, WantBreaks= NULL, RDSigFig = 3) {

  MyRhos <- unlist(ListUseRho);
  MySigmas = unlist(ListSigmaNoise);
  
  MyRS = paste(MyRhos, MySigmas, sep="o"); 
  OORS = paste(OrderedRhos, OrderedSigmas, sep="o"); 
  IMatch <- match(OORS, MyRS);
  OrderedRhos = OrderedRhos[!is.na(IMatch)]
  OrderedSigmas = OrderedSigmas[!is.na(IMatch)]
  NNeedColumns <- list();
  for (ii in 1:length(NeedColumns)) {
    if (!is.na(IMatch[ii])) {
      NNeedColumns[[length(NNeedColumns)+1]] <- NeedColumns[[ii]];
    }
  }
  if (length(NNeedColumns) >= 1) {
    NeedColumns = NNeedColumns;
  }
  IMatch = IMatch[!is.na(IMatch)];
  
  
  LLen <- length(unlist(NeedColumns));
  MeanGo = matrix(0, length(WantLines), LLen);
  SDGo = MeanGo + 0.0;
  InLen <- 0;
  for (ii in 1:length(IMatch)) {
     if (is.character(NeedColumns[[ii]])) {
        IM <- match(NeedColumns[[ii]], colnames(ListMeanSum[[ii]]));
     }
     matchLines <- match(WantLines, rownames(ListMeanSum[[ii]]));
     MeanGo[(1:length(matchLines))[!is.na(matchLines)], 
       InLen + 1:length(IM)] <- ListMeanSum[[ii]][matchLines[!is.na(matchLines)],
       IM];
     SDGo[(1:length(matchLines))[!is.na(matchLines)], 
       InLen + 1:length(IM)] <- ListSDSum[[ii]][matchLines[!is.na(matchLines)],
       IM]; 
     InLen <- InLen + length(IM);
  }
   
  rownames(MeanGo) = WantLines;  rownames(SDGo) = WantLines;
  colnames(MeanGo) = unlist(NeedColumns); colnames(SDGo) = unlist(NeedColumns);
  
  if (!is.null(WantBreaks)) {
    if (length(WantBreaks) == length(WantLines)) {
      if (all(WantBreaks %in% c(0,1))) {
         WantBreaks = (1:length(WantLines))[WantBreaks == 1];
      }
    }
  }
  rd0T <- matrix("", NROW(MeanGo), NCOL(MeanGo));
  for (ii in 1:NROW(MeanGo)) {
    for (jj in 1:NCOL(MeanGo)) {
      rd0T[ii,jj] = paste(rd0(MeanGo[ii,jj], SigFig = RDSigFig),
        "(", rd0(SDGo, SigFig=RDSigFig), ")", sep="");
    }
  }
  eval(parse(text=SetGText("MeanGo", "globalenv()", S=1)));
  eval(parse(text=SetGText("SDGo", "globalenv()", S=1)));  
  RMT  = paste("\\begin{table}[htbp] 
    \\begin{centring}
    \\begin{tabular}{r|", paste(rep("r", NCOL(MeanGo)), collapse="&"),
      "|} \\\\
    ", sep="");
  for (ii in 1:length(NeedColumns)) {
    RMT <- paste(RMT, " & 
      \\multicolumn{", length(NeedColumns[[ii]]), "}{|c|}",
      "{", "\\sigma=", round(ListSigmaNoise[[ii]],2), "}", sep="");
  }
  RMT = paste(RMT, " \\\\
    ", sep="");
  eval(parse(text=GetG0Text("TimeString", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("ListP", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("ListK", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("ListN", "globalenv()", S=1)));
  pUse <- mean(unlist(ListP));
  kUse <- mean(unlist(ListK));
  nUse <- mean(unlist(ListN));

  
  Matchers <- rbind(c("l2error", "$L_{2}$"), c("Type1","T1"),
  c("Type2", "T2"), c("time:elapsed","time")); 
  Matchers[Matchers[,1] == "time:elapsed",2] <- TimeString;
  for (ii in 1:length(NeedColumns)) {
    RMT <- paste(RMT, " & ",
      paste(Matchers[match(NeedColumns[[ii]], Matchers[,1]),2], collapse = " & ")
    )
  }
  RMT <- paste(RMT, " \\\\
   ", sep="");
  NameLines <- AllKnownFormat[match(WantLines, AllKnownNameFunctions)];
  for (ii in 1:NCOL(MeanGo)) {
    if (ii %in% WantBreaks) {
      RMT <- paste( RMT, " \\hline 
        ", sep="")
    }
    RMT <- paste(RMT, NameLines[ii], " & ",
      paste(rd0T[ii,], collapse=" & "), " \\\\ 
      ", sep="");
  }
  if (kUse == 6) {
    kTreatment <- "true non-zero Beta are (1,1,1,-1,-1,-1).";
  } else {
    kTreatment <- paste("number of non-zero Beta are ", kUse, sep="");
  }
  RMT <- paste(RMT <- " \\hline 
  \\end{tabular}
  \\caption{Table of competing estimator performance where n=", nUse,
  ", p=", pUse, ", and ", kTreatment, "}
  \\end{centering}
  \\label{TabRhos", paste(tSeq(unlist(ListUseRho)), collapse=""),
    paste(tSeq(unlist(ListSigmaNoise)), collapse=""), "}
  \\end{table} ", sep="");
  eval(parse(text=SetGText("RMT", "globalenv()", S=1)));
  return(RMT);  
}

GroupInfoForTables <- function(OnDirs = NULL, TSEDir = "", SecondsMinutesHours = "Seconds", ChooseRho = NULL) {
  if (is.null(OnDirs)) {
    return;
  }
  ListMeanSum <- list();  ListSDSum <- list();
  ListMyMean <- list();  ListMySDSum <- list();
  ListMeanMIP <- list();  ListSDMIP <- list();
  ListMyMIP <- list();  ListMySDMIP <- list();
  ListAMatSD <- list();  ListAMatNeed <- list();
  ListSigmaNoise <- list();  ListUseRho <- list();
  ListN <- list();  ListP <- list(); ListK <- list();
  ListPG <- list();  ListKG <- list();
  GoL <- 0;
  for (ii in 1:length(OnDirs)) {
    ART <- -1;
    MyT <- "ART <- GetInfoForTables(OnDirs[ii], TSEDir, SecondsMinutesHours);";
    try(eval(parse(text=MyT)));
    if (ART != -1) {
    if (!is.null(ChooseRho) && is.numeric(ChooseRho) && length(ChooseRho) == 1) {
      eval(parse(text=GetG0Text("UseRho", "globalenv()", S=1)));
      if (UseRho == ChooseRho) {
        DoGo = TRUE;
      } else { DoGo <- FALSE; }
    } else { DoGo == TRUE; }
    if (DoGo == TRUE) {
    GoL <- GoL + 1;
    eval(parse(text=GetG0Text("MeanSum", "globalenv()", S=1)));
    ListMeanSum[[GoL]] <- MeanSum
    eval(parse(text=GetG0Text("SDSum", "globalenv()", S=1)));
    ListSDSum[[GoL]] <- SDSum
    GetInfoForTables(OnDirs[ii], TSEDir, SecondsMinutesHours);
    eval(parse(text=GetG0Text("MyMean", "globalenv()", S=1)));
    ListMyMean[[GoL]] <- MyMean
    eval(parse(text=GetG0Text("MySDSum", "globalenv()", S=1)));
    ListMySDSum[[GoL]] <- MySDSum
    
    eval(parse(text=GetG0Text("pSimulation", "globalenv()", S=1)));
    ListP[[ii]] <- pSimulation;
    eval(parse(text=GetG0Text("kSimulation", "globalenv()", S=1)));
    ListK[[ii]] <- kSimulation;
    eval(parse(text=GetG0Text("nSimulation", "globalenv()", S=1)));
    ListP[[ii]] <- nSimulation;
    eval(parse(text=GetG0Text("kGSimulation", "globalenv()", S=1)));
    ListKG[[ii]] <- kGSimulation;
    eval(parse(text=GetG0Text("pGSimulation", "globalenv()", S=1)));
    ListPG[[ii]] <- pGSimulation;

    eval(parse(text=GetG0Text("MeanMIP", "globalenv()", S=1)));
    ListMeanMIP[[GoL]] <- MeanMIP
    eval(parse(text=GetG0Text("SDMIP", "globalenv()", S=1)));
    ListSDMIP[[GoL]] <- SDMIP
    
    eval(parse(text=GetG0Text("MyMIP", "globalenv()", S=1)));
    ListMyMIP[[GoL]] <- MyMIP
    eval(parse(text=GetG0Text("MySDMIP", "globalenv()", S=1)));
    ListMySDMIP[[GoL]] <- MySDMIP

    eval(parse(text=GetG0Text("AMatSD", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("AMatNeed", "globalenv()", S=1)));
    ListAMatSD[[GoL]] <- AMatSD;
    ListAMatNeed[[GoL]] <- AMatNeed;
    eval(parse(text=GetG0Text("SigmaNoise", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("UseRho", "globalenv()", S=1)));  
    ListSigmaNoise[[GoL]] <- SigmaNoise;
    ListUseRho[[GoL]] <- UseRho;
    } 
    }   
  }
  eval(parse(text=SetGText("ListMeanSum", "globalenv()", S=1)));
  eval(parse(text=SetGText("ListSDSum", "globalenv()", S=1)));  
  eval(parse(text=SetGText("ListMyMean", "globalenv()", S=1)));
  eval(parse(text=SetGText("ListMySDSum", "globalenv()", S=1)));  
  eval(parse(text=SetGText("ListMeanMIP", "globalenv()", S=1)));
  eval(parse(text=SetGText("ListSDMIP", "globalenv()", S=1)));  

  eval(parse(text=SetGText("ListMyMIP", "globalenv()", S=1)));
  eval(parse(text=SetGText("ListMySDMIP", "globalenv()", S=1)));    
  eval(parse(text=SetGText("ListAMatSD", "globalenv()", S=1)));
  eval(parse(text=SetGText("ListAMatNeed", "globalenv()", S=1)));    
    
  eval(parse(text=SetGText("ListSigmaNoise", "globalenv()", S=1)));
  eval(parse(text=SetGText("ListUseRho", "globalenv()", S=1)));    

  eval(parse(text=SetGText("ListN", "globalenv()", S=1)));   
  eval(parse(text=SetGText("ListP", "globalenv()", S=1)));   
  eval(parse(text=SetGText("ListK", "globalenv()", S=1)));   
  eval(parse(text=SetGText("ListKG", "globalenv()", S=1)));   
  eval(parse(text=SetGText("ListPG", "globalenv()", S=1)));   
}
GetInfoForTables <- function(OnDir="", TSEDir =  "", SecondsMinutesHours="Seconds") {
if (is.character(SecondsMinutesHours)) {
  if (SecondsMinutesHours %in% c("Seconds", "seconds")) {
    SecondsMinutesHours <- 1;
  } else if  (SecondsMinutesHours %in% c("Minutes", "minutes")) {
    SecondsMinutesHours <- 60;
  } else if (SecondsMinutesHours %in% c("Hours", "hours")) {
    SecondsMinutesHours <- 60 * 60;
  } else if (SecondsMinutesHours %in% c("Days", "days")) {
    SecondsMinutesHours <- 60 * 60 * 24;
  }
  if (is.numeric(SecondsMinutesHours)) {
    if (SecondsMinutesHours == 2) {
      SecondsMinutesHours = 60;
    } else if (SecondsMinutesHours==3) {
      SecondsMinutesHours = 60 * 60
    } else if (SecondsMinutesHours == 4) {
      SecondsMinutesHours = 60 * 60 * 24;
    }
  }
} 

if (TSEDir == "") {
  try(library(AlanDirectories)); try(SetSaveHomes());
  try(TSEDir <-  paste(DropHome, "//TwoSimExample//", sep=""));
}
if (OnDir == "") {
  print("ERROR: GeTInfoForTables, OnDir is not a directory");
  return;
}

##OnDir <- "N1000p1ef05k6sig2piSV3sigSV3Rhop9";
MyFiles <- unlist(list.files(paste(TSEDir, "//", OnDir, sep="")));


SFiles <- MyFiles[substr(MyFiles,1,nchar("S")) == "S"];
MIPFiles <- MyFiles[substr(MyFiles,1,nchar("MIP")) == "MIP"];
CIFiles <- MyFiles[substr(MyFiles,1,nchar("CI")) == "CI"];

if (length(SFiles) <= 0) {
   return(-1);
}

MySKnown <- substr(SFiles, 1+nchar("S"), nchar(SFiles)-nchar("Summary.RData"));
MyCIKnown <- substr(CIFiles, 1+nchar("CI"), nchar(CIFiles)-nchar("Summary.RData"));
MyMIPKnown <- substr(MIPFiles, 1+nchar("MIP"), nchar(MIPFiles)-nchar("Summary.RData"));

SigmaNoise <- NULL;  rm(SigmaNoise);  UseRho <- NULL; rm(UseRho);
load(paste(TSEDir, "//", OnDir, "//", SFiles[1], sep=""));
nSimulation <- -1;
pSimulation <- -1;
if (exists("n") && !is.null(n) && length(n) == 1 && is.numeric(n)) {             
  nSimulation <- n;  eval(parse(text=SetGText("nSimulation", "globalenv()", S=1)));
}  else {
  IrD <- (strsplit(OnDir, "N")[[1]])[2];
  IrD <- (strsplit(IrD, "p")[[1]])[1];
  AV <- (strsplit(IrD, "e"));
  if (length(AV[[1]]) == 1) {
    try(nSimulation <- as.numeric(AV[[1]][1]));  
    eval(parse(text=SetGText("nSimulation", "globalenv()", S=1)));
  } else if (length(AV[[1]]) == 2) {
    Ones = NULL;  Tens = NULL;
    try(Ones <- as.numeric(AV[[1]][1]));
    try(Tens <- as.numeric(AV[[1]])[2]);
    if (is.numeric(Ones) && is.numeric(Tens)) {
      n <- Ones * 10^(Tens);  nSimulation <- n;
      eval(parse(text=SetGText("nSimulation", "globalenv()", S=1)));
    }
  }
}
nSimulation <- -1;
if (exists("p") && !is.null(p) && length(p) == 1 && is.numeric(p)) {
  pSimulation <- p;  eval(parse(text=SetGText("pSimulation", "globalenv()", S=1)));
}  else {
  IrD <- (strsplit(OnDir, "p")[[1]])[2];
  IrD <- (strsplit(IrD, "k")[[1]])[1];
  AV <- (strsplit(IrD, "e"));
  if (length(AV[[1]]) == 1) {
    try(pSimulation <- as.numeric(AV[[1]][1]));  
    eval(parse(text=SetGText("pSimulation", "globalenv()", S=1)));
  } else if (length(AV[[1]]) == 2) {
    Ones = NULL;  Tens = NULL;
    try(Ones <- as.numeric(AV[[1]][1]));
    try(Tens <- as.numeric(AV[[1]])[2]);
    if (is.numeric(Ones) && is.numeric(Tens)) {
      p <- Ones * 10^(Tens);  pSimulation <- p;
      eval(parse(text=SetGText("pSimulation", "globalenv()", S=1)));
    }
  }
}
kSimulation <- -1;
if (exists("k") && !is.null(k) && length(k) == 1 && is.numeric(k)) {
  kSimulation <- k;  eval(parse(text=SetGText("kSimulation")));
}  else {
  IrD <- (strsplit(OnDir, "k")[[1]])[2];
  IrD <- (strsplit(IrD, "sig")[[1]])[1];
  AV <- (strsplit(IrD, "e"));
  if (length(AV[[1]]) == 1) {
    try(kSimulation <- as.numeric(AV[[1]][1]));  
    eval(parse(text=SetGText("kSimulation", "globalenv()", S=1)));
  } else if (length(AV[[1]]) == 2) {
    Ones = NULL;  Tens = NULL;
    try(Ones <- as.numeric(AV[[1]][1]));
    try(Tens <- as.numeric(AV[[1]])[2]);
    if (is.numeric(Ones) && is.numeric(Tens)) {
      k <- Ones * 10^(Tens);  kSimulation <- k;
      eval(parse(text=SetGText("kSimulation", "globalenv()", S=1)));
    }
  }
}
if (exists("sigma") && exists("CovRho")) {
  SigmaNoise <- sigma;  UseRho <- CovRho;
  eval(parse(text=SetGText("sigma", "globalenv()", S=1)));  
  eval(parse(text=SetGText("CovRho", "globalenv()", S=1)));
  eval(parse(text=SetGText("SigmaNoise", "globalenv()", S=1)));  
  eval(parse(text=SetGText("UseRho", "globalenv()", S=1)));
} else if (exists("SigmaNoise") || !exists("SigmaNoise")) {
   IrD <- strsplit(OnDir, "sig")[[1]];   try(IrD <- IrD[2]);
   try(IrD <- strsplit(IrD, "pi")[[1]]);
   LS <- IrD[1];  try(IS <- strsplit(IrD, "p")[[1]]);
   if (length(IS) <= 1) {
     try(SigmaNoise <- as.numeric(IrD[1]));
   } else {
     try(AFT <- as.numeric(IS[2]));
     if (log(AFT,10) == round(log(AFT,10))) {
        AFT <- AFT / 10^(log(AFT,10)+1);
     } else {
        AFT <- AFT / 10^(ceiling(log(AFT,10)));
     }
     try(SigmaNoise <- as.numeric(IS[1]) + AFT);
   }
   eval(parse(text=SetGText("SigmaNoise", "globalenv()", S=1)));
   IrD <- strsplit(OnDir, "Rho")[[1]];  IrD <- IrD[2];
   IS <- strsplit(IrD, "p")[[1]];
   if (length(IS) <= 1) {
     try(UseRho <- as.numeric(IrD[1]));
   } else {
     try(AFT <- as.numeric(IS[2]));
     if (log(AFT,10) == round(log(AFT,10))) {
        AFT <- AFT / 10^(log(AFT,10)+1);
     } else {
        AFT <- AFT / 10^(ceiling(log(AFT,10)));
     }
     if (IS[1] == "") { IS[1] = "0"; }
     try(UseRho <- as.numeric(IS[1]) + AFT);
   }
   eval(parse(text=SetGText("UseRho", "globalenv()", S=1)));
}

MeanSum <- matrix(NA, length(SFiles), NCOL(MyOutMatrix)+1);
SDSum <- MeanSum+0.0;
rownames(MeanSum) <- MySKnown;
rownames(SDSum) <- MySKnown;
if (length(colnames(MyOutMatrix)) != length(MyOutMatrix)) {
  RN <- colnames(MyOutMatrix);
} else {
  RN <- names(MyOutMatrix);
}
    colnames(MeanSum) <- c(RN, "N");
    colnames(SDSum) <- colnames(MeanSum);
n <- NULL;  p <- NULL; rm(n);  rm(p);
load(paste(TSEDir, "//", OnDir, "//", SFiles[1], sep=""));
if (exists("n") && !is.null(n) && length(n) == 1 && is.numeric(n)) {
  nSimulation <- n+ 0.0;  
  eval(parse(text=SetGText("nSimulation", "globalenv()", S=1)));
}
if (exists("p") && !is.null(p) && length(p) == 1 && is.numeric(p)) {
  pSimulation <- p+ 0.0;  
  eval(parse(text=SetGText("pSimulation", "globalenv()", S=1)));
}
for (ii in 1:length(SFiles)) {
  load(paste(TSEDir, "//", OnDir, "//", SFiles[ii], sep=""));
  UseOutMatrix <- MyOutMatrix[MyOutMatrix[,3] != -2,];
  if (length(UseOutMatrix) >= 1) {
  if (NCOL(UseOutMatrix) == 1 ||NROW(UseOutMatrix) == 1) {
  SDSum[ii,1:length(UseOutMatrix)] <-   NA;
  SDSum[ii,NCOL(SDSum)] <- 1; 
  MeanSum[ii,1:length(UseOutMatrix)] <- UseOutMatrix;
  MeanSum[ii,NCOL(MeanSum)] <- 1;   
  } else {
  MeanSum[ii,1:NCOL(UseOutMatrix)] <- colMeans(UseOutMatrix, na.rm=TRUE);
  MeanSum[ii,NCOL(MeanSum)] <- NROW(UseOutMatrix);
  SDSum[ii,1:NCOL(UseOutMatrix)] <-   apply(UseOutMatrix, 2, sd);
  SDSum[ii,NCOL(SDSum)] <- NROW(UseOutMatrix);
  }
  }
  UseOutMatrix <- NULL;   MyOutMatrix <- NULL;
}



eval(parse(text=GetG0Text("FirstKeepColumns", envir="globalenv()", S=1)));
if (!exists("FirstKeepColumns") || is.null(FirstKeepColumns) || length(FirstKeepColumns) <= 0) {
FirstKeepColumns <- c("piAM2Lasso", "piAM2LassoNI", "2LassoCV", "ENetFixed", "ENetFixedNI", 
  "ENetCV", "SCADminConvex", "SCADMinLRY", 
  "LarsCp",  "MCPCVMR", "MCPMinR",
  "GenerateBayesSpikeAutoMean", "GenerateBayesSpikeInfoPrior",
    "GenerateBayesSpikeInfoPriorNI", 
  "IshwaranSpike", "BayesLasso",   "BayesVarSel")
}

MyInterestingRows <- match(FirstKeepColumns, rownames(MeanSum));
KeepCols <- c(2,3,4,8);
MyMean <- MeanSum[MyInterestingRows, KeepCols];
rownames(MyMean) <- FirstKeepColumns;
MyMean[,4] = MyMean[,4] / SecondsMinutesHours;
MySDSum <- SDSum[MyInterestingRows, KeepCols];
rownames(MySDSum) <- FirstKeepColumns
MySDSum[,4] = MySDSum[,4] / SecondsMinutesHours;
eval(parse(text=SetGText("MeanSum", "globalenv()", S=1)));
eval(parse(text=SetGText("SDSum", "globalenv()", S=1)));

TimeString <- "time-sec"
if (SecondsMinutesHours == 60) {
  TimeString <- "time-min";
} else if (SecondsMinutesHours == 60 * 60) {
  TimeString <- "time-hour";
} else if (SecondsMinutesHours == 60 * 60 * 24) {
  TimeString <- "time-day";
}
eval(parse(text=SetGText("TimeString", "globalenv()", S=1)));
eval(parse(text=SetGText("MyMean", "globalenv()", S=1)));
eval(parse(text=SetGText("MySDSum", "globalenv()", S=1)));

load(paste(TSEDir, "//", OnDir, "//", MIPFiles[1], sep=""));


Needs <- c("SphereLoss", "AUC", "MarginalHellinger",   "RTIF");   
MeanMIP <- matrix(NA, length(MIPFiles), 4);
SDMIP <- MeanMIP+0.0;  NSum = MeanMIP + 0.0;
rownames(MeanMIP) <- MyMIPKnown;
rownames(SDMIP) <- MyMIPKnown;
rownames(NSum) <- MyMIPKnown;
for (ii in 1:length(MIPFiles)) {
  load(paste(TSEDir, "//", OnDir, "//", MIPFiles[ii], sep=""));
  for (jj in 1:4) {
    MyText <- paste("NNeed <- AllMIPList$", Needs[jj], sep="");
    eval(parse(text=MyText));
    MeanMIP[ii,jj] <- mean(NNeed[!is.na(NNeed)]); 
    NSum[ii,jj] <- length(NNeed[!is.na(NNeed)]);
    SDMIP[ii,jj] <- sd(NNeed[!is.na(NNeed)]);
  }

}
colnames(MeanMIP) = Needs
colnames(SDMIP) = Needs
colnames(NSum) = Needs
MeanMIP[,colnames(MeanMIP)=="SphereLoss"] <- (p-
   MeanMIP[,colnames(MeanMIP)=="SphereLoss"])/p;
SDMIP[,colnames(SDMIP)=="SphereLoss"] <- (
   SDMIP[,colnames(SDMIP)=="SphereLoss"])/p;
eval(parse(text=SetGText("MeanMIP", "globalenv()", S=1)));
eval(parse(text=SetGText("SDMIP", "globalenv()", S=1)));


MyInterestingRows <- match(FirstKeepColumns, rownames(MeanSum));

MIPInterestingRows <- match(FirstKeepColumns, rownames(MeanMIP));
KeepCols <- c(1,2,3);
MyMIP <- MeanMIP[MyInterestingRows, KeepCols];
rownames(MyMIP) <- FirstKeepColumns;
MySDMIP <- SDMIP[MyInterestingRows, KeepCols];
rownames(MySDMIP) <- FirstKeepColumns
eval(parse(text=SetGText("MySDMIP", "globalenv()", S=1)));
eval(parse(text=SetGText("MySDMIP", "globalenv()", S=1)));

AMatNeed <- matrix(NA, length(FirstKeepColumns), 4 + 3 + 1);
AMatSD <- AMatNeed + 0.0;  
  
MyInterestingRows <- match(FirstKeepColumns, rownames(MeanSum));
KeepCols <- c(2,3,4,8);
AMatNeed[,1:4] <- MeanSum[MyInterestingRows, KeepCols];
rownames(AMatNeed) <- FirstKeepColumns;
colnames(AMatNeed) <- c(colnames(MeanSum)[KeepCols], colnames(MeanMIP)[1:3], "n");
AMatNeed[,4] = AMatNeed[,4] / SecondsMinutesHours;
AMatNeed[,8] <- MeanSum[MyInterestingRows,9];
AMatNeed[,5:7] <-MeanMIP[MIPInterestingRows, c(1,2,3)]
AMatSD[,1:4] <- SDSum[MyInterestingRows, KeepCols];
AMatSD[,5:7] <- SDMIP[MIPInterestingRows,c(1,2,3)];
AMatSD[,8] <- MeanSum[MyInterestingRows,9];
rownames(AMatSD) <- FirstKeepColumns
AMatSD[,4] = SDSum[MyInterestingRows,4] / SecondsMinutesHours;
colnames(AMatSD) <- colnames(AMatNeed);
eval(parse(text=SetGText("AMatSD", "globalenv()", S=1)));
eval(parse(text=SetGText("AMatNeed", "globalenv()", S=1)));
return(1);
}


AllKnownNameFunctions <- c("2LassoCV", "2LassoSpreadNI", "2LassoSpread", "B2Lasso", "ENetFixedNI",
"ENetFixed", "ENetCV", "FDLimitLasso", "FDLimitLassoNI",
"GenerateBayesSpikeAutoMean",   
"GenerateBayesSpikeAutoPrior", 
"GenerateBayesSpikeInfoMeanNI", "GenerateBayesSpikeInfoMean", "GenerateBayesSpikeInfoPriorNI",
"GenerateBayesSpikeInfoPrior", "IshwaranSpike", "LarsCpNI", "LarsCp", "LarsFixedNI", "LarsFixed",
"LassoLYNI", "LassoLY", "LassoStoddenNI", "LassoStodden", "LimitLassoNI", "LimitLasso", 
"piALimitLassoNI", "piALimitLasso", "piAM2LassoNI", "piAM2Lasso", "SCAD4xStoddenNI", "SCAD4xStodden",
"SCADminConvex", "SCADMinLRY", "HorseShoe", "BayesLasso", "MCPCVMR", "MCPMinR", "BayesVarSel");

AllKnownFormat <- c("2-Lasso CV", "2-Lasso - $\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}$ = k-Noise/p",
  "2-Lasso - $\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}$ = k-True/p",
  "B2Lasso", "Elastic Net, fix k-true", "Elastic Net, fix k-Noise", "Elastic Net, CV",
  "FD-LimLasso, k-True", "FD LimLasso, k-Noise",
  "GB (1,p) mean", "GB (1,p) median", "GB (k-Noise,p) mean", 
  "GB (k-True,p) mean", "GB (k-Noise,p)  median", "GB (k-True,p) median",
  "Ishwaran Spike\\&Slab", "Lars Cp (sig-Noise)", "Lars Cp (sig)",
  "Lars fix k-Noise", "Lars fix k", "Lasso LY k-noise", "Lasso LY",
  "Lasso Stodden k-noise", "Lasso Stodden k-True", "LimLasso $\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}$ = k-Noise/p",
  "LimLasso $\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}$ = k-True/p",
  "LimLasso (k-Noise,p)", "LimLasso (k-True,p)",
  "2-Lasso (k-Noise,p)", "2Lasso (k-True,p)", "SCAD 4xStodden(k-Noise)", "SCAD 4xStodden(kTrue)",
  "SCAD minConvex", "SCAD LRY", "HorseShoe", "BayesLasso", "MCP CVMR", "MCP MinR", "BayesVarSel");  
