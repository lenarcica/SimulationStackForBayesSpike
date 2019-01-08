
DeclareUniqueSimulationIdentifier <- function(ExperimentName = "", jobii = 0,
  WorkingRow = 0,n,p,k,sigma,LogitRegression=0) {
  if (ExperimentName == "") {
    eval(parse(text=GetG0Text("ExperimentName", S=1)));
    if (ExperimentName == 0) { ExperimentName = "" }
  }
  if (jobii[1] == 0) {
    eval(parse(text=GetG0Text("jobii", S=1))); 
  }
  if (WorkingRow == 0) {
     eval(parse(text=GetG0Text("WorkingRow", S=1))); 
  }
  ARealTime <- GetRealTime();
  UniqueSimulationIdentifier = paste(ExperimentName,"J", jobii, "R", WorkingRow,
   "n", n, "p", p, "k", k, "s", tSeq(round(sigma,1)),
   "d", ARealTime$JustDays, "s", tSeq(round(ARealTime$JustSeconds)),  sep="");
  if (LogitRegression == 1) {
  UniqueSimulationIdentifier = paste(ExperimentName,"J", jobii, "R", WorkingRow,
   "n", n, "p", p, "k", k, "s", tSeq(round(sigma,1)),
   "d", ARealTime$JustDays, "s", tSeq(round(ARealTime$JustSeconds)), "L1", sep="");    
  }
  return(UniqueSimulationIdentifier);
}

if (!exists("SimCor")) {
  SimCor <- function(p, rho) {
    RetVector <- rnorm(p);  ACons <- sqrt(1-rho^2)
    for (tt in 2:p) {
      RetVector[tt] = RetVector[tt-1] * rho + RetVector[tt] * ACons;
    }
    return(SimCor);
  }
}

DeclareUniqueGroupSimulationIdentifier <- function(ExperimentName = "", jobii = 0,
  WorkingRow = 0,n,pGroups, GroupSize, kGroups,sigma, LogitRegression = 0) {
  if (ExperimentName == "") {
    eval(parse(text=GetG0Text("ExperimentName", S=1)));
    if (ExperimentName == 0) { ExperimentName = "" }
  }
  if (jobii[1] == 0) {
    eval(parse(text=GetG0Text("jobii", S=1))); 
  }
  if (WorkingRow == 0) {
     eval(parse(text=GetG0Text("WorkingRow", S=1))); 
  }
  ARealTime <- GetRealTime();
  UniqueSimulationIdentifier = paste(ExperimentName,"GrJ", jobii, "R", WorkingRow,
   "n", n, "pG", pGroups, "Gs", GroupSize, "kGroups", kGroups, "s", tSeq(round(sigma,1)),
      "d", ARealTime$JustDays, "s", tSeq(round(ARealTime$JustSeconds)),  sep="");
  if (LogitRegression == 1) {
   UniqueSimulationIdentifier = paste(ExperimentName,"GrJ", jobii, "R", WorkingRow,
   "n", n, "pG", pGroups, "Gs", GroupSize, "kGroups", kGroups, "s", tSeq(round(sigma,1)),
      "d", ARealTime$JustDays, "s", tSeq(round(ARealTime$JustSeconds)), "L1", sep="");   
    
  }
  return(UniqueSimulationIdentifier);
}
AZeroOut <- function(ANum=0, maxNum=1000) {
  TrySetupText <- "
  if (is.null(ANum)) {
    print(paste(\"AZeroOut: Well you gave us NULL ANum!\", sep=\"\"));
    BreakMe1 = \"-2 = -1;\"
    eval(parse(text=BreakMe1));
  } 
  ";
  try(eval(parse(text=TrySetupText)));
  TrySetupText <- "
  if (ANum < 0) {
    print(paste(\"AZeroOut: No we do not Zero out when ANum = \", ANum, sep=\"\"));
    flush.console();
    BreakMe1 = \"-2 = -1;\"
    eval(parse(text=BreakMe1));
  }
  ";
  try(eval(parse(text=TrySetupText)));
  TrySetupText <- "
  if (!is.numeric(ANum)) {
    print(paste(\"AZeroOut: No, ANum = \", ANum, \" is not a numeric! \", sep=\"\"));flush.console();
    BreakMe1a = \"-2=-1;\"
    eval(parse(text=BreakMe1a));
  }
  " ;
  try(eval(parse(text=TrySetupText)));
  TrySetupText <- "
  if (!is.numeric(maxNum)) {
    print(paste(\"AZeroOut: No, maxNum = \", maxNum, \" is not a numeric! \", sep=\"\"));flush.console();
    BreakMe1b = \"-2=-1;\"
    eval(parse(text=BreakMe1a));  
  }
  if (maxNum <= 0) {
    print(paste(\"AZeroOut: No, maxNum = \", maxNum, \" is less than zero! \", sep=\"\")); flush.console();
    BreakMe2 = \"-2 = -1;\"
    eval(parse(text=BreakMe2));
  }
  ANum <- ANum[1];  maxNum <- maxNum[1];
  ";
  try(eval(parse(text=TrySetupText)));
  lmaxNum <- ceiling(log(maxNum, 10))-1;
  
  if (ANum == 0) { return(paste(rep("0", lmaxNum), collapse="")); }
  if (ANum > maxNum) { return(10^(lmaxNum+1)-1); }
  lANum <- floor(log(ANum, 10));
  if (lANum > lmaxNum) {
    print(paste("AZeroOut: What is Wrong: you have ANum = ", ANum, " but maxNum = ", maxNum, sep="")); flush.console();
    -2 = -1;
  }
  if (lANum == lmaxNum) { return(paste(ANum , sep="")); }
  return(paste(paste(rep("0", lmaxNum-lANum), collapse=""), ANum, sep=""));
}

DeclareUniqueProcessIdentifier <- function(NInstance=0, TimeRun =-1) {
  if (TimeRun == -1) {
    ADate <- unlist(strsplit(as.character(Sys.Date()), "-"));
    ADate <- ADate[ADate != ""];
    if (is.null(ADate) || length(ADate) < 3) {
      ADate <- c("6", "66", "6666");
    }  else {
      ADate <- c("6", "66", "6666");
      try(ADate <- paste(ADate[1], AZeroOut(as.numeric(ADate[2]),100), AZeroOut(as.numeric(ADate[3]),100), sep=""));
    }
    ADT <- unlist(strsplit(date(), " "));
    ADT <- ADT[ADT != ""];
    GGZ <- strsplit(ADT, ":");
    lGZ <- rep(0, length(GGZ));
    for (ii in 1:length(GGZ)) {
      lGZ[ii] <- length(GGZ[[ii]]);
    } 
    AOn <- (1:length(lGZ))[lGZ >= 3];
    ##ABT <- unlist(strsplit(ADT[4], ":"));
    ABT <- GGZ[[AOn]];
    ABT <- ABT[ABT != ""];
    NewDate <- paste(ADate, AZeroOut(as.numeric(ABT[1]), 100),
      AZeroOut(as.numeric(ABT[2]), 100), AZeroOut(as.numeric(ABT[3]), 100), sep="")
    tt = proc.time();
    tt = paste( tt[1:2], collapse="");
    TimeRun = tSeq(tt); 
    for (ii in 1:5) { TimeRun = tSeq(TimeRun);}
    TimeRun <- paste(NewDate, "PP", TimeRun, sep="");
  }
  RR1 <- AZeroOut(sample(0:999, size=1),1000);
  RR2 <- AZeroOut(sample(0:999, size=1),1000);
  RR3 <- AZeroOut(sample(0:999, size=1),1000);
  if (!exists("NInstance") || is.null(NInstance) || is.character(NInstance) ||
    is.numeric(NInstance)  && NInstance <= 0) {
    NInstance <- 999;  
  }
  AZN <- "FFFFF";
  try(AZN <- AZeroOut(NInstance, 100000));
  if (AZN == "FFFFF") {
     print("Hey: DeclareUniqueProcessIdentifier, we still can't transcribe NInstance!");
     flush.console();
     try(print(paste("NInstance is ", NInstance, sep=""))); flush.console();
     print("We are going with default NInstance is FFFFF");
  }
  UniqueProcessIdentifier = paste("N", AZN[1], "T", TimeRun[1], 
    "RR1", RR1[1], "RR2", RR2[1], "RR3", RR3[1], sep="");
  if (length(UniqueProcessIdentifier) > 1) {
    UniqueProcessIdentifier <- UniqueProcessIdentifier[1];
  }
  eval(parse(text = SetGText("UniqueProcessIdentifier","globalenv()", S=1)));
  return(UniqueProcessIdentifier);
}


##############################################################
##   SimMeData()
##    This function simulates a sparse linear regression dataset
##    The covariates will be correlated.
##    Noise is Sigma noise level.  
##    BRealVecs is vector of first set of non-zero values
##    ppCov,sigmaFactCov are parameters for correlating covariates.
##
##
SimMeData <- function(n = -999, p = -999, k = -999, sigma = -999,
     NLen = -999, kLen = -999, kActiveLen = -999,  ElseGetX = NULL, 
     SigmaNoise = -999, GenerateBetaVec = DefaultGenerateBetaVec, 
     CorrelationXmatrix = DefaultCorrelationXmatrix, tNoiseDF = DefaultNoiseDF,
     LogitRegression = FALSE, LogitFactor = 1, Beta0 = 0, ExperimentName="", jobii= 0, 
     WorkingRow = 0, AlreadyLocked=TRUE, MeanCenterXX = TRUE,
     UniqueProcessIdentifier = "", ISample = 1,
     DontSave=FALSE, InOrder=FALSE, dfRobit = -1) {
  if (ExperimentName == "") {
    eval(parse(text=GetG0Text("ExperimentName", S=1)));
    if (ExperimentName == 0) { ExperimentName = "" }
  }
  if (jobii[1] == 0) {
    eval(parse(text=GetG0Text("jobii", S=1))); 
  }
  if (WorkingRow == 0) {
     eval(parse(text=GetG0Text("WorkingRow", S=1))); 
  }
 UniqueSimulationIdentifier = DeclareUniqueSimulationIdentifier(ExperimentName = ExperimentName, 
   jobii = jobii, WorkingRow = WorkingRow,n,p,k,sigma, as.integer(LogitRegression));
 if (kLen == -999) { if (p > 0) {kLen = p;} else{ kLen = 20;} }
 if (NLen == -999) { if (n > 0) {NLen = n;} else{ NLen = round(kLen * .8);} }
 if (!exists("k")) {
   print("k Doesn't exist!"); flush.console();  
 }                   
 if (kActiveLen == -999) {
    if (k > 0) {kActiveLen = k} else { kActiveLen = round(kLen * .4); }  
 }
 if (SigmaNoise == -999) {
   if (sigma > 0) { SigmaNoise = sigma } else { SigmaNoise <- .5; }
 }
  BetasRealVD <- GenerateBetaVec(kLen, kActiveLen, InOrder=InOrder);
  BetasReal <- BetasRealVD[[1]];
  BBsReal <- BetasReal; BBsReal[abs(BetasReal) > 0] <- 1;

  DGG <- CorrelationXmatrix(4);
  
  NameX <- NULL;
  if (!is.null(ElseGetX) && NROW(ElseGetX) == NLen && NCOL(ElseGetX) == kLen) {
    XX <- ElseGetX;
    NameX <- colnames(ElseGetX); 
  } else if (DGG$type== "CorrelationXmatrix1" && kLen > 250)  {
    XX <- matrix( rnorm(NLen*kLen,0,1), kLen, NLen );
    Rem <- sqrt(1-DGG$PosCorr^2);
    for (ii in 2:kLen) {
      XX[ii,] <- XX[ii-1,] * DGG$PosCorr + Rem * XX[ii,]; 
    }
  } else {
    IIMax <- 30;  iiSim = 1;
    while( iiSim <= IIMax ) {
      CorrRT <- CorrelationXmatrix(kLen);
      vD <- eigen(CorrRT[[1]])$values;
      if (min(vD) > 0) {
        break;
      }
      iiSim = iiSim+ 1;
      if (iiSim >= IIMax) {
        print(paste("SimMeData, eroor could not get positive definite",
          " Cor matrix in ", IIMax, " tries. ", sep="") );
        return(-1);
      }
    }
    CorrMat = CorrRT[[1]];
    LCovX <- t(chol(CorrMat));
    XX <- LCovX %*% matrix( rnorm( NLen * kLen, 0, 1), kLen, NLen );
  }
  if (MeanCenterXX == TRUE) {
    XX <- XX - rowMeans(XX);
  }
  XX <- t(XX);
  if (LogitRegression == TRUE && LogitFactor > 1) {
    BetasReal <- BetasReal * LogitFactor;
  }
   
  SigmaSqNoise = SigmaNoise^2;
  if (tNoiseDF <= 0) {
    YY <-    XX %*% BetasReal + SigmaNoise * rnorm(NLen, 0,1 );
  } else {
    YY <-   XX %*% BetasReal + SigmaNoise * rnorm(NLen, 0, 1) / 
      sqrt( rchisq(NLen, tNoiseDF) / tNoiseDF );
  }
  if (MeanCenterXX == TRUE) {
    YY <- YY - mean(YY);
  }
  BADUniqueProcessIdentifier <- 1;
  ATryGo <- "
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
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
    eval(parse(text=GetG0Text("UniqueProcessIdentifier")));
  }
  UniqueProcessThatCreatedMeIdentifier =  UniqueProcessIdentifier;
  if (dfRobit >= 0.0) {
    Beta0 * rnorm(1,0,1) * SigmaNoise;
    YY = XX %*% BetasReal + Beta0;
    if (dfRobit == 0.0) {
      pYY <- pnorm(YY/SigmaNoise,0,1);
    } else {
      pYY = pt(YY/SigmaNoise, df=dfRobit)
    }
    YY = rbinom( length(pYY), 1, pYY); 
  } else if (LogitRegression == TRUE) {
    Beta0 * rnorm(1,0,1) * SigmaNoise;
    YY = XX %*% BetasReal + Beta0;
    pYY = exp(YY) / ( 1 + exp(YY))
    YY = rbinom( length(pYY), 1, pYY); 
    ##XX = cbind(rep(1, length(YY)), XX);
  }

    eval(parse(text=GetG0Text("piSV", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("sigSV", "globalenv()", S=1)));
    kActiveNoise = round(length(BetasReal[BetasReal!=0]) + sqrt(NCOL(XX)) * piSV * rnorm(1, 0, 1))
    if (kActiveNoise <= 0) { kActiveNoise <- 1; }
    if (kActiveNoise >= length(BetasReal)) { kActiveNoise <- length(BetasReal)-1; }
    if (kActiveNoise >= NROW(XX)) { kActiveNoise <- NROW(XX)-1; }

	 puse <- k / p;
	 if (puse > 1) { puse = (p-1)/p;}
	 if (puse < 0) {puse = 1/p;}
	 puseNoise <- kActiveNoise / p; 
	 
	 ##eval(parse(text=GetG0Text("piSV", S=1)));

    if (sigSV > 0) { 
      SigmaNoisesim <- sqrt(SigmaNoise^2 * rchisq(1, sigSV)/sigSV)
    } else {
       SigmaNoisesim <- (SigmaNoise)
    }
    eval(parse(text=GetG0Text("piSV", "globalenv()", S=1)));
    eval(parse(text=GetG0Text("sigSV", "globalenv()", S=1)));
    
    eval(parse(text=GetGText("PriorStrengthFunction")));
    eval(parse(text=GetGText("PriorStrength")));
    eval(parse(text=GetGText("DefaultPriorStrengthFunction")));
    if (is.null(PriorStrength)) {
      if (!is.null(PriorStrengthFunction)) {
         PriorStrength <- PriorsStrengthFunction(n=n,p=p, pGroups=NULL); 
      } else {
         PriorStrength <- DefaultPriorStrengthFunction(n=n,p=p, pGroups=NULL)
      }
    }   
    eval(parse(text=GetGText("SigmaPriorStrengthFunction")));
    eval(parse(text=GetGText("SigmaPriorStrength")));
    eval(parse(text=GetGText("DefaultSigmaPriorStrengthFunction")));
    if (is.null(SigmaPriorStrength)) {
      if (!is.null(SigmaPriorStrengthFunction)) {
         SigmaPriorStrength <- SigmaPriorsStrengthFunction(n=n,p=p, pGroups=pGroups); 
      } else {
         SigmaPriorStrength <- DefaultSigmaPriorStrengthFunction(n=n,p=p, pGroups=pGroups)
      }
    } 
    if (!is.numeric(SigmaPriorStrength) || length(SigmaPriorStrength) != 1) {
      SigmaPriorStrength = 0;  
    }    
    SMSNoise = list(XX=XX, YY=YY, puse = puseNoise, 
          BetasReal = BetasReal, BBsReal = BBsReal,
          n=n,p=p,k=k,   ElseGetX=ElseGetX,
          NLen = NLen, kLen = kLen, kActiveLen = kActiveNoise,
          SigmaSqNoise = SigmaNoisesim^2,  SigmaNoise = SigmaNoisesim,
          SigmaNoiseReal = SigmaNoise, kActiveRealLen = k,
          BetasName = GetBetaName(GenerateBetaVec),
          CorrName = GetCorrName(CorrelationXmatrix),
          LogitRegression=LogitRegression, dfRobit = dfRobit,
          tNoiseDF = tNoiseDF, SigmaPriorStrength = SigmaPriorStrength,
          UniqueSimulationIdentifier=UniqueSimulationIdentifier,
          UniqueProcessThatCreatedMeIdentifier = UniqueProcessThatCreatedMeIdentifier,
          piSV=piSV, sigSV=sigSV, PriorStrength=PriorStrength);      
   MyRet = list(XX=XX, YY=YY, puse = puse, 
          BetasReal = BetasReal, BBsReal = BBsReal,
          n=n,p=p,k=k,  NameX = NameX,
          NLen = NLen, kLen = kLen, kActiveLen = k,
          SigmaSqNoise = SigmaNoise^2,  SigmaNoise = SigmaNoise,
          SigmaNoiseReal = SigmaNoise, kActiveRealLen = k,
          BetasName = GetBetaName(GenerateBetaVec),
          CorrName = GetCorrName(CorrelationXmatrix),
          LogitRegression=LogitRegression, dfRobit = dfRobit,
          tNoiseDF = tNoiseDF,  SigmaPriorStrength = SigmaPriorStrength,
          UniqueSimulationIdentifier=UniqueSimulationIdentifier,
          UniqueProcessThatCreatedMeIdentifier = UniqueProcessThatCreatedMeIdentifier,
          SMSNoise = SMSNoise, piSV=piSV, sigSV=sigSV, PriorStrength=PriorStrength);
   if (!exists("DontSave") || 
     is.null(DontSave) || !is.logical(DontSave) || DontSave == FALSE) {
     SaveCurrentNewSimulation(MyRet, UniqueSimulationIdentifier=UniqueSimulationIdentifier, 
       AlreadyLocked=AlreadyLocked, 
       UniqueProcessIdentifier = UniqueProcessIdentifier, 
       ISample = ISample);
   }
   return(MyRet);
}

DoubleSimulation <- function(SMS, AlreadyLocked, DontSave=FALSE) {
  
  if (is.null(SMS$puseGroups)) {
      ##attach(SMS$SMSNoise);
      SMN <- SMS$SMSNoise;
      NewSMSoise <-   list(XX=cbind(SMN$XX,SMN$XX), YY=SMN$YY, puse = SMN$puse, 
          BetasReal = c(SMN$BetasReal, SMN$BetasReal)/2, BBsReal = c(SMN$BBsReal, SMN$BBsReal),
          BetasPartReal = c(SMN$BetasReal, rep(0, length(SMN$BetasReal))),
          BetasOtherReal = c(rep(0, length(SMS$BetasReal)), SMN$BetasReal),
          n=SMN$n,p=NCOL(SMN$XX)*2,k=SMN$k*2,  ElseGetX = paste(SMN$ElseGetX, "x2", sep=""),
          NLen = SMN$NLen, kLen = SMN$kLen*2, kActiveLen = SMN$k*2,
          SigmaSqNoise = SMN$SigmaNoise^2,  SigmaNoise = SMN$SigmaNoise,
          SigmaNoiseReal = SMN$SigmaNoise, kActiveRealLen = SMN$k*2,
          BetasName = SMN$BetasName,
          CorrName = SMN$CorrName,
          LogitRegression=SMN$LogitRegression, dfRobit = SMN$dfRobit,
          tNoiseDF = SMN$tNoiseDF,  SigmaPriorStrength = SMN$SigmaPriorStrength,
          UniqueSimulationIdentifier=paste("DD", SMN$UniqueSimulationIdentifier, sep=""),
          UniqueProcessThatCreatedMeIdentifier = SMN$UniqueProcessThatCreatedMeIdentifier,
          SMSNoise = NULL, piSV=SMN$piSV, sigSV=SMN$sigSV, PriorStrength=SMN$PriorStrength);   
      ##attach(SMS);
      MyRet = list(XX=cbind(SMS$XX,SMS$XX), YY=SMS$YY, puse = SMS$puse, 
          BetasReal = c(SMS$BetasReal, SMS$BetasReal)/2, BBsReal = c(SMS$BBsReal, SMS$BBsReal),
          n=SMS$n,p=NCOL(SMS$XX)*2,k=SMS$k*2,  NameX = paste(SMS$NameX, "x2", sep=""),
          NLen = SMS$NLen, kLen = SMS$kLen*2, kActiveLen = SMS$k*2,
          BetasPartReal = c(SMS$BetasReal, rep(0, length(SMS$BetasReal))),
          BetasOtherReal = c(rep(0, length(SMS$BetasReal)), SMS$BetasReal),
          SigmaSqNoise = SMS$SigmaNoise^2,  SigmaNoise = SMS$SigmaNoise,
          SigmaNoiseReal = SMS$SigmaNoise, kActiveRealLen = SMS$k*2,
          BetasName = SMS$BetasName,
          CorrName = SMS$CorrName,
          LogitRegression=SMS$LogitRegression, dfRobit = SMS$dfRobit,
          tNoiseDF = SMS$tNoiseDF,  SigmaPriorStrength = SMS$SigmaPriorStrength,
          UniqueSimulationIdentifier=paste("DD", SMS$UniqueSimulationIdentifier, sep=""),
          UniqueProcessThatCreatedMeIdentifier = SMS$UniqueProcessThatCreatedMeIdentifier,
          SMSNoise = SMS$NewSMSoise, piSV=SMS$piSV, sigSV=SMS$sigSV, PriorStrength=SMS$PriorStrength);  
        
      if (!exists("DontSave") || 
       is.null(DontSave) || !is.logical(DontSave) || DontSave == FALSE) {
       try(SaveCurrentNewSimulation(MyRet, UniqueSimulationIdentifier=MyRet$UniqueSimulationIdentifier, 
         AlreadyLocked=AlreadyLocked, UniqueProcessIdentifier = MyRet$UniqueProcessIdentifier, 
         ISample = ISample));
      }
      return(MyRet); 
  }
  ##attach(SMS$SMSNoise);
  SMN <- SMS$SMSNoise
  NewSMSNoise = list(XX=cbind(SMN$XX,SMN$XX), YY=SMN$YY, puse = SMN$puse*2,  puseGroups=SMN$puseGroups*2,
          BetasReal = c(SMN$BetasReal, SMN$BetasReal)/2, BBsReal = c(SMN$BBsReal, SMN$BBsReal),
          NLen = SMN$NLen, kLen = SMN$kLen*2, kActiveLen = 2*SMN$kActiveLen,
          p = NCOL(SMN$XX)*2, n = length(SMN$YY), kGroups=SMN$kGroups*2,
          pGroups=SMN$pGroups*2,  ElseGetX=paste(SMN$ElseGetX, "x2", sep=""),
          BetasPartReal = c(SMN$BetasReal, rep(0, length(SMN$BetasReal))),
          BetasOtherReal = c(rep(0, length(SMN$BetasReal)), SMN$BetasReal),
          SigmaSqNoise = SMN$SigmaNoise^2,  SigmaNoise = SMN$SigmaNoise,
          SigmaNoiseReal = SMN$SigmaNoise, kActiveRealLen = SMN$kActiveLen*2,
          BetasName = SMN$BetasName,
          CorrName = SMN$CorrName,
          LogitRegression=SMN$LogitRegression,  dfRobit=SMN$dfRobit,
          tNoiseDF = SMN$tNoiseDF, SigmaPriorStrength = SMN$SigmaPriorStrength,
          UniqueSimulationIdentifier=paste("DD", SMN$UniqueSimulationIdentifier, sep=""),
          UniqueProcessThatCreatedMeIdentifier = SMN$UniqueProcessThatCreatedMeIdentifier,
          FirstGroupIndex = SMN$FirstGroupIndex, EndGroupIndices=SMN$EndGroupIndices,
          indexIndices=SMN$indexIndices, TrueOnGroups=SMN$TrueOnGroups, SMSNoise = NULL, 
          piSV=SMN$piSV,sigSV=SMN$sigSV, PriorStrength = SMN$PriorStrength, Beta0=SMN$Beta0);
    ##attach(SMS);      
          
    MyRet = list(XX=cbind(SMS$XX,SMS$XX), YY=SMS$YY, puse = SMS$puse*2,  puseGroups=SMS$puseGroups*2,
          BetasReal = c(SMS$BetasReal, SMS$BetasReal)/2, BBsReal = c(SMS$BBsReal, SMS$BBsReal),
          NLen = SMS$NLen, kLen = SMS$kLen*2, kActiveLen = SMS$kActiveLen*2,
          p = NCOL(SMS$XX)*2, n = length(SMS$YY), kGroups=SMS$kGroups*2,
          BetasPartReal = c(SMS$BetasReal, rep(0, length(SMS$BetasReal))),
          BetasOtherReal = c(rep(0, length(SMS$BetasReal)), SMS$BetasReal),
          pGroups=SMS$pGroups*2,  NameX=paste(SMS$NameX, "x2", sep=""),
          SigmaSqNoise = SMS$SigmaNoise^2,  SigmaNoise = SMS$SigmaNoise,
          SigmaNoiseReal = SMS$SigmaNoise, kActiveRealLen = SMS$kActiveLen*2,
          BetasName = SMS$BetasName,
          CorrName = SMS$CorrName,
          LogitRegression=SMS$LogitRegression,  dfRobit=SMS$dfRobit,
          tNoiseDF = SMS$tNoiseDF, SigmaPriorStrength = SMS$SigmaPriorStrength,
          UniqueSimulationIdentifier=paste("DD", SMS$UniqueSimulationIdentifier, sep=""),
          UniqueProcessThatCreatedMeIdentifier = SMS$UniqueProcessThatCreatedMeIdentifier,
          FirstGroupIndex = SMS$FirstGroupIndex, EndGroupIndices=SMS$EndGroupIndices,
          indexIndices=SMS$indexIndices, TrueOnGroups=SMS$TrueOnGroups, SMSNoise = SMS$NewSMSNoise, 
          piSV=SMS$piSV,sigSV=SMS$sigSV, PriorStrength = SMS$PriorStrength, Beta0=SMS$Beta0);
    if (!exists("DontSave") || 
       is.null(DontSave) || !is.logical(DontSave) || DontSave == FALSE) {
       try(SaveCurrentNewSimulation(MyRet, UniqueSimulationIdentifier=MyRet$UniqueSimulationIdentifier, 
         AlreadyLocked=AlreadyLocked, UniqueProcessIdentifier = MyRet$UniqueProcessIdentifier, 
         ISample = ISample));
    }
    return(MyRet);
}



##############################################################
##   SimMeData()
##    This function simulates a sparse linear regression dataset
##    The covariates will be correlated.
##    Noise is Sigma noise level.  
##    BRealVecs is vector of first set of non-zero values
##    ppCov,sigmaFactCov are parameters for correlating covariates.
##
##
SimGroupData <- function(n = -999, pGroups = -999, GroupSize=5, kGroups = -999, sigma = -999,
     SigmaNoise = -999, tNoiseDF = -1, GenerateBetaGroupVec = DefaultGenerateBetaGroupVec, 
     CorrelationXmatrix = DefaultCorrelationXmatrix, ElseGetX = NULL,  ElseTauEndList = NULL,
     LogitRegression = FALSE, LogitFactor = 1, Beta0 = 0, ExperimentName="", jobii= 0, 
     WorkingRow = 0, AlreadyLocked=TRUE, MeanCenterXX = TRUE,
     UniqueProcessIdentifier = "", ISample = 1, GroupsSumToZero=0,
     DontSave=FALSE, InOrder=FALSE, dfRobit=-1) {
  if (!exists("ExperimentName") || ExperimentName == "") {
    eval(parse(text=GetG0Text("ExperimentName", S=1)));
    if (ExperimentName == 0) { ExperimentName = "" }
  }
  if (jobii[1] == 0) {
    eval(parse(text=GetG0Text("jobii", S=1))); 
  }
  if (!exists("WorkingRow") || WorkingRow == 0) {
     eval(parse(text=GetG0Text("WorkingRow", S=1))); 
  }
  if (!exists("ElseTauEndList")) {ElseTauEndList = NULL; }
  if (!exists("ElseGetX")) { ElseGetX = NULL; }
  if (!exists("pGroups") || pGroups <= 0) { 
    eval(parse(text=GetG0Text("pGroups", S=1)));
    if (pGroups <= 0) { pGroups <- 10; }
  }
  if (!exists("GroupSize") || GroupSize <= 0) { 
    eval(parse(text=GetG0Text("GroupSize", S=1)));
    if (GroupSize <= 0) { GroupSize <- 5; }
  }
  if (!exists("kGroups") || kGroups <= 0) { 
    eval(parse(text=GetG0Text("kGroups", S=1)));
    if (pGroups <= 0) { pGroups <- 2; }
  }
  if (!exists("tNoiseDF") || tNoiseDF <= 0) {
    eval(parse(text=GetG0Text("tNoiseDF", S=1))); 
  }
  if (!exists("GroupsSumToZero") || 
    is.null(GroupsSumToZero) || GroupsSumToZero == 0) {
    eval(parse(text=GetG0Text("GroupsSumToZero", S=1)));
    if (is.null(GroupsSumToZero) || GroupsSumToZero == 0) {
      GroupsSumToZero <- TRUE;     
    } 
  }
  
  
  BADUniqueProcessIdentifier <- 1;
  ATryGo <- "
  if (!exists(\"UniqueProcessIdentifier\")) {
    AFilePrint(\"UhOh:HitASimulation: UniqueProcessIdentifier does not exist!!\");
    UniqueProcessIdentifier <- 1;
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
    eval(parse(text=GetG0Text("UniqueProcessIdentifier", "globalenv()", S=1)));
  }
  UniqueProcessThatCreatedMeIdentifier =  UniqueProcessIdentifier;
  if (pGroups <= 0) {
    print(paste("TwoSimR5: SimGroupData supplied pGroups = ", pGroups, ", nope!  Go supply real positive number!", sep=""));
    flush.console();
    return(-999);
  }
  if (GroupSize <= 0) {
    print(paste("TwoSimR5: SimGroupData: Nope, GroupSize = ", GroupSize, " nope!, Go supply real positive number!", sep=""));
    flush.console();
    return(-999);
  }
  if (!is.null(ElseTauEndList) && !ElseGetX && length(ElseTauEndList) >= 1) {
    if (max(ElseTauEndList) == NCOL(ElseGetX)) {
      p = NCOL(ElseGetX);
      pGroups = length(ElseTauEndList);  GroupSize = max(ElseTauEndList[2:length(ElseTauEndList)] - ElseTauEndList[1:(length(ElseTauEndList)-1)]);
    } else {
      ElseTauEndList = NULL;  ElseGetX = NULL;
      p = pGroups * GroupSize;
      k = kGroups * GroupSize;      
    }
  } else {
    ElseTauEndList = NULL;  ElseGetX = NULL;
    p = pGroups * GroupSize;
    k = kGroups * GroupSize;
  }
 UniqueSimulationIdentifier = DeclareUniqueGroupSimulationIdentifier(
   ExperimentName = ExperimentName, 
   jobii = jobii, WorkingRow = WorkingRow,n=n,pGroups = pGroups, 
   GroupSize= GroupSize, kGroups= kGroups,sigma=sigma, as.integer(LogitRegression));
 
 if (SigmaNoise == -999) {
   if (sigma > 0) { SigmaNoise = sigma } else { SigmaNoise <- .5; }
 }

  BetasRealVD <- GenerateBetaGroupVec(pGroups=pGroups, kGroups=kGroups, GroupSize=GroupSize, 
    GroupsSumToZero=GroupsSumToZero, InOrder=InOrder, ElseTauEndList=ElseTauEndList);
  BetasReal <- BetasRealVD[[1]];
  BBsReal <- BetasReal; BBsReal[abs(BetasReal) > 0] <- 1;
  
  FirstGroupIndex = BetasRealVD$FirstGroupIndex;
  EndGroupIndices = BetasRealVD$EndGroupIndices;
  indexIndices <- rep(0,p); OnI =  FirstGroupIndex;
  TrueOnGroups <- rep(0, pGroups);
  
  OnI <- 1;
  for (ii in 1:length(EndGroupIndices)) { 
    indexIndices[OnI:EndGroupIndices[ii]] <- ii;

    if (max(BetasReal[OnI:EndGroupIndices[ii]]) > 0) {
      TrueOnGroups[ii] <- 1;
    }
    OnI <- EndGroupIndices[ii]+1;
  }

  DGG <- CorrelationXmatrix(4);
 
  kLen <- length(BetasReal);  NLen = n; 
  NameX <- NULL; 
  if (!is.null(ElseGetX) && NROW(ElseGetX) == n && NCOL(ElseGetX) == kLen) {
    XX <- ElseGetX;
    NameX <- colnames(ElseGetX);
    tauEndList <- ElseTauEndList;
  } else if (DGG$type== "CorrelationXmatrix1" && p > 250)  {
    XX <- matrix( rnorm(n*p,0,1), p, n );
    Rem <- sqrt(1- DGG$PosCorr^2);
    for (ii in 2:p) {
      XX[ii,] <- XX[ii-1,] * DGG$PosCorr + XX[ii,] * Rem
      ##XX[1:p,ii] = SimCor(p, DGG$PosCorr);
    }
  } else {
    IIMax <- 30;  iiSim = 1;
    while( iiSim <= IIMax ) {
      CorrRT <- CorrelationXmatrix(kLen);
      vD <- eigen(CorrRT[[1]])$values;
      if (min(vD) > 0) {
        break;
      }
      iiSim = iiSim+ 1;
      if (iiSim >= IIMax) {
        print(paste("SimMeData, eroor could not get positive definite",
          " Cor matrix in ", IIMax, " tries. ", sep="") );
        return(-1);
      }
    }
    CorrMat = CorrRT[[1]];
    LCovX <- t(chol(CorrMat));
    XX <- LCovX %*% matrix( rnorm( n * p, 0, 1), p, n );
  }
  if (MeanCenterXX == TRUE) {
    XX <- XX - rowMeans(XX);
  }
  XX <- t(XX);
  if (LogitRegression == TRUE && LogitFactor > 1) {
    BetasReal <- BetasReal * LogitFactor;
  } 
  SigmaSqNoise = SigmaNoise^2;
  if (tNoiseDF <= 0) {
    YY <-    XX %*% BetasReal + SigmaNoise * rnorm(NLen, 0,1 );
  } else {
    YY <-   XX %*% BetasReal + SigmaNoise * rnorm(NLen, 0, 1) / 
      sqrt( rchisq(NLen, tNoiseDF) / tNoiseDF );
  }
  if (MeanCenterXX == TRUE) {
    YY <- YY - mean(YY);
  }
   
  if (!exists("LogitRegression")) {
     LogitRegression <- FALSE;
  }
  if (dfRobit >= 0.0) {
    YY = XX %*% BetasReal + Beta0;
    if (dfRobit == 0.0) {
      pYY = pnorm(YY / SigmaNoise)
    } else {
      pYY = pt(YY / SigmaNoise, dfRobit)
    }
    YY = rbinom( length(pYY), 1, pYY); 
  }  else if (LogitRegression == TRUE) {
    ##Beta0 * rnorm(1,0,1) * SigmaNoise;
    YY = XX %*% BetasReal + Beta0;
    pYY = exp(YY) / ( 1 + exp(YY)) ;
    pYY[YY >= 10] = 1 - exp(-YY[YY >= 10]);
    pYY[YY < -10] = exp(YY[YY < -10]);
    YY = rbinom( length(pYY), 1, pYY); 
    ##XX = cbind(rep(1, length(YY)), XX);
  }
  eval(parse(text=GetG0Text("piSV")));
  kActiveNoise = round(length(BetasReal[BetasReal!=0]) + sqrt(NCOL(XX)) * piSV * rnorm(1, 0, 1))
  if (kActiveNoise <= 0) { kActiveNoise <- 1; }
  if (kActiveNoise >= length(BetasReal)) { kActiveNoise <- length(BetasReal)-1; }
  if (kActiveNoise >= NROW(XX)) { kActiveNoise <- NROW(XX)-1; }

  kGroupsNoise = round(kGroups + sqrt(pGroups) * piSV * rnorm(1, 0, 1))
  if (kGroupsNoise <= 0) { kGroupsNoise <- 1; }
  if (kGroupsNoise >= pGroups) { kGroupsNoise <- pGroups-1; }
  if (kGroupsNoise >= floor(NROW(XX) /GroupSize)  ) { kGroupsNoise <- floor(NROW(XX) / GroupSize)-1 ; }

   puseGroups <- kGroups / pGroups;
   puse <- length(BetasReal[BetasReal != 0.0]) / length(BetasReal);
	 if (puse > 1) { puse = (length(BetasReal)-1)/length(BetasReal);}
	 if (puse < 0) {puse = 1/length(BetasReal);}
     if (puseGroups > 1) { puseGroups <- (pGroups-1)/pGroups; }
     if (puseGroups < 0) { puseGroups <- 1 / pGroups; }
   puseNoise <- kActiveNoise/NCOL(XX);
   puseGroupsNoise <- kGroupsNoise / pGroups;
    
	 
	 ##eval(parse(text=GetG0Text("piSV", S=1)));
   eval(parse(text=GetG0Text("piSV", "globalenv()", S=1)));
   eval(parse(text=GetG0Text("sigSV", "globalenv()", S=1)));
   if (sigSV > 0) { 
     SigmaNoisesim <- sqrt(SigmaNoise^2 * rchisq(1, sigSV)/sigSV)
   } else {
      SigmaNoisesim <- (SigmaNoise)
   }

    eval(parse(text=GetGText("PriorStrengthFunction")));
    eval(parse(text=GetGText("PriorStrength")));
    eval(parse(text=GetGText("DefaultPriorStrengthFunction")));
    if (is.null(PriorStrength)) {
      if (!is.null(PriorStrengthFunction)) {
         PriorStrength <- PriorsStrengthFunction(n=n,p=p, pGroups=pGroups); 
      } else {
         PriorStrength <- DefaultPriorStrengthFunction(n=n,p=p, pGroups=pGroups)
      }
    } 
    if (!is.numeric(PriorStrength) || length(PriorStrength) != 1) {
      PriorStrength = 0;  
    }
    eval(parse(text=GetGText("SigmaPriorStrengthFunction")));
    eval(parse(text=GetGText("SigmaPriorStrength")));
    eval(parse(text=GetGText("DefaultSigmaPriorStrengthFunction")));
    if (is.null(SigmaPriorStrength)) {
      if (!is.null(SigmaPriorStrengthFunction)) {
         SigmaPriorStrength <- SigmaPriorsStrengthFunction(n=n,p=p, pGroups=pGroups); 
      } else {
         SigmaPriorStrength <- DefaultSigmaPriorStrengthFunction(n=n,p=p, pGroups=pGroups)
      }
    } 
    if (!is.numeric(SigmaPriorStrength) || length(SigmaPriorStrength) != 1) {
      SigmaPriorStrength = 0;  
    }

     
   SMSNoise <-  list(XX=XX, YY=YY, puse = puseNoise,
      puseGroups = puseGroupsNoise, ElseGetX = ElseGetX, 
          BetasReal = BetasReal, BBsReal = BBsReal,
          NLen = NLen, kLen = kLen, kActiveLen = kActiveNoise,
          p = NCOL(XX), n=length(YY), kGroups=kGroupsNoise,
          SigmaSqNoise = SigmaNoisesim^2,  SigmaNoise = SigmaNoisesim,
          SigmaNoiseReal = SigmaNoise, kActiveRealLen = kActiveLen,
          BetasName = GetBetaName(GenerateBetaVec),
          CorrName = GetCorrName(CorrelationXmatrix),
          LogitRegression=LogitRegression,  dfRobit=dfRobit,
          tNoiseDF = tNoiseDF,  SigmaPriorStrength = SigmaPriorStrength,
          UniqueSimulationIdentifier=UniqueSimulationIdentifier,
          UniqueProcessThatCreatedMeIdentifier = UniqueProcessThatCreatedMeIdentifier,
          FirstGroupIndex = FirstGroupIndex, EndGroupIndices=EndGroupIndices,
          indexIndices=indexIndices, TrueOnGroups=TrueOnGroups, piSV=piSV, sigSV=sigSV,
          PriorStrength = PriorStrength);
   MyRet = list(XX=XX, YY=YY, puse = puse,  puseGroups=puseGroups,
          BetasReal = BetasReal, BBsReal = BBsReal,
          NLen = NLen, kLen = kLen, kActiveLen = kActiveLen,
          p = NCOL(XX), n = length(YY), kGroups=kGroups,
          pGroups=pGroups,  NameX=NameX,
          SigmaSqNoise = SigmaNoise^2,  SigmaNoise = SigmaNoise,
          SigmaNoiseReal = SigmaNoise, kActiveRealLen = kActiveLen,
          BetasName = GetBetaName(GenerateBetaVec),
          CorrName = GetCorrName(CorrelationXmatrix),
          LogitRegression=LogitRegression,  dfRobit=dfRobit,
          tNoiseDF = tNoiseDF, SigmaPriorStrength = SigmaPriorStrength,
          UniqueSimulationIdentifier=UniqueSimulationIdentifier,
          UniqueProcessThatCreatedMeIdentifier = UniqueProcessThatCreatedMeIdentifier,
          FirstGroupIndex = FirstGroupIndex, EndGroupIndices=EndGroupIndices,
          indexIndices=indexIndices, TrueOnGroups=TrueOnGroups, SMSNoise = SMSNoise, 
          piSV=piSV,sigSV=sigSV, PriorStrength = PriorStrength, Beta0=Beta0);
   
   if (!exists("DontSave") || 
     is.null(DontSave) || !is.logical(DontSave) || DontSave == FALSE) {
     try(SaveCurrentNewSimulation(MyRet, UniqueSimulationIdentifier=UniqueSimulationIdentifier, 
       AlreadyLocked=AlreadyLocked, UniqueProcessIdentifier = UniqueProcessIdentifier, 
       ISample = ISample));
   }
   return(MyRet);
}


###############################################################################
##  UpdateSimValues() :
##
##    Uldates Simulation parameters from the OneVV2
##
##  n,p,k,sigma are set as well as piSV/sigSV parameters 
##    which determine innacuracy of "k supplied" to algorithms
##
UpdateSimValues <- function() {
  eval(parse(text=GetG0Text("jobii", S=1)));
  eval(parse(text=GetG0Text("OneVV2", S=1)));
  eval(parse(text=GetG0Text("OnSimType", "globalenv()", S=1)));
  if (!is.null(OnSimType) && OnSimType == "Group") {
     UpdateGroupSimValues();
  }
  if (jobii == 0) { jobii = 1; }
  AFilePrint(paste("Onjobii = ", jobii, sep=""));
  AFilePrint(paste("   OneVV2 = ", paste(OneVV2[jobii,], collapse=", "), sep=""));
  ##print(paste("dimOnevv2 = c(", paste(dim(OneVV2), collapse=", "),
  ##  sep=""));
  ##print(OneVV2);
  eval(parse(text=GetG0Text("OneVV2", S=1)));
  eval(parse(text=GetG0Text("n", S=1)));
  eval(parse(text=GetG0Text("p", S=1)));
  eval(parse(text=GetG0Text("k", S=1)));
  eval(parse(text=GetG0Text("NN", S=1)));
  eval(parse(text=GetG0Text("NR", S=1)));
  eval(parse(text=GetG0Text("Gimme", S=1)));  
  eval(parse(text=GetG0Text("sigma", S=1)));
  n = OneVV2[jobii,1];   NN = n;
  p = OneVV2[jobii,2];   NR = p;
  k = OneVV2[jobii,3];   Gimme = k;
  sigma = OneVV2[jobii,4];
  eval(parse(text=SetGText("OneVV2",S=1)));
  eval(parse(text=SetGText("n",S=1)));
  eval(parse(text=SetGText("p",S=1)));
  eval(parse(text=SetGText("k",S=1)));
  eval(parse(text=SetGText("NN",S=1)));  
  eval(parse(text=SetGText("NR",S=1)));  
  eval(parse(text=SetGText("Gimme",S=1)));
  eval(parse(text=SetGText("sigma",S=1)));  
  if (dim(OneVV2)[2] >= 6) {
    ##piSV = OneVV2[jobii,5];  sigSV = OneVV2[jobii,6];
    eval(parse(text=GetG0Text("piSV", S=1)));  eval(parse(text=GetG0Text("sigSV", S=1)));
    piSV <- OneVV2[jobii,5];  sigSV <- OneVV2[jobii,6];
    eval(parse(text=LockGText("piSV",S=1))); eval(parse(text=LockGText("sigSV",S=1))); 
  }
}


UpdateGroupSimValues <- function() {
  eval(parse(text=GetG0Text("jobii", S=1)));
  eval(parse(text=GetG0Text("OneVV2", S=1)));
  eval(parse(text=GetG0Text("OnSimType", "globalenv()", S=1)));
  if (jobii == 0) { jobii = 1; }
  AFilePrint(paste("Onjobii = ", jobii, sep=""));
  AFilePrint(paste("   OneVV2 = ", paste(OneVV2[jobii,], collapse=", "), sep=""));
  ##print(paste("dimOnevv2 = c(", paste(dim(OneVV2), collapse=", "),
  ##  sep=""));
  ##print(OneVV2);
  eval(parse(text=GetG0Text("OneVV2", S=1)));
  eval(parse(text=GetG0Text("n", S=1)));
  eval(parse(text=GetG0Text("p", S=1)));
  eval(parse(text=GetG0Text("k", S=1)));
  eval(parse(text=GetG0Text("pGroups", S=1)));
  eval(parse(text=GetG0Text("kGroups", S=1)));
  eval(parse(text=GetG0Text("GroupSize", S=1)));
  eval(parse(text=GetG0Text("NN", S=1)));
  eval(parse(text=GetG0Text("NR", S=1)));
  eval(parse(text=GetG0Text("Gimme", S=1)));  
  eval(parse(text=GetG0Text("sigma", S=1)));
  eval(parse(text=GetG0Text("TDfNoise", S=1)));
  eval(parse(text=GetG0Text("sigSV", S=1)));
  eval(parse(text=GetG0Text("piSV", S=1)));
  n = OneVV2[jobii,1];   NN = n;
  pGroups = OneVV2[jobii,2];
  GroupSize = OneVV2[jobii,3];
  kGroups = OneVV2[jobii,4];
  p = pGroups*GroupSize;   NR = p;
  k = kGroups*GroupSize;   Gimme = k;
  sigma = OneVV2[jobii,4];
  if (NCOL(OneVV2) == 7) {
    try(TDfNoise <- OneVV2[jobii,7]);
  }
  if (NCOL(OneVV2) == 5) {
    try(piSV <- OneVV2[jobii, 5]);
  }
  if (NCOL(OneVV2) == 6) {
    try(sigSV <- OneVV2[jobii,6]);
  }
  eval(parse(text=SetGText("OneVV2",S=1)));
  eval(parse(text=SetGText("n",S=1)));
  eval(parse(text=SetGText("p",S=1)));
  eval(parse(text=SetGText("k",S=1)));
  eval(parse(text=SetGText("pGroups",S=1)));
  eval(parse(text=SetGText("kGroups",S=1)));
  eval(parse(text=SetGText("GroupSize", S=1)));
  eval(parse(text=SetGText("NN",S=1)));  
  eval(parse(text=SetGText("NR",S=1)));  
  eval(parse(text=SetGText("Gimme",S=1)));
  eval(parse(text=SetGText("sigma",S=1))); 
  eval(parse(text=SetGText("piSV", S=1)));
  eval(parse(text=SetGText("sigSV", S=1)));
  eval(parse(text=SetGText("TDfNoise", S=1))); 
  if (dim(OneVV2)[2] >= 6) {
    ##piSV = OneVV2[jobii,5];  sigSV = OneVV2[jobii,6];
    eval(parse(text=GetG0Text("piSV", S=1)));  eval(parse(text=GetG0Text("sigSV", S=1)));
    piSV <- OneVV2[jobii,5];  sigSV <- OneVV2[jobii,6];
    eval(parse(text=LockGText("piSV",S=1))); eval(parse(text=LockGText("sigSV",S=1))); 
  }
}