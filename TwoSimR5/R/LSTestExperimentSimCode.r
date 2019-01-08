TwoSimLoadIn <- function() {
  ##MyColNames = c("PiD", "l1error", "l2error", "Type1", "Type2", "L.zero", "time1", "time2", "time3") 
  MyColNames = c("l1error", "l2error", "Type1", "Type2", "L.zero", "time1", "time2", "time3") 
  MyPermaColNames <- MyColNames;
  eval(parse(text=SetGText("MyColNames", "globalenv()", S=1)));
  eval(parse(text=SetGText("MyPermaColNames", "globalenv()", S=1)));
}
SimMeData22 <- function(n, p, k, sigma) {	 
	# Utility functions to normalize with L_p norm
	norm <- function(v, p=2) { sum(abs(v)^p)^(1/p) }
	normalize <- function (v, p=2) { v / norm(v, p) }

	beta <- rep(0, p)
	nonzeros <- sample(1:p, size=k)
	if (simBetas == 0)  {
	   beta[nonzeros] <- runif(length(nonzeros), min=1, max=100)
	} else {
     beta[nonzeros]  <- simBetas * (2* rbinom(length(nonzeros), 1, .5) - 1);
  }
	BetasReal <- matrix(beta, p, 1)
	BBsReal <- BetasReal; BBsReal[abs(BetasReal) > 0] <- 1;

	puse <- k / p;
  
	# Normalized matrix of predictors
	XX <- normalize(matrix(rnorm(n), n, 1))
	for (i in 1:(p-1)) {
	C <- normalize(matrix(rnorm(n), n, 1))
        XX <- cbind(XX, C)
	}
  
  if (simSigma >= 0)  {
     SigmaSqNoise = simSigma^2;
     sigma = simSigma
  } else {
     SigmaSqNoise = sigma^2
  }

	z <- rnorm(n, mean=0, sd=sigma)
	
	YY <- XX %*% BetasReal + z;
   
	MyRet = list(XX=XX, YY=YY, puse = puse, 
         BetasReal = BetasReal, BBsReal = BBsReal,
         Nlen = n, kLen = p, kActiveLen = k,
         SigmaSqNoise = sigma^2,  SigmaNoise = sigma);

	return(MyRet);
}



# Normalizes L_p norm
lnorm <- function(x,p) {
sum(abs(x)^p)^(1/p)}


TableCompareSaverR5$methods(
   AssessMIPFit = function(SMS, FittedOb,
     UniqueProcessIdentifier,
     UniqueSimulationIdentifier, 
     NameFunction, 
     Hit = "+",
     quoteMore = " - ", 
     verbose = 0, DontRecord=FALSE, AlreadyLocked=FALSE) {
     if (is.null(FittedOb$MIPReport)) {
       AFilePrint(paste("AssessMIPFit: Hey, we tried to do MIP fit for ", NameFunction, " but got NULL!"));
     } 
     MyTestT <- " 
     if (is.logical(verbose) && verbose == FALSE) {
       verbose = FALSE;
     } else if (is.logical(verbose) && verbose == TRUE) { verbose = 2; 
     } else if (!is.numeric(verbose)) {
       printf(\"AssessMIPFit: Hey, this isn't a good verbose;\"); verbose = 2;
     } else {
       verbose = .self$verbose-2;
     } ";
     try(eval(parse(text=MyTestT)));
     if (verbose >= 1) {
       AFilePrint(paste("AssessMIPFit; Starting for ", NameFunction, ", Hit = ", Hit, sep=""));
       flush.console();
     }
    if (is.null(FittedOb$MIPReport) || length(FittedOb$MIPReport) <= 0) {
      AFilePrint(paste("AssessMIPFit: No MIPReport Given for ", NameFunction, sep="")); flush.console();
      return(-1);
    }
    Type1 = NULL; Type2=NULL; SphereLoss = NULL; AUC=NULL; MarginalHellinger=NULL;
    try(Type1 <- Type1Function(FittedOb$MIPReport, SMS$BetasReal));
    try(Type2 <- Type2Function(FittedOb$MIPReport, SMS$BetasReal));
    try(SphereLoss <- SphereLossFunction(FittedOb$MIPReport, SMS$BetasReal));
    try(AUC <- AUCFunction(FittedOb$MIPReport, SMS$BetasReal));
    try(MarginalHellinger <- MarginalHellingerFunction(FittedOb$MIPReport, SMS$BetasReal));
    
    Type1Random=NULL; Type2Random=NULL; SphereLossRandom=NULL;
    AUCRandom=NULL; MarginalHellingerRandom=NULL; 
    if (!is.null(SMS$TrueOnGroups) && length(SMS$TrueOnGroups) >= 1) {
      if (!is.null(FittedOb$MIPRandom) && length(FittedOb$MIPRandom) == length(SMS$TrueOnGroups)) {
      try(Type1Random <- Type1Function(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(Type2Random <- Type2Function(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(SphereLossRandom <- SphereLossFunction(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(AUCRandom <- AUCFunction(FittedOb$MIPRandom, SMS$TrueOnGroups));
      try(MarginalHellingerRandom <- MarginalHellingerFunction(FittedOb$MIPRandom, SMS$TrueOnGroups));
      }
    }
  AllMIPReturn <- list(Type1=Type1, Type2=Type2, SphereLoss=SphereLoss,
  AUC=AUC, MarginalHellinger=MarginalHellinger, Type1Random=Type1Random,
  Type2Random=Type2Random, SphereLossRandom=SphereLossRandom,
  AUCRandom=AUCRandom,MarginalHellingerRandom = MarginalHellingerRandom)
  if (is.logical(DontRecord) && DontRecord == TRUE) {
    return(AllMIPReturn);
  }
  eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
  if (verbose >= 2) {
      AFilePrint(paste("----  AssessMIPFit[", UniqueSimulationIdentifier, 
        "] --- CurrentSmallContainDir = ", CurrentSmallContainDir, " for NameFunction = ", NameFunction, sep=""));
      flush.console();
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
    AFilePrint("AssessMIPFit: UnqiueProcessIdentifier is NULL!"); flush.console();  
  }
  if (verbose >= 2) {
    AFilePrint(paste("  --- AssedMIPFit, ", NameFunction, " looking for HMS. ", sep="")); flush.console();
  }
  FunctionDir <- paste(CurrentSmallContainDir, "//", NameFunction, sep="");
  try(dir.create(FunctionDir, showWarnings=FALSE,recursive=TRUE));
  try(Oldwd <- getwd());
  try(setwd(CurrentSmallContainDir));
  while(LockMeIn(verbose=as.numeric(verbose), quoteMore=quoteMore, LFDir = NameFunction)==FALSE){
    Sys.sleep(runif(1,0,4));
  }
  GoalFile <- paste(NameFunction, "//", "MIP", UniqueSimulationIdentifier, ".RData", sep="");
  if (any(unlist(list.files(FunctionDir)) == 
    paste("MIP", UniqueSimulationIdentifier, ".RData", sep=""))) {
    AFilePrint(paste("Woah, AssessedMIP, Function = ", NameFunction, " UniqueSimulationIdentifier = ",
      UniqueSimulationIdentifier, ", We already have in directory!  We're saving over anyway!", sep=""));
    unlink(GoalFile);
  }
  save(AllMIPReturn=AllMIPReturn, NameFunction=NameFunction,
     Hit=Hit,  FunctionDir = FunctionDir, file=GoalFile)
  if (FALSE) {
    save(AllMIPReturn=AllMIPReturn, NameFunction=NameFunction,
     Hit=Hit,  FunctionDir = FunctionDir, file=paste(NameFunction, "//ATryMIP.RData", sep=""));
  }
  try(setwd(CurrentSmallContainDir));
  try(UnLockMe(verbose=verbose, quoteMore=paste(quoteMore, " - AssedMIP - ", NameFunction, sep=""), 
    LFDir = NameFunction));
  try(setwd(Oldwd));
  if (verbose >= 3) {
    AFilePrint(paste("  --- AssessedMIP: Finish ", NameFunction, " length = ", length(AllMIPReturn), sep="")); flush.console();
  }
  return(AllMIPReturn);
});
TableCompareSaverR5$methods(
   AssessCIFit = function(SMS, FittedOb,
     UniqueProcessIdentifier,
     UniqueSimulationIdentifier, 
     NameFunction, 
     Hit = "+",
     quoteMore = " - ", 
     verbose = 0, DontRecord=FALSE, AlreadyLocked=FALSE) {
     if (is.null(FittedOb$CIEst)) {
       AFilePrint(paste("AssessCIFit: Hey, we tried to do CI fit for ", NameFunction, " but got NULL!"));
     } 
     MyTestT <- " 
     if (is.logical(verbose) && verbose == FALSE) {
       verbose = FALSE;
     } else if (is.logical(verbose) && verbose == TRUE) { verbose = 2; 
     } else if (!is.numeric(verbose)) {
       printf(\"AssessCIFit: Hey, this isn't a good verbose;\"); verbose = 2;
     } else {
       verbose = .self$verbose-2;
     } ";
     try(eval(parse(text=MyTestT)));
     if (verbose >= 1) {
       AFilePrint(paste("AssessCIFit; Starting for ", NameFunction, ", Hit = ", Hit, sep=""));
       flush.console();
     }
     Betas <- SMS$BetasReal;
     IndexOnBetas <- (1:length(Betas))[Betas!=0.0];
     IndexOffBetas <- (1:length(Betas))[Betas==0.0];
     IndexSmallBetas <- (1:length(Betas))[Betas!=0.0 &
       Betas^2 <= SMS$SigmaNoiseReal / length(SMS$NLen)];
     CIQuantiles <- FittedOb$CIQuantiles; 
     
     if (CIQuantiles[1] == .5) {
       AK <- 1; 
     } else { AK <- 0;}
     QuantileCov <- rep(0, length(CIQuantiles)/2);
     ArAt <- AK;
     for (jj in 1:length(QuantileCov)) {
       QuantileCov[jj] <- CIQuantiles[ArAt + 2] - CIQuantiles[ArAt+1];
       ArAt <- ArAt +2;
     }
     CITime = FittedOb$CITime;
     
     AllCIReturn <- list();
     for (tjt in 1:length(FittedOb$CIEst)) {
       if (any(CIQuantiles == .5)) {
         MedianAbsDiff <-sum(abs(Betas - FittedOb$CIEst[[tjt]][,CIQuantiles == .5]))
         MedianOnAbsDiff <- sum(abs(Betas[IndexOnBetas] - FittedOb$CIEst[[tjt]][IndexOnBetas,CIQuantiles == .5]))
         MedianOffAbsDiff <- sum(abs(Betas[IndexOffBetas] - FittedOb$CIEst[[tjt]][IndexOffBetas,CIQuantiles == .5]))
         MedianSqDiff <- sum((Betas - FittedOb$CIEst[[tjt]][,CIQuantiles == .5])^2)
         MedianOnSqDiff <- sum((Betas[IndexOnBetas] - FittedOb$CIEst[[tjt]][IndexOnBetas,CIQuantiles == .5])^2)
         MedianOffSqDiff <- sum((Betas[IndexOffBetas] - FittedOb$CIEst[[tjt]][IndexOffBetas,CIQuantiles == .5])^2)
         if (length(IndexSmallBetas) >= 1) {
           MedianSmallAbsDiff <- sum(abs(Betas[IndexSmallBetas] - FittedOb$CIEst[[tjt]][IndexSmallBetas,CIQuantiles == .5]))
           MedianSmallSqDiff <- sum((Betas[IndexSmallBetas] - FittedOb$CIEst[[tjt]][IndexSmallBetas,CIQuantiles == .5])^2)            
         } else {
           MeanSmallAbsDiff <- NULL; MeanSmallSqDiff = NULL;  
           MedianSmallAbsDiff <- NULL; MedianSmallSqDiff = NULL;         
         }          
       } else {
         MedianAbsDiff = NULL;  MedianOnAbsDiff = NULL;
         MedianOffAbsDiff = NULL;  MedianSmallAbsDiff = NULL;
         MedianSqDiff = NULL;  MedianOnSqDiff = NULL;
         MedianOffSqDiff = NULL;  MedianSmallSqDiff = NULL;
       }
       ArAt <- AK;
       CoverageQuantiles <- rep(0, length(QuantileCov));
       CoverageOnQuantiles <- rep(0, length(QuantileCov));
       CoverageOffQuantiles <- rep(0, length(QuantileCov));
       CoverageSmallQuantiles <- rep(0, length(QuantileCov));              
       MeanWidthQuantiles <- rep(0, length(QuantileCov));
       MeanOnWidthQuantiles <- rep(0, length(QuantileCov));
       MeanOffWidthQuantiles <- rep(0, length(QuantileCov));
       MeanSmallWidthQuantiles <- rep(0, length(QuantileCov));
       CIConservativeSignLoss <- -1;
       CILiberalSignLoss <- -1;
                     
       MeanSqWidthQuantiles <- rep(0, length(QuantileCov));
       MeanSqOnWidthQuantiles <- rep(0, length(QuantileCov));
       MeanSqOffWidthQuantiles <- rep(0, length(QuantileCov));
       MeanSqSmallWidthQuantiles <- rep(0, length(QuantileCov)); 
       
       CINormalLoss <- rep(0, length(QuantileCov));
       CIUnifLoss <- rep(0, length(QuantileCov));
       CIOnNormalLoss <- rep(0, length(QuantileCov));
       CIOnUnifLoss <- rep(0, length(QuantileCov));
       CIOffNormalLoss <- rep(0, length(QuantileCov));
       CIOffUnifLoss <- rep(0, length(QuantileCov));
       CISmallNormalLoss <- rep(0, length(QuantileCov));
       CISmallUnifLoss <- rep(0, length(QuantileCov));
                          
       for (jj in 1:length(QuantileCov)) {
         CoverageQuantiles[jj] <- sum( FittedOb$CIEst[[tjt]][,ArAt+1] <=
            Betas & FittedOb$CIEst[[tjt]][,ArAt+2] >= Betas) /length(Betas);
         CoverageOnQuantiles[jj] <- sum( FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+1] <=
            Betas[IndexOnBetas] & FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+2] >= Betas[IndexOnBetas])/
            length(IndexOnBetas);
         CoverageOffQuantiles[jj] <- sum( FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+1] <=
            Betas[IndexOffBetas] & FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+2] >= Betas[IndexOffBetas])/
            length(IndexOffBetas);
         CINormal <- CINewLoss(FittedOb$CIEst[[tjt]][, ArAt+1], 
           FittedOb$CIEst[[tjt]][, ArAt+2], Betas, QuantileCov[jj]);
         if (length(CIConservativeSignLoss) == 1 && CIConservativeSignLoss[1] == -1) {
           CIConservativeSignLoss <- CICalcConservativeSignLoss(FittedOb$CIEst[[tjt]][, ArAt+1],
            FittedOb$CIEst[[tjt]][, ArAt+2], Betas);
         } else {
           CIConservativeSignLoss <- c(CIConservativeSignLoss, 
             CICalcConservativeSignLoss(FittedOb$CIEst[[tjt]][, ArAt+1],
              FittedOb$CIEst[[tjt]][, ArAt+2], Betas));                  
         }
         if (length(CILiberalSignLoss) == 1 && CILiberalSignLoss[1] == -1) {
           CILiberalSignLoss <- CICalcLiberalSignLoss(FittedOb$CIEst[[tjt]][, ArAt+1],
            FittedOb$CIEst[[tjt]][, ArAt+2], Betas);
         } else {
           CILiberalSignLoss <- c(CILiberalSignLoss, 
             CICalcLiberalSignLoss(FittedOb$CIEst[[tjt]][, ArAt+1],
              FittedOb$CIEst[[tjt]][, ArAt+2], Betas));                  
         }
         CINormalLoss[jj] <- sum(CINormal);
         CIOnNormalLoss[jj] <- sum(CINormal[IndexOnBetas]);
         CIOffNormalLoss[jj] <- sum(CINormal[IndexOffBetas]);
         try(CIUnif <- CIUnifLoss(FittedOb$CIEst[[tjt]][, ArAt+1], 
           FittedOb$CIEst[[tjt]][, ArAt+2], Betas, QuantileCov[jj]));
         CIUnifLoss[jj] <- sum(CIUnif);
         CIOnUnifLoss[jj] <- sum(CIUnif[IndexOnBetas]);
         CIOffUnifLoss[jj] <- sum(CIUnif[IndexOffBetas]);
        
         MeanWidthQuantiles[jj] <- sum(mean(FittedOb$CIEst[[tjt]][,ArAt+2]-
           FittedOb$CIEst[[tjt]][,ArAt+1]));
         MeanOnWidthQuantiles[jj] <- sum(mean(FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+1]));           
         MeanOffWidthQuantiles[jj] <- sum(mean(FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+1]));
         MeanSqWidthQuantiles[jj] <- sum(mean((FittedOb$CIEst[[tjt]][,ArAt+2]-
           FittedOb$CIEst[[tjt]][,ArAt+1])^2));
         MeanSqOnWidthQuantiles[jj] <- sum(mean((FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexOnBetas,ArAt+1])^2));           
         MeanSqOffWidthQuantiles[jj] <- sum(mean((FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexOffBetas,ArAt+1])^2));           


         if (length(IndexSmallBetas) >= 1) {
           CoverageSmallQuantiles[jj] <- sum( FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+1] <=
            Betas[IndexSmallBetas] & FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+2] >= Betas[IndexSmallBetas])/
            length(IndexSmallBetas);
           MeanSmallWidthQuantiles[jj] <- sum(mean(FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+1]));                 
           MeanSqSmallWidthQuantiles[jj] <- sum(mean((FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+2]-
           FittedOb$CIEst[[tjt]][IndexSmallBetas,ArAt+1])^2)); 
           try(CISmallNormalLoss[jj] <- sum(CINormal[IndexSmallBetas]));
           try(CISmallUnifLoss[jj] <- sum(CIUnif[IndexSmallBetas]));

         } else {
           CoverageSmallQuantiles <- NULL; MeanSqSmallWidthQuantiles = NULL;
           MeanSmallWidthQuantiles = NULL;    
           CISmallNormalLoss <- NULL;  CISmallUnifLoss <- NULL;      
         } 
         ArAt <- ArAt+2;           
       }
       APiece <- list(
        NameOfInterval = names(FittedOb$CIEst)[tjt],
        NameFunction = NameFunction, Hit=Hit,
        CIQuantiles = FittedOb$CIQuantiles,
        MedianAbsDiff = MedianAbsDiff,
        MedianOnAbsDiff = MedianOnAbsDiff,
        MedianOffAbsDiff = MedianOffAbsDiff, MedianSqDiff = MedianSqDiff,
        MedianOnSqDiff = MedianOnSqDiff, MedianOffSqDiff = MedianOffSqDiff,
        MedianSmallAbsDiff = MedianSmallAbsDiff,  MedianSmallSqDiff = MedianSmallSqDiff,   
         CoverageQuantiles = CoverageQuantiles,
         MeanWidthQuantiles = MeanWidthQuantiles,
         MeanSqWidthQuantiles = MeanSqWidthQuantiles,
         CoverageOnQuantiles =CoverageOnQuantiles,
         MeanOnWidthQuantiles = MeanOnWidthQuantiles,
         MeanSqOnWidthQuantiles = MeanSqOnWidthQuantiles,
         CoverageOffQuantiles = CoverageOffQuantiles,
         MeanOffWidthQuantiles = MeanOffWidthQuantiles,
         MeanSqOffWidthQuantiles = MeanSqOffWidthQuantiles,
         CoverageSmallQuantiles = CoverageSmallQuantiles,
         MeanSmallWidthQuantiles = MeanSmallWidthQuantiles,
         MeanSqSmallWidthQuantiles = MeanSqSmallWidthQuantiles,
         CINormalLoss =  CINormalLoss,
         CIUnifLoss = CIUnifLoss,
         CIOnNormalLoss =  CIOnNormalLoss,
         CIOnUnifLoss = CIOnUnifLoss,
         CIOffNormalLoss = CIOffNormalLoss, 
         CIOffUnifLoss = CIOffUnifLoss,
         CISmallNormalLoss = CISmallNormalLoss,
         CIConservativeSignLoss = CIConservativeSignLoss,
         CILiberalSignLoss = CILiberalSignLoss,
         CISmallUnifLoss = CISmallUnifLoss);
       AllCIReturn[[tjt]] <- APiece;
   }
   
  if (is.logical(DontRecord) && DontRecord == TRUE) {
    return(AllCIReturn);
  }
  eval(parse(text=GetG0Text("CurrentSmallContainDir", S=1)));
  if (verbose >= 2) {
      AFilePrint(paste("----  AssessCIFit[", UniqueSimulationIdentifier, 
        "] --- CurrentSmallContainDir = ", CurrentSmallContainDir, " for NameFunction = ", NameFunction, sep=""));
      flush.console();
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
    AFilePrint("AssessCIFit: UnqiueProcessIdentifier is NULL!"); flush.console();  
  }
  if (verbose >= 2) {
    AFilePrint(paste("  --- AssedCIFit, ", NameFunction, " looking for HMS. ", sep="")); flush.console();
  }
  FunctionDir <- paste(CurrentSmallContainDir, "//", NameFunction, sep="");
  try(dir.create(FunctionDir, showWarnings=FALSE,recursive=TRUE));
  try(Oldwd <- getwd());
  try(setwd(CurrentSmallContainDir));
  while(LockMeIn(verbose=as.numeric(verbose), quoteMore=quoteMore, LFDir = NameFunction)==FALSE){
    try(Sys.sleep(runif(1,0,4)));
  }
  GoalFile <- paste(NameFunction, "//", "CIT", UniqueSimulationIdentifier, ".RData", sep="");
  if (any(unlist(list.files(NameFunction)) == 
    paste("CIT", UniqueSimulationIdentifier, ".RData", sep=""))) {
    AFilePrint(paste("Woah, AssessedCI, Function = ", NameFunction, " UniqueSimulationIdentifier = ",
      UniqueSimulationIdentifier, ", We already have in directory!", sep=""));
    unlink(GoalFile);
  }
  CIEst = NULL;
  if (!exists("FittedOb") || is.null(FittedOb)) {
    AFilePrint("Why did FittedOb AssesCI? that not happen what?"); flush.console(); 
  } else if (is.null(FittedOb$CIEst)) {
    AFilePrint("What: Why is FittedOb$CIEst null?"); flush.console();
  } else {
    CIEst <- FittedOb$CIEst;
  }
  save(AllCIReturn=AllCIReturn, NameFunction=NameFunction, Betas=Betas, CIQuantiles=CIQuantiles,
    CIEst = CIEst, Hit=Hit, CITime=CITime, file=GoalFile)
  UnLockMe(verbose=verbose, quoteMore=paste(quoteMore, " - AssedCI - ", NameFunction, sep=""), LFDir = NameFunction);
  try(setwd(Oldwd));
  if (verbose >= 3) {
    AFilePrint(paste("  --- AssessedCI: Finish ", NameFunction, " length = ", length(AllCIReturn), sep="")); flush.console();
  }
  return(AllCIReturn);
  }
);

TableCompareSaverR5$methods(
  RunOneSimExperiment =  function(verbose = -1, DoCI = DefaultDoCI, ForceFunction=NULL, 
    ForceAFail=FALSE,...) {
  try(setwd(.self$OriginalOldwd));
  try(eval(parse(text=AOverallPaste(".self$ListOfParams"))), silent=TRUE);
  #Use custom version of Alan's SimMeData
  if (!exists("ForceFunction")) { ForceFunction <- NULL; }
  if (exists("verbose") && !is.null(verbose) && is.numeric(verbose) &&
	  length(verbose) == 1 && verbose == -1) {
	  MyT <- "verbose <- .self$verbose;"
	  try(eval(parse(text=MyT)), silent=TRUE);
  } else if (!exists("verbose") || is.null(verbose)  || !is.numeric(verbose)) {
    eval(parse(text=GetG0Text("verbose")));
    if (is.null(verbose) || !is.numeric(verbose) || verbose[1] == 0) {
      MyT <- "verbose <- .self$verbose";
      try(eval(parse(text=MyT)), silent=TRUE);
    } else if (is.logical(verbose)) {
      if (verbose == TRUE) {
        if (.self$verbose >= 1) {
        } else {
          try(.self$verbose <- as.integer(1));
        }
      }  else {
        try(.self$verbose <- as.integer(0));
      }
    } else if (is.numeric(verbose) && verbose[1] == -1) {
      MyT <- "verbose <- .self$verbose";
      try(eval(parse(text=MyT)), silent=TRUE);
    } else {
      try(.self$verbose <- verbose);
    }
  } else {
    MyT <- "verbose <- .self$verbose";
    try(eval(parse(text=MyT)), silent=TRUE);
  }


  ##SS <-LookForNewSimulation(verbose = .self$verbose);
  SS <- LookToGetNewSimulation(verbose = .self$verbose, 
    quoteMore = "LookToGetNewSimulation : ", ForceFunction = ForceFunction,
    TCS = .self) 
  if (is.null(SS)) { 
    AFilePrint("LookToGetNewSimulation Returns Null, must have gone poorly!");
    flush.console(); 
    try(setwd(.self$OriginalOldwd));
    return(-1); 
  } else if (is.numeric(SS) && SS[1] == 111) {
     AFilePrint("LookToGetNewSimulation Returns 111, there was a full success!")
     try(setwd(.self$OriginalOldwd));
     return(1);
  } else if (is.numeric(SS) && SS[1] <= -1) {
    AFilePrint("LookForNewSimulationReturns a -1, bad sign!"); flush.console();
    try(setwd(.self$OriginalOldwd));
    return(-1);
  } else if (verbose == TRUE || (is.numeric(verbose) && verbose > 0)) {
    AFilePrint("RunOneSimExperiment: We got a new Simulation"); flush.console(); 
  }
  
  SMS = SS$SMS;
  MyNameFunction = SS$MyNameFunction;
  IFunction = SS$OnNumFunction;
  eval(parse(text=SetGText("IFunction", S=1)));
  eval(parse(text=SetGText("MyNameFunction", S=1)));
  eval(parse(text=SetGText("SMS", S=1)));  
  if (is.null(SMS$UniqueSimulationIdentifier)) {
    AFilePrint("RunExperimentRun:  Error: SMS$UniqueSimulationIdentifier is NULL! "); flush.console();
    try(setwd(.self$OriginalOldwd));
    return(-666);
  }
  if (verbose > 0) {
    AFilePrint(paste("We're going to run function: ", .self$OnNameFunctions[IFunction], sep=""));
    flush.console();
  }
  
  RTI = unlist(.self$OnUseFunctions[IFunction]);
  eval(parse(text=SetGText("SMS", S=1)));
    ##AFilePrint(paste(" Doing Analysis NameFunctions[", ii, "] = ",
    ##   NameFunctions[ii], sep="")); 
  if (p > 1000) {
    AFilePrint(paste("Doing Analysis NameFunctions[", IFunction, "] = ",
        .self$OnNameFunctions[IFunction], sep=""));  flush.console();
  }
  limlas = NULL;
  try(limlas <- (RTI[[1]])(SMS));
  eval(parse(text=SetGText("limlas", S=1)));
  if (verbose > 0) {
    AFilePrint(paste("We just ran function: ", .self$OnNameFunctions[IFunction], sep=""));
    flush.console();
  }
  if (ForceAFail == TRUE) {
    AFilePrint(paste("We ran function ", .self$OnNameFunctions[IFunction], 
      " but deliberately call it a failure. "));
    limlas <- NULL;
  }
  if (is.null(limlas)) {
    eval(parse(text=GetG0Text("NameFunctions")));
    AFilePrint(paste("We've got a limlas error, we are on ii = ", IFunction, " and NameFunctions = ",
      .self$OnNameFunctions[IFunction], sep="")); flush.console();
    beta.hat = SMS$BetasReal * 0 - 666;
    l1error = -2;
    l2error = -2;
    Type1 = -2; Type2 = -2; 
    L.zero <- -2 
    time = c(NA,NA,NA);
    Hit = "D"
    FittedOb = NULL;
  } else if (!is.null(limlas)) {
    FittedOb = limlas;
    beta = SMS$BetasReal
  	beta.hat = limlas$BetaFit
  	if (IFunction == 1) {
      BetaStart = beta.hat;
      eval(parse(text=SetGText("BetaStart", S=1)));
    }
  	
  	# l1 and l2 errors
  	l1error <- lnorm(beta.hat - beta, 1) / lnorm(beta, 1)
  	l2error <- lnorm(beta.hat - beta, 2) / lnorm(beta, 2)
    if (beta.hat[1] == -999) {
      l2error = -999;
      l1error = -999;
    }
    # Type1 and Type2 errors
    Type1 <- sum(((SMS$BBsReal - limlas$BBHatsFit)==-1)*1)
    Type2 <- sum(((-1*SMS$BBsReal + limlas$BBHatsFit)==-1)*1)
  	
    # L_0 Error
    L.zero <- 0.5*Type1 + 0.5*Type2  	
    if (SMS$BBsReal[1] <= -1) {
      Type1 = -2;
      Type2 = -2;
      L.zero = -2;
    }
  

  	if (length(limlas$FitTime) >= 3) {
     	time <- as.numeric(limlas$FitTime[1:3])
    } else {
      time <- c(NA, NA, NA);
    }
    eval(parse(text=GetG0Text("CurrentSmallContainDir")));
    if (p > 5000) {
      paste("    Onjobii = ", jobii, 
        " Finished ", IFunction, " = ", .self$OnNameFunctions[IFunction], 
        " in time ", sum(time[1:2]), sep="");
      dir.create(CurrentSmallContainDir, showWarnings = FALSE)
      FilesInDir = unlist(list.files(CurrentSmallContainDir));
      if (any(FilesInDir == "Progress.txt")) {
        fp <- file(paste(CurrentSmallContainDir, "//Progress.txt", sep=""), "a")
      }   else {
        fp <- file(paste(CurrentSmallContainDir, "//Progress.txt", sep=""), "w")
      }
      writeLines(text=paste("Onjobii = ", jobii, 
                         " Finished ", IFunction, " = ", .self$OnNameFunctions[IFunction], 
                         " in time ", sum(time[1:2]), sep=""),
                   con = fp);
      close(fp);
    }
   
   Hit = "+"
   }
    if (!is.numeric(l1error) || is.na(l1error) || !is.finite(l1error) || length(l1error) != 1) { l1error = NA; }
    if (!is.numeric(l2error) || is.na(l2error) || !is.finite(l2error) || length(l2error) != 1) { l2error = NA; }
    if (!is.numeric(Type1) || is.na(Type1) || !is.finite(Type1) || length(Type1) != 1) { Type1 = NA; }
    if (!is.numeric(Type2) || is.na(Type2) || !is.finite(Type2) || length(Type2) != 1) { Type2 = NA; }
    if (!is.numeric(L.zero) || is.na(L.zero) || !is.finite(L.zero) || length(L.zero) != 1) { L.zero = NA; }        
   ## SimVec = c(l1error, l2error, Type1, Type2, L.zero, as.numeric(time[1:3]));
   MyV = c(l1error, l2error, Type1, Type2, L.zero, time);
   if (!is.numeric(IFunction) || IFunction > length(.self$OnNameFunction) ||
     IFunction <= 0) {
     AFilePrint(paste("Oh No, IFunction = ", IFunction, " this won't work.!", sep=""));  
     flush.console();
   }
   TryToAddMyVToRDataTable <- "
     AT <- FALSE; 
     AddMyVToRDataTable(OnFunction = .self$OnNameFunctions[IFunction],
       RDataTableOutName = .self$RDataTableOutName, MyV=MyV, 
       PiD=SMS$UniqueSimulationIdentifier, verbose=verbose, AlreadyLocked=FALSE);
     AT <- TRUE;
   ";
   try(eval(parse(text=TryToAddMyVToRDataTable)));
   if (AT == FALSE) {
     AFilePrint("LSTExperientSimCode: We tried to AddMyVToRDataTable but failed.");
   }

   AHit = NULL;
   if (IFunction > length(.self$OnNameFunctions)) {
     AFilePrint(paste("LSTExperimentSimCode: Whats this, IFunction = ", IFunction,
       " but length OnNameFunctions is ", length(.self$OnNameFunctions),
       sep="")); flush.console();
   }
   if (verbose >= 1) {
     AFilePrint(paste("LSTExperimentSimCode: RunOneSimExperiment Now Run HitAFunctionForSimulation",
       .self$OnNameFunctions[IFunction], sep="")); flush.console();
   }
   if (FALSE) {
     try(AHit <- HitASimulation(UniqueProcessIdentifier = .self$UniqueProcessIdentifier,
       UniqueSimulationIdentifier = SMS$UniqueSimulationIdentifier, 
       NameFunction = .self$OnNameFunctions[IFunction], 
       Hit = Hit,
       quoteMore = paste(quoteMore, " - ", sep=""), 
       verbose = max(verbose -1, verbose)));
   }
   if (ForceAFail == TRUE) {
     AFilePrint(paste("Note We are deliberately choosing to count ",
       .self$OnNameFunctions[IFunction], " now as a fail.", sep=""));
      Hit <- "D";
   }
   AHit <- NULL;
   try(AHit <- HitAFunctionForSimulation(UniqueProcessIdentifier = .self$UniqueProcessIdentifier,
     UniqueSimulationIdentifier = SMS$UniqueSimulationIdentifier, 
     NameFunction = .self$OnNameFunctions[IFunction], 
     Hit = Hit,
     quoteMore = paste(quoteMore, " - ", sep=""), 
     verbose = max(verbose -1, verbose)));

   if (verbose >= 1) {
     AFilePrint(paste("LSTExperimentSimCode: RunOneSimExperiment We have finished HitASimulation on ",
       .self$OnNameFunctions[IFunction], " getting AHit returned as : ", AHit, sep="")); flush.console();
   }
   try(setwd(.self$OriginalOldwd));
   if (!exists("AHit") || is.null(AHit) || as.character(AHit) == "FailWrite") {
      AFilePrint(paste("We had a  Write Fail at  ", .self$OnNameFunctions[IFunction], sep=""));
      eval(parse(text=SetGText("limlas", "globalenv()", S=1)));
      eval(parse(text=SetGText("MyV", "globalenv()", S=1)));
      flush.console(); return(-1);
   } else if (is.numeric(AHit) && (AHit == -1)) {
     AFilePrint(paste("We Have a HitSimulation NumericFail at ", .self$OnNameFunctions[IFunction], sep=""));
     flush.console(); return(-1);
   }
   if (verbose >= 1) {
     AFilePrint(paste("LSTExperimentSimCode: Try Read a new FileL, ", sep=""));flush.console();
     AFilePrint(paste(" IFunction=", IFunction,
       ", NameFunction = ", .self$OnNameFunctions[IFunction], sep="")); flush.console();
      AFilePrint(paste(" --- currentcontaindir = ", CurrentSmallContainDir, sep="")); 
        flush.console();
      AFilePrint(paste(" --- TableOutName = ", .self$TableOutName, sep="")); flush.console();
      AFilePrint(paste(" AlreadyLocked = TRUE", sep="")); flush.console();
   }                       
   FileL = NULL;
   try(FileL <- ReadOneFile(OnNameFunction = .self$OnNameFunctions[IFunction], 
     AContainDir = CurrentSmallContainDir, 
     RDataTableOutName = .self$RDataTableOutName,
     TableOutName = .self$TableOutName, 
     verbose = max(verbose-1,0),
     AlreadyLocked = FALSE));
   if (verbose >= 1) {
     AFilePrint(paste("LSTExperimentSimCode: Finished ReadOneFile", sep="")); flush.console();
   }
   if (is.null(FileL)) {
     AFilePrint("RunOneSimExperiment: Error FileL returned NULL ");  flush.console();
     AFilePrint(paste("  Current Function = ", .self$OnNameFunctions[IFunction], sep="")); flush.console();
     return(-1);
   }
   if (is.null(FileL$FileCon) || length(FileL$FileCon)<=1) {
     AFilePrint("RunOneSimExperiment: Error FileL$FileCon is effed up returned NULL ");  flush.console();
     eval(parse(text=SetGText("FileL")));
     AFilePrint(paste("  Current Function was ", .self$OnNameFunctions[IFunction], sep="")); flush.console();
     ATText <- "
     OnNameFunction <- .self$OnNameFunctions[IFunction]; 
     ASmallContainDir <- CurrentSmallContainDir;
     RDataTableOutName <- .self$RDataTableOutName;
     ";
     try(eval(parse(text=ATText)));
     eval(parse(text=SetGText("OnNameFunction", "globalenv()", S=1)));
     eval(parse(text=SetGText("ASmallContainDir", "globalenv()", S=1)));
     eval(parse(text=SetGText("RDataTableOutName", "globalenv()", S=1)));
     try(eval(parse(text=SetGText("AlreadyLocked", "globalenv()", S=1))));
     return(-1);
   }
   if (verbose >= 1) {
     ##AFilePrint("RunOneSimExperiment: Set new part .self$FileCons, FileLenS!"); flush.console();
     AFilePrint("RunOneSimExperiment: Ignoring FileCons Setup"); flush.console();
   }
   if (FALSE && IFunction > length(.self$FileCons)) {
     AFilePrint(paste("RunOneSimExperiment: I See a weird issue, IFunction = ",
       IFunction, " but length(.self$FileCons) = ", length(.self$FileCons),
       ", what will happen? ", sep="")); flush.console();
   }
   if (FALSE) {
   ## We're separating here the running of experiments with the running of analyses
   ##  
   ## The SummarySimulations List measures completed tasks
   ##
   ## Eventually we'll recover Give FileCons properties for TableCompareSaver
   ##  However, the TableCompareSaver will not be storing results duirng the
   ##  Experimental phase.  We need to write post software that accomplishes this.
   ##
   TestSetFileL <- "Fun = -1; 
   try(.self$FileCons[[IFunction]] <- -1);
   if (!is.null(FileL$FileCon)) {
     try(.self$FileCons[[IFunction]] <- FileL$FileCon);
   }
   try(.self$FileMyDirectories[IFunction] <- FileL$FileMyDirectory);
   try(.self$FileUniqueProcessIds[[IFunction]] <- FileL$FileUniqueProcessIds);
   try(.self$FileUniqueSimulationIds[[IFunction]] <- FileL$FileMySimulationIds);
   Fun = 1";  
   try(eval(parse(text=TestSetFileL)));
   if (Fun == -1) {
     AFilePrint("ERROR ERROR ERROR ERROR ERROR ERROR ERROR");  flush.console();
     AFilePrint(paste("ERROR: RunOne Sim Experiment.", sep="")); flush.console();
     AFilePrint(paste("ERROR In setting FileCons to IFunction = ", IFunction, "!"));
     AFilePrint(paste("length(FileCons) = ", length(.self$FileCons), sep=""));
     AFilePrint("Current FileL$FileCon is "); flush.console();
     AFilePrint(FileL$FileCon);
     AFilePrint(".self$FileCons is "); flush.console();
     AFilePrint(.self$FileCons); flush.console();
     AFilePrint("Figure Out why this failed?"); flush.console();
     MYFAILFLAG = NULL;
     rm(MYFAILFLAG);
     NOTGOOD = MYFAILFLAG;
   }
   if (IFunction > length(.self$FileLenS)) {
     AFilePrint(paste("RunOneSimExperiment: I See a weird issue, IFunction = ",
       IFunction, " but length(.self$FileLenS) = ", length(.self$FileLenS),
       ", what will happen? ", sep="")); flush.console();
   }
   try(.self$FileLenS[IFunction] <- as.integer(FileL$FileSLen));
   }
   if (verbose >= 3) {
     try(AFilePrint(paste("LSTExperimentSimCode: RunOneSimExperiment created limlas", sep=""))); flush.console();
   }
   if (AHit == "+" && is.logical(DoCI) && DoCI == TRUE  && !is.null(FittedOb$CIEst) ) {
     try(.self$AssessCIFit(SMS = SMS, FittedOb = FittedOb,
     UniqueProcessIdentifier = .self$UniqueProcessIdentifier,
     UniqueSimulationIdentifier = SMS$UniqueSimulationIdentifier, 
     NameFunction = .self$OnNameFunctions[IFunction], 
     Hit = Hit,
     quoteMore = paste(quoteMore, " - ", sep=""), 
     verbose = max(verbose -1, verbose), AlreadyLocked=FALSE)); 
   }
   
   if (AHit == "+" && !is.null(FittedOb$MIPReport) ) {
     try(.self$AssessMIPFit(SMS = SMS, FittedOb = FittedOb,
     UniqueProcessIdentifier = .self$UniqueProcessIdentifier,
     UniqueSimulationIdentifier = SMS$UniqueSimulationIdentifier, 
     NameFunction = .self$OnNameFunctions[IFunction], 
     Hit = Hit,
     quoteMore = paste(quoteMore, " - ", sep=""), 
     verbose = max(verbose -1, verbose), AlreadyLocked=FALSE)); 
   }
   
   ## Nothing to unlock here
   ##UnLockMe(verbose=verbose, quoteMore=quoteMore, LFDir = MDir);
   if (verbose >= 3) {
     AFilePrint("RunOneSimExperiment:  All Finished, return(1)!"); flush.console();
   } 
   try(setwd(.self$ASmallContainDir)); 
   try(UnLockMyProcess(verbose=0, quoteMore="UnLock SumPerformance on NoBS", LFDir = "SumPerformance",
        NoSpaceLFDir =FALSE, SumLock="No", LockIn="No"));
   try(UnLockMyProcess(verbose=0, quoteMore="UnLock SumPerformance on NoBS", LFDir = .self$OnNameFunctions[IFunction],
        NoSpaceLFDir =FALSE, SumLock="No", LockIn="No"));
   if (any(names(limlas) == "GibbsSS")) {
     try(limlas$GibbsSS$TBSR5$unsave(TRUE));
     try(limlas$GibbsSS <- NULL);
   }
   if (any(names(limlas) == "BackGibbsSS")) {
     try(limlas$BackGibbsSS$TBSR5$unsave(TRUE));
     try(limlas$BackGibbsSS <- NULL);
   }
   try(setwd(.self$OriginalOldwd));
   return(1);
   if (!is.null(limlas$CI)) {
    if (length(IntAssess)  == 1 && is.numeric(IntAssess) && IntAssess[1] == -1) {
      IntAssess = list();
      eval(parse(text=GetGText("LenIntAssess", S=1)));
      LenIntAsses = 1;
      eval(parse(text=SetGText("LenIntAssess", S=1)));
      IntAssess[[1]] = list(iiint = IFunction, nameF = .self$OnNameFunctions[IFunction],
      AssessedCI = AssessCI(limlas$CI, beta));     
    } else {
      eval(parse(text=GetG0Text("LenIntAssess", S=1)));
      if (LenIntAssess != length(IntAssess)) {
        AFilePrint(paste("RunOneSimExperiment: Weird, LenIntAsses = ",
          LenIntAssess, " but length(IntAssess) = ", length(IntAssess),
          ", wonder what will happen?", sep="")); flush.console();
      }
      LenIntAssess =  LenIntAssess+1;
      IntAssess[[LenIntAssess]] = list(iiint = IFunction, nameF=.self$OnNameFunctions[IFunction],
      AssessedCI = AssessCI(limlas$CI, beta));
      eval(parse(text=SetGText("LenIntAssess",S=1)));
    }
  }
  Retter = list(LS=NULL, IntAssess=IntAssess);   
  return(Retter);
});
