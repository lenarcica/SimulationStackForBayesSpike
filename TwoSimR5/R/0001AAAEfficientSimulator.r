###############################################################################
##  EfficientSimulator  -- Alan Lenarcic
##
##    01-12-10 - 2019
##
##   There are many estimators in R for sparse fixed model selection, and
##   this code is intended to try and unify the input formats across a number
##   of estimators.  All of these estimators were current as of the last
##   time I had access to UNC Killdevil in Summer 2017.  However,
##   not all estimator results were used in Lenarcic and Valdar paper, mostly
##   because better or more canonical examples produced similar results.
##   Clearly, this covers only a tiny subset of available options, but
##   a purpose of this was to focus on Bayesian and Penalized estimators
##   which had a sparse setting, and to look specifically in Type 1 Type 2 
##   performance when the Beta versus noise setting was .25 up to 2.5.
##   All of these would be R packages with implementations sufficient for single
##   thread analyses, and no need for CUDA or multi-thread to compete. At the
##   time of the TwoSimR5 development (2010-2016), most methods fell into this category.
##
##   While it is hopeful that code like these simulations could help other
##    projects doing simulation performance tests in some other region of
##    the p,n,k,sigma,Beta,Cov(X), etc. choices every group likely has
##    a different region of interest and specialty, and TwoSimR5 seems
##    to demonstrate a role for all of these estimators, depending on the
##    location of noise space one belongs to, as well as the amount of prior
##    information one can come into the experiment holding.
##
##    Helper File for 221 project
##
##    This file corrects issues in TwoLasso package "SimForMe2" and represents
##      a fair testing situation for producing Stodden plots comparing TwoLasso
##      estimators to LARS estimates
##
##    For the most part you will only be using 
##      "GenerateLarsFixed" and "GenerateLimitLasso"
##      estimates and plotting to graphs.  However the function
##                       "SimulateCompare()"
##      gives a method for generating a table rating performance of
##      all of the LimitLasso estimators.
##    
##   SimMeData is a function at the end that generates data to rate LimitLasso
##     and TwoLasso against each other
##
##   Key inputs are a function "CorrelationXmatrix" which generates a 
##     Correlation matrix.  Choose one of CorrelationXmatrix1 or
##     CorrleationXmatrix2, and consider reasonable inputs.
##     DefaultCorrelationXmatrix is set yere to CorrelationXmatrix1
##     with an input of PosCorr = .3
##
##   Generating the Betas vector comes from a function declared
##     "GenerateBetaVec".  Choose one of the GenerateBetaVec functions
##     from GenerateBetaVec1, GenerateBetaVec2, GenerateBetaVec3
##     default is GenerateBetaVec2
##
##   An issue is that I've discovered some programming issue in underlying 
##     C code.  Though the main algorithms run as I intended them to,
##     I wrote my own helper version of the LARS algorithm in the function
##     LarsCC2() which is in desperate need of Lapack optimizing.  I'll be
##     working on an update of the underlying package C++ code for the next
##     few weeks.  So depending on how the preliminary results come out on this
##     assignment, the best runs will get into the next package update,
##     possibly getting in a publication.
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
BetaStart <<- 0;


  ## Covers the possible R packages that should be installed to get best use out of TwoSim.
## source("c://Stat//2008Summer//LarsProject//Code//EfficientSimulator.R");
## source("C://Stat//2008Summer//LarsProject//TwoLassoPackage//ablocker_ps2//jtermFDL//jtermFDL//EfficientSimulator.R")
   try(library(lars, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
   try(library(corpcor, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
   ##try(library(TwoLassoCpp, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);  
     ## If TwoLasso has installed correctly.
   try(library(ncvreg, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);  ##SCAD is located Here
   try(library(elasticnet, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE); ## elastic net enet is here
   try(library(glmpath, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);  ## more elastic net
   try(library(spikeslab, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);  ## Ishwaran Spike
   try(library(corpcor, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
   try(library(monomvn, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
   

   
   ALT3 <- 10;  ALT2 <- 10 * 11 / 2;
   
   
   LARSSEEKMETHOD <<- 1;
   STANDARDFLAG <<- 0;
   tNoiseDF <<- 0; DefaultNoiseDF <<- 0;

   #######################################################################
   ###   SimKActiveLen;  SimSigmaNoise;
   ###
   ### To Give False Input to simulations, we give
   ###
   piSV = 0;
   SimKActiveLen <- function (kActiveLen, NLen, kLen)  {
            min(c(NLen-1, kLen -1, max(round(kActiveLen + rnorm(1, 0, sqrt(kLen)*piSV)), 1)));
   }
   sigSV = 0;
    SimSigmaNoise <- function (SigmaNoise) {
       if (sigSV > 0) {
         return(sqrt(SigmaNoise^2 * rchisq(1, sigSV) /sigSV));
       } else{
         return(SigmaNoise);
       }
    }


StableHPDinterval <- function(MCMCThing, prob=.05) {
  library(coda);
  if (length(MCMCThing) == 0) {
    if (NCOL(MCMCThing) >= 1) {
      AR <- (matrix(0, NCOL(MCMCThing),2));
      try(rownames(AR) <- colnames(MCMCThing));
      return(AR);
    }
    return(c(0,0));
  }
  ART <- NULL;
  if (!is.mcmc(MCMCThing) && !is.list(MCMCThing)) {
    try(MCMCThing <- suppressWarnings(as.mcmc(MCMCThing)), silent=TRUE);
  }
  try(ART <- HPDinterval(MCMCThing, prob=prob), silent=TRUE);
  if (is.null(ART)) {
    ART <- matrix(0, NCOL(MCMCThing), 2);
    rownames(ART) <- colnames(MCMCThing);
    for (jj in 1:NCOL(MCMCThing)) {
      BBT <- NULL;
      try(BBT <- HPDinterval(MCMCThing[,jj], prob=prob), silent=TRUE);
      if (!is.null(BBT) && !all(is.na(BBT)) && length(BBT) == 2) {
        ART[jj,] <- BBT;
      }
    }
  } else {
    return(ART);
  }
  return(ART);
}



   
CorrelationXmatrix1 <- function(kLen, PosCorr) {
  if (abs(PosCorr) > 1) {
      AFilePrint("CorrleationXmatrix1, PosCorr is too large, doesn't work\n");
  }
  if (PosCorr == 0) {
        list(
    corrX = diag(kLen),
    type="CorrelationXmatrix1", PosCorr = PosCorr )
  }
  return( list(
    corrX = PosCorr^( abs( matrix(rep((1:kLen),kLen),kLen,kLen) - 
                          t(matrix(rep((1:kLen),kLen), kLen,kLen)  ))  ),
    type="CorrelationXmatrix1", PosCorr = PosCorr )
    );

}
SinerX <- function(Inpts, a,b) {
   rt <-  sin(a * Inpts) / ( a * Inpts);
}
CorrelationXmatrix2 <- function(kLen, a,b) {
  if (abs(PosCor) > 1) {
      AFilePrint("CorrleationXmatrix1, PosCorr is too large, doesn't work\n");
  }
  return( list(
     corrX = SinerX( abs( matrix(rep((1:kLen),kLen),kLen,kLen) - 
                          t(matrix(rep((1:kLen),kLen), kLen,kLen)  )) , a,b  ),
     type="CorrelationXmatrix2", a=a, b=b
     )
  );
}

SStandardizeXX   <- function(XX, meanP = FALSE, SumToOne = TRUE) {
  if (meanP == FALSE && SumToOne == TRUE) {
     if (length(dim(XX)) == 2) {
        XX <- t( ( t(XX) / sqrt(colSums(XX^2) ))) 
     }   else {
        XX <- XX / sqrt(sum(XX^2));
     }
     return(XX);
  }
  if (meanP == FALSE)  {
       if ( length(dim(XX)) == 2) {
             XX <- t( ( t(XX) / apply(XX,2,sd) ));
       }  else {
         if (is.matrix(XX)) {
            XX <-  (XX) / as.vector(XX);
         } else {
           XX <- (XX) / sd(as.vector(XX));
         }        
       }
       return(XX);  
  }
  if (SumToOne == TRUE) {
     if (length(dim(XX)) == 2) {
        XX = t( t(XX) - meanC(XX));
        XX <- t( ( t(XX) / ( apply(XX,2,sd) * (length(XX[,1]) -1) ))) 
     }   else {
        XX = XX - mean(XX);
        XX <- XX / sqrt(sum(XX^2));
     }
     return(XX);  
  }
   return(StandardizeXX(XX));
}


## GenerateBetaVec1 
  GenerateBetaVec1 <- function(kLen = 100, kAwant = 6) {
     if (kAwant > kLen) { AFilePrint("GenerateBetaVec1 Entry Error");
                          return(-1);
     }
     return(list(
      BetaVec =
        sample(
           c(rep(1, kAwant), rep(0, NROn-kAwant))
               , NROn, replace=FALSE
           )
       , type="GenerateBetaVec1",
         kLen = kLen, kAwant = kAwant
             )
        );  
  }
  GenerateBetaVec2 <- function(kLen = 100, kAwant = 6) {
     if (kAwant > kLen) { AFilePrint("GenerateBetaVec2 Entry Error");
                          return(-1);
     }
     return(list(
        BetaVec=sample(
           c(rep(1, floor(kAwant/2)),
             rep(-1, floor(kAwant/2)),
             rep(1, kAwant - 2 * floor(kAwant/2)), 
             rep(0, kLen-kAwant) )
             , kLen, replace=FALSE )
        , type="GenerateBetaVec2",
           kLen = kLen, kAwant = kAwant
             )             
           ); 
  }
  GenerateBetaVec2a <- function(kLen = 100, kAwant = 6) {
     if (kAwant > kLen) { AFilePrint("GenerateBetaVec2 Entry Error");
                          return(-1);
     }
     fLL = floor(kAwant/2);
     if (kAwant %% 2 == 1) {
        BetaVec = c( t( cbind(rep(1,fLL),rep(-1,fLL)) ), 
                      1, rep(0, kLen-kAwant)
                   )
     } else {
        BetaVec = c( t( cbind(rep(1,fLL),rep(-1,fLL)) ), 
                       rep(0, kLen-kAwant)
                   )     
     }
     return(list(
        BetaVec=BetaVec,
           type="GenerateBetaVec2a",
           kLen = kLen, kAwant = kAwant
             )             
           ); 
  }  
  
    GenerateBetaVec2b <- function(kLen = 100, kAwant = 6) {
     if (kAwant > kLen) { AFilePrint("GenerateBetaVec2 Entry Error");
                          return(-1);
     }
     fLL = floor(kAwant/2);
     ILLoc <- sample(1:kLen, kAwant, replace=FALSE);
     if (kAwant %% 2 == 1) {
        BetaVec = rep(0, kLen);
        BetaVec[ILLoc[1: (kAwant %/% 2 + 1)]] = 1;
                BetaVec[ILLoc[((kAwant %/% 2 + 2)):kAwant]] = -1;   
     } else {
        BetaVec = rep(0, kLen);
        BetaVec[ILLoc[1: (kAwant %/% 2 )]] = 1;
                BetaVec[ILLoc[((kAwant %/% 2) + 1):kAwant]] = -1;      
     }
     return(list(
        BetaVec=BetaVec,
           type="GenerateBetaVec2b",
           kLen = kLen, kAwant = kAwant
             )             
           ); 
  } 
  
GenerateBetaVec3 <- function(kLen = 100, kAwant = 10, MN = 100) {
  if (kAwant > kLen) { AFilePrint("GenerateBetaVec3 Entry Error");
                       return(-1);
  }
  return(list(
    BetaVec= sample( c( runif(kAwant, -MN,MN), rep(0, NROn-kAwant) ),
      kLen, replace=FALSE )), 
      type="GenerateBetaVec3", kLen = kLen, kAwant = kAwant, MN = MN);  
  }

GetBetaName <- function(GenerateBetaVec) {
  BetaVD <- GenerateBetaVec(4,1);
    BetaName <- BetaVD$type;
    if (length(BetaVD) > 4) {
      for (ii in 5:length(BetaVD)) {
        BetaName <- paste(BetaName, names(BetaVD)[ii],
                                    tSeq(BetaVD[[ii]]), sep="" );
      }
    }
    return(BetaName);
}
GetCorrName <- function(CorrelationXmatrix) {
   CorrVD <- CorrelationXmatrix(4);
   CorrName <- CorrVD$type
      for (ii in 3:length(CorrVD)) {
        CorrName <- paste(CorrName, names(CorrVD)[ii],
                                    tSeq(CorrVD[[ii]])
                     ,sep="" );
      }
    return(CorrName);
}
if (!exists("tSeq")) {
  tSeq <- function (Number) 
{
    CN <- sub("e", "e", as.character(Number))
    CN <- sub("-", "n", as.character(CN))
    CN <- sub("\\+", "f", as.character(CN))
    CN <- sub("\\.", "p", as.character(CN))
    CN <- sub("0p", "p", as.character(CN))
    return(CN)
}
}

#########################################################################################
## FileNamingClenture
##
##   Naming standard for Files
##
FileNamingClenture <- function(GenerateBetaVec, CorrelationXmatrix,
      NLen, kLen, kActiveLen, SigmaNoise, NumReps, LARSSeekFlag) {
  BetaName <- GetBetaName(GenerateBetaVec);
  CorrName <- GetCorrName(CorrelationXmatrix);
	  NC <- paste(BetaName, CorrName, "NLen", NLen,
	         "kLen", kLen, "kALen", kActiveLen,
	         "LSF", tSeq(LARSSeekFlag),
	         "SQN", tSeq(SigmaNoise),
	         "Reps", NumReps, 
	         ".csv", sep="");
  return(NC);
} 

DirNamingClenture <- function(GenerateBetaVec, CorrelationXmatrix,
      NLen, kLen, kActiveLen, SigmaNoise, NumReps, LARSSeekFlag) {
  BetaName <- GetBetaName(GenerateBetaVec);
  CorrName <- GetCorrName(CorrelationXmatrix);
	  NC <- paste(BetaName, CorrName, "NLen", NLen,
	         "kLen", kLen, "kALen", kActiveLen,
	         "LSF", tSeq(LARSSeekFlag),
	         "SQN", tSeq(SigmaNoise),
	         "Reps", NumReps, sep="");
  return(NC);
}   

###############################################################################
FileNameGenerator <- function( DIRECTORY, 
       GenerateBetaVec, CorrelationXmatrix,
      NLen, kLen, kActiveLen, SigmaNoise, NumReps, LARSSeekFlag ) {
   FileName <-  FileNamingClenture(GenerateBetaVec, CorrelationXmatrix,
      NLen, kLen, kActiveLen, SigmaNoise, NumReps, LARSSeekFlag);
  paste(DIRECTORY, "//", FileName, sep="");
}

  
  DefaultGenerateBetaVec <- GenerateBetaVec2;
  DefaultCorrelationXmatrix <- function(kLen) {
           CorrelationXmatrix1(kLen, .25);
  }
  NLen = 20; kLen = 10; SigmaNoise = .5; NumReps = 500;
  kActiveLen = 4;
  LARSSeekFlag = 3;
  musize =1; SigmaEta = 2; SigmaBar = -999;
  FileName = "c://Stat//Playfile.csv";
  PrintOutFlags = 1;
  GenerateBetaVec = DefaultGenerateBetaVec
  CorrelationXmatrix = DefaultCorrelationXmatrix;
  
EstimatorList = c("LarsFixed", "LarsCp", "LassoQLY",
                  "LimitRidge",
                  "TwoLasso9x", "LimitLasso", 
                  "piALimitLasso", "FDLimitLasso",
                  "XLMMedian", "XLLimitLasso");
  InputTable <- NULL;
  ## If InputTable has been read before, retrieve with:
  ## InputTable <- read.table(paste(FileName,"a",sep=""), sep=",", header=TRUE);
##
SimulateCompare <- function(NLen, kLen, SigmaNoise, NumReps, 
                            kActiveLen,
                            LARSSeekFlag = 3,
                            GenerateBetaVec = DefaultGenerateBetaVec,
                            CorrelationXmatrix = DefaultCorrelationXmatrix,
                            musize = 10, SigmaEta = 2, SigmaBar =-999, 
                            FileName = NULL, PrintOutFlags = 0
                           ) 
{ 
   PPract <- kActiveLen / kLen;
      
   ReturnMatrix <- MakeReturnMatrix(InputTable, EstimatorList, NumReps);
   ii4Start <- min((1:length(ReturnMatrix[,1]))[ ReturnMatrix[,1] == 0 ] );
   
   
      ## Just a TestSim        
	     SMS <- SimMeData(NLen = NLen, kLen = kLen, kActiveLen = kActiveLen,
       SigmaNoise = SigmaNoise, GenerateBetaVec = GenerateBetaVec, 
       CorrelationXmatrix = CorrelationXmatrix);
      if (SigmaEta > 0) {
         SigmaBar <- SigmaNoise^2;
      }  
      
   ## Start Loop                 
for (ii4 in ii4Start:NumReps) {
        if (PrintOutFlags >= 1) {
          AFilePrint(paste("ii4 = ", ii4, ", starting beginning"));
          flush.console();
        }
       ##Simulate the Data
	     SMS <- SimMeData(NLen = NLen, kLen = kLen, kActiveLen = kActiveLen,
       SigmaNoise = SigmaNoise, GenerateBetaVec = GenerateBetaVec, 
       CorrelationXmatrix = CorrelationXmatrix)

	     if (SigmaBar == -999) {
	       SigmaBarS <- (.5 * sd(as.vector(SMS$YY)))^2
	     } else {
	       SigmaBarS <- SigmaBar;
	     }

      ## Generate Lars fit with our own implementation to check
      HitLars <- LarsCC2(SMS$XX,SMS$YY, BetStart = -1, lambda = 1/ BetaBarMin * 
	                log((1-SMS$puse)/SMS$puse), BetOlds = -1,
	                     StandardFlag = 1); 	     
	  FitL <- list();
	  ##1 Generate the Lars Fit
    FitL[[1]]<- FitLarsFixed <- GenerateLarsFixed(SMS);
    FitL[[2]]<-  FitLarsCp <- GenerateLarsCp(SMS);
    FitL[[3]]<-  FitLassoQLY <- GenerateLassoQLY(SMS);
 	      if (PrintOutFlags >= 2) {  
	        AFilePrint(paste("ii4 = ", ii4, ", starting with EMRIDGE"));
	          flush.console();
	      }
    FitL[[4]]<-  FitEMRidge <- GenerateEMRidge(SMS, LARSSeekFlag);	 
    FitL[[5]]<- FitTwoLasso9X <- GenerateTwoLasso9X(SMS, LARSSeekFlag); 
        if (PrintOutFlags >= 2) {  
	        AFilePrint(paste("ii4 = ", ii4, ", starting with LimitLasso"));
	          flush.console();
	      }      
   FitL[[6]]<-   FitLimitLasso <- GenerateLimitLasso(SMS, LARSSeekFlag);
        if (PrintOutFlags >= 2) {  
	        AFilePrint(paste("ii4 = ", ii4, ", starting with piALimitLasso"));
	          flush.console();
	      } 
   FitL[[7]]<-   FitpiALimitLasso <- GeneratepiALimitLasso(SMS, 
     LARSSeekFlag, musize, SigmaEta, SigmaBarS);  
	      if (PrintOutFlags >= 2) {  
	        AFilePrint(paste("ii4 = ", ii4, ", starting with FDLimitLasso"));
	          flush.console();
	      }	                                            
   FitL[[8]]<-   FitFDLimitLasso <- GenerateFDLimitLasso(SMS
                                            , LARSSeekFlag);       
   FitL[[9]]<-   FitXLMMedian   <- GenerateXLMMedian(SMS);                                       
   FitL[[10]]<-   FitXLLimitLasso   <- GenerateXLLimitLasso(SMS,FitLimitLasso); 
   
   Header = 8;
   ReturnMatrix[ii4,1:Header] <- c( SMS$BetasName, SMS$CorrName, SMS$NLen,
            SMS$kLen, SMS$kActiveLen, SMS$SigmaNoise, NumReps,
            LARSSeekFlag)
      for (iiFit in 1:length(FitL)) {
         Type1 <- length(SMS$BBsReal[SMS$BBsReal - FitL[[iiFit]]$BBHatsFit 
                                 < 0]);
         Type2 <- length(SMS$BBsReal[SMS$BBsReal - FitL[[iiFit]]$BBHatsFit 
                                 > 0]);
         Distance <- sum( (SMS$BetasReal - FitL[[iiFit]]$BetaFit)^2 );  
         Time <- FitL[[iiFit]]$FitTime                             
         ReturnMatrix[ii4, Header  + iiFit + 
             (0:5) * length(FitL)] <- c(Type1, Type2, Distance, Time[1:3]);
          
      } 
    write.table( file = paste(FileName, "a", sep=""), 
        x=ReturnMatrix[1:ii4,], sep=", ", eol="\n",
        na ="NA", row.names=FALSE,col.names=TRUE);
   } 
      write.table( file = paste(FileName, sep=""), 
        x=ReturnMatrix, sep=", ", eol="\n",
        na ="NA", row.names=FALSE,col.names=TRUE);  
  return(ReturnMatrix);
}   

MakeReturnMatrix <- function(InputTable, EstimatorList, NumReps) {
    Header = 8;
    ReturnMatrix <- matrix(0, NumReps, 
                           Header+length(EstimatorList) * (2 + 1 + 3));
    colnames <- c("BetaGen", "CorrGen", "NLen", "kLen", "kActiveLen", 
                   "SigmaNoise", "NumReps", "LARSSeekFlag",
                  paste( "Type1", EstimatorList,sep=""),
                  paste( "Type2", EstimatorList,sep=""),
                  paste( "SumDeltaBetaSq", EstimatorList, sep=""),
                  paste( "t1", EstimatorList, sep="" ),
                  paste( "t2", EstimatorList, sep="" ),
                  paste( "t3", EstimatorList, sep="" )
                  );  
    colnames(ReturnMatrix) <- colnames        
    ## If there is an input table for data, puts data inside it.
    if (is.null(InputTable)) {  
      } else if ( length(InputTable[1,]) != length(ReturnMatrix[1,])) {
         AFilePrint(paste("SimForMe2: InputTable Not shaped like Return Matrix,",
                     "  dim(InputTable) = (",
                      dim(InputTable)[1], ", ", dim(InputTable)[2], ")", sep=""));
         AFilePrint(paste("FileName: ", FileName, sep=""));      
      } else {
         LIT <- length(InputTable[InputTable[,1] != 0,1]);
         ReturnMatrix[1:LIT,] <- 
             unlist(InputTable[1:LIT,]);
         ii4Start <- LIT +1;
         AFilePrint(paste("SimForMe2: InputTable for has nonzero length = ", 
                                          ii4Start-1, " of ", NumReps, sep=""));
         AFilePrint(paste("FileName: ", FileName, sep=""));
    } 
    return(ReturnMatrix);                           
}

################################################################################
### Generate Fits
###
###   These functions are designed to take SMS simulations and then generate
###    fits using R and C++ code based upon Lasso and TwoLasso estimators
###
###   An issue with the TwoLasso fits is that they start off with a relatively
###    slow implementation of LARS, which will be bad when kActiveLen >= 50
###    (I was stupid, didn't know about Lapack, so LarsCC2 uses self coded
###     row-reduction based matrix inversion).  If I can fix LarsCC2 in time
###     for the assignment it will be a miracle, but I will try.
###

MinMin  = .005;
GenerateLarsFixed <- function(SMS, DoCI = DefaultDoCI,...) {
  kLen = SMS$kLen;
  T1 <- proc.time();
	library(lars, warn.conflicts = FALSE, quietly=TRUE);
	if (!is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
	if (kLen > 500) {
	  HitLars2Ob <- lars(SMS$XX, SMS$YY, use.Gram=FALSE);
	} else {
 	  HitLars2Ob <- lars(SMS$XX, SMS$YY);     
  }
  } else {
     HitLars2Ob = glmpath(SMS$XX, SMS$YY,
      nopenalty.subset = NULL, family = binomial,
      weight = rep(1, n), offset = rep(0, n), lambda2 = 1e-5,
      max.steps = 10*min(n, m), max.norm = 100*m,
               min.lambda = (if (m >= n) 1e-6 else 0), 
      max.vars = min(length(SMS$YY)-1, length(SMS$XX[1,])),
      max.arclength = Inf, frac.arclength = 1, add.newvars = 1,
      bshoot.threshold = 0.1, relax.lambda = 1e-8,
      standardize = TRUE, function.precision = 3e-13,
      eps = .Machine$double.eps, trace = FALSE);
     HitLars2Ob$beta = HitLars2Ob$b.predictor;
  }
  
	if (is.null(HitLars2Ob$beta)) {
    Lars2Fit <- rep(0, kLen);
    Lars2Fit2 = Lars2Fit;    
  } else {
	  if (length(HitLars2Ob$beta[,1]) < SMS$kActiveLen) {
	    ITO = length(HitLars2Ob$beta[,1]);
	  } else {
	    if (kLen > 1000) {
	      ITO <- LarsGiveMeClose(HitLars2Ob, IntClose = round(SMS$kActiveLen), 
          MinMin = .00001); 
      } else {
 	      ITO <- LarsGiveMeClose(HitLars2Ob, IntClose = round(SMS$kActiveLen), 
          MinMin = .00001);          
      } 
    }
    if (ITO <= 0) {ITO = 1;}
    if (ITO > length(HitLars2Ob$beta[,1]))  { 
      ITO =length(HitLars2Ob$beta[,1]);
    } 
	  Lars2Fit <- HitLars2Ob$beta[ITO,];
	  Lars2Fit[abs(Lars2Fit) > MinMin] <- 1;
	  Lars2Fit[abs(Lars2Fit) <= MinMin] <- 0;  
	  Lars2Fit2 <- Lars2Fit;
    BetaSM <- 0;
	  if (sum(Lars2Fit) > 0 || length(SMS$XX[1,Lars2Fit == 1]) > 0) { 
	    XXSM <-  SMS$XX[,Lars2Fit == 1];
	    XXSXX <- t(XXSM) %*% XXSM;
	    if (any(is.null(XXSXX)) || 
        length(dim(XXSXX)) != 2 || dim(XXSXX)[1] != dim(XXSXX)[2] || 
	      dim(XXSXX)[1] < 1 || any(is.nan(XXSXX))) { 
	      Lars2Fit2 <- Lars2Fit; BetaSM <- 0;
	    } else if (length(XXSXX) == 1) {
	      BetaSM <- 1 / XXSXX * sum(XXSM * SMS$YY);
	    } else if ( !is.null(dim(XXSXX)) && dim(XXSXX)[1] > 1 && 
        dim(XXSXX)[1] == dim(XXSXX)[2] && det(XXSXX) <= 0.00000001 ) {
	      library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
	      BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;
	    } else if (!is.null(dim(XXSXX)) && det(XXSXX) > .00000001) {
		    LS <- svd(XXSXX);
		    if (max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
          BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;  
        } else {
		      BetaSM <- try(solve(XXSXX) %*% t(XXSM) %*% SMS$YY);
		    }		        
	    } else {
 	      library(corpcor, warn.conflicts = FALSE, quietly=TRUE);
	      BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;         
      }
	      Lars2Fit2[Lars2Fit != 0 ] <- BetaSM;
	    } else {
	      Lars2Fit2 <- Lars2Fit;  BetaSM <- 0;
	    }
    }
	T2 <- proc.time();
	Lars2FitTime <- T2 - T1;
	return(list( type="LarsFixed", BetaFit = Lars2Fit2, BBHatsFit = Lars2Fit,
	  FitTime = Lars2FitTime, OtherBetas = HitLars2Ob$beta[ITO,],
      HitLars2Ob=HitLars2Ob, ITO=ITO, Lars2Fit=Lars2Fit,
        Lars2Fit2=Lars2Fit2, BetaSM=BetaSM, MIPReport = Lars2Fit) );                  
}

################################################################################
### Generate ENetFixed
###
###   Gives ENet the right parameters to choose optimal shrinkage lambda2
###
MinMin  = .005;
GenerateENetFixed <- function(SMS, DoCI = DefaultDoCI,...) {
  kLen = SMS$kLen;
  ##AFilePrint("ENetFixed Run"); flush.console();
  T1 <- proc.time();
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
	library(elasticnet, warn.conflicts = FALSE, quietly=TRUE);
  } else {
    library(glmpath, warn.conflicts=FALSE, quietly=TRUE)
  }
  VarData <- var(SMS$BetasReal[SMS$BetasReal != 0]) / var(SMS$YY);
	lambda2 = (1 / VarData) / sqrt(length(SMS$YY));
	lambda2 = 1/VarData;
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE)  {
	if (kLen > 500) {
	  ##HitLars2Ob <- lars(SMS$XX, SMS$YY, use.Gram=FALSE);
	  HitENET2Ob <-  enet(SMS$XX, SMS$YY, lambda = as.numeric(lambda2),
      max.steps= SMS$kActiveLen+1, normalize=TRUE, 
      intercept=FALSE, trace = FALSE, eps = .Machine$double.eps);
	} else {
 	  ##HitLars2Ob <- lars(SMS$XX, SMS$YY);   
   	HitENET2Ob <-  enet(SMS$XX, SMS$YY, lambda = as.numeric(lambda2),
      max.steps= 2*SMS$kActiveLen+1, normalize=TRUE, intercept=FALSE, 
      trace = FALSE, eps = .Machine$double.eps);         
  }
  HitENET2Ob$beta = HitENET2Ob$beta.pure;
 
  if (is.null(HitENET2Ob$beta)) {
    ENET2Fit <- rep(0, kLen);
    ENET2Fit2 = ENET2Fit;    
  }     else {
  	if (length(HitENET2Ob$beta[,1]) < SMS$kActiveLen) {
  	  ITO = length(HitENET2Ob$beta[,1]);
  	} else {
  	  if (kLen > 1000) {
  	    ITO <- LarsGiveMeClose(HitENET2Ob, 
          IntClose = round(SMS$kActiveLen) , MinMin = .00001); 
      } else {
 	      ITO <- LarsGiveMeClose(HitENET2Ob, 
          IntClose = round(SMS$kActiveLen) , MinMin = .00001);          
      } 
    }
    if (ITO <= 0) {ITO = 1;}
    if (ITO > length(HitENET2Ob$beta[,1]))  { 
      ITO =length(HitENET2Ob$beta[,1]);
    } 
	  ENET2Fit <-  rep(0, kLen);
    ENET2Fit[ HitENET2Ob$allset[HitENET2Ob$beta[ITO,] != 0]] = 1;   
    #HitENET2Ob$beta[ITO,];
	  #ENET2Fit[abs(ENET2Fit) > MinMin] <- 1;
	  #ENET2Fit[abs(ENET2Fit) <= MinMin] <- 0;  
	  ENET2Fit2 <- ENET2Fit;
    BetaSM <- 0;
	  if (sum(ENET2Fit) > 0 || length(SMS$XX[1,ENET2Fit == 1]) > 0) { 
	    XXSM <-  SMS$XX[,ENET2Fit == 1];
	    XXSXX <- t(XXSM) %*% XXSM;
	    if (any(is.null(XXSXX)) || 
        length(dim(XXSXX)) != 2 || dim(XXSXX)[1] != dim(XXSXX)[2] || 
	      dim(XXSXX)[1] < 1 || any(is.nan(XXSXX))) { 
	      ENET2Fit2 <- ENET2Fit; BetaSM <- 0;
	    } else if (length(XXSXX) == 1) {
	      BetaSM <- 1 / XXSXX * sum(XXSM * SMS$YY);
	    } else if ( !is.null(dim(XXSXX)) && dim(XXSXX)[1] > 1 && 
        dim(XXSXX)[1] == dim(XXSXX)[2] && det(XXSXX) <= 0.00000001 ) {
	      library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
	      BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;
	    } else if (!is.null(dim(XXSXX)) && det(XXSXX) > .00000001) {
		    LS <- svd(XXSXX);
		    if (max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
          BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;  
        } else {
		      BetaSM <- try(solve(XXSXX) %*% t(XXSM) %*% SMS$YY);
		    }		        
	    } else {
 	      library(corpcor, warn.conflicts = FALSE, quietly=TRUE);
	      BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;         
      }
	    ENET2Fit2[ENET2Fit != 0 ] <- BetaSM;
	  } else {
	    ENET2Fit2 <- ENET2Fit;  BetaSM <- 0;
	  }
   }
   } else {
	  ##HitLars2Ob <- lars(SMS$XX, SMS$YY, use.Gram=FALSE);
	  HitENET2Ob <-  glmpath(SMS$XX, SMS$YY, 
       weight = rep(1, n), offset = rep(0, n), lambda2 = 1e-5,
      max.steps= 2*SMS$kActiveLen+1, standardize = FALSE, 
       trace = FALSE, eps = .Machine$double.eps);
      HitENET2Ob$beta <- HitENET2Ob$b.predictor[,2:NCOL(HitENET2Ob$b.predictor)];
      if (length(HitENET2Ob$beta[,1]) < SMS$kActiveLen) {
  	    ITO = length(HitENET2Ob$beta[,1]);
  	  } else {
        iIND <- HitENET2Ob$beta
        iIND[iIND != 0.0] <- 1;
        iS <- rowSums(iIND);
        ITO <- sort(abs(iS-SMS$kActiveLen), index=TRUE)$ix[1];	
      }
      if (ITO <= 0) {ITO = 1;}
      if (ITO > length(HitENET2Ob$beta[,1]))  { 
        ITO =length(HitENET2Ob$beta[,1]);
      } 
	  ENET2Fit <-  rep(0, kLen);
      ENET2Fit[ HitENET2Ob$beta[ITO,] != 0] = 1;  
      MYGLM <- glm(SMS$YY~SMS$XX[,ENET2Fit == 1], family = binomial(link="logit"))
      ENET2Fit2 = 0 * ENET2Fit;
      ENET2Fit2[ENET2Fit != 0.0] <- MYGLM$coefficients[2:length(MYGLM$coefficients)]; 
  }
	 T2 <- proc.time();
	 ENET2FitTime <- T2 - T1;
	 ##AFilePrint("Finished RunningENetFixed"); flush.console();
	 return(list( type="ENetFixed", BetaFit = ENET2Fit2, BBHatsFit = ENET2Fit,
	   FitTime =ENET2FitTime, OtherBetas = HitENET2Ob$beta[ITO,],
       HitENET2Ob=HitENET2Ob, ITO=ITO, ENET2Fit=ENET2Fit,
       ENET2Fit2=ENET2Fit2, BetaSM=BetaSM, MIPReport = ENET2Fit) );                  
}


################################################################################
### Generate HorseShoeOld
###
###   Gives HorseShoe Median estimate
###
MinMin  = .005;
GenerateHorseShoe <- function(SMS, tauSq = 1, 
  CountSamp = DefaultHorseShoeSamples, DoCI = DefaultDoCI,...) {
  AFilePrint("GenerateHorseShoe starting ");
  
   library(monomvn, warn.conflicts=FALSE, quietly=TRUE);
  t1 = proc.time();
  sdY = sd(as.vector(SMS$YY));
  sdX = apply(SMS$XX,2,sd);
  meanX = colMeans(SMS$XX); meanY = mean(SMS$YY);
  mX = t( (t(SMS$XX)- meanX) / sdX);
  mY = (SMS$YY-meanY)/sdY;
  RB = blasso(X = mX, y =mY, T = 1000, thin = NULL, RJ = TRUE, M = NULL,
       beta = NULL, lambda2 = 1, s2 = var(mY-mean(mY)),
       case = "hs",
       mprior = 0, rd = NULL, ab = NULL, theta=0, rao.s2 = TRUE,
       normalize = TRUE, verb = 0);
  PointEst = colMeans(RB$beta) * sdY / sdX;
  BBOnEst = PointEst;  BBOnEst[abs(PointEst) > .005] = 1;  BBOnEst[abs(PointEst) < .005] = 0;
  t2 = proc.time();
  if (DoCI == TRUE) {
    library(coda);
    BetaChains <- as.mcmc(sdY * t(RB$beta) / sdX);
    CITime1 = proc.time();
    eval(parse(text=GetG0Text("GetConfidenceIntervals")));
    ASG <- summary(BetaChains, quantiles = GetConfidenceIntervals);
    BetaSymmetricIntervals <- ASG[[2]];
    colnames(BetaSymmetricIntervals) <- GetConfidenceIntervals;
    rownames(BetaSymmetricIntervals) <- paste("Beta", 1:NROW(BetaSymmetricIntervals), sep="");
    BetaHPDIntervals <- BetaSymmetricIntervals * 0;
    if (GetConfidenceIntervals[1] == .5) {
      AK <- 1;  BetaHPDIntervals[,1] <- BetaSymmetricIntervals[,1];
    } else {AK <- 0;}
    NIntervals <- (length(GetConfidenceIntervals) - AK)/2;  AtAk <- AK;
    for (jj in 1:NIntervals) {
      CIP <- GetConfidenceIntervals[AtAk +2] - GetConfidenceIntervals[AtAk +1];
      try(HI <- TwoSimR5:::StableHPDinterval(BetaChains, CIP), silent=TRUE);
      try(BetaHPDIntervals[,AtAk+1:2] <- HI);
      AtAk <- AtAk +2;
    }
    CIEst <- list(BetaSymmetricIntervals=BetaSymmetricIntervals,
      BetaHPDIntervals=BetaHPDIntervals);
    CIQuantiles =  GetConfidenceIntervals;
    CITime2 = proc.time();
    try(CIEst <- CleanCIEst(CIEst));
    CITime = CITime2-t1;
  } else {CIEst = NULL; CIQuantiles = NULL;  CITime = NULL; }
  
  TimeBayesLasso <- t2 - t1;                      
	return(list( type="HoreShoeMonoMVN",
	                BetaFit = PointEst,
                  BBHatsFit = BBOnEst,
	                FitTime =TimeBayesLasso,
	                OtherBetas = PointEst,
                  EMLARSObject = RB,
                  CIEst=CIEst, CITime=CITime, CIQuantiles=CIQuantiles, MIPReport = BBOnEst) );
                  
  ##  A Test Code:                 
    X = matrix(rnorm(10*100),100,10); Y = as.vector(X %*% c(1,1,0,0,0,0,0,0,0,-1) + .5*rnorm(100));
    sdX = apply(X,2,sd); sdY = sd(as.vector(Y));
    mX =  (X-colMeans(X))/ apply(X,2,sd);  mY = (Y-mean(Y))/sd(as.vector(Y));
    RB = blasso(X = mX, y = mY, T = 1000, thin = NULL, RJ = TRUE, M = NULL,
       beta = NULL, lambda2 = 1, s2 = var(mY-mean(mY)),
       case = c("default"),
       mprior = 0, rd = NULL, ab = NULL, theta=0, rao.s2 = TRUE,
       normalize = TRUE, verb = 1)
   PointEst = sdY       
   
   
  kLen = SMS$kLen;
  library(TwoLassoCpp);
  AFilePrint("Starting SampleGibbsHorseShoe ");
  ##AFilePrint("ENetFixed Run"); flush.console();
  T1 <- proc.time();
  GibbHS<- TwoLassoCpp:::SampleGibbsHorseShoe(X=SMS$X, Y=SMS$Y, 
    tauSq, LengthSamples=CountSamp, SigmaSq=SMS$SigmaSqNoise, printFlag=10);
	HL <- GibbHS$RecordBetajs[,round(CountSamp/4):CountSamp];
	MakeAns = HL[,1];
  for (ii in 1:length(HL[,1])) {
    MakeAns[ii] = median(HL[ii,]);
  }
  HSAns <- MakeAns;
  HSAns[abs(HSAns) < MinMin] = 0;
  ActiveHS = HSAns;
  ActiveHS[ActiveHS != 0] = 1;
	     
	T2 <- proc.time();
	
	if (DoCI == TRUE) {
    library(coda);
    BetaChains <- as.mcmc(t(GibbHS$RecordBetajs));
    CITime1 = proc.time();
    eval(parse(text=GetG0Text("GetConfidenceIntervals")));
    ASG <- summary(BetaChains, quantiles = GetConfidenceIntervals);
    BetaSymmetricIntervals <- ASG[[2]];
    colnames(BetaSymmetricIntervals) <- GetConfidenceIntervals;
    rownames(BetaSymmetricIntervals) <- paste("Beta", 1:NROW(BetaSymmetricIntervals), sep="");
    BetaHPDIntervals <- BetaSymmetricIntervals * 0;
    if (GetConfidenceIntervals[1] == .5) {
      AK <- 1;  BetaHPDIntervals[,1] <- BetaSymmetricIntervals[,1];
    } else {AK <- 0;}
    NIntervals <- (length(GetConfidenceIntervals) - AK)/2;  AtAk <- AK;
    for (jj in 1:NIntervals) {
      CIP <- GetConfidenceIntervals[AtAk +2] - GetConfidenceIntervals[AtAk +1];
      try(HI <- TwoSimR5:::StableHPDinterval(BetaChains, CIP), silent=TRUE);
      try(BetaHPDIntervals[,AtAk+1:2] <- HI);
      AtAk <- AtAk +2;
    }
    CIEst <- list(BetaSymmetricIntervals=BetaSymmetricIntervals,
      BetaHPDIntervals=BetaHPDIntervals);
    CIQuantiles =  GetConfidenceIntervals;
    try(CIEst <- CleanCIEst(CIEst));
    CITime2 = proc.time();
    CITime = CITime2-CITime1;
  } else {CIEst = NULL; CIQuantiles = NULL;  CITime = NULL; }
	HSFitTime <- T2 - T1;
	##AFilePrint("Finished RunningENetFixed"); flush.console();
	return(list( type="HorseShoe", BetaFit = HSAns, BBHatsFit = ActiveHS,
	  FitTime =HSFitTime, OtherBetas = MakeAns, HitOb = HL,
    CIEst=CIEst, CIQuantiles=CIQuantiles, CITime=CITime, MIPReport = ActiveHS) );                  
}

################################################################################
### Generate HorseShoeOld
###
###   Gives HorseShoe Median estimate
###
MinMin  = .005;
GenerateHorseShoeOld <- function(SMS, tauSq = 1, 
  CountSamp = DefaultHorseShoeSamples, DoCI = DefaultDoCI,...) {
  AFilePrint("GenerateHorseShoe starting ");
  kLen = SMS$kLen;
  library(TwoLassoCpp);
  AFilePrint("Starting SampleGibbsHorseShoe ");
  ##AFilePrint("ENetFixed Run"); flush.console();
  T1 <- proc.time();
  GibbHS<- TwoLassoCpp:::SampleGibbsHorseShoe(X=SMS$X, Y=SMS$Y, 
    tauSq, LengthSamples=CountSamp, SigmaSq=SMS$SigmaSqNoise, printFlag=10);
	HL <- GibbHS$RecordBetajs[,round(CountSamp/4):CountSamp];
	MakeAns = HL[,1];
  for (ii in 1:length(HL[,1])) {
    MakeAns[ii] = median(HL[ii,]);
  }
  HSAns <- MakeAns;
  HSAns[abs(HSAns) < MinMin] = 0;
  ActiveHS = HSAns;
  ActiveHS[ActiveHS != 0] = 1;
	     
	T2 <- proc.time();
	
	if (DoCI == TRUE) {
    library(coda);
    BetaChains <- as.mcmc(t(GibbHS$RecordBetajs));
    CITime1 = proc.time();
    eval(parse(text=GetG0Text("GetConfidenceIntervals")));
    ASG <- summary(BetaChains, quantiles = GetConfidenceIntervals);
    BetaSymmetricIntervals <- ASG[[2]];
    colnames(BetaSymmetricIntervals) <- GetConfidenceIntervals;
    rownames(BetaSymmetricIntervals) <- paste("Beta", 1:NROW(BetaSymmetricIntervals), sep="");
    BetaHPDIntervals <- BetaSymmetricIntervals * 0;
    if (GetConfidenceIntervals[1] == .5) {
      AK <- 1;  BetaHPDIntervals[,1] <- BetaSymmetricIntervals[,1];
    } else {AK <- 0;}
    NIntervals <- (length(GetConfidenceIntervals) - AK)/2;  AtAk <- AK;
    for (jj in 1:NIntervals) {
      CIP <- GetConfidenceIntervals[AtAk +2] - GetConfidenceIntervals[AtAk +1];
      try(HI <- TwoSimR5:::StableHPDinterval(BetaChains, CIP), silent=TRUE);
      try(BetaHPDIntervals[,AtAk+1:2] <- HI);
      AtAk <- AtAk +2;
    }
    CIEst <- list(BetaSymmetricIntervals=BetaSymmetricIntervals,
      BetaHPDIntervals=BetaHPDIntervals);
    CIQuantiles =  GetConfidenceIntervals;
    try(CIEst <- CleanCIEst(CIEst));
    CITime2 = proc.time();
    CITime = CITime2-CITime1;
  } else {CIEst = NULL; CIQuantiles = NULL;  CITime = NULL; }
	HSFitTime <- T2 - T1;
	##AFilePrint("Finished RunningENetFixed"); flush.console();
	return(list( type="HorseShoe", BetaFit = HSAns, BBHatsFit = ActiveHS,
	  FitTime =HSFitTime, OtherBetas = MakeAns, HitOb = HL,
    CIEst=CIEst, CIQuantiles=CIQuantiles, CITime=CITime, MIPReport = ActiveHS) );                  
}

################################################################################
### Generate SpikeAndSlab
###
###   Gives Manaual "Spike and Slab" Using our R code (Old, we probably should do it)
###
MinMin  = .005;
GenerateSpikeAndSlab <- function(SMS, tauSqA = 1, CountSamp = 2000, DoCI = DefaultDoCI,
  ...) {
  kLen = SMS$kLen;
  ##AFilePrint("ENetFixed Run"); flush.console();
  T1 <- proc.time();
  GibbsSS <- SampleGibbsSpikeSlab(SMS$X, SMS$Y, 
    tauASq = tauSqA, piA=SMS$puse, LengthSamples = CountSamp,
    SigmaSq = SMS$SigmaSqNoise, printFlag = 0)
	HL <- GibbsSS$RecordBetajs[,round(CountSamp/4):CountSamp];
	MakeAns = HL[,1];
  for (ii in 1:length(HL[,1])) {
    MakeAns[ii] = median(HL[ii,]);
  }
  HSAns <- MakeAns;
  HSAns[abs(HSAns) < MinMin] = 0;
  ActiveHS = HSAns;
  ActiveHS[ActiveHS != 0] = 1;
  T2 <- proc.time();
	HSFitTime <- T2 - T1;
  library(coda);
 ##AS <- summary(as.mcmc(t(GibbsSS$RecordBjs)));
  MIPReport <- rowMeans(GibbsSS$RecordBjs);
  CIQuantiles = NULL;  CITime = NULL;  CIEst = NULL;
  if (is.logical(DoCI) && DoCI == TRUE) {
     eval(parse(text=GetG0Text("GetConfidenceIntervals", "globalenv()", S=1)));
     if (length(GetConfidenceIntervals) > 0 && GetConfidenceIntervals[1] != 0) {
       CITime1 <- proc.time();
       LongM <- as.mcmc(t(GibbsSS$RecordBetajs));
       CIQuantiles = GetConfidenceIntervals;
       SymmetricMatrix <- matrix(0, NCOL(LongM), length(CIQuantiles));
       UnshrunkSymmetricMatrix <- SymmetricMatrix * 0.0;
       for (jj in 1:NCOL(LongM)) {
         try(SymmetricMatrix[jj,] <- quantile(LongM[,jj], CIQuantiles));
         ATT <- NULL;
         try(ATT <- LongM[LongM[,jj] != 0.0,jj]);
         if (length(ATT) >= 1) {
           try(UnshrunkSymmetricMatrix[jj,] <- quantile(ATT, CIQuantiles));
         }
       }
       colnames(SymmetricMatrix) <- CIQuantiles;
       colnames(UnshrunkSymmetricMatrix) <- CIQuantiles;
       try(rownames(SymmetricMatrix) <- paste("Beta:", 1:NCOL(LongM), sep=""));
       try(rownames(UnshrunkSymmetricMatrix) <- paste("Beta:", 1:NCOL(LongM), sep=""));
       HPDMatrix <- SymmetricMatrix*0;
       HPDUnshrunk <- UnshrunkSymmetricMatrix*0;        IO = 0;
       if (CIQuantiles[1] == .5) {
         HPDMatrix[,1] <- SymmetricMatrix[,1];
         HPDUnshrunk[,1] <- UnshrunkSymmetricMatrix[,1]; IO <- IO +1;
       }
       while(IO +2 < length(CIQuantiles)) {
         Down <- IO +1;  Up <- IO+2;
         AQ <- CIQuantiles[Up]-CIQuantiles[Down];
         try(APP <- TwoSimR5:::StableHPDinterval(LongM, AQ), silent=TRUE);
         HPDMatrix[,c(Down,Up)] <- APP;
         for (jj in 1:NCOL(LongM)) {
           ATT <- NULL;
           try(ATT <- LongM[LongM[,jj] != 0.0 & !is.na(LongM[,jj]),jj])
           if (length(ATT) >= 1) {
             try(AT <- TwoSimR5:::StableHPDinterval(as.mcmc(ATT), AQ), silent=TRUE)
             try(HPDUnshrunk[jj,c(Down,Up)]  <- as.numeric(AT));
           }
         }
         IO <- IO + 2;
       }
       CIEst <- list(BetaSymmetricQuantiles = SymmetricMatrix, BetaSymmetricUnshrinkQuantiles =
         UnshrunkSymmetricMatrix, BetaHPDQuantiles = HPDMatrix, BetaHPDUnshrinkQuantiles = HPDUnshrunk);
       CITime2 <- proc.time();
       CITime <- CITime2[1:3]-CITime1[1:3];
     } else {CIEst <- NULL; CITime <- NULL; CIQuantiles=NULL; }
  }
  
	##AFilePrint("Finished RunningENetFixed"); flush.console();
	return(list( type="SpikeAndSlab", BetaFit = HSAns, BBHatsFit = ActiveHS,
	  FitTime =HSFitTime, OtherBetas = MakeAns, HitOb = HL, MIPReport = MIPReport,
    CIQuantiles=CIQuantiles, CITime=CITime, CIEst=CIEst) );                  
}
IntroGenerateBayesSpikeText <- function() {
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
  if (!exists(\"AutoInfo\")) { AutoInfo = \"Auto\"; }
  if (!exists(\"DoLogitPostPreProb\")) { DoLogitPostPreProb = 1; }
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE)  {
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
  if (!is.null(TemperatureList)) {
  } else if (NumTemp < 2) {
    TemperatureList = NULL; EEMergeEvery = NULL;
  }  else {
     NumTemp <- round(NumTemp);
     TemperatureList <- (2^(1:(NumTemp-1)));
  }
  if (TestAtGoal == TRUE) {
    if (!is.null(SMS$BetasPartReal) && length(SMS$BetasPartReal) == SMS$p) {
       BetaStart <- SMS$BetasPartReal;
    }  else {
      print(paste(\"Hey: TestAtGoal is TRUE: I'm thinking you'd rather that SMS$BetasPartReal was set up! = \", SMS$p, sep=\"\"));
      flush.console();
        BetaStart <- SMS$BetasReal;
    }
    NoNoiseBetaStart <- TRUE;
  }  else {
    BetaStart <- NULL; NoNoiseBetaStart <- FALSE;
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
  } else if (AutoInfo == \"AutoSqrt\") {
    if (kLen >= 1100 && PriorStrength > 0) {
      PiAPrior = c( PriorStrength * 1/sqrt(SMS$kLen) + 1, sqrt(PriorStrength) +
        1);
      SigmaPrior = c(SigmaPriorStrength, .2*var(SMS$Y));
    } else if (SMS$kLen >= 1000 && PriorStrength <= 0) {
      PiAPrior = c(2, sqrt(SMS$kLen)+1);
      SigmaPrior = c(SMS$n, .2*var(SMS$Y));
    } else if (PriorStrength > 0) {
      PiAPrior = c(PriorStrength * 1 / sqrt(SMS$kLen)+1, PriorStrength +1);
      SigmaPrior = c(SigmaPriorStrength, SMS$SigmaSqNoise);
    } else {
      SigmaPrior = c(-1,-1); PiAPrior = c(-1,-1);
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
################################################################################
### GenerateBayesSpikeErrorPrior
###
###   Gives BayesSpikeEstimate Median estimate
###
MinMin  = .005;
GenerateBayesSpike <- function(SMS, tauSqA = 1, CountSamp = 1000, 
  NumChains = 3, PriorStrength = 0,  burnin = 25,
  CutOff = .1, DoCI = DefaultDoCI, DoMedian = TRUE, 
  SigmaPriorStrength = 1, AutoInfo="Auto", DoLogitPostPreProb = 1,
  NumTemp = 1, TestAtGoal = FALSE, TemperatureList = NULL, EEMergeEvery = 10, RevertTemperatureEvery = -1, ...) {
  library(BayesSpike);
  kLen = SMS$kLen;
  ##AFilePrint("ENetFixed Run"); flush.console();
  T1 <- proc.time();
  if (!exists("NumChains")) { NumChains = 3; }
  if (!exists("CountSamp")) { CountSamp = 1000; }
  if (!exists("PriorStrength")) { PriorStrength = 0; }
  if (!exists("tauSqA")) { tauSqA = 1; }
  if (!exists("CutOff")) { CutOff = .1; }
  if (!exists("AutoInfo")) { AutoInfo = "Auto"; }
  if (!exists("DoLogitPostPreProb")) { DoLogitPostPreProb = 1; }
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE)  {
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
    AutoInfo = "Auto"
  } else if (is.logical(AutoInfo) && AutoInfo == FALSE) {
    AutoInfo = "Info"
  } else if (is.character(AutoInfo) && !(AutoInfo %in% c("Auto", "Info"))) {
    print(paste("Error AUTOINFO NOT SET!", sep=""));
  } else if (is.character(AutoInfo) && AutoInfo %in% c("Auto", "Info")) {
  } else {
    AutoInfo = "Auto";
  }
  if (!is.null(TemperatureList) && length(TemperatureList) >= 2) {
    TemperatureList <- sort(TemperatureList, decreasing=TRUE);
    NumTemp <- length(TemperatureList);
  } else if (NumTemp < 2) {
    TemperatureList = NULL; EEMergeEvery = NULL;
  }  else if (NumTemp >= 2) {
    NumTemp <- round(NumTemp);
    TemperatureList <- (2^((1:(NumTemp-1))/NumTemp));
  }  else {
    TemperatureList = NULL; EEMergeEvery = NULL;
  }
  if (TestAtGoal == TRUE) {
    if (!is.null(SMS$BetasPartReal) && length(SMS$BetasPartReal) == SMS$p) {
       BetaStart <- SMS$BetasPartReal;
    }  else {
      print(paste("Hey: TestAtGoal is TRUE: I'm thinking you'd rather that SMS$BetasPartReal was set up! = ", SMS$p, sep=""));
      flush.console();
        BetaStart <- SMS$BetasReal;
    }
    NoNoiseBetaStart <- TRUE;
  }  else {
    BetaStart <- NULL; NoNoiseBetaStart <- FALSE;
  }
  if (AutoInfo == "Info") {  
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
  } else if (AutoInfo == "AutoSqrt") {
    if (kLen >= 1100 && PriorStrength > 0) {
      PiAPrior = c( PriorStrength * 1/sqrt(SMS$kLen) + 1, sqrt(PriorStrength) +
        1);
      SigmaPrior = c(SigmaPriorStrength, .2*var(SMS$Y));
    } else if (SMS$kLen >= 1000 && PriorStrength <= 0) {
      PiAPrior = c(2, sqrt(SMS$kLen)+1);
      SigmaPrior = c(SMS$n, .2*var(SMS$Y));
    } else if (PriorStrength > 0) {
      PiAPrior = c(PriorStrength * 1 / sqrt(SMS$kLen)+1, PriorStrength +1);
      SigmaPrior = c(SigmaPriorStrength, SMS$SigmaSqNoise);
    } else {
      SigmaPrior = c(-1,-1); PiAPrior = c(-1,-1);
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
  eval(parse(text=GetG0Text("BSSaveDir")));
  eval(parse(text=GetG0Text("UniqueProcessIdentifier")));
  if (is.character(BSSaveDir) && (nchar(BSSaveDir) < nchar(UniqueProcessIdentifier) ||
    substr(BSSaveDir, nchar(BSSaveDir) - nchar(UniqueProcessIdentifier)+1,
    nchar(BSSaveDir)) != UniqueProcessIdentifier)) {
    BSSaveDir <- paste(BSSaveDir, "//", UniqueProcessIdentifier, sep="");      
  }
  if (kLen <= 250 || is.null(BSSaveDir) || is.numeric(BSSaveDir) || !is.character(BSSaveDir) ||
    BSSaveDir %in% c("NoSave", "NOSAVE", "NoSAVE", "nosave")) {
    GibbsSS <- BayesSpike:::BayesSpikeRegression(X=SMS$X, Y=SMS$Y, 
      tauASq = tauSqA, PiAStart=SMS$puse, MaxGibbsIters = CountSamp,
      NumChains=NumChains,
      SigmaSq = SMS$SigmaSqNoise, Verbose = 1, tauEndList = NULL, PiAPrior = PiAPrior,
      SigmaPrior = SigmaPrior, DoSave=FALSE,
      DoRecord=c(1,1,0,0,0,1,1), IndexFirstRandomEffect = -1, StartRunProbVector = 5,
      DoLogitPreProb = DoLogitPreProb, DoLogitPostProb = DoLogitPostProb,
      dfRobit = dfRobit, PreRunMaxGibbsIters = PreRunMaxGibbsIters, 
      AlterWeightFlag = AlterWeightFlag, burnin = burnin, TemperatureList=TemperatureList,
      EEMergeEvery = EEMergeEvery, RevertTemperatureEvery = RevertTemperatureEvery,
      BetaStart=BetaStart);
     T2 <- proc.time();
     ##Keep = GibbsSS$PctBetas;
     MIPReport = GibbsSS$MIP;
     if (is.null(MIPReport)) {
       print(paste("BayesSpikeErrorPrior: we got null MIP for MIPReport kLen = ", kLen, sep=""));
       flush.console();
       MIPReport <- rep(.5, kLen);
     }
     SubSetCoords <- (1:length(MIPReport))[MIPReport >= .5];
     GibbsSS$SubSetCoords <- SubSetCoords;
     if (length(SubSetCoords) >= 1 && DoLogitPostProb == TRUE) {
       LSCL <- GibbsSS$TBSR5$LoadSubCodaList(StartIter = 1, EndIter = GibbsSS$TBSR5$MaxGibbsIters);
       ALongBetaMatrix <- AStackCoda(GibbsSS$SubCodaList, GibbsSS$burnin);
        try(WeightCodas <- GibbsSS$TBSR5$ABayesSpikeCL$AlterWeightCodaList)
        try(StackWeightCodas <- AStackCoda(WeightCodas, GibbsSS$burnin));
        if (sd(as.vector(StackWeightCodas)) >= 1.5) {
          StackWeightCodas[StackWeightCodas > quantile(StackWeightCodas,.99)] <- quantile(StackWeightCodas, .99);
         StackWeightCodas <- StackWeightCodas / sd(as.vector(StackWeightCodas)) * 1.5;
        }
        ReWeight <- exp(StackWeightCodas-max(StackWeightCodas));
        ReWeight <- ReWeight/ sum(ReWeight);
        MyQuant <- QuantileMe(ReWeight, ALongBetaMatrix, QuantileLocs = .5, Verbose = 0);
       HSAns <- rep(0, GibbsSS$p);
       HSAns[SubSetCoords] <- MyQuant;
       MakeAns <- rep(0, GibbsSS$p);
       if (length(SubSetCoords) == 1) {
         MakeAns[SubSetCoords] <- sum(as.vector(ReWeight) * ALongBetaMatrix);
       } else {
         MakeAns[SubSetCoords] <- colSums(as.vector(ReWeight) * ALongBetaMatrix); 
       }
     }  else {
     	HL <- GibbsSS$CodaList
     	NHL <- list();
	    for (ii in 1:length(HL)) {
        NHL[[ii]] <- as.mcmc(HL[[ii]][round(NROW(HL[[ii]])/4):NROW(HL[[ii]]),]);
      }
     try(NHL <- as.mcmc.list(NHL));
     MySum <- summary(NHL);
     MakeAns <- MySum[[2]][1:GibbsSS$p,3];
     KeepIters <- (1:NCOL(SMS$X))[MySum[[2]][1:GibbsSS$p,3] != 0.0];
     HSAns <- MySum[[2]][1:GibbsSS$p,3]; 
     MakeAns <- MySum[[1]][1:GibbsSS$p,1];
     }
  } else {
    GibbsSS <- BayesSpike:::BayesSpikeRegression(X=SMS$X, Y=SMS$Y, 
      tauASq = tauSqA, PiAStart=SMS$puse,InitKKs=6, MaxGibbsIters = CountSamp,
      NumChains=NumChains,
      SigmaSq = SMS$SigmaSqNoise, Verbose = 0, tauEndList = NULL, PiAPrior = PiAPrior,
      SigmaPrior = SigmaPrior, DoSave=TRUE, DoSaveTBSR5=FALSE, DoRecord = c(0,0,0,0,0,0,0),
      SaveDir = BSSaveDir, IndexFirstRandomEffect = -1,
      StartRunProbVector = 5,
      DoLogitPreProb = DoLogitPreProb, DoLogitPostProb = DoLogitPostProb,
      dfRobit = dfRobit, PreRunMaxGibbsIters = PreRunMaxGibbsIters, 
      AlterWeightFlag = AlterWeightFlag, burnin=burnin, TemperatureList=TemperatureList,
      EEMergeEvery = EEMergeEvery, RevertTemperatureEvery = RevertTemperatureEvery,
      BetaStart=BetaStart); 
     T2 <- proc.time();
     HL <- NULL;
     Keep = GibbsSS$MIP;  ## PctBetas
     MIPReport = Keep;
     if (is.null(MIPReport)) {
       print("BayesSpikeError Prior: we didn't get MIP out of GibbsSS!");
       Keep <- rep(.5, NCOL(SMS$X));
       MIPReport <- Keep;
     }
     SubSetCoords <- (1:length(Keep))[Keep >= .5];
     GibbsSS$SubSetCoords <- SubSetCoords;

     if (length(SubSetCoords) >= 1 && DoLogitPostProb == TRUE) {
       LSCL <- GibbsSS$TBSR5$LoadSubCodaList(StartIter = 1, EndIter = GibbsSS$TBSR5$MaxGibbsIters);
       ALongBetaMatrix <- AStackCoda(GibbsSS$SubCodaList, GibbsSS$burnin);
        try(WeightCodas <- GibbsSS$AlterWeightCodaList)
        try(StackWeightCodas <- AStackCoda(WeightCodas, GibbsSS$burnin));
        if (sd(as.vector(StackWeightCodas)) >= 1.5) {
          StackWeightCodas[StackWeightCodas > quantile(StackWeightCodas,.99)] <- quantile(StackWeightCodas, .99);
         StackWeightCodas <- StackWeightCodas / sd(as.vector(StackWeightCodas)) * 1.5;
        }
        ReWeight <- exp(StackWeightCodas-max(StackWeightCodas));
        ReWeight <- ReWeight/ sum(ReWeight);
        MyQuant <- QuantileMe(ReWeight, ALongBetaMatrix, QuantileLocs = .5, Verbose = 0);
       HSAns <- rep(0, GibbsSS$p);
       HSAns[SubSetCoords] <- MyQuant;
       MakeAns <- rep(0, GibbsSS$p);
       if (length(SubSetCoords) == 1) {
         MakeAns[SubSetCoords] <- sum(as.vector(ReWeight) * ALongBetaMatrix);
       } else {
         MakeAns[SubSetCoords] <- colSums(as.vector(ReWeight) * ALongBetaMatrix); 
       }
     }  else {
       ACoda2 <- GibbsSS$SubCodaList;
       MySumP <- summary(ACoda2);
       KeepIters <- SubSetCoords[MySumP[[2]][,3] != 0.0];
       HSAns <- rep(0, GibbsSS$p);
       HSAns[SubSetCoords] <- MySumP[[2]][,3];
       MakeAns <- rep(0, GibbsSS$p);
       MakeAns[SubSetCoords] <- MySumP[[1]][,1];
     }
  }
   

  if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
     MakeAns <- MakeAns[2:length(MakeAns)];
     HSAns <- HSAns[2:length(HSAns)];
     MIPReport <- MIPReport[2:length(MIPReport)];
  }
  GibbsSS$TBSR5$unsave(TRUE);
  
  ActiveHS = HSAns;
  ActiveHS[ActiveHS != 0] = 1;

	HSFitTime <- T2 - T1;
    
    OtherHitOb <- NULL;	
	if (DoCI == TRUE) {
	  eval(parse(text=GetG0Text("GetConfidenceIntervals")));

	  CITime1 <- proc.time();
      NGibbsSS <- BayesSpike:::BayesSpikeRegression(X=SMS$X, Y=SMS$Y, 
      tauASq = tauSqA, PiAStart=SMS$puse, MaxGibbsIters = CountSamp,
      NumChains=NumChains, 
      SigmaSq = SMS$SigmaSqNoise, Verbose = 0, tauEndList = NULL, PiAPrior = PiAPrior,
      SigmaPrior = SigmaPrior, DoSave=TRUE, DoSaveTBSR5=FALSE, DoRecord = c(0,0,0,0,0,0,0),
      SaveDir = BSSaveDir, IndexFirstRandomEffect = -1, DoLongCI = TRUE, StartRunProbVector=5,
      DoLogitPreProb = DoLogitPreProb, DoLogitPostProb = DoLogitPostProb,
      dfRobit = dfRobit, PreRunMaxGibbsIters = PreRunMaxGibbsIters, 
      AlterWeightFlag = AlterWeightFlag, burnin = burnin, TemperatureList=TemperatureList,
      EEMergeEvery = EEMergeEvery,  RevertTemperatureEvery = RevertTemperatureEvery,
      BetaStart=BetaStart); 
      NGibbsSS$TBSR5$KeepPosteriorQuantiles <- GetConfidenceIntervals;
      try(NGibbsSS$burnin <- 25);
	  CIEst <- list();
	   try(CIEst[[1]] <- NGibbsSS$TBSR5$BetaSymmetricQuantiles);
	   try(CIEst[[2]] <- NGibbsSS$TBSR5$BetaSymmetricUnshrinkQuantiles);
	   try(CIEst[[3]] <- NGibbsSS$TBSR5$BetaHPDQuantiles);
	   try(CIEst[[4]] <- NGibbsSS$TBSR5$BetaHPDUnshrinkQuantiles);  
	   try(CIEst[[5]] <- NGibbsSS$TBSR5$BetaSymmetricLongQuantiles);
	   try(CIEst[[6]] <- NGibbsSS$TBSR5$BetaHPDLongQuantiles); 
    try(CIEst <- CleanCIEst(CIEst))   	   
    CITime2 <- proc.time();
      names(CIEst) <- c("BetaSymmetricQuantiles", "BetaSymmetricUnshrinkQuantiles",
        "BetaHPDQuantiles", "BetaHPDUnshrinkQuantiles", 
        "BetaSymmetricLongQuantiles", "BetaHPDLongQuantiles");
      CIQuantiles <- NGibbsSS$TBSR5$KeepPosteriorQuantiles;
    CITime <- CITime2-CITime1;
    OtherHitOb <- GibbsSS;  GibbsSS <- NGibbsSS;
    NGibbsSS$TBSR5$unsave(TRUE);
    if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
      for (ii in 1:length(CIEst)) {
        try(CIEst[[ii]] <- CIEst[[ii]][2:NROW(CIEst[[ii]]),]);
        try(rownames(CIEst[[ii]]) <- paste("Beta:", 1:NROW(CIEst[[ii]]), sep=""));
      }
    }
  } else {
      CIEst <- NULL;   CIQuantiles <- NULL;  CITime <- NULL;
      OtherHitOb <- NULL;
  }  
	##AFilePrint("Finished BayesSpike"); flush.console();
	return(list( type=paste("BayesSpike", AutoInfo), BetaFit = HSAns, BBHatsFit = ActiveHS,
	  FitTime =HSFitTime, OtherBetas = MakeAns, HitOb = GibbsSS,
      OtherHitOb = OtherHitOb, GibbsSS = GibbsSS, BackGibbsSS = OtherHitOb,
      CIEst=CIEst, CIQuantiles=CIQuantiles, CITime=CITime,MIPReport=MIPReport,
      BetaStart=BetaStart, TemperatureList=TemperatureList,
      DoMedian=DoMedian) );                  
}


################################################################################
### GenerateBayesSpikeArtPrior
###
###   Gives BayesSpike Learning Prior Median estimate
###
MinMin  = .005;
GenerateBayesSpikeNo <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  CutOff = .1, NumChains=3, DoCI = DefaultDoCI, DoMedian = TRUE, 
  NumTemp = 1, TestAtGoal = FALSE,  EEMergeEvery = 10, RevertTemperatureEvery = -1, ...) {
    
  if (!exists("PriorStrength")) { PriorStrength <- 0; }
  if (!exists("CutOff")) { CutOff <- .1; }
  if (!exists("NumChains")) { NumChains <- 3; }
  if (!exists("DoCI")) { DoCI <- DefaultDoCI; }
  if (!exists("CountSamp")) { CountSamp <- 2000; }
  if (!exists("tauSqA")) { tauSqA <- 1; }
  library(BayesSpike);
  kLen = SMS$kLen;
  ##AFilePrint("BayesSpikeArtPrior Run"); flush.console();
  if (NCOL(SMS$XX) == 1) {
    X <- SMS$XX - mean(SMS$XX);
  } else {
    X <- t(t(SMS$XX) - colMeans(SMS$XX));
  }
  Y <- SMS$YY - mean(SMS$YY);
  T1 <- proc.time();
  Prob = .08 + .42 * min(c(length(SMS$Y) / NCOL(SMS$X),1));

  if (NumTemp < 2) {
    TemperatureList = NULL;
  }  else {
     NumTemp <- round(NumTemp);
     TemperatureList <- (2^(1:(NumTemp-1)));
  }
  if (TestAtGoal == TRUE) {
    if (!is.null(SMS$BetasPartReal) && length(SMS$BetasPartReal) == SMS$p) {
       BetaStart <- SMS$BetasPartReal;
    }  else {
      print(paste("Hey: I'm sure you'd rather that SMS$BetasPartReal was set up! = ", SMS$p, sep=""));
      flush.console();
        BetaStart <- SMS$BetasReal;
    }
    NoNoiseBetaStart <- TRUE;
  }  else {
    BetaStart <- NULL; NoNoiseBetaStart <- FALSE;
  }
  if (kLen >= 1100 && PriorStrength > 0) {
    PiAPrior = c(PriorStrength * 1/SMS$kLen + 1, PriorStrength +
      1);
    SigmaPrior = c(PriorStrength, SMS$SigmaSqNoise);
  } else if (SMS$kLen >= 1000 && PriorStrength <= 0) {
    PiAPrior = c(2, SMS$kLen+1);
    SigmaPrior = c(SMS$n, SMS$SigmaSqNoise);
  } else if (PriorStrength > 0) {
    PiAPrior = c(PriorStrength * 1 / SMS$kLen+1, PriorStrength +1);
    SigmaPrior = c(PriorStrength, SMS$SigmaSqNoise);
  } else {
    SigmaPrior = c(-1,-1); PiAPrior = c(-1,-1);
  }
  eval(parse(text=GetG0Text("BSSaveDir")));
  eval(parse(text=GetG0Text("UniqueProcessIdentifier")));
  if (is.character(BSSaveDir) && (nchar(BSSaveDir) < nchar(UniqueProcessIdentifier) ||
    substr(BSSaveDir, nchar(BSSaveDir) - nchar(UniqueProcessIdentifier)+1,
    nchar(BSSaveDir)) != UniqueProcessIdentifier)) {
    BSSaveDir <- paste(BSSaveDir, "//", UniqueProcessIdentifier, sep="");      
  }
  if (kLen <= 250 || is.null(BSSaveDir) || is.numeric(BSSaveDir) || !is.character(BSSaveDir) ) {
    GibbsSS <- BayesSpike:::BayesSpikeRegression(X=X, Y=Y, 
      tauASq = tauSqA, PiAStart=SMS$puse, MaxGibbsIters = CountSamp,
      SigmaSq = SMS$SigmaSqNoise, Verbose = 0, tauEndList = NULL, PiAPrior = PiAPrior,
      SigmaPrior = SigmaPrior, DoSave=FALSE, SaveDir = NULL, 
      DoRecord=c(1,0,1,1,0,0,0), IndexFirstRandomEffect = -1, StartRunProbVector = 5,
      TemperatureList=TemperatureList, BetaStart=BetaStart, RevertTemperatureEvery = RevertTemperatureEvery, NoNoiseBetaStart=NoNoiseBetaStart);
     T2 <- proc.time();
    ## Keep = GibbsSS$PctBetas;
    MIPReport <- NULL;
    try(MIPReport <- GibbsSS$MIP);
    if (is.null(MIPReport)) {
      print("BayesSpikeAutoPrior: sorry MIPReport is null."); flush.console();
      MIPReport <- rep(.5, kLen);
    }
     HL <- GibbsSS$CodaList;
     NHL <- list();
	    for (ii in 1:length(HL)) {
        NHL[[ii]] <- as.mcmc(HL[[ii]][round(NROW(HL[[ii]])/4):NROW(HL[[ii]]),]);
      }
     try(NHL <- as.mcmc.list(NHL));
     MySum <- summary(NHL);
     MakeAns <- MySum[[2]][1:GibbsSS$p,3];
     KeepIters <- (1:NCOL(SMS$X))[MySum[[2]][1:GibbsSS$p,3] != 0.0];
     HSAns <- MySum[[2]][1:GibbsSS$p,3]; 
     MakeAns <- MySum[[1]][1:GibbsSS$p,1];
  } else {
    GibbsSS <- BayesSpike:::BayesSpikeRegression(X=X, Y=Y, 
      tauASq = tauSqA, PiAStart=SMS$puse, MaxGibbsIters = CountSamp,
      SigmaSq = SMS$SigmaSqNoise, Verbose = 0, InitKKs=6, tauEndList = NULL, PiAPrior = PiAPrior,
      SigmaPrior = SigmaPrior, DoSave=TRUE, DoSaveTBSR5=FALSE,
      SaveDir = BSSaveDir, IndexFirstRandomEffect = -1, StartRunProbVector = 5,
      TemperatureList=TemperatureList, RevertTemperatureEvery = RevertTemperatureEvery, BetaStart=BetaStart, NoNoiseBetaStart=NoNoiseBetaStart); 
     T2 <- proc.time();
     print("Done Running BayesSpikeRegresion for large p first time"); flush.console();
     MIPReport <- NULL;
     HL <- NULL;
     try(MIPReport <- GibbsSS$MIP);
     if (is.null(MIPReport)) {
       print("---------------------------------------------------------------");
       print("---BayesSpikeAutoPrior: Big Bad issue why is MIPReport NULLL?")
       print("--- BayesSpikeAutoPrior: with kLen big, MIPReport is null."); flush.console();
       print("--- Why did this fail?");
       MIPReport <- rep(.5, kLen);
       Keep <- rep(.5, kLen);
       eval(parse(text=SetGText("GibbsSS", "globalenv()", S=1)));
     } else {
      print("--- BayesSpikeAutoPrior: Supposedly good MIP");
      if (all(abs(MIPReport - .5) < .001)) {
        print("--- BayesSpikeAutoPrior: Issue, all of MIP is .5 for some reason?"); flush.console();
      }
       Keep = GibbsSS$MIP;     ##PctBetas;
       MIPReport = Keep;
     }
     eval(parse(text=SetGText("Keep", "globalenv()", S=1)));
     eval(parse(text=SetGText("MIPReport", "globalenv()", S=1)));
     SubSetCoords <- (1:length(Keep))[Keep >= .5];
     if (length(SubSetCoords) <= 0 && TestAtGoal == TRUE) {
      ## Test accuracy when TestAtGoal is TRUE;
      SubSetCoords <- sort(Keep, decreasing=TRUE, index=TRUE)$ix[1:12];
     }
     GibbsSS$SubSetCoords <- SubSetCoords;
     eval(parse(text=SetGText("SubSetCoords", S=1)));
     if (length(SubSetCoords) <= 0) {
       print("Hey: After all that SubSetCoords is still length zero!");
     }
     ACoda2 <- GibbsSS$SubCodaList;
     if (is.null(ACoda2) || length(ACoda2) <= 0) {
       print("Hey: ACoda2 was returned as a null, the rest of this won't work!"); flush.console();
     }
     MySumP <- summary(ACoda2);
     if (length(MySumP) <= 1) {
       print("Hey: MySumP does not have 2 dimensions!");
       eval(parse(text=SetGText("MySumP",S=1)));
       eval(parse(text=SetGText("ACoda2",S=1)));
     }
     try(KeepIters <- SubSetCoords[MySumP[[2]][,3] != 0.0]);
     HSAns <- rep(0, GibbsSS$p);
     try(HSAns[SubSetCoords] <- MySumP[[2]][,3]);
     MakeAns <- rep(0, GibbsSS$p);
     try(MakeAns[SubSetCoords] <- MySumP[[1]][,1]);
  }
 
  
	##HL <- GibbsSS$RecordBetajs[,round(CountSamp/4):CountSamp];
  HSAns[abs(HSAns) < MinMin] = 0;
  ActiveHS = HSAns;
  ActiveHS[ActiveHS != 0] = 1;

	HSFitTime <- T2 - T1;
	if (DoCI == TRUE) {
	  eval(parse(text=GetG0Text("GetConfidenceIntervals")));

	  CITime1 <- proc.time();
      NGibbsSS <- BayesSpike:::BayesSpikeRegression(X=X, Y=Y, 
        tauASq = tauSqA, PiAStart=SMS$puse, MaxGibbsIters = CountSamp,
        SigmaSq = SMS$SigmaSqNoise, Verbose = 0, tauEndList = NULL, PiAPrior = PiAPrior,
        SigmaPrior = SigmaPrior, DoSave=TRUE, DoSaveTBSR5=FALSE,
        SaveDir = BSSaveDir, IndexFirstRandomEffect = -1, DoLongCI = TRUE,
        TemperatureList=TemperatureList, BetaStart=BetaStart, NoNoiseBetaStart); 
      NGibbsSS$TBSR5$KeepPosteriorQuantiles <- GetConfidenceIntervals;
	  CIEst <- list();
	   try(CIEst[[1]] <- NGibbsSS$TBSR5$BetaSymmetricQuantiles);
	   if (SMS$p >= 10000) {
         print("Got BetaSymmetric Quantiles"); flush.console();
       }
       try(CIEst[[2]] <- NGibbsSS$TBSR5$BetaSymmetricUnshrinkQuantiles);
	   if (SMS$p >= 10000) {
         print("Got Beta UnshrinkSymmetric Quantiles"); flush.console();
       }
	   try(CIEst[[3]] <- NGibbsSS$TBSR5$BetaHPDQuantiles);
       if (SMS$p >= 10000) {
         print("Got Beta HPD Quantiles"); flush.console();
       }
	   try(CIEst[[4]] <- NGibbsSS$TBSR5$BetaHPDUnshrinkQuantiles);    
	   	if (SMS$p >= 10000) {
         print("Got Beta Unshrink Quantiles"); flush.console();
       }
       try(CIEst[[5]] <- NGibbsSS$TBSR5$BetaSymmetricLongQuantiles);
	   try(CIEst[[6]] <- NGibbsSS$TBSR5$BetaHPDLongQuantiles);  	   
	   if (SMS$p >= 10000) {
         print("Got Beta Long Quantiles"); flush.console();
       }
       CITime2 <- proc.time();
       names(CIEst) <- c("BetaSymmetricQuantiles", "BetaSymmetricUnshrinkQuantiles",
         "BetaHPDQuantiles", "BetaHPDUnshrinkQuantiles",
         "BetaSymmetricLongQuantiles",
         "BetaHPDLongQuantiles");
    CIQuantiles <- NGibbsSS$TBSR5$KeepPosteriorQuantiles;
    CITime <- CITime2-CITime1;
    try(CIEst <- CleanCIEst(CIEst));
    BackGibbsSS <- GibbsSS;  GibbsSS <- NGibbsSS;
  } else {
      CIEst <- NULL;   CIQuantiles <- NULL;  CITime <- NULL;
      BackGibbsSS <- NULL;
  }  
	##AFilePrint("Finished RunningENetFixed"); flush.console();
	return(list( type="BayesSpike", BetaFit = HSAns, BBHatsFit = ActiveHS,
	  FitTime =HSFitTime, OtherBetas = MakeAns, HitOb = GibbsSS,
      OtherHitOb <- BackGibbsSS,  GibbsSS = GibbsSS, BackGibbsSS=BackGibbsSS,
      CIEst=CIEst, CIQuantiles=CIQuantiles, CITime = CITime, MIPReport = MIPReport, Keep=Keep,
      DoMedian=DoMedian, TemperatureList=TemperatureList, BetaStart=BetaStart, NoNoiseBetaStart=NoNoiseBetaStart) );                  
}







MinMin  = .005;
GenerateISpike <- function(SMS, tauSqA = 1, CountSamp = 2000, PriorStrength = 0,
  CutOff = .1) {
  library(spikeslab, warn.conflicts=FALSE, quietly=TRUE)
  kLen = SMS$kLen;
  ##library(spikeslab)
  ##AFilePrint("ENetFixed Run"); flush.console();
  T1 <- proc.time();
  
  if (dim(SMS$X)[2] > dim(SMS$X)[1]) {
    bigp.smalln = TRUE;
  } else { bigp.smalln = FALSE; }
  
  MyS <- spikeslab(SMS$Y~SMS$X, data = NULL, x = SMS$X, y = SMS$Y, n.iter1 = 500, 
    n.iter2 = 500, mse = TRUE, bigp.smalln = bigp.smalln, bigp.smalln.factor = 1, 
    screen = bigp.smalln, r.effects = NULL, max.var = 500, 
    center = FALSE, intercept = FALSE, fast = TRUE, beta.blocks = 5, 
    verbose = FALSE, ntree = 300, seed = NULL) 

  HL = MyS$bma;
  HL2 = HL;  HL2[abs(HL) < MinMin] = 0;
  HL2[HL2 != 0] = 1;
  
  T2 <- proc.time();
	HSFitTime <- T2 - T1;
	##AFilePrint("Finished RunningENetFixed"); flush.console();
	return(list( type="SpikeAndSlab", BetaFit = HL, BBHatsFit = HL2,
	  FitTime =HSFitTime, OtherBetas = HL, HitOb = MyS, MIPReport = HL2) );                  
}



################################################################################
### Generate ENetFixed
###
###   Gives ENet the right parameters to choose optimal shrinkage lambda2
###
MinMin  = .005;
GenerateENetCVMin <- function(SMS, DoCI = DefaultDoCI, ...) {
  ##AFilePrint("RunningGenerateENetCVMin");      flush.console();
  library(elasticnet, warn.conflicts = FALSE, quietly=TRUE);
  kLen = SMS$kLen;
  T1 <- proc.time();
  if (is.null(SMS$LogitRegression)  || SMS$LogitRegression ==FALSE) {
	library(elasticnet, warn.conflicts = FALSE, quietly=TRUE);
	VarData <- var(SMS$BetasReal[SMS$BetasReal != 0]) / var(SMS$YY);
	## lambda2 = (1 / VarData) / sqrt(length(SMS$YY));
	lambda2 = 1/VarData;
	lambda2Vals <- c(lambda2 * .01, lambda2 * .1, lambda2, lambda2 * 10, lambda2 * 100, lambda2 * 1000, 1);
	MDD <- list();
	for (ii in 1:length(lambda2Vals)) {
    MDD[[ii]] = cv.enet(SMS$XX, SMS$YY, K=10, lambda=lambda2Vals[ii], mode="step", 
      s=1:min(c(SMS$kActiveLen *2, SMS$kLen-1, length(SMS$YY)-1)), normalize=TRUE, 
      intercept=FALSE, trace = FALSE, eps = .00001);
  }
	Minerizer <- rep(0, length(lambda2Vals));
	for (ii in 1:length(lambda2Vals)) {
    Minerizer[ii] = min(MDD[[ii]]$cv);
  }
  MyWant <- sort(Minerizer, index=TRUE)$ix[1];
  UseAns <- sort(MDD[[MyWant]]$cv, index=TRUE)$ix[1];
  HitENET2Ob <-  enet(SMS$XX, SMS$YY, lambda = lambda2Vals[MyWant],
    max.steps= max(MDD[[MyWant]]$s[UseAns],2), normalize=TRUE, 
    intercept=FALSE, trace = FALSE, eps = .00001);
    HitENET2Ob$beta = HitENET2Ob$beta.pure; 
  } else if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
	  library(glmnet, warn.conflicts = FALSE, quietly=TRUE);
	  ##VarData <- var(SMS$BetasReal[SMS$BetasReal != 0]) / var(SMS$YY);
	  ## lambda2 = (1 / VarData) / sqrt(length(SMS$YY));
	  ##alpha2 = 1/VarData;
	  ##alpha2Vals <- c(alpha2 * .01, alpha2 * .1, alpha2, 
    ##  alpha2 * 10, alpha2 * 100, alpha2 * 1000, 1);
    alpha2Vals <- c(.011, .1, .25, .5, .75, .9, .99, 1);
	  MDD <- list();
	  for (ii in 1:length(alpha2Vals)) {
	    MDD[[ii]] = list();
      try(MDD[[ii]] <- cv.glmnet(SMS$XX, SMS$YY, #K=10, 
        family = "binomial", alpha=alpha2Vals[ii], #mode="step", 
        s=1:min(c(SMS$kActiveLen *2, SMS$kLen-1, length(SMS$YY)-1)), 
        #normalize=TRUE, intercept=FALSE, trace = FALSE, eps = .00001
      ));
      MDD[[ii]]$cv = MDD[[ii]]$cvm;
    }
    MaxINF =  99999999999999999
	  Minerizer <- rep(MaxINF , length(alpha2Vals));
	  for (ii in 1:length(alpha2Vals)) {
	    if (!is.null(MDD[[ii]]$cv))  {
        Minerizer[ii] = min(MDD[[ii]]$cv,na.rm = TRUE);
      }
    }
    MyWant <- sort(Minerizer, index=TRUE)$ix[1];
    if (min(Minerizer) == MaxINF) {
      AFilePrint("AllMinierizers were to effing big")
      HitENET2Ob = list();
      HitENET2Ob$beta = rep(0, kLen);
      UseAns = 1;
      BetaSM = NULL;
    }  else {
      ##UseAns <- sort(MDD[[MyWant]]$cv, index=TRUE)$ix[1];
      AFilePrint("Run One More HitENET2Obv");
      HitENET2ObCV <-  cv.glmnet(SMS$XX, SMS$YY, family = "binomial", 
        lambda = NULL, alpha = alpha2Vals[MyWant]
        ##max.steps= max(MDD[[MyWant]]$s[UseAns],2), normalize=TRUE, 
        ##intercept=FALSE, trace = FALSE, eps = .00001
        );  
      MyST = sort(HitENET2ObCV$cvm, index=TRUE)$ix[1];
      AFilePrint(paste("MyST = ", MyST, sep="")); flush.console();
      HitENET2Ob <- glmnet(SMS$XX, SMS$YY, family = "binomial",
        lambda = HitENET2ObCV$lambda[MyST], alpha= alpha2Vals[MyWant]);
      AFilePrint(paste("Using MyST = ", MyST," , lambda = ", 
        HitENET2ObCV$lambda[MyST], "we get HitENET2Ob$beta = ", sep=""));
        flush.console();
      AFilePrint(paste("names(HitENET2Ob) = ", paste(names(HitENET2Ob), collapse=", "),
       sep=""));  flush.console();
      AFilePrint(paste("Is HitENET2Ob$beta  vector?  ", 
        is.vector(HitENET2Ob$beta), sep="")); flush.console();
      AFilePrint(paste(as.numeric(HitENET2Ob$beta), collapse=", "));flush.console();
      ##HitENET2Ob$beta = HitENET2Ob$beta;
    }

    ENET2Fit1 <-  as.vector(as.numeric(HitENET2Ob$beta));
    ENET2Fit = ENET2Fit1;
    ENET2Fit[ ENET2Fit1 != 0 ] = 1;
    ENET2Fit2 = ENET2Fit1;
    if (length(ENET2Fit1) != kLen)  {
      AFilePrint(paste("Error, ENET2Fit1 is length ", length(ENET2Fit1),
        " but not kLen = ", kLen, sep="")); flush.console();
    }
    if (sum(ENET2Fit) == 0) {
      BetaSM = NULL;
	    T2 <- proc.time();
	    ENET2FitTime <- T2 - T1;    
      return(list( type="ENetCVMin", BetaFit = ENET2Fit2, BBHatsFit = ENET2Fit,
	    FitTime =ENET2FitTime, OtherBetas = ENET2Fit1,
      HitENET2Ob=HitENET2Ob, ITO=MyST, ENET2Fit=ENET2Fit, ENET2Fit2=ENET2Fit2,
	    BetaSM=BetaSM ) );      
    }
    BetaSM = ENET2Fit1[ENET2Fit == 1];
    if (sum(ENET2Fit) < length(SMS$YY)) {
      AFilePrint("Running last glm"); flush.console();
      AFilePrint(paste("ENET2Fit = ", paste(ENET2Fit, collapse=", "), sep=""));
      glm2 = NULL;
      try(glm2 <- glm(SMS$YY~SMS$XX[,ENET2Fit == 1], family=binomial(link="logit")));
      if (!is.null(glm2)) { BetaSM = coefficients(glm2); }
      if (!is.null(coefficients(glm2))) {
        BetaSM = coefficients(glm2);
        BetaSM = BetaSM[2:length(BetaSM)];
      }  else {
        AFilePrint("That BetaSM is Null, so we're in trouble")
      }
      ENET2Fit2 = ENET2Fit;
      AFilePrint("Putting in BetaSM"); flush.console();
      ENET2Fit2[ENET2Fit == 1] = BetaSM;     
    }
    
	  T2 <- proc.time();
	  ENET2FitTime <- T2 - T1;    
    return(list( type="ENetCVMin", BetaFit = ENET2Fit2, BBHatsFit = ENET2Fit,
	  FitTime =ENET2FitTime, OtherBetas = ENET2Fit1,
    HitENET2Ob=HitENET2Ob, ITO=MyST, ENET2Fit=ENET2Fit, ENET2Fit2=ENET2Fit2,
	  BetaSM=BetaSM, MIPReport = ENET2Fit) );  
  }
  if (UseAns == 1) {UseAns = 2;}
  ITO = UseAns;

  ENET2Fit <-  rep(0, kLen);
  ENET2Fit[ HitENET2Ob$allset[HitENET2Ob$beta[UseAns,] != 0]] = 1;   
  ##ENET2Fit <- HitENET2Ob$beta[UseAns,];
	##ENET2Fit[abs(ENET2Fit) > MinMin] <- 1;
	##ENET2Fit[abs(ENET2Fit) <= MinMin] <- 0;  
	ENET2Fit2 <- ENET2Fit;
  BetaSM <- 0;
	if (sum(ENET2Fit) > 0 || length(SMS$XX[1,ENET2Fit == 1]) > 0) { 
	  XXSM <-  SMS$XX[,ENET2Fit == 1];
	  XXSXX <- t(XXSM) %*% XXSM;
	  if (any(is.null(XXSXX)) || 
      length(dim(XXSXX)) != 2 || dim(XXSXX)[1] != dim(XXSXX)[2] || 
	    dim(XXSXX)[1] < 1 || any(is.nan(XXSXX))) { 
	    ENET2Fit2 <- ENET2Fit; BetaSM <- 0;
	  } else if (length(XXSXX) == 1) {
	    BetaSM <- 1 / XXSXX * sum(XXSM * SMS$YY);
	  } else if ( !is.null(dim(XXSXX)) && dim(XXSXX)[1] > 1 && 
      dim(XXSXX)[1] == dim(XXSXX)[2] && det(XXSXX) <= 0.00000001 ) {
	    library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
	    BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;
	  } else if (!is.null(dim(XXSXX)) && det(XXSXX) > .00000001) {
		  LS <- svd(XXSXX);
		  if (max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
        BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;  
      } else {
		    BetaSM <- try(solve(XXSXX) %*% t(XXSM) %*% SMS$YY);
		  }		        
	  } else {
 	    library(corpcor, warn.conflicts = FALSE, quietly=TRUE);
	    BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;         
    }
	  ENET2Fit2[ENET2Fit != 0 ] <- BetaSM;
	} else {
	  ENET2Fit2 <- ENET2Fit;  BetaSM <- 0;
	}

	T2 <- proc.time();
	ENET2FitTime <- T2 - T1;
	##   AFilePrint("Finished RunningGenerateENetCVMin");
	##   flush.console();
	return(list( type="ENetCVMin", BetaFit = ENET2Fit2, BBHatsFit = ENET2Fit,
	  FitTime =ENET2FitTime, OtherBetas = HitENET2Ob$beta[ITO,],
    HitENET2Ob=HitENET2Ob, ITO=ITO, ENET2Fit=ENET2Fit, ENET2Fit2=ENET2Fit2,
	  BetaSM=BetaSM, MIPReport = ENET2Fit ) );                  
}      

################################################################################
##  This generates the LarsCp fit.
##
##  Note, for Logit regression, we'll use BIC minimizer
##
GenerateLarsCp <- function(SMS, DoCI = DefaultDoCI,...) {
  library(lars);
  kLen = SMS$kLen;
  n = NULL; m = NULL;
  try(n <- length(SMS$YY)); 
  try(m <- length(SMS$XX[1,]));
  if (is.null(n)) {
    AFilePrint("LarsCp, N is Null!"); flush.console();
    return(NULL);
  }
  if (is.null(m)) {
    AFilePrint("LarsCp, m is NULL!");  flush.console();
    return(NULL);
  }
	T1 <- proc.time(); 
	##AFilePrint("Starting the LarsCp"); flush.console();
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE)  { 
    if (kLen > 1000) { 
	    HitLars2Ob <- lars(SMS$XX, SMS$YY, use.Gram=FALSE); 
	  } else {
      HitLars2Ob <- lars(SMS$XX, SMS$YY, use.Gram=TRUE);     
    }
    Cp <- HitLars2Ob$Cp;
    if (any(is.nan(Cp)) || any(is.na(Cp)) ) {
		  Cp <-  vector("numeric", length(HitLars2Ob$beta[,1]));
	    for (jj in 1:length(Cp)) {
	      Cp[jj] <- ( sum ( (SMS$YY - SMS$XX %*% as.vector(HitLars2Ob$beta[jj,])
          )^2 ) / SMS$SigmaSqNoise  - length(SMS$YY) +
          2 * length(HitLars2Ob$beta[jj, HitLars2Ob$beta[jj,] != 0]));
	    }      
	  }
	  BetaSM <- 0;
	  if (any(is.nan(Cp)) || any(is.na(Cp))) {
      LarsCpFitAll <- HitLars2Ob$beta[1,] * 0 - 999;
	    LarsCpFitOnes <- LarsCpFitAll; 
	    LarsCpFitMP <- LarsCpFitAll;	      
	  } else if (length(Cp) != 1 && Cp[1] != -999) {
	    STL <- sort(Cp, index=TRUE);
		  LarsCpFitAll <- HitLars2Ob$beta[STL$ix[1],];
		  LarsCpFitOnes <- LarsCpFitAll;
		  LarsCpFitOnes[abs(LarsCpFitAll) > MinMin] <- 1;
		  LarsCpFitOnes[abs(LarsCpFitAll) <= MinMin] <- 0;
		  LarsCpFitMP <- LarsCpFitOnes;
		  if (sum(abs(LarsCpFitOnes)) == 0) {
        LarsCpFitMP = 0 * LarsCpFitAll;
      }  else if (sum(LarsCpFitOnes) > 0 || 
        length(SMS$XX[1, LarsCpFitOnes == 1]) > 0)  { 
		    XXSM <- SMS$XX[, LarsCpFitOnes == 1];
		    MyS = NULL;
		    try(MyS <- lm(SMS$YY~XXSM-1));
		    if (!is.null(MyS)) {
          BetaSM = coefficients(MyS);
        } else {
          XXSXX <- t(XXSM) %*% XXSM;
          try(BetaSM <- solve(XXSXX) %*% t(XXSM) %*% SMS$YY);
          if (is.null(BetaSM)) {
            library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
            try(BetaSM <- pseudoinverse(XXSXX) %*% t(XXSM) %*% SMS$YY);
          }
          if (is.null(BetaSM)) {
            BetaSM = rep(0, length(XXSM[1,]));
          }
        }
        LarsCpFitMP = 0 * LarsCpFitOnes;
        LarsCpFitMP[LarsCpFitOnes == 1] = BetaSM;
      }
    }
  } else {
   library(glmpath, warn.conflicts=FALSE, quietly=TRUE)
   n <- length(SMS$YY); m <-length(SMS$XX[1,]);
   HitLars2Ob = NULL;
   ##AFilePrint("LarsCp: About to run glmpath"); flush.console();
   n <- length(SMS$YY); m <- NCOL(SMS$XX);
   try(HitLars2Ob <- glmpath(SMS$XX, SMS$YY,
      nopenalty.subset = NULL, family = binomial,
      weight = rep(1, n), offset = rep(0, n), lambda2 = 1e-5,
      max.steps = 10*min(n, m), max.norm = 100*m,
               min.lambda = (if (m >= n) 1e-6 else 0), 
      max.vars = min(length(SMS$YY)-1, length(SMS$XX[1,])),
      max.arclength = Inf, frac.arclength = 1, add.newvars = 1,
      bshoot.threshold = 0.1, relax.lambda = 1e-8,
      standardize = TRUE,
      eps = .Machine$double.eps, trace = FALSE));
   if (is.null(HitLars2Ob)) {
     ##AFilePrint("glmpath fails for finding LarsCp on this animal");
     ##flush.console();
     LarsCpFitOnes = rep(0, m);
     LarsCpFitMP = rep(0, m);
     LarsCpFitAll = rep(0,m); 
	  T2 <- proc.time();
	  LarsCpT <- T2 - T1;
	  return(list( type="LarsCp", BetaFit = LarsCpFitMP, BBHatsFit = LarsCpFitOnes,
	    FitTime =LarsCpT, OtherBetas = LarsCpFitAll, MIPReport = LarsCpFitOnes) );      
   } 
   ##AFilePrint("LarsCp: Made it to ask for bic"); flush.console();  
   Cp <- HitLars2Ob$bic;
   if (is.null(Cp)) {
     LarsCpFitOnes = rep(0, m);
     LarsCpFitMP = rep(0, m);
     LarsCpFitAll = rep(0,m); 
     T2 <- proc.time();
     ##AFilePrint("glmpath didn't create nonnull Cp?"); flush.console();
	   LarsCpT <- T2 - T1;
	   return(list( type="LarsCp", BetaFit = LarsCpFitMP, BBHatsFit = LarsCpFitOnes,
	     FitTime =LarsCpT, OtherBetas = LarsCpFitAll, MIPReport = LarsCpFitOnes) );   
   } else {
     ##AFilePrint("LarsCp: Extracting Betas"); flush.console();
     HitLars2Ob$beta = HitLars2Ob$b.predictor[,2:NCOL(HitLars2Ob$b.predictor)];
     MyX = sort(Cp, index=TRUE)$ix[1];
     if (is.null(HitLars2Ob$beta)) {
      LarsCpFitMP =  rep(0,m);
      LarsCpFitOnes = rep(0,m);
      LarsCpFitAll = rep(0, m);
      AFilePrint("HitLars2Ob$beta was null, can't work!");
      T2 <- proc.time();
      AFilePrint("glmpath didn't create nonnull Cp?");  flush.console();
	    LarsCpT <- T2 - T1;
	    return(list( type="LarsCp", BetaFit = LarsCpFitMP, BBHatsFit = LarsCpFitOnes,
	     FitTime =LarsCpT, OtherBetas = LarsCpFitAll, MIPReport = LarsCpFitOnes) );     
     }
     ##AFilePrint(paste("LarsCp: Extracting UseX, MyX = ", MyX, sep="")); flush.console();
     ##AFilePrint(paste("dimbeta = c(", 
     ##  paste(dim(HitLars2Ob$beta), collapse=", "), ")", sep="")); flush.console();
     ##AFilePrint(paste("HitBeta = c(",
     ##  paste(HitLars2Ob$beta[MyX,], collapse=", "), ")", sep="")); flush.console();
     UseX = (1:m)[HitLars2Ob$beta[MyX,]!= 0];
     ##AFilePrint(paste("UseX = c(",
     ##  paste(UseX, collapse=", "), ")", " and m = ", m, sep="")); flush.console();     
     LarsCpFitAll = HitLars2Ob$beta[MyX,UseX];
     LarsCpFitMP =  LarsCpFitAll;
     LarsCpFitOnes = rep(0,m);
     LarsCpFitOnes[UseX] = 1;
     AFilePrint(paste("LarsCp: length of UseX = ", UseX, sep="")); flush.console();     
     if (length(UseX) > 0) {
       LarsCpFitOnes = rep(0,m);
       LarsCpFitOnes[UseX] = 1;
       MyGlm = NULL;
       try(MyGlm <- glm(SMS$YY~SMS$XX[,LarsCpFitOnes == 1], family=binomial(link="logit")));
       if (!is.null(MyGlm)) {
         BetaSM = NULL;
         try(BetaSM <- coefficients(MyGlm)[2:length(MyGlm$coefficients)]);
         if (is.null(BetaSM)) {
           T2 <- proc.time();
           AFilePrint("BetaSM did not return")
	         LarsCpT <- T2 - T1;
	         return(list( type="LarsCp", BetaFit = LarsCpFitMP, BBHatsFit = LarsCpFitOnes,
	          FitTime =LarsCpT, OtherBetas = LarsCpFitAll, MIPReport = LarsCpFitOnes) );      
         }
         LarsCpFitMP = rep(0,m);
         if (any(is.na(UseX))) {
           AFilePrint("We've got to deal with na UseX"); flush.console();
           try(
             LarsCpFitMP[UseX[!is.na(UseX)]] <- 
               BetaSM[(1:length(UseX[!is.na(UseX)]))] );
           AFilePrint("We've dealt with na UseX"); flush.console();
         } else {
           LarsCpFitMP[UseX] = BetaSM;
           AFilePrint("A successful fill of LarsCpFitMP"); flush.console();
         }
         T2 <- proc.time();
      	 LarsCpT <- T2 - T1;
	       return(list( type="LarsCp", BetaFit = LarsCpFitMP, BBHatsFit = LarsCpFitOnes,
	        FitTime =LarsCpT, OtherBetas = LarsCpFitAll, MIPReport = LarsCpFitOnes) );    
       }
     } else {
       T2 <- proc.time();
       LarsCpFitOnes = rep(0,m);
       LarsCpFitAll = rep(0,m);
       LarsCpFitMP = rep(0,m);
       AFilePrint("Length UseX == 0"); flush.console();
	     LarsCpT <- T2 - T1;
	     return(list( type="LarsCp", BetaFit = LarsCpFitMP, BBHatsFit = LarsCpFitOnes,
	       FitTime =LarsCpT, OtherBetas = LarsCpFitAll, MIPReport = LarsCpFitOnes) );        
     
     }
   }
  }

	T2 <- proc.time();
	LarsCpT <- T2 - T1;
	return(list( type="LarsCp", BetaFit = LarsCpFitMP, BBHatsFit = LarsCpFitOnes,
	  FitTime =LarsCpT, OtherBetas = LarsCpFitAll, MIPReport = LarsCpFitOnes) );         
}


GenerateLassoQLY <- function(SMS, DoCI = DefaultDoCI,...) {
  library(lars);
	    ##AFilePrint(paste("ii4 = ", ii4, ", starting with lars"));
	     T1 <- proc.time();
	     BetaSM <- 0;
        if (length(SMS$XX[1,]) > 1000 && length(SMS$XX[1,]) > length(SMS$XX[,1]) 
           ) { 
	      HitLars2Ob <- lars(SMS$XX, SMS$YY, use.Gram=FALSE); 
	      } else {
        HitLars2Ob <- lars(SMS$XX, SMS$YY, use.Gram=TRUE);     
        }
		   ##   HitLars2Ob <- lars(SMS$XX, SMS$YY);
		       q <- SMS$puse;
		       lambdaLook <- 1 * (1-q) / q  * 2 * sqrt( 2 * SMS$SigmaSqNoise);
		       LTS <- max(HitLars2Ob$lambda);
		       if (LTS < lambdaLook) { 
		            Lars4Fit <- HitLars2Ob$beta[1,] * 0;
		            Lars4Fit4A <- HitLars2Ob$beta[1,] * 0;
		       } else { 
		           ITO22 <- (1:length(HitLars2Ob$lambda))[HitLars2Ob$lambda > lambdaLook];
		           ITOOn <- (1:length(HitLars2Ob$lambda))[HitLars2Ob$lambda == lambdaLook];
			       ITO33 <- (1:length(HitLars2Ob$lambda))[HitLars2Ob$lambda < lambdaLook];
			       minLambda <- min(HitLars2Ob$lambda);
			       maxLambda <- max(HitLars2Ob$lambda);
			       if (minLambda > lambdaLook) {
			          Lars4Fit <- rep(1, length(HitLars2Ob$beta[1,]));
			          Lars4Fit4A <- HitLars2Ob$beta[ length(HitLars2Ob$beta[,1]),];
			       } else if (maxLambda < lambdaLook) {
			          Lars4Fit <- rep(0, length(HitLars2Ob$beta[1,]));
			          Lars4Fit4A <- rep(0, length(HitLars2Ob$beta[1,]));
			       } else if ( length(HitLars2Ob$lambda[HitLars2Ob$lambda > 0]) <= 0) {
			          Lars4Fit <- HitLars2Ob$beta[1,] * 0;
			          Lars4Fit4A <- rep(0, length(HitLars2Ob$beta[1,]));
			       } else if (length(ITOOn) > 0) {
			          Lars4Fit <- HitLars2Ob$beta[ ITOOn[1] ];
			          Lars4Fit[abs(Lars4Fit) >= MinMin] <- 1;
			          Lars4Fit[abs(Lars4Fit) < MinMin] <- 0;
			          Lars4Fit4A <- HitLars2Ob$beta[ITOOn[1],];
			       } else if (length(ITO33) >= 1) {
			         Lars4Fit <- HitLars2Ob$beta[ min(ITO33),];
			         Lars4Fit[abs(Lars4Fit) >= MinMin] <- 1;
			         Lars4Fit[abs(Lars4Fit) < MinMin] <- 0;	
			         Lars4Fit4A <- HitLars2Ob$beta[ min(ITO33),];		       			          
			       } else if (max(ITO22)+1 <= length(HitLars2Ob$beta[,1]) ) {
			         Lars4Fit <- HitLars2Ob$beta[ max(ITO22)+1,];
			         Lars4Fit[abs(Lars4Fit) >= MinMin] <- 1;
			         Lars4Fit[abs(Lars4Fit) < MinMin] <- 0;
			         Lars4Fit4A <- HitLars2Ob$beta[max(ITO22),];
			       } else {
			          AFilePrint("I Can't figure out what's going");
			          Lars4Fit <- rep(-999, length(HitLars2Ob$beta[1,]));
			       }
			   }
		       Lars4Fit4 <- Lars4Fit;
		           if (sum(Lars4Fit) > 0) { 
			             XXSM <- SMS$XX[, Lars4Fit == 1];
			             XXSXX <- t(XXSM) %*% XXSM;
			             if (any(is.null(XXSXX)) || length(dim(XXSXX)) != 2 || dim(XXSXX)[1] != dim(XXSXX)[2] || 
	                           dim(XXSXX)[1] < 1 || any(is.nan(XXSXX))) { 
	                           Lars4Fit4 <- Lars4Fit;
	                           BetaSM <- 0;
	                    } else if (length(XXSXX) == 1) {
	                           BetaSM <- 1 / XXSXX * sum(XXSM * SMS$YY);
		                } else if (!is.null(dim(XXSXX)) && dim(XXSXX)[1] > 1 && dim(XXSXX)[1] == dim(XXSXX)[2] && 
		                            det(XXSXX) <= 0.00000001 ) {
		                 library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
		                 BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;
		                } else if (!is.null(dim(XXSXX)) && det(XXSXX) > .00000001){
		                  BetaSM <- solve(t(XXSM) %*% XXSM) %*% t(XXSM) %*% SMS$YY;
		                } else {
  		                library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
		                  BetaSM <- pseudoinverse( XXSXX) %*% t(XXSM) %*% SMS$YY;                  
                    }
			             Lars4Fit4[Lars4Fit != 0 ] <- BetaSM; 	       
		            } else {
		               Lars4Fit4 <- Lars4Fit;
		            }
           T2 <- proc.time();
           TLars4 <- T2 - T1;	
	         return(list( type="LassoQLY",
	                BetaFit = Lars4Fit4, BBHatsFit = Lars4Fit,
	                FitTime =TLars4,
	                OtherBetas = Lars4Fit4A, MIPReport = Lars4Fit) );                
           
}
GenerateEMRidge <- function(SMS, LARSSeekFlag, DoCI = DefaultDoCI,...) {
	       ##USENOISE <- sqrt(SMS$SigmaSqNoise);  ##USENOISE <- SMS$SigmaSqNoise
	       USENOISE <- SMS$SigmaSqNoise;
   library(TwoLassoCpp);
        T1 <- proc.time();
   if (LARSSeekFlag == 0) {
	      StlambdaA = .5/ (2*USENOISE);
	      StlambdaD = .5 / (2*USENOISE);
      
	  } else if (LARSSeekFlag >= 1) {
	      MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)
	      rnd <- round(SMS$puse * kLen * LARSSeekFlag); rnd = min( kLen, rnd);
	      ITO <- MyLarsGiveMeClose(MyFit, 
	            IntClose = rnd, MinMin = 0.000000001);   
	     StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	     StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	     ##StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE);
	     ##StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE);
	  } else if (LARSSeekFlag < 0) {
	     rnd <- round(SMS$kActiveLen);
       StlambdaA = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;
       StlambdaD = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;  
    } 
      tauAsqStart <- 1 / StlambdaA ;
      tauDsqStart <- 2 / StlambdaD;
      if (!is.finite(tauAsqStart)) {
         tauAsqStart = .25 * sd(as.vector(SMS$YY)) / min(apply(SMS$XX,2,sd));
      }
      if (!is.finite(tauDsqStart)) {
         tauDsqStart = .3 * sd(as.vector(SMS$YY)) / min(apply(SMS$XX,2,sd));
      }
        ##AFilePrint(paste("ii4 = ", ii4, ", starting with EMRIDGE"));
        
	    EMCFit = EMRIDGE(yys=SMS$YY, xxs=SMS$XX, n1 = 100, n2 = 4, 
	          ppiuse = SMS$puse,
	     tauDsqStart = tauDsqStart, tauAsqStart = tauAsqStart, 
	     tauDMult = 2^(-1/8), tauAMult = 2^(1/8), 
	     tauAMultStop = 50, BStartBetas = -1, 
	     SigmaSqNoise = USENOISE, EFlag = FALSE,
	     logp = log(SMS$puse), log1mp = log(1-SMS$puse), 
	     SigmaNoiseInv = 1/USENOISE, NOPRINT = FALSE);

	  if (!is.null(EMCFit) && !is.null(EMCFit$FailureFlag) &&
            is.numeric(EMCFit$FailureFlag) && EMCFit$FailureFlag ==1 ) {
	     AFilePrint("EMCFit Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
         FailToFitEMOb[[OnFails]] <<- EMCFit;	     
	  }                
	   
	  EMCFitBBs <- as.vector(EMCFit$BBHatsFinal);
	  EMCFitBBs[abs(EMCFit$BBHatsFinal) <= MinMin] = 0;
	  EMCFitBBs[abs(EMCFit$BBHatsFinal) > MinMin] = 1;
	   T2 <- proc.time();
	     EMFitT <- T2 - T1;
	   ##AFilePrint(paste("ii4 = ", ii4, ", starting with HitBatWingFast"));
	      if (PrintOutFlags >= 2) {  
	        AFilePrint(paste("ii4 = ", ii4, ", starting with HitBatWingFast"));
	          flush.console();
	      }	  
		 return(list( type="EMRidge",
	                BetaFit = EMCFit$BBHatsFinal, BBHatsFit = EMCFitBBs,
	                FitTime =EMFitT,
	                OtherBetas = EMCFit$BBHatsFinal,
                  EMObject = EMCFit, MIPReport = EMCFit$MIP) );  
}
GenerateTwoLasso9X <- function(SMS, LARSSeekFlag, DoCI = DefaultDoCI,...) {
    library(TwoLassoCpp);
    USENOISE <- SMS$SigmaSqNoise; kLen <- SMS$kLen;
	  T1 <- proc.time();
	  if (LARSSeekFlag == 0) {
	      StlambdaA = .5/ (2*USENOISE);
	      StlambdaD = .5 / (2*USENOISE);     
	  } else if (LARSSeekFlag >= 1) {
	      MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)
	      rnd <- round(SMS$puse * kLen * LARSSeekFlag); rnd = min( kLen, rnd);
	      ITO <- MyLarsGiveMeClose(MyFit, 
	            IntClose = rnd, MinMin = 0.000000001);   
	     StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	     StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	     ##StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE);
	     ##StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE);
	  } else if (LARSSeekFlag < 0) {
	     rnd <- round(SMS$kActiveLen);
       StlambdaA = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;
       StlambdaD = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;  
    } 
    if (is.null(SMS$LogitRegresion) || SMS$LogitRegression == FALSE) {
	    HitBatWingFast <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, 
	      FixKa = -100, sigmaNoiseSq = USENOISE,
	      RatWant = .2, StlambdaD=StlambdaD, StlambdaA=StlambdaA, 
	      lambdaDMultC = 3^(1/2), lambdaAMultC = 3^(-1/4),
	      lambdaAmultStop = 5, TotalRuns = 5, 
	      NumEMConv = 4, MultEMCons = .99, NumCDOConv = 40, CDOEpsilon = .000001);
    } else if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
      LambdaAK = StlambdaA * lambdaAMultC^(0:4);
      LambdaDK = StlambdaD * lambdaDMultC^(0:4);
      OrderSeq = rep(10,5);
      
       if (SMS$n < SMS$p) {
        LogitNoise = SMS$n/SMS$p;
      } else { LogitNoise = 1.0; } 
      HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = LogitNoise, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = c(5,1,1), Verbose = 0, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = -1, DoLogit=TRUE, TargetMinimumBeta = .75, LogitNoise = n/p);  
      ##HitBatWingFast = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
      ##  StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
      ##  LambdaAK = LambdaAK,LambdaDK = LambdaDK,OrderSeq = OrderSeq,
      ##  StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
      ##  MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
      ##  TotalRuns = 5, MaxLambdaASteps = 50,
      ##  RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
      ##  MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = -1)  
    }
	  if (!is.null(HitBatWingFast$FailureFlag) && 
      length(HitBatWingFast$Failureflag) >= 1 &&
      HitBatWingFast$FailureFlag[1] == 1 ) {
	     AFilePrint("HitBatWingFast Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
	     FailToFitEMOb[[OnFails]] <<-HitBatWingFast;
	  }  	
    if (length(HitBatWingFast$ReturnBetas) > 1) {
	  LarsFitBBsFast <- as.vector(HitBatWingFast$ReturnBetas);
	  LarsFitBBsFast[abs(HitBatWingFast$ReturnBetas) <= MinMin] = 0;
	  LarsFitBBsFast[abs(HitBatWingFast$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBsFast = SMSReal * 0 - 10;
	  }
	  T2 <- proc.time();
	     EMLARSFastT <- T2 - T1;                      
		 return(list( type="TwoLasso9X",
	                BetaFit = HitBatWingFast$ReturnBetas,
                  BBHatsFit = LarsFitBBsFast,
	                FitTime =EMLARSFastT,
	                OtherBetas = HitBatWingFast$ReturnBetas,
                  EMLARSObject = HitBatWingFast, MIPReport = HitBatWingFast$BBHatsFinal) );
}


GenerateSCAD <- function(SMS, LARSSeekFlag, DoConvexMin=1, DoCI = DefaultDoCI,...) {
    library(ncvreg, warn.conflicts=FALSE, quietly=TRUE);
     T1 = proc.time();
	     ##AFilePrint(paste("ii4 = ", ii4, ", starting with longer EMLARS"));
    USENOISE <- SMS$SigmaSqNoise;	  kLen <- SMS$kLen;  
    rnd <- round(SMS$puse * kLen);  rnd = min( kLen, rnd);
    rnd <- SMS$kActiveLen;
    if (LARSSeekFlag < 0) {
       UseGamma = sqrt(length(SMS$YY)) * sqrt(USENOISE) / rnd / abs(LARSSeekFlag);
    } else {
       UseGamma = sqrt(length(SMS$YY)) * sqrt(USENOISE) / rnd;
    }
      ##sXX = StandardizeXX(SMS$XX);
      sXX = SMS$XX;
      ## sYY = StandardizeXX(SMS$YY);
      sYY = SMS$YY;
      ###sdX = sd(SMS$XX);  ##sdY = sd(SMS$YY);
    
    if (length(sYY) == 0) {
		 return(list( type="SCADLARSeek",
	                BetaFit = rep(0,kLen),
                  BBHatsFit = rep(0,kLen),
	                FitTime =c(0,0,0,0,0),
	                OtherBetas = rep(0,kLen),
	                LARSSeekFlag = LARSSeekFlag,
                  EMLARSObject = NULL, MIPReport = BBHatsFit) );	       
    }
    ##if (length(BetaStart) > 1 && length(BetaStart) == kLen) {
    ##  HitSCAD <- scad.onestep(BetaStart,y=sYY,x = sXX,UseGamma);
     ## HitLasso <- CoordinateDescent(sXX,sYY, OnGamma = UseGamma,
     ##    NumCDOConv = 10000,  CDOEpsilon = .00000001, OnBeta = BetaStart); 
    ##} else {
    ##  HitSCAD <- scad.onestep(rep(0, length(sXX[1,]) ),y=sYY,x = sXX,UseGamma);
     ## HitLasso <- CoordinateDescent(sXX,sYY, OnGamma = UseGamma,
     ##    NumCDOConv = 10000,  CDOEpsilon = .00000001 );     
    ##}
    
      if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE ) {
       fit2 <- ncvreg(X = sXX, y= sYY, gamma = 3.7, penalty="SCAD");
      } else {
       fit2 <- ncvreg(X = sXX, y= sYY, gamma = 3.7, penalty="SCAD", 
         family = "binomial");      
      }
      if (DoConvexMin == 1) {
        HitSCAD = fit2$beta[2:length(fit2$beta[,1]), fit2$convex.min];
      } else if (DoConvexMin==0) {
        LMin <- sort( abs(fit2$lambda - UseGamma), index=TRUE)$ix[1];
        HitSCAD = fit2$beta[2:length(fit2$beta[,1]), LMin]    
      } else {
         ## WLT Method  (Wang, Hansheng and Li, Runze and Tsai, Chih-Ling, 2007)
         sigmaSqF2 <- fit2$lambda * 0;
         dfF2 <- fit2$lambda * 0;        
         for (ii in 1:length(fit2$lambda)) {
             BetaNot0 <- fit2$beta[2:length(fit2$beta[,1]),ii];
             Errors <- sYY -fit2$beta[1,ii]- sXX %*% BetaNot0;
             IIX <- (1:length(BetaNot0))[BetaNot0 != 0];
             if (length(IIX) > 0) {
                     BetaD <- BetaNot0[IIX];
                     Pen = rep(0, length(BetaD));
                     Pen[BetaD < fit2$lambda[ii]] <- fit2$lambda[ii];
                     Pen[BetaD >= fit2$lambda[ii] & BetaD < (fit2$gamma * fit2$lambda[ii]) ] <-
                       (fit2$gamma * fit2$lambda[ii] - 
                           BetaD[ BetaD >= fit2$lambda[ii] & BetaD < (fit2$gamma * fit2$lambda[ii]) ])  /
                        (fit2$gamma  - 1);
                     SigmaLambda <-matrix(0, length(BetaD), length(BetaD));
                     diag(SigmaLambda) = ( Pen / abs(BetaD) )
                     DFlambda = sum(diag( sXX[,IIX] %*% solve(
                       t(sXX[,IIX]) %*% sXX[,IIX] +
                        n * SigmaLambda ) %*% t(sXX[,IIX])
                          ) );
             } else {
                  DFlambda = 0;
             }
             dfF2[ii] <-DFlambda
             sigmaSqF2[ii] = sum (Errors^2) / n;
          }
          MinLambdaii <- sort( log(sigmaSqF2) + dfF2 * log(n)/n,index=TRUE)$ix[1];
          HitSCAD = fit2$beta[2:length(fit2$beta[,1]), MinLambdaii]    
      }
    HitSCAD = as.numeric(HitSCAD);
    if (length(HitSCAD) != length(sXX[1,]) || is.null(HitSCAD) ||
        any(!is.finite(HitSCAD)) ) {
      BBHatsFit = sXX[1,] * 0;
      HitSCAD = BBHatsFit;    
    }  else {
      BBHatsFit = HitSCAD;
      BBHatsFit[HitSCAD != 0] = 1; 
    } 
      
    if (sum(BBHatsFit) >= 1 && !is.null(SMS$LogitRegression) && 
      (SMS$LogitRegression == TRUE || SMS$LogitRegression == 1)) {
      glm2 = NULL;
      try(glm2 <- glm(SMS$YY~SMS$XX[,BBHatsFit == 1], family=binomial(link="logit")));
      if (!is.null(glm2)) { BetaSM = coefficients(glm2); }
      if (!is.null(coefficients(glm2))) {
        BetaSM = coefficients(glm2);
        BetaSM = BetaSM[2:length(BetaSM)];
      }  else {
        AFilePrint("That BetaSM is Null, so we're in trouble")
      }  
      BetaFit <- rep(0, NCOL(SMS$XX));
      BetaFit[BBHatsFit == 1] <- BetaSM;
	} else  if (sum(BBHatsFit) >= 1) {
		             XXSM <- SMS$XX[, BBHatsFit == 1];
		             XXSXXST <- t(XXSM) %*% XXSM;
		             if (any(is.null(XXSXXST)) || length(dim(XXSXXST)) != 2 
		                || dim(XXSXXST)[1] != dim(XXSXXST)[2] || 
		                dim(XXSXXST)[1] < 1 || any(is.nan(XXSXXST))) { 
		                BTF1 <- HitSCAD[BBHatsFit == 1];                
		             } else if (length(XXSXXST) == 1) {
		                BTF1 <- 1 / XXSXXST * sum(XXSM * SMS$YY);
		             } else if (det(XXSXXST) <= .00000001) {
 		                library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                     ##   "Determinant is Samll, Going with corpcor", sep=""));
		                BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, passed corpcor", sep=""));                                  
		             } else if (!is.null(dim(XXSXXST)) && dim(XXSXXST)[1] > 1 && 
		                        dim(XXSXXST)[1] == dim(XXSXXST)[2] && det(XXSXXST) <= 0.00000001 ) {
 		                library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, Going with corpcor", sep=""));
		                BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, passed corpcor", sep=""));   
		             } else if (!is.null(dim(XXSXXST)) && det(XXSXXST) > .0001 ) {
		                LS <- svd(XXSXXST);
		                if (max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
                         BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
                    }   else {
		                     BTF1 <- try(solve(XXSXXST) %*% t(XXSM) %*% SMS$YY);
		                }
		             }  else {
 		                library(corpcor,warn.conflicts=FALSE, quietly=TRUE);
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, Going with corpcor", sep=""));
		                BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, passed corpcor", sep=""));                   
                 }
      BetaFit = BBHatsFit;
      BetaFit[BBHatsFit == 1] = BTF1;     
	  } else {
	    BetaFit = BBHatsFit * 0;
	  }
	   T2 <- proc.time();
	   StoddenFit1T <- T2 - T1;
		 return(list( type="SCAD",
	                BetaFit = BetaFit,
                  BBHatsFit = BBHatsFit,
	                FitTime =StoddenFit1T,
	                OtherBetas = HitSCAD,
                  EMLARSObject = NULL, MIPReport = BBHatsFit) );	   
}



GenerateMCP <- function(SMS, LARSSeekFlag, DoConvexMin=1, DoCI = DefaultDoCI,...) {
    library(ncvreg, warn.conflicts=FALSE, quietly=TRUE);
     T1 = proc.time();
	     ##AFilePrint(paste("ii4 = ", ii4, ", starting with longer EMLARS"));
    USENOISE <- SMS$SigmaSqNoise;	  kLen <- SMS$kLen;  
    rnd <- round(SMS$puse * kLen);  rnd = min( kLen, rnd);
    rnd <- SMS$kActiveLen;
    if (LARSSeekFlag < 0) {
       UseGamma = sqrt(length(SMS$YY)) * sqrt(USENOISE) / rnd / abs(LARSSeekFlag);
    } else {
       UseGamma = sqrt(length(SMS$YY)) * sqrt(USENOISE) / rnd;
    }
      ##sXX = StandardizeXX(SMS$XX);
      sXX = SMS$XX;
      ## sYY = StandardizeXX(SMS$YY);
      sYY = SMS$YY;
      ###sdX = sd(SMS$XX);  ##sdY = sd(SMS$YY);
    
    if (length(sYY) == 0) {
		 return(list( type="MCPLARSeek",
	                BetaFit = rep(0,kLen),
                  BBHatsFit = rep(0,kLen),
	                FitTime =c(0,0,0,0,0),
	                OtherBetas = rep(0,kLen),
	                LARSSeekFlag = LARSSeekFlag,
                  EMLARSObject = NULL, MIPReport = BBHatsFit) );	       
    }
    ##if (length(BetaStart) > 1 && length(BetaStart) == kLen) {
    ##  HitSCAD <- scad.onestep(BetaStart,y=sYY,x = sXX,UseGamma);
     ## HitLasso <- CoordinateDescent(sXX,sYY, OnGamma = UseGamma,
     ##    NumCDOConv = 10000,  CDOEpsilon = .00000001, OnBeta = BetaStart); 
    ##} else {
    ##  HitSCAD <- scad.onestep(rep(0, length(sXX[1,]) ),y=sYY,x = sXX,UseGamma);
     ## HitLasso <- CoordinateDescent(sXX,sYY, OnGamma = UseGamma,
     ##    NumCDOConv = 10000,  CDOEpsilon = .00000001 );     
    ##}
    
      if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE ) {
       fit2 <- ncvreg(X = sXX, y= sYY, gamma = 3.7, penalty="MCP");
      } else {
       fit2 <- ncvreg(X = sXX, y= sYY, gamma = 3.7, penalty="MCP", 
         family = "binomial");      
      }
      if (DoConvexMin == 1) {
        HitSCAD = fit2$beta[2:length(fit2$beta[,1]), fit2$convex.min];
      } else if (DoConvexMin==0) {
        LMin <- sort( abs(fit2$lambda - UseGamma), index=TRUE)$ix[1];
        HitSCAD = fit2$beta[2:length(fit2$beta[,1]), LMin]    
      } else {
         sigmaSqF2 <- fit2$lambda * 0;
         dfF2 <- fit2$lambda * 0;        
         for (ii in 1:length(fit2$lambda)) {
             BetaNot0 <- fit2$beta[2:length(fit2$beta[,1]),ii];
             Errors <- sYY -fit2$beta[1,ii]- sXX %*% BetaNot0;
             IIX <- (1:length(BetaNot0))[BetaNot0 != 0];
             if (length(IIX) > 0) {
                     BetaD <- BetaNot0[IIX];
                     Pen = rep(0, length(BetaD));
                     Pen[BetaD < fit2$lambda[ii]] <- fit2$lambda[ii];
                     Pen[BetaD >= fit2$lambda[ii] & BetaD < (fit2$gamma * fit2$lambda[ii]) ] <-
                       (fit2$gamma * fit2$lambda[ii] - 
                           BetaD[ BetaD >= fit2$lambda[ii] & BetaD < (fit2$gamma * fit2$lambda[ii]) ])  /
                        (fit2$gamma  - 1);
                     SigmaLambda <-matrix(0, length(BetaD), length(BetaD));
                     diag(SigmaLambda) = ( Pen / abs(BetaD) )
                     DFlambda = sum(diag( sXX[,IIX] %*% solve(
                       t(sXX[,IIX]) %*% sXX[,IIX] +
                        n * SigmaLambda ) %*% t(sXX[,IIX])
                          ) );
             } else {
                  DFlambda = 0;
             }
             dfF2[ii] <-DFlambda
             sigmaSqF2[ii] = sum (Errors^2) / n;
          }
          MinLambdaii <- sort( log(sigmaSqF2) + dfF2 * log(n)/n,index=TRUE)$ix[1];
          HitSCAD = fit2$beta[2:length(fit2$beta[,1]), MinLambdaii]    
      }
    HitSCAD = as.numeric(HitSCAD);
    if (length(HitSCAD) != length(sXX[1,]) || is.null(HitSCAD) ||
        any(!is.finite(HitSCAD)) ) {
      BBHatsFit = sXX[1,] * 0;
      HitSCAD = BBHatsFit;    
    }  else {
      BBHatsFit = HitSCAD;
      BBHatsFit[HitSCAD != 0] = 1; 
    } 
      
    if (sum(BBHatsFit) >= 1 && !is.null(SMS$LogitRegression) && 
      (SMS$LogitRegression == TRUE || SMS$LogitRegression == 1)) {
      glm2 = NULL;
      try(glm2 <- glm(SMS$YY~SMS$XX[,BBHatsFit == 1], family=binomial(link="logit")));
      if (!is.null(glm2)) { BetaSM = coefficients(glm2); }
      if (!is.null(coefficients(glm2))) {
        BetaSM = coefficients(glm2);
        BetaSM = BetaSM[2:length(BetaSM)];
      }  else {
        AFilePrint("That BetaSM is Null, so we're in trouble")
      }  
      BetaFit <- rep(0, NCOL(SMS$XX));
      BetaFit[BBHatsFit == 1] <- BetaSM;
	} else  if (sum(BBHatsFit) >= 1) {
		             XXSM <- SMS$XX[, BBHatsFit == 1];
		             XXSXXST <- t(XXSM) %*% XXSM;
		             if (any(is.null(XXSXXST)) || length(dim(XXSXXST)) != 2 
		                || dim(XXSXXST)[1] != dim(XXSXXST)[2] || 
		                dim(XXSXXST)[1] < 1 || any(is.nan(XXSXXST))) { 
		                BTF1 <- HitSCAD[BBHatsFit == 1];                
		             } else if (length(XXSXXST) == 1) {
		                BTF1 <- 1 / XXSXXST * sum(XXSM * SMS$YY);
		             } else if (det(XXSXXST) <= .00000001) {
 		                library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                     ##   "Determinant is Samll, Going with corpcor", sep=""));
		                BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, passed corpcor", sep=""));                                  
		             } else if (!is.null(dim(XXSXXST)) && dim(XXSXXST)[1] > 1 && 
		                        dim(XXSXXST)[1] == dim(XXSXXST)[2] && det(XXSXXST) <= 0.00000001 ) {
 		                library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, Going with corpcor", sep=""));
		                BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, passed corpcor", sep=""));   
		             } else if (!is.null(dim(XXSXXST)) && det(XXSXXST) > .0001 ) {
		                LS <- svd(XXSXXST);
		                if (max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
                         BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
                    }   else {
		                     BTF1 <- try(solve(XXSXXST) %*% t(XXSM) %*% SMS$YY);
		                }
		             }  else {
 		                library(corpcor,warn.conflicts=FALSE, quietly=TRUE);
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, Going with corpcor", sep=""));
		                BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		                ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
                    ##    "Determinant is Samll, passed corpcor", sep=""));                   
                 }
      BetaFit = BBHatsFit;
      BetaFit[BBHatsFit == 1] = BTF1;     
	  } else {
	    BetaFit = BBHatsFit * 0;
	  }
	   T2 <- proc.time();
	   StoddenFit1T <- T2 - T1;
		 return(list( type="MCP",
	                BetaFit = BetaFit,
                  BBHatsFit = BBHatsFit,
	                FitTime =StoddenFit1T,
	                OtherBetas = HitSCAD,
                  EMLARSObject = NULL, MIPReport = BBHatsFit) );	   
}


GenerateTwoLassoSpread <- function(SMS, DoCI = DefaultDoCI,...) {
    library(TwoLassoCpp);
    USENOISE <- SMS$SigmaSqNoise; kLen <- SMS$kLen;
    AFilePrint("0001AAEfficientSimulator:  Running 2 Lasso Spread.")
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
	      LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) * sqrt( 2 * SDA^2 / USENOISE)
        LambdaASeq = c( 1.1 * LambdaAS1, .05 * LambdaAS1);
        if (LambdaASeq[1] >= LambdaDSeq[1]) {
           LambdaASeq[1] = LambdaDSeq[1] + LambdaASeq[1];
           LambdaDSeq[1] = LambdaASeq[1] - LambdaDSeq[1];
           LambdaASeq[1] = LambdaASeq[1] - LambdaDSeq[1];
        }
        LambdaASeq = c(LambdaASeq[1], LambdaASeq);
        LambdaDSeq = c(LambdaDSeq[1], LambdaDSeq);        
	      OrderSeq = c(20, 4,1);
	  if (SMS$LogitRegression == FALSE) { 
     HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = FALSE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF);   
	   ## HitBatWingFast <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, 
	   ##   FixKa = -100, sigmaNoiseSq = USENOISE,
	   ##   RatWant = .2, StlambdaD=LambdaDSeq[1], 
	   ##   StlambdaA = LambdaASeq[1], LambdaDK = LambdaDSeq,
     ##   LambdaAK = LambdaASeq, OrderSeq = OrderSeq, 
     ##   TotalRuns = 5, NumEMConv = 4, MultEMCons = .99, NumCDOConv = 40, 
     ##   CDOEpsilon = .000001, TDFNu = SMS$tNoiseDF);
    } else if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
      ##LambdaAK = StlambdaA * lambdaAMultC^(0:4);
      ##LambdaDK = StlambdaD * lambdaDMultC^(0:4);
      ##OrderSeq = rep(10,5);
      
      ##HitBatWingFast = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
      ##  StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
      ##  LambdaAK = LambdaASeq,LambdaDK = LambdaDSeq,OrderSeq = OrderSeq,
      ##  StLambdaA = LambdaASeq[1], StLambdaD = LambdaDSeq[1], 
      ##  MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
      ##  TotalRuns = 5, MaxLambdaASteps = 50,
      ##  RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
      ##  MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = -1,
      ##  LogitNoise=LogitNoise) 
      if (SMS$n < SMS$p) {
        LogitNoise = SMS$n/SMS$p;
      } else { LogitNoise = 1.0; } 
      HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = LogitNoise, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = c(5,1,1), Verbose = 0, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = -1, DoLogit=TRUE, TargetMinimumBeta = .75, LogitNoise = n/p);  
    }

	  if (!is.null(HitBatWingFast$FailureFlag) && 
      length(HitBatWingFast$FailureFlag) >= 1 &&
      HitBatWingFast$FailureFlag == 1 ) {
       eval(parse(text=GetG0Text("OnFails", "globalenv()", S=1)));
       eval(parse(text=GetG0Text("FailToFitSim", "globalenv()", S=1)));
       eval(parse(text=GetG0Text("FailToFitEMOb", "globalenv()", S=1)));
	     AFilePrint("HitBatWingFast Failed "); 
       OnFails <- OnFails+1;
	     FailToFitSim[[OnFails]] <- SMS;
	     FailToFitEMOb[[OnFails]] <-HitBatWingFast;
       eval(parse(text=SetGText("OnFails", "globalenv()", S=1)));
       eval(parse(text=SetGText("FailToFitSim", "globalenv()", S=1)));
       eval(parse(text=SetGText("FailToFitEMOb", "globalenv()", S=1)));
	  }  	
    if (length(HitBatWingFast$FinalBeta) > 1) {
	  LarsFitBBsFast <- as.vector(HitBatWingFast$FinalBeta);
	  LarsFitBBsFast[abs(HitBatWingFast$FinalBeta) <= MinMin] = 0;
	  LarsFitBBsFast[abs(HitBatWingFast$FinalBeta) > MinMin] = 1;  
	  } else {
	    LarsFitBBsFast = SMS$BetasReal * 0 - 10;
	  }
	  T2 <- proc.time();
	     EMLARSFastT <- T2 - T1; 
    MIPReport <- NULL;
    if (!is.null(HitBatWingFast$TLS$RecBBOn1)) {
      try(MIPReport <- HitBatWingFast$TLS$RecBBOn1[,2]);
    } else {
      try(MIPReport <- HitBatWingFast$TLS$BBOn1);
    }
    CIEst = NULL;  CITime = NULL;  
    if (DoCI == TRUE) {
      CITime1 <- proc.time();
     if (is.null(SMS$LogitRegresson) || SMS$LogitRegression == FALSE) { 
       HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = 0, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, DoLogit=FALSE); 
      } else {
        if (SMS$n < SMS$p) {
          LogitNoise = SMS$n/SMS$p;
        } else { LogitNoise = 1.0; } 
          HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
            SigmaSq = LogitNoise, LambdaAK = LambdaASeq,
            LambdaDK = LambdaDSeq, OrderSeq = c(5,1,1), Verbose = 0, 
            InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
            TDFNu = -1, DoLogit=TRUE, TargetMinimumBeta = .75, LogitNoise = n/p);          
      }
      MIPReport <- NULL;
      if (!is.null(HitBatWingFast$TLS$RecBBOn1)) {
        try(MIPReport <- HitBatWingFast$TLS$RecBBOn1[,2]);
      } else {
        try(MIPReport <- HitBatWingFast$TLS$BBOn1);
      }
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      HitBatWingFast$TLS$ConfidenceQuantiles <- GetConfidenceIntervals;
      HitBatWingFast$TLS$LambdaIndexForConfidenceIntervals <- length(LambdaDSeq)-2
      CITime2 <- proc.time();  CITime <- CITime2-CITime1;
      CIEst <- list();
      MIPMIP <-  HitBatWingFast$TLS$RecBBOn1[, length(LambdaDSeq)-1];
      CIEst[[1]]<- HitBatWingFast$TLS$ConfidenceMatrix;
      CIEst[[2]]<- HitBatWingFast$TLS$UnshrunkConfidenceMatrix;
      CIEst[[3]] <- HitBatWingFast$TLS$HPDMatrix;
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix");
      try(CIEst <- CleanCIEst(CIEst));
      CIQuantiles <- HitBatWingFast$TLS$ConfidenceQuantiles;
    } else {
      CIEst <- NULL;   CIQuantiles  <- NULL;  CITime = NULL;
      MIP = MIPMIP;
    }                                                                                                             
		 return(list( type="TwoLassoSpread",
	                BetaFit = HitBatWingFast$FinalBeta,
                  BBHatsFit = LarsFitBBsFast,
	                FitTime =EMLARSFastT,
	                OtherBetas = HitBatWingFast$TLS$Beta,
                  TwoLassoObject = HitBatWingFast,CIEst=CIEst, 
                  CIQuantiles = CIQuantiles , CITime = CITime,
                  MIP = MIPMIP, MIPReport = MIPReport) );
}

GenerateR2Lasso <- function(SMS, R = .5, DoCI = DefaultDoCI, ...) {
  library(TwoLassoCpp);
  USENOISE <- SMS$SigmaSqNoise; kLen <- SMS$kLen;
	T1 <- proc.time();
	SDD <- sqrt(min(apply(SMS$XX,2,sd))^2 * (length(SMS$XX[,1])-1) );
	LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD, 
    sqrt(3.7 / USENOISE) * SDD *  max( (2*length(SMS$XX[,1])), 200));
  ## LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD
  ##                  , 
  ##                sqrt(3.7 / USENOISE) * SDD *  200
  ##                    );                           
	SABS <- min( abs(SMS$BetasReal[SMS$BetasReal != 0])  );
	SDA <- SDD / sqrt(length(SMS$XX[,1])-1);
	LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) * sqrt( 2 * SDA^2 / USENOISE)
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
    SigmaVec = c(USENOISE * NN^R, USENOISE, USENOISE);
    LambdaASeq = c(LambdaASeq[1] * NN^(.5 * R), LambdaASeq[1], 
      LambdaASeq[length(LambdaASeq)]);
    LambdaDSeq = c(LambdaDSeq[1] * NN^(-.5 * R), LambdaDSeq[1],
      LambdaDSeq[length(LambdaDSeq)] );
  } else {
    SigmaVec = -1;
  }   		      
	if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
    HitBatWingFast <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, 
	    FixKa = -100, sigmaNoiseSq = USENOISE, RatWant = .2, 
      StlambdaD=LambdaDSeq[1], StlambdaA = LambdaASeq[1], LambdaDK = LambdaDSeq,
      LambdaAK = LambdaASeq, OrderSeq = OrderSeq, TotalRuns = 5, 
	    NumEMConv = 4, MultEMCons = .99, NumCDOConv = 40, CDOEpsilon = .000001, 
      TDFNu = SMS$tNoiseDF, SigmaVec= SigmaVec);
  } else if (!is.null(SMS$LogitRegression) && SMS$LogitRegression == TRUE) {
    ##LambdaAK = StlambdaA * lambdaAMultC^(0:4);
    ##  OrderSeq = rep(10,5);
    ##HitBatWingFast = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
    ##  StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
    ##  LambdaAK = LambdaASeq,LambdaDK = LambdaDSeq,OrderSeq = OrderSeq,
    ##  StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
    ##  MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
    ##  TotalRuns = 5, MaxLambdaASteps = 50,
    ##  RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
    ##  MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = -1);
      if (SMS$n < SMS$p) {
        LogitNoise = SMS$n/SMS$p;
      } else { LogitNoise = 1.0; } 
      HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = LogitNoise, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = c(5,1,1), Verbose = 0, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = -1, DoLogit=TRUE, TargetMinimumBeta = .75, LogitNoise = n/p);  
  }
	  if (!is.null(HitBatWingFast$FailureFlag) &&
      length(HitBatWingFast$FailureFlag) >= 1 &&
      HitBatWingFast$FailureFlag ==1 ) {
	     AFilePrint("HitBatWingFast Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
	     FailToFitEMOb[[OnFails]] <<-HitBatWingFast;
	  }  	
    if (!is.null(HitBatWingFast$BBHatsAll) &&
      length(dim(HitBatWingFast$BBHatsAll)) == 2) {
      MIPReport = HitBatWingFast$BBHatsAll[2,];
    } else {
      MIPReport = HitBatWingFast$BBHatsFinal
    }
    if (length(HitBatWingFast$ReturnBetas) > 1) {
	  LarsFitBBsFast <- as.vector(HitBatWingFast$ReturnBetas);
	  LarsFitBBsFast[abs(HitBatWingFast$ReturnBetas) <= MinMin] = 0;
	  LarsFitBBsFast[abs(HitBatWingFast$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBsFast = SMS$BetasReal * 0 - 10;
	  }
	  T2 <- proc.time();
	     EMLARSFastT <- T2 - T1;                      
		 return(list( type="TwoLassoSpread",
	                BetaFit = HitBatWingFast$ReturnBetas,
                  BBHatsFit = LarsFitBBsFast,
	                FitTime =EMLARSFastT,
	                OtherBetas = HitBatWingFast$ReturnBetas,
                  EMLARSObject = HitBatWingFast, MIPReport=MIPReport) );
}

GenerateBayesVarSel <- function(SMS, DoCI = DefaultDoCI,...) {
  library(BayesVarSel);
MyData <- cbind(SMS$YY, SMS$XX);
colnames(MyData) = c("y", paste("x:", 1:NCOL(SMS$XX), sep=""));
MyData <- as.data.frame(MyData);
t1 = proc.time();
Model.Bvs<- GibbsBvs(formula="y~.", data=MyData, prior.betas = "Robust",
         prior.models = "Constant", n.iter=1500, init.model = "Null",
         n.burnin = 50, time.test = FALSE)
         
MySum <- summary(Model.Bvs)
t2 = proc.time();

MIPMIP <- as.numeric(as.character(MySum[[1]][,1]));
ATI <- (1:length(MIPMIP))[MySum[[1]][,2] == "*"];
BBHatsFit <- rep(0, length(MIPMIP));
BBHatsFit[ATI] <- 1;
MeanVar <- rep(0, length(MIPMIP));
aLM <- NULL;
try(aLM <- lm(SMS$YY~SMS$XX[,ATI]));
if (!is.null(aLM)) {
try(MeanVar[ATI] <- aLM$coef[2:length(aLM$coef)])
}  

  TimeBayesVarSel <- t2 - t1; 
  CIEst = NULL; CITime = NULL; CIQuantiles = NULL;
  MIPReport = MIPMIP;                     
	return(list( type="BayesVarSel",
	    BetaFit = MeanVar,
        BBHatsFit = BBHatsFit,
	    FitTime =TimeBayesVarSel,
	    OtherBetas = MeanVar,
        EMLARSObject = Model.Bvs,
        CIEst=CIEst, CITime=CITime, CIQuantiles=CIQuantiles, MIPReport = MIPReport) );
    
}


GenerateBayesLasso <- function(SMS, DoCI = DefaultDoCI,...) {
  library(monomvn, warn.conflicts=FALSE, quietly=TRUE);
  t1 = proc.time();
  sdY = sd(as.vector(SMS$YY));
  sdX = apply(SMS$XX,2,sd);
  meanX = colMeans(SMS$XX); meanY = mean(SMS$YY);
  mX = t( (t(SMS$XX)- meanX) / sdX);
  mY = (SMS$YY-meanY)/sdY;
  RB = blasso(X = mX, y =mY, T = 1000, thin = NULL, RJ = TRUE, M = NULL,
       beta = NULL, lambda2 = 1, s2 = var(mY-mean(mY)),
       case = "default",
       mprior = 0, rd = NULL, ab = NULL, theta=0, rao.s2 = TRUE,
       normalize = TRUE, verb = 0);
  PointEst = colMeans(RB$beta) * sdY / sdX;
  BBOnEst = PointEst;  BBOnEst[abs(PointEst) > .005] = 1;  BBOnEst[abs(PointEst) < .005] = 0;
  t2 = proc.time();
  if (DoCI == TRUE) {
    library(coda);
    BetaChains <- as.mcmc(sdY * t(RB$beta) / sdX);
    CITime1 = proc.time();
    eval(parse(text=GetG0Text("GetConfidenceIntervals")));
    ASG <- summary(BetaChains, quantiles = GetConfidenceIntervals);
    BetaSymmetricIntervals <- ASG[[2]];
    colnames(BetaSymmetricIntervals) <- GetConfidenceIntervals;
    rownames(BetaSymmetricIntervals) <- paste("Beta", 1:NROW(BetaSymmetricIntervals), sep="");
    BetaHPDIntervals <- BetaSymmetricIntervals * 0;
    if (GetConfidenceIntervals[1] == .5) {
      AK <- 1;  BetaHPDIntervals[,1] <- BetaSymmetricIntervals[,1];
    } else {AK <- 0;}
    NIntervals <- (length(GetConfidenceIntervals) - AK)/2;  AtAk <- AK;
    for (jj in 1:NIntervals) {
      CIP <- GetConfidenceIntervals[AtAk +2] - GetConfidenceIntervals[AtAk +1];
      try(HI <- TwoSimR5:::StableHPDinterval(BetaChains, CIP), silent=TRUE);
      try(BetaHPDIntervals[,AtAk+1:2] <- HI);
      AtAk <- AtAk +2;
    }
    CIEst <- list(BetaSymmetricIntervals=BetaSymmetricIntervals,
      BetaHPDIntervals=BetaHPDIntervals);
    CIQuantiles =  GetConfidenceIntervals;
    CITime2 = proc.time();
    try(CIEst <- CleanCIEst(CIEst));
    CITime = CITime2-t1;
  } else {CIEst = NULL; CIQuantiles = NULL;  CITime = NULL; }
  
  TimeBayesLasso <- t2 - t1;                      
	return(list( type="BayesLasso",
	                BetaFit = PointEst,
                  BBHatsFit = BBOnEst,
	                FitTime =TimeBayesLasso,
	                OtherBetas = PointEst,
                  EMLARSObject = RB,
                  CIEst=CIEst, CITime=CITime, CIQuantiles=CIQuantiles, MIPReport = BBOnEst) );
                  
  ##  A Test Code:                 
    X = matrix(rnorm(10*100),100,10); Y = as.vector(X %*% c(1,1,0,0,0,0,0,0,0,-1) + .5*rnorm(100));
    sdX = apply(X,2,sd); sdY = sd(as.vector(Y));
    mX =  (X-colMeans(X))/ apply(X,2,sd);  mY = (Y-mean(Y))/sd(as.vector(Y));
    RB = blasso(X = mX, y = mY, T = 1000, thin = NULL, RJ = TRUE, M = NULL,
       beta = NULL, lambda2 = 1, s2 = var(mY-mean(mY)),
       case = c("default"),
       mprior = 0, rd = NULL, ab = NULL, theta=0, rao.s2 = TRUE,
       normalize = TRUE, verb = 1)
   PointEst = sdY       
       
}


Generate2LassoCV <- function(SMS, DoCI = DefaultDoCI,...) {
  library(TwoLassoCpp);
t1 = proc.time();
 SABS <- min( abs(SMS$BetasReal[SMS$BetasReal != 0])  );
 if (SMS$LogitRegression == TRUE) {
    DoLogit = TRUE;
 } else {
    DoLogit = FALSE;
 }
 CFoldOut <- CFoldValidateTwoLassoCpp(SMS$XX,SMS$YY, DoLogit=DoLogit, nPi=20, nSigma=1, cFolds =5,
   OrderSeq = c(20,4,1), SABS = .05, MaxCauchy = 100, CauchyEpsilon = .00001,
   Verbose = 2)
   if (length(CFoldOut$ReturnBetas) > 1) {
     ReturnBetas = CFoldOut$ReturnBetas;
   } else if (length(CFoldOut$FinalBeta) > 1) {
     ReturnBetas = CFoldOut$FinalBeta;   
   }
   BBOnEst = CFoldOut$BBOn1;  BBOnEst[ BBOnEst < .005]  = 0; 
   BBOnEst[BBOnEst > .005] = 1;
t2 = proc.time();
  if (!is.null(CFoldOut$RecBBOn1) && length(dim(CFoldOut$RecBBOn1)) == 2 &&
    NROW(CFoldOut$RecBBOn1) >= 1 && NCOL(CFoldOut$RecBBOn1) >= 1) {
    MIPReport = CFoldOut$RecBBOn1[,2];
  } else {
    MIPReport = CFoldOut$BBOn1;
  }

    if (DoCI == TRUE) {
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      if (!is.null(CFoldOut$TLS)) {
      CITime1 <- proc.time();
      ## try(CFoldOut$TLS$ConfidenceQuantiles <- c(.5, .25, .75, .05,.95, .025, .975, .01, .99, .005, .995));
      try(CFoldOut$TLS$ConfidenceQuantiles <- GetConfidenceIntervals);
      try(CFoldOut$TLS$LambdaIndexForConfidenceIntervals <- length(CFoldOut$TLS$LambdaDK)-1);
      CIEst <- list();
      try(CIEst[[1]]<- CFoldOut$TLS$ConfidenceMatrix);
      try(CIEst[[2]]<- CFoldOut$TLS$UnshrunkConfidenceMatrix);
      try(CIEst[[3]] <- CFoldOut$TLS$HPDMatrix);
      try(names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix"));
      CITime2 <- proc.time();
      try(CIQuantiles<- CFoldOut$TLS$ConfidenceQuantiles);      
      } else{
      CITime1 <- proc.time();
      ## try(CFoldOut$ConfidenceQuantiles <- c(.5, .25, .75, .05,.95, .025, .975, .01, .99, .005, .995));
      try(CFoldOut$ConfidenceQuantiles <- GetConfidenceIntervals);
      try(CFoldOut$LambdaIndexForConfidenceIntervals <- length(CFoldOut$LambdaDK)-1);
      CIEst <- list();
      try(CIEst[[1]]<- CFoldOut$ConfidenceMatrix);
      try(CIEst[[2]]<- CFoldOut$UnshrunkConfidenceMatrix);
      CITime2 <- proc.time();
      try(names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix"));
      try(CIQuantiles<- CFoldOut$ConfidenceQuantiles);
      }
      CITime <- CITime2 - CITime1;
      CIEst <- CleanCIEst(CIEst);
    } else {
      CIEst <- NULL;  CITime <- NULL;  CIQuantiles <- NULL;
      CITime1 = NULL;  CITime2 = NULL;
    }     
  
  Time2LassoCV <- t2 - t1;                      
	return(list( type="TwoLassoCV",
	                BetaFit = ReturnBetas,
                  BBHatsFit = BBOnEst,
	                FitTime =Time2LassoCV,
	                OtherBetas = ReturnBetas,
                  EMLARSObject = CFoldOut,
                  CIEst = CIEst, CITime = CITime, CIQuantiles=CIQuantiles, 
                  CITime1 = CITime1, CITime2=CITime2, MIPReport = MIPReport, DoLogit=DoLogit) );
       
}

CleanCIEst <- function(CIEst) {
  for (ii in 1:length(CIEst)) {
    if (any(CIEst[[ii]] == -999)) {
      for (jj in 1:NROW(CIEst[[ii]])) {
        if (any(CIEst[[ii]][jj,] == -999)) {
          for (kk in 1:NCOL(CIEst[[ii]])) {
            if (CIEst[[ii]][jj,kk] == -999) {
              if ( as.numeric(colnames(CIEst[[ii]])[kk]) > .5) {
                ABD = CIEst[[ii]][jj, as.numeric(colnames(CIEst[[ii]])) >= .5];
                ABD = ABD[ABD != -999];
                ABD = ABD[!is.na(ABD)];
                if (length(ABD) >= 1) {
                  ATT <- max(ABD);
                  CIEst[[ii]][jj,kk] <- ATT;
                } else {
                  CIEst[[ii]][jj,kk] <- 0;
                }
              } else if (as.numeric(colnames(CIEst[[ii]])[kk]) < .5) {
                ABD = CIEst[[ii]][jj, as.numeric(colnames(CIEst[[ii]])) <= .5];
                ABD = ABD[ABD != -999];
                ABD = ABD[!is.na(ABD)];
                if (length(ABD) >= 1) {
                  ATT <- min(ABD);
                  CIEst[[ii]][jj,kk] <- ATT;
                } else {
                  CIEst[[ii]][jj,kk] <- 0;
                }              
              }
            }
          }
        }
      }
    }
  }
  return(CIEst);
}


GenerateB2Lasso <- function(SMS, R = .5, DoCI = DefaultDoCI,...) {
  library(TwoLassoCpp);
  USENOISE <- SMS$SigmaSqNoise; kLen <- SMS$kLen;
  UseNoise = USENOISE;
	T1 <- proc.time();
	SDD <- sqrt(min(apply(SMS$XX,2,sd))^2 * (length(SMS$XX[,1])-1) );
	p = NCOL(SMS$XX);
  LambdaA0 = sqrt(2+.01 / UseNoise);
  LambdaD0 = 2 / sqrt(UseNoise);
  n = length(SMS$YY);
  LambdaDSeq = c( sqrt(n) * LambdaD0, n^(1.5) * LambdaD0);
  LambdaASeq = c( sqrt(n) * LambdaA0, 1/sqrt(n) * LambdaA0);
  OrderSeq = c(20,1);
  InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)), 5);
      
	if (SMS$LogitRegression == FALSE) {
	  HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
      SigmaSq = USENOISE, LambdaAK = LambdaASeq,
      LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
      InitKKs = InitKKs, RecordFlag = FALSE, HoldOn = FALSE, 
      TDFNu = SMS$tNoiseDF);
    ##HitBatWingFast <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, 
	  ##  FixKa = -100, sigmaNoiseSq = USENOISE, RatWant = .2, 
    ##  StlambdaD=LambdaDSeq[1], StlambdaA = LambdaASeq[1], LambdaDK = LambdaDSeq,
    ##  LambdaAK = LambdaASeq, OrderSeq = OrderSeq  , TotalRuns = 5, 
	  ##  NumEMConv = 4, MultEMCons = .99, NumCDOConv = 40, CDOEpsilon = .000001, 
    ##  TDFNu = SMS$tNoiseDF, SigmaVec= SigmaVec);
  } else if (SMS$LogitRegression == TRUE) {
    ##LambdaAK = StlambdaA * lambdaAMultC^(0:4);
    ##  OrderSeq = rep(10,5);
    if (FALSE) {
    HitBatWingFast = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
      StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
      LambdaAK = LambdaASeq,LambdaDK = LambdaDSeq,OrderSeq = OrderSeq,
      StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
      MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
      TotalRuns = 5, MaxLambdaASteps = 50,
      RecordBetaFlag = TRUE, InitKKs=max(abs(SMS$puse * NCOL(SMS$XX)), 5), PrintFlag = 0,
      MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = -1);
    }
    ##HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
    ##  SigmaSq = USENOISE, LambdaAK = LambdaASeq,
    ##  LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
    ##  InitKKs = InitKKs, RecordFlag = FALSE, HoldOn = FALSE, 
    ##  TDFNu = SMS$tNoiseDF, DoLogit=TRUE);
     if (SMS$n < SMS$p) {
        LogitNoise = SMS$n/SMS$p;
      } else { LogitNoise = 1.0; } 
      HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = LogitNoise, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = c(5,1,1), Verbose = 0, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = -1, DoLogit=TRUE, TargetMinimumBeta = .75, LogitNoise = n/p);  
  }
	  if (!is.null(HitBatWingFast$FailureFlag) && 
      length(HitBatWingFast$FailureFlag) >= 1 &&
      HitBatWingFast$FailureFlag ==1 ) {
	     AFilePrint("HitBatWingFast Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
	     FailToFitEMOb[[OnFails]] <<-HitBatWingFast;
	  } 
    if (length(HitBatWingFast$FinalBeta) > 1) {
	  EstBeta = HitBatWingFast$FinalBeta;
    } else if (length(HitBatWingFast$ReturnBetas) > 1) {
    EstBeta = HitBatWingFast$ReturnBetas; 
	  } else {
  	EstBeta = SMS$BetasReal * 0;
	  } 	
    if (length(HitBatWingFast$FinalBeta) > 1) {
	  LarsFitBBsFast <- as.vector(HitBatWingFast$FinalBeta);
	  LarsFitBBsFast[abs(HitBatWingFast$FinalBeta) <= MinMin] = 0;
	  LarsFitBBsFast[abs(HitBatWingFast$FinalBeta) > MinMin] = 1;  
    } else if (length(HitBatWingFast$ReturnBetas) > 1) {
	  LarsFitBBsFast <- as.vector(HitBatWingFast$ReturnBetas);
	  LarsFitBBsFast[abs(HitBatWingFast$ReturnBetas) <= MinMin] = 0;
	  LarsFitBBsFast[abs(HitBatWingFast$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBsFast = SMS$BetasReal * 0 - 10;
	  }
    MIPReport <- NULL;
    if (!is.null(HitBatWingFast$RecBBOn1) && length(dim(HitBatWingFast$RecBBOn1)) == 2) {
      try(MIPReport <- HitBatWingFast$RecBBOn1[,2]);
    } else {
      try(MIPReport <- HitBatWingFast$BBOn1);
    }
    if (DoCI == TRUE) {
      CITime1 <- proc.time();
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      HitBatWingFast$TLS$ConfidenceQuantiles <- GetConfidenceIntervals;
      HitBatWingFast$TLS$LambdaIndexForConfidenceIntervals <- length(LambdaDSeq)-1;
      CITime2 <- proc.time();
      CIEst <- list();
      CIEst[[1]]<- HitBatWingFast$TLS$ConfidenceMatrix;
      CIEst[[2]]<- HitBatWingFast$TLS$UnshrunkConfidenceMatrix;
      CIEst[[3]] <- HitBatWingFast$TLS$HPDMatrix;
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix");
      CIQuantiles<- HitBatWingFast$TLS$ConfidenceQuantiles;
      CITime <- CITime2-CITime1;
      try(CIEst <- CleanCIEst(CIEst));
    } else {
      CIEst <- NULL;   CIQuantiles<- NULL;  CITime = NULL;
    }  
	  T2 <- proc.time();
	     EMLARSFastT <- T2 - T1;                      
		 return(list( type="B2Lasso",
	                BetaFit = EstBeta,
                  BBHatsFit = LarsFitBBsFast,
	                FitTime =EMLARSFastT,
	                OtherBetas = HitBatWingFast$ReturnBetas,
                  EMLARSObject = HitBatWingFast, CIEst=CIEst, CIQuantiles=CIQuantiles,
                  CITime=CITime, MIPReport = MIPReport) );
}







GenerateTrueBeta <- function(SMS, DoCI = DefaultDoCI,...) {
  BBHatsFit = SMS$BBsReal
   
  if ((is.logical(DoCI) && DoCI == TRUE) ||
    (is.numeric(DoCI) && DoCI >= 1)) {
    try(eval(parse(text=GetG0Text("GetConfidenceIntervals", "globalenv()", S=1))));
    if (!is.null(GetConfidenceIntervals) && GetConfidenceIntervals[1] != 0) {
      ConfidenceQuantiles <- GetConfidenceIntervals; 
    } else {
      ConfidenceQuantiles <- NULL; 
    }
  } else {
    ConfidenceQuantiles <- NULL;
  } 
  CITime <- NULL; CIQuantiles <- NULL; CIEst <- NULL;   
       
  T1 <- proc.time();
	if(sum(BBHatsFit) >= 1) {
		XXSM <- SMS$XX[, BBHatsFit == 1];
		XXSXXST <- t(XXSM) %*% XXSM;
		if (any(is.null(XXSXXST)) || length(dim(XXSXXST)) != 2 
		  || dim(XXSXXST)[1] != dim(XXSXXST)[2] || 
		  dim(XXSXXST)[1] < 1 || any(is.nan(XXSXXST))) { 
		  BTF1 <- HitLasso$FinalBeta[BBHatsFit == 1];                
		} else if (length(XXSXXST) == 1) {
		  BTF1 <- 1 / XXSXXST * sum(XXSM * SMS$YY);
      if (!is.null(ConfidenceQuantiles)) {
        CITime1 <- proc.time();
        sSig <- (sum(SMS$YY^2) - sum( sum(XXSM * SMS$YY) * BTF1 )) / (length(SMS$YY)-2);
        try(CIOne <- BTF1 + sqrt(sSig/as.numeric(XXSXXST)) * qnorm(ConfidenceQuantiles));
        try(names(CIOne) <- ConfidenceQuantiles);
        CIReal <- matrix(0, NCOL(SMS$XX), length(ConfidenceQuantiles));
        try(CIReal[BBHatsFit==1,] <- CIOne);
        try(colnames(CIReal) <- ConfidenceQuantiles);
        try(rownames(CIReal) <- paste("Beta:", 1:NCOL(SMS$XX), sep=""));
        CIEst <- list();
        CIEst[[1]] <- CIReal;  names(CIEst) <- "OracleConfidence"
        CITime2 <- proc.time();
        CITime <- CITime2[1:3] - CITime1[1:3];
        CIQuantiles <- ConfidenceQuantiles;
      }
		} else if (det(XXSXXST) <= .00000001) {
 		  library(corpcor,warn.conflicts=FALSE, quietly=TRUE);
 		##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
    ##   "Determinant is Samll, Going with corpcor", sep=""));
		BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
    ##    "Determinant is Samll, passed corpcor", sep=""));                                  
		} else if (!is.null(dim(XXSXXST)) && dim(XXSXXST)[1] > 1 && 
		  dim(XXSXXST)[1] == dim(XXSXXST)[2] && det(XXSXXST) <= 0.00000001 ) {
 		  library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
 		##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
    ##    "Determinant is Samll, Going with corpcor", sep=""));
		  BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
    ##    "Determinant is Samll, passed corpcor", sep=""));   
		} else if (!is.null(dim(XXSXXST)) && det(XXSXXST) > .0001 ) {
		  LS <- svd(XXSXXST);
		  if (max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
        BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
      } else {
		    try(BTF1 <- solve(XXSXXST) %*% t(XXSM) %*% SMS$YY);
		  }
       CITime1 <- proc.time();
      sSig <- (sum(SMS$YY^2) - (
        as.vector(t(XXSM) %*% SMS$YY) * BTF1)) / (length(SMS$YY) - 1 - NCOL(XXSM));
      sMI <- NULL;
		  if (max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
        sMI <- pseudoinverse( XXSXXST);
      } else {
		    try(sMI <- solve(XXSXXST));
		  }
      try(CIOne <- matrix(rep(BTF1, length(ConfidenceQuantiles)), length(BTF1), 
        length(ConfidenceQuantiles)) + 
        sqrt(as.numeric(sSig) * diag(sMI))  %*% t(qnorm(ConfidenceQuantiles)));
        try(colnames(CIOne) <- ConfidenceQuantiles);
      CIReal <- matrix(0, NCOL(SMS$XX), NCOL(CIOne));
      CIReal[BBHatsFit == 1,] <- CIOne;
      try(rownames(CIReal) <- paste("Beta:", 1:NCOL(SMS$XX), sep=""));
      try(colnames(CIReal) <- ConfidenceQuantiles);
        CIEst <- list();
        CIEst[[1]] <- CIReal;  names(CIEst) <- "OracleConfidence"
        CITime2 <- proc.time();
        CITime <- CITime2[1:3] - CITime1[1:3];
        CIQuantiles <- ConfidenceQuantiles;
    }  else {
 		  library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
 		  ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
      ##    "Determinant is Samll, Going with corpcor", sep=""));
		  BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		  ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
      ##    "Determinant is Samll, passed corpcor", sep=""));                   
    }
    BetaFit = BBHatsFit;
    BetaFit[BBHatsFit == 1] = BTF1;     
	} else {
	  BetaFit = BBHatsFit * 0;
	}
	T2 <- proc.time();
	FitTime <- T2 - T1;
  MIPReport <- BBHatsFit;
	return(list( type="TrueBeta", BetaFit = BetaFit, BBHatsFit = BBHatsFit,
	  FitTime =FitTime, EMLARSObject = NULL, MIPReport = MIPReport,
    CIEst =CIEst, CIQuantiles=CIQuantiles, CITime = CITime) );		     
}


GenerateFailure <- function(SMS) {
  BBHatsFit = SMS$BBsReal  * 0 * -2;
  BetaFit = SMS$BBsReal * 0 - 999;
       
	FitTime <- c(-999,-999,-999)
	return(list( type="LassoFailure", BetaFit = BetaFit, BBHatsFit = BBHatsFit,
	  FitTime =FitTime, EMLARSObject = NULL, MIPReport = NULL, CIEst = NULL, CIQuantiles=NULL, CITime = NULL) );		   
	   
}

GenerateLassoStodden <- function(SMS, LARSSeekFlag, DoCI = DefaultDoCI, ...) {
  library(TwoLassoCpp);
  LARSSeekFlag = 1;
  T1 = proc.time();
	##AFilePrint(paste("ii4 = ", ii4, ", starting with longer EMLARS"));
  USENOISE <- SMS$SigmaSqNoise;	  kLen <- SMS$kLen;  
  rnd <- round(SMS$puse * kLen);  rnd = min( kLen, rnd);
  rnd <- SMS$kActiveLen;
  if (LARSSeekFlag < 0) {
    UseGamma = sqrt(length(SMS$YY)) * sqrt(USENOISE) / rnd * abs(LARSSeekFlag);
  } else {
    UseGamma = sqrt(length(SMS$YY)) * sqrt(USENOISE) / rnd;
  }
  ##sXX = StandardizeXX(SMS$XX);
  sXX = SMS$XX;
  ## sYY = StandardizeXX(SMS$YY);
  sYY = SMS$YY;
  ###sdX = sd(SMS$XX);  ##sdY = sd(SMS$YY);
    
  if (length(sYY) == 0) {
		return(list( type="LassoStodden", BetaFit = rep(0,kLen),
      BBHatsFit = rep(0,kLen), FitTime =c(0,0,0,0,0), OtherBetas = rep(0,kLen),
        EMLARSObject = NULL, MIPReport = rep(0, kLen) ) );	       
  }
  if (!is.null(SMS$LogitRegression)  || SMS$LogitRegression == FALSE) {
  if (!exists("BetaStart")) { BetaStart = 0; }
  if (length(BetaStart) > 1 && length(BetaStart) == kLen) {
    HitLasso <- CoordinateDescent(sXX,sYY, OnGamma = UseGamma,
      NumCDOConv = 10000,  CDOEpsilon = .00000001, OnBeta = BetaStart); 
  } else {
    HitLasso <- CoordinateDescent(sXX,sYY, OnGamma = UseGamma,
      NumCDOConv = 10000,  CDOEpsilon = .00000001 );     
  }
  } else if (SMS$LogitRegression == TRUE) {
    HitLasso <- glmnet(sXX, sYY,
      family="binomial", weights, offset=NULL,
      alpha = 1, ##nlambda = 100, lambda.min = ifelse(nobs<nvars,0.05,0.0001), 
      lambda=UseGamma*2,
      standardize = TRUE, thresh = 1e-04,  dfmax = nvars + 1, 
      pmax = min(dfmax * 1.2, nvars),
      exclude, penalty.factor = rep(1, nvars), maxit=100, 
      HessianExact = FALSE, type = c("covariance", "naive"))
    HitLasso$beta = HitLasso$b.predictor;
  }
  BBHatsFit = HitLasso$FinalBeta;
  BBHatsFit[HitLasso$FinalBeta != 0] = 1;  
      
	if (sum(BBHatsFit) >= 1) {
		XXSM <- SMS$XX[, BBHatsFit == 1];
		XXSXXST <- t(XXSM) %*% XXSM;
		if (any(is.null(XXSXXST)) || length(dim(XXSXXST)) != 2 
		  || dim(XXSXXST)[1] != dim(XXSXXST)[2] || 
		  dim(XXSXXST)[1] < 1 || any(is.nan(XXSXXST))) { 
		  BTF1 <- HitLasso$FinalBeta[BBHatsFit == 1];                
		} else if (length(XXSXXST) == 1) {
		  BTF1 <- 1 / XXSXXST * sum(XXSM * SMS$YY);
		} else if (det(XXSXXST) <= .00000001) {
 		  library(corpcor,warn.conflicts=FALSE, quietly=TRUE);
 		  ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
      ##   "Determinant is Samll, Going with corpcor", sep=""));
		  BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		  ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
      ##    "Determinant is Samll, passed corpcor", sep=""));                                  
		} else if (!is.null(dim(XXSXXST)) && dim(XXSXXST)[1] > 1 && 
		  dim(XXSXXST)[1] == dim(XXSXXST)[2] && det(XXSXXST) <= 0.00000001 ) {
 		  library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
 		  ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
      ##    "Determinant is Samll, Going with corpcor", sep=""));
		  BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		  ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
      ##    "Determinant is Samll, passed corpcor", sep=""));   
		} else if (!is.null(dim(XXSXXST)) && det(XXSXXST) > .0001 ) {
		  LS <- svd(XXSXXST);
		  if (max(LS$d) < 0  || min(LS$d) / max(LS$d) < .000001) {
        BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
      } else {
		    BTF1 <- try(solve(XXSXXST) %*% t(XXSM) %*% SMS$YY);
		  }
		} else {
 		  library(corpcor, warn.conflicts=FALSE, quietly=TRUE);
 		  ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
      ##    "Determinant is Samll, Going with corpcor", sep=""));
		  BTF1 <- pseudoinverse( XXSXXST) %*% t(XXSM) %*% SMS$YY;  
 		  ##AFilePrint(paste("On Member ",  min(FileLenS) + njj, 
      ##    "Determinant is Samll, passed corpcor", sep=""));                   
    }
    BetaFit = BBHatsFit;
    BetaFit[BBHatsFit == 1] = BTF1;     
	} else {
	  BetaFit = BBHatsFit * 0;
	}
	T2 <- proc.time();
	StoddenFit1T <- T2 - T1;
	return(list( type="LassoStodden", BetaFit = BetaFit, BBHatsFit = BBHatsFit,
	  FitTime =StoddenFit1T, OtherBetas = HitLasso$FinalBeta,
    MIPReport = BBHatsFit,
    EMLARSObject = HitLasso) );	   
}

###############################################################################
###  GenerateLimitLasso()
###
###    A function for using EM2Lasso to fit data Y~X
###
###    Note the features here in case you believe specifiction is appropriate
###     for your own simulations.
###
###     USENOISE is set to SMS$SigmaSqNoise which is the TRUE variance
###
###     "pi_A" or ppiuse in EM2Lasso is given SMS$puse = kActiveLen / kLen
###       this is the TRUE proportion of active factors in the model.
###    
###      "LARSSeekFlag" if it is set to 0, one chooses to specify
###           StlambdaA = StlambdaD =Stlambda = .5 / (2 * USENOISE)
###      Essentially to use a small value for lambda based upon lack of
###       knowledge of the distribution.
###
###     If "LARSSeekFlag >= 1", then one first uses a LARS algorithm to find
###       a "lambda" parameter that contains kActiveLen * LARSSeekFlag active
###       factors.  So one stars with a fraction (but larger than the fraction
###       one desires) of the covariates active.
###
####    Smoothing constants for limit lasso here are 2^(1/8), 2^(-1/8), rather
####     slow, but a pain in the but to play with.  TotalRuns=100 so that
####    Final lambdaA and lambdaD values for "Limit-Lasso" estimate are
####      StlambdaA * 2^(-100/8) and StlambdaD * 2^(100/8).
GenerateLimitLasso <- function(SMS, LARSSeekFlag, Ignorant = FALSE, DoCI = DefaultDoCI, ...) {
	     ##AFilePrint(paste("ii4 = ", ii4, ", starting with longer EMLARS"));
  library(TwoLassoCpp);
  USENOISE <- SMS$SigmaSqNoise;	  kLen <- SMS$kLen; 
  MP <- 1;
  if (Ignorant == FALSE) {
    ppiuse = SMS$puse
  } else {
    ppiuse =  1/NCOL(SMS$XX);
    USENOISE = 1;
  }   
	T1 <- proc.time();
	if (LARSSeekFlag == 0) {
	  StlambdaA = .5/ (2*USENOISE);
	  StlambdaD = .5 / (2*USENOISE);
  } else if (LARSSeekFlag >= 1) {     
	  if (LARSSEEKMETHOD == 1)   {
      rnd <- round(SMS$puse * kLen);  rnd = min( kLen, rnd);
      rnd <- SMS$kActiveLen;
      MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)	  
	    ITO <- MyLarsGiveMeClose(MyFit, IntClose = rnd, MinMin = 0.000000001);   
	    StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	    StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);     
    } else {
      rnd1 <- round(SMS$kActiveLen * LARSSeekFlag);  rnd1 = min( kLen, rnd1);
      MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)	  
	    ITO <- MyLarsGiveMeClose(MyFit, IntClose = rnd1, MinMin = 0.000000001);   
	      StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE);
	      StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE);              
    }     
	  ##StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE);
	  ##StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE);
	}	else if (LARSSeekFlag < 0) {
	  rnd <- round(SMS$puse * kLen);  rnd = min( kLen, rnd);
	  rnd <- SMS$kActiveLen;
    StlambdaA = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;
    StlambdaD = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;  
  } 
  
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
    TotalRuns = 100;  lambdaAmultStop = 50;
    LambdaAK = StlambdaA * (2^(-1/8))^(0:TotalRuns);
    LambdaAK[0:TotalRuns > lambdaAmultStop]  <- LambdaAK[lambdaAmultStop];
    LambdaDK = StlambdaD * (2^(1/8))^(0:TotalRuns);
    OrderSeq <- c(8, rep(4, TotalRuns));
    
   HitLL <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = 1, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF); 
       
    ##HitLL <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = ppiuse, 
	  ##  sigmaNoiseSq = USENOISE, RatWant = .2, StlambdaD=StlambdaD, 
	  ##  StlambdaA=StlambdaA, lambdaDMultC = 2^(1/8), lambdaAMultC = 2^(-1/8),
	  ##  lambdaAmultStop = lambdaAmultStop, TotalRuns = TotalRuns, NumEMConv = 4, 
    ##  MultEMCons = .99, NumCDOConv = 40, CDOEpsilon = .000001,
	  ##  FixKa = -100, TDFNu = SMS$tNoiseDF);
	} else if (SMS$LogitRegression == TRUE) {
	  TotalRuns = 100;
    LambdaAK = StlambdaA * (2^(-1/8))^(0:TotalRuns);
    LambdaDK = StlambdaD * (2^(1/8))^(0:TotalRuns);
    OrderSeq = rep(10,TotalRuns);
   HitLL <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = 1, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, DoLogit = TRUE); 
    ##HitLL = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
    ##  StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
    ##  LambdaAK = LambdaAK,LambdaDK = LambdaDK,OrderSeq = OrderSeq,
    ##  StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
    ##  MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
    ##  TotalRuns = 5, MaxLambdaASteps = 50,
    ##  RecordBetaFlag = TRUE, InitKKs=max(abs(SMS$puse * NCOL(SMS$XX)),5), PrintFlag = 0,
    ##  MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = -1)   
  }
	if (!is.null(HitLL$FailureFlag) && length(HitLL$FailureFlag) >= 1 &&
    HitLL$FailureFlag ==1 ) {
	  AFilePrint(paste("HitLL: LimitLasso: normal version, SMS$puse = ",
      SMS$puse, " Failed ", sep="")); 
    OnFails <<- OnFails+1;     
	}  	  
  
    MP <- sort(abs(1/USENOISE - .25*(LambdaDK-LambdaAK)^2), index=TRUE)$ix[1];
    MIPReport <- NULL;
    HitLL$MP <- MP;
    if (!is.null(HitLL$TLS$RecBBOn1) && !is.null(HitLL$TLS$RecBBOn1)) {
      try(MIPReport <- HitLL$TLS$RecBBOn1[,MP]);
    } else if (!is.null(HitLL$BBHatsAll) && length(dim(HitLL$BBHatsAll)) == 2 &&
      NROW(HitLL$BBHatsAll) >= 1) {
      try(MIPReport <- HitLL$BBHatsAll[MP,]);
    } else {
      try(MIPReport <- HitLL$BBHatsFinal);
    }  
      
    if (!is.null(HitLL$FinalBeta) && length(HitLL$FinalBeta) >= 1) {
		  AFitBeta <- as.vector(HitLL$FinalBeta);
	  } else if (length(HitLL$ReturnBetas) > 1) {
		  AFitBeta <- as.vector(HitLL$ReturnBetas);
		}
	                    
    if (!is.null(HitLL$FinalBeta) && length(HitLL$FinalBeta) >= 1) {
		  LarsFitBBs1 <- as.vector(HitLL$FinalBeta);
		  LarsFitBBs1[abs(HitLL$FinalBeta) <= MinMin] = 0;
		  LarsFitBBs1[abs(HitLL$FinalBeta) > MinMin] = 1; 
	  } else if (length(HitLL$ReturnBetas) > 1) {
		  LarsFitBBs1 <- as.vector(HitLL$ReturnBetas);
		  LarsFitBBs1[abs(HitLL$ReturnBetas) <= MinMin] = 0;
		  LarsFitBBs1[abs(HitLL$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBs1 = SMS$BetaReal * 0 - 10;
	  }
	   T2 <- proc.time();
	   LarsFit1T <- T2 - T1;
     if ((is.logical(DoCI) && DoCI == TRUE) ||
      (is.numeric(DoCI) && DoCI >= 1)) {
      CITime1 <- proc.time();
      if (SMS$LogitRegression == FALSE) { 
       HitLL <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
         SigmaSq = USENOISE, LambdaAK = LambdaAK,
         LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
         InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = 1, HoldOn = FALSE, 
         TDFNu = SMS$tNoiseDF, DoLogit=FALSE); 
       
    ##HitLL <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = ppiuse, 
	  ##  sigmaNoiseSq = USENOISE, RatWant = .2, StlambdaD=StlambdaD, 
	  ##  StlambdaA=StlambdaA, lambdaDMultC = 2^(1/8), lambdaAMultC = 2^(-1/8),
	  ##  lambdaAmultStop = lambdaAmultStop, TotalRuns = TotalRuns, NumEMConv = 4, 
    ##  MultEMCons = .99, NumCDOConv = 40, CDOEpsilon = .000001,
	  ##  FixKa = -100, TDFNu = SMS$tNoiseDF);
	} else if (SMS$LogitRegression == TRUE) {
	  TotalRuns = 100;
      LambdaAK = StlambdaA * (2^(-1/8))^(0:TotalRuns);
      LambdaDK = StlambdaD * (2^(1/8))^(0:TotalRuns);
      OrderSeq = rep(10,TotalRuns);
      HitLL <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = 1, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, DoLogit = TRUE); 
    } 
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      HitLL$TLS$ConfidenceQuantiles <- GetConfidenceIntervals;
      HitLL$TLS$LambdaIndexForConfidenceIntervals <- length(HitLL$TLS$LambdaDK)-2
      CITime2 <- proc.time();  CITime <- CITime2-CITime1;
      CIEst <- list();
      CIEst[[1]]<- HitLL$TLS$ConfidenceMatrix;
      CIEst[[2]]<- HitLL$TLS$UnshrunkConfidenceMatrix;
      CIEst[[3]] <- HitLL$TLS$HPDMatrix;
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix");
      try(CIEst <- CleanCIEst(CIEst));
      CIQuantiles <- HitLL$TLS$ConfidenceQuantiles;
    } else {
      CIEst <- NULL;   CIQuantiles  <- NULL;  CITime = NULL;
    }            
		 return(list( type="LimitLasso",
	                BetaFit = AFitBeta,
                  BBHatsFit = LarsFitBBs1,
	                FitTime =LarsFit1T,
	                OtherBetas = HitLL$ReturnBetas,
                  EMLARSObject = HitLL, MIPReport = MIPReport,
                  CIEst=CIEst,CIQuantiles=CIQuantiles,CITime=CITime) );	   
}

###############################################################################
###  GeneratepiALimitLasso()
###
###    A function for using EM2PILasso to fit data Y~X
###
###    Note the features here in case you believe specifiction is appropriate
###     for your own simulations.
###
###    Difference between EM2Lasso/GenerateLimitLasso is that this algorithm
###     uses priors to arrive at estimates for pi_A and SigmaSq rather than
###     accepting a fixed number
###
###    SigmaEta and SigmaBarS are inverse-chisquare input parameters for SigmaSq
###     For simulations, inputing "SigmaBarS = SMS$SigmaSqNoise"
###      , aka, an unbiased prior for the unknown noise distribution seems acceptable
###     SigmaEta is number of degrees of freedom for the prior.  So that it is
###      informative, having it be, say "n/2", where n is the amount of data
###      seems like one acceptable choice.  
###     SigmaEta = 5 or SigmaEta = 10 are two other posibilities that ensure a proper
###      prior without too myuch informative power
###
###    EM2PILasso ordinarily takes two Beta parameters Beta(  m_1, m_2 )
###     This version of the algorithm suggests setting musize = m_1 + m_2
###     so that one has essentially "musize" worth of data, and then has an unbiased
###     prior     m_1 / (m_1 + m_2) =   True pi_A
###    
###      "LARSSeekFlag" if it is set to 0, one chooses to specify
###           StlambdaA = StlambdaD =Stlambda = .5 / (2 * USENOISE)
###      Essentially to use a small value for lambda based upon lack of
###       knowledge of the distribution.
###
###     If "LARSSeekFlag >= 1", then one first uses a LARS algorithm to find
###       a "lambda" parameter that contains kActiveLen * LARSSeekFlag active
###       factors.  So one stars with a fraction (but larger than the fraction
###       one desires) of the covariates active.
###
####    Smoothing constants for limit lasso here are 2^(1/8), 2^(-1/8), rather
####     slow, but a pain in the but to play with.  TotalRuns=100 so that
####    Final lambdaA and lambdaD values for "Limit-Lasso" estimate are
GeneratepiALimitLasso <- function(SMS, LARSSeekFlag, 
  musize, SigmaEta, SigmaBarS, DoCI = DefaultDoCI, ...) {
  library(TwoLassoCpp);
  USENOISE <- SMS$SigmaSqNoise;	  kLen <- SMS$kLen;  
  T1 <- proc.time();	  
	if (LARSSeekFlag == 0) {
	  StlambdaA = .5/ (2*USENOISE);
	  StlambdaD = .5 / (2*USENOISE);
  } else if (LARSSeekFlag >= 1) {
    if (LARSSEEKMETHOD == 1)   {
      rnd <- round(SMS$puse * kLen);  rnd = min( kLen, rnd);
      rnd <- SMS$kActiveLen;
      MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)	  
	    ITO <- MyLarsGiveMeClose(MyFit, IntClose = rnd, MinMin = 0.000000001);   
	    StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	    StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);     
  } else {
    rnd1 <- round(SMS$kActiveLen * LARSSeekFlag);  rnd1 = min( kLen, rnd1);
    MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)	  
	  ITO <- MyLarsGiveMeClose(MyFit, IntClose = rnd1, MinMin = 0.000000001);   
	  StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE);
	  StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE);              
  }     
 } else if (LARSSeekFlag < 0) {
	 rnd <- round(SMS$kActiveLen);  rnd = min( kLen, rnd);
   StlambdaA = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;
   StlambdaD = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;  
 }  
                                      
 if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
   TotalRuns = 100;   lambdaAmultStop = 50;
   LambdaAK = StlambdaA * (2^(-1/8))^(0:TotalRuns);
   LambdaDK = StlambdaD * (2^(1/8))^(0:TotalRuns);
   OrderSeq = c(20, rep(4,TotalRuns-1),20);   
   HitBatWingPI <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = as.integer(1), HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, PiAPrior = c(musize*SMS$puse, musize*(1-SMS$puse)),
       SigmaPrior = c(SigmaBarS, SigmaEta) ); 
	 ##HitBatWingPI <- EMPI2Lasso(xxs=SMS$XX, yys=SMS$YY, pPpiuseSt = SMS$puse, 
	 ##  sigmaNoiseSq = USENOISE, RatWant = .2, StlambdaD= StlambdaD, 
	 ##  StlambdaA=StlambdaA, lambdaDMultC = 2^(1/8), lambdaAMultC = 2^(-1/8),
	 ##  lambdaAmultStop = lambdaAmultStop, TotalRuns = TotalRuns, NumEMConv = 4, MultEMCons = .99,
	 ##  m1 = musize * SMS$puse, m2 = musize * (1-SMS$puse),
	 ##  SigmaEta = SigmaEta, SigmaBarSq = SigmaBarS,
	 ##  StandardFlag = 0, NumCDOConv = 150, CDOEpsilon = .000001, 
   ##  TDFNu = SMS$tNoiseDF);	
 } else if (SMS$LogitRegression == TRUE) {
   LambdaAK = StlambdaA * (2^(-1/8))^(0:TotalRuns);
   LambdaDK = StlambdaD * (2^(1/8))^(0:TotalRuns);
   OrderSeq = c(20, rep(4,TotalRuns-1),20);  
   HitBatWingPI <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = 1, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, PiAPrior = c(musize*SMS$puse, musize*(1-SMS$puse)),
       SigmaPrior = c(SigmaBarS, SigmaEta), DoLogit = TRUE); 
    ##HitBatWingPI = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
    ##  StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
    ##  LambdaAK = LambdaAK,LambdaDK = LambdaDK,OrderSeq = OrderSeq,
    ##  StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
    ##  MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
    ##  TotalRuns = 5, MaxLambdaASteps = 50,
    ##  RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
    ##  MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = 
    ##   c(musize * SMS$puse, musize * (1-SMS$puse))  ); 
 }
 if (!is.null(HitBatWingPI) && !is.null(HitBatWingPI$FailureFlag) &&
   is.numeric(HitBatWingPI$FailureFlag) && length(HitBatWingPI$FailureFlag) == 1 &&
    HitBatWingPI$FailureFlag ==1 ) {
	 AFilePrint("HitBatWingPI Failed "); OnFails <<- OnFails+1;
	 FailToFitSim[[OnFails]] <<- SMS;
   FailToFitEMOb[[OnFails]] <<- HitBatWingPI;	    	     
 }  
 if (length(HitBatWingPI) == 1) {
	 BatWingFitBBsPI = SMS$BetaReal * 0 - 10;
   BetaFit = SMS$BetaReal*0;
 } else if (length(HitBatWingPI$FinalBeta) > 1) {
   BetaFit = HitBatWingPI$FinalBeta;
 	 BatWingFitBBsPI <- as.vector(HitBatWingPI$FinalBeta);
	 BatWingFitBBsPI[abs(HitBatWingPI$FinalBeta) <= MinMin] = 0;
	 BatWingFitBBsPI[abs(HitBatWingPI$FinalBeta) > MinMin] = 1; 
 } else if ( length(HitBatWingPI$ReturnBetas) > 1) {
   BetaFit = as.vector(HitBatWingPI$ReturnBetas);
	 BatWingFitBBsPI <- as.vector(HitBatWingPI$ReturnBetas);
	 BatWingFitBBsPI[abs(HitBatWingPI$ReturnBetas) <= MinMin] = 0;
	 BatWingFitBBsPI[abs(HitBatWingPI$ReturnBetas) > MinMin] = 1;     
 } else {
   BetaFit = SMS$BetaReal*0;
	 BatWingFitBBsPI = SMS$BetaReal * 0 - 999;
 }
	 if ((length(BatWingFitBBsPI)==1 && is.nan(BatWingFitBBsPI)) 
       || is.nan(BatWingFitBBsPI[1])) {
	   BatWingFitBBsPI = (1:length(BatWingFitBBsPI)) * 0 - 999;
	 }
   if (length(HitBatWingPI$FinalBeta) > 1) {
 	   LarsFitBBsPI <- as.vector(HitBatWingPI$FinalBeta);
	   LarsFitBBsPI[abs(HitBatWingPI$FinalBeta) <= MinMin] = 0;
	   LarsFitBBsPI[abs(HitBatWingPI$FinalBeta) > MinMin] = 1; 
	 } else if (length(HitBatWingPI$ReturnBetas) > 1) {
		 LarsFitBBsPI <- as.vector(HitBatWingPI$ReturnBetas);
		 LarsFitBBsPI[abs(HitBatWingPI$ReturnBetas) <= MinMin] = 0;
		 LarsFitBBsPI[abs(HitBatWingPI$ReturnBetas) > MinMin] = 1;  
	 } else {
	   LarsFitBBsPI = SMSReal * 0 - 10;
	 }
	 T2 <- proc.time();
   
    LambdaAK = StlambdaA * (2^(-1/8))^(0:TotalRuns);
    LambdaDK = StlambdaD * (2^(1/8))^(0:TotalRuns);
    LambdaAK[0:TotalRuns >= lambdaAmultStop] <- LambdaAK[lambdaAmultStop];
    MP <- sort(abs(USENOISE-.25*(LambdaDK-LambdaAK)^2), index=TRUE)$ix[1];
   MIPReport <- NULL;
   if (!is.null(HitBatWingPI$RecBBOn1) && length(dim(HitBatWingPI$RecBBOn1)) == 2 &&
     NROW(HitBatWingPI$RecBBOn1) >= 1) {
     try(MIPReport <- HitBatWingPI$RecBBOn1[,MP]);
   } else if (!is.null(HitBatWingPI$BBHatsAll) && length(dim(HitBatWingPI$BBHatsAll)) == 2 &&
     NROW(HitBatWingPI$BBHatsAll) >= 1) {
     try(MIPReport <- HitBatWingPI$BBHatsAll[MP,]);   
   } else if (!is.null(HitBatWingPI$BBHatsFinal)) {
     try(MIPReport <- HitBatWingPI$BBHatsFinal);
   } else if (!is.null(HitBatWingPI$BBOn1)) {
     try(MIPReport <- HitBatWingPI$BBOn1);
   }
   if ((is.logical(DoCI) && DoCI == TRUE) ||
      (is.numeric(DoCI) && DoCI >= 1)) {
      CITime1 <- proc.time();
      OOSeq <- c(8, rep(3, length(LambdaAK)-2),8);
      if (SMS$LogitRegression == FALSE) { 
       HitBatWingPI <-  TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = 1, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior = c(SigmaBarS, SigmaEta), PiAPrior = c(musize*SMS$puse, musize*(1-SMS$puse)), DoLogit=FALSE); 
      } else {
        if (SMS$n < SMS$p) { LogitNoise <- SMS$n / SMS$p; } else {
          LogitNoise = 1.0; 
        }
        HitBatWingPI <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = 1, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, PiAPrior = c(musize*SMS$puse, musize*(1-SMS$puse)),
       SigmaPrior = c(SigmaBarS, SigmaEta), DoLogit = TRUE, LogitNoise=LogitNoise); 
      }
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      HitBatWingPI$TLS$ConfidenceQuantiles <- GetConfidenceIntervals;
      HitBatWingPI$TLS$LambdaIndexForConfidenceIntervals <- length(HitBatWingPI$TLS$LambdaDK)-2
      CITime2 <- proc.time();  CITime <- CITime2-CITime1;
      CIEst <- list();
      CIEst[[1]]<- HitBatWingPI$TLS$ConfidenceMatrix;
      CIEst[[2]]<- HitBatWingPI$TLS$UnshrunkConfidenceMatrix;
      CIEst[[3]] <- HitBatWingPI$TLS$HPDMatrix;
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix");
      try(CIEst <- CleanCIEst(CIEst));
      CIQuantiles <- HitBatWingPI$TLS$ConfidenceQuantiles;
    } else {
      CIEst <- NULL;   CIQuantiles  <- NULL;  CITime = NULL;
    }            
    OrderSeq = c(20, rep(4,98),20,20,20);
	 LarsFitPIT <- T2 - T1;
	 return(list( type="piALimitLasso", BetaFit = BetaFit,
     BBHatsFit = LarsFitBBsPI, FitTime =LarsFitPIT,
	   OtherBetas = HitBatWingPI$ReturnBetas, EMLARSObject = HitBatWingPI,
     MIPReport = MIPReport, CIEst=CIEst, CITime=CITime, CIQuantiles=CIQuantiles) );
}

GeneratepiAM2Lasso <- function(SMS, LARSSeekFlag, 
  musize, SigmaEta, SigmaBarS, R = 0, DoCI = DefaultDoCI, ...) {
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
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
	  HitBatWingPI <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = FALSE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior = SigmaPrior, PiAPrior = PiAPrior);   
  } else if (SMS$LogitRegression == TRUE) {
    if (FALSE) {
    HitBatWingPI = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
      StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
      LambdaAK = LambdaASeq,LambdaDK = LambdaDSeq,OrderSeq = OrderSeq,
      StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
      MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
      TotalRuns = 5, MaxLambdaASteps = 50,
      RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
      MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = 
       c(musize *SMS$puse, musize * (1-SMS$puse)) );
    }
    if (SMS$n < SMS$p) {
      LogitNoise = SMS$n / SMS$p;
    } else {
      LogitNoise = 1.0;
    }
    HitBatWingPI <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = FALSE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior = SigmaPrior, PiAPrior = PiAPrior, DoLogit=TRUE,
       LogitNoise=LogitNoise); 
  }	
   MIPReport <- NULL;
   if (!is.null(HitBatWingPI$RecBBOn1) && length(dim(HitBatWingPI$RecBBOn1)) == 2
     && NCOL(HitBatWingPI$RecBBOn1) >= 1) {
      try(MIPReport <- HitBatWingPI$RecBBOn1[,2]);
   } else if (!is.null(HitBatWingPI$BBHatsAll) && length(dim(HitBatWingPI$BBHatsAll)) == 2
     && NCOL(HitBatWingPI$BBHatsAll) >= 1) {
     try(MIPReport <- HitBatWingPI$BBHatsAll[2,]);
   } else if (!is.null(HitBatWingPI$BBHatsFinal)) {
     try(MIPReport <- HitBatWingPI$BBHatsFinal);
   } else {
     try(MIPReport <- HitBatWingPI$BBOn1);
   }
	  if (!is.null(HitBatWingPI) && !is.null(HitBatWingPI$FailureFlag) &&
            is.numeric(HitBatWingPI$FailureFlag) && length(HitBatWingPI$FailureFlag) > 1 &&
            HitBatWingPI$FailureFlag ==1 ) {
	     AFilePrint("HitBatWingPI Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
         FailToFitEMOb[[OnFails]] <<- HitBatWingPI;	    	     
	  }  
	  if (length(HitBatWingPI) == 1) {
	    BatWingFitBBsPI = SMSReal * 0 - 10;
      } else if ( length(HitBatWingPI$FinalBeta) > 1) {
		  BatWingFitBBsPI <- as.vector(HitBatWingPI$FinalBeta);
		  BatWingFitBBsPI[abs(HitBatWingPI$FinalBeta) <= MinMin] = 0;
		  BatWingFitBBsPI[abs(HitBatWingPI$FinalBeta) > MinMin] = 1;  
	  } else if ( length(HitBatWingPI$ReturnBetas) > 1) {
		  BatWingFitBBsPI <- as.vector(HitBatWingPI$ReturnBetas);
		  BatWingFitBBsPI[abs(HitBatWingPI$ReturnBetas) <= MinMin] = 0;
		  BatWingFitBBsPI[abs(HitBatWingPI$ReturnBetas) > MinMin] = 1;     
	  } else {
	    BatWingFitBBsPI = SMSReal * 0 - 999;
	  }
	  if ((length(BatWingFitBBsPI)==1 && is.nan(BatWingFitBBsPI)) || is.nan(BatWingFitBBsPI[1])) {
	      BatWingFitBBsPI = (1:length(BatWingFitBBsPI)) * 0 - 999;
	  }
    if (length(HitBatWingPI$FinalBeta) > 1) {
		  BetaFit <- as.vector(HitBatWingPI$FinalBeta);
	  } else if (length(HitBatWingPI$ReturnBetas) > 1) {
		  BetaFit <- as.vector(HitBatWingPI$ReturnBetas);	
	  } else {
	    BetaFit = SMSReal * 0 - 10;
	  }
    if (length(HitBatWingPI$FinalBeta) > 1) {
		  LarsFitBBsPI <- as.vector(HitBatWingPI$FinalBeta);
		  LarsFitBBsPI[abs(HitBatWingPI$FinalBeta) <= MinMin] = 0;
		  LarsFitBBsPI[abs(HitBatWingPI$FinalBeta) > MinMin] = 1;  
	  } else if (length(HitBatWingPI$ReturnBetas) > 1) {
		  LarsFitBBsPI <- as.vector(HitBatWingPI$ReturnBetas);
		  LarsFitBBsPI[abs(HitBatWingPI$ReturnBetas) <= MinMin] = 0;
		  LarsFitBBsPI[abs(HitBatWingPI$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBsPI = SMSReal * 0 - 10;
	  }
	  T2 <- proc.time();
	  LarsFitPIT <- T2 - T1;
    if ((is.logical(DoCI) && DoCI == TRUE) ||
      (is.numeric(DoCI) && DoCI >= 1)) {
      CITime1 <- proc.time();
      if (SMS$LogitRegression == FALSE) { 
       HitBatWingPI <-  TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior = SigmaPrior, PiAPrior = PiAPrior, DoLogit=FALSE); 
      }  else  if (SMS$LogitRegression == TRUE) {
       if (SMS$n < SMS$p) {
         LogitNoise <- SMS$n/SMS$p;
       } else {
         LogitNoise <- 1.0;
       }
       HitBatWingPI <-  TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior = SigmaPrior, PiAPrior = PiAPrior, DoLogit=TRUE,
       LogitNoise=LogitNoise); 
      } 
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      HitBatWingPI$TLS$ConfidenceQuantiles <- GetConfidenceIntervals;
      HitBatWingPI$TLS$LambdaIndexForConfidenceIntervals <- length(HitBatWingPI$TLS$LambdaDK)-2
      CITime2 <- proc.time();  CITime <- CITime2-CITime1;
      CIEst <- list();
      CIEst[[1]]<- HitBatWingPI$TLS$ConfidenceMatrix;
      CIEst[[2]]<- HitBatWingPI$TLS$UnshrunkConfidenceMatrix;
      CIEst[[3]] <- HitBatWingPI$TLS$HPDMatrix;
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix");
      try(CIEst <- CleanCIEst(CIEst));
      CIQuantiles <- HitBatWingPI$TLS$ConfidenceQuantiles;
    } else {
      CIEst <- NULL;   CIQuantiles  <- NULL;  CITime = NULL;
      MIP = MIPMIP;
    }            
	    return(list( type="piALimitLasso",
	                BetaFit = BetaFit,
                  BBHatsFit = LarsFitBBsPI,
	                FitTime =LarsFitPIT,
	                OtherBetas = HitBatWingPI$ReturnBetas,
                  EMLARSObject = HitBatWingPI, MIPReport = MIPReport,
                  CIEst=CIEst, CITime=CITime,
                  CIQuantiles=CIQuantiles) );
}

GeneratepiAM2Lasso <- function(SMS, LARSSeekFlag, 
  musize, SigmaEta, SigmaBarS, R = 0, DoCI = DefaultDoCI,...) {
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
  SigmaPrior <- c(SigmaBarS, SigmaEta);
  PiAPrior <- musize * c(SMS$puse, 1- SMS$puse);
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
  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
	  HitBatWingPI <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = FALSE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior = SigmaPrior, PiAPrior = PiAPrior);
  } else if (SMS$LogitRegression == TRUE) {
      if (SMS$n < SMS$p) {
        LogitNoise <- SMS$n / SMS$p;
      } else { LogitNoise = 1.0; }
	  HitBatWingPI <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = FALSE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior = SigmaPrior, PiAPrior = PiAPrior,
       DoLogit = TRUE, LogitNoise=LogitNoise);
    ##HitBatWingPI = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
    #3  StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
    ##  LambdaAK = LambdaASeq,LambdaDK = LambdaDSeq,OrderSeq = OrderSeq,
    ##  StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
    ##  MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
    ##  TotalRuns = 5, MaxLambdaASteps = 50,
    ##  RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
    ##  MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = -1, PriorPi = 
    ##   c(musize *SMS$puse, musize * (1-SMS$puse)) );     
  }	
  MIPReport <- NULL;
   if (!is.null(HitBatWingPI$RecBBOn1) && length(dim(HitBatWingPI$RecBBOn1)) == 2 &&
     NCOL(HitBatWingPI$RecBBOn1) >=  1) {
     try(MIPReport <- HitBatWingPI$RecBBOn1[,2]);
   } else if (!is.null(HitBatWingPI$BBHatsAll) && length(dim(HitBatWingPI$BBHatsAll)) == 2 &&
     NROW(HitBatWingPI$BBHatsAll) >=  1) {
     try(MIPReport <- HitBatWingPI$BBHatsAll[2,]);
   } else if (!is.null(HitBatWingPI$BBHatsFinal)) {
     try(MIPReport <- HitBatWingPI$BBHatsFinal);
   } else {
     try(MIPReport <- HitBatWingPI$BBOn1)
   }
	  if (!is.null(HitBatWingPI) && !is.null(HitBatWingPI$FailureFlag) &&
            is.numeric(HitBatWingPI$FailureFlag) && length(HitBatWingPI$FailureFlag) >= 1 &&
            HitBatWingPI$FailureFlag ==1 ) {
	     AFilePrint("HitBatWingPI Failed "); OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
         FailToFitEMOb[[OnFails]] <<- HitBatWingPI;	    	     
	  }  
	  if (length(HitBatWingPI) == 1) {
	    BatWingFitBBsPI = SMSReal * 0 - 10;
	  } else if ( length(HitBatWingPI$FinalBeta) > 1) {
		  BatWingFitBBsPI <- as.vector(HitBatWingPI$FinalBeta);
		  BatWingFitBBsPI[abs(HitBatWingPI$FinalBeta) <= MinMin] = 0;
		  BatWingFitBBsPI[abs(HitBatWingPI$FinalBeta) > MinMin] = 1;     
	  } else if ( length(HitBatWingPI$ReturnBetas) > 1) {
		  BatWingFitBBsPI <- as.vector(HitBatWingPI$ReturnBetas);
		  BatWingFitBBsPI[abs(HitBatWingPI$ReturnBetas) <= MinMin] = 0;
		  BatWingFitBBsPI[abs(HitBatWingPI$ReturnBetas) > MinMin] = 1;     
	  } else {
	    BatWingFitBBsPI = SMSReal * 0 - 999;
	  }
	  if ((length(BatWingFitBBsPI)==1 && is.nan(BatWingFitBBsPI)) || is.nan(BatWingFitBBsPI[1])) {
	      BatWingFitBBsPI = (1:length(BatWingFitBBsPI)) * 0 - 999;
	  }
	  if (length(HitBatWingPI$FinalBeta) > 1) {
		  LarsFitBBsPI <- as.vector(HitBatWingPI$FinalBeta);
		  LarsFitBBsPI[abs(HitBatWingPI$FinalBeta) <= MinMin] = 0;
		  LarsFitBBsPI[abs(HitBatWingPI$FinalBeta) > MinMin] = 1;  
	  } else if (length(HitBatWingPI$ReturnBetas) > 1) {
		  LarsFitBBsPI <- as.vector(HitBatWingPI$ReturnBetas);
		  LarsFitBBsPI[abs(HitBatWingPI$ReturnBetas) <= MinMin] = 0;
		  LarsFitBBsPI[abs(HitBatWingPI$ReturnBetas) > MinMin] = 1;  
	  } else {
	    LarsFitBBsPI = SMSReal * 0 - 10;
	  }
	  if (length(HitBatWingPI$FinalBeta) > 1) {
		  BetaFit <- as.vector(HitBatWingPI$FinalBeta);
	  } else if (length(HitBatWingPI$ReturnBetas) > 1) {
		  BetaFit <- as.vector(HitBatWingPI$ReturnBetas);
	  } else {
	    BetaFit = SMSReal * 0 - 10;
	  }
    if ((is.logical(DoCI) && DoCI == TRUE) ||
      (is.numeric(DoCI) && DoCI >= 1)) {
      CITime1 <- proc.time();
      if (SMS$LogitRegression == FALSE) { 
       HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior=SigmaPrior, PiAPrior=PiAPrior, DoLogit=FALSE); 
      } else {
       if (SMS$n < SMS$p) { LogitNoise = SMS$n/SMS$p} else {LogitNoise = 1.0; }
       HitBatWingFast <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = SMS$puse, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = max(abs(SMS$puse * NCOL(SMS$XX)),5), RecordFlag = TRUE, HoldOn = FALSE, 
       TDFNu = SMS$tNoiseDF, SigmaPrior=SigmaPrior, PiAPrior=PiAPrior, DoLogit=TRUE,
       LogitNoise=LogitNoise);         
      }
      MIPReport <- NULL;
      if (!is.null(HitBatWingFast$TLS$RecBBOn1)) {
        try(MIPReport <- HitBatWingFast$TLS$RecBBOn1[,2]);
      } else {
        try(MIPReport <- HitBatWingFast$TLS$BBOn1);
      }
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      HitBatWingFast$TLS$ConfidenceQuantiles <- GetConfidenceIntervals;
      HitBatWingFast$TLS$LambdaIndexForConfidenceIntervals <- length(LambdaDSeq)-2
      CITime2 <- proc.time();  CITime <- CITime2-CITime1;
      CIEst <- list();
      MIPMIP <-  HitBatWingFast$TLS$RecBBOn1[, length(LambdaDSeq)-1];
      CIEst[[1]]<- HitBatWingFast$TLS$ConfidenceMatrix;
      CIEst[[2]]<- HitBatWingFast$TLS$UnshrunkConfidenceMatrix;
      CIEst[[3]] <- HitBatWingFast$TLS$HPDMatrix;
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix");
      try(CIEst <- CleanCIEst(CIEst));
      CIQuantiles <- HitBatWingFast$TLS$ConfidenceQuantiles;
    } else {
      CIEst <- NULL;   CIQuantiles  <- NULL;  CITime = NULL;
      MIP = MIPMIP;
    }            
	    T2 <- proc.time();
	    LarsFitPIT <- T2 - T1;
	    		 return(list( type="piAM2Lasso",
	                BetaFit = BetaFit,
                  BBHatsFit = LarsFitBBsPI,
	                FitTime =LarsFitPIT,
	                OtherBetas = HitBatWingPI$ReturnBetas,
                  EMLARSObject = HitBatWingPI, MIPReport=MIPReport,
                  CIEst = CIEst, CIQuantiles=CIQuantiles, CITime = CITime) );
}


###############################################################################
###  GenerateFDLimitLasso()
###
###    A function for using Fermi-Dirac mean-field based EM2Lasso to fit data Y~X
###
###    This type of 2Lasso fits an active set size, this is parameter FixKA sent
###     to the EM2Lasso function.  This function merely counts the number of true
###     non-zero parameters in "SMS" and gives that to EM2Lasso
###
###    Noise parameter is also taken from true value in SMS
###
###
####    Smoothing constants for limit lasso here are 2^(1/8), 2^(-1/8), rather
####     slow, but a pain in the but to play with.  TotalRuns=100 so that
####    Final lambdaA and lambdaD values for "Limit-Lasso" estimate are
GenerateFDLimitLasso <- function(SMS, LARSSeekFlag, DoCI = DefaultDoCI,...) {
  library(TwoLassoCpp)
  USENOISE <- SMS$SigmaSqNoise;	  kLen <- SMS$kLen;  
	T1 <- proc.time();	  
  if (LARSSeekFlag == 0) {
	  StlambdaA = .5/ (2*USENOISE);
	  StlambdaD = .5 / (2*USENOISE);
  } else if (LARSSeekFlag >= 1) {
    if (LARSSEEKMETHOD == 1)   {
      rnd <- round(SMS$puse * kLen);  rnd = min( kLen, rnd);
      rnd <- SMS$kActiveLen;
      MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)	  
	    ITO <- MyLarsGiveMeClose(MyFit, IntClose = rnd, MinMin = 0.000000001);   
	    StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);
	    StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE * LARSSeekFlag);     
    } else {
      rnd1 <- round(kActiveLen * LARSSeekFlag);  rnd1 = min( kLen, rnd1);
      MyFit <- LarsCC2(xxs = SMS$XX, yys = SMS$YY, lambda = 0.000000001)	  
	    ITO <- MyLarsGiveMeClose(MyFit, IntClose = rnd1, MinMin = 0.000000001);   
	    StlambdaA =  MyFit$Lambda[ITO] / (2*USENOISE);
	    StlambdaD = MyFit$Lambda[ITO] / (2*USENOISE);              
    }     
	} else if (LARSSeekFlag < 0) {
	  rnd <- round(SMS$puse * kLen);   rnd = min( kLen, rnd);
	  rnd <- SMS$kActiveLen;
    StlambdaA = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;
    StlambdaD = abs(LARSSeekFlag) * sqrt(USENOISE) / rnd;  
  }           
	FixKa =  SMS$kActiveLen;
	if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
    TotalRuns = 100;  lambdaAmultStop = 50;
    LambdaDK <- StlambdaD/10 * (2^(1/8))^(0:TotalRuns);
    LambdaAK <- StlambdaA/10 * (2^(-1/8))^(0:TotalRuns);
    LambdaAK[0:TotalRuns > lambdaAmultStop] <- LambdaAK[lambdaAmultStop];  
    OrderSeq = c(8,rep(4, TotalRuns-1),8);
    
    HitBatWingFD <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = -SMS$puse, 
      FixKa = FixKa, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = FixKa, RecordFlag = FALSE, HoldOn = FALSE, 
       CauchyEpsilon = .00001, MaxCauchy = 150,
       TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = NULL);
  } else if (SMS$LogitRegression == TRUE) {
    TotalRuns = 100;
    LambdaAK = StlambdaA/10 * (2^(-1/8))^(0:TotalRuns);
    LambdaDK = StlambdaD/10 * (2^(1/8))^(0:TotalRuns);
    OrderSeq = c(8,rep(4,TotalRuns-1), 8);
    if (FALSE) {
    HitBatWingFD = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
      StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
      LambdaAK = LambdaAK,LambdaDK = LambdaDK,OrderSeq = OrderSeq,
      StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
      MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
      TotalRuns = 5, MaxLambdaASteps = 50,
      RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
      MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = FixKa, PriorPi = -1)    
    }
    HitBatWingFD <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = -SMS$puse, 
      FixKa = FixKa, 
       SigmaSq = USENOISE, LambdaAK = LambdaAK,
       LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = FixKa, RecordFlag = FALSE, HoldOn = FALSE, 
       CauchyEpsilon = .00001, MaxCauchy = 150,
       TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = NULL, DoLogit = TRUE);
  }
  if (is.null(HitBatWingFD)) {
    AFilePrint("FDLimitLasso: Fit Hit BatWing FD is NULL! ");  flush.console();
  } else if (is.null(HitBatWingFD$FinalBeta)) {
    AFilePrint("FDLimitLasso: Fit FinalBeta is Null!");flush.console();
  } else if (length(HitBatWingFD$FinalBeta) != kLen) {
    AFilePrint("FDLimitLasso: Fit Final Beta is bad length"); flush.console();
  } else {
    BetaFit <- HitBatWingFD$FinalBeta;
  }
  MIPReport <- NULL;
  try(MP <- sort(abs(USENOISE-.25*(LambdaDK-LambdaAK)^2), index=TRUE)$ix[1]);
  if (!is.null(HitBatWingFD$RecBBOn1) && length(dim(HitBatWingFD$RecBBOn1)) == 2 &&
    NCOL(HitBatWingFD$RecBBOn1) >= 1) {
    try(MIPReport <- HitBatWingFD$RecBBOn1[,MP]);
  } else if (!is.null(HitBatWingFD$BBHatsAll) && length(dim(HitBatWingFD$BBHatsAll)) == 2 &&
    NROW(HitBatWingFD$BBHatsAll) >= 1) {
    try(MIPReport <- HitBatWingFast$BBHatsAll[MP, ]);
  } else if (!is.null(HitBatWingFD$BBHatsFinal)) {
    try(MIPReport <- HitBatWingFD$BBHatsFinal);
  } else {
    try(MIPReport <- HitBatWingFD$BBOn1);
  }
  if (!is.null(HitBatWingFD) && !is.null(HitBatWingFD$FailureFlag) &&
          is.numeric(HitBatWingFD$FailureFlag) &&
          length(HitBatWingFD$FailureFlag) && HitBatWingFD$FailureFlag ==1 ) {
	     AFilePrint(paste("HitBatWingFD Failed, we tried FixKa = ", HitBatWingFD$FixKa, sep="")); 
	     OnFails <<- OnFails+1;
	     FailToFitSim[[OnFails]] <<- SMS;
         FailToFitEMOb[[OnFails]] <<- HitBatWingFD;
  }
  if (length(BetaFit) == kLen) {
      BBHatsFit = BetaFit;
      BBHatsFit[BetaFit!= 0] = 1;      
  } else if (length(HitBatWingFD$ReturnBetas) == kLen) {
           BetaFit = HitBatWingFD$ReturnBetas;
           BBHatsFit = BetaFit;
           BBHatsFit[BetaFit!= 0] = 1;
  }  else {
       BetaFit = rep(0, kLen);
       BBHatsFit = rep(0, kLen);
  }       
                  
  if (length(BetaFit) > 1) {
   LarsFitBBsFD <- as.vector(BetaFit);
   LarsFitBBsFD[abs(BetaFit) <= MinMin] = 0;
   LarsFitBBsFD[abs(BetaFit) > MinMin] = 1;  
  } else {
    LarsFitBBsFD = SMS$BetasReal * 0 - 10;
  }
  T2 <- proc.time();
  LarsFitFDT <- T2 - T1;	   
 ##AFilePrint(paste("ii4 = ", ii4, ", starting with HitBatWingFD"));
  if ((is.logical(DoCI) && DoCI == TRUE) ||
      (is.numeric(DoCI) && DoCI >= 1)) {
      CITime1 <- proc.time();
      OOSeq <- c(8, rep(3, length(LambdaAK)-2),8)
    if (SMS$LogitRegression == FALSE) { 
       HitBatWingFD <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = -SMS$puse, 
        FixKa = FixKa, 
        SigmaSq = USENOISE, LambdaAK = LambdaAK,
        LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
        InitKKs = FixKa, RecordFlag = 1, HoldOn = FALSE, 
        CauchyEpsilon = .00001, MaxCauchy = 150,
        TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = NULL, DoLogit=FALSE);
    } else {
       HitBatWingFD <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = -SMS$puse, 
        FixKa = FixKa, 
        SigmaSq = USENOISE, LambdaAK = LambdaAK,
        LambdaDK = LambdaDK, OrderSeq = OrderSeq, Verbose = -1, 
        InitKKs = FixKa, RecordFlag = 1, HoldOn = FALSE, 
        CauchyEpsilon = .00001, MaxCauchy = 150,
        TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = NULL, DoLogit = TRUE);
      }
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      HitBatWingFD$TLS$ConfidenceQuantiles <- GetConfidenceIntervals;
      HitBatWingFD$TLS$LambdaIndexForConfidenceIntervals <- length(HitBatWingFD$TLS$LambdaDK)-2
      CITime2 <- proc.time();  CITime <- CITime2-CITime1;
      CIEst <- list();
      CIEst[[1]]<- HitBatWingFD$TLS$ConfidenceMatrix;
      CIEst[[2]]<- HitBatWingFD$TLS$UnshrunkConfidenceMatrix;
      CIEst[[3]] <- HitBatWingFD$TLS$HPDMatrix;
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix");
      try(CIEst <- CleanCIEst(CIEst));
      CIQuantiles <- HitBatWingFD$TLS$ConfidenceQuantiles;
    } else {
      CIEst <- NULL;   CIQuantiles  <- NULL;  CITime = NULL;
    }  
 return(list( type="FDLimitLasso", BetaFit = BetaFit,
   BBHatsFit = LarsFitBBsFD, FitTime =LarsFitFDT, OtherBetas = BetaFit,
   EMLARSObject = HitBatWingFD, MIPReport = MIPReport,
   CIEst=CIEst, CIQuantiles=CIQuantiles, CITime=CITime) );        	
}

GenerateFDM2Lasso <- function(SMS, LARSSeekFlag, R=0, DoCI = DefaultDoCI,...) {
  library(TwoLassoCpp);
  USENOISE <- SMS$SigmaSqNoise;	  kLen <- SMS$kLen;  
	T1 <- proc.time();
	SDD <- sqrt(min(apply(SMS$XX,2,sd))^2 * (length(SMS$XX[,1])-1) );
	LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD, 
    sqrt(3.7 / USENOISE) * SDD *  max( (2*length(SMS$XX[,1])), 200));
                             
	SABS <- min( abs(SMS$BetasReal[SMS$BetasReal != 0])  );
	SDA <- SDD / sqrt(length(SMS$XX[,1])-1);
  LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) * sqrt( 2 * SDA^2 / USENOISE)
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
    SigmaVec = c(USENOISE * NN^R, USENOISE, USENOISE);
    LambdaASeq = c(LambdaASeq[1] * NN^(.5 * R), LambdaASeq[1], 
      LambdaASeq[length(LambdaASeq)]);
    LambdaDSeq = c(LambdaDSeq[1] * NN^(-.5 * R), LambdaDSeq[1],
      LambdaDSeq[length(LambdaDSeq)] );
  } else {
    SigmaVec = -1;
  }          
	FixKa =  SMS$kActiveLen;

  if (is.null(SMS$LogitRegression) || SMS$LogitRegression == FALSE) {
    if (FALSE) {
	  HitBatWingFD <- EM2Lasso(xxs=SMS$XX, yys=SMS$YY, ppiuse = -SMS$puse, 
	    FixKa = FixKa, sigmaNoiseSq = USENOISE, RatWant = .2, 
	    LambdaAK=LambdaASeq,LambdaDK = LambdaDSeq,
      OrderSeq=OrderSeq,  TotalRuns = 100, 
	    NumEMConv = 4, MultEMCons = .99, NumCDOConv = 150, CDOEpsilon = .000001,
      TDFNu = SMS$tNoiseDF, SigmaVec=SigmaVec);
    }
    HitBatWingFD <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = -SMS$puse, 
      FixKa = FixKa, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = FixKa, RecordFlag = FALSE, HoldOn = FALSE, 
       CauchyEpsilon = .00001, MaxCauchy = 150,
       TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = NULL);
  } else if (SMS$LogitRegression == TRUE) {
    if (FALSE) {
    HitBatWingFD = GLMG2Lasso(X= SMS$XX, y01=SMS$YY, 
      StartBeta=-999, StartBeta0=0, piA =SMS$puse, sigmaSq=2,
      LambdaAK = LambdaASeq, LambdaDK = LambdaDSeq, OrderSeq = OrderSeq,
      StLambdaA = StlambdaA, StLambdaD = StlambdaD, 
      MultLambdaA = lambdaAMultC,MultLambdaD=lambdaDMultC,
      TotalRuns = 5, MaxLambdaASteps = 50,
      RecordBetaFlag = TRUE, InitKKs=5, PrintFlag = 0,
      MaxCDOEpsilon=.000001, MaxCDOLoops=80, FixKa = FixKa, PriorPi = -1)  
    }
    if (SMS$n < SMS$p) {
      LogitNoise <- SMS$n/SMS$p;
    } else { LogitNoise = 1.0; }
    HitBatWingFD <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = -SMS$puse, 
      FixKa = FixKa, 
       SigmaSq = USENOISE, LambdaAK = LambdaASeq,
       LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
       InitKKs = FixKa, RecordFlag = FALSE, HoldOn = FALSE, 
       CauchyEpsilon = .00001, MaxCauchy = 150,
       TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = NULL, DoLogit=TRUE,
       LogitNoise=LogitNoise);
  }
  if (!is.null(HitBatWingFD) && !is.null(HitBatWingFD$FailureFlag) &&
    is.numeric(HitBatWingFD$FailureFlag) && HitBatWingFD$FailureFlag ==1 ) {
	  AFilePrint(paste("HitBatWingFD Failed, we tried FixKa = ", HitBatWingFD$FixKa, sep="")); 
	  OnFails <<- OnFails+1;
	  FailToFitSim[[OnFails]] <<- SMS;
    FailToFitEMOb[[OnFails]] <<- HitBatWingFD;
  }
    if (length(HitBatWingFD$FinalBeta) == kLen) {
      BetaFit = HitBatWingFD$FinalBeta;
      BBHatsFit = BetaFit;
      BBHatsFit[BetaFit!= 0] = 1;    
    } else if (length(HitBatWingFD$ReturnBetas) == kLen) {
      BetaFit = HitBatWingFD$ReturnBetas;
      BBHatsFit = BetaFit;
      BBHatsFit[BetaFit!= 0] = 1;
    } else {
      BetaFit = rep(0, kLen);
      BBHatsFit = rep(0, kLen);
    }       
                    
	if (length(BetaFit) > 1) {
		LarsFitBBsFD <- as.vector(BetaFit);
		LarsFitBBsFD[abs(BetaFit) <= MinMin] = 0;
		LarsFitBBsFD[abs(BetaFit) > MinMin] = 1;  
	} else {
	  LarsFitBBsFD = SMS$BetasReal * 0 - 10;
	}
  T2 <- proc.time();
  MIPReport <- NULL;
  if (!is.null(HitBatWingFD$TLS$RecBBOn1) && length(dim(HitBatWingFD$TLS$RecBBOn1)) == 2 &&
    NROW(HitBatWingFD$TLS$RecBBOn1) >= 1) {
    try(MIPReport <- HitBatWingFD$TLS$RecBBOn1[,2]);
  } else if (!is.null(HitBatWingFD$RecBBOn1) && length(dim(HitBatWingFD$RecBBOn1)) == 2 &&
    NROW(HitBatWingFD$RecBBOn1) >= 1) {
    try(MIPReport <- HitBatWingFD$RecBBOn1[,2]);
  } else   if (!is.null(HitBatWingFD$BBHatsAll) && length(dim(HitBatWingFD$BBHatsAll)) == 2 &&
    NROW(HitBatWingFD$BBHatsAll) >= 1)  {
    try(MIPReport <- HitBatWingFD$BBHatsAll[2,]);
  } else if (!is.null(HitBatWingFD$BBHatsFinal)) {
    try(MIPReport <- HitBatWingFD$BBHatsFinal);    
  } else {
    try(MIPReport <- HitBatWingFD$BBOn1);
  }
	LarsFitFDT <- T2 - T1;	  
  if ((is.logical(DoCI) && DoCI == TRUE) ||
      (is.numeric(DoCI) && DoCI >= 1)) {
      CITime1 <- proc.time();
	##AFilePrint(paste("ii4 = ", ii4, ", starting with HitBatWingFD"));
    if (SMS$LogitRegression == FALSE) { 
       HitBatWingFD <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = -SMS$puse, 
        FixKa = FixKa, 
        SigmaSq = USENOISE, LambdaAK = LambdaASeq,
        LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
        InitKKs = FixKa, RecordFlag = 1, HoldOn = FALSE, 
        CauchyEpsilon = .00001, MaxCauchy = 150,
        TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = NULL, DoLogit=FALSE);
    } else {
      if (SMS$n < SMS$p) {
        LogitNoise <- SMS$n/SMS$p
      } else {
        LogitNoise <- 1.0;
      }
       HitBatWingFD <- TwoLassoCpp(X=SMS$XX,Y=SMS$YY, PiA = -SMS$puse, 
        FixKa = FixKa, 
        SigmaSq = USENOISE, LambdaAK = LambdaASeq,
        LambdaDK = LambdaDSeq, OrderSeq = OrderSeq, Verbose = -1, 
        InitKKs = FixKa, RecordFlag = 1, HoldOn = FALSE, 
        CauchyEpsilon = .00001, MaxCauchy = 150,
        TDFNu = SMS$tNoiseDF, SigmaPrior = NULL, PiAPrior = NULL, DoLogit = TRUE,
        LogitNoise=LogitNoise);
      }
      eval(parse(text=GetG0Text("GetConfidenceIntervals")));
      HitBatWingFD$TLS$ConfidenceQuantiles <- GetConfidenceIntervals;
      HitBatWingFD$TLS$LambdaIndexForConfidenceIntervals <- length(HitBatWingFD$TLS$LambdaDK)-2
      CITime2 <- proc.time();  CITime <- CITime2-CITime1;
      CIEst <- list();
      CIEst[[1]]<- HitBatWingFD$TLS$ConfidenceMatrix;
      CIEst[[2]]<- HitBatWingFD$TLS$UnshrunkConfidenceMatrix;
      CIEst[[3]] <- HitBatWingFD$TLS$HPDMatrix;
      names(CIEst) <- c("ConfidenceMatrix", "UnshrunkConfidenceMatrix", "HPDMatrix");
      try(CIEst <- CleanCIEst(CIEst));
      CIQuantiles <- HitBatWingFD$TLS$ConfidenceQuantiles;
    } else {
      CIEst <- NULL;   CIQuantiles  <- NULL;  CITime = NULL;
    } 
	return(list( type="FDM2Lasso", BetaFit = BetaFit,
    BBHatsFit = LarsFitBBsFD, FitTime =LarsFitFDT,
	  OtherBetas = BetaFit, EMLARSObject = HitBatWingFD, MIPReport = MIPReport,
      CIEst=CIEst, CIQuantiles=CIQuantiles, CITime=CITime) );        	
}


###############################################################################
###  GenerateXLMMedian()
###
###    A function for using Marginal Median inference
###
###    This is not our algorithm, but one taken from Xiao-Li Meng's
###     "Who cares if it is a white or a black cat" example
###
###    Theoretically, the prior is modelled as a mixture of a continuous prior and
###     a delta function.  Because this posterior is hard to integrate, only the
###     marginals are inspected and a pseudo-inverse input is used.
###
###    This version only takes SMS object as input, so it steals true piA value
###       noise value, et all from the data
GenerateXLMMedian <- function(SMS, Interval = Interval, DoCI = DefaultDoCI, ...) {
    library(TwoLassoCpp);
	T1 <- proc.time();
	if (!exists("DoCI") || is.null(DoCI)) {
    eval(parse(text=GetG0Text("DoCI", S=1)));
  }
  if ((is.logical(DoCI) && DoCI[1] == TRUE) || 
    (is.numeric(DoCI) && DoCI[1] >= 1)) {
    eval(parse(text=GetG0Text("GetConfidenceIntervals", S=1)));
    ConfidenceQuantiles <- GetConfidenceIntervals;  
  }
	XLFit <- XLOperation(SMS$XX, SMS$YY, SMS$puse, SMS$SigmaSqNoise, 
	  TauOther = 60, Interval=Interval, ConfidenceQuantiles=ConfidenceQuantiles);
	XLFitBBs <- XLFit$BBHatsFinal;
	XLFitBBs[abs(XLFit$BBHatsFinal) <= MinMin] = 0;
	XLFitBBs[abs(XLFit$BBHatsFinal) > MinMin] = 1;   
	XLFitBBs <- as.vector(XLFitBBs);
	XLFitBetas <- XLFit$BBHatsFinal * XLFit$ReturnBetas;
	            
	T2 <- proc.time();

	XLFitPIT <- T2 - T1;
	return(list( type="XLMMedian", BetaFit = XLFitBetas,
    BBHatsFit = XLFitBBs, FitTime = XLFit$PointTime,
	  OtherBetas = XLFitBetas, XLFitObject = XLFit, CIEst = XLFit$CIEst,
    CITime = XLFit$CITime, ConfidenceQuantiles = XLFit$ConfidenceQuantiles,
    CIQuantiles = XLFit$ConfidenceQuantiles,
    MIPReport = XLFit$PosteriorProbB1) );        		     
}

###############################################################################
###  GenerateXLLimitLasso()
###
###    A function for using Marginal Median inference on a Limit-Lasso inference
###
###    Here we demonstrate that Limit Lasso is a better input than pseudo estimate
###      used in XLMMedian
###
###     To use, first fit data using GenerateLimitLasso();  Then input that
###      object as FitLimitLasso to this algorithm
### 
###    This algorithm tends to do better than LimitLasso itself.  
###
GenerateXLLimitLasso <- function(SMS, FitLimitLasso, Interval = .95, DoCI = DefaultDoCI) {
  library(TwoLassoCpp);
	if (!exists("DoCI") || is.null(DoCI)) {
    eval(parse(text=GetG0Text("DoCI", S=1)));
  }
  if ((is.logical(DoCI) && DoCI[1] == TRUE) || 
    (is.numeric(DoCI) && DoCI[1] >= 1)) {
    eval(parse(text=GetG0Text("GetConfidenceIntervals", S=1)));
    ConfidenceQuantiles <- GetConfidenceIntervals;  
  }
  OtherBetas = NULL;
  if (!is.null(FitLimitLasso$BetaFit)) {
    OtherBetas = FitLimitLasso$BetaFit;
  } else if (!is.null(FitLimitLasso$FinalBeta)) {
    OtherBetas = FitLimitLasso$FinalBeta;
  } else if (!is.null(FitLimitLasso$ReturnBetas)) {
    OtherBetas = FitLimitLasso$ReturnBetas;
  } else if (!is.null(FitLimitLasso$TwoLassoObject) &&
    !is.null(FitLimitLasso$TwoLassoObject$FinalBeta)) {
    OtherBetas = FitLimitLasso$TwoLassoObject$FinalBeta;
  } else if (!is.null(FitLimitLasso$TLS) &&
    !is.null(FitLimitLasso$TLS$Beta)) {
    OtherBetas = FitLimitLasso$TLS$Beta  
  } else if (!is.null(FitLimitLasso$TwoLassoObject$TLS) &&
    !is.null(FitLimitLasso$TwoLassoObject$TLS$Beta)) {
    OtherBetas =  FitLimitLasso$TwoLassoObject$TLS$Beta; 
  }    
	     T1 <- proc.time();
	     XLPostFit <- XLPostAnalysis(
            BetaProposed = FitLimitLasso$BetaFit, 
	          xxs=SMS$XX, yys=SMS$YY, ppiuse = SMS$puse, sigmaNoiseSq = SMS$SigmaSqNoise, 
	         tausqA = max(60, .5 / max(SMS$SigmaSqNoise,1) * sum(SMS$XX[,1])), 
           ConfidenceQuantiles=ConfidenceQuantiles);
	     XLPostFitBBs <- XLPostFit$BBHatsFinal;
	     XLPostFitBBs[abs(XLPostFit$BBHatsFinal) <= MinMin] = 0;
	     XLPostFitBBs[abs(XLPostFit$BBHatsFinal) > MinMin] = 1;  
	     XLPostFitBBs <- as.vector(XLPostFitBBs);	     
	     XLPostFitBetas <- XLPostFit$ReturnBetas; 
       

       	               
	     
	     try(FitTime <- XLPostFit$PointTime[1:3])
	     try(FitTime <- XLPostFit$PointTime[1:3] + FitLimitLasso$FitTime[1:3]);
	     T2 <- proc.time();
	     XLPostFitPIT <- T2 - T1;	   
	    		 return(list( type="XLMMedian",
	                BetaFit = XLPostFitBetas,
                  BBHatsFit = XLPostFitBBs,
	                FitTime = FitTime,
	                OtherBetas = FitLimitLasso$FinalBeta,
                  XLFitObject = XLPostFit, CIEst = XLPostFit$CIEst,
                  CITime = XLPostFit$CITime, ConfidenceQuantiles=ConfidenceQuantiles,
                  CIQuantiles = ConfidenceQuantiles,
                  MIPReport = XLPostFit$PosteriorProbB1) );                
}


  colPermaList <- c("red", "pink",  "purple", "orange", "blue", "green", "cyan", 
         "magenta", "lightblue",
         "aquamarine", "brown")
  colPermaList2 <- c("red", "pink",  "salmon", "green", "cyan", "lightblue",
               "magenta", "aquamarine", "turquoise", "purple", "blue", "orange", "darkblue", "brown",
               colors()[107], colors()[551], colors()[120], colors()[567], 
               colors()[624], colors()[652], colors()[655],
               colors()[477],colors()[461],colors()[598],colors()[639], colors()[179], colors()[25],
               colors()[125],
               colors()[23]);
  AlphPermaList <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V","W","X","Y","Z","1","2", "3");      
  pchListDefault = c(15, 15, 15, 16, 16, 16, 18, 16, 16, 15, 16, 16,16, 15,16,16,16, 15, 15, 18, 18,16,16,16,16,24,24,16,24,24); 
  UnderColor <- c("black", "black",  "black", "black", "black", "black",
               "black", "black", "black", "black", "black", "black", "white", "black", "white", "white", "white", "white", "white",
               "black", "black","white","white","white","white", "white", "black", "black",
               "white", "white");
  UnderColor[c(26,27)]= "black" 
  PointStartA = 52;
    PointDiv1 = 11;
    PointDiv2 = 22; DB2 = 8;
  PointStartB = 52 
    PointDiv3 = 11;
    PointDiv4 = 22;   DB4 = 8;                    
######################################################################
##   PrintMyPlotterV2  
##
##   Useful in printing out data in Points distance format.
##
##
##
## 
PrintMyPlotterV3<- function(jobsL, OneVV2, PrMeVec, TMM = FALSE,
  colList, AlphList, xlims = -999, ylims = -999,pchList = pchListDefault,
  UnderCol = UnderColor) {
##PnamesV2 <- c("Sim Real", "Lasso-Fixed 7", "Lasso-Cp", "Lasso Min-Yuan", 
##              "Limit-Ridge", "2-Lasso-Fast", "Limit-2-Lasso", "Limit-2PI-Lasso", "Posterior Median",
##              "Limit-Lasso FD);
  PnV2 <- c("Lasso-Fixed", "Lasso-Cp", "Lasso w=1 LY", 
              "Limit-Ridge", "2-Lasso-9X", "Limit-2-Lasso", "Limit-2PI-Lasso", "Marg Median",
              "Limit-Lasso FD", "MargM Lim-Lasso");
  PnV3 <- c("Lasso-\nFixed", "Lasso-\nCp", "Lasso \nw=1 LY", 
              "Limit-\nRidge", "2-Lasso-\n9X", "Limit-2\n-Lasso", "Limit-2PI\n-Lasso", "Marg\n Median",
              "Limit-\nLasso FD", "MargM Lim\n-Lasso");   
  PnV4 <- c("Lasso-\nFixed", "Lars-\nCp", "Lasso-\nStodden",
    "Limit-\nLasso", "Limit-2PI\n-Lasso", "Lim-\nLasso FD",
    "Marg\n Median", "MargM Lim\n-Lasso", "Ignorant \n Limit-Lasso", 
     "Lasso \nw=1 LY", "2-Lasso\n9X", "Limit-Ridge", "M2Lasso", "SCAD-minConvex", 
     "M2PiLasso", "M2FDLasso", "M2MMedian", "SCAD-4xStodden", "SCAD-minWLT", 
     "ENet-Fixed", "ENet-CVmin",
     "R2Lasso", "R2piLasso", "R2FDLasso", "R2XLMLasso", "HorseShoe", "B2Lasso",
     "Bayes Lasso"); 
                             
   cto = 0;
   ALT2 <- 10 * 11 / 2;
   if ( is.vector(OneVV2)) {
        iin <- 1; jjn <- 1;
        par(mfrow=c(1,1));
          OneVV2 <- t( matrix( c(OneVV2, OneVV2), length(OneVV2), 2));
   } else if ( length(OneVV2[,1]) == 12) {
      iin <- 3;
      jjn <- 4;
       matlayout <- rbind(
     c( 24,15,25 ,16,26 ,17,27),
     c(1,2,3,4,5,6,7),
     c(28 ,18,29 ,19,30 ,20,31),
     c(8,9,10,11,12,13,14),
     c(32 ,21,33 ,22,34 ,23,35)   );
nf <- layout(matlayout, c(1,.02, 1,.02, 1,.02, 1), c(1,.02, 1, .02,1), TRUE)
       matlayout <- rbind(
     c(1, 25,16,26 ,17,27 ,18,28),
     c(1, 2,3,4,5,6,7,8),
     c(1, 29 ,19,30 ,20,31 ,21,32),
     c(1, 9,10,11,12,13,14,15),
     c(1, 33 ,22,34 ,23,35 ,24,36)   );
nf <- layout(matlayout, c(.5, 1,.02, 1,.02, 1,.02, 1), c(1,.02, 1, .02,1), TRUE)
  par(plt=c(0,1,0,1));
  plot(1,1, xlab="", ylab="", axes=FALSE, type="n", ylim=c(0,1), xlim=c(0,1));
      nNames <- PnV4[PrMeVec];
      ylen <- length(AlphList) + 4;
      yspots <- (ylen - .5  - 1:length(AlphList))/ylen;
    for (ii in 1:length(AlphList)) {
      points(x=.09, y=yspots[ii], pch=pchList[ii], 
          cex=(PointStartA - ii) / PointDiv1, col=colList[ii]);
      points(x=.092, y=yspots[ii]+.000, pch=AlphList[ii], col=UnderCol[ii], cex=(PointStartA - ii - DB2) / PointDiv2);
      text( x=.55, y=yspots[ii]-.001, labels=paste(AlphList[ii], "-", nNames[ii]) );
    }
    text(x=.5, y=yspots[length(AlphList)]-1/ylen, labels=("X-axis = Type II\nTrue Missed"));
    text(x=.5, y=yspots[length(AlphList)]-2.5/ylen, labels=("Y-axis = Type I\nFalse Pos"));    
##layout.show(nf)
    par(plt=c(0,1,0,1));
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8), type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
     par(plt=c(0,1,0,1));
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8), type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       par(plt=c(0,1,0,1)); 
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");      
      par(plt=c(.13,.975,.155,.83))
   } else if ( length(OneVV2[,1]) == 8) {
      iin <- 2;
      jjn <- 4;
      par(mfrow=c(2,4));
   } else if ( length(OneVV2[,1]) == 4) {
      iin <- 2;
      jjn <- 2;
      par(mfrow=c(2,2));
   } else if (length(OneVV2[,1]) == 6) {
      iin <- 2;
      jjn <- 3;
      matlayout <- rbind(
     c(1, 11,7,12 ,8,13),
     c(1, 2,3,4,5,6),
     c(1, 14 ,9,15 ,10,16)   );
     nf <- layout(matlayout, c(.5, 1,.02, 1,.02, 1), c(1,.02, 1), TRUE) 
     par(plt=c(0,1,0,1));
     plot(1,1, xlab="", ylab="", axes=FALSE, type="n", ylim=c(0,1), xlim=c(0,1));
      nNames <- PnV4[PrMeVec];
      ylen <- length(AlphList) + 4;
      yspots <- (ylen - .5  - 1:length(AlphList))/ylen;
    for (ii in 1:length(AlphList)) {
      points(x=.09, y=yspots[ii], pch=pchList[ii], cex =(PointStartA-ii)/ PointDiv1, col=colList[ii]);
      points(x=.092, y=yspots[ii]+.000, pch=AlphList[ii] , cex=(PointStartA - ii - DB2) / PointDiv2,
         col=UnderCol[ii]);
      text( x=.55, y=yspots[ii]-.001, labels=paste(AlphList[ii], "-", nNames[ii]) );
    }
    text(x=.5, y=yspots[length(AlphList)]-1/ylen, labels=("X-axis = Type II\nTrue Missed"));
    text(x=.5, y=yspots[length(AlphList)]-2.5/ylen, labels=("Y-axis = Type I\nFalse Pos"));  
     par(plt=c(0,1,0,1));
    for (ii in 2:6)    {
      plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main=""); 
    } 
    for (ii in 7:10) {
     plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main=""); 
    }  
     par(plt=c(.13,.975,.155,.83))               
   } else if ( length(OneVV2[,1]) == 9) {
      iin <- 3;
      jjn <- 3;
      par(mfrow=c(3,3));
             matlayout <- rbind(
     c( 24,15,25 ,16,26 ,17,27),
     c(1,2,3,4,5,6,7),
     c(28 ,18,29 ,19,30 ,20,31),
     c(8,9,10,11,12,13,14),
     c(32 ,21,33 ,22,34 ,23,35)   );
     
nf <- layout(matlayout, c(1,.02, 1,.02, 1,.02, 1), c(1,.02, 1, .02,1), TRUE)
       matlayout <- rbind(
     c(1, 18,12,19 ,13,20),
     c(1, 2,3,4,5,6),
     c(1, 21 ,14,22 ,15,23),
     c(1, 7,8,9,10,11),
     c(1, 24 ,16,25 ,17,26)   );
nf <- layout(matlayout, c(.5, 1,.02, 1,.02, 1), c(1,.02, 1, .02,1), TRUE)
  par(plt=c(0,1,0,1));
  plot(1,1, xlab="", ylab="", axes=FALSE, type="n", ylim=c(0,1), xlim=c(0,1));
      nNames <- PnV4[PrMeVec];
      ylen <- length(AlphList) + 4;
      yspots <- (ylen - .5  - 1:length(AlphList))/ylen;
    for (ii in 1:length(AlphList)) {
      points(x=.09, y=yspots[ii], pch=pchList[ii], cex =(PointStartA-ii)/ PointDiv1, col=colList[ii]);
      points(x=.092, y=yspots[ii]+.000, pch=AlphList[ii] , cex=(PointStartA - ii - DB2) / PointDiv2,
         col=UnderCol[ii]);
      text( x=.55, y=yspots[ii]-.001, labels=paste(AlphList[ii], "-", nNames[ii]) );
    }
    text(x=.5, y=yspots[length(AlphList)]-1/ylen, labels=("X-axis = Type II\nTrue Missed"));
    text(x=.5, y=yspots[length(AlphList)]-2.5/ylen, labels=("Y-axis = Type I\nFalse Pos"));    
##layout.show(nf)
    par(plt=c(0,1,0,1));
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8), type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    ##plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    ##plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
     par(plt=c(0,1,0,1));
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8), type="l", lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       lines(x=c(0,0), y=c(0,1), lwd=3, col="black");
    plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    ##plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
    ##plot(x=c(0,1), y=c(0,0), xlim=c(.2,.8), ylim=c(.2,.8),type="l",  lwd=4, col="black", axes=FALSE, xlab="", ylab="", main="");
       par(plt=c(0,1,0,1)); 
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      ##plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
      ##    type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      ##plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
      ##    type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
          type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");  
      ##plot(x=c(0,0), y=c(0,1), xlim=c(.2,.8), ylim=c(.2,.8),
      ##    type="l", lwd=8, col="black", axes=FALSE, xlab="", ylab="", main="");      
      par(plt=c(.13,.975,.155,.83))
   } else if ( length(OneVV2[,1]) == 16) {
      iin <- 4; jjn <- 4;
      par(mfrow=c(4,4));
   } else {
      return;
   }   
   AFilePrint("Ready To Put Down limits");
   flush.console();
       for (ii in 1:iin) {
          for (jj in 1:jjn) {
                      cto <- cto +1;
            if (length(xlims) == 1 && xlims == -999) {
               xlimss = c(0,7);
            } else if (length(xlims) == 1) { 
               xlimss = c(0, xlims);
            } else if (length(xlims) == 2) {
               xlimss = c(xlims[1], xlims[2]);
            } else if (length(xlims) < length(OneVV2[,1])) {
               xlimss = c(0,7);
            } else if (length(xlims) == length(OneVV2[,1])) {
               xlimss = c(0, xlims[  cto]);
            } else if (length(xlims) < 2 * length(OneVV2[,1])) {
               xlimss = c(0,7);
            } else if (length(xlims) == 2 * length(OneVV2[,1])) {
               xlimss = xlims[  cto ,];
            } else { 
               xlimss = c(0,7);
            }
            JBM = jobsL[[cto]]
            if (length(ylims) ==1 && ylims == -999) {
               ylimss =  c(0, min(c(6, max(OneVV2[cto,2] +.5,1))))
               LL = length(JBM[1,]) / NumDetails;
               
               MN = 6;
               for (ii in 1:length(PrMeVec)) {
                   TV = JBM[, 3 + NumDetails * (PrMeVec[ii]-1)];
                   TV = as.numeric(TV);
                   TV = TV[!is.na(TV) & TV >= 0 & !is.nan(TV)];  TV2 = TV;
                   if (PrMeVec[ii] == 2 || PrMeVec[ii]==3) {
                      ## AFilePrint(paste("mean(TV) for PrMeVec[",ii,"] is ",
                      ##        mean(TV), sep=""));
                      ## TV2 = JBM[, 3 + NumDetails * (PrMeVec[ii]-1)];
                      ## TV2 = TV2[!is.na(TV2) & TV2 >= 0];    
                       ##AFilePrint(paste("mean(TV2) for PrMeVec[",ii,"] is ",
                       ##       mean(TV2), sep=""));   
                       ##AFilePrint(paste(" That Is cto = ", cto, sep=""));
                       ##AFilePrint(paste(" n = ", OneVV2[cto,1], "; p = ", OneVV2[cto,2],
                       ##            "; k = ", OneVV2[cto,3], "; sigma= ",
                       ##            OneVV2[cto,4], "; CovRho = ", CovRho, sep=""));                                              
                   }
                   if (length(TV2) <= 0) { MTV = MN}  else {
                   MTV <- mean(TV2[is.na(TV2) == FALSE]); }
                   if (length(TV2) >= 1) {
                      MN = max( MTV, MN);
                   }   else {
                       AFilePrint(paste("This Is Weird, TV[", ii, "] or  ",
                          "NameFunctions[", ii, "] = ",
                          NameFunctions[ii], "has no data \n\n\n*******\n\n", sep="") );
                   }
               }
               if (length(MN) == 0 || is.null(MN) || is.na(MN)) {
                  MN = 10;
               }
               if (length(MN) > 1) { MN = MN[1];}
               if (MN > 10) {
                   ylimss = c(-1,-MN);
                   AFilePrint(paste("We Got a Big MN = ", MN, sep=""));
               }  else if (MN > 6) {
                    AFilePrint(paste("We moved MN up at least = ", MN, sep=""));
                    ylimss = c(0,MN);
               }  else {
                  ylimss = c(0,MN);
              }
            } else if (length(ylims) == 1) { 
               ylimss = c(0, ylims);
            } else if (length(ylims) == 2) {
               ylimss = c(ylims[1], ylims[2]);
            } else if (length(ylims) < length(OneVV2[,1])) {
               ylimss =  c(0, min(c(6, max(OneVV[cto,2] +.5,1))))
            } else if (length(ylims) == length(OneVV2[,1])) {
               ylimss = c(0, ylims[  cto + jj]);
            } else if (length(ylims) < 2 * length(OneVV2[,1])) {
               ylimss =  c(0, min(c(6, max(OneVV[cto,2] +.5,1))))
            } else if (length(ylims) == 2 * length(OneVV2[,1])) {
               ylimss = ylims[  cto ,];
            } else { 
               ylimss =  c(0, min(c(10, max(OneVV2[cto,2] +.5,1))))
            }            
            
            ##AFilePrint(paste("Now Printing SubSPlot[[", cto, "]]", sep=""));
            ##flush.console();

            SubSPlot <- jobsL[[cto]]
            ##if (max(colMeans(SubSPlot[, 3 + (PrMeVec-1) *NumDetails])) > max(ylimss)) {
            ##    ylimss[2] = max(colMeans(SubSPlot[, 3 + (PrMeVec-1) *NumDetails]));
            ##}
            ##plot(x=c(0,1), y=c(0,1), type="n", 
            ##       main=paste("n = ", OneVV[cto,4], ", kppa = ", OneVV[cto,5], 
            ##                  ", sig^2 = ", OneVV[cto,6], sep=""), ylab="", xlab="",
            ##       xlim=xlimss, ylim = ylimss
            ##     );
            if (ylimss[1] >= 0) {
              plot(x=c(0,1), y=c(0,1), type="n",, ylab="", xlab="",
                   xlim=c(-.1,xlimss[2]), ylim = ylimss
                 );            
            } else { 
              ATM = max(abs(ylimss[2])); 
              YL = c(.05, .1, .5, 1, 5, 10, 20, 40, 70, 100,200,500,1000,5000, 10000);
              YL = c(0, .5,1,5,10,20,35,50,75,100, 200, 500, 1000, 5000, 10000)
              ATYL = min(YL[YL > ATM]);
              ##plot(x=c(0,1), y=c(0,1), type="n",, ylab="", xlab="",
              ##         xlim=c(-.1, xlimss[2]), ylim=c(-2,log(ATYL,2) ), axes=FALSE );
              plot(x=c(0,1), y=c(0,1), type="n",, ylab="", xlab="",
                       xlim=c(-.1, xlimss[2]), ylim=c(0,sqrt(ATYL) ), axes=FALSE );              
              axis(1)
              ##axis(2, at=log(YL,2), labels=YL)
              axis(2, at=sqrt(YL), labels=YL);
           }
            title(main=substitute(list(n,p,sigma)== group("(",list(a,x,y),")"),              
                           list(a=OneVV2[cto,1], x=OneVV2[cto,2], y=OneVV2[cto,4])))   
           ## title(main=substitute(list(n,kappa,sigma)== group("(",list(a, x,y),")"),              
           ##                list(a=OneVV2[cto,1], x=OneVV2[cto,2], y=OneVV2[cto,4])))                                
            for (tt in 1:length(PrMeVec)) {
                if (TMM == FALSE) {
                   PrV1 <-as.numeric(SubSPlot[, 3 + (PrMeVec[tt]-1) * NumDetails]);
                   PrV1 <- PrV1[!is.na(PrV1) & PrV1 >= 0 & !is.nan(PrV1) ];
                   PrV2 <-as.numeric(SubSPlot[, 4 + (PrMeVec[tt]-1) * NumDetails]);
                   PrV2 <- PrV2[!is.na(PrV2) & PrV2 >= 0 & !is.nan(PrV2) ];                   
                } else {
                   PrV1 <-as.numeric(SubSPlot[, 3 + (PrMeVec[tt]-1) * NumDetails]);
                   PrV1[is.na(PrV1) | PrV1 < 0 | is.nan(PrV1)] <- max(PrV1[!is.na(PrV1) & PrV1 >= 0 ]);
                   PrV2 <-as.numeric(SubSPlot[, 4 + (PrMeVec[tt]-1) * NumDetails]);
                   PrV2[is.na(PrV2) | PrV2 < 0 | is.nan(PrV2)]  <- max(PrV1[!is.na(PrV2) & PrV2 >= 0 ]);                                
                }
                 PrV1 = as.numeric(PrV1);  PrV2 = as.numeric(PrV2);  
                if (!exists("PrV1")) {
                   PrV1 = 0;
                   AFilePrint("PrV1 Doesn't Exist?");
                }
                if (!exists("PrV2")) {
                   PrV2 = 0;
                   AFilePrint("PrV2 Doesn't Exist?");
                }                
                if (any(is.na(PrV1))) {
                   PrV1[is.na(PrV1)] = 0;
                }
                if (any(is.na(PrV2))) {
                   AFilePrint(paste("PrV2 Has NAs"));
                   AFilePrint(paste("PrV2 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV2), " Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));
                   PrV1[is.na(PrV2)] = 0;
                } 
                if (any(!is.finite(PrV2))) {
                   AFilePrint(paste("PrV2 has Infs"));
                   AFilePrint(paste("PrV2 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV2), " Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));                
                   PrV2[!is.finite(PrV2)] = 0;
                }
                if (any(!is.finite(PrV1))) {
                   AFilePrint(paste("PrV1 has Infs"));
                   AFilePrint(paste("PrV1 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV2), " Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));                
                   PrV1[!is.finite(PrV1)] = 0;
                }  
                if (any(is.nan(PrV2))) {
                   AFilePrint(paste("PrV2 has NANS"));
                   AFilePrint(paste("PrV2 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV2), " Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));                
                   PrV2[is.nan(PrV2)] = 0;
                }
                if (any(is.nan(PrV1))) {
                   AFilePrint(paste("PrV1 Has NANS"));
                   AFilePrint(paste("PrV1 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV2), " Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));                
                   PrV1[is.nan(PrV1)] = 0;
                }                               
                if (is.null(PrV2)) {
                   PrV2 = 0;
                }
                if (is.null(PrV1)) {
                   PrV1 = 0;
                } 
                           
                if (length(PrV2) <= 0) { 
                   AFilePrint(paste("PrV2 is Highly Flawed"));
                   AFilePrint(paste("PrV2 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV2), " Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep="")); 
                   PrV2 = 0;
                }
                if (length(PrV1) <= 0) { 
                   AFilePrint(paste("PrV1 is Highly Flawed"));
                   AFilePrint(paste("PrV1 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV2), " Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep="")); 
                   PrV1 = 0;
                }              
                if (mean(PrV2) < 0 ) {
                    AFilePrint(paste("PrV2 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV2), " Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));
                }
                if (mean(PrV1) < 0 ) {
                    AFilePrint(paste("PrV1 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is ", mean(PrV1), ", Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));
                }
                if (mean(PrV1) > abs(ylimss[2]) ) {
                    ##AFilePrint(paste("PrV1 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                    ##       " is ", mean(PrV1), " > ", abs(xlimss[2]),
                    ##       ", Name = ", NameFunctions[PrMeVec[tt]], 
                    ##       "F Me", sep=""));
                }
                if (length(PrV2) < 0 || mean(PrV2) > abs(xlimss[2]) ) {
                    ##AFilePrint(paste("PrV2 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                    ##       " is ", mean(PrV2), " > ", abs(ylimss[2]),
                    ##       ", Name = ", NameFunctions[PrMeVec[tt]], 
                    ##       "F Me", sep=""));
                }                
                if (length(PrV2) < 0 || is.na(mean(PrV2)) ) {
                    AFilePrint(paste("PrV2 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is NA, Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));
                }
                if (length(PrV1) < 0 || is.na(mean(PrV1)) ) {
                    AFilePrint(paste("PrV1 for PrMeVec[",tt, "] = ", PrMeVec[tt],
                           " is NA, Name = ", NameFunctions[PrMeVec[tt]], 
                           "F Me", sep=""));
                }                
                if (ylimss[1] >= 0)     {
                  points( x=mean(PrV2), y= mean(PrV1), col= colList[tt],
                           pch=pchList[tt], lwd=(PointStartB-tt)/PointDiv3, cex =(PointStartB-tt)/PointDiv3 ); 
                  points( x=mean(PrV2)+.01, y= mean(PrV1)-.01, col=UnderCol[tt],
                           pch=AlphList[tt], cex=(PointStartB  - tt - DB4) / PointDiv4);      
                } else {
                 ## points( x=mean(PrV2), y= log(mean(PrV1),2), col= colList[tt],
                 ##          pch=pchList[tt], lwd=(PointStartB-tt)/PointDiv3, cex =(PointStartB-tt)/PointDiv3 ); 
                 ## points( x=mean(PrV2)+.01, y= log(mean(PrV1),2)-.01, col= "black",
                 ##          pch=AlphList[tt], cex=(PointStartB  - tt - DB4) / PointDiv4);      
                 points( x=mean(PrV2), y= sqrt(mean(PrV1)), col= colList[tt],
                           pch=pchList[tt], lwd=(PointStartB-tt)/PointDiv3, cex =(PointStartB-tt)/PointDiv3 ); 
                 points( x=mean(PrV2)+.01, y= sqrt(mean(PrV1))-.01, col= UnderCol[tt],
                           pch=AlphList[tt], cex=(PointStartB  - tt - DB4) / PointDiv4);                                       
                }                                 
            }
        }
     }
  return(1);
 }
 
 
 AssessCI <- function(CIs, RealBetas) {
     TotalCover <- length(RealBetas[  CIs[,1] <= RealBetas & CIs[,2] >= RealBetas]);
     MeanLength <-  mean( CIs[,2] - CIs[,1]);
     sdLength <- sd(as.vector(CIs[,2]) - as.vector(CIs[,1]));
     MedianLength <- median(CIs[,2] - CIs[,1]);
     
     TotalActive <- length(RealBetas[RealBetas != 0]);
     TotalCoverActive<- length(RealBetas[ RealBetas != 0 & CIs[,1] <= RealBetas &
                                        CIs[,2] >= RealBetas]);
     MeanLengthActive <- mean(CIs[RealBetas!=0,2] - CIs[RealBetas!=0,1]);
     sdLengthActive <- sd(as.vector(CIs[RealBetas!=0,2]) - 
       as.vector(CIs[RealBetas!=0,1]));
     MedianLengthActive <- median(CIs[RealBetas!=0,2] - CIs[RealBetas!=0,1]);

     TotalDeactive <- length(RealBetas) - TotalActive;
     TotalCoverDeActive<- length(RealBetas[ RealBetas == 0 & CIs[,1] <= RealBetas &
                                        CIs[,2] >= RealBetas]);
     MeanLengthDeActive <- mean(CIs[RealBetas==0,2] - CIs[RealBetas==0,1]);
     sdLengthDeActive <- sd(as.vector(CIs[RealBetas==0,2]) - 
       as.vector(CIs[RealBetas==0,1]));
     MedianLengthDeActive <- median(CIs[RealBetas==0,2] - CIs[RealBetas==0,1]);
     
     return(c( TotalCover / length(RealBetas), 
                 MeanLength, sdLength, MedianLength,
               TotalCoverActive / TotalActive,
                 MeanLengthActive, sdLengthActive, MedianLengthActive,
               TotalCoverDeActive / TotalDeactive,
                 MeanLengthDeActive, sdLengthDeActive,
                 MedianLengthDeActive));
 }


z <- FALSE;