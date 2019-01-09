################################################################################
## 0012 CollectSaveResults.r
##   Collect Simulation results saved to many directories and files into
##    easier collated file
##
##     (c) 2009-2019 Alan Lenarcic
##     The code was written to support simulations for Lenarcic and Valdar methods
##      paper.
##
##   CollectSaveResults is largely concerned with "post estimation" stage of the algorithm
##    While simulations have already been run and performance assessed, due to LSF threading
##    they have been saved to many different files.  CollectSaveResults collates
##    this information.
##
##   The TwoSimR5 approach to conducint simulations is
##     1. Decide upon which N_M estimators to use in a problem, install functions
##      through EfficientSimulator and declare parameter settings with DeclareAllTstFunctions
##     2. Define a simulation (n,p,k, sigma, covarianceX, size of Beta...)
##     2. Simulate N_S 500-1000 indepenent simulations from t
##     3. Open N_M * (N_S/N_divisor) threads on Killdevil to attempt to solve each simulation with each estimator
##        (Killdevil requires 1 minute or more processes with 1 hour kill time, so some threads will solve many problems,
##           so each thread gets to solve multiple randomly picked estimator sand simulation problems)
##        Each solution must:
##          a. fit the estimator
##          b. Assess the L1, L2, etc error of the estimator
##          c. Assess Credibility or Model Inclusion performance if applciable
##          d. Save results on this estimator somewhere to remind program to not rerun.
##     4. Rerun threads when an estimator may have failed due to a timeout or memory error thus every
##       estimator gets at least two chances to complete a run.  Most time they don't need any, but sometimes
##       a problem is begun right before a 1 hour timeout.  As problems get larger certain R packages hit
##       memory limits, and this gives us a chance to double check.
##     5.Collect the results of all of these simulations by estimator
##     6.Produce Latex usable tables of results or other interesting RData.   
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
CheckArgFunction <- function() {
  eval(parse(text=GetG0Text("Myargs")));
  eval(parse(text=GetG0Text("NameFunctions")));
  if (length(Myargs) == 1) {
    LTO = 0;
    print(paste("CheckArgFunction: No Arg given for function")); flush.console();
  } else if (length(Myargs) >= 2) {
    print(paste("Working with Myargs[2] = ", Myargs[2], sep=""));    flush.console();
    try(LTO <- round(as.numeric(Myargs[2])));
    if (LTO >= 1 && LTO <= length(NameFunctions)) {
      OnNameFunction <- NameFunctions[LTO];
    }
  } else {
    print("CheckArgFunction: You didn't supply args.  We will default.");
    flush.console();
    LTO <- 0;
  }
  if (LTO < 0 || LTO > length(NameFunctions)) {
    LTO = -999;  
    print("LTO is rejectable we quit"); flush.console();
    quit();
  } else if (LTO == 0) {
    LTO = 0;
    OnFunctionNums = 1:length(NameFunctions);
    RunFunctions <- NameFunctions;
  } else {
    OnFunctionNums = LTO; 
    RunFunctions <- NameFunctions[LTO];
  }
  eval(parse(text=SetGText("RunFunctions", "globalenv()", S=1)));
  eval(parse(text=SetGText("OnFunctionNums", "globalenv()", S=1)));  
  eval(parse(text=SetGText("LTO", "globalenv()", S=1))); 
  return();
}

CollapseOnFunction <- function(OnFunction) {
  MyS <- 0;
  Oldwd <- getwd();
  eval(parse(text=SetGText("Oldwd", "globalenv()", S=1)));
  MyT <- "
    setwd(ASmallContainDir);
    MyS = 1;
  ";
  try(eval(parse(text=MyT)));
  if (MyS == 0) {
    print(paste("CollapseOnFunction: Woah, no ASmallContainDir = ", ASmallContainDir, sep=""));
    flush.console(); return(-999);
  }
  if (!any(unlist(list.files()) == OnFunction)) {
    print(paste("CollapseOnFunction: Woah, no OnFunction = ", OnFunction, " in ", ASmallContainDir, sep=""));
    flush.console(); return(-999);
  }
  MyFiles <- list.files(OnFunction);
  if (length(MyFiles) <= 0) {
    print(paste("CollapseOnFunction: Woah, No files fo OnFunction = ", OnFunction, " in ", ASmallContainDir, sep=""));
    flush.console(); return(-999);
  }
  SMyFiles <- MyFiles[substr(MyFiles, 1, nchar(paste("S", OnFunction, sep=""))) ==  paste("S", OnFunction, sep="")];
  if (length(SMyFiles) <= 0) {
    print(paste("CollapseOnFunction: Woah, No completed SFiles fo OnFunction = ", OnFunction, " in ", ASmallContainDir, 
      " but length(MyFiles) = ", length(MyFiles), sep=""));
    flush.console(); return(-999);    
  }

  BetaStart <- NULL; NoNoiseBetaStart <- NULL; TemperatureList <- NULL;
  load(paste(OnFunction, "//", SMyFiles[1], sep=""));
  if (!exists("MyV") || is.null(MyV) || length(MyV) <= 0) {
    print(paste("CollapseOnFunction: Woah, MyV is not good for this file ", SMyFiles[1], sep=""));
    flush.console();
    return(-999);
  }
  MyOutMatrix <- matrix(NA, length(SMyFiles), length(MyV));
  ProcessVector <- rep(NA, length(SMyFiles));
  SimulationVector <- rep(NA, length(SMyFiles));
  AlternateTimeVector <- rep(NA, length(SMyFiles));
  
  TemperatureListMatrix <- NULL;
  if (!is.null(TemperatureList)) {
    TemperatureListMatrix <- matrix(NA, length(SMyFiles), length(TemperatureList));
  }
  AKS <- c("FinishesTime", "FunctionName", "MyV", "NameFunction", 
    "NowToSaveTime", "RunProcessIdentifier", "SimulationIdentifier",
    "StartsTime");
  RMAKS <- function() {
     Mytxt <- "";
     for (kk in 1:length(AKS)) { 
        Mytxt <- paste(Mytxt,   "
          ", AKS[kk], " <- NULL; rm(",  AKS[kk], ");", sep="")            
     }
     return(Mytxt);
  }          

  eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
  for (iOni in 1:length(SMyFiles)) {
    if (verbose >= 2 || (verbose >=1 && iOni %% 10 == 0)) {
      print(paste("Loading SFile ", iOni, "/", length(SMyFiles), " ",
        SMyFiles[iOni], sep="")); flush.console();
    }
    try(eval(parse(text=RMAKS())));
    try(load(paste(OnFunction, "//", SMyFiles[iOni], sep="")));
    if (exists("MyV") && !is.null(MyV) && length(MyV) >= 1) {
      try(MyOutMatrix[iOni, ] <- MyV[1:NCOL(MyOutMatrix)]);  
    }
    if (exists("RunProcessIdentifier")) {
      try(ProcessVector[iOni] <- RunProcessIdentifier);
    }
    if (exists("SimulationIdentifier")) {
      try(SimulationVector[iOni] <- SimulationIdentifier);   
    }
    if (exists("FinishesTime") && exists("StartsTime")) {
      try(SecDiff <- FinishesTime$JustSeconds  - StartsTime$JustSeconds);
      try(DayDiff <- 24* 60 * 60 * (FinishesTime$JustDays - StartsTime$JustDays));
      try(AlternateTimeVector[iOni] <-  DayDiff+SecDiff);
    }
    if (!is.null(TemperatureListMatrix)) {
      try(TemperatureListMatrix[iOni,] <- TemperatureList);
    }
  } 
  try(setwd(ASmallContainDir));
  try(dir.create("CollectSimulationsSummary", recursive=TRUE, showWarnings=FALSE));
  
  eval(parse(text=GetG0Text("n", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("p", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("k", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("pGroups", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("kGroups", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("sigma", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("piAPriorData", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("CovRho", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("sigmaPriorData", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("piSV", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("sigSV", "globalenv()", S=1)));  
  
  SigmaNoise = sigma;
  print(paste("-- 012CollectSaveResuls.r: Now to Save S", OnFunction, "Summary.RData", sep="")); flush.console();
  try(colnames(MyOutMatrix)  <- c("l1error", "l2error", "Type1", "Type2", "L.zero", "time:user",
    "time:system", "time:elapsed"));
  try(eval(parse(text=RMAKS())));    
  try(save(AlternateTimeVector = AlternateTimeVector, SMyFiles=SMyFiles,
    MyOutMatrix=MyOutMatrix, ProcessVector=ProcessVector, sigma = sigma,
    SimulationVector=SimulationVector, sigma=sigma, n=n, p=p, k=k,
    pGroups=pGroups, kGroups=kGroups, piAPriorData=piAPriorData, 
    CovRho=CovRho,  sigmaPriorData=sigmaPriorData, piSV=piSV, sigSV=sigSV,
    TemperatureListMatrix=TemperatureListMatrix,
    file=paste("CollectSimulationsSummary//",
      "S", OnFunction, "Summary.RData", sep="")));
      
  CIMyFiles <- MyFiles[substr(MyFiles, 1, nchar(paste("CI", OnFunction, sep=""))) ==  paste("CI", OnFunction, sep="")];
  if (length(CIMyFiles) <= 0) {
    if (verbose >= 1) {
    print(paste("CollapseOnFunction: Woah, No completed CIFiles fo OnFunction = ", OnFunction, " in ", ASmallContainDir, 
      " but length(MyFiles) = ", length(MyFiles), sep=""));
    flush.console();
    }    
  } else {
  load(paste(OnFunction, "//", CIMyFiles[1], sep=""));
  if (!exists("AllCIReturn") || is.null(AllCIReturn) || length(AllCIReturn) <= 0) {
    print(paste("CollapseOnFunction: Woah, MyV is not good for this file ", CIMyFiles[1], sep=""));
    flush.console();
    return(-999);
  }
  AllCIList <- list(); 
  for (idi in 1:length(AllCIReturn)) {
    NameCIList <- names(AllCIReturn[[idi]]);
    NewThing <- list();
    for (kk in 1:length(names(AllCIReturn[[idi]]))) {
      if (is.null(AllCIReturn[[idi]][[kk]]) || length(AllCIReturn[[idi]][[kk]]) <= 1) {
        NewThing[[kk]] <- rep(NA, length(CIMyFiles));
      } else if (length(AllCIReturn[[idi]][[kk]]) >= 2) {
        NewThing[[kk]] <- matrix(NA, length(CIMyFiles), length(AllCIReturn[[idi]][[kk]]));
      }
      names(NewThing)[kk] <- NameCIList[kk];
    }
    AllCIList[[idi]] <- NewThing;
    names(AllCIList)[idi] <- AllCIReturn[[idi]]$NameOfInterval;
  }
  BetaMatrix <- matrix(NA, length(CIMyFiles), length(Betas));
  CITimeMatrix <- matrix(NA, length(CIMyFiles), length(CITime));
  CISavedQuantiles <- CIQuantiles;
  HitVector <- rep(NA, length(CIMyFiles));
  SimulationVector <- rep(NA, length(CIMyFiles));
  ProcessVector <- rep(NA, length(CIMyFiles));
    
  AKS <- c("AllCIReturn", "Betas",  "CIEst", "CIQuantiles", "CITime", "FunctionName",       
    "Hit",  "NameFunction", "RunProcessIdentifier", "SimulationIdentifier");
  RMAKS <- function() {
     Mytxt <- "";
     for (kk in 1:length(AKS)) { 
        Mytxt <- paste(Mytxt,   "
          ", AKS[kk], " <- NULL; rm(",  AKS[kk], ");", sep="")            
     }
     return(Mytxt);
  }          

  eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
  for (iOni in 1:length(CIMyFiles)) {
    if (verbose >= 2 || (verbose >=1 && iOni %% 10 == 0)) {
      print(paste("Loading CIFile ", iOni, "/", length(CIMyFiles), " ",
        CIMyFiles[iOni], sep="")); flush.console();
    }
    try(eval(parse(text=RMAKS())));
    try(load(paste(OnFunction, "//", CIMyFiles[iOni], sep="")));
    if (exists("AllCIReturn")) {
      MyNList <- rep("", length(AllCIReturn));
      for (idi in 1:length(AllCIReturn)) {
        MyNList[idi] <- AllCIReturn[[idi]]$NameOfInterval
      }
      try(names(AllCIList) <- MyNList);
      for (idi in 1:length(AllCIReturn)) {
        NameCIList <- names(AllCIReturn[[idi]]);
        for (kk in 1:length(NameCIList)) {
          if (is.null(AllCIReturn[[idi]][[kk]]) || length(AllCIReturn[[idi]][[kk]]) <= 0) {
            try(AllCIList[[idi]][[kk]][iOni] <- NA);
          } else if (length(AllCIReturn[[idi]][[kk]]) == 1) {
            try(AllCIList[[idi]][[kk]][iOni] <- AllCIReturn[[idi]][[kk]]);            
          } else if (length(AllCIReturn[[idi]][[kk]]) >= 2) {
            try(AllCIList[[idi]][[kk]][iOni,] <- 
              AllCIReturn[[idi]][[kk]][1:NCOL(AllCIList[[idi]][[kk]])]);
          }
        }
        ##try(names(AllCIReturn)[idi] <- AllCIReturn[[idi]]$NameOfInterval
      }
    }
    if (exists("Betas") && !is.null(Betas) && length(Betas) >= 1) {
      try(BetaMatrix[iOni, ] <- Betas[1:NCOL(BetaMatrix)]);  
    }
    if (exists("CITime") && !is.null(CITime) && length(CITime) >= 1) {
      try(CITimeMatrix[iOni, ] <- CITime[1:NCOL(CITimeMatrix)]);  
    }
    if (exists("Hit")) {
      try(HitVector[iOni] <- Hit);
    }
    if (exists("RunProcessIdentifier")) {
      try(ProcessVector[iOni] <- RunProcessIdentifier);
    }
    if (exists("SimulationIdentifier")) {
      try(SimulationVector[iOni] <- SimulationIdentifier);   
    }
  } 
  try(setwd(ASmallContainDir));
  try(eval(parse(text=RMAKS()))); 
  print(paste("--0012CollectSaveResults.r  About to save CI", OnFunction, "Summary.RData", sep=""));
  flush.console();   
  try(save(AllCIList=AllCIList,
    CIMyFiles=CIMyFiles, CITimeMatrix=CITimeMatrix, 
    CISavedQuantiles=CISavedQuantiles,
    BetaMatrix=BetaMatrix, HitVector=HitVector,
    ProcessVector=ProcessVector,
    SimulationVector=SimulationVector, sigma=sigma, n=n, p=p, k=k,
    pGroups=pGroups, kGroups=kGroups, piAPriorData=piAPriorData, 
    CovRho=CovRho,  sigmaPriorData=sigmaPriorData, piSV=piSV, sigSV=sigSV,
    TemperatureListMatrix=TemperatureListMatrix,
    file=paste("CollectSimulationsSummary//",
      "CI", OnFunction, "Summary.RData", sep="")));
  }

  MIPMyFiles <- MyFiles[substr(MyFiles, 1, nchar(paste("MIP", OnFunction, sep=""))) ==  paste("MIP", OnFunction, sep="")];
  if (length(MIPMyFiles) <= 0) {
    if (verbose >= 1) {
    print(paste("CollapseOnFunction: Woah, No completed MIPFiles fo OnFunction = ", OnFunction, " in ", ASmallContainDir, 
      " but length(MyFiles) = ", length(MyFiles), sep=""));
    flush.console();
    }    
  } else {
  load(paste(OnFunction, "//", MIPMyFiles[1], sep=""));
  if (!exists("AllMIPReturn") || is.null(AllMIPReturn) || length(AllMIPReturn) <= 0) {
    print(paste("CollapseOnFunction: Woah, MyV is not good for this file ", MIPMyFiles[1], sep=""));
    flush.console();
    return(-999);
  }
  AllMIPList <- list(); 
  NameMIPList <- names(AllMIPReturn);
    NewThing <- list();
    for (kk in 1:length(NameMIPList)) {
      if (is.null(AllMIPReturn[[kk]]) || length(AllMIPReturn[[kk]]) <= 1) {
        NewThing[[kk]] <- rep(NA, length(MIPMyFiles));
      } else if (length(AllMIPReturn[[kk]]) >= 2) {
        ANN <- names(AllMIPReturn)[kk];
        aLLN <- names(AllMIPReturn[[kk]]);
        MyText <- paste(
        "NewThing[[kk]] <- matrix(NA, length(MIPMyFiles), length(AllMIPReturn[[kk]]),
          dimnames=list(Files=(1:length(MIPMyFiles)),  ", ANN, " = aLLN))", sep="");
        try(eval(parse(text=MyText)));
      }
      ##names(NewThing)[kk] <- NameMIPList[kk];
    }
    names(NewThing) <- NameMIPList;
    AllMIPList <- NewThing;
  SimulationVector <- rep(NA, length(MIPMyFiles));
  ProcessVector <- rep(NA, length(MIPMyFiles));
    
  AKS <- c( "AllMIPReturn", "FunctionName", "NameFunction",        
   "RunProcessIdentifier", "SimulationIdentifier");
  RMAKS <- function() {
     Mytxt <- "";
     for (kk in 1:length(AKS)) { 
        Mytxt <- paste(Mytxt,   "
          ", AKS[kk], " <- NULL; rm(",  AKS[kk], ");", sep="")            
     }
     return(Mytxt);
  }          

  eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
  for (iOni in 1:length(MIPMyFiles)) {
    if (verbose >= 2 || (verbose >=1 && iOni %% 10 == 0)) {
      print(paste("Loading MIPFile ", iOni, "/", length(MIPMyFiles), " ",
        MIPMyFiles[iOni], sep="")); flush.console();
    }
    try(eval(parse(text=RMAKS())));
    try(load(paste(OnFunction, "//", MIPMyFiles[iOni], sep="")));
    if (exists("AllMIPReturn")) {
        NameMIPList <- names(AllMIPReturn);
        for (kk in 1:length(NameMIPList)) {
          if (is.null(AllMIPReturn[[kk]]) || length(AllMIPReturn[[kk]]) <= 0) {
            try(AllMIPList[[kk]][iOni] <- NA);
          } else if (length(AllMIPReturn[[kk]]) == 1) {
            try(AllMIPList[[kk]][iOni] <- AllMIPReturn[[kk]]);            
          } else if (length(AllMIPReturn[[kk]]) >= 2) {
            try(AllMIPList[[kk]][iOni,] <- 
              AllMIPReturn[[kk]][1:NCOL(AllMIPList[[kk]])]);
          }
        }
    }
    if (exists("RunProcessIdentifier")) {
      try(ProcessVector[iOni] <- RunProcessIdentifier);
    }
    if (exists("SimulationIdentifier")) {
      try(SimulationVector[iOni] <- SimulationIdentifier);   
    }
  } 
  try(setwd(ASmallContainDir));
  try(eval(parse(text=RMAKS())));    
  try(save(AllMIPList=AllMIPList,
    MIPMyFiles=MIPMyFiles,
    ProcessVector=ProcessVector,
    SimulationVector=SimulationVector, sigma=sigma, n=n, p=p, k=k,
    pGroups=pGroups, kGroups=kGroups, piAPriorData=piAPriorData, 
    CovRho=CovRho,  sigmaPriorData=sigmaPriorData, piSV=piSV, sigSV=sigSV,
    TemperatureListMatrix=TemperatureListMatrix,
    file=paste("CollectSimulationsSummary//",
      "MIP", OnFunction, "Summary.RData", sep="")));
  }


  BetaCIMyFiles <- MyFiles[substr(MyFiles, 1, nchar(paste("BetaCI", OnFunction, sep=""))) ==  paste("BetaCI", OnFunction, sep="")];
  
  MatrixVars <- c("NameOfInterval", "NameFunction", "Hit",
    "CIQuantiles", "MedianAbsDiffs", "MedianSqDiffs", "TableBetas", "UniqueBetas", "iMIPReport", "BetaStart");
  ArrayVars <-  c("iCoverageQuantiles", "iMeanWidthQuantiles", "iMeanSqWidthQuantiles", "iCINormalLoss", "iCIUnifLoss", "FittedCIs")


  if (length(BetaCIMyFiles) <= 0) {
    if (verbose >= 1) {
    print(paste("CollapseOnFunction: Woah, No completed BetaCIFiles fo OnFunction = ", OnFunction, " in ", ASmallContainDir, 
      " but length(MyFiles) = ", length(MyFiles), ", so we will not do the Betas", sep=""));
    flush.console();
    }    
  } else {
  load(paste(OnFunction, "//", BetaCIMyFiles[1], sep=""));
  if ("iMIPReport" %in% names(AllCIReturn[[1]])) {
    
  } else {
    MatrixVars <- MatrixVars[MatrixVars != "iMIPReport"];
  }
  if (!exists("AllCIReturn") || is.null(AllCIReturn) || length(AllCIReturn) <= 0) {
    print(paste("CollapseOnFunction: Woah, MyV is not good for this file ", CIMyFiles[1], sep=""));
    flush.console();
    return(-999);
  }
  AllCIList <- list(); 
  for (idi in 1:length(AllCIReturn)) {
    NameCIList <- names(AllCIReturn)[idi];
    NewThing <- list();
    OnGo <- 1;
    for (kk in 1:length(MatrixVars)) {
      MyS <- paste("if (is.null(AllCIReturn[[idi]]$", MatrixVars[kk], ")|| 
        length(AllCIReturn[[idi]]$", MatrixVars[kk], ") <= 1) {
        NewThing[[OnGo]] <- rep(NA, length(BetaCIMyFiles));
      } else if (length(AllCIReturn[[idi]]$", MatrixVars[kk], ") >= 2) {
        ColMatVars = names(AllCIReturn[[idi]]$", MatrixVars[kk], ");
        if (is.null(ColMatVars)) {
          ColMatVars <- (1:length(AllCIReturn[[idi]]$", MatrixVars[kk], "));
        }
        if (MatrixVars[kk] == \"iMIPReport\") {
          ColMatVars <- AllCIReturn[[idi]]$UniqueBetas;
        }
        NewThing[[OnGo]] <- matrix(NA, length(BetaCIMyFiles), length(AllCIReturn[[idi]]$", MatrixVars[kk], "),
          dimnames=list(Files=1:length(BetaCIMyFiles), ", MatrixVars[kk], " = ColMatVars));
      }
      names(NewThing)[OnGo] <- MatrixVars[kk];
      OnGo <- OnGo + 1;
      ", sep="");
      eval(parse(text=MyS));
    }
    for (kk in 1:length(ArrayVars)) {
      Dim1Name <- dimnames(AllCIReturn[[idi]][[ArrayVars[kk]]])
      NName <- ArrayVars[kk];
      MyS <- paste("if (is.null(AllCIReturn[[idi]]$", ArrayVars[kk], ") || 
        length(AllCIReturn[[idi]]$", ArrayVars[kk], ") <= 1) {
        MyObject <- rep(NA, length(BetaCIMyFiles));
      } else if (!is.null(dim(AllCIReturn[[idi]]$", ArrayVars[kk], "))) {
        MyObject <- array(data = NA, dim=c(length(BetaCIMyFiles), NROW(AllCIReturn[[idi]]$", ArrayVars[kk], "),
          NCOL(AllCIReturn[[idi]]$", ArrayVars[kk], ")), dimnames=list(Files=(1:length(BetaCIMyFiles)), 
          ", names(Dim1Name)[1], "=Dim1Name[[1]], ", names(Dim1Name)[2], "=Dim1Name[[2]]));
      }  else if (length(AllCIReturn[[idi]]$", ArrayVars[kk], ") >= 2) {
        MyObject <- matrix(NA, length(BetaCIMyFiles), length(AllCIReturn[[idi]]$", ArrayVars[kk], "),
          dimnames=list(Files=(1:length(BetaCIMyFiles)), ", NName,    " = colnames(AllCIReturn[[idi]]$", ArrayVars[kk], ")));
      }
      NewThing[[OnGo]] <- MyObject;
      names(NewThing)[OnGo] <- ArrayVars[kk];
   
      ", sep="");
      eval(parse(text=MyS));
      OnGo <- OnGo + 1;
    }
    AllCIList[[idi]] <- NewThing;
    names(AllCIList)[idi] <- AllCIReturn[[idi]]$NameOfInterval;
  }
  if (exists("Betas")) {
    BetaMatrix <- matrix(NA, length(BetaCIMyFiles), length(Betas));
  } else { BetaMatrix <- NULL; }
  if (exists("BetasOn")) {
     BetasOnMatrix <- matrix(NA, length(BetaCIMyFiles), length(BetasOn));
  }  else { BetasOnMatrix <- NULL; }
  if (exists("IndexOnBetaStart")) {
    IndexOnBetaStartMatrix <- matrix(NA, length(BetaCIMyFiles), length(IndexOnBetaStart));
  }  else {
    IndexOnBetaStartMatrix <- NULL;
  }
  if (exists("BetaStartOn")) {
    BetaStartOnMatrix <- matrix(NA, length(BetaCIMyFiles), length(BetaStartOn));
  }  else {BetaStartOnMatrix <- NULL; }
  if (exists("IndexOnBetas")) {
     IndexOnBetasMatrix <- matrix(NA, length(BetaCIMyFiles), length(IndexOnBetas));
  }  else { IndexOnBetasMatrix <- NULL; }
  if (exists("FitBetaOn")) {
     FitBetaOnMatrix <- matrix(NA, length(BetaCIMyFiles), length(FitBetaOn));
  }  else { FitBetaOnMatrix <- NULL; }
  if (exists("MIPBetaOn")) {
     MIPBetaOnMatrix <- matrix(NA, length(BetaCIMyFiles), length(MIPBetaOn));
  }  else { MIPBetaOnMatrix <- NULL; }
  if (exists("MeanMIPBetaOff")) {
     MeanMIPBetaOffArray <- rep(NA, length(BetaCIMyFiles));
  }  else { MeanMIPBetaOffArray <- NULL; }
  
                   
  CITimeMatrix <- matrix(NA, length(BetaCIMyFiles), length(CITime));
  CISavedQuantiles <- CIQuantiles;
  HitVector <- rep(NA, length(BetaCIMyFiles));
  SimulationVector <- rep(NA, length(BetaCIMyFiles));
  ProcessVector <- rep(NA, length(BetaCIMyFiles));
    
  AKS <- c("AllCIReturn", "Betas",  "CIEst", "CIQuantiles", "CITime", "FunctionName",       
    "Hit",  "NameFunction", "RunProcessIdentifier", "SimulationIdentifier");
  RMAKS <- function() {
     Mytxt <- "";
     for (kk in 1:length(AKS)) { 
        Mytxt <- paste(Mytxt,   "
          ", AKS[kk], " <- NULL; rm(",  AKS[kk], ");", sep="")            
     }
     return(Mytxt);
  }          

  eval(parse(text=GetG0Text("verbose", "globalenv()", S=1)));
  for (iOni in 1:length(BetaCIMyFiles)) {
    if (verbose >= 2 || (verbose >=1 && iOni %% 10 == 0)) {
      print(paste("Loading CIFile ", iOni, "/", length(BetaCIMyFiles), " ",
        BetaCIMyFiles[iOni], sep="")); flush.console();
    }
    try(eval(parse(text=RMAKS())));
    try(load(paste(OnFunction, "//", BetaCIMyFiles[iOni], sep="")));
    if (exists("AllCIReturn")) {
      for (idi in 1:length(AllCIReturn)) {
        NameCIList <- names(AllCIReturn[[idi]]);
        try(names(AllCIList)[idi] <- AllCIReturn[[idi]]$NameOfInterval);
        LN <- names(AllCIList[[idi]]);
        for (kk in 1:length(AllCIList[[idi]])) {
         MyS <- paste("
          if (is.null(AllCIReturn[[idi]]$", LN[kk], ") || length(AllCIReturn[[idi]]$", LN[kk], ") <= 0) {
            try(AllCIList[[idi]][[kk]][iOni] <- NA);
          } else if (!is.null(dim(AllCIReturn[[idi]]$", LN[kk], " )) &&
            NROW(AllCIReturn[[idi]]$", LN[kk], ") > 1 &&
            NCOL(AllCIReturn[[idi]]$",LN[kk], ") > 1) {
            try(AllCIList[[idi]][[kk]][iOni,,] <- AllCIReturn[[idi]]$", LN[kk], ");
          } else if (length(AllCIReturn[[idi]]$", LN[kk], ") == 1) {
            try(AllCIList[[idi]][[kk]][iOni] <- AllCIReturn[[idi]]$", LN[kk], ");            
          } else if (length(AllCIReturn[[idi]]$", LN[kk], ") >= 2) {
            try(AllCIList[[idi]][[kk]][iOni,] <- 
              AllCIReturn[[idi]]$", LN[kk], "[1:NCOL(AllCIList[[idi]][[kk]])]);
          }
          ", sep="");
          eval(parse(text=MyS));
        }
      }
    }
    if (exists("Betas") && !is.null(Betas) && length(Betas) >= 1) {
      try(BetaMatrix[iOni, ] <- Betas[1:NCOL(BetaMatrix)]);  
    }
    if (exists("IndexOnBetaStart") && !is.null(IndexOnBetaStart) && length(IndexOnBetaStart) >= 1) {
      try(IndexOnBetaStartMatrix[iOni,] <- IndexOnBetaStart);
    }
    if (exists("BetaStartOn") && !is.null(BetaStartOn) && length(BetaStartOn) >= 1) {
      try(BetaStartOnMatrix[iOni,] <- BetaStartOn);
    }
    if (exists("BetasOn") && !is.null(BetasOn) && length(BetasOn) == NCOL(BetasOnMatrix)) {
      try(BetasOnMatrix[iOni, ] <- BetasOn);  
    }
    if (exists("IndexOnBetas") && !is.null(IndexOnBetas) && length(IndexOnBetas) == NCOL(IndexOnBetasMatrix)) {
      try(IndexOnBetasMatrix[iOni, ] <- IndexOnBetas);  
    }
    if (exists("FitBetaOn") && !is.null(FitBetaOn) && length(FitBetaOn) == NCOL(FitBetaOnMatrix)) {
      try(FitBetaOnMatrix[iOni, ] <- FitBetaOn);  
    }
    if (exists("MIPBetaOn") && !is.null(MIPBetaOn) && length(MIPBetaOn) == NCOL(MIPBetaOnMatrix)) {
      try(MIPBetaOnMatrix[iOni, ] <- MIPBetaOn);  
    }
    if (exists("MeanMIPBetaOff") && !is.null(MeanMIPBetaOff) && length(MeanMIPBetaOff) == 1) {
      try(MeanMIPBetaOffArray[iOni] <- MeanMIPBetaOff);
    }
    if (exists("CITime") && !is.null(CITime) && length(CITime) >= 1) {
      try(CITimeMatrix[iOni, ] <- CITime[1:NCOL(CITimeMatrix)]);  
    }
    if (exists("Hit")) {
      try(HitVector[iOni] <- Hit);
    }
    if (exists("RunProcessIdentifier")) {
      try(ProcessVector[iOni] <- RunProcessIdentifier);
    }
    if (exists("SimulationIdentifier")) {
      try(SimulationVector[iOni] <- SimulationIdentifier);   
    }
  } 
  try(setwd(ASmallContainDir));
  try(eval(parse(text=RMAKS())));    
  try(save(AllCIList=AllCIList,
    MIPBetaOnMatrix = MIPBetaOnMatrix, FitBetaOnMatrix=FitBetaOnMatrix,
    IndexOnBetasMatrix = IndexOnBetasMatrix, BetasOnMatrix=BetasOnMatrix,
    MeanMIPBetaOffArray = MeanMIPBetaOffArray,  
    IndexOnBetaStartMatrix = IndexOnBetaStartMatrix,
    BetaStartOnMatrix = BetaStartOnMatrix,
    BetaCIMyFiles=BetaCIMyFiles, CITimeMatrix=CITimeMatrix, 
    CISavedQuantiles=CISavedQuantiles,
    BetaMatrix=BetaMatrix, HitVector=HitVector,
    ProcessVector=ProcessVector,
    SimulationVector=SimulationVector, sigma=sigma, n=n, p=p, k=k,
    pGroups=pGroups, kGroups=kGroups, piAPriorData=piAPriorData, 
    CovRho=CovRho,  sigmaPriorData=sigmaPriorData, piSV=piSV, sigSV=sigSV,
    TemperatureListMatrix=TemperatureListMatrix,
    file=paste("CollectSimulationsSummary//",
      "BetaCI", OnFunction, "Summary.RData", sep="")));
  }


  try(setwd(Oldwd))    
}
