  try(library(methods, warn.conflicts = FALSE, quietly=TRUE));
  try(library(R.methodsS3, warn.conflicts=FALSE, quietly=TRUE)); 
  try(library(R.oo, warn.conflicts=FALSE, quietly=TRUE));

MAXDIMXKeep  = 4;

if (FALSE) {
  library(AlanDirectories); SetSaveHomes();
  source(paste(DropHome, "//TwoLassoPackage//TwoLassoCpp//R//TwoLassoSexp.R", sep=""))
}
TwoLassoR5 = setRefClass("TwoLassoR5",
      fields = list(
     X = "matrix", Y="vector", XtX = "matrix",
     .XtY = "vector", XTY = function(AA = NULL) {
       if (!is.null(AA)) {
         if (length(AA) == .self$pCoef) {
           try(.self$.XtY <- AA);  return;
         } else {
           print(paste("TwoLassoR5: Attempt to set XtY in Error, ",
             "Give vector length: ", .self$pCoef)); return(-1);
         }
       }
       return(.self$.XtY)
     },
     XtY = function(AA = NULL) {
       if (!is.null(AA)) {
         if (length(AA) == .self$pCoef) {
           .self$.XtY <- AA;  return;
         } else {
           print(paste("TwoLassoR5: Attempt to set XtY in Error, ",
             "Give vector length: ", .self$pCoef)); return(-1);
         }
       }
       return(.self$.XtY)
     },
     OnBeta = "vector", BBOn1 = "vector",
     ReturnBetas = "vector", StartBeta = "vector",
     InverseGammaConstant="numeric",
     TDFNu="numeric",
     StandardFlag = "numeric",
     CauchyEpsilon="numeric", MaxCauchy="numeric",
     pCoef = "numeric", nSample = "numeric",
     .tt1 = "integer", .tt2 = "integer", 
     tt1 = function(att1=NULL) {
       if (!is.null(att1) && is.numeric(att1) && att1 >= 1) {
         try(.self$.tt1 <- as.integer(att1-1));
         return;
       }
       return(.self$.tt1+1);
     },
     tt2 = function(att2 = NULL) {
       if (!is.null(att2) && is.numeric(att2) && att2 >= 1) {
         try(.self$.tt2 <- as.integer(att2-1));
         return;
       }
       return(.self$.tt2+1);
     },
     WLSWeights= "numeric", 
     .Verbose="numeric",
     Verbose = function(AT = NULL) {
       if (is.null(AT)) { return(.self$.Verbose); }
       try(.self$.Verbose <- as.numeric(AT));
       if (!is.null(.self$.TCSPointer)) {
         try(.Call("setVerbose", .self$.TCSPointer, as.integer(AT)));
       }
       return;
     },
     LambdaAK = "vector", 
     LambdaDK="vector", OrderSeq= "vector",
     PiA = "numeric", SigmaSq = "numeric",
     PiAVector= "vector", SigmaVector="vector", PiAPrior="vector",
     SigmaPrior="vector", InitKKs="numeric",       
     sdX="vector", sdY = "numeric",
     GroupSexp="numeric", OnGammas = "numeric",
     RecordFlag = "numeric",
     RecBBOn1 = "matrix", USOnBeta = "vector",
     RecOnBeta = "matrix", USRecOnBeta = "vector",
     TCpp = "ANY", FailureFlag = "numeric", 
     RecordBetaCVFinish = "matrix",  StartBetaMatrix = "matrix",
     SigmaVectorInputs ="vector", PiAVectorInputs = "vector",
     LambdaAKInputs="vector", LambdaDKInputs = "vector",
     L2ShrinkagePrior = "vector", L2ShrinkageRecords = "vector",
     YSq = "numeric", .TCSPointer = "ANY", 
     DoLogit = "logical", sDoLogit = function(iIn) {
       if (.self$DoLogit == TRUE)  {return(1); }
       return(0);
     },
     ConfidenceMatrix = function(aN = NULL) {
       if (!is.null(aN)) {
         print("getConfidenceMatrix: cannot be set, only get!"); return(NULL);
       }
       if (is.null(.self$.TCSPointer) || is.numeric(.self$.TCSPointer)) {
         print("getConfidenceMatrix: Sorry only operable if TCSPointer is set!");
         return(NULL);
       }
       AT <- NULL;
       try(AT <- .Call("getConfidenceMatrix", .self$.TCSPointer));
       if (is.null(AT)) {
         print("Sorry, TCSPointer: Did not get ConfidenceMatrix"); flush.console();
         return(NULL);
       }
       return(AT);
     }, 
     CredibilityMatrix = function(aN = NULL) {
       if (!is.null(aN)) {
         print("getConfidenceMatrix: cannot be set, only get!"); return(NULL);
       }
       if (is.null(.self$.TCSPointer) || is.numeric(.self$.TCSPointer)) {
         print("getConfidenceMatrix: Sorry only operable if TCSPointer is set!");
         return(NULL);
       }
       AT <- NULL;
       try(AT <- .Call("getConfidenceMatrix", .self$.TCSPointer));
       if (is.null(AT)) {
         print("Sorry, TCSPointer: Did not get ConfidenceMatrix"); flush.console();
         return(NULL);
       }
       return(AT);
     }, 
     UnshrunkConfidenceMatrix = function(aN = NULL) {
       if (!is.null(aN)) {
         print("getUnshrunkConfidenceMatrix: cannot be set, only get!"); return(NULL);
       }
       if (is.null(.self$.TCSPointer) || is.numeric(.self$.TCSPointer)) {
         print("getUnshrunkConfidenceMatrix: Sorry only operable if TCSPointer is set!");
         return(NULL);
       }
       AT <- NULL;
       try(AT <- .Call("getUnshrunkConfidenceMatrix", .self$.TCSPointer));
       if (is.null(AT)) {
         print("Sorry, (TwoLassoR5) TCSPointer: Did not get UnshrunkConfidenceMatrix"); flush.console();
         return(NULL);
       }
       return(AT);
     },
     UnshrunkCredibilityMatrix = function(aN = NULL) {
       if (!is.null(aN)) {
         print("getUnshrunkConfidenceMatrix: cannot be set, only get!"); return(NULL);
       }
       if (is.null(.self$.TCSPointer) || is.numeric(.self$.TCSPointer)) {
         print("getUnshrunkConfidenceMatrix: Sorry only operable if TCSPointer is set!");
         return(NULL);
       }
       AT <- NULL;
       try(AT <- .Call("getUnshrunkConfidenceMatrix", .self$.TCSPointer));
       if (is.null(AT)) {
         print("Sorry, (TwoLassoR5) TCSPointer: Did not get UnshrunkConfidenceMatrix"); flush.console();
         return(NULL);
       }
       return(AT);
     },
     ConfidenceQuantiles = function(aN = NULL) {
       if (!is.null(aN)) {
         if (is.null(.self$.TCSPointer)) {
           print("SetupConfidenceQuantiles(TwoLassoR5): can set, if TCS pointer is NULL"); return(NULL);
         }
         AT <- NULL;
         try(AT <- .Call("SetupConfidenceQuantiles", .self$.TCSPointer, aN));
         if (is.null(AT) || (is.numeric(AT) && AT[1] < 0)) {
           print("setConfiddenceQuantiles(TwoLassoR5): set error. "); return(NULL);
         }
         return(AT);
       }
       if (is.null(.self$.TCSPointer) || is.numeric(.self$.TCSPointer)) {
         print("setConfidenceQuantiles: Sorry only operable if TCSPointer is set!");
         return(NULL);
       }
       AT <- NULL;
       try(AT <- .Call("getConfidenceQuantiles", .self$.TCSPointer));
       if (is.null(AT)) {
         print("Sorry, (TwoLassoR5) TCSPointer: Did not get ConfidenceQuantiles"); flush.console();
         return(NULL);
       }
       return(AT);
     },
     LambdaIndexForConfidenceIntervals = function(aN = NULL) {
       if (!is.null(aN)) {
         if (is.null(.self$.TCSPointer)) {
           print("setLambdaIndexForConfidenceIntervals(TwoLassoR5): can set, if TCS pointer is NULL"); return(NULL);
         }
         AT <- NULL;
         try(AT <- .Call("setLambdaIndexForConfidenceIntervals", .self$.TCSPointer, aN));
         if (is.null(AT) || (is.numeric(AT) && AT[1] < 0)) {
           print("setLambdaIndexForConfidenceIntervals(TwoLassoR5): set error. "); return(NULL);
         }
         return(AT);
       }
       if (is.null(.self$.TCSPointer) || is.numeric(.self$.TCSPointer)) {
         print("setLambdaIndexForConfidenceIntervals: Sorry only operable if TCSPointer is set!");
         return(NULL);
       }
       AT <- NULL;
       try(AT <- .Call("getLambdaIndexForConfidenceIntervals", .self$.TCSPointer));
       if (is.null(AT)) {
         print("Sorry, (TwoLassoR5) TCSPointer: Did not get LambdaIndexForConfidenceIntervals"); flush.console();
         return(NULL);
       }
       return(AT);
     },
     MyCFold = "ANY",
     RunFlag = "numeric", FixKa = "integer"), 
     methods = list(nSample = function() { length(Y) },
       getNames = function()  {
         return(c("Y", "X", "XtX", "XtY", "InverseGammaConstant",
           "TDFNu", "nSample", "StartBeta", "RecordFlag",
           "StandardFlag", "SigmaSq", "PiA", "RecBBOn1",
           "RecOnBeta", "USOnBeta", "USRecOnBeta", 
           "PiAVector", "SigmaAVector"));
       },
       initialize = function(Y=NULL,X=NULL, XtX = NULL, XtY = NULL,
       InverseGammaConstant = 1, TDFNu = -1, nSample = -1, 
       DoLogit = FALSE,
       StartBeta = NULL, RecordFlag = 0, StandardFlag = 0, 
       SigmaSq = 1, PiA = .5,
       LambdaAK = rep(1,2), LambdaDK = rep(1,2), OrderSeq = rep(3,2), 
       LambdaAKInputs = NULL, LambdaDKInputs = NULL, 
       RunFlag = 0,
       GroupSEXP = NULL, 
       CauchyEpsilon = .0001, MaxCauchy = 40,
       PiAVector = c(.5,.5), SigmaPrior = NULL, PiAPrior = NULL,  Verbose = 0, 
       InitKKs = 20, PiAVectorInputs = NULL, SigmaVectorInputs = NULL,
       StartBetaMatrix = NULL, FixKa = -100,
       MakeCpp = FALSE, YSq = NULL, L2ShrinkagePrior = NULL, 
       RecordL2Shrinkage = FALSE,
       ...) {
      try(.self$MyCFold <- NULL);
      
      try(.self$StandardFlag <- as.numeric(StandardFlag));
      try(.self$MaxCauchy <- as.numeric(MaxCauchy));  
      try(.self$CauchyEpsilon <- as.numeric(CauchyEpsilon));
      try(.self$RunFlag <- as.numeric(RunFlag));
      if (length(SigmaPrior) <= 1) { try(.self$SigmaPrior <- NULL);
      } else { .self$SigmaPrior <- SigmaPrior[1:2]; }
      if (length(PiAPrior) <= 1) { try(.self$PiAPrior <- NULL); 
      } else {try(.self$PiAPrior <- PiAPrior[1:2]);}
      try(.self$FailureFlag <- 0);
      try(.self$Verbose <- as.integer(Verbose) );
      try(.self$InitKKs <- as.integer(InitKKs)); 
      if (Verbose > 0) {
        print("Started TwoLassoR5");flush.console(); 
      }       
      if (!is.null(PiAPrior) && length(PiAPrior) >= 2) {
        .self$PiAPrior <- PiAPrior[1:2];
      }  else { try(.self$PiAPrior <- NULL); }
 
      if (Verbose > 0) {
        print("TwoLassoR5: Filled PiAPrior, setting XX")
      }
      if (Verbose > 0) {
        print(paste("TwoLassoR5: TDFNu set to ", TDFNu, sep="")); flush.console();
      }
      try(.self$DoLogit <- as.logical(DoLogit));
      .self$TDFNu <- TDFNu;
      try(.self$.TCSPointer <- NULL);
      if (TDFNu <= 0 && !is.null(X) && dim(X)[2] < MAXDIMXKeep) {
        .self$sdX <- apply(X,2,sd); 
        if (StandardFlag == 0) {
          if (Verbose > 0) {
            print("TwoLassoR5: Converting to XtX, XtY"); flush.console();
          }
          .self$XtX <<- t(X) %*% X;
          .self$.XtY <- as.vector(t(X) %*% Y);
          eval(parse(text=SetGText("nSample", "globalenv()", S=1)));
          .self$nSample <- length(Y);
          eval(parse(text=SetGText("nSample", "globalenv()", S=1)));          
          .self$pCoef <- length(X[1,]);
          .self$X <- X; .self$Y <- Y; .self$YSq <- sum(Y^2);
        } else {
          if (Verbose > 0) {
            print("Standardize, X, Y, for permanent use"); flush.console();
          }
          .self$X <- StandardizeXX(X);
          .self$Y <- StandardizeXX(Y);
          .self$YSq <- 1; 
        }
      } else if (TDFNu  <= 0 && !is.null(XtX) && !is.null(XtY) &&
        !is.null(dim(XtX)) && length(dim(XtX)) == 2 
        && sum(abs(dim(XtX))) > 0) {
        if (nSample < 0  && is.null(Y)) {
          tryCatch("Set TwoLassoR5 Error, nSample < 0, but no Y given");
        }
        if (Verbose > 0) {
          print("TwoLassoR5: About to use XtX, XtY"); flush.console(); 
          print(paste(" -- Dim XtX = (", 
            paste(dim(XtX), collapse=", "), ") and ",
            " length(XtY) = ", length(XtY), sep="")); flush.console();
        }
        .self$sdX <- diag(XtX); .self$sdY <- 1;
        if (StandardFlag == 0) {
          .self$XtX <- XtX;  .self$.XtY <- XtY; .self$nSample <- nSample;
        } else {
          .self$XtX <- XtX / sqrt( diag(XtX) %*% t(diag(XtX)));  
          .self$.XtY <- XtY / sqrt(diag(XtX));
        }
        pCoef <<- length(XtX[1,]);      
        YSq <<- YSq;
      } else {
        if (Verbose > 0) {
          print("TwoLassoR5: Default X, Y, else"); flush.console();
        }
        X <<- X; Y <<- Y; YSq <<- sum(Y^2);
        pCoef <<- length(X[1,]);  .self$nSample <- length(Y);
      }
      if (Verbose > 0) {
        print("TwoLassoR5: Set X,Y, XtY"); flush.console();
      }
      if (StandardFlag > 0) {
        USRecOnBeta <<- rep(0, pCoef);
      }
      if (Verbose > 0) {
        print("TwoLassoR5: Started TDFNu");flush.console(); 
      } 
      if (TDFNu > 0) {
        if (is.null(X)) {
          tryCatch("Set TwoLassoR5: We don't do T Noise without X, Y");
        }
        WLSWeights <<- rep(1, pCoef);
      }
      print("TwoLassoR5:  Setting InverseGammaConsant"); flush.console();
      if (!is.null(InverseGammaConstant)) {
        InverseGammaConstant <<- InverseGammaConstant;
      }
      ##if (!is.null(TDFNu)) { TDFNu <<- TDFNu;  }
      TCpp <<- NULL;
      print(" Testing FailureFlag"); flush.console();
      FailureFlag <<- 0;  
      ##.self$lock(FailureFlag, pCoefs,nSample);
      ##pCoefs$lock();  nSample$lock();
      print(" Filling StartBeta")
      if (!is.null(StartBeta) && length(StartBeta) == pCoef) {
        OnBeta <<- StartBeta;
      }  else { OnBeta <<- rep(0, pCoef); }
      ##print(" Finishing Initialize")

      if (Verbose > 0) {
        print("Setting LambdaDK"); flush.console();
      }
      LambdaDK <<- LambdaDK;
      if (length(LambdaAK) < length(LambdaDK))  {
        LambdaAK <<- c(LambdaAK, rep(LambdaAK[length(LambaAK)], length(LambdaDK)- length(LambdaAK)));
      }  else if (length(LambdaAK) > length(LambdaDK)) {
        LambdaAK <<- LambdaAK[1:length(LambdaDK)];
      } else { LambdaAK <<- LambdaAK; }
      if (length(OrderSeq) < length(LambdaDK))  {
        OrderSeq <<- c(OrderSeq, rep(OrderSeq[length(OrderSeq)], length(LambdaDK)- length(OrderSeq)));
      }  else if (length(OrderSeq) > length(LambdaDK)) {
        .self$OrderSeq <- OrderSeq[1:length(LambdaDK)];
      } else { .self$OrderSeq <- OrderSeq; }
      if (Verbose > 0) {
        print("Setting PiAVector"); flush.console();
      }
      
      try(.self$FixKa <- as.integer(FixKa));

      if (!exists("PiAVector") ) {
        if ( (!exists("PiAVector")) || !is.null(PiA)) {
          try( .self$PiAVector <- rep(PiA[1],length(OrderSeq)) );
          try(.self$PiA <- as.numeric(PiA[1]));  
        } else {
          try(.self$PiAVector <- as.numeric( rep(.5, length(OrderSeq)) ) );
          try(.self$PiA <- as.numeric(.5) );
        }
      } else if (!exists("PiAVector") || is.null(PiAVector) || length(PiAVector) == 0) {
        if (!is.null(PiA)) {
          try(.self$PiAVector <- rep(PiA[1],length(OrderSeq)) );
          try(.self$PiA <- as.numeric(PiA[1]));
        } else {
          try(.self$PiAVector <- as.numeric(rep(.5, length(OrderSeq))));
          try(.self$PiA <- as.numeric(.5)); 
        }
      } else { 
        if (length(PiAVector) < length(OrderSeq)) {
          try(.self$PiAVector <- c(PiAVector, rep(PiAVector[length(PiAVector)],
            length(OrderSeq) -length(PiAVector))) );
          try(.self$PiA <- as.numeric(PiAVector[1]) );
        } else {
           try(.self$PiAVector <- as.numeric(PiAVector[1:length(OrderSeq)]));
           try(.self$PiA <- as.numeric(PiAVector[1])); 
        }
      }
      if (Verbose > 0)  {
        print("Filling BBOn1"); flush.console();
      }
      BBOn1 <<- rep(PiAVector[1], pCoef);
      if (Verbose > 0) {
        print("Setting SigmaVector"); flush.console();
      }
      if (is.null(SigmaVector) || length(SigmaVector) == 0) {
        if (!exists("SigmaSq")) { SigmaSq <<- 1; } else {
          try(.self$SigmaSq  <- as.numeric(SigmaSq));
        }
        if (!is.null(SigmaSq)) {
          SigmaVector <<- rep(.self$SigmaSq[1],length(OrderSeq));
        } else {
          SigmaVector <<- rep(1, length(OrderSeq));
        }
      } else { 
        if (length(SigmaVector) < length(OrderSeq)) {
          SigmaVector <<- c(SigmaVector, rep(SigmaVector[length(SigmaVector)],
            length(OrderSeq) -length(SigmaVector)));
        } else {
           SigmaVector <<- SigmaVector[1:length(OrderSeq)]; 
        }
        try(.self$SigmaSq <- as.numeric(SigmaVector[1]) );
      }
      if (Verbose > 0) {
        print("Setting LambdaAKInputs, PiAVectorInputs, SigmaVector"); flush.console();
      }
      if (length(PiAVectorInputs) <= 1 && length(SigmaVectorInputs) > 1) {
        if (is.null(PiAVectorInputs) || length(PiAVectorInputs) < 1) {
          MyT <- "PiAVectorInputs <- rep(.self$PiA, length(SigmaVectorInputs))";
        } else if (length(PiAVectorInputs) == 1) {
          MyT <- "PiAVectorInputs <- rep(PiAVectorInputs[1], length(SigmaVectorInputs))";
        }
        eval(parse(text=MyT));
      } else if (  length(SigmaVectorInputs) > 1 && 
        length(PiAVectorInputs) < length(SigmaVectorInputs) ) {
          MyT <- "PiAVectorInputs <- c(PiAVectorInputs,  
            rep(PiAVectorInputs[length(PiAVectorInputs)],
              length(SigmaVectorInputs) - length(PiAVectorInputs)) )";
        eval(parse(text=MyT));
      }
      if (length(SigmaVectorInputs) <= 1 && length(SigmaVectorInputs) > 1) {
        if (is.null(SigmaVectorInputs) || length(SigmaVectorInputs) < 1) {
          MyT <- "SigmaVectorInputs <- rep(.self$SigmaSq, length(PiAVectorInputs))";
        } else if (length(SigmaVectorInputs) == 1) {
          MyT <- "SigmaVectorInputs <- rep(SigmaVectorInputs[1], length(PiAVectorInputs))";
        }
        eval(parse(text=MyT));
      }  else if (  length(PiAVectorInputs) > 1 && 
        length(SigmaVectorInputs) < length(PiAVectorInputs) ) {
          MyT <- "SigmaVectorInputs <- c(SigmaVectorInputs,  
            rep(SigmaVectorInputs[length(SigmaVectorInputs)],
              length(PiAVectorInputs) - length(SigmaVectorInputs)) )";
        eval(parse(text=MyT));
      }
      if (!is.null(PiAVectorInputs) && length(PiAVectorInputs) >= 1) {
        try(.self$PiAVectorInputs <- PiAVectorInputs);
      } else {
        try(.self$PiAVectorInputs <- vector("numeric",0));
      }
      if (!is.null(SigmaVectorInputs) && length(SigmaVectorInputs) >= 1) {
        try(.self$SigmaVectorInputs <- SigmaVectorInputs);      
      } else {
        try(.self$SigmaVectorInputs <- vector("numeric", 0));    
      }

      .self$.TCSPointer <- NULL;
      if (length(PiAVectorInputs) != length(SigmaVectorInputs)) {
        tryCatch("PiAVectorInputs and SimgaVectorInputs need to be relengthed");
      }
      if (!is.null(LambdaDKInputs) && length(LambdaDKInputs) > 0) { .self$LambdaDKInputs <- LambdaDKInputs; }
      if (!is.null(LambdaAKInputs) && length(LambdaAKInputs) > 0) { .self$LambdaAKInputs <- LambdaAKInputs; }
      if (!is.null(PiAVectorInputs)  && length(PiAVectorInputs) > 0) { 
        .self$RecordBetaCVFinish <- matrix(0, length(PiAVectorInputs), pCoef); 
      }
      if (!is.null(StartBetaMatrix)) { .self$StartBetaMatrix <- StartBetaMatrix; }

      if (Verbose > 0) {
        print(paste("Dim RecordBetaCVFinish = (", 
          paste(dim(RecordBetaCVFinish), collapse=", "), ")", sep=""));
        print(paste(" length SigmaInputs = ", length(.self$SigmaVectorInputs),
          "  and lengthPiAInputs = ", length(.self$PiAVectorInputs), sep=""));
        flush.console();
      }
      MaxLen = length(.self$PiAVectorInputs);
      if (length(.self$SigmaVectorInputs) > MaxLen) {
        MaxLen =  length(.self$SigmaVectorInputs);
      }
      if (Verbose > 0)  {
        print(paste("  MaxLen = ", MaxLen, sep="")); flush.console();
      }
      if (MaxLen > 0) {
        RecordBetaCVFinish <<- matrix(0, pCoef, MaxLen);      
      } 

      
      if (Verbose > 0) {
        print("L2ShrinkagePrior Setup"); flush.console();
      }
      if (!is.null(L2ShrinkagePrior)) {
        L2ShrinkagePrior <<- L2ShrinkagePrior;
        L2ShrinkageRecords <<- rep(1.1, length(OrderSeq));      
      }
      OnGammas <<-  BBOn1 * LambdaAK[1] + (1.0 - BBOn1) * LambdaDK[1]; 
      
      if (Verbose >1) {
        print("Set Sigma, PiA")
      }
              
      if (Verbose > 0) {
        print("Setting RecordFlag info"); flush.console();
      }

      try(.self$tt1 <- as.integer(0));  try(.self$tt2 <- as.integer(0));
      if (MakeCpp == TRUE) {
        if (Verbose > 0) {
          print("TwoLassoR5: Setting and creating TwoLassoSexp Cpp class");
            flush.console();
        }
        TwoLassoSexp <- modTwoLassoSexpCL$TwoLassoSexp;
        ListLambda <- list(.self$LambdaAK, .self$LambdaDK, .self$OrderSeq);
        if (!is.null(XtX)  && sDoLogit == 0) {
          ListFlags = c(.self$Verbose, -n, .self$PiA, .self$FixKa, OnGamma1, 
            .self$Sigma, .self$InitKKs, RecordFlag,
            CauchyEpsilon,  MaxCauchy, .self$sDoLogit);
          ATLS <- new( TwoLassoSexp, -n, .self$XtX, .self$XtY, 
            .self$Verbose, .self$InitKKs, ListLambda)   
          if (!is.null(YSq)) {
            ATLS$GiveYYSq(YSq);
          }         
        } else {
          ListFlags = c(.self$Verbose, n, .self$PiA, FixKa, OnGamma1, 
            .self$Sigma, InitKKs, RecordFlag,
            CauchyEpsilon,  MaxCauchy, .self$sDoLogit);
          ATLS <- new( TwoLassoSexp, n, .self$X, .self$Y, 
            .self$Verbose, .self$InitKKs, ListLambda)   
          ATLS$SetupSumYYSq();                     
        }
        if (ATLS$SuccessFlag <= 0) {
          print(paste("Hey TwoLassoSexp Setup, ATLS SuccessFlag is ", ATLS$SuccessFlag,
            " we had a failure!", sep=""));
          eval(parse(text=SetGText("ATLS", envir="globalenv()", S=1)));
          return(ATLS);
        }
        if ( (!is.numeric(RecordFlag) && RecordFlag == TRUE) || 
          (is.numeric(RecordFlag) && RecordFlag > 0)) {
          .self$RecBBOn1 <- ATLS$sRecBBOn1;
          .self$RecOnBeta <- ATLS$sRecOnBeta;
          if (StandardFlag > 0) {     
             .self$USRecOnBeta <- ATLS$sRecOnBeta;        
          }
        }
        ##if (Verbose > 1) {
        ##    print("ATLS: Setup Lambda"); flush.console();
        ##}
       
        ##ATLS$SetupLambda(.self$SLambdaAK, .self$DLambdaAK, .self$OrderSeq);
        if (!is.null(.self$PiAPrior)) {
          if (Verbose > 1) {
            print("ATLS: Setup PiA Prior"); flush.console();
          }
          ATLS$SetupPiAPrior(.self$PiAPrior);
        }
        if (!is.null(.self$SigmaPrior)) {
          if (Verbose > 1) {
            print("ATLS: Setup Sigma Prior"); flush.console();
          }
          ATLS$SetupSigmaPrior(.self$SigmaPrior);
        }
        if (Verbose > 1) { 
          print("ATLS: Setup Caucby"); flush.console();
        }
        eval(parse(text=SetGText("ATLS", envir="globalenv()", S=1)));
        ATLS$SetupCauchy(.sef$MaxCauchy, .self$CauchyEpsilon)
        if (Verbose > 1) { 
          print("ATLS: Setup tt"); flush.console();
        }
        ATLS$Setuptt(.self$.tt1, .self$.tt2)
        
        if (!is.null(.self$RecBBOn1)) {
          if (Verbose > 1) {
            print("Setup RecBBOn1 "); flush.console();     
          }
          ATLS$SetupRecords(.self$RecOnBeta, .self$RecBBOn1);
        }
        if (Verbose > 1) { 
          print("ATLS: Setup L2 Shinkage"); flush.console();
        }
        ATLS$SetupL2Shrinkage(L2ShrinkagePrior, L2ShrinkageRecords);
        if (.self$FixKa > 0) {
          ATLS$FixKa = .self$FixKa;
        }
        TCpp <<- ATLS;    
      }
      .self
    }
  )
);
  
  
TwoLassoR5$lock("pCoef", "nSample", "X", "Y", "XtX",
  "LambdaAK", "LambdaDK", "OrderSeq", "TDFNu", "LambdaDKInputs", "LambdaAKInputs")   



setConstructorS3("TwoLassoOb", function(
  RetBBOn1 = NULL, nSample = 0, pCoefs =0,
  BBOn1=NULL, RecOnBeta=NULL, RecBBOn1=NULL, USBeta=NULL, 
  USRecOnBeta=NULL, X=NULL, Y=NULL, XtX = NULL, XtY = NULL,
  tt1 = NULL, tt2 = NULL,
  LambdaAK = NULL, LambdaDK=NULL, OrderSeq=NULL,
  PiAVector=NULL, SigmaVector=NULL, PiAPrior=NULL,
  SigmaPrior=NULL, InitKKs=NULL, 
  OnGammas = NULL, StandardFlag = 0,
  InverseGammaConstant=-1,
  CauchyEpsilon=.00001, MaxCauchy=200,
  TDFNu=-1, WLSWeights=NULL, Verbose=0,
  GroupSexp=NULL, OnBeta = 0, sdX = NULL, sdY = NULL,
  LambdaDKInputs = NULL, LambdaAKInputs = NULL,
  PiAVectorInputs=NULL, SigmaVectorInputs=NULL,
  TCSPointer = NULL, FailureFlag = 0, 
  L2ShrinkagePrior = NULL, L2ShrinkageRecords = NULL, Type = "TwoLasso", ...) {
  if (Verbose > 0) {
    print(paste("TwoLassoOb: Creating with Verbose = ", Verbose, sep=""));
    flush.console();
  }
  
  ListOfPrivateData <- c("ReturnBetas", "USOnBeta", 
    "USReturnBetas", "OnBeta", "BBOn1",
    "OnGammas", "X", "Y", "XtX", "XtY",
    "tt1", "tt2", "LambdaAK", "LambdaDK", "OrderSeq",
    "sdX", "sdY",
    "PiAPrior", "SigmaPrior", "PiAVector", "SigmaVector",
    "GroupSexp", "RecOnBeta", "RecBBOn1", "FailureFlag",
    "Type","WLSWeights", "Verbose", "pCoefs", "nSample",
    "MaxCauchy", "CauchyEpsilon", "LambdaAKInputs", "LambdaDKInputs",
    "PiAVectorInputs", "SigmaVectorInputs", "L2ShrinkagePrior",
    "L2ShrinkageRecords"
    )
  FailureFlag = NULL;
  RT =  extend(Object(), "TwoLassoOb",     
     InverseGammaConstant=InverseGammaConstant,
     TDFNu=TDFNu, ListOfPrivateData = ListOfPrivateData,
     StandardFlag = StandardFlag,
     .CauchyEpsilon=CauchyEpsilon, .MaxCauchy=MaxCauchy,
     .nSample = nSample, .pCoefs = pCoefs,
     .tt1 = tt1, .tt2 = tt2, .WLSWeights=WLSWeights, .Verbose=Verbose,
     .LambdaAK = LambdaAK, .LambdaDK=LambdaDK, .OrderSeq=OrderSeq,
     .PiAVector=PiAVector, .SigmaVector=SigmaVector, .PiAPrior=PiAPrior,
     .SigmaPrior=SigmaPrior, InitKKs=InitKKs,       
     .sdX=sdX, .sdY = sdY,
     .GroupSexp=GroupSexp, .OnGammas = OnGammas,
     .X=X, .Y=Y, .XtX = XtX, .XtY = XtY,
     .OnBeta = OnBeta, .BBOn1=BBOn1,
     .RecOnBeta=RecOnBeta, .RecBBOn1=RecBBOn1, .USOnBeta = NULL,
     .USRecOnBeta = NULL,
     .TCSPointer = TCSPointer, .FailureFlag = FailureFlag, .Type = Type,
     .ReturnBetas = NULL, ReturnBetas = NULL, .USRecOnBeta = NULL,
     .LambdaAKInputs = LambdaAKInputs, .LambdaDKInputs = LambdaDKInputs,
     .SigmaVectorInputs =SigmaVectorInputs, .PiAVectorInputs = PiAVectorInputs,
     .L2ShrinkagePrior = L2ShrinkagePrior, .L2ShrinkageRecords = NULL)
  return(RT); 
});

setMethodS3("getMaxCauchy", "TwoLassoOb", function(this,...) {
  return(this$.MaxCauchy);
});
setMethodS3("getCauchyEpsilon", "TwoLassoOb", function(this,...) {
  return(this$.CauchyEpsilon);
});
setMethodS3("getTCSPointer", "TwoLassoOb", function(this,...) {
  print("No: TCSPointer is not available."); flush.console();
  if (is.null(.TCSPointer[1])) {
    print("  Anyway: TCSPointer is Null "); flush.console();
  }
  return(-1);
});
setMethodS3("getType", "TwoLassoOb", function(this,...) {
  return(this$.Type);
});
setMethodS3("getL2ShrinkagePrior", "TwoLassoOb", function(this,...) {
  return(this$.L2ShrinkagePrior);
});
setMethodS3("getL2ShrinkageRecords", "TwoLassoOb", function(this,...) {
  return(this$.L2ShrinkageRecords);
});
setMethodS3("getPCoefs", "TwoLassoOb", function(this,...) {
   return(this$.pCoefs[1]);
});
setMethodS3("getFailureFlag", "TwoLassoOb", function(this,...) {
   return(this$.FailureFlag[1]);
});
setMethodS3("getNSample", "TwoLassoOb", function(this,...) {
   return(this$.nSample[1]);
});
setMethodS3("getTt1", "TwoLassoOb", function(this,...) {
   return(this$.tt1[1:length(this$.tt1)]);
});
setMethodS3("getTt2", "TwoLassoOb", function(this,...) {
   return(this$.tt2[1:length(this$.tt2)]);
});
setMethodS3("getLambdaAK", "TwoLassoOb", function(this,...) {
   return(this$.LambdaAK[1:length(this$.LambdaAK)]);
});
setMethodS3("getLambdaDK", "TwoLassoOb", function(this,...) {
   return(this$.LambdaDK);
});
setMethodS3("getOrderSeq", "TwoLassoOb", function(this,...) {
   return(this$.OrderSeq);
});
setMethodS3("getPiAPrior", "TwoLassoOb", function(this,...) {
   return(this$.PiAPrior);
});
setMethodS3("getPiAVector", "TwoLassoOb", function(this,...) {
   return(this$.PiAVector);
});
setMethodS3("getSigmaVector", "TwoLassoOb", function(this,...) {
   return(this$.SigmaVector);
});
setMethodS3("getSigmaPrior", "TwoLassoOb", function(this,...) {
   return(this$.SigmaPrior);
});
setMethodS3("getSdX", "TwoLassoOb", function(this,...) {
   return(this$.sdX);
});
setMethodS3("getSdY", "TwoLassoOb", function(this,...) {
   return(this$.sdY);
});
setMethodS3("getGroupSexp", "TwoLassoOb", function(this,...) {
   return(this$.GroupSexp);
});
setMethodS3("getOnGammas", "TwoLassoOb", function(this,...) {
   return(this$.OnGammas);
});
setMethodS3("getOnBeta", "TwoLassoOb", function(this,...) {
   return(this$.OnBeta[1:length(this$.OnBeta)]);
});
setMethodS3("getRecBBOn1", "TwoLassoOb", function(this,...) {
   if (is.null(this$.RecBBOn1)) { return(NULL); }
   return(matrix(as.numeric(this$.RecBBOn1),
     length(this$.RecBBOn1[,1]),  length(this$.RecBBOn1[1,])
   ) );
});
setMethodS3("getBBOn1", "TwoLassoOb", function(this,...) {
   return(this$.BBOn1[1:length(this$.BBOn1)]);
});
setMethodS3("getX", "TwoLassoOb", function(this,...) {
   if (is.null(this$.X)) { return(NULL); }
   return(this$.X);
});
setMethodS3("getPiAVectorInputs", "TwoLassoOb", function(this,...) {
   return(matrix(as.numeric(this$.PiAVectorInputs),
     length(this$.PiAVectorInputs)[,1],
     length(this$.PiAVectorInputs[1,])  ));
});
setMethodS3("getSigmaVectorInputs", "TwoLassoOb", function(this,...) {
   return(matrix(as.numeric(this$.SigmaVectorInputs),
     length(this$.SigmaVectorInputs)[,1],
     length(this$.SigmaVectorInputs[1,])  ));
});

setMethodS3("getY", "TwoLassoOb", function(this,...) {
   if (is.null(this$.Y)) { return(NULL); }
   return(this$.Y);
});
setMethodS3("getXtX", "TwoLassoOb", function(this,...) {
   if (is.null(this$.XtX)) { return(NULL); }
   return(this$.XtX);
});
setMethodS3("getXtY", "TwoLassoOb", function(this,...) {
   if (is.null(this$.XtY)) { return(NULL); }
   return(this$.XtY[1:length(this$.XtY)]);
});
setMethodS3("getUSOnBeta", "TwoLassoOb", function(this,...) {
   return(this$.OnBeta[1:length(this$.OnBeta)]);
});
setMethodS3("getReturnBetas", "TwoLassoOb", function(this,...) {
   if (is.null(this$.OnBeta)) { return(NULL); }
   if (this$StandardFlag == TRUE && !is.null(this$.sdX) &&
     !is.null(this$.sdY) ) {
     return( this$.sdY * this$.OnBeta / this$.sdX );  
   }
   return(as.vector(as.numeric(this$.OnBeta)) );
});
setMethodS3("getUSRecOnBeta", "TwoLassoOb", function(this,...) {
   if (is.null(this$.RecOnBeta)) { return(NULL); }
   return(matrix(as.numeric(this$.RecOnBeta),
     length(this$.RecOnBeta[,1]),  length(this$.RecOnBeta[1,]) )
    );
});
setMethodS3("getRecOnBeta", "TwoLassoOb", function(this,...) {
   if (is.null(this$.RecOnBeta)) { return(NULL); }
   if (this$StandardFlag == TRUE && !is.null(this$.sdX) &&
     !is.null(this$.sdY)) {
     return(matrix(as.numeric( this$.sdY * 1/this$.sdX * this$.RecOnBeta),
       length(this$.RecOnBeta[,1]),  length(this$.RecOnBeta[1,]) )
      );     
   }
   return(matrix(as.numeric(this$.RecOnBeta),
     length(this$.RecOnBeta[,1]),  length(this$.RecOnBeta[1,]) )
    );
});
##################################################################
## TwoLasso: Basic Coding of the EM2Lasso Algorithm
##   EM2Lasso is an algorithm solving the Two-Lasso problem for 
##   a number of specified lambda values
##    (specified by multiplying initial values by a multiplicative constant)
##
##   Does not estimate pi_A or sigma
##
##   By Default uses CoordinateDescent Code from CoordinateDescent.cc
##   However, It can use LARS algorithm in stead from LarsC.cc
Two2LassoCpp <- function(X = -1, Y = -1, nSample = -1, 
  PiA = .5, SigmaSq = -999, 
  LambdaDK = -999, LambdaAK = -999, OrderSeq = -999,
  StLambdaD = -999, StLambdaA = -999, lambdaDMultC = 2^(2/8), 
  lambdaAMultC = 2^(-2/8), lambdaAmultStop = 20, TotalRuns = 20, 
  MaxCauchy = 20, CauchyEpsilon = .001, StartBeta = -999, StandardFlag = 0,
  XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
  InitKKs = -5, WLSWeights = -1, TDFNu = -1,
  Verbose = -1, SigmaVec = -1,
  RecordFlag = 0,
  PiAPrior = c(-1.5,-1.5), m1 = -1, m2 = -1, SigmaPrior = c(-1,-1),
  SigmaBar = -1, SigmaDf = -1, HoldOn = FALSE, Groupers = NULL,
  RunFlag = 1, L2ShrinkagePrior = c(-999,-999), 
  RecordL2Shrinkage = FALSE, DoCrossValidate = FALSE,...) {
  if (Verbose > 0) {
    print("Two2LassoCpp:: Start"); flush.console();
  } 
  try(library(methods, warn.conflicts = FALSE, quietly=TRUE));
  try(library(R.methodsS3, warn.conflicts=FALSE, quietly=TRUE));  
  try(library(R.oo, warn.conflicts=FALSE, quietly=TRUE)); 
  
  GroupSexp = NULL; TCSPointer = NULL;
  if(!is.loaded("LarsConvergencyShell") ) {
    LoadTwoCC();
  } 
  if (is.null(LambdaAK) || LambdaAK[1] == -999) {
    if (StLambdaA < 0 & PiA > 0 && PiA < 1.0) {
       StLambdaA =  .6 * log((1-PiA)/PiA);      
    } else if (StLambdaA < 0 && FixKa > 0) {
      pDtry <- FixKa / pCoefS; 
      StLambdaA <-  .6 * log((1-pDtry)/pDtry);
    } else if (StLambdaA < 0) {
      print("Why did you not give us StLambdaA");
    }
  } 
  if (LambdaDK[1] == -999) {
    if (StLambdaD < 0 & PiA > 0 && PiA < 1.0) {
      StLambdaD =  .6 * log((1-PiA)/PiA);      
    } else if (StLambdaD < 0 && FixKa > 0) {
      pDtry <- FixKa / pCoefs; 
      StLambdaD =  .6 * log((1-pDtry)/pDtry);
    } else if (StLambdaD < 0) {
      print("Why did you not give us StlambdaD");
    }
  }    

  if (PiA   >= 1) {PiA = .99;}
  if (PiA <=0) {PiA = .01;}
  if (is.na(PiA) || is.nan(PiA)) { PiA = .5;}

  if (is.matrix(X) == FALSE || length(X) == 1 || X == -1) {
    if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == -1) {
      print("TwoLasso, Error no proper X given, XTX False"); return;  
    }
    if ( (is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) ||
      length(XTY) == 1 || XTY == -1) {
      print("TwoLasso, Error no proper X given, XTY False"); return;
    }
    if (length(nSample) > 1 || nSample == -1) {
      print("TwoLasso, Error, no proper NLen given");
      return(-1);
    }
    XTXFlag = 1;
    XtY = as.double(XtY);   p = length(XtY);
    TnSample = -nSample; pCoefs = length(XTX[1,]);
    WLSWeights = -1;  TDFNu = -1;
    sdX = sqrt(diag(XtX)); sdY = NULL;
    if (Verbose > 0) {
      print(paste("Two2Lasso: We're going to supply XtX, XtY with TnSample = ", TnSample, sep=""));
      flush.console();
    }
  } else if (( is.vector(Y) == FALSE && is.matrix(Y) == FALSE) || 
    length(Y) == 1 || Y == -1) {
    print("TwoLasso, Error no proper Y given"); 
    print(" Y is : "); print( Y ); return(-1);
  } else {
    XTXFlag = 0;
    if (length(Y) > 1000 && length(X[1,]) < length(Y) &&
      (length(WLSWeights) != length(Y) ||
        length(TDFNu) == 1 && TDFNu > 0  ) ) {
      XTXFlag = 1;
      XtX = as.double( t(X)%*% X );
      XtY = as.double( t(X) %*% Y );  p = length(XtY);
      TnSample = -length(Y); pCoefs = length(X[1,]); 
      sdX = apply(X,2,sd); sdY = sd(as.vector(Y));    
      if (Verbose > 0) {
        print(paste("Two2Lasso: Y too long, We're going to supply XtX, XtY with TnSample = ", TnSample, sep=""));
        flush.console();
      }
    }  else {
      TnSample = length(Y); pCoefs = length(X[1,]); p = pCoefs;
      if (length(WLSWeights) != length(Y)) {
        WLSWeights = -1;
      }
      if (pCoefs > 10000) {
        InitKKs = round(pCoefs *PiA);
        if (InitKKs > pCoefs) {InitKKs = pCoefs;}
      }
      if (InitKKs < 0) {
        InitKKs = min(50, length(X[1,]));
      }
      if (length(TDFNu) == 1 && TDFNu > 0) {
        WLSWeights = rep(1, length(Y))
      } else {
        TDFNu = -1;
      }
      sdX = apply(X,2,sd); sdY = sd(as.vector(Y));
      if (Verbose > 0) {
        print(paste("Two2Lasso: We're going to supply X, Y with TnSample = ", TnSample, sep=""));
        flush.console();
      }   
    }      
  }   
  nSample = abs(TnSample);
   if (InitKKs < 0 && pCoefs > 500) {
     InitKKs = min(round(pCoefs * 1.5 * PiA), pCoefs-1);
   } 
          
   if (SigmaSq < 0) {
     SigmaSq = (.4 * sd (Y) / mean(apply(X,2,sd)))^2;
   }

   if (InverseGammaConstant < 0) {
     print(paste("InverseGammaConstant is less than zero is : ",
       InverseGammaConstant, sep=""));
     return(-1);
   }
   if (length(StartBeta) < pCoefs) {
     StartBeta = (1:pCoefs) * 0;
   }  

   sdYYs <- sd(as.vector(Y));
   sdXXs <- apply(X,2,sd);  
   StandardFlag = as.integer(StandardFlag) 
   if (StandardFlag == 1) {
     Y = StandardizeXX(Y);
     X = StandardizeXX(X);
     SigmaSq = SigmaSq  / sdYYs^2;
     StLambdaD = StLambdaD * sdYYs;
     StLambdaA = StLambdaA * sdYYs;
   }
   if (InitKKs > pCoefs) {
     InitKKs = pCoefs -1;
   } 

   OnBeta = as.vector(matrix(as.numeric(StartBeta), 1, length(StartBeta)));
   BBOn1 = OnBeta * 0;
   if (is.vector(LambdaDK) == FALSE || length(LambdaDK) == 1 || LambdaDK == -999) {
     LambdaDK = StLambdaDK * LambdaDMultiplier^(0:(TotalRuns-1));
     LambdaAK = StLambdaDK * LambdaAMultiplier^(0:(TotalRuns-1));
     LambdaAStop = round(LambdaAStop);
     LambdaAK[(1:TotalRuns) > LambdaAStop] = LambdaAK[LambdaAStop]
   } else {
     TotalRuns = min(length(LambdaDK), length(LambdaAK));
     LambdaDK= LambdaDK[1:TotalRuns];  LambdaAK = LambdaAK[1:TotalRuns];
   }
   if (is.vector(OrderSeq) == TRUE && length(OrderSeq) >= TotalRuns) {
     OrderSeq = OrderSeq[1:TotalRuns];
   } else if (length(OrderSeq) == 1 && OrderSeq == -999) {
     OrderSeq <- rep(8, TotalRuns);
   }
   
   OnGammas = rep(0,pCoefs) + PiA * LambdaAK[1] + (1-PiA) * LambdaDK[1];
   if (InitKKs > pCoefs) {
     InitKKs = pCoefs -1;
   } 
   if (length(SigmaVec) < TotalRuns) {
     SigmaVec = rep(SigmaSq, TotalRuns);
   }
   if (is.null(Verbose)) { Verbose = -1; }
   if (length(Verbose) != 1) { Verbose = -1;}
   if (Verbose > 3) {
     print(paste("LambdaAK = ", paste(LambdaAK, collapse=", "), " and SigmaVec = ", 
     paste(SigmaVec, collapse=", "), sep=""));  flush.console();
   }
   if (RecordFlag > 0) {
     RecBBOn1 = matrix(0, length(BBOn1), TotalRuns);
     RecOnBeta = matrix(0, length(OnBeta), TotalRuns);
   } else {
     RecBBOn1 = -1;
     RecOnBeta = -1;
   }
   PiAVector = rep(PiA, TotalRuns);  SigmaVector = rep(SigmaSq, TotalRuns);
   if (m1 > 0 && m2 > 0) {PiAPrior = c(m1,m2);}
   if (SigmaBar > 0 && SigmaDf > 0) {SigmaPrior = c(SigmaBar, SigmaDf); }
   
      if (is.null(Verbose)) { Verbose = -1; }
   if (length(Verbose) != 1) { Verbose = -1;}
   if (Verbose > 3) {
     print(paste("Dimension of X will be: (", 
       paste(dim(X), collapse=", "), ")", sep="")); flush.console();
   }
   GroupSEXP = NULL;
   tt1 = 0; tt2 = rep(0, length(LambdaAK));
   
   ##if(!is.null(L2ShrinkagePrior) && 
   ##  length(L2ShrinkagePrior) == 2 && L2ShrinkagePrior[1] >= -2) {  
   ##  if (RecordL2Shrinkage == TRUE) {
   ##    L2ShrinkageRecords = rep(0, length(LambdaAK));
   ##  } else { L2ShrinkageRecords = NULL;}
   ##} else {L2ShrinkageRecords = NULL;  L2ShrinkagePrior = NULL;}
   
   MakeCpp = FALSE;
   if (DoCrossValidate == TRUE) {
     MakeCpp = TRUE;
   }
   RT = TwoLassoR5$new(
     X=X, Y=Y, XtX = XtX, XtY = XtY,
     LambdaAK = LambdaAK, LambdaDK=LambdaDK, OrderSeq=OrderSeq,
     PiA=PiA, SigmaSq=SigmaSq,
     PiAVector=PiAVector, SigmaVector=SigmaVector, PiAPrior=PiAPrior,
     SigmaPrior=SigmaPrior, InitKKs=InitKKs, 
     OnGammas = OnGammas, StartBeta = StartBeta, 
     InverseGammaConstant=InverseGammaConstant,
     CauchyEpsilon=CauchyEpsilon, MaxCauchy=MaxCauchy,
     TDFNu=TDFNu, WLSWeights=WLSWeights, Verbose=Verbose,
     GroupSexp=GroupSexp, RecordFlag = RecordFlag,
     FailureFlag=FailureFlag, StandardFlag = StandardFlag, 
     L2ShrinkagePrior = L2ShrinkagePrior,
     RunFlag = RunFlag, MakeCpp = MakeCpp);
   if (!is.null(RT$XtX)) { 
     TnSample = -RT$nSample;
   } 
     
   if (HoldOn == FALSE) {
     if (Verbose > 0) {
       print("HoldOn is False: We'll do TwoLassoRegression Shell");
       flush.console();
     }
     TCpp =.Call("TwoLassoRegressionShell", nSample = TnSample, pCoef = RT$pCoef,
       RT$Y, RT$X, RT$XtY, RT$XtX, RT$tt1, RT$tt2, 
       RT$LambdaAK, RT$LambdaDK, RT$OrderSeq,
       RT$BBOn1, RT$OnBeta, RT$RecBBOn1, RT$RecOnBeta,
       RT$OnGammas, RT$PiAVector, RT$SigmaVector,
       RT$PiAPrior, RT$SigmaPrior,                                     
       RT$InitKKs, RT$InverseGammaConstant,
       RT$CauchyEpsilon, RT$MaxCauchy,
       RT$TDFNu, RT$WLSWeights, RT$L2ShrinkagePrior,
       RT$L2ShrinkageRecords, RT$Verbose,
       RT$GroupSexp);  TLS = NULL;  TCSPointer = NULL;   TCSPointer = TCpp;
    
   } else {
     if (Verbose > 0) {
       print("HoldOn is False: We'll do TwoLassoCreateTwoLassoSexpShell");
       flush.console();
     }
    
     TCpp = .Call("TwoLassoCreateTwoLassoSexpShell", nSample = TnSample, 
       pCoef = RT$pCoef, RT$Y, RT$X, RT$XtY, RT$XtX, 
       RT$tt1, RT$tt2, 
       RT$LambdaAK, RT$LambdaDK, RT$OrderSeq, 
       RT$BBOn1, RT$OnBeta, 
       RT$RecBBOn1, RT$RecOnBeta,
       RT$OnGammas, RT$PiAVector, RT$SigmaVector,
       RT$PiAPrior, RT$SigmaPrior,
       RT$InitKKs, RT$InverseGammaConstant,
       RT$CauchyEpsilon, RT$MaxCauchy,
       RT$TDFNu, RT$WLSWeights, L2ShrinkagePrior = RT$L2ShrinkagePrior,
       L2ShrinkageRecords = RT$L2ShrinkageRecords, Verbose = RT$Verbose,
       RT$GroupSexp, RunFlag = RT$RunFlag);
     TCSPointer = TCpp;
     RT$.TCSPointer = TCSPointer;
   }
   TCSPointer = TCpp;
   if (Verbose > 0) {
     print("TwoLassoCpp: Finished .Call Code"); flush.console();
   }
   RT$ReturnBetas = as.vector(matrix(as.numeric(RT$OnBeta), length(RT$OnBeta),1));
   RT$USBeta = as.vector(matrix(as.numeric(RT$OnBeta), length(RT$OnBeta),1)); 
     
   if (RT$StandardFlag > 0) {
     RT$ReturnBetas = sdYYs * RT$USBeta / sdXXs
   } 
   if (Verbose > 0) {
      print("TwoLassoCpp: Create RT"); flush.console();
   }

   FailureFlag = 0;
    if (nSample == -999) {
      FailureFlag = 1;
    }

   RT$TCpp = TCpp;  RT$FailureFlag = FailureFlag;

   if (Verbose > 0) {
      print("TwoLassoCpp: SetUSReturnBetas"); flush.console();
   }
   if (Verbose > 0) {
     print("Two2LassoCpp: Returning Home"); flush.console();
   }
   return(RT);              
} 

##################################################################
## TwoLasso: Basic Coding of the EM2Lasso Algorithm
##   EM2Lasso is an algorithm solving the Two-Lasso problem for 
##   a number of specified lambda values
##    (specified by multiplying initial values by a multiplicative constant)
##
##   Does not estimate pi_A or sigma
##
##   By Default uses CoordinateDescent Code from CoordinateDescent.cc
##   However, It can use LARS algorithm in stead from LarsC.cc
.Two2LassoCppMiniValidate <- function(X = -1, Y = -1, nSample = -1, 
  PiAVectorInputs = -1, SigmaVectorInputs = -1,
  LambdaDK = -999, LambdaAK = -999, LambdaDKInputs = -999,
  LambdaAKInputs = -999, OrderSeq = -999,
  StLambdaD = -999, StLambdaA = -999, lambsdaDMultC = 2^(2/8), 
  lambdaAMultC = 2^(-2/8), lambdaAmultStop = 20, TotalRuns = 20, 
  MaxCauchy = 20, CauchyEpsilon = .001, StartBeta = -999, StandardFlag = 0,
  XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
  InitKKs = -5, WLSWeights = -1, TDFNu = -1,
  Verbose = -1, 
  RecordFlag = 0, L2ShrinkagePrior = NULL,
  mStrength = -1, SigmaPriorStrength = -1, Groupers = NULL,
  StartBetaMatrix = NULL, Two2LassoCppMiniValidate =FALSE, 
  RecordL2Shrinkage = FALSE,
  RunFlag = TRUE,...) {
  eval(parse(text=MyExistsF("nSample", -1)));
  if (!exists("LambdaAK") || is.null(LambdaAK)) { LambdaAK = -999; }
  if (!exists("LambdaDK") || is.null(LambdaDK)) { LambdaAK = -999; }
  eval(parse(text=MyExistsF("TotalRuns", 20)));
  eval(parse(text=MyExistsF("MaxCauchy", 20)));
  eval(parse(text=MyExistsF("lambdaAmultStop", 20)));
  eval(parse(text=MyExistsF("TotalRuns", 20)));
  eval(parse(text=MyExistsF("TDFNu", -1)));
  eval(parse(text=MyExistsF("WLSWeights", -1)));
  eval(parse(text=MyExistsF("FixKa", -100)));
  eval(parse(text=MyExistsF("InitKKs", -5)));
  eval(parse(text=MyExistsF("RecordFlag", 0)));
  eval(parse(text=MyExistsF("mStrength", -1)));
  eval(parse(text=MyExistsF("SigmaPriorStrength", -1)));
  eval(parse(text=MyExistsF("Groupers", "NULL")));
  eval(parse(text=MyExistsF("StandardFlag", 0)));
  eval(parse(text=MyExistsF("RunFlag", "TRUE")));
  eval(parse(text=MyExistsF("XtX", -1)));
  eval(parse(text=MyExistsF("XtY", -5)));
  eval(parse(text=MyExistsF("InverseGammaConstant", 1.0)));
  
  try(library(methods, warn.conflicts = FALSE, quietly=TRUE));
  try(library(R.methodsS3, warn.conflicts=FALSE, quietly=TRUE));      
  try(library(R.oo)); 

  if (exists("SampleInputs") && !exists("PiAVectorInputs")
     && !exists("SigmaVectorInputs") ) {
    PiAVectorInputs = SampleInputs[,1];
    SigmaVectorInputs = SampleInputs[,2];   
  }  
  if (length(PiAVectorInputs) < 2 || 
    length(PiAVectorInputs) != length(SigmaVectorInputs)) {
    print("Two2LassoMiniValidate: Give me a string of Pi, Sigma Inputs!");
    flush.console();
    return(-1);  
  }
  if (Verbose >= -1) {
     print("TwoLassoMiniValidate: Starting. "); flush.console();
  }
  TCSPointer = NULL;

  if (mStrength > 0) {
    PiAPrior = c(mStrength, mStrength)
  } else {
    PiAPrior = c(-1,-1);
  }
  if (SigmaPriorStrength>0) {
    SigmaPrior = c(SigmaVectorInputs[1], SigmaPriorStrength);
  } else {
    SigmaPrior = c(-1,-1);
  }
  GroupSexp = NULL; TCSPointer = NULL;
  if(!is.loaded("LarsConvergencyShell") ) {
    LoadTwoCC();
  } 
  if (is.null(LambdaAK) || LambdaAK[1] == -999) {
    if (is.null(LambdaAKInputs) || length(LambdaAKInputs) <= 1) {
       print("Two2LassoMiniValidate: Give Me a LambdaAK or LambdaAK Inputs!  Return-1 error"); return(-1);
    } else {
      LambdaAK <- LambdaAKInputs[,1];
    }
  } 
  if (is.null(LambdaDK) || LambdaDK[1] == -999) {
    if (is.null(LambdaDKInputs) || length(LambdaDKInputs) <= 1) {
       print("Two2LassoMiniValidate: Give Me a LambdaDK or LambdaDK Inputs! Return-1 error"); return(-1);
    } else {
      LambdaDK <- LambdaDKInputs[,1];
    }
  }    


  if (is.matrix(X) == FALSE || length(X) == 1 || X == -1) {
    if (is.matrix(XTX) == FALSE || length(XTX) == 1 || XTX == -1) {
      print("Two2LassoMiniValidate, Error no proper X given, XTX False"); return;  
    }
    if ( (is.matrix(XTY) == FALSE && is.vector(XTY) == FALSE) ||
      length(XTY) == 1 || XTY == -1) {
      print("Two2LassoMiniValidate, Error no proper X given, XTY False"); return;
    }
    if (length(nSample) > 1 || nSample == -1) {
      print("Two2LassoMiniValidate, Error, no proper NLen given");
      return(-1);
    }
    XTXFlag = 1;
    XtY = as.double(XtY);
    TnSample = -nSample; pCoefs = length(XTX[1,]);
    WLSWeights = -1;  TDFNu = -1;
    sdX = sqrt(diag(XtX)); sdY = NULL;  InitKKs = -1;
  } else if (( is.vector(Y) == FALSE && is.matrix(Y) == FALSE) || 
    length(Y) == 1 || Y == -1) {
    print("TwoLasso, Error no proper Y given"); 
    print(" Y is : "); print( Y ); return(-1);
  } else {
    XTXFlag = 0;
    if (length(Y) > 1000 && length(X[1,]) < length(Y) &&
      (length(WLSWeights) != length(Y) ||
        length(TDFNu) == 1 && TDFNu > 0  ) ) {
      XTXFlag = 1;
      XtX = as.double( t(X)%*% X );
      XtY = as.double( t(X) %*% Y );
      TnSample = -length(Y); pCoefs = length(X[1,]); 
      sdX = apply(X,2,sd); sdY = sd(as.vector(Y));    
    }  else {
      TnSample = length(Y); pCoefs = length(X[1,]);
      if (length(WLSWeights) != length(Y)) {
        WLSWeights = -1;
      }
      if (pCoefs > 10000) {
        InitKKs = round(pCoefs *PiAVectorInputs[1]);
        if (InitKKs > pCoefs) {InitKKs = pCoefs;}
      }
      if (length(TDFNu) == 1 && TDFNu > 0) {
        WLSWeights = rep(1, length(Y))
      } else {
        TDFNu = -1;
      }
      if (InitKKs < 0) {
        InitKKs = min(50, length(X[1,]));
      }
      sdX = apply(X,2,sd); sdY = sd(as.vector(Y));
    }      
  }   
  nSample = abs(TnSample);
   if (InitKKs < 0 && pCoefs > 500) {
     InitKKs = min(round(pCoefs * 1.5 * PiAVectorInputs[1]), pCoefs-1);
   } 

   if (InverseGammaConstant < 0) {
     print(paste("InverseGammaConstant is less than zero is : ",
       InverseGammaConstant, sep=""));
     return(-1);
   }
   if (!exists("StartBeta") || is.null(StartBeta) || length(StartBeta) < pCoefs) {
     StartBeta = (1:pCoefs) * 0;
   }  

   sdYYs <- sd(as.vector(Y));
   sdXXs <- apply(X,2,sd);   
   if ((is.logical(StandardFlag) && StandardFlag == TRUE) ||
     (is.numeric(StandardFlag) && StandardFlag >= 1)) {
     Y = StandardizeXX(Y);
     X = StandardizeXX(X);
     SigmaVectorInputs = SigmaVectorInputs  / sdYYs^2;
     StLambdaD = StLambdaD * sdYYs;
     LambdaDK = LambdaDK * sdYYs;
     StLambdaA = StLambdaA * sdYYs;
     LambdaAK = LambdaAK * sdYYs;
   }
   if (InitKKs > pCoefs) {
     InitKKs = pCoefs -1;
   } 

   OnBeta = as.vector(matrix(as.numeric(StartBeta), 1, length(StartBeta)));
   BBOn1 = OnBeta * 0;
   if (is.vector(LambdaDK) == FALSE || length(LambdaDK) == 1 || LambdaDK == -999) {
     LambdaDK = StLambdaDK * LambdaDMultiplier^(0:(TotalRuns-1));
     LambdaAK = StLambdaDK * LambdaAMultiplier^(0:(TotalRuns-1));
     LambdaAStop = round(LambdaAStop);
     LambdaAK[(1:TotalRuns) > LambdaAStop] = LambdaAK[LambdaAStop]
   } else {
     TotalRuns = min(length(LambdaDK), length(LambdaAK));
     LambdaDK= LambdaDK[1:TotalRuns];  LambdaAK = LambdaAK[1:TotalRuns];
   }
   if (is.vector(OrderSeq) == TRUE && length(OrderSeq) >= TotalRuns) {
     OrderSeq = OrderSeq[1:TotalRuns];
   } else if (length(OrderSeq) == 1 && OrderSeq == -999) {
     OrderSeq <- rep(8, TotalRuns);
   }
   PiAVector = rep(PiAVectorInputs[1], length(OrderSeq));
   SigmaVector = rep(SigmaVectorInputs[1], length(OrderSeq));
   SigmaSq = SigmaVector[1];
   
   OnGammas = rep(0,p) + PiAVectorInputs[1] * LambdaAK[1] + 
     (1-PiAVectorInputs[1]) * LambdaDK[1];
   if (InitKKs > pCoefs) {
     InitKKs = pCoefs -1;
   } 

   if (is.null(Verbose)) { Verbose = -1; }
   if (length(Verbose) != 1) { Verbose = -1;}
   if (Verbose > 3) {
     print(paste("LambdaAK = ", paste(LambdaAK, collapse=", "), 
     " and SigmaVectorInputs = ", 
     paste(SigmaVectorInputs, collapse=", "), sep=""));  flush.console();
   }
  ## RecBBOn1 = matrix(0, length(BBOn1), length(PiAVectorInputs));
  ## RecOnBeta = matrix(0, length(OnBeta), length(PiAVectorInputs));
   
   if (Verbose > 3) {
     print(paste("Dimension of X will be: (", 
       paste(dim(X), collapse=", "), sep="")); flush.console();
   }
   GroupSEXP = NULL;
   tt1 = 0; tt2 = rep(0, length(LambdaAK));
   
   if (Verbose > 1) {
     print(paste("Dimension of XtX will be ", 
       paste(dim(XtX), collapse=", "))); flush.console();
     print(paste("Dimension of XtY will be ", 
       paste(dim(XtX), collapse=", "))); flush.console();  
     print(paste("Thus InitKKs = ", InitKKs, sep="")); flush.console();   
   }

   RT = TwoLassoR5$new(
     X=X, Y=Y, XtX = XtX, XtY = XtY,
     LambdaAK = LambdaAK, LambdaDK=LambdaDK, OrderSeq=OrderSeq,
     PiA = PiA, SigmaSq = SigmaSq,
     PiAVector=PiAVector, SigmaVector=SigmaVector, PiAPrior=PiAPrior,
     SigmaPrior=SigmaPrior, InitKKs=InitKKs, 
     OnGammas = OnGammas, StartBeta = StartBeta, 
     InverseGammaConstant=InverseGammaConstant,
     CauchyEpsilon=CauchyEpsilon, MaxCauchy=MaxCauchy,
     TDFNu=TDFNu, WLSWeights=WLSWeights, Verbose=Verbose,
     GroupSexp=GroupSexp, RecordFlag =RecordFlag, 
     FailureFlag=FailureFlag, StandardFlag = as.integer(StandardFlag), 
     L2ShrinkagePrior = L2ShrinkagePrior,
     LambdaDKInputs = LambdaDKInputs, LambdaAKInputs = LambdaAKInputs,
     PiAVectorInputs=PiAVectorInputs, SigmaVectorInputs=SigmaVectorInputs, 
     StartBetaMatrix=StartBetaMatrix, FixKa = FixKa, 
     RunFlag = RunFlag, DoCrossValidate=TRUE);
  RT$RecBBOn1 = matrix(0, length(RT$OnBeta), length(PiAVectorInputs));
  RT$RecOnBeta = matrix(0, length(RT$OnBeta), length(PiAVectorInputs));
  if (Verbose > 0) {
    print(paste("Before move in Dim of RCVF = (", 
      paste(RT$RecordBetaCVFinish, collapse=", "), ")", sep=""));
    flush.console();
  }   
  if (Verbose > 0) {
    print("HoldOn is False: We'll do TwoLassoMiniValidate Shell (TwoLassoSexp.R): Running TwoLassoCrossValidateShell:");
    flush.console();     ;
  }
  TCSPointer = NULL;
  MyTryText <- "Fun = -1; 
  if (length(RT$PiA) <= 0) { RT$PiA <- .5; }
  TCSPointer =.Call(\"TwoLassoCrossValidateShell\", RT$Verbose,
    nSample = TnSample, pCoefs = RT$pCoef,
    RT$Y, RT$X, RT$.XtY, RT$XtX, RT$tt1, RT$tt2, 
    RT$LambdaAK, RT$LambdaDK, 
    RT$PiAVectorInputs, RT$SigmaVectorInputs,
    RT$LambdaAKInputs, RT$LambdaDKInputs,
    RT$StartBetaMatrix, RT$RecordBetaCVFinish,
    RT$OrderSeq, RT$BBOn1, RT$OnBeta, RT$RecBBOn1, RT$RecOnBeta,
    RT$OnGammas, RT$PiA, RT$SigmaSq, RT$PiAPrior, RT$SigmaPrior,                                     
    RT$InitKKs, RT$InverseGammaConstant,
    RT$CauchyEpsilon, RT$MaxCauchy,
    RT$TDFNu, RT$WLSWeights, RT$L2ShrinkagePrior, RT$L2ShrinkageRecords,
    RT$GroupSexp);  
  Fun = 1; "
  try(eval(parse(text=MyTryText)));
  try(eval(parse(text=SetGText("TCSPointer", envir="globalenv()", S=1))));
  try(RT$.TCSPointer <- TCSPointer);
  if (Verbose >= 2) {
    print("Just finished TwoLassoCrossValidate (TwoLassoSexp.R)!");
  }
  if (Fun == -1) {
    print(" **************************************** ");
    print(" *** TwoLasso ERROR ERROR ERROR ERROR");
    print(" *** After Error, Setting RT to Global and trigger error. "); 
    flush.console();
    eval(parse(text=SetGText("RT", "globalenv()", S=1)));
    eval(parse(text=SetGText("TnSample", "globalenv()", S=1)));
    eval(parse(text=SetGText("nSample", "globalenv()", S=1)));
    BADFAIL = NULL;  rm(BADFAIL);
    WHYFAIL = BADFAIL;
  }
   
  ##RT$RecOnBeta <- RT$OnBeta;
  SuccessFlag = NULL;
  try(SuccessFlag <- .Call("WhatIsSuccessFlag", TCSPointer));
  if (is.null(SuccessFlag)) {
    print("TwoLasso: Woah, SuccessFlag out of TCSPointer was NULL!"); flush.console();
    return(TCSPointer);
  }
  if (as.integer(SuccessFlag)  <= 0) {
    print(paste("TwoLasso: Bad Error: SuccessFlag was ", SuccessFlag, " at the end! \n"));
    flush.console();
  }
   
  try(RT$.TCSPointer <- TCSPointer);  
 
   if (Verbose > 0) {
     print("TwoLasso: Finished .Call Code"); flush.console();
   }
   RT$USBeta = as.vector(matrix(as.numeric(RT$OnBeta), length(OnBeta),1)); 
   ##RT$RecOnBeta <- RT$USBeta;  
   if (StandardFlag > 0) {
     RT$ReturnBetas = sdYYs * RT$USBeta / sdXXs
     ##RT$RecOnBeta <- RT$ReturnBetas;
   } 
   if (Verbose > 0) {
      print("TwoLasso: Create RT"); flush.console();
   }
   ##RT$RetBeta = as.vector(matrix(as.numeric(OnBeta), length(OnBeta),1));
   FailureFlag = 0;
    if (nSample == -999) {
      FailureFlag = 1;
    }

   try(RT$.TCSPointer <- TCSPointer);  try(RT$.FailureFlag <- FailureFlag);

   if (Verbose > 0) {
      print("TwoLasso: SetUSReturnBetas"); flush.console();
   }
   if (Verbose > 0) {
     print("Two2Lasso: Returning Home"); flush.console();
   }
   return(RT);              
} 


TwoLassoR5$methods(
  TwoLassoRefit = function(RunRegress = FALSE, NewOnBeta =-1, NewSigma=-1, NewPiA= -1,
  NewLambdaAK = -1, NewLambdaDK =-1, NewOrderSeq = -1,
  NewPiAPrior = -1, NewSigmaPrior = -1,
  NewVerbose = -25,...
  ) {
  if (!is.null(NewVerbose)) {
    .self$Verbose = NewVerbose[1];
  }
  if (.self$Verbose > 0) {
    print("TwoLassoRefit: Starting"); flush.console();
  }
  if (is.null(.self$OnBeta)) {
     print("TwoLassoRefit: this does not have an OnBeta");
     return(-1);
  }
  if (is.null(.self$.TCSPointer)) {
    print("TwoLassoRefit: this does not have Pointer to Object");
  }
  if (NewOnBeta != -1) {
    if (.self$Verbose > 0) {
      print("TwoLassoRefit: Changing OnBeta!"); flush.console();
    }
    if (length(NewOnBeta) != length(.self$OnBeta)) {
      print("TwoLassoRefit: NewOnBeta is not long enough");
      flush.console();
    } else {
      .Call("AlterVector",.self$OnBeta, length(.self$OnBeta), NewOnBeta);
    }
  }
  Zero = 0;
  if (NewSigma[1] > 0) {
    if (.self$Verbose > 0) {
      print("TwoLassoRefit: Changing Sigma!"); flush.console();
    }
    if (length(NewSigma) == length(.self$SigmaVector)) {
      .Call("AlterVector",.self$SigmaVector, Zero, 
         NewSigma);
      ##.self$SigmaVector[1:length(.self$SigmaVector)] = 
      ##  NewSigma[1:length(.self$SigmaVector)];
    } else {
      .Call("AlterVector",.self$SigmaVector, Zero,
         rep(NewSigma[1], length(.self$SigmaVector)));
      ##.self$SigmaVector[1:length(.self$SigmaVector)] = 
      ##  NewSigma[1];    
    }
  }
  if (NewPiA[1] > 0) {
    if (.self$Verbose > 0) {
      print("TwoLassoRefit: Changing NewPiA!"); flush.console();
    }
    if (length(NewPiA[1]) == length(.self$PiAVector)) {
      .Call("AlterVector",.self$PiAVector, Zero, 
         NewPiA);
    } else {
      .Call("AlterVector", .self$PiAVector, Zero, 
         rep(NewPiA[1], length(.self$PiAVector) ) );  
    }
  }
  if (length(NewLambdaAK) == length(.self$LambdaAK)) {
    if (.self$Verbose > 0) {
      print("TwoLassoRefit: Changing LambdAK!"); flush.console();
    }
    .Call("AlterVector",.self$LambdaAK, Zero, 
       NewLambdaAK);
  } 
  if (length(NewLambdaDK) == length(.self$LambdaDK)) {
    if (.self$Verbose > 0) {
      print("TwoLassoRefit: Changing LambdDK!"); flush.console();
    }
    .Call("AlterVector",.self$LambdaDK, Zero, 
       NewLambdaDK);
  } 
  if (length(NewOrderSeq) == length(.self$OrderSeq)) {
    if (.self$Verbose > 0) {
      print("TwoLassoRefit: Changing OrderSeq!"); flush.console();
    }
    .Call("AlterVector",.self$OrderSeq, Zero, 
       NewOrderSeq);
  }

  if (.self$Verbose > 0) {
     print("TwoLassoReRegress: R To RunRegress?"); flush.console();
  }
  SuccessFlag = NULL;
  try(SuccessFlag <- .Call("WhatIsSuccessFlag", .self$.TCSPointer));
  if (is.null(SuccessFlag)) {
    print("TwoLasso: Woah, SuccessFlag out of TCSPointer was NULL!\n");
    return(TCSPointer);
  }
  if (as.integer(SuccessFlag)  <= 0) {
    print(paste("TwoLasso: Bad Error: SuccessFlag was ", SuccessFlag, " at the end! \n"));
    flush.console();
  }
  if (RunRegress == TRUE) {
    DoRunRegress = 3;
    if (.self$Verbose >= 1) {
      print("Running TwoLassoReRegressShell! "); flush.console();
    }
    if (is.null(NewVerbose)) { NewVerbose = 2; }
    if (is.null(NewSigmaPrior)) { NewSigmaPrior = c(2.0,2.0); }
    LoopAgain = .Call("TwoLassoReRegressShell", .self$.TCSPointer,
      DoRunRegress, 
      NewPiAPrior=NewPiAPrior, NewSigmaPrior=NewSigmaPrior, NewVerbose=NewVerbose);
  } else if (.self$Verbose > -2) {
    print("TwoLassoReFit: You must confirm: RunRegress=TRUE to run");
    flush.console();
  }
  if (.self$Verbose > 0) {
     print("TwoLassoReFit: R Finished"); flush.console();
  }
  return(.self);
});


###############################################################################
## M2Items
##
##  Derive key quantities useful for M2Lasso
M2Items <- function(X, sigma, OrderSeq = c(20,4,1), SABS = .05,
  Verbose = 0) {
    if (is.null(X) || length(X) == 1) { return(NULL); }
    p <-length(X[1,]);  pCoefs <- p;
  if (length(sigma)  == 1) {
    USENOISE <- sigma; pCoefs <- p;
    SDD <- sqrt(min(apply(X,2,sd))^2 * (length(X[,1])-1) );
    LambdaDSeq = c(sqrt(3.7 / USENOISE) * SDD, 
      sqrt(3.7 / USENOISE) * SDD *  max( (2*length(X[,1])), 200) );
                            
    SDA <- SDD / sqrt(length(X[,1])-1);
    LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISE )) * sqrt( 2 * SDA / USENOISE)
    LambdaASeq = c( 1.1 * LambdaAS1, .05 * LambdaAS1);
    if (LambdaASeq[1] >= LambdaDSeq[2]) {
      LambdaASeq[1] = LambdaDSeq[1] + LambdaASeq[1];
      LambdaDSeq[1] = LambdaASeq[1] - LambdaDSeq[1];
      LambdaASeq[1] = LambdaASeq[1] - LambdaDSeq[1];
    }
    LambdaASeq = c(LambdaASeq[1], LambdaASeq);
    LambdaDSeq = c(LambdaDSeq[1], LambdaDSeq);
      OrderSeq = c(20, 4,1);
    if (Verbose > 1) {
      print("M2STuff: Returning option 1"); flush.console();
    }
    return(list(LambdaASeq=LambdaASeq, LambdaDSeq=LambdaDSeq, OrderSeq= OrderSeq));
  } else {
    USENOISEVector <- sigma; pCoefs <- p;
    SDD <- sqrt(min(apply(X,2,sd))^2 * (length(X[,1])-1) );
    LambdaDKInputs = rbind(sqrt(3.7 / USENOISEVector) * SDD, 
      sqrt(3.7 / USENOISEVector) * SDD *  max( (2*length(X[,1])), 200) );
                            
    SDA <- SDD / sqrt(length(X[,1])-1);
    LambdaAS1 <-  exp( - SABS^2 * SDA^2 / (2 * USENOISEVector )) * sqrt( 2 * SDA / USENOISEVector)
    LambdaAKInputs = rbind( 1.1 * LambdaAS1, 1/pCoefs * LambdaAS1);

      OrderSeq = c(24, 1);
    if (Verbose > 1) {
      print("M2STuff: Returning option 2"); flush.console();
    }
    return(list(LambdaAKInputs=LambdaAKInputs, 
      LambdaDKInputs=LambdaDKInputs, 
      LambdaASeq = LambdaAKInputs[,1], LambdaDSeq = LambdaDKInputs[,1],
      OrderSeq= OrderSeq));  
  
  }
}

M2Lasso <- function(X = -1, Y = -1, nSample = -1, PiA = -.5, SigmaSq = -999, 
  OrderSeq = c(25,10,1),
  lambdaDMultC = 2^(2/8), 
  lambdaAMultC = 2^(-2/8), lambdaAmultStop = 20, TotalRuns = 20, 
  MaxCauchy = 80, CauchyEpsilon = .001, StartBeta = -999, StandardFlag = 0,
  XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
  InitKKs = -5, WLSWeights = -1, dfTNoise = -1,
  Verbose = -1, SigmaVec = -1,
  RecordFlag = 0, PiAPrior = c(-1.5,-1.5), m1 = -1, m2 = -1, SigmaPrior = c(-1,-1),
  SigmaBar = -1, SigmaDf = -1, HoldOn = FALSE, Groupers = NULL,
  RunFlag = 1, L2ShrinkagePrior = c(-999,-999), RecordL2Shrinkage = FALSE,
  SABS = .05) {
  if (is.null(X) || length(X) == 1 && X == -1) {return(NULL);}
  pCoefs <- length(X[1,]);

  if (SigmaSq < 0 && SigmaPrior[1] < 0) {
    SigmaSq = .25 * (sd(as.vector(Y)))^2;
    SigmaPrior = c(sqrt(length(Y)), SigmaSq);
  }
  if (PiA < 0 && PiAPrior[1] < 0) {
    PiA = .5;
    PiAPrior = c(sqrt(n), sqrt(n));
  }
  M2Stuff <- M2Items(X, sigma=SigmaSq, 
      OrderSeq = OrderSeq, SABS = SABS);
  if (Verbose > 0) {
    print(paste("M2Lasso, sequences derived are: "));
    print(paste("LambdaASeq: ", paste(M2Stuff$LambdaASeq, collapse=", "), ""));
    print(paste("LambdaDSeq: ", paste(M2Stuff$LambdaDSeq, collapse=", "), ""));
    print(paste("OrderSeq: ", paste(M2Stuff$OrderSeq, collapse=", "), "")); 
  }
  NewLambdaAK = M2Stuff$LambdaASeq;  NewLambdaDK = M2Stuff$LambdaDSeq;
    NewOrderSeq = M2Stuff$OrderSeq;   
  FinalM2Out  <<- Two2LassoCpp(X = X, Y = Y, 
    PiA = PiA, SigmaSq = SigmaSq, 
    LambdaDK = NewLambdaDK, LambdaAK = NewLambdaAK, OrderSeq = NewOrderSeq,
    MaxCauchy = MaxCauchy, CauchyEpsilon = CauchyEpsilon, 
    StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
    XtX = XtX, XtY = XtY, InverseGammaConstant = 1, FixKa = -100, 
    InitKKs = InitKKs, WLSWeights = WLSWeights, TDFNu = dfTNoise,
    Verbose = Verbose, RecordFlag = RecordFlag, HoldOn = TRUE, 
    Groupers = NULL,  L2ShrinkagePrior = L2ShrinkagePrior, 
    RecordL2Shrinkage = RecordL2Shrinkage, PiAPrior = PiAPrior, SigmaPrior = SigmaPrior,
    RunFlag = 1);
  return(FinalM2Out);
}

