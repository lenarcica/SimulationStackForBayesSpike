## CrossValidatorTwoLasso;
FastCFoldValidate2Lasso <- function(X,Y, DoLogit=FALSE, nPi=20, nSigma=10, cFolds =5,
  OrderSeq = c(20,4,1), SABS=.05, MaxCauchy = 100, CauchyEpsilon = .00001,
  Verbose = 0) {
  p = length(X[1,]);  n = length(Y);
  FastCFoldValidate2Lasso = Verbose;
  if (Verbose > 0) {
    print("FastCFoldValidate2Lasso: Starting"); flush.console();
  }
  if (all(Y %in% c(0,1))) {
    if (DoLogit==FALSE) {
      print("CrossValidatorTwoLasso.r::FastCFoldValidate2Lasso() why is all Y in (0,1) but DoLogit=FALSE?");
    }
  }

  SampleInputs = TwoLassoSampleInputs(X,Y,nPi, nSigma, DoLogit=DoLogit);
  M2Stuff <- M2Items(X, sigma=SampleInputs[1,2], 
    OrderSeq = OrderSeq, SABS = SABS)
  LambdaAK = M2Stuff$LambdaASeq;  LambdaDK = M2Stuff$LambdaDSeq;
  OrderSeq = M2Stuff$OrderSeq;
  try(eval(parse(text=GetG0Text(".HiddenTwoLassoOb", "globalenv()"))), silent=TRUE);
  .HiddenTwoLassoOb <- Two2LassoCpp(X = X, Y = Y, DoLogit=DoLogit,
    PiA = SampleInputs[1,1], SigmaSq = SampleInputs[1,2], 
    LambdaDK = LambdaDK, LambdaAK = LambdaAK, OrderSeq = OrderSeq,
    MaxCauchy = MaxCauchy * .5, CauchyEpsilon = CauchyEpsilon * sqrt(p), 
    StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
    XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
    InitKKs = -5, WLSWeights = -1, TDFNu = -1,
    Verbose = -1, RecordFlag = 0, HoldOn = TRUE, Groupers = NULL,
    RunFlag = 1)
   try(eval(parse(text=LockGText(".HiddenTwoLassoOb", "globalenv()"))), silent=TRUE);
   
  eval(parse(text=GetG0Text("TwoLassoStayFitFunction", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("TwoLassoStayFitFunction", "TWOLASSONAMESPACE", S=1))); 
  TwoLassoStayFitFunction <- function(X,Y, DoLogit=FALSE, SampleInput,StartBeta=-1, StandardFlag = FALSE,...) {
    if (FastCFoldValidate2Lasso > 2) {
      print(paste("FastCFoldValidate2Lasso: TwoLassoStayFitRun SI: (",
        paste(SampleInput, collapse=", "),")", sep="")); flush.console();
    }
    M2Stuff <- M2Items(.HiddenTwoLassoOb$X, sigma=SampleInput[2], 
      OrderSeq = OrderSeq, SABS = SABS, Verbose = Verbose);
    NewLambdaAK = M2Stuff$LambdaASeq;  NewLambdaDK = M2Stuff$LambdaDSeq;
    NewOrderSeq = M2Stuff$NewOrderSeq;
    HideAgain = .HiddenTwoLassoOb$TwoLassoRefit(RunRegress=TRUE, 
      NewSigma = SampleInput[2], NewPiA = SampleInput[1],
      NewLambdaAK = NewLambdaAK, NewLambdaDK = NewLambdaDK, 
      NewOrderSeq = NewOrderSeq, StandardFlag = StandardFlag);
    return(HideAgain$ReturnBetas);
  }
  eval(parse(text=LockGText("TwoLassoStayFitFunction", "globalenv()", S=1)));
  try(eval(parse(text=LockGText("TwoLassoStayFitFunction", "TWOLASSONAMESPACE", S=1))), silent=TRUE);  
  RunOnAChange <- function(XTrain, YTrain, DoLogit=FALSE) {
    if (FastCFoldValidate2Lasso > 1) {
      print("FastCFoldValidate2Lasso: RunOnAChange Executing"); flush.console();
    }
    M2Stuff <- M2Items(XTrain, sigma=SampleInputs[1,2], 
      OrderSeq = OrderSeq, SABS = SABS)
    LambdaAK = M2Stuff$LambdaASeq;  LambdaDK = M2Stuff$LambdaDSeq;
    OrderSeq = M2Stuff$OrderSeq;
    if (!exists("DoLogit")) {
      if (all(YTrain %in% c(0,1))) {
         DoLogit = TRUE;
      } else { DoLogit = FALSE; }
    }  else if (all(YTrain %in% c(0,1))) {
      if (DoLogit==FALSE) {
        print("CrossValidatorTwoLasso.r::RunOnAChange(): Why is DoLogit False, if all Y in 0,1?"); flush.console();
      }
    }
    try(eval(parse(text=GetG0Text(".HiddentTwoLassoOb", "gloablenv()"))), silent=TRUE);
    .HiddenTwoLassoOb <- Two2LassoCpp(X = XTrain, Y = YTrain, 
      DoLogit=DoLogit,
      PiA = SampleInputs[1,1], SigmaSq = SampleInputs[1,2], 
      LambdaDK = LambdaDK, LambdaAK = LambdaAK, OrderSeq = OrderSeq,
      MaxCauchy = MaxCauchy, CauchyEpsilon = CauchyEpsilon, 
      StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
      XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
      InitKKs = -5, WLSWeights = -1, TDFNu = -1,
      Verbose = -1, RecordFlag = 0, HoldOn = TRUE, Groupers = NULL,
      RunFlag = 1);
    try(eval(parse(text=LockGText(".HiddenTwoLassoOb", "globalenv()", S=1))), silent=TRUE);
  }
  
  MyCFold = CFoldValidate(X,Y, DoLogit=DoLogit, TwoLassoStayFitFunction,SampleInputs,
    cFolds = cFolds, RunOnAChange = RunOnAChange,CFoldVerbose=Verbose -2);
  try(eval(parse(text=SetGText("MyCFold", "globalenv()", S=1))))   
  M2Stuff <- M2Items(X, sigma=MyCFold$BestSampleInputs[2], 
      OrderSeq = OrderSeq, SABS = SABS);
  NewLambdaAK = M2Stuff$LambdaASeq;  NewLambdaDK = M2Stuff$LambdaDSeq;
    NewOrderSeq = M2Stuff$OrderSeq; 
  
   eval(parse(text=GetG0Text("FinalM2Out", "globalenv()")));    
  FinalM2Out  <- Two2LassoCpp(X = X, Y = Y, DoLogit=DoLogit, 
    PiA = MyCFold$BestSampleInputs[1], SigmaSq = MyCFold$BestSampleInputs[2], 
    LambdaDK = NewLambdaDK, LambdaAK = NewLambdaAK, OrderSeq = NewOrderSeq,
    MaxCauchy = MaxCauchy, CauchyEpsilon = CauchyEpsilon, 
    StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
    XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
    InitKKs = -5, WLSWeights = -1, TDFNu = -1,
    Verbose = -1, RecordFlag = 0, HoldOn = TRUE, Groupers = NULL,
    RunFlag = 1);
  eval(parse(text=SetGText("FinalM2Out", "globalenv()", S=1)));
  try(FinalM2Out$MyCFold <- MyCFold);
  return(FinalM2Out);
}

MyExistsF <- function(ATT, AN) {
  MyT <- paste("
    if (!exists(\"", ATT, "\")) { 
      ", ATT, " = ", AN, ";
    }", sep="");
  return(MyT);
}
## CrossValidatorTwoLasso;
#############################################################################
## FastestCFoldValidate2Lasso 
##
##   A function running cross validation
##
FastestCFoldValidate2Lasso <- function(X,Y, DoLogit=FALSE, nPi=20, nSigma=10, cFolds =5,
  OrderSeq = c(20,1), SABS=.05, MaxCauchy = 100, CauchyEpsilon = .00001,
  Verbose = 0, RFastestCFoldValidate2Lasso = NULL,
  L2ShrinkagePrior = NULL,  RecordL2Shrinkage = FALSE,
  RandomStartBetas = 2) {
  eval(parse(text=MyExistsF("nPi", 20)));
  eval(parse(text=MyExistsF("nSigma", 10)));
  eval(parse(text=MyExistsF("cFolds", 5)));
  eval(parse(text=MyExistsF("OrderSeq", "c(20,1)")));
  eval(parse(text=MyExistsF("SABS", .05)));
  eval(parse(text=MyExistsF("MaxCauchy", 100)));
  eval(parse(text=MyExistsF("CauchyEpsilon", .00001)));
  eval(parse(text=MyExistsF("L2ShrinkagePrior", "NULL")));
  eval(parse(text=MyExistsF("RandomStartBetas", 2)));
  eval(parse(text=MyExistsF("RecordL2Shrinkage", "FALSE")));
  if (all(Y %in% c(0,1))) {
     if (DoLogit==FALSE) {
        print("CrossValidatorTwoLasso.r::FastestCFoldValidate2Lasso() why is all Y in (0,1) but DoLogit is FALSE?"); flush.console();
     }
  }
  
  p = length(X[1,]);  n = length(Y);

  if (!exists("RFastestCFoldValidate2Lasso") || 
    is.null(RFastestCFoldValidate2Lasso) ||
    (is.numeric(RFastestCFoldValidate2Lasso) && 
      RFastestCFoldValidate2Lasso < Verbose) ||
    (is.logical(RFastestCFoldValidate2Lasso) && 
      RFastestCFoldValidate2Lasso == FALSE && Verbose > 0) ) {
    RFastestCFoldValidate2Lasso = Verbose; 
    eval(parse(text=SetGText("RFastestCFoldValidate2Lasso", "globalenv()", "S=1")));   
  }
  
  if (!is.null(RFastestCFoldValidate2Lasso) && 
    ((is.numeric(RFastestCFoldValidate2Lasso) && RFastestCFoldValidate2Lasso > 0) ||
      (is.logical(RFastestCFoldValidate2Lasso) && RFastestCFoldValidate2Lasso == TRUE) ) ) {
    print("FastestCFoldValidate2Lasso: Starting"); flush.console();
    print(paste(" RFastestCFoldValidate2Lasso = ", 
      RFastestCFoldValidate2Lasso, sep="")); flush.console();
  }
  
  if (all(Y %in% c(0.0,1.0))) {
    if (DoLogit == FALSE) {
      print("FastestCFoldValidate2Lasso: Don't you want to do a Logit based? "); flush.console();
    }
  }
  if (Verbose > 2) {
    print("FastestCFoldValidate2Lasso: Starting SampleInputs"); flush.console();
  }

  SampleInputs = TwoLassoSampleInputs(X,Y,nPi, nSigma, DoLogit=DoLogit);

  if (Verbose > 2) {
    print("FastestCFoldValidate2Lasso: Creating M2Stuff"); flush.console();
  }
  try(eval(parse(text=GetG0Text(".M2Stuff", "globalenv()",S=1))),silent=TRUE);
  .M2Stuff = M2Items(X, sigma=SampleInputs[,2], 
    OrderSeq = OrderSeq, SABS = SABS, Verbose = Verbose);
  eval(parse(text=SetGText(".M2Stuff", "globalenv()", S=1)));
  if (Verbose > 1) {
    print(".M2Stuff - gotten, pasteto global"); flush.console();
  }
  eval(parse(text=SetGText(".M2Stuff", "globalenv()")));
  LambdaAKInputs = .M2Stuff$LambdaAKInputs;
  LambdaDKInputs = .M2Stuff$LambdaDKInputs;   
  LambdaAK = .M2Stuff$LambdaASeq;
  LambdaDK = .M2Stuff$LambdaDSeq; 
  if (is.numeric(RandomStartBetas) && RandomStartBetas > 0) {
     MyXeigen = NULL;
     try(MyXeigen <- eigen(t(X) %*% X));
     if (is.null(MyXeigen)) {
       StartBetaMatrix = matrix(rnorm(p * ceiling(RandomStartBetas)),
         p, ceiling(RandomStartBetas));
     } else {
       Di = MyXeigen$values;
       Di[Di != 0] = 1/ Di[Di != 0];
       try(Di[Di < 0] <- 0);
       try(Di[is.na(Di)] <- 0);
       PointMean <- t(MyXeigen$vectors) %*% diag(Di) %*% 
         MyXeigen$vectors  %*% t(X) %*% Y;
       SigmaEst =  sum( (Y - X %*% PointMean)^2) / (length(Y)-1);
       StartBetaMatrix = matrix(rnorm(p * ceiling(RandomStartBetas)),
         p, ceiling(RandomStartBetas));
       StartBetaMatrix = matrix(rep(PointMean, RandomStartBetas),
         length(PointMean), RandomStartBetas) +
         sqrt(SigmaEst) * (t(MyXeigen$vectors) %*% 
           sqrt(diag(Di))) %*% MyXeigen$vectors %*% StartBetaMatrix;
       eval(parse(text=SetGText("StartBetaMatrix", "globalenv()", S=1)));
     } 
  } else {
    StartBetaMatrix = NULL;
    eval(parse(text=SetGText("StartBetaMatrix", "globalenv()", S=1)));
  }
  if (Verbose > 0) {
    print("Got StartBetaMatrix"); flush.console();
  }
  if (Verbose > 0) {
    print("FastestCFold2Lasso:  checking RFastestCFoldValidate2Lasso"); 
    flush.console();
  }   
  if (Verbose > 1) {
    print(paste("LoadedSampleInputs, length(LambdaAK) == ", LambdaAK, sep=""));
    flush.console();
  }
  if (!is.null(RFastestCFoldValidate2Lasso) && 
    ((is.numeric(RFastestCFoldValidate2Lasso) && RFastestCFoldValidate2Lasso > 0) ||
      (is.logical(RFastestCFoldValidate2Lasso) && RFastestCFoldValidate2Lasso == TRUE))) {
    print("FastestCFoldValidate2Lasso: Inputs Blank"); flush.console();
    print(paste("RFastestCFoldValidate2Lasso = ",
      RFastestCFoldValidate2Lasso, sep="")); flush.console();
  } 
    
  if (!is.null(RFastestCFoldValidate2Lasso) && 
    (RFastestCFoldValidate2Lasso > 2)) {
    if (!exists("ii")) { ii = 0; }
      print(paste("Getting Models for  = ", ii, 
        sep="")); flush.console();
  } 

    
  if (!is.null(RFastestCFoldValidate2Lasso) && 
    (RFastestCFoldValidate2Lasso > 0 ||
      RFastestCFoldValidate2Lasso == TRUE)) {
    print("RFastestCFoldValidate2Lasso: Inputs Declared"); flush.console();
  }  
  OrderSeq = .M2Stuff$OrderSeq;
  try(eval(parse(text=GetG0Text(".FastestCFoldValidate2Lasso",S=1))), silent=TRUE);
  .FastestCFoldValidate2Lasso <- RFastestCFoldValidate2Lasso;
  try(eval(parse(text=LockGText(".FastestCFoldValidate2Lasso", S=1))), silent=TRUE);
  
  if (Verbose > 1) {
    print("Setting TwoLassoStayFitFunction"); flush.console();
  }
  try(eval(parse(text=GetG0Text("TwoLassoStayFitFunction", "globalenv()", S=1))), silent=TRUE);
  TwoLassoStayFitFunction <- function(X,Y,DoLogit=FALSE, SampleInput,StartBeta=-1,
    Verbose = 0, RFastestCFoldValidate2Lasso = .FastestCFoldValidate2Lasso, StandardFlag=FALSE, ...) {
    if (!exists("RFastestCFoldValidate2Lasso") || 
      is.null(RFastestCFoldValidate2Lasso) ) {
      RFastestCFoldValidate2Lasso = Verbose;
      eval(parse(text=SetGText("RFastestCFoldValidate2Lasso", S=1)));    
    }  else if (!exists(".FastestCFoldValidate2Lasso") || 
      is.null(.FastestCFoldValidate2Lasso) ) {
      .RFastestCFoldValidate2Lasso = RFastestCFoldValidate2Lasso;    
      eval(parse(text=SetGText(".FastestCFoldValidate2Lasso", S=1)));
    }
    if (!is.null(RFastestCFoldValidate2Lasso) && RFastestCFoldValidate2Lasso > 2) {
      print(paste("RFastestCFoldValidate2Lasso: TwoLassoStayFitRun SI: (",
        paste(SampleInput, collapse=", "),")", sep="")); flush.console();
    }
    eval(parse(text=MyExistsF("OrderSeq", "c(20,1)")));
    eval(parse(text=GetG0Text(".M2Stuff", "globalenv()", S=1)));
    LambdaAKInputs = .M2Stuff$LambdaAKInputs;
    LambdaDKInputs = .M2Stuff$LambdaDKInputs;   
    LambdaAK = .M2Stuff$LambdaASeq;
    LambdaDK = .M2Stuff$LambdaDSeq; 
   
   if (all(Y %in% c(0.0,1.0))) {
     if (DoLogit==FALSE) {
      print("CrossValidatorTwoLasso.r::TwoLassoStayFitFunction() why is DoLogit FALSE?"); flush.console();
     }
   }
   if (is.null(LambdaDKInputs) || length(LambdaDKInputs) <= 0) {
     print("TwoLassoStayFitFunction;  Uh oh, LambdaDKInputs is Null.");
     flush.console();
   }
   if (is.null(LambdaAKInputs) || length(LambdaAKInputs) <= 0) {
     print("TwoLassoStayFitFunction;  Uh oh, LambdaAKInputs is Null.");
     flush.console();
   }
   if (Verbose > 2) {
     print("Launching Two2LassoMiniValidate")
   }
   try(eval(parse(text=GetG0Text(".HiddenTwoLassoOb", "globalenv()"))), silent=TRUE);  
   .HiddenTwoLassoOb <- TwoLassoCpp:::.Two2LassoCppMiniValidate(X = X, Y = Y, 
      DoLogit=DoLogit,
      PiAVectorInputs = SampleInputs[,1], 
      SigmaVectorInputs = SampleInputs[,2], 
      LambdaDK = LambdaDK, LambdaAK = LambdaAK, 
      LambdaDKInputs = LambdaDKInputs, LambdaAKInputs = LambdaAKInputs,
      OrderSeq = OrderSeq,  MaxCauchy = MaxCauchy * .5, 
      CauchyEpsilon = CauchyEpsilon * sqrt(p), 
      StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
      XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
      InitKKs = -5, WLSWeights = -1, TDFNu = -1,
      Verbose = RFastestCFoldValidate2Lasso-3, RecordFlag = 0,  Groupers = NULL,
      L2ShrinkagePrior = L2ShrinkagePrior, 
      RecordL2Shrinkage = RecordL2Shrinkage,StartBetaMatrix=StartBetaMatrix,
      RunFlag = TRUE,  ...)
    try(eval(parse(text=LockGText(".HiddenTwoLassoOb", "globalenv()"))), silent=TRUE);
    return(.HiddenTwoLassoOb$RecordBetaCVFinish);
  }
  try(eval(parse(text=LockGText("TwoLassoStayFitFunction", "globalenv()"))), silent=TRUE);
  if(Verbose > 0) {
    print("Declaring RunOnAChange"); flush.console();
  }
  try(parse(eval(text=GetG0Text("RunOnAChange", "globalenv()",S=1))), silent=TRUE);
  RunOnAChange <- function(XTrain, YTrain,...) {  
  }
  try(parse(eval(text=SetGText("RunOnAChange", "gloablenv()"))), silent=TRUE);
  
  if (Verbose > 0) {
    print("RFastestCFoldValidate2Lasso: Run CFoldValidate")
  }
  try(eval(parse(text=SetGText("SampleInputs", "globalenv()", S=1))));
  try(eval(parse(text=LockGText("MyCFold", "globalenv()", S=1))), silent=TRUE);
  MyCFold = CFoldValidate(X,Y,DoLogit=DoLogit,TwoLassoStayFitFunction,SampleInputs,
    cFolds = cFolds, RunOnAChange = RunOnAChange,
    CFoldVerbose=RFastestCFoldValidate2Lasso -2,
    FitEmAll = TRUE, DoLogit=DoLogit);
  try(parse(eval(text=SetGText("MyCFold", "gloablenv()"))), silent=TRUE);    
 
  M2Stuff <- M2Items(X, sigma=MyCFold$BestSampleInputs[2], 
      OrderSeq = OrderSeq, SABS = SABS);
  NewLambdaAK = M2Stuff$LambdaASeq;  NewLambdaDK = M2Stuff$LambdaDSeq;
    NewOrderSeq = M2Stuff$OrderSeq;  
  if (Verbose > 0) {
    print("RFastestCFoldValidate2Lasso: SetFinalM2Out"); flush.console();  
  } 
  try(eval(parse(text=GetG0Text("FinalM2Out", "globalenv()", S=1))), silent=TRUE);
  FinalM2Out  <- Two2LassoCpp(X = X, Y = Y, DoLogit=DoLogit,
    PiA = MyCFold$BestSampleInputs[1], SigmaSq = MyCFold$BestSampleInputs[2], 
    LambdaDK = NewLambdaDK, LambdaAK = NewLambdaAK, OrderSeq = NewOrderSeq,
    MaxCauchy = MaxCauchy, CauchyEpsilon = CauchyEpsilon, 
    StartBeta = rep(0, length(X[1,])), StandardFlag = 0,
    XtX = -1, XtY = -1, InverseGammaConstant = 1, FixKa = -100, 
    InitKKs = -5, WLSWeights = -1, TDFNu = -1,
    Verbose = RFastestCFoldValidate2Lasso, RecordFlag = 0, HoldOn = TRUE, 
    Groupers = NULL,
    L2ShrinkagePrior = L2ShrinkagePrior, RecordL2Shrinkage = RecordL2Shrinkage,     
    RunFlag = 1);
  FinalM2Out$MyCFold = MyCFold;
  try(eval(parse(text=SetGText("FinalM2Out", "globalenv()"))), silent=TRUE);
  if (Verbose > 0) {
    print("RFastestCFoldValidate2Lasso: All done, FinalM2Out"); flush.console();
  }
  return(FinalM2Out);
}


CFoldValidate2Lasso <- function(X,Y, nPi = 20, nSigma = 10, cFolds = 5,
  Verbose = 0, DoLogit = FALSE) {
  FitFunction = TwoLassoFitFunction;
  SampleInputs = TwoLassoSampleInputs(X,Y,nPi, nSigma, DoLogit=DoLogit)
  MyCFold = CFoldValidate(X,Y,FitFunction, SampleInputs, cFolds = cFolds, DoLogit=DoLogit);
  
  FinalM2Out = M2Lasso(xxs = X, yys=Y, DoLogit=DoLogit,
    ppiuse = MyCFold$BestSampleInputs[1], 
    sigmaNoiseSq = MyCFold$BestSampleInputs[2]);
  FinalM2Out$MyCFold = MyCFold;
  return(FinalM2Out);

}

##############################################################
## TwoLassoFitFunction
##
##  Given X,Y, SampleInput (which is PiA, SigmaSq)
##  this orders M2Lasso to simply return the output Betas.
TwoLassoFitFunction <- function(X,Y, DoLogit=FALSE, SampleInput, StartBeta = -1) {
  PiAInput = SampleInput[1];  SigmaSqInput = SampleInput[2];
  Out <- M2Lasso(xxs=X, yys=Y, DoLogit=DoLogit, ppiuse = PiAInput,
   sigmaNoiseSq = SigmaSqInput, StartBeta = StartBeta);
  return(Out$ReturnBetas);
}
TwoLassoSampleInputs <- function(X,Y, nPi=1, nSigma=1, DoLogit=FALSE) {
  p = length(X[1,]); n = length(Y); 
  if (DoLogit==TRUE) { nPi = nPi * nSigma;  nSigma = 1; }
  if (p >= 10000) {
    LogitScale <- log(500) * (1:nPi - nPi/2) / (nPi/2);
   ##if (n < p) {
   ##  LogitScale = LogitScale + log(n/p)/(1-log(n/p));
   ##}
   Pis = exp(LogitScale - log(p)) / ( 1 + exp(LogitScale-log(p))); 
   if (n < p) {
     Pis = Pis * (n/p);
   }   
  } else {
  LogitScale <- (1:nPi - nPi/2) / (nPi/2)^(1/3);
  ##if (n < p) {
  ##  LogitScale = LogitScale + log(n/p)/(1-log(n/p));
  ##}
  Pis = exp(LogitScale) / ( 1 + exp(LogitScale));
  if (n < p) {
    Pis = Pis * (n /p);
  }
  }

  if (DoLogit==TRUE) {
    SigmaValues=1; 
    if (n < p) { SigmaValues = n/p; }
  } else if (p >= 10000) {
    SigmaBig = var(Y);
    SigmaValues = SigmaBig * 2^(-1 -4*((1:nSigma)) / (nSigma) )     
  } else {
    SigmaBig = var(Y);
    SigmaValues = SigmaBig * 2^(-1 -((nSigma:1)) / (nSigma)^(1/3) )
    if (2*n < p) {
      SigmaValues = SigmaValues * 2^( log((2*n)/p) )
    }
  }
  SampleSet = cbind( rep(Pis, nSigma), rep(SigmaValues, each=nPi) );
  return(SampleSet);
} 

CFoldValidateInitiateText <- function() {
  return("
  if (!exists(\"CFoldVerbose\")) { CFoldVerbose = 0; }
    if (NCOL(X) >= 50000) {
    CFoldVerbose = CFoldVerbose + 3;
  }
  if (CFoldVerbose > 0)  {
    print(paste(\"CrossValidatorTwoLasso.r:::CFoldValidate[Verb=\", CFoldVerbose, \"]: Running\", sep=\"\")); flush.console();
  }
  eval(parse(text=TwoLassoCpp:::MyExistsF(\"cFolds\", 5)));
  eval(parse(text=TwoLassoCpp:::MyExistsF(\"CFoldVerbose\", 2)));
  eval(parse(text=TwoLassoCpp:::MyExistsF(\"FitEmAll\", FALSE)));    
  eval(parse(text=TwoLassoCpp:::MyExistsF(\"DoRandom\", FALSE))); 
  if (!exists(\"FitFunction\")) {
    eval(parse(text=GetG0Text(\"TwoLassoStayFitFunction\", \"globalenv()\",S=1)));
    FitFunction =  TwoLassoStayFitFunction;
  }


  if (CFoldVerbose > 0) {
    print(paste(\"CrossValidatorTwoLasso.r::CFoldValidate: Setting Y to center. \", sep=\"\")); flush.console();
  }
  Y = Y - mean(Y); X = t( t(X) - colMeans(X) ); #
  if (CFoldVerbose > 0) {
    print(paste(\"CrossValidatorTwoLasso.r::CFoldValidate: succesfully set Y to center. \", sep=\"\")); flush.console();
  }
  ### gotta be mean standardized here because of the cross validation
  RandomSetOrder = sample(1:length(Y), size = length(Y), replace=FALSE);
  Brks = ceiling( (1:(cFolds)) * (length(Y) / cFolds));
  SetsOrdered = list();   SetsOrdered[[1]] = 1:(Brks[1]) 
  SetsOrdered[[1]] = 1:Brks[1];
  if (length(Brks) >= 2) {
    for (ii in 2:(length(Brks))) {
      SetsOrdered[[ii]] = (Brks[ii-1]+1):(Brks[ii]);
      SetsOrdered[[ii]] = SetsOrdered[[ii]][SetsOrdered[[ii]] >= 1 & SetsOrdered[[ii]] <= length(Y)]   
    }
  }
  NewSetsOrdered= list();
  for (ii in 1:length(SetsOrdered)) {
    if (length(SetsOrdered[[ii]]) > 0) {
      NewSetsOrdered[[length(NewSetsOrdered)+1]] <- SetsOrdered[[ii]];
    }
  }
  if (length(NewSetsOrdered) == 0) {
    print(paste(\"CrossValidatorTwoLasso.r: ERROR No we cannot fold \", cFolds, 
      \" times for length Y = \", length(Y), sep=\"\")); flush.console();
    return(-1);
  }
  SetsOrdered = NewSetsOrdered;
  
  p = NCOL(X);
  
  if (CFoldVerbose > 0) {
    print(paste(\"CrossValidatorTwoLasso.r::CFoldValidate(): Now Now Allocate MSMatrix. \", sep=\"\")); flush.console();
  }
  MSMatrix <- NULL;
  try(MSMatrix <- matrix(0, NROW(SampleInputs), cFolds));
  try(MatrixNonZero <- matrix(0, NROW(SampleInputs), cFolds));
  if (is.null(MSMatrix) || length(MSMatrix) != NROW(SampleInputs) * cFolds) {
    print(paste(\"CrossValidatorTwoLasso.r::CFoldValidate() failure to allocate MSMatrix[\",
      NROW(SampleInpust), \", \", CFolds, \"]\", sep=\"\"));
    print(\"ERROR: In MSMatrix allocation\");
    return(-999);
  }
  StartBeta = -1;
  EveryFitBetas <- list();
  ");  
}
CFoldValidateInFoldStart <- function() {
  return("
  if (!exists(\"Foldtt\")) { Foldtt = 1; }
 if (CFoldVerbose > 0)  {
      print(paste(\"CrossValidatorTwoLasso.r::CFoldValidate()[Verb=\", CFoldVerbose, 
        \"]: Doing Foldtt = \", Foldtt, \"/\", cFolds, sep=\"\"));
      flush.console();
    }
    if (DoRandom == TRUE) {
      TestSet = RandomSetOrder[SetsOrdered[[Foldtt]]];
    } else {
      TestSet = SetsOrdered[[Foldtt]];
    }
    if (CFoldVerbose > 0)  {
      print(paste(\"CrossValidatorTwoLasso.r::CFoldValidate()[Verb=\", CFoldVerbose, 
        \"]: Setting TrainSet, Testyy for Foldtt =  \", Foldtt, sep=\"\"));
      flush.console();
    }
    TrainSet = (1:length(Y))[ !( (1:length(Y)) %in%  TestSet) ];
    TestY = Y[TestSet];  TrainY = Y[TrainSet];
    TestX = X[TestSet,];  TrainX = X[TrainSet,];
    if (!is.null(RunOnAChange)) {
      RunOnAChange(TrainX, TrainY);
    }  
  ");
}
CFoldValidate <- function(X,Y,FitFunction, SampleInputs, cFolds = 5,
  DoRandom = FALSE, RunOnAChange = NULL,CFoldVerbose = 0, FitEmAll = FALSE,
  DoLogit = FALSE) {
  eval(parse(text=CFoldValidateInitiateText()));
  for (Foldtt in 1:cFolds) {
    eval(parse(text=CFoldValidateInFoldStart()));
    if (FitEmAll == TRUE) {
       if (CFoldVerbose > 0)  {
         print(paste("CrossValidatorTwoLasso.r::CFoldValidate()[Verb=", CFoldVerbose, 
           "]: Running FitFunction for Fit BetasAll  ", Foldtt, sep=""));
         flush.console();
       }
       FitBetasAll = FitFunction(TrainX, TrainY, SampleInputs, StartBeta =StartBeta,
         StandardFlag=FALSE, TestX=TestX)
       eval(parse(text=SetGText("FitBetasAll", "globalenv()", S=1)));
       eval(parse(text=SetGText("TestX", "globalenv()", S=1)));
       eval(parse(text=SetGText("TrainX", "globalenv()", S=1)));
       eval(parse(text=SetGText("SampleInputs", "globalenv()", S=1)));
       eval(parse(text=SetGText("StartBeta", "globalenv()", S=1)));
       eval(parse(text=SetGText("FitFunction", "globalenv()", S=1)));
       eval(parse(text=SetGText("TestY", "globalenv()", S=1)));
       if (is.null(FitBetasAll) || !is.numeric(FitBetasAll)) {
         print("Oh No, FitBetasAll is Not a vector!"); flush.console();
       }
       if (NROW(FitBetasAll) != NCOL(TrainX)) {
         print(paste("CrossValidatorTwoLasso.r::FitBetasAll():: Oh Oh Oh FitBetasAll length ", length(FitBetasAll),
           "  but NCOL(TrainX) = ", NCOL(TrainX), sep="")); flush.console();
       }
       if (NROW(FitBetasAll) != NCOL(TestX)) {
         print(paste("Oh Oh Oh FitBetasAll length ", length(FitBetasAll),
           "  but NCOL(TrainX) = ", NCOL(TrainX), sep="")); flush.console();
       }
       if (is.list(FitBetasAll)) {
         EstYs <- FitBetasAll$FitY;
         FitBetasAll <- FitBetasAll$FitBeta;
       } else {
         EstYs = TestX %*% FitBetasAll;
       }
       if (DoLogit == FALSE) {
         MSMatrix[,Foldtt] = colMeans( (EstYs - TestY)^2 );
       } else {
         PPY <- EstYs -log(1.0+exp(EstYs));
         PPY[EstYs > 10] = EstYs[EstYs > 10];
         PPY[EstYs < -10] = EstYs[EstYs < -10];
         MPPY <- log(1.0-exp(PPY));
         MPPY[EstYs > 10] <- -EstYs[EstYs > 10]
         MPPY[EstYs < -10] <- -EstYs[EstYs < -10];
         MSMatrix[,Foldtt] =  colSums(PPY*TestY + MPPY*(1-TestY));
       }
       StartBeta = FitBetasAll[,1]
       EveryFitBetas[[Foldtt]] <- FitBetasAll;
       for (kk in 1:NCOL(FitBetasAll)) {
         MatrixNonZero[kk,Foldtt] = length(FitBetasAll[FitBetasAll[,kk] != 0.0,kk]);
       }
       if (CFoldVerbose > 0)  {
         print(paste("CrossValidatorTwoLasso.r::CFoldValidate(): Foldtt =  ", Foldtt,
           " MSMatrix[1,Foldtt] = ", MSMatrix[1,Foldtt], "done setting EveryFitBetas.", sep=""));
         flush.console();
       }
    } else {
     ## FitBetasAll <- matrix(0, NCOL(TestX), length(SampleInputs[,1]));  
     ## for (Inputii in 1:length(SampleInputs[,1])) {
        if (CFoldVerbose > 2)  {
          print(paste("CFoldValidate[Verb=", CFoldVerbose, "]: Doing Foldtt = ", Foldtt, 
            ", Inputii = ", Inputii, sep=""));
          flush.console();
        }
       if (CFoldVerbose > 0)  {
         print(paste("CFoldValidate[Verb=", CFoldVerbose, 
           "]: Foldtt =  ", Foldtt, " Run FitFunction on TrainX, TrainY.", sep=""));
         flush.console();
       }
        FitBetas = FitFunction(TrainX, TrainY, SampleInputs, StartBeta =StartBeta,
          StandardFlag = FALSE,)
        if (DoLogit == TRUE) {
          EstY = TestX %*% matrix( FitBetas, p,1);
          MSMatrix[, Foldtt] = colMeans((EstY-TestY)^2);
        } else {
         PPY <- EstYs -log(1.0+exp(EstYs));
         PPY[EstYs > 10] = EstYs[EstYs > 10];
         PPY[EstYs < -10] = EstYs[EstYs < -10];
         MPPY <- log(1.0-exp(PPY));
         MPPY[EstYs > 10] <- -EstYs[EstYs > 10]
         MPPY[EstYs < -10] <- -EstYs[EstYs < -10];
         MSMatrix[,Foldtt] =  colSums(PPY*TestY + MPPY*(1-TestY));
        }      
        MatrixNonZero[1,Foldtt] = length(FitBetas[FitBetas != 0.0]);
        StartBeta = FitBetas;
        ##FitBetasAll[, Inputii] <- FitBetas;
      ##}
      EveryFitBetas[[Foldtt]] <- FitBetas;    
       if (CFoldVerbose >= 2)  {
         print(paste("CFoldValidate: Foldtt =  ", Foldtt, " SSq[1] was ", 
           MSMatrix[1, Foldtt],
           "again done setting EveryFitBetas.", sep=""));
         flush.console();
       }
    }
  }
  MSV = rowMeans(MSMatrix);
  BestSampleInputs = SampleInputs[ sort(MSV, index=TRUE)$ix[1], ];
  if (CFoldVerbose >= 1) {
    print(paste("CrossValidatorTwoLasso.r::CFoldValidate: We have finished all folds. ", 
      sep="")); 
    flush.console();
  }
  return(list(BestSampleInputs = BestSampleInputs, MSV = MSV, 
    MSMatrix = MSMatrix, DoRandom = DoRandom, SetsOrdered = SetsOrdered,
      RandomSetOrder = RandomSetOrder, 
      SampleInputs = SampleInputs, EveryFitBetas = EveryFitBetas,
      MatrixNonZero = MatrixNonZero, DoLogit=DoLogit));
  ##RealSets = Set
}