#############################################################################
##  2009-2011 Alan Lenarcic
##
##     Implementation of this alternate estimator as requested by a suggestion
##  from outside request, did not quite achieve what it was hoped to.
##  Using Coordinate descent against an integrated penalty proved to require
##  a bit more intuition than the Bayesian EM version used in other penalties.
##  It is believed that a continuous implementation of a psy function
##  would make the work on this estimator easier
##
##    In any case, this estimator had multiple free parameters, and such
##  the method to conduct cross validation did not have success.
##  Cross validation likely requires more experience from using east estimator
##  to better search relevant space for anestimator.
##
## 
SaveDataOnce <- function(OnLambda, OnGammaSq) {
   MyDirs = unlist(searchpaths());
   if (any( substr(MyDirs, nchar(MyDirs) - nchar("TwoLasso")+1, nchar(MyDirs))
     == "TwoLasso")) {
     PathOfPackage = MyDirs[
       substr(MyDirs, nchar(MyDirs) - nchar("TwoLasso")+1, nchar(MyDirs))
     == "TwoLasso" ];    
   }
   WholePathData = paste(PathOfPackage, "//data", sep="");
   dir.create(WholePathData, showWarnings=FALSE);
   NameFile = NameNEGRData(OnLambda, OnGammaSq);
   ListFiles = list.files(WholePathData);
   if (any(ListFiles == NameFile)) {
      load( paste(WholePathData, "//", NameFile, sep=""), .GlobalEnv);
      if (!exists("DensBeta") || is.null(DensBeta)) {
         DensBeta = SIntegral / 
           ( .5 * sum( SIntegral[2:(length(SIntegral)-1)] * 
             (Beta[3:(length(SIntegral))] - Beta[1:(length(SIntegral)-2)] ) ) +
            SIntegral[1] *(Beta[2] - Beta[1]) +
            SIntegral[length(SIntegral)] * (Beta[length(SIntegral)] -
              Beta[length(SIntegral)-1])
           )
      }
      return(list(Beta=Beta, QuadraticPen = QuadraticPen,
        SIntegral=SIntegral, SInvSum = SInvSum,
        SSecondInvSum = SSecondInvSum, DensBeta = DensBeta));
   } else {
      AllBetaValues =  .5 * exp((-400:300)/50);
      LessOneBetaValues = AllBetaValues[AllBetaValues <= 1];
      MoreOneBetaValues = AllBetaValues[AllBetaValues > 1];
      
      LogMaxBetaSq = log( max(AllBetaValues)^2 ) / log(10);
      NPsi = 49000;
      PsiValues = c((1:9)/10000,(1:999)/1000, 10^( 5 * LogMaxBetaSq * (0:NPsi) / NPsi) );
      logOutInvPsiValues = -log(min(PsiValues)) + (0:50)/10
      On2 = MakeQuadraticPen(OnLambda, OnGammaSq, VerboseInt = 25,
        Beta = MoreOneBetaValues, PsiValues = PsiValues,
        logOutInvPsiValues = logOutInvPsiValues);
        
      NPsi = 50000;
      PsiValues = c((1:9)/10000,(1:999)/1000, 1+1000 * (1:NPsi) / (NPsi) );
      DoubleFactor = 2*abs( min(log(AllBetaValues)));
      logOutInvPsiValues = -log( min(PsiValues)   ) + 
        ( DoubleFactor + 4 + log( min(PsiValues) )) *(0:NPsi) / NPsi;
      On1 = MakeQuadraticPen(OnLambda, OnGammaSq, VerboseInt = 25,
        Beta = LessOneBetaValues, PsiValues = PsiValues,
        logOutInvPsiValues = logOutInvPsiValues );
      Beta = c(On1$Beta, On2$Beta);
      QuadraticPen = c(On1$QuadraticPen, On2$QuadraticPen);
      SInvSum = c(On1$SInvSum, On2$SInvSum);  
      SSecondInvSum = c(On1$SSecondInvSum, On2$SSecondInvSum);
      SIntegral = c(On1$SIntegral, On2$SIntegral);
      DensBeta = SIntegral / 
        ( .5 * sum( SIntegral[2:(length(SIntegral)-1)] * 
          (Beta[3:(length(SIntegral))] - Beta[1:(length(SIntegral)-2)] ) ) +
          SIntegral[1] *(Beta[2] - Beta[1]) +
          SIntegral[length(SIntegral)] * (Beta[length(SIntegral)] -
           Beta[length(SIntegral)-1])
        )
      save(Beta=Beta, QuadraticPen=QuadraticPen, 
        SIntegral=SIntegral, SInvSum = SInvSum,
        SSecondInvSum = SSecondInvSum, 
        file=paste(WholePathData, "//", NameFile, sep="") ) 
      return(list(Beta=Beta, QuadraticPen = QuadraticPen,
        SIntegral=SIntegral, SInvSum = SInvSum,
        SSecondInvSum = SSecondInvSum)); 
   }
}


NameNEGRData <- function(OnLambda, OnGammaSq) {
 paste("OnL", tSeq(OnLambda), "OnG", tSeq(OnGammaSq), ".Rdata", sep="")
}


MakeQuadraticPen <- function(OnLambda, OnGammaSq, VerboseInt = 100,
  Beta = .1 * exp((-450:450)/40), PsiValues = (1:100000) / 100,
  logOutInvPsiValues = log(100) + (0:10000)/10) {   
   QuadraticPen = Beta * 0;
   SIntegral = Beta * 0; SSecondInvSum = Beta*0; SInvSum = Beta*0;
   
   .Call("NEGLassoMakeTableShell", QuadraticPen, PsiValues, logOutInvPsiValues,
     Beta, OnLambda, OnGammaSq, SIntegral,  SInvSum,
     SSecondInvSum, VerboseInt) 
   return(list(QuadraticPen = QuadraticPen, Beta = Beta,
     SIntegral = SIntegral, SInvSum = SInvSum, SSecondInvSum = SSecondInvSum,
     OnLambda = OnLambda, OnGammaSq = OnGammaSq,
     PsiValues=PsiValues, logOutInvPsiValues=logOutInvPsiValues));
} 


# On = MakeInvPsiOfBeta(OnLambda = .5, OnGammaSq = 1, VerboseInt= 100);
# plot( On$InvPsiValues~ On$Beta, type="l");
# plot( log(On$InvPsiValues) ~ log(On$Beta), type="l");

SuperNEGLasso <- function(X,Y, OnLambda = .5, OnGammaSq = 1,
  OnBeta = -1, NRandomStarts = 5,
  OnGammas = -999, InitKKs = -5,
  NumCDOConv = -999, CDOEpsilon = -999, 
  NumEMSteps = 10, EMCauchy=.00001,
  TotalLoops = 100, MaxEpsilon=.00001,
  StartBeta = NULL, Verbose = 0) {
  SVDX = NULL;
  try(SVDX <- svd(X));
  if (is.null(SVDX)) {
    print("SuperNEG: Cannot work with NULL svd of X"); flush.console();
    return(NULL);
  }
  PenaltyInfo = SaveDataOnce(OnLambda, OnGammaSq);
  DensBeta = PenaltyInfo$DensBeta;
  if (Verbose > 0) {
    print("SuperNEGLasso: Starting "); flush.console();
  }    
  A = SVDX$v;
  library(corpcor);
  MLLSBeta <- NULL;
    try(MLLSBeta <- pseudoinverse( t(X) %*% X )  %*% t(X) %*% Y);
  if (Verbose > 0) {
    print("SuperNEGLasso: Got MLLSBeta! "); flush.console();
  } 
  if (is.null(MLLSBeta)) {
    try( MLLSBeta <- (FitNEGLasso(X=X, Y=Y, OnLambda=OnLambda, 
      OnGammaSq= OnGammaSq, OnBeta = -1, OnGammas= OnGammas, 
      InitKKs = InitKKs, NumCDOConv = NumCDOConv, CDOEpsilon = CDOEpsilon, 
      NumEMSteps=NumEMSteps, EMCauchy = EMCauchy, 
      TotalLoops = TotalLoops, MaxEpsilon = MaxEpsilon,
      StartBeta = StartBeta, Verbose=Verbose-1))$ReturnBetas ); 
    if (is.null(MLLSBeta)) {
      print("SuperNEG: MLSSbeta is Still Null \n"); 
      flush.console();return(NULL);
    }
  }
  MyFits = list();
  Scores = rep(0, NRandomStarts);
  ReturnBetas = matrix(0, NRandomStarts, length(X[1,]) );
  for (ii in 1:NRandomStarts) {
    Z = rnorm( length(A[,1]) );
    TryBeta = MLLSBeta + Z - A %*% t(A) %*% Z;
    RT1 = NULL;
    if (Verbose > 0) {
      print(paste("SuperNEGLasso: Fitting ",ii, sep="")); flush.console();
    }
    try( RT1 <- FitNEGLasso(X=X, Y=Y, OnLambda=OnLambda, OnGammaSq=OnGammaSq,
      StartBeta = TryBeta, OnGammas=OnGammas, InitKKs=InitKKs, 
      NumCDOConv=NumCDOConv, CDOEpsilon=CDOEpsilon, 
      NumEMSteps, EMCauchy=EMCauchy, TotalLoops=TotalLoops, MaxEpsilon=MaxEpsilon,
      Verbose = Verbose -1) ); 
    if (is.null(RT1) || length(RT1) == 1) {
       print("SuperNEGLasso: Failed to get a fit at this one");
    } else {
      if (Verbose > 0) {
        print(paste("SuperNEGLasso: Succesful fit for ii = ", ii, sep=""));
        flush.console();
      }
      MapBackReturnBetas = .MapBackBeta(PenaltyInfo$Beta, abs(RT1$ReturnBetas))
      Score = .5 * sum( (Y - X %*% RT1$ReturnBetas)^2 ) -
        sum(log(DensBeta[MapBackReturnBetas]));
      RT1$Score = Score;
      ReturnBetas[ii,] = RT1$ReturnBetas;
      Scores[ii] = Score;
      if (Verbose > 0) {
        print(paste("SuperNEGLasso: Fit was ", Score, sep=""));
        print(paste(" CurrentScores: ", paste(Scores, collapse=", "), sep=""));
        flush.console();
      }
    }
    MyFits[[ii]] = RT1;
  }
  if (Verbose > 0) {
    print("SuperNEGLasso: AllDone  "); flush.console();
  }
  return(list(ReturnBetas = ReturnBetas, Scores=Scores, AllFits=MyFits));  
}

FitNEGLasso <-function(X, Y, OnLambda = .5, OnGammaSq = 1,
  OnGammas = -999, InitKKs = -5,
  NumCDOConv = -999, CDOEpsilon = -999, 
  NumEMSteps = 10, EMCauchy=.00001,
  TotalLoops = 100, MaxEpsilon=.00001,
  StartBeta = NULL, Verbose = 0) {
  if (length(X) == 1) {
    print("FitNEGLasso: X is not good"); return; 
  }
  if (length(Y) == 1) {print("FitNEGLasso: Y is not Good")}
  PenaltyInfo = SaveDataOnce(OnLambda, OnGammaSq);
  StartPenii = .MapBackBeta(PenaltyInfo$Beta, 1);
  StartPenalty = PenaltyInfo$QuadraticPen[StartPenii] * .5;
  p = length(X[1,]);
  if (Verbose > 0) {
    print("NEGLasso: Starting to run first CDO");
    flush.console();
  }
  
  if (Verbose > 1) {PrintFlag = Verbose -1;} else {PrintFlag = -1;}
  ##  We find LaPaplace distribution with variance matching Psi
  ##  Since 2/Lambda = 1/ Gamma = Psi, this is LaPlace with 1/Psi
  FirstCDO = NULL;
  if (is.null(StartBeta)  || length(StartBeta) != p ) {
    FirstCDO = CoordinateDescent(xx = X, yy = Y, 
      OnBeta = rep(0, length(X[1,])), 
      OnGammas = rep(StartPenalty, p),
      InitKKs=InitKKs, NumCDOConv=NumCDOConv, CDOEpsilon=CDOEpsilon,
      TotalLoops=TotalLoops, MaxEpsilon=MaxEpsilon, Verbose = PrintFlag);
     OnBeta = FirstCDO$ReturnBetas;
  } else {
    OnBeta = StartBeta;   
  }
  LinBeta = PenaltyInfo$QuadraticPen[
    .MapBackBeta(PenaltyInfo$Beta,abs(OnBeta)) ]   
  for (jj in 1:NumEMSteps) {
    NewCDO = CoordinateDescent(xx = X, yy = Y, 
      OnBeta = OnBeta, 
      OnGammas = LinBeta,
      InitKKs=InitKKs, NumCDOConv=NumCDOConv, CDOEpsilon=CDOEpsilon,
      TotalLoops=TotalLoops, MaxEpsilon=MaxEpsilon, Verbose = PrintFlag);  
    NewBeta = NewCDO$ReturnBetas;
    if (sum(abs(NewBeta-OnBeta)) <= EMCauchy) {
      OnBeta = NewBeta;
      break;
    }
    OnBeta = NewBeta;
    LinBeta = PenaltyInfo$QuadraticPen[
      .MapBackBeta(PenaltyInfo$Beta,abs(OnBeta)) ];
  }
  LinBeta = PenaltyInfo$QuadraticPen[
      .MapBackBeta(PenaltyInfo$Beta,abs(OnBeta)) ];
  RetList = list( ReturnBetas = OnBeta, LinearPenalty = LinBeta,
    CDOOut = NewCDO, FirstCDO = FirstCDO);
  return(RetList);
}
##CoordinateDescent(xx = -1, yy = -1, XTX = -1, XTY = -1, NLen = -1,
##   TotalLoops = -999, MaxEpsilon = -999, OnBeta = -1, OnGamma = -1, OnLambda = -1,
##   RecordBetasFlag = FALSE, OnGammas = -999, InitKKs = -5,
##   NumCDOConv = -999, CDOEpsilon = -999, WLSWeights = -1)


.MapBackBeta <- function(BetaList, FindBeta) {
  LinBeta = log(BetaList);
  logFindBeta = log(FindBeta);
  ##LinBeta = LinBeta[1] + (0:( length(LinBeta)-1)) * ( 
  ##  LinBeta[length(LinBeta)] -LinBeta[1] ) / (length(LinBeta)-1);
  iiFind = round( (logFindBeta - LinBeta[1] ) * (length(LinBeta)-1) /
     ( LinBeta[length(LinBeta)] -LinBeta[1] ) + 1 );
  iiFind[iiFind < 1] = 1;
  iiFind[iiFind > length(LinBeta)] = length(LinBeta);
  return(iiFind); 
}
