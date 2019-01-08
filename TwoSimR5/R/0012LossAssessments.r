## 0012LossAssessments.r
##
##
##  Functions for Type 1, Hellinger, Sphere, AUC loss to meausre performance.
##

################################################################################
##  Sphere Loss
##
##    Measure the sphere loss between true (0,1,1,0,0,0,1,1...)
##  vector of true effects and the MIP vector (.1,.95,.8,.35,...)
##  that encodes certainty.  Ideal sphere distance is perfect 0,1 bets.
##
##  Let Sphere loss be 
##    sum_i   p_i^(1{i in model})*(1-p_i)^(1{i not in model})/ sqrt( p_i^2+(1-p_i^2))
##  If p_i was perfect zero or one, denominator for each weight is 1.  
##  If we were perfectly right, the total sphere score would be "p"
##  If we were perfectly wrong, total sphere score would be "0".
##  
##  If guess .75 for on and get it right, the denominator is .79
##  Thus we score a .948, which is pretty close to 1 for this score.
##  If we guess .75 for on and get it wrong (it is not in the model), 
##  we score a .314, which is much
##  better than zero.  When we vote .5, we will always score .707.
##
##  Thus a score of less .707 * p suggests less performance that just 
##  arbitrarily guessing.
SphereLossFunction <- function(MIPVec, BetaTrue) {
  Denom <- sqrt( MIPVec^2 + (1-MIPVec)^2);
   if (length(BetaTrue) < length(MIPVec) &&
    sum(abs(round(BetaTrue) - BetaTrue)) < .0001) {
    NotBetaTrue <- (1:length(MIPVec))[!( (1:length(MIPVec)) %in% BetaTrue) ];
    return( sum(MIPVec[BetaTrue]/Denom[BetaTrue]) +
      sum((1-MIPVec[NotBetaTrue])/ Denom[NotBetaTrue]) );
  } else if (length(BetaTrue )  != length(MIPVec)) {
    AFilePrint("AUC: Cant calculate because BetaTrue and MIPVec are different lengths!");
    return(-1);
  }
  BetaTrue <- (1:length(BetaTrue))[BetaTrue != 0.0];
  NotBetaTrue <- (1:length(MIPVec))[!( (1:length(MIPVec)) %in% BetaTrue) ];
    return( sum(MIPVec[BetaTrue]/Denom[BetaTrue]) +
      sum((1-MIPVec[NotBetaTrue])/ Denom[NotBetaTrue]) );
  return(Asum);
} 

##  AUCFunction()
##
##    AUC is "Area under Curve" of a ROC curve
##
##  This sorts the data and totals how many true effects enter the model
##  before the remainder leaves.
##
##  AUC is area under curve.
##
##  If MIP - model inclusion is given we order from largest inclusion to
##   smallest inclusion.  
##
##  AST be the order of iterices as they enter the model.
##  Let RKK be the positions that the the non-zero coordinates of BetaTrue enter the model
##
##   For instance if MIP is c(.9,.3,.6,.4,.5)
##    And the sort order is c(1, 3, 5, 4, 2)
##    Thus the position order is c(1,5,2,4,3)
##    And if Beta was       c(1, 1, 0, 0, 1)
##    Then  RKK is          c(1, 5, 3)
##
##    AUC is integrated number of True positives observed before false positives
##        We get (2 + 0+  1)/(2*3) or 50% AUC
##
##   MIPVec = c(.3,.35,.2,.7,.3,.5,.9,.2,.1,.85)
##   Sort Order is  c(7, 10, 4, 6, 2, 1, 5, 3, 8, 9)
##   Position Order is c(6, 5, 8, 3, 7, 4, 1, 9, 10, 2)
##   And if Beta was   c(1, 0, 0, 1,1,0,1,0,1,1)
##   RKK is           c(1, 2, 3, 6, 7, 10)
##   IBN <-  (1:length(MIPVec)
##     We get   (4 +4 + 4+ 2 + 2) / (4*6)
RankTotalIntegrateFunction <- function(MIPVec, BetaTrue) {
  if (length(BetaTrue) < length(MIPVec) &&
    sum(abs(round(BetaTrue) - BetaTrue)) < .0001) {
    AST <- sort(MIPVec, index=TRUE, decreasing=TRUE)$ix;
    POST <- sort(AST, index=TRUE)$ix;
    RKK <- sort(POST[BetaTrue]);
    Asum <- sum((length(MIPVec) - (RKK-1)))/
      (length(BetaTrue) * (length(MIPVec)));
    return(Asum);
  } else if (length(BetaTrue )  != length(MIPVec)) {
    AFilePrint("AUC: Cant calculate because BetaTrue and MIPVec are different lengths!");
    return(-1);
  }
  BetaTrue <- (1:length(BetaTrue))[BetaTrue != 0.0];
  AST <- sort(MIPVec, index=TRUE, decreasing=TRUE)$ix;
  POST <- sort(AST, index=TRUE)$ix;
  RKK <- sort(POST[BetaTrue]);
  Asum <- sum((length(MIPVec) - (RKK-1)))/
      (length(BetaTrue) * (length(MIPVec)));
  return(Asum);
} 

AUCFunction <- function(MIPVec, BetaTrue) {
  if (length(BetaTrue) < length(MIPVec) &&
    sum(abs(round(BetaTrue) - BetaTrue)) < .0001) {
    AST <- sort(MIPVec, index=TRUE, decreasing=TRUE)$ix;
    POST <- sort(AST, index=TRUE)$ix;
    RKK <- sort(POST[BetaTrue]);
    IBN <- (1:length(MIPVec))[!(1:length(MIPVec) %in% RKK)]
    MyR <- rep(0, length(RKK));
    for (ii in 1:length(RKK)) {
      MyR[ii] <- length(IBN[IBN > RKK[ii]]);
    }
    Asum <- sum(MyR )/
      (length(BetaTrue) * (length(MIPVec)-length(BetaTrue)));
    return(Asum);
  } else if (length(BetaTrue )  != length(MIPVec)) {
    AFilePrint("AUC: Cant calculate because BetaTrue and MIPVec are different lengths!");
    return(-1);
  }
  BetaTrue <- (1:length(BetaTrue))[BetaTrue != 0.0];
  AST <- sort(MIPVec, index=TRUE, decreasing=TRUE)$ix;
  POST <- sort(AST, index=TRUE)$ix;
  RKK <- sort(POST[BetaTrue]);
  IBN <- (1:length(MIPVec))[!(1:length(MIPVec) %in% RKK)]
  MyR <- rep(0, length(RKK));
  for (ii in 1:length(RKK)) {
    MyR[ii] <- length(IBN[IBN > RKK[ii]]);    
  }
  Asum <- sum( MyR )/
      (length(BetaTrue) * (length(MIPVec)-length(BetaTrue)));
  return(Asum);
} 

AUCFunctionOld <- function(MIPVec, BetaTrue) {
  if (length(BetaTrue) < length(MIPVec) &&
    sum(abs(round(BetaTrue) - BetaTrue)) < .0001) {
    AST <- sort(MIPVec, index=TRUE, decreasing=TRUE)$ix;
    RKK <- sort(AST[BetaTrue]);
    Asum <- sum(length(MIPVec) - RKK - (length(BetaTrue)-1:length(RKK)) )/
      (length(BetaTrue) * (length(MIPVec)-length(BetaTrue)));
    return(Asum);
  } else if (length(BetaTrue )  != length(MIPVec)) {
    AFilePrint("AUC: Cant calculate because BetaTrue and MIPVec are different lengths!");
    return(-1);
  }
  BetaTrue <- (1:length(BetaTrue))[BetaTrue != 0.0];
  AST <- sort(MIPVec, index=TRUE, decreasing=TRUE)$ix;
  POST <- sort(AST, index=TRUE)$ix;
  RKK <- sort(AST[BetaTrue]);
  Asum <- sum(length(MIPVec) - RKK - (length(BetaTrue)-1:length(RKK)) )/
      (length(BetaTrue) * (length(MIPVec)-length(BetaTrue)));
  return(Asum);
} 

#################################################################################
##  Type 1 Error
##     
##    How many false effects do we accept as True?
Type1Function <- function(MIPVec, BetaTrue) {
  if (length(BetaTrue) < length(MIPVec) &&
    sum(abs(round(BetaTrue) - BetaTrue)) < .0001) {
    NotBetaTrue <- (1:length(MIPVec))[!(1:length(MIPVec) %in% BetaTrue)];
    RT <- length(NotBetaTrue[MIPVec[NotBetaTrue] >= .5]);
    return(RT);
  } else if (length(BetaTrue )  != length(MIPVec)) {
    AFilePrint("AUC: Cant calculate because BetaTrue and MIPVec are different lengths!");
    return(-1);
  }
  RT <- length(MIPVec[BetaTrue == 0.0 & MIPVec >= .5]);
  return(RT);
} 

################################################################################
##  Type 2 Error
##
##   How many true effects do we ignore as False?
Type2Function <- function(MIPVec, BetaTrue) {
  if (length(BetaTrue) < length(MIPVec) &&
    sum(abs(round(BetaTrue) - BetaTrue)) < .0001) {
    RT <- length(BetaTrue[MIPVec[BetaTrue] < .5]);
    return(RT);
  } else if (length(BetaTrue )  != length(MIPVec)) {
    AFilePrint("AUC: Cant calculate because BetaTrue and MIPVec are different lengths!");
    return(-1);
  }
  RT <- length(MIPVec[BetaTrue != 0.0 & MIPVec < .5]);
  return(RT);
} 

############################################################################
##  Hellinger Distance
##
##    This is a distance that can be taken between absolutes and probabilities
## (there is no such thing as a Kullbeck Leiber distance from the truth.)
##
##  We report
##  sum_i  (1-sqrt(p_i))^(Beta_i != 0) + (1-sqrt(1-p_i))^(Beta_i == 0)
##
##  Note that if we report perfect 0's and 1's and get them correctly,
##   we will get a Hellinger distance of 0.0.
##  If we report incorrect 1's and 0's, getting them all false,
##   we get a maximum Hellinger distance of p.
##  If we report 50% for all probabilities, we score a perfect p*.293.
##  Reporting .75 and getting an effect correct generates a distance of .134
##  Reporting .75 and getting an effect false generates a distance of .5.
##  This means if .75 is a real posterior score our expected socre is .225
## 
MarginalHellingerFunction <- function(MIPVec, BetaTrue) {
  if (length(BetaTrue) < length(MIPVec) &&
    sum(abs(round(BetaTrue) - BetaTrue)) < .0001) {
    FFit <- (sum(1-sqrt(MIPVec[BetaTrue])) + 
      sum(1-sqrt(1-MIPVec[!(1:length(BetaTrue)  %in% BetaTrue)])))/length(MIPVec);
    return(FFit);
  } else if (length(BetaTrue )  != length(MIPVec)) {
    AFilePrint("AUC: Cant calculate because BetaTrue and MIPVec are different lengths!");
    return(-1);
  }
  BetaTrue <- (1:length(BetaTrue))[BetaTrue != 0.0];
  FFit <- (sum(1-sqrt(MIPVec[BetaTrue])) + 
    sum(1-sqrt(1-MIPVec[!(1:length(BetaTrue)  %in% BetaTrue)])))/length(MIPVec);
  return(FFit);

}  
