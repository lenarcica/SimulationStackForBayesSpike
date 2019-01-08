## CINewLoss.R

## bxglm cjbyo izcky cpalh
##
##Int (A-X)^2 exp((-X-mu)/sigma^2) dx
##Int (A-mu+mu-X)^2 exp((-X-mu)/sigma^2) dx
##Int (A-mu+mu-X)^2 exp((-X-mu)^2/2sigma^2) dx/sqrt(2pi sigma^2)
## Int (A-mu)^2*exp((-X-mu)/sigma^2 )dx/sqrt(2pi sigma^2) = (A-mu)^2 * Prob
##  -2 (A-mu)  Int (X-mu) exp( -(X-mu)^2/2sigma^2)dx/sqrt(2pi sigma^2)
## Int (X-mu)^2 exp(-(X-mu)^2/2sigma^2)dx/sqrt(2pi sigma^2);
## Int_B^A sqrt(2 sigma^2 v) exp(-v) dv / sqrt(2pi);
## sqrt(2 sigma^2) Gamma(3/2) Int_b^a v^{3/2-1} exp(-v)dv
if (FALSE) {
  DownInt <- c(0,0,1,0,-1,1);
  UpInt <- c(.5,0,1.5,.5,-.5,1);
  Truth <- c(.25, .1, 1.25, .35, -.25,1)
  Prob <- .95;
  
  Oni <-1 ;
  Up<-UpInt[Oni]; Down <- DownInt[Oni];
  Alpha <- 1-Prob;
  ABD <- qnorm(Prob + Alpha/2) - qnorm(Alpha/2);
  OurSigma <- (Up - Down) / ABD;
  Mu <- (Down+Up)/2;
  XXX <- Down+0:50000/50000*(Up-Down);
  sum((XXX-Mu)^2 * exp(-(XXX-Mu)^2/(2*OurSigma^2)) /(sqrt(2*pi*OurSigma^2))) * (XXX[2]-XXX[1]);
  AA = Truth[Oni];
  sum((XXX-AA)^2 * exp(-(XXX-Mu)^2/(2*OurSigma^2)) /(sqrt(2*pi*OurSigma^2))) * (XXX[2]-XXX[1]);
  
  bbb = (Up-Mu)^2 / (2*OurSigma^2);
  aaa = (Down-Mu)^2 / (2*OurSigma^2);
  VVV <- 0 + 0:10000/10000 * bbb;
  sum( (2*OurSigma^2)*sqrt(VVV) * exp(-VVV) / sqrt(pi)) * (VVV[2]-VVV[1])
  A3 <- 1/sqrt(pi) * OurSigma^2 * gamma(3/2) *
    (sign(Up-Mu)*pgamma(bbb, 3/2) - sign(Down-Mu)*pgamma(aaa,3/2));
  
  XXX <- 0:10000/10000*2;
  sum( XXX^(3/2-1)* exp(-XXX) ) * (XXX[2]-XXX[1]);
  pgamma(max(XXX),3/2)  * gamma(3/2);
  
  A = -1;
  XXX <- Down+0:50000/50000*(Up-Down);
  (A-Mu)^2 * sum( exp(-(XXX-Mu)^2/(2*OurSigma^2)) / sqrt(2* pi *OurSigma^2)) * (XXX[2]-XXX[1]);

  XXX <- Down+0:50000/50000*(Up-Down);
  -2 * sum( (XXX-Mu) * exp(-(XXX-Mu)^2/(2*OurSigma^2))) * (XXX[2] - XXX[1]);
  A2 <- -2  * OurSigma^2 *  (exp(-aaa) + exp(-bbb));   
}
CICalcConservativeSignLoss <- function(DownInt, UpInt, Truth) {
  return(length(DownInt[
    (DownInt > 0 & UpInt > 0 & Truth < 0) |
    (DownInt < 0 & UpInt < 0 & Truth > 0)]));

}
CICalcLiberalSignLoss <- function(DownInt, UpInt, Truth) {
  return(length(DownInt[
    (DownInt > 0 & UpInt > 0 & Truth <= 0) |
    (DownInt < 0 & UpInt < 0 & Truth >= 0)]));

}

CINewLoss <- function(DownInt, UpInt, Truth, Prob) {
  Alpha <- 1-Prob;
  ABD <- qnorm(Prob + Alpha/2) - qnorm(Alpha/2);
  OurSigma <- (UpInt - DownInt) / ABD;
  Mus <- (UpInt+DownInt)/2;
  A1 <- Prob * (Truth - Mus)^2;
    aaa <- (DownInt-Mus)^2 / (2*OurSigma^2);
    bbb <- (UpInt-Mus)^2 / (2*OurSigma^2);
  A2 <- -2 * ( Truth- Mus) * OurSigma^2 *  (exp(-aaa) - exp(-bbb));
  A3 <- 1/sqrt(pi) * OurSigma^2 * gamma(3/2) *
    (sign(UpInt-Mus)*pgamma(bbb, 3/2) - sign(DownInt-Mus)*pgamma(aaa,3/2));
  PD <- A1+A2+A3;
  if (length(DownInt) == 1) {
    if (DownInt[1] == UpInt[1]) {
       PD <- Prob * (Truth[1] - DownInt[1])^2;
    }
  } else {
    if (length(DownInt[DownInt==UpInt]) >= 1) {
      PD[DownInt == UpInt] <- Prob * (Truth[DownInt == UpInt]-DownInt[DownInt == UpInt])^2;
    }
  }
  PD <- PD / Prob^2;
  return(PD);    
}
CIUnifLoss <- function(DownInt, UpInt, Truth, Prob) {
  ## int_a^b (T-x)^2 dx = -(T-b)^2/3+(T-a)^3 /3
  RD <- Prob*( (UpInt-Truth)^3/3-(DownInt-Truth)^3/3 ) / ((UpInt-DownInt));
  if (length(DownInt)==1) {
    if (DownInt==UpInt[1]) {
      RD <- (DownInt - Truth)^2 * Prob;
      return(RD);
    }
    RD = RD / Prob^2;
    return(RD);
  }
  if (length(DownInt[DownInt==UpInt]) >= 1) {
    RD[DownInt==UpInt] <- Prob*(DownInt[DownInt==UpInt] - Truth[DownInt==UpInt])^2;
  }
  RD <- RD / Prob^2;
  return(RD);
  ##XX <-DownInt + (0:10000)/10000*(UpInt-DownInt);
  ##sum((XX-Truth)^2 * (XX[2]-XX[1]));
}