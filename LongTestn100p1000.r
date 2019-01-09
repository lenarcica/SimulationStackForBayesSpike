################################################################################
## LongTestn100p1000.r
##    (c) 2009-2019 Alan Lenarcic
##
##  An example input file to be "source(LongTestn100p1000.r)" on loading libary(TwoSimR5)
##  This just gives parameter settings for a n=100, p=1000 model.
##
NRList=c(1000, 1000, 1000)
NNList=c(100, 100, 100);
GimmeList=c(6); 
SigmaList = c(.5);

NNList = c(25,25,25,25, 25, 15, 25, 35, 45);
NRList = c(10,20,30,40, 50, 50, 50, 50, 100, 100, 100, 100);
GimmeList = c(6,6,6,6,6,6,6,6,6,6,6,6);
SigmaList = c(.5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5);


SkipLoad= 0;
Interval = .95;

MyPrintRoutines = "PrintRoutinesNormal.R"

##NNList = c(25,25,25,25, 25, 15, 25, 35, 45);
##NRList = c(20,30,40,40, 50, 50, 50, 50, 100, 100, 100, 100);
##GimmeList = c(6,6,6,6,6,6,6,6,6,6,6,6);
##SigmaList = c(.5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5);

ADDDIRCONT = "Test2013April";
n <- 100;  p <- 1000;  sigma <- .5;  k <- 6;

## ADDIRCONT = "Square2";
##NNList = c(25,25,25,15,25,35, 25,50,50);
##NRList = c(10,20,30, 50,50,50, 100,100,100);
##GimmeList=c(6,6,6,6,6,6,6,6,6);
##SigmaList = c(.5, .5, .5, .5, .5, .5, .5, .5, .75);

LARSSEEKMETHOD = 1;
LARSSeekFlag = -.8 ;
   piAPriorData = 1;
   piAPriorData = .1;
   sigmaPriorData = .1;
Interval = .95;


MaxRunProcTime <- .66 * 60 * 60 * 24

##  ADDDIRCONT = "SquareN";
##   ADDDIRCONT = "Square1";
##NNList = c(25,25,25,15,25,25, 25,50,50);
##NRList = c(20,30,40, 50,50,100, 100,100,100);
##GimmeList=c(6,6,6,6,6,6,6,6,6);
##SigmaList = c(.5, .5, .5, .5, .5, .75, .25, .5, .75);

##NumReps = 100000;
NumReps = 800;
CovRho = .2;
DirectorForImage = "NineSquareJMLR";
DirectorForImage = paste(DirectorForImage, "NumReps", NumReps, "CovRho", tSeq(CovRho), sep="");
  DefaultGenerateBetaVec  = GenerateBetaVec2;  BM = 5;

  DefaultCorrelationXmatrix <- function(kLen) {
      CorrelationXmatrix1(kLen, CovRho);
  }

if(!exists("Myargs")) {
  eval(parse(text=GetG0Text("Myargs", S=1)));
  Myargs = commandArgs(TRUE); 
  if (length(Myargs) < 2) {
    Myargs = c(1,1); 
  }
  eval(parse(text=SetGText("Myargs", S=1)));
}
  
##NRList = c(20,20,20);
##NNList = c(30,40,50);
##GimmeList=c(55,55,55);

#argMe =  as.numeric(args[2]);
##print(paste(" argMe is ", argMe, sep=""));
##NRii = ceiling(argMe  / (length(NNList) * length(GimmeList) * length(SigmaList) ) );
##if (NRii > 0 && NRii <= length(NRList)) {
##   NR = NRList[NRii];
##   p = NR;
##}  else {
##   quit();
##}
##ArgResid =   argMe -length(NNList) * length(GimmeList)  *length(SigmaList) * (NRii -1)
##NNii = ceiling(  ( ArgResid  ) /
##              (length(GimmeList) *length(SigmaList) ) );
##if (NNii > 0 && NNii <= length(NNList)) {
##    NN = NNList[NNii];
##    n = NN;
##}  else {
##    quit();
##}
##
##ArgResid2 = ArgResid - length(GimmeList) *length(SigmaList) * (NNii-1);
##GLii = ceiling( (ArgResid2) / length(SigmaList) );
##if (GLii > 0 && GLii <= length(GimmeList)) {
##     Gimme = GimmeList[GLii];
##     k = Gimme;
##}    else {
##    quit();
##}
##
##ArgResid3 = ArgResid2 - length(SigmaList) * (GLii-1);
##SGii = ArgResid3;
##if (SGii > 0 && SGii <= length(SigmaList)) {
##    sigma = SigmaList[SGii];
##}  else {
##   quit();
##}

   #######################################################################
   ###   SimKActiveLen;  SimSigmaNoise;
   ###
   ### To Give False Input to simulations, we give
   ###
   eval(parse(text=GetG0Text("piSV", S=1)));
   piSV = 3;
   SimKActiveLen <- function (kActiveLen, NLen, kLen)  {
           OneNum = round( kActiveLen + sqrt(kLen) * piSV * rnorm(1,0,1) )
           if (OneNum >= kLen) { OneNum = kLen -1; }
           if (OneNum >= NLen) { OneNum = NLen -1; }
           if (OneNum <= 0) { OneNum = 1; }
           return(OneNum);
   }
   eval(parse(text=SetGText("piSV", S=1)));
   eval(parse(text=GetG0Text("sigSV", S=1)));
   sigSV=3;
    SimSigmaNoise <- function (SigmaNoise) {
       if (sigSV > 0) {
         ##print("SuckSuckSuckSuck") ;
         return(sqrt(SigmaNoise^2 * rchisq(1, sigSV) /sigSV));
       } else{
         return(SigmaNoise);
       }
    }
    eval(parse(text=SetGText("sigSV", S=1)));   
   
   
   #######################################################################
   ###   SimKActiveLen;  SimSigmaNoise;
   ###
   ### To Give False Input to simulations, we give
   ###
   ##SimKActiveLen <- function (kActiveLen)  {
   ##            kActiveLen + rnorm(1, 0, sqrt(kActiveLen)/2);
   ##}
   ##SimSigmaNoise <- function (SigmaNoise) {
   ##      sqrt(SigmaNoise * rchisq(10) / 10);
   ##}
   eval(parse(text=GetG0Text("DoCI", S=1)));
   DoCI = TRUE;
   eval(parse(text=SetGText("DoCI", S=1)));
  



