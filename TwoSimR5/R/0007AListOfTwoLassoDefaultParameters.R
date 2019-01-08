
   
SetEfficientSimulatorDefaults <- function() {   
   eval(parse(text=GetG0Text("TWOSIMENVIRONMENT", S=1)));
   if (is.numeric(TWOSIMENVIRONMENT) && TWOSIMENVIRONMENT == 0) {
     TWOSIMENVIRONMENT = environment();
     eval(parse(text=SetGText("TWOSIMENVIRONMENT", envir="globalenv()", S=1)));
   }
   
   MySetter <- function(AText, AVal) {
     eval(parse(text=GetG0Text("TWOSIMENVIRONMENT", S=1, envir="globalenv()")));
     eval(parse(text=GetG0Text(AText, S=1, envir="globalenv()")));
     eval(parse(text=paste(AText, " = ", AVal, sep="")));
     eval(parse(text=SetGText(AText, S=1, envir= "globalenv()")));
     eval(parse(text=SetGText(AText, S=1, envir= "TWOSIMENVIRONMENT")));     
   }
   MySetter("OnSimType", "\"Default\"");
   MySetter("ALT3", 10);
   MySetter("ALT2", 10*11/2);
   MySetter("LARSSEEKMETHOD", 1);
   MySetter("STANDARDFLAG", 0);
   MySetter("tNoiseDF", 0);
   MySetter("DefaultNoiseDF", 0);
   
   MySetter("piSV",0);  MySetter("sigSV", 0);
   #######################################################################
   ###   SimKActiveLen;  SimSigmaNoise;
   ###
   ### To Give False Input to simulations, we give
   ###
   eval(parse(text=GetG0Text("SimKActiveLen", S=1, envir="globalenv()")));
   SimKActiveLen <- function (kActiveLen, NLen, kLen)  {
            min(c(NLen-1, kLen -1, max(round(kActiveLen + rnorm(1, 0, sqrt(kLen)*piSV)), 1)));
   }
   eval(parse(text=SetGText("SimKActiveLen", S=1, envir="globalenv()")));
   eval(parse(text=SetGText("SimKActiveLen", S=1, envir="TWOSIMENVIRONMENT"))); 
   
   eval(parse(text=GetG0Text("OnSimType", S=1)));
   if (OnSimType == "Group" || OnSimType == "GroupRobit" || OnSimType == "GroupLogit") {
     MySetter("pGroups", 10);
     MySetter("kGroups", 2);
     MySetter("GroupSize",5);
     MySetter("tNoiseDF", -1);
     MySetter("MeanCenterXX", TRUE)
     MySetter("Beta0", .15);
     MySetter("GroupsSumToZero", TRUE);
   }
   MySetter("NLen", 20);
   MySetter("kLen", 10);
   MySetter("LARSSeekFlag", 3);
   MySetter("kActiveLen", 4);
   MySetter("NumReps", 500);
   MySetter("SigmaNoise", .5);
   MySetter("musize", 1); 
   MySetter("SigmaEta", 2);
   MySetter("SigmaBar", -999); 
   
   eval(parse(text=GetG0Text("DefaultConfidenceIntervals", S=1, envir="globalenv()")));
   DefaultConfidenceIntervals = c(.5, .25, .75, .05,.95, .025, .975, .01, .99, .005, .995);  
   eval(parse(text=SetGText("DefaultConfidenceIntervals", S=1, envir="globalenv()")));
   eval(parse(text=SetGText("DefaultConfidenceIntervals", S=1, envir="TWOSIMENVIRONMENT")));
   eval(parse(text=GetG0Text("DefaultConfidenceIntervals")));
   GetConfidenceIntervals = DefaultConfidenceIntervals;
   
   TestConfidenceIntervals <- function(ToTest, AStringHelp="") {
     if (ToTest[1] == .5) {
       AK <- 1;
     } else {AK = 0;}
     AtAK  = AK;
     if (length(ToTest) <= 1) {
       print(paste(AStringHelp, "TestConfidenceIntervals;, ",
         "Bad User Supplied Confidence Intervals (GetConfidenceIntervals) have length ", 
         length(ToTest), ".", sep="")); return(-1);
     }
     NtoTest <- (length(ToTest)-AK) /2;
     if (length(ToTest) != NtoTest * 2 + AK) {
       print(paste(AStringHelp, "TestConfidenceIntervals;, ",
         "Bad User Supplied Confidence Intervals (GetConfidenceIntervals) have length ", 
         length(ToTest), ".", sep=""));   return(-1);
     }  
     for (ii in 1:length(NtoTest)) {
       if(ToTest[AtAK+1] <= 0.0 || ToTest[AtAK+1] >= 1.0) {
          print(paste(AStringHelp, "TestConfidenceIntervals;, ",
         "Bad User Supplied Confidence Intervals: (", 
         paste(ToTest, collapse=", "),")", ".", sep=""));   return(-1);
       }
       if(ToTest[AtAK+2] <= 0.0 || ToTest[AtAK+2] >= 1.0) {
          print(paste(AStringHelp, "TestConfidenceIntervals;, ",
         "Bad User Supplied Confidence Intervals: (", 
         paste(ToTest, collapse=", "),")", ".", sep=""));   return(-1);
       }
        if(ToTest[AtAK+1] >=  ToTest[AtAK+2]) {
          print(paste(AStringHelp, "TestConfidenceIntervals;, ",
         "Bad User Supplied Confidence Intervals: (", 
         paste(ToTest, collapse=", "),")", ".", sep=""));   return(-1);
       }
       AtAK <- AtAK+2;
     }  
     return(1);
   }
   eval(parse(text=SetGText("TestConfidenceIntervals", S=1, envir="globalenv()")));
   eval(parse(text=SetGText("TestConfidenceIntervals", S=1, envir="TWOSIMENVIRONMENT")));
   TestConfidenceIntervals(GetConfidenceIntervals);  
   eval(parse(text=SetGText("GetConfidenceIntervals", S=1, envir="globalenv()")));
   eval(parse(text=SetGText("GetConfidenceIntervals", S=1, envir="TWOSIMENVIRONMENT")));
   MySetter("DefaultDoCI", TRUE);
   
   eval(parse(text=GetG0Text("FileName", S=1, envir="globalenv()")));
   FileName <- "c://Stat//Playfile.csv";
   eval(parse(text=SetGText("FileName", S=1, envir="TWOSIMENVIRONMENT"))); 
   eval(parse(text=SetGText("FileName", S=1, envir="globalenv()")));
   
     
   MySetter("PrintOutFlags", 1);

}