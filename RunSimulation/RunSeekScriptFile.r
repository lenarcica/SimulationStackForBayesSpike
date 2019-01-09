###############################################################################
##  RunSeekScriptFile
## 
##  This function runs one thread with "RunAProcessExperiment()" as the major
##   function that runs a statistical estimator against a simulation.
##
## 
## R --vanilla --args GroupLongTestn100p1000.r 1 1 1 1 1 < RunSeekScriptFile.r
## R --vanilla --args PracticeTestn100p100.r 1 1 1 1 1
library(TwoSimR5);
TwoSimR5:::ScriptGetArgs();
CheckUniqueProcessIdentifier();    ## In 2LassoSaveRoutines.$

##try(library(TwoLassoCpp, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);

TwoSimR5:::ScriptGetArgs();

verbose = 1;
Practice = 1;
source(paste(as.character(argeese), sep=""));

##MyID <- as.numeric(Myargs[1]);
TwoSimR5:::MakeDirSaves();

ABStart <- proc.time();

ALTFiles <- unlist(list.files(paste(ASmallContainDir, "//Submissions", sep="")));
TotalSlaves <- length(ALTFiles);
WantFile <- paste("Request", ZZMyID, ".RData", sep="");
if (WantFile %in% ALTFiles) {
  setwd(ASmallContainDir);
  try(load(paste("Submissions//", WantFile, sep="")));
  setwd(OriginalOldwd);
} else {
  print(paste("Oh No No RData file for MyID = ", MyID, sep=""));
  q();
}

ABB <- min((1:NROW(SubMitList))[SubMitList[,3] %in% c("0", 0, " 0")])

SMSName <- SubMitList[ABB,1];
FunctionName  <- SubMitList[ABB,2]; 
ABack <- TwoSimR5:::RunAProcessExperiment(SMSName=SubMitList[ABB,1], 
  FunctionName =SubMitList[ABB,2], verbose = verbose, DoCI = DefaultDoCI)
ABFinish <- proc.time();
if (ABFinish[3] - ABStart[3] <= 5) {
  Sys.sleep(5- ABFinish[3]+ABStart[3])
}

SubMitList[ABB,3] <- as.character(ABack);
setwd(ASmallContainDir);
try(save(SubMitList=SubMitList, ASmallContainDir=ASmallContainDir, 
       ALargeContainDir=ALargeContainDir, 
       CurrentLargeContainDir=CurrentLargeContainDir, 
       CurrentSmallContainDir=CurrentSmallContainDir, file=paste("Submissions//", WantFile, sep="")));
q();
