###############################################################################
## RunMakeScriptFiles.r
##
##   This R process is run by a file like "RunScriptPracticeMake.lsf"
## It is run to generate 200-1000 simulations from a parameter settings
##  provided by a settings file like "GroupLongTestn100p1000.r" or "LongTestn100p1000.r"
## It creates directories for all of the statistical estimators that will be trialed by
## this estimator. 
##  Afterwards RunSeekScriptFile will run "RunAProcessExperiment()" against the estimators*simulation combinations.
##

## R --vanilla --args GroupLongTestn100p1000.r 1 1 1 1 1 < RunMakeScriptFiles.r
## R --vanilla --args GroupLongTestn100p1000.r 1 1 1 1 1 
Oldwd <- getwd();
library(TwoSimR5);    library(AlanDirectories);

TwoSimR5:::ScriptGetArgs();
CheckUniqueProcessIdentifier();    ## In 2LassoSaveRoutines.$


verbose = 1;
Practice = 1;
if (length(list.files("c:/Stat")) >= 1) {
  argeese <- paste(DropHome, "//TwoLassoPackage//Groupn100p1000/GroupLongTestn100p1000.r", sep="");
}
setwd(Oldwd);
source(paste(as.character(argeese), sep=""));


TwoSimR5:::MakeDirSaves();
NowDelete=TRUE;
AssessFirstTime(ALargeContainDir, ASmallContainDir, TargetTotal = TargetTotal, TotalThralls=TotalThralls)
