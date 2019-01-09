###############################################################################
##  TwoSimR5CompellingElements
##
##   (c) 2009-2014 Alan Lenarcic
##    The following attempts to Record and describe the names of important elements of TwoSimR5
##
##    Does not appear to do much work
##

## OriginalDir:  Overall Directory where sims will be saved to
##
##  On X230: this directory is "c:/Stat//TwoSimR5/Saves/";
##
## ADDDIRCONT: Should be set in an input table
##  One Default is "PaperOfficialJuly082012"
##
## MethodDir: attached to OriginalDir//ADDDIRCONT//CR_CovRho_//MethodDir
##   Where files are saved for one partjob in a directory.
##   Where CovRho is correlation of design matrix.
##
## AllDir: Sets the globaldirectory to save a simulation: A Test:
##  AllDir <- "c:/Stat//TwoSimR5/Saves///DefaultADDDIRCONT//CR0//LSM1LSF1piAP0sPD0"
##
## CurrentContainDir: Adds another sequence describing compelling choices of 
##  Betasimulation, N, k, p, values.
##
## GenerateBetaVec2CorrelationXmatrix1PosCorrp2NLen25kLen20kALen6LSFnp8SQNp5Reps800
##



if (FALSE) {
  OnNameFunctions=NameFunctions;   
  TargetLength=TargetTotalSimulations; NotBS=TRUE;
}
# Table CompareSaverGenerate


## LSTestPrintPictureTable.R
##    -- File which runs scripts
##

##  ListOfParams  essential list of n,k,p,sigma to use for a simulation.


##  FileCons
##     FileCons is a list of length # of analyses in a group
##   For every Algorithm, a table is generated of completed simulations
##   FileCons is a list of these tables
##
##  FileLenS (the length of each table in FileCons)
##

###########################################################
## MySummaryPerformanceTable
##
##  A table recording finished analyses.


## MaxDuplicateTrials
##
##  How many machines could be working on the same problem?
##  2 is a good default
##
##
## InitiateDuplicate
##  Initiates a Duplicate if previous attempt hasn't returned in this time
##  24 * 60 * 60 /2;
##   Half a day
##
## CallFailure
##  Designates a Failure if most recent previous attempt has gone over this time
##  24 * 60 * 60 * 2
##     Two days listed as a failure

##MySummaryListPerformance:
##  List Of number of simulations
##  Each Simulation is a list for each function

z <- FALSE;
try(rm(z));