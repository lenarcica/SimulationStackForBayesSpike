################################################################################
## 0011AddMyVToRTable.r
##   Manipulate an RData to save results of a simulation for TwoSimR5 package
##
##     (c) 2009-2019 Alan Lenarcic
##     The code was written to support simulations for Lenarcic and Valdar methods
##      paper.
##
##   Here "MyV" is the result of a simulation, usually L1, L2, Type1, Type2 Error
##    and runtime estimates.  Since the LSF server UNC Killdevil is the
##    designed target for deployment, and since any given R package or estimator
##    could fail or run out of memory or time, a try-catch framework is
##    necessary to save results of simulations in a controlled matter to disk.
##
##

## LICENSE INFO: R CODE
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/
#
#  Note, TwoSimR5 code is predominantly built around running existing 
#  selector impleentations which exist as R packages, most of which have
#  a GNU license.
AddMyVToRDataTable <- function(OnFunction = OnFunction,
       RDataTableOutName, MyV,
       PiD, AlreadyLocked=FALSE, verbose=0) {
  AFilePrint(paste("   --  AddMyVToRDataTable putting on OnFunction = ", OnFunction, sep=""));
  eval(parse(text=GetG0Text("CurrentSmallContainDir")));
  MDir <- paste(CurrentSmallContainDir, "//", OnFunction, sep="");
  dir.create(paste(CurrentSmallContainDir, "//", OnFunction, sep=""), showWarnings=FALSE, recursive=TRUE);
  quoteMore = paste(" FinishedSimulation: 0011AddMyVToRTable.r: ", IFunction, " : ", OnFunction, sep=""); flush.console();
  Oldwd <- getwd();
   if (AlreadyLocked==FALSE) {
     TryLockText <- "
       DidLock = FALSE; 
       try(setwd(CurrentSmallContainDir)); 
       while(LockMeIn(verbose=verbose, quoteMore=quoteMore, 
         LFDir = OnFunction, MoreLock = \"\", LockIn=\"No\", SumLock=\"No\")==FALSE) {
          DidLock <- FALSE;  Sys.sleep(4);
       } 
       try(setwd(Oldwd));
       DidLock <- TRUE;
     ";
     try(eval(parse(text=TryLockText)));  
     if (DidLock == FALSE) {
        AFilePrint("***************************************************************");
        AFilePrint("** AddMyVToRDataTable: we tried to perform Lock but failed.");
        AFilePrint("** -- Was anything given to the quoteMore?  ");
        AFilePrint("**  -- Did we have to perform a lock?  ");
        AFilePrint(paste("** -- our OnFunction was ", OnFunction, sep=""));
        AFilePrint(paste("** -- our CurrentSmallContainDir = ", CurrentSmallContainDir, sep=""));
        try(eval(parse(text=SetGText("OnFunction"))));
     }   
   }

   ListFiles <- unlist(list.files(MDir));
   RDD <- NULL;  MyColNames = NULL;  MyPiDs = NULL;  MyUniqueProcesses=NULL;
   TryToOpenText <- "
     SuccessTry <- FALSE;
     if (length(ListFiles) >= 1 && any(ListFiles == RDataTableOutName)) {
       try(setwd(CurrentSmallContainDir));
       try(load(paste(OnFunction, \"//\", RDataTableOutName, sep=\"\")));
       try(setwd(Oldwd));
     }
     SuccessTry <- TRUE;
   ";
   try(eval(parse(text=TryToOpenText)));
   if (SuccessTry == FALSE) {
     AFilePrint(paste("AddMyVToRDataTable: Something went weird in try to Load MDir = ", MDir, sep=""));
   }
   
   eval(parse(text=GetG0Text("UniqueProcessIdentifier")));
   if (is.null(RDD)) {
      AFilePrint(paste("   --  AddMyVToRDataTable putting on OnFunction = ", OnFunction, ": Null RDD ", sep=""));
      eval(parse(text=GetG0Text("MyPermaColNames", "globalenv()", S=1)));
      MyPiDs <- PiD;
      RDD <- matrix(MyV, 1, NROW(MyV));
      try(colnames(RDD) <- MyPermaColNames);
      MyColNames <- MyPermaColNames;
      try(rownames(RDD) <- PiD);
      MyUniqueProcesses <- UniqueProcessIdentifier;
      try(save(RDD=RDD, MyPiDs=MyPiDs, MyUniqueProcesses=MyUniqueProcesses, MyColNames=MyColNames,
        RDataTableOutName = RDataTableOutName, MDir=MDir, 
        file=paste(MDir, "//", RDataTableOutName, sep="")));
      if (AlreadyLocked==FALSE) {
        try(setwd(CurrentSmallContainDir));
        UnLockMe(verbose=verbose, 
          quoteMore=paste(quoteMore, " - AddMyVToRDataTable - ", OnFunction, sep=""),
          LFDir = OnFunction);
        try(setwd(Oldwd));
      } 
      return(1);
   }
   NewPiDs <- c(MyPiDs, PiD);
   NewUniqueProcesses <- c(MyUniqueProcesses, UniqueProcessIdentifier);
   if (is.matrix(RDD)) {
      if (NCOL(RDD) !=  length(MyColNames)) {
        AErrorPrint(paste("Hey AddMyVToRDataTable: NCOL(RDD) = ", NCOL(RDD), " but length(MyColNames) = ", 
          length(MyColNames), sep=""));
        tryCatch("AddMyVToRDataTable: NCOLs of RDD is Bad!");
      }
      if (NCOL(RDD) !=  length(MyV)) {
        AErrorPrint(paste("Hey AddMyVToRDataTable: NCOL(RDD) = ", NCOL(RDD), " but length(MyV) = ", 
          length(MyV), sep=""));
        tryCatch("AddMyVToRDataTable: NCOLs of RDD is Bad!");
      }    
      NewRDD <- rbind(RDD, MyV);
      try(rownames(NewRDD) <- NewPiDs);
      try(colnames(NewRDD) <- MyColNames);
   } else if (is.vector(RDD)) {
      if (length(RDD) !=  length(MyColNames)) {
        AErrorPrint(paste("Hey AddMyVToRDataTable: length(RDD) = ", length(RDD), " but length(MyColNames) = ", 
          length(MyColNames), sep=""));
        tryCatch("AddMyVToRDataTable: NCOLs of RDD vector is Bad!");
      }
      if (length(RDD) !=  length(MyV)) {
        AErrorPrint(paste("Hey AddMyVToRDataTable: length(RDD) = ", length(RDD), " but length(MyV) = ", 
          length(MyV), sep=""));
        tryCatch("AddMyVToRDataTable: Length of RDD  vector is Bad!");
      }    
      NewRDD <- rbind(RDD, MyV);
      try(rownames(NewRDD) <- NewPiDs);
      try(colnames(NewRDD) <- MyColNames);
   }
   ASort <- sort(NewPiDs, index=TRUE)$ix;
   MyPiDs <- NewPiDs[ASort];
   MyUniqueProcesses <- NewUniqueProcesses[ASort];
   RDD <- NewRDD[ASort,];
   try(save(RDD=RDD, MyPiDs=MyPiDs, MyUniqueProcesses=MyUniqueProcesses, MyColNames=MyColNames,
     RDataTableOutName = RDataTableOutName, MDir=MDir, 
     file=paste(MDir, "//", RDataTableOutName, sep="")));
   if (AlreadyLocked==FALSE) {
      try(setwd(CurrentSmallContainDir));
      UnLockMe(verbose=verbose, 
        quoteMore=paste(quoteMore, " - AddMyVToRDataTable - ", OnFunction, sep=""),
        LFDir = OnFunction);
      try(setwd(Oldwd));
   }     
   return(1);
}