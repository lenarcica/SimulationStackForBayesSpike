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