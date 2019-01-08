###############################################################################
### Acpp.r
###
###  Alan Lenarcic 07/27/2014
###
###
###  Implementation of functions necessary to run Acpp cloaking of a gc protected list.

### LockG, SetG, checking what happens:
if (!exists("SetAcppText")) {
SetAcppText <- function(AText) {
  if (!is.character(AText)) {
    print("SetAcppText: deep error, AText is not character");
    flush.console();
  }
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ACPPNAMESPACE), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ACPPNAMESPACE ) 
    }
    assign( \"", AText, "\", ", AText, ", ACPPNAMESPACE);
    ", sep="");
}
}
if (!exists("SetGText")) {
SetGText <- function(AText, envir = "globalenv()", S=0) {
  if (!is.character(AText)) {
    print("SetGText: big error, AText is not char"); flush.console();
  }
  if (S==0) {
  return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, " ); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    ", sep=""));
  } else {
   return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, " ); 
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    ", sep="")); 
  }
}
}

if (!exists("LockACPPText")) {
LockACPPText <- function(AText) {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ACPPNAMESPACE), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", RSIMNAMESPACE ) 
    }
    try(unlockBinding( \"", AText, "\", ACPPNAMESPACE ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ACPPNAMESPACE);
    lockBinding( \"", AText, "\", ACPPNAMESPACE);
    ", sep="");
}
}



if (!exists("LockGText")){
LockGText <- function(AText, envir="globalenv()", S=1)  {
  if (S==0) {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=FALSE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=FALSE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep=""));
  } else {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep=""));   
  }
}
}


if (!exists("GetG0Text")) {
GetG0Text <- function(AText, envir = "globalenv()", S=0) {
  if (S== 0) {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=FALSE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));
  } else {
    return(paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\",  ", envir, " ), silent=TRUE);
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE);
    }
    try(unlockBinding( \"", AText, "\", ", envir, " ), silent=TRUE); 
    if (exists(\"", AText, "\", ", envir, " )) {
    ", AText, " <- get(\"", AText, "\", ", envir, " );
    } else {
      ", AText, " = 0;
    }
    ", sep=""));   
  }
}
}

g=1;
eval(parse(text=SetGText("g", envir="globalenv()", S=1)));
  
## Hey a dude has got to somehow make directory loading automatic for his computer!
.onLoadAcpp <- function(libname, pkgname) {
  ## load the module and store it in our namespace
  ##ACPPNAMESPACE <- environment()
  ## The OnLoad creates an ACPPNAMESPACE 
   AP = paste("TryO = NULL; try( TryO <- bindingIsActive(\"", "ACPPNAMESPACE", "\", ", "environment()", "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", "  ) 
    }
    try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", "  ), silent=TRUE); 
    try(unlockBinding( \"", "ACPPSAVESPACE", "\", ", "environment()", "  ), silent=TRUE); 
    try(ACPPNAMESPACE <- environment(), silent=TRUE);
    try(ACPPNAMESPACE <- .Call(\"GetMyPackageNamespace\"));
    try(ACPPSAVESPACE <- new.env());
    assign( \"", "ACPPNAMESPACE", "\", ", "ACPPNAMESPACE", ", ", "environment()", " );
    assign( \"", "ACPPSAVESPACE", "\", ", "ACPPSAVESPACE", ", ", "environment()", " );
    try(lockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPSAVESPACE", "\", ", "environment()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPNAMESPACE", "\", ", "ACPPNAMESPACE", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPSAVESPACE", "\", ", "ACPPNAMESPACE", " ), silent=TRUE);
    ", sep="");
  try(eval(parse(text=AP)), silent=TRUE);
  
 if (!exists("ACPPNAMESPACE", globalenv())) {
   AP = paste("TryO = NULL; try( TryO <- bindingIsActive(\"", "ACPPNAMESPACE", "\", ", "environment()", "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", "  ), silent=TRUE) 
    }
    try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", "  ), silent=TRUE); 
    try(unlockBinding( \"", "ACPPSAVESPACE", "\", ", "environment()", "  ), silent=TRUE); 
    try(ACPPNAMESPACE <- environment(), silent=TRUE);
    try(ACPPNAMESPACE <- .Call(\"GetMyPackageNamespace\"));
    try(assign( \"", "ACPPNAMESPACE", "\", ", "RSIMNAMESPACE", ", ", "environment()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPSAVESPACE", "\", ", "ACPPNAMESPACE", " ), silent=TRUE);
    try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "globalenv()", "  ), silent=TRUE); 
    try(unlockBinding( \"", "ACPPSAVESPACE", "\", ", "globalenv()", "  ), silent=TRUE); 
    try(lockBinding( \"", "ACPPNAMESPACE", "\", ", "globalenv()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPSAVESPACE", "\", ", "globalenv()", " ), silent=TRUE);
    ", sep="");
    try(eval(parse(text=AP)), silent=TRUE);
 }
 if (!exists("LockGText")) {
 LockGText <- function(AText, envir="globalenv()", S=0)  {
  if (S == 1) {
   paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep="");  
  } else {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep="");
  }
  }
  }
  eval(parse(text=LockGText("LockGText", envir="globalenv()", S=1)));
  try(eval(parse(text=LockGText("LockGText", envir="ACPPNAMESPACE", S=1))));
  try(eval(parse(text=LockGText("ACPPNAMESPACE", envir="globalenv()", S=1))));
  try(eval(parse(text=LockGText("ACPPSAVESPACE", envir="globalenv()", S=1))));
  try(eval(parse(text=LockGText("ACPPSAVESPACE", envir="ACPPNAMESPACE", S=1))));
    
  Verbose = 0; eval(parse(text=SetGText("Verbose", S=1)));
  

  ##eval(parse(text=GetG0Text("KillVerboseChains",S=1)));
  ##KillVerboseChains = TRUE;
  ##eval(parse(text=LockGText("KillVerboseChains",S=1)));
  try(AlanDirectories:::SetSaveHomes(), silent=TRUE);
}


.onAttachAcpp <- function(libname, pkgname) {
  ## load the module and store it in our namespace
  ##ACPPNAMESPACE <- environment()
   AP = paste("TryO = NULL; try( TryO <- bindingIsActive(\"", "ACPPNAMESPACE", "\", ", "environment()", "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", "  ) 
    }
    try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", "  ), silent=TRUE); 
    try(unlockBinding( \"", "ACPPSAVESPACE", "\", ", "environment()", "  ), silent=TRUE); 
    try(ACPPNAMESPACE <- environment());    
    try(ACPPNAMESPACE <- .Call(\"GetMyPackageNamespace\"));
    if (exists(\"ACPPSAVESPACE\", ACPPNAMESPACE)) {
      ACPPSAVESPACE <- get(\"ACPPSAVESPACE\", ACPPNAMESPACE);
    } else { ACPPSAVESPACE <- NULL; }
    if (is.null(ACPPSAVESPACE) || !is.environment(ACPPSAVESPACE)) {
      ACPPSAVESPACE <- new.env();
    }
    assign( \"", "ACPPNAMESPACE", "\", ", "ACPPNAMESPACE", ", ", "environment()", " );
    assign( \"", "ACPPSAVESPACE", "\", ", "ACPPSAVESPACE", ", ", "environment()", " );
    assign( \"", "ACPPSAVESPACE", "\", ", "ACPPSAVESPACE", ", ", "ACPPNAMESPACE", " );
    try(lockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPSAVESPACE", "\", ", "ACPPNAMESPACE", " ), silent=TRUE);
    ", sep="");
  try(eval(parse(text=AP)), silent=TRUE);
  
 if (!exists("ACPPNAMESPACE", globalenv())) {
   AP = paste("TryO = NULL; try( TryO <- bindingIsActive(\"", "ACPPNAMESPACE", "\", ", "environment()", "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", "  ), silent=TRUE) 
    }
    try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", "  ), silent=TRUE); 
    try(unlockBinding( \"", "ACPPSAVESPACE", "\", ", "environment()", "  ), silent=TRUE); 
    try(ACPPNAMESPACE <- environment(), silent=TRUE);
    try(ACPPNAMESPACE <- .Call(\"GetMyPackageNamespace\"));
    try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "ACPPNAMESPACE", "  ), silent=TRUE); 
    try(unlockBinding( \"", "ACPPSAVESPACE", "\", ", "ACPPNAMESPACE", "  ), silent=TRUE); 
    try(assign( \"", "ACPPNAMESPACE", "\", ", "ACPPNAMESPACE", ", ", "environment()", " ), silent=TRUE);
    try(assign( \"", "ACPPSAVESPACE", "\", ", "ACPPSAVESPACE", ", ", "environment()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPNAMESPACE", "\", ", "environment()", " ), silent=TRUE);
    try(unlockBinding( \"", "ACPPNAMESPACE", "\", ", "globalenv()", "  ), silent=TRUE); 
    try(unlockBinding( \"", "ACPPSAVESPACE", "\", ", "globalenv()", "  ), silent=TRUE); 
    try(lockBinding( \"", "ACPPNAMESPACE", "\", ", "globalenv()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPSAVESPACE", "\", ", "globalenv()", " ), silent=TRUE);
    try(lockBinding( \"", "ACPPSAVESPACE", "\", ", "ACPPNAMESPACE", " ), silent=TRUE);
    ", sep="");
    try(eval(parse(text=AP)), silent=TRUE);
  }

if (!exists("LockGText")) {
LockGText <- function(AText, envir="globalenv()", S=0)  {
  if (S == 1) {
   paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep="");  
  } else {
  paste("TryO = NULL; try( TryO <- bindingIsActive(\"", AText, "\", ", envir, "), silent=TRUE)
    if (!is.null(TryO) && (TryO == TRUE || TryO == FALSE)) {
      unlockBinding( \"", AText, "\", ", envir, "  ) 
    }
    try(unlockBinding( \"", AText, "\", ", envir, "  ), silent=TRUE); 
    assign( \"", AText, "\", ", AText, ", ", envir, " );
    lockBinding( \"", AText, "\", ", envir, " );
    ", sep="");
  }
}
}
  eval(parse(text=LockGText("LockGText", S=1)));
  try(eval(parse(text=LockGText("LockGText", "ACPPNAMESPACE", S=1))));
  
  Verbose = 0; eval(parse(text=SetGText("Verbose", S=1)));
  

  ##eval(parse(text=GetG0Text("KillVerboseChains",S=1)));
  ##KillVerboseChains = TRUE;
  ##eval(parse(text=LockGText("KillVerboseChains",S=1)));
  try(AlanDirectories:::SetSaveHomes(), silent=TRUE);
}

if (!exists(".onAttach")) {
  .onAttach = .onAttachAcpp;
}
   




.CreateList <- function() {
  try(eval(parse(text=GetG0Text("ACPPNAMESPACE", envir="globalenv()", S=1))));
  try(ACPPNAMESPACE <- .Call("GetMyPackageNamespace"));
  try(eval(parse(text=GetG0Text("ACPPSAVESPACE", envir="ACPPNAMESPACE", S=1))));
  ##try(print(paste(".CreateList, ACPPSAVESPACE = ", sep="")));
  flush.console();
  try(print(ACPPSAVESPACE));  flush.console();
  try(eval(parse(text=GetG0Text(".AcppList", envir="ACPPSAVESPACE", S=1))));
  if (is.null(.AcppList) || length(.AcppList) <= 0 ||
    (is.numeric(.AcppList) && length(.AcppList) == 1 && .AcppList[1] == 0)) {
    eval(parse(text=GetG0Text(".AcppList", "ACPPSAVESPACE", S=1)));
    print(" -- ACPP: Creating .AcppList for first time. "); flush.console();       
    .AcppList <- vector("list", 5)
    try(eval(parse(text=SetGText(".AcppList", "ACPPSAVESPACE", S=1))));
    eval(parse(text=LockGText(".AcppList", envir="ACPPSAVESPACE", S=1)));
    print(paste(" -- ACPP: After creation of .Acpp, length of .AcppList is ", length(.AcppList), sep=""));
    flush.console();
  }
  ##eval(parse(text=LockGText("ACPPSAVESPACE", "ACPPNAMESPACE", S=1)));
  eval(parse(text=LockGText("ACPPNAMESPACE", envir="globalenv()", S=1)));
}

.DoubleAList <- function() {
  eval(parse(text=GetG0Text("ACPPNAMESPACE", envir="globalenv()", S=1)));
  try(ACPPNAMESPACE <- .Call("GetMyPackageNamespace"));  
  eval(parse(text=GetG0Text("ACPPSAVESPACE", envir="ACPPNAMESPACE", S=1)));
  eval(parse(text=GetG0Text(".AcppList", envir="ACPPSAVESPACE", S=1)));
  try(print(paste(".DoubleAList, ACPPSAVESPACE = ", sep="")));
  try(print(ACPPSAVESPACE)); flush.console();    
  if (is.null(.AcppList) || length(.AcppList) <= 0 ||
    (is.numeric(.AcppList) && length(.AcppList) == 1 && .AcppList[1] == 0)) {
      .CreateList(); 
  }
  NewACppList <- .Call("DoubleUpaList", .AcppList);
  .AcppList <- NewACppList;  NewACppList <- NULL;
  eval(parse(text=LockGText(".AcppList", envir="ACPPSAVESPACE", S=1)));
  eval(parse(text=LockGText("ACPPSAVESPACE", envir="ACPPNAMESPACE", S=1)));  
  eval(parse(text=LockGText("ACPPNAMESPACE", envir="globalenv()", S=1)));  
}

.InsertIntoAList <- function(AOb) {
  print(" -- Running .InsertIntoAList, Start. "); flush.console();
  eval(parse(text=GetG0Text("ACPPNAMESPACE", envir="globalenv()", S=1)));
  try(ACPPNAMESPACE <- .Call("GetMyPackageNamespace"));  
  try(eval(parse(text=GetG0Text("ACPPSAVESPACE", envir="ACPPNAMESPACE", S=1))));
  if (is.null(ACPPNAMESPACE) || (is.numeric(ACPPNAMESPACE) && ACPPNAMESPACE[1] == 0)) {
    print("  -- Have to set ACPPNAMESPACE again. "); flush.console();
    try(ACPPNAMESPACE <- .Call("GetMyPackageNamespace"));  
    eval(parse(text=SetGText("ACPPNAMESPACE", envir="globalenv()", S=1)));
    eval(parse(text=SetGText("ACPPNAMESPACE", envir="ACPPNAMESPACE", S=1))); 
  }
  if (is.null(ACPPSAVESPACE) || (is.numeric(ACPPSAVESPACE) && ACPPSAVESPACE[1] == 0)) {
    print("  -- Insert into a list Have to set ACPPSAVESPACE again. "); flush.console();
    ACPPSAVESPACE <- new.env();  
    try(eval(parse(text=SetGText("ACPPSAVESPACE", envir="ACPPNAMESPACE", S=1))));
    eval(parse(text=SetGText("ACPPSAVESPACE", envir="globalenv()", S=1))); 
  }
  print(" --- .InsertIntoAList, getting .AcppList"); flush.console();
  try(eval(parse(text=GetG0Text(".AcppList", envir="ACPPSAVESPACE", S=1))));
    
  print(" -- Running InsertInto A List "); flush.console();
  
  if (is.null(.AcppList) || length(.AcppList) <= 0 ||
    (is.numeric(.AcppList) && length(.AcppList) == 1 && .AcppList[1] == 0)) {
      ##print(" -- ACpp: InsertIntoAList: AcppList was NULL, creating list. "); flush.console();
      .CreateList();
      ##print(" -- ACpp: InsertIntoAList: Completed .CreateList"); flush.console();

  } else {
    ##print("-- Acpp: InsertIntoAList: presumably, .AcppList was not zero, here it is"); flush.console();
    ##try(print(.AcppList)); flush.console();
  }   
  try(eval(parse(text=GetG0Text(".AcppList", envir="ACPPSAVESPACE", S=1))));
  print(" -- ACpp: About to lock AOb into .AcppList");  flush.console();
  print(" -- Here is AOb: ");  print(AOb); flush.console();
  if (!exists(".AcppList")) {
    print(" -- .AccpList: We don't know where .AcppList is?  why is this?"); flush.console();
  }  else {
    ##print("We supposedly believe .AcppList exists. "); flush.console();
  }
  if (length(.AcppList) == 1 && is.numeric(.AcppList) && .AcppList[1] == 0) {
    print(" -- InsertIntoAList: Big ERror, .AcppList is 0 which means it's not where we want it!");
    flush.console();
  }
  print(paste(" -- Now length .AcppList is ", length(.AcppList), sep="")); flush.console();
  print(paste(" --  this is .AcppList "));  print(.AcppList); flush.console();
  print(paste(" -- About to Lock in .AcppList", sep="")); flush.console();
  ID <- -2;
  try(ID <- .Call("LockMeIn", .AcppList, AOb));
  if (ID == -2) {
    print("--------------------------------------------------------");
    print("-- InsertIntoAList, error, during LockMeIn.             ");
    print(paste("--  We believe length .AcppList is ", length(.AcppList))); flush.console();
    print("-- Why did this fail!");
    print("-- Saving .AcppList to GlobalEnv as Back .AcppList");
    BackAcppList <- .AcppList;
    eval(parse(text=SetGText("BackAcppList", envir="globalenv()", S=1)));
    BackAOb <- AOb;
    eval(parse(text=SetGText("BackAOb", "globalenv()", S=1)));
    print("-- Let's trigger InsertIntoAList Error. "); flush.console();
    eval(parse(text="1=2"));
  } else if (ID == -1) {
     print(paste(" -- InsertIntoAList, we are Doubling up AList length ", length(.AcppList), sep="")); flush.console();
     print(paste(" -- This could be very interesting. ")); flush.console();
     .AcppList <- .Call("DoubleUpaList", .AcppList);
     if (is.null(.AcppList) || length(.AcppList) <= 1) {
       print("-----------------------------------------------------------------");
       print("--- Error: DoubleUpaList failed!  Acpp is Broken"); flush.console();
       print("--- Memory is out and so the program will fail. "); flush.console();
       eval(parse(text="1=2"));
     }
     ID <- .Call("LockMeIn", .AcppList, AOb);
     ##print(paste("-- Here after DoubleUpAList and insert, AcppList is length ",
     ##  length(.AcppList), sep="")); flush.console();
     ##print("-- Here is the whole ACppList"); flush.console();
     ##print(.AcppList); flush.console();
     if (ID == -1) {
        print("----------------------------------------------------------------");
        print("--- InsertIntoAList, error, we ran LockMeIn but ID is still negative!");
        print("--- Sign that error occured. "); flush.console();
     }
  }  
  eval(parse(text=LockGText(".AcppList", "ACPPSAVESPACE", S=1)));
  ##eval(parse(text=LockGText("ACPPSAVESPACE", "ACPPNAMESPACE", S=1)));
  eval(parse(text=LockGText("ACPPNAMESPACE", "globalenv()", S=1)));
  return(ID);
}

.GetAcppList <- function() { 
  eval(parse(text=GetG0Text("ACPPNAMESPACE", "globalenv()", S=1)));
  eval(parse(text=GetG0Text("ACPPSAVESPACE", "ACPPNAMESPACE", S=1)));
  eval(parse(text=GetG0Text(".AcppList", "ACPPSAVESPACE", S=1)));
  try(eval(parse(text=LockGText("ACPPNAMESPACE", "globalenv()", S=1))));
  try(eval(parse(text=LockGText(".AcppList", "ACPPSAVESPACE",S=1))));
  return(.AcppList);
    
}

.GetACPPSAVESPACE <- function() { 
  eval(parse(text=GetG0Text("ACPPNAMESPACE", "globalenv()", S=1)));
  ##print(" -- We are running .GetACPPSAVESPACE "); flush.console();
  if (!exists("ACPPSAVESPACE", ACPPNAMESPACE)) {
    print("Acpp::Acpp.r: .GetACPPSAVESPACE: Hey, ACPPSAVESPACE not in ACPPNAMESPACE");
    flush.console();
  }
  eval(parse(text=GetG0Text("ACPPSAVESPACE", "ACPPNAMESPACE", S=1)));
  try(eval(parse(text=LockGText("ACPPNAMESPACE", "globalenv()", S=1))));
  return(ACPPSAVESPACE);
    
}


.DeleteFromAList <- function(AnInt) {
  if (AnInt < 0) { return(-1); }
  ##print(paste("--- .DeleteFromAList: R running. AInt=", AnInt, sep="")); flush.console();
  ##eval(parse(text=GetG0Text("ACPPNAMESPACE", "globalenv()", S=1)));  
  try(ACPPNAMESPACE <- .Call("GetMyPackageNamespace"));
  eval(parse(text=GetG0Text("ACPPSAVESPACE", "ACPPNAMESPACE", S=1)));  
  eval(parse(text=GetG0Text(".AcppList", "ACPPSAVESPACE", S=1)));    
  if (AnInt >= length(.AcppList)) { return(-2); }
  ##print(".DeleteFromAList: Now Running KillFromaList"); flush.console();
  .Call("KillFromaList", .AcppList, AnInt) 
  ##print(paste("--- .DeleteFromAList: Successfully ran .Call(KillFromaList, ", AnInt,")", sep="")); flush.console();
  eval(parse(text=SetGText(".AcppList", "ACPPSAVESPACE", S=1)));
  ##try(eval(parse(text=SetGText("ACPPSAVESPACE", "ACPPNAMESPACE", S=1))));
  eval(parse(text=SetGText("ACPPNAMESPACE", "globalenv()", S=1)));
  return(1); 
}


###############################################################################
## Old School linker code adapted from Rcpp(0.10.6):RcppLpath.R
##
##

#### HACK from Rcpp

## make sure system.file returns an absolute path
Acpp.system.file <- function(...){
    tools::file_path_as_absolute( base::system.file( ..., package = "Acpp" ) )
}

## identifies if the default linking on the platform should be static
## or dynamic. Currently only linux uses dynamic linking by default
## although it works fine on mac osx as well
staticLinking <- function() {
    ! grepl( "^linux", R.version$os )
}

## Use R's internal knowledge of path settings to find the lib/ directory
## plus optinally an arch-specific directory on system building multi-arch
AcppLdPath <- function() {
    if (nzchar(.Platform$r_arch)) {	## eg amd64, ia64, mips
        path <- Acpp.system.file("lib",.Platform$r_arch)
    } else {
        path <- Acpp.system.file("lib")
    }
    path
}

## Provide linker flags -- i.e. -L/path/to/libRcpp -- as well as an
## optional rpath call needed to tell the Linux dynamic linker about the
## location.  This is not needed on OS X where we encode this as library
## built time (see src/Makevars) or Windows where we use a static library
## Updated Jan 2010:  We now default to static linking but allow the use
##                    of rpath on Linux if static==FALSE has been chosen
##                    Note that this is probably being called from LdFlags()
AcppLdFlags <- function(static=staticLinking()) {
    rcppdir <- RcppLdPath()
    if (static) {                               # static is default on Windows and OS X
        flags <- paste(rcppdir, "/libAcpp.a", sep="")
        if (.Platform$OS.type=="windows") {
            flags <- asBuildPath(flags)
        }
    } else {					# else for dynamic linking
        flags <- paste("-L", rcppdir, " -lAcpp", sep="") # baseline setting
        if ((.Platform$OS.type == "unix") &&    # on Linux, we can use rpath to encode path
            (length(grep("^linux",R.version$os)))) {
            flags <- paste(flags, " -Wl,-rpath,", rcppdir, sep="")
        }
    }
    invisible(flags)
}

# indicates if Rcpp was compiled with GCC >= 4.3
canUseCXX0X <- function() .Call( "canUseCXX0X", PACKAGE = "Rcpp" )

## Provide compiler flags -- i.e. -I/path/to/Rcpp.h
RcppCxxFlags <- function(cxx0x=FALSE) {
    # path <- RcppLdPath()
    path <- Rcpp.system.file( "include" )
    if (.Platform$OS.type=="windows") {
        path <- asBuildPath(path)
    }
    paste("-I", path, if( cxx0x && canUseCXX0X() ) " -std=c++0x" else "", sep="")
}

## Shorter names, and call cat() directly
## CxxFlags defaults to no using c++0x extensions are these are considered non-portable
CxxFlags <- function(cxx0x=FALSE) {
    cat(RcppCxxFlags(cxx0x=cxx0x))
}
## LdFlags defaults to static linking on the non-Linux platforms Windows and OS X
LdFlags <- function(static=staticLinking()) {
    cat(RcppLdFlags(static=static))
}

print("Acpp.r: Complete Loading"); flush.console();
