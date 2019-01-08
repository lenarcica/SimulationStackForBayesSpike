GenerateDefaultSimulationFunctions <- function() {
 eval(parse(text=GetG0Text("TWOSIMENVIRONMENT", S=1)));
   if (is.numeric(TWOSIMENVIRONMENT) && TWOSIMENVIRONMENT == 0) {
     TWOSIMENVIRONMENT = environment();
     eval(parse(text=SetGText("TWOSIMENVIRONMENT", envir="globalenv()", S=1)));
   }
   
   eval(parse(text=GetG0Text("SimSigmaNoise", envir="globalenv()", S=1)));
   SimSigmaNoise <- function (SigmaNoise) {
       if (sigSV > 0) {
         return(sqrt(SigmaNoise^2 * rchisq(1, sigSV) /sigSV));
       } else{
         return(SigmaNoise);
       }
    }
    eval(parse(text=SetGText("SimSigmaNoise", envir="globalenv()", S=1)));
    eval(parse(text=SetGText("SimSigmaNoise", envir="TWOSIMENVIRONMENT", S=1)));    


  ## Create One form of Positive definite correlation data matrix
  eval(parse(text=GetG0Text("CorrelationXmatrix1", envir="globalenv()", S=1)));
  CorrelationXmatrix1 <- function(kLen, PosCorr) {
    if (abs(PosCorr) > 1) {
       print("CorrleationXmatrix1, PosCorr is too large, doesn't work\n");
   }
   if (PosCorr == 0) {
         list(
     corrX = diag(kLen),
     type="CorrelationXmatrix1", PosCorr = PosCorr )
   }
   return( list(
      corrX = PosCorr^( abs( matrix(rep((1:kLen),kLen),kLen,kLen) - 
                          t(matrix(rep((1:kLen),kLen), kLen,kLen)  ))  ),
      type="CorrelationXmatrix1", PosCorr = PosCorr )
   );
  }
  eval(parse(text=SetGText("CorrelationXmatrix1", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("CorrelationXmatrix1", envir="TWOSIMENVIRONMENT", S=1)));   

  eval(parse(text=GetG0Text("SinerX", envir="globalenv()", S=1)));
  SinerX <- function(Inpts, a,b) {
   rt <-  sin(a * Inpts) / ( a * Inpts);
  }
  eval(parse(text=SetGText("SinerX", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("SinerX", envir="TWOSIMENVIRONMENT", S=1)));   

  eval(parse(text=GetG0Text("CorrelationXmatrix2", envir="globalenv()", S=1)));
  CorrelationXmatrix2 <- function(kLen, a,b) {
  if (abs(PosCor) > 1) {
      print("CorrleationXmatrix1, PosCorr is too large, doesn't work\n");
  }
  return( list(
     corrX = SinerX( abs( matrix(rep((1:kLen),kLen),kLen,kLen) - 
                          t(matrix(rep((1:kLen),kLen), kLen,kLen)  )) , a,b  ),
     type="CorrelationXmatrix2", a=a, b=b
     )
  );
  }
  eval(parse(text=SetGText("CorrelationXmatrix2", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("CorrelationXmatrix2", envir="TWOSIMENVIRONMENT", S=1)));   


  SetBetaVecs();
  
  eval(parse(text=GetG0Text("GenerateBetaVec2", S=1, envir="globalenv()")));
  eval(parse(text=GetG0Text("DefaultGenerateBetaVec", S=1, envir="globalenv()")));
  eval(parse(text=GetG0Text("GenerateBetaVec", S=1, envir="globalenv()")));
  DefaultGenerateBetaVec <- GenerateBetaVec2;
  GenerateBetaVec = DefaultGenerateBetaVec
  eval(parse(text=SetGText("DefaultGenerateBetaVec", S=1, envir="globalenv()")));
  eval(parse(text=SetGText("DefaultGenerateBetaVec", envir="TWOSIMENVIRONMENT", S=1)));  
  eval(parse(text=SetGText("GenerateBetaVec", S=1, envir="globalenv()")));
  eval(parse(text=SetGText("GenerateBetaVec", envir="TWOSIMENVIRONMENT", S=1)));
   
  eval(parse(text=GetG0Text("DefaultCorrelationXmatrix", S=1, envir="globalenv()"))); 
  eval(parse(text=GetG0Text("CorrelationXmatrix", S=1, envir="globalenv()"))); 
  DefaultCorrelationXmatrix <- function(kLen) {
           CorrelationXmatrix1(kLen, .25);
  }
  CorrelationXmatrix = DefaultCorrelationXmatrix;
  eval(parse(text=SetGText("DefaultCorrelationXmatrix", S=1, envir="globalenv()")));
  eval(parse(text=SetGText("DefaultCorrelationXmatrix", envir="TWOSIMENVIRONMENT", S=1)));  
  eval(parse(text=SetGText("CorrelationXmatrix", S=1, envir="globalenv()")));
  eval(parse(text=SetGText("CorrelationXmatrix", envir="TWOSIMENVIRONMENT", S=1)));  
  
 
  SetBetaVecNameFunctions();
  SetBetaGroupVecs();
  
  eval(parse(text=GetG0Text("GenerateBetaGroupVec2", S=1, envir="globalenv()")));
  eval(parse(text=GetG0Text("DefaultGenerateBetaGroupVec", S=1, envir="globalenv()")));
  eval(parse(text=GetG0Text("GenerateBetaGroupVec", S=1, envir="globalenv()")));
  DefaultGenerateBetaGroupVec <- GenerateBetaGroupVec2;
  GenerateBetaGroupVec = DefaultGenerateBetaGroupVec
  eval(parse(text=SetGText("DefaultGenerateBetaGroupVec", S=1, envir="globalenv()")));
  eval(parse(text=SetGText("DefaultGenerateBetaGroupVec", envir="TWOSIMENVIRONMENT", S=1)));  
  eval(parse(text=SetGText("GenerateBetaGroupVec", S=1, envir="globalenv()")));
  eval(parse(text=SetGText("GenerateBetaGroupVec", envir="TWOSIMENVIRONMENT", S=1)));
}

############################################################################
## SetBetaVecs();
##
##  These are default ways to simulate Beta Vectors
##
SetBetaVecs <- function() {
  eval(parse(text=GetG0Text("TWOSIMENVIRONMENT", S=1)));

## GenerateBetaVec1 
  GenerateBetaVec1 <- function(kLen = 100, kAwant = 6, InOrder = FALSE,...) {
     ## function GenerateBetaVec1: In 0006GenerateDefaultSimulationFunctions.R
     if (kAwant > kLen) { print("GenerateBetaVec1 Entry Error");
                          return(-1);
     }
     BetaVec = c(rep(1, kAwant), rep(0, NROn-kAwant))
     if (InOrder==FALSE) {
       BetaVec = sample(BetaVec, NROn, replace=FALSE);
     }
     return(list(
      BetaVec =BetaVec,        
         type="GenerateBetaVec1",
         kLen = kLen, kAwant = kAwant, InOrder=InOrder
             )
        );  
  }
 eval(parse(text=SetGText("GenerateBetaVec1", envir="globalenv()", S=1)));
 eval(parse(text=SetGText("GenerateBetaVec1", envir="TWOSIMENVIRONMENT", S=1)));  
  
  GenerateBetaVec2 <- function(kLen = 100, kAwant = 6, InOrder=FALSE,...) {
     ## function GenerateBetaVec2: In 0006GenerateDefaultSimulationFunctions.R
     if (kAwant > kLen) { print("GenerateBetaVec2 Entry Error");
                          return(-1);
     }
     BetaVec =c(rep(1, floor(kAwant/2)),
             rep(-1, floor(kAwant/2)),
             rep(1, kAwant - 2 * floor(kAwant/2)), 
             rep(0, kLen-kAwant))
     if (InOrder==FALSE) {
       BetaVec = sample(BetaVec, kLen, replace=FALSE);
     }
     return(list(
        BetaVec=BetaVec,
         type="GenerateBetaVec2",
           kLen = kLen, kAwant = kAwant, InOrder=InOrder
             )             
           ); 
  }
  eval(parse(text=SetGText("GenerateBetaVec2", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaVec2", envir="TWOSIMENVIRONMENT", S=1)));   
  
  GenerateBetaVec2a <- function(kLen = 100, kAwant = 6, InOrder=FALSE,...) {
     ## function GenerateBetaVec2a: In 0006GenerateDefaultSimulationFunctions.R
     if (kAwant > kLen) { print("GenerateBetaVec2 Entry Error");
                          return(-1);
     }
     fLL = floor(kAwant/2);
     if (kAwant %% 2 == 1) {
        BetaVec = c( t( cbind(rep(1,fLL),rep(-1,fLL)) ), 
                      1, rep(0, kLen-kAwant)
                   )
     } else {
        BetaVec = c( t( cbind(rep(1,fLL),rep(-1,fLL)) ), 
                       rep(0, kLen-kAwant)
                   );     
     }
     if (InOrder==FALSE) {
       BetaVec = sample(BetaVec, kLen, replace=FALSE);
     }
     return(list(
        BetaVec=BetaVec,
           type="GenerateBetaVec2a",
           kLen = kLen, kAwant = kAwant, InOrder=InOrder)             
           ); 
  } 
  eval(parse(text=SetGText("GenerateBetaVec2a", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaVec2a", envir="TWOSIMENVIRONMENT", S=1)));    
  
  GenerateBetaVec2b <- function(kLen = 100, kAwant = 6, InOrder=FALSE,...) {
     ## function GenerateBetaVec2b: In 0006GenerateDefaultSimulationFunctions.R
     if (kAwant > kLen) { print("GenerateBetaVec2 Entry Error");
                          return(-1);
     }
     fLL = floor(kAwant/2);
     ILLoc <- sample(1:kLen, kAwant, replace=FALSE);
     if (kAwant %% 2 == 1) {
        BetaVec = rep(0, kLen);
        BetaVec[ILLoc[1: (kAwant %/% 2 + 1)]] = 1;
                BetaVec[ILLoc[((kAwant %/% 2 + 2)):kAwant]] = -1;   
     } else {
        BetaVec = rep(0, kLen);
        BetaVec[ILLoc[1: (kAwant %/% 2 )]] = 1;
                BetaVec[ILLoc[((kAwant %/% 2) + 1):kAwant]] = -1;      
     }
     if (InOrder==FALSE) {
        BetaVec = sample(BetaVec, kLen, replace=FALSE)
     }
     return(list(
        BetaVec=BetaVec,
           type="GenerateBetaVec2b",
           kLen = kLen, kAwant = kAwant, InOrder=InOrder
             )             
           ); 
  } 
  eval(parse(text=SetGText("GenerateBetaVec2b", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaVec2b", envir="TWOSIMENVIRONMENT", S=1)));   
  
  GenerateBetaVec3 <- function(kLen = 100, kAwant = 10, MN = 100, InOrder=FALSE,...) {
    if (kAwant > kLen) { print("GenerateBetaVec3 Entry Error");
                       return(-1);
    }
    BetaVec = c( runif(kAwant, -MN,MN), rep(0, kLen-kAwant) )
     if (InOrder==FALSE) {
       BetaVec = sample(BetaVec, kLen, replace=FALSE);
     }
  return(list(
    BetaVec= sample(,
      kLen, replace=FALSE )), 
      type="GenerateBetaVec3", kLen = kLen, kAwant = kAwant, MN = MN, InOrder=InOrder);  
  }
  eval(parse(text=SetGText("GenerateBetaVec3", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaVec3", envir="TWOSIMENVIRONMENT", S=1)));   
  

}


############################################################################
## SetBetaVecs();
##
##  These are default ways to simulate Beta Vectors
##
SetBetaGroupVecs <- function() {
  eval(parse(text=GetG0Text("TWOSIMENVIRONMENT", S=1)));

## GenerateBetaVec1 
  eval(parse(text=GetG0Text("GenerateBetaGroupVec1", envir="globalenv()", S=1)));
  GenerateBetaGroupVec1 <- function(pGroups = pGroups, GroupSize=GroupSize, 
    kGroups = kGroups, GroupsSumToZero=TRUE, InOrder=FALSE, ElseTauEndList=NULL, ...) {
     ## function GenerateGroupBetaVec1: In 0006GenerateDefaultSimulationFunctions.R
     if (!is.null(ElseTauEndList)) {
        pGroups = length(ElseTauEndList);
     }
     if (kGroups > pGroups) { print("text Entry Error");
                          return(-1);
     }
     if (InOrder==FALSE) {
       ASampleGroups <- sample( c(rep(1, kGroups), rep(0, pGroups-kGroups)),
         pGroups, replace=FALSE );
     } else {
      ASampleGroups = c(rep(1, kGroups), rep(0, pGroups-kGroups));  
     }
    if (!is.null(ElseTauEndList)) {
      MyBeta <- rep(0, max(ElseTauEndList));
      ASG <- (1:length(ASampleGroups))[ASampleGroups == 1];
      for (ii in 1:length(ASG)) {
        if (ASG[ii] == 1) {
          GroupSize = ElseTauEndList[ASG[ii]]; St = 1;
        } else {
          St = ElseTauEndList[ASG[ii]-1]+1;
          GroupSize = ElseTauEndList[ASG[ii]] - ElseTauEndList[ASG[ii]-1];
        }
        if ((GroupSize %% 2) == 1)  {
          GSS <- c(rep(1, GroupSize %/% 2),
           rep(-1, GroupSize %/% 2),0);
        } else {
          GSS <- c(rep(1, GroupSize %/% 2),
            rep(-1, GroupSize %/% 2));
        }
        AM <- sample(GSS, replace=FALSE); AM = AM-mean(AM);
        MyBeta[St:ElseTauEndList[ASG[ii]]] <- AM;
      }
    } else {
     MyBeta <- rep(0, GroupSize*pGroups);
     if ((GroupSize %% 2) == 1)  {
       GSS <- c(rep(1, GroupSize %/% 2),
         rep(-1, GroupSize %/% 2),0);
     } else {
       GSS <- c(rep(1, GroupSize %/% 2),
         rep(-1, GroupSize %/% 2));
     }
     for (ii in 1:pGroups) {
       if (ASampleGroups[ii] == 1) {
         AM <- sample(GSS, replace=FALSE); AM = AM - mean(AM);
         MyBeta[(ii-1)*GroupSize + 1:GroupSize] <-  AM
       }
     } 
    }
     return(list(
       BetaVec = MyBeta, ASampleGroups=ASampleGroups,
       type="GenerateBetaGroupVec1",
       pGroups = pGroups, kGroups = kGroups,
       GroupSize=GroupSize, p = pGroups*kGroups, k=GroupSize*kGroups,
       FirstGroupIndex=1, EndGroupIndices = (1:pGroups)*GroupSize, InOrder=InOrder
      ));  
  }
  eval(parse(text=SetGText("GenerateBetaGroupVec1", 
    envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaGroupVec1", 
    envir="TWOSIMENVIRONMENT", S=1)));  
 
  eval(parse(text=GetG0Text("GenerateBetaGroupVec2", envir="globalenv()", S=1))); 
  GenerateBetaGroupVec2 <- function(pGroups = pGroups, GroupSize=GroupSize, 
    kGroups = kGroups, GroupsSumToZero=FALSE, InOrder=FALSE,ElseTauEndList = NULL, ...) {
    ## function GenerateBetaVec2: In 0006GenerateDefaultSimulationFunctions.R

    if (!is.null(ElseTauEndList)) {
        pGroups = length(ElseTauEndList);
    }
     if (kGroups > pGroups) { print("GenerateBetaGroupVec2 Entry Error");
                          return(-1);
    }
    if (InOrder==FALSE) {
      ASampleGroups <- sample( c(rep(1, kGroups), rep(0, pGroups-kGroups)),
         pGroups, replace=FALSE );
    } else {
      ASampleGroups <- c(rep(1, kGroups), rep(0, pGroups-kGroups));         
    }
    if (!is.null(ElseTauEndList)) {
      MyBeta <- rep(0, max(ElseTauEndList));
      ASG <- (1:length(ASampleGroups))[ASampleGroups == 1];
      p = max(ElseTauEndList);
      for (ii in 1:length(ASG)) {
        if (ASG[ii] == 1) {
          GroupSize = ElseTauEndList[ASG[ii]]; St = 1;
        } else {
          St = ElseTauEndList[ASG[ii]-1]+1;
          GroupSize = ElseTauEndList[ASG[ii]] - ElseTauEndList[ASG[ii]-1];
        }
        if ((GroupSize %% 2) == 1)  {
          GSS <- c(rep(1, GroupSize %/% 2),
           rep(-1, GroupSize %/% 2),0);
        } else {
          GSS <- c(rep(1, GroupSize %/% 2),
            rep(-1, GroupSize %/% 2));
        }
        AM <- sample(GSS, replace=FALSE); AM = AM-mean(AM);
        MyBeta[St:ElseTauEndList[ASG[ii]]] <- AM;
      }
    } else {
    MyBeta <- rep(0, GroupSize*pGroups);
     if ((GroupSize %% 2) == 1)  {
       GSS <- c(rep(1, GroupSize %/% 2),
         rep(-1, GroupSize %/% 2),0);
     } else {
       GSS <- c(rep(1, GroupSize %/% 2),
         rep(-1, GroupSize %/% 2));
     }
     for (ii in 1:pGroups) {
       if (ASampleGroups[ii] == 1) {
         AM <- sample(GSS, replace=FALSE)
         AM = AM-mean(AM);
         MyBeta[(ii-1)*GroupSize + 1:GroupSize] <-  AM
       }
     } 
    }
     return(list(
       BetaVec = MyBeta,type="GenerateBetaGroupVec2",
         pGroups = pGroups, kGroups = kGroups, GroupSize=GroupSize,
         p=length(MyBeta), k=kGroups*GroupSize,
       FirstGroupIndex=1, EndGroupIndices = (1:pGroups)*GroupSize, InOrder=InOrder
         )); 
  }
  eval(parse(text=SetGText("GenerateBetaGroupVec2", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaGroupVec2", envir="TWOSIMENVIRONMENT", S=1)));   
  
  eval(parse(text=GetG0Text("GenerateBetaGroupVec2a", envir="globalenv()", S=1)));
  GenerateBetaGroupVec2a <- function(pGroups = pGroups, GroupSize=GroupSize, 
    kGroups = kGroups, Spread = 1, GroupsSumToZero=FALSE,
    InOrder=FALSE,ElseTauEndList = NULL, ...) {
     ## function GenerateBetaVec2a: In 0006GenerateDefaultSimulationFunctions.R
    if (!is.null(ElseTauEndList)) {
        pGroups = length(ElseTauEndList);
    }
     if (kGroups > pGroups) { print("GenerateBetaVec2a Entry Error");
                          return(-1);
     }
    if (InOrder==FALSE) {
      ASampleGroups <- sample( c(rep(1, kGroups), rep(0, pGroups-kGroups)),
         pGroups, replace=FALSE );
    } else {
      ASampleGroups <- c(rep(1, kGroups), rep(0, pGroups-kGroups))  
    }
                  k <- 0;
  if (!is.null(ElseTauEndList)) {
      MyBeta <- rep(0, max(ElseTauEndList));
      ASG <- (1:length(ASampleGroups))[ASampleGroups == 1];
      for (ii in 1:length(ASG)) {
        if (ASG[ii] == 1) {
          GroupSize = ElseTauEndList[ASG[ii]]; St = 1;
        } else {
          St = ElseTauEndList[ASG[ii]-1]+1;
          GroupSize = ElseTauEndList[ASG[ii]] - ElseTauEndList[ASG[ii]-1];
        }
        AM <- 2*rbinom(GroupSize, 1, .5)-1;
        AM = AM-mean(AM);   AM = AM * Spread;
        MyBeta[St:ElseTauEndList[ASG[ii]]] <- AM; k <- k + length(AM);
      }
    } else {
    MyBeta <- rep(0, GroupSize*pGroups);
     for (ii in 1:pGroups) {
       if (ASampleGroups[ii] == 1) {
         AM <- 2*rbinom(GroupSize, 1, .5)-1;
         AM = AM-mean(AM);   AM = AM * Spread;
         MyBeta[(ii-1)*GroupSize + 1:GroupSize] <-  AM
       }
     } 
     k <- GroupSize * kGroups;
    }
     return(list(
       BetaVec = MyBeta,type="GenerateBetaGroupVec2a",
         pGroups = pGroups, kGroups = kGroups, GroupSize=GroupSize,
         p=length(MyBeta), k=k, Spread=Spread,
       FirstGroupIndex=1, EndGroupIndices = (1:pGroups)*GroupSize,
       InOrder=InOrder)); 
  } 
  eval(parse(text=SetGText("GenerateBetaGroupVec2a", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaGroupVec2a", envir="TWOSIMENVIRONMENT", S=1)));    
 
  eval(parse(text=GetG0Text("GenerateBetaGroupVec2b", envir="globalenv()", S=1))); 
  GenerateBetaGroupVec2b <- function(pGroups = pGroups, GroupSize=GroupSize, 
    kGroups = kGroups, Spread=1, GroupsSumToZero=FALSE, InOrder=FALSE, ElseTauEndList = NULL, ...) {
     ## function GenerateBetaVec2b: In 0006GenerateDefaultSimulationFunctions.R
     if (!is.null(ElseTauEndList)) {
       pGroups = length(ElseTauEndList);
     }
     if (kGroups > pGroups) { print("GenerateBetaVec2b Entry Error");
       return(-1);
     }
     if (InOrder==FALSE) {
     ASampleGroups <- sample( c(rep(1, kGroups), rep(0, pGroups-kGroups)),
       pGroups, replace=FALSE );
     } else {
      ASampleGroups <- c(rep(1, kGroups), rep(0, pGroups-kGroups));
     }
     ## function GenerateBetaVec2a: In 0006GenerateDefaultSimulationFunctions.R
    if (!is.null(ElseTauEndList)) {
        pGroups = length(ElseTauEndList);
    }
     if (kGroups > pGroups) { print("GenerateBetaVec2a Entry Error");
                          return(-1);
     }
    if (InOrder==FALSE) {
      ASampleGroups <- sample( c(rep(1, kGroups), rep(0, pGroups-kGroups)),
         pGroups, replace=FALSE );
    } else {
      ASampleGroups <- c(rep(1, kGroups), rep(0, pGroups-kGroups))  
    }
  
  k <- 0;
  if (!is.null(ElseTauEndList)) {
      MyBeta <- rep(0, max(ElseTauEndList));
      ASG <- (1:length(ASampleGroups))[ASampleGroups == 1];
      for (ii in 1:length(ASG)) {
        if (ASG[ii] == 1) {
          GroupSize = ElseTauEndList[ASG[ii]]; St = 1;
        } else {
          St = ElseTauEndList[ASG[ii]-1]+1;
          GroupSize = ElseTauEndList[ASG[ii]] - ElseTauEndList[ASG[ii]-1];
        }
        k <- k + GroupSize;
        AM <- Spread* rnorm(GroupSize, 0, 1);
         if (GroupsSumToZero == TRUE) { AM = AM-mean(AM);  }
        MyBeta[St:ElseTauEndList[ASG[ii]]] <- AM;
      }
    } else {
     MyBeta <- rep(0, GroupSize*pGroups);
     for (ii in 1:pGroups) {
       if (ASampleGroups[ii] == 1) {
         AM <- Spread* rnorm(GroupSize, 0, 1);
         if (GroupsSumToZero == TRUE) { AM = AM-mean(AM);  }
         MyBeta[(ii-1)*GroupSize + 1:GroupSize] <-  AM
       }
     } 
     k <- kGroups * GroupSize;
     }
     return(list(
       BetaVec = MyBeta,type="GenerateBetaGroupVec2b",
         pGroups = pGroups, kGroups = kGroups, GroupSize=GroupSize,
         p=length(MyBeta), k=k,
       FirstGroupIndex=1, EndGroupIndices = (1:pGroups)*GroupSize, InOrder=InOrder
         )); 
  } 
  eval(parse(text=SetGText("GenerateBetaGroupVec2b", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaGroupVec2b", envir="TWOSIMENVIRONMENT", S=1)));   
 
  eval(parse(text=GetG0Text("GenerateBetaGroupVec3", envir="globalenv()", S=1))); 
  GenerateBetaGroupVec3 <- function(pGroups = pGroups, GroupSize=GroupSize, 
    kGroups = kGroups, Spread=1, GroupsSumToZero=FALSE, InOrder=FALSE,ElseTauEndList = NULL, ...) {
     ## function GenerateBetaVec2b: In 0006GenerateDefaultSimulationFunctions.R
     if (!is.null(ElseTauEndList)) {
       pGroups = length(ElseTauEndList);
     }
    if (kGroups > pGroups) { print("GenerateBetaVec3 Entry Error");
                       return(-1);
    }
    if (InOrder==FALSE) {
      ASampleGroups <- sample( c(rep(1, kGroups), rep(0, pGroups-kGroups)),
        pGroups, replace=FALSE );
    } else {
      ASampleGroups <- c(rep(1, kGroups), rep(0, pGroups-kGroups))   
    }
    k <- 0;
    if (!is.null(ElseTauEndList)) {
      MyBeta <- rep(0, max(ElseTauEndList));
      ASG <- (1:length(ASampleGroups))[ASampleGroups == 1];
      for (ii in 1:length(ASG)) {
        if (ASG[ii] == 1) {
          GroupSize = ElseTauEndList[ASG[ii]]; St = 1;
        } else {
          St = ElseTauEndList[ASG[ii]-1]+1;
          GroupSize = ElseTauEndList[ASG[ii]] - ElseTauEndList[ASG[ii]-1];
        }
        k <- k + GroupSize;
        AM <- Spread* rnorm(GroupSize, 0, 1);
         if (GroupsSumToZero==TRUE) { AM = AM-mean(AM); }
        MyBeta[St:ElseTauEndList[ASG[ii]]] <- AM;
      }
    } else {
     MyBeta <- rep(0, GroupSize*pGroups);
     for (ii in 1:pGroups) {
       if (ASampleGroups[ii] == 1) {
         AM <- Spread* rnorm(GroupSize, 0, 1);
         if (GroupsSumToZero==TRUE) { AM = AM-mean(AM); }
         MyBeta[(ii-1)*GroupSize + 1:GroupSize] <-  AM
       }
     } 
     k <- GroupSize * kGroups;
    }
     return(list(
       BetaVec = MyBeta,type="GenerateBetaGroupVec3",
         pGroups = pGroups, kGroups = kGroups, GroupSize=GroupSize,
         p=length(MyBeta), k=k,
       FirstGroupIndex=1, EndGroupIndices = (1:pGroups)*GroupSize, InOrder=InOrder
         )); 
  }
  eval(parse(text=SetGText("GenerateBetaGroupVec3", envir="globalenv()", S=1)));
  eval(parse(text=SetGText("GenerateBetaGroupVec3", envir="TWOSIMENVIRONMENT", S=1)));   
  

}

SetBetaVecNameFunctions <- function() {

eval(parse(text=GetG0Text("GetBetaName", envir="globalenv()", S=1)));
GetBetaName <- function(GenerateBetaVec) {
  BetaVD <- GenerateBetaVec(4,1);
    BetaName <- BetaVD$type;
    if (length(BetaVD) > 4) {
      for (ii in 5:length(BetaVD)) {
        BetaName <- paste(BetaName, names(BetaVD)[ii],
                                    tSeq(BetaVD[[ii]]), sep="" );
      }
    }
    return(BetaName);
}  
eval(parse(text=SetGText("GetBetaName", envir="globalenv()", S=1)));
eval(parse(text=SetGText("GetBetaName", envir="TWOSIMENVIRONMENT", S=1)));   
 
eval(parse(text=GetG0Text("GetCorrName", envir="globalenv()", S=1))); 
GetCorrName <- function(CorrelationXmatrix) {
   CorrVD <- CorrelationXmatrix(4);
   CorrName <- CorrVD$type
      for (ii in 3:length(CorrVD)) {
        CorrName <- paste(CorrName, names(CorrVD)[ii],
                                    tSeq(CorrVD[[ii]])
                     ,sep="" );
      }
    return(CorrName);
}
eval(parse(text=SetGText("GetCorrName", envir="globalenv()", S=1)));
eval(parse(text=SetGText("GetCorrName", envir="TWOSIMENVIRONMENT", S=1)));   
  

}