#######################################################################################
#######################################################################################
###  This code is for making Latex demonstration table summaries of simulation output
###
###

tSeq <- function (Number) 
{
    CN <- sub("e", "e", as.character(Number))
    CN <- sub("-", "n", as.character(CN))
    CN <- sub("\\+", "f", as.character(CN))
    CN <- sub("\\.", "p", as.character(CN))
    CN <- sub("0p", "p", as.character(CN))
    return(CN)
}

SetAllComparisonRCodeFunctions <- function() {

#########################################################################################
###
###   These are Latex titles, row and column headers
###
EstimatorNames<-  c("Lasso Fixed", "Lars Cp", "Lasso Lin and Yuan", "Limit Ridge", "Quick Two Lasso", 
     "Limit Lasso", "Marginal Median");
FunctionPlot <- paste(" $ \\begin{array} {c} ", "\\mbox{\\footnotesize{\\# II}} \\\\ \\hline",
                    "\\mbox{\\footnotesize{\\# I}} \\\\ \\hline", 
                    "\\mbox{\\footnotesize{$\\sum \\delta_{\\mbox{\\tiny{$\\beta$}}}^2$}} \\\\ \\hline",
                    ##"\\mbox{\\footnotesize{\\% Perf}} \\\\ \\hline ", 
                    "\\mbox{\\footnotesize{Run}}",
                    "\\end{array} $ ", sep="");
EstimatorColNames<-  c(
      paste("\\begin{array}{c} \\mbox{LARS} \\\\",
         " \\mbox{Fixed $\\kappa_{\\mbox{\\tiny{$\\mathcal{A}$}}}$",
         " } \\end{array} ", sep=""),
      paste("\\begin{array}{c} \\mbox{LARS} \\\\",
         " \\mbox{$C_{\\mbox{\\tiny{p}}}$",
         "} \\end{array}", sep=""),
      paste("\\begin{array}{c} \\mbox{\\small{L}\\footnotesize{asso $w$=1}}",
            " \\\\ \\mbox{\\small{L}\\footnotesize{in \\& }\\small{Y}\\foootnotesize{uan}}} ",
            " \\end{array}", sep=""), 
      paste( "\\begin{array}{c}",
         " \\mbox{\\small{L}\\footnotesize{im}",
         "\\small{R}\\footnotesize{idge}} \\\\",
         " \\mbox{$\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}$",
         " \\small{K}\\footnotesize{nown}} \\end{array}", sep=""), 
      paste("\\begin{array}{c} \\mbox{\\small{T}\\footnotesize{wo}",
         "\\small{L}\\footnotesize{asso}}\\\\",
         " \\mbox{$\\times$ 9} \\end{array}", sep=""), 
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{L}\\footnotesize{im}\\small{L}",
         "\\footnotesize{asso}} \\\\",
         " \\mbox{$\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}$",
         " \\small{K}\\footnotesize{nown}}",
         " \\end{array}", sep=""),
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{L}\\footnotesize{im}",
         "\\small{L}\\footnotesize{asso}} \\\\",
         " \\mbox{$\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}$",
         " Est.} \\end{array}", sep=""),
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{P}\\footnotesize{sd}-\\small{M}",
         "\\footnotesize{arg}} \\\\ \\mbox{\\small{M}\\footnotesize{edian}} ",
         "\\end{array}", sep=""),
      paste("\\begin{array}{c} \\mbox{Fermi-D} \\\\ ",
         "\\mbox{\\small{L}\\footnotesize{im}\\small{L}",
         "\\footnotesize{asso}} \\end{array}", sep=""),
      paste("\\begin{array}{c} \\mbox{\\small{M}\\footnotesize{arg}",
         "\\small{M}\\footnotesize{edian}} \\\\ ",
         "\\mbox{\\small{L}\\footnotesize{im}\\small{L}\\footnotesize{asso}} \\end{array}", sep="")
   );
TopPlot <- c(" \\begin{array}{c} \\mbox{Mean Type II}  \\\\ \\mbox{Real Factors Missed} \\\\ \\end{array} ",
  " \\begin{array}{c} \\mbox{Mean Type I}  \\\\ \\mbox{Real Factors Missed} \\\\ \\end{array} ",
  " \\begin{array}{c} \\mbox{\\% True Model} \\\\ \\end{array}",
  " \\begin{array}{c} \\mbox{SD Type II}  \\\\ \\mbox{Real Factors Missed} \\\\ \\end{array} ",
  " \\begin{array}{c} \\mbox{SD Type I}  \\\\ \\mbox{Real Factors Missed} \\\\ \\end{array} ",
  " \\begin{array}{c} \\sum \\left( \\hat{\\beta}_{j} - \\beta_{j-\\mbox{\\tiny{TRUE}}} \\right)^2 \\\\ \\end{array}",
  " \\begin{array}{c} \\mbox{Computation Time} \\\\ \\mbox{User (sec)}  \\end{array} ",
  " \\begin{array}{c} \\mbox{Computation Time} \\\\ \\mbox{Computer (sec)}  \\end{array} ",
  " \\begin{array}{c} \\mbox{Computation Time} \\\\ \\mbox{Total (sec)}  \\end{array} "
  )
EstimatorColNames2 <- paste( "$ ", EstimatorColNames, " $", sep="");
  eval(parse(text=SetGText("EstimatorColNames2",S=1)));
  eval(parse(text=SetGText("TopPlot",S=1)));
  eval(parse(text=SetGText("EstimatorNames",S=1)));
  eval(parse(text=SetGText("FunctionPlot",S=1)));
  eval(parse(text=SetGText("EstimatorColNames",S=1)));
}
#####################################################################################
###  rd0 is a function for Latex formatting of numbers to reduce their space occupied in tables
###
###
rd0 <- function(RoundNumber, SigFig = 3) {
  
if (length(RoundNumber) == 1) {
  if (is.na(RoundNumber)) {
    return("-");
  } else if (RoundNumber < 1.0 && round(RoundNumber,2) == 1.0) {
    return("1.");
  } else if (RoundNumber >= .001 && RoundNumber < .01) {
    MSSplit <- unlist(strsplit(as.character(round(RoundNumber,SigFig+1)), "\\."));
     if (length(MSSplit) == 1) {  
          MSSplit <- c(MSSplit,"0");
     }
     return( paste( ".", 
       (MSSplit)[2], sep=""));
  } else if (RoundNumber >= .01 && RoundNumber < 1) {
    MSSplit <- unlist(strsplit(as.character(round(RoundNumber,SigFig)), "\\."));
    if (length(MSSplit) == 1) {  
          MSSplit <- c(MSSplit,"0");
    }
     return( paste( ".", 
       (MSSplit)[2], sep=""));
  } else if (RoundNumber >= 100) {
     L2 <- floor(log(RoundNumber,10));
     MSSplit <- unlist(strsplit(as.character(round(RoundNumber/10^(L2),SigFig-1)), "\\."));
     if (length(MSSplit) == 1) {
        MSSplit <- c( MSSplit,"0");
     }
      return( paste("\\mbox{", MSSplit[1],".","{\\footnotesize", MSSplit[2],
            "e", L2, "", "\\normalsize}}", sep="") );
  } else if (RoundNumber >= 10) {
     MSSplit <- unlist(strsplit(as.character(round(RoundNumber,SigFig-2)), "\\."));
     if (length(MSSplit) == 1) {  
          MSSplit <- c(MSSplit,"0");
       }
    return( paste("\\mbox{", MSSplit[1],".","{\\footnotesize", MSSplit[2],"\\normalsize}}", sep="") );
  }  else if (RoundNumber >= 1 && RoundNumber < 10) {
    MSSplit <- unlist(strsplit(as.character(round(RoundNumber,SigFig-1)), "\\."));
     if (length(MSSplit) == 1) {  
          MSSplit <- c(MSSplit,"0");
       }    
    return( paste("\\mbox{", MSSplit[1],".","{\\footnotesize", MSSplit[2],"\\normalsize}}", sep="") );    
  } else if (RoundNumber > 0 && RoundNumber < .01) {
     L2 <- floor(log(RoundNumber,10));
     MSSplit <- unlist(strsplit(as.character(round(RoundNumber/10^(L2),SigFig-1)), "\\."));
     if (length(MSSplit) == 1) {
        MSSplit <- c( MSSplit,"0");
     }
      return( paste("\\mbox{", MSSplit[1],".","{\\footnotesize", MSSplit[2],
            "e", L2, "", "\\normalsize}}", sep="") ); 
  } else if (RoundNumber == 0) {
     return("\\mbox{0.\\footnotesize{0}}");  
  } else {
     return(as.character(round(RoundNumber,2)));
    }
  } else {
    RTV <- RoundNumber;
    for (ii in 1:length(RoundNumber)) {
       RTV[ii] = rd0(RoundNumber[ii], SigFig=SigFig);
    }
    return(RTV);
    RTV[RoundNumber >= 0 & RoundNumber < 1] <- 
      paste( ".", 
    (unlist(strsplit(as.character(RoundNumber[RoundNumber >= 0 & RoundNumber < 1]), "\\."))[2]), 
       sep="")
     MSSplit <- unlist(strsplit(as.character(round(RoundNumber[RoundNumber >= 1 & RoundNumber < 10],2)), "\\."));
     RTV[RoundNumber >= 1 & RoundNumber < 10] <- paste("\\mbox{", MSSplit[1],".","{\\footnotesize", MSSplit[2],"\\normalsize}}", sep="");       
     MSSplit <- unlist(strsplit(as.character(round(RoundNumber[RoundNumber >= 10],1)), "\\."));
     RTV[RoundNumber >= 10] <- paste("\\mbox{", MSSplit[1],".","{\\footnotesize", MSSplit[2],"\\normalsize}}", sep="");
     return(RTV);
  }
}
    ##############################
    ## BoxArangement  : 
    ##   Mean Type II  (sd Type II)  
    ##   Mean Type 1 (sd Type I) 
    ##       Mean sum ( hat beta j - beta j True )^2 (sd sum)
    ##   Computer Time 
###############################################################################
##   MySaveFileName ()
##
##    Based upon characteristics of table, picks a title for Latex file to save
##    
##
##
##    
    MySaveFileName <- function(OneVV, KPAm, NCount, PrMeVec, LL = FALSE) {
       STD <- LoadSavedTableDirectory();
       if (LL== TRUE) {
         My = "L" } else { My = "" } ;
       name <- paste(STD,"/","OutputTable", My, "KP", KPAm, "CNT", NCount,
                     "TB", paste(PrMeVec, collapse=""), 
                     "mNN", tSeq(min(OneVV[,4])), "MNN", tSeq(max(OneVV[,4])),
                     "mKP", tSeq(min(OneVV[,5])), "MKP", tSeq(max(OneVV[,5])),
                     "msig", tSeq(min(OneVV[,6])), "Msig", tSeq(max(OneVV[,6])),
                     ".tex", sep="");
       return(name);
    }
    
#############################################################################
##   DoAllTheSaving <- function(GFMAA, OneVV, KPAm, PrMeVec, NCount)
##       
##     Saves the table created by these function.  Call this function with
##       GFMAA: simulation out put, OneVV matrix of columns requested
##       KPAm is a statement of what the size of active set was before doing study
##       PrMeVec: which of the 10 types of simulation estimators to use
##       NCount:  How many N was the sample size per parameter set.
    DoAllTheSaving <- function(GFMAA, OneVV, KPAm, PrMeVec, NCount, rndit=2) {
       OnMyFileName <- MySaveFileName(OneVV, KPAm, NCount, PrMeVec);
       TableGot <- CreatePrintTable(GFMAA, OneVV, KPAm, PrMeVec);
       BiggerTable <- matrix(0, length(TableGot$MyPrintTB[,1]) +1,
                                length(TableGot$MyPrintTB[1,]) + 2);
       BiggerTable[1, 3:length(BiggerTable[1,])] <- EstimatorColNames2[PrMeVec];   
       BiggerTable[2:length(BiggerTable[,1]),1] <- TableGot$RowsNames;
       BiggerTable[2:length(BiggerTable[,1]),2] <- rep( FunctionPlot, length(TableGot$MyPrintTB[,1]));
       BiggerTable[2:length(BiggerTable[,1]), 3:length(BiggerTable[1,]) ] <- TableGot$MyPrintTB;
       BiggerTable[1,1] = ""; BiggerTable[1,2] = "";
       if (KPAm == 6) {
          BiggerTable[1,1] = "$ \\beta_{\\mbox{\\tiny{1:6}}} =   \\begin{array}{c} (   1,-1, 1,  \\\\   -1, 1, -1  ) \\end{array}$";
            TDCaption <- paste("$ \\beta_{\\mbox{\\tiny{1:6}}}", 
                                  " = \\left( 1,-1,1,-1,1,-1 \right)$",
                                " and $\\sigma = ", round(max(OneVV[,6]),rndit),"$", sep="");
       } else if (KPAm == 4) {
          BiggerTable[1,1] = "$ \\beta_{\\mbox{\\tiny{1:4}}} = \\left( 4,3,-2.5,1 \\right)$$\\mbox{ }$";      
            TDCaption <- paste("$ \\beta_{\\mbox{\\tiny{1:4}}}", 
                                  " = \\left( 4,3,-2.5,1 \right)$",
                                " and $\\sigma = ", round(max(OneVV[,6]),rndit), "$", sep="");          
       }
       ArrayColumns <- paste("{@{\\extracolsep{-1.25mm}}|c@{\\hspace{-.5mm}}|c@{\\hspace{-.5mm}}|",
          paste(rep( "@{\\hspace{-.5mm}}|c", length(BiggerTable[1,])-1), collapse=""),
          "@{\\hspace{-.5mm}}|}", sep="");
       StartF <- paste(" \\begin{tabular} ", ArrayColumns,  " \\hline ", sep="");
       MyF <- file(OnMyFileName, open="wt", blocking=FALSE );
        writeLines(StartF, con=MyF);
        close(MyF);
        write.table(x=BiggerTable, file=OnMyFileName, append=TRUE,
            sep = " & \n", eol=" \\\\ \\hline \\hline \n", na="NA", quote=FALSE,
            row.names=FALSE, col.names=FALSE,);
        ##  open(MyF, "at");
        MyF <- file(OnMyFileName, open="at", blocking=FALSE );
        writeLines(" \\end{tabular} \n", con=MyF);
        close(MyF);
    } 
    
#############################################################################
##   DoAllTheSavingL <- function(GFMAA, OneVV, KPAm, PrMeVec, NCount)
##       
##   Same as DoAllTheSaving() except, this makes a Latex "longtable"
##     Saves the table created by these function.  Call this function with
##       GFMAA: simulation out put, OneVV matrix of columns requested
##       KPAm is a statement of what the size of active set was before doing study
##       PrMeVec: which of the 10 types of simulation estimators to use
##       NCount:  How many N was the sample size per parameter set.    
    DoAllTheSavingL <- function(GFMAA, OneVV, KPAm, PrMeVec, NCount, rndit=2) {
       OnMyFileName <- MySaveFileName(OneVV, KPAm, NCount, PrMeVec, LL = TRUE);
       TableGot <- CreatePrintTable(GFMAA, OneVV, KPAm, PrMeVec);
       BiggerTable <- matrix(0, length(TableGot$MyPrintTB[,1]) +1,
                                length(TableGot$MyPrintTB[1,]) + 2);
       BiggerTable[1, 3:length(BiggerTable[1,])] <- EstimatorColNames2[PrMeVec];   
       BiggerTable[2:length(BiggerTable[,1]),1] <- TableGot$RowsNames;
       BiggerTable[2:length(BiggerTable[,1]),2] <- rep( FunctionPlot, length(TableGot$MyPrintTB[,1]));
       BiggerTable[2:length(BiggerTable[,1]), 3:length(BiggerTable[1,]) ] <- TableGot$MyPrintTB;
       BiggerTable[1,1] = ""; BiggerTable[1,2] = "";
       if (KPAm == 6) {
          BiggerTable[1,1] = "$ \\beta_{\\mbox{\\tiny{1:6}}} =   \\begin{array}{c} (   1,-1, 1,  \\\\   -1, 1, -1  ) \\end{array}$";
            TDCaption <- paste("$ \\beta_{\\mbox{\\tiny{1:6}}}", 
                                  " = \\left( 1,-1,1,-1,1,-1 \\right)$",
                                " and $\\sigma = ", round(max(OneVV[,6]),rndit), "$", sep="");
       } else if (KPAm == 4) {
          BiggerTable[1,1] = "$ \\beta_{\\mbox{\\tiny{1:4}}} = \\left( 4,3,-2.5,1 \\right)$";      
            TDCaption <- paste("$ \\beta_{\\mbox{\\tiny{1:4}}}", 
                                  " = \\left( 4,3,-2.5,1 \\right)$$\\mbox{ }$",
                                " and $\\sigma = ", round(max(OneVV[,6]),rndit), "$", sep="");          
       }
       ArrayColumns <- paste("{@{\\extracolsep{-1.25mm}}|c@{\\hspace{-.5mm}}|c@{\\hspace{-.5mm}}|",
          paste(rep( "@{\\hspace{-.5mm}}|c", length(BiggerTable[1,])-1), collapse=""),
          "@{\\hspace{-.5mm}}|}", sep="");
       StartF <- paste(" \\begin{longtable} ", ArrayColumns,  " \\hline ", sep="");
       MyF <- file(OnMyFileName, open="wt", blocking=FALSE );
        writeLines(StartF, con=MyF);
        close(MyF);
        write.table(x=BiggerTable, file=OnMyFileName, append=TRUE,
            sep = " & \n", eol=" \\\\ \\hline \\hline \n", na="NA", quote=FALSE,
            row.names=FALSE, col.names=FALSE,);
        ##  open(MyF, "at");
        MyF <- file(OnMyFileName, open="at", blocking=FALSE );
        writeLines(paste(" \\caption{", TDCaption, "}", sep=""), con=MyF);
           tabnameS <- unlist(strsplit(OnMyFileName, "/"));
           tabnameS <- tabnameS[length(tabnameS)];
           tabnameS <- unlist(strsplit(tabnameS, "\\."));
           tabnameS <- tabnameS[1];
        writeLines(paste(" \\label{tabl:", tabnameS, "}", sep=""), con=MyF); 
        writeLines(" \\end{longtable} \n", con=MyF);
        close(MyF);
    }   
    
#############################################################################
##   CreatePrintTable <- function(GFMAA, OneVV, KPAm, PrMeVec) 
##       
##     Helper function to DoAllTheSaving, gets the numbers for the table
    CreatePrintTable <- function(GFMAA, OneVV, KPAm, PrMeVec) {
        PrintTB <- matrix(0, length(OneVV[,1]), length(PrMeVec)) 
        for (cto in 1:length(OneVV[,1])) {
           PrintTB[cto,] <- SubTable(GFMAA, OneVV, PrMeVec, cto, rndit = 2);
        }
        MyPrintTB <- PrintTB;
        for (ii in 1:length(PrintTB[,1])) { 
          for (jj in 1:length(PrintTB[1,])) {
            MyPrintTB[ii,jj] <- paste(" $ ", PrintTB[ii,jj], " $ ", sep="");
          }
        }
        RowsNames <- paste( " $ ", SubRows(OneVV, KPAm, rndit = 2), " $ ", sep="");
        ColsNames <- EstimatorColNames2[PrMeVec];
        RetMakeMe <- list(MyPrintTB = MyPrintTB, RowsNames=RowsNames, ColsNames = ColsNames);
        return(RetMakeMe);
    }
 
#############################################################################
##   SubRows <- function(OneVV, KPAm, rndit = 2)
##       
##     Helper function to CreatePrintTable, creates Latex string explaining 
##      characteristics of individual sim.      
    SubRows <- function(OneVV, KPAm, rndit = 2) {
      rt <- paste(" \\begin{array}{c}  P_{\\mbox{\\tiny{xcor}}} = ", 
             ".", unlist(strsplit(as.character(round(OneVV[,2], rndit)), "\\."))[2],
                 " \\mbox{ , } \\xi = ", 
             ".", unlist(strsplit(as.character(round(OneVV[,3], rndit)), "\\."))[2],                 
                 " \\\\",
                 " \\kappa_{\\mbox{\\tiny{$\\mathcal{A}$}}} = ", KPAm, 
                 " \\mbox{ , } \\sigma = ", OneVV[, 6], "\\\\",
                 " n = ", OneVV[,4], "\\mbox{ , } ", 
                 " \\kappa = ", OneVV[,5], "\\end{array} ", sep="");
      return(rt);
    }
#############################################################################
##  SubTable <- function(GFMAA, OneVV, PrMeVec, cto, rndit =2, TMM = FALSE) 
##       
##     Helper function to CreatePrintTable, creates rows for Latex table 
##       
    SubTable <- function(GFMAA, OneVV, PrMeVec, cto, rndit =2, TMM = FALSE) {    
         MeanTII  = PrMeVec *0; SdTII = PrMeVec * 0;
         MeanTI = PrMeVec * 0; SdTI = PrMeVec * 0;
         MeanBetaSq = PrMeVec * 0; SdBetaSq = PrMeVec * 0;
         PercentPer = PrMeVec * 0;
         CompTU  = PrMeVec * 0;  SdCompTU = PrMeVec * 0;
         CompTC = PrMeVec * 0; SdCompTC = PrMeVec * 0;

            ALTA <- 10 * 11 / 2;  ALTB <- 10;
            SubSPlot <- GFMAA[ GFMAA[,1] == OneVV[cto,1] & GFMAA[,2] == OneVV[cto,2] & 
                   GFMAA[,3] == OneVV[cto,3] & GFMAA[,4] == OneVV[cto,4] & 
                   GFMAA[,5] == OneVV[cto,5] & GFMAA[,6] == OneVV[cto,6], ];
            if (length(SubSPlot) == 0) {
               print("SubTable, cannot get any for OneVV = ");
               print(OneVV[cto,]);
               return(0);
            }
            for (tt in 1:length(PrMeVec)) {
                if (TMM == FALSE) {
                   PrV1 <-SubSPlot[, 7 + PrMeVec[tt]];
                   PrV1 <- PrV1[!is.na(PrV1) & PrV1 >= 0 ];
                   PrV2 <-SubSPlot[, 7 + ALTA + PrMeVec[tt]];
                   PrV2 <- PrV2[!is.na(PrV2) & PrV2 >= 0 ];   
                   PrV3 <- SubSPlot[, 7 + ALTA*2 + PrMeVec[tt]]; 
                   PrV3 <- PrV3[!is.na(PrV3) & PrV3 >= 0 ];  
                   PrV4 <- SubSPlot[, 7 + ALTA*3 + PrMeVec[tt]]; 
                   PrV4 <- PrV4[!is.na(PrV4) & PrV4 >= 0 ];   
                   PrV5 <- SubSPlot[, 7 + ALTA*3 + ALTB*2 + PrMeVec[tt]]; 
                   PrV5 <- PrV5[!is.na(PrV5) & PrV5 >= 0 ];                                                   
                } else {
                   PrV1 <-SubSPlot[, 7 + PrMeVec[tt]];
                   PrV1[is.na(PrV1) | PrV1 < 0 ] <- max(PrV1[!is.na(PrV1) & PrV1 >= 0 ]);
                   PrV2 <-SubSPlot[, 7 + ALTA + PrMeVec[tt]];
                   PrV2[is.na(PrV2) | PrV2 < 0]  <- max(PrV1[!is.na(PrV2) & PrV2 >= 0 ]);  
                   PrV3 <- SubSPlot[, 7 + ALTA*2 + PrMeVec[tt]]; 
                   PrV3[is.na(PrV3) | PrV3 < 0] <- max(PrV3[!is.na(PrV3) & PrV3 >= 0 ]);   
                   PrV4 <- SubSPlot[, 7 + ALTA*3 + PrMeVec[tt]]; 
                   PrV4[is.na(PrV4) | PrV4 < 0] <- max(PrV4[!is.na(PrV4) & PrV4 >= 0 ]); 
                   PrV5 <- SubSPlot[, 7 + ALTA*3 + ALTB*2 + PrMeVec[tt]]; 
                   PrV5[is.na(PrV5) | PrV5 < 0] <- max(PrV5[!is.na(PrV5) & PrV5 >= 0 ]);                                                                                            
                }
                MeanTII[tt] <- mean(PrV1);  SdTII[tt] <- sd(as.vector(PrV1));
                MeanTI[tt] <- mean(PrV2); SdTI[tt] <- sd(as.vector(PrV2));
                MeanBetaSq[tt] <- mean(PrV3); SdBetaSq[tt] <- sd(as.vector(PrV3));
                PercentPer[tt] <- length( PrV1[PrV1 ==0 & PrV2 == 0] ) /
                                  length(SubSPlot[,1]);
                CompTU[tt] <- mean(PrV4); CompTC[tt] <- mean(PrV5);
                SdCompTU[tt] <- sd(as.vector(PrV4)); 
                SdCompTC[tt] <- sd(as.vector(PrV5));                                  
            }
            rt <- WhatGoesEachBoxAA(tt=0, MeanTII, SdTII, MeanTI, SdTI, MeanBetaSq, SdBetaSq,
                PercentPer, CompTU, SdCompTU, CompTC, SdCompTC, rndit);
            return(rt);
     }
#############################################################################
##  WhatGoesEachBoxAA <- function(tt =0, MeanTII, SdTII,...
##       
##     Helper function to CreatePrintTable, use to input all sumary 
##     statistics one desires for simulation
##        
     WhatGoesEachBoxAA <- function(tt =0, MeanTII, SdTII, MeanTI, SdTI, MeanBetaSq, SdBetaSq,
            PercentPer, CompTU, SdCompTU, CompTC, SdCompTC, rndit=2) {
          SdComptTU = SdCompTU;
       if (tt <= 0) {
     ##   rt <- paste(" \\begin{array}{c} \\hline ",
      rt <- paste(" \\begin{array}{c} ",
           rd0(round(MeanTII,rndit)), "\\mbox{ (", 
                   rd0(round(SdTII,rndit)), ")} \\\\ \\hline ",
        rd0(round(MeanTI, rndit)), "\\mbox{ (", rd0(round(SdTI,rndit)), ") } \\\\ \\hline",
        rd0(round(MeanBetaSq, rndit)), "\\mbox{ (", rd0(round(SdBetaSq, rndit)), ") } \\\\ \\hline ",
        ##round(PercentPer, rndit), "\\mbox{\\%} \\\\ \\hline ",
        rd0(round(CompTC, rndit)), " \\\\ ",
        ##round(CompTC, rndit), "\\mbox{s (", round(SdComptTU, rndit), ") } \\\\",
        " \\end{array} ", sep="");
       } else {
       ##  rt <- paste(" \\begin{array}{|c|c|} \\hline ",
       rt <- paste(" \\begin{array}{c} ",
           rd0(round(MeanTII[tt],rndit)), "\\mbox{ (", 
                   rd0(round(SdTII[tt],rndit)), ")} \\\\ \\hline ",
        rd0(round(MeanTI[tt], rndit)), "\\mbox{ (", rd0(round(SdTI[tt],rndit)), ") } \\\\ \\hline",
        rd0(round(MeanBetaSq[tt], rndit)), "\\mbox{ (", rd0(round(SdBetaSq[tt], rndit)), ") } \\\\ \\hline",
        ##round(PercentPer[tt]*100, rndit), "\\mbox{\\%} \\\\ \\hline ",
        rd0(round(CompTC, rndit)), " \\\\ ",
        ##rd0(round(CompTC[tt], rndit)), "\\mbox{s (", rd0(round(SdComptTU[tt], rndit)), ") } \\\\",
        " \\end{array} ", sep="");      
       }
       return(rt)
    }
#############################################################################
##  WhatGoesEachBoxBB <- function(tt =0, MeanTII, SdTII, MeanTI,,...
##       
##     Helper function to CreatePrintTable, use to input all sumary 
##     statistics one desires for simulation
##      
    WhatGoesEachBoxBB <- function(tt =0, MeanTII, SdTII, MeanTI, SdTI, MeanBetaSq, SdBetaSq,
            PercentPer, CompTU, SdCompTU, CompTC, SdCompTC, rndit=2) {
       if (tt <= 0) {
     ##   rt <- paste(" \\begin{array}{|c|c|} \\hline ",
      rt <- paste(" \\begin{array}{c|c} \\hline ",
           rd0(round(MeanTII,rndit)), "\\mbox{ (", 
                   rd0(round(SdTII,rndit)), ")} & ",
        rd0(round(MeanTI, rndit)), "\\mbox{ (", rd0(round(SdTI,rndit)), ") } \\\\ \\hline",
        rd0(round(MeanBetaSq, rndit)), "\\mbox{ (", rd0(round(SdBetaSq, rndit)), ") } & ",
        rd0(round(PercentPer, rndit)), "\\mbox{\\%} \\\\ \\hline ",
        rd0(round(CompTU, rndit)), "\\mbox{sc (", rd0(round(SdCompTU, rndit)), ") } & ",
        rd0(round(CompTC, rndit)), "\\mbox{sc (", rd0(round(SdCompTU, rndit)), ") } \\\\ \\hline",
        " \\end{array} ", sep="");
       } else {
       ##  rt <- paste(" \\begin{array}{|c|c|} \\hline ",
       rt <- paste(" \\begin{array}{c|c} \\hline ",
           rd0(round(MeanTII[tt],rndit)), "\\mbox{ (", 
                   rd0(round(SdTII[tt],rndit)), ")} & ",
        rd0(round(MeanTI[tt], rndit)), "\\mbox{ (", rd0(round(SdTI[tt],rndit)), ") } \\\\ \\hline",
        rd0(round(MeanBetaSq[tt], rndit)), "\\mbox{ (", rd0(round(SdBetaSq[tt], rndit)), ") } & ",
        rd0(round(PercentPer[tt]*100, rndit)), "\\mbox{\\%} \\\\ \\hline ",
        rd0(round(CompTU[tt], rndit)), "\\mbox{sc (", rd0(round(SdCompTU[tt], rndit)), ") } & ",
        rd0(round(CompTC[tt], rndit)), "\\mbox{sc (", rd0(round(SdCompTU[tt], rndit)), ") } \\\\ \\hline",
        " \\end{array} ", sep="");      
       }
       return(rt)
    }
#############################################################################
##  LoadSavedTableDirectory <- function()
##       
##     Tries to identify directory to save Latex Tables into.

##      
   LoadSavedTableDirectory <- function() {
   FilesInDIR <- unlist(list.files(.Library));
     if (any(FilesInDIR == "PrintTables")) {
         PathMeR <- MakePathMe();  
         FilesInDIR <- unlist(list.files(PathMeR));
         if (any(FilesInDIR == "PrintTables")) {
           PathMeR = paste(PathMeR, "PrintTables", sep="");
         } else {
           PathMeR = paste(PathMeR, "PrintTables", sep="");
           dir.create(PathMeR, showWarnings = FALSE, recursive = FALSE);
         }
     } else {
         PathMeR = "c://Stat//2008Summer//LarsProject//code//PrintTables//"
     }
     SavedOutPutDirectory = PathMeR;
     return(SavedOutPutDirectory);
}
#######################################################################################
#######################################################################################