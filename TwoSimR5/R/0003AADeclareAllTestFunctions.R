###############################################################################
##  0003AADeclareAllTestFunctions.r  -- Alan Lenarcic
##
##    01-09-2019
##
##   DeclareAllTestFunctions
##
##    Helps install Test functions into the package with setting curried
##  into the answer.
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

try(library(lars, warn.conflicts = FALSE, quietly=TRUE), silent=TRUE);
try(library(corpcor, warn.conflicts = FALSE, quietly=TRUE), silent=TRUE);
try(library(TwoLassoCpp, warn.conflicts = FALSE, quietly=TRUE), silent=TRUE);
try(library(methods, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
try(library(R.methodsS3, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
try(library(R.oo, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);
try(library(Rcpp, warn.conflicts=FALSE, quietly=TRUE), silent=TRUE);

DeclareAllDefaultFunctions <- function() {
  DefaultAllUseFunctions <- c(GenerateLarsFixed, GenerateLarsCp, GenerateLassoStoddenR,
  GenerateLimitLassoR, GeneratepiALimitLassoR, 
  GenerateFDLimitLassoR,  GenerateXLMMedianR,                 
  GenerateXLLimitLassoR, GenerateIgnorantLimitLassoR,
  GenerateLassoQLYR, GenerateTwoLasso9XR, GenerateEMRidgeR,
  GenerateTwoLassoSpreadR, GenerateSCADCVMR,
  GeneratepiAM2LassoR,GenerateFDM2LassoR,
  GenerateXLM2LassoR,GenerateSCADFourVR, GenerateSCADMinR,
  GenerateENetFixedR, GenerateENetCVMinR,
  GenerateHorseShoeR, GenerateSpikeAndSlabR, GenerateB2LassoR, GenerateBayesLassoR,
  Generate2LassoCVR, GenerateISpike,
  GenerateBayesSpikeInfoR, GenerateBayesSpikeAutoR,
  GenerateLarsFixedNoiseR, GenerateLarsCpNoiseR, GenerateLassoStoddenNoiseR,
  GenerateLimitLassoNoiseR, GeneratepiALimitLassoNoiseR, 
  GenerateFDLimitLassoNoiseR,  GenerateXLMMedianNoiseR,                 
  GenerateXLLimitLassoNoiseR,
  GenerateLassoQLYNoiseR, GenerateTwoLasso9XNoiseR, GenerateEMRidgeNoiseR,
  GenerateTwoLassoSpreadNoiseR,
  GeneratepiAM2LassoNoiseR,GenerateFDM2LassoNoiseR,
  GenerateXLM2LassoNoiseR,  GenerateSCADFourVR,
  GenerateENetFixedNoiseR, 
  GenerateBayesSpikeInfoNoiseR, GenerateBayesSpikeAutoMeanR,
  GenerateBayesSpikeInfoMeanR, GenerateBayesSpikeInfoMeanNoiseR,
  GenerateMCPCVMR, GenerateMCPMinR,
  GenerateBayesSpikeAutoLogitSlowR,
  GenerateBayesSpikeInfoLogitSlowR,
  GenerateBayesVarSelR,
  GenerateHorseShoeOldR,
  GenerateBayesSpikeAutoSQRTR,
  GenerateBayesSpikeAutoSQRTMeanR,
  GenerateBayesSpikeAutoTemp1R,
  GenerateBayesSpikeAutoTemp2R,
  GenerateBayesSpikeAutoTemp3R,
  GenerateBayesSpikeAutoTemp4R,
  GenerateBayesSpikeAutoTemp2c10R,
  GenerateBayesSpikeAutoTemp2c05R,
  GenerateBayesSpikeAutoTemp2c20R,
  GenerateBayesSpikeAutoTemp3c05R,
  GenerateBayesSpikeAutoTemp3c10R,
  GenerateBayesSpikeAutoTemp3c20R
  );                     
DefaultAllNameFunctions <- c("LarsFixed", "LarsCp", "LassoStodden",
    "LimitLasso", "piALimitLasso", "FDLimitLasso",
    "XLMMedian", "XLLimitLasso", "IgnorantLimitLasso", "LassoLY", "TwoLasso9X", 
    "LimitRidge", "2LassoSpread", "SCADminConvex", "piAM2Lasso", "FDM2Lasso", "XLM2Lasso",
    "SCAD4xStodden", "SCADMinLRY", "ENetFixed", "ENetCV",
    "HorseShoe",
    "OurFirstSpike", "B2Lasso", "BayesLasso", "2LassoCV", "IshwaranSpike",
    "GenerateBayesSpikeInfoPrior", "GenerateBayesSpikeAutoPrior",
    "LarsFixedNI", "LarsCpNI", "LassoStoddenNI",
    "LimitLassoNI", "piALimitLassoNI", "FDLimitLassoNI",
    "XLMMedianNI", "XLLimitLassoNI",
    "LassoLYNI", "TwoLasso9XNI", 
    "LimitRidgeNI", "2LassoSpreadNI",
    "piAM2LassoNI", "FDM2LassoNI", "XLM2LassoNI", "SCAD4xStoddenNI", 
    "ENetFixedNI",
    "GenerateBayesSpikeInfoPriorNI", "GenerateBayesSpikeAutoMean",
    "GenerateBayesSpikeInfoMean",     "GenerateBayesSpikeInfoMeanNI",
    "MCPCVMR", "MCPMinR",
    "GenerateBayesSpikeAutoLogitSlow", "GenerateBayesSpikeInfoLogitSlow",
    "BayesVarSel", "HorseShoeAL",
  "GenerateBayesSpikeAutoSQRTR",
  "GenerateBayesSpikeAutoSQRTMeanR",
  "GenerateBayesSpikeAutoTemp1R",
  "GenerateBayesSpikeAutoTemp2R",
  "GenerateBayesSpikeAutoTemp3R",
  "GenerateBayesSpikeAutoTemp4R",
  "GenerateBayesSpikeAutoTemp2c10R",
  "GenerateBayesSpikeAutoTemp2c20R",
  "GenerateBayesSpikeAutoTemp2c05R",
  "GenerateBayesSpikeAutoTemp3c10R",
  "GenerateBayesSpikeAutoTemp3c20R",
  "GenerateBayesSpikeAutoTemp3c05R" );   


eval(parse(text=GetG0Text("TWOSIMENVIRONMENT", S=1)));
eval(parse(text=GetG0Text("DefaultAllGroupFunctions", S=1))); 
DefaultAllGroupFunctions <- c(
  GenerateGroupBayesSpikeAutoR, 
  GenerateGroupBayesSpikeInfoR, GenerateGroupBayesSpikeInfoNoiseR,
  GenerateGroupBayesSpikeAutoMeanR,GenerateGroupBayesSpikeInfoMeanR,
  GenerateGroupBayesSpikeInfoMeanNoiseR,
  GenerateGroupTwoLassoSpread,
  GenerateGroupLassoReg,
  GenerateGroupGrpReg, GenerateStandGLReg, GenerateSGLReg,
  GenerateGroupTwoLassoSpreadNI,
  GenerateSGLRegNI,
  GenerateGroupBayesSpikeAutoLogitSlowR, GenerateGroupBayesSpikeInfoLogitSlowR
);
eval(parse(text=SetGText("DefaultAllGroupFunctions", S=1)));
try(eval(parse(text=SetGText("DefaultAllGroupFunctions", 
  envir="TWOSIMENVIRONMENT", S=1))), silent=TRUE); 
DefaultAllGroupNameFunctions <- c(
  "GroupBayesSpikeAuto", "GroupBayesSpikeInfo",
  "GroupBayesSpikeSpikeInfoNoise", "GroupBayesSpikeAutoMean",
  "GroupBayesSpikeInfoMean", "GroupBayesSpikeInfoMeanNoise",
  "GenerateGroupTwoLassoSpread", 
  "GenerateGroupLassoReg", "GenerateGroupGrpReg",
  "GenerateStandGLReg", "GenerateSGLReg",
  "GenerateGroupTwoLassoSpreadNI",
  "GenerateSGLRegNI",
  "GenerateGroupBayesSpikeAutoLogitSlowR", "GenerateGroupBayesSpikeInfoLogitSlowR"
);
eval(parse(text=SetGText("DefaultAllGroupNameFunctions", S=1)));
try(eval(parse(text=SetGText("DefaultAllGroupNameFunctions", 
  envir="TWOSIMENVIRONMENT", S=1))),silent=TRUE);  
eval(parse(text=GetG0Text("DefaultPriorStrengthFunction", "globalenv()", S=1)));
DefaultPriorStrengthFunction <- function(n=100, p=100, pGroups = NULL, ...) {
  if (!is.null(pGroups) && is.numeric(pGroups) && pGroups[1] >= 0) {
    return(pGroups);
  }
  return(p);
}
eval(parse(text=SetGText("DefaultPriorStrengthFunction", "globalenv()", S=1)));
eval(parse(text=GetG0Text("DefaultSigmaPriorStrengthFunction", "globalenv()", S=1)));
DefaultSigmaPriorStrengthFunction <- function(n=100, p=100, pGroups = NULL, ...) {
  return(n);
}
eval(parse(text=SetGText("DefaultSigmaPriorStrengthFunction", "globalenv()", S=1)));
eval(parse(text=SetGText("DefaultAllGroupFunctions", "globalenv()", S=1)));
eval(parse(text=SetGText("DefaultAllUseFunctions", "globalenv()", S=1)));
eval(parse(text=SetGText("DefaultAllGroupNameFunctions", "globalenv()", S=1)));
eval(parse(text=SetGText("DefaultAllNameFunctions", "globalenv()", S=1)));
    
}
DeclareAllTestFunctions <- function(SubSetInt = 0, SubSetNameFunctions = 0) {


 DeclareAllDefaultFunctions();               



eval(parse(text=GetG0Text("DefaultAllGroupNameFunctions", "globalenv()", S=1)));
eval(parse(text=GetG0Text("DefaultAllNameFunctions", "globalenv()", S=1)));
eval(parse(text=GetG0Text("DefaultAllGroupFunctions", "globalenv()", S=1)));
eval(parse(text=GetG0Text("DefaultAllUseFunctions", "globalenv()", S=1)));
AllNameFunctions <- DefaultAllNameFunctions;
AllGroupNameFunctions <- DefaultAllGroupNameFunctions;
AllGroupFunctions <- DefaultAllGroupFunctions;
AllUseFunctions <- DefaultAllUseFunctions;

if (length(SubSetInt) > 1 || SubSetInt != 0) {
  
} else if (length(SubSetNameFunctions) > 1 || SubSetNameFunctions != 0) {
  SubSetInt = match(SubSetNameFunctions, AllNameFunctions);
  SubSetInt=SubSetInt[!is.na(SubSetInt)];
} else {
  SubSetInt = (1:length(AllNameFunctions)); 
}
NameFunctions = AllNameFunctions[SubSetInt];
UseFunctions = AllUseFunctions[SubSetInt];
eval(parse(text=SetGText("UseFunctions",S=1)));
eval(parse(text=SetGText("NameFunctions",S=1)));   


eval(parse(text=GetG0Text("OnSimType", "globalenv()", S=1)));
if (!is.null(OnSimType) && is.character(OnSimType) &&
  OnSimType %in% c("Group", "GroupRobit", "GroupLogit")) {
   NameFunctions = DefaultAllGroupFunctions;
   UseFunctions = DefaultAllGroupNameFunctions;
   eval(parse(text=SetGText("UseFunctions",S=1)));
   eval(parse(text=SetGText("NameFunctions",S=1)));   
}    
                  
PnV4 <- c("Lasso-\nFixed", "Lars-\nCp", "Lasso-\nStodden",
    "Limit-\nLasso", "Limit-2PI\n-Lasso", "Lim-\nLasso FD",
    "Marg\n Median", "MargM Lim\n-Lasso", "Ignorant \n Limit-Lasso", 
     "Lasso \nw=1 LY", "2-Lasso\n9X", "Limit-Ridge", "M2Lasso", "SCAD-minConvex", 
     "piAM2Lasso", "FDM2Lasso", "XLM2Lasso",
    "SCAD-4xStodden", "SCAD-MinWLT", "ENet-Fixed", "ENet-CVmin",
    "R2Lasso", "R2PiLasso", "R2FDLasso", "R2XLMLasso", "HorseShoe", "Our First Spike",
    "B2Lasso", "Bayes Lasso", "Ishwaran Spike", "BayesSpikeE", "BayesSpikeA", "MCP-CV", "MCP Min R")[SubSetInt];
eval(parse(text=SetGText("PnV4", S=1)));
EstimatorNames12 <- c("Lasso Fixed", "Lars Cp", "Lasso-Stodden",
  "Limit-Lasso", "Limit-2PILasso", "FDLimit-Lasso",
  "Marginal Median", "MM Limit-Lasso", "Ignorant LimLasso",
  "Lasso Lin and Yuan", "2Lasso 9X",  "Limit-Ridge", 
  "M2Lasso", "SCAD-minConvex", 
  "piA-M2Lasso", "FDM2Lasso", "XLM2Lasso",
  "SCAD-4xStodden", "SCAD-MinWLT","ENet-Fixed", "ENet-CVmin",
  "R2Lasso", "R2PiLasso", "R2FDLasso", "R2XLMLasso", 
  "HorseShoe", "Our First Spike",
  "B2Lasso", "Bayes Lasso", "Ishwaran Spike", "ByesSpikeE", "BayesSpikeA", "MCP-CV", "MCP Min R")[SubSetInt];
LeftNames12 <- c("Lars-Fixed", "Lars-$C_{p}$", "Lasso-Stodden",
  "L2Lasso", "$\\hatpia$-L2Lasso", "FDLimit-Lasso",
  "Marginal Median", "MM Limit-Lasso", "Ignorant LimLasso",
  "Lasso Lin and Yuan", "2Lasso 9X",  "Limit-Ridge", 
  "M2Lasso", "SCAD-minConvex", 
  "$\\hatpia$-M2Lasso", "FDM2Lasso", "XLM2Lasso",
  "SCAD-4xStodden", "SCAD-MinGIC","ENet-Fixed", "ENet-CVmin",
  "R2Lasso", "R2$\\hatpia$Lasso", "R2FDLasso", "R2XLMLasso", 
  "HorseShoe", "Spike and Slab",
  "B2Lasso", "Bayes Lasso", "Ishwaran Spike", "BayesSpike E-Prior", 
  "BayesSpike A-Prior")[SubSetInt];
eval(parse(text=SetGText("EstimatorNames12",S=1)));
eval(parse(text=SetGText("LeftNames12",S=1)));
  eval(parse(text=GetG0Text("piAPriorData")));
  if (is.null(piAPriorData) || (!is.numeric(piAPriorData) && piAPriorData == "")) {
    piAPriorData = 0; 
  }
  eval(parse(text=GetG0Text("piSV")));
  if (piSV > 0) {
    ka = "\\kappa{\\mbox{\\tiny{F}}}"
    piaKnown = "\\small{W}\\footnotesize{rong}";
    pia =  "\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}"
    piaW = "\\pi_{\\mbox{\\tiny{F}}}";
    ka = "k";
    kap = "p";
  }  else {
    ka = "\\kappa_{\\mbox{\\tiny{$\\mathcal{A}$}}}"
    piaKnown = "\\small{K}\\footnotesize{nown}";
    piaW =  "\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}"
    pia = "\\pi_{\\mbox{\\tiny{$\\mathcal{A}$}}}"
    kap = "p";
    ka = "k"
  }
  eval(parse(text=SetGText("piSV",S=1)));
FunctionPlot12 <- paste(" $ \\begin{array} {c} ", "\\mbox{\\footnotesize{\\# II}} \\\\ \\hline",
                    "\\mbox{\\footnotesize{\\# I}} \\\\ \\hline", 
                    "\\mbox{\\footnotesize{$\\sum \\delta_{\\mbox{\\tiny{$\\beta$}}}^2$}} \\\\ \\hline",
                    ##"\\mbox{\\footnotesize{\\% Perf}} \\\\ \\hline ", 
                    "\\mbox{\\footnotesize{Run}}",
                    "\\end{array} $ ", sep="");
FunctionPlot3 <- c( "$L_{2}$", "Type 1", "Type 2", "\\emph{time}");
eval(parse(text=SetGText("FunctionPlot12",S=1)));
eval(parse(text=SetGText("FunctionPlot3",S=1)));
EstimatorColNames12<-  c(
      paste("\\begin{array}{c} \\mbox{LARS} \\\\",
         " \\mbox{Fixed $", ka, "$",
         " } \\end{array} ", sep=""),
      paste("\\begin{array}{c} \\mbox{LARS} \\\\",
         " \\mbox{$C_{\\mbox{\\tiny{p}}}$",
         "} \\end{array}", sep=""),
      paste("\\begin{array}{c} \\mbox{Lasso} \\\\",
         " \\mbox{Stodden}",
         " \\end{array}", sep=""),         
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{L}\\footnotesize{im}\\small{L}",
         "\\footnotesize{asso}} \\\\",
         " \\mbox{$", pia, "$ ",
         piaKnown,  "}",
         " \\end{array}", sep=""),
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{L}\\footnotesize{im}",
         "\\small{L}\\footnotesize{asso}} \\\\",
         " \\mbox{$", pia, "$",
         " Est. $n\\times$", round(piAPriorData,3), 
         "} \\end{array}", sep=""),
      paste("\\begin{array}{c} \\mbox{Fermi-D} \\\\ ",
         "\\mbox{\\small{L}\\footnotesize{im}\\small{L}",
         "\\footnotesize{asso}} \\end{array}", sep=""),         
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{P}\\footnotesize{sd}-\\small{M}",
         "\\footnotesize{arg}} \\\\ \\mbox{\\small{M}\\footnotesize{edian}} ",
         "\\end{array}", sep=""),
      paste("\\begin{array}{c} \\mbox{\\small{M}\\footnotesize{arg}",
         "\\small{M}\\footnotesize{edian}} \\\\ ",
         "\\mbox{\\small{L}\\footnotesize{im}\\small{L}\\footnotesize{asso}} \\end{array}", sep=""),
      paste("\\begin{array}{c} \\mbox{\\small{T}\\footnotesize{wo}",
         "\\small{I}\\footnotesize{gnorant}}\\\\",
         " \\mbox{Lim-Lasso} \\end{array}", sep=""), 
      paste("\\begin{array}{c} \\mbox{\\small{L}\\footnotesize{asso $w$=1}}",
            " \\\\ \\mbox{\\small{L}\\footnotesize{in \\& }\\small{Y}\\foootnotesize{uan}}} ",
            " \\end{array}", sep=""),           
      paste("\\begin{array}{c} \\mbox{\\small{T}\\footnotesize{wo}",
         "\\small{L}\\footnotesize{asso}}\\\\",
         " \\mbox{$\\times$ 9} \\end{array}", sep=""),   
      paste( "\\begin{array}{c}",
         " \\mbox{\\small{L}\\footnotesize{im}",
         "\\small{R}\\footnotesize{idge}} \\\\",
         " \\mbox{$", pia, "$",
         " ", piaKnown, "} \\end{array}", sep=""), 
      paste( "\\begin{array}{c}",
         " \\mbox{\\small{M}\\footnotesize{2Lasso}} \\\\",
         " \\mbox{$", pia, "$",
         " ", piaKnown, "} \\end{array}", sep=""),
      paste( "\\begin{array}{c}",
         " \\mbox{\\small{SCAD}} \\\\",
         " \\mbox{min Convex}",
         " \\end{array}", sep=""),  
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{M}\\footnotesize{2Lasso}}",
         "\\\\",
         " \\mbox{$", pia, "$",
         " Est. $n\\times$", round(piAPriorData,3), 
         "} \\end{array}", sep=""),  
     paste("\\begin{array}{c} \\mbox{Fermi-D} \\\\ ",
         "\\mbox{\\small{M}\\footnotesize{2Lasso}}",
         " \\end{array}", sep=""),  
     paste("\\begin{array}{c} \\mbox{\\small{M}\\footnotesize{arg}",
         "\\small{M}\\footnotesize{edian}} \\\\ ",
         "\\mbox{\\small{M}\\footnotesize{2Lasso}} \\end{array}", sep=""),
       paste( "\\begin{array}{c}",
         " \\mbox{\\small{SCAD}} \\\\",
         " \\mbox{$\\lambda = 4 \\sigma / \\hat{k}$}",
         " \\end{array}", sep=""), 
           paste( "\\begin{array}{c}",
         " \\mbox{\\small{SCAD}} \\\\",
         " \\mbox{\\small{min W.L.T.}}",
         " \\end{array}", sep=""),
         paste( "\\begin{array}{c}",
         " \\mbox{\\small{ENET}} \\\\",
         " \\mbox{Fixed $", ka, "$",
         " } \\end{array} ", sep=""),
           paste( "\\begin{array}{c}",
         " \\mbox{\\small{ENet}} \\\\",
         " \\mbox{\\small{min CV}}",
         " \\end{array}", sep=""),
            paste( "\\begin{array}{c}",
         " \\mbox{\\small{R}\\footnotesize{2Lasso}} \\\\",
         " \\mbox{$", pia, "$",
         " ", piaKnown, "} \\end{array}", sep=""),
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{R}\\footnotesize{2Lasso}}",
         "\\\\",
         " \\mbox{$", pia, "$",
         " Est. $n\\times$", round(piAPriorData,3), 
         "} \\end{array}", sep=""),  
     paste("\\begin{array}{c} \\mbox{Fermi-D} \\\\ ",
         "\\mbox{\\small{R}\\footnotesize{2Lasso}}",
         " \\end{array}", sep=""),  
     paste("\\begin{array}{c} \\mbox{\\small{M}\\footnotesize{arg}",
         "\\small{M}\\footnotesize{edian}} \\\\ ",
         "\\mbox{\\small{R}\\footnotesize{2Lasso}} \\end{array}", sep=""), 
     paste("\\begin{array}{c} \\mbox{\\small{H}\\footnotesize{orse}} \\\\ ", 
         "\\mbox{\\small{S}\\footnotesize{hoe}} \\end{array}",sep=""),
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{O}\\footnotesize{ur }\\small{S}",
         "\\footnotesize{pike}} \\\\",
         " \\mbox{$", pia, "$ ",
         piaKnown,  "}",
         " \\end{array}", sep=""),
      paste( "\\begin{array}{c}",
         " \\mbox{\\small{B}\\footnotesize{2Lasso}} \\\\",
         " \\mbox{$", pia, "$",
         " ", piaKnown, "} \\end{array}", sep=""),                                                                                                 
      paste( "\\begin{array}{c}",
         " \\mbox{\\small{B}\\footnotesize{ayes}} \\\\",
         " \\mbox{\\small{L}\\footnotesize{asso}} \\end{array}", sep=""),
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{I}\\footnotesize{shwaran }\\small{S}",
         "\\footnotesize{pike}} \\\\",
         " \\mbox{$", pia, "$ ",
         piaKnown,  "}",
         " \\end{array}", sep=""),
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{B}\\footnotesize{ayes}\\small{S}",
         "\\footnotesize{pike}} \\\\",
         " \\mbox{Errored Prior}",
         " \\end{array}", sep=""),
      paste("\\begin{array}{c} ",
         "\\mbox{\\small{B}\\footnotesize{ayes}\\small{S}",
         "\\footnotesize{pike}} \\\\",
         " \\mbox{Artifice Prior}",
         " \\end{array}", sep="")  
   )[SubSetInt];
  eval(parse(text=SetGText("EstimatorColNames12",S=1)));
}
   