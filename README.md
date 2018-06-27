# TetradAnalysis
# R functions for Confirmatory Tetrad Anlysis.

"TetradAnalysisNoRandom" function allows you to run Confirmatory Tetrad Analysis from Bollen 1993.

"NestedOrNot" function allows you to test two models are nested in terms of tetrads or not.

# Instructions to use this package:
Release it in a folder, e.g. "D:/TetradAnalysisRpackage"
The use R codes below to install package:

install.packages("devtools")

library("devtools")

devtools::install_github("klutometis/roxygen")

library(roxygen2)

setwd("D:/TetradAnalysisRpackage/ConfirmatoryTetradAnalysis")

setwd("..")

install("ConfirmatoryTetradAnalysis")

install.packages("lavaan")

library("lavaan")

TetradAnalysisNoRandom

NestedOrNot

?TetradAnalysisNoRandom

?NestedOrNot
