#' NestedOrNot Function
#'
#' "NestedOrNot " is used to know two models are nested in terms of tetrads or not (Bollen, 1993).
#' @param "Data1", An optional data frame containing the observed variables used in the first model.
#' @param "Data2", An optional data frame containing the observed variables used in the second model.
#' @param "Model1" A description of the user-specified first model. The model is described using the lavaan model synta(Rosseel, 2012).
#' @param "Model2" A description of the user-specified second model. The model is described using the lavaan model synta(Rosseel, 2012).
#' @param "Size1" Size of the first data sample.
#' @param "Size2" Size of the second data sample.

#' @keywords Nested or not in Confirmatory Tetrad Analysis

#' @return "Model1COV",    Implied variance-covariance matrix from the first model.
#' @return "Model2COV",    Implied variance-covariance matrix from the second model.
#' @return "NRVT1",     Nonredundant vanishing tetrad setting implied by the first model.
#' @return "NRVT2",     Nonredundant vanishing tetrad setting implied by the second model.
#' @return "NRVT_Num1",     Nonredundant vanishing tetrad number for the first model.
#' @return "NRVT_Num1",     Nonredundant vanishing tetrad number for the second model.
#' @return "NEST",     Two models are nested or not nested in terms of tetrads.


#' @references Bollen, K. A., & Ting, K. F. (1993). Confirmatory tetrad analysis. Sociological methodology, 147-175.
#' @references Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. Journal of Statistical Software, 48(2), 1-36.

#' @export
#' @examples
#' ## Input sample variance-covariance matrix "(Sampledata) ".
#' data<-c(2.18, 0.63, 0.50, 0.22, 0.39, 0.42, 0.63, 1.44, 0.50, 0.12, 0.11, 0.15, 0.50, 0.50, 1.71, 0.87, 0.59, 0.48, 0.22, 0.12, 0.87, 1.46, 0.54, 0.42, 
#'	0.39, 0.11, 0.59, 0.54, 1.12, 0.23, 0.42, 0.15, 0.48, 0.42, 0.23, 1.42)
#' Sampledata <- matrix(data, nrow = 6, ncol = 6, byrow = TRUE)
#' colnames(Sampledata) <- c('x1','x2','x3','x4','x5','x6')
#'
#' ## Sample size
#' Size<-100
#'
#' ## Input model
#' FactorModel1<-'
#' Y1 =~ x1 + x2 + x3
#' Y2 =~ x4 + x5 + x6
#' Y2 ~ Y1
#' x3 ~ x4
#' '
#'
#'FactorModel2<-'
#'Y1 =~ x1 + x2 + x3
#'Y2 =~ x4 + x5 + x6
#'Y2 ~ Y1
#'x4 ~~ x3
#''
#'
#' ## Use "NestedOrNot " to run Confirmatory Tetrad Analysis
#' NestedOrNot(Sampledata,FactorModel1,Size,Sampledata,FactorModel2,Size)



##################################################################################
################################## Nested or not Function #############################
##################################################################################
### Information needed  ###
### 1. Two sample covariance matrices "Samplecov1" and "Samplecov2". ###
### 2. Two SEM models "FactorModel1" and "FactorModel2". ###
### 3. Two sample size "Size1" and "Size2". ###

NestedOrNot<-function(Samplecov1,FactorModel1,Size1,Samplecov2,FactorModel2,Size2){
	Model1<-TetradAnalysisNoRandom(Samplecov1,FactorModel1,Size1)
	Model2<-TetradAnalysisNoRandom(Samplecov2,FactorModel2,Size2)
	EmpiricalTetrads1<-Model1$EmpiricalTetrads
	EmpiricalTetrads2<-Model2$EmpiricalTetrads

	if (nrow(EmpiricalTetrads1)<=nrow(EmpiricalTetrads2)){		# Check which model has more vanishing tetrads after empirical method. #
	if (all(EmpiricalTetrads1[,1] %in% EmpiricalTetrads2[,1])==TRUE){	# If all vanishing tetrads for one model are all exist in another model, then they are nested. # 
		nest<-"Two models are nested"}else{
		nest<-"Two models are not nested"}
	}else{
		if (all(EmpiricalTetrads2[,1] %in% EmpiricalTetrads1[,1])==TRUE){
		nest<-"Two models are nested"}else{
		nest<-"Two models are not nested"}
	}
	results<-list(Model1COV=Model1$Modelcov,
			  Model2COV=Model2$Modelcov,
			  NRVT1=Model1$NRVT,
			  NRVT2=Model2$NRVT,
			  NRVT_Num1=Model1$NRVT_Num,		
			  NRVT_Num2=Model2$NRVT_Num,
			  NEST=nest)
### Output nonredundant tetrad setting, nonredundant tetrad number and nested or not status for these two models. ###

	return(results)
}
########### NestedOrNot Function End#######################################################
