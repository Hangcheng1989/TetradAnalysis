#' NestedOrNot Function
#'
#' This function allows you to know two SEM models are nested in terms of tetrads or not.
#' @param "Samplecov1", "Samplecov2", "FactorModel1", "FactorModel2", "Size1" and "Size2" needed for this function.
#' @keywords Nested or not in Confirmatory Tetrad Analysis
#' @export
#' @examples


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







