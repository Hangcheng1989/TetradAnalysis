#' Confirmatory Tetrad Analysis Function
#'
#' This function allows you to run Confirmatory Tetrad Analysis from Bollen 1993.
#' @param "Samplecov", "FactorModel", and "Size"
#' @keywords Confirmatory Tetrad Analysis
#' @export
#' @examples
#' TetradAnalysisNoRandom()


##################################################################################
#################################### Tetrad Function #################################
##################################################################################
TetradAnalysisNoRandom<-function(Samplecov,FactorModel,Size){

### "sweep_all" function is necessary to find Nonredundant vanishing tetrads from empirical vanishing tetrads ###
sweep<-function(A,k) {
	B<-matrix(0,nrow=nrow(A),ncol=ncol(A))
	if (abs(A[k,k])<2.2*10^(-10)){
		B=A
		for (i in 1:nrow(B)){
			B[i,k]=0
			B[k,i]=0
		}
	}
	
	else if (A[k,k] !=0) {
		for (i in 1:nrow(A)) {
			for (j in 1: ncol(A)) {
			if (i== k & j==k){
			  	B[i,j]=1/A[k,k]
				}
			else if (i == k & j != k) {
				B[i,j] = A[k,j]/A[k,k]
				}
			else if (i != k & j == k) {
				B[i,j] = -A[i,k]/A[k,k]
				}
			else if (i != k & j != k) {
				B[i,j] = A[i,j] - ( A[i,k]*A[k,j] / A[k,k])
				}
		 	}
		}
	}
	return(B)
}
sweep_all<-function(A) {
	B<-matrix(0,nrow=nrow(A),ncol=ncol(A))
	for (i in 1: nrow(A)) {
		B=sweep(A,i)
		A=B
	}
	return(B)
}

##################################################################################
### Tetrad Function Part 1: Use Empirical Method to find model implied tetrad #####################
##################################################################################

### Information needed for Part 1 ###
### 1. Sample covariance matrix "Samplecov". ###
### 2. SEM model need to be tested "FactorModel". ###
### 3. Sample size "Size". ###

### Use function "ExampleTetrad" to get all model implied vanishing tetrads using empirical method. ###
EmpiricalMethod<-function(FactorModel,Samplecov){
### Package "lavaan" is used to get model implied covariance matrix "Modelcov", the synatax used to sepecify model is the same as syntax in "lavaan". ###
	library(lavaan)
	FactorOut<-cfa(FactorModel, sample.cov=Samplecov, sample.nobs =Size)
	Modelcov<-data.frame(lavTech(FactorOut, "cov.ov"))*Size/(Size-1)
### List all tetrads in matrix "Tetrad", when there are N observed variables in inputted model. ###
	N<-ncol(Modelcov)
	x<-combn(1:N,4)
	TetradNum<-choose(N, 4)
	x1<-t(x)
	x2<-cbind(x[1,],x[4,],x[2,],x[3,])
	x3<-cbind(x[1,],x[3,],x[4,],x[2,])
	y<-rbind(x1,x2,x3)
	rownames(y)<-rep(1:TetradNum,3)
	Y<-y[order(as.numeric(rownames(y))),]
	Tetrad<-Y[,1]*1000+Y[,2]*100+Y[,3]*10+Y[,4]
### Calculate all tetrad values based on model implied covariance matrix "Modelcov". ###
### Save all tetrad values in a vector named "T", if the absolute tetrad value is less than 0.001, we assume it is a vanishing tetrad. ###
	T<-matrix(1,nc=1,nr=nrow(Y))
	for(i in 1:nrow(Y)){
	T[i,1]<-round(Modelcov[Y[i,2],Y[i,1]]*Modelcov[Y[i,4],Y[i,3]]
		  -Modelcov[Y[i,3],Y[i,1]]*Modelcov[Y[i,4],Y[i,2]],2)
	}
	colnames(T) <- c("ImpliedValue")
	implied<-matrix(,nrow=nrow(T))
	colnames(implied) <- c("Implied")
	for(i in 1:nrow(T)){
	if(abs(T[i,1])>0.001) {implied[i,1]<-0} else{implied[i,1]<-1}
	}
	Tvalue<-matrix(1,nc=4,nr=nrow(Y))
	for (i in 1:nrow(Y)){
	Tvalue[i,1]<-min(Y[i,1],Y[i,2])*10+max(Y[i,1],Y[i,2])
	Tvalue[i,2]<-min(Y[i,3],Y[i,4])*10+max(Y[i,3],Y[i,4])
	Tvalue[i,3]<-min(Y[i,1],Y[i,3])*10+max(Y[i,1],Y[i,3])
	Tvalue[i,4]<-min(Y[i,2],Y[i,4])*10+max(Y[i,2],Y[i,4])
	}
	colnames(Tvalue) <- c("1","2","3","4")
	EmpiricalMethod_list<-list(EmpiricalMethod=cbind(Tetrad,Tvalue,implied,T), Modelcov=Modelcov)
	return(EmpiricalMethod_list)
}

##################################################################################
########## Tetrad Function Part 2: Find Nonredundant vanishing tetrads ########################
##################################################################################

### Information needed for Part 2 ###
### 1. Model implied vanishing tetrads "ModelTetrad". ###
ExampleTetrad<-EmpiricalMethod(FactorModel,Samplecov)				# Use "EmpiricalMethod" function. #
ModelTetrad<-ExampleTetrad$EmpiricalMethod[ExampleTetrad$EmpiricalMethod[,6]==1,]	# Implied=1 means tetrad vanished in model. #
EmpiricalTetrads<-ModelTetrad[,1,drop=FALSE]					# Model implied tetrad names. #
### 2. All covariances in the model. ###
COV<-ModelTetrad[,2:5,drop=FALSE]
### 3. Model implied covariance matrix "Modelcov". ###
Modelcov<-ExampleTetrad$Modelcov

##################################################################################
### Use function "NRVT_Fun" to identify model implied nonredundant vanishing tetrads. ###
NRVT_Fun<-function(COV,ModelTetrad,Modelcov) {
 	if (nrow(COV)==1){			# If there is only one model implied vanishing tetrad, identify is unnecessary. #
 		NRVT<-ModelTetrad[1,1]
		NRVT_Num<-1
	}

	else if (nrow(COV)>1) {		# If there is more than one model implied vanishing tetrad, identify is necessary. #
		uc<-unique(as.vector(COV))	# Find all unique covariances among vanishing tetrads. #
		SigmaSS<-matrix(1,nc=length(uc),nr=length(uc))
		for (i in 1:length(uc)){
			for (j in 1:length(uc)) {
			e<-uc[i]%/%10
			f<-uc[i]%%10
			g<-uc[j]%/%10
			h<-uc[j]%%10
			efgh<-Modelcov[e,g]*Modelcov[f,h]+Modelcov[e,h]*Modelcov[f,g]
			SigmaSS[i,j]<-efgh
			}
		}
### Get model implied SigmaTT based on model implied covariance matrix "Modelcov". ###
		dtaudsigmaCOV<-matrix(0,ncol=length(uc),nrow=nrow(COV))
		for (i in 1:nrow(dtaudsigmaCOV)){
			for(j in 1:ncol(dtaudsigmaCOV)){
			if (COV[i,1]==uc[j]) {dtaudsigmaCOV[i,j]=Modelcov[COV[i,2]%/%10,COV[i,2]%%10]} 
			if (COV[i,2]==uc[j]) {dtaudsigmaCOV[i,j]=Modelcov[COV[i,1]%/%10,COV[i,1]%%10]} 
			if (COV[i,3]==uc[j]) {dtaudsigmaCOV[i,j]=-Modelcov[COV[i,4]%/%10,COV[i,4]%%10]} 
			if (COV[i,4]==uc[j]) {dtaudsigmaCOV[i,j]=-Modelcov[COV[i,3]%/%10,COV[i,3]%%10]} 
			}
		}
### use sweep operator to identify set of nonredundant vanishing tetrads. ###
		SigmaTT<-dtaudsigmaCOV %*% SigmaSS %*% t(dtaudsigmaCOV)
		SweepResults<-sweep_all(SigmaTT)
		rowsum<-rowSums(SweepResults)
		nrvt <- matrix(0, ncol=1, nrow=nrow(COV))
		NRVT_Num<-0
		for (i in 1: nrow(COV)) {
			if ( rowsum[i] != 0 ) {			# If rowsum==0, that tetrad is redundant.#
			nrvt[i,] = ModelTetrad[i,1]
			NRVT_Num<-NRVT_Num+1
				}
			}
	NRVT<-nrvt[nrvt[,1]!=0,]
	
	}
	NRVT_list<-list(NRVT=NRVT, NRVT_Num=NRVT_Num, uc=uc,Modelcov=Modelcov)
	return(NRVT_list)
}

NRVT_Results<-NRVT_Fun(COV,ModelTetrad,Modelcov)
NRVT_Num<-NRVT_Results$NRVT_Num
NRVT<-NRVT_Results$NRVT
UniqueCov<-NRVT_Results$uc

##################################################################################
########## Tetrad Function Part 3: Individual Test and Multivariate Test  #########################
##################################################################################

### Information needed for Part 3 ###
### 1. Model implied nonredundant vanishing tetrads "NRVT". ###
NRVT_Results<-NRVT_Fun(COV,ModelTetrad,Modelcov)	# Use "NRVT_Fun" function. #
NRVT<-NRVT_Results$NRVT
### 2. Model implied nonredundant vanishing tetrad number "NRVT_Num" ###
NRVT_Num<-NRVT_Results$NRVT_Num
### 3. All unique covariances in vanishing tetrads "uc". ###
uc<-NRVT_Results$uc

##################################################################################
######### Individual test for every nonredundant tetrad #####################################
Tetrad<-matrix(1,nc=1,nr=NRVT_Num)
AVAR<-matrix(1,nc=1,nr=NRVT_Num)
Teststat<-matrix(1,nc=1,nr=NRVT_Num)
Pvalue<-matrix(1,nc=1,nr=NRVT_Num)
Fun<-t(t(NRVT))				# Get all nonredundant tetrads. #
NRVT_COV<-matrix(0, nrow=NRVT_Num,ncol=4)	
for (k in 1:NRVT_Num){
		e<-Fun[k,1]%/%1000
		f<-Fun[k,1]%/%100%%10
		g<-Fun[k,1]%/%10%%10
		h<-Fun[k,1]%%10
	sigma<-c(min(e,f)*10+max(e,f),min(g,h)*10+max(g,h),min(e,g)*10+max(e,g),min(f,h)*10+max(f,h))
	NRVT_COV[k,1]<-sigma[1]
	NRVT_COV[k,2]<-sigma[2]
	NRVT_COV[k,3]<-sigma[3]
	NRVT_COV[k,4]<-sigma[4]
	SigmaSS<-matrix(1,nc=length(sigma),nr=length(sigma))	# Get SigmaSS based on Sample covariance matrix "Samplecov". #
		for (i in 1:length(sigma)){
			for (j in 1:length(sigma)) {
			ee<-sigma[i]%/%10
			ff<-sigma[i]%%10
			gg<-sigma[j]%/%10
			hh<-sigma[j]%%10
			efgh<-Samplecov[ee,gg]*Samplecov[ff,hh]+Samplecov[ee,hh]*Samplecov[ff,gg]
		SigmaSS[i,j]<-efgh
			}
		}
dt<-matrix(c(Samplecov[g,h],Samplecov[e,f],-Samplecov[f,h],-Samplecov[e,g]),ncol=4,nrow=1)
VAR<-(1/Size)* dt %*% SigmaSS %*% t(dt)
T<-Samplecov[e,f]*Samplecov[g,h]-Samplecov[e,g]*Samplecov[f,h]		# Calculate T values for individual tests. #
stat<-T^2/VAR
Chipvalue<-1-pchisq(stat,df=1)						# Calculate P-values for individual tests. #
	Tetrad[k,1]<-round(T,3)
	AVAR[k,1]<-round(VAR,3)
	Teststat[k,1]<-round(stat,3)
	Pvalue[k,1]<-round(Chipvalue,3)
}
Result<-cbind(Fun,Tetrad,AVAR,Teststat,Pvalue)
colnames(Result) <- c("Model Implied Tetrad","t","AVAR","TestStatistic","P-value")	# All results for individual tests. #

######### Multivariate test for overall model fit #################################################
### Get SigmaSS based on sample covariance matrix "Samplecov". ###
	SigmaSS<-matrix(1,nc=length(uc),nr=length(uc))	
		for (i in 1:length(uc)){
			for (j in 1:length(uc)) {
			e<-uc[i]%/%10
			f<-uc[i]%%10
			g<-uc[j]%/%10
			h<-uc[j]%%10
			efgh<-Samplecov[e,g]*Samplecov[f,h]+Samplecov[e,h]*Samplecov[f,g]
			SigmaSS[i,j]<-efgh
				}
			}

### Get SigmaTT, use dtau/dsigma and based on sample covariance matrix "Samplecov". ###
	dtaudsigmaCOV<-matrix(0,ncol=length(uc),nrow=NRVT_Num)
		for (i in 1:nrow(dtaudsigmaCOV)){
			for(j in 1:ncol(dtaudsigmaCOV)){
			if (NRVT_COV[i,1]==uc[j]) {dtaudsigmaCOV[i,j]=Samplecov[NRVT_COV[i,2]%/%10,NRVT_COV[i,2]%%10]} 
			if (NRVT_COV[i,2]==uc[j]) {dtaudsigmaCOV[i,j]=Samplecov[NRVT_COV[i,1]%/%10,NRVT_COV[i,1]%%10]} 
			if (NRVT_COV[i,3]==uc[j]) {dtaudsigmaCOV[i,j]=-Samplecov[NRVT_COV[i,4]%/%10,NRVT_COV[i,4]%%10]} 
			if (NRVT_COV[i,4]==uc[j]) {dtaudsigmaCOV[i,j]=-Samplecov[NRVT_COV[i,3]%/%10,NRVT_COV[i,3]%%10]} 
			}
		}
SigmaTT<-dtaudsigmaCOV %*% SigmaSS %*% t(dtaudsigmaCOV)

t<-t(t(Result[,2]))				# Get T values for all individual tests, use function "Result". #
T<-Size*t(t) %*% solve(SigmaTT) %*% t		# T value for multivariate test. #
MultiPvalue<-1-pchisq(T,df=NRVT_Num)		# P-value for multivariate test. #
colnames(T) <- c("TestStatistic For Multivariate Test")
colnames(MultiPvalue) <- c("P-value For Multivariate Test")

##################################################################################
########## Tetrad Function Part 4: Output all analysis results #################################
##################################################################################

results<-list(T=T, MultiPvalue=MultiPvalue,NRVT=NRVT,NRVT_Num=NRVT_Num,
		  EmpiricalTetrads=EmpiricalTetrads,Modelcov=Modelcov)
return(results)
}

########### Tetrad Function End#######################################################
