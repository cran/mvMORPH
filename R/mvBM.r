################################################################################
##                                                                            ##
##                               mvMORPH: mvBM                                ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor                                          ##
##                                                                            ##
################################################################################


mvBM<-function(tree,data,error=NULL,reg=NULL,model=c("BMM","BM1"),constraint=FALSE,simmap.tree=TRUE,scale.height=FALSE,method=c("L-BFGS-B","Nelder-Mead","subplex"),random.start=FALSE,control=list(maxit=20000),mod="ER",pseudoinverse=FALSE,diagnostic=TRUE,echo=TRUE){

#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}
# select default model
model<-model[1]
# Use Moore-Penrose inverse
if(pseudoinverse!=TRUE){ pseudoinverse<-solve }
##------------------------Create VCV matrix-----------------------------------##
if(simmap.tree==FALSE & is.null(reg)){
model<-"BM1"
cat("No selective regimes specified for 'reg' argument, only a BM1 model could be estimated","\n")
}
if(simmap.tree==TRUE){

if(scale.height==TRUE){
maxHeight<-max(nodeHeights(tree))
tree$edge.length<-tree$edge.length/maxHeight
tree$mapped.edge<-tree$mapped.edge/maxHeight
}

# Compute vcv for SIMMAP tree
	C1<-vcv.phylo(tree)
	if(!is.null(rownames(data))) { 
	 if(any(tree$tip.label==rownames(data))){
  C1<-C1[rownames(data),rownames(data)]
  }else if(echo==TRUE){
  
  cat("row names of the data matrix must match tip names of your phylogeny!","\n")
  }}else if(echo==TRUE){
	cat("species in the matrix are assumed to be in the same order as in the phylogeny, otherwise specify rownames of 'data'","\n")
	}
		if(model=="BMM"){
	multi.tre<-list()
  class(multi.tre)<-"multiPhylo"
	C2<-array(dim=c(nrow(C1),ncol(C1),ncol(tree$mapped.edge)))
	for(i in 1:ncol(tree$mapped.edge)){
		multi.tre[[i]]<-tree
		multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
		multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
		temp<-vcv.phylo(multi.tre[[i]])
		if(any(tree$tip.label==rownames(data))) { 
		C2[,,i]<-temp[rownames(data),rownames(data)]
		}else{
		C2[,,i]<-temp
		}
	}
 }
}else{

# Compute vcv for nodes-mapped tree
  C1<-vcv.phylo(tree)
  	if(!is.null(rownames(data))) { 
  	 if(any(tree$tip.label==rownames(data))){
  C1<-C1[rownames(data),rownames(data)]
  }else{
  cat("row names of the data matrix must match tip names of your phylogeny!","\n")
  }  
  }else if(echo==TRUE){
	cat("species in the matrix are assumed to be in the same order as in the phylogeny, otherwise specify rownames of 'data'","\n")
	}
	 	if(model=="BMM"){
	multi.tre<-list(); class(multi.tre)<-"multiPhylo"
	# number of regimes taken from OUwie (Beaulieu et al. 2012)
	regimes<-rerootingMethod(tree,reg,model=mod)
  state<-numeric(nrow(regimes$marginal.anc))
  # Use marginal ancestral state
  for(i in 1:nrow(regimes$marginal.anc)){
  state[i]<-which.max(regimes$marginal.anc[i,][])
  }
  tree$node.label<-state
	tot.states<-factor(c(tree$node.label,as.numeric(as.factor(reg))))
	k<-length(levels(tot.states))
	int.states<-factor(tree$node.label)
	tree$node.label=as.numeric(int.states)
	      #Obtain root state and internal node labels
				root.state<-tree$node.label[1]
				int.state<-tree$node.label[-1]
				n=max(tree$edge[,1])
				#New tree matrix to be used for subsetting regimes
				edges=cbind(c(1:(n-1)),tree$edge,nodeHeights(tree))
				if(scale.height==TRUE){
					edges[,4:5]<-edges[,4:5]/max(nodeHeights(tree))
				}
				edges=edges[sort.list(edges[,3]),]

				mm<-c(as.numeric(as.factor(reg)),int.state)
				mm<-as.numeric(mm)
				regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
				#Generates an indicator matrix from the regime vector
				for (i in 1:length(mm)) {
					regime[i,mm[i]] <- 1
				}
				#Finishes the edges matrix
				edges=cbind(edges,regime)
				edges=edges[sort.list(edges[,1]),]

	C2<-array(dim=c(nrow(C1),ncol(C1),k))
	for(i in 1:k){
		multi.tre[[i]]<-tree

		multi.tre[[i]]$edge.length<-edges[,5+i]*tree$edge.length

		temp<-vcv.phylo(multi.tre[[i]])
			if(!is.null(rownames(data))) { 
		C2[,,i]<-temp[rownames(data),rownames(data)]
		}else{
		C2[,,i]<-temp
		}
	} 
  }
	}
##------------------------Parameters------------------------------------------##
# number of species (tip)
n<-dim(data)[1]
# number of variables
k<-dim(data)[2]
# number of selective regimes
if(model!="BM1"){
p<-dim(C2)[3]
}else{ p<-1 } # for calculating free parameters in the LRT test

# bind data to a vector
  if(!is.matrix(data)){
data<-as.matrix(data)}
# method for the optimizer
method<-method[1]
# Taken from phytools
D <- matrix(0, n * k, k)
    for (i in 1:(n * k)){
     for (j in 1:k){
      if ((j - 1) * n < i && i <= j * n){D[i, j] = 1 }
      }
      }
# bind error to a vector
if(!is.null(error)){error<-as.vector(error)}

##------------------LogLikelihood function for multiple rates per traits------##

lik.Mult<-function(param,dat,C,D,index.mat,sig,error){ ##
  n<-dim(dat)[1] # nombre d'individus (tip)
  k<-dim(dat)[2] # nombre de variables
  p<-dim(C)[3] # correspond au nombre de sous arbres i.e. le nombre de régimes

  # matrice de taux d'évolutions
 	sig[] <- c(param)[index.mat]

  #calcul de la matrice de variance pour les différents mappings
  sig3D<-array(dim=c(k,k,p))

  for(i in 1:p){   sig3D[,,i]<-sym.par(sig[i,])} # construction d'une matrice positive definite

  # multivariate vcv matrix
  V<-matrix(0,n*k,n*k)
  for(i in 1:p){
    Cmat<-C[,,i]
    Sig<-sig3D[,,i]
    V<-V+kronecker(Sig,Cmat)
  }
  # add measurement error
  if(!is.null(error)){
			diag(V)<-diag(V)+error
		}
  #Inversion de la matrice de covariance ajustée pour les différentes branches
  dat<-as.vector(dat)
  a <- pseudoinverse(t(D) %*% pseudoinverse(V) %*% D) %*% (t(D) %*% pseudoinverse(V) %*%dat)
  ##Likelihood
  logL <- -t(dat - D %*% a) %*% pseudoinverse(V) %*% (dat - D %*% a)/2 - n * k * log(2 * pi)/2 - determinant(V)$modulus[1]/2

  list(loglik=-logL, ancstate=a)
  
}
##---------------------Loglik BM1---------------------------------------------##
lik.BM1<-function(param,dat,C,D,error){ ##
  n<-dim(dat)[1] # nombre d'individus (tip)
  k<-dim(dat)[2] # nombre de variables
  # matrice de taux d'évolutions
 	sig<-sym.par(param) 
  # multivariate vcv matrix
  V<-kronecker(sig,C)
  # add measurement error
  if(!is.null(error)){
			diag(V)<-diag(V)+error
		}
  #Inversion de la matrice de covariance ajustée pour les différentes branches
  dat<-as.vector(dat)
  a <- pseudoinverse(t(D) %*% pseudoinverse(V) %*% D) %*% (t(D) %*% pseudoinverse(V) %*%dat)
  ##Likelihood
  logL <- -t(dat - D %*% a) %*% pseudoinverse(V) %*% (dat - D %*% a)/2 - n * k * log(2 * pi)/2 - determinant(V)$modulus[1]/2
  list(loglik=-logL, ancstate=a)
}
if(model=="BMM"){
##---------------------Optimization BMM---------------------------------------##
# number of parameters
npar=(k*(k+1)/2)
# sigma matrix
sig<-matrix(1,p,npar)

# index matrix of rates
index.mat<-matrix(1:length(sig),p,npar,byrow=TRUE)

# initial values for the optimizer
sig1<-estVarBM(as.matrix(data),C1)
sig1<-sym.unpar(sig1)

starting<-NULL
for(i in 1:p){
  starting<-c(starting,sig1)
}
if(random.start==TRUE){startval<-runif(length(starting))}else{startval<-rep(1,length(starting))}
# Optimizer
if(method!="subplex"){
estim<-optim(par=starting*startval,fn=function (par) { lik.Mult(param=par, dat=data, C=C2, D=D, index.mat=index.mat,sig=sig, error=error)$loglik },control=control,hessian=TRUE,method=method)   #mettre les options maxit et method dans le menu
}else{
estim<-subplex(par=starting*startval,fn=function (par){lik.Mult(param=par,dat=data,C=C2,D=D,index.mat=index.mat,sig=sig,error=error)$loglik},control=control,hessian=TRUE)   #mettre les options maxit et method dans le menu
}

}else if(model=="BM1"){
##---------------------Optimization BM1---------------------------------------##
# number of parameters
npar=(k*(k+1)/2)
# initial values for the optimizer
sig1<-estVarBM(as.matrix(data),C1)
sig1<-sym.unpar(sig1)
if(random.start==TRUE){startval<-runif(length(sig1))}else{startval<-rep(1,length(sig1))}
# Optimizer
if(method!="subplex"){
estim<-optim(par=sig1*startval,fn=function(par){lik.BM1(param=par,dat=data,C=C1,D=D,error=error)$loglik},control=control,hessian=TRUE,method=method) 
}else{
estim<-subplex(par=sig1*startval,fn=function(par){lik.BM1(param=par,dat=data,C=C1,D=D,error=error)$loglik},control=control,hessian=TRUE) 
} 
}
##-----------------Summarizing results----------------------------------------##
if(model=="BMM"){
matResults<-matrix(1,p,npar)
matResults[]<-c(estim$par)[index.mat]
resultList<-array(dim = c(k, k, p))
states=vector()
 for(i in 1:p){
 resultList[,,i]<-sym.par(matResults[i,])
  states[i]<-multi.tre[[i]]$state
}
dimnames(resultList)<-list(colnames(data), colnames(data), states)

#ancestral states estimates
anc<-lik.Mult(param=estim$par,dat=data,C=C2,D=D,index.mat=index.mat,sig=sig,error=error)$ancstate
}else if(model=="BM1"){
 resultList<-sym.par(estim$par)
 colnames(resultList)<-colnames(data)
 rownames(resultList)<-colnames(data)
 #ancestral states estimates
anc<-lik.BM1(param=estim$par,dat=data,C=C1,D=D,error=error)$ancstate
}
# LogLikelihood
LL<--estim$value
# models parameters
if(model=="BMM"){
nparam=k+length(unique(index.mat))  #k+(p*k) = p for each regimes, k for each rates, k for each ancestral states   or(k+length(unique(index.mat))?
}else if(model=="BM1"){
nparam=k+length(estim$par)        #k+k= k for each rates and k for each ancestral states
}
# AIC
AIC<--2*LL+2*nparam
# AIC corrected
AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) #Hurvich et Tsai, 1989
##---------------------Diagnostics--------------------------------------------##

if(estim$convergence==0 & diagnostic==TRUE){  
cat("\n","successful convergence of the optimizer","\n") 
}else if(estim$convergence==1 & diagnostic==TRUE){  
cat("\n","maximum limit iteration has been reached, please consider increase maxit","\n")
}else if(diagnostic==TRUE){  
cat("\n","convergence of the optimizer has not been reached, try simpler model","\n") 
}
# Hessian eigen decomposition to check the derivatives
hess<-eigen(estim$hessian)$values
if(any(hess<0)){
hess.value<-1
if(diagnostic==TRUE){
cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
}else{
hess.value<-0
if(diagnostic==TRUE){
cat("a reliable solution has been reached","\n")}
}

##-------------------Print results--------------------------------------------##
if(echo==TRUE){
cat("\n")
cat("Summary results for multiple rates",model,"model","\n")
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat(nparam,"parameters")
cat("\n")
cat("Estimated rates matrix","\n")
print(resultList)
cat("\n")
cat("Estimated ancestral state","\n")
cat(anc)
cat("\n")
}

if(constraint==TRUE){
##-------------------Choleski constrained model-------------------------------##
##--------LogLikelihood function constrained for a unique rate per traits-----##

lik.Mult.Const<-function(param,dat,C,D,index.mat,sig,error){ ##
  n<-dim(dat)[1] # nombre d'individus (tip)
  k<-dim(dat)[2] # nombre de variables
  p<-dim(C)[3] # correspond au nombre de sous arbres i.e. le nombre de régimes

  # matrice de taux d'évolutions
 	sig[] <- c(param)[index.mat]

  #calcul de la matrice de variance pour les différents mappings
  sig3D<-array(dim=c(k,k,p))

  for(i in 1:p){
  low.chol<-build.chol(sig[i,],k)
  R<-low.chol%*%t(low.chol)
  sig3D[,,i]<-R
  }

  # multivariate vcv matrix
  V<-matrix(0,n*k,n*k)
  for(i in 1:p){
    Cmat<-C[,,i]
    Sig<-sig3D[,,i]
    V<-V+kronecker(Sig,Cmat)
  }
  # add measurement error
  if(!is.null(error)){
			diag(V)<-diag(V)+error
		}
  
  dat<-as.vector(dat)
  #Inversion de la matrice de covariance ajustée pour les différentes branches
  a <- pseudoinverse(t(D) %*% pseudoinverse(V) %*% D) %*% (t(D) %*% pseudoinverse(V) %*%dat)
  ##Likelihood
  logL <- -t(dat - D %*% a) %*% pseudoinverse(V) %*% (dat - D %*% a)/2 - n * k * log(2 * pi)/2 - determinant(V)$modulus[1]/2
  list(loglik=-logL, ancstate=a)
}
##---------------------LogLik BM1 Constraints---------------------------------##
lik.BM1.cons<-function(param,dat,C,D,error){ ##
  n<-dim(dat)[1] # nombre d'individus (tip)
  k<-dim(dat)[2] # nombre de variables
  # matrice de taux d'évolutions
 	low.chol<-build.chol(param,k)
  sig<-low.chol%*%t(low.chol)
  # multivariate vcv matrix
  V<-kronecker(sig,C)
  # add measurement error
  if(!is.null(error)){
			diag(V)<-diag(V)+error
		}
  #Inversion de la matrice de covariance ajustée pour les différentes branches
  dat<-as.vector(dat)
  a <- pseudoinverse(t(D) %*% pseudoinverse(V) %*% D) %*% (t(D) %*% pseudoinverse(V) %*%dat)
  ##Likelihood
  logL <- -t(dat - D %*% a) %*% pseudoinverse(V) %*% (dat - D %*% a)/2 - n * k * log(2 * pi)/2 - determinant(V)$modulus[1]/2
  list(loglik=-logL, ancstate=a)
}

if(model=="BMM"){
##---------------------Optimization-------------------------------------------##

# number of parameters for the constrained model
npar.cons=(k*(k-1)/2)+1
# sigma matrix
sigc<-matrix(1,p,npar.cons)
# index matrix
index.mat.cons<-matrix(1:length(sigc),p,npar.cons,byrow=TRUE)
# initial values for the optimizer
sig1<-estVarBM(as.matrix(data),C1)
# starting values following Adams (2012)
sigma.mn<-mean(diag(sig1))
R.offd<-rep(0,(k*(k-1)/2))
# # même chose mais pour chaque regimes dans la matrice index
valstart=c(sigma.mn,R.offd)
starting.cons<-NULL
for(i in 1:p){
  starting.cons<-c(starting.cons,valstart)
}
if(random.start==TRUE){startval<-runif(length(starting.cons))}else{startval<-rep(1,length(starting.cons))}
# Optimizer for constrained model
if(method!="subplex"){
estim.const<-optim(par=starting.cons*startval,fn=function(par){lik.Mult.Const(param=par,dat=data,C=C2,D=D,index.mat=index.mat.cons,sig=sigc,error=error)$loglik},control=control,hessian=TRUE,method=method)
}else{
estim.const<-subplex(par=starting.cons*startval,fn=function(par){lik.Mult.Const(param=par,dat=data,C=C2,D=D,index.mat=index.mat.cons,sig=sigc,error=error)$loglik},control=control,hessian=TRUE)
}
}else if(model=="BM1"){
##---------------------Optimization BM1 constrained---------------------------##
#Same as Adams 2012
# initial values for the optimizer
sig1<-estVarBM(as.matrix(data),C1)
# starting values following Adams (2012)
sigma.mn<-mean(diag(sig1))
R.offd<-rep(0,(k*(k-1)/2))
# # même chose mais pour chaque regimes dans la matrice index
valstart=c(sigma.mn,R.offd)
if(random.start==TRUE){startval<-runif(length(valstart))}else{startval<-rep(1,length(valstart))}
# Optimizer for constrained model
if(method!="subplex"){
estim.const<-optim(par=valstart*startval,fn=function(par){lik.BM1.cons(param=par,dat=data,C=C1,D=D,error=error)$loglik},control=control,hessian=TRUE,method=method)
}else{
estim.const<-subplex(par=valstart*startval,fn=function(par){lik.BM1.cons(param=par,dat=data,C=C1,D=D,error=error)$loglik},control=control,hessian=TRUE)
}
}
##---------------------Diagnostics--------------------------------------------##

if(estim.const$convergence==0 & diagnostic==TRUE){  
cat("\n","successful convergence of the optimizer for the constrained model","\n") 
}else if(estim.const$convergence==1 & diagnostic==TRUE){
cat("\n","maximum limit iteration has been reached, please consider increase maxit","\n")
}else if(diagnostic==TRUE){
cat("\n","convergence of the optimizer has not been reached, try simpler constrained model","\n")
}
# Hessian eigen decomposition to check the derivatives
hess.cons<-eigen(estim.const$hessian)$values
if(any(hess.cons<0)){
hess.val.cons<-1
if(diagnostic==TRUE){
cat("unreliable solution has been reached, check hessian eigenvectors or try simpler constrained model","\n") }
}else{
hess.val.cons<-0
if(diagnostic==TRUE){
cat("a reliable solution has been reached","\n")}
}

##-----------------Summarizing results----------------------------------------##

# LogLikelihood
LLc<--estim.const$value

if(model=="BMM"){
nparamc=k+length(unique(index.mat.cons)) #p+k = p for each regimes, k for each ancestral states, 1 rates   #
}else if(model=="BM1"){
nparamc=k+length(estim.const$par)   #1+k = k for each ancestral states, 1 for rates
}
# maximum likelihood estimates of rates matrix
if(model=="BMM"){
matResultsC<-matrix(1,p,npar.cons)
matResultsC[]<-c(estim.const$par)[index.mat.cons]
resultListC<-array(dim = c(k, k, p))
statesC=vector()
 for(i in 1:p){
 low.chol<-build.chol(matResultsC[i,],k)
 R<-low.chol%*%t(low.chol)
 resultListC[,,i]<-R
  statesC[i]<-multi.tre[[i]]$state
}
dimnames(resultListC)<-list(colnames(data), colnames(data), statesC)

# ancestral states estimates
anc.c<-lik.Mult.Const(param=estim.const$par,dat=data,C=C2,D=D,index.mat=index.mat.cons,sig=sigc,error=error)$ancstate     
}else if(model=="BM1"){
 low.chol<-build.chol(estim.const$par,k)
 resultListC<-low.chol%*%t(low.chol)
 colnames(resultListC)<-colnames(data)
 rownames(resultListC)<-colnames(data)
# ancestral states estimates
anc.c<-lik.BM1.cons(param=estim.const$par,dat=data,C=C1,D=D,error=error)$ancstate
}
# AIC
Cons.AIC<--2*LLc+2*nparamc
# AIC corrected
Cons.AICc<-Cons.AIC+((2*nparamc*(nparamc+1))/(n-nparamc-1)) #Hurvich et Tsai, 1995

##-------------------Print results--------------------------------------------##
if(echo==TRUE){
cat("\n")
cat("Summary results for the constrained",model,"model","\n")
cat("LogLikelihood:","\t",LLc,"\n")
cat("AIC:","\t",Cons.AIC,"\n")
cat("AICc:","\t",Cons.AICc,"\n")
cat(nparamc,"parameters")
cat("\n")
cat("Estimated rates matrix","\n")
print(resultListC)
cat("\n")
cat("Estimated ancestral state","\n")
cat(anc.c,"\n")
cat("\n")
} 
##-------------------LRT comparison of the models-----------------------------##

LRT<-(2*((LL-LLc)))
#nvariables-1 degrees of freedom
LRT.prob<-pchisq(LRT,(k-1)*p,lower.tail=FALSE)  
if(echo==TRUE){ 
cat("LRT p-value:",LRT.prob,"\n")
}
}# end of constrained model

##-------------------Store results--------------------------------------------##

if(constraint==TRUE){
 results<-list(LogLik.m=LL, AIC.mult=AIC, AICc.mult=AICc, rates.m=resultList, anc=anc, LogLik.cons=LLc, AIC.cons=Cons.AIC, AICc.cons=Cons.AICc, rates.cons=resultListC, anc.cons=anc.c, LRT=LRT, pval=LRT.prob, convergence=estim$convergence, convergence.const=estim.const$convergence, hess.values=hess.value, hess.values.cons=hess.val.cons)
}else{
 results<-list(LogLik.m=LL, AIC.mult=AIC, AICc.mult=AICc, rates.m=resultList, anc=anc,convergence=estim$convergence, hess.values=hess.value)
}
#End
}