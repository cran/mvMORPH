################################################################################
##                                                                            ##
##                               mvMORPH: mvOU                                ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, subplex                                 ##
##                                                                            ##
################################################################################

mvOU<-function(tree,data,error=NULL,sigma=NULL,alpha=NULL,model=c("OUM","OU1"),scale.height=FALSE,diagnostic=TRUE,method=c("L-BFGS-B","Nelder-Mead","subplex"),pseudoinverse=FALSE,echo=TRUE,control=list(maxit=20000)){

#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}
# Use Moore-Penrose inverse
if(pseudoinverse!=TRUE){ pseudoinverse<-solve }
# Check the order of the dataset and the phylogeny 

	if(!is.null(rownames(data))) { 
   if(any(tree$tip.label==rownames(data))){
	data<-data[tree$tip.label,] }else if(echo==TRUE){
  cat("row names of the data matrix must match tip names of your phylogeny!","\n")
	}}else if(echo==TRUE){
	cat("species in the matrix are assumed to be in the same order as in the phylogeny, otherwise specify rownames of 'data'","\n")
	}

##-------------------------Calculation of parameters--------------------------##
  # number of species
n<-length(tree$tip.label)

  # choose model
model=model[1]
  # choose method for the optimizer
method=method[1]

  # max node height for standardisation
mb<-max(nodeHeights(tree))
  # regimes number
  if(model=="OUM"){
k<-length(colnames(tree$mapped.edge))
}else{ k<-1 }
  # number of traits
p<-ncol(data)
  # number of submatrix
par<-p*p
  # ancestor times calculation
mt<-mTime(tree,scale.height)
  # bind data to a vector
  if(is.matrix(data)){
dat<-as.vector(data) }else{ dat<-as.vector(as.matrix(data))}
  # bind error to a vector
if(!is.null(error)){error<-as.vector(error)}


  # initial alpha and sigma matrix if not provided
if(is.null(sigma)){sigma=sym.unpar(diag(1,p))}
if(is.null(alpha)){alpha=sym.unpar(diag(1,p))}
nalpha<-length(alpha)
nsigma<-length(sigma)

##-----------------------Precalculate regime indexation-----------------------##
  # root to tip lineage indexation
root2tip <- .Call("seq_root2tipM", tree$edge, n, tree$Nnode)
# Si OU1 sur un objet 'phylo'
if(model=="OU1"){
if(scale.height==TRUE){
valLineage<-sapply(1:n,function(z){rev(unlist(
sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$edge.length[val]<-tree$edge.length[val]/mb},simplify=FALSE)))
} ,simplify=FALSE)}else{
valLineage<-sapply(1:n,function(z){rev(unlist(
sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$edge.length[val]<-tree$edge.length[val]},simplify=FALSE)))
} ,simplify=FALSE)
}
}else{
# Données de temps par régimes et par branches
if(scale.height==TRUE){
valLineage<-sapply(1:n,function(z){rev(unlist(
sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$maps[[val]]<-tree$maps[[val]]/mb},simplify=FALSE)))
} ,simplify=FALSE)
}else{
valLineage<-sapply(1:n,function(z){rev(unlist(
sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$maps[[val]]<-tree$maps[[val]]},simplify=FALSE)))
} ,simplify=FALSE)
}
}

# Indexer les régimes
if(model=="OUM"){ 
# Indexing factors
facInd<-factor(colnames(tree$mapped.edge))
indice<-lapply(1:n,function(z){rev(unlist(
lapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); factor(names(tree$maps[[val]]),levels=facInd)})))
})
}else if(model=="OU1"){
indice<-lapply(1:n,function(z){ as.factor(rep(1,length(valLineage[[z]])))})
}
# Liste avec dummy matrix
indiceA<-indiceReg(n,indice)
listReg<-sapply(1:n,function(x){sapply(1:p,function(db){regimeList(indiceA[[x]],k=k)},simplify=FALSE)},simplify=FALSE)

# mapped epochs
epochs<-sapply(1:n,function(x){lineage<-as.numeric(c(cumsum(valLineage[[x]])[length(valLineage[[x]])],(cumsum(valLineage[[x]])[length(valLineage[[x]])]-cumsum(valLineage[[x]])))); lineage[which(abs(lineage)<1e-15)]<-0; lineage },simplify=FALSE)



##-----------------------Likelihood Calculation-------------------------------##

	devianc<-function(alpha,sigma,dat,error,mt){

		eig<-eigen(alpha)
		N<-length(dat)# 

    V<-.Call("simmap_covar",nterm=as.integer(n),bt=mt$mDist,lambda=eig$values,S=eig$vectors,sigma.sq=sigma)
		#Transformed C code from ouch
		W<-.Call("simmap_weights",nterm=as.integer(n), epochs=epochs,lambda=eig$values,S=eig$vectors,beta=listReg)
  
		
		if (any(is.nan(diag(V))) || any(is.infinite(diag(V)))) return(1000000)

		if(!is.null(error)){
			diag(V)<-diag(V)+error
		}

		theta<-Inf
		try(theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%dat, silent=TRUE)
		if(any(theta==Inf)){
			return(10000000)
		}

		DET<-determinant(V, logarithm=TRUE)

		logl<--.5*(t(W%*%theta-dat)%*%pseudoinverse(V)%*%(W%*%theta-dat))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))  
		if(is.infinite(logl)){
			return(10000000)
		}
		return(-logl)
	
	}
##----------------------Likelihood optimization-------------------------------##
if(method!="subplex"){
estim <- optim(
                   par=c(alpha,sigma),
                   fn = function (par) {
                     devianc(
                               alpha=sym.par(par[seq_len(nalpha)]),
                               sigma=sym.par(par[nalpha+seq_len(nsigma)]),
                               error=error,
                               dat=dat,
                               mt=mt
                               )
                   },
                   gr=NULL,
                   hessian=TRUE,
                   method=method,
                   control=control
                   )
}else{
estim <- subplex(
                   par=c(alpha,sigma),
                   fn = function (par) {
                     devianc(
                               alpha=sym.par(par[seq_len(nalpha)]),
                               sigma=sym.par(par[nalpha+seq_len(nsigma)]),
                               error=error,
                               dat=dat,
                               mt=mt
                               )
                   },
                   hessian=TRUE,
                   control=control
                   )
}
                   
##---------------------theta estimation---------------------------------------##

est.theta<-function(estimML){

    alpha=sym.par(estimML[seq_len(nalpha)])
    sigma=sym.par(estimML[nalpha+seq_len(nsigma)])
  	N<-length(dat) 
    eig<-eigen(alpha)

    V<-.Call("simmap_covar",nterm=as.integer(n),bt=mt$mDist,lambda=eig$values,S=eig$vectors,sigma.sq=sigma)
   	if(!is.null(error)){
			diag(V)<-diag(V)+error
		}
		#Transformed C code from ouch
		W<-.Call("simmap_weights",nterm=as.integer(n), epochs=epochs,lambda=eig$values,S=eig$vectors,beta=listReg)
		theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%dat
		
	  theta
		}
res.theta<-est.theta(estim$par)

##---------------------Diagnostics--------------------------------------------##
hess<-eigen(estim$hessian)$values

if(estim$convergence==0 & diagnostic==TRUE){ 
cat("successful convergence of the optimizer","\n")
}else if(estim$convergence==1 & diagnostic==TRUE){  
cat("maximum limit iteration has been reached, please consider increase maxit","\n")
}else if(diagnostic==TRUE){
cat("convergence of the optimizer has not been reached, try simpler model","\n")
}

if(any(hess<0)){
hess.val<-1  
if(diagnostic==TRUE){
cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
}else{
hess.val<-0
if(diagnostic==TRUE){
cat("a reliable solution has been reached","\n")}
}

##-------------------Summarize Results----------------------------------------##
LL<- -estim$value 
#nparam=nalpha+nsigma+(k*p)   # for all parameters  (optional)
#or for free parameters
nparam=(2*p)+(k*p)
# maximum likelihood estimates of alpha and sigma
estim.alpha<-estim$par[seq_len(nalpha)]
estim.sigma<-estim$par[nalpha+seq_len(nsigma)]
# estimates standards errors
estim.se<-sqrt(diag(pseudoinverse(estim$hessian)))
estim.se.alpha<-estim.se[seq_len(nalpha)]
estim.se.sigma<-estim.se[nalpha+seq_len(nsigma)]
# alpha matrix
alpha.mat<-sym.par(estim.alpha)
# sigma matrix
sigma.mat<-sym.par(estim.sigma)
# AIC
AIC<- -2*LL+2*nparam
# AIC corrected
AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) #Hurvich et Tsai, 1995
# matrix of estimated theta values
theta.mat<-matrix(res.theta,k)
if(model=="OUM"){
rownames(theta.mat)<-colnames(tree$mapped.edge)}else{
rownames(theta.mat)<-"OU1"}
colnames(theta.mat)<-colnames(data)

##-------------------Print results--------------------------------------------##
if(echo==TRUE){
cat("\n")
cat("Summary results","\n")
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat("\n")
cat("Estimated theta values","\n")
print(theta.mat)
cat("\n")
cat("ML alpha values","\n")
print(alpha.mat)
cat("\n")
cat("ML sigma values","\n")
print(sigma.mat)
}

##------------------List results----------------------------------------------##
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha.mat=alpha.mat, sigma.mat=sigma.mat, alpha=estim.alpha, sigma=estim.sigma, alpha.se=estim.se.alpha, sigma.se=estim.se.sigma, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val ) 


}
