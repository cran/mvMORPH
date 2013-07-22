################################################################################
##                                                                            ##
##                               mvMORPH: mvEB                                ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, geiger                                  ##
##                                                                            ##
################################################################################

mvEB<-function(tree,data,error=NULL,low=-3, up=0, tol=c(0.00000001,Inf), scale.height=FALSE,control=list(maxit=20000),pseudoinverse=FALSE, diagnostic=TRUE,echo=TRUE){

# Use Moore-Penrose inverse
if(pseudoinverse!=TRUE){ pseudoinverse<-solve }
# scale height of the tree
if(scale.height==TRUE){
maxHeight<-max(nodeHeights(tree))
tree$edge.length<-tree$edge.length/maxHeight
}
# number of traits
p<-ncol(data)
if(is.null(p)){
p<-1
}
# bind error to a vector
if(!is.null(error)){error<-as.vector(error)}
# number of specimens
n<-length(tree$tip.label)
# bind error to a vector
if(!is.null(error)){error<-as.vector(error)}
# number of parameters
npar<-(p*(p+1)/2)
# compute the vcv for starting values
C1<-vcv.phylo(tree)
# initial values for the optimizer
sig1<-estVarBM(as.matrix(data),C1)
sig1<-sym.unpar(sig1)
# initial value for exponent parameter
ebval<-runif(1,low,up)
# taken from phytools
D <- matrix(0, n * p, p)
    for (i in 1:(n * p)){
     for (j in 1:p){
      if ((j - 1) * n < i && i <= j * n){D[i, j] = 1 }
      }
      }
##--------------Maximum Likelihood estimation of EB---------------------------##  
## modified function LikLambda from Freckleton 2012
likEB <- function(dat, error, tree, a, sig) { 
  treeEB<-ebTree(tree, a=a)
	Vmat<-vcv.phylo(treeEB)
  n<-dim(dat)[1] # nombre d'individus (tip)
  k<-dim(dat)[2] # nombre de variables

  # multivariate vcv matrix
  V<-kronecker(sig,Vmat)
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
##---------------Optimizing function------------------------------------------##
	lower = c(low,rep(tol[1], npar))
	upper = c(up,rep(tol[2], npar))
estim<-NA	
try( estim <- optim(
                   par=c(ebval,sig1),
                   fn = function (par) {
                     likEB(
                               a=par[1],
                               error=error,
                               dat=as.matrix(data),
                               tree=tree,
                               sig=sym.par(par[1+seq_len(npar)])
                               )$loglik
                   },
                   gr=NULL,
                   hessian=TRUE,
                   method ="L-BFGS-B",
                   lower = lower, 
                   upper = upper,
                   control=control
                   ),TRUE )
# Try with pseudo-inverse if any NA
 if(any(is.na(estim))==TRUE){
  lower = c(low,rep(-Inf, npar))
 	upper = c(up,rep(Inf, npar))
 	pseudoinverse<-pseudoinverse
   multiTry(estim <- optim(
                   par=c(ebval,sig1)*runif(length(sig1)+1),
                   fn = function (par) {
                     likEB(
                               a=par[1],
                               error=error,
                               dat=as.matrix(data),
                               tree=tree,
                               sig=sym.par(par[1+seq_len(npar)])
                               )$loglik
                   },
                   gr=NULL,
                   hessian=TRUE,
                   method ="L-BFGS-B",
                   lower = lower, 
                   upper = upper,
                   control=control
                   ) )
 }  
 # Try with new starting values              
 if(any(eigen(estim$hessian)$values<0) | estim$convergence!=0){
    	lower = c(low,rep(-Inf, npar))
    	upper = c(up,rep(Inf, npar))
    estim <- optim(
                   par=c(ebval,sig1)*runif(length(sig1)+1),
                   fn = function (par) {
                     likEB(
                               a=par[1],
                               error=error,
                               dat=as.matrix(data),
                               tree=tree,
                               sig=sym.par(par[1+seq_len(npar)])
                               )$loglik
                   },
                   gr=NULL,
                   hessian=TRUE,
                   method ="L-BFGS-B",
                   lower = lower, 
                   upper = upper,
                   control=control
                   )
    }

##-----------------Summarizing results----------------------------------------##
# Rate matrix
resultList<-sym.par(estim$par[1+seq_len(npar)])
colnames(resultList)<-colnames(data)
rownames(resultList)<-colnames(data)
#ancestral states estimates
anc<-likEB(a=estim$par[1],dat=as.matrix(data),tree=tree,sig=sym.par(estim$par[1+seq_len(npar)]),error=error)$ancstate
# rate parameter
r=estim$par[1]
# LogLikelihood
LL<--estim$value
# models parameters
nparam=1+(2*p)  #p ancestral states, p rates, 1 EB 'g' parameter
# AIC
AIC<--2*LL+2*nparam
# AIC corrected
AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) #Hurvich et Tsai, 1995
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
cat("Summary results for Early Burst or ACDC model","\n")
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat("Rate change:","\t",r,"\n")
cat("\n")
cat("Estimated rates matrix","\n")
print(resultList)
cat("\n")
cat("Estimated ancestral state","\n")
cat(anc)
cat("\n")
}

##-------------------Store results--------------------------------------------##
 
results<-list(LogLik.m=LL, AIC=AIC, AICc=AICc, r=r, rates.m=resultList, anc=anc, convergence=estim$convergence, hess.values=hess.value)

}