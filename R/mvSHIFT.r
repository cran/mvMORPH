################################################################################
##                                                                            ##
##                               mvMORPH: mvSHIFT                             ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, subplex                                 ##
##                                                                            ##
################################################################################

mvSHIFT<-function(tree,data,age=NULL,error=NULL,sigma=NULL,alpha=NULL,sig=NULL,model=c("ER","RR","EC","SR"),scale.height=FALSE,diagnostic=TRUE,method=c("L-BFGS-B","Nelder-Mead","subplex"),pseudoinverse=FALSE,echo=TRUE,control=list(maxit=20000)){

#set age shift with make.era from phytools if the age is provided
if(!is.null(age)){
 tot<-max(nodeHeights(tree))
 limit=tot-age
 tree<-make.era.map(tree,  c(0,limit))
}
#Shift rate model
if(model=="SR"){
mvBM(tree,data,model="BMM",scale.height=scale.height,diagnostic=diagnostic,method=method,pseudoinverse=pseudoinverse,echo=echo,control=control,error=error) 
}else{
# Use Moore-Penrose inverse
if(pseudoinverse!=TRUE){ pseudoinverse<-solve }
#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}


#set parameters for ecological constraint model
if(model=="EC" || model=="constraint"){
before<-2
after<-1
model<-"RR"
nmod<-"Ecological Constraint"
}else{
before<-1
after<-2
if(model=="RR"){
nmod<-"Release and Radiate"
}else if(model=="ER"){
nmod<-"Ecological Release"
}
}

#scale tree
if(scale.height==TRUE){
maxHeight<-max(nodeHeights(tree))
tree$edge.length<-tree$edge.length/maxHeight
tree$mapped.edge<-tree$mapped.edge/maxHeight
}

  #multiple vcv matrix
multi.tre<-list()
  class(multi.tre)<-"multiPhylo"
	vcvList<-list()
	for(i in 1:ncol(tree$mapped.edge)){
		multi.tre[[i]]<-tree
		multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
		multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
		temp<-vcv.phylo(multi.tre[[i]])
		if(any(tree$tip.label==rownames(data))) { 
		vcvList[[i]]<-temp[rownames(data),rownames(data)]
		}else{
		vcvList[[i]]<-temp
		}
	}

  # number of species (tip)
n<-dim(data)[1]
  # number of variables
p<-dim(data)[2]
  # choose model
model=model[1]
  # choose method for the optimizer
method=method[1]
  # initial alpha and sigma matrix if not provided
if(is.null(sigma)){sigma=sym.unpar(diag(1,p))}
if(is.null(alpha)){alpha=sym.unpar(diag(1,p))}
if(is.null(sig)){sig=sym.unpar(diag(1,p))}
nalpha<-length(alpha)
nsigma<-length(sigma)
  # conditions
if(model=="release" || model=="ER"){
nsig<-nsigma
}else if(model=="radiate" || model=="RR"){
nsig<-nsigma*2
}
  
  # bind data to a vector
  if(is.matrix(data)){
dat<-as.vector(data) }else{ dat<-as.vector(as.matrix(data))}
  # method for the optimizer
method<-method[1]
  # Taken from phytools
W <- matrix(0, n * p, p)
    for (i in 1:(n * p)){
     for (j in 1:p){
      if ((j - 1) * n < i && i <= j * n){W[i, j] = 1 }
      }
      }
  # bind error to a vector
if(!is.null(error)){error<-as.vector(error)}

##------------------------Likelihood function---------------------------------##
# release and radiate
lik.shift<-function(alpha,sigma,sig,dat,error,vcvList){ 
  
  #calcul de la matrice de variance pour le modèle OU + BM
  eig<-eigen(alpha)
  Vou<-.Call("simmap_covar",nterm=as.integer(n),bt=vcvList[[before]],lambda=eig$values,S=eig$vectors,sigma.sq=sigma)
  V<-Vou+kronecker(sig,vcvList[[after]])
  
  # add measurement error
  if(!is.null(error)){
			diag(V)<-diag(V)+error
		}

  dat<-as.vector(dat)
  ## ancestral states
  a <- pseudoinverse(t(W) %*% pseudoinverse(V) %*% W) %*% (t(W) %*% pseudoinverse(V) %*%dat)
  ##Likelihood
  logL <- -t(dat - W %*% a) %*% pseudoinverse(V) %*% (dat - W %*% a)/2 - n * p * log(2 * pi)/2 - determinant(V)$modulus[1]/2

  return(-logL)
}


# rate shift
#mvBM(tree,data,scale.height=scale.height,maxit=maxit,method=method,random.start=FALSE,diagnostic=TRUE,echo=TRUE)

##----------------------Optimization------------------------------------------##
# starting values
if(model=="RR" || model=="radiate"){
starting<-c(alpha,sigma,sig)
}else if(model=="ER" || model=="release"){
starting<-c(alpha,sigma)
}

if(method!="subplex"){
estim <- optim(
                   par=starting,
                   fn = function (par) {
                     lik.shift(
                               alpha=sym.par(par[seq_len(nalpha)]),
                               sigma=sym.par(par[nalpha+seq_len(nsigma)]),
                               sig=sym.par(par[nsig+seq_len(nsigma)]),
                               error=error,
                               dat=dat,
                               vcvList=vcvList
                               )
                   },
                   gr=NULL,
                   hessian=TRUE,
                   method=method,
                   control=control
                   )
}else{
estim <- subplex(
                   par=starting,
                   fn = function (par) {
                     lik.shift(
                               alpha=sym.par(par[seq_len(nalpha)]),
                               sigma=sym.par(par[nalpha+seq_len(nsigma)]),
                               sig=sym.par(par[nsig+seq_len(nsigma)]),
                               error=error,
                               dat=dat,
                               vcvList=vcvList
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
    sig=sym.par(estimML[nsig+seq_len(nsigma)])
  	N<-length(dat) 
    eig<-eigen(alpha)

    Vou<-.Call("simmap_covar",nterm=as.integer(n),bt=vcvList[[before]],lambda=eig$values,S=eig$vectors,sigma.sq=sigma)
    V<-Vou+kronecker(sig,vcvList[[after]])
  
   # add measurement error
   if(!is.null(error)){
			diag(V)<-diag(V)+error
		}
    
		theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%(t(W)%*%pseudoinverse(V)%*%dat)
		
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

# free parameters
if(model=="release" || model=="ER"){
nparam=(2*p)+p
}else if(model=="radiate" || model=="RR"){
nparam=(3*p)+p
}
# maximum likelihood estimates of alpha and sigma
estim.alpha<-estim$par[seq_len(nalpha)]
estim.sigma<-estim$par[nalpha+seq_len(nsigma)]
if(model=="RR"){
estim.sig<-estim$par[nsig+seq_len(nsigma)]
}
# estimates standards errors
estim.se<-sqrt(diag(pseudoinverse(estim$hessian)))
estim.se.alpha<-estim.se[seq_len(nalpha)]
estim.se.sigma<-estim.se[nalpha+seq_len(nsigma)]
if(model=="RR"){
estim.se.sig<-estim.se[nsig+seq_len(nsigma)]
}
# alpha matrix
alpha.mat<-sym.par(estim.alpha)
# sigma matrix
sigma.mat<-sym.par(estim.sigma)
# sig matrix
if(model=="RR"){
sig.mat<-sym.par(estim.sig)
}
# AIC
AIC<- -2*LL+2*nparam
# AIC corrected
AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) #Hurvich et Tsai, 1995
# matrix of estimated theta values
theta.mat<-matrix(res.theta,1)

rownames(theta.mat)<-"theta"
colnames(theta.mat)<-colnames(data)

##-------------------Print results--------------------------------------------##
if(echo==TRUE){
cat("\n")
cat("Summary results for the",nmod,"model","\n")
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
if(model=="RR" || model=="radiate"){
cat("\n")
cat("ML sigma radiation values","\n")
print(sig.mat)
}
}

##------------------List results----------------------------------------------##
if(model=="ER" || model=="release"){
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha.mat=alpha.mat, sigma.mat=sigma.mat, alpha=estim.alpha, sigma=estim.sigma, alpha.se=estim.se.alpha, sigma.se=estim.se.sigma, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val )
}else if(model=="RR" || model=="radiate"){
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha.mat=alpha.mat, sigma.mat=sigma.mat, sig.mat=sig.mat, alpha=estim.alpha, sigma=estim.sigma, alpha.se=estim.se.alpha, sigma.se=estim.se.sigma, sig.se=estim.se.sig, convergence=estim$convergence, hessian=estim$hessian, hess.values=hess.val )
} 

}
} #fin de la fonction