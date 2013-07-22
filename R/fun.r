
##------------------------Fonctions necessaires-------------------------------##

# Calcul d'un matrice positive semi-definite (Choleski transform) from OUCH by A. King
sym.par<-function (x) {
  nchar <- floor(sqrt(2*length(x)))
  if (nchar*(nchar+1)!=2*length(x)) {
    stop("a symmetric matrix is parameterized by a triangular number of parameters",call.=FALSE)
  }
  y <- matrix(0,nchar,nchar)
  y[lower.tri(y,diag=TRUE)] <- x
  y%*%t(y)
}

# Inverse of Sym.par, return a vector (from OUCH) by A. King
sym.unpar <- function (x) {
  y <- t(chol(x))
  y[lower.tri(y,diag=TRUE)]
}

  
# Compute matrix time to common ancestor
mTime<-function(phy,scale.height){
  if(is.ultrametric(phy)){
   if(scale.height==TRUE){
   mdist<-vcv.phylo(phy)/max(vcv.phylo(phy))
   }else{
    mdist<-vcv.phylo(phy)}
    mSpdist<-rep(max(mdist),length(phy$tip.label))
  }else{
   vcv<-vcv.phylo(phy)
   if(scale.height==TRUE){
   vstand<-vcv/max(vcv)
   }else{vstand<-vcv}
   mSpdist<-diag(vstand)
   mdist<-vstand
  }
list(mSpDist=mSpdist, mDist=mdist)
}

# Binary regime coding
regimeList<-function(mm,k){
nReg=length(mm)
regime <- matrix(0,nrow=nReg,ncol=k)
for(i in 1:nReg){
 regime[i,mm[i]]<-1
 }
return(regime)
}
# Set root regime
indiceReg<-function(n,indice){
for(i in 1:n){
  val=length(indice[[i]])
  indice[[i]][val+1]<-indice[[i]][val]
  }
  return(indice)
}



##-----------------Function used in mvBM--------------------------------------##

# Compute phylogenetically corrected mean vector (Freckleton, 2012: Methods in Ecology and Evolution)          
estMean <- function(y, V) {
	iV <- solve(V, tol = .Machine$double.eps)
	xdum <- matrix(rep(1, dim(y)[1] ) , nrow = dim(y)[1])
	xVix <- crossprod(xdum, iV %*% xdum)
	xViy <- crossprod(xdum, iV %*% y)
	mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy # This is a bad thing!
	return(mu)
	}

# Compute phylogenetically corrected variance matrix for traits (Freckleton, 2012)            
estVarBM <- function(x, V) {
	mu <- estMean(x, V)
	iV <- solve(V, tol = .Machine$double.eps)
	e <- x
	for( i in 1:length(mu) ) e[,i] <- e[,i]  - mu[i]
	s2 <- crossprod(e, iV%*%e)
	n <- dim(x)[1]
	k <- dim(x)[2]
	return(s2 / n)
	}

# Constrainded Choleski decomposition (Adams, 2013: Systematic Biology)                   
build.chol<-function(b,p){
 c.mat<-matrix(0,nrow=p,ncol=p)
 c.mat[lower.tri(c.mat)] <- b[-1]
 c.mat[p,p]<-exp(b[1])
 c.mat[1,1]<-sqrt(sum((c.mat[p,])^2))
 if(p>2){
 for (i in 2:(p-1)){
 c.mat[i,i]<-ifelse((c.mat[1,1]^2-sum((c.mat[i,])^2))>0,sqrt(c.mat[1,1]^2-sum((c.mat[i,])^2)), 0)
 }}
 return(c.mat)
}

##----------------------Functions_for_EB--------------------------------------##

# Modified function from old version of Geiger package 1.3 for tree transformations
ebTree<-function (phy, endRate=NULL, a=NULL) 
{
    if(is.ultrametric(phy)) {	
    	times <- branching.times(phy);
    } 
    if(!is.ultrametric(phy)) {
    	times <- BranchingTimesFossil(phy)[1:phy$Nnode];
    }
   
 
 	if(a==0) return(phy)   
    names(times) <- (as.numeric(names(times)))
    for (i in 1:length(phy$edge.length)) {
        bl <- phy$edge.length[i]
        age = times[which(names(times) == phy$edge[i, 1])]
        t1 = max(times) - age
        t2 = t1+bl
        phy$edge.length[i] = (exp(a*t2)-exp(a*t1))/(a)
    }
    phy
}

#While loop to escape error message , maxIter fixed at 10
multiTry <- function (expr, silent = TRUE, maxIter = 10,
    quotedExpr = substitute(expr),  envir = parent.frame())
{
    while ((maxIter <- maxIter - 1) >= 0 && inherits(tmp <- try(eval(quotedExpr,
        envir = envir), silent = silent), "try-error")) {
    }
    if (maxIter < 0) {
        stop("reached maxIter iterations, last error message is ",
            as.character(tmp))
    }
    tmp
}

# Branching time for non-ultrametric Tree, taken from Slater 2013: Methods in Ecology and Evolution (more efficient than the previous version)
BranchingTimesFossil <- function (phy, tol = .Machine$double.eps^0.5) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    phy2 <- phy
    phy <- new2old.phylo(phy)
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(as.numeric(phy$edge[, 2]) == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    bt <- abs(xx - max(xx));
	
	for(i in 1:length(bt)) {
		
		if(bt[i]<.Machine$double.eps^0.5) bt[i] <- 0; 	}
	
	names(bt) <- c(seq(nb.tip+1, nb.tip+nb.node), phy$tip.label)
	
	
	return(bt);
}
