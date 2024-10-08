\name{mvBM}
\alias{mvBM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Brownian Motion models of continuous traits evolution
%%  ~~function to do ... ~~
}
\description{
This function allows the fitting of multivariate multiple rates of evolution under a Brownian Motion model. This function can also fit constrained models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvBM(tree, data, error = NULL, model = c("BMM", "BM1"), 
    param = list(constraint = FALSE, smean = TRUE, trend=FALSE), 
    method = c("rpf", "pic", "sparse", "inverse", "pseudoinverse"),
    scale.height = FALSE, optimization = c("L-BFGS-B", "Nelder-Mead", "subplex"),
    control = list(maxit = 20000), precalc = NULL, diagnostic = TRUE, echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
Phylogenetic tree in SIMMAP format by default. A "phylo" object can also be used with the "BM1" model.
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns (preferentially with names and in the same order than in the tree). NA values are allowed with the "rpf", "inverse", and "pseudoinverse" methods.
%%     ~~Describe \code{data} here~~
}
  \item{error}{
Matrix or data frame with species in rows and continuous trait sampling variance (squared standard errors) in columns.
%%     ~~Describe \code{error} here~~
}
  \item{model}{
"BMM" for multi-rate and multi-selective regimes, and "BM1" for a unique rate of evolution per trait.
%%     ~~Describe \code{model} here~~
}
  \item{param}{
List of arguments to be passed to the function. See details.
%%     ~~Describe \code{param} here~~
}
  \item{method}{
Choose between "rpf", "sparse", "inverse", "pseudoinverse", or "pic" for log-likelihood computation during the fitting process. See details.
%%     ~~Describe \code{method} here~~
}
  \item{scale.height}{
Whether the tree should be scaled to unit length or not.
%%     ~~Describe \code{scale.height} here~~
}
  \item{optimization}{
Methods used by the optimization routines (see ?optim and ?subplex for details). The "fixed" method returns the log-likelihood function only.
%%     ~~Describe \code{optimization} here~~
}
  \item{control}{
Max. bound for the number of iteration of the optimizer; other options can be fixed in the list (see ?optim or ?subplex).
%%     ~~Describe \code{control} here~~
}
  \item{precalc}{
Optional. Precalculation of fixed parameters. See ?mvmorph.Precalc.
%%     ~~Describe \code{precalc} here~~
}
  \item{diagnostic}{
Whether the diagnostics of convergence should be returned or not.
%%     ~~Describe \code{diagnostic} here~~
}
  \item{echo}{
Whether the results must be returned or not.
%%     ~~Describe \code{echo} here~~
}
}
\details{
The mvBM function fits a homogeneous multivariate Brownian Motion (BM) process:
\deqn{dX(t) = \Sigma^{1/2}dW(t)}

With possibly multiple rates (\eqn{\Sigma_i}) in different parts ("i" selective regimes) of the tree (see O'Meara et al., 2006; Revell and Collar, 2009, see also details of the implementation in Clavel et al. 2015). Note that the function uses the non-censored approach of O'Meara et al. (2006) by default (i.e., a common ancestral state is assumed for the different regimes), but it is possible to specify multiple ancestral states (i.e., one for each regime) through the "smean" parameter (smean=FALSE) in the "param" list.

The "method" argument allows the user to try different algorithms for computing the log-likelihood. The \code{"rpf"} and \code{"sparse"} methods use fast GLS algorithms based on factorization for avoiding the computation of the inverse of the variance-covariance matrix and its determinant involved in the log-likelihood estimation. The \code{"inverse"} approach uses the "stable" standard explicit computation of the inverse and determinant of the matrix and is therefore slower. The \code{"pseudoinverse"} method uses a generalized inverse that is safer for matrix near singularity but highly time consuming. The \code{"pic"} method uses a very fast algorithm based on independent contrasts. It should be used with strictly dichotomic trees (i.e., no polytomies) and is currently not available for the multivariate "BMM" model. See ?mvLL for more details on these computational methods.


The \bold{"param"} \code{list} arguments:


\bold{"constraint"} - The "constraint" argument in the "param" list allows the user to compute the joint likelihood for each trait by assuming they evolved independently (\bold{ constraint="diagonal"}, or \bold{ constraint="equaldiagonal"}). If \bold{ constraint="equal"}, the sigma values are constrained to be the same for each studied trait using the constrained Cholesky decomposition proposed by Adams (2013) or a separation strategy based on spherical parameterization (when p>2) because of an unstable behavior observed for the constrained Cholesky (Clavel et al. 2015).

This approach is extended here to the multi-rate case by specifying that the rates must be the same in different parts of the tree (common selective regime). It's also possible to constraint the rate matrices in the "BMM" model to share the same eigen-vectors (\code{constraint="shared"}); the same variance but different covariances (\code{ constraint="variance"}); the same correlation but different variances (\code{ constraint="correlation"}); or to fit a model with different but proportional rates matrices (\code{constraint="proportional"}).

Finally, user-defined constrained models can be specified through a numeric matrix (square and symmetric) with integer values taken as indices of the parameters. For instance, for three traits:

\code{constraint=matrix(c(1,3,3,3,2,3,3,3,2),3)}. 

Covariances constrained to be zero are introduced by NA values, e.g., 

\code{constraint=matrix(c(1,4,4,4,2,NA,4,NA,3),3)}.

Difference between two nested fitted models can be assessed using the "LRT" function. See example below and ?LRT.


\bold{ "decomp"} - For the general case (unconstrained models), the sigma matrix is parameterized by various methods to ensure its positive definiteness (Pinheiro and Bates, 1996). These methods are the "cholesky", "eigen+", and "spherical" parameterizations.


\bold{"smean"} - Default set to TRUE. If FALSE, the ancestral state for each selective regime is estimated (e.g., Thomas et al., 2006).


\bold{"trend"} - Default set to FALSE. If TRUE, the ancestral state is allowed to drift linearly with time. This model is identifiable only with non-ultrametric trees. Note that it is possible to provide a vector of integer indices to constrain the estimated trends (see the vignettes).


\bold{"sigma"} - Starting values for the likelihood estimation. By default the theoretical expected values are used as starting values for the likelihood optimization (for measurement errors, multiple rates,...). The user can specify starting values with a list() object for the "BMM" model (e.g., two objects in the list for a two-regime analysis), or a simple vector of values for the "BM1" model. The parameterization is done using various factorizations for symmetric matrices (e.g., for the "decomp" argument; Pinheiro & Bates, 1996). Thus, you should provide p*(p+1)/2 values, with p the number of traits (e.g., random numbers or the values from the cholesky factor of a symmetric positive definite sigma matrix; see example below). If a constrained model is used, the number of starting values is (p*(p-1)/2)+1.

If no selective regime is specified the function works only with the model "BM1".


N.B.: Mapping of ancestral states can be done using the "make.simmap", "make.era.map" or "paintSubTree" functions from the "phytools" package.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik }{The log-likelihood of the optimal model.}
\item{AIC }{Akaike Information Criterion for the optimal model.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated ancestral states.}
\item{sigma }{Evolutionary rate matrix for each selective regime.}
\item{convergence }{Convergence status of the optimizing function; "0" indicates convergence (See ?optim for details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached. (See ?mvOU).}
\item{param }{List of model fit parameters (optimization, method, model, number of parameters...).}
\item{llik }{The log-likelihood function evaluated in the model fit "$llik(par, root.mle=TRUE)".}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Adams D.C. 2013. Comparing evolutionary rates for different phenotypic traits on a phylogeny using likelihood. Syst. Biol. 62:181-192.

Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods Ecol. Evol. 6(11):1311-1319.

O'Meara B.C., Ane C., Sanderson M.J., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution. Evolution. 60:922-933.

Revell L.J. 2012. phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3:217-223.

Revell L.J., Collar D.C. 2009. Phylogenetic analysis of the evolutionary correlation using likelihood. Evolution. 63:1090-1100.

Thomas G.H., Freckleton R.P., Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proc. R. Soc. B. 273:1619-1624.

%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}
\note{
The "pic" method is not yet implemented for the multivariate "BMM" model.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvgls}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{mvOUTS}}
\code{\link{mvRWTS}}
\code{\link{mvSIM}}
\code{\link{LRT}}
\code{\link{optim}}
\code{\link[phytools]{brownie.lite}}
\code{\link[phytools]{evol.vcv}}
\code{\link[phytools]{make.simmap}}
\code{\link[phytools]{make.era.map}}
\code{\link[phytools]{paintSubTree}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Simulate the traits
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(0,0)
data<-mvSIM(tree, param=list(sigma=sigma, ntraits=2, theta=theta,
            names_traits=c("head.size","mouth.size")), model="BM1", nsim=1)

## Fitting the models
# BMM - Analysis with multiple rates
mvBM(tree, data)

# BM1 - Analysis with a unique rate matrix
fit1<-mvBM(tree, data, model="BM1", method="pic")

# BM1 constrained
fit2<-mvBM(tree, data, model="BM1", method="pic", param=list(constraint="equal"))

# Comparison with LRT test
LRT(fit1,fit2)

# Random starting values
mvBM(tree, data, model="BMM", method="sparse", param=list(sigma=list(runif(3), runif(3))))

# Specified starting values (from the Cholesky factor)
chol_factor<-chol(sigma)
starting_values<-chol_factor[upper.tri(chol_factor,TRUE)]
mvBM(tree, data, model="BMM", method="sparse",
    param=list( sigma=list(starting_values, starting_values)))


# Multiple mean
mvBM(tree, data, model="BMM", method="sparse", param=list(smean=FALSE))


# Introduce some missing cases (NA values)
data2<-data
data2[8,2]<-NA
data2[25,1]<-NA

mvBM(tree, data2, model="BM1")


## FAST FOR THE UNIVARIATE CASE!!
\donttest{
 set.seed(14)
 tree2<-pbtree(n=5416) # Number of Mammal species
# Setting the regime states of tip species
 sta<-as.vector(c(rep("group_1",2000),rep("group_2",3416))); names(sta)<-tree2$tip.label

# Making the simmap tree with mapped states
 tree2<-make.simmap(tree2,sta , model="ER", nsim=1)
 col<-c("blue","orange"); names(col)<-c("Group_1","Group_2")
 plotSimmap(tree2,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Simulate a trait evolving by brownian motion on the tree
 trait<-rTraitCont(tree2)

# Fitting the models
 mvBM(tree2, trait, model="BMM", method="pic")
 mvBM(tree2, trait, model="BM1", method="pic")
}

}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Brownian Motion }
\keyword{ Evolutionary rates }
\keyword{ Cholesky constraint }% __ONLY ONE__ keyword per line
