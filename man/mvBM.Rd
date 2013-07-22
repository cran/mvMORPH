\name{mvBM}
\alias{mvBM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Brownian Motion models of continuous traits evolution
%%  ~~function to do ... ~~
}
\description{
This function allows the fitting of multivariate multiple rates of evolution. This function also allows to fit constrained models. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvBM(tree, data, error = NULL, reg = NULL, model = c("BMM", "BM1"), 
 constraint = FALSE, simmap.tree = TRUE, scale.height = FALSE, 
 method = c("L-BFGS-B","Nelder-Mead","subplex"), random.start=FALSE,
 control=list(maxit=20000),mod="ER",pseudoinverse=FALSE, diagnostic = TRUE,
 echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
  Phylogenetic tree in SIMMAP format by default.
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
  Matrix or data frame with species in rows and continuous traits in columns (preferentially with names and in the same order than in the tree).
%%     ~~Describe \code{data} here~~
}
  \item{error}{
   Matrix or data frame with species in rows and continuous traits standard error (squared) in columns.
%%     ~~Describe \code{error} here~~
}
  \item{reg}{
  Vector of current species regimes. Use this argument only if you use non-SIMMAP tree.
%%     ~~Describe \code{reg} here~~
}
  \item{model}{
  "BMM" for a multi-rate and multi-selective regimes, and "BM1" for a unique rate of evolution per trait.
%%     ~~Describe \code{model} here~~
}
  \item{constraint}{
  If TRUE, the constrained model will be also calculated (see details).
%%     ~~Describe \code{constraint} here~~
}
  \item{simmap.tree}{
  Set to FALSE if you specify the reg argument. In this case, ancestral states are estimated from the posterior probabilities of loglikelihood estimates using the rerooting Method of the "phytools" package.
%%     ~~Describe \code{simmap.tree} here~~
}
  \item{scale.height}{
  Scale the tree to relative length.
%%     ~~Describe \code{scaleHeight} here~~
}
  \item{method}{
  Methods used by the optimization function. (see ?optim and ?subplex for details).
%%     ~~Describe \code{method} here~~
}
\item{random.start}{
  Whether the starting values of the optimizer must be estimated or random. (assessing convergence with random starting values may be helpful).
%%     ~~Describe \code{method} here~~
}
  \item{control}{
  Max. bound for the number of iteration of the optimizer; other options can be fixed on the list (see ?optim or ?subplex).
%%     ~~Describe \code{control} here~~
}
  \item{mod}{
  Rates of transition for ancestral states reconstruction if simmap.tree is set to FALSE (see ?ace for details). 
%%     ~~Describe \code{mod} here~~
}
 \item{pseudoinverse}{
  Whether Moore-Penrose pseudoinverse should be used in calculation (slower).
%%     ~~Describe \code{pseudoinverse} here~~
}
  \item{diagnostic}{
  Whether the diagnostic on convergence and reliability of the estimates must be returned.
%%     ~~Describe \code{diagnostic} here~~
}
  \item{echo}{
  Whether the summary statistics must be returned (or just stored)
%%     ~~Describe \code{echo} here~~
}
}
\details{
The mvBM function fit a multivariate Brownian Motion (BM) process, with unique or multiple BM rates (see O'Meara et al., 2006; and Revell and Collar, 2009). Note that the function use the non-censored approach of O'Meara et al. (2006).
Constrained models force the sigma values for each studied trait to be the same using the constrained Cholesky decomposition proposed by Adams (2013), which allows considering the covariation between traits.
 This approach is extended here to the multi-rate case, by specifying that the rates must be the same in different parts of the tree (common selective regime).
 If the tree is in "phylo" format, the user may specify the selective regime for each tip species (in the same order as in the phylogeny) with the argument "reg". If no selective regimes are specified the function works only with the model "BM1".
 
 N.B. Mapping of ancestral states can be done using the "make.simmap", "make.era.map" or "paintSubTree" functions from the "phytools" package.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ 
\item{LogLik.m }{The log-likelihood of the optimal model (note that LogLik.cons return the Log-likelihood for the constrained model if the "constraint" argument is turned on).}
\item{AIC.mult }{Akaike Information Criterion for the optimal model (AIC.cons for the constrained model)}
\item{AICc.mult }{Sample size-corrected AIC (AICc.cons for the constrained model)}
\item{rates.m }{Evolutionary rates matrix for each selective regimes (rates.cons for the constrained model)}
\item{anc }{Estimated ancestral states (anc.cons for the constrained model)}
\item{convergence }{Convergence status of the optimizing function; "0" mean convergence (convergence.cons for the constrained model) (See ?optim for details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" mean that a reliable estimate has been reached (hess.values.cons for the constrained model). See ?mvOU for details.}
\item{LRT }{If the constraint argument is "TRUE", the LRT (Log-ratio test) score between the observed and constrained model is returned}
\item{pval }{If the constraint argument is "TRUE", the p-value of the LRT test (comparison with Chi-square distribution) is returned}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Adams D.C. 2013. Comparing evolutionary rates for different phenotypic traits on a phylogeny using likelihood. Syst. Biol. 62:181-192. 

O'Meara B.C., Ane C., Sanderson M.J., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution. Evolution. 60:922-933. 

Revell L.J. 2012. phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3:217-223. 

Revell L.J., Collar D.C. 2009. Phylogenetic analysis of the evolutionary correlation using likelihood. Evolution. 63:1090-1100. 




%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{optim}}
\code{\link{brownie.lite}}
\code{\link{evol.vcv}}
\code{\link{make.simmap}}
\code{\link{make.era.map}}
\code{\link{paintSubTree}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Toy Exemple
  set.seed(123)
  # Generating a random tree
  tree<-pbtree(n=25)

  # Setting the regime states of tip species
  sta<-as.vector(c(rep("Forest",10),rep("Savannah",15))); names(sta)<-tree$tip.label

  # Making the simmap tree with mapped states
  tree<-make.simmap(tree,sta , model="ER", nsim=1)
  col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

  # Plot of the phylogeny for illustration
  plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

  # 2 Random traits evolving along the phylogeny
  data<-data.frame(head.size=rTraitCont(tree), mouth.size=rTraitCont(tree))

  # Names of the species
  rownames(data)<-tree$tip.label

## Run the analysis!!
  # BM model with multiple regimes
  mvBM(tree,data)

  # BM model with a single rate per trait
  mvBM(tree,data, model="BM1")

  # Comparing evolutionary rates between traits with a constrained BM model (see details)
  # see Adams (2013)
  mvBM(tree,data, model="BM1", constraint=TRUE)
  
  mvBM(tree,data, model="BMM", constraint=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Brownian Motion }
\keyword{ Evolutionary rates }
\keyword{ Cholesky constraint }% __ONLY ONE__ keyword per line
