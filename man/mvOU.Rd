\name{mvOU}
\alias{mvOU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Ornstein-Uhlenbeck model of continuous traits evolution
%%  ~~function to do ... ~~
}
\description{
This function allows the fitting of a multivariate OU model while allowing a given tree branch to be subdivided by multiple selective regimes using SIMMAP-like mapping of ancestral states. Species measurement error or dispersion can also be included in the model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvOU(tree, data, error = NULL, sigma = NULL, alpha = NULL, 
 model = c("OUM", "OU1"), scale.height=FALSE, diagnostic = TRUE, 
 method = c("L-BFGS-B","Nelder-Mead","subplex"), pseudoinverse=FALSE, 
 echo = TRUE, control=list(maxit=20000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
  Phylogenetic tree with mapped ancestral states in SIMMAP format. 
  (See make.simmap function from phytools package). A "phylo" object can be used with model "OU1".
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
  Matrix or data frame with species in rows and continuous traits in columns.
%%     ~~Describe \code{data} here~~
}
  \item{error}{
  Matrix or data frame with species in rows and continuous traits standard error (squared) in columns.
  
%%     ~~Describe \code{error} here~~
}
  \item{sigma}{
  Starting values of the sigma matrix before optimization (optional). Note that if you use "constraint" instead of NULL, sigma covariances matrix off-diagonal are omitted. You can also use a list argument to provide starting values for each traits (see detail).
%%     ~~Describe \code{sigma} here~~
}
  \item{alpha}{
  Starting values of the alpha matrix before optimization (optional). Note that if you use "constraint" instead of NULL, sigma covariances matrix off-diagonal are omitted. You can also use a list argument to provide starting values for each traits (see detail).
%%     ~~Describe \code{alpha} here~~
}
  \item{model}{
  Choose between "OUM" for a multiple selective regime model, or "OU1" for a unique selective regime for the whole tree.
%%     ~~Describe \code{model} here~~
}
 \item{scale.height}{
  Whether the tree should be scaled to length 1.
%%     ~~Describe \code{diagnostic} here~~
}
  \item{diagnostic}{
  Whether the diagnostics of convergence should be returned.
%%     ~~Describe \code{diagnostic} here~~
}
 \item{method}{
  Methods used by the optimization function. (see ?optim and ?subplex for details).
%%     ~~Describe \code{method} here~~
}
 \item{pseudoinverse}{
  Whether Moore-Penrose pseudoinverse should be used in calculation (slower)
%%     ~~Describe \code{pseudoinverse} here~~
}
  \item{echo}{
  Should the results be returned
%%     ~~Describe \code{echo} here~~
}
   \item{control}{
  Max. bound for the number of iteration of the optimizer; other options can be fixed on the list (see ?optim or ?subplex).
%%     ~~Describe \code{control} here~~
}
}
\details{
The mvOU function fit a multivariate model of evolution according to a Ornstein-Uhlenbeck process. This function use a modified version of the OUCH package (Butler & King, 2004) and extends it by allowing the user to incorporate measurement error and use SIMMAP-like mapping of ancestral states. SIMMAP mapping allows one to assign parts of branchs to different selective regimes, and allows testing for changes in trait variance which are not synchronous with the species divergence event.

Mapping of ancestral states can be done using the "make.simmap", "make.era.map" or "paintSubTree" functions from the "phytools" package.

It is possible to test for the significance of the off-diagonal sigma and alpha matrix in the full model by making comparison with a constrained model (using sigma="constraint", or alpha="constraint"). You can also provide starting values for the constrained model. For instance, for two traits use sigma=list("constraint", c(0.5,0.5))  (or alpha=list("constraint", c(0.5,0.5))).

The hessian is the matrix of second order partial derivatives of the likelihood function with respect to the maximum likelihood parameter values. This matrix provide a measure of the steepness of the likelihood surface in the region of the optimum. The eigen-decomposition of the hessian matrix returned by the optimizing function allow to assess the reliability of the fit of the model (even if the optimizer have converged).
When the optimization function does not converge on a stable result, the user may consider increasing the "maxit" argument in the "control" option, or try a simpler model with less parameters to estimate. Changing the starting values ("alpha" and "sigma" options) as well as the optimizing function ("method" options) may help sometimes (e.g., alpha=runif(3) for two traits analysis with random starting values (i.e., the lower triangular alpha matrix)).

Note that Bartoszek et al.,(2012) proposed the mvSLOUCH package dedicated to multivariate Ornstein-Uhlenbeck processes, which allows fitting models with a randomly evolving predictor variables. They also provide a detailled mathematical description of the multivariate OU model.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik}{Log-Likelihood of the optimized model.}
\item{AIC}{Akaike Information Criterion calculated for the best model.}
\item{AICc}{Akaike Information Criterion corrected for small sample size.}
\item{theta.mat}{Matrix of estimated theta values for each trait and selective regimes.}
\item{alpha.mat}{Matrix of estimated alpha values (strength of selection) for the studied traits (diagonal).}
\item{sigma.mat}{Matrix of estimated sigma values (drift) for the studied traits (diagonal).}
\item{alpha}{alpha values estimated by the optimizing function.}
\item{sigma}{sigma values estimated by the optimizing function.}
\item{alpha.se}{standard errors of the alpha estimated values.}
\item{sigma.se}{standard errors of the sigma estimated values.}
\item{convergence}{Convergence status of the optimizing function. See ?optim for more details.}
\item{hessian}{Hessian matrix (see details).}
\item{hess.value}{If value is 0, this mean that the eigen-values of the hessian matrix are all positives. Value of 1 means that the optimizing function may have converged to a non-reliable estimation.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bartoszek K., Pienaar J., Mostad P., Andersson S., Hansen T.F. 2012. A phylogenetic comparative method for studying multivariate adaptation. J. Theor. Biol. 314:204-215. 

Beaulieu J.M., Jhwueng D.-C., Boettiger C., O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution. 66:2369-2389. 

Butler M.A., King A.A. 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. Am. Nat. 164:683-695. 


%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}
\note{

This function use a modified version of the C code from the "OUCH" package build by Aaron King, as well as a C code part from the "ape" package by Emmanuel Paradis.  I kindly thank those authors for sharing their sources.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvBM}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{optim}}
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

## Run the analysis!
  # OU model with unique optima
  mvOU(tree,data, model="OU1", scale.height=TRUE)
  
  # Analysis with multiple optima!
  mvOU(tree,data, scale.height=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ornstein Uhlenbeck }
\keyword{ Hessian }
\keyword{ OU }% __ONLY ONE__ keyword per line
