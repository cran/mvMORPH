\name{mvSHIFT}
\alias{mvSHIFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate change in mode of continuous traits evolution 
%%  ~~function to do ... ~~
}
\description{
This function allow the fitting of different models of evolution after a fixed point. This allows fitting model of change in mode of evolution following a given event.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvSHIFT(tree, data, age = NULL, error = NULL, sigma = NULL, alpha = NULL, 
 sig = NULL, model = c("ER", "RR", "EC", "SR"), scale.height = FALSE, 
 diagnostic = TRUE, method = c("L-BFGS-B", "Nelder-Mead", "subplex"), 
 pseudoinverse=FALSE, echo = TRUE, control=list(maxit=20000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
   Phylogenetic tree with shift mapped. (See "make.era.map" function from "phytools" package). A "phylo" object can be used if the "age" argument is provided.
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
   Matrix or data frame with species in rows and continuous traits in columns
%%     ~~Describe \code{data} here~~
}
  \item{age}{
   Age at which the shift in mode of evolution occur (in unit of time of the provided tree)
%%     ~~Describe \code{age} here~~
}
  \item{error}{
   Matrix or data frame with species in rows and continuous traits standard error (squared) in columns
  
%%     ~~Describe \code{error} here~~
}
  \item{sigma}{
   Starting values of the sigma matrix of the OU process before optimization.(optional)
%%     ~~Describe \code{sigma} here~~
}
  \item{alpha}{
   Starting values of the alpha matrix of the OU process before optimization.(optional)
%%     ~~Describe \code{alpha} here~~
}
  \item{sig}{
   Starting values of the sigma matrix of the BM process after or before the time shift (in "RR" or "EC" models, see details).(optional)
%%     ~~Describe \code{sig} here~~
}
  \item{model}{
   Choose between "RR" for ecological release and radiate model, "ER" for ecological release model, "EC" for ecologically constrained model, or "SR" for shift rate model (see details).
%%     ~~Describe \code{model} here~~
}
  \item{scale.height}{
   Whether the tree should be scaled to length 1
%%     ~~Describe \code{scale.height} here~~
}
  \item{diagnostic}{
  Whether the diagnostics of convergences should be returned
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
The mvSHIFT function fit a shift in mode or rate of evolution at a fixed point in time, as previously proposed by some authors (O'Meara et al. 2006; O'Meara, 2012; Slater, 2013). Shift in mode of evolution could be mapped on a modified "phylo" object using the "make.era.map" function from the "phytools" package.
Note that only one shift is allowed by the current version. The age of the shift could be otherwise directly provided in the function by the "age" argument in unit of times of the tree.

The function allows to fit model of "ecological release" and "ecological release and radiate" following Slater (2013), as well as a model of constrained ecology "EC" (e.g., after invasion of a competitive species in a given ecosystem) where traits are constrained in a Ornstein-Uhlenbeck process after a fixed point in time. The "SR" model allows fitting different (brownian) rates before and after the shift point (note that this model could also be fitted using the mvBM function).
The models "RR" for "radiate and release" or "ER" for "ecological release", fit an Ornstein-Uhlenbeck process before the fixed point while the drift parameter is not constrained after this point ("ecological release") or is allowed to vary ("release and radiate"). (see Slater, 2013).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik}{Log-Likelihood of the optimized model.}
\item{AIC}{Akaike Information Criterion calculated for the best model.}
\item{AICc}{Akaike Information Criterion corrected for small sample size.}
\item{theta.mat}{Matrix of estimated theta values for each traits and selective regimes (ancestral states).}
\item{alpha.mat}{Matrix of estimated alpha values (strength of selection) for studied traits (diagonal).}
\item{sigma.mat}{Matrix of estimated sigma values (drift) for studied traits (diagonal).}
\item{sig.mat}{Matrix of estimated sigma values for the BM process after the shift (diagonal, only in the "RR" model).}
\item{alpha}{alpha values estimated by the optimizing function.}
\item{sigma}{sigma values estimated by the optimizing function.}
\item{alpha.se}{standard errors of the alpha estimated values.}
\item{sigma.se}{standard errors of the sigma estimated values.}
\item{sig.se}{standard errors of the sigma estimated values for the BM process after the shift (only in "RR" model).}
\item{convergence}{Convergence statut of the optimizing function. See ?optim for more details.}
\item{hessian}{Hessian matrix (see "mvOU" details).}
\item{hess.value}{If value is 0, this mean that the eigen-values of the hessian matrix are all positives. Value of 1 means that the optimizing function may have converged to a non-reliable estimation.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
O'Meara B.C. 2012. Evolutionary inferences from phylogenies: a review of methods. Annu. Rev. Ecol. Evol. Syst. 43:267-285. 

O'Meara B.C., Ane C., Sanderson M.J., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution. Evolution. 60:922-933. 

Slater G.J. 2013. Phylogenetic evidence for a shift in the mode of mammalian body size evolution at the Cretaceous-Palaeogene boundary. Methods Ecol. Evol. 4:734-744. 

%% ~put references to the literature/web site here ~
}
\author{

Julien Clavel
%%  ~~who you are~~
}
\note{

Changes in rate of evolution and optima can also be fitted using the mvBM and mvOU functions using a 'make.era.map' transformed tree.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvOU}}
\code{\link{mvBM}}
\code{\link{mvEB}}
\code{\link{optim}}
\code{\link{subplex}}
\code{\link{paintSubTree}}
\code{\link{make.era.map}}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Toy Exemple
  set.seed(123)
  # Generating a random tree
  tree<-pbtree(n=25)

  # Making the simmap tree with the shift at a fixed point in time
  tot<-max(nodeHeights(tree))
  age=tot-0.3    # The shift occured 0.3 Ma ago
  tree<-make.era.map(tree,c(0,age))

  # Plot of the phylogeny for illustration
  plotSimmap(tree,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

  # 2 Random traits evolving along the phylogeny
  data<-data.frame(head.size=rTraitCont(tree), mouth.size=rTraitCont(tree))

  # Names of the species
  rownames(data)<-tree$tip.label

## Run the analysis!

  # Note that different rates before and after the shift could be fitted as:
  mvBM(tree,data)
  # or
  mvSHIFT(tree,data,0.3,model="SR")
  
  
  # Not run
  # Analysis with ecological release
  # mvSHIFT(tree,data,model="ER")

  # Analysis with ecological release and radiate
  # mvSHIFT(tree,data,model="RR")
  
 
   
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ornstein Uhlenbeck }
\keyword{ Shifts }
\keyword{ Brownian Motion }
\keyword{ Evolutionary rates }
\keyword{ RR }
\keyword{ EC }
\keyword{ SR }
\keyword{ ER }% __ONLY ONE__ keyword per line
