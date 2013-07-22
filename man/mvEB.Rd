\name{mvEB}
\alias{mvEB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Early Burst model of continuous traits evolution
%%  ~~function to do ... ~~
}
\description{
This function allows the fitting of multivariate Early Burst (EB) or ACDC model of evolution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvEB(tree, data, error = NULL, low = -3, up = 0, tol=c(0.00000001,Inf),
 scale.height = FALSE, control=list(maxit=20000), pseudoinverse=FALSE,
 diagnostic = TRUE, echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
   Phylogenetic tree (phylo object).
%%     ~~Describe \code{tree} here~~
}
   \item{data}{
   Matrix or data frame with species in rows and continuous traits in columns (preferentially with names and in the same order than in the tree).
%%     ~~Describe \code{data} here~~
}
  \item{error}{
   Matrix or data frame with species in rows and continuous traits standard error in columns.
%%     ~~Describe \code{error} here~~
}
  \item{low}{
   Lower bound of the exponent of the relationship between rate and time in the exponential change model.
%%     ~~Describe \code{low} here~~
}
  \item{up}{
   Upper (maximum) bound of the exponent of the relationship between rate and time in the exponential change model. Default value is zero for a decelerating model (see details).
%%     ~~Describe \code{up} here~~
}
  \item{tol}{
   Optimizer bounds for the search of the evolutionary rate matrix.
%%     ~~Describe \code{tol} here~~
}
  \item{scale.height}{
   Scale the tree to relative length.
%%     ~~Describe \code{scale.height} here~~
}
  \item{control}{
  Max. bound for the number of iteration of the optimizer; other options can be fixed on the list (see ?optim or ?subplex).
%%     ~~Describe \code{control} here~~
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
   Whether the summary statistics must be returned (or just stored).
%%     ~~Describe \code{echo} here~~
}
}
\details{
The Early Burst model (Harmon et al. 2010) is a special case of the ACDC model of Blomberg et al. (2003). Using an upper bound larger than zero transform the model to the accelerating rates of character evolution of Blomberg et al. (2003). The lower bound for decelerating rates (as assumed in Early Burst) could be fixed as log(min.rate) / T, where T is the depth of the tree and min.rate is the minimum rate that could be assumed for the model (as proposed by Slater and Pennell (in press))
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik }{The log-likelihood of the optimal model.}
\item{AIC }{Akaike Information Criterion for the optimal model.}
\item{AICc }{Sample size-corrected AIC.}
\item{r }{Exponent rate.}
\item{rates.m }{Evolutionary character rates matrix.}
\item{anc }{Estimated ancestral states.}
\item{convergence }{Convergence status of the optimizing function; "0" mean convergence (convergence.cons for the constrained model) (See ?optim for details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" mean that a reliable estimate has been reached (hess.values.cons for the constrained model). See ?mvOU for details.}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Blomberg S.P., Garland T.J., Ives A.R. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution. 57:717-745.

Harmon L.J., Losos J.B., Davies J.T., Gillespie R.G., Gittleman J.L., Jennings B.W., Kozak K.H., McPeek M.A., Moreno-Roark F., Near T.J., Purvis A., Ricklefs R.E., Schluter D., Schulte II J.A., Seehausen O., Sidlauskas B.L., Torres-Carvajal O., Weir J.T., Mooers A.O. 2010. Early bursts of body size and shape evolution are rare in comparative data. Evolution. 64:2385-2396.

Slater G.J., Pennell M. 2013. Robust regression and posterior predictive simulation increase power to detect early bursts of trait evolution. Syst. Biol. in press. 



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
\code{\link{mvBM}}
\code{\link{mvSHIFT}}
\code{\link{optim}}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Toy Exemple
  set.seed(123)
  
  # Generating a random tree
  tree<-pbtree(n=25)
  tree2<-tree
  
  # Transform branch length according to an early burst 
  # Not run!
  #library(geiger)
  #tree2=exponentialchangeTree(tree, a=-1)

  # 2 Random traits evolving along the phylogeny
  data<-data.frame(head.size=rTraitCont(tree2), mouth.size=rTraitCont(tree2))

## Run the analysis!!
  # EB model with two traits
  mvEB(tree,data)

  # Comparison with BM model
  mvBM(tree,data,model="BM1",simmap.tree=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Early Burst }
\keyword{ ACDC }
\keyword{ EB }
\keyword{ Evolutionary rates }% __ONLY ONE__ keyword per line
