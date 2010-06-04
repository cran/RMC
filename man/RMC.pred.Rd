\name{RMC.pred}
\alias{RMC.pred}
\alias{RMC.pred}
\title{Prediction of local area stationary distribution for an arbitrary number of covariate combinations.}
\description{Predicts the "local area stationary distribution" (with standard errors), and the average of all predictions (with standard errors). The latter prediction and the standard error can be taken to be an areal prediction if the number of points in the area is large and the covariance is not too high.}
\usage{
RMC.pred( fit, fit2=NULL, pts)
}
\arguments{
\item{ fit}{an object resulting from an RMC fit. This model must relate to the set of observations whose marginal distribution is required.}
\item{ fit2}{an object resulting from an RMC fit whose outcomes are covariates for the outcomes in the object fit. See Foster et al. (2009) for details on the process.}
\item{ pts}{a data matrix whose covariates must match those in the RMC objects fit and fit2. This matrix defines the points where the predictions are to occur.}
}
\details{
Predictions are made at each specified point by calculating the stationary distribution of the transition matrix for that combination of covariates. The variance, due to parameter uncertainty, of these prediction is obtained using the delta approximation. For predictions of a univariate model these derivatives are found using automatic differentiation (see Griewank 2001) implemented using the CppAD tool (Bell 2007). For predictions of a bivariate model the derivatives are found using finite differences and hence this procedure can run a little slow (well, a lot slow). If this causes problems then a more sophisticated implementation can be considered.

Global (or areal if the points come from a contiguous region of space) predictions are also calculated as the simple average of the predicted points. The variance of this average can be used as a prediction variance if the number of prediction points is large and the covariance between them is not too severe.

Predictions for univariate chained data are obtained by specifying a single model to the "fit" argument (keeping the "fit2" argument NULL). Predictions for the second variable of bivariate chained data are obtained by specifying the target univariate model in the "fit" argument, and the non-target univariate model in "fit2". This methodology was developed to predict fauna types where geomorphology (also chained data) was used as a covariate in the fauna model.

Details are given in Foster et al. 2009.
}
\value{
\item{ area}{the average of the predictions. Can be used as an areal prediction}
\item{ pts}{point predictions for each of the rows in the prediction matrix}
\item{ vcov}{variance and covariance matrix of the areal predictions. This is only calculated if the parameter estimate's variance for fit (and fit2 if not null) are present}
}
\references{
Bell BM. 2007. CppAD: a package for C++ algorithmic differentiation, COIN-OR. http://www.coin-or.org/CppAD/, version 2007/02/07.

Foster, S.D., Bravington, M.V., Williams, A., Althaus, F, Laslett, G.M., and Kloser, R{}.J. (2008) Analysis and prediction of faunal distributions from video and multi-beam sonar data using Markov models. \emph{Environmetrics}, \bold{20}: 541-560.

Griewank A. (2001) Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation. \emph{SIAM}. Philadelphia.
}
\author{Scott D. Foster}
\seealso{\code{\link{RMC.mod}} to estimate the Markov model, \code{\link{MVfill}} to impute any missing values in a particular sub-set of covariates, and \code{\link{sim.chain}} to simulate chained data.}
\examples{
#fit model to non-stationary data including all covariates
fm.est1 <- RMC.mod( states=dataEG2[,2], chain.id=dataEG2[,1], X=dataEG2[,3:4], vcov=TRUE)	#estimate the model
#perform predictions
pred1 <- RMC.pred( fit=fm.est1, fit2=NULL, pts=dataEG2[,3:4])
tmp <- cbind( table( dataEG2[,"state"]) / nrow( dataEG2), pred1$area, sqrt( diag( pred1$vcov)))
colnames( tmp) <- c( "observed", "predicted", "se")
rownames( tmp) <- paste( "cat",1:fm.est1$stuff$n.cats, sep="_")
print( tmp)
####Simulate and predict from bivariate non-stationary chained data
n.cats1 <- fm.est1$stuff$n.cats; n.cats2 <- 5; n.covars <- 2; n.covars2 <- n.covars*n.cats1 + n.covars
#fill in missing values (if any) for the previous outcomes that are now covariates
level1Probs <- MVfill( fm.est1, states=dataEG2[,2], chain.id=dataEG2[,1], X=dataEG2[,3:4])
#setting up design matrix -- note the order, it is important for RMC.pred
X2 <- matrix( rep( dataEG2[,3:4], n.cats1+1), nrow=nrow( dataEG2))
for( ii in 1:n.cats1)
  X2[,ii*n.covars+1:n.covars] <- X2[,ii*n.covars+1:n.covars] * rep( level1Probs[,ii], n.covars)
colnames(X2) <- paste( c( "const", "rand"), rep(c("",1:n.cats1), each=2), sep="")
#specify parameter values, note that parameterisation gives state 1 as a reference state for all categories. This is for convenience only.
gamma <- matrix( rnorm( n.covars2*n.cats2, sd=0.5), nrow=n.covars2, ncol=n.cats2)	#initial design matrix
for( ii in 1:n.cats2)
  gamma[sample( 2:n.covars2, sample(n.covars2-1:4, 1)), ii] <- 0
beta <- matrix( rnorm( n.covars2*n.cats2, sd=0.5), nrow=n.covars2, ncol=n.cats2)
beta[n.covars+1:n.covars,] <- 0
beta[,n.cats2] <- 0
#simulate chains
chains2 <- sim.chain( n.chains=5, n.obs=rep( 1000, 5), n.cats=n.cats2, n.covars=n.covars2, beta=beta, gamma=gamma, X=X2)	#simulate the chained categorical data
#specify RMC model to be estimated
my.phi.id <- ifelse( gamma!=0, 1, 0)	#model controlling matrix
my.pi.id <- apply( beta, FUN=function(x)\{if( any( x!=0)) 1 else 0\}, MARG=1)	#model controlling matrix
#fit model
fm.est2 <- RMC.mod( states=chains2[,2], chain.id=chains2[,1], X=X2, phi.id=my.phi.id, pi.id=my.pi.id, vcov=TRUE)	#estimate the model
#perform predictions
pred2cond <- RMC.pred( fit=fm.est2, fit2=NULL, pts=X2)
pred2marg <- RMC.pred( fit=fm.est2, fit2=fm.est1, pts=dataEG2[,3:4])
#check against empirical value
tmp <- cbind( table( chains2[,"state"]) / nrow( chains2), pred2cond$area, sqrt( diag( pred2cond$vcov)), pred2marg$area, sqrt( diag( pred2marg$vcov)))
colnames( tmp) <- c("observed","conditional prediction", "conditional se", "marginal prediction", "marginal se")
rownames( tmp) <- paste( "category", 1:n.cats2, sep="_")
print( tmp)
}
\keyword{misc}
