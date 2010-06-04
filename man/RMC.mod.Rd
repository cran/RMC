\name{RMC.mod}
\alias{RMC.mod}
\alias{RMC.mod}
\title{Estimation of categorical discrete-time non-stationary Markov chain models with simple parameterisation.}
\description{Estimation of categorical Markov chain models whose parameterisation is based on a simple reversible Markov model and that can be extended to non-stationary cases. The model is parameterised by two vectors of parameters: one describing the probability of moving from each state (phi) and the other describing the probability of moving into each state given that a movement will occur (pi). Non-stationary models are incorporated by letting each of these vectors depend on covariates.}
\usage{
RMC.mod( states, chain.id=NULL, X=NULL, phi.id=NULL, pi.id=NULL, vcov=FALSE, inits=NULL, contr=list( maxit=1000, epsg=1e-8, epsf=1e-8, epsx=1e-8), quiet=FALSE, penalty=0)
}
\arguments{
\item{ states}{observed ordered chained data. If there are multiple chains then chains are stacked on top of each other. Argument must be supplied}
\item{ chain.id}{vector (length matches states) of identifiers for the individual chains. If NULL then it is assumed that all observations form a single chain.}
\item{ X}{design matrix (covariates) for the two vectors of probabilities. If NULL then X is assumed to contain an intercept term only. If not NULL then model will depend on phi.id and pi.id matrices (see below). X must be of dimensions nrow(X)=length(states) and ncol(X)=number of covariates. Typically will be created with a call to model.matrix}
\item{ phi.id}{indicator matrix of zeros and ones showing which covariates to include in the model for which element of phi (zero means not included and one means included). Each element of phi corresponds to the probability of moving from an observed state. phi.id must be of dimensions nrow(phi.id)=ncol(X)=number of covariates and ncol(phi.id)=number of states. If NULL then all covariates are included. Covariates are included via a logistic model for each element of phi}
\item{ pi.id}{indicator vector of zeros and ones showing which covariates to include in the model for all elements of pi (zero means not included and one means included). Each element of pi correspond to the probability of moving to that state given that a movement will occur. pi.id must have length equal to the number of covariates and indicates if that covariate is included in the model for pi. Covariates are included via the additive logistic transformation (Aitchison 1982)}
\item{ vcov}{boolean indicating if the variance matrix of the parameter estimates should be calculated. TRUE indicates that it is calculated}
\item{ inits}{initial values for the parameters. Must be of appropriate length and ordered as phi parameters and the pi parameters. If NULL then initial values are assumed to be zero. The ordering of this vector is:phi parameters for category 1, category 2, etc followed by pi parameters for transformed category 1, transformed category 2, etc.}
\item{ contr}{list containing control values for the optimisation procedure. maxit specifies the maximum number of iterations before optimisation is stopped. epsg, epsf and epsx give the stopping tolerances for gradients, relative function and estimates respectively}
\item{ quiet}{boolean indicating if any output is wanted. TRUE indicates that output is generated}
\item{ penalty}{experimental argument for an optional quadratic penalty on the parameters. A non-zero value indicates that the sum of the squared parameters must be less than or equal to the value. A value of zero indicates no penalty and is the default.}
}
\details{
The observed chained categorical data (in argument states) is modelled according to that described in Foster et al (2009). The Markov process is assumed to be parameterised by two vectors, phi and pi. The phi parameters indicate the probability of moving from each state and the pi probabilities prescribe the probability of moving to each state given that a move will occur. This process is reversible if the parameters do not change within a chain. The probabilities are allowed to vary within a chain by specifying these two vectors of probabilities as functions of covariates (possibly index number).

Since the model has simple form then the stationary distribution is known (up to normalisation constant) and hence, the (log-)likelihood is calculated exactly.

Optimisation is performed using a quasi-Newton method implemented in the LBFGS code from the ALGLIB website (see references). First derivatives for the optimisation are obtained using automatic differentiation (Griewank 2001) using the CppAD tool for C++ (Bell 2007). This saves an awful lot of mucking around with derivative free methods and increases speed. If you do not already use automatic differentiation then you may want to look into it.
}
\value{
\item{ Upon successful completion the function returns}{}
\item{ pars}{the parameter estimates ordered as phi parameters and then pi parameters. The ordering of this vector is:phi parameters for category 1, category 2, etc followed by pi parameters for transformed category 1, transformed category 2, etc.}
\item{ like}{the maximised log-likelihood}
\item{ scores}{the gradients calculated at the estimates. Ordered to match the pars vector}
\item{ vcov}{the variance matrix of the estimates if vcov==TRUE and NULL if vcov==FALSE}
\item{ conv}{the convergence code from the quasi-Newton optimiser}
\item{ time}{the time taken to perform the fit}
\item{ niter}{the number of iterations required by the optimiser}
\item{ stuff}{quite literal:stuff used for model specification and optimisation. Generally not of use to the user}
}
\references{
Aitchison J. (1982) The statistical analysis of compositional data. \emph{The Journal of the Royal Statistical Society-series B} \bold{44}: 139-177.

ALGLIB http://www.alglib.net/ (accessed June 2008)

Bell BM. 2007. CppAD: a package for C++ algorithmic differentiation, COIN-OR. http://www.coin-or.org/CppAD/, version 2007/02/07.

Griewank A. (2001) Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation. \emph{SIAM}. Philadelphia.

Foster, S.D., Bravington, M.V., Williams, A., Althaus, F, Laslett, G.M., and Kloser, R{}.J. (2008) Analysis and prediction of faunal distributions from video and multi-beam sonar data using Markov models. \emph{Environmetrics}, \bold{20}: 541-560.
}
\author{Scott D. Foster}
\seealso{\code{\link{RMC.pred}} for predicting the stationary distribution at arbitrary combinations of covariates. \code{\link{diagnos}} and \code{\link{diagnos.envel}} for graphical diagnostic methods for models of class RMC.}
\examples{
#estimate a model for the stationary example data, dataEG1
fm.est1 <- RMC.mod( states=dataEG1[,2], chain.id=dataEG1[,1], X=dataEG1[,3])
#estimate a model for the non-stationary example data, dataEG2
fm.est2 <- RMC.mod( states=dataEG2[,2], chain.id=dataEG2[,1], X=dataEG2[,-(1:2)])
}
\keyword{misc}
