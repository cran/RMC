\name{diagnos}
\alias{diagnos}
\alias{diagnos}
\alias{diagnos.envel}
\title{Calculation of Markov residuals for discrete-time non-stationary Markov models with simple parameterisation.}
\description{Calculates patch and movement residuals for Markov models with a simple parameterisation. The models themselves are for categorical Markov processes that are usefully described by models whose parameterisation is based on a simple reversible Markov model and that can be extended to non-stationary cases. Non-stationary models are incorporated by letting the transition matrix vary with covariates. Simulation envelopes are created using diagnos.envel.}
\usage{
diagnos(obs.states, chain.id, X=NULL, fit)
diagnos.envel( obs.states, chain.id, X=NULL, fit, perc=c( 0.025, 0.975), B=100, contr=list( print.iter=50))
}
\arguments{
\item{ obs.states}{observed ordered chained data. If there are multiple chains then chains are stacked on top of each other. Argument must be supplied}
\item{ chain.id}{vector (length matches states) of identifiers for the individual chains. If NULL then it is assumed that all observations form a single chain.}
\item{ X}{design matrix as passed to the model fitting routine RMC.mod. If NULL then a matrix with 1 column full of ones is assumed.}
\item{ fit}{a fitted model formed by a call to the estimation function RMC.mod. Must match up with the X argument.}
\item{ perc}{the percentiles of the simulations to take for the simulation envelope}
\item{ B}{the number of simulations for the simulation envelope}
\item{ contr}{list describing control parameters for the function. Currently consists of a single value for how often printing is to be performed}
}
\details{
Calculates the patch and movement residuals for the given Markov model. The methods are given in Foster and Bravington (2009). The patch residuals assess the model's compatibility with the data by inspecting the probabilities of observing each fully observed patch. The movement residuals assess the model's ability to describe each of the movement categories in the transition matrix.

Usage of diagnos.envel will produce, in addition to the functionality of diagnos, simulation envelopes
}
\value{
\item{ *for diagnos* A list with the following elements}{}
\item{ patch}{a list with number of elements equal to the number of states. Each element contains the patch residuals for each state}
\item{ movement}{a square matrix containing the movement residuals from and to each state}
\item{ njumps}{a square matrix containing the number of jumps from each state to each other state}
\item{ *for diagnos.envel* a list with the following elements}{}
\item{ patch}{a list with elements equal to the number of states. Each list element contains a matrix with observed patch residuals, expected patch quantiles and, lower and upper simulation envelope bounds. All values are given on uniform and normal deviate scales}
\item{ movement}{a matrix with observed movement residuals, expected movement residuals and, lower and upper simulation envelope bounds. All values are given on uniform and normal deviate scales}
\item{ njumps}{a list with length equal to the number of states. Each element contains the number of observed and simulated movements from each a particular state to each other state}
}
\references{Foster, S.D. and Bravington, M.V. (2009) Graphical Diagnostics for Markov Models for Categorical Data. \emph{Journal of Computational and Graphical Statistics}, to appear.}
\author{Scott D. Foster}
\seealso{\code{\link{RMC.mod}} for estimation of the Markov models and \code{\link{hrplot}} to plot the calculated residuals from diagnos.envel.}
\examples{
#estimate a model for stationary example data, dataEG1
fm.est <- RMC.mod( states=dataEG1[,2], chain.id=dataEG1[,1], X=dataEG1[,3])
#calculate residuals
res <- diagnos( dataEG1[,2], dataEG1[,1], X=dataEG1[,3, drop=FALSE], fit=fm.est)
res.envel <- diagnos.envel( dataEG1[,2], dataEG1[,1], X=dataEG1[,3,drop=FALSE], fit=fm.est, B=25)
}
\keyword{misc}
