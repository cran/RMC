\name{hrplot}
\alias{hrplot}
\alias{hrplot}
\title{Plot horizontal residual plots for a RMC model.}
\description{Plot movement or a particular state's patch residuals for objects created by the estimation function RMC.mod.}
\usage{
hrplot( envel.obj, ylim, ...)
}
\arguments{
\item{ envel.obj}{a matrix of residuals and their upper and lower simulation limits. The easiest method to obtain this object is via a call to the diagnos.envel function. The object to be passed is an element of the resulting list. If the matrix has a column called "labels" then these are used as plot marks in the plot (useful for movement residuals).}
\item{ ylim}{a numeric vector of length 2 giving the minimum and maximum values for the y-axis. This parameter is very similar to that described in par but has slightly different implementation (very slight).}
\item{ ...}{other graphical parameters to be passed to the plot function.}
}
\details{This function is simple but care will need to be exercised in its use. See the example for best instruction and see Foster and Bravington (2009) for details about what the residuals are.}
\references{Foster, S.D. and Bravington, M.V. (2009) Graphical Diagnostics for Markov Models for Categorical Data. \emph{Journal of Computational and Graphical Statistics}, to appear.}
\author{Scott D. Foster}
\seealso{\code{\link{RMC}} to estimate the Markov model, and \code{\link{diagnos}} and \code{\link{diagnos.envel}} to calculate residuals and simulation envelopes.}
\examples{
#produces Figure 1 of Foster and Bravington (2009)
#fit RMC model to dataEG1
fm.est <- RMC.mod( states=dataEG1[,2], chain.id=dataEG1[,1], X=dataEG1[,3])
#calculate the simulation envelope (many fewer samples than Foster and Bravington (2009)
B<-50
temp.est <- diagnos.envel( obs.states=dataEG1[,2], chain.id=dataEG1[,1], X=dataEG1[,3,drop=FALSE], fit=fm.est, B=B)
#plot the patch residuals for the second state and the movement residuals
par(mfrow=c(1,2))
hrplot( temp.est[["patch"]][[2]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Stationary Data", pch=20)
hrplot( temp.est$movement, ylab="Movement Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Stationary Data", pch=20)
}
\keyword{misc}
