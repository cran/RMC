\name{simRandWalk}
\alias{simRandWalk}
\alias{simRandWalk}
\title{Simulate a continuous auto-regressive process.}
\description{Simulates an autoregressive process for a specified number of sets of observations. The first observation in each chain is drawn from an independent normal and subsequent observations are drawn from another normal with mean equal to the last observation.}
\usage{
simRandWalk( nc=5, ni=rep( 1000, nc), init.var=1, seq.var=0.1)
}
\arguments{
\item{ nc}{the number of chains to simulate}
\item{ ni}{the length of each chain (must be a vector of length nc)}
\item{ init.var}{the variance to use for the initial random number generation (distribution will have mean zero)}
\item{ seq.var}{the variance to use in the sequential simulation}
}
\details{The autoregressive process is simulated for each of nc chains. The first element of each chain is simulated from a normal with zero mean and variance init.var. The subsequent values are simulated from a random draw from a normal with mean equal to the previous observation and variance seq.var.}
\value{A matrix with sum( ni) rows and 2 columns. The first column has elements 1:nc and indicates the chain the the observation belongs to. The second column contains the random values}
\author{Scott D. Foster}
\examples{
ni <- c( 30, 300, 3000)
simDat <- simRandWalk( nc=3, ni=ni, init.var=1, seq.var=0.1)
par( mfrow=c( 1, 3))
plot( 1:ni[1], simDat[1:ni[1],2], type='b', pch=20, ylab="Random Variable", xlab="Index", main="Chain 1")
plot( 1:ni[2], simDat[ni[1]+1:ni[2],2], type='b', pch=20, ylab="Random Variable", xlab="Index", main="Chain 2")
plot( 1:ni[3], simDat[sum(ni[1:2])+1:ni[3],2], type='b', pch=20, ylab="Random Variable", xlab="Index", main="Chain 3")
}
\keyword{misc}
