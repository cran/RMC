\name{sim.chain}
\alias{sim.chain}
\alias{sim.chain}
\title{Simulate Markov chain data from a Markov model.}
\description{Simulate chained data from a simple reversible Markov model (see Foster et al 2009 for details). Simulates stationary and non-stationary data, the later is formed by defining the transition matrix as a function of covariates.}
\usage{
sim.chain( n.chains = 5, n.obs = rep( 100, n.chains), n.cats = 3, n.covars = 1, beta = NULL, gamma = NULL, X = NULL)
}
\arguments{
\item{ n.chains}{scalar specifying the number of chains to simulate. The default value is 5.}
\item{ n.obs}{vector of length n.chains. Indicates the number of observations per chain. Default is rep( 100, n.chains) for 5 chains of 100 observations.}
\item{ n.cats}{scalar specifying the number of categories in the chain. Default is 3 states.}
\item{ n.covars}{The number of covariates that affect the transition matrix. The constant must be considered to be one of these covariates. Default is n.covars=1 for constant term only (stationary chain).}
\item{ beta}{elements of the matrix of coefficients for the probability of moving into each state, used to partially define the transition matrices (also need gamma). These values are transformed via the additive logistic transform. The dimension of the beta matrix must be nrow=n.covars and ncol=n.cats with rows indexing the covariates and columns indexing the transformed probabilities (transformed from a simplex with the additive logistic transform). The first column of this matrix must be zero, reflecting that the last category's transition probability is constrained.}
\item{ gamma}{elements of the matrix of coefficients for the probability of moving from any particular state, used to partially define the transition matrices (also need beta). These parameters are transformed using the logistic transform. The dimension of the gamma matrix must be nrow=n.covars and ncol=n.cats with rows indexing covariates and columns indexing states.}
\item{ X}{The design matrix for the covariates. If NULL (the default) then the design matrix is filled with a constant and random (normal) variables. The design matrix must be numeric with all factors suitably converted into dummy variables and so on.}
}
\details{The observed chained categorical model is defined according to Foster et al (2009). The Markov process is assumed to be parameterised by two vectors, phi and pi. The phi parameters indicate the probability of moving from each state and the pi probabilities prescribe the probability of moving to each state given that a move will occur. This process is reversible if the parameters do not change within a chain. The probabilities are allowed to vary within a chain by specifying these two vectors of probabilities as functions of covariates (possibly index number).}
\value{The function returns a matrix with columns for the chain identification (values from 1 to n.chains), the simulated chained data, and a column for each of the scaled X matrix.}
\references{Foster, S.D., Bravington, M.V., Williams, A., Althaus, F, Laslett, G.M., and Kloser, R{}.J. (2008) Analysis and prediction of faunal distributions from video and multi-beam sonar data using Markov models. \emph{Environmetrics}, \bold{20}: 541-560.}
\author{Scott D. Foster}
\examples{
#simulate 3 chains with 4 states of length 1500. There is one (constant) covariate. The beta parameters are for the pi probabilities and the gamma parameters are for the phi probabilities. These are the same parameters used to generate dataEG1
chain <- sim.chain( n.chains=3, n.obs=rep( 1500, 3), n.cats=4, n.covars=1, beta=c(0,0.3,-0.3,0), gamma=c(0.5,0.2,1,0))
#Plot the first 100 observations of each chain
par( mfrow=c(3,1))
for( ii in 1:3)\{
 tmpDat <- chain[chain[,"chain"]==ii,]
 plot( 1:100, tmpDat[1:100,"state"], main=paste( "Chain", ii), ylab="State", xlab="Index", type='b', pch=20)
 abline(h=1:4, lty=2)
\}
}
\keyword{misc}
