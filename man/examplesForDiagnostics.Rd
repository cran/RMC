\name{examplesForDiagnostics}
\alias{examplesForDiagnostics}
\alias{examplesForDiagnostics}
\title{The examples section of this help file provides code to produce the plots in Section 5.1 of Foster and Bravington (2009).}
\description{This help file provides code to generate illustrative, simulated examples for different types of data/model discrepancies. The examples are those given in Section 5.1 of Foster and Bravington (2009). The code is provided so that any interested reader of that paper can, if they want, reproduce the results.}
\references{Foster, S.D. and Bravington, M.V. (2009) Graphical Diagnostics for Markov Models for Categorical Data. \emph{Journal of Computational and Graphical Statistics}, to appear.}
\author{Scott D. Foster}
\seealso{\code{\link{RMC.mod}} for estimation of the Markov model, \code{\link{diagnos}} and \code{\link{diagnos.envel}} for calculation of residuals and simulation envelopes, \code{\link{hrplot}} for plotting the residuals, and \code{\link{sim.chain}} for simulating chained data.}
\examples{
## common parameters for the examples
nc <- 5
n <- 1000
n.cats <- 4
B<-50	#number of simulations for simulated envelopes. The paper uses 1000 but this can be pretty slow
####################################################
####	Example of stationary well-fitting models	####
####################################################
set.seed(21)
#simulating data -- will be the same as dataEG1
chain <- sim.chain( n.chains=nc, n.obs=rep( n, nc), n.cats=n.cats, n.covars=1, beta=c(0,0.3,-0.3,0), gamma=c(0.5,0.2,1,0))
#plotting start of first chain as an example (Figure 1 of paper)
m <- 100
plot( 1:m, head( chain[chain[,"chain"]==2,"state"],m), type='b', pch=19, main="Start of Example Chain", ylab="State", xlab="Index", axes=FALSE)
abline(h=c(1:n.cats), lty=3, col=grey(0.5))
axis(1)
axis( 2, 1:n.cats, 1:n.cats)
box()
#fitting the model
fm.est <- RMC.mod( states=chain[,2], chain.id=chain[,1], X=chain[,3])
#defining true model
fm <- fm.est
fm$pars <- c( 0.5,0.2,1,0,0.3,-0.3,0)
#generating simulation envelope
temp.est <- diagnos.envel( obs.states=chain[,2], chain.id=chain[,1], X=chain[,3,drop=FALSE], fit=fm.est, B=B)
#plotting patch residuals (Figure 2 of paper)
par( mfrow=c( 1,2))
my.cat <- 2
hrplot( temp.est[["patch"]][[my.cat]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Example Data, state 2", pch=19)
my.cat <- 3
hrplot( temp.est[["patch"]][[my.cat]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Example Data, state 3", pch=19)
#plotting movement residuals (Figure 3 of paper)
hrplot( temp.est$movement, ylab="Movement Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Example Data", pch=20)
#############################################################
####	Example of non-stationary good well-fitting models	####
#############################################################
set.seed( 10)
#simulating data -- will be the same as dataEG2
X <- cbind( rep( 1, nc*n), simRandWalk( nc=nc, ni=rep( n, nc), init.var=1, seq.var=0.1)[,-1])
colnames( X) <- c( "const", "rand2")
n.covs <- 2
gpar <- matrix( c( c( -0.6, 0), c( 1.1, 1.5), c( 0.2, 0), c( -1.25, -1.5)), nrow=n.covs, ncol=n.cats)
bpar <- matrix( c( rep( 0, n.covs), c( -0.9, 0.5), c( 0.8, -0.4), c( -0.4, -0.7)), nrow=n.covs, ncol=n.cats)
chain.ns <- sim.chain( n.chains=nc, n.obs=rep( n, nc), n.cats=n.cats, n.covars=n.covs, beta=bpar, gamma=gpar, X=X)
#setting up model
my.phi.id <- ifelse( gpar!=0, 1, 0)
my.pi.id <- apply( bpar, FUN=function(x)\{if( any( x!=0)) 1 else 0\}, MARG=1)
#fitting the model
fm.est1 <- RMC.mod( states=chain.ns[,2], chain.id=chain.ns[,1], X=chain.ns[,3:4], phi.id=my.phi.id, pi.id=my.pi.id)
#generating simulation envelope
temp1 <- diagnos.envel( chain.id=chain.ns[,1], obs.states=chain.ns[,2], X=chain.ns[,3:4], fit=fm.est1, B=B)
#plotting residuals (Figure 4 of paper)
par(mfrow=c(1,3))
my.cat <- 2
hrplot( temp1[["patch"]][[my.cat]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Non-Stationary Data, state 2", pch=20)
my.cat <- 3
hrplot( temp1[["patch"]][[my.cat]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Non-Stationary Data, state 3", pch=20)
hrplot( temp1$movement, ylab="Movement Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Non-Stationary Data", pch=20)
##################################
####		Adding  outliers		####
##################################
set.seed( 25)
#simulating data -- will be the same as dataEG3patch and dataEG3movement
chain.orig <- sim.chain( n.chains=nc, n.obs=rep( n, nc), n.cats=n.cats, n.covars=1, beta=c(0,0.3,-0.3,0), gamma=c(0.5,0.2,1,0))
chain1 <- chain2 <- chain.orig
chain1[301:320,"state"] <- 3
ids <- sample( setdiff( 1:(nc*n), seq( from = n, to = n*nc, by = n) ) , 100)
chain2[ ids,"state"] <- 3
chain2[ ids+1, "state"] <- 4
#fit the models
fm1 <- RMC.mod( states=chain1[,2], chain.id=chain1[,1], X=chain1[,3])
fm2 <- RMC.mod( states=chain2[,2], chain.id=chain2[,1], X=chain2[,3])
#generate simulation envelopes
temp1 <- diagnos.envel( chain.id=chain1[,1], obs.state=chain1[,2], X=chain1[,3,drop=FALSE], fm1, B=B)
temp2 <- diagnos.envel( chain.id=chain2[,1], obs.state=chain2[,2], X=chain2[,3,drop=FALSE], fm2, B=B)
#plotting residuals (Figure 5 of paper)
par(mfrow=c(1,2))
hrplot( temp1[["patch"]][[3]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Data with Patch Outlier", pch=20)
hrplot( temp2$movement, ylab="Movement Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Data with Repetition Outlier", pch=20)
##################################
####	Omitting Covariates		####
##################################
set.seed(13)
#simulating data -- will be the same as dataEG4
X <- cbind( rep( 1, nc*n), simRandWalk( nc=nc, ni=rep( n, nc), init.var=1, seq.var=0.1)[,-1])
colnames( X) <- c( "const", "rand2")
n.covs <- 2
gpar <- matrix( c( c( -0.6, 0), c( 1.1, 1.5), c( 0.2, 0), c( -1.25, -1.5)), nrow=n.covs, ncol=n.cats)
bpar <- matrix( c( rep( 0, n.covs), c( -0.9, 0.5), c( 0.8, -0.4), c( -0.4, -0.7)), nrow=n.covs, ncol=n.cats)
chain <- sim.chain( n.chains=nc, n.obs=rep( n, nc), n.cats=n.cats, n.covars=n.covs, beta=bpar, gamma=gpar, X=X)
#setting up model
my.phi.id <- ifelse( gpar!=0, 1, 0)
my.pi.id <- apply( bpar, FUN=function(x)\{if( any( x!=0)) 1 else 0\}, MARG=1)
#fit the models (correct and incorrect)
fm.est <- RMC.mod( states=chain[,2], chain.id=chain[,1], X=chain[,3,drop=FALSE])
fm.est1 <- RMC.mod( states=chain[,2], chain.id=chain[,1], X=chain[,3:4], phi.id=my.phi.id, pi.id=my.pi.id)
#generate simulation envelopes
temp <- diagnos.envel( chain.id=chain[,1], obs.states=chain[,2], X=chain[,3,drop=FALSE], fit=fm.est, B=B)
temp1 <- diagnos.envel( chain.id=chain[,1], obs.states=chain[,2], X=chain[,-(1:2)], fit=fm.est1, B=B)
#plotting residuals
par( mfrow=c( 2, 3))
my.set <- c( 1, 4)
hrplot( temp[["patch"]][[my.set[1]]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Mean Only Model", pch=20)
hrplot( temp[["patch"]][[my.set[2]]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Mean Only Model", pch=20)
hrplot( temp$movement, ylab="Movement Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Mean Only Model", pch=20)
	
hrplot( temp1[["patch"]][[my.set[1]]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Correct Model", pch=20)
hrplot( temp1[["patch"]][[my.set[2]]], ylab="Patch Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Correct Model", pch=20)
hrplot( temp1$movement, ylab="Movement Residuals - Normal Quantiles", xlab="Normal Quantiles", main="Correct Model", pch=20)
}
\keyword{misc}
