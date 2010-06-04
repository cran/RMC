\name{MVfill}
\alias{MVfill}
\alias{MVfill}
\title{Fill in missing values via a single imputation from the fitted model.}
\description{Inserts the fitted probability of observing the chained data used in estimation for the model (fm). To be used when the outcome variable is used as a covariate at subsequent stages of an analysis. The fitted probability is conditional on the nearest observed value in the chain.}
\usage{
MVfill( fm, states=NULL, chain.id=NULL, X=NULL)
}
\arguments{
\item{ fm}{a fitted model for the outcome that is required to be filled. Must be the result of a call to \code{\link{RMC.mod}}}
\item{ states}{the outcome vector used to estimate the fitted model. Must be sequentially ordered within chains}
\item{ chain.id}{a vector indicating which states belong to which chains}
\item{ X}{the design matrix for the model used to create fm. Column ordering must match that provided to fm and row ordering must match that in states (and chain.id)}
}
\value{\item{ a matrix with number of columns equal to the number of categories of the outcome variable. Each column contains the observed value of that variable (if available) or the predicted probability of that observed variable (if no observation). Note that missing values at the ends of the chains are not imputed and are removed prior to estimation from RMC.mod.}{}}
\references{Foster, S.D., Bravington, M.V., Williams, A., Althaus, F, Laslett, G.M., and Kloser, R{}.J. (2008) Analysis and prediction of faunal distributions from video and multi-beam sonar data using Markov models. \emph{Environmetrics}, \bold{20}: 541-560.}
\author{Scott D. Foster}
\examples{
#form data with NAs in the outcomes
dataEG2.NA <- dataEG2
dataEG2.NA[sample(1:nrow( dataEG2.NA), 1000),"state"] <- NA	#specify lots of NAs
#fit a model using dataEG2.NA
fm.est2.NA <- RMC.mod( states=dataEG2.NA[,2], chain.id=dataEG2.NA[,1], X=dataEG2.NA[,-(1:2)])
#impute the missing values. The result will be a matrix
imputedData <- MVfill( fm.est2.NA, states=dataEG2.NA[,2], chain.id=dataEG2.NA[,1], X=dataEG2.NA[,3:4])
print( head( imputedData, 20))
}
\keyword{misc}
