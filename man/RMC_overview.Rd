\name{RMC}
\alias{RMC}
\title{Package Description}
\description{
RMC is a package that fits and predicts reversible Markov models with a particular parameterisation described in Foster et al (2009). The core work-horse of the estimation routines is the function \code{\link{RMC.mod}}, but also see \code{\link{MVfill}} for (single) imputing of associated chained covariates, and \code{\link{RMC.pred}} for prediction of the local area stationary distribution.

Also contained in the RMC package is a bunch of methods that provide graphical diagnostics for this class of models (see Foster and Bravington 2010). This is performed through the functions \code{\link{diagnos}} for calculation of residuals, \code{\link{diagnos.envel}} for calculation of residuals and simulation envelopes, and \code{\link{hrplot}} for subsequent plotting.
}
\references{
Foster, S.D., Bravington, M.V., Williams, A., Althaus, F, Laslett, G.M., and Kloser, R{}.J. (2008) Analysis and prediction of faunal distributions from video and multi-beam sonar data using Markov models. \emph{Environmetrics}, \bold{20}: 541-560.

Foster, S.D. and Bravington, M.V. (2009) Graphical Diagnostics for Markov Models for Categorical Data. \emph{Journal of Computational and Graphical Statistics}, to appear.
}
\author{Scott D. Foster}
\seealso{\code{\link{diagnos}}, \code{\link{diagnos.envel}}, \code{\link{examplesForDiagnostics}}, \code{\link{hrplot}}, \code{\link{MVfill}}, \code{\link{RMC.mod}}, \code{\link{RMC.pred}}.}
\keyword{misc}
