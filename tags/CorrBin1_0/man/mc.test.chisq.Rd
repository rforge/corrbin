\name{mc.test.chisq}
\alias{mc.test.chisq}
\title{Test the assumption of reproducibility}
\description{
 \code{mc.test.chisq} tests whether thee assumption of reproducibility is violated 
 in the data.}
\usage{
mc.test.chisq(cbdata)
}
\arguments{
  \item{cbdata}{a \code{\link{CBData}} object}
}
\details{
The assumption of marginal compatibility (AKA interpretability) implies that the marginal probability
of response does not depend on clustersize. Stefanescu and Turnbull (2003), and Pang and Kuk (2007) developed
a Cochran-Armitage type test for trend in the marginal probability of success as a function of the clustersize.
\code{mc.test.chisq} implements a generalization of that test extending it to multiple treatment groups. }
\value{
A list with the following components:
  \item{overall.chi }{the test statistic; sum of the statistics for each group}
  \item{overall.p }{p-value of the test}
  \item{individual}{a list of the results of the test applied to each group separately:
  \itemize{
  \item{chi.sq }{the test statistic for the group}
  \item{p }{p-value for the group}}}
}
\references{ 
Stefanescu, C. & Turnbull, B. W. (2003) Likelihood inference for exchangeable binary data with varying cluster sizes. \emph{Biometrics}, 59, 18-24

Pang, Z. & Kuk, A. (2007) Test of marginal compatibility and smoothing methods for exchangeable binary data with unequal cluster sizes. \emph{Biometrics}, 63, 218-227
}
\author{Aniko Szabo}
\seealso{\code{\link{CBData}} for constructing a CBData object, \code{\link{mc.est}} for estimating the distribution
under marginal compatibility. }
\examples{
data(shelltox)
mc.test.chisq(shelltox)
}
\keyword{htest}
