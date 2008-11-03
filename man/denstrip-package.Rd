\name{denstrip-package}
\alias{denstrip-package}
\docType{package}
\title{
  Overview of the denstrip package
}
\description{
  Graphical methods for compactly 
  illustrating and comparing distributions, particularly distributions arising from parameter
  estimation or prediction.
}
\details{
\code{\link{denstrip}} implements the \emph{density strip} for
illustrating a single univariate distribution.  The darkness of the
density strip at a point is proportional to the density at that point. 
A shortcut function \code{\link{denstrip.normal}} draws the strip for the
given normal distribution. 

\code{\link{densregion}} implements the \emph{density region}, which
illustrates the uncertainty surrounding a continuously-varying quantity
as a two-dimensional shaded region with darkness proportional to the
density.  There are shortcut functions \code{\link{densregion.normal}}
and \code{\link{densregion.survfit}} for computing and drawing the
region for normally-distributed predictions and survival curves,
respectively.

\code{\link{sectioned.density}} implements the \emph{sectioned density plots}
of Cohen and Cohen (2006).  These illustrate distributions using
occlusion and varying shading.  They were developed for
summarising data, but can also be used for illustrating known
distributions.

\code{\link{vwstrip}} can be used to draw \emph{varying-width strips} to
illustrate distributions, in a similar manner to the \emph{violin plot} for
summarising data. The width of the strip is proportional to the density.
A shortcut function \code{\link{vwstrip.normal}} draws the strip for the
given normal distribution. 

\code{\link{bpstrip}} adapts the \emph{box-percentile plot} to illustrate a
distribution instead of observed data.  This strip has width
proportional to the probability of a more extreme point.

\code{\link{cistrip}} implements the popular point and line figure for
illustrating point and interval estimates, for example from multiple
regression. 

These methods are discussed in more detail by Jackson (2008).

Each function is designed to add a graphic to an existing set of plot
axes. The plots can be added to either base graphics or \pkg{lattice}
panels.

}
\author{
Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
}
\references{
  Jackson, C. H. (2008) Displaying uncertainty with shading. \emph{The American Statistician}, 62(4):340-347.
  \url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip.pdf}
}
\keyword{ package }
