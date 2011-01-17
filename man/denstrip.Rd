\name{denstrip}
\alias{denstrip}
\alias{panel.denstrip}
\title{Density strips}
\description{
  The density strip illustrates a univariate distribution as a shaded
  rectangular strip, whose darkness at a point is proportional to the
  probability density.  The strip is darkest at the maximum density and
  fades into the background at the minimum density.  It may be used to generalise
  the common point-and-line drawing of a point and interval estimate,
  by representing the entire posterior or predictive distribution of the
  estimate.  This function adds a density strip to an existing plot. 
}
\usage{
denstrip(x, dens, at, width, horiz=TRUE, colmax, colmin="white",
         scale=1, gamma=1, ticks=NULL, tlen=1.5, twd, mticks=NULL,
         mlen=1.5, mwd, lattice=FALSE, \dots)
panel.denstrip(...)
}
\arguments{
  \item{x}{ Either the vector of points at which the density is
    evaluated (if \code{dens} supplied), or a sample from the distribution
    (if \code{dens} not supplied). }
  \item{dens}{ Density at \code{x}. If \code{dens} is not supplied,
    the density of the sample \code{x} is estimated by kernel density
    estimation, using \code{\link[stats]{density}(x,\dots)}.}
  \item{at}{ Position of the centre of the strip on the y-axis (if
    \code{horiz=TRUE}) or the x-axis (if \code{horiz=FALSE}).  }
  \item{width}{ Thickness of the strip, that is, the length of its
    shorter dimension.  Defaults to 1/30 of the axis range. }
  \item{horiz}{ Draw the strip horizontally (\code{TRUE}) or vertically (\code{FALSE}).}
  \item{colmax}{ Colour at the maximum density, either as a built-in R
    colour name (one of \code{\link{colors}()}) or an RGB hex
    value. Defaults to \code{par("fg")} which is normally 
    \code{"black"}, or \code{"#000000"}. Or in \pkg{lattice}, defaults   
    to \code{trellis.par.get("add.line")$col}.}
  \item{colmin}{ Colour to shade the minimum density, likewise. 
    Defaults to "white".  If this is set to \code{"transparent"},
    and the current graphics device supports transparency (see \code{\link{rgb}}), then
    overlapping strips will merge smoothly. }
  \item{scale}{ Proportion of \code{colmax} to shade the maximum
    density, for example \code{scale=0.5} with \code{colmax="black"} for
    a mid-grey colour.  }
  \item{gamma}{ Gamma correction to apply to the colour palette.
    The default of 1 should give an approximate perception of
    darkness proportional to density, but this may need to be adjusted
    for different displays. Values of \code{gamma} greater than 1
    produce colours weighted towards the lighter end, and values of
    between 0 and 1 produce darker colours.
  }
  \item{ticks}{ Vector of \code{x}-positions on the strip to draw tick marks, or \code{NULL} for no ticks. }
  \item{tlen}{ Length of these tick marks relative to the strip width. }
  \item{twd}{ Line thickness of these marks (defaults to
    \code{\link{par}("lwd")}, or in \pkg{lattice}, to

    \code{\link{trellis.par.get}("add.line")$lwd*2}.). }
  \item{mticks}{ \code{x}-position to draw a thicker tick mark or tick
    marks (for example, at the mean or median). }
  \item{mlen}{ Length of this mark relative to the strip width. }
  \item{mwd}{ Line thickness of this mark (defaults to
    \code{\link{par}("lwd")*2}, or in \pkg{lattice}, to    
    \code{\link{trellis.par.get}("add.line")$lwd*2}.). }
  \item{lattice}{Set this to \code{TRUE} to make \code{\link{denstrip}}
    a \pkg{lattice} panel function instead of a base graphics function. \cr
    \code{panel.denstrip(x,...)} is equivalent to
    \code{denstrip(x, lattice=TRUE, ...)}.}
  \item{\dots}{Additional arguments supplied to \code{\link[stats]{density}(x,\dots)}, if
    the density is being estimated.  For example, \code{bw} to change
    the bandwidth of the kernel. }
}
\references{ Jackson, C. H.  (2008) Displaying uncertainty with
  shading. \emph{The American Statistician}, 62(4):340-347. 
  \url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip.pdf}

  An add-on which enables the WinBUGS 1.4 software for Bayesian analysis to
  draw density strips is available from
\url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip_wbpatch.txt}. Open
  this file in WinBUGS and select \code{Tools->Decode->Decode All}.  

}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\seealso{ \code{\link{denstrip.legend}}, \code{\link{densregion}}.}
\examples{
## Illustrate a known standard normal distribution
## Various settings to change the look of the plot

x <- seq(-4, 4, length=10000)
dens <- dnorm(x)
plot(x, xlim=c(-5, 5), ylim=c(-5, 5), xlab="x", ylab="x", type="n")
denstrip(x, dens, at=0) # default width 
denstrip(x, dens, width=0.5, at=0) 
denstrip(x, dens, at=-4, ticks=c(-2, 0, 2)) 
denstrip(x, dens, at=-3, ticks=c(-2, 2), mticks=0) 
denstrip(x, dens, at=-2, ticks=c(-2, 2), mticks=0, mlen=3, 
         mwd=4, colmax="#55AABB") 
denstrip(x, dens, at=1, ticks=c(-2, 2), tlen=3, twd=3) 
denstrip(x, dens, at=-4, ticks=c(-2, 2), mticks=0, colmax="darkgreen", 
         horiz=FALSE)
x <- rnorm(1000) # Estimate the density
denstrip(x, width=0.2, at=-3, ticks=c(-2, 2), mticks=0, colmax="darkgreen",
         horiz=FALSE)
denstrip(x, at=2, width=0.5, gamma=2.2) 
denstrip(x, at=3, width=0.5, gamma=1/2.2) 

### Specifying colour of minimum density 
par(bg="lightyellow")
plot(x, xlim=c(-5, 5), ylim=c(-5, 5), xlab="x", ylab="x", type="n")
x <- seq(-4, 4, length=10000)
dens <- dnorm(x)
## Equivalent ways of drawing same distribution 
denstrip(x, dens, at=-1, ticks=c(-2, 2), mticks=0, colmax="darkmagenta")
denstrip(x, dens, at=-2, ticks=c(-2, 2), mticks=0, colmax="darkmagenta",
         colmin="lightyellow")
## ...though the next only works if graphics device supports transparency 
denstrip(x, dens, at=-3, ticks=c(-2, 2), mticks=0, colmax="darkmagenta",
         colmin="transparent")
denstrip(x, dens, at=-4, ticks=c(-2, 2), mticks=0, colmax="#8B008B", colmin="white")

## Alternative to density regions (\link{densregion.survfit}) for
## survival curves - a series of vertical density strips with no
## interpolation

library(survival)
fit <- survfit(Surv(time, status) ~ 1, data=aml, conf.type="log-log")
plot(fit, col=0)
lse <- (log(-log(fit$surv)) - log(-log(fit$upper)))/qnorm(0.975)
n <- length(fit$time)
lstrip <- fit$time - (fit$time-c(0,fit$time[1:(n-1)])) / 2
rstrip <- fit$time + (c(fit$time[2:n], fit$time[n])-fit$time) / 2
for (i in 1:n) { 
    y <- exp(-exp(qnorm(seq(0,1,length=1000)[-c(1,1000)], 
                        log(-log(fit$surv))[i], lse[i])))
    z <- dnorm(log(-log(y)), log(-log(fit$surv))[i], lse[i])
    denstrip(y, z, at=(lstrip[i]+rstrip[i])/2,
                 width=rstrip[i]-lstrip[i],
                 horiz=FALSE, colmax="darkred")
}
par(new=TRUE)
plot(fit, lwd=2)

## Use for lattice graphics (first example from help(xyplot))

library(lattice)
Depth <- equal.count(quakes$depth, number=8, overlap=.1)
xyplot(lat ~ long | Depth, data = quakes,
       panel = function(x, y) { 
           panel.xyplot(x, y)
           panel.denstrip(x, horiz=TRUE, at=-10, ticks=mean(x))
           panel.denstrip(y, horiz=FALSE, at=165, ticks=mean(y))
       } 
       )

## Lattice example data: heights of singing voice types

bwplot(voice.part ~ height, data=singer, xlab="Height (inches)",
       panel=panel.violin, xlim=c(50,80))
bwplot(voice.part ~ height, data=singer, xlab="Height (inches)",
       panel = function(x, y) {
           xlist <- split(x, factor(y))
           for (i in seq(along=xlist))
               panel.denstrip(x=xlist[[i]], at=i)
       },
       xlim=c(50,80)
       )
}



\keyword{aplot}
\keyword{color}
