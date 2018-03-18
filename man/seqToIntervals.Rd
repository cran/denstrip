\name{seqToIntervals}
\alias{seqToIntervals}
\title{Find contiguous sequences in a vector of integers}
\description{
  Get all sequences of contiguous values in a vector of integers.
}
\usage{
seqToIntervals(x)
}
\arguments{
  \item{x}{A vector of integers, for example, representing
    indices. \code{x} is coerced to integer, sorted, and unique values
    extracted, if necessary, before finding the contiguous sequences. }
}
\value{
  A matrix with one row for each sequence, and two columns containing
  the start and the end of the sequence, respectively.
}
\seealso{
  \code{\link{sectioned.density}}
}
\examples{
seqToIntervals(1:10)  # [1 10]
seqToIntervals(c(1:10, 15:18, 20))  # [1 10; 15 18; 20 20]
# vectorised, so efficient for large vectors x
seqToIntervals(sample(1:1000000, size=999990)) 
}
\author{Chris Jackson <chris.jackson@mrc-bsu.cam.ac.uk>. Thanks to Tobias Verbeke for the suggestion.}
\keyword{manip}
\keyword{arith}
