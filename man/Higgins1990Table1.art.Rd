% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Higgins1990Table1.art}
\alias{Higgins1990Table1.art}
\title{Aligned Rank Transformed Version of Higgins1990Table1}
\format{
A data frame with 36 observations on the following 10 variables.
\describe{
     \item{Subject}{a factor with levels \code{"s1"} .. \code{"s36"}}
     \item{Row}{a factor with levels \code{"r1"} .. \code{"r3"}}
     \item{Column}{a factor with levels \code{"c1"} .. \code{"c3"}}
     \item{Response}{a numeric vector}
     \item{aligned.Response..for.Row}{a numeric vector}
     \item{aligned.Response..for.Column}{a numeric vector}
     \item{aligned.Response..for.Row.Column}{a numeric vector}
     \item{ART.Response..for.Row}{a numeric vector}
     \item{ART.Response..for.Column}{a numeric vector}
     \item{ART.Response..for.Row.Column}{a numeric vector}
}
}
\source{
Wobbrock, J. O., Findlater, L., Gergle, D., and Higgins, J. J.
\emph{ARTool}. \url{https://depts.washington.edu/acelab/proj/art/}.
}
\description{
The ART version of \code{\link{Higgins1990Table1}} as produced by the
original ARTool, used to test the correctness of \code{\link{art}} output.
}
\seealso{
\code{\link{Higgins1990Table1}}, \code{\link{art}}.
}
\keyword{datasets}
\keyword{internal}
