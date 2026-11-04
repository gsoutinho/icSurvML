\name{brca}
\alias{brca}
\docType{data}
\title{Breast Cancer Cosmesis Data}
\description{
A dataset containing interval-censored times to cosmetic deterioration (breast
retraction) for breast cancer patients undergoing two different treatments.
}
\usage{
brca
}
\format{
A data frame with 94 rows and 4 columns:
\describe{
  \item{left}{The left endpoint of the censoring interval (in months).}
  \item{right}{The right endpoint of the censoring interval (in months). Uses \code{Inf} if right-censored.}
  \item{treatment}{A factor indicating the treatment group: \code{"Rad"} (Radiotherapy only) or \code{"RadChem"} (Radiotherapy and Chemotherapy).}
  \item{right2}{A supplementary time variable often used in specialized interval-censored frameworks.}
}
}
\source{
Finkelstein, D.M. and Wolfe, R.A. (1985). A semiparametric model for 
regression analysis of interval-censored failure time data. \emph{Biometrics}, \strong{41}, 933-945.
}
\keyword{datasets}