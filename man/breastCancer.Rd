\name{breastCancer}
\alias{breastCancer}
\docType{data}
\title{Semi-Synthetic Interval-Censored Breast Cancer Data}
\description{
A semi-synthetic dataset based on the German Breast Cancer Study Group (GBSG) 
trial of node-positive breast cancer patients. The data follow an illness-death 
model with three states: "disease-free" (state 1), "recurrence" (intermediate 
transient state, state 2), and "death" (state 3). The first gap time (transition 
from disease-free to recurrence) was reformulated as interval-censored, to mimic 
the periodic clinical follow-up visits at which recurrence is detected. Times to 
death or censoring are given as exact (possibly right-censored) times.
}
\usage{
breastCancer
}
\format{
A data frame with 686 rows and 22 columns:
\describe{
  \item{age}{Age of the patient in years.}
  \item{menopause}{Menopausal status: 1 = premenopausal, 2 = postmenopausal.}
  \item{hormone}{Hormonal therapy: 1 = no, 2 = yes.}
  \item{size}{Tumour size, in millimetres.}
  \item{grade}{Tumour grade (1, 2 or 3).}
  \item{nodes}{Number of positive lymph nodes.}
  \item{prog_recp}{Progesterone receptor status (fmol).}
  \item{estrg_recp}{Estrogen receptor status (fmol).}
  \item{left_rec}{Left endpoint of the censoring interval for recurrence (in days).}
  \item{right_rec}{Right endpoint of the censoring interval for recurrence (in days). \code{NA} if the recurrence time is right-censored.}
  \item{rectime}{Recorded time to recurrence, or censoring time if no recurrence was observed (in days).}
  \item{censrec}{Recurrence indicator: 1 if recurrence occurred, 0 otherwise.}
  \item{survtime}{Total survival time since surgery (in days). These times are exact or right-censored.}
  \item{censdead}{Death indicator: 1 if death occurred, 0 otherwise.}
  \item{t23}{Time between recurrence and death or censoring (in days); 0 for patients without recurrence.}
  \item{X1to1}{Indicator of remaining in the disease-free state until the end of follow-up.}
  \item{X2to3}{Indicator of a transition from recurrence to death.}
  \item{X2to2}{Indicator of remaining in the recurrence state until the end of follow-up.}
  \item{right_rec_inf}{Same as \code{right_rec}, with \code{Inf} in place of \code{NA}.}
  \item{Z.mid}{Midpoint of the recurrence interval, used for midpoint imputation. For right-censored observations, it equals the censoring time.}
  \item{left}{Left endpoint of the recurrence interval used in the analysis (in days).}
  \item{right}{Right endpoint of the recurrence interval used in the analysis (in days); \code{Inf} if right-censored.}
}
}
\source{
Soutinho, G. et al. (2026). Proposal of transition probability estimators 
based on a landmarking strategy extending de Uña-Álvarez and Meira-Machado (2015). 
Based on the German Breast Cancer Study Group data (\code{gbcsCS} in the 
\pkg{condSURV} package), originally reported in Schumacher, M. et al. (1994), 
\emph{Journal of Clinical Oncology}, \strong{12}(10), 2086--2093.
}
\keyword{datasets}