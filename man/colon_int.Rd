\name{colon_int}
\alias{colon_int}
\docType{data}
\title{Synthetic Interval-Censored Colon Cancer Data}
\description{
A synthetic dataset inspired by a well-known clinical trial on adjuvant chemotherapy 
for colon cancer (Moertel et al., 1990). The data simulates an illness-death model 
with three states: "disease-free", "recurrence" (intermediate transient state), and "death". 
To replicate real-world oncology protocols, the first gap time (transition from state 1 to 2, 
representing recurrence) was reformulated as interval-censored based on a front-loaded 
schedule of clinical follow-up visits. Times for total survival or transition to death 
(states 2-3 or 1-3) are given as exact, uncensored times.
}
\usage{
colon_int
}
\format{
A data frame with 929 rows and 8 columns:
\describe{
  \item{time1}{The true or original time for the first transition (recurrence).}
  \item{left.rec}{The left endpoint of the censoring interval for recurrence (in days).}
  \item{right.rec}{The right endpoint of the censoring interval for recurrence (in days). Uses \code{NA} or \code{Inf} if right-censored.}
  \item{Stime}{The total survival time from the start of the study (in days). These times are exact and not censored.}
  \item{event1}{An indicator variable for the first transition (recurrence): 1 if occurred, 0 otherwise.}
  \item{event}{The final survival status indicator: 1 if death occurred, 0 otherwise.}
  \item{time1.mean}{The mean or midpoint of the interval-censored recurrence time, often used in specialized estimation frameworks.}
  \item{id}{A unique identification number for each patient.}
  \item{rx}{Treatment factor: Levamisole (Lev), 5-FU + Levamisole (Lev+5FU), or Observation (Obs).}
  \item{sex}{Gender indicator: 1 for male, 0 for female.}
  \item{age}{Age of the patient in years at the start of the study.}
  \item{obstruct}{Obstruction of colon by tumor indicator: 1 if yes, 0 if no.}
  \item{perfor}{Perforation of colon indicator: 1 if yes, 0 if no.}
  \item{adhere}{Adherence to neighboring organs indicator: 1 if yes, 0 if no.}
  \item{nodes}{Number of lymph nodes with detectable cancer.}
  \item{differ}{Differentiation of tumor factor: 1 = well, 2 = moderate, 3 = poor.}
  \item{extent}{Extent of local spread factor: 1 = submucosa, 2 = muscle, 3 = serosa, 4 = contiguous structures.}
  \item{surg}{Time from surgery to randomization: 1 = short, 0 = long.}
  \item{node4}{More than 4 positive lymph nodes indicator: 1 if yes, 0 if no.}
}
}
\source{
Soutinho, G. et al. (2026). Proposal of transition probability estimators 
based on a landmarking strategy extending de Uña-Álvarez and Meira-Machado (2015). 
Inspired by the original clinical trial data from Moertel, C.G. et al. (1990), 
\emph{New England Journal of Medicine}, \strong{322}(6), 352-358.
}
\keyword{datasets}