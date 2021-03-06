% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reestim_cond.R
\name{reestimate_pwr}
\alias{reestimate_pwr}
\title{Sample size re-estimation method}
\usage{
reestimate_pwr(alpha, D, d_k, beta0, RT, scenario)
}
\arguments{
\item{alpha}{Numeric values}

\item{D}{Numeric values}

\item{d_k}{Numeric values}

\item{beta0}{Numeric values}

\item{RT}{Numeric values}

\item{scenario}{Values should be between 0 and 1}
}
\value{

}
\description{
Sample size re-estimation using the Conditional Power and Predicitve power approach
}
\examples{
reestimate_pwr(alpha = 0.05,D = 82,d_k = 53,beta0 = 0.459,RT = 1.78,scenario = 1)
cond_result <- reestimate_pwr(alpha = 0.05,D = 82,d_k = 53,beta0 = 0.459,RT = 1.78,scenario = 1)
cond_result
event_cond <- search.best.n.fixed(nstart = 150,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
event_trend <- search.best.n.trend(nstart = 60,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
event_pred <- search.best.n.pred(nstart = 60,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
required_cond <- (event_cond*2 - cond_result$d_k)
required_trend <- (event_trend*2 - cond_result$d_k)
required_pred <- (event_pred*2 - cond_result$d_k)
}
