#' @title  Predicitve Power approach to sample size recalculation
#'
#' @param nstart Numeric values
#' @param nmax Numeric values
#' @param D Numeric values
#' @param beta0 Numeric values
#' @param scenario Numeric values between 0 and 4
#'
#' @return
#' @export
#'
#' @examples
#' search.best.n.pred(nstart = 60,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
#'
search.best.n.pred <- function(nstart,nmax,D,beta0,scenario){
  if (D <= 0 ) stop("Error:The value D is not allowed. Total event should be greater than 0")
  if (nstart <= 0 ) stop("Error:The value nstart is not allowed. nstart should be greater than 0")
  if (nmax <= 0 ) stop("Error:The value nmax is not allowed. nmax should be greater than 0")
  if (beta0 <= 0 | beta0 >= 10)   stop("Error:The value beta0 is not allowed, The Beta0 should be between 0 and 10")
  if (scenario <= 0 | scenario >= 5)   stop("Error:The value scenario is not allowed, The scenario should be between 1 and 4")

  y <- c(); k <- 0
  n1 <- nstart
  #nmax <- 2*n1
  n3 <- nmax
  n2 <- ceiling(mean(c(n1, n3)))
  y[n1] <- reestimate_pwr(alpha=0.05,d_k = n1,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$pred_power
  y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$pred_power
  y[n3] <- reestimate_pwr(alpha=0.05,d_k = n3,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$pred_power
  repeat{
    if(y[n3] < 0.90) {
      message('WARNING: C. power ', 0.90, ' can not be achieved. Increase nmax.')
      return(n3)
      break
    }

    if(all(y[c(n1, n2, n3)] >= 0.90)){
      return(n1)
      break

    } else if(all(y[c(n2, n3)] >= 0.90)){
      n3 <- n2
      n2 <- ceiling(mean(c(n1, n3)))
      y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$pred_power

    } else if(y[n3] >= 0.90){
      n1 <- n2
      n2 <- ceiling(mean(c(n1, n3)))
      y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$pred_power
    }

    k <- k + 1
    if(k > log2(nmax - nstart)) {
      return(n2)
      break
    }
  }
}



