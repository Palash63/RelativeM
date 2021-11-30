
#' @title Sample size re-estimation method
#' @description Sample size re-estimation using the Conditional Power and Predicitve power approach
#'
#' @param alpha Numeric values
#' @param D Numeric values
#' @param d_k Numeric values
#' @param beta0 Numeric values
#' @param RT Numeric values
#' @param scenario Values should be between 0 and 1
#'
#' @return
#' @export
#'
#' @examples
#' reestimate_pwr(alpha = 0.05,D = 82,d_k = 53,beta0 = 0.459,RT = 1.78,scenario = 1)
#' cond_result <- reestimate_pwr(alpha = 0.05,D = 82,d_k = 53,beta0 = 0.459,RT = 1.78,scenario = 1)
#' cond_result
#' event_cond <- search.best.n.fixed(nstart = 150,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
#' event_trend <- search.best.n.trend(nstart = 60,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
#' event_pred <- search.best.n.pred(nstart = 60,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
#' required_cond <- (event_cond*2 - cond_result$d_k)
#' required_trend <- (event_trend*2 - cond_result$d_k)
#' required_pred <- (event_pred*2 - cond_result$d_k)


reestimate_pwr <- function(alpha,D,d_k,beta0,RT,scenario){
  if (alpha <= 0 | alpha >= 1)   stop("Error:The value alpha is not allowed, The alpha should be between 0 and 1")
  if (D <= 0 ) stop("Error:The value D is not allowed. Total event should be greater than 0")
  if (d_k <= 0 ) stop("Error:The value d_k is not allowed. Interim should be greater than 0")
  if (RT <= 0 ) stop("Error:The value RT is not allowed. RT should be greater than 0")
  if (beta0 <= 0 | beta0 >= 10)   stop("Error:The value beta0 is not allowed, The Beta0 should be between 0 and 10")
  if (scenario <= 0 | scenario >= 5)   stop("Error:The value scenario is not allowed, The scenario should be between 1 and 4")

  D <- D/2
  d_k <- d_k/2
  z_alpha <-  qnorm(alpha,lower.tail = F)
  theta <- log(RT)
  if(scenario == 1){
  sample_fixed <- sample_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,
                              RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)

  sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,
                            RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)

  } else if (scenario == 2){
    sample_fixed <- sample_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 2.00,
                                RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)

    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 2.00,
                              RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)

  } else if (scenario == 3){
    sample_fixed <- sample_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.50,
                                RT_P2 = 1.667,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)

    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.50,
                              RT_P2 = 1.667,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)

  } else if (scenario == 4){
    sample_fixed <- sample_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.667,
                                RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)

    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.667,
                              RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)

  }
  i_k <- 1/ (1/(d_k*(sample_int$beta0)^2) + 1/(d_k*(sample_int$beta1)^2) )
  I_k <- 1/ (1/(D*(sample_fixed$beta0)^2) + 1/(D*(sample_fixed$beta1)^2) )
  if(I_k >= i_k) {
    z_k <- log(sample_int$RT_mid)/sample_int$pool_sd
    #z_k = 0.10 # Lower Z_k is better: In this case
    # Upper One sided conditional power
    lambda = i_k/I_k
    b.int = z_k*sqrt(lambda)
    # CP: under current trend
    theta_trend = b.int/lambda
    theta  = log(1.78)
    cond_power <- 1- pnorm(( z_k*sqrt(i_k) - z_alpha*sqrt(I_k) + theta*(I_k-i_k))/(sqrt(I_k-i_k)))
    cond_power_trend <- 1- pnorm(( z_k*sqrt(i_k) - z_alpha*sqrt(I_k) + theta_trend*(I_k-i_k))/(sqrt(I_k-i_k)))
    Futility.cond = 1- cond_power
    pred_power <- 1-pnorm(((z_k*sqrt(I_k)) - z_alpha*sqrt(i_k))/sqrt(I_k - i_k))
    res <- data.frame(2*d_k,cond_power, cond_power_trend,pred_power,Futility.cond)
    names(res) <- c("d_k","cond_power", "cond_power_trend","pred_power","Futility.cond")
    return(res)
  } else {
    stop("I_k became less than i_k i.e You need to choose lower event")
  }

}


#'  Conditional Power with fixed effect  approach to sample size recalculation
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
#' search.best.n.fixed(nstart = 60,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
#'

search.best.n.fixed <- function(nstart,nmax,D,beta0,scenario){
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
  #scenario = scenario
  y[n1] <- reestimate_pwr(alpha=0.05,d_k = n1,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power
  y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power
  y[n3] <- reestimate_pwr(alpha=0.05,d_k = n3,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power
  repeat{
    if(y[n3] < 0.80) {
      message('WARNING: C. power ', 0.80, ' can not be achieved. Increase nmax.')
      return(n3)
      break
    }

    if(all(y[c(n1, n2, n3)] >= 0.80)){
      return(n1)
      break

    } else if(all(y[c(n2, n3)] >= 0.80)){
      n3 <- n2
      n2 <- ceiling(mean(c(n1, n3)))
      y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power

    } else if(y[n3] >= 0.80){
      n1 <- n2
      n2 <- ceiling(mean(c(n1, n3)))
      y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power
    }

    k <- k + 1
    if(k > log2(nmax - nstart)) {
      return(n2)
      break
    }
  }
}



#'  Conditional Power with Current Trend  approach to sample size recalculation
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
#' search.best.n.trend(nstart = 60,nmax = 300,D = 82,beta0 = 0.459,scenario = 1)
#'
search.best.n.trend <- function(nstart,nmax,D,beta0,scenario){
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
  y[n1] <- reestimate_pwr(alpha=0.05,d_k = n1,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power_trend
  y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power_trend
  y[n3] <- reestimate_pwr(alpha=0.05,d_k = n3,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power_trend
  repeat{
    if(y[n3] < 0.80) {
      message('WARNING: C. power ', 0.80, ' can not be achieved. Increase nmax.')
      return(n3)
      break
    }

    if(all(y[c(n1, n2, n3)] >= 0.80)){
      return(n1)
      break

    } else if(all(y[c(n2, n3)] >= 0.80)){
      n3 <- n2
      n2 <- ceiling(mean(c(n1, n3)))
      y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power_trend

    } else if(y[n3] >= 0.80){
      n1 <- n2
      n2 <- ceiling(mean(c(n1, n3)))
      y[n2] <- reestimate_pwr(alpha=0.05,d_k = n2,D = D,beta0 = beta0,RT = 1.78,scenario = scenario)$cond_power_trend
    }

    k <- k + 1
    if(k > log2(nmax - nstart)) {
      return(n2)
      break
    }
  }
}


