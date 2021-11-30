#' @title Generate Time to event data using the framework of Relative time method
#' @description Generate Survival data for Sample size re-estimation method using the framework of Relative time method
#' @param a numeric variable
#' @param f numeric variable
#' @param beta0 numeric variable
#' @param median_c numeric variable
#' @param p1 numeric variable
#' @param p2 numeric variable
#' @param RT_P1 numeric variable
#' @param RT_P2 numeric variable
#' @param scenario numeric variable between 1 and 4
#'
#' @return
#' @export
#'
#' @examples
#' Data(a = 12,f = 6,beta0 = 0.25,median_c = 4,p1 = 0.10,p2 = 0.90,RT_P1 = 1.52,RT_P2 = 1.98,scenario = 1)
#' set.seed(222)
#' interim <- Data(a = 12,f = 6,beta0 = 0.25,median_c = 4,p1 = 0.10,p2 = 0.90,RT_P1 = 1.52,RT_P2 = 1.98,scenario = 1)
#' sum(interim$new.event)
#' \dontrun{
#' set.seed(222)
#' interim <- Data(a = 12,f = 6,beta0 = 0.25,median_c = 4,p1 = 0.10,p2 = 0.90,RT_P1 = 1.52,RT_P2 = 1.98,scenario = 1)
#' sum(interim$new.event)
#' #' }
#'

Data <- function(a,f,beta0,median_c,p1,p2,RT_P1,RT_P2,scenario){

  if (a <= 0 ) stop("Error:The value a for ACCRUAL TIME (Months) is not allowed. The ACCRUAL TIME (Months) should be greater than 0")
  if (f <= 0 ) stop("Error:The value f for follow-up time (Months) is not allowed. The follow-up time (Months) should be greater than 0")

  if (beta0 <= 0 | beta0 >= 10)   stop("Error:The value beta0 is not allowed, The Beta0 should be between 0 and 10")
  if (median_c <= 0 | median_c >= 20) stop("Error:The value for MEDIAN TIME (months) of Control arm is not allowed. The MEDIAN should be greater than 0")

  if (p1 <= 0 ) stop("Error:The value p1 is not allowed. p1 should be greater than 0")
  if (p2 <= 0 ) stop("Error:The value p2 is not allowed. p2 should be greater than 0")

  if (p1 == 0.10 & scenario %in% c(3,4)) stop("Error:For this p1 value, Only scenario 1 and 2 applicable")
  if (p1 == 0.25 & scenario %in% c(1,2)) stop("Error:For this p1 value, Only scenario 3 and 4 applicable")

  if (RT_P1 <= 0 ) stop("Error:The value RT_P1 is not allowed. RT_P1 should be greater than 0")
  if (RT_P2 <= 0 ) stop("Error:The value RT_P2 is not allowed. RT_P2 should be greater than 0")

  if (RT_P1 == 1.52 & scenario %in% c(2,3,4)) stop("Error:For this RT_P1 value, Only scenario 1 applicable")
  if (RT_P1 == 2.00 & scenario %in% c(1,3,4)) stop("Error:For this RT_P1 value, Only scenario 2 applicable")

  if (RT_P1 == 1.50 & scenario %in% c(1,2,4)) stop("Error:For this RT_P1 value, Only scenario 3 applicable")
  if (RT_P1 == 1.667 & scenario %in% c(1,2,3)) stop("Error:For this RT_P1 value, Only scenario  4 applicable")

  if (scenario <= 0 | scenario >= 5)   stop("Error:The value scenario is not allowed, The scenario should be between 1 and 4")

  if (scenario == 1) {
  # Secnerio 1
  sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,
                            RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  } else if (scenario == 2){
  # # Secnerio 2
  sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 2.00,
                            RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  } else if (scenario == 3){
  # Secnerio 3
  sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.50,
                              RT_P2 = 1.667,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  } else if (scenario == 4){
  # Secnerio 4
  sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.667,
                              RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  }

  sigma0 = 1/beta0
  theta0 = median_c/(log(2)^sigma0)
  mu0 = log(theta0)
  t_p1 = qgamma(p1,1)
  t_p2 = qgamma(p2,1)
  g_p1 = log(t_p1)
  g_p2 = log(t_p2)
  mu1 = (log(RT_P1)*g_p2 - log(RT_P2)*g_p1)/(g_p2 - g_p1) + mu0
  sigma1 = (log(RT_P1) - log(RT_P2))/ (g_p1 - g_p2) +sigma0
  beta1 = 1/sigma1
  theta1 = exp(mu1)

  # Control  data
  n_admin_length <- (ceiling(sample_final$n_admin_arm0))
  u1 <- runif(n_admin_length,0,1)
  y = theta0*log(1/u1)^sigma0
  event <- rep(1,n_admin_length)
  arm <- rep("C",n_admin_length)
  a_time <- a*runif(n_admin_length,0,1)
  new.y <- ifelse(y+a_time > a+f, (a+f) - a_time,y)
  new.event <- ifelse(y+a_time > a+f, 0,1)
  data.h0 <- data.frame(new.y,new.event,arm)

  # Treatment  data
  n_admin_length <- (ceiling(sample_final$n_admin_arm1))
  #set.seed(11145)
  u2 <- runif(n_admin_length,0,1)
  y2 = theta1*log(1/u2)^sigma1
  event <- rep(1,n_admin_length)
  arm <- rep("T",n_admin_length)
  a_time.2 <- a*runif(n_admin_length,0,1)
  new.y <- ifelse(y2+a_time.2 > a+f, (a+f) - a_time.2,y2)
  new.event <- ifelse(y2+a_time.2 > a+f, 0,1)
  data.h1 <- data.frame(new.y,new.event,arm)
  data <- rbind(data.h0,data.h1)
  return(data)
}
