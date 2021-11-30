#' @title  Generate Interim survival data based on Known Weibull control shape parameter
#'
#' @param beta_c numeric variable
#' @param beta_t numeric variable
#' @param median_c numeric variable
#' @param median_t numeric variable
#' @return Generate Interim survival data
#' @export
#'
#' @examples
#' sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = 0.25,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
#'
#' sample_final
#'
#' sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = 0.25,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
#' sample_int
#'
#' data_int <- interim_data(beta_c = 0.25,beta_t =0.2447551, median_c = 4,median_t = 6)
#'
#' \dontrun{
#' sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = 0.25,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
#'
#' sample_final
#'
#' sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = 0.25,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
#' sample_int
#'
#' data_int <- interim_data(beta_c = 0.25,beta_t =0.2447551, median_c = 4,median_t = 6)
#'
#' }
#'


interim_data <- function(beta_c,beta_t, median_c,median_t){
  if (median_c <= 0 | median_c >= 20) stop("Error:The value for MEDIAN TIME (months) of Control arm is not allowed. The MEDIAN should be greater than 0")
  if (median_t <= 0 | median_t >= 20) stop("Error:The value for MEDIAN TIME (months) of Treatment arm is not allowed. The MEDIAN should be greater than 0")
  if (beta_c <= 0 | beta_c >= 10)   stop("Error:The value beta is not allowed, The Control shape parameter should be between 0 and 10")
  if (beta_t <= 0 | beta_t >= 10)   stop("Error:The value beta is not allowed, The Treatment shape parameter should be between 0 and 10")

  theta_c <- function(beta_c){
    if(beta_c < 0){stop("positive shape parameter is required.")}
    value <- median_c/(log(2)^(1/beta_c)) # 4 is the median survival in control arm
    return(value)
  }

  data_c <- function(n,shape,scale){
    if(censor == 1){
      set.seed(123)
      time <- rweibull(n,shape = beta_c,scale = theta_c(beta_c))
      status <- rep(1,n)
      data.frame(time = time,status = status)
    } else {
      set.seed(111)
      time_old <- rweibull(n,shape = beta_c,scale = theta_c(beta_c))
      censor_theta <- theta_c(beta_c)*(censor/(1-censor))^(1/beta_c)
      set.seed(100)
      censor_dis <- rweibull(n,shape = beta_c,scale = censor_theta)
      time <- pmin(time_old,censor_dis)
      status <- ifelse(time == time_old,1,0)
      data.frame(time = time,status = status)
    }
    }
  censor = (sample_int$d_true_0*0.80) #event rate
  n = ceiling(sample_final$n_admin_arm0)
  #beta_c <- 0.75 # Assume follows true beta0 = 0.75 and event rate = 0.8
  theta_crt = theta_c(beta_c)
  x <- data_c(n,beta_c,theta_c(beta_c)) # around 54 events
  x$group <- rep("Control",length(x$time))
  #sum(x$status)
  # Treatment arm data generation
  # Lets try with median treatment effect is 13

  theta_t <- function(beta_t){
    if(beta_t < 0){stop("positive shape parameter is required.")}
    value <- median_t/(log(2)^(1/beta_t))
    return(value)
  }

  data_t <- function(n,shape,scale){
    if(censor == 1){
      set.seed(12345)
      time <- rweibull(n,shape = beta_t,scale = theta_t(beta_t))
      status <- rep(1,n)
      data.frame(time = time,status = status)
    } else {
      set.seed(999)
      time_old <- rweibull(n,shape = beta_t,scale = theta_t(beta_t))
      censor_theta <- theta_t(beta_t)*(censor/(1-censor))^(1/beta_t)
      set.seed(1999)
      censor_dis <- rweibull(n,shape = beta_t,scale = censor_theta)
      time <- pmin(time_old,censor_dis)
      status <- ifelse(time == time_old,1,0)
      data.frame(time = time,status = status)
    }
    }

  censor = (sample_int$d_true_1*0.80)
  n = ceiling(sample_final$n_admin_arm0)
  #beta_t <- 0.704
  theta_trt = theta_t(beta_t)
  y <- data_t(n,beta_t,theta(beta_t))
  y$group <- rep("Trt",length(y$time))
  #sum(y$status)

  data <- rbind(x,y)
  return(data)
}
