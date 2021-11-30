#' @title Operating characteristics of the Relative time method
#' @description This function calculate the Type I error rate based on this method and associated Power under various scenarios.
#' @param num_sims numeric variable
#' @param a numeric variable
#' @param f numeric variable
#' @param Hyp Binary variable (0/1); 0 indicate the Type I error rate and 1 Indicate the associated Power
#' @param p1 numeric variable
#' @param p2 numeric variable
#' @param beta0 numeric variable
#' @param scenario numeric variable between 1 and 4
#'
#' @return
#' @import  plyr
#' @export
#'
#' @examples
#' num_sims <- 10000
#' sims <- plyr::ldply(1:num_sims, type.1.err, a = 12, f=12, Hyp=0,p1=0.1, p2=0.9,beta0 = 0.25,scenario = 1)
#' result <- apply(sims, MARGIN = 2 , FUN=mean, na.rm = TRUE)
#' result
#' \dontrun{
#' num_sims <- 10000
#' sims <- plyr::ldply(1:num_sims, type.1.err, a = 12, f=12, Hyp=0,p1=0.1, p2=0.9,beta0 = 0.25,scenario = 1)
#' result <- apply(sims, MARGIN = 2 , FUN=mean, na.rm = TRUE)
#' result
#' }
#'

type.1.err <- function(num_sims,a,f, Hyp,p1,p2,beta0,scenario){
  if (num_sims <= 0 ) stop("Error:The value num_sims is not allowed. The number of simulation should be greater than 0")
  if (a <= 0 ) stop("Error:The value a for ACCRUAL TIME (Months) is not allowed. The ACCRUAL TIME (Months) should be greater than 0")
  if (f <= 0 ) stop("Error:The value f for follow-up time (Months) is not allowed. The follow-up time (Months) should be greater than 0")
  if (p1 <= 0 ) stop("Error:The value p1 is not allowed. p1 should be greater than 0")
  if (p1 == 0.10 & scenario %in% c(3,4)) stop("Error:For this p1 value, Only scenario 3 and 4 applicable")
  if (p1 == 0.25 & scenario %in% c(1,2)) stop("Error:For this p1 value, Only scenario 1 and 2 applicable")
  if (p2 <= 0 ) stop("Error:The value p2 is not allowed. p2 should be greater than 0")
  if (beta0 <= 0 | beta0 >= 10)   stop("Error:The value beta0 is not allowed, The Beta0 should be between 0 and 10")
  if (scenario <= 0 | scenario >= 5)   stop("Error:The value scenario is not allowed, The scenario should be between 1 and 4")

  if (scenario == 1){
    # scenario : 1
    sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,
                                RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
  } else if(scenario == 2){
    #scenario : 2
    sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 2.00,
                                RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
  } else if(scenario == 3){
    # scenario : 3
    sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.50,
                                RT_P2 = 1.667,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
  } else if(scenario == 4){
    # # scenario : 4
    sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.667,
                                RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
  }

  # Under null they are the same Theta 1, theta0, beta1 and beta0
  if(Hyp == 0){
    theta1 <- sample_final$theta0
    sigma1 <- sample_final$sigma0
    theta0 <- sample_final$theta0
    sigma0 <- sample_final$sigma0
    n_admin_length <- (ceiling(sample_final$n_admin_arm0))
    #set.seed(123)
    u1 <- runif(n_admin_length,0,1)
    y = theta0*log(1/u1)^sigma0
    event <- rep(1,n_admin_length)
    arm <- rep("C",n_admin_length)
    a_time <- a*runif(n_admin_length,0,1)
    data.c <- data.frame(y,event,arm,a_time)
    new.y <- ifelse(y+a_time > a+f, (a+f) - a_time,y)
    new.event <- ifelse(y+a_time > a+f, 0,1)
    data.h0 <- data.frame(new.y,new.event,arm)
    n_admin_length <- (ceiling(sample_final$n_admin_arm1))

    #set.seed(111)
    u2 <- runif(n_admin_length,0,1)
    y2 = theta1*log(1/u2)^sigma1
    event <- rep(1,n_admin_length)
    arm <- rep("T",n_admin_length)
    a_time.2 <- a*runif(n_admin_length,0,1)
    new.y <- ifelse(y2+a_time.2 > a+f, (a+f) - a_time.2,y2)
    new.event <- ifelse(y2+a_time.2 > a+f, 0,1)
    data.h1 <- data.frame(new.y,new.event,arm)

  } else {
    theta1 <- sample_final$theta1
    sigma1 <- sample_final$sigma1
    theta0 <- sample_final$theta0
    sigma0 <- sample_final$sigma0
    n_admin_length <- (ceiling(sample_final$n_admin_arm0))

    #set.seed(12345)
    u1 <- runif(n_admin_length,0,1)
    y = theta0*log(1/u1)^sigma0
    event <- rep(1,n_admin_length)
    arm <- rep("C",n_admin_length)
    a_time <- a*runif(n_admin_length,0,1)
    new.y <- ifelse(y+a_time > a+f, (a+f) - a_time,y)
    new.event <- ifelse(y+a_time > a+f, 0,1)
    data.h0 <- data.frame(new.y,new.event,arm)
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
  }

  # Conrol arm : shape and scale parameter from the interim data
  fit <- survival::survreg(Surv(new.y,new.event) ~ 1, dist="weibull",data.h0)
  #summary(fit)
  Est_theta0 <- unname(exp(summary(fit)$coefficients[1])) # theta0: Scale
  MyBeta0 <- 1/summary(fit)$scale #beta0 # Shape
  Sig0 <- summary(fit)$scale
  wei_sd_c <- sqrt(summary(fit)$var[1])


  # Treatment arm
  fit1 <- survival::survreg(Surv(new.y,new.event) ~ 1, dist="weibull",data.h1)
  #summary(fit)
  Est_theta1 <- unname(exp(summary(fit1)$coefficients[1])) # theta1
  MyBeta1 <- 1/summary(fit1)$scale # beta1
  Sig1 <- summary(fit1)$scale
  wei_sd_t <- sqrt(summary(fit1)$var[1])
  #median_t = Theta1*((log(2))**Sigma1)


  # Now goal is to calculate the RT_mid from the interim data
  Theta_ratio = Est_theta1/Est_theta0
  p1_ratio = (log(1/p1))**(Sig1 - Sig0)
  p2_ratio = (log(1/p2))**(Sig1 - Sig0);
  p_mid_ratio = (log(2/(p1 + p2)))**(Sig1 - Sig0)

  exact_0 <- sample_final$d_true_0
  exact_1 <- sample_final$d_true_1
  exact_d <- sample_final$d_exact

  n_event_C <- sample_final$n_evt_arm0
  SD_pooled <- sample_final$pool_sd
  n0 <- sample_final$n_evt_arm0
  Check_mid <- sample_final$RT_mid

  Admin_ratio = (exact_0/exact_1)**(Sig1 - Sig0)

  TestStat1 = log(Theta_ratio*p1_ratio*Admin_ratio)
  TestStat2 = log(Theta_ratio*p2_ratio*Admin_ratio)
  TestStat_mid = log(Theta_ratio*p_mid_ratio*Admin_ratio)

  z1 = sqrt(n_event_C/exact_d)*TestStat1/SD_pooled
  z2 = sqrt(n_event_C/exact_d)*TestStat2/SD_pooled
  zmid = sqrt(n_event_C)*TestStat_mid/SD_pooled

  z1_count <- ifelse(z1 > 1.645,1,0)
  z2_count <- ifelse(z2 > 1.645,1,0)
  zmid_count <- ifelse(zmid > 1.645,1,0)

  Est_Omega0 = log((Sig0/Sig1)*(Est_theta0**(1/Sig0))/(Est_theta1**(1/Sig1)))
  Beta_diff = (1/Sig1) - (1/Sig0)

  Est_RT_midval = Theta_ratio*p_mid_ratio
  Est_LL_RT_midval = exp(log(Est_RT_midval) - 1.645*(SD_pooled/sqrt(n0)))
  Est_UL_RT_midval = exp(log(Est_RT_midval) + 1.645*(SD_pooled/sqrt(n0)))
  RelBias = (Est_RT_midval - Check_mid)/Check_mid

  RT_Coverage <- ifelse( (Check_mid > Est_LL_RT_midval) & (Check_mid < Est_UL_RT_midval),1,0)

  #return(c(zmid_count,RelBias,RT_Coverage))
  if(Hyp == 0){
    results <- c(zmid_count)
    names(results) <- "Type I error"
  } else if (Hyp == 1) {
    results <- c(zmid_count,RelBias,RT_Coverage)
  names(results) <- c("Power","Relative Bias","RT_Coverage")
  }
  return(results)
}


