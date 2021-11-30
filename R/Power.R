#' Title Stochastic curtailment test calculation
#'
#' @param alpha numeric variable
#' @param beta_c numeric variable
#' @param scenario numeric variable
#'
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
#' result = Interim_test(alpha = 0.05*0.90,beta_c = 0.25,scenario =1)
#' result
#'
#' \dontrun{
#' sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = 0.25,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
#'
#' sample_final
#'
#' sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = 0.25,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
#' sample_int
#'
#' data_int <- interim_data(beta_c = 0.25,beta_t =0.245, median_c = 4,median_t = 6) # Secnerio 2: They are different : RT_P1 = 1.52 AND RT_P2 = 1.98
#'
#' result = Interim_test(alpha = 0.05*0.90,beta_c = 0.25,scenario =1)
#' result
#' }
#'
#'
#'
Interim_test <- function(alpha,beta_c,scenario){

  if (alpha <= 0 | alpha >= 1)   stop("Error:The value alpha is not allowed, Alpha should be between 0 and 1")
  if (beta_c <= 0 | beta_c >= 10)   stop("Error:The value beta is not allowed, The Beta should be between 0 and 10")
  #if (beta_t <= 0 | beta_t >= 10)   stop("Error:The value beta is not allowed, The Beta0 should be between 0 and 10")
  if (scenario == 1){
  sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = beta_c,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta_c,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
  data_int <- interim_data(beta_c = beta_c,beta_t = sample_final$beta1, median_c = 4,median_t = 6) # Secnerio 2: They are different : RT_P1 = 1.52 AND RT_P2 = 1.98
  } else if (scenario == 2){
    sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = beta_c,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 2.00,RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta_c,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 2.00,RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
    data_int <- interim_data(beta_c = beta_c,beta_t = sample_final$beta1, median_c = 4,median_t = 6) # Secnerio 2: They are different : RT_P1 = 1.52 AND RT_P2 = 1.98

  }else if (scenario == 3){
    sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = beta_c,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.50,RT_P2 = 1.667,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta_c,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.50,RT_P2 = 1.667,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
    data_int <- interim_data(beta_c = beta_c,beta_t = sample_final$beta1, median_c = 4,median_t = 6) # Secnerio 2: They are different : RT_P1 = 1.52 AND RT_P2 = 1.98

  } else if (scenario == 4){
    sample_final <- sample_size(alpha = 0.05,sides = 1,beta0  = beta_c,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.667,RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta_c,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.667,RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 6,p_event=0.8)
    data_int <- interim_data(beta_c = beta_c,beta_t = sample_final$beta1, median_c = 4,median_t = 6) # Secnerio 2: They are different : RT_P1 = 1.52 AND RT_P2 = 1.98

  }
  #' sample_int
  #cat("\n Conditional power to reject null hypothesis H0: theta = 0")
  z_alpha <-  qnorm(alpha,lower.tail = F)
  beta_t = sample_final$beta1
  # Control arm : shape and scale parameter from the interim data
  library(survival)
  datax <- subset(data_int,group == "Control")
  fit <- survreg(Surv(time,status) ~ 1, dist="weibull",datax)
  #summary(fit)
  wei_scale_c <- unname(exp(summary(fit)$coefficients[1])) # theta0: Scale
  wei_shape_c <- 1/summary(fit)$scale #beta0 # Shape
  wei_sd_c <- sqrt(summary(fit)$var[1])


  # Treatment arm
  datay <- subset(data_int,group == "Trt")
  fit <- survreg(Surv(time,status) ~ 1, dist="weibull",datay)
  #summary(fit)

  wei_scale_t <- unname(exp(summary(fit)$coefficients[1])) # theta1
  wei_shape_t <- 1/summary(fit)$scale # beta1
  wei_sd_t <- sqrt(summary(fit)$var[1])
  #median_t = Theta1*((log(2))**Sigma1)


  # Now goal is to calculate the RT_mid from the interim data
  r=1

  mu_d = log(sample_final$RT_mid)
  SD = sqrt((1/beta_c^2 + 1/(r*beta_t^2))/sample_final$n_evt_arm0)
  #SD = sample_final$pool_sd
  #I = 1/(1/sample_final$n_evt_arm0 + 1/sample_final$n_evt_arm0) #23 FOE sECNERIO 2
  I = 1/(1/(beta_c*sample_final$n_evt_arm0) + 1/(beta_t*sample_final$n_evt_arm0)) #23 FOE sECNERIO 2
  #I
  I_k = I


  p_k <- 0.50 # Median survival
  d_k1 <- sum(datax$status) # event in control arm
  d_k2 <- sum(datay$status) # event in treatment arm
  d_k <- sum(data_int$status)
  D <- ceiling(2*sample_final$n_evt_arm0)
  RT_k <- unname(log((wei_scale_t/wei_scale_c)*(log(1/(1-p_k)))^(1/wei_shape_t - 1/wei_shape_c) )) # Interim stage mu
  sd_k <- sqrt(((1/r*wei_shape_t^2) + (1/wei_shape_c^2))/d_k) # Interim stage at k std 3 d_k = 73;r=1
  #sd_k <- sqrt(wei_sd_c^2+wei_sd_t^2)
  #sd_k <- sqrt(wei_sd_c^2/d_k1+wei_sd_t^2/d_k2)   # I am changing here
  #i_k <- 1/sd_k^2 # interim information time
  #i_k <- 1/ (wei_sd_c^2/d_k1 + wei_sd_t^2/d_k2)
  #i_k = 1/(1/d_k1 + 1/d_k2) # 12.98
  i_k = 1/(1/(beta_c*d_k1) + 1/(beta_t*d_k2)) # 12.98
  #I = 1/(1/70 + 1/70)
  t = i_k/I_k # or t = d_k/D # I can vary the information time/ Information fraction # 0.56
  #t = d_k/D

  z_k <- (RT_k- mu_d)/sd_k # Interim stage test statistics z_k = -0.968 FOR SECNERIO 2
  z_relative  = z_k
  theta = log(mu_d)

  #CP_n0 <- 1 - pnorm((- (z_k*sqrt(i_k)) + z_alpha*sqrt(I_k))/(sqrt(I_k-i_k)))
  CP_n0 <- pnorm((- (z_k*sqrt(i_k)) - z_alpha*sqrt(I_k))/(sqrt(I_k-i_k))) # 0.0823

  CP_A <- pnorm((- (z_k*sqrt(i_k)) - z_alpha*sqrt(I_k)-theta*(I_k-i_k))/(sqrt(I_k-i_k)))  # Lower One sided based on alternative hypothesis : 0.409

  #CP_A <- 1 - pnorm((- (z_k*sqrt(i_k)) + z_alpha*sqrt(I_k) - log(mu_d)*(I_k-i_k))/(sqrt(I_k-i_k)))
  CP_A_F = 1- CP_A # Futility under Design

  # Under current trend # Changed it under the Current Trend
  theta_T <- z_k*sqrt(t) # Trend: A review of methods for futility stopping based on conditional power: Lachin # -0.727
  CP_T <- pnorm((- (z_k*sqrt(i_k)) - z_alpha*sqrt(I_k)-theta_T*(I_k-i_k))/(sqrt(I_k-i_k)))  # Upper One sided FOR design/protocol Theta :0.819



  # # Adaptive condiitonal power
  C <- 1 # or 2.326 #
  CP_ADAP <- pnorm((- (z_k*sqrt(i_k)) - z_alpha*sqrt(I_k)-(z_k+C*sd_k)*(I_k-i_k))/(sqrt(I_k-i_k)))
  #1- pnorm(((z_alpha*sqrt(I) - z_k*sqrt(i_k) - ((z_k+ C)*(I- i_k))/sqrt(i_k) ))/(sqrt(I-i_k)))  # Upper One sided FOR design/protocol Theta
  CP_ADAP_F = 1- CP_ADAP # Futility under Design


  # # Predictive power with uniform prior
  #
  # PP_U <- pnorm((-(z_k*sqrt(I_k)) - z_alpha*sqrt(i_k))/sqrt(I_k - i_k)) # 0.343
  # PP_FU <- 1- PP_U

  # Predictive power with uniform prior
  # alpha = 0.10
  # z_alpha <-  qnorm(alpha,lower.tail = F)
  PP_U <- pnorm((-(z_k*sqrt(I_k)) - z_alpha*sqrt(i_k))/sqrt(I_k - i_k))
  PP_FU <- 1- PP_U

  # Posterior predicitive probability under uniform prior
  BPP_U <- pnorm((-(z_k*sqrt(I_k)) - z_alpha*sqrt(i_k)- theta*sqrt(i_k) )/sqrt(I_k - i_k))
  #pnorm(( z_k*sqrt(I) - qnorm(0.95,lower.tail = FALSE)*sqrt(i_k)  - mu_d*sqrt(i_k) )/(sqrt(I-i_k)))
  BPP_FU <- 1- BPP_U


  # Predictive power with Optimistic prior (Dignam 1998 paper)
  sigma <- sqrt(sd_k^2*d_k) # common sigma var(log(HR_mid)) = sigma^2/d_k
  #mu_0 <-  0.601 # theta = 0.58 # prior
  sd_0 <- mu_d/qnorm(alpha,lower.tail = F) # Diagnam (1998) page 8 : p(mu_I > mu_E) = c (0.05) => 0.95 = p(mu_I < mu_E ) or
  # pnomr(0.95) = mu_E-0/sd_0 || Frayer (1997 tutorial)
  d_0 <- ceiling((sigma^2/(sd_0^2))) # Enthusianstic prior # 39 EVENTS
  mu_0 <- mu_d
  part_1 <- sqrt((d_0 * (D-d_k))/((d_0 + d_k)*(d_0+D))) * (sqrt(d_0)*mu_0)/sigma
  part_2 <- sqrt(((D-d_k) * (d_0 + D))/((D - d_k)*(d_0+d_k))) * (sqrt(d_k)*RT_k)/sigma
  part_3 <- sqrt((D*(d_0+d_k))/((D-d_k)*(d_0+D)) )*qnorm(alpha,lower.tail = F)
  PP_E <-  pnorm(-part_3 + part_1  + part_2  ) #
  PP_FE <- 1- PP_E

  # Weakly informative prior
  mu_w <- 0
  part_1 <- sqrt((d_0 * (D-d_k))/((d_0 + d_k)*(d_0+D))) * (sqrt(d_0)*mu_w)/sigma
  part_2 <- sqrt(((D-d_k) * (d_0 + D))/((D - d_k)*(d_0+d_k))) * (sqrt(d_k)*RT_k)/sigma
  part_3 <- sqrt((D*(d_0+d_k))/((D-d_k)*(d_0+D)) )*qnorm(alpha,lower.tail = F)
  PP_W <-  pnorm(-part_3 + part_1  + part_2  ) #
  PP_FW <- 1- PP_W



  # Posterior Predictive probability wiith Optimistic prior
  # eta = 0.95
  delta_I = 0.50
  part_4 <- (((d_0*mu_0 + d_k*RT_k)/(d_0 + d_k))- delta_I) * (sqrt(d_0+d_k)/sigma) # Change the value here
  #part_4 <- (d_0*mu_0 + d_k*log(RT_Q))/(sqrt(d_0+d_k)*sigma) # if sigma is 1.51 and delta_I = 0
  part_5 <- sqrt((d_0+D)/(D-d_k))
  part_6 <- sqrt((d_0 +d_k)/(D-d_k))*qnorm(0.05,lower.tail = FALSE)
  BPP_E <-   pnorm(part_6 + part_4*part_5 )
  PPP_FE <- 1- BPP_E


  # WEAKLY informative

  mu_w <- 0
  part_4 <- (((d_0*mu_w + d_k*RT_k)/(d_0 + d_k))- delta_I ) * (sqrt(d_0+d_k)/sigma) # Change the value here
  #part_4 <- (d_0*mu_0 + d_k*log(RT_Q))/(sqrt(d_0+d_k)*sigma) # if sigma is 1.51 and delta_I = 0
  part_5 <- sqrt((d_0+D)/(D-d_k))
  part_6 <- sqrt((d_0 +d_k)/(D-d_k))*qnorm(0.05,lower.tail = FALSE)
  BPP_W <-  pnorm(part_6 + part_4*part_5 )
  PPP_FW <- 1- BPP_W


  # # Uniform prior
  # d_0 <- 0
  # mu_w <- 0
  # part_4 <- (((d_0*mu_w + d_k*RT_k)/(d_0 + d_k))- delta_I ) * (sqrt(d_0+d_k)/sigma) # Change the value here
  # #part_4 <- (d_0*mu_0 + d_k*log(RT_Q))/(sqrt(d_0+d_k)*sigma) # if sigma is 1.51 and delta_I = 0
  # part_5 <- sqrt((d_0+D)/(D-d_k))
  # part_6 <- sqrt((d_0 +d_k)/(D-d_k))*qnorm(0.05,lower.tail = FALSE)
  # BPP_U <-  pnorm(part_6 + part_4*part_5 )
  # PPP_Fu <- 1- BPP_U

  # Simulation based approach
  # Try to make similar simulation code for the predicitve prbability paper for enthusiastic prior

  #d_0 = 19
  #mu_0 = 0.58
  #d_k = 84
  #RT_k = 0.37
  #sigma = 1.46

  #D = 140

  prob1 <- c(NA)
  prob2 <- c(NA)
  #data <- c(NA)
  set.seed(111)
  for (i in 1: 10000){
    prob1[i] <- 1- pnorm(RT_k, mean = ((d_0*mu_0 + d_k*(RT_k))/(d_0 + d_k)),
                         sd = sqrt(sigma^2/(d_0+d_k)), lower.tail = TRUE)
    mean(prob1)


    k = ( delta_I *(d_0 +D) - qnorm(0.05,lower.tail = FALSE) * sqrt(sigma^2*(d_0+D)) -  (d_0*mu_0 + d_k*RT_k) )/(D-d_k)
    mean = ((d_0*mu_0 + d_k*RT_k)/(d_0 + d_k))
    sd = sqrt(sigma^2/(d_0+d_k)+ (sigma^2/(D-d_k)))
    data<- suppressWarnings(rnorm(10000,mean = mean,
                                  sd = sd))

    prob2[i] <- 1- pnorm(k, mean = mean(data),
                         sd = sd(data), lower.tail = TRUE)

  }



  BPP_E_sim <- mean(prob2,na.rm = T)

  # Simulation based approach
  # Try to make similar simulation code for the predicitve prbability paper for enthusiastic prior

  #d_0 = 19
  #mu_0 = 0.58
  #d_k = 84
  #RT_k = 0.37
  #sigma = 1.46

  #D = 140

  prob1_w <- c(NA)
  prob2_w <- c(NA)
  #data_w<- c(NA)
  set.seed(111)
  for (i in 1: 10000){
    prob1_w[i] <- 1- pnorm(RT_k, mean = ((d_0*mu_w + d_k*(RT_k))/(d_0 + d_k)),
                           sd = sqrt(sigma^2/(d_0+d_k)), lower.tail = TRUE)
    mean(prob1_w)


    k = ( delta_I *(d_0 +D) - qnorm(0.05,lower.tail = FALSE) * sqrt(sigma^2*(d_0+D)) -  (d_0*mu_w + d_k*RT_k) )/(D-d_k)

    mean = ((d_0*mu_w + d_k*RT_k)/(d_0 + d_k))
    sd = sqrt(sigma^2/(d_0+d_k)+ (sigma^2/(D-d_k)))

    data_w <- suppressWarnings(rnorm(10000,mean = mean,sd = sd))

    prob2_w[i] <- 1- pnorm(k, mean = mean(data_w),
                           sd = sd(data_w), lower.tail = TRUE)

  }



  BPP_W_sim <- mean(prob2_w,na.rm = T)

  #Uniform prior: Working but not going to use this
  # d_0 = 0.00001
  # prob1_u <- c(NA)
  # prob2_u <- c(NA)
  # #data_w<- c(NA)
  # set.seed(123)
  # for (i in 1: 10000){
  #
  #   prob1_u[i] <- 1- pnorm(RT_k, mean = ((d_0*mu_w + d_k*(RT_k))/(d_0 + d_k)),
  #                          sd = sqrt(sigma^2/(d_0+d_k)), lower.tail = TRUE)
  #   mean(prob1_u)
  #
  #
  #   k = ( delta_I *(d_0 +D) - qnorm(0.05,lower.tail = FALSE) * sqrt(sigma^2*(d_0+D)) -  (d_0*mu_w + d_k*RT_k) )/(D-d_k)
  #
  #   mean = ((d_0*mu_w + d_k*RT_k)/(d_0 + d_k))
  #   sd = sqrt(sigma^2/(d_0+d_k)+ (sigma^2/(D-d_k)))
  #
  #   data_u <- suppressWarnings(rnorm(10000,mean = mean,sd = sd))
  #
  #   prob2_u[i] <- 1- pnorm(k, mean = mean(data_u),
  #                          sd = sd(data_u), lower.tail = TRUE)
  #
  # }

  #
  #
  # BPP_u_sim <- mean(prob2_u,na.rm = T)
  #
  #
  return(data.frame(beta_t,d_k1,d_k2,D,i_k,I_k,t,z_relative,RT_k,CP_n0,CP_A,CP_T,CP_ADAP,PP_U,PP_E,PP_W,BPP_U,BPP_E,BPP_E_sim,BPP_W,BPP_W_sim,sigma,sd_k,sd_0,SD,d_0,mu_d))

}
