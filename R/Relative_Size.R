
#' @title Calculate Sample size using the Relative Time method under Non-Proportional Hazard and Non-Proportional Time scenarios
#' @param alpha numeric variable
#' @param sides binary variable
#' @param beta0 numeric variable
#' @param median_c numeric variable
#' @param p1 numeric variable
#' @param p2 numeric variable
#' @param p_user numeric variable
#' @param RT_P1 numeric variable
#' @param RT_P2 numeric variable
#' @param qmin numeric variable
#' @param qmax numeric variable
#' @param Power numeric variable
#' @param r numeric variable
#' @param a numeric variable
#' @param f numeric variable
#' @param p_event numeric variable
#'
#' @return Sample size Calculation using the Relative Time method proposed by Phadnis and Mayo (2021).
#' @import ggplot2 survival
#' @export
#'
#' @examples
#' Relative_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
#' \dontrun{
#' Relative_size(alpha = 0.05,sides = 1,beta0  = 1,median_c = 4,p1 = 0.1,p2 = 0.9,RT_P1 = 1.52,RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
#' }
#'
Relative_size <- function(alpha,sides,beta0,median_c,p1,p2,p_user=NULL,RT_P1,RT_P2,qmin,qmax,Power,r,a,f,p_event){

  if (alpha <= 0 | alpha >= 1)    stop("Error:The value alpha is not allowed, The type I error should be between 0 and 1")
  if (sides != 1 & sides != 2)    stop("Error:The value sides is not allowed. Sides should be either 1 or 2")
  if (beta0 <= 0 | beta0 >= 10)   stop("Error:The value beta0 is not allowed, The Weibull shape paramter should be between 0 and 10")
  if (median_c <= 0 )             stop("Error:The value for MEDIAN TIME (Months) OF CONTROL ARM is not allowed. The MEDIAN should be > 0")
  if (p1 <= 0 | p1 >= 1)          stop("Error:The value p1 is not allowed, p1 should be between 0 and 1")
  if (p2 <= 0 | p2 >= 1)          stop("Error:The value p2 is not allowed, p2 should be between 0 and 1")
  if (RT_P1 <= 1 | RT_P1 >= 100)  stop("Error:The value RT_P1 is not allowed, RT_P1 should be between 1 and 100")
  if (RT_P2 <= 1 | RT_P2 >= 100)  stop("Error:The value RT_P2 is not allowed, RT_P2 should be between 1 and 100")
  if (qmin <= 0 | qmin >= p1)     stop("Error:The value qmin is not allowed, qmin should be between 1 and p1")
  if (qmax >= 1 | qmax <= p2)     stop("Error:The value qmax is not allowed, qmax should be between p2 and 1")
  if (r <= 0 | r >= 10)           stop("Error:The value r for ALLOCATION RATIO r is not allowed. The ALLOCATION RATIO r should between 0 and 10")
  if (a <= 0 )                    stop("Error:The value a for ACCRUAL TIME is not allowed. The ACCRUAL TIME should be greater than 0")
  if (f <= 0 )                    stop("Error:The value f for FOLLOW-UP TIME is not allowed. The FOLLOW-UP TIME should be greater than 0")
  if (p_event <= 0 | p_event >= 1) stop("Error:The value p_event for EVENT RATE is not allowed. The EVENT RATE should between 0 and 1")

  # Actual Code Start Here

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
  median_t = theta1*((log(2))^sigma1)
  #mid_p = (p1+p2)/2
  if(!is.null(p_user)){
    mid_p = p_user
  } else{
    mid_p = (p1+p2)/2
  }
  # set.seed(111)
  t_mid = qgamma(mid_p,1)
  g_mid = log(t_mid)
  RT_mid = exp((mu1 - mu0) + g_mid*(sigma1 - sigma0))
  RT_midd = round(RT_mid,0.0001)

  g_min = log(qgamma(qmin,1))
  g_max = log(qgamma(qmax,1))
  RT_min = exp((mu1 - mu0) + g_min*(sigma1 - sigma0))
  RT_max = exp((mu1 - mu0) + g_max*(sigma1 - sigma0))


  arm0_sd = 1/beta0
  arm1_sd = (1/beta1)*sqrt(1/r)
  pool_sd = sqrt(arm0_sd^2 + arm1_sd^2 )

  p_alpha = 1-alpha/sides
  z_crit = qnorm(p_alpha,lower.tail = T)
  z_pow = qnorm(Power,lower.tail = T)


  log_HR_t0 = log((beta1/beta0)*(theta0^(beta0))/(theta1^(beta1)))
  t_avg = (median_c+median_t)/2
  HR_t_avg_RT = exp(log_HR_t0)*((t_avg)^(beta1 - beta0))

  n_evt_arm0 = ( (((z_crit + z_pow)*pool_sd)/log(RT_mid))^2 ) # Chnage made here
  n_evt_arm0
  n_evt_arm1 = n_evt_arm0*r # Change made here
  n_evt_arm1

  s_f_0 = exp(-((f/theta0)^beta0))
  s_a_half_0  = exp(-(((f+a*0.5)/theta0)^beta0))
  s_a_f_0  = exp( - (((a+f)/theta0)^beta0))
  d_0 = 1 - ((s_f_0 + 4*s_a_half_0 +s_a_f_0)/6) # EVENT RATE FOR THE CONTROL ARM
  I1_0 = ((a+f)/a) * (s_f_0 - s_a_f_0)
  h_0 = (1/beta0)+1
  g_0 = gamma(h_0)
  m1_0 = ((a+f)/theta0)^beta0
  m2_0 = (f/theta0)^beta0
  IG1_0 = pgamma(m1_0,h_0)
  IG2_0 = pgamma(m2_0,h_0)
  I2_0 = (theta0/a)*g_0*(IG1_0 - IG2_0)
  d_true_0 = I1_0 - I2_0 + 1- s_f_0

  s_f_1 = exp(-((f/theta1)^beta1))
  s_a_half_1  = exp(-(((f+a*0.5)/theta1)^beta1))
  s_a_f_1  = exp( - (((a+f)/theta1)^beta1))
  d_1 = 1 - ((s_f_1 + 4*s_a_half_1 +s_a_f_1)/6)
  I1_1 = ((a+f)/a) * (s_f_1 - s_a_f_1)
  h_1 = (1/beta1)+1
  g_1 = gamma(h_1)
  m1_1 = ((a+f)/theta1)^beta1
  m2_1 = (f/theta1)^beta1
  IG1_1 = pgamma(m1_1,h_1)
  IG2_1 = pgamma(m2_1,h_1)
  I2_1 = (theta1/a)*g_1*(IG1_1 - IG2_1)
  d_true_1 = I1_1 - I2_1 + 1- s_f_1 # EVENT RATE FOR THE TREATMENT ARM

  d_approx = (d_0 +d_1)/2
  d_exact = (d_true_0 +d_true_1)/2
  n_admin_arm0 = ( n_evt_arm0/d_exact )
  #n_admin_arm0 # SAMPLE SIZE UNADJUSTED FOR 1-RHO (DROP OUT DUE TO LOSS)
  n_admin_arm1 = n_admin_arm0*r

  N_arm0 = ceiling(n_admin_arm0/(p_event)) # 1- drop out rate
  N_arm1 = N_arm0*r


  # The internal macro from SAS code

  x = matrix(data=NA, nrow=3, ncol=2)
  y = matrix(data=NA, nrow=3, ncol=1)
  for (i in 1:3) {
    x[i,1] = 1
  }

  x[1,2] = g_p1 # low_g
  x[2,2] = g_p2 #high_g
  y[1,1] = log(RT_P1)
  y[2,1] = log(RT_P2)

  # Check_min = RT_min
  if ( (RT_min <= 1) & (RT_P1 < RT_P2)){
    x[3,2] = log(qgamma(qmin,1))
    y[3,1] = log(1.01)
  }

  if ( (RT_max <= 1) & (RT_P1 > RT_P2)){
    x[3,2] = log(qgamma(qmax,1))
    y[3,1] = log(1.01)
  }


  if ((RT_min > 1) & (RT_P1 < RT_P2)){
    x[3,2] = log(qgamma(qmin,1))
    y[3,1] = log(RT_min)
  }

  if ((RT_max > 1) & (RT_P1 > RT_P2)){
    x[3,2] = log(qgamma(qmax,1))
    y[3,1] = log(RT_max)
  }

  param = solve(t(x)%*%x)%*%t(x)%*%y
  sigma0 = 1/beta0
  theta0 = median_c/(log(2)^sigma0)
  mu0 = log(theta0)
  sigma_new = param[2,1] +sigma0
  mu_new = param[1,1] +mu0

  RT_new_p1 <- NULL
  RT_new_p2 <- NULL

  if ((RT_min <= 1) & (RT_P1 < RT_P2)){
    RT_new_p1 = ceiling(100*exp(mu_new - mu0 + g_p1*(sigma_new - sigma0)))/100
    RT_new_p2 = floor(100*exp(mu_new - mu0 + g_p2*(sigma_new - sigma0)))/100
  }


  if ((RT_max <= 1) & (RT_P1 > RT_P2)){
    RT_new_p1 = floor(100*exp(mu_new - mu0 + g_p1*(sigma_new - sigma0)))/100
    RT_new_p2 = ceiling(100*exp(mu_new - mu0 + g_p2*(sigma_new - sigma0)))/100
  }

  if ((RT_min <= 1) & (RT_P1 < RT_P2)){
    cat("=====================================================================================\n")
    cat("Sample Size Calculation Using the Concept of Relative Time \n")
    cat("=====================================================================================\n")

    cat("This combination of input for p1, p2, RT_p1 and RT_p2 results in survival curves that cross and hence not allowed.","\n")
    cat("Try again with one of the following options:","\n")
    cat("[1] Choose smaller value for p1 for your choice of RT_p1","\n")
    cat("[2] Choose larger value of RT_p1 for your choice of p1","\n")
    cat("[3] Choose larger value for p2 for your choice of RT_p2","\n")
    cat("[4] Choose smaller value for RT_p2 for your choice of p2","\n")
    cat("[5] Choose a larger value for q_min","\n")
    cat("Recommendation: If not wanting to change values of p1,p2 and q_min,
      consider using RT_p1 =", RT_new_p1, "and RT_p2 =", RT_new_p2, "as your input design values for the effect size.","\n")
  } else if((RT_max <= 1) & (RT_P1 > RT_P2)) {
    cat("=====================================================================================\n")
    cat("Sample Size Calculation Using the Concept of Relative Time \n")
    cat("=====================================================================================\n")

    cat("This combination of input for p1, p2, RT_p1 and RT_p2 results in survival curves that cross and hence not allowed.","\n")
    cat("Try again with one of the following options:","\n")
    cat("[1] Choose larger value for p1 for your choice of RT_p1","\n")
    cat("[2] Choose smaller value of RT_p1 for your choice of p1","\n")
    cat("[3] Choose larger value for p2 for your choice of RT_p2","\n")
    cat("[4] Choose larger value for RT_p2 for your choice of p2","\n")
    cat("[5] Choose a larger value for q_max","\n")
    cat("Recommendation: If not wanting to change values of p1,p2 and q_max,
      consider using RT_p1 =", RT_new_p1, "and RT_p2 =", RT_new_p2, "as your input design values for the effect size.","\n")
  } else {


    #/* The code for displaying the text messages when the two survival curves cross ENDS here */
    cat("=====================================================================================\n")
    cat("Sample Size Calculation Using the Concept of Relative Time method \n")
    cat("=====================================================================================\n")
    writeLines("\n")
    cat("The effect size used for the sample size calculations is based on a TIME RATIO of", RT_mid, "calculated at p = ",mid_p,"\n" )
    cat("A sample size of", N_arm0, "in the CONTROL ARM", "and", N_arm0, "in the TREATMENT ARM is needed to detect this
      effect size with", Power, "POWER using a", sides, "SIDED test with TYPE I ERROR set at", alpha,"\n" )
    cat("The EVENT RATE for this study is", p_event, ",ACCRUAL TIME is", a, "months , FOLLOW-UP TIME is", f, "months, and ALLOCATION RATIO is", r,"\n")

    writeLines("\n")
    cat("=====================================================================================\n")
    cat("Sample Size Calculation Statistics:  \n")
    cat("=====================================================================================\n")
    return(data.frame(sigma0,beta0,theta0,sigma1,beta1,theta1,RT_min, RT_mid ,RT_max,pool_sd,d_true_0,d_true_1,d_approx,
                      d_exact,n_evt_arm0,n_admin_arm0,n_admin_arm1,
                      N_arm0,N_arm1))
  }
}

