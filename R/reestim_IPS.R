#' @title Calculate the event using the IPS approach

#' @description Here, we used the IPS approach to re-estimate the sample size based on the relative time method
#'
#' @param beta0 numeric variable
#' @param Interim.event numeric variable
#' @param fixed_sample numeric variable
#' @param scenario numeric variable between 1 and 4
#'
#' @return
#' @export
#'
#' @examples
#' event(beta0 = 0.462,Interim.event = 53,fixed_sample = 92,scenario = 1)

event <- function(beta0,Interim.event,fixed_sample,scenario){
  if (Interim.event <= 0 ) stop("Error:The value Interim.event is not allowed. Interim event should be greater than 0")
  if (fixed_sample <= 0 ) stop("Error:The value fixed_sample is not allowed. fixed_sample should be greater than 0")

  if (beta0 <= 0 | beta0 >= 10)   stop("Error:The value beta0 is not allowed, The Beta0 should be between 0 and 10")
  if (scenario <= 0 | scenario >= 5)   stop("Error:The value scenario is not allowed, The scenario should be between 1 and 4")

  beta = beta0
  # Sec 1
  if (scenario == 1){
    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta,median_c = 4,p1 = 0.10,p2 = 0.90,RT_P1 = 1.52,
                              RT_P2 = 1.98,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  } else if (scenario == 2){
    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.10,p2 = 0.90,RT_P1 = 2.00,
                              RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  } else if (scenario == 3){
    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.50,
                              RT_P2 = 1.667,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  } else if (scenario == 4){
    sample_int <- sample_size(alpha = 0.05,sides = 1,beta0  = beta0,median_c = 4,p1 = 0.25,p2 = 0.75,RT_P1 = 1.667,
                              RT_P2 = 1.50,qmin = 0.001,qmax = 0.999,Power = 0.8,r = 1,a = 12,f = 12,p_event=0.8)
  }
  beta1 <- sample_int$beta1
  event <- ceiling(2*(sample_int$n_evt_arm0))
  event_req <- event - Interim.event
  event_req <- ifelse(event_req<0,"No new event required",event_req)
  sample_size <- ceiling(2*(sample_int$n_admin_arm0)) - fixed_sample
  sample_size <- ifelse(sample_size<0,"No new sample size required",sample_size)
  res <- c(beta,beta1,event_req,event,sample_size)
  names(res) <- c("beta","beta1","event_req","event","sample_size")
  return(res)

}





