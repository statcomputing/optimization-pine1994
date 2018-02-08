observed_sample <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
                     3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
start_value <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38, mean(observed_sample))
alpha <- c(1, 0.64, 0.25)
#a represents thita
loglikelihood <- function(a){
  n <- length(observed_sample)
  l_a <- -1*(n*log(pi) + sum(log(1 + (a-observed_sample)^2)))
  return(l_a)
}

#graph the log-likelihood
a_value <- seq(-10, 60, by = 0.05)
log_likelihood <- sapply(a_value, loglikelihood)
plot(a_value, log_likelihood, xlab = "theta", ylab = "log_likelihood", type = "l")

#use Newton-Raphson method to find the MLE for theta
#parameter a of the function is the start value of theta
deriv_llf <- function(theta){
  deriv <- 0
  for(i in 1:length(observed_sample))
    deriv <- deriv + (-2) * ((theta - observed_sample[i])/(1 + (theta - observed_sample[i])^2))
  return(deriv)
}

hessian_llf <- function(theta){
  hess <- 0
  for(i in 1:length(observed_sample))
    hess <- hess + (-2) * (sum((1 - (theta - observed_sample[i])^2)/(1 + (theta - observed_sample[i])^2)^2))
  return(hess)
}
##question 1_b---------------------------------
nr_method <- function(a){
  theta <- array()
  theta[1] <- a
  #initialize the h_t equal to 1 in order to begin the loop 
  h_t <- 1
  i <- 1
  #we set the max loop number be 100 and the upper bound of error be 0.0001
  while(abs(h_t) > 0.0001 & i<100){
    h_t <- (-1)*deriv_llf(theta[i])/hessian_llf(theta[i])
    theta[i+1] <- theta[i] + h_t
    i <- i+1
  }
  return(theta[i])
}
# a <- nr_method(7)
# a
nr_method_value <- sapply(start_value, nr_method)
nr_method_value

##question 1_c---------------------------------
fixpoint_method <- function(a,alpha){
  theta <- array()
  theta[1] <- a
  h_t <- 1
  i <- 1
  while(abs(h_t) > 0.0001 & i<1000){
    h_t <- alpha*deriv_llf(theta[i])
    theta[i+1] <- theta[i] + h_t
    i <- i+1
  }
  return(theta[i])
}

a1 <- fixpoint_method_value <- sapply(start_value,fixpoint_method, alpha = alpha[1])
a2 <- fixpoint_method_value <- sapply(start_value,fixpoint_method, alpha = alpha[2])
a3 <- fixpoint_method_value <- sapply(start_value,fixpoint_method, alpha = alpha[3])

##question 1_d-------------------------------
fisherscore_method <- function(a){
  theta <- array()
  theta[1] <- a
  #initialize the h_t equal to 1 in order to begin the loop 
  h_t <- 1
  i <- 1
  #we set the max loop number be 100 and the upper bound of error be 0.0001
  while(abs(h_t) > 0.01 & i<1000){
    h_t <- deriv_llf(theta[i])/(length(observed_sample) / 2)
    theta[i+1] <- theta[i] + h_t
    i <- i+1
  }
  return(theta[i])
}

fisherscore_method_value <- sapply(start_value, fisherscore_method)
#refine the value by newton raphson method
refine_value <- sapply(fisherscore_method_value,nr_method)
refine_value
