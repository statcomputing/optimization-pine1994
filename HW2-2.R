observed_sample <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
       2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)
n <- length(observed_sample)
##question 2-a----------------------------
loglikelihood <- function(a){
  l_a <- 0
  for(i in 1:n)
    l_a <- l_a + log(1/(2*pi)-cos(observed_sample[i]-a)/(2*pi))
  return(l_a)
}


a_value <- seq(-pi, pi, by = 0.02)
log_likelihood <- sapply(a_value, loglikelihood)
plot(a_value, log_likelihood, xlab = "theta", ylab = "log_likelihood", type = "l")

sample_mean <- mean(observed_sample)

cond_expected <- function(a){
  expected_value <- pi+sin(a)
  return(expected_value)
} 
# a_value <- seq(-pi, pi, by = 0.02)
# as <- sapply(a_value, cond_expected)
# plot(a_value, as, xlab = "theta", ylab = "log_likelihood", type = "l")

##question 2-b---------------------------------
moment_theta <- function(a){
  pi+sin(a)-sample_mean
}

r1 <- uniroot(moment_theta,c(-pi, pi/2))$root
r2 <- uniroot(moment_theta,c(pi/2, pi))$root

##question 2-c--------------------------------
start_value <- c(r1,r2)

deriv_llf <- function(a){
  deriv <- sum(sin(observed_sample-a)/(1-cos(observed_sample-a)))
  return(deriv)
}

hessian_llf <- function(a){
  hess <-(sum(1/(1-cos(observed_sample-a))))
  return(hess)
}

nr_method <- function(a){
  theta <- array()
  theta[1] <- a
  #initialize the h_t equal to 1 in order to begin the loop 
  h_t <- 1
  i <- 1
  #we set the max loop number be 100 and the upper bound of error be 0.0001
  while(abs(h_t) > 0.0001 & i<1000){
    h_t <- (-1)*deriv_llf(theta[i])/hessian_llf(theta[i])
    theta[i+1] <- theta[i] + h_t
    i <- i+1
  }
  return(theta[i])
}

nr_method_value <- sapply(start_value, nr_method)
nr_method_value

##question 2-d--------------------------------

nr_method(2.7)
nr_method(-2.7)

##question 2-e--------------------------------

start_value <- seq(-pi,pi,length.out = 200)
nr_method_value <- sapply(start_value, nr_method)
nr_value <- round(nr_method_value,6)
freq <- as.data.frame(table(nr_value))
interval_index <- array()
interval_index[1] <- freq$Freq[1]
for(i in 2:length(freq$Freq)){
  interval_index[i] <- freq$Freq[i]+interval_index[i-1]
}
freq[1]
internal <- c()
for(i in 1:length(interval_index)){
  internal[i] <- start_value[interval_index[i]]
}
#internal
table_converge <- cbind(freq[1],internal)
table_converge
# nr_value_all <- unique(nr_value)
# rank_nrvalue <- array()
# rank_nrvalue[1] <- -pi
















