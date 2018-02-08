rm(list = ls())
beetles <- data.frame( days = c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154),
                       beetles = c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024))
N0 <- beetles$beetles[1]
n <- length(beetles$days)
A_calculate <- function(k,r){
  A <- matrix(0,nrow = n, ncol = 2)
  for(i in 1:n){
    A[i,1] <- (N0^2-N0^2*exp(-1*r*beetles$days[i]))/((N0+(k-N0)*exp(-1*r*beetles$days[i]))^2)
    A[i,2] <- beetles$days[i]*(k-N0)*exp(-1*r*beetles$days[i])*k*N0/((N0+(k-N0)*exp(-1*r*beetles$days[i]))^2)
  }
  return(A)
}

Z_calculate <- function(k,r){
  Z <- matrix(0,nrow = n, ncol = 1)
  for(i in 1:n){ 
    Z[i,1] <- beetles$beetles[i] - k*N0/(N0+(k-N0)*exp(-1*r*beetles$days[i]))
  }
  return(Z)
}

h_t_calculate <- function(k,r){
 A <- A_calculate(k,r)
 Z <- Z_calculate(k,r)
 h_t <- solve(t(A) %*% A) %*% t(A) %*% Z
 return(h_t)
}

##question 3-a-----------------------------
gauss_newton_method <- function(k0,r0){
  h_t <- c(1000,0.1)
  error <- sum(abs(h_t))
  error_modified_a <- 2
  error_modified_b <- 1
  i <- 1
  while(abs(error_modified_b - error_modified_a) > 0.00000001 & i<1000){
    error_modified_a <- error_modified_b
    h_t <- h_t_calculate(k0,r0)
    k0 <- k0+h_t[1]
    r0 <- r0+h_t[2]
    error <- sum(abs(h_t))
    error_modified_b <- sum(abs(h_t))
    i <- i+1
    #print(h_t)
  }
  theta <- c(k0,r0,i)  
  return(theta)
}

a <- gauss_newton_method(1200,0.5)

##question 3-b--------------------------------
library(rgl)
squared_error_calculate <- function(k,r){
  squared_error <- 0
  for(i in 1:n){
    squared_error <- squared_error + (k*N0/(N0 + (k-N0)*exp(-1*r*beetles$days[i]))-beetles$beetles[i]) ^ 2
  }
  return(squared_error)
}
r_value <- seq(0.01,0.2,length.out = 100) 
r_value_complete <- rep(r_value, times = rep(100,100))

k_value <- seq(800,1300,length.out = 100)
k_value_complete <- rep(k_value,100)

N <- c()
for (i in 1:(length(r_value)*length(k_value))){
    N[i] <- squared_error_calculate(k_value_complete[i], r_value_complete[i])
}
#View(N)
plot3d(x = k_value_complete,y = r_value_complete, z = N, col="pink", size=5)
library(MASS)

k_value <- seq(100,1500,length.out = 100)
r_value <- seq(0.01,0.5,length.out = 100)

N <- matrix(0, nrow = 100, ncol = 100)
for(i in 1:length(k_value)){
  for(j in 1:length(r_value)){  
    N[i,j] <- squared_error_calculate(k_value[i],r_value[j])
  }
}
contour(k_value,r_value, N, col = 'red', drawlabel=FALSE, main="Density estimation :cont Plot")


##question 3-c-----------------------------


l_k_derivative <- function(k,r,sigma){
  l_k_derivative_value <- 0
  for(i in 1:n){
    temp <- (N0^2-N0^2*exp(-1*r*beetles$days[i]))/((N0+(k-N0)*exp(-1*r*beetles$days[i]))^2)
    l_k_derivative_value <- l_k_derivative_value + (1/(sigma*sigma))*(log(beetles$beetles[i])-log(k*N0/(N0+(k-N0)*exp(-1*r*beetles$days[i]))))/(k*N0/(N0+(k-N0)*exp(-1*r*beetles$days[i])))*temp
  }
  return(l_k_derivative_value)
}
#l_r_derivative(1000,0.1,1)
l_r_derivative <- function(k,r,sigma){
  l_r_derivative_value <- 0
  for(i in 1:n){
    temp <- beetles$days[i]*(k-N0)*exp(-1*r*beetles$days[i])*k*N0/((N0+(k-N0)*exp(-1*r*beetles$days[i]))^2)
    l_r_derivative_value <- l_r_derivative_value + (1/(sigma*sigma))*(log(beetles$beetles[i])-log(k*N0/(N0+(k-N0)*exp(-1*r*beetles$days[i]))))/(k*N0/(N0+(k-N0)*exp(-1*r*beetles$days[i])))*temp
  }
  return(l_r_derivative_value)
}

l_sigma_derivative <- function(k,r,sigma){
  l_sigma_derivative_value <- 0
  for(i in 1:n){
    temp <- k*N0/(N0+(k-N0)*exp(-1*r*beetles$days[i]))
    l_sigma_derivative_value <- l_sigma_derivative_value +((log(beetles$beetles[i])-log(temp))^2)/(sigma^3)
  }
  l_sigma_derivative_value <- l_sigma_derivative_value - n/(2*sigma)
  return(l_sigma_derivative_value)
}
#l_sigma_derivative(1000,0.1,1)

l_derivative <- function(k,r,sigma){
  l_derivative_value <- c()
  l_derivative_value[1] <- l_k_derivative(k,r,sigma)
  l_derivative_value[2] <- l_r_derivative(k,r,sigma) 
  l_derivative_value[3] <- l_sigma_derivative(k,r,sigma) 
  return(l_derivative_value)
}
#a <- l_derivative(1000,0.1,1)

fixpoint_method <- function(k0,r0,sigma0){
  h_t <- c(1000,1,1)
  error <- sum(abs(h_t))
  i <- 1
  error_list <- c()
  while(error>0.001 & i<3000){
    h_t <- l_derivative(k0,r0,sigma0)
    #result <- c(k0,r0,sigma0,i)
    #print(result)
    k0 <- k0 + h_t[1] * 1000#rnorm(1, 5000000, 1000000)
    r0 <- r0+h_t[2] * 0.0003#rnorm(1, 0.0005, 1)
    sigma0 <- sigma0+h_t[3]* 0.0003#rnorm(1, 0.05, 0.01)
    error <- sum(abs(h_t))
    error_list[i] <-error 
    i <- i+1

    #print(error)
  }
  result <- c(k0,r0,sigma0,i,error_list)
  #print(error)
  return(result)
}
a <- fixpoint_method(800,0.15,0.1)
a_error <- c()
for(i in 6:length(a)){
  a_error[i-5] <- a[i]
}
var(a_error)
b <- fixpoint_method(1200,0.15,0.1)
b
# 820.3800335
# 0.1926401
# 0.9109554
# 672.0907892 
# 0.4004609   
# 0.8286135





