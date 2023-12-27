pop <- function(N, contam){
  
  # Intializing values of covariates
  set.seed(1)
  x1 <- rnorm(N, mean = 0, sd = 1)
  set.seed(2)
  x2 <- rnorm(N, mean = 0, sd = 1)
  
  # Initializing error terms
  # Combine two error term distributions reshuffle using SRSWOR
  set.seed(3)
  error <- sample(c(rnorm((1-contam)*N, mean = 0, sd = 5), 
                 rnorm((contam)*N, mean = 0, sd = 25)))
  
  # Initializing response values (with beta_0 = 1, beta_1 = 1, and beta_2 = 1)
  y <- 1 + x1 + x2 + error
  popn <- list(y = y, x1 = x1, x2 = x2, error = error)
  return(popn)
}

n <- 100
contam <- 0.5

popn <- pop(n, contam)

one <- as.matrix(rep(1, n))
y <- as.matrix(popn$y)
x <- as.matrix(cbind(rep(1, n), popn$x1, popn$x2))
beta <- as.matrix(c(0, 0, 0))


f_prime <- function(beta, n){
  
  a <- as.matrix(c(0, 0, 0))
  
  for(i in 1:n){
    z <- tanh(x[i,1]*beta[1,] + x[i,2]*beta[2,] + x[i,3]*beta[3,] - y[i,])*as.matrix(x[i, ])
    a <- a + z
  }
  return(a)
}


f_prime_2 <- function(beta, n){
  
  # Initialize 
  a <- as.matrix(cbind(rep(0, 3), rep(0, 3), rep(0, 3)))
  
  for(i in 1:n){
    cross <- as.matrix(x[i, ])%*%t(as.matrix(x[i, ]))
    z <- (1/cosh(x[i,1]*beta[1,] + x[i,2]*beta[2,] + x[i,3]*beta[3,] - y[i,]))^2*cross
    a <- a + z
  }
  return(a)
}

newton_fx <- function(beta, iter){
  prev <- beta
  k <- 1
  for(k in 1:iter){
    guess <- prev - solve(f_prime_2(prev, n))%*%f_prime(prev, n)
    prev <- guess
  }
  return(guess)
}


