library(qgam)
library(mgcv)
library(qgam)
library(glmnet)
library(readr)

#-------------------------------------------------------------------------------

#Section 2.2

x <- seq(1, 10, length.out=100)
y1 <- x + rnorm(100)
y2 <- sin(x) + rnorm(100)

#linear effect of simple linear regression on Y1
linear <- glm(y1 ~ x - 1)
linear_mse <- sum(linear$residuals^2) / length(x)
linear_mse

#Figure 2.1
plot(linear)
plot(x, y1, xlab=expression(X), ylab=expression(Y1))
lines(x, linear$coefficients * x, lwd=2, col="blue")

#linear effect of simple linear regression on Y2
nonlin <- glm(y2 ~ x - 1)
nonlin_mse <- sum(nonlin$residuals^2) / length(x)
nonlin_mse

#Figure 2.1
plot(nonlin)
plot(x, y2, xlab=expression(X), ylab=expression(Y1))
lines(x, x * nonlin$coefficients, lwd=2, col="blue")

#-------------------------------------------------------------------------------

#Subsection 2.3.1

X <- cbind(x, x^2, x^3)
cubic <- glm(y2 ~ X)
cubic_mse <- sum(cubic$residuals^2) / length(x)
cubic_mse

plot(cubic)

#Figure 2.2
plot(x, y2, xlab = "X", ylab = "Y2")
lines(x, t(cubic$coefficients) %*% t(cbind(1, X)), lwd = 2, col="blue")

#Figure 2.3
plot(x, rep(1,length(x)) * cubic$coefficients[1], type = "l", col="blue", xlab = "x", ylab = expression(beta[0]))
plot(x, x * cubic$coefficients[2], type = "l", col="blue", xlab = "x", ylab = expression(beta[1] * x))
plot(x, x^2 * cubic$coefficients[3], type = "l", col="blue", xlab = "x", ylab = expression(beta[2] * x^2))
plot(x, x^3 * cubic$coefficients[4], type = "l", col="blue", xlab = "x", ylab = expression(beta[3] * x^3))

#-------------------------------------------------------------------------------

#Subsection 2.3.2

#store x and y2 in dataframe format
df <- data.frame( X = x , Y2 = y2 )
#fit cubic spline
fitSpline <- gam(Y2 ~ s(X , bs ="cr"), data = df)

#Figure 2.4
plot( fitSpline, col="blue", xlab="X", ylab="Y2", ylim=c(min(y2), max(y2)))
points (x , y2 )

#mean squared error of model
fitSpline_mse <- sum(fitSpline$residuals^2) / length(x)
fitSpline_mse

#-------------------------------------------------------------------------------

#Subsection 2.3.4

#load dataset from loon package
library(loon.data)
data("diabetes")

#define matrix to store 10 sets of 10 predictor coefficients
beta <- matrix(0, 10, 10)
for (i in 1:10){
  #initialise the hyperparameter
  gamma <- 0.01 * i
  #train cubic regression spline model using hyperparameter
  fit <- gam(RelativeWeight ~ s(FastingPlasmaGlucose, k = 10, bs = "cr"),
              data = diabetes, scale = gamma )
  #store the predictor coefficients
  beta [i,] <- fit$coefficients
}

#Figure 2.5
plot(seq(1, 10), rep(0, 10), type="b", xlab=expression(gamma), ylab=expression(beta), ylim=c(-0.10, 0.10))
for (i in 2:10){
  lines(seq(1, 10), beta[,i], type="b", col=i)
}

#-------------------------------------------------------------------------------

#Section 3.2

df <- data.frame( X = x , Y2 = y2 )
#fit cubic regression spline onto data
fitSpline <- gam(Y2 ~ s(X , bs ="cr"), data = df)
#output the fitted values
y2pred <- predict(fitSpline)

#compute the 0.95-th conditional quantile for each fitted values
cond_quan <- rep(0, length(y2pred))
for (i in 1:length(y2pred)){
  cond_quan[i] <- qnorm(0.95, y2pred[i], 1)
}

#Figure 3.1
plot( fitSpline, col="blue", xlab="X", ylab="Y2", ylim=c(min(y2), max(y2)))
points(x, y2)
lines(x, cond_quan, col = "red")

#-------------------------------------------------------------------------------

#Subsection 3.2.1

#load the dataset
#may consider changing the location of file
X2021 <- read_delim("C:/Users/Fuaddin/Documents/2021.csv", delim = ";", 
                    escape_double = FALSE, trim_ws = TRUE)

#data manipulation process
data <- na.omit(data.frame(y=X2021$O3, x=X2021$Temperature))
data <- data[data$y > 0,]
n <- dim(data)[1]
#observe the data
head(data)

#Figure 3.2
plot(data$x, data$y, ylab = "Temperature", xlab = "O3")

x <- data$x
y <- data$y
#create 20 cubic spline basis transformations of the observations x
sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data, knots=NULL)[[1]]
X <- sm$X

#fit the cubic regression spline model
fit <- gam(y ~ X - 1, data=data)

#Figure 3.3
plot(data$x, data$y, xlab='Temperature', ylab='O3', col='grey')
lines(data$x[order(data$x)], X[order(data$x),]%*%fit$coefficients, lwd=2, col='blue')

#calculating the parameter for the quantile estimation
yPred <- predict(fit)
condQuan <- rep(0, length(yPred))
sigma <- sqrt(sum((data$y - yPred)^2) / length(yPred))
#estimate the quantile using parametric method
for (i in 1:length(yPred)){
  condQuan[i] <- qnorm(0.95, yPred[i], sigma)}
plot(data$x, data$y, xlab='Temperature', ylab='O3', col='grey', ylim=c(-2, 150))
lines(data$x[order(data$x)], condQuan[order(data$x)], lwd=2, col='red')

#-------------------------------------------------------------------------------

#Subsection 3.3.1

#significance level
a <- 0.05

#calculate the statistics and store as p_hat
w <- rep (0 , length (y))
for (i in 1:length(y)){
  if (y[i] < condQuan[i]){
    w[i] <- 1
  }
}
p_hat <- sum(w) / length(y)

#calculate the 0.95 wilson score interval
z <- qnorm(1 - a/2, 0, 1)
r1 <- p_hat + z^2 / (2 * n)
r2 <- z * sqrt((p_hat * (1 - p_hat)) / n + z^2 / (4 * n^2))
r3 <- 1 + z^2 / n
CI <- c((r1 - r2) / r3, (r1 + r2) / r3)

print(paste("The confidence interval for p is (",CI[1],", ",CI[2],")"))

#-------------------------------------------------------------------------------

#Subsection 3.4.1

tau <- 0.95
x <- data$x
y <- data$y
sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data, knots=NULL)[[1]]
X <- sm$X
n <- length(y)
iBeta <- fit$coefficients

# takes inputs quantile , basis transformed covariates , response , parameter
loss_pinball <- function(tau, X, y, par){
  # define the quantile estimate mu
  mu <- X %*% par
  # define z as difference between response and estimate
  z <- y - mu
  p <- rep(0, length(y))
  # compute pinball loss in terms of z
  p[z >= 0] <- tau * z[z >= 0]
  p[z < 0] <- (tau - 1) * z[z < 0]
  return(sum(p))
}

# takes inputs quantile , basis transformed covariates , response , parameter
grad_pinball <- function(tau, X, y, par){
  # define the quantile estimate mu
  mu <- X %*% par
  # define z as difference between response and estimate
  z <- y - mu
  g <- rep(0, length(y))
  # compute derivative of pinball loss in terms of z
  g[z >= 0] <- tau
  g[z < 0] <- tau - 1
  return(g)
}

# takes inputs quantile , basis transformed covariates , response , parameter
negGrad_pinball <- function(tau, X, y, par){
  # find the directional derivative for each components of X
  tmp <- lapply(1:length(y),
  function(ii){
      # use negative sign gradient since we want to minimise the function
      a <- -grad_pinball(tau, X[ii,], y[ii], par)
      return( a * X[ii, ])})
  out <- Reduce("+", tmp)
  return(out)
}

fit_pinball <- optim(par = iBeta, loss_pinball, gr = negGrad_pinball, method =
                       "BFGS", tau = tau, X = X, y = data$y)

#Figure 3.5
plot(x, y, xlab = 'Temperature', ylab = 'O3', col = 'grey', ylim=c(-2, 150))
lines(x[order(x)], X[order(x),] %*% fit_pinball$par, col = 'blue', lwd = 2)
lines(data$x[order(data$x)], condQuan[order(data$x)], col='red', lwd = 2)
legend(-5, 130, legend=c("Parametric model", "Non-parametric model"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

#-------------------------------------------------------------------------------

#Subsection 3.4.2

#compute the statistic and Wilson score interval for hypothesis test
a <- 0.05
w <- rep(0, length(y))
for (i in 1:length(y)){
  if (y[i] <= X[i,] %*% fit_pinball$par){
    w[i] <- 1
  }
}

p_hat <- sum(w) / length(y)

z <- qnorm(1 - a, 0, 1)
r1 <- p_hat + z^2 / (2 * n)
r2 <- z * sqrt((p_hat * (1 - p_hat)) / n + z^2 / (4 * n^2))
r3 <- 1 + z^2 / n
CI <- c((r1 - r2) / r3, (r1 + r2) / r3)

print(paste("The confidence interval for p is (",CI[1],", ",CI[2],")"))


# takes input dataframe and number of bins to make
create_bin <- function (data, k){
  x <- data$x
  range <- seq(min (x), max(x))
  # set k number of percentiles to partition data
  percentile <- seq(0 , 100 , 100/k)
  # set the highest percentile to 100
  percentile[length(percentile)] <- 100
  # find k number of evenly spaced percentile using quantile
  breaks <- quantile (range , probs = percentile/100)
  # for the first percentile we minus the minimum observations by one
  # since cut consider open interval to define the lower bound of bin
  breaks [1] <- min(x) - 1
  # the last percentile be the maximum observations
  breaks[length(breaks)] <- max(x)
  # classify the observations into k bins using cut function
  bin <- cut(x, breaks = breaks)
  # append the bin ranges into the original dataframe
  data$bins <- bin
  return(data)
}

# takes input dataframe and predictor coefficient
quantile_bin <- function(data, sm, par){
  y <- data$y
  x <- data$x
  # initialised the bins for observations in dataframe data
  blocks <- sort(unique(data$bins))
  k <- length(blocks)
  # initialised k length of zero vector to store the estimate of p
  p <- rep(0, k)
  for(i in 1:length(blocks)){
    # collect observations from the same bin
    ranged <- data[data$bins == blocks[i],]
    sum_p <- 0
    for(j in 1:dim(ranged)[1]){
      # calculate the predicted conditional quantile of data in the bin
      pred <- PredictMat(sm, data.frame(x = ranged$x[j]))%*%par
      # if observation of response is smaller than predicted quantile
      if(ranged$y[j] < pred){
        # increment the count for sum_z
        sum_p <- sum_p + 1
        }
      }
    # define estimate of p as as sum_z divided by number of data in bin
    p[i] <- sum_p / dim(ranged)[1]
    }
  return (p)
}

# takes input estimate of p , dataframe and significance level
ci <- function(qb_estim, data, a = 0.05){
  k <- length(qb_estim)
  # initialised the bins for observations in dataframe data
  blocks <- sort(unique(data$bins))
  # define the 1 - a/2 quantile of standard normal distribution
  z <- qnorm(1 - a/2, 0, 1)
  # create dataframe to store the estimate of p , number of observations ,
  # ranges of bins , lower bound interval , upper bound interval
  df <- data.frame(quantile_estimate = qb_estim, number = rep(0, k),
                    bins = blocks, ci_low = rep (0, k), ci_up = rep(0, k))
  for(i in 1:k){
    # compute the wilson score interval for every estimate of p
    p <- qb_estim[i]
    n <- dim(data[data$bins == blocks[i],])[1]
    r1 <- p + z^2 / (2 * n)
    r2 <- z * sqrt ((p * (1 - p)) / n + z^2 / (4 * n^2))
    r3 <- 1 + z^2 / n
    low <- (r1 - r2) / r3
    up <- (r1 + r2) / r3
    
    # store the number of observations in bin
    df[i,] $number <- n
    # store the lower bound of wilson score
    df[i,] $ci_low <- low
    # store the upper bound of wilson score
    df[i,] $ci_up <- up
    }
  return (df)
}

# takes input of wilson score interval and quantile
ci_bin <- function(ci_estim, tau){
  k <- dim(ci_estim)[1]
  range <- seq(1, k)
  # plot the estimate of p for each bins
  plot(range, ci_estim$quantile_estimate, ylim = c(min(ci_estim$ci_low), max(ci_estim$ci_up)), xlab = "Bins", ylab = "probability")
  abline(h = tau, col = "red")
  
  for (i in 1:k){
    # plot the wilson score interval for each bins
    arrows(range[i], ci_estim[i,]$ci_low, range[i], ci_estim[i,]$ci_up, angle = 90, code = 3, length = 0.1)
  }
  return (ci_estim)
}

#partition our data into 5 bins
data <- create_bin(data, 5)
# define the basis functions in terms of univariate x
x <- data$x
sm <- smoothCon(s(x, k = 20, bs = "cr"), data = data , knots = NULL)[[1]]
# calculate the statistics in each bin
qb_pinball <- quantile_bin(data, sm, fit_pinball$par)
# compute Wilson score interval for p
df_pinball <- ci(qb_pinball, data, 0.05)
# plot the results
ci_bin(df_pinball, tau)

#-------------------------------------------------------------------------------

#Section 3.5

#define index for training samples
train <- sample(1:dim(data)[1], size = 0.7 * length(x), replace = FALSE)

#train-test split
data_train <- data[train,]
data_test <- data[-train,]

x_train <- data_train$x
y_train <- data_train$y

#initialised test set
x_test <- data_test$x
y_test <- data_test$y

#tranform x into cubic spline basis term
sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data_train, knots=NULL)[[1]]
X_train <- sm$X 

tau <- 0.95
#initialiser for predictor coefficients
iBeta <- fit$coefficients

#train model using training data
train_pinball <- optim(par = iBeta, loss_pinball, gr = negGrad_pinball,
                       method = "BFGS", tau = tau, X = X_train, y = y_train)

#to transform univariate in test set into cubic spline terms
X_test <- PredictMat(sm, data.frame(x=x_test))
#calculate the pinball loss between fitted values and actual response in test set
pinball_non <- loss_pinball(tau, X_test, y_test, train_pinball$par) / length(y_test)
#return the values
pinball_non

#Figure 3.7
plot(x_test, y_test, xlab="Temperature", ylab="O3", col = "grey")
lines(x_test[order(x_test)], X_test[order(x_test),] %*% fit_pinball$par, type = 'l', col = "blue", lwd = 2)
lines(x_test[order(x_test)], X_test[order(x_test),] %*% train_pinball$par, type = 'l', col = "red", lwd = 2)
legend(-5, 130, legend=c("Training data", "Full data"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

#set the number of bin samples to 1 since we want to check for whole test set
data_single_test <- create_bin(data_test, 1)
#quantile_bin function to compute the statistics
#ci function to compute the Wilson score interval
ci(quantile_bin(data_single_test, sm, train_pinball$par), data_single_test, 0.05)

#hypothesis test for bin samples
qbt_pinball <- quantile_bin(data_test, sm, train_pinball$par)
dft_pinball <- ci(qbt_pinball, data_test, 0.05)
#Table 3.1 and Figure 3.8
ci_bin(dft_pinball, tau)


#-------------------------------------------------------------------------------

#Subsection 3.6.1

y <- rnorm(1000, 0, 1)
tau <- 0.95
true_quantile <- qnorm(tau, 0, 1)
pinball <- rep(0, length(y))

for(i in 1:length(y)){
  z <- y[i] - true_quantile
  if(z < 0){
    pinball[i] <- -(1 - tau) * z
  }
  else{
    pinball[i] <- tau * z
  }
}

sig <- 3
lam <- 0.1
log_elf <- rep(0, length(y))
for(i in 1:length(y)){
  z <- (y[i] - true_quantile) / sig
  log_elf[i] <- (1-tau) *  - lam * log1pexp( z / lam )
}

#Figure 3.10
hist(y, probability = TRUE, main = "Y~N(0,1)")
lines(y[order(y)], pinball[order(y)], lwd = 2, col = 'blue')
lines(y[order(y)], - log_elf[order(y)] * 10, lwd = 2, col = 'red')
abline(v = true_quantile, lty = 2, lwd = 2, col = 'darkgreen')
legend(-3, 0.3, legend=c("Pinball loss", "Negative log ELF", "0.95-th quantile"),
       col=c("blue", "red", "darkgreen"), lty = 1:1:2, cex = 0.7)

#-------------------------------------------------------------------------------

#Subsection 3.6.3

#computes ELF density and its derivatives w.r.t mu
dlf <- function(y, tau, mu, lam, log = FALSE, deriv = 0)
{
  #calculate the log of ELF density function 
  sig <- 1 #fixed sigma to 1
  z <- (y - mu) / sig
  out <- (1-tau) * z - lam * log1pexp( z / lam ) - 
    log( sig * lam * beta(lam*(1-tau), tau*lam) )
  #if not log return the ELF density function
  if( !log ) out <- exp(out)
  #if deriv is positive compute the derivative of the log ELF density
  if( deriv > 0 )
  {
    out <- list("d" = out)
    dl <- dlogis(y, mu, lam*sig)
    pl <- plogis(y, mu, lam*sig)
    out$D <- sum((pl - (1-tau)) / sig)
    #if deriv is more than one compute the second derivative  
    if(deriv > 1)
    {
      out$D2 <- sum( - dl / sig )
    }
  }
  return(out)
}


#takes input hyperparameter, predictor coefficient, positive semidefinite
penalty <- function(k, par, S = sm$S[[1]], deriv = 0)
{
  #compute the penalty component in the generalised pinball loss
  w <- k * t(par) %*% S %*% par
  w <- list("w" = w)
  #if deriv is positive compute the derivative of the penalty component
  if (deriv > 0)
  {
    w$D <- 2 * k * t(par) %*% S
  }
  return(w)
}


#takes input predictor coefficient, quantile, weight, covariates, 
#response observation, semi positive definite matrix and hyperparameter
penalised_negllkFun <- function(par, tau, lam, X, y, S = sm$S[[1]], k = 0)
{
  #compute the conditional quantile estimate
  mu <- X %*% par
  #compute the semi-parametric objective function
  out <- - sum( dlf(y = y, tau = tau, mu = mu, lam = lam, log = TRUE)) 
  + penalty(k, par, S)$w
  return(out)
}

#takes input predictor coefficient, quantile, weight, covariates, 
#response observation, semi positive definite matrix and hyperparameter
penalised_negGrad <- function(par, tau, lam, X, y, S = sm$S[[1]], k = 0)
{
  #compute the conditional quantile estimate
  mu <- X %*% par
  #compute the descending gradient of semi-parametric objective function
  tmp <- lapply(1:length(mu),
                function(ii){
                  a <- - dlf(y = y[ii], tau = tau, 
                             mu = mu[ii], lam = lam, log = TRUE, deriv = 1)$D 
                  
                  return( a * X[ii, ] + penalty(k, par, S, deriv = 1)$D )})
  out <- Reduce("+", tmp)
  return(out)
}

#-------------------------------------------------------------------------------

#Algorithm 1 (normal distribution)

# to compute the negative log likelihood function for normal assumption
negative_log_likelihood_normal <- function(params, y){
  mu <- params[1]
  sigma <- params[2]
  # Calculate negative log-likelihood
  neg_log <- -sum(dnorm(y, mean = mu, sd = sigma, log = TRUE))
  
  return(neg_log)
}

# Function to estimate ML parameters of normal distribution using optim
estimate_normal_mle <- function(y){
  # Initial guess for parameters
  par_0 <- c(mu = 1, sigma = 1)
  # Minimize negative log-likelihood using optim
  result <- optim(par = par_0, fn = negative_log_likelihood_normal, y = y, method = "L-BFGS-B")
  # Extract estimated parameters
  par_mle <- result$par
  
  return(par_mle)
}

# to compute the weight parameter of the ELF loss function by normal assumption
lambda_normal <- function(tau, y){
  par <- estimate_normal_mle(y)
  n <- length(y)
  # a <- 0.5(mean(y) - mean(log(y)))
  mu_tau <- qnorm(tau, par[1], par[2])
  f_mu0 <- dnorm(mu_tau, par[1], par[2])
  f_mu1 <- dnorm(mu_tau + 1e-10, par[1], par[2])
  #using numerical method to estimate parameter of distribution in the pdf 
  f_gau_prime <- (f_mu0 - f_mu1) / 1e-10
  out <- (9 / (n * pi^4) * f_mu0 / f_gau_prime^2 ) ^ (1/3)
  
  return(out)
}

#-------------------------------------------------------------------------------

#Algorithm 1 (gamma distribution)

# to compute the negative log likelihood function for gamma assumption
negative_log_likelihood_gamma <- function(params, y) {
  a <- params[1]
  b <- params[2]
  # Calculate negative log-likelihood
  neg_log <- -sum(dgamma(y, shape = a, rate = b, log = TRUE))
  
  return(neg_log)
}

# Function to estimate ML parameters of gamma distribution using optim
estimate_gamma_mle <- function(y) {
  # Initial guess for parameters
  par_0 <- c(a = 10, b = 10)
  # Minimize negative log-likelihood using optim
  result <- optim(par = par_0, fn = negative_log_likelihood_gamma, y = y, method = "L-BFGS-B", lower=c(0.1,0.1))
  # Extract estimated parameters
  par_mle <- result$par
  
  return(par_mle)
}

# to compute the weight parameter of the ELF loss function by gamma assumption
lambda_gamma <- function(tau, y){
  par <- estimate_gamma_mle(y)
  n <- length(y)
  mu_tau <- qgamma(tau, par[1], par[2])
  f_mu0 <- dgamma(mu_tau, par[1], par[2])
  f_mu1 <- dgamma(mu_tau + 1e-10, par[1], par[2])
  #using numerical method to estimate parameter of distribution in the pdf 
  f_gau_prime <- (f_mu0 - f_mu1) / 1e-10
  out <- (9 / (n * pi^4) * f_mu0 / f_gau_prime^2 ) ^ (1/3)
  
  return(out)
}

#-------------------------------------------------------------------------------

#Subsection 3.6.4

#Algorithm 2 (gamma distribution)

#refer to pseudo code of Algorithm 2
#takes input covariates, response, initial predictors, hyperparameters candidates
cvFive_gamma <- function(X, y, tau, par, hp_can){
  #line 1
  cv <- rep(1e10, length(hp_can))
  #determine 5-fold cv
  loop <- seq(1, dim(X)[1], floor(dim(X)[1]/5))
  loop <- append(loop, dim(X)[1])
  #line 2
  for(k in 1:length(hp_can)){
    sum_cv <- 0
    #line 3
    for (i in 1:(length(loop) - 1)){
      index <- seq(-loop[i], -loop[i + 1])
      #line 4
      l_gamma <- lambda_gamma(tau, y[index])
      #line 7
      fit <- optim(par = par, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau = tau, 
                    lam = l_gamma, X = X[index, ], y = y[index], k = hp_can[k])
      #line 8
      sum_cv <- sum_cv - sum(dlf(y[-index], tau, X[-index, ] %*% fit$par, l_gamma, log = TRUE)) / length(y[-index])
    }
    #line 9
    cv[k] <- sum_cv / length(loop)
  }
  #line 10
  return(cv)
}


#Algorithm 2 (normal distribution)

#takes input covariates, response, initial predictors, hyperparameters candidates
cvFive_normal <- function(X, y, tau, par, hp_can){
  cv <- rep(1e10, length(hp_can))
  loop <- seq(1, dim(X)[1], floor(dim(X)[1]/5))
  loop <- append(loop, dim(X)[1])
  for(k in 1:length(hp_can)){
    sum_cv <- 0
    for (i in 1:(length(loop) - 1)){
      index <- seq(-loop[i], -loop[i + 1])
      l_normal <- lambda_normal(tau, y[index])
      fit <- optim(par = par, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau = tau, 
                   lam = l_normal, X = X[index, ], y = y[index], k = hp_can[k])
      sum_cv <- sum_cv - sum(dlf(y[-index], tau, X[-index, ] %*% fit$par, l_normal, log = TRUE)) / length(y[-index])
    }
    cv[k] <- sum_cv / length(loop)
  }
  return(cv)
}

#-------------------------------------------------------------------------------

#initialised important values
tau <- 0.95
y <- data$y
x <- data$x
sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data, knots=NULL)[[1]]
X <- sm$X
#tuning hyperparameter candidates
hp_can <- c(5*1e-6, 1e-5, 5*1e-5, 1e-4, 5*1e-4, 1e-3, 5*1e-3, 1e-2)

#five-fold cv to choose hyperparameter for gamma distribution
cv_gamma <- cvFive_gamma(X, y, tau, fit$coefficients, hp_can)
#five-fold cv to choose hyperparameter for normal distribution
cv_normal <- cvFive_normal(X, y, tau, fit$coefficients, hp_can)
total_cv <- c(cv_gamma, cv_normal)

#Figure 3.11
plot(log(hp_can), cv_normal, type = "b", xlab = "log(gamma)", ylab = "CV", col="red", ylim = c(min(total_cv), max(total_cv)))
lines(log(hp_can), cv_gamma, type = "b", col="blue")
legend(-12, 4.078, legend=c("Normal Distribution", "Gamma Distribution"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

#compute the weight parameter for gamma assumption
l_gamma <- lambda_gamma(tau, y)
#compute the weight parameter for normal assumption
l_normal <- lambda_normal(tau, y)

#extract the best tuning hyperparameter
best_hp_gamma <- hp_can[which.min(cv_gamma)]
#fit semi-parametric model using the hyperparameter
best_gamma <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau = tau,
                    lam = l_gamma, X = X, y = y, k = best_hp_gamma)
best_hp_normal <- hp_can[which.min(cv_normal)]
best_normal <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau =tau,
                     lam = l_normal, X = X, y = y, k = best_hp_normal)

#Figure 3.12
plot(x[order(x)], y[order(x)], xlab="Temperature", ylab="O3", col = "grey")
#plotting the quantile estimate
lines(x[order(x)], X[order(x),]%*%best_gamma$par, type = 'l', col = "red", lwd=2)
lines(x[order(x)], X[order(x),]%*%best_normal$par, type = 'l', col = "blue", lwd=2)
legend(-5, 120, legend=c("Gamma Distribution", "Normal Distribution"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

#-------------------------------------------------------------------------------

#Subsection 3.6.5

data_single <- create_bin(data, 1)
#best_gamma$par is the semi-parametric model coefficients for gamma assumption
ci(quantile_bin(data_single, sm, best_gamma$par), data_single, 0.05)
#best_normal$par is the semi-parametric model coefficients for normal assumption
ci(quantile_bin(data_single, sm, best_normal$par), data_single, 0.05)

#hypothesis test using gamma assumption
gamma_best_qb <- quantile_bin(data, sm, best_gamma$par)
df_gamma <- ci(gamma_best_qb, data, 0.05)
#Figure 3.13(b)
ci_bin(df_gamma, tau)

#hypothesis test using normal assumption
normal_best_qb <- quantile_bin(data, sm, best_normal$par)
df_normal <- ci(normal_best_qb, data, 0.05)
#Figure 3.13(a)
ci_bin(df_normal, tau)

#-------------------------------------------------------------------------------

#Subsection 3.6.6

#train-test split data
train <- sample(1:dim(data)[1], size = 0.7 * length(x), replace = FALSE)

data_train <- data[train,]
data_test <- data[-train,]
data_single_test <- data_single[-train,]

x_train <- data_train$x
y_train <- data_train$y

x_test <- data_test$x
y_test <- data_test$y

sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data_train, knots=NULL)[[1]]
X_train <- sm$X 

#compute weight parameter for both distribution assumptions
l_gamma <- lambda_gamma(tau, y_train)
l_normal <- lambda_normal(tau, y_train)

#train semi-parametric model using training data
train_gamma <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau = tau,
                    lam = l_gamma, X = X_train, y = y_train, k = best_hp_gamma)
train_normal <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau =tau,
                     lam = l_normal, X = X_train, y = y_train, k = best_hp_normal)

#Figure 3.14
X_test <- PredictMat(sm, data.frame(x=x_test))
plot(x_test[order(x_test)], y_test[order(x_test)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
#plotting the quantile estimate
lines(x_test[order(x_test)], X_test[order(x_test),]%*%train_gamma$par, type = 'l', col = "red", lwd=2)
lines(x_test[order(x_test)], X_test[order(x_test),]%*%train_normal$par, type = 'l', col = "blue", lwd=2)
legend(-5, 120, legend=c("Gamma Distribution", "Normal Distribution"),
       col=c("red", "blue"), lty=1:1, cex=0.8)


#train_gamma$par is the semi-parametric model coefficients for gamma assumption
ci(quantile_bin(data_single_test, sm, train_gamma$par), data_single_test, 0.05)
#train_normal$par is the semi-parametric model coefficients for normal assumption
ci(quantile_bin(data_single_test, sm, train_normal$par), data_single_test, 0.05)

#hypothesis test for bin samples
gamma_train_qb <- quantile_bin(data_test, sm, train_gamma$par)
df_test_gamma <- ci(gamma_train_qb, data_test, 0.05)
ci_bin(df_test_gamma, tau)

normal_train_qb <- quantile_bin(data_test, sm, train_normal$par)
df_test_normal <- ci(normal_train_qb, data_test, 0.05)
ci_bin(df_test_normal, tau)

#total pinball loss values
pinball_gamma <- loss_pinball(tau, X_test, y_test, train_gamma$par) / length(y_test)
pinball_gamma
pinball_normal <- loss_pinball(tau, X_test, y_test, train_normal$par) / length(y_test)
pinball_normal

#find the minimum cross validation score
min(cv_gamma)
min(cv_normal)

#-------------------------------------------------------------------------------

#Appendix B

#Simulation for 0.90-th semi-parametric quantile regression model
tau1 <- 0.9
y <- data$y
x <- data$x
sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data, knots=NULL)[[1]]
X <- sm$X
hp_can1 <- c(5*1e-5, 1e-4, 5*1e-4, 1e-3, 5*1e-3, 1e-2, 5e-2, 1e-1)
loop <- seq(1, dim(X)[1], floor(dim(X)[1]/5))
loop <- append(loop, dim(X)[1])

cv_gamma1 <- cvFive_gamma(X, y, tau1, fit$coefficients, hp_can1)
cv_normal1 <- cvFive_normal(X, y, tau1, fit$coefficients, hp_can1)
total_cv1 <- c(cv_gamma1, cv_normal1)

#Figure B.1
plot(log(hp_can1), cv_normal1, type = "b", xlab = "log(gamma)", ylab = "CV", col="red", ylim = c(min(total_cv1), max(total_cv1)))
lines(log(hp_can1), cv_gamma1, type = "b", col="blue")
legend(-5, 6.90, legend=c("Normal Distribution", "Gamma Distribution"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

best_hp_gamma1 <- hp_can1[which.min(cv_gamma1)]
best_hp_gamma1
best_hp_normal1 <- hp_can1[which.min(cv_normal1)]
best_hp_normal1

set.seed(1)
train1 <- sample(1:dim(data)[1], size = 0.7 * length(x), replace = FALSE)

data_train1 <- data[train1,]
data_test1 <- data[-train1,]
data_single_test1 <- data_single[-train1,]

x_train1 <- data_train1$x
y_train1 <- data_train1$y

x_test1 <- data_test1$x
y_test1 <- data_test1$y

sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data_train1, knots=NULL)[[1]]
X_train1 <- sm$X 

l_gamma1 <- lambda_gamma(tau1, y_train1)
l_normal1 <- lambda_normal(tau1, y_train1)


train_gamma1 <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau = tau1,
                     lam = l_gamma1, X = X_train1, y = y_train1, k = best_hp_gamma1)
train_normal1 <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau =tau1,
                      lam = l_normal1, X = X_train1, y = y_train1, k = best_hp_normal1)

#Figure B.1
X_test1 <- PredictMat(sm, data.frame(x=x_test1))
plot(x_test1[order(x_test1)], y_test1[order(x_test1)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
#plotting the quantile estimate
lines(x_test1[order(x_test1)], X_test1[order(x_test1),]%*%train_gamma1$par, type = 'l', col = "red", lwd=2)
legend(-5, 120, legend="Gamma Distribution", col="red", lty=1, cex=0.8)

plot(x_test1[order(x_test1)], y_test1[order(x_test1)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
lines(x_test1[order(x_test1)], X_test1[order(x_test1),]%*%train_normal1$par, type = 'l', col = "blue", lwd=2)
legend(-5, 120, legend="Normal Distribution", col="blue", lty=1, cex=0.8)

#Figure B.2
gamma_train_qb1 <- quantile_bin(data_test1, sm, train_gamma1$par)
df_test_gamma1 <- ci(gamma_train_qb1, data_test1, 0.05)
ci_bin(df_test_gamma1, tau1)

#Figure B.2
normal_train_qb1 <- quantile_bin(data_test1, sm, train_normal1$par)
df_test_normal1 <- ci(normal_train_qb1, data_test1, 0.05)
ci_bin(df_test_normal1, tau1)


l_gamma1
l_normal1

best_hp_gamma1
best_hp_normal1

#train_gamma$par is the semi-parametric model coefficients for gamma assumption
ci(quantile_bin(data_single_test1, sm, train_gamma1$par), data_single_test1, 0.05)
#train_normal$par is the semi-parametric model coefficients for normal assumption
ci(quantile_bin(data_single_test1, sm, train_normal1$par), data_single_test1, 0.05)

pinball_gamma1 <- loss_pinball(tau1, X_test1, y_test1, train_gamma1$par) / length(y_test1)
pinball_gamma1
pinball_normal1 <- loss_pinball(tau1, X_test1, y_test1, train_normal1$par) / length(y_test1)
pinball_normal1

min(cv_gamma1)
min(cv_normal1)
#-------------------------------------------------------------------------------

#Simulation for 0.85-th semi-parametric quantile regression model

tau2 <- 0.85
y <- data$y
x <- data$x
sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data, knots=NULL)[[1]]
X <- sm$X
hp_can2 <- c(5*1e-5, 1e-4, 5*1e-4, 1e-3, 5*1e-3, 1e-2, 5e-2, 1e-1)
loop <- seq(1, dim(X)[1], floor(dim(X)[1]/5))
loop <- append(loop, dim(X)[1])

cv_gamma2 <- cvFive_gamma(X, y, tau2, fit$coefficients, hp_can2)
cv_normal2 <- cvFive_normal(X, y, tau2, fit$coefficients, hp_can2)
total_cv2 <- c(cv_gamma2, cv_normal2)

#Figure B.3
plot(log(hp_can2), cv_normal2, type = "b", xlab = "log(gamma)", ylab = "CV", col="red", ylim = c(min(total_cv2), max(total_cv2)))
lines(log(hp_can2), cv_gamma2, type = "b", col="blue")
legend(-10, 5.42, legend=c("Normal Distribution", "Gamma Distribution"),
       col=c("red", "blue"), lty=1:1, cex=0.8)


best_hp_gamma2 <- hp_can2[which.min(cv_gamma2)]
best_hp_normal2 <- hp_can2[which.min(cv_normal2)]

set.seed(2)
train2 <- sample(1:dim(data)[1], size = 0.7 * length(x), replace = FALSE)

data_train2 <- data[train2,]
data_test2 <- data[-train2,]
data_single_test2 <- data_single[-train2,]

x_train2 <- data_train2$x
y_train2 <- data_train2$y

x_test2 <- data_test2$x
y_test2 <- data_test2$y

sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data_train2, knots=NULL)[[1]]
X_train2 <- sm$X 

l_gamma2 <- lambda_gamma(tau2, y_train2)
l_normal2 <- lambda_normal(tau2, y_train2)

train_gamma2 <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau = tau2,
                      lam = l_gamma2, X = X_train2, y = y_train2, k = best_hp_gamma2)
train_normal2 <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau =tau2,
                       lam = l_normal2, X = X_train2, y = y_train2, k = best_hp_normal2)

#Figure B.3
X_test2 <- PredictMat(sm, data.frame(x=x_test2))
plot(x_test2[order(x_test2)], y_test2[order(x_test2)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
#plotting the quantile estimate
lines(x_test2[order(x_test2)], X_test2[order(x_test2),]%*%train_gamma2$par, type = 'l', col = "red", lwd=2)
legend(-5, 120, legend="Gamma Distribution", col="red", lty=1, cex=0.8)

plot(x_test2[order(x_test2)], y_test2[order(x_test2)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
lines(x_test2[order(x_test2)], X_test2[order(x_test2),]%*%train_normal2$par, type = 'l', col = "blue", lwd=2)
legend(-5, 120, legend="Normal Distribution", col="blue", lty=1, cex=0.8)

#Figure B.4
gamma_train_qb2 <- quantile_bin(data_test2, sm, train_gamma2$par)
df_test_gamma2 <- ci(gamma_train_qb2, data_test2, 0.05)
ci_bin(df_test_gamma2, tau2)

#Figure B.4
normal_train_qb2 <- quantile_bin(data_test2, sm, train_normal2$par)
df_test_normal2 <- ci(normal_train_qb2, data_test2, 0.05)
ci_bin(df_test_normal2, tau2)


l_gamma2
l_normal2

best_hp_gamma2
best_hp_normal2

#train_gamma$par is the semi-parametric model coefficients for gamma assumption
ci(quantile_bin(data_single_test2, sm, train_gamma2$par), data_single_test2, 0.05)
#train_normal$par is the semi-parametric model coefficients for normal assumption
ci(quantile_bin(data_single_test2, sm, train_normal2$par), data_single_test2, 0.05)

pinball_gamma2 <- loss_pinball(tau2, X_test2, y_test2, train_gamma2$par) / length(y_test2)
pinball_gamma2
pinball_normal2 <- loss_pinball(tau2, X_test2, y_test2, train_normal2$par) / length(y_test2)
pinball_normal2

min(cv_gamma2)
min(cv_normal2)
#-------------------------------------------------------------------------------

#Simulation for 0.80-th semi-parametric quantile regression model

tau3 <- 0.8
y <- data$y
x <- data$x
sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data, knots=NULL)[[1]]
X <- sm$X
hp_can3 <- c(5*1e-5, 1e-4, 5*1e-4, 1e-3, 5*1e-3, 1e-2, 5e-2, 1e-1)
loop <- seq(1, dim(X)[1], floor(dim(X)[1]/5))
loop <- append(loop, dim(X)[1])

cv_gamma3 <- cvFive_gamma(X, y, tau3, fit$coefficients, hp_can3)
cv_normal3 <- cvFive_normal(X, y, tau3, fit$coefficients, hp_can3)
total_cv3 <- c(cv_gamma3, cv_normal3)

#Figure B.5
plot(log(hp_can3), cv_normal3, type = "b", xlab = "log(gamma)", ylab = "CV", col="red", ylim = c(min(total_cv3), max(total_cv3)))
lines(log(hp_can3), cv_gamma3, type = "b", col="blue")
legend(-10, 5.9, legend=c("Normal Distribution", "Gamma Distribution"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

best_hp_gamma3 <- hp_can3[which.min(cv_gamma3)]
best_hp_normal3 <- hp_can3[which.min(cv_normal3)]

set.seed(3)
train3 <- sample(1:dim(data)[1], size = 0.7 * length(x), replace = FALSE)

data_train3 <- data[train3,]
data_test3 <- data[-train3,]
data_single_test3 <- data_single[-train3,]

x_train3 <- data_train3$x
y_train3 <- data_train3$y

x_test3 <- data_test3$x
y_test3 <- data_test3$y

sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data_train3, knots=NULL)[[1]]
X_train3 <- sm$X 

l_gamma3 <- lambda_gamma(tau3, y_train3)
l_normal3 <- lambda_normal(tau3, y_train3)

train_gamma3 <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau = tau3,
                      lam = l_gamma3, X = X_train3, y = y_train3, k = best_hp_gamma3)
train_normal3 <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau =tau3,
                       lam = l_normal3, X = X_train3, y = y_train3, k = best_hp_normal3)

#Figure B.5
X_test3 <- PredictMat(sm, data.frame(x=x_test3))
plot(x_test3[order(x_test3)], y_test3[order(x_test3)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
#plotting the quantile estimate
lines(x_test3[order(x_test3)], X_test3[order(x_test3),]%*%train_gamma3$par, type = 'l', col = "red", lwd=2)
legend(-5, 120, legend="Gamma Distribution", col="red", lty=1, cex=0.8)

plot(x_test3[order(x_test3)], y_test3[order(x_test3)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
lines(x_test3[order(x_test3)], X_test3[order(x_test3),]%*%train_normal3$par, type = 'l', col = "blue", lwd=2)
legend(-5, 120, legend="Normal Distribution", col="blue", lty=1, cex=0.8)

#Figure B.6
gamma_train_qb3 <- quantile_bin(data_test3, sm, train_gamma3$par)
df_test_gamma3 <- ci(gamma_train_qb3, data_test3, 0.05)
ci_bin(df_test_gamma3, tau3)

#Figure B.6
normal_train_qb3 <- quantile_bin(data_test3, sm, train_normal3$par)
df_test_normal3 <- ci(normal_train_qb3, data_test3, 0.05)
ci_bin(df_test_normal3, tau3)


l_gamma3
l_normal3

best_hp_gamma3
best_hp_normal3

#train_gamma$par is the semi-parametric model coefficients for gamma assumption
ci(quantile_bin(data_single_test3, sm, train_gamma3$par), data_single_test3, 0.05)
#train_normal$par is the semi-parametric model coefficients for normal assumption
ci(quantile_bin(data_single_test3, sm, train_normal3$par), data_single_test3, 0.05)

pinball_gamma3 <- loss_pinball(tau3, X_test3, y_test3, train_gamma3$par) / length(y_test3)
pinball_gamma3
pinball_normal3 <- loss_pinball(tau3, X_test3, y_test3, train_normal3$par) / length(y_test3)
pinball_normal3

min(cv_gamma3)
min(cv_normal3)

#-------------------------------------------------------------------------------

#Simulation for 0.75-th semi-parametric quantile regression model

tau4 <- 0.75
y <- data$y
x <- data$x
sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data, knots=NULL)[[1]]
X <- sm$X
hp_can4 <- c(5*1e-5, 1e-4, 5*1e-4, 1e-3, 5*1e-3, 1e-2, 5e-2, 1e-1)
loop <- seq(1, dim(X)[1], floor(dim(X)[1]/5))
loop <- append(loop, dim(X)[1])

cv_gamma4 <- cvFive_gamma(X, y, tau4, fit$coefficients, hp_can4)
cv_normal4 <- cvFive_normal(X, y, tau4, fit$coefficients, hp_can4)
total_cv4 <- c(cv_gamma4, cv_normal4)

#Figure B.7
plot(log(hp_can4), cv_normal4, type = "b", xlab = "log(gamma)", ylab = "CV", col="red", ylim = c(min(total_cv4), max(total_cv4)))
lines(log(hp_can4), cv_gamma4, type = "b", col="blue")
legend(-10, 6.70, legend=c("Normal Distribution", "Gamma Distribution"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

best_hp_gamma4 <- hp_can4[which.min(cv_gamma4)]
best_hp_normal4 <- hp_can4[which.min(cv_normal4)]

set.seed(4)
train4 <- sample(1:dim(data)[1], size = 0.7 * length(x), replace = FALSE)

data_train4 <- data[train4,]
data_test4 <- data[-train4,]
data_single_test4 <- data_single[-train4,]

x_train4 <- data_train4$x
y_train4 <- data_train4$y

x_test4 <- data_test4$x
y_test4 <- data_test4$y

sm <- smoothCon(s(x, k = 20, bs = "cr"), data=data_train4, knots=NULL)[[1]]
X_train4 <- sm$X 

l_gamma4 <- lambda_gamma(tau4, y_train4)
l_normal4 <- lambda_normal(tau4, y_train4)

train_gamma4 <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau = tau4,
                      lam = l_gamma4, X = X_train4, y = y_train4, k = best_hp_gamma4)
train_normal4 <- optim(par = fit$coefficients, penalised_negllkFun, gr = penalised_negGrad, method = "BFGS", tau =tau4,
                       lam = l_normal4, X = X_train4, y = y_train4, k = best_hp_normal4)

#Figure B.7
X_test4 <- PredictMat(sm, data.frame(x=x_test4))
plot(x_test4[order(x_test4)], y_test4[order(x_test4)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
#plotting the quantile estimate
lines(x_test4[order(x_test4)], X_test4[order(x_test4),]%*%train_gamma4$par, type = 'l', col = "red", lwd=2)
legend(-5, 120, legend="Gamma Distribution", col="red", lty=1, cex=0.8)

plot(x_test4[order(x_test4)], y_test4[order(x_test4)], xlab="Temperature", ylab="O3", col = "grey", ylim=c(-2, 140))
lines(x_test4[order(x_test4)], X_test4[order(x_test4),]%*%train_normal4$par, type = 'l', col = "blue", lwd=2)
legend(-5, 120, legend="Normal Distribution", col="blue", lty=1, cex=0.8)

#Figure B.8
gamma_train_qb4 <- quantile_bin(data_test4, sm, train_gamma4$par)
df_test_gamma4 <- ci(gamma_train_qb4, data_test4, 0.05)
ci_bin(df_test_gamma4, tau4)

#Figure B.8
normal_train_qb4 <- quantile_bin(data_test4, sm, train_normal4$par)
df_test_normal4 <- ci(normal_train_qb4, data_test4, 0.05)
ci_bin(df_test_normal4, tau4)


l_gamma4
l_normal4

best_hp_gamma4
best_hp_normal4

#train_gamma$par is the semi-parametric model coefficients for gamma assumption
ci(quantile_bin(data_single_test4, sm, train_gamma4$par), data_single_test4, 0.05)
#train_normal$par is the semi-parametric model coefficients for normal assumption
ci(quantile_bin(data_single_test4, sm, train_normal4$par), data_single_test4, 0.05)

pinball_gamma4 <- loss_pinball(tau4, X_test4, y_test4, train_gamma4$par) / length(y_test4)
pinball_gamma4
pinball_normal4 <- loss_pinball(tau4, X_test4, y_test4, train_normal4$par) / length(y_test4)
pinball_normal4

min(cv_gamma4)
min(cv_normal4)