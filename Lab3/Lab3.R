library(ggplot2)
library(mvtnorm)
library(rstan)
set.seed(12345)

# 1 a)

# load data. We need the natural log of the data
data <- log(readRDS("Precipitation.rds"))

# initial values
mu_0 <- mean(data)
nu_0 <- length(data) - 1
sigma_0_2 <- var(data)
tau_0_2 <- sigma_0_2

# implement samplers using full conditional posteriors from Lecture 7, slide 16
sigma_2_sampler <- function(data, nu_0, sigma_0_2, mu){
  N <- length(data)
  nu_n <- N - 1
  s_2 <- (nu_0 * sigma_0_2 + sum((data - mu)^2)) / (N + nu_0)
  x <- rchisq(1, nu_n)
  random_draws_inv_chisquare <- nu_n * s_2 / x
  return (random_draws_inv_chisquare)
}

mu_sampler <- function(data, sigma_n_2, tau_0_2, mu_0){
  n = length(data)
  part1 = n / sigma_n_2
  part2 = part1 + 1 / tau_0_2
  w = part1 / part2
  mu_n <- w * mean(data) + (1-w) * mu_0
  tau_n_2 = 1 / part2
  x <- rnorm(1, mean = mu_n, sd = tau_n_2)
  return(x)
}

# Gibbs sampling
gibbs <- function(nDraws, data, mu_0, nu_0, sigma_0_2, tau_0_2){
  # create empty matrix for draws
  gibbsDraws <- matrix(0,nDraws,3)
  # initial values
  mu <- mu_0
  nu <- nu_0
  sigma_2 <- sigma_0_2
  data_and_samples <- data
  
  # repeat: draw mu, draw sigma^2, draw from joint posterior
  for (i in 1:nDraws){
    
    # Update sigma_2 given previous mu
    sigma_2 <- sigma_2_sampler(data_and_samples, nu_0, sigma_0_2, mu)
    gibbsDraws[i,1] <- sigma_2
    
    # Update mu given new sigma_2
    mu <- mu_sampler(data_and_samples, sigma_2, tau_0_2, mu_0)
    gibbsDraws[i,2] <- mu
    
    # draw from joint posterior
    gibbsDraw <- rnorm(1,mean = mu, sd = sigma_2)
    data_and_samples <- c(data_and_samples, gibbsDraw)
    gibbsDraws[i,3] <- gibbsDraw
  }
  # gibbsDraws is a matrix with nDraws rows and three columns.
  # the columns represent:
  # [1] sigma^2 draws, [2] mu draws, [3] joint posterior draws
  return(gibbsDraws)
}

# draw from gibbs sampler
draws <- gibbs(1000, data, mu_0, nu_0, sigma_0_2, tau_0_2)

plot_df <- data.frame(
  sigma2 = draws[,1],
  mu = draws[,2]
)

ggplot(plot_df) +
  geom_path(aes(x = sigma2, y = mu))


# autocorrelations of sigma^2 draws
a_sigma_2 <- acf(draws[,1])
# calculate Inefficiency Factor
IF_sigma_2 <- 1+2*sum(a_sigma_2$acf[-1])

# autocorrelations of mu draws
a_mu <- acf(draws[,2])
# calculate Inefficiency Factor
IF_mu <- 1+2*sum(a_mu$acf[-1])

# 1 b)
plot_df1 <- data.frame(
  data = exp(data)
)
plot_df_2 <- data.frame(
  draws = exp(draws[,3])
)

ggplot()+
  geom_histogram(data = plot_df1, bins=100, aes(x = data, y=..density..),
                 color ="black", fill="gray") +
  geom_density(data = plot_df_2, aes(x = draws), alpha=.2, fill="#FF6666") +
  xlim(0,100)


# 2 a)
df <- read.table("eBayNumberOfBidderData.dat", header=TRUE)
y <- df[,1]
x <- df[,2:10]

#Store some trivial information for future calculations
parameters <- colnames(x)
x_col = ncol(x)
x_row = nrow(x)

#Change X according to description (Ignore Const)
x_glm <- df[,3:10]

#Use glmmodel with poisson
glmModel <- glm(nBids ~ PowerSeller + VerifyID + Sealed + Minblem + MajBlem + 
                  LargNeg + LogBook + MinBidShare, data = df, family = poisson())


# 2 b)

LogPostPoisson <- function(betas,y,x,mu,Sigma,n){
  #Loglikelihood
  logLik <- sum((y*as.matrix(x)%*%as.matrix(betas))-exp(as.matrix(x)%*%as.matrix(betas)))
  #Prior
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  return(logLik + logPrior)
}

# Select the initial values for beta
initVal <- matrix(0,x_col,1)

sigma_prior <- 100*solve(t(x)%*%as.matrix(x))
mu_prior <- as.matrix(rep(0,x_col))


#Get the optimal beta values
OptimRes <- optim(initVal,LogPostPoisson,gr=NULL,y,x,mu_prior,sigma_prior,n=x_row,method=c("BFGS"),
                  control=list(fnscale=-1),hessian=TRUE)


#Posterior Mode
posterior_mode <- OptimRes$par

#Negative Hessian
negative_hessian = -OptimRes$hessian

#Inverse negative hessian
inverse_negative_hessian = solve(negative_hessian)

#Number of draws
draws = 10000

beta_values <- rmvnorm(draws, mean = posterior_mode, sigma = inverse_negative_hessian)


# 2 c)

#LogPostFunc with starting betas
LogPostFunc <- function(betas, x){
  #Create initial values for prior
  x_col = ncol(x)
  Sigma <- 100*solve(t(x)%*%as.matrix(x))
  mu <- as.matrix(rep(0,x_col))
  #Loglikelihood
  logLik <- sum((y*as.matrix(x)%*%as.matrix(betas))-exp(as.matrix(x)%*%as.matrix(betas)))
  #Prior
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  return(logLik + logPrior)
}
#Metropolis random walk
RWMSampler <- function(draws, c, omega, x,theta, LogPostFunc){
  acc_counter <- 0
  #Create empty matrix for storing thetas
  theta_matrix <- matrix(theta, draws, nrow(omega))
  #For loop to make new draw/theta_value
  for (index in 2:draws){
    #Create new theta according to description
    theta <- as.vector(rmvnorm(n = 1 , mean = theta_matrix[index-1,], sigma = c*omega))
    #Use new theta in the posterior density function
    y <- LogPostFunc(betas = theta, x)
    #Create random value between 0-1
    u <- runif(1)
    #Calculcate alpha according to the slides
    alpha <- exp(y - LogPostFunc(theta_matrix[index-1,], x))
    #If alpha is larger than the random value, store new theta
    if (u <= alpha){
      theta_matrix[index,] <- theta
      acc_counter <- acc_counter + 1
    }
    #Else reuse the theta from earlier draw
    else{
      theta_matrix[index,] <- theta_matrix[index-1,]
    }
  }
  #print(acc_counter/draws)
  return(theta_matrix)
  
}
#Tuning parameter c and number of draws, 
c = 0.2
draws = 1000
#Sample from Metropolis Random walk
beta_values_RWM <- as.data.frame(RWMSampler(draws, c, inverse_negative_hessian,x,0,  LogPostFunc))

#Create dataframe and plot all variables
colnames(beta_values_RWM) <- parameters
beta_values_RWM["iteration"] <- c(1:draws)

#Const
ggplot(beta_values_RWM, aes(x = iteration, y=Const))+
  geom_line()+
  labs(title = "Const ", x = "Iteration", y = "Const Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))


#PowerSeller
ggplot(beta_values_RWM, aes(x = iteration, y=PowerSeller))+
  geom_line()+
  labs(title = "PowerSeller", x = "Iteration", y = "PowerSeller Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(beta_values_RWM, aes(x = iteration, y=PowerSeller))+
  geom_line()+
  labs(title = "PowerSeller", x = "Iteration", y = "PowerSeller Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#VerifyID
ggplot(beta_values_RWM, aes(x = iteration, y=VerifyID))+
  geom_line()+
  labs(title = "VerifyID", x = "Iteration", y = "VerifyID Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#Sealed
ggplot(beta_values_RWM, aes(x = iteration, y=VerifyID))+
  geom_line()+
  labs(title = "VerifyID", x = "Iteration", y = "VerifyID Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#Minblem
ggplot(beta_values_RWM, aes(x = iteration, y=Minblem))+
  geom_line()+
  labs(title = "Minblem", x = "Iteration", y = "Minblem Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#MajBlem
ggplot(beta_values_RWM, aes(x = iteration, y=MajBlem))+
  geom_line()+
  labs(title = "MajBlem", x = "Iteration", y = "MajBlem Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#LargNeg
ggplot(beta_values_RWM, aes(x = iteration, y=LargNeg))+
  geom_line()+
  labs(title = "LargNeg", x = "Iteration", y = "LargNeg Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#LogBook
ggplot(beta_values_RWM, aes(x = iteration, y=LogBook))+
  geom_line()+
  labs(title = "LogBook", x = "Iteration", y = "LogBook Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#MinBidShare
ggplot(beta_values_RWM, aes(x = iteration, y=MinBidShare))+
  geom_line()+
  labs(title = "MinBidShare", x = "Iteration", y = "MinBidShare Value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))


# 2 d)

#Const, Power seller, VerifyID, Sealed, MinBlem, MajBlem, LargNeg, LogBook, MinBidShare
x_pred <- as.matrix(c(1, 1, 0, 1, 0, 1, 0, 1.2, 0.8))
#Beta values after convergence
beta_values <- colMeans(beta_values_RWM[500:draws,1:9])

#Create lambda (for rpois) with new beta values
lambda <- exp(t(x_pred)*beta_values)

#Numer of draws
number_of_draws = 1000

#Get prediction with the 
results <- rpois(number_of_draws, lambda)

#Create dataframe and plot it
results_df <- data.frame(results = results,
                         index = c(1:number_of_draws))
colnames(results_df) <- c("Prediction", "Iteration")

ggplot(results_df, aes(x=Prediction))+
  geom_histogram(bins=30)+
  labs(title = "Predictive distribution.", x = "Bids", y = "Frequency")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#Probability of no bidders
Prob_no_bidders <- length(which(results==0))/number_of_draws


# 3 a)

sampler <- function(T, mu, phi, sigma2){
  x <- rep(mu,T)
  error <- rnorm(T, 0,sigma2)
  for(t in 2:T){
    x[t] <- mu + phi * (x[t-1] - mu) + error[t]
  }
  return(x)
}

mu = 13
sigma2 = 3
phi = 0.5
T=300

samples <- sampler(T, mu, phi, sigma2)

plot(samples)
lines(samples)


# 3 b)

phi = 0.2
samples_02 <- sampler(T, mu, phi, sigma2)
phi = 0.95
samples_095 <- sampler(T, mu, phi, sigma2)

StanModel = '
data {
  int<lower=0> T; // Number of observations
  vector[T] x;
}
parameters {
  real mu;
  real sigma2;
  real phi;
}
model {
  for(t in 2:T){
    x[t] ~ normal(mu + phi * (x[t-1] - mu), sigma2);
  } 
}'

data <- list(T=T, x=samples_02)
warmup <- 1000
niter <- 2000
fit_x <- stan(model_code=StanModel,data=data, warmup=warmup,iter=niter,chains=4,seed=12345)



data <- list(T=T, x=samples_095)
warmup <- 1000
niter <- 2000
fit_y <- stan(model_code=StanModel,data=data, warmup=warmup,iter=niter,chains=4,seed=12345)

# Print the fitted model
print(fit_x,digits_summary=3)

# Print the fitted model
print(fit_y,digits_summary=3)

#Plot x dataset
traceplot(fit_x)
pairs(fit_x, pars = c("mu","phi"))

#Plot y dataset
traceplot(fit_y)
pairs(fit_y, pars = c("mu","phi"))