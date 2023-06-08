library(readxl)
library(tibble)
library(dplyr)
library(mvtnorm)
library(ggplot2)

set.seed(12345)

# 1 ####
# a)

# load data
lin_temp <- read_xlsx("./Linkoping2022.xlsx")
lin_temp <- lin_temp %>% add_column(day = 1:365/365)

# function for inverse chi square draw
# inputs:
# draws: wanted number of draws
# nu: degrees of freedom
# sigma2: standard deviation
# returns:
# vector of random draws from inverse chi square
inv_chi_squared_draw <- function(draws, nu, sigma2){
  x <- rchisq(draws, nu)
  random_draws_inv_chisquare <- nu * sigma2 / x
  return (random_draws_inv_chisquare)
}

# function to calculate polynomial regression
# inputs:
# x: x values
# beta: the beta parameters
# returns:
# data.frame containing:
# x0: the first summand of the polynom
# x1: the second summand of the polynom
# x2: the third summand of the polynom
# day: the input x
# y_reg: a polynomial regression based on the input x
calc_pol_reg <- function(x, beta){
  n <- length(x)
  degree <- length(beta) - 1
  # xp is the matrix used in polynomial regression. It contains
  # vectors of 1, x, x^2, ...
  xp <- matrix(data = 1,nrow=n,ncol=1)
  xp <- cbind(xp,x)
  for(i in 2:degree){
    new_x <- x^i
    xp <- cbind(xp, new_x)
  }
  # for some reason it has to be converted into a matrix again
  xp <- as.matrix(xp)
  # the return data frame
  ret_df <- data.frame(
    x0 = beta[1],
    x1 = xp[,2] * beta[2],
    x2 = xp[,3] * beta[3],
    day = x,
    # this is the actual regression part
    y_reg = xp%*%t(beta)
  )
  return(ret_df)
}
# function to generate curves
# inputs:
# n: number of curves to generate
# x: the x values that the regression is based on
# mu0: vector of means
# omega0: matrix of omega 0
# nu0: degrees of freedom of inverse chi sqare
# sigma02: standard deviation of inverse chi square
# returns:
# data.frame containing:
# x: the input x
# y_reg_1 to y_reg_n: n calculated polynomial regressions
generate_curves <- function(n,x,mu0,omega0,nu0,sigma02){
  ret_df <- data.frame(
    x <- x
  )
  # vector to store names of y_reg_i columns
  df_names <- c()
  # repeat the same process n times
  for (i in 1:n){
    # draw sigma2
    sigma2 <- inv_chi_squared_draw(1, nu0, sigma02)
    # draw beta values
    beta <- rmvnorm(1, mean = mu0, sigma = sigma2 * solve(omega0))
    # calculate polynomial regression
    pol_reg <- calc_pol_reg(x, beta)
    # add to return data.frame
    ret_df <- cbind(ret_df, pol_reg$y_reg)
    # create name for this column in return data.frame
    df_names <- c(df_names,paste0("y_reg_",i))
  }
  # set names
  names(ret_df) <- c("x",df_names)
  return(ret_df)
}

# tweaked parameters
# first mu controls y-axis offset
# second mu controls steepnes of linear part
# third mu controls square part
mu0 <- matrix(c(-12,115,-110),ncol = 1)
omega0 <- diag(0.01, 3)
# as last time, nu is chosen as (n-1)
nu0 <- dim(lin_temp)[1] - 1
# way smaller sigma02
sigma02 <- 0.2

# generate 10 curves to see if they fit
curves <- generate_curves(10,lin_temp$day, mu0,omega0,nu0,sigma02)
# get the names of these curves to plot them in a loop later
y_names <- names(curves)[2:dim(curves)[2]]
# add real temperature values to data frame 
plot_df <- cbind(curves, temp = lin_temp$temp)
# plot real temperature values
curve_plot <- ggplot(plot_df) +
  geom_point(aes(x = x, y = temp))
# plot generated curves in a loop
for(col in 1:length(y_names)){
  curve_plot <- curve_plot + geom_line(aes_string(x="x",y=y_names[col]))
}
# change labels
curve_plot <- curve_plot + xlab("day") + ylab("temperature")
# draw plot
curve_plot


# b)

# number of samples (days)
n <- dim(lin_temp)[1]
# matrix containing vectors of 1, x, x^2
x_mat <- as.matrix(cbind(
  rep(1,n),
  lin_temp$day,
  lin_temp$day ^ 2))

# temperature of samples
y_mat <- as.matrix(lin_temp$temp)
# save x transposed x, because that will be needed multiple times
xtx <- t(x_mat) %*% x_mat

# calculate betahat according to lecture 5, slide 4
betahat <- solve(xtx) %*% t(x_mat) %*% y_mat
# calculate parameters for posterior
mu_n <- solve(xtx + omega0) %*% (xtx %*% betahat + omega0 %*% mu0)
omega_n <- xtx + omega0
nu_n <- nu0 + n
nu_n_sigma_n_2 <- nu0 * sigma02 + (t(y_mat) %*% y_mat + t(mu0) %*% omega0 %*%
                                     mu0 - t(mu_n) %*% omega_n %*% mu_n)
sigma_n_2 <- as.vector(nu_n_sigma_n_2 / nu_n)
sigma2_post <- inv_chi_squared_draw(1, nu_n, sigma_n_2)

# function to get betas and sigma2 values
# inputs:
# n: number of samples
# mu: vector of means
# omega: matrix of omega
# nu: degrees of freedom of inverse chi sqare
# sigma2: standard deviation of inverse chi square
# returns:
# data.frame containing:
# sample_num: sample numbers from 1 to n
# sigma2: n values for sigma2
# beta0: n values for beta0
# beta1: n values for beta1
# beta2: n values for beta2
get_params <- function(n,mu,omega,nu,sigma2){
  beta0_vec <- c()
  beta1_vec <- c()
  beta2_vec <- c()
  sigma2_vec <- c()
  # draw n random samples
  for (i in 1:n){
    # sigma2 draw from inv chi sqare
    sigma2_vec <- c(sigma2_vec, inv_chi_squared_draw(1, nu, sigma2))
    # beta draws from multivariate normal
    beta <- rmvnorm(1, mean = mu, sigma = sigma2_vec[i] * solve(omega))
    beta0_vec <- c(beta0_vec,beta[1])
    beta1_vec <- c(beta1_vec,beta[2])
    beta2_vec <- c(beta2_vec,beta[3])
  }
  # create return data.frame
  ret_df <- data.frame(
    sample_num = 1:n,
    sigma2 = sigma2_vec,
    beta0 = beta0_vec,
    beta1 = beta1_vec,
    beta2 = beta2_vec
  )
  return(ret_df)
}

# set n
n = 10000
# get n draws of sigma and betas
plot_df <- get_params(n,mu_n,omega_n,nu_n,sigma_n_2)

# create histograms
ggplot(plot_df, aes(x=sigma2))+
  geom_histogram(color ="black", fill="gray",bins = 30)+
  labs(title = expression ("Histogram of Sigma"^2), x = expression ("Sigma"^2), y = "Frequency")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_df, aes(x=beta0))+
  geom_histogram(color ="black", fill="gray",bins = 30)+
  labs(title = expression ("Histogram of Beta"[0]), x = expression ("Beta"[0]), y = "Frequency")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_df, aes(x=beta1))+
  geom_histogram(color ="black", fill="gray",bins = 30)+
  labs(title = expression ("Histogram of Beta"[1]), x = expression ("Beta"[1]), y = "Frequency")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_df, aes(x=beta2))+
  geom_histogram(color ="black", fill="gray",bins = 30)+
  labs(title = expression ("Histogram of Beta"[2]), x = expression ("Beta"[2]), y = "Frequency")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))
# 3 rows (the betas), n columns (the number of beta samples we took)
beta_mat <- matrix(
  c(plot_df$beta0,plot_df$beta1,plot_df$beta2), nrow = 3, ncol = n, byrow = TRUE
)
# y_mat has 365 rows (the days in the year),
# n columns (the number of beta samples we took)
y_mat <- x_mat %*% beta_mat

# generate mean, lower 5% and upper 5% boundaries from samples for each day
# of the year
mean_vec <- c()
lower_vec <- c()
upper_vec <- c()
# for each day of the year
for(row in 1:dim(y_mat)[1]){
  # get lower and upper indices
  lower_index <- floor(0.05 * n + 1)
  upper_index <- floor(0.95 * n)
  # calculate mean of values
  values <- y_mat[row,]
  mean_vec <- c(mean_vec, mean(values))
  # calculate lower and upper boundaries from sorted array
  sorted_val <- sort(values)
  lower_vec <- c(lower_vec, sorted_val[lower_index])
  upper_vec <- c(upper_vec, sorted_val[upper_index])
}

# put everything in df
result_df <- data.frame(
  day = 1:365,
  temp = lin_temp$temp,
  means = mean_vec,
  low = lower_vec,
  up = upper_vec
)

# plot results
ggplot(result_df) +
  geom_point(aes(x = day, y = temp))+
  geom_line(aes(x=day,y=means, color = "Mean"))+
  geom_line(aes(x=day,y=low, color = "90 % Probability Intervals"))+
  geom_line(aes(x=day,y=up, color = "90 % Probability Intervals"))+
  labs(title = "Scatterplot with probability intervals", x = "Day of the Year", y = "Temperature")+
  scale_color_manual(name="Definition of lines",
                     values = c("blue", "red"))+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5)) 


# c)

# x-values for maximum
x_tilde <- - plot_df$beta1 / (2* plot_df$beta2) * 365

# structure in data frame
plot_df <- data.frame(
  x = x_tilde
)

# plot
ggplot(plot_df, aes(x=x))+
  geom_histogram(color ="black", fill="gray",bins = 30)+
  labs(title = expression ("Histogram of Beta"[2]), x = expression ("Beta"[2]), y = "Frequency")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))



#2A

#Read in data and get X and y
df <- read.table("WomenAtWork.dat", header=TRUE)

y <- as.matrix(df["Work"])
X <- as.matrix(df[1:nrow(df), 2:ncol(df)])

#Get parameters and number of rows and columns
parameters <- colnames(X)
x_col = ncol(X)
x_row = nrow(X)

#Copied from given example
#log posterior for the logistic
LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  return(logLik + logPrior)
}

# Select the initial values for beta
initVal <- matrix(0,x_col,1)

#Create matrix for Mu with zeros
mu <- as.matrix(rep(0,x_col)) # Prior mean vector

#Calculate Sigma according to the description
tau <- 2
Sigma <-  (tau**2) *diag(x_col)


OptimRes <- optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),
                  control=list(fnscale=-1),hessian=TRUE)

#Posterior Mode
posterior_mode <- OptimRes$par
#print('The posterior mode is:')
#print(posterior_mode)

#Negative Hessian
negative_hessian = -OptimRes$hessian
#print("The Negative Hessian is:")
#print(negative_hessian)

#Inverse negative hessian
inverse_negative_hessian = solve(negative_hessian)
#print("The Inverse Negative Hessian is:")
#print(inverse_negative_hessian)

#Compute an approximate 95% equal tail posterior probability interval for the regression 
#coefficient to the variable NSmallChild.
draws = 10000
beta_values <- rmvnorm(draws, mean = posterior_mode, sigma = inverse_negative_hessian)
beta_n_small_child <- beta_values[1:draws, 6]


#Comparison of estimation results (of the posterior means) to the maximum likelihood estimates.
glmModel <- glm(Work ~ 0 + ., data = df, family = binomial)
#print("Using GLMmodel following coefficients are obtained")
#print(glmModel$coefficients)

#print("With comparison to the drawn beta_values")
#print(glmModel$coefficients-colMeans(beta_values))


sorted_beta_n_small_child <- sort(beta_n_small_child)
lower <- sorted_beta_n_small_child[251]
higher <- sorted_beta_n_small_child[9750]


df_2a <- data.frame(
  x <- 1:draws,
  y <- beta_n_small_child
)

ggplot(df_2a, aes(x=y))+
  geom_histogram(bins=30, aes(y=..density..), color ="black", fill="gray")+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(aes(xintercept = lower, color = "Lower bound"))+
  geom_vline(aes(xintercept = higher, color = "Upper bound"))+
  scale_color_manual(name="Definition of lines",
                     values = c("blue", "red"))+
  labs(title = "Posterior probability interval for NSmallChild with 95% equal tail ", 
       x = "NSmallChild coefficient", 
       y = "Density")+
  theme_light()





#2B
#Draw beta and calculate logistic regression model with own values as x.
simulate_posterior_draws <- function(X, posterior_mode, inverse_negative_hessian){
  betas <- rmvnorm(1, mean = posterior_mode, sigma = inverse_negative_hessian)
  linPred <- X%*%t(betas);
  post_log_reg <- exp(linPred)/(1 + exp(linPred))
  #1-probability to get the y=0, else we got y=1
  return (1-post_log_reg)
}

# x corresponds to a 40-year-old woman, with two children
#(4 and 7 years old), 11 years of education, 7 years of experience, and a husband with an 
#income of 18.
#"Constant"    "HusbandInc"  "EducYears"   "ExpYears"    "Age" "NSmallChild" "NBigChild" 
x_2b <- matrix(data = c(1, 18, 11, 7, 40, 1, 1), ncol=7, nrow=1)


#Draw 10000 samples and store it.
draws = 10000
pos_pred_distribution <- c()

for (n in 1:draws){
  pos_pred_distribution <- c(pos_pred_distribution, 
                             simulate_posterior_draws(x_2b, 
                                                      posterior_mode, 
                                                      inverse_negative_hessian))
}

#Create a dataframe and plot the results
df_2b <- data.frame(
  x <- 1:draws,
  y <- pos_pred_distribution
)

ggplot(df_2b, aes(x=y))+
  geom_histogram(aes(y=..density..), color ="black", fill="gray", bins=30)+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(title = "Posterior predictive distribution of Pr(y = 0|x)", x = "Pr(y = 0|x)", y = "Density")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))




#2C

#Draw beta and calculate logistic regression model with own values as x. In addition predict
#number of women (13) that are not working through binomial distribution
simulate_posterior_draws_multiple <- function(X, posterior_mode, inverse_negative_hessian){
  betas <- rmvnorm(1, mean = posterior_mode, sigma = inverse_negative_hessian)
  linPred <- X%*%t(betas );
  post_log_reg <- exp(linPred)/(1 + exp(linPred))
  #Size = 13 --> Number of women
  prob <- rbinom(n = 1, size = 13, prob = 1-post_log_reg)
  return (prob)
}


#Draw 10000 samples and store it.
draws = 10000
pos_pred_distribution_2c <- c()

for (n in 1:draws){
  pos_pred_distribution_2c <- c(pos_pred_distribution_2c, 
                                simulate_posterior_draws_multiple(x_2b, 
                                                                  posterior_mode, 
                                                                  inverse_negative_hessian))
}


#Create a dataframe and plot the results
df_2c <- data.frame(
  x <- 1:draws,
  y <- pos_pred_distribution_2c
)


ggplot(df_2c, aes(x=y))+
  geom_histogram(aes(y=..density..), color ="black", fill="gray", bins=30)+
  labs(title = "Posterior predictive distribution of Pr(y = 0|x) for 13 women", 
       x = "Number of women", y = "Density")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))