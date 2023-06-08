library(ggplot2)

set.seed(12345)

#1####
#a)

alpha0 <- 8
beta0 <- 8
s <- 22
f <- 70 - s

new_alpha <- alpha0 + s
new_beta <- beta0 + f

post <- rbeta(10000,new_alpha, new_beta) # shape1 = alpha, shape2 = beta

n_seq <- seq(10, 10000, 10)
means <- c()
sds <- c()

for (n in n_seq)
{
  means <- c(means, mean(post[1:n]))
  sds <- c(sds,sd(post[1:n]))
}

result_df <- data.frame(
  n = n_seq,
  dist_mean = means,
  dist_sd = sds
)

actual_mean <- new_alpha / (new_alpha + new_beta)
actual_sd <- sqrt((new_alpha * new_beta) / 
                    ((new_alpha + new_beta)**2 * (new_alpha + new_beta + 1)))

ggplot(result_df, aes(x = n, y = dist_mean, color = "Sample mean")) +
  geom_line() +
  geom_hline(aes(yintercept = actual_mean, color = "Posterior mean"))+
  scale_color_manual(name="Definition of lines",
                     values = c("red", "black"))+
  labs(title = "Line plot of the mean", x = "Sample size", y = "Mean")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))



ggplot(result_df, aes(x = n, y = dist_sd, color = "Sample standard deviation")) +
  geom_line() +
  geom_hline(aes(yintercept = actual_sd, color = "Posterior standard deviation")) +
  scale_color_manual(name="Definition of lines",
                     values = c("red", "black"))+
  labs(title = "Line plot of the standard deviation", x = "Sample size",
       y = "Standard deviation")+
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))

#b)

sampled_prob <- sum(post>0.3) / length(post)
# shape1 = alpha, shape2 = beta
actual_prob <- pbeta(0.3,new_alpha, new_beta,lower.tail=FALSE)

#c)

sigma <- post / (1-post)

df_1c <- data.frame(
  x = 1:length(sigma),
  y = sigma
)

ggplot(df_1c, aes(x=y))+
  geom_histogram(bins=30, aes(y=..density..), color ="black", fill="gray")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(title = "Posterior distribution of Phi", x = "Phi", y = "Density")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#2####

#a)

samples <- c(33,24,48,32,55,74,23,17)
smean <- mean(samples)
n <- length(samples)
mu <- 3.6
tau2 <- sum((log(samples) - mu)^2) / n

x <- rchisq(10000, n-1)
sigma2 <- (n-1) * tau2 / x

df_2a <- data.frame(
  x = 1:length(sigma2),
  y = sigma2
)

ggplot(df_2a, aes(x=y))+
  geom_histogram(bins=50, aes(y=..density..), color ="black", fill="gray")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(title = "Posterior distribution", x = expression ("Sigma"^2), y = "Density")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#b)

probs <- pnorm(sqrt(sigma2) / sqrt(2),mean=0,sd=1)
G <- 2 * probs - 1

df_2b <- data.frame(
  x = 1:length(G),
  y = G
)

ggplot(df_2b, aes(x=y))+
  geom_histogram(bins=30, aes(y=..density..), color ="black", fill="gray")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(title = "Posterior distribution of Gini coefficients",
       x = "Gini Coefficient", y = "Density")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#c)

sorted_G <- sort(G)
lower <- sorted_G[251]
higher <- sorted_G[9750]

ggplot(df_2b, aes(x=y))+
  geom_density()+
  geom_vline(aes(xintercept = lower, color = "Lower bound"))+
  geom_vline(aes(xintercept = higher, color = "Upper bound"))+
  scale_color_manual(name="Definition of lines",
                     values = c("blue", "red"))+
  labs(title = "Posterior distribution of Gini coefficients with equal-tailed interval",
       x = "Gini Coefficient", y = "Density")+
  theme_light()


# d)

dens_G <- density(G, n = 10000)
# probability mass * n
prob_mass <- sum(dens_G$y)
# get 95% of mass
threshold_mass <- 0.95 * prob_mass

# sort in descending order of y
desc_order <- order(dens_G$y,decreasing = TRUE)
sorted_y <- dens_G$y[desc_order]
sorted_x <- dens_G$x[desc_order]

sum <- 0
i <- 1
# in our case there are always only two x-value boundaries, so this way of
# getting "lower" and "higher" works. If there are more than two x-values, this
# method will not be sufficient
lower <- 1
higher <- 0
# sum up values until threshold_mass is reached
while(sum < threshold_mass){
  sum <- sum + sorted_y[i]
  if(sorted_x[i] < lower){
    lower <- sorted_x[i]
  }
  if(sorted_x[i] > higher){
    higher <- sorted_x[i]
  }
  i <- i + 1
}

# this is the accumulated mass we got
print(sum)
# this is the y value at which we want to place our threshold
print(sorted_y[i])
# these are the x values at that threshold
print(lower)
print(higher)

result_df <- data.frame(
  x <- dens_G$x,
  y <- dens_G$y
)

ggplot(result_df, aes(x = x, y = y, color = "Gini coefficients")) +
  geom_line() +
  geom_hline(aes(yintercept = sorted_y[i], color = "Highest Posterior Density
Interval"))+
  scale_color_manual(name="Definition of lines",
                     values = c("black", "red"))+
  labs(title = "Posterior distribution of Gini coefficients with HPDI",
       x = "Gini coefficient", y = "Density")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#3####
#a)

I0 <- function(kappa){
  return(besselI(kappa, 0))
}

f <- function(kappa, lambda, y, mu, normConst){
  res1 <- exp(kappa * (-lambda + sum(cos(y-mu))))
  res2 <- 1 / (I0(kappa)^n)
  return(res1 * res2 * normConst)
}

obs <- c(-2.79, 2.33, 1.83, -2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.54)

kappa <- seq(0,10,0.01)
lambda <- 0.5
n <- length(obs)
mu <- 2.4

posterior <- f(kappa, lambda, obs, mu, 1)
# approximate normalizing constant
nConst <- sum(posterior * (max(kappa) - min(kappa)) / (length(kappa)))
post_normalized <- posterior/nConst
df_3a <- data.frame(
  x = kappa,
  y = post_normalized
)  

ggplot(df_3a, aes(x = x, y = y)) +
  geom_line(color = "black") +
  labs(title = "Normalized posterior distribution of Kappa",
       x = "Kappa", y = "Posterior value")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))

#b)
index <- which.max(post_normalized)
mode <- index / length(kappa) * (max(kappa) - min(kappa))
