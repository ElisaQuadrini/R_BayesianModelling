library(readxl)
library(coda)
library(mvtnorm)

data <- read_xlsx("dataset_10sites_clean.xlsx")
mod_data = data[c(1:175), c(2,3,4,5,9)]

mod_data$Material = as.factor(mod_data$Material)
mod_data$Method = as.factor(mod_data$Method)
mod_data$SiteName = as.factor(mod_data$SiteName)

X <- model.matrix(Error ~ ., data = mod_data)# Design matrix

X[, "Age"] <- scale(X[, "Age"])

head(X)
dim(X)

# ----- INPUT -----

# Design matrix X and response y

n <- nrow(X)
p <- ncol(X)

y <- mod_data$Error

# Priors
mu_beta <- rep(0, p)
Sigma_beta <- diag(100, p)
a <- 1
b <- 0.01
delta <- 100

# MCMC settings
S <- 10000
beta_curr <- as.vector(lm(y ~ X - 1)$coefficients)
alpha_curr <- 2

samples <- matrix(NA, nrow = S, ncol = p+1)
colnames(samples) <- c(paste0("beta", 1:p), "alpha")

# Initialization structures for the adaptive proposal
m_j <- rep(0, p)
s2_j <- rep(1, p)  # varianze iniziali di fallback


# ----- LOG POSTERIOR FUNCTION -----

log_posterior <- function(beta, alpha, y, X, mu_beta, Sigma_beta, a, b) {
  mu <- as.vector(X %*% beta)
  if (any(mu <= 0) || alpha <= 0) {
      cat("Invalid proposal: some mu <= 0 or alpha <= 0\n")
      return(-Inf)
    }
  
  loglik <- sum(dgamma(y, shape = alpha, rate = alpha / mu, log = TRUE))
  logprior_beta <- dmvnorm(beta, mean = mu_beta, sigma = Sigma_beta, log = TRUE)
  logprior_alpha <- dgamma(alpha, shape = a, rate = b, log = TRUE)
  
  return(loglik + logprior_beta + logprior_alpha)
}

# ----- MCMC LOOP -----
set.seed(42)

for (s in 1:S) {
  
  # Proposal for each component of beta 
  for (j in 1:p) {
    
    if (s > 1) {
      m_j[j] <- mean(samples[1:(s - 1), j])
      s2_j[j] <- var(samples[1:(s - 1), j])
      if (is.na(s2_j[j]) || s2_j[j] == 0) s2_j[j] <- 1e-6
      beta_j_prop <- rnorm(1, mean = m_j[j], sd = sqrt(s2_j[j]))
    } else {
      beta_j_prop <- rnorm(1, mean = beta_curr[j], sd = sqrt(0.01))
    }
    
    beta_prop <- beta_curr
    beta_prop[j] <- beta_j_prop
    
    log_post_curr <- log_posterior(beta_curr, alpha_curr, y, X, mu_beta, Sigma_beta, a, b)
    log_post_prop <- log_posterior(beta_prop, alpha_curr, y, X, mu_beta, Sigma_beta, a, b)
    
    log_q_ratio <- dnorm(beta_curr[j], mean = m_j[j], sd = sqrt(s2_j[j]), log = TRUE) -
      dnorm(beta_j_prop, mean = m_j[j], sd = sqrt(s2_j[j]), log = TRUE)
    
    log_r <- log_post_prop - log_post_curr + log_q_ratio
    
    if (is.finite(log_r) && log(runif(1)) < log_r) { beta_curr[j] <- beta_j_prop }
  }
  
  # Proposal for alpha
  alpha_prop <- rgamma(1, shape = delta, rate = delta / alpha_curr)
  
  log_post_curr <- log_posterior(beta_curr, alpha_curr, y, X, mu_beta, Sigma_beta, a, b)
  log_post_prop <- log_posterior(beta_curr, alpha_prop, y, X, mu_beta, Sigma_beta, a, b)
  
  log_q_ratio <- dgamma(alpha_curr, shape = delta, rate = delta / alpha_prop, log = TRUE) -
    dgamma(alpha_prop, shape = delta, rate = delta / alpha_curr, log = TRUE)
  
  log_r_alpha <- log_post_prop - log_post_curr + log_q_ratio
  
  if (is.finite(log_r_alpha) && log(runif(1)) < log_r_alpha) {
    alpha_curr <- alpha_prop
  }
  
  # ----- OUTPUT -----
  samples[s, ] <- c(beta_curr, alpha_curr)
}




### CONVERGENCE CHECKS

thin <- 10
samples_thin <- samples[seq(1, nrow(samples), by = thin), ]



par(mfrow = c(2,2))
for(j in 1:4){
  acf(samples_thin[,j])
}


par(mfrow = c(3,2), mar = c(3,3,2,1))
for(j in 5:8){
  acf(samples_thin[,j])
}



# ----- PREDICTIVE CHECKS -----
library(bayesplot)

# Estrai i campioni dalla matrice samples
post.beta <- samples[, 1:p]       # S x p
post.alpha <- samples[, p + 1]    # vettore di lunghezza S

# Boxplot e istogrammi dei coefficienti
par(mfrow = c(1, 1))
boxplot(post.beta, main = "Posterior distribution of beta coefficients")

par(mfrow = c(3, 3))
for (i in 1:p) {
  hist(post.beta[, i], main = paste0("Beta_", i), xlab = "", ylab = "")
}
hist(post.alpha, main = "Posterior alpha", xlab = "alpha")

# Predictive simulation 

n_iter <- nrow(post.beta)
n_obs <- nrow(X)

# Calcolo delle medie predette mu = X * beta per ogni iterazione
mu <- X %*% t(post.beta)  # dimensione: n_obs x n_iter

# Matrice per y_pred simulati
y_pred <- matrix(NA, nrow = n_obs, ncol = n_iter)

# Simulazioni dalla distribuzione Gamma(alpha, rate = alpha / mu)
for (s in 1:n_iter) {
  mu_s <- mu[, s]
  alpha_s <- post.alpha[s]
  rate_s <- alpha_s / mu_s
  y_pred[, s] <- rgamma(n_obs, shape = alpha_s, rate = rate_s)
}

# Trasponi: y_pred deve avere (S x N)
y_pred <- t(y_pred)

# Thinning
y_pred <- y_pred[seq(1, nrow(y_pred), by = 10), ]

# Plot posterior predictive checks 

# Se `y` è un dataframe (es. mod_data$y), adattalo:
ppc_dens_overlay(y, y_pred)
ppc_hist(y, y_pred[1:5, ])



# ----- Monte Carlo evaluation through Posterior Predictive model checking -----

t_rep <- numeric(nrow(y_pred))  # statistica su ogni replica di un sample, es. media

for (s in 1:nrow(y_pred)) {
  t_rep[s] <- mean(y_pred[s, ])
}

# Ora t_rep ~ distribuzione predittiva di t(Y)
par(mfrow = c(1,1))
hist(t_rep, breaks = 50, main = "Posterior Predictive Distribution of t(Y)", xlab = "t(Y)", col = "lightblue")
abline(v = mean(y), col= "red", lwd = 2)




