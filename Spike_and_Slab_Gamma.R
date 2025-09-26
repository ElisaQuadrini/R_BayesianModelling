

library(readxl)
library(coda)
library(mvtnorm)

data <- read_xlsx("dataset_10sites_clean.xlsx")
mod_data = data[c(1:175), c(2,3,4,5,9)]

mod_data$Material = as.factor(mod_data$Material)
# mod_data$Method <- factor(mod_data$Method, levels = c("Radiometric", "AMS","Unknown"))
mod_data$Method = as.factor(mod_data$Method)
mod_data$SiteName = as.factor(mod_data$SiteName)

X <- model.matrix(Error ~ ., data = mod_data)# Design matrix

X[, "Age"] <- scale(X[, "Age"])

head(X)
dim(X)

p <- ncol(X)
n <- nrow(mod_data)

y <- mod_data$Error





# ----- MODEL SELECTION VIA SPIKE AND SLAB -----

# --------- INPUT: dati e dimensioni -------------------
n <- nrow(X)
p <- ncol(X)
S <- 10000

# --------- PRIORS -------------------------------------
beta0      <- rep(0, p)        # prior mean for beta
sigma2_0   <- rep(10, p)      # prior variance for beta (slab)
w          <- 0.5              # prior inclusion prob for gamma_j
a <- 1                         # shape for alpha ~ Gamma(a,b)
b <- 0.01                      # rate

# --------- INITIAL VALUES -----------------------------
beta_curr  <- rep(0, p)
gamma_curr <- rep(0, p)        # all predictors initially included
alpha_curr <- 2

# --------- MCMC STORAGE -------------------------------
samples_beta  <- matrix(NA, nrow = S, ncol = p)
samples_gamma <- matrix(NA, nrow = S, ncol = p)
samples_alpha <- numeric(S)

# --------- ADAPTIVE PROPOSAL --------------------------
m_j  <- rep(0, p)
s2_j <- rep(1, p)

# --------- LOG-POSTERIOR FUNCTION ---------------------
log_posterior <- function(beta, gamma, alpha, y, X,
                          beta0, sigma2_0, w, a, b) {
  
  # Imposta beta coerente con gamma (0 se escluso)
  beta_eff <- beta
  beta_eff[gamma == 0] <- 0
  
  mu <- as.vector(X %*% beta_eff)
  if (any(mu <= 0) || alpha <= 0) return(-Inf)
  
  loglik <- sum(dgamma(y, shape = alpha, rate = alpha / mu, log = TRUE))
  
  logprior_beta <- sum(sapply(1:length(beta), function(j) {
    if (gamma[j] == 1) {
      dnorm(beta[j], mean = beta0[j], sd = sqrt(sigma2_0[j]), log = TRUE)
    } else {
      0
    }
  }))
  
  logprior_gamma <- sum(dbinom(gamma, size = 1, prob = w, log = TRUE))
  logprior_alpha <- dgamma(alpha, shape = a, rate = b, log = TRUE)
  
  return(loglik + logprior_beta + logprior_gamma + logprior_alpha)
}




# --------- MCMC LOOP ----------------------------------
set.seed(42)
for (s in 1:S) {
  
  # ---- Step 1: Update gamma_j via Gibbs --------------
  for (j in 1:p) {
    gamma_prop <- gamma_curr
    gamma_prop[j] <- 1 - gamma_curr[j]  # flip inclusion
    
    beta_prop <- beta_curr
    # Se proponi gamma_j = 0 ⇒ beta_j = 0
    if (gamma_prop[j] == 0) beta_prop[j] <- 0
    
    log_r <- log_posterior(beta_prop, gamma_prop, alpha_curr, y, X,
                           beta0, sigma2_0, w, a, b) -
      log_posterior(beta_curr, gamma_curr, alpha_curr, y, X,
                    beta0, sigma2_0, w, a, b)
    
    if (is.finite(log_r) && log(runif(1)) < log_r) {
      gamma_curr <- gamma_prop
      beta_curr <- beta_prop
    }
  }
  
  
  # --------- STEP 2: beta_j via MH ---------------
  for (j in 1:p) {
    
    
    if (gamma_curr[j] == 1) {
      
      beta_prop <- beta_curr
      
      # Adaptive proposal
      if (s > 1) {
        m_j[j] <- mean(samples_beta[1:(s - 1), j])
        s2_j[j] <- var(samples_beta[1:(s - 1), j])
        if (is.na(s2_j[j]) || s2_j[j] == 0) s2_j[j] <- 1e-6
        sd_prop <- sqrt(s2_j[j])
      } else {
        sd_prop <- 0.1
      }
      
      beta_prop[j] <- rnorm(1, mean = beta_curr[j], sd = sd_prop)
      
      
      # Calcola log-posterior
      log_post_prop <- log_posterior(beta_prop, gamma_curr, alpha_curr, y, X,
                                     beta0, sigma2_0, w, a, b)
      log_post_curr <- log_posterior(beta_curr, gamma_curr, alpha_curr, y, X,
                                     beta0, sigma2_0, w, a, b)
      log_r <- log_post_prop - log_post_curr
      
      
      # MH accettazione
      if (is.finite(log_r) && log(runif(1)) < log_r) {
        beta_curr[j] <- beta_prop[j]
      } 
  }
  
  
  
  # ---- Step 3: Update alpha via MH -------------------
  log_alpha_curr <- log(alpha_curr)
  log_alpha_prop <- rnorm(1, mean = log_alpha_curr, sd = 0.1)
  alpha_prop     <- exp(log_alpha_prop)
  
  log_r_alpha <- log_posterior(beta_curr, gamma_curr, alpha_prop, y, X,
                               beta0, sigma2_0, w, a, b) -
    log_posterior(beta_curr, gamma_curr, alpha_curr, y, X,
                  beta0, sigma2_0, w, a, b) +
    log_alpha_prop - log_alpha_curr
  
  if (is.finite(log_r_alpha) && log(runif(1)) < log_r_alpha) {
    alpha_curr <- alpha_prop
  }
  
  # ---- Save samples ----------------------------------
  samples_beta[s, ]  <- beta_curr
  samples_gamma[s, ] <- gamma_curr
  samples_alpha[s]   <- alpha_curr
  }
}

# --------- Posterior Inclusion Probabilities ----------
post_inclusion_prob <- colMeans(samples_gamma)
names(post_inclusion_prob) <- paste0("X", 1:p)

# --------- Plot ---------------------------------------
barplot(post_inclusion_prob,
        col = "steelblue", ylim = c(0, 1),
        ylab = expression(hat(p)[j]),
        main = "Posterior Inclusion Probabilities",
        las = 2)

# --------- Selected Variables -------------------------
selected <- which(post_inclusion_prob > 0.9)
cat("Variabili selezionate:", paste0("X", selected), "\n")








library(R2jags)

jags_data <- list(
  y = y,
  X = X,
  n = n,
  p = p,
  w = 0.5,         # inclusion prior
  sigma2 = 10,     # slab variance
  a0 = 1,          # prior shape for alpha
  b0 = 1           # prior rate for alpha
)

model_string <- "
model {
  for (i in 1:n) {
    mu_raw[i] <- inprod(X[i, ], beta[])
    mu[i] <- max(mu_raw[i], 1e-3)   
    y[i] ~ dgamma(alpha, alpha / mu[i])
  }

  for (j in 1:p) {
    gamma[j] ~ dbern(w)
    theta[j] ~ dnorm(0, 1 / sigma2)
    beta[j] <- gamma[j] * theta[j]
  }

  alpha ~ dgamma(a0, b0)
}
"

params <- c("beta", "gamma", "alpha")

inits <- function() {
  flag <- TRUE
  
  while (flag) {
    gamma <- rbinom(p, 1, 0.5)
    theta <- rnorm(p, 0, 1)
    beta <- gamma * theta
    mu <- as.vector(X %*% beta)
    
    flag <- all(mu > 0)
  }
  
  return(list(gamma = gamma, theta = theta, alpha = 1))
}


model_out <- jags(
  data = jags_data,
  inits = inits,
  parameters.to.save = params,
  model.file = textConnection(model_string),
  n.chains = 2,
  n.iter = 5000,
  n.burnin = 1000,
  n.thin = 2
)

print(model_out)


post_gamma <- as.matrix(model_out$BUGSoutput$sims.list$gamma)
inclusion_prob <- colMeans(post_gamma)
inclusion_prob

col <- colnames(X)

barplot(inclusion_prob,
        col = "steelblue",
        ylim = c(0, 1),
        ylab = expression(hat(p)[j]),
        main = "Posterior Inclusion Probabilities",
        names.arg = col,
        las = 2)  # ruota le etichette per maggiore leggibilità

# --------- Selected Variables -------------------------
selected <- which(inclusion_prob > 0.9)
cat("Selected variables:", paste0(" ", col(selected)), "\n")





# ----- Best model ------

mod_best <- mod_data[, c(1,2,4,5)]
mod_best$Charcoal <- ifelse(mod_data$Material == "charcoal", 1, 0)
mod_best$Bone <- ifelse(mod_data$Material == "bone", 1, 0)


mod_best$Method = as.factor(mod_best$Method)
mod_best$SiteName = as.factor(mod_best$SiteName)

X <- model.matrix(Error ~ ., data = mod_best)# Design matrix

X[, "Age"] <- scale(X[, "Age"])

head(X)
dim(X)

p <- ncol(X)
n <- nrow(mod_data)

y <- mod_data$Error


#### Con metropolis 


# ----- INPUT -----

# Design matrix X and response y
# Make sure to define or load them here
# X <- ...
# y <- ...

n <- nrow(X)
p <- ncol(X)

# Priors
mu_beta <- rep(0, p)
Sigma_beta <- diag(100, p)
a <- 1
b <- 0.01
delta <- 40

# sigma_log_alpha <- 0.1         # std dev on log-alpha proposal

# MCMC settings
S <- 10000
beta_curr <- as.vector(lm(y ~ X - 1)$coefficients)
alpha_curr <- 2

samples <- matrix(NA, nrow = S, ncol = p+1)
colnames(samples) <- c(paste0("beta", 1:p), "alpha")
accepted <- 0

# Inizializza strutture per proposal adattiva
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
  
  beta_prop <- beta_curr
  
  # Adattamento: media e varianza fino all’iterazione precedente
  if (s > 1) {
    for (j in 1:p) {
      m_j[j] <- mean(samples[1:(s - 1), j])
      s2_j[j] <- var(samples[1:(s - 1), j])
      if (is.na(s2_j[j]) || s2_j[j] == 0) s2_j[j] <- 1e-6
      beta_prop[j] <- rnorm(1, mean = m_j[j], sd = sqrt(s2_j[j]))
    }
  } else {
    beta_prop <- as.vector(rmvnorm(1, mean = beta_curr, sigma = diag(0.01, p)))
  }
  
  ### --- Proposal per alpha: Gamma centrata in alpha_curr ---
  alpha_prop <- rgamma(1, shape = delta, rate = delta / alpha_curr)
  
  # Calcolo log-posterior
  log_post_curr <- log_posterior(beta_curr, alpha_curr, y, X, mu_beta, Sigma_beta, a, b)
  log_post_prop <- log_posterior(beta_prop, alpha_prop, y, X, mu_beta, Sigma_beta, a, b)
  
  ### --- Log proposal ratio per MH ---
  # q(curr | prop) - q(prop | curr)
  log_q_ratio <- dgamma(alpha_curr, shape = delta, rate = delta / alpha_prop, log = TRUE) -
    dgamma(alpha_prop, shape = delta, rate = delta / alpha_curr, log = TRUE)
  
  # MH ratio
  log_r <- log_post_prop - log_post_curr + log_q_ratio
  
  # Accept/reject
  if (is.finite(log_r) && log(runif(1)) < log_r) {
    beta_curr <- beta_prop
    alpha_curr <- alpha_prop
    accepted <- accepted + 1
  }
  
  samples[s, ] <- c(beta_curr, alpha_curr)
}


# ----- OUTPUT -----
cat("Acceptance rate:", accepted / S, "\n")
summary <- apply(samples, 2, mean)
print(summary)



# ----- OUTPUT -----
cat("Acceptance rate:", accepted / S, "\n")
summary <- apply(samples, 2, mean)
print(summary)




### CONVERGENCE CHECKS
par(mfrow = c(3,2), mar = c(3,3,2,1), bty = "l")
for(j in 1:4){
  plot(samples[1:100,j], type = "l", main = paste0("Beta_",j-1),
       ylab = "")
} 
par(mfrow = c(3,2), mar = c(3,3,2,1), bty = "l")
for(j in 5:7){
  if (j == 7)
  {
    plot(samples[1:100,j], type = "l", main = "alpha",
         ylab = "")
  }
  else {
    plot(samples[1:100,j], type = "l", main = paste0("Beta_",j-1),
         ylab = "")
  }
  
}

par(mfrow = c(2,2))
for(j in 1:4){
  acf(samples[,j])
}


par(mfrow = c(3,2), mar = c(3,3,2,1))
for(j in 5:7){
  acf(samples[,j])
}


ess = apply(samples[,1:7], 2, effectiveSize) # Effective sample size for each beta
ess

thin <- 10
samples_thin <- samples[seq(1, nrow(samples), by = thin), ]

ess = apply(samples_thin[,1:7], 2, effectiveSize) # Effective sample size for each beta
ess


par(mfrow = c(2,2))
for(j in 1:4){
  acf(samples_thin[,j])
}


par(mfrow = c(3,2), mar = c(3,3,2,1))
for(j in 5:7){
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

# mean 
par(mfrow = c(1,1))
hist(t_rep, breaks = 50, main = "Posterior Predictive Distribution of the Mean", xlab = "t(Y)", col = "lightblue", freq= FALSE)
abline(v = mean(y), col= "red", lwd = 2)

t2_rep <- numeric(nrow(y_pred))
for (s in 1:nrow(y_pred)) {
  t2_rep[s] <- var(y_pred[s,])
}
# variance
par(mfrow = c(1,1))
hist(t2_rep, breaks = 50, main = "Posterior Predictive Distribution of the Variance", xlab = "t(Y)", col = "lightblue", freq= FALSE)
abline(v = var(y), col= "red", lwd = 2)








