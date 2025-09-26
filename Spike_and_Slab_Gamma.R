

library(readxl)
library(coda)
library(mvtnorm)
library(R2jags)

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


# ----- MODEL SELECTION VIA SPIKE AND SLAB ------

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
        las = 2)  

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

