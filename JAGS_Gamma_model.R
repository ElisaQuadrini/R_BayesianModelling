## JAGS model
library(R2jags)
library(readxl)
library(coda)

data <- read_xlsx("dataset_10sites_clean.xlsx")
mod_data = data[c(1:175), c(2,3,4,5,9)]

mod_data$Material = as.factor(mod_data$Material)
mod_data$Method = as.factor(mod_data$Method)
mod_data$SiteName = as.factor(mod_data$SiteName)

X <- model.matrix(Error ~ . - 1, data = mod_data) # Design matrix
X <- cbind(1, scale(X)) # Standardize the design matrix, add the intercept
head(X)
dim(X)

sig2_beta = 100
a0 = 1
b0 = 1
# Step 1: Define the model 
model_code = function() {
  # Priors
  for(j in 1:p){
    beta[j] ~ dnorm(beta0, prec_beta)
  }
  a ~ dgamma(a0, b0)

  # Likelihood
  # for(i in 1:n){
  #  log(mu[i]) <- inprod(X[i, ], beta)  # calcola mu[i] come funzione di X[i,] e beta
  #   y[i] ~ dgamma(a, a / exp(mu[i]))
  # }
  
  
  # Likelihood
  for(i in 1:n){
    mu[i] <- inprod(X[i, ], beta)    # SENZA log
    y[i] ~ dgamma(a, a / mu[i])
  }
  
}

y = mod_data$Error

# Step 2: specify the data  
jags_data = list(y = c(y), X = X, n = nrow(X), p = ncol(X), a0 = 1, b0=1,
                 beta0 = 0, prec_beta = 1/sig2_beta)

# Step 3: parameters definition
model1_parameters = c("beta","a")

# Step 4: initial values (not mandatory)
p = ncol(X)
init_values = function(){list(beta = rep(0, p))}
Niter <- 30000 # Number of retained samples
burn_in <- 30000 # Burn-in period

# Step 4: run the model 
out_model1 = jags(data = jags_data , inits = init_values, 
                  parameters.to.save = model1_parameters,
                  model.file = model_code, 
                  n.chains = 1, n.iter = Niter+burn_in, n.burnin = burn_in,
                  n.thin = 10)



out_model1$BUGSoutput
head(out_model1$BUGSoutput$sims.array) # sims.array is a 3D array with iterations, number of chains and parameter to be estimated
dim(out_model1$BUGSoutput$sims.array)

# Diagnostic
dim(out_model1$BUGSoutput$sims.matrix) # sims.matrix is a matrix with as rows the iterations*numberofchains, as columns the parameters
head(out_model1$BUGSoutput$sims.matrix)

sim = out_model1$BUGSoutput$sims.matrix
par(mfrow = c(3,2), mar = c(3,3,2,1), bty = "l")
for(j in 1:4){
  if (j == 1) {
    plot(sim[1:100,j], type = "l", main = paste0("a"),
         ylab = "")
  }
  plot(sim[1:100,j], type = "l", main = paste0("Beta_",j-1),
       ylab = "")
}
par(mfrow = c(3,2), mar = c(3,3,2,1), bty = "l")
for(j in 5:10){
  plot(sim[1:100,j], type = "l", main = paste0("Beta_",j-1),
       ylab = "")
}

par(mfrow = c(2,2))
for(j in 1:4){
  acf(sim[,j])
}


par(mfrow = c(3,2), mar = c(3,3,2,1))
for(j in 5:10){
  acf(sim[,j])
}


ess = apply(sim[,1:10], 2, effectiveSize) # Effective sample size for each beta
ess

# Acceptance rate: non posso farla perche non sto usando un metropolis hastings ma solo un gibbs sampling  
1 - rejectionRate(as.mcmc(sim[,1:7]))
rejectionRate(as.mcmc(sim[,1:7]))


######################################
library(R2WinBUGS)
library(rjags)

model_code <- "
model {
  for(j in 1:p){
    beta[j] ~ dnorm(beta0, prec_beta)
  }
  a ~ dgamma(a0, b0)
  
  for(i in 1:n){
    log(mu[i]) <- inprod(X[i, ], beta)
    y[i] ~ dgamma(a, a / exp(mu[i]))
  }
}
"

# Passa la stringa al modello
mod <- jags.model(textConnection(model_code),
                  data = jags_data,
                  inits = init_values,
                  n.chains = 1)

# Mostra i sampler usati
list.samplers(mod)

# Esegui il burn-in
update(mod, 1000)

# Simula i parametri
samples <- coda.samples(mod, 
                        variable.names = model1_parameters, 
                        n.iter = 10000)
samples
###########################################

# Predictive checks --------------------------------------------------------

post.a = out_model1$BUGSoutput$sims.matrix[, 1]
post.beta = out_model1$BUGSoutput$sims.matrix[, 2:10]

par(mfrow = c(1, 1))
boxplot(post.beta, main = "Posterior distribution of beta coefficients")

for(i in 1:9) {
  hist(post.beta[,i])
}

hist(post.a)

# simulate data from the posterior predictive 

# post.beta: matrice (n_iterazioni x p)
# post.a: vettore lungo n_iterazioni
# X: matrice (n_osservazioni x p)

n_iter <- nrow(post.beta)
n_obs <- nrow(X)

# Calcola mu per ogni iterazione e osservazione
eta = X %*% t(post.beta)   # (n_obs x n_iter)
dim(eta)
mu = exp(eta)              # media gamma per ogni iterazione e osservazione

# Matrice per salvare le simulazioni
y_pred = matrix(NA, nrow = n_obs, ncol = n_iter)

# Simula i dati
for(s in 1:n_iter){
  y_pred[, s] = rgamma(n_obs, shape = post.a[s], rate = post.a[s] / mu[, s])
}


dim(y_pred)
y_pred = t(y_pred) # it must be a matrix with S rows and N columns 

y_pred = y_pred[seq(1,nrow(y_pred),by = 10),] # thinning
dim(y_pred)


# A cool library for some plots 
library(bayesplot)
# ypred must be S by N matrix of draws from the posterior (or prior) predictive distribution. 
# The number of rows, S, is the size of the posterior (or prior) sample used to generate 
# ypred.
# The number of columns, N is the number of predicted observations (length(y)).

# Confronto tra dati osservati y e prime 50 simulazioni predittive
ppc_dens_overlay(y,y_pred[1:50,])
ppc_hist(y, y_pred[1:5, ]) 

ones_pred = apply(y_pred,1,function(i) length(which(i == 1)))
summary(ones_pred)
obs = length(which(y == 1))
obs
ppc_stat(y, ypred, stat = function(i) length(which(i == 1))) + theme(aspect.ratio = 1)









