# Radiocarbon Measurement: the use of Bayesian techniques 

## PROJECT_Quadrini pdf file
This project applies Bayesian modeling to radiocarbon (C14) measurements using data from two archaeological sites. A Bayesian Gamma regression with an identity link was implemented to investigate factors influencing measurement error. Posterior inference was performed via Metropolis–Hastings and spike-and-slab variable selection, with diagnostics and predictive checks confirming model reliability.

### Source dataset
*dataset_10sites_clean.xlsx*: radiocarbon data from two archaeological sites (France and Jordan), with 175 observations and four predictors (Age, Material, Method, SiteName). This is a cleaned subset coming from an original dataset, [p3k14c_2022.06](https://www.p3k14c.org/download/), a global-scale archaeological radiocarbon database composed of 180,070 radiocarbon dates, coming from samples of sites from all over the world.

### Project Workflow  
- **Exploratory analysis**: Checked the distribution of the response variable (*Error* of radiocarbon age estimation), which was positive and right-skewed, suggesting the use of a Gamma model.  
- **Model specification**: Defined a Bayesian Gamma regression model with identity link, parameterized by shape parameter α and mean µ depending on covariates.  
- **Prior distributions**: Assigned a Normal prior for regression coefficients β (weakly informative) and a Gamma prior for α.  
- **Posterior formulation**: Derived the joint posterior distribution \( p(\beta, \alpha | y) \propto p(y|\beta,\alpha)p(\beta)p(\alpha) \).  
- **Model implementation**:  
  - Used a **Metropolis-Hastings algorithm** with adaptive proposals to approximate the posterior distribution of β and α.  
  - Initialized β with OLS estimates and α = 2, applying thinning to improve convergence.  
- **Variable selection**: Applied a **spike-and-slab prior** (via JAGS) to assess variable inclusion probabilities, discarding less relevant predictors (e.g., “Material: wood”).  
- **Posterior inference**: Obtained concentrated posterior distributions with narrow credible intervals, indicating stable estimates.  
- **Diagnostics**:  
  - Verified convergence through traceplots, autocorrelation, Effective Sample Size (ESS), and Geweke test.  
  - Performed posterior predictive checks (mean, variance, density overlays) to assess model fit against observed data.  
- **Results and interpretation**: Found that predictors such as Age, Material, and Method influence radiocarbon measurement error, improving precision of C14 age estimates. However, the model cannot reduce uncertainty in *calibrated* calendar ages, which depend on external calibration curves (e.g., IntCal20).  
- **Extra**: Provided a Bayesian Directed Acyclic Graph (DAG) representation of the radiocarbon model, illustrating in the reality the links between true calendar age, radiocarbon age, lab error, and calibration curve uncertainty.


### Modello_Gamma_MH file
This is the R script used to perform the Model Implementation step described in the project workflow.  
The design matrix is built from predictors such as Age, Material, Method, and SiteName, with Age standardized. Weakly informative priors are set: a multivariate Normal prior for the regression coefficients β and a Gamma prior for the shape parameter α. Posterior inference is performed using a custom Metropolis–Hastings MCMC algorithm with adaptive proposals for each coefficient and for α, generating 10,000 samples. Convergence is checked via autocorrelation plots and thinning, while posterior inference is summarized through boxplots and histograms of β and α. Finally, posterior predictive checks are performed by simulating new data from the fitted Gamma model and comparing them to the observed distribution using density overlays, histograms, and Monte Carlo evaluation of summary statistics.

### Spike_and_Slab_Gamma file 
The script implements a Bayesian Gamma regression with spike-and-slab priors for automatic variable selection (see the workflow). 
The model is specified in JAGS: each regression coefficient β is linked to an inclusion indicator γ, where γ ~ Bernoulli(w). If γ = 0, the coefficient is excluded (spike); if γ = 1, it comes from a Normal slab prior. Priors are set as before and posterior inference is performed via MCMC with 2 chains and thinning. Posterior inclusion probabilities are computed for each covariate and visualized in a barplot: variables with high posterior inclusion (probability > 0.9) are selected as relevant predictors. A reduced “best model” dataset is then built using the selected variables, with dummy coding for material types (charcoal, bone). The design matrix is recomputed for this final model, preparing the ground for subsequent Bayesian regression analysis with the reduced set of predictors.

### Required R Packages  
`readxl`,`coda`,`mvtnorm`,`R2jags` 

