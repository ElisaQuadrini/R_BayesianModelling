# Radiocarbon Measurement: the use of Bayesian techniques 

## PROJECT_Quadrini pdf file
This project applies Bayesian modeling to radiocarbon (C14) measurements using data from two archaeological sites. A Bayesian Gamma regression with an identity link was implemented to investigate factors influencing measurement error. Posterior inference was performed via Metropolis–Hastings and spike-and-slab variable selection, with diagnostics and predictive checks confirming model reliability.

### Source dataset
*dataset_10sites_clean.xlsx*: radiocarbon data from two archaeological sites (France and Jordan), with 175 observations and four predictors (Age, Material, Method, SiteName). This is a cleaned subset coming from an original dataset, [p3k14c_2022.06](https://www.p3k14c.org/download/)  

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
- **Extra**: Provided a Bayesian Directed Acyclic Graph (DAG) representation of the radiocarbon model, illustrating links between true calendar age, radiocarbon age, lab error, and calibration curve uncertainty.  
