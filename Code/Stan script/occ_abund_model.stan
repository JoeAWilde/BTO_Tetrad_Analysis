data {
  int<lower=0> N;
  
  int<lower=1> N_tets;
  array[N] int tets;

  int<lower=1> N_counties;
  array[N] int county;
  
  int<lower=0> X;
  matrix[N, X] preds;
  
  array[N] int occ;
  array[N] real<lower=0, upper=1> abund;

  matrix[N_tets, N_tets] D; // Euclidean distance matrix between all tetrads
}

parameters {

  // --- GP Parameters ---
  real<lower=0> sigma_gp; // Marginal standard deviation
  real<lower=0> rho;      // Length scale (decay rate)
  vector[N_tets] eta_spatial_std; // Standard normal deviates


  // --- Fixed Effects ---
  real a_occ;
  real a_abund;
  vector[X] beta;
  real p_det_eta;
  
  // --- Hierarchical Precision Parameter ---
  // Allowing precision to vary by county addresses heteroscedasticity
  vector<lower=0>[N_counties] phi_county; 
}

transformed parameters {
  // --- 1. Compute Spatial Random Effects (eta_spatial) ---
  // Generates the spatially correlated tetrad-level random effects
  matrix[N_tets, N_tets] Sigma;
  
  // Cholesky decomposition for stable MVN sampling
  matrix[N_tets, N_tets] L; 
  vector[N_tets] eta_spatial;

  for (i in 1:(N_tets - 1)) {
    Sigma[i, i] = square(sigma_gp) + 1e-9; 
    for (j in (i + 1):N_tets) {
      // Squared Exponential (SE) kernel: sigma_gp^2 * exp(-D_ij^2 / (2 * rho^2))
      Sigma[i, j] = square(sigma_gp) * exp(-square(D[i, j]) / (2.0 * square(rho)));
      Sigma[j, i] = Sigma[i, j];
    }
  }
  Sigma[N_tets, N_tets] = square(sigma_gp) + 1e-9;
  
  L = cholesky_decompose(Sigma);
  eta_spatial = L * eta_spatial_std; // The shared spatial effect
}

model {
  // --- Priors ---
  
  // GP Priors (Semi-informative)
  rho ~ inv_gamma(5, 5);
  sigma_gp ~ student_t(3, 0, 1);
  eta_spatial_std ~ normal(0, 1);
  
  // Fixed Effects Priors
  a_occ ~ normal(0, 1.5);
  a_abund ~ normal(0, 1.5);
  beta ~ normal(0, 1);
  p_det_eta ~ normal(0, 1);
  
  // Hierarchical Precision Prior: Use an exponential for strictly positive phi
  phi_county ~ exponential(0.1); 

  // --- Likelihood (The Joint ZAB Model) ---
  
  real p_det;      
  real p_occ;      
  real abund_eta;  
  real abund_mu;   
  int c_idx; // County index for current observation

  p_det = inv_logit(p_det_eta);
  
  for (n in 1:N) {
    c_idx = county[n];
    
    // Occupancy (psi): Depends on GP and Occupancy Intercept
    // Logit(psi) = a_occ + eta_spatial[tets[n]]
    p_occ = inv_logit(a_occ + eta_spatial[tets[n]]);
    
    // Abundance Mean (theta): Depends on GP, Abundance Intercept, and Covariates
    // Logit(theta) = a_abund + eta_spatial[tets[n]] + X * beta
    abund_eta = a_abund + eta_spatial[tets[n]] + dot_product(preds[n], beta);
    abund_mu = inv_logit(abund_eta);
    
    // Likelihood
    if(occ[n] == 1) { // If recorded as occupied (abund > 0)
      // L(y>0) = psi * p * f(abund | theta, phi_county)
      target += log(p_occ) + log(p_det) + 
                beta_lpdf(abund[n] | abund_mu * phi_county[c_idx], 
                                     (1.0 - abund_mu) * phi_county[c_idx]);
    } else { // If not recorded as occupied (abund = 0)
      // L(y=0) = (psi * (1-p)) + (1-psi) -> Must use log_sum_exp
      target += log_sum_exp(
        log(p_occ) + log(1.0 - p_det), // Occupied but not detected
        log(1.0 - p_occ)               // Not occupied
      );
    }
  }
}

generated quantities {
  array[N_counties] real phi_county_inv;
  
  // Calculate the inverse-logit of the detection parameter
  real p_det_median = inv_logit(p_det_eta);

  // You can now extract the precision for each county
  for(k in 1:N_counties) {
    phi_county_inv[k] = 1.0 / phi_county[k]; // A common way to report variance proxy
  }
  vector[N] occ_sim;
  
  real p_occ_sim;
  real abund_eta_sim;
  real abund_mu_sim;
  real p_det_sim;

  for(n in 1:N) {
    p_det_sim = inv_logit(p_det_eta);
    
    p_occ_sim = inv_logit(a_occ + eta_spatial[tets[n]]);
    
    abund_eta_sim = a_abund + eta_spatial[tets[n]] + dot_product(preds[n], beta);
    abund_mu_sim = inv_logit(abund_eta_sim);
    
    if(bernoulli_rng(p_occ_sim) == 1) {
      if(bernoulli_rng(p_det_sim) == 1) {
        occ_sim[n] = beta_rng(abund_mu_sim * phi_county[county[n]], (1 - abund_mu_sim) * phi_county[county[n]]);
      } else {
        occ_sim[n] = 0;
      }
    } else {
      occ_sim[n] = 0;
    }
  }
}

