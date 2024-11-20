data {
  int<lower=0> N;
  
  int<lower=1> N_tets;
  array[N] int tets;
  
  int<lower=0> X;
  matrix[N, X] preds;
  
  array[N] int occ;
  array[N] real<lower=0, upper=1> abund;
}

parameters {
  vector[N_tets] a_occ_tets_std;
  
  real a_occ_mu;
  real<lower=0> sd_occ_tets;
  
  vector[N_tets] a_abund_tets_std;
  
  real a_abund_mu;
  real<lower=0> sd_abund_tets;
  
  vector[X] beta;
  
  real p_det_eta;
  
  real<lower=0> phi;
}

transformed parameters {
  vector[N_tets] a_occ_tets = a_occ_mu + (sd_occ_tets * a_occ_tets_std);
  vector[N_tets] a_abund_tets = a_abund_mu + (sd_abund_tets * a_abund_tets_std);
}
model {
  a_occ_tets_std ~ normal(0, 1);
  a_abund_tets_std ~ normal(0, 1);
  
  a_occ_mu ~ normal(0, 1);
  a_abund_mu ~ normal(0, 1);
  
  sd_occ_tets ~ exponential(1);
  sd_abund_tets ~ exponential(1);
  
  beta ~ normal(0, 1);
  
  p_det_eta ~ normal(0.5, 1);
  
  phi ~ exponential(1);
  
  //Create temporary variables for likelihood calcs.
  real p_occ;
  
  real abund_eta;
  real abund_mu;
  
  real p_det;
  
  for (n in 1:N) {
    p_det = inv_logit(p_det_eta);
    
    p_occ = inv_logit(a_occ_tets[tets[n]]);
    abund_eta = a_abund_tets[tets[n]];
    
    for(x in 1:X) {
      abund_eta += beta[x] * preds[n, x];
    }
    abund_mu = inv_logit(abund_eta);
    
    if(occ[n] == 1) { //if recorded as occupied
      target += log(exp(bernoulli_lpmf(1 | p_occ)) * exp(bernoulli_lpmf(1 | p_det)) * exp(beta_lpdf(abund[n] | abund_mu * phi, (1 - abund_mu) * phi ))); // Occupied, detected
    } else { //if not recorded as occupied
      target += log(exp(bernoulli_lpmf(1 | p_occ)) * exp(bernoulli_lpmf(0 | p_det))); // Occupied but not detected
      
      target += bernoulli_lpmf(0 | p_occ); //Not occupied
    }
  }
}

generated quantities {
  vector[N] occ_sim;
  
  real p_occ_sim;
  real abund_eta_sim;
  real abund_mu_sim;
  real p_det_sim;

  for(n in 1:N) {
    p_det_sim = inv_logit(p_det_eta);
    
    p_occ_sim = inv_logit(a_occ_tets[tets[n]]);
    
    abund_eta_sim = a_abund_tets[tets[n]];
    for(x in 1:X) {
      abund_eta_sim += beta[x] * preds[n, x];
    }
    abund_mu_sim = inv_logit(abund_eta_sim);
    
    if(bernoulli_rng(p_occ_sim) == 1) {
      if(bernoulli_rng(p_det_sim) == 1) {
        occ_sim[n] = beta_rng(abund_mu_sim * phi, (1 - abund_mu_sim) * phi);
      } else {
        occ_sim[n] = 0;
      }
    } else {
      occ_sim[n] = 0;
    }
  }
}

