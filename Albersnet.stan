
# Single-run Binomial STAN Model ####


data {
  
  
  // number of observations
  
  int<lower=1> N;
  
  
  // outcome variable
  
  int<lower=0> virus_shared[N];
  
  
  // main effects
  
<<<<<<< HEAD
  real<lower = 0> space[N];
  
  real<lower = 0> phylo[N];
  
  real<lower = 0> space_phylo[N];
=======
>>>>>>> cad80240d1896d4a6ecb76a0d17a74815a1d3f19
  
  // variables related to varying effects
  
  int<lower=1> N_species;
  int<lower=1, upper=N_species> species1[N];
  int<lower=1, upper=N_species> species2[N];
  
  real d_cites_s[N_species];
  int<lower=0> domestic[N_species];
}

///////////////////////////////////////////////////////////////////////////////
  
  
  parameters {
    
    
    // parameters related to main effects
    
    real mu_alpha;
    
<<<<<<< HEAD
    real beta_space; 
    
    real beta_phylo; 
    
    real beta_inter; 
=======
>>>>>>> cad80240d1896d4a6ecb76a0d17a74815a1d3f19
    
    // parameters related to varying effects
    
    vector[N_species] alpha_species;
    
    real beta_d_cites_s;
    
    real beta_domestic;
    
    real<lower=0> sigma;
  }

///////////////////////////////////////////////////////////////////////////////
  
  
  transformed parameters {
    
    
    vector[N_species] varying_effects_predictor;
    
    vector[N] alpha;
    
    
    // calculate varying effects predictor
    
    for (i in 1:N_species) {
      
      varying_effects_predictor[i] = 
        
        beta_d_cites_s*d_cites_s[i] + beta_domestic*domestic[i];
    }
    
    
    // calculate overall linear predictor
    
    for (i in 1:N) {
      
      alpha[i] = 
        
<<<<<<< HEAD
        mu_alpha + 
        beta_space*space[i] + beta_phylo*phylo[i] + beta_inter*space_phylo[i] + 
        alpha_species[species1[i]] + alpha_species[species2[i]];
=======
        mu_alpha + alpha_species[species1[i]] + alpha_species[species2[i]];
>>>>>>> cad80240d1896d4a6ecb76a0d17a74815a1d3f19
    }
  }

///////////////////////////////////////////////////////////////////////////////
  
  
  model {
    
    
    // priors for main effects
    
    mu_alpha ~ normal(0, 5);
    
<<<<<<< HEAD
    beta_space ~ normal(0, 5); 
    
    beta_phylo ~ normal(0, 5); 
    
    beta_inter ~ normal(0, 5); 
=======
>>>>>>> cad80240d1896d4a6ecb76a0d17a74815a1d3f19
    
    
    // priors for varying effects
    
    alpha_species ~ normal(varying_effects_predictor, sigma);
    
    beta_d_cites_s ~ normal(0, 1);
    
    beta_domestic ~ normal(0, 1);
    
    sigma ~ cauchy(0, 1);
    
    
    // likelihood
    
    virus_shared ~ bernoulli_logit(alpha);
  }