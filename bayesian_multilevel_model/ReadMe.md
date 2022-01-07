# Bayesian Multilevel Model <!-- omit in toc -->

This folder contains Stan Model files and R Scripts for running a Bayesian Multilevel Model estimating fuel consumption for households based on primary cooking fuel and city ward.

- [1. Data Sources](#1-data-sources)
- [2. Multilevel Modelling Overview](#2-multilevel-modelling-overview)
- [3. Model Specification](#3-model_specifications)
- [4. MCMC Sampling with Stan](#4-mcmc-sampling-with-stan)

## 1. Data Sources

> The Bayesian Multilevel Models use a statewide NSS dataset alongside a synthetic population of households for the city of interest with energy and socio-economic attributes.

The NSS Consumer Surveys are carried out regularly by the Ministry of Statistics and Programme Implementation, and microdata is available online (http://microdata.gov.in/nada43/index.php/home).

The synthetic population generated using teh script in this repository draws individual level instances (households) from the Indian Human Development Survey microdata. Information about this survey as well as the data can be accessed online (https://ihds.umd.edu/)

In addition a ward map as a geojson or shapefile will be needed to produce the city ward graph used to calculate the ICAR coefficient based on neighbouring wards. 

## 2. Multilevel Modelling Overview

This aim of this model is to estimate distributions for mean LPG and Firewood consumption and Residual Biomass use for households in a city within a a given state, based on primary cooking fuel choice and city ward.

The procedure carried out by the R script in this folder performs the following steps:

1. Load synthetic population dataset for given city;
2. Load and process NSS Consumer Survey Data for state of interest;
3. Normalise input data;
4. Pass the processed NEED and EPC dataframes to RStan for MCMC sampling;
5. Transform energy consumption distributions back to kWh (reverse normalisation).

## 3. Stan Multilevel Model


### 3.1 Model Versions

There are four versions of the model which incorprate different combinations of the model components.

- Base Model: Household level effects only (basic hierarchical regression with coefficients varying by household typology)
- Random Effects Model: Base model + Random  Effect coefficient varying by city ward.
- ICAR Model: Base model + ICAR (Autoregressive Spatial Coefficient) varying by city ward.

Based on results in our study we found that the Random Effects and BYM (ICAR plus Random Effects Coefficient) versions performed best when tested against real world data.

### 3.2 Model Inputs

Note that all inputs must be in numerical format and should be normalised prior to passing to Stan.

```stan
data {
  int<lower = 0> N; // No. Households
  int<lower = 0> M; // No. Households For Pred
  int<lower = 2> K; // No. Predictors CC
  int<lower = 2> KRB; // No. of Predictors RB
  int<lower = 2> KV; // No. of Predictors 
  int<lower = 0> J; // No. Latent cook choice groups
  int<lower = 0> W; // No. of Wards in City
  int<lower = 1, upper = J> cook_choice[N]; // Primary Cooking Fuel
  int<lower = 1, upper = J> cook_choice_test[M]; // Primary Cooking Fuel
  int<lower = 1, upper = W> ward[M];
  
  int<lower=0> nodes;
  int<lower=0> N_edges;
  int<lower=1, upper=nodes> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nodes> node2[N_edges];  // and node1[i] < node2[i]
  
  vector[N] y; // Vector of outputs (LPG)
  vector[N] z; // Vector of outputs (FW)
  vector[M] y_p; // Vector of previous data (LPG)
  vector[M] z_p; // Vector of previous data (FW)
  vector[N] HH_size; // No. Persons in Household
  vector[N] HH_exp; // Household Expenditure (URP)
  matrix[N,K] x;
  vector[M] HH_size_test; // No. Persons in Household
  vector[M] HH_exp_test; // Household Expenditure (URP)
  real sigma_z_p;
  real sigma_y_p;
  int rb_p[M]; // Vector of historical RB
  int rb[N]; // Vector of NSS RB
  int v_p[M]; // Vector of historical fan ownership
  int v[N]; //
  vector[N] maj_relig;
  vector[N] rcard;
  vector[N] owned;
  vector[N] self_emp;
  vector[N] stsc;
  vector[M] maj_relig_test;
  vector[M] rcard_test;
  vector[M] owned_test;
  vector[M] self_emp_test;
  vector[M] stsc_test;
  matrix[M, K] x_test; // Vector of inputs
}

```

### 3.3 Model Parameters

This model makes use of non-centered parameterisations to improve sampling efficiency, and resolve issues arising from 'funnels' in the parameter-space (described further in the Stan User Manual: https://mc-stan.org/docs/2_19/stan-users-guide/reparameterization-section.html).

In summary we say that the distribution for a coefficient such as `a ~ N(mu_a,sigma_a)` can be represented as `a = mu_a + sigma_a*yay` where `mu_E` is a location parameter, sigma_E is a scale parameter, and yay is a standard normal distribution (`N(0,1)`). This reprarameterisation will allow for more efficient sampling in cases with little data as explained by Betancourt and Girolami (2013).

The coefficients for the model equation (e.g. `a`, `b1`, and `b2`) are defined in the transformed parameter block by a location parameter `mu`, scale parameter `sigma`, and a corresponding standard normal distribition. The later parameters are defined in the parameter block and are actually what Stan is sampling.

In addition we also set a prior on the scale and location hyperparameters using the historical fuel consumption data from the IHDS embedded in the synthetic population. By specifcying a model for the prior data with 'historical' coefficients (`pay`, `peta_1`, `peta_2`, etc.), and relating these 'historical' coefficients to the hyperparameters for the actual model coefficients (i.e. `peta_1 = mu_b1 + sigma_b1 * zee`) we use the historical embedded fuel consumption data to set a prior on the hyperparameters which we can then update with the more recent NSS consumer survey data.

```stan
transformed parameters{
  vector[M] peta_1;
  vector[M] peta_2;
  vector[M] pay;
  
  vector[M] pee_1;
  vector[M] pee_2;
  vector[M] pdee;
  
  vector[M] rb_hat_test;
  vector[N] rb_hat;
    
  vector[J] a;
  vector[J] b1;
  vector[J] b2;
  
  vector[J] d;
  vector[J] c1;
  vector[J] c2;
    
  vector[W] f;
  vector[W] g;
  vector[W] h;
  
  real mu_beta; // Slope 2
  real mu_delta; // Slope 2
  real<lower = 0> sigma_beta;
  real<lower = 0> sigma_delta;
  
  matrix[N, J] x_beta = x * beta;
  
  a = mu_a + sigma_a*yay;
  b1 = mu_b1 + sigma_b1*yee_1;
  b2 = mu_b2 + sigma_b2*yee_2;
  
  d = mu_d + sigma_d*yde;
  c1 = mu_c1 + sigma_c1*yet_1;
  c2 = mu_c2 + sigma_c2*yet_2;
  
  f = mu_ff + sigma_ff*faa;
  g = mu_g + sigma_g*gaa;
  h = mu_h + sigma_h*haa;
  
  for(p in 1:M){
    peta_1[p]=mu_b1+sigma_b1*zee[p];
    peta_2[p]=mu_b2+sigma_b2*zee[p];
    pay[p]=mu_a+sigma_a*zee[p];
    pee_1[p]=mu_c1+sigma_c1*zee[p];
    pee_2[p]=mu_c2+sigma_c2*zee[p];
    pdee[p]=mu_d+sigma_d*zee[p];
    
    rb_hat_test[p] = delta[1] + delta[2] * HH_exp_test[p]
                    + delta[3] * maj_relig_test[p] + delta[4] *rcard_test[p]
                    + delta[5] * owned_test[p] + delta[6] *self_emp_test[p];
 
  }
  
    for(i in 1:N){
    rb_hat[i] = delta[1] + delta[2] * HH_exp[i]
                    + delta[3] * maj_relig[i] + delta[4] *rcard[i]
                    + delta[5] * owned[i] + delta[6] *self_emp[i];
                  
    }
    
    mu_beta = 1*mu_beta_raw;
    sigma_beta = 1*sigma_beta_raw;
    mu_delta = 1*mu_delta_raw;
    sigma_delta = 1*sigma_delta_raw;
}
```

### 3.4 Model

Recall that the actual model for fuel consumption takes the form of `Fuel ~ normal(mu_fuel, sigma_fuel)` where `mu_fuel[p] = a[cook_choice[p]] + b1[cook_choice[p]].*HH_size + b2[cook_choice[p]].*HH_exp + icar[ward[p]] + rand_eff[ward[p]]`.

The first step is to define the hierarchical prior on the coefficient hyperparameters using the historical embedded data in the syunthetic population (our prior). This is also used to estimate the ICAR and Random Effect coefficients, as we do not have any data available that could be used to estimate spatial differnces within the city other than in the synthetic population. One caveat is that the data embeded in the synthetic population is several years out of date and does not include any information about primary cooking fuel choice.

A categorical logit is used to estimate likely primary fuel choice for the synthetic population, and model coefficient for each fuel of interest are estimated using the more recent NSS data.

Using these distributions fuel consumption estimates can be sampled in the generated quantiaties block of the Stan model for each synthetic household based on its estimated primary cooking fuel and the city ward it is located in.

```stan
model {
  
  to_vector(beta) ~ normal(mu_beta, sigma_beta);
  delta ~ normal(mu_delta, sigma_delta);
  mu_beta_raw ~ std_normal();
  sigma_beta_raw ~ std_normal();
  mu_delta_raw ~ std_normal();
  sigma_delta_raw ~ std_normal();

  // hierarchical priors
  zee ~ std_normal();
  faa ~ std_normal();
  gaa ~ std_normal();
  haa ~ std_normal();

  for(p in 1:M){
     z_p[p] ~ normal(pdee[p] + pee_1[p].*HH_size_test[p] + pee_2[p].*HH_exp_test[p]  + phi_f[ward[p]] + f[ward[p]], sigma_z_p);
     y_p[p] ~ normal(pay[p] + peta_1[p].*HH_size_test[p] + peta_2[p].*HH_exp_test[p] + phi_l[ward[p]] + g[ward[p]], sigma_y_p);
     rb_p[p] ~ bernoulli_logit(rb_hat_test[p] + phi_rb[ward[p]] + h[ward[p]]);
  }
  
  phi_f ~ icar_normal_lpdf(nodes, node1, node2);
  phi_l ~ icar_normal_lpdf(nodes, node1, node2);
  phi_rb ~ icar_normal_lpdf(nodes, node1, node2);
  
  // soft sum-to-zero constraint on phi
  // more efficient than mean(phi) ~ normal(0, 0.001)
  sum(phi_f) ~ normal(0, 0.001 * N);
  sum(phi_l) ~ normal(0, 0.001 * N);
  sum(phi_rb) ~ normal(0, 0.001 * N);

  yay~std_normal();
  yee_1~std_normal();
  yee_2~std_normal();
  yde~std_normal();
  yet_1~std_normal();
  yet_2~std_normal();
  
  sigma_f ~ std_normal();
  sigma ~ std_normal();
  
  z ~ normal(d[cook_choice]+ c1[cook_choice].*HH_size + c2[cook_choice].*HH_exp, sigma_f);
  y ~ normal(a[cook_choice] + b1[cook_choice].*HH_size + b2[cook_choice].*HH_exp, sigma);
  
  for(i in 1:N){
    cook_choice[i] ~ categorical_logit((x_beta[i])');
    rb[i] ~ bernoulli_logit(rb_hat[i]);
  }
    
}
```

## 4. MCMC Sampling with Stan

### 4.1 Stan Settings

We recommend following RStan recommended settings:

```r
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
```

In addition the following Stan mcmc sampling settings used are:

```r
iter=4000, 
warmup=1000,
chains=4,
control = list(max_treedepth = 10,adapt_delta = 0.9)
```

### 4.2 Sampling Energy Consumption for Synthetic Households

This is done in the Stan model file using the generated quantitaties block.

### 4.3. Troubleshooting

The Stan documentation provides some excellent discussion and examples concerning common errors and warning messages, however there are a few specific issues encountered when using this. 

1. In some cities there will be virtually no primary users of a particular fuel in the NSS data (this is most commonly the case with kerosene whose sale has been restricted in some states/cities) - this can lead to a problem in sampling as there is no data to use for estimating fuel consumption for that group of households. We found it helpful to add a few instances to the state dataset from the national dataset of households representing the group. This will not significantly impact the fuel consumption estimated for the synthetic population, but will ensure that the sampler has some data to use in estimating distributions for each group.  
2. For large cities (with over 10,000 households) it can take over a day to run one of these models using 16GB RAM and and i7 4-core processor (and may require running the chains in parallel due to memory contraints). For this reason it may be useful to reduce the size of the synthetic population to 1/10th or even 1/20th of actual size - we did not find any significant loss in output accuracy as a result of doing this.

## 5. Notes on Data

- Outputs will need to be transformed back to kWh values as the inputs for the Stan model were all normalised.
