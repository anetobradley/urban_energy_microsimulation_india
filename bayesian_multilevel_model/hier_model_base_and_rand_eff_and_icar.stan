// iheus_neapolitan_hier_BYM.stan
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
/* 
* 
*/
functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]);
 }
}

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
  //matrix[N, KRB] x_rb;
  //matrix[M, KRB] x_rb_test;
  matrix[M, K] x_test; // Vector of inputs
}

parameters {
  matrix[K,J] beta;
  vector[KRB+1] delta;
  vector[KV+1] epsilon;
  vector[J] yay;
  vector[J] yee_1;
  vector[J] yee_2;
  vector[J] yde;
  vector[J] yet_1;
  vector[J] yet_2;
  //vector[J] a;
  //vector[J] b1;
  //vector[J] b2;
  //vector[J] d;
  //vector[J] c1
  vector[W] faa;
  vector[W] gaa;
  vector[W] haa;
  vector[W] laa;
  //vector[W] f;
  //vector[W] g;
  //vector[W] h;
  //vector[J] c2;
  real mu_a; // Slope 1
  real mu_b1; // Slope 1
  real mu_b2; // Slope 2
  real mu_beta_raw; // Slope 2
  real mu_delta_raw; // Slope 2
  real mu_epsilon_raw; // Slope 2
  real mu_d; // Slope 1
  real mu_c1; // Slope 2
  real mu_c2; // Slope 2
  real mu_ff; // Slope 1
  real mu_g; // Slope 1
  real mu_h; // Slope 2
  real mu_l; // Slope 2
  real<lower = 0> sigma_c1;
  real<lower = 0> sigma_c2;
  real<lower = 0> sigma_b1;
  real<lower = 0> sigma_b2;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_d;
  real<lower = 0> sigma_beta_raw;
  real<lower = 0> sigma_delta_raw;
  real<lower = 0> sigma_epsilon_raw;
  real<lower = 0> sigma;
  real<lower = 0> sigma_f;  
  real<lower = 0> sigma_ff;
  real<lower = 0> sigma_g;
  real<lower = 0> sigma_h;
  real<lower = 0> sigma_l;
  vector[M] zee;
  vector[nodes] phi_f; // spatial effects
  vector[nodes] phi_l; // spatial effects
  vector[nodes] phi_rb; // spatial effects
  vector[nodes] phi_app; // spatial effects
}

transformed parameters{
  vector[M] peta_1;
  vector[M] peta_2;
  vector[M] pay;
  
  vector[M] pee_1;
  vector[M] pee_2;
  vector[M] pdee;
  
  vector[M] rb_hat_test;
  vector[N] rb_hat;
  
  vector[M] app_hat_test;
  vector[N] app_hat;
    
  vector[J] a;
  vector[J] b1;
  vector[J] b2;
  
  vector[J] d;
  vector[J] c1;
  vector[J] c2;
    
  vector[W] f;
  vector[W] g;
  vector[W] h;
  vector[W] l;
  
  real mu_beta; // Slope 2
  real mu_delta; // Slope 2
  real mu_epsilon; // Slope 2
  real<lower = 0> sigma_beta;
  real<lower = 0> sigma_delta;
  real<lower = 0> sigma_epsilon;
  
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
  l = mu_l + sigma_l*haa;
  
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
    
    app_hat_test[p] = epsilon[1] + epsilon[2] * HH_exp_test[p]
                    + epsilon[3] * stsc_test[p] + epsilon[4] *HH_size_test[p]
                    + epsilon[5] *self_emp_test[p];
  }
  
    for(i in 1:N){
    rb_hat[i] = delta[1] + delta[2] * HH_exp[i]
                    + delta[3] * maj_relig[i] + delta[4] *rcard[i]
                    + delta[5] * owned[i] + delta[6] *self_emp[i];
                  
    app_hat[i] = epsilon[1] + epsilon[2] * HH_exp[i]
                    + epsilon[3] * stsc[i] + epsilon[4] *HH_size[i]
                    + epsilon[5] *self_emp[i];
    }
    mu_beta = 1*mu_beta_raw;
    sigma_beta = 1*sigma_beta_raw;
    mu_delta = 1*mu_delta_raw;
    sigma_delta = 1*sigma_delta_raw;
    mu_epsilon = 1*mu_epsilon_raw;
    sigma_epsilon = 1*sigma_epsilon_raw;
}

model {
  
  to_vector(beta) ~ normal(mu_beta, sigma_beta);
  delta ~ normal(mu_delta, sigma_delta);
  epsilon ~ normal(mu_epsilon, sigma_epsilon);
  mu_beta_raw ~ std_normal();
  sigma_beta_raw ~ std_normal();
  mu_delta_raw ~ std_normal();
  sigma_delta_raw ~ std_normal();
  mu_epsilon_raw ~ std_normal();
  sigma_epsilon_raw ~ std_normal();
  
  
  // priors
  //mu_a ~ normal(0,5);
  //mu_b1 ~ normal(0,5);
  //mu_b2 ~ normal(0,5);
  //mu_d ~ normal(0,2);
  //mu_c1 ~ normal(0,2);
  //mu_c2 ~ normal(0,5);
  
  // hierarchical priors
  
  zee~std_normal();
  faa ~ std_normal();
  gaa ~ std_normal();
  haa ~ std_normal();
  laa ~ std_normal();
  
  for(p in 1:M){
     z_p[p] ~ normal(pdee[p] + pee_1[p].*HH_size_test[p] + pee_2[p].*HH_exp_test[p]  + phi_f[ward[p]] + f[ward[p]], sigma_z_p);
     y_p[p] ~ normal(pay[p] + peta_1[p].*HH_size_test[p] + peta_2[p].*HH_exp_test[p] + phi_l[ward[p]] + g[ward[p]], sigma_y_p);
     rb_p[p] ~ bernoulli_logit(rb_hat_test[p] + phi_rb[ward[p]] + h[ward[p]]);
     v_p[p] ~ bernoulli_logit(app_hat_test[p] + l[ward[p]] + phi_app[ward[p]]);
  }
  
  phi_f ~ icar_normal_lpdf(nodes, node1, node2);
  phi_l ~ icar_normal_lpdf(nodes, node1, node2);
  phi_rb ~ icar_normal_lpdf(nodes, node1, node2);
  phi_app ~ icar_normal_lpdf(nodes, node1, node2);
  // soft sum-to-zero constraint on phi
  // more efficient than mean(phi) ~ normal(0, 0.001)
  sum(phi_f) ~ normal(0, 0.001 * N);
  sum(phi_l) ~ normal(0, 0.001 * N);
  sum(phi_rb) ~ normal(0, 0.001 * N);
  sum(phi_app) ~ normal(0, 0.001 * N);
  
  //a ~ normal(mu_a, sigma_a);
  //b1 ~ normal(mu_b1, sigma_b1);
  //b2 ~ normal(mu_b2, sigma_b2);
  //d ~ normal(mu_d, sigma_d);
  //c1 ~ normal(mu_c1, sigma_c1);
  //c2 ~ normal(mu_c2, sigma_c2);

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
    v[i] ~ bernoulli_logit(app_hat[i]);
  }
    
}

generated quantities{
  vector[N] y_rhat;
  vector[N] z_rhat;
  int cook_rhat[N];
  
  vector[M] rb_test;
  vector[M] app_test;
  vector[M] y_rhat_test;
  vector[M] z_rhat_test;
  int cook_rhat_test[M];
  
  for(i in 1:N){
    cook_rhat[i] = categorical_logit_rng((x[i]*beta)'); // Posterior predictive
    z_rhat[i] = normal_rng(d[cook_rhat[i]]+ c1[cook_rhat[i]]*HH_size[i] + c2[cook_rhat[i]]*HH_exp[i], sigma_f);//+ c2[cook_rhat[i]]*HH_exp[i], sigma_f);
    y_rhat[i] = normal_rng(a[cook_rhat[i]] + b1[cook_rhat[i]]*HH_size[i] + b2[cook_rhat[i]]*HH_exp[i], sigma);
    }
  
  for(m in 1:M){
    rb_test[m] = inv_logit(rb_hat_test[m]+ phi_rb[ward[m]] + h[ward[m]]);
    app_test[m] = inv_logit(app_hat_test[m] + phi_app[ward[m]]  + l[ward[m]]);
    cook_rhat_test[m] = categorical_logit_rng((x_test[m]*beta)'); // Posterior predictive
    z_rhat_test[m] = normal_rng(d[cook_rhat_test[m]]+ c1[cook_rhat_test[m]]*HH_size_test[m] + c2[cook_rhat_test[m]]*HH_exp_test[m] + phi_f[ward[m]] + f[ward[m]], sigma_f );//+ c2[cook_rhat_test[m]]*HH_exp_test[m], sigma_f);
    y_rhat_test[m] = normal_rng(a[cook_rhat_test[m]] + b1[cook_rhat_test[m]]*HH_size_test[m] + b2[cook_rhat_test[m]]*HH_exp_test[m] + phi_l[ward[m]] + g[ward[m]], sigma);
    }
  
  }
