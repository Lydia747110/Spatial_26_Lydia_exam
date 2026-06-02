functions {

  matrix cov_GPL2(matrix D, real alpha_sq, real rho_sq, real delta) {

    int N = rows(D);
    matrix[N, N] K;

    for (i in 1:(N-1)) {
      K[i, i] = alpha_sq + delta;

      for (j in (i + 1):N) {
        K[i, j] = alpha_sq * exp(-rho_sq * square(D[i,j]));
        K[j, i] = K[i, j];
      }
    }
    K[N,N] = alpha_sq + delta;

    return K;
  }

}

data {

  int<lower=1> N_regions;
  int<lower=2> N_parties;

  array[N_regions, N_parties] int votes;

  matrix[N_regions, N_regions] Dmat;
  vector[N_parties-1] alpha_mu;

}

parameters {

  // baseline logits
  vector[N_parties - 1] alpha;

  // GP hyperparameters
  real<lower=0> rho_sq;
  real<lower=0> eta_sq;

  // latent GP values
  matrix[N_regions, N_parties - 1] f_raw;

}

transformed parameters {

  matrix[N_regions, N_parties] eta;
  array[N_regions] simplex[N_parties] p;

  {
    matrix[N_regions, N_regions] K;
    matrix[N_regions, N_regions] L_K;

    K = cov_GPL2(Dmat, eta_sq, rho_sq, 1e-2);

    L_K = cholesky_decompose(K);

    for (k in 1:(N_parties - 1)) {

      vector[N_regions] f;

      f = L_K * col(f_raw, k);

      eta[,k] = alpha[k] + f;
    }

    // reference category
    eta[,N_parties] = rep_vector(0, N_regions);

    // softmax probabilities
    for (r in 1:N_regions) {
      p[r] = softmax(to_vector(eta[r]));
    }
  }

}

model {

  // priors
  alpha ~ normal(alpha_mu, 0.5);

  rho_sq ~ exponential(0.8);
  eta_sq ~ exponential(2);

  to_vector(f_raw) ~ normal(0, 1);

  // likelihood
  for (r in 1:N_regions) {
    votes[r] ~ multinomial(p[r]);
  }
}
