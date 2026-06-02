// make a function for spatial autocorrelation/gaussian process
functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}
data{
     int N_regions;
    array[N_regions] int n_voters;
    array[N_regions] int V;
     vector[N_regions] k;
     vector[N_regions] theta;
     matrix[N_regions,N_regions] D;
}
parameters{
     vector[N_regions] z_b;
     vector[N_regions] z_p;
     real<lower=0> etasq_b;
     real<lower=0> rhosq_b;
     real b_base; //average effect of distance on voting probability
     
     real<lower=0> etasq_p;
     real<lower=0> rhosq_p;
     real p_base_logit; //average probability
     
}
transformed parameters{
     vector[N_regions] b; //total magnitude of distance effect on voting
     vector[N_regions] log_b;
     vector[N_regions] b_place; //area specific effect
     matrix[N_regions,N_regions] L_SIGMAb;
     matrix[N_regions,N_regions] SIGMAb;
     vector[N_regions] p_true; //actual proportion of people who voted
     vector[N_regions] p_max; //hypothetical max if distance=0
     vector[N_regions] p_max_logit;
     vector[N_regions] p_place_logit; //area specific probabilities
     matrix[N_regions,N_regions] L_SIGMAp;
     matrix[N_regions,N_regions] SIGMAp;
     
     //construct covariance matrix for gaussian process
    SIGMAb = cov_GPL2(D, etasq_b, rhosq_b, 0.01); //0.01 is arbitrary, no effect on code
    L_SIGMAb = cholesky_decompose(SIGMAb);
    
    b_place = L_SIGMAb * z_b;
    SIGMAp = cov_GPL2(D, etasq_p, rhosq_p, 0.01);
    L_SIGMAp = cholesky_decompose(SIGMAp);
    
    p_place_logit = L_SIGMAp * z_p;
    
    for ( i in 1:N_regions ) {
        log_b[i] = b_base + b_place[i]; //construct b
        b[i] = exp(log_b[i]);
        p_max_logit[i] = p_base_logit + p_place_logit[i]; //construct p_true and p_max
        p_max[i] = inv_logit(p_max_logit[i]);
        p_true[i] = p_max[i] * (1 + theta[i] * b[i])^(-k[i]);
    }
}

model{

    //priors 
    b_base ~ normal( -3.5 , 0.3 );
    rhosq_b ~ exponential( 0.5 );
    etasq_b ~ exponential( 2 );
    z_b ~ normal( 0 , 1 );
    
    p_base_logit ~ normal( 2.2 , 0.1 );
    rhosq_p ~ exponential( 0.5 );
    etasq_p ~ exponential( 3 );
    z_p ~ normal( 0 , 1 );
    
    
    
    V ~ binomial( n_voters , p_true );
}
