data {
  int<lower=0> N;
  int<lower=0> B1;
  int<lower=0> B2;
  int<lower=0> Q;
  vector[N] Y;
  int<lower=0> var1[N];
  int<lower=0> var2[N];
  real mmu;
  real smu;
  real a;
  real b;
  real mtheta;
  real stheta;
}

parameters {
  real muall;
  vector[B1] b1;
  vector[B2] b2;
  vector[Q] lambda_raw;
  vector<lower=0>[Q] sthetaBeta1;
  vector<lower=0>[Q] sthetaBeta2;
  vector[B1] thetaBeta1[Q];
  vector[B2] thetaBeta2[Q];
  real<lower=0> sy;
}

transformed parameters {
  vector[Q] lambda = sort_asc(lambda_raw);
  vector[N] mu;
  vector[N] inte;
  vector[B1] beta1[Q];
  vector[B2] beta2[Q];
  for (n in 1:N) {
    if (Q == 1) {
      inte[n] = lambda[1] * beta1[var1[n],1] * beta2[var2[n],1];
    } else {
      inte[n] = dot_product(lambda, beta1[var1[n]] .* beta2[var2[n]]);
    }
    mu[n] = muall + b1[var1[n]] + b2[var2[n]] + inte[n];
  }
  for (q in 1:Q){
    {
      vector[B1] thetaBeta1New = thetaBeta1[q] - mean(thetaBeta1[q]);
      beta1[q] = thetaBeta1New / (sqrt(sum(square(thetaBeta1New)))+0.000001);
    }
    {
      vector[B2] thetaBeta2New = thetaBeta2[q] - mean(thetaBeta2[q]);
      beta2[q] = thetaBeta2New / (sqrt(sum(square(thetaBeta2New)))+0.000001);
    }
  }
}

model {
  muall ~ normal(mmu, smu);
  b1 ~ normal(0, 1);
  b2 ~ normal(0, 1);
  for (q in 1:Q){
    thetaBeta1[q] ~ normal(mtheta, sthetaBeta1[q]);
    thetaBeta2[q] ~ normal(mtheta, sthetaBeta2[q]);
    lambda_raw[q] ~ normal(0, 1);
  }
  sy ~ gamma(a, b);
  Y ~ normal(mu, sy);
}
