data {
  // number of observations
  int<lower=0> N;
  int<lower=0> I;
  int<lower=0> J;
  int<lower=0> K;
  int Q;
  // response
  vector[N] y;
  // index of genotype
  int genotype[I*J*K];
  // index of environment
  int environment[I*J*K];
  // index of time
  int time[I*J*K];
  // hyperparameters
  real mmu;
  real smu;
  real stheta;
  // priors on sy
  real a;
  real b;
}

parameters {
  real muall;
  vector[I] gen;
  vector[J] env;
  vector[K] t;
  vector[Q] lambda_raw;
  matrix[I,Q] thetaG;
  matrix[J, Q] thetaD;
  matrix[K, Q] thetaP;
  real sg;
  real se;
  real st;
  real slambda;
  real<lower=0> sy;
}

transformed parameters{
  vector[Q] lambda;

  matrix[I,Q] gamma;
  matrix[J, Q] delta;
  matrix[K, Q] phi;

  matrix[I, Q] thetaGNew;
  vector[Q] sqrtThetaG;
  vector[Q] mG;

  matrix[J, Q] thetaDNew;
  vector[Q] sqrtThetaD;
  vector[Q] mD;

  matrix[K, Q] thetaPNew;
  vector[Q] sqrtThetaP;
  vector[Q] mP;

   for(q in 1:Q){
    mG[q] = sum(thetaG[1:I,q])/I;
    for(i in 1:I){
    thetaGNew[i,q] = (thetaG[i,q] - mG[q])^2;
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:I,q]) + 0.000001));
    for(i in 1:I){
      gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q];
    }
  }

  for(q in 1:Q){
    mD[q] = sum(thetaD[1:J,q])/J;
    for(j in 1:J){
    thetaDNew[j,q] = (thetaD[j,q] - mD[q])^2;
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:J,q])));
    for(j in 1:J){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q];
    }
  }

  for(q in 1:Q){
    mP[q] = sum(thetaP[1:K,q])/K;
    for(k in 1:K){
    thetaPNew[k,q] = (thetaP[k,q] - mP[q])^2;
    }
    sqrtThetaP[q] = sqrt(1/(sum(thetaPNew[1:J,q])));
    for(k in 1:K){
      phi[k,q] = thetaPNew[k,q]*sqrtThetaP[q];
    }
  }

  lambda = sort_desc(lambda_raw);

}

model {
  vector[N] mu;
  vector[N] blin;

  for (n in 1:N) {
    blin[n] = 0;
   }

  for (n in 1:N) {
    for (q in 1:Q) {
      blin[n] = blin[n] + (lambda[q] * gamma[genotype[n], q] * delta[environment[n], q] * phi[time[n], q]);
    }
  }

  for (n in 1:N) {
    y[n] ~ normal(mu[n], sy);
    mu[n] = muall + gen[genotype[n]] + env[environment[n]] + t[time[n]] + blin[n];
   }

  // Priors
  // Prior on grand mean
   muall ~ normal(mmu, smu^-2);

  // Prior on genotype effect
  for(i in 1:I) {
  gen[i] ~ normal(0, sg^-2);
  }

  // Prior on environment effect
  for(j in 1:J) {
  env[j] ~ normal(0, se^-2);
  }

  // Prior on time effect
  for(k in 1:K) {
  t[k] ~ normal(0, st^-2);
  }


  // Priors on gamma
  for(q in 1:Q){
    for(i in 1:I){
      thetaG[i,q] ~ normal(0,stheta);
    }
  }

   // Priors on delta
   for(q in 1:Q){
    for(j in 1:J){
      thetaD[j,q] ~ normal(0,stheta);
    }
  }

   // Priors on phi
   for(q in 1:Q){
    for(k in 1:K){
      thetaP[k,q] ~ normal(0,stheta);
    }
  }

  // Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ normal(0, slambda^-2)T[0,];
  }


  // Priors
  sg ~ student_t(0, 1, 1)T[0,];
  se ~ student_t(0, 1, 1)T[0,];
  st ~ student_t(0, 1, 1)T[0,];
  slambda ~ student_t(0, 1, 1)T[0,];


  // Prior on residual standard deviation
   sy ~ gamma(a, b); // inverse of sy
}

