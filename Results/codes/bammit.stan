data {
  // number of observations
  int<lower=0> N;
  int<lower=0> B1;
  int<lower=0> B2;
  int Q;
  // response
  vector[N] y;
  // index of genotype
  int var1[N];
  // index of environment
  int var2[N];
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
  vector[B1] gen;
  vector[B2] env;
  vector[Q] lambda_raw;
  matrix[B1,Q] thetaB1;
  matrix[B2, Q] thetaB2;
  real sg;
  real se;
  real slambda;
  real<lower=0> sy;
}

// transformed parameters{
//   vector[Q] lambda;
//
//   matrix[B1,Q] beta1;
//   matrix[B2, Q] beta2;
//
//   matrix[B1, Q] thetaB1New;
//   vector[Q] sqrtThetaB1;
//   vector[Q] mB1;
//
//   matrix[B2, Q] thetaB2New;
//   vector[Q] sqrtThetaB2;
//   vector[Q] mB2;
//
//    for(q in 1:Q){
//     mB1[q] = sum(thetaB1[1:B1,q])/B1;
//     for(i in 1:B1){
//     thetaB1New[i,q] = (thetaB1[i,q] - mB1[q])^2;
//     }
//     sqrtThetaB1[q] = sqrt(1/(sum(thetaB1New[1:B1,q])));
//     for(i in 1:B1){
//       beta1[i,q] = thetaB1New[i,q]*sqrtThetaB1[q];
//     }
//   }
//
//   for(q in 1:Q){
//     mB2[q] = sum(thetaB2[1:B2,q])/B2;
//     for(j in 1:B2){
//     thetaB2New[j,q] = (thetaB2[j,q] - mB2[q])^2;
//     }
//     sqrtThetaB2[q] = sqrt(1/(sum(thetaB2New[1:B2,q])));
//     for(j in 1:B2){
//       beta2[j,q] = thetaB2New[j,q]*sqrtThetaB2[q];
//     }
//   }
//
//   lambda = sort_desc(lambda_raw);
//
// }

model {
  vector[N] mu;
  // vector[N] blin;
  //
  // for (n in 1:N) {
  //   blin[n] = 0;
  //  }
  //
  // for (n in 1:N) {
  //   for (q in 1:Q) {
  //     blin[n] = blin[n] + (lambda[q] * beta1[var1[n], q] * beta2[var2[n], q]);
  //   }
  // }

  for (n in 1:N) {
    y[n] ~ normal(mu[n], sy);
    mu[n] = muall + gen[var1[n]]; //+ env[var2[n]]; // + blin[n];
   }

  // Priors
  // Prior on grand mean
   muall ~ normal(mmu, smu^-1);

  // Prior on genotype effect
  for(i in 1:B1) {
  gen[i] ~ normal(0, sg^-2);
  }

  // Prior on environment effect
  for(j in 1:B2) {
  env[j] ~ normal(0, se^-2);
  }
  //
  // // Priors on gamma
  // for(q in 1:Q){
  //   for(i in 1:B1){
  //     thetaB1[i,q] ~ normal(0,stheta);
  //   }
  // }

   // Priors on delta
  //  for(q in 1:Q){
  //   for(j in 1:B2){
  //     thetaB2[j,q] ~ normal(0,stheta);
  //   }
  // }
  //
  //
  // // Prior on eigenvalues
  // for(q in 1:Q) {
  //   lambda_raw[q] ~ normal(0, slambda^-2)T[0,];
  // }
  //
  //
  // Priors
  sg ~ student_t(0, 1, 1)T[0,];
  se ~ student_t(0, 1, 1)T[0,];
  // slambda ~ student_t(0, 1, 1)T[0,];


  // Prior on residual standard deviation
   sy ~ gamma(a, b); // inverse of sy
}

