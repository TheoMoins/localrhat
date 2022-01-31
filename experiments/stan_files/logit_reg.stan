data {
  int<lower=1> D;
  int<lower=0> N;
  real<lower=0> tau;
  int<lower=0, upper=1> y[N];
  real x[N,D];
}

parameters {
  real bet[D];
}

model {
  for (d in 1:D)
    bet[d] ~ normal(0, tau);    

  for (n in 1:N)
    y[n] ~ bernoulli(inv_logit(dot_product(x[n], bet)));
}
