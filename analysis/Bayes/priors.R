
n=10000
sigma_mean = rnorm(n,10,10);
sigma_sd = abs(rnorm(n,0,1));
sigma = rlnorm(n,sigma_mean,sigma_sd)

plot(density(sigma))


n=10000
sigma_mean = rnorm(n,50,10);
sigma_sd = abs(rnorm(n,0,10));
sigma = rnorm(n,sigma_mean,sigma_sd)

plot(density(sigma))

pnorm((76.5788-100)/4.00138)