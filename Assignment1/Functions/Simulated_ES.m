function ES_sim = Simulated_ES(nobs, a, b, c, d, xi, seed)
X = stabgen(nobs, a, b, d, c, seed);
q = quantile(X, xi);
Plo = X(X < q);
ES_sim = mean(Plo);
end
