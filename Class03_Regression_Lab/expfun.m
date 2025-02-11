function [y] = expfun(beta, x)
y = beta(1) * exp(-beta(2) * x);
