function [p] = gaussian(x, mu, sigma)
p = 1 / (sigma * sqrt(2*pi)) * exp(-1/2 * ((x - mu)/sigma).^2);