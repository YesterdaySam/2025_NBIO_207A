
function F = gaussian(x, xdata);
% Fit gauusian distribution of length x centred at mu and of width sigma 

mu = x(1); sigma = x(2);
F = (1/sigma*(sqrt(2*pi)))*exp(-(xdata-mu).^2/(2*sigma^2));


% 