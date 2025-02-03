
function F = gaussian(x, mu, sigma);
% generate gauusian distribution of length x centred at mu and of width sigma 

F = (1/(sigma*sqrt(2*pi)))*exp(-(x-mu).^2/(2*sigma^2));

%