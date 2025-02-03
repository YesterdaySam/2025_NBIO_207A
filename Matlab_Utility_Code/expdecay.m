
function F = expdecay(x,xdata)
% x = model coefficients
% xdata = vector of time points


F = x(1)*exp(-xdata/x(2)); 