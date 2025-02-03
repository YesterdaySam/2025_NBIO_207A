
function F = sigmoid(x,xdata)
% x = model coefficients
% xdata = vector of voltage points


F = 1./(1+exp(-(xdata-x(1))/x(2))); 