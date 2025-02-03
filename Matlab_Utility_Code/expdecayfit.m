
function x = expdecayfit(ydata,A,e)

xdata = 1:length(ydata);
x0=[A e];

x = lsqcurvefit(@expdecay,x0,xdata,ydata);