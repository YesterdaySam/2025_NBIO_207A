
function x = mygaussianfit(x0,xdata,ydata)

x = lsqcurvefit(@mygaussian,x0,xdata,ydata);