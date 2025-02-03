function x = sigmoidfit(xdata,ydata)

x0=[-30,30];
x = lsqcurvefit(@sigmoid,x0,xdata,ydata);