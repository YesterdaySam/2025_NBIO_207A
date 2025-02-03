function [p,littleflag] = myfcdf(x,v1,v2);
%FCDF   F cumulative distribution function.
%   P = FCDF(X,V1,V2) returns the F cumulative distribution function
%   with V1 and V2 degrees of freedom at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2007/07/27 16:45:00 $
%**************************************************************************
%**************************************************************************
%    Kluged to return non-one number under conditions in which p
%    would be approximated by 1.  This allows this function to
%    return p values very close to 1. Work done by:
%              A. Katrin Schenk (1,2) and Evan H. Goulding (3)
%              (1)UCSF Department of Neurology 
%              (2)Sloan-Swartz Center for Theoretical Neuroscience
%              (3)UCSF Department of Psychiatry
%              University of California, San Francisco
%              Rock Hall
%              1550 4th st.
%              San Francisco, CA 94501
%              schenk@phy.ucsf.edu, evan.goulding@ucsf.edu
%              May 20, 2007.
%***************************************************************************
%***************************************************************************

if nargin < 3, 
    error('Requires three input arguments.'); 
end

[errorcode x v1 v2] = distchck(3,x,v1,v2);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

%   Initialize P to zero.
p = zeros(size(x));

t = (v1 <= 0 | v2 <= 0 | isnan(x) | isnan(v1) | isnan(v2));
p(t) = NaN;
s = (x==Inf) & ~t;
if any(s)
   p(s) = 1;
   t = t | s;
end

% Compute P when X > 0.
k = find(x > 0 & ~t & isfinite(v1) & isfinite(v2));
if any(k), 
% use A&S formula 26.6.2 to relate to incomplete beta function 
    % Also use 26.5.2 to avoid cancellation by subtracting from 1
    xx = x(k)./(x(k) + v2(k)./v1(k));
    %[p(k),littleflag] = mybetainc(xx, v1(k)/2, v2(k)/2);
    [p(k)] = betainc(xx, v1(k)/2, v2(k)/2); littleflag = 0;
    
end

if any(~isfinite(v1(:)) | ~isfinite(v2(:)))
   k = find(x > 0 & ~t & isfinite(v1) & ~isfinite(v2) & v2>0);
   if any(k)
      p(k) = chi2cdf(v1(k).*x(k),v1(k));
   end
   k = find(x > 0 & ~t & ~isfinite(v1) & v1>0 & isfinite(v2));
   if any(k)
      p(k) = 1 - chi2cdf(v2(k)./x(k),v2(k));
   end
   k = find(x > 0 & ~t & ~isfinite(v1) & v1>0 & ~isfinite(v2) & v2>0);
   if any(k)
      p(k) = (x(k)>=1);
   end
end
