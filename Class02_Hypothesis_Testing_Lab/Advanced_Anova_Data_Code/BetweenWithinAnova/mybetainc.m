function [y,littleflag] = mybetainc(x,a,b)
%BETAINC Incomplete beta function.
%   Y = BETAINC(X,Z,W) computes the incomplete beta function for
%   corresponding elements of X, Z, and W.  The elements of X must be
%   in the closed interval [0,1].  The arguments X, Z and W must all
%   be the same size (or any of them can be scalar).
%
%   The incomplete beta function is defined as
%
%     I_x(z,b) = 1./BETA(z,w) .* 
%                 integral from 0 to x of t.^(z-1) .* (1-t).^(w-1) dt
%
%   See also BETA, BETALN.

%   Ref: Abramowitz & Stegun, Handbook of Mathematical Functions, sec. 26.5,
%   especially 26.5.8, 26.5.20 and 26.5.21.
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
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2007/07/27 16:45:00 $
littleflag=0;
if nargin < 3
   error('Requires three input arguments.')
elseif any(x(:) < 0 | x(:) > 1 | isnan(x(:))) | ~isreal(x)
   error('X must be in the interval [0,1].')
elseif any(a(:) < 0 | isnan(a(:))) | ~isreal(a)
   error('Z must be nonnegative.')
elseif any(b(:) < 0 | isnan(b(:))) | ~isreal(b)
   error('W must be nonnegative.')
end

try
  % Preallocate y (using the size rules for plus)
  y = x + a + b; 
  y(:) = 0;
catch
  error('X, Z and W must all the same size (or any of them can be scalar).')
end

if ~isempty(y)
  bt = exp(gammaln(a+b)-gammaln(a)-gammaln(b) + ...
       a.*log(x+(x==0)) + b.*log(1-x+(x==1)));

  k = find(x < (a+1) ./ (a+b+2));
  if ~isempty(k)
     if length(x) == 1, xk = x; else, xk = x(k); end
     if length(a) == 1, ak = a; else, ak = a(k); end
     if length(b) == 1, bk = b; else, bk = b(k); end
     y(k) = bt(k) .* betacore(xk,ak,bk) ./ ak;
  end

  k = find(x >= (a+1) ./ (a+b+2));
  if ~isempty(k)
     if length(x) == 1, xk = x; else, xk = x(k); end
     if length(a) == 1, ak = a; else, ak = a(k); end
     if length(b) == 1, bk = b; else, bk = b(k); end
     minusterm = bt(k) .* betacore(1-xk,bk,ak) ./ bk;
     if minusterm<1e-10
       littleflag = 1;
       y(k) = minusterm;
       % Note that the true incomplete beta function is 1 minus the
       % y(k) as define right above this line!  we kluged this for
       % use with p values
     else
       y(k) = 1-minusterm;
     end
  end

  k = find(isnan(y));
  if ~isempty(k)
     % Continued fraction in betacore failed, use approximations.
     if length(x) == 1, xk = x; else, xk = x(k); end
     if length(a) == 1, ak = a; else, ak = a(k); end
     if length(b) == 1, bk = b; else, bk = b(k); end
     w1 = (bk*xk).^(1/3);
     w2 = (ak*(1-xk)).^(1/3);
     y(k) = 0.5*erfc(-3/sqrt(2)*((1-1/(9*bk))*w1-(1-1/(9*ak))*w2)./ ...
        sqrt(w1.^2/bk+w2.^2/ak));
     k = find((ak+bk-1).*(1-xk) < 0.8);
     if ~isempty(k)
        if length(x) == 1, xk = x; else, xk = xk(k); end
        if length(a) == 1, ak = a; else, ak = ak(k); end
        if length(b) == 1, bk = b; else, bk = bk(k); end
        s = 0.5*((ak+bk-1)*(3-xk)-(bk-1)).*(1-xk);
        y(k) = 1-gammainc(s,bk);
     end
  end

  k = find(x == 0);
  if length(k) > 0
     if length(x) == 1, y(:) = 0; else, y(k) = 0; end
  end

  k = find(x == 1);
  if length(k) > 0
     if length(x) == 1, y(:) = 1; else, y(k) = 1; end
  end
end
