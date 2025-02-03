% [sig] = WATSONU(a1, a2, pval)
%         Given two distributions on a circle, uses the Watson U^2 test to
%         determine the whether the differences between the CDFs are 
%         significant at the specified p value.
%         pval should be .05, .02, .01, .005, .002, or .001
%         Assumes that 0 <= a1, a2 <= 2pi
%         From _Circular Statistics in Biology_ by Edward Batshelet, 1981, Academic
%         Press, pg. 115
function [sig] = watsonu(a1, a2, pval)

% load the table of the estimates of u
load /home/loren/matlab/stats/watsonu.mat;

% make sure a1 and a2 are column vectors 
if (size(a1,1) > 1)
	a1 = a1';
end
if (size(a2,1) > 1)
	a2 = a2';
end

n = length(a1);
m = length(a2);
if (n > m)
	% switch the two variables
	tmp = a1;
	a1 = a2;
	a2 = tmp;
	tmp = n;
	n = m;
	m = tmp;
end

% check n and m for excessively low values
if (min([n m]) < 5)
	sig = 0;
end

% sort a1 and a2 and put in a 1 / nvalues in the second column
a1 = [sort(a1) ; ones(1,n) / n]';
a2 = [sort(a2) ; ones(1,m) / m]';

% add in the points from the other distribution
tmpa1 = a1;
tmp = zeros(m,2);
tmp(:,1) = a2(:,1);
a1 = sortrows([a1 ; tmp], 1);

tmp = zeros(n,2);
tmp(:,1) = tmpa1(:,1);
a2 = sortrows([a2 ; tmp], 1);


% get the cdfs of the two distributions from 0 to 2pi with indeces
a1cdf = cumsum(a1(:,2));
a2cdf = cumsum(a2(:,2));

% get the differences in the cdfs
d = abs(a1cdf - a2cdf);

% calculate U^2
N = n + m;
Usqr = (n * m / N^2) * (sum(d.^2) - 1/N * (sum(d))^2);

% get the relevant column of the U^2 table
alphacol = find(watsonUalphas == pval);
if (isempty(alphacol))
	error('Error: WATSONU pval invalid');
end

% look up the value in the U^2 table
if ((n >= 12) | (m >= 12))
	uval = watsonUsqr(alphacol, 11, 1);
else
	uval = watsonUsqr(alphacol, n, m);
end


if (Usqr > uval)
	sig = 1;
else
	sig = 0;
end
