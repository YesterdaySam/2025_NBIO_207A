%function [p] = gaussian(x, mu, sigma)
%p = 1 / (sigma * sqrt(2*pi)) * exp(-1/2 * ((x - mu)/sigma).^2);


% 2) Gaussian distributions 
% ----------------------------
x = -6:.0012:6;
f = gaussian(x, 0, 1);
figure(100); hold on; plot(f);
%a. >1 stdev from mean
a = find((x > 1) | (x < -1));
p = sum(f(a))/sum(f)
0.3172
%b. >2 stdev from mean
a = find((x > 2) | (x < -2));
p = sum(f(a))/sum(f)
0.0455
%c. >3 stdev from mean
a = find((x > 3) | (x < -3));
p = sum(f(a))/sum(f)
0.0027



% 3) Gaussian CDF
%-------------
c = cumsum(f) ./ sum(f);
figure(100); hold on; plot(c, 'g');
s = [];
r = rand(1000,1);
for i = 1:1000
    xr(i) = x(min(find(c > r(i))));
end
%a > 1 stdev from mean
length(find((xr < -1) | (xr > 1))) / 1000
0.3070
%b. > 2 stdev from mean
length(find((xr < -2) | (xr > 2))) / 1000
0.049
%c. > 3 stdev from mean
length(find((xr < -3) | (xr > 3))) / 1000
0.0040






% 4) Distance between gaussian and binomial
% -------------------------------------------
pvals = 0.05:0.05:.5;
xvalues = 1:50;
n = 50;
gdist = zeros(length(pvals),1);
figure; hold on;
for pind = 1:length(pvals)
    subplot(5,2,pind);hold on;
    p1 = pvals(pind);
    p0 = 1-p1;
    % calculate the mean and stdev for the gaussian
    mn = n * p1;
    sigma = (n * p1 * p0) ^ 0.5;
    for xind = 1:length(xvalues)
        x = xvalues(xind);
        binomial_dist(xind) = nchoosek(n,x) * p1^x * p0^(n-x);
    end
    gaussian_dist = gaussian(xvalues, mn, sigma);
    plot(binomial_dist,'k');
    plot(gaussian_dist, 'r');
    gdist(pind) = sum(sqrt((gaussian_dist - binomial_dist) .^2));
end
plot(pvals, gdist);
title('Distanct between Gaussian and Binomial')
xlabel('p(1)');
ylabel('Distance');





% 5) Draws from distribution
% --------------------
%distribution
figure; hold on; 
mu = 10;
sigma = 5;
ndraws = 100000;
y = normrnd(mu, sigma, ndraws, 1);
[counts bins] = hist(y, 100);
bar(bins, counts./sum(counts));

%a. generate 1000 draws of five elements
nsamp = 5;
ind = ceil(rand(nsamp, 1000) * ndraws);
means = mean(y(ind));
mn = repmat(means, nsamp, 1);
stdn = sqrt(sum((y(ind) - mn).^2)/nsamp);
stdn1 = sqrt(sum((y(ind) - mn).^2)/(nsamp - 1));
sten = stdn / sqrt(nsamp);
sten1 = stdn1 / sqrt(nsamp);

% comparison of standard deviation estimates
mean(stdn)
4.25
mean(stdn1)
4.75
% stdn1 is closer (but still an underestimate
%b. standard error of the mean
std(means)
2.24
mean(sten)
1.90
mean(sten1)
2.12


%sten1 is closer

% pdf and cdf of std error estimates
% --------------------------------
figure; hold on;
subplot(1,2,1)
[counts b] = hist(sten1, 50);
counts = counts ./ sum(counts);
bar(b, counts);
xlabel('Standard Error');
ylabel('Probability');
%cdf
sten1c = sort(sten1);
c = cumsum(sten1c);
c = c ./ c(end);
subplot(1,2,2)
plot(sten1c, c)
xlabel('Standard Error');
ylabel('Probability');