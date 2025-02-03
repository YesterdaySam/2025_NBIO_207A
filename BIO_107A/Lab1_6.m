

% Lab 1.6
% ---------
% Distributions - estimates

sample1 = generate_random_data(20,'normal',0,1); 
sample2 = generate_random_data(20,'normal',2,1); % 
[X1,Y1] = cumhist(sample1, [-10 10], 0.1);
[X2,Y2] = cumhist(sample2, [-10 10], 0.1);
figure;
plot(X1,Y1,'r-');
hold on
plot(X2,Y2,'b-');
xlabel('values of 2 samples of size 20');
ylabel('Percent of data');

[h,pvalue] = kstest2(sample1,sample2);


% Mean values - distr

S = [ 0.5 7 10 4 3 5 ]; % a data set with 6 elements
number_of_elements = length(S);  % length(S) returns the number of elements of S
R = randperm(number_of_elements); % returns a random permutation
R(1:3); % the first 3 elements of the random permutation
newsample = S(R(1:3));  % the new sample of 3 points from the data set


cd data
squirrel_births = load('squirrel_births.txt','-ascii');
cd ..
[samplemeans_squirrelbirths, truemean_squirrelbirths] = simulate_random_sampling(squirrel_births, 30, 1000);
title('The sample mean of squirrel birth days for 1000 sampling experiments');

mean(samplemeans_squirrelbirths), 
truemean_squirrelbirths
sample30_prc = prctile(squirrel_births, 95)

simulate_random_sampling(squirrel_births, 10, 1000);
title('Sample means of squirrel birth days for 1000 experiments, 10 samples');

sample10_prc = prctile(squirrel_births, 95)

simulate_random_sampling(squirrel_births, 30, 1000);
title('Sample means of squirrel birth days for 1000 experiments, 30 samples');

simulate_random_sampling(squirrel_births, 60, 1000);
title('Sample means of squirrel birth days for 1000 experiments, 60 samples');

prctile(squirrel_births, 95)
sample60_prc = prctile(squirrel_births, 95)

% Shape of distribution
% Uniform distribution

cd data
human_births = load('human_births.txt','-ascii');
cd ..
[samplemeans_humanbirths, truemeans_humanbirths]  = simulate_random_sampling(human_births, 30, 1000);
title('Sample means of birth days for 1000 sampling experiments, 30 samples');
mean(samplemeans_humanbirths), 
truemeans_humanbirths


cd data
heads_coinflips = load('heads_coinflips.txt','-ascii');
cd ..
[samplemeans_coinflips, truemeans_coinflips] = simulate_random_sampling(heads_coinflips, 30, 1000);
title('Sample means of number of heads for 1000 sampling experiments, 30 samples');
mean(samplemeans_coinflips), 
truemeans_coinflips




% Normal distribution and Central Limit Theorem
simulate_random_sampling(squirrel_births, 30, 1000);
title('Sample means of squirrel birth days for 1000 experiments, 30 samples');

R = randperm(length(squirrel_births));
S_30 = squirrel_births(R(1:30));

Sm = mean(S_30);
S_standarddeviation = std(S_30);  % std calculates the standard deviation
Std_error = S_standarddeviation/sqrt(30);

% Distribution of sample means by central limit theorem
figure;
X = 1:200;
dX = 1;
Nus = dX*exp(-(power(X-Sm,2)/(2*power(Std_error,2))))/sqrt(2*pi*power(Std_error,2));
cumulative = cumsum(Nus);
hold on
plot(X,Nus,'g-');
plot(X,cumulative,'g--');
prctile(cumulative,95)

% Differences between 2 means. T-distribution
R2 = randperm(length(squirrel_births));
S2_30 = squirrel_births(R2(1:30));

[h,pvalue] = ttest2(S_30,S2_30);
pvalue,

%
%T-test of 2 sample means recipe
%Question answered: Are the means of 2 sample distributions S1 and S2 identical?  Technically, we assume that S1 and S2 are samples of quantities that are normally distributed, however the test has been shown to still work for most data that has some central tendency and lacks hard thresholds or edges.

%Statistic calculated:  The T statistic (formula in Baldi).

%How the statistic is distributed when the null hypothesis is true (that is, how it is distributed when the 2 samples have identical means): As T2_CDF(N1, N2, T_statistic), where N1 and N2 are the number of samples of S1 and S2, respectively

%How to perform this test in Matlab: If S1 and S2 are variables with N1 and N2 random samples, then one can perform the K-S test with [h,pvalue] = ttest2(S1,S2); If pvalue is less than the critical value alpha (such as 0.05), then we can reject the null hypothesis and say with 1-alpha confidence (such as 95% confidence) that S1 and S2 have different means.




