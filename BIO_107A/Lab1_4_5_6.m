% Lab 1.4
% ---------
% Cumulative histogram
cd data;

sample1 = load('change_mineral_breastfeeding.txt','-ascii'); 
sample2 = load('change_mineral_other.txt','-ascii');
cd ..

sort([sample1; sample2]);
unique(sample1);
unique([sample1; sample2]);

bin_edges = unique([sample1; sample2]);
bin_counts1 = histc(sample1, bin_edges);
bin_counts2 = histc(sample2, bin_edges);

bin_edges = unique([sample1; sample2]);
bin_counts1 = histc(sample1, bin_edges);
bin_counts2 = histc(sample2, bin_edges);

figure;
bar(bin_edges, bin_counts1);
title('Histogram 1');
figure;
bar(bin_edges, bin_counts2);
title('Histogram 2');

A = [ 1 2 3 4 5];
cumsum(A);

cumsum1 = cumsum(bin_counts1) / sum(bin_counts1);
cumsum2 = cumsum(bin_counts2) / sum(bin_counts2);

figure;
plot(bin_edges,cumsum1,'b-');
hold on;
plot(bin_edges,cumsum2,'r-');
xlabel('Sample values');
ylabel('Cumulative sum of samples');

% Distance between curves

A = [ 1 2 3 4 5];
B = [ 9 8 7 6 5];
B-A;

abs(-1);
abs(B-A);
abs(A-B);

[maxvalueA,maxlocationA] = max(A);
[maxvalueB,maxlocationB] = max(B);
[maxdiff,maxdiff_location] = max(abs(B-A));






% Lab 1.5
% ---------
% Statistical Inference

sample1 = generate_random_data(20,'normal',0,1); 
sample2 = generate_random_data(20,'normal',0,1);
[X1,Y1] = cumhist(sample1,[-10 10],0.1);
[X2,Y2] = cumhist(sample2,[-10 10],0.1);
figure(90);
plot(X1,Y1,'bo-');
hold on;
plot(X2,Y2,'ro-');
xlabel('Value');
ylabel('Percent of samples');
axis([-2 2 0 100]);

[maxdiff,maxdiff_loc,Xval,sample1CDF,sample2CDF] = cumulative_hist_diff(sample1,sample2);
figure(90);
Xvalues = [ Xval(maxdiff_loc) Xval(maxdiff_loc) ];
Yvalues = 100 * [ sample1CDF(maxdiff_loc) sample2CDF(maxdiff_loc) ];
plot(Xvalues,Yvalues,'k-','linewidth',2);
maxdiff,
maxdiff_o=maxdiff


% while loops
i = 1;
while (i<10)
    i;
    i = i + 1;
end

md = [];
i = 1;
while (i<1001),
    sample1 = generate_random_data(20,'normal',0,1);
    sample2 = generate_random_data(20,'normal',0,1);
    maxdiff = cumulative_hist_diff(sample1,sample2);
    md(i) = maxdiff;
    i = i + 1;
end;

size(md);
[Xmd,Ymd] = cumhist(md,[0 1],0.01);
figure(100);
plot(Xmd,Ymd,'k-');
xlabel('Max diff b/w 2 samples of size 20');
ylabel('Percent of data');


% 2 distributions - differences
sample1 = generate_random_data(20,'normal',0,1); 
sample2 = generate_random_data(20,'normal',2,1); % Note different distribution
[X1,Y1] = cumhist(sample1, [-10 10], 0.1);
[X2,Y2] = cumhist(sample2, [-10 10], 0.1);
maxdiff = cumulative_hist_diff(sample1,sample2);
maxdiff_n = maxdiff

figure(101);
plot(X1,Y1,'r-');
hold on
plot(X2,Y2,'b-');
xlabel('values of 2 samples of size 20');
ylabel('Percent of data');
figure(100);
hold on;
plot([maxdiff maxdiff],[0 100],'g-'); % we'll plot a green vertical bar


% Likelihood that 2 samples of data are same (or different)
% ---------------------------------------------------------
%ks test

x = 0:0.01:1; % create the X points
n1 = 20;
n2 = 20;
math_cdf = ks2_cdf(n1, n2, x);
math_cdf_percent = 100*math_cdf; % convert to percent rather than fraction
figure(100);
hold on;
plot(x,math_cdf_percent,'k--'); % black dashed line

% Hypothesis Testing

%Ho and H1






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

[h,pvalue] = kstest2(sample1,sample2)


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

samplemeans_squirrelbirths, 
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
samplemeans_humanbirths, 
truemeans_humanbirths


cd data
heads_coinflips = load('heads_coinflips.txt','-ascii');
cd ..
[samplemeans_coinflips, truemeans_coinflips] = simulate_random_sampling(heads_coinflips, 30, 1000);
title('Sample means of number of heads for 1000 sampling experiments, 30 samples');
samplemeans_coinflips, 
truemeans_coinflips




% Normal distribution and Central Limit Theorem
simulate_random_sampling(squirrel_births, 30, 1000);
title('Sample means of squirrel birth days for 1000 experiments, 30 samples');

R = randperm(length(squirrel_births));
S_30 = squirrel_births(R(1:30));

Sm = mean(S_30);
S_standarddeviation = std(S_30);  % std calculates the standard deviation
Std_error = S_standarddeviation/sqrt(30);

figure;
X = 1:200;
dX = 1;
Nus = dX*exp(-(power(X-Sm,2)/(2*power(Std_error,2))))/sqrt(2*pi*power(Std_error,2));
cumulative = cumsum(Nus);
hold on
plot(X,Nus,'g-');
plot(X,cumulative,'g--');


% Differences between 2 means
R2 = randperm(length(squirrel_births));
S2_30 = squirrel_births(R2(1:30));

[h,pvalue] = ttest2(S_30,S2_30);
pvalue,





