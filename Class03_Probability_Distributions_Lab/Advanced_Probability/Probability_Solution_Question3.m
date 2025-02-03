%Probability Problem 3
% ----------------

load Probability_Data1.mat
s1=responses;
% A. Plot a histogram with 1 count bins 
% ---------------------------------
% note that my version of hist may be different from the normal matlab version, so
 % the bins you should specify would be 0:1:100 and 0:5:100 
stimA = 1:2:length(s1); 
stimB = 2:2:length(s1); 
[na1 x1] = hist(s1(stimA),.5:1:100); 
[nb1 x1] = hist(s1(stimB),.5:1:100);
figure; hold on; 
plot(x1, nb1, 'b'); hold on 
plot(x1+.1, na1, 'r'); 
set(gca, 'XLim', [20 80], 'FontSize', 18);


% B. Plot a histogram with 5 count bins 
% ---------------------------------
stimA = 1:2:length(s1); 
stimB = 2:2:length(s1); 
[na5 x5] = hist(s1(stimA),2.5:5:100); 
[nb5 x5] = hist(s1(stimB),2.5:5:100); 
figure; hold on;
plot(x5, nb5, 'b'); 
hold on 
plot(x5+.5, na5, 'r'); 
set(gca, 'XLim', [20 80], 'FontSize', 18); 
legend('Stimulus A', 'Stimulus B');

% C. pdfs and cdfs. 
% ---------------------------------
figure 
% we'll use four subplots, one for each plot 
subplot(2,2,1) 
%to get the pdf for stimulus A responses, we divide na1 by the total number of 
%trials, which is the sum of na1 
pdfa1 = na1./sum(na1); 
plot(x1, pdfa1, 'k-', 'LineWidth', 2); 
set(gca, 'FontSize', 14); 
h = xlabel('# Spikes'); 
set(h, 'FontSize', 16); 
h = ylabel('Probability'); 
set(h, 'FontSize', 16); 
h = title('pdf of stimA responses'); 
set(h, 'FontSize', 16); 
%to get the cdf for stimulus A responses, we take the cumulative sum of na1 and 
%divide by the total sum: 
subplot(2,2,2) 
cdfa1 = cumsum(na1); 
cdfa1 = cdfa1 ./ cdfa1(end); 
% note that you could also use ecdf function 
plot(x1, cdfa1, 'k-', 'LineWidth', 2); 
set(gca, 'FontSize', 14); 
h = xlabel('# Spikes'); 
set(h, 'FontSize', 16); 
h = ylabel('Cumulative probability'); 
set(h, 'FontSize', 16); 
h = title('cdf of stimA responses'); 
set(h, 'FontSize', 16); 
subplot(2,2,3) 
pdfb1 = nb1./sum(nb1); 
plot(x1, pdfb1, 'k-', 'LineWidth', 2); 
set(gca, 'FontSize', 14); 
h = xlabel('# Spikes'); 
set(h, 'FontSize', 16); 
h = ylabel('Probability'); 
set(h, 'FontSize', 16); 
h = title('pdf of stimB responses'); 
set(h, 'FontSize', 16); 
subplot(2,2,4) 
cdfb1 = cumsum(nb1); 
cdfb1 = cdfb1 ./ cdfb1(end); 
plot(x1, cdfb1, 'k-', 'LineWidth', 2); 
set(gca, 'FontSize', 14); 
h = xlabel('# Spikes'); 
set(h, 'FontSize', 16); 
h = ylabel('Cumulative probability'); 
set(h, 'FontSize', 16); 
h = title('cdf of stimB responses'); 
set(h, 'FontSize', 16);




% D. Probability of stimulus A given 44 spikes , 1 count bins
% ---------------------------------
% normalize the histograms to be probabilities
na1 = na1./sum(na1);
nb1 = nb1./sum(nb1);
% find the bin number for a count of 44
[m ba] = min(abs(x1 - 44));
% p (A | 44 spikes) = P(44 spikes | A) P(A) / P(44 spikes)
p = na1(ba) * .5 / (na1(ba) * .5 + nb1(ba) * .5)
0.2500


% E. Probability of stimulus A given 44 spikes , 5 count bins
% ---------------------------------
% normalize the histograms to be probabilities
na5 = na5./sum(na5);
nb5 = nb5./sum(nb5);
% find the bin number for a count of 44
[m ba] = min(abs(x5 - 44));
% p (A | 44 spikes) = P(44 spikes | A) P(A) / P(44 spikes)
p = na5(ba) * .5 / (na5(ba) * .5 + nb5(ba) * .5)
0.3953

%E is likely to be closer because it appears to smooth over noise due to sampling.

