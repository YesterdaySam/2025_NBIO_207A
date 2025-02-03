clc
clear all;
close all;
%%
% load the data
load('ValleyTime_SWRinhpairs.mat')
%%
% take a look at the data
bin_width = 20; % ms
bin_edges = -500:bin_width:500;

awakecounts = histcounts(awake,bin_edges); 
awakeprc = awakecounts./sum(awakecounts).*100;% 0-100%
awakepairnum = length(awake);

sleepcounts = histcounts(sleep,bin_edges); 
sleepprc = sleepcounts./sum(sleepcounts).*100;% 0-100%
sleeppairnum = length(sleep);

% plot histograms
bin_centers = bin_edges(1:end-1) + bin_width/2;
figure('position',[100 1000 500 300],'color','w'),
bar(bin_centers,awakeprc,'cyan');
hold on
stairs(bin_edges,[sleepprc,sleepprc(end)],'k','linewidth',2)
plot([0 0],[0, max(max(awakeprc),max(sleepprc))],'r--')
ylabel('Percentage' ,'fontsize',14);
xlabel('Valley time (ms)','fontsize',14);
legend(['Awake group, N=' int2str(awakepairnum)],['Sleep group, N=' int2str(sleeppairnum)]);
box off;

% see the kurtosis of the original distributions
awake_kurtosis = kurtosis(awake) % normal distribution as 0
sleep_kurtosis = kurtosis(sleep)
%%
% simple bootstrapping
runs = 1000;% 1000 times
resamplenum = min(awakepairnum,sleeppairnum);% use the minimal number of pairs to match samples
awake_resample_k = zeros(1,runs);
sleep_resample_k = zeros(1,runs);
for i = 1:runs
    awake_resample = datasample(awake,resamplenum);
    sleep_resample = datasample(sleep,resamplenum);
    awake_resample_k(i) = kurtosis(awake_resample)-3;
    sleep_resample_k(i) = kurtosis(sleep_resample)-3;
end

pval_kurtosis = mean(awake_resample_k < sleep_resample_k)% bootstrappped p-value

%%
% plot histograms of bootstrapped results
minedge = min(min(awake_resample_k),min(sleep_resample_k));
maxedge = max(max(awake_resample_k),max(sleep_resample_k));
bin_width = (maxedge-minedge)/40; % 40 bins
bin_edges = minedge:bin_width:maxedge;

awake_kcounts = histcounts(awake_resample_k,bin_edges); 
awake_kprc = awake_kcounts./sum(awake_kcounts).*100;% 0-100%

sleep_kcounts = histcounts(sleep_resample_k,bin_edges); 
sleep_kprc = sleep_kcounts./sum(sleep_kcounts).*100;% 0-100%

%%
% plot histogram of the bootstrap samples
bin_centers = bin_edges(1:end-1) + bin_width/2;
figure('position',[700 600 500 300],'color','w'),
bar(bin_centers,awake_kprc,'r');
hold on
bar(bin_centers,sleep_kprc,'cyan')
plot([awake_kurtosis,awake_kurtosis],[0, max(awake_kprc)],'k','linewidth',2)
plot([sleep_kurtosis,sleep_kurtosis],[0, max(sleep_kprc)],'b','linewidth',2)

ylabel('Percentage' ,'fontsize',14);
xlabel('Kurtosis','fontsize',14);
title(['Histogram of the bootstrap samples (N=' int2str(runs),')'],'fontsize',14)
legend('Awake','Sleep');
box off;


