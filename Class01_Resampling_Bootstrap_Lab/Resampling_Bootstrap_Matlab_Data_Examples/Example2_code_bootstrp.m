clc
clear all
close all
%%
load('SWRcofiring_Spatialcorr_day1_vs_day45.mat')
%%
% simple bootstrapping
runs = 100000;% 100000 times
resamplenum = 246;% use the minimal number of pairs to match samples
day1_resample_corr = zeros(1,runs);
day45_resample_corr = zeros(1,runs);
for i = 1:runs
    [day1_resample,idx1] = datasample(SWR_cofiring_day1,resamplenum);
    day1_resample_corr(i) = corr(day1_resample,spatial_corr_day1(idx1));
    [day45_resample,idx2] = datasample(SWR_cofiring_day45,resamplenum);
    day45_resample_corr(i) = corr(day45_resample,spatial_corr_day45(idx2));
end

pval_corr = mean(day1_resample_corr < day45_resample_corr)% bootstrappped p-value
