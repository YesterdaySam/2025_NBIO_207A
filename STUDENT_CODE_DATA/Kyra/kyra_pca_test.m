%%
cd('F:\Research\Code\2025_NBIO_207A\STUDENT_CODE_DATA\Kyra')

load('KASData.mat')

dat = ForClass;
% Organized into ramp1, ramp2, oa1, oa2 = 4 conditions / experiment
% Each condition has 2 variables Temp, Freq sampled N times where N =
% bursts

% Reorganize into iterable struct
ds(1).ramp(1).T = dat.Ramp1{1}.Temp;
ds(2).ramp(1).T = dat.Ramp1{2}.Temp;
ds(1).ramp(2).T = dat.Ramp2{1}.Temp;
ds(2).ramp(2).T = dat.Ramp2{2}.Temp;
ds(1).ramp(1).F = dat.Ramp1{1}.Freq;
ds(2).ramp(1).F = dat.Ramp1{2}.Freq;
ds(1).ramp(2).F = dat.Ramp2{1}.Freq;
ds(2).ramp(2).F = dat.Ramp2{2}.Freq;
ds(1).ramp(3).T = dat.OA1{1}.Temp;
ds(2).ramp(3).T = dat.OA1{2}.Temp;
ds(1).ramp(4).T = dat.OA2{1}.Temp;
ds(2).ramp(4).T = dat.OA2{2}.Temp;
ds(1).ramp(3).F = dat.OA1{1}.Freq;
ds(2).ramp(3).F = dat.OA1{2}.Freq;
ds(1).ramp(4).F = dat.OA2{1}.Freq;
ds(2).ramp(4).F = dat.OA2{2}.Freq;

%% Plot bursts across temperature by frequency for 1 experiment
expmt = 2;

cmaphot = hot(6);
cmapcool = cool(4);
figure; hold on 
plot(dat.Ramp1{expmt}.Temp, dat.Ramp1{expmt}.Freq,'.','Color',cmaphot(1,:))
plot(dat.Ramp2{expmt}.Temp, dat.Ramp2{expmt}.Freq,'.','Color',cmaphot(2,:))
plot(dat.OA1{expmt}.Temp, dat.OA1{expmt}.Freq,'.','Color',cmapcool(1,:))
plot(dat.OA2{expmt}.Temp, dat.OA2{expmt}.Freq,'.','Color',cmapcool(2,:))
xlabel('Temp')
ylabel('Frequency')

%% PCA on 1 experiment of ramp only
expmt = 2;

nbursts1 = length(dat.Ramp1{expmt}.Temp);
nbursts2 = length(dat.Ramp2{expmt}.Temp);
r1inds = 1:nbursts1;
r2inds = nbursts1+1:(nbursts1+nbursts2);

combineRamp = [dat.Ramp1{expmt}.Temp, dat.Ramp1{expmt}.Freq; dat.Ramp2{expmt}.Temp, dat.Ramp2{expmt}.Freq];

[COEFF, score, LATENT, TSQUARED, EXPLAINED, MU] = pca(combineRamp);

%% Repeat plotting on PCAs

figure; hold on
plot(score(r1inds,1),'.','Color',cmaphot(1,:))
plot(score(r2inds,1),'.','Color',cmaphot(2,:))

%% PCA on many experiments of ramp and OA
expmt = 1;
ct = 0;
cmapCombi = [cmaphot(1:2,:); cmapcool(1:2,:)];
combineDat = [];
sinds = [];

for j = 1:2
    expmt = j;
    for i = 1:4
        combineDat = [combineDat; ds(expmt).ramp(i).T, ds(expmt).ramp(i).F];
        nBursts = length(ds(expmt).ramp(i).T);
        sinds(j,i).ramp = (1:nBursts) + ct;
        ct = ct + nBursts;
    end
end

% Confirm data reorginzation successful
% figure; hold on
% for i = 1:4
%     plot(ds(expmt).ramp(i).T, ds(expmt).ramp(i).F, '.')
% end

[COEFF, score, LATENT, TSQUARED, EXPLAINED, MU] = pca(combineDat);

figure; hold on

for j = 1:2
    expmt = j;
    for i = 1:4
        plot(score(sinds(j,i).ramp,1),'.','Color',cmapCombi(i,:))
    end
end

%% Linear model on PC1
slopes = [];
for j = 1:2
    expmt = j;
    for i = 1:4
        xvals = 1:length(ds(expmt).ramp(i).T);
        tmpmdl = fitlm(xvals,score(sinds(j,i).ramp,1));
        slopes(j,i) = tmpmdl.Coefficients.Estimate(2);
        ps(j,i) = tmpmdl.Coefficients.pValue(2);
    end
end






