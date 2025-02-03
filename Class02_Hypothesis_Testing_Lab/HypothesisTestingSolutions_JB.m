% 03: Hypothesis testing problem set
% Jordan Breffle
% 2/25/2021

gitHubFolder = 'C:\Users\jtbre\Documents\GitHub\';


%% Question 1

%% 1a %
alpha = 0.05;

g1mean = 20; % group 1 mean
g2mean = 30;

g1STD = 10;
g2STD = 10;

n = 10;     % number of mice in each group
N = 1000;   % Number of simulated experiments

simPvals = zeros(size(N));
for i = 1:N
    g1 = (randn(1, n)*g1STD)+g1mean ;
    g2 = (randn(1, n)*g2STD)+g2mean ;
    
    [H, P] = ttest2(g1, g2, 'alpha', alpha);
    simPvals(i) = P;
end

figure; histogram(simPvals, 'Normalization', 'Probability')
xlabel('p-value'); ylabel('Probability')

% Approximately 56%
TruePosPerc = sum(simPvals<alpha)/N*100


%% 1bi %
alpha = 0.05;

g1mean = 25; % group 1 mean
g2mean = 25;

g1STD = 10;
g2STD = 10;

n = 10;     % number of mice in each group
N = 1000;   % Number of simulated experiments

simPvals = zeros(size(N));
for i = 1:N
    g1 = (randn(1, n)*g1STD)+g1mean ;
    g2 = (randn(1, n)*g2STD)+g2mean ;
    
    [H, P] = ttest2(g1, g2, 'alpha', alpha);
    simPvals(i) = P;
end

figure; histogram(simPvals, 'Normalization', 'Probability')
xlabel('p-value'); ylabel('Probability')

% Approximately 5%
FalsePosPerc = sum(simPvals<alpha)/N*100



%% 1bii %
g1mean = 25;
g2mean = 25;

g1STD = 10;
g2STD = 10;

n = 10;
N = 1000;

alpha_vec = 0.005:0.005:0.10;
simHvalMat = zeros(numel(alpha_vec), N);

% increment through each alpha in alpha_vec, for each alpha do the same simulation as in 1bi
for j = 1:numel(alpha_vec)
    alpha = alpha_vec(j);

    simPvals = zeros(size(N));
    for i = 1:N
        g1 = (randn(1, n)*g1STD)+g1mean ;
        g2 = (randn(1, n)*g2STD)+g2mean ;
        
        [H, P] = ttest2(g1, g2, 'alpha', alpha);
        simHvalMat(j, i) = H;
    end
    
end

figure; hold on;
plot([0, max(alpha_vec)*1.15], [0, max(alpha_vec)*1.15], 'k--')
errorbar(alpha_vec, mean(simHvalMat, 2), std(simHvalMat, [], 2)/sqrt(N))
xlabel('t-test \alpha'); ylabel('Pr(False Positive)')

xlim([0, max(alpha_vec)*1.15])
ylim([0, max(alpha_vec)*1.15])

% the false positive rate is the same as the alpha


%% 1c %
alpha = 0.05;

gmean = 25;
gSTD = 10;

n = 10;
N = 1000;

simHvals = zeros(size(N));
for i = 1:N
    
    % generate data for the four groups
    g1 = (randn(1, n)*gSTD)+gmean ;
    g2 = (randn(1, n)*gSTD)+gmean ;
    g3 = (randn(1, n)*gSTD)+gmean ;
    g4 = (randn(1, n)*gSTD)+gmean ;
    
    % isPos=1 if any group is different from any of the others
    isPos = ttest2(g1, g2, 'alpha', alpha) | ttest2(g1, g3, 'alpha', alpha) | ttest2(g1, g4, 'alpha', alpha) | ...
                ttest2(g2, g3, 'alpha', alpha) | ttest2(g2, g4, 'alpha', alpha) | ttest2(g3, g4, 'alpha', alpha) ;
    simHvals(i) = isPos;
end

figure; histogram(simHvals, 'Normalization', 'Probability')
xlabel('Had false-positive'); ylabel('Probability')

% Approximately 20%
FalsePosPerc = sum(simHvals)/N*100


%% Question 2
clear g1 g2 g3 g4

part = 'A'; % Choose 'A', 'B', or 'C'

switch part
    case 'A'
        % Only significant group effect
        g1{1} = normrnd(50 * ones(1,10), 3 * ones(1,10));
        g1{2} = normrnd(50 * ones(1,10), 3 * ones(1,10));
        g1{3} = normrnd(50 * ones(1,10), 3 * ones(1,10));

        g2{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g2{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g2{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));

        g3{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g3{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g3{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));
    case 'B'
        % Only significant time effect
        g1{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g1{2} = normrnd(50 * ones(1,10), 3 * ones(1,10));
        g1{3} = normrnd(500 * ones(1,10), 3 * ones(1,10));

        g2{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g2{2} = normrnd(50 * ones(1,10), 3 * ones(1,10));
        g2{3} = normrnd(500 * ones(1,10), 3 * ones(1,10));

        g3{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g3{2} = normrnd(50 * ones(1,10), 3 * ones(1,10));
        g3{3} = normrnd(500 * ones(1,10), 3 * ones(1,10));
    case 'C'
        % Only interaction effect
        g1{1} = normrnd(10 * ones(1,10), 3 * ones(1,10));
        g1{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g1{3} = normrnd(0 * ones(1,10), 3 * ones(1,10));

        g2{1} = normrnd(0 * ones(1,10), 3 * ones(1,10));
        g2{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g2{3} = normrnd(10 * ones(1,10), 3 * ones(1,10));

        g3{1} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g3{2} = normrnd(5 * ones(1,10), 3 * ones(1,10));
        g3{3} = normrnd(5 * ones(1,10), 3 * ones(1,10));
end

y = [];
group = [];
time = [];
for i = 1:3
    y = [y g1{i}];
    group = [group ones(size(g1{i}))];
    time = [time i*ones(size(g1{i}))];
end
for i = 1:3
    y = [y g2{i}];
    group = [group 2*ones(size(g2{i}))];
    time = [time i*ones(size(g2{i}))];
end
for i = 1:3
    y = [y g3{i}];
    group = [group 3*ones(size(g3{i}))];
    time = [time i*ones(size(g3{i}))];
end

[p atab stats] = anovan(y, {group time}, 3, 3, strvcat('group', 'time'));
comp = multcompare(stats, .05, 'off',[] ,[] ,[1 2]);

disp(['Question 2, part ', part, ', p-vals: ', regexprep(num2str(round(p, 4)'),'\s+',', ')])


%% Question 3

load([gitHubFolder, 'NBIO207A_Spring2021\03_Hypothesis_Testing\ps3_3.mat'])
figure; hold on
bins = 0:5:50;
histogram(group1, bins)
histogram(group2, bins)
legend({'Group 1', 'Group 2'})
xlabel('Escape Latency (s)'); ylabel('Count (mice)')
title('Question 3')

disp('Question 3:')
pVar = vartestn([group1', group2'], 'Display', 'off')

if pVar<0.05
    [H, P] = ttest2(group1, group2, 'Vartype','unequal')
else
    [H, P] = ttest2(group1, group2)
end

std(group1)
std(group2)

mean(group1)
mean(group2)

% A) An unpaired 2-sample t-test, with unpooled variance should be used
% B) They are statistically significantly different
% C) Group 2 is faster at the water maze than group 1


%% Question 4

load([gitHubFolder, 'NBIO207A_Spring2021\03_Hypothesis_Testing\ps3_4.mat'])
figure; hold on
bins = 0:5:50;
histogram(group1, bins)
histogram(group2, bins)
legend({'Group 1', 'Group 2'})
xlabel('Escape Latency (s)'); ylabel('Count (mice)')
title('Question 4')

disp('Question 4:')
pVar = vartestn([group1', group2'], 'Display', 'off')

if pVar<0.05
    [H, P] = ttest2(group1, group2, 'Vartype','unequal')
else
    [H, P] = ttest2(group1, group2)
end

std(group1)
std(group2)

mean(group1)
mean(group2)

% A) An unpaired 2-sample t-test, with pooled variance should be used
% B) They are statistically significantly different
% C) Group 2 is faster at the water maze than group 1




%% Question 5

load([gitHubFolder, 'NBIO207A_Spring2021\03_Hypothesis_Testing\ps3_5.mat'])
figure; hold on
bins = 0:5:50;
histogram(group1, bins)
histogram(group2, bins)
legend({'Group 1', 'Group 2'})
xlabel('Escape Latency (s)'); ylabel('Count (mice)')
title('Question 5')

disp('Question 5:')
pVar = vartestn([group1', group2'], 'Display', 'off')

if pVar<0.05
    [H, P] = ttest2(group1, group2, 'Vartype','unequal')
else
    [H, P] = ttest2(group1, group2)
end
%[H, P] = ttest2(group1, group2)
% They are different


std(group1)
std(group2)

mean(group1)
mean(group2)

% A) An unpaired 2-sample t-test, with unpooled variance should be used
% B) They are statistically significantly different
% C) Group 2 is faster at the water maze than group 1


%% Question 6

% spike times across 50 trials for two different cells
% a stimulus is presented at t=0.5 seconds

load([gitHubFolder, 'NBIO207A_Spring2021\03_Hypothesis_Testing\ps3_6.mat'])
nTrials = numel(spikes1);
stimT = 0.5; % time of stimulation (seconds)


%% 6a
spikes1PreFR = cellfun(@(x)sum(x<stimT), spikes1)./stimT;
spikes2PReFR = cellfun(@(x)sum(x<stimT), spikes2)./stimT;

c1PreMean = mean(spikes1PreFR)
c1PreVar = var(spikes1PreFR)

c2PreMean = mean(spikes2PReFR)
c2PreVar = var(spikes2PReFR)


%% 6b

binSize = 0.1;
tMax = 2;
nBins = tMax/binSize;

zscoreMat1 = zeros(nBins, nTrials);
zscoreMat2 = zeros(nBins, nTrials);
FrMat1 = zeros(nBins, nTrials);
FrMat2 = zeros(nBins, nTrials);
for i = 1:nBins
    binFR1 = [cellfun(@(x)sum(x<(i*binSize) & x>((i-1)*binSize)), spikes1)./binSize]';
    binZscores1 = (binFR1-c1PreMean) / c1PreVar;
    zscoreMat1(i,:) = binZscores1;
    FrMat1(i,:) = binFR1;

    binFR2 = [cellfun(@(x)sum(x<(i*binSize) & x>((i-1)*binSize)), spikes2)./binSize]';
    binZscores2 = (binFR2-c2PreMean) / c2PreVar;
    zscoreMat2(i,:) = binZscores2;
    FrMat2(i,:) = binFR2;
end

tBins = ((1:nBins)*binSize) - stimT - (binSize/2);

figure; hold on;
plot([-0.5, 1.5], [0, 0], 'k--')
errorbar(tBins, mean(zscoreMat1, 2), std(zscoreMat1, [], 2)/sqrt(nTrials), 'o')
errorbar(tBins, mean(zscoreMat2, 2), std(zscoreMat2, [], 2)/sqrt(nTrials), 'o')
legend({'0 z-score', 'Cell 1', 'Cell 2'})
xlabel('Time (s)'); ylabel('z-score')
ylim([-5, 55])

%% 6c

figure; hold on;
plot([-0.5, 1.5], [0, 0], 'k--')
errorbar(tBins, mean(FrMat1, 2), std(FrMat1, [], 2)/sqrt(nTrials), 'o')
errorbar(tBins, mean(FrMat2, 2), std(FrMat2, [], 2)/sqrt(nTrials), 'o')
legend({'0 Hz', 'Cell 1', 'Cell 2'})
xlabel('Time (s)'); ylabel('Firing rate (Hz)')
ylim([-5, 50])

% Both cells have similar firing rates in the post-stimulus presentation
% period, but because cell 2 has a much lower baseline firing rate it has a
% much higher z-score. 
%
% The problem in interpretation is that we need to determine whether
% absolute response matters in our experiment or whether change relative to 
% baseline is what matters. 
% 
% The difference in baseline firing may be unimportant for our question of
% interest, in that case the z-score comparison would not be meaningful. 



