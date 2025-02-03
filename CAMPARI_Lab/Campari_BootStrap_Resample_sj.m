% shuffled controls for CaMPARI2 experiments
% choose set number of cells at random from each slice
% shuffle slices and randomly split into two groups and measure difference
% between
% compare to known slice differences: Dep slices vs Control slices

function [pval_means, sig_ranksums] = Campari_BootStrap_Resample(CAMPRAT_MD)


% First Specify which animals to include (or which experimental datasets to include)
drop = input('Do you want to analyze the drop or rebound timepoint? Drop = 1, Rebound = 0   ');
if ~drop
    hm4di = input('Do you want to analyze rebound WITH BF ACh inhibition? YES = 1, NO = 0   ');
end
% These are all just different experimental groups of animals to choose
% from
if isfield(CAMPRAT_MD, 'expt_groups')
    if contains(CAMPRAT_MD.expt_groups, 'scn1')
        if drop 
            animals = [1, 2, 3, 4]; % scn1 mice imaged on MD3
        else
            animals = [6, 7, 9, 10]; % scn1 mice imaged on MD6
            %animals = [11, 12, 13]; % scn1 mice imaged on MOCK MD3
        end
    else
        if drop 
            animals = [1, 2, 5, 6, 12, 13]; % animals to analyze - imaged on MD3
        else
            if hm4di
                animals = [14, 15, 18, 19, 20, 21, 23]; % animals analyzed on MD5 after Hm4Di and CNO on MD3-4
                %animals = [14, 15, 18, 19, 20, 21]; % animals analyzed on MD5 after Hm4Di and CNO on MD3-4
            else
                animals = [7, 8, 9, 10, 11]; % animals to analyze - imaged on MD6
                %animals = [22, 25, 26, 27, 28]; % animals imaged on MD5, CNO on MD3+4 but NO DREADDS
                %animals = [32, 33, 34, 35, 36, 37]; % animals imaged with NO MD, but 2.5 days of CNO in drinking water
                %animals = [38, 41]; %control animals - no MD or CNO or anything, new recovery method for slicing
                %animals = [39, 40]; % control animals - no MD but 2.5 days CNO in drinking water, plus new recovery method for slicing
            end
        end
    end
else
    if drop 
        animals = [1, 2, 5, 6, 12, 13]; % animals to analyze - imaged on MD3
    else
        if hm4di
            %animals = [14, 15, 18, 19, 20, 21]; % animals analyzed on MD5 after Hm4Di and CNO on MD3-4
            animals = [14, 15, 18, 19, 20, 21, 23]; % animals analyzed on MD5 after Hm4Di and CNO on MD3-4
        else
            animals = [7, 8, 9, 10, 11]; % animals to analyze - imaged on MD6
            %animals = [22, 25, 26, 27, 28]; % animals imaged on MD5, CNO on MD3+4 but NO DREADDS
            %animals = [32, 33, 34, 35, 36, 37]; % animals imaged with NO MD, but 2.5 days of CNO in drinking water
            %animals = [38, 41]; %control animals - no MD or CNO or anything, new recovery method for slicing
            %animals = [39, 40]; % control animals - no MD but 2.5 days CNO in drinking water, plus new recovery method for slicing
        end
    end
end



% get data of interest - from experimental group of interest
DATA = CAMPRAT_MD.CELLS(ismember(CAMPRAT_MD.CELLS(:,3), animals), :);
for ii = 1:length(DATA(:,1)) % unique slice codes now
    DATA(ii,4) = str2double([num2str(DATA(ii,3)), num2str(DATA(ii,4))]);
end

% Get cell data by slice (variable to be shuffled)
all_slices = unique([DATA(:,4)]);
bySlice = struct;
dep_slices = [];
con_slices = [];
all_numCells = []; % number of cells in each slice



% SJ: Things to do
%------------------
% Don't subset cells in first go. 
% Get a Test Statistic for comparing dep and con cumulative distributions
% Max distance between two distributions is a good test statistic used in
% ks tests

% In resmapling, we will generate a distribition of this test statistic
% parameter for shuffled distributions

% So, shuffle slice identity, and 
% a) randomly draw same total number of dep and con cells as original. 
% b) Or choose a number of cells to draw. such as mean(N_depcells, N_concells)).
% Get max distance between two suffled distributions for each run.
% Compare if test statistic is signficantly higher than 95% of resampled max distance distribution 



dep_cells_rg_byslice=[]; con_cells_rg_byslice=[]; % Gather dep and con cells
N_dep_cells_perslice=[]; N_con_cells_perslice=[]; % Count N_cells per slice 

for ss = 1:length(all_slices)
    bySlice.(['slice' num2str(all_slices(ss))]) = DATA(DATA(:,4)==all_slices(ss),:);
    all_numCells = [all_numCells, length(bySlice.(['slice' num2str(all_slices(ss))])(:,1))];
    if bySlice.(['slice' num2str(all_slices(ss))])(1, 5) == 1
        dep_slices = [dep_slices, all_slices(ss)];
        N_dep_cells_perslice=[N_dep_cells_perslice,length(bySlice.(['slice' num2str(all_slices(ss))])(:,1))];
        dep_cells_rg_byslice=[dep_cells_rg_byslice;bySlice.(['slice' num2str(all_slices(ss))])(:,2)];
    else
        con_slices = [con_slices, all_slices(ss)];
        N_con_cells_perslice=[N_con_cells_perslice,length(bySlice.(['slice' num2str(all_slices(ss))])(:,1))];
        con_cells_rg_byslice=[con_cells_rg_byslice;bySlice.(['slice' num2str(all_slices(ss))])(:,2)];
    end
end

% Total cells in con and dep slices
N_con_cells = sum(N_con_cells_perslice);
N_dep_cells = sum(N_dep_cells_perslice);

% Normal hypothesis test
[h,p_kstest] = kstest2(con_cells_rg_byslice, dep_cells_rg_byslice)

% Normal plot
figure; hold on
c_color = [0,0,0]; d_color = [1,0,1];
[hc,stats_con] = cdfplot(con_cells_rg_byslice);
hc.Color = c_color;
hc.LineWidth = 3;
hc.DisplayName = ['Control', ' n = '  num2str(N_con_cells)];

[hd,stats_dep] = cdfplot(dep_cells_rg_byslice);
hd.Color = d_color;
hd.LineWidth = 3;
hd.DisplayName = ['Deprived', ' n = '  num2str(N_dep_cells)];

set(gca,'xscale','log');
xlim([10e-2, 10e0])


% CDF maximum distance plot - get max distance between distributions and
% plot
% -------------------------
[maxdiff, maxdiff_location, Xvalues, con_CDF, dep_CDF] = ...
    cumulative_hist_diff(con_cells_rg_byslice,dep_cells_rg_byslice);
figure;
plot(Xvalues,con_CDF,'b-');
hold on;
plot(Xvalues,dep_CDF,'r-');

% Plot maxdiff and location - maxdiff is the test statistic
% --------------------------------------
X_loc = [ Xvalues(maxdiff_location) Xvalues(maxdiff_location) ];
Y_loc = [ con_CDF(maxdiff_location) dep_CDF(maxdiff_location)];
plot(X_loc, Y_loc,'k-','linewidth',2);

% Other test statistics
%-------------------
meandiff = mean(con_cells_rg_byslice) - mean(dep_cells_rg_byslice);


% Resample over the control and deprived cell rg values
% -----------------------------------------


% 1) Use same number of dep and con cells as original
% -----------------------------------------
all_rg = [con_cells_rg_byslice; dep_cells_rg_byslice];
N_all_cells = N_con_cells+N_dep_cells;
nruns = 10000;
store_maxdist = []; store_meandist = [];
for i = 1:nruns
    r = randperm(N_all_cells);
    dummy_con_rg = all_rg(r(1:N_con_cells));
    dummy_dep_rg = all_rg(r(N_con_cells+1:end)); % should be N_dep_cells
    [dmaxdiff, dmaxdiff_location, dXvalues, dcon_CDF, ddep_CDF] = ...
    cumulative_hist_diff(dummy_con_rg,dummy_dep_rg);
    dmean_diff = mean(dummy_con_rg) - mean(dummy_dep_rg);

    store_maxdist = [store_maxdist,dmaxdiff];
    store_meandist = [store_meandist,dmean_diff ];

end

% Plot distribution of dummy max_distances
lower = min(store_maxdist); higher = max(store_maxdist);
[N_hist,Edges] = histcounts(store_maxdist);
figure; hold on;
bar(Edges(1:end-1), N_hist)
plot(Edges(1:end-1), N_hist,'k--','LineWidth',2)
plot(maxdiff*ones(1,max(N_hist)),[1:max(N_hist)],'k-','LineWidth',4);

resamp_pval1 = length(find(store_maxdist>=maxdiff))/length(store_maxdist)



% Plot distribution of dummy max_distances
lower = min(store_meandist); higher = max(store_meandist);
[N_hist,Edges] = histcounts(store_meandist);
figure; hold on;
bar(Edges(1:end-1), N_hist)
plot(Edges(1:end-1), N_hist,'k--','LineWidth',2)
plot(meandiff*ones(1,max(N_hist)),[1:max(N_hist)],'k-','LineWidth',4);

resamp_pval_mean1 = length(find(store_meandist>=meandiff))/length(store_meandist)




% 2) Use mean number of cells between con vs dep
% -----------------------------------------
N_mean_cells = mean([N_con_cells,N_dep_cells]); %968-total is 1936
nruns = 10000;
store_maxdist = [];
for i = 1:nruns
    r = randperm(N_all_cells);
    dummy_con_rg = all_rg(r(1:N_mean_cells));
    dummy_dep_rg = all_rg(r(N_mean_cells+1:end));
    [dmaxdiff, dmaxdiff_location, dXvalues, dcon_CDF, ddep_CDF] = ...
    cumulative_hist_diff(dummy_con_rg,dummy_dep_rg);
    store_maxdist = [store_maxdist,dmaxdiff];
end

% Plot distribution of dummy max_distances
lower = min(store_maxdist); higher = max(store_maxdist);
[N_hist,Edges] = histcounts(store_maxdist);
figure; hold on;
bar(Edges(1:end-1), N_hist)
plot(Edges(1:end-1), N_hist,'k--','LineWidth',2)
plot(maxdiff*ones(1,max(N_hist)),[1:max(N_hist)],'k-','LineWidth',4);

resamp_pval2 = length(find(store_maxdist>=maxdiff))/length(store_maxdist)


keyboard;







% ---------------------------------------------------------------------
% Juliet's code

% numCells = min(all_numCells); % number of cells to analyze from each slice, to control for differences
% 
% % plot actual deprived vs control?
% 
% % get cell r/g data from deprived and control slices
% dep_cells = zeros(1, numCells*length(dep_slices));
% con_cells = zeros(1, numCells*length(con_slices));
% for d = 1:length(dep_slices)
%     if d == 1
%         start_in = 1;
%     else
%         start_in = ((d-1)*numCells)+1;
%     end
%     dep_cells(start_in:start_in+numCells-1) = datasample(bySlice.(['slice' num2str(dep_slices(d))])(:, 2)', numCells, 'Replace', false);
% end
% for c = 1:length(con_slices)
%     if c == 1
%         start_in = 1;
%     else
%         start_in = ((c-1)*numCells)+1;
%     end
%     con_cells(start_in:start_in+numCells-1) = datasample(bySlice.(['slice' num2str(con_slices(c))])(:, 2)', numCells, 'Replace', false);
% end
% 
% disp ('pval for cellNum controlled - difference in deprived vs control groups:')
% %pval = ranksum(DATA(DATA(:, 5) == 1, 2), DATA(DATA(:, 5) == 0 , 2)) % redo this with only taking numCells from each slice
% pval = ranksum(dep_cells, con_cells) 
% 
% % simple bootstrapping
% runs = 1000; % 1000 times
% resamplenum = min(length(dep_slices),length(con_slices)); % use the minimal number of slices to match samples
% ranksum_resample = zeros(1, runs);
% dep_resample_m = zeros(1, runs);
% con_resample_m = zeros(1, runs);
% for i = 1:runs
%     dep_resample = datasample(all_slices,resamplenum);
%     con_resample = datasample(all_slices,resamplenum);
%     
%     dep_cells_resample = zeros(1, numCells*resamplenum);
%     con_cells_resample = zeros(1, numCells*resamplenum);
%     for dd = 1:length(resamplenum)
%         if dd == 1
%             start_in = 1;
%         else
%             start_in = ((dd-1)*numCells)+1;
%         end
%         dep_cells_resample(start_in:start_in+numCells-1) = datasample(bySlice.(['slice' num2str(dep_resample(dd))])(:, 2)', numCells, 'Replace', false);
%         con_cells_resample(start_in:start_in+numCells-1) = datasample(bySlice.(['slice' num2str(con_resample(dd))])(:, 2)', numCells, 'Replace', false);
%     end
%     
%     % ks test for dep cells vs con cells?
%     % wilcoxon for dep cells vs con cells
%     ranksum_resample(i) = ranksum(dep_cells_resample, con_cells_resample);
%     dep_resample_m(i) = mean(dep_cells_resample);
%     con_resample_m(i) = mean(con_cells_resample);
% end
% 
% pval_means = mean(dep_resample_m < con_resample_m); % bootstrappped p-value
% sig_ranksums = (sum(find(ranksum_resample < 0.05))/runs)*100; % percentage of reshuffled slices that were still significantly different (not directional)
% end
