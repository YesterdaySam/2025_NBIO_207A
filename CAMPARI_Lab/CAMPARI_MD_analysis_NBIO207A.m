

function CAMPARI_MD_analysis(CAMPARI)

data_by_animal = 0; % set this to one ONLY if you want to plot data from each animal in desired group separately; otherwise set to 0

anim_age = [5,  2,  1,  7,  8,  6,  9,  13, 10, 11, 12, 14, 15, 17, 18, 19, 20, 21, 23; ...
    29, 26, 26, 32, 32, 29, 29, 28, 29, 30, 28, 30, 30, 31, 28, 28, 29, 29, 29];


age_color = colormap(jet(7));
%c_color = [.5,0,.5]; % input in scalar!!!
%c_color = [0,0,1];
c_color = [0,0,0];
%d_color = [.8,.6,1]; % input in scalar!!!
%d_color = [.6,.8,1];
d_color = [1,0,1];

drop = input('Do you want to analyze the drop or rebound timepoint? Drop = 1, Rebound = 0   ');
if ~drop
    hm4di = input('Do you want to analyze rebound WITH BF ACh inhibition? YES = 1, NO = 0   ');
end

if drop
    animals = [1, 2, 5, 6, 12, 13]; % animals to analyze - imaged on MD3 (EARLY MD/DROP time point)
else
    if hm4di
        animals = [14, 15, 18, 19, 20, 21, 23]; % animals analyzed on MD5 (LATE MD) after Hm4Di and CNO on MD3-4
    else
        animals = [7, 8, 9, 10, 11]; % animals to analyze - imaged on MD6 (LATE MD/REBOUND time point)
    end
end

% this section is just for some extra analysis for the hm4di experimental
% group
animals_CurveDiffs = NaN(1,length(animals));
animals_MedianDiffs = NaN(1, length(animals));
Hm4Di_score = [14, 15, 18, 19, 20, 21, 23; ...
    4416339802, 6264129831, 3210807479, 2184855441, 2941503892, 2190614446, 0; ... % average of whole image raw integrated densities
    761091392, 1079530903, 553335577.2, 376527790.8, 506925053.9, 377520270.8, 0; ... % average of whole image integrated density
    4223.6348, 5986.564043, 3069.540125, 2087.809143, 2811.61125, 2094.8145, 0; ... % average of whole image mean greyscale values
    1320.1, 2859, 999.9, 491.4, 992.7, 582.9, 753.3; ... % average of area of hm4di-positive projections in V1 images using Sydney's analysis with DEFiNE macro
    0.839, 1.084, 0.486, 0.479, 0.478, 0.749, 0.469; ... % total ratio of RFP/CHAT in BF of analyzed slices (from Sydney)
    0.502, 0.971, 0.636, 0.140, 0.757, 0.793, 0.671; ... % ratio of HDB only RFP/Chat in BF of analyzed slices (from Sydney)
    1.068, 1.627, 0.406, 0.531, 0.134, 0.235, 0.231; ... % ratio of VDB only RFP/Chat in BF of analyzed slices (from Sydney)
    1.52, 0.912, 0.259, 1.541, 0.405, 1.052, 0.394]; % ratio of above HDB (MCPO etc) only RFP/Chat in BF of analyzed slices (from Sydney)


rg_cells = CAMPARI.CELLS(ismember(CAMPARI.CELLS(:,3), animals), 2);
rg_deprived = CAMPARI.CELLS(ismember(CAMPARI.CELLS(:,3), animals), 5);
cellSizes = CAMPARI.CELLS(ismember(CAMPARI.CELLS(:,3), animals), 9);

rg_dep = rg_cells(rg_deprived==1);
rg_con = rg_cells(rg_deprived==0);
[h,p_kstest] = kstest2(rg_dep, rg_con)

% plot cumulative distributions of cells from each slice for each
% animal, to see how dep vs control slice vary within animals
for rr = 1:length(animals)
    %figure(200+rr)
    %hold on
    this_anim = CAMPARI.CELLS(CAMPARI.CELLS(:,3) == animals(rr), :);
    these_slices = unique(this_anim(:,4));
    dep_slices = [];
    con_slices = [];
    for ss = 1:length(these_slices)
        these_rg = this_anim(this_anim(:,4) == these_slices(ss), 2);
        %figure(200+rr)
        %[hi,stats_i] = cdfplot(these_rg);
        %hi.LineWidth = 2;
        if these_slices(ss) >= 20 % deprived slices always start with a 2, control slices always start with a 1
            dep_slices = [dep_slices; these_rg];
            %hi.Color = d_color;
            figure(100)
            hold on
            [hdd, stats_i] = cdfplot(these_rg);
            age = find(anim_age(1,:)==animals(rr));
            hdd.Color = age_color(anim_age(2,age)-25,:);
            hdd.LineWidth = 2;
            
        else
            con_slices = [con_slices; these_rg];
            %hi.Color = c_color;
            figure(99)
            hold on
            [hcc, stats_i] = cdfplot(these_rg);
            age = find(anim_age(1,:)==animals(rr));
            hcc.Color = age_color(anim_age(2,age)-25,:);
            hcc.LineWidth = 2;
        end
        %if rem(ss, 2)
        %    hi.LineStyle = '-';
        %else
        %    hi.LineStyle = ':';
        %end
        %hi.DisplayName = ['Slice Num: ', num2str(these_slices(ss)), ' n = '  num2str(length(these_rg))];
    end
    
    if drop
        myfigname = ['MD3 timepoint Animal ' num2str(animals(rr))];
    else
        myfigname = ['MD6 timepoint Animal ' num2str(animals(rr))];
    end
    %figure(200+rr)
    %title(myfigname)
    %legend('Location', 'southeast')
    %ylabel ('Proportion of Cells')
    %xlabel ('R/G Ratio')
    %set(gca,'xscale','log');
    %xlim([10e-2, 10e0])
    
    if data_by_animal
        % combine control slices and deprived slices per animal and plot by
        % animal
        if ~isempty(dep_slices) && ~isempty(con_slices)
            
            % This is plotting cdfs of each animal on its own plot
            figure(600+rr)
            
            hold on
            [hii,stats_ii] = cdfplot(dep_slices);
            hii.LineWidth = 2;
            hii.Color = d_color;
            hii.DisplayName = ['Deprived', ' n = '  num2str(length(dep_slices))];
            [hiii,stats_iii] = cdfplot(con_slices);
            hiii.LineWidth = 2;
            hiii.Color = c_color;
            hiii.DisplayName = ['Control', ' n = '  num2str(length(con_slices))];
            legend('Location', 'southeast')
            ylabel ('Proportion of Cells')
            xlabel ('R/G Ratio')
            my_anim_figname = [myfigname, ' BY ANIMAL'];
            title(my_anim_figname)
            set(gca,'xscale','log');
            xlim([10e-2, 10e0])
            
            
            % find difference in integral of deprived vs control cumulative distributions from 10th
            % percentile to 90th percentile (inverted - so with r/g ratios on y
            % axis and percentage of cells on x axis)
            x_dep = hii.XData; % dep data
            y_dep = hii.YData;
            x_con = hiii.XData; % control data
            y_con = hiii.YData;
            lims_dep = [find(y_dep>=.1,1), find(y_dep<=.9,1,'last')];
            lims_con = [find(y_con>=.1,1), find(y_con<=.9,1,'last')];
            i_con = trapz(y_con(lims_con(1):lims_con(2)), x_con(lims_con(1):lims_con(2)));
            i_dep = trapz(y_dep(lims_dep(1):lims_dep(2)), x_dep(lims_dep(1):lims_dep(2)));
            Animal_Difference = i_con-i_dep;
            fprintf('Difference between Animal %d control and depived curves: \n %d  \n', animals(rr), Animal_Difference)
            animals_CurveDiffs(rr) = Animal_Difference;
            animals_MedianDiffs(rr) = (median(x_dep)/median(x_con))*100;
        end
    end
end

if data_by_animal == 1
    
    if exist('hm4di', 'var') && hm4di == 1
        % plot correlation of difference between control and deprived
        % distributions and level of Hm4Di infected projections in V1m
        xs = NaN(1,length(animals));
        ys = NaN(1,length(animals));
        for aa = 1:length(animals)
            hm4di_anim = find(Hm4Di_score(1,:)==animals(aa));
            if ~isempty(hm4di_anim)
                xs(aa) = Hm4Di_score(7,hm4di_anim);
                %ys(aa) = animals_CurveDiffs(aa);
                ys(aa) = animals_MedianDiffs(aa);
            end
        end
        figure()
        scatter(xs,ys, 100, 'k', 'filled')
        %xlabel ('Amount of Hm4Di in V1')
        %ylabel ('Control minus Deprived Distributions')
        xlabel('Amount of Hm4di in HDB')
        ylabel('Median deprived R/G as % of Control')
        set (gca, 'box', 'off', 'fontsize', 18)
        disp ('correlation between Hm4di expression and difference btwn control and deprived hemisphere:')
        [rho, pval] = corr(xs', ys')
    end
    
end

% plot all control slices cumulative distributions and all deprived slices
% cumulative distributions together by AGE
figure(99)
title('Control Slices')
ylabel ('Proportion of Cells')
xlabel ('R/G Ratio')
hold on
set (gca, 'box', 'off', 'fontsize', 18)
colormap(jet(7))
%caxis([0,6])
cbar = colorbar('Ticks', 0:(1/6):1, 'TickLabels', {'26', '27', '28', '29', '30', '31', '32'});
cbar.Label.String = 'Age at Slicing';
set(gca,'xscale','log');
xlim([10e-2, 10e0])

figure(100)
title('Deprived Slices')
ylabel ('Proportion of Cells')
xlabel ('R/G Ratio')
hold on
set (gca, 'box', 'off', 'fontsize', 18)
colormap(jet(7))
%caxis([0,6])
cbar = colorbar('Ticks', 0:(1/6):1, 'TickLabels', {'26', '27', '28', '29', '30', '31', '32'});
cbar.Label.String = 'Age at Slicing';
set(gca,'xscale','log');
xlim([10e-2, 10e0])

% make cumulative distribution plots of all deprived and control cells for
% this experimental group
figure(900)
hold on

[hc,stats_con] = cdfplot(rg_con);
hc.Color = c_color;
%hc.Color = [.5, .5, .5];
hc.LineWidth = 3;
%hc.LineStyle = '--';
hc.DisplayName = ['Control', ' n = '  num2str(length(rg_con))];

[hd,stats_dep] = cdfplot(rg_dep);
hd.Color = d_color;
%hd.Color = [.8, .6, .8];
hd.LineWidth = 3;
%hd.LineStyle = '--';
hd.DisplayName = ['Deprived', ' n = '  num2str(length(rg_dep))];

set(gca,'xscale','log');
xlim([10e-2, 10e0])


if drop
    title('MD3');
else
    if hm4di
        title ('MD5, BF ACh Inhibition on MD3-4')
    else
        title ('MD6')
    end
end

legend ('Location', 'southeast')
ylabel ('Proportion of Cells')
xlabel ('R/G Ratio')
set (gca, 'box', 'off', 'fontsize', 18)

% now prepare for boxplots / violin plots

cells_deprived = [];
cells_deprived_L5 = [];
cells_deprived_L2 = [];

cells_control = [];
cells_control_L5 = [];
cells_control_L2 = [];

for i = 1:length(CAMPARI.CELLS(:,1)) % loop through every single cell with r/g ratio
    if ismember(CAMPARI.CELLS(i,3), animals) % only analyze selected animals's data
        
        if CAMPARI.CELLS(i,5) == 1
            cells_deprived = [cells_deprived; CAMPARI.CELLS(i,2)];
            if CAMPARI.CELLS(i,8) == 2
                cells_deprived_L2 = [cells_deprived_L2; CAMPARI.CELLS(i,2)];
            elseif CAMPARI.CELLS(i,8) ==  5
                cells_deprived_L5 = [cells_deprived_L5;  CAMPARI.CELLS(i,2)];
            end
            
        else
            cells_control = [cells_control; CAMPARI.CELLS(i,2)];
            if CAMPARI.CELLS(i,8) == 2
                cells_control_L2 = [cells_control_L2; CAMPARI.CELLS(i,2)];
            elseif CAMPARI.CELLS(i,8) ==  5
                cells_control_L5 = [cells_control_L5;  CAMPARI.CELLS(i,2)];
            end
        end
    end
    
end

mean_deprived = mean(cells_deprived);
median_deprived = median(cells_deprived);
mean_control = mean(cells_control);
median_control = median(cells_control);
%fprintf('MEAN deprived r/g = %d ; MEDIAN deprived r/g = %d \n', mean_deprived, median_deprived);
%fprintf('MEAN control r/g = %d ; MEDIAN control r/g = %d \n', mean_control, median_control);
fprintf('MEAN deprived r/g = %d percent of control ; MEDIAN deprived r/g = %d percent of control \n', (mean_deprived/mean_control)*100, (median_deprived/median_control)*100);

% Need groups to be consistent length for plotting
if length(cells_control) > length(cells_deprived)
    add = length(cells_control) - length(cells_deprived);
    cells_deprived(end+1:end+add) = NaN;
elseif length(cells_control) < length(cells_deprived)
    add = length(cells_deprived) - length(cells_control);
    cells_control(end+1:end+add) = NaN;
end

Lay_max = max([length(cells_control_L2), length(cells_control_L5), length(cells_deprived_L2), length(cells_deprived_L5)]);
if  length(cells_control_L2) < Lay_max
    add_lay = Lay_max - length(cells_control_L2);
    cells_control_L2(end+1:end+add_lay) = NaN;
end
if length(cells_control_L5) < Lay_max
    add_lay = Lay_max - length(cells_control_L5);
    cells_control_L5(end+1:end+add_lay) = NaN;
end
if length(cells_deprived_L2) < Lay_max
    add_lay = Lay_max - length(cells_deprived_L2);
    cells_deprived_L2(end+1:end+add_lay) = NaN;
end
if length(cells_deprived_L5) < Lay_max
    add_lay = Lay_max - length(cells_deprived_L5);
    cells_deprived_L5(end+1:end+add_lay) = NaN;
end

r_g_values = [cells_control, cells_deprived];

xlabels = {'Control Hem', 'Deprived Hem'};

r_g_values_byLayer = [cells_control_L2, cells_deprived_L2, cells_control_L5, cells_deprived_L5];
xlabels_layer  = {'Control L2/3', 'Deprived L2/3', 'Control L5/6', 'Deprived  L5/6'};

p_ranksum = ranksum(cells_control, cells_deprived)
figure()
% facecolors for Nick's paper: [.5, .5, .5; .3, .3, 1]
if strcmp(c_color, 'k') || isequal(c_color, [0,0,0])
    c_color_violin = [.5, .5, .5];
else
    c_color_violin = c_color;
end
[hv, lv] = violin(log10(r_g_values),'xlabel', xlabels,'facecolor',[c_color_violin; d_color],'edgecolor','k',...
    'bw',0.3, 'mc', '','medc', 'k');
set(lv,'visible','off')
ylabel('log10(R / G ratio)')
set (gca, 'box', 'off', 'fontsize', 25)
if drop
    title ('MD 3')
else
    if hm4di
        title('V1m Campari2 After 4 days MD and BF ACh inhibition')
    else
        title ('MD 6')
    end
end

% plot by putative layer
if ~isempty(cells_control_L2) && ~isempty (cells_control_L5)
    
    figure()
    [hv, lv] = violin(log10(r_g_values_byLayer),'xlabel',xlabels_layer,'facecolor',[c_color_violin; d_color; c_color_violin; d_color],'edgecolor','k',...
        'bw',0.3, 'mc', '','medc', 'k');
    set(lv,'visible','off')
    ylabel('log10(R / G ratio)')
    set (gca, 'box', 'off', 'fontsize', 25)
    if drop
        title ('MD 3')
    else
        if hm4di
            title('V1m Campari2 After 4 days MD and BF ACh inhibition')
        else
            title ('MD 6')
        end
    end
    
    p_wilcoxon_L2 = ranksum(cells_control_L2, cells_deprived_L2)
    p_wilcoxon_L5 = ranksum(cells_control_L5, cells_deprived_L5)
    
    % plot cumulative distributions of cells by putative layer
    figure()
    hold on
    [hc2,stats_con_2] = cdfplot(cells_control_L2);
    hc2.Color = c_color;
    hc2.LineWidth = 3;
    hc2.DisplayName = ['L2/3 Control', ' n = '  num2str(length(cells_control_L2))];
    [hd2,stats_dep_2] = cdfplot(cells_deprived_L2);
    hd2.Color = d_color;
    hd2.LineWidth = 3;
    hd2.DisplayName = ['L2/3 Deprived', ' n = '  num2str(length(cells_deprived_L2))];
    [hc5,stats_con_5] = cdfplot(cells_control_L5);
    hc5.Color = c_color;
    hc5.LineWidth = 3;
    hc5.LineStyle = ':';
    hc5.DisplayName = ['L5/6 Control', ' n = '  num2str(length(cells_control_L5))];
    [hd5,stats_dep_5] = cdfplot(cells_deprived_L5);
    hd5.Color = d_color;
    hd5.LineWidth = 3;
    hd5.LineStyle = ':';
    hd5.DisplayName = ['L5/6 Deprived', ' n = '  num2str(length(cells_deprived_L5))];
    
    if drop
        title ('Distribution of Cells on MD3 (by layer)')
    else
        if hm4di
            title ('Distribution of Cells on MD5 after BF ACh Inhibiton (by Layer)')
        else
            title ('Distribution of Cells on MD6 (by layer)')
        end
    end
    
    legend ('Location', 'southeast')
    ylabel ('Proportion of Cells')
    xlabel ('R/G Ratio')
    set (gca, 'box', 'off', 'fontsize', 18)
    set(gca,'xscale','log');
    xlim([10e-2, 10e0])
    
    mean_deprived_L2 = nanmean(cells_deprived_L2);
    median_deprived_L2 = nanmedian(cells_deprived_L2);
    std_deprived_L2 = nanstd(cells_deprived_L2)
    mean_control_L2 = nanmean(cells_control_L2);
    median_control_L2 = nanmedian(cells_control_L2);
    std_control_L2 = nanstd(cells_control_L2)
    mean_deprived_L5 = nanmean(cells_deprived_L5);
    median_deprived_L5 = nanmedian(cells_deprived_L5);
    std_deprived_L5 = nanstd(cells_deprived_L5)
    mean_control_L5 = nanmean(cells_control_L5);
    median_control_L5 = nanmedian(cells_control_L5);
    std_control_L5 = nanstd(cells_control_L5)
    %fprintf('MEAN deprived r/g = %d ; MEDIAN deprived r/g = %d \n', mean_deprived, median_deprived);
    %fprintf('MEAN control r/g = %d ; MEDIAN control r/g = %d \n', mean_control, median_control);
    fprintf('MEAN L2 deprived r/g = %d percent of control ; MEDIAN L2 deprived r/g = %d percent of control \n', (mean_deprived_L2/mean_control_L2)*100, (median_deprived_L2/median_control_L2)*100);
    fprintf('MEAN L5 deprived r/g = %d percent of control ; MEDIAN L5 deprived r/g = %d percent of control \n', (mean_deprived_L5/mean_control_L5)*100, (median_deprived_L5/median_control_L5)*100);
    
end


end