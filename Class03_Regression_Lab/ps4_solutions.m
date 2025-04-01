%% 1
% x = stim contrast for each trial; y1 = passive condition FR; y2 = active
% condition FR

load('ps4_problem1.mat')

figure; plot(x,y1,'kx')
hold on; plot(x,y2,'bo')

[~,p_tt] = ttest2(y1,y2);
[p_rs] = ranksum(y1,y2);

% They are not significantly different mean FR responses in the 2
% conditions

[B1,BINT1,R1,RINT1,STATS1] = regress(y1',x',0.05)
[B2,BINT2,R2,RINT2,STATS2] = regress(y2',x',0.05)











