
% 1
load ps4_problem1.mat
if ~lillietest(y1) | ~lillietest(y2)
p = ranksum(y1,y2)
end
x_n = [ones(91,1) x'];
[b,bint,r,rint,stats] = regress(y1',x_n);
stats(1)


%a. means not significantly different, look at correlation
[c1 p1] = corrcoef(x,y1);
[c2 p2] = corrcoef(x,y2);
% both correlations significant, reasonable to do
% regression analysis
x = [x ; ones(size(x))];
[b1 bint1 r rint stats1] = regress(y1', x', .05);
[b2 bint2 r rint stats2] = regress(y2', x', .05);

b1, b2
bint1, bint2
% the slopes are significantly different, but the intercepts are not. Thus,
% the two groups are different


%2
close all
clear
load ps4_problem2.mat
x = repmat([1:40],10, 1)';
g1 = g1';
g2 = g2';
figure(1);
plot(x(:), g1(:), 'ro', x(:), g2(:),'bo');
legend('g1', 'g2')
xlabel('Trial number')
ylabel('Escape latencies')

mg1 = mean(g1,1);
mg2 = mean(g2,1);
err1 = std(g1)/sqrt(10);
err2 = std(g2)/sqrt(10);
% figure;
% errorbar(1:40, mg1, err1, 'ro'); hold on
% errorbar(1:40, mg2, err2, 'bo');
% legend('g1', 'g2')
% xlabel('Trial number')
% ylabel('Escape latencies')


%2c
% %c. anova
% y = [];
% group = [];
% time = [];
% for i = 1:10
%     y = [y g1(i,:)];
%     group = [group ones(size(g1(i,:)))];
%     time = [time 1:40];
%     y = [y g2(i,:)];
%     group = [group 2*ones(size(g1(i,:)))];
%     time = [time 1:40];
% end
% % group = [ones(size(g1)); 2*ones(size(g2))];
% % time = repmat(1:40, 2*size(g1,1),1);
% % y = reshape(y,[],1);
% % group =reshape(group,[],1);
% % time = reshape(time,[],1);
% [p atab stats] = anovan(y, {time group}, 3, 3, strvcat('group', 'time'));
% [comp m h] = multcompare(stats, .05, 'off','tukey-kramer',[] ,[1 2]);
% % the interesting differences are the groups with the same time, which are
% % group indeces separated by 40
% sig = find((sign(comp(:,3)) == sign(comp(:,5))) & (comp(:,1) == comp(:,2) - 40));
% % sig is empty, so none of the times are different




% 2d. Linear Regression
load ps4_problem2.mat
x = repmat([1:40], [10,1]);
x = reshape(x,[],1);
x = [ones(size(x)) x];
g1= reshape(g1,[],1);
g2 = reshape(g2,[],1);
[b1,bint1,r1,rint1,stats1] = regress(g1,x);
[b2,bint2,r2,rint2,stats2] = regress(g2,x);
y1 = b1(1) + b1(2)*[0:45];
y1_ci1 = bint1(1,1) + bint1(2,1)*[0:45];
y1_ci2 = bint1(1,2) + bint1(2,2)*[0:45];
y2 = b2(1) + b2(2)*[0:45];
y2_ci1 = bint2(1,1) + bint2(2,1)*[0:45];
y2_ci2 = bint2(1,2) + bint2(2,2)*[0:45];

figure(1); hold on
plot([0:45], y1, 'r');
fill( [0:45 45:-1:0], [y1_ci1 fliplr(y1_ci2)],'r'); alpha(.2)
plot([0:45], y2, 'b');
fill( [0:45 45:-1:0], [y2_ci1 fliplr(y2_ci2)],'b'); alpha(.2)

legend('g1', 'g2')
xlabel('Trial number')
ylabel('Escape latencies')



% 2e. Non-Linear Regression
load ps4_problem2.mat

% x = [];
% for i = 1:40
% x = [x ; ones(10,1) * i ones(10,1)];
% end
% % % Non-linear regression, y = Aexp(-Bx)
% %function [y] = expfn(beta, x)
% %y = beta(1) * exp(-beta(2) * x);
% %end
% x = x(:,1);
% g1row = reshape(g1, 400, 1);
% g2row = reshape(g2, 400, 1);
% [beta1 r j] = nlinfit(x, g1row, @expfn, [50 1]);
% bint1 = nlparci(beta1, r, j);
% [beta2 r j] = nlinfit(x, g2row, @expfn, [50 1]);
% bint2 = nlparci(beta2, r, j);



x = repmat([1:40], [10,1]);
x = reshape(x,[],1);
g1= reshape(g1,[],1);
g2 = reshape(g2,[],1);
[b1,rnl1,j1,covb1,mse1] = nlinfit(x,g1,@expfun,[50 1/50]);
ci1 = nlparci(b1,rnl1,'jacobian',j1);
y1_ci1 = expfun(ci1(:,1), x);
y1_ci2 = expfun(ci1(:,2), x);
[b2,rnl2,j2,covb2,mse2] = nlinfit(x,g2,@expfun,[50 1/50]);
ci2 = nlparci(b2,rnl2,'jacobian',j2);
y2_ci1 = expfun(ci2(:,1), x);
y2_ci2 = expfun(ci2(:,2), x);


b1, b2

% overlapping intercepts but different time constants
% plot data

load ps4_problem2.mat
mg1 = mean(g1,1);
mg2 = mean(g2,1);
err1 = std(g1)/sqrt(10);
err2 = std(g2)/sqrt(10);

figure;
errorbar(1:40, mg1, err1, 'ro'); hold on
errorbar(1:40, mg2, err2, 'bo');
legend('g1', 'g2')
xlabel('Trial number')
ylabel('Escape latencies')

hold on
plot(1:40, expfun(b1,1:40),'r')
plot(x,y1_ci1, 'r--')
plot(x,y1_ci2, 'r--')
% fill( [x fliplr(x)], [y1_ci1 fliplr(y1_ci2)],'r'); alpha(.2)
plot(1:40, expfun(b2,1:40),'b')
plot(x,y2_ci1, 'b--')
plot(x,y2_ci2, 'b--')
% fill( [x fliplr(x)], [y2_ci1 fliplr(y2_ci2)],'b'); alpha(.2)


%Looking at the sum square error, the nonlinear fit is better
sum(rnl1.^2)
sum(rnl2.^2)
sum(r1.^2)
sum(r2.^2)


%3.
clear
load ps4_problem3.mat
edges = [0:0.01:0.300];
spk1 = [];

spk2 = [];
for iTrial = 1:size(spiketimes1,1)
n = histc(spiketimes1{iTrial,1},edges);
spk1 = [spk1; n];
n = histc(spiketimes2{iTrial,1},edges);
spk2 = [spk2; n];
end
mspk1 = mean(spk1);
sepk1 = std(spk1)/sqrt(50);
mspk2 = mean(spk2);
sepk2 = std(spk2)/sqrt(50);
figure
errorbar(edges,mspk1,sepk1,'Color','r')
hold on
errorbar(edges,mspk2,sepk2,'Color','b')
xlabel('Time (ms)')
ylabel('Mean Spike Count')

%3b
spk1rate = sum(spk1,2)/0.300;
spk2rate = sum(spk2,2)/0.300;
totalSpks = [spk1rate; spk1rate];
y = [zeros(49,1); ones(49,1)];
y = [y ones(98,1)];
[b, dev, stats] = glmfit(totalSpks, y, 'binomial');
stats(1).p

%3c
spk1rate = sum(spk1,1);
spk2rate = sum(spk2,1);

