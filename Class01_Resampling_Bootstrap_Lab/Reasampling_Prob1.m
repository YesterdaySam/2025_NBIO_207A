
[X,Y]=GetSample(10); % the data provided


%1a

 N=10;
 numsims=100000; % simulate many times
 d=nan(numsims,1);
% 
% for i=1:numsims
% Xboot=X(ceil(N*rand(N,1))); % the bootstrap / sampling from data w/ replacement
% Yboot=Y(ceil(N*rand(N,1)));
% d(i)=mean(Xboot)-mean(Yboot);
% end
% figure
% hist(d,40);
% title('1a');



%1b

% for i=1:numsims
% paired_indices=ceil(N*rand(N,1));
% Xboot=X(paired_indices);
% Yboot=Y(paired_indices);
% d(i)=mean(Xboot)-mean(Yboot);
% end
% figure
% hist(d,40);
% title('1b');
% 
% d_stdest=std(d) % ex. 0.15
% lower_CI=prctile(d,5) % ex. -0.95
% median=prctile(d,50) % ex. -0.70
% upper_CI=prctile(d,95) % ex. -0.4457


%2a

% d_data=mean(X)-mean(Y) % ex. -0.70
% dperm=nan(numsims,1);
% Z=[X;Y];
% for i=1:numsims
% Zperm=Z(randperm(length(Z)));
% dperm(i)=mean(Zperm(1:length(X)))-mean(Zperm((length(X)+1):length(Z)));
% end
% figure
% hist(dperm,40);
% title('2a');
% 
% lower_CI=prctile(dperm,5) % ex. -1.10
% median=prctile(dperm,50) % ex. -0.00
% upper_CI=prctile(dperm,95) %


%2b

d_data_2b=mean(X)-mean(Y) % ex. -0.70
dperm=nan(numsims,1);
Z=[X Y]; % here concatenate into columns to visualize "pairing"
Zperm=[];
for i=1:numsims
for j=1:10 % flipping columns randomly
if rand>0.5
Zperm(j,:)=Z(j,:);
else
Zperm(j,:)=fliplr(Z(j,:));
end
end
dperm(i)=mean(Zperm(:,1))-mean(Zperm(:,2));
end
figure
hist(dperm,40);
title('2b');
lower_CI=prctile(dperm,5) % ex. -0.46
median=prctile(dperm,50) % ex. 0.00
upper_CI=prctile(dperm,95) % ex. 0.44




