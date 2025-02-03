
n=1000;
R_store=[];
len = length(spatial_corr_day1);

for i=1:n

    curr_spatial_corr = spatial_corr_day1(randperm(len));  % shuffle
    [R_curr] = corrcoef(curr_spatial_corr, SWR_cofiring_day1);
    R_store(i)=R_curr(1,2);
    
    %figure; hold on; plot(curr_spatial_corr, SWR_cofiring_day1,'g.');
    %keyboard;
    

end
mean(R_store)
conf = prctile(R_store, 95)

edge_val = [min(R_store):0.01:max(R_store)];

hxc = histc(R_store, edge_val);
figure; hold on; bar(edge_val, hxc);



hx=hist(R_store,[min(R_store):0.01:max(R_store)]);
figure; hold on; bar(hx);

