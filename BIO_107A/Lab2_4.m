% Lab2.4
% Dealing with a lot of data

help tumorfit  % read the help
help tumorplot % read the help
cd data
tumor_data = load('tumor_example.txt','-ascii');
cd ..
[change, rate,a,b,c,d] = tumorfit(tumor_data,5);
figure;
tumorplot(tumor_data,a,b,c,d);

% Dealing with whole data sets, with results from individual experiments stored in different directories

cd data/tumor_study

ls placebo/exp001    % refers to exp001 relative to the current directory
ls placebo/exp001/tumor_data.txt
td = load('placebo/exp001/tumor_data.txt','-ascii');

% String variables
%myfilename = 'placebo/exp001/tumor_data.txt'
%td2 = load(myfilename,'-ascii'); 

cd ..
[change_d,rate_d,a_d,b_d,c_d,d_d] = analyze_tumors('tumor_study/drug');
[change_p,rate_p,a_p,b_p,c_p,d_p] = analyze_tumors('tumor_study/placebo');

[change_d,rate_d,a_d,b_d,c_d,d_d]=analyze_tumors_plot('tumor_study/drug',1);
[change_p,rate_p,a_p,b_p,c_p,d_p]=analyze_tumors_plot('tumor_study/placebo',1);

