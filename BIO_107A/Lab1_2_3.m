%Lab 1.2
% -------

A = [ 4 8 6 5 7 3 0 9 12];

%Accessing the elements of a variable
%We use parenthesis () to access the elements of an array or matrix variable:

myvariable = 5  % fine
my_variable = 5 % fine, underscores are fine
%5th_variable = 5 % does this work?
%my+var = 5  % does this work? 
%my var = 5  % does this work?

A(1)  % the first element
A(end) % the last element
A(5) % the 5th element

A(1,2)  % the 2nd column of the 1st row
A(1,end) % the last column of the 1st row
A(1,:) % the entire 1st row; the colon (:) inside parentheses means "everything"

%Functions in Matlab

mean_a = mean(A)

median_a = median(A)

percentile_a_15 = prctile(A,15)

%Definition
%[output1, output2, ...] = function_name(input1, input2, ...)
help mean
lookfor percentile 

% generates "bell curve" data with mean 1 and standard deviation 1
% multiplying the output by 100 changes the mean to 100 and
% standard deviation to 100
mydata1 = 100*generate_random_data(10,'normal',1,1); 

% generates "bell curve" data with mean 0 and standard deviation 1
mydata2 = 100*generate_random_data(10,'normal',0,1); 

mydata1
mydata2
mean(mydata1)
mean(mydata2)

% Figures and Plotting (simple)

f = figure;
plot(1,mydata1,'og')
hold on
plot(2,mydata2,'sb');

myax = axis
axis([0 3 -300 500])
xlabel('Class number');
ylabel('Change in student happiness');
title('Change in student happiness in Class 1 vs. Class 2');

% Character strings
mystring = 'This is a test'
mystring(1)
mystring(end)
mystring(5)

% generates "bell curve" data with mean 10 and standard deviation 1
mydata3 = 100*generate_random_data(10,'normal',10,1); 
% generates "bell curve" data with mean 0 and standard deviation 1
mydata4 = 100*generate_random_data(10,'normal',0,1); 

f2 = figure;
plot(1,mydata3,'og')
hold on
plot(2,mydata4,'sb');
axis([0 3 -300 1500])
xlabel('Class number');
ylabel('Change in student happiness');
title('Change in student happiness in Class 1 vs. Class 2');

figure(f) % brings the 1st figure to the front; any drawing would now happen here
title(['This is the first figure']);
figure(f2); % brings the second figure to the front
title(['This is the second figure']);

figure(1)
title(['This is still the first figure']);

%Lab 1.3
% -------


% generates "bell curve" data with mean 1 and standard deviation 1
mydata1 = 100*generate_random_data(10,'normal',1,1); 
% generates "bell curve" data with mean 0 and standard deviation 1
mydata2 = 100*generate_random_data(10,'normal',0,1); 

mn1 = mean(mydata1);
mn2 = mean(mydata2);


figure;
bar(1,mn1,'w');
hold on; % leave the current plot on the figure
bar(2,mn2,'w');
plot(1,mydata1,'og')
plot(2,mydata2,'sb')
axis([0 3 -300 500])
xlabel('Class number');
ylabel('Change in student happiness');
title('Change in student happiness in Class 1 vs. Class 2');

% Histogram
figure;
hist(mydata1,10); % creates a histogram with 10 bins
xlabel('Change in happiness after Class 1');
ylabel('Number of occurrences');

mydata1000 = 100*generate_random_data(1000,'normal',1,1);
figure;
hist(mydata1000,100); % creates a histogram with 100 bins
xlabel('Change in happiness after Class 1');
ylabel('Number of occurrences');


% Bins
figure;
hist(mydata1,100); % creates a histogram with 100 bins
xlabel('Change in happiness after Class 1');
ylabel('Number of occurrences');

%Histc
bin_edges = [ 0 1 2 3 4 5 ];
data = [ 0.5 0.5 1.5 1.5 3.5 3.5 3.7 4 4 ];
N = histc(data, bin_edges); % get the counts for each bin
N,  % look at the counts
bin_centers = bin_edges + 0.5; 
figure;
bar(bin_centers,N);
xlabel('Value');
ylabel('Number of occurrences');

% Colon operator
A = [ 1 2 3 4 5 ];
A = 1:5;

B1=5:-1:1 ;
B2=0:0.1:1;

C1 = [ 0 1 2 3 7 8 9 10];
C2 = [ 0 1 2 3 4 5 6 7];

C1(1);
C1(2);
(C1(1)+C1(2))/2;

C1(1:end-1);
C1(2:end);
(C1(1:end-1)+C1(2:end))/2;
(C2(1:end-1)+C2(2:end))/2;

% Script M files
bin_edges = [ -1000:100:1000 ];
N = histc(mydata1000, bin_edges); % get the counts for each bin
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2; 
  % histc actually returns an extra bin in N containing values that
  % are exactly bin_edges(end); we don't want this so let's drop this bin
N = N(1:end-1);
figure;
bar(bin_centers,N);
xlabel('Change in happiness after Class 1');
ylabel('Number of occurrences');


% Scripts and Functions - histbins.m

%Cumulative Histogram

Y = [ 0 0:0.1:100 100]; % we'll compute this in 0.1% steps
X = [min(mydata1)-1 prctile(mydata1, 0:0.1:100) max(mydata1)+1];
figure;
plot(X,Y,'k-'); % plot a black line
ylabel('Fraction of participants');
xlabel('Change in happiness after Class 1');







