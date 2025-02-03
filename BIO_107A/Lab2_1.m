% Describing data with fits (models)

% Introduction

%Scatterplots
% -----------
cd data;
florida = load('manatee_powerboat.txt','-ascii');
cd ..

florida,  % see the data
figure;
plot(florida(:,2),florida(:,3),'ko'); % plot with black circles
xlabel('Powerboats registered (thousands)');
ylabel('Manatee deaths from powerboat collision');
box off; % turn off the ugly extra box

% Covariance
% ----------
x1 = rand(20,1);  % creates an array of random numbers distributed between 0 and 1
y1 = x1; % y1 is x1
x2 = rand(20,1);
y2 = rand(20,1);
x3 = rand(20,1);
y3 = -x3;

figure;
subplot(2,2,1);
plot(x1,y1,'ko'); 
ylabel('y1'); xlabel('x1');
subplot(2,2,2);
plot(x2,y2,'ko'); 
ylabel('y2'); xlabel('x2');
subplot(2,2,3);
plot(x3,y3,'ko'); 
ylabel('y3'); xlabel('x3');

cov(x2,y2)
var(x2)
var(y2)

C1 = cov(x1,y1)
C2 = cov(x2,y2)
C3 = cov(x3,y3)


% Correlation
r_1 = C1(1,2)/sqrt(C1(1,1)*C1(2,2))
r_2 = C2(1,2)/sqrt(C2(1,1)*C2(2,2))
r_3 = C3(1,2)/sqrt(C3(1,1)*C3(2,2))


corrcoef(x1,y1)
corrcoef(x2,y2)
corrcoef(x3,y3)
corrcoef(florida(:,2),florida(:,3))

[R,p] = corrcoef(florida(:,2),florida(:,3))


% Correlation and Causation

% Linear fits, aka regression
% --------------------------

x = florida(:,2);
y = florida(:,3);

figure;
plot(x,y, 'bo');
xlabel('Powerboats registered (thousands)');
ylabel('Manatee deaths from powerboat collision');
box off;

[R,p] = corrcoef(x,y)
C = cov(x,y);

m = R(1,2)*sqrt(C(2,2)/C(1,1))
b = mean(y) - m*mean(x)

X = [ 400 1100]; % chosen to match our X-axis of interest
Y = m * X + b
hold on
plot(X,Y,'k--');

Y_fit = m * x + b;
error = y-Y_fit;

[1 2 3].*[4 5 6];

squarederror = error.*error;% OR, error.^2 is el-by-el raising to 2nd power
meanerr = mean(squarederror)

m1 = m + 0.05; b1 = b + 1;
Y_fit1 = m1 * x + b1;
meanerr2 = mean((y-Y_fit1).^2);

plot([400 1100],[400 1100]*m1+b1,'g--');


%Non-linear models
%-----------------

cd data
load census  % loads variables 'cdate' and 'pop'
cd ..
figure;
plot(cdate,pop,'o')
ylabel('Population');
xlabel('Census year');
hold on;

s = fitoptions('Method','NonlinearLeastSquares',... % use squared error
               'Lower',[0,0],...   % lower bounds for [a b]
               'Upper',[Inf,max(cdate)],... % upper bounds for [a b]
               'Startpoint',[1 1]); % startpoint, search will start at a=1,b=1


f = fittype('a*(x-b)^2','options',s);
[c,gof] = fit(cdate,pop,f)

x_values = 1750:2000;
y_values = c(x_values); % this returns c evaluated with best a, b
hold on;
plot(x_values,y_values,'b--');


