% Lab 2.2

cd data
armdata = load('armdata.txt','-ascii');
cd ..

% Fit a gaussian
figure;
plot(armdata(:,1),armdata(:,2),'k');
xlabel('Probe tap position (in mm)');
ylabel('Nerve impulses per second')

methods('fittype')

gauss = fittype('a+b*exp(-((x-c).^2)/((2*d^2)))');
fo = fitoptions(gauss),

fo.StartPoint = [0; 0.5; 0; 0.5]; % guess a=0, b=0.5, c=0; d=0.5
gauss = setoptions(gauss,fo); % install the options in the fittype
[nerve,nervegof] = fit(armdata(:,1),armdata(:,2),gauss)
hold on
h = plot(armdata(:,1),nerve(armdata(:,1)),'k--');

% Watch fit happen, as parameters changed
[nerve,nervegof,errorovertime]=watchfithappen(armdata(:,1),armdata(:,2),gauss,30);

% Plot error over time
figure(100);
plot(errorovertime);
xlabel('Fit iteration number');
ylabel('Sum of squared error between fit and data');

% Value of c near 6
gauss2 = gauss;
fo = fitoptions(gauss2);
fo.StartPoint = [0; 0.5; 6; 0.5];
gauss2 = setoptions(gauss2,fo);

[nerve2,nervegof2,errortime2]=watchfithappen(armdata(:,1),armdata(:,2),gauss2,30);

figure(100);
hold on;
plot(errortime2,'k-');
legend('Bad starting position','Good starting position');

% Local error minima vs global minimum error
% ----------------------
gauss3 = gauss;
fo = fitoptions(gauss3);
[maxvalue,location] = max(armdata(:,2)); % maxvalue will be guess for b
position_guess = armdata(location,1); % guess parameter c
width_guess = 1; % arbitrarily pick a width guess, parameter d
fo.StartPoint = [mean(armdata(:,2)); maxvalue; position_guess; width_guess];

fo,

fo.Lower = [-maxvalue; 0; min(armdata(:,1)); -Inf];
fo.Upper = [maxvalue; Inf; max(armdata(:,1)); Inf];
gauss3 = setoptions(gauss3,fo);

[nerve3,nervegof3,errortime3]=watchfithappen(armdata(:,1),armdata(:,2),gauss3,30);

figure(100);
plot(errortime3,'g-');
legend('Bad starting position','Good starting position','Great starting position');


%GIGO

% Situation 1a = function-data mismatch
% -------

x = 1:0.1:10;
y = sin(x);
s = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0],'Upper',[Inf max(x)],'StartPoint',[1 1])
f = fittype('a*(x-b)^2','options',s);
[c,gof]=fit(x',y',f)   % need to make sure x and y are columns
figure;
hold on;
plot(x,y,'k');
xlabel('x'); ylabel('y');
plot(x,c(x),'b--');


% Situation 1b
% --------

cd data
morearmdata = load('morearmdata.txt','-ascii');
cd ..
gauss4 = gauss;
fo = fitoptions(gauss4);
[maxvalue,location] = max(morearmdata(:,2)); % maxvalue will be guess for b
position_guess = morearmdata(location,1); % guess parameter c
width_guess = 1; % arbitrarily pick a width guess, parameter d
fo.StartPoint = [mean(morearmdata(:,2)); maxvalue; position_guess; width_guess];
fo.Lower = [-maxvalue; 0; min(morearmdata(:,1)); -Inf];
fo.Upper = [maxvalue; Inf; max(morearmdata(:,1)); Inf];
gauss4 = setoptions(gauss4,fo);
[nerve4,ngof4,et4]=watchfithappen(morearmdata(:,1),morearmdata(:,2),gauss4,30);

% Omit noisy data points
5*std(morearmdata(:,2))
maxvalue,

% Situation 2 - Functipon-data match, but not enough data points

armdata_low = armdata(1:200:end,:); % sample every 2mm
gauss5 = gauss;
fo = fitoptions(gauss5);
[maxvalue,location] = max(armdata_low(:,2)); % maxvalue will be guess for b
position_guess = armdata_low(location,1); % guess parameter c
width_guess = 1; % arbitrarily pick a width guess, parameter d
fo.StartPoint = [mean(armdata_low(:,2)); maxvalue; position_guess; width_guess];
fo.Lower = [-maxvalue; 0; min(armdata_low(:,1)); -Inf];
fo.Upper = [maxvalue; Inf; max(armdata_low(:,1)); Inf];
gauss5 = setoptions(gauss5,fo);
[nerve5,ngof5,et5]=watchfithappen(armdata_low(:,1),armdata_low(:,2),gauss5,30);
hold on
plot(armdata(:,1),armdata(:,2),'m');

% Number of parameters and model selection
% -----------------------------------------

load census

% Exponential fit
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[Inf,max(cdate)],...
               'Startpoint',[1 1]);
f = fittype('a*(x-b)^2','options',s);
[c,gof] = fit(cdate,pop,f)
figure;
plot(cdate,pop,'o'); xlabel('Year'); ylabel('Population, millions');
x_values = 1750:2000;
y_values = c(x_values); % this returns c evaluated with best a, b
hold on;
plot(x_values,y_values,'b--');

% 3rd order polynomial fit
s2 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0 0 0],...
               'Upper',[Inf,max(cdate),Inf,max(cdate)],...
               'Startpoint',[1 1 1 1]);
f2 = fittype('a*(x-b)^2+c*(x-d)^3','options',s2);
[c2,gof2] = fit(cdate,pop,f2)
x_values = 1750:2000;
y2_values = c2(x_values); % this returns c evaluated with best a, b
hold on;
plot(x_values,y2_values,'m--');


% Nested F-test. Degrees of freedom
% ------------------------
df_more = length(pop) - 4; % let df be number of data points - number parameters
df_less = length(pop) - 2; % same
SSE_moreP = gof2.sse;
SSE_lessP = gof.sse;
changeinerror = ((SSE_lessP - SSE_moreP)/SSE_moreP);
changeindegreesoffreedom = ((df_less-df_more)/df_more);
F = changeinerror / changeindegreesoffreedom;

% We can obtain a P value for the nested F test by evaluating the cumulative density 
%function for F, which is provided by the Matlab function fcdf (see help fcdf):

P_nestedF = 1-fcdf(F,df_less-df_more,df_more),

% Model Selection imp.

% Number of parameters and overfitting
% ------------------------------






