%% ridge
openExample('stats/RidgeRegressionExample')

load acetylene
plotmatrix([x1 x2 x3])

X = [x1 x2 x3];
D = x2fx(X,'interaction');
D(:,1) = []; % No constant term
k = 0:1e-5:5e-3;
B = ridge(y,D,k);

figure
plot(k,B,'LineWidth',2)
ylim([-100 100])
grid on 
xlabel('Ridge Parameter') 
ylabel('Standardized Coefficient') 
title('Ridge Trace') 
legend('x1','x2','x3','x1x2','x1x3','x2x3')

%% lasso
openExample('stats/RemoveRedundantPredictorsWithCVFitsExample')

rng default % For reproducibility
X = randn(100,5);
weights = [0;2;0;-3;0]; % Only two nonzero coefficients
y = X*weights + randn(100,1)*0.1; % Small added noise
[B,FitInfo] = lasso(X,y,'CV',10,'PredictorNames',{'x1','x2','x3','x4','x5'});

idxLambdaMinMSE = FitInfo.IndexMinMSE;
minMSEModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinMSE)~=0)

idxLambda1SE = FitInfo.Index1SE;
sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)

lassoPlot(B,FitInfo,'PlotType','CV');
legend('show') % Show legend

%% lasso- using acetylene
openExample('stats/LassoPlotWithCrossValidatedFitsExample')

load acetylene

X = [x1 x2 x3];
D = x2fx(X,'interaction');
D(:,1) = []; % No constant term

rng default % For reproducibility 
[B,FitInfo] = lasso(D,y,'CV',10);

lassoPlot(B,FitInfo,'PlotType','CV');
legend('show') % Show legend

%% elastic net
openExample('stats/LassoAndElasticNetWithCrossValidationExample')
