function AnovaTable = BetweenWithinAnova(data,between,within,subjects,sstype,names,printon,alpha)

%   Between- and Within- Subject (Repeated Measures) Variables Analysis of Variance Test.
%   Used to analyze designs that have a combination of between and within-subject variables.
%   Uses either Type I/II or III sum of squares.
%   
%   Syntax:     function AnovaTable =
%                 BetweenWithinAnova(data,between,within,subjects,names,printon,alpha)
%
%      
% Inputs:
%     data        - dependant variable (numeric) data vector (ndata X 1)
%     between     - between subject grouping variable (numeric) (ndata X 1)
%     within      - within subject grouping variable (numeric) (ndata X 1)
%     subjects    - subject grouping variable (numeric) (ndata X 1)   
%     sstype      - either 1 or 3. Use 3 for unbalanced experiments
%     names       - a cell of length 2 for names of factors between
%                   and within.
%     printon   - a flag to turn on and off printing of
%                   output.
%     alpha       - significance level (default = 0.05).
%
% Outputs:
%     AnovaTable  - Structure containing the anova results
%                 - prints complete Analysis of Variance Table to screen when
%                   printon = 'yes'
%
%
%**************************************************************************
%**************************************************************************
%    COMPLETELY REWORKED, INCLUDING ADDING ABILITY TO DO TYPE III
%    SUM OF SQUARES AND FULLY DEBUGGED AND CHECKED AGAINST SAS
%    AND SPSS OUTPUT BY:
%              A. Katrin Schenk (1,2) and Evan H. Goulding (3)
%              (1)UCSF Department of Neurology 
%              (2)Sloan-Swartz Center for Theoretical Neuroscience
%              (3)UCSF Department of Psychiatry
%              University of California, San Francisco
%              Rock Hall
%              1550 4th st.
%              San Francisco, CA 94501
%              schenk@phy.ucsf.edu, evan.goulding@ucsf.edu
%              May 20, 2007.
%***************************************************************************
%***************************************************************************
%    Based minimally on original code by
%               A. Trujillo-Ortiz, R. Hernandez-Walls and R.A. Trujillo-Perez
%               Universidad Autonoma de Baja California
%               Ensenada, Baja California
%               Mexico.


%    References:
%    Lane, M. D. http://davidmlane.com/hyperstat/questions/Chapter_14.html
%    Huck, S. W. (2000), Reading Statistics and Research. 3rd. ed. 
%             New-York:Allyn&Bacon/Longman Pub. Chapter 16.
%

format long
% Set defaults.
if nargin < 5,
  sstype = 1;
  names    = {'BetweenFactor','WithinFactor'};
  printon  = 'no';
  alpha = .05;
elseif nargin < 6,
  names    = {'BetweenFactor','WithinFactor'};
  printon  = 'no';
  alpha = .05;
elseif nargin < 7
  
  printon  = 'no';
  alpha = .05;
elseif nargin < 8
  
  alpha = .05;

end

if nargin < 4, 
   error('Requires at least four input arguments.');
   return;
end;

% Get data length
Nd = length(data);

% Get the number of levels for between (b) and within subject variables (w) and the number of subjects (s).
ub = unique(between);
uw = unique(within);
us = unique(subjects);

b = length(ub); % how many between subj levels
w = length(uw); % how many within subj levels
s = length(us); % how many subjects

% Print info about levels and subjects to screen.
if strcmp(printon,'yes')
    fprintf(['The number of ' names{1,1} ' (between subject) levels is:%2i\n\n'], b);
    fprintf(['The number of ' names{1,2} ' (within subject) levels is:%2i\n\n'], w);
    fprintf('The number of subjects is:%2i\n\n', s);
end

% total degrees of freedom (number of data points -1)
dfTO = Nd-1;  
% degrees of freedom for between subjects
dfB = b-1;  
% degrees of freedom for within subjects
dfW = w-1;  
% degrees of freedom of the between-error
dfEB = s-b; 
% degrees of freedom of the between x within
dfBW = dfB*dfW;  
% degrees of freedom of the within subjects error.
dfEBW = dfEB*dfW;  
  
% Get sum of squares for independent variable for between-subjects.
B = [];
for i = 1:b
  Xe = find(between==ub(i));
  tmp = sum(data(Xe))^2/length(data(Xe));
  B = [B,tmp];
end
% Get between subjects error SS
S = [];
for k = 1:s
  Xe = find(subjects==us(k));
  tmp = sum(data(Xe))^2/length(data(Xe));
  S = [S,tmp];
end;
BW = [];
for i = 1:b
  for j = 1:w
    Xe = find((between==ub(i)) & (within==uw(j)));
    tmp = sum(data(Xe))^2/length(data(Xe));
    BW = [BW,tmp];
  end
end

if sstype==1
  % correction term
  C = (sum(data))^2/Nd;
  
  % total sum of squares
  SSTO = sum(data.^2)-C;  
  
  % sum of squares for the between subject independent variable
  SSB = sum(B)-C;
 
  % Get sum of squares for the within-subjects independent variable
  W = [];
  for j = 1:w
    Xe = find(within==uw(j));
    tmp = sum(data(Xe))^2/length(data(Xe));
    W =[W,tmp];
  end
  
  % sum of squares for within subjects
  SSW = sum(W)-C;
  
  % sum of squares of the between x within
  SSBW = sum(BW)-sum(B)-sum(W)+C;  
  
elseif sstype==3
  % Make into column vectors if not already..
  data = data(:);
  subjects = subjects(:);
  within = within(:);
  between = between(:);

  % sort according to subject variable for ease
  [ss,ib]=unique([subjects within],'rows');
  swithin = within(ib);
  sdata = data(ib);
  sbetween= between(ib);
  ssubjects = subjects(ib);
  
  % Make x1
  x1 = ones(Nd,1);
  
  % Make x2, between subjects column, must sum to zero
  parity = 1;
  for bind = 1:b
    bis = find(sbetween==ub(bind));
    x2(bis,:) = parity;   
    parity = -parity;
  end 
  if sum(x2)>1e-6
    disp('bewteen groups column does not sum to zero, groups unbalanced!')
    
  end
  
  % Make within subjects columns, each w X w
  % block should sum to zero...
  rw = eye(w-1);
  rw(end+1,1:w-1) = -1*ones(1,w-1);
  xw = repmat(rw,s,1);
  
  if sum(sum(xw))>1e-6
    disp('within groups column does not sum to zero!')
    keyboard
  end
  
  % Make between X within subjects columns, each w X w
  % block should sum to zero...
  rb =repmat(x2,1,w-1); 
  xbw = rb.*xw;
  if sum(sum(xbw))>1e-6
    disp('within X between groups column does not sum to zero!')
    keyboard
  end
  X = [x1 x2 xw xbw];

  fullSS = sdata'*X*inv(X'*X)*X'*sdata;
  Xnobetween = X(:,1:1);
  Xnobetween = [Xnobetween X(:,3:end)];
  nobetweenSS = sdata'*Xnobetween*inv(Xnobetween'*Xnobetween)* ...
      Xnobetween'*sdata;
  Xnowithin = X(:,1:2);
  Xnowithin = [Xnowithin X(:,3+w-1:end)];
  nowithinSS = sdata'*Xnowithin*inv(Xnowithin'*Xnowithin)*Xnowithin'* ...
      sdata;
  Xnobw = X(:,1:2+w-1);
  nobwSS = sdata'*Xnobw*inv(Xnobw'*Xnobw)*Xnobw'*sdata;
  
  SSTO = fullSS;
  SSB = fullSS-nobetweenSS;
  SSW = fullSS-nowithinSS;
  SSBW = fullSS-nobwSS;

end

% mean square for between subjects
MSB = SSB/dfB;  

% mean square for within subjects
MSW = SSW/dfW;  

% sum of squares of the between-error
SSEB = sum(S)-sum(B);

% mean square for the between-error
MSEB = SSEB/dfEB;  

% mean square for the between x within
MSBW = SSBW/dfBW;  

% Sum of squares for the within subjects error.
SSEBW = sum(data.^2)-sum(BW)-sum(S)+sum(B);


% mean square for the within subjects error.
MSEBW = SSEBW/dfEBW;   

% F-statistics calculation
F1 = MSB/MSEB;
F2 = MSW/MSEBW;                             
F3 = MSBW/MSEBW;

% degrees of freedom re-definition
v1 = dfB;
v2 = dfEB;
v3 = dfW;
v4 = dfBW;
v5 = dfEBW;
v6 = dfTO;

%Probability associated with the F-statistics (P values).
[fcdf1,lf1]=myfcdf(F1,v1,v2);
if lf1==0
P1 = 1 -fcdf1;
else
  P1 = fcdf1;
end
[fcdf2,lf2]=myfcdf(F2,v3,v5);
if lf2==0
P2 = 1 - fcdf2;
else
  P2 = fcdf2;
end
[fcdf3,lf3]=myfcdf(F3,v4,v5);
if lf3==0
  P3 = 1 - fcdf3;
else
  P3 = fcdf3;
end

if strcmp(printon,'yes')
  disp('Between- and Within- Subject Variables Analysis of Variance Table.')
  fprintf('---------------------------------------------------------------------------\n');
  disp('SOV                    SS          df           MS             F        P')
  fprintf('---------------------------------------------------------------------------\n');
  fprintf([names{1,1} '              %11.3f%10i%15.3f%14.3f %f\n\n'],SSB,v1,MSB,F1,P1);
  fprintf(['Error(' names{1,1} ')       %11.3f%10i%15.3f\n\n'],SSEB,v2,MSEB);
  fprintf([names{1,2} '               %11.3f%10i%15.3f%14.3f %f\n\n'],SSW,v3,MSW,F2,P2);
  fprintf([names{1,1} ' x ' names{1,2} '        %11.3f%10i%15.3f%14.3f %f\n\n'],SSBW,v4,MSBW,F3,P3);
  fprintf(['Error(' names{1,1} ' x ' names{1,2} ') %11.3f%10i%15.3f\n\n'],SSEBW,v5,MSEBW);
  fprintf('Total             %11.3f%10i\n\n',SSTO,v6);
  fprintf('---------------------------------------------------------------------------\n');
end

AnovaTable.SStype = sstype;
AnovaTable.VarNames       = names;
AnovaTable.Between_SS         = SSB;
AnovaTable.Between_df         = v1;
AnovaTable.Between_MS         = MSB;
AnovaTable.Between_F          = F1;
AnovaTable.Between_P          = P1;
AnovaTable.Between_Err_SS     = SSEB;
AnovaTable.Between_Err_df     = v2;
AnovaTable.Between_Err_MS     = MSEB;
AnovaTable.Within_SS         = SSW;
AnovaTable.Within_df         = v3;
AnovaTable.Within_MS         = MSW;
AnovaTable.Within_F          = F2;
AnovaTable.Within_P          = P2;
AnovaTable.BetweenxWithin_SS     = SSBW;
AnovaTable.BetweenxWithin_df     = v4;
AnovaTable.BetweenxWithin_MS     = MSBW;
AnovaTable.BetweenxWithin_F      = F3;
AnovaTable.BetweenxWithin_P      = P3;
AnovaTable.BetweenxWithin_Err_SS = SSEBW;
AnovaTable.BetweenxWithin_Err_df = v5;
AnovaTable.BetweenxWithin_Err_MS = MSEBW;
AnovaTable.Total_SS       = SSTO;
AnovaTable.Total_df       = v6;

