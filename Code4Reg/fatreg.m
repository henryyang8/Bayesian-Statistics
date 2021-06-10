clear all
close all


 disp('Multiple Linear Regression: Body Fat Example')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;
load('C:\STAT\Reg\Regdat\fat.dat')
%load('C:\Courses\bmestatu\fat.dat')
%load  'fat.dat'
     casen = fat(:,1); %case number
broz = fat(:,2); %dependent variable Y
     siri = fat(:,3); %function of densi
     densi = fat(:,4); %invasive (makes people wet!)
age = fat(:,5); %below are predictors
weight = fat(:,6);   
height = fat(:,7);   
adiposi = fat(:,8); 
    % ffwei is functionally dependent on density//brozek and siri indices.
neck = fat(:,10); 
chest = fat(:,11);
abdomen  = fat(:,12);   
hip = fat(:,13); 
thigh = fat(:,14);
knee = fat(:,15);   
ankle = fat(:,16);
biceps = fat(:,17);
forearm = fat(:,18);     
wrist = fat(:,19);
vecones = ones(size(broz)); % necessary for intercept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[betahat,Ibeta,res,Ires,stats] = regress(broz, [vecones abdomen hip ankle biceps])
 
%
 disp('=================================================================')
 disp('             p = 15, all 14 variables + intercept')
 disp('=================================================================')
 

X =[vecones age weight height adiposi  neck chest abdomen ...
    hip thigh knee ankle biceps forearm wrist];
% X is design matrix, n x p where n is the number of subjects and p is
% the number of parameters, or number of predictors + 1.
varnames = ['intercept=1 ' 'age=2 ' 'weight=3 ' ...
    'height=4 ' 'adiposi=5 ' 'neck=6 ' 'chest=7 ' 'abdomen=8 '...
    'hip=9' 'thigh=10 ' 'knee=11 ' 'ankle=12 ' 'biceps=13 ' 'forearm=14 ' 'wrist=15 '];
[n, p] = size(X)
% coefficients in multivariate linear regression
b = inv(X' * X) * X' * broz 
% hat matrix for the predictions and residuals: yhat=H y.
H = X * inv(X' * X) * X';
brozhat = H * broz;
res = broz - brozhat;
% ANOVA TABLE-----------------------------------
%SUMS OF SQUARES
 II = eye(n); %identity matrix
 JJ = ones(n); %matrix of ones, n x n
 SST = broz' * (II - 1/n * JJ) * broz
 SSR = broz' * (H - 1/n * JJ) * broz
 SSE =  broz' * (II - H) * broz
 % of course SST = SSR + SSE-------------
 MSR = SSR/(p-1)
 MSE = SSE/(n-p)
 F = MSR/MSE
 pval = 1-fcdf(F, p-1, n-p)
 %---------------------------------------
 Rsq = 1 - SSE/SST
 Rsqadj = 1 - (n-1)/(n-p) * SSE/SST
 s = sqrt(MSE)
 %----------------------------
 covmatb = MSE * inv(X' * X); % covariances among b's
 sb = sqrt(diag(covmatb)); %standard deviations for each coefficient
 tb = b./sb; %t values for testing and confidence intervals: t_{n-p} cut-points.
 
 pvals = 2 * tcdf(-abs(tb), n-p); %p-vals for test that population regression
 % coefficients are equal to 0. Used for rough ellimination of variables.
 
 % take all predictors (no intercept)
 XX =[age weight height adiposi neck chest abdomen ...
    hip thigh knee ankle biceps forearm wrist];
 corr(XX); %correlation among predictors (multicolinearity)
 VIFs =diag(inv(corr(XX)))  ; %diagonal of inverse correlation matrix -> VIF's
 
 
disp('---------------------------------------------')
disp('     var#       t        pval        VIF     ') 
disp('---------------------------------------------')
 [(1:p)'  tb pvals [0; VIFs] ]
 varnames
 %varnames = ['intercept=1 ' 'age=2 ' 'weight=3 ' ...
 %  'height=4 ' 'adiposi=5 ' 'neck=6 ' 'chest=7 ' 'abdomen=8 '...
 %  'hip=9' 'thigh=10 ' 'knee=11 ' 'ankle=12 ' 'biceps=13 ' ...
 %  'forearm=14 ' 'wrist=15 '];

%% 
  
 disp('================================================')
 disp('             p = 5, 4 vars + intercept                ')
 disp('================================================')
 

 X2 =[ones(size(weight)) weight   abdomen forearm wrist ];
 [n, p2] = size(X2)
 b2 = inv(X2' * X2) * X2' * broz
 H2 = X2 * inv(X2' * X2) * X2';
 brozhat2 = H2 * broz;
 res2 = broz - brozhat2;
 SST = broz' * (II - 1/n * JJ) * broz;
 SSR2 = broz' * (H2 - 1/n * JJ) * broz
 SSE2 =  broz' * (II - H2) * broz
 MSR2 = SSR2/(p2-1)
 MSE2 = SSE2/(n-p2)
 F2 = MSR2/MSE2
 pval2 = 1-fcdf(F2, p2-2, n-p2)
 %----------------------------
 Rsq2 = 1 - SSE2/SST
 Rsqadj2 = 1 - (n-1)/(n-p2) * SSE2/SST
 s2 = sqrt(MSE2)
 %----------------------------
 covmatb2 = MSE2 * inv(X2' * X2); % covariances among b's
 sb2 = sqrt(diag(covmatb2));
 tb2 = b2./sb2;
 pvals2 = 2 * tcdf(-abs(tb2), n-p2);
 varnames2 = ['weight=1  '   'abdomen=2 '  'forearm=3 ' 'wrist=4']
 
X2no1 =[weight   abdomen forearm wrist ];
corr(X2no1);
VIFs2 = diag(inv(corr(X2no1)));

disp('---------------------------------------------')
disp('     var#       t        pval        VIF     ') 
disp('---------------------------------------------')
 [(1:p2)'  tb2 pvals2 [0; VIFs2] ]

 %%
 X =[age weight height adiposi  neck chest abdomen ...
    hip thigh knee ankle biceps forearm wrist];
vecones = ones(size(age));
varnames = ['age=1 ' 'weight=2 ' ...
    'height=3 ' 'adiposi=4 ' 'neck=5 ' 'chest=6 ' 'abdomen=7 '...
    'hip=8 ' 'thigh=9  ' 'knee=10 ' 'ankle=11 ' 'biceps=12 ' 'forearm=13 ' 'wrist=14 '];
 [B,BINT,R,RINT,STATS] = regress(broz, [vecones X])
%%
stepwise(X, broz)

