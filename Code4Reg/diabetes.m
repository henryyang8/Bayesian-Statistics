 clear all
 close all force

 disp('Simple Linear Regression: Insulin and Acidity in Blood')
 lw = 2; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 9;
  %First column: Age
  %Second Column: Deficit in Alcaline in Blood [x]
  %Third Column: Log(C_peptide) level [y]
dia = [ 5.2,-8.1,4.8 ; ...
8.8,-16.1,4.1 ; ...
10.5,-.9,5.2 ; ...
10.6,-7.8,5.5 ; ...
10.4,-29,5 ; ...
1.8,-19.2,3.4 ; ...
12.7,-18.9,3.4 ; ...
15.6,-10.6,4.9 ; ...
5.8,-2.8,5.6 ; ...
1.9,-25,3.7 ; ...
2.2,-3.1,3.9 ; ...
4.8,-7.8,4.5 ; ...
7.9,-13.9,4.8 ; ...
5.2,-4.5,4.9 ; ...
.9,-11.6,3 ; ...
11.8,-2.1,4.6 ; ...
7.9,-2,4.8 ; ...
11.5,-9,5.5 ; ...
10.6,-11.2,4.5 ; ...
8.5,-.2,5.3 ; ...
11.1,-6.1,4.7 ; ...
12.8,-1,6.6 ; ...
11.3,-3.6,5.1 ; ...
1,-8.2,3.9 ; ...
14.5,-.5,5.7 ; ...
11.9,-2,5.1 ; ...
8.1,-1.6,5.2 ; ...
13.8,-11.9,3.7 ; ...
15.5,-.7,4.9 ; ...
9.8,-1.2,4.8 ; ...
11,-14.3,4.4 ; ...
12.4,-.8,5.2 ; ...
11.1,-16.8,5.1 ; ...
5.1,-5.1,4.6 ; ...
4.8,-9.5,3.9 ; ...
4.2,-17,5.1 ; ...
6.9,-3.3,5.1 ; ...
13.2,-.7,6 ; ...
9.9,-3.3,4.9 ; ...
12.5,-13.6,4.1 ; ...
13.2,-1.9,4.6 ; ...
8.9,-10,4.9 ; ...
10.8,-13.5,5.1];
%===================================
x = dia(:,2); 
y = dia(:,3); 
n = length(x);




p = 2; %number of parameters (beta0, beta1)
% beta0 estimated by b0; beta1 estimated by b1
% Sums of Squares
xbar=mean(x); ybar = mean(y);
SXX = (x-xbar)' * (x-xbar); % = sum( (x - mean(x)).^2 )
SYY = (y-ybar)' * (y-ybar); % = sum( (y - mean(y)).^2 )
SXY = (x-xbar)' * (y-ybar); % = sum( (x - mean(x)).* (y - mean(y)) )

% estimators of coefficients beta1 and beta0
b1 = SXY/SXX  
    %b1 = 0.0494
b0 = mean(y) - b1 * mean(x)
    %b0 = 5.1494
% predictions
yhat = b0 + b1 * x;
%residuals
res = y - yhat;
% ANOVA Identity
SST = sum(  (y - mean(y)).^2 ) %SST = 21.8070, which is SYY 
SSR = sum(  (yhat - mean(y)).^2 ) %SSR = 5.2079
SSE = sum(  (y - yhat).^2 ) %SSE = 16.5990, which is sum(res.^2)
% forming F and test of adequacy of linear regression
MSR = SSR/(p - 1)  %MSR = 5.2079
MSE = SSE/(n - p)  %MSE = 0.4049,   sigma2hat
F = MSR/MSE  %F=12.8637
pvalue = 1-fcdf(F, p-1, n-p)   %pvalue = 8.8412e-004
%H_0: regression has beta1=0, that is, no need for linear fit
%  Other measures of goodness of fit
R2 = SSR/SST     %R2 = 0.2388
R2adj = 1 - (n-1)/(n-p)* SSE/SST  %R2adj = 0.2203
s = sqrt(MSE)   %s = 0.6363

%  Standard error of coefficient estimators
sb1 = s/sqrt(SXX)   %sb1 = 0.0138
sb0 = s * sqrt(1/n + (mean(y))^2/SXX )  %sb0 = 0.1170
% are the coefficients equal to 0?
t1 = b1/sb1   % t1 = 3.5866
pb1 = 2 * (1-tcdf(abs(t1),n-p) )  %pb1 = 8.8412e-004
t0 = b0/sb0  % t0 = 44.0013
pb0 = 2 * (1-tcdf(abs(t0),n-p) )  %pb1 =  0

% predicting y for the new observation x, CI and PI
newx = -22; %alcaline defficiency -22
ypred = b0 + b1 * newx   %ypred =   4.0618
sym = s * sqrt(1/n + (mean(x) - newx)^2/SXX ) %sym = 0.2142, s for mean y 
syp = s * sqrt(1 + 1/n + (mean(x) - newx)^2/SXX ) 
                                          %syp = 0.6714, s for predicted y
%intervals CI and PI
alpha = 0.05;
%mean response interval
lbym = ypred - tinv(1-alpha/2, n-p) * sym;
rbym = ypred + tinv(1-alpha/2, n-p) * sym; 
% prediction interval
lbyp = ypred - tinv(1-alpha/2, n-p) * syp;
rbyp = ypred + tinv(1-alpha/2, n-p) * syp;
%print the intervals
[lbym rbym]   %[3.6293    4.4943]
[lbyp rbyp]   %[2.7059    5.4176]


figure(2)
% scatterplot of (x,y)'s and regression
plot(x, y,  'ko', 'MarkerSize',12,  'MarkerFaceColor','g')
hold on
plot(x, yhat, 'r-', 'linewidth', lw)
xlabel('Alcaline Defficiency')
ylabel('C-Peptide Level (log units)')
hold off
%print -depsc 'C:\Springer\Reg\Regeps\diabetes1.eps'

figure(3)
%residuals against x
stem(x, res) 

figure(4)
%residuals against y
stem(y, res) 

figure(5)
%residuals against yhat
stem(yhat, res) 

%=====================================

CoeffTable0 = dataset({['Const';'PredX'   ],'Predictor'},{[b0;b1],'Coefs'},{[sb0;sb1],'StdErr'}, ...
                     {[t0; t1],'tStat'},{[pb0; pb1],'pVals'})
                 
fprintf('\n');
fprintf('%10s','s','R2','adjR2');
fprintf('\n')
fprintf('%10.4f',s,R2,R2adj);
fprintf('\n')


fprintf('\n')
fprintf('Regression ANOVA Table');
fprintf('\n\n')

fprintf('%6s','Source');
fprintf('%10s','df','SS','MS','F','P');
fprintf('\n')

fprintf('%6s','Regr');
fprintf('%10.4f',p-1,SSR,MSR,F,pvalue);
fprintf('\n')

fprintf('%6s','Resid');
fprintf('%10.4f',n-p,SSE,MSE);
fprintf('\n')

fprintf('%6s','Total');
fprintf('%10.4f',n-1,SST);
fprintf('\n')

%%
%===========BUILT IN MATLAB FUNCTION====================
XX = [ones(size(x)) x];
[B,BINT,R,RINT,STATS] = regress(y, XX)
figure(1)
   rcoplot(R, RINT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XXX = [ones(size(x)) dia(:,1) x];
[B,BINT,R,RINT,STATS] = regress(y, XXX )
figure(11)
   rcoplot(R, RINT)
   statsr = regstats(y, dia(:,1:2) )
