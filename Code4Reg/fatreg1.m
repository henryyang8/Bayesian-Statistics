clear all
close all

 disp('Linear Regression: Body Fat Example One Variable')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;

load  'C:\STAT\Reg\Regdat\fat.dat'
y = fat(:,2);   %brozek index -- y
x  = fat(:,12);   %abdomen   -- x
%
n = length(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2; %number of parameters (beta0, beta1)
% beta0 estimated by b0; beta1 estimated by b1
% Abdomen measurement is "x" and the response brox is "y".
% Sums of Squares
SXX = sum( (x - mean(x)).^2 )    %2.9185e+04
SYY = sum(  (y - mean(y)).^2 )   %1.5079e+04
SXY = sum( (x - mean(x)).* (y - mean(y)) )  %1.7070e+04

% 
% estimators of coefficients beta1 and beta0
b1 = SXY/SXX    %0.5849
b0 = mean(y) - b1 * mean(x)   %-35.1966
% predictive equation (regression equation)
yhat = b0 + b1 * x;
%residuals
res = y - yhat;
% ANOVA Identity
SST = sum(  (y - mean(y)).^2 ) %the same as SYY=1.5079e+04
SSR = sum(  (yhat - mean(y)).^2 )   %9.9841e+03
SSE = sum(  (y - yhat).^2 ) %which is sum(res.^2)=5.0949e+03
% forming F and test for adequacy of linear relationship
MSR = SSR/(p - 1)   %9.9841e+03
MSE = SSE/(n - p)  %should be sigma2hat = 20.3797
F = MSR/MSE     %489.9029
pvalue = 1-fcdf(F, p-1, n-p) %H_0: regression has beta1=0, no need for 
% linear fit, pvalue = 0
%  Other measures of goodness of fit
R2 = SSR/SST  %0.6621
R2adj = 1 - (n-1)/(n-p)* SSE/SST  %0.6608
s = sqrt(MSE)  %4.5144
%  Standard error of coefficient estimators
sb1 = s/sqrt(SXX)   %0.0264
sb0 = s * sqrt(1/n + (mean(x))^2/SXX )  % 2.4623

% are the coefficients equal to 0?
t1 = b1/sb1  %22.1337
pb1 = 2 * (1-tcdf(abs(t1),n-p) )  %0
t0 = b0/sb0  % -14.2942
pb0 = 2 * (1-tcdf(abs(t0),n-p) )  %0
% 
% slope H0: beta1=1/2 vs H1: beta1 > 1/2
tt1 = (b1 - 0.5)/sb1  %3.2125
ptt1 = 1-tcdf(tt1, n-p)  %7.4436e-04
% intercept  H0: beta0 = -30 vs H1: beta0 < -30
tt0 =(b0 - (-30))/sb0   %-2.1105 
ptt0 = tcdf(tt0, n-p)  %0.0179
% CIs
[b0 - tinv(1-0.05/2, n-p)*sb0, b0 + tinv(1-0.05/2, n-p)*sb0] %-40.0461  -30.3471
[b1 - tinv(1-0.05/2, n-p)*sb1, b1 + tinv(1-0.05/2, n-p)*sb1]  %0.5328    0.6369
%
newx = 110; %   New person with abdomen x=110
ypred = b0 + b1 * newx    % 29.1414
sym = s * sqrt(1/n + (mean(x) - newx)^2/SXX ) %s for  mean yhat =   0.5416
syp = s * sqrt(1 + 1/n + (mean(x) - newx)^2/SXX ) %s for prediction yhat = 4.5468
%intervals CI and PI
alpha = 0.05;
%mean response interval
lbym = ypred - tinv(1-alpha/2, n-p) * sym;
rbym = ypred + tinv(1-alpha/2, n-p) * sym;
% prediction interval
lbyp = ypred - tinv(1-alpha/2, n-p) * syp;
rbyp = ypred + tinv(1-alpha/2, n-p) * syp;
%print the intervals
[lbym rbym]     %28.0746   30.2081
[lbyp rbyp]     %20.1865   38.0962
%
figure(1)
% scatterplot of (x,y)'s and regression line
plot(x, y, 'ro')
hold on
plot(x, yhat, 'k-')
hold off
xlabel('x'); ylabel('y hat')

figure(2)
%residuals against x
stem(x, res) 
xlabel('x'); ylabel('residual')

figure(3)
%residuals against order
stem((1:n), res) 
xlabel('order number'); ylabel('residual')


figure(4)
%residuals against predicted value brozhat
plot(yhat, res,'o') 
xlabel('yhat'); ylabel('residual')

figure(5)
%residuals against predicted value brozhat
plot(y, res,'o') 
xlabel('y'); ylabel('residual')

figure(6)
%residuals against predicted value brozhat
qqplot( res) 


%% 
% Univariate Regression via a built in MATLAB  m-file


vecones = ones(size(y));
X = [vecones x];
[bx,bintx,rx,rintx,statsx] = regress(y, X);

bx       %beta estimators, here intercept b0 and slope b1
bintx  % CI's for b0 and b1, default 95%,  to get a 100*(1-ALPHA)% level CI, use regress(Y,X,ALPHA) .
rx        %residuals
rintx  %CIs for the residuals (useful for outlier detection), see figure 4 below.
statsx    % 1. the R-square statistic, 2. the F statistic, 3. p value
               % for the full model, and 4. an estimator of the error variance.

%  bx =
%  -35.1966
%     0.5849
%
%  bintx =
%   -40.0461  -30.3471
%     0.5328    0.6369
%              
%  rx = 
%    -2.0361
%    -6.4493
%     8.3847
%     . . . 
%    -1.7187
%     1.2472
%     2.4360
% 
% 
% rintx =
%   -10.9154    6.8432
%   -15.2900    2.3914
%    -0.4414   17.2108
%    . . . 
%    -10.5523    7.1149
%    -7.6310   10.1254
%    -6.4111   11.2830
% 
% 
% statsx =
%     0.6621  489.9029    0.0000   20.3797
%               
figure(4)
E = (rintx(1:n,2) - rintx(1:n,1))/2;
errorbar((1:n),rx(1:n), E)
hold on
plot((1:n),rx(1:n),'ro','markersize',5,'MarkerFaceColor','g') 
plot((1:n),1.96 * sqrt(statsx(4))*ones(1,n), 'r:','LineWidth',3)
plot((1:n),-1.96 * sqrt(statsx(4))*ones(1,n), 'r:','LineWidth',3)
axis tight


