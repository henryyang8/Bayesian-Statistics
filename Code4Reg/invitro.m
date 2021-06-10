 disp('Linear Regression: In Virto Fertiloization Succes Rate by Age of Mother')
 lw = 2; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 7;
 age = (24:46)'; % AGE
 successrate =  [38.7  38.6  38.9  41.4  39.7  41.1  38.7  37.6  36.3  36.9  35.7  33.8 ...
       33.2  30.1  27.8  22.7  21.3  15.4  11.2  9.2  5.4  3.0  1.6]'; %SUCCESS RATE   
 x = (33:46)'; %  33 <= AGE <= 46
 y =  [ 36.9  35.7  33.8  33.2  30.1  27.8  22.7  21.3 ...
      15.4  11.2  9.2  5.4  3.0  1.6]'; %SUCCESS RATE for 33 <= AGE <= 46
   
vecones = ones(size(x));
n = length(vecones);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2; %number of parameters (beta0, beta1)
% beta0 estimated by b0; beta1 estimated by b1
%sums of squares
% Age is "x" and the response success rate is "y".
% Sums of Squares
xbar=mean(x); ybar = mean(y);
SXX = (x-xbar)' * (x-xbar); % = sum( (x - mean(x)).^2 )
SYY = (y-ybar)' * (y-ybar); % = sum(  (y - mean(y)).^2 )
SXY = (x-xbar)' * (y-ybar); % = sum( (x - mean(x)).* (y - mean(y)) )

% estimators of coefficients beta1 and beta0
b1 = SXY/SXX  
    %b1 = -3.0226
b0 = mean(y) - b1 * mean(x)
    %b0 = 139.9156
% predictions
yhat = b0 + b1 * x;
%residuals
res = y - yhat;
% ANOVA Identity
SST = sum(  (y - mean(y)).^2 ) %SST = 2.1208e+003, which is SYY 
SSR = sum(  (yhat - mean(y)).^2 ) %SSR = 2.0785e+003
SSE = sum(  (y - yhat).^2 ) %SSE = 42.2470, which is sum(res.^2)
% forming F and test of adequacy of linear regression
MSR = SSR/(p - 1)  %MSR = 2.0785e+003
MSE = SSE/(n - p)  %MSE = 3.5206,   sigma2hat
F = MSR/MSE  %F=590.3900
pvalue = 1-fcdf(F, p-1, n-p)   %pvalue = 1.4219e-011
%H_0: regression has beta1=0, that is, no need for linear fit
%  Other measures of goodness of fit
R2 = SSR/SST     %R2 = 0.9801
R2adj = 1 - (n-1)/(n-p)* SSE/SST  %R2adj = 0.9784
s = sqrt(MSE)   %s = 1.8763

%  Standard error of coefficient estimators
sb1 = s/sqrt(SXX)   %sb1 =  0.1244
sb0 = s * sqrt(1/n + (mean(y))^2/SXX )  %sb0 = 2.6016

% are the coefficients equal to 0?
t1 = b1/sb1   % t1 = -24.2979
pb1 = 2 * (1-tcdf(abs(t1),n-p) )  %pb1 = 1.4219e-011
t0 = b0/sb0  % t0 = 53.7800
pb1 = 2 * (1-tcdf(abs(t0),n-p) )  %pb1 =  1.1102e-015

% predicting y for the new observation x, CI and PI
newx = 40.5; %AGE = 40 1/2
ypred = b0 + b1 * newx   %ypred = 17.4988
sym = s * sqrt(1/n + (mean(x) - newx)^2/SXX ) %sym = 0.5167 s for mean y 
syp = s * sqrt(1 + 1/n + (mean(x) - newx)^2/SXX ) 
                                          %syp = 1.9462s for predicted y
%intervals CI and PI
alpha = 0.05;
%mean response interval
lbym = ypred - tinv(1-alpha/2, n-p) * sym;
rbym = ypred + tinv(1-alpha/2, n-p) * sym; 
% prediction interval
lbyp = ypred - tinv(1-alpha/2, n-p) * syp;
rbyp = ypred + tinv(1-alpha/2, n-p) * syp;
%print the intervals
[lbym rbym]   %[16.3731   18.6245]
[lbyp rbyp]   %[13.2585   21.7391]


figure(1)
% scatterplot of (x,y)'s and regression
plot(age, successrate, 'ko','markersize', msize)
hold on
plot(x, y, 'k-', 'linewidth', lw)
plot(x, yhat, 'r:', 'linewidth', lw)
hold off

figure(2)
%residuals against x
stem(x, res) 

figure(3)
%residuals against y
stem(y, res) 

figure(4)
%residuals against yhat
stem(yhat, res) 

% built in matlab's black-box m-file
XX = [vecones x x.^2 x.^3];
[bx,bintx,rx,rintx,statsx] = regress(y, XX)
 rcoplot(rx, rintx)
