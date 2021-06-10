 
clear all
close all

 disp('Diabetes Mellitus: C-Peptide')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;
 
 
Deficit =[-8.1 -16.1 -0.9 -7.8 -29.0 -19.2 -18.9 -10.6 -2.8...
 -25.0 -3.1 -7.8 -13.9 -4.5 -11.6 -2.1 -2.0 -9.0 -11.2 -0.2...
 -6.1 -1 -3.6 -8.2 -0.5 -2.0 -1.6 -11.9 -0.7 -1.2 -14.3 -0.8... 
 -16.8 -5.1 -9.5 -17.0 -3.3 -0.7 -3.3 -13.6 -1.9 -10.0 -13.5];

C_peptide =[ 4.8 4.1 5.2 5.5 5 3.4 3.4 4.9 5.6 3.7 3.9 ... 
 4.5 4.8 4.9 3.0 4.6 4.8 5.5 4.5 5.3 4.7 6.6 5.1 3.9 ...
 5.7 5.1 5.2 3.7 4.9 4.8 4.4 5.2 5.1 4.6 3.9 5.1 5.1 ...
 6.0 4.9 4.1 4.6 4.9 5.1];

%  Age=[ 5.2      8.8    10.5   10.6   10.4     1.8  12.7    15.6     5.8 ...   
%             1.9      2.2      4.8    7.9     5.2      0.9   11.8     7.9      11.5  ...  
%           10.6      8.5     11.1    12.8     11.3      1.0   14.5    11.9     8.1  ... 
%           13.8    15.5      9.8   11.0    12.4     11.1  5.1      4.8      4.2 ...
%             6.9      13.2      9.9  12.5    13.2     8.9   10.8  ]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2; %number of parameters (beta0, beta1)
%"Deficit" measurement is "x",  "C_peptide" is "y".
x = Deficit' ;  %column vector
y = C_peptide' ; %column vector
n = length(x);
% Sums of Squares
SXX = sum( (x - mean(x)).^2 )  %SXX=2.1310e+003
SYY = sum(  (y - mean(y)).^2 )  %SYY=21.807
SXY = sum( (x - mean(x)).* (y - mean(y)) ) %SXY=105.3477
% estimators of coefficients beta1 and beta0
b1 = SXY/SXX   %0.0494
b0 = mean(y) - b1 * mean(x)  %5.1494
% predictions
yhat = b0 + b1 * x; 
%residuals
res = y  - yhat;
% ANOVA Identity
SST = sum(  (y - mean(y)).^2 ) %which is SYY
SSR = sum(  (yhat - mean(y)).^2 ) %5.2079
SSE = sum(  (y - yhat).^2 ) %=sum(res.^2), 16.599
% forming F and test of adequacy of linear regression
MSR = SSR/(p - 1)  %5.2079
MSE = SSE/(n - p)  %estimator of variance, 0.4049
s = sqrt(MSE) %0.6363
F = MSR/MSE %12.8637
pvalue = 1-fcdf(F, p-1, n-p) 
%testing H_0: regression has beta1=0, 
%that is no need for linear fit, p-val = 0.00088412
%  Other measures of goodness of fit
R2 = SSR/SST  %0.2388
  % also R2 = 1- SSE/SST
  % this is corr(x,y) squared: 0.4887^2.
R2adj = 1 - (n-1)/(n-p)* SSE/SST %0.2203
%  Standard deviations of coefficient estimators
sb1 = s/sqrt(SXX) %0.0138
sb0 = s * sqrt(1/n + (mean(x))^2/SXX ) %0.1484

% are the coefficients equal to 0?
t1 = b1/sb1  %3.5866
pb1 = 2 * (1-tcdf(abs(t1),n-p) ) %0.00088412
t0 = b0/sb0  %34.6927
pb1 = 2 * (1-tcdf(abs(t0),n-p) ) %0
% 
% predicting y for the new observation x, CI and PI
newx = -5; %Deficit = -5
y_newx = b0 + b1 * newx  % 4.9022
sym = s * sqrt(1/n + (mean(x) - newx)^2/SXX ) 
  %st.dev. for mean response, sym = 0.1063
syp = s * sqrt(1 + 1/n + (mean(x) - newx)^2/SXX ) 
  %st.dev. for the prediction syp = 0.6451
alpha = 0.05;
%mean response interval
lbym = y_newx - tinv(1-alpha/2, n-p) * sym;
rbym = y_newx + tinv(1-alpha/2, n-p) * sym;
% prediction interval
lbyp = y_newx - tinv(1-alpha/2, n-p) * syp;
rbyp = y_newx + tinv(1-alpha/2, n-p) * syp;
%print the intervals
[lbym rbym]   % 4.6875    5.1168
[lbyp rbyp]     % 3.5994    6.2050

figure(1)
% scatterplot of (x,y)'s and regression
% scatterplot of (x,y)'s and regression
plot(x, y,  'ko', 'MarkerSize',12,  'MarkerFaceColor','g')
hold on
plot(x, yhat, 'r-', 'linewidth', lw)
xlabel('Alcaline Defficiency')
ylabel('C-Peptide Level (log units)')
hold off
figure(2)
%residuals against x
stem(x, res) 
xlabel('Alcaline Defficiency')
ylabel('Residuals')
figure(3)
%residuals against yhat
stem(yhat, res) 
xlabel('yhat')
ylabel('Residuals')



figure(4)
histn(res, -3, 0.4, 3)
hold on
plot((-3:0.1:3),normpdf( (-3:0.1:3), 0, 0.6363),'r-','linewidth',3)
xlabel('Residuals')
% built in matlab's black-box m-file
%%
XX = [ones(size(y)) x];
[bx,bintx,rx,rintx,statsx] = regress(y, XX)
