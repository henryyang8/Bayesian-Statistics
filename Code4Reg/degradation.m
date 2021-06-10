 disp('Scaffold Degradation')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 16;
 msize = 10;

degradation =[...
0  42.520	 ;...	
0  71.590    ;...	
0  40.063	 ;...	
1  68.397	 ;...	
1  53.527	 ;...	
1  40.676	 ;...
2  21.724	 ;...	
2  35.032	 ;...	
2  56.687	 ;...	
3  44.029	 ;...	
3  45.058	 ;...	
3  27.579	 ;...	
4  44.929	 ;...	
4  29.348	 ;...	
4  37.259	 ;...	
5  28.625	 ;...	
5  25.956	 ;...	
5  20.179	 ;...	
6  14.994	 ;...	
6  9.051	 ;...	
6  14.923	 ;...
7  2.692	 ;...	
7  15.688	 ;...	
7   3.420	 ];
%---------------------------
day   = degradation(:,1);
mod1  = degradation(:,2);   


close all; plot(day,  mod1, 'o','MarkerSize',msize,...
'MarkerFaceColor','g','MarkerEdgeColor','k')
[b1 bint1 res1 resint1 stats1] = regress(  mod1, [ones(size(day)) day]);
hold on
plot( [-1:8], b1(1)+ b1(2)*[-1:8], 'r-', 'LineWidth',2)
xlabel('day'); ylabel('Modulus at f=1')
axis tight
%print -depsc 'C:\Springer\Reg\Regeps\finalp.eps'


n = length(mod1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2; %number of parameters (beta0, beta1)
% beta0 estimated by b0; beta1 estimated by b1
%sums of squares
% 'Abdomen measurement is "x" and the response brox is "y".
% Sums of Squares
SXX = sum( (day - mean(day)).^2 )  %126
SYY = sum(  (mod1 - mean(mod1)).^2 ) %8.0763e+003
SXY = sum( (day - mean(day)).* (mod1 - mean(mod1)) ) % -832.9750

% estimators of coefficients beta1 and beta0
b1 = SXY/SXX   % -6.6109
b0 = mean(mod1) - b1 * mean(day) % 56.2193
% predictions
mod1hat = b0 + b1 * day;
%residuals
res = mod1 - mod1hat;
% ANOVA Identity
SST = sum(  (mod1 - mean(mod1)).^2 ) %8.0763e+003 which is SYY
SSR = sum(  (mod1hat - mean(mod1)).^2 ) %5.5067e+003
SSE = sum(  (mod1 - mod1hat).^2 ) %2.5695e+003 which is sum(res.^2)
% forming F and test of adequacy of linear regression
MSR = SSR/(p - 1) %5.5067e+003
MSE = SSE/(n - p)  %116.7973 should be sigma2hat
F = MSR/MSE   %47.1477
pvalue = 1-fcdf(F, p-1, n-p) %6.7600e-007 H_0: regression has beta1=0, no need for 
% linear fit
%  Other measures of goodness of fit
R2 = SSR/SST %0.6818
R2adj = 1 - (n-1)/(n-p)* SSE/SST  %0.6674
s = sqrt(MSE)  %10.8073
%  Standard error of coefficient estimators
sb1 = s/sqrt(SXX)  %0.9628
sb0 = s * sqrt(1/n + (mean(mod1))^2/SXX ) %31.9264

% are the coefficients equal to 0?
t1 = b1/sb1  %-6.8664
pb1 = 2 * (1-tcdf(abs(t1),n-p) )  %6.7600e-007
t0 = b0/sb0 %1.7609
pb0 = 2 * (1-tcdf(abs(t0),n-p) ) %0.0922

% predicting y for the new observation x, CI and PI
newx = 5.5; 
mod1pred = b0 + b1 * newx  %19.8593
sym = s * sqrt(1/n + (mean(day) - newx)^2/SXX ) % 2.9282,  s for y mean
syp = s * sqrt(1 + 1/n + (mean(day) - newx)^2/SXX ) %11.1969, s for y prediction
%intervals CI and PI
alpha = 0.05;
%mean response interval
lbym = mod1pred - tinv(1-alpha/2, n-p) * sym;
rbym = mod1pred + tinv(1-alpha/2, n-p) * sym;
% prediction interval
lbyp = mod1pred - tinv(1-alpha/2, n-p) * syp;
rbyp = mod1pred + tinv(1-alpha/2, n-p) * syp;
%print the intervals
[lbym rbym] %  13.7865   25.9320
[lbyp rbyp] %-3.3618   43.0803

% PINK
% (a) Use the moduli for frequency $f= 1.$  Write down linear equation model
% {\tt  mod1 = b0 + b1 * day}, where $b_0$ and $b_1$ are estimators of the population
% intercept and slope. What is $R^2$ for your regression. [56.2193 - 6.6109
% day], 0.6818.
% 
% (b) Test the hypothesis that population intercept is equal to 100,
% versus the alternative that it is smaller than 100.
% 
tpink = (b0 - 100)/sb0   %-1.3713
ppink = tcdf(tpink, n-p) %0.0921
% (c) Find 96 \% Confidence Interval for the population slope.
% 
[b1 - norminv(0.98) * sb1, b1 + norminv(0.98) * sb1]  % -8.5882   -4.6336
% (d) At time {\tt day = 5.5} find the  prediction of the modulus.
% What is the standard deviation for this predicted value?
%   19.8593, 11.1969


