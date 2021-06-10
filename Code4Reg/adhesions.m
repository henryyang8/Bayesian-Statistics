%%
%----------------------------------------------------------

% Postoperative adhesions are formed after surgical cardiac and great vessels procedures, 
% as part of the healing process. The scar tissue makes the re-entry complex,
% increasing the rate of iatrogenic lesions.
% 
% Currently, as reoperations represent 10 to 20\% of heart surgeries,
% various methods have been investigated to prevent or decrease the
% severity of postoperative adhesions.
% 
% The surgical time spent in the adhesiolysis procedure and 
% amounts of sharp dissection are new powerful measurement
% tools to evaluate the outcomes of methods used to prevent 
% pericardial adhesions, as reported by Lopes et al 2009.
%The authors found linear relationship between the logarithm of the
% amount of sharp dissection [lasd] and severity score [sesco].

% Jackson Brandão Lopes;
% Luís Alberto Oliveira Dallan;
% Luiz Felipe Pinho Moreira;
% Mário Castro Carreiro;
% Flávia Luana Barbosa Rodrigues; 
% Pedro de Castro Mendes; 
% Noedir Antônio Groppo Stol
% New quantitative variables to measure postoperative 
% pericardial adhesions. Useful tools in experimental research.
% Acta Cir. Bras. vol.24 no.2 São Paulo Mar./Apr. 2009

%data from Lopes
shdiss =[14 108	39  311	24	112	104 ...
382	42	74	67 145	21	93	75 ...
381	36	36	73  239	35	69];

lasd = log(shdiss);

% % data from Lopes, not used here
% dtime =[ 7	8.13	7.9 ...
% 38.62	7.68	8.47	9.33 ...
% 40.3	8.1	9.18	11.12 ...
% 20.83	6.32	11.97	8.98 ...
% 42.17	7.1	8.28	11.67 ...
% 27.67	6.75	9.65];
     
%data from Lopes: Adhesion Score
sesco  =[6	12	9 ...
18	7	12	12 ...
17	8	11	14 ...
15	7	12	13 ...
18	9	9	10 ...
16	7	10];


%plot(shdiss, lasd, 'o')       

%              
% shdissmy = [ 12   26   28   36   37  40 ...
%                    49  83  94  96  99  99 105 120 ...
%                    131 134  135 156 251 332  411 412]
%                
% lasdmy = [ 2.4849    3.2581    3.3322    3.5835 ...
%     3.6109    3.6889    3.8918   4.4188 ...
%     4.5433    4.5643    4.5951    4.5951 ...
%     4.6540    4.7875    4.8752    4.8978 ...
%     4.9053    5.0499    5.5255    5.8051 ...
%     6.0186    6.0210];                      
% sescomy =    [6   7   7  7  9  9  8  14 13  10 ...
%         10 10 11 12 12 12  12  15 16  18 17 18]

figure(1)
plot(lasd,  sesco, 'o', 'MarkerSize',14)

 [B,BINT,R,RINT,STATS] =regress( sesco', [ones(size(lasd')) lasd'] )

 figure(2)
 plot(sesco, R,'*')

  disp('Linear Regression: ')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;


x = lasd';
y = sesco';
n = length(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2; %number of parameters (beta0, beta1)
% beta0 estimated by b0; beta1 estimated by b1
%sums of squares
% 'Abdomen measurement is "x" and the response brox is "y".
% Sums of Squares
SXX = sum( (x - mean(x)).^2 )   %17.9017
SYY = sum(  (y - mean(y)).^2 )  %279.5000
SXY = sum( (x - mean(x)).* (y - mean(y)) )  %65.8276

% estimators of coefficients beta1 and beta0
b1 = SXY/SXX   %3.6772
b0 = mean(y) - b1 * mean(x)  %-5.0651
% predictions
yhat = b0 + b1 * x;
%residuals
res = y - yhat;
% ANOVA Identity
SST = sum(  (y - mean(y)).^2 ) %which is SYY = 279.5000
SSR = sum(  (yhat - mean(y)).^2 )  % 242.0592
SSE = sum(  (y - yhat).^2 ) %which is sum(res.^2) = 37.4408
% forming F and test of adequacy of linear regression
MSR = SSR/(p - 1)  % 242.0592
MSE = SSE/(n - p)  %should be sigma2hat,     1.8720 
F = MSR/MSE   %129.3023
pvalue = 1-fcdf(F, p-1, n-p) %H_0: regression has beta1=0, no need for 
% linear fit    pval= 3.4983e-010
%  Other measures of goodness of fit
R2 = SSR/SST  %0.8660
R2adj = 1 - (n-1)/(n-p)* SSE/SST  %0.8593
s = sqrt(MSE)   % 1.3682
%  Standard error of coefficient estimators
sb1 = s/sqrt(SXX)  % 0.3234
sb0 = s * sqrt(1/n + (mean(x))^2/SXX )   %1.4857   %%%3.7303

% are the coefficients equal to 0?
t1 = b1/sb1    %11.3711
pb1 = 2 * (1-tcdf(abs(t1),n-p) )  %  3.4983e-010
t0 = b0/sb0  % -3.4093
pb0 = 2 * (1-tcdf(abs(t0),n-p) )    % 0.0028 %%% 0.1896

% predicting y for the new observation x, CI and PI
newx = 4; %
ypred = b0 + b1 * newx   % 9.6436
sym = s * sqrt(1/n + (mean(x) - newx)^2/SXX ) %s for y mean   0.3343
syp = s * sqrt(1 + 1/n + (mean(x) - newx)^2/SXX ) %s for y prediction   1.4085

%intervals CI and PI
alpha = 0.05;
%mean response interval
lbym = ypred - tinv(1-alpha/2, n-p) * sym;
rbym = ypred + tinv(1-alpha/2, n-p) * sym;
% prediction interval
lbyp = ypred - tinv(1-alpha/2, n-p) * syp;
rbyp = ypred + tinv(1-alpha/2, n-p) * syp;
%print the intervals
[lbym rbym]  % 8.9463   10.3409
[lbyp rbyp]  % 6.7055   12.5816

%--------------------------------------
%YELLOW
[b0, b1]
R2
t=(b1 - 3)/sb1
pval = 1- tcdf(t, n-2)
crit = tinv(0.95, n-2)

[b0 - tinv(0.975, n-2)*sb0, b0 + tinv(0.975, n-2)*sb0]

ym = b0 + b1 * 4
[ym - tinv(0.995, n-2) * sym,  ym + tinv(0.995, n-2) * sym]

% ans =   -5.0651    3.6772
% R2 = 0.8660
% t =2.0940
% pval =0.0246
% crit =1.7247
%( [-12.8463    2.7161])                [-8.1642   -1.9660]
% ym = 9.6436
% 8.6924   10.5947

