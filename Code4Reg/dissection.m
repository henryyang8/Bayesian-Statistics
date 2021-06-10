close all
clear all


 disp('Regression: Sharp Dissection')
 lw = 2.5; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;

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
% The authors found linear relationship between the logarithm of the
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

% lasd = log(shdiss);
lasd= [ 2.6391    4.6821    3.6636    5.7398    3.1781    4.7185 ...
    4.6444    5.9454    3.7377    4.3041    4.2047    4.9767...
    3.0445    4.5326    4.3175    5.9428    3.5835    3.5835 ...
    4.2905    5.4765    3.5553    4.2341];


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
plot(lasd,  sesco, 'ko', 'MarkerSize',14,  'MarkerFaceColor','g')
hold on

 [B,BINT,R,RINT,STATS] =regress( sesco', [ones(size(lasd')) lasd'] )

 figure(2)
 plot(sesco, R, 'ko', 'MarkerSize',10,  'MarkerFaceColor','g')
hold on
plot([0 20],[0,0],'LineWidth',2)
  disp('Linear Regression: ')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;
figure(3)
qqplot(R)

x = lasd';
y = sesco';
n = length(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2; %number of parameters (beta0, beta1)
% beta0 estimated by b0; beta1 estimated by b1
%sums of squares
% 'Abdomen measurement is "x" and the response brox is "y".
% Sums of Squares
SXX = sum( (x - mean(x)).^2 )   %17.7180
SYY = sum(  (y - mean(y)).^2 )  %283.4545
SXY = sum( (x - mean(x)).* (y - mean(y)) )  %67.6813

% estimators of coefficients beta1 and beta0
b1 = SXY/SXX   %3.8199
b0 = mean(y) - b1 * mean(x)  %-5.0396
figure(1)
plot( 2:0.1:7, b0 + b1 .* (2:0.1:7), 'r-','LineWidth', 3)
legend('(Log Dissection, Adhesion Score)','Linear Fit',2)
xlabel('Log Dissection')
ylabel('Adhesion Severity Score')
% predictions
yhat = b0 + b1 * x;
%residuals
res = y - yhat;
% ANOVA Identity
SST = sum(  (y - mean(y)).^2 ) %which is SYY = 283.4545
SSR = sum(  (yhat - mean(y)).^2 )  % 258.5367
SSE = sum(  (y - yhat).^2 ) %which is sum(res.^2) = 24.9179
% forming F and test of adequacy of linear regression
MSR = SSR/(p - 1)  % 258.5367
MSE = SSE/(n - p)  %should be sigma2hat,     1.2459
F = MSR/MSE   %207.5109
pvalue = 1-fcdf(F, p-1, n-p) %H_0: regression has beta1=0, no need for 
% linear fit    pval=  5.0626e-012
%  Other measures of goodness of fit
R2 = SSR/SST  %0.9121
R2adj = 1 - (n-1)/(n-p)* SSE/SST  %0.9077
s = sqrt(MSE)   % 1.1162
%  Standard error of coefficient estimators
sb1 = s/sqrt(SXX)  % 0.2652
sb0 = s * sqrt(1/n + (mean(x))^2/SXX )   %1.1695

% are the coefficients equal to 0?
t1 = b1/sb1    %14.4052
pb1 = 2 * (1-tcdf(abs(t1),n-p) )  %    5.0626e-012
t0 = b0/sb0  %  -4.3093
pb0 = 2 * (1-tcdf(abs(t0),n-p) )    %  3.4137e-004

% predicting y for the new observation x, CI and PI
newx = 4; % corresponding to amount exp(4)=54.5982
ypred = b0 + b1 * newx   %10.24
sym = s * sqrt(1/n + (mean(x) - newx)^2/SXX ) %s for y mean  0.2525
syp = s * sqrt(1 + 1/n + (mean(x) - newx)^2/SXX ) %s for y prediction   1.1444

%intervals CI and PI
alpha = 0.05;
%mean response interval
lbym = ypred - tinv(1-alpha/2, n-p) * sym;
rbym = ypred + tinv(1-alpha/2, n-p) * sym;
% prediction interval
lbyp = ypred - tinv(1-alpha/2, n-p) * syp;
rbyp = ypred + tinv(1-alpha/2, n-p) * syp;
%print the intervals
[lbym rbym]  % 9.7134   10.7666
[lbyp rbyp]  % 7.8528   12.6272

%--------------------------------------
%YELLOW FINAL SPRING 2010
[b0, b1]
R2
t=(b1 - 3)/sb1
pval = 1- tcdf(t, n-2)
crit = tinv(0.95, n-2)

[b0 - tinv(0.975, n-2)*sb0, b0 + tinv(0.975, n-2)*sb0]

ym = b0 + b1 * 4
[ym - tinv(0.995, n-2) * sym,  ym + tinv(0.995, n-2) * sym]


% ans = -5.0396    3.8199
% R2 = 0.9121
% t =3.0920
% pval = 0.0029
% crit =1.7247
% ans = -7.4791   -2.6001
% ym =10.2400
% ans = 9.5216   10.9584



