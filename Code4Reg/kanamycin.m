% Kanamycin
% Miller (1980) describes a project and provides data on
% assessing the precision of non-invasive measuring of kanamycin
% concentration in neonates.
% 
% Premature babies are susceptible to infections and kanamycin (an aminoglycoside)
% is used for treatment of sepsis. Because kanamycin is ineffective at low doses and
% potentially harmful at high doses, it is necessary to constantly monitor its levels
% in a premature baby's body during the treatment.
% The standard procedure for measuring serum kanamycin levels is to take
% blood samples from a  heel. Unfortunately, due to frequent blood sampling,
% neonates are left with badly bruised  heels.
% 
% Kanamycin is routinely administered through an umbilical catheter.
% An alternative  procedure for measuring serum kanamycin would be to
% reverse flow in the catheter and draw blood sample from it.
% The concern about this non-invasive method is that the blood drawn 
% from the point close to infusion may have elevated level compared to blood samples from a 
% more distant points in the body.
% 
% In a carefully designed experimental setup, blood samples from 20 babies are
% obtained simultaneously from an umbilical catheter and a heel venapuncture (using heelstick).
% If the agreement is satisfactory,   physicians would be willing to use
% the catheter values instead of heelstick values.


clear all
close all
 disp('Kanamycin in Premature Babies')
 lw = 2.5; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;

%Baby Heelstick Catheter 
kmy = [...
 1 23.0 25.2;  2 33.2 26.0;  3 16.6 16.3;  4 26.3 27.2;...
 5 20.0 23.2;  6 20.0 18.1;  7 20.6 22.2;  8 18.9 17.2;...
 9 17.8 18.8; 10 20.0 16.4; 11 26.4 24.8; 12 21.8 26.8;...
13 14.9 15.4; 14 17.4 14.9; 15 20.0 18.1; 16 13.2 16.3;...
17 28.4 31.3; 18 25.9 31.2; 19 18.9 18.0; 20 13.8 15.6];

heel = kmy(:,2);   %response
cath = kmy(:,3);   %predictor
%
n = length(heel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 2; %number of parameters (beta0, beta1)
% beta0 estimated by b0; beta1 estimated by b1
% Sums of Squares
SXX = sum( (cath - mean(cath)).^2 ) %553.2900
SYY = sum( (heel - mean(heel)).^2 ) %494.5095
SXY = sum( (cath - mean(cath)).* (heel - mean(heel)) ) %435.4350

% estimators of coefficients beta1 and beta0
b1 = SXY/SXX %0.7870
b0 = mean(heel) - b1 * mean(cath) %4.2101

% predictions
heelhat = b0 + b1 * cath;
%residuals
res = heel - heelhat;
% ANOVA Identity
SST = sum(  (heel - mean(heel)).^2 ) %494.5095, which is SYY
SSR = sum(  (heelhat - mean(heel)).^2 ) %342.6840
SSE = sum(  (heel - heelhat).^2 ) %151.8255, which is sum(res.^2)
% forming F and test of adequacy of linear regression
MSR = SSR/(p - 1) %342.6840
MSE = SSE/(n - p)  %8.4347, should be sigma2hat
F = MSR/MSE  %40.6276
pvalue = 1-fcdf(F, p-1, n-p) %5.2876e-006
%  Other measures of goodness of fit
R2 = SSR/SST  %0.6930
R2adj = 1 - (n-1)/(n-p)* SSE/SST  %0.6759
s = sqrt(MSE)  %2.9043
%  Standard error of coefficient estimators
sb1 = s/sqrt(SXX)   %0.1235
sb0 = s * sqrt(1/n + (mean(cath))^2/SXX )   %2.6909

% are the coefficients equal to 0?
t1 = b1/sb1  %6.3740
pb1 = 2 * (1-tcdf(abs(t1),n-p) ) %5.2876e-006
t0 = b0/sb0  %1.5646
pb0 = 2 * (1-tcdf(abs(t0),n-p) )   %0.1351
%
% test the slope H0: beta1 = 1 vs H1: beta1 < 1
t = (b1 - 1)/sb1  %-1.7252
pval = tcdf(t, n-p) %0.0508
% test the slope H0: beta1 = 0.5 vs H1: beta1 > 0.5
t = (b1 - 0.5)/sb1  %2.3244
pval = 1-tcdf(t, n-p) %0.0160
% find 88% CI for the slope
[b1 - tinv(0.94, n-p)*sb1, b1 + tinv(0.94, n-p)*sb1]
   %0.5855    0.9885
%
% test the intercept H0: beta0 = 9 vs H1: beta0 < 9
t = (b0 - 9)/sb0    -1.7800  
pval = tcdf(t, n-p)   %0.0460
% test the intercept H0: beta0 = -5 vs H1: beta0 ~= -5
t = (b0 - (-5))/sb0     %3.4227
pval = 2 * tcdf(-abs(t), n-p)  %0.0030
% find 97% CI for the intercept
[b0 - tinv(0.985, n-p)*sb0, b0 + tinv(0.985, n-p)*sb0]
   %-2.1302   10.5504

% predicting y for the new observation x, CI and PI
newx = 20; %cath = 20
heelpred = b0 + b1 * newx  %19.9500
%s for mean resp
sym = s * sqrt(1/n + (mean(cath) - newx)^2/SXX ) %0.6648
%s for predicted resp
syp = s * sqrt(1 + 1/n + (mean(cath) - newx)^2/SXX ) %2.9794  
%intervals CI and PI
alpha = 0.05;
%mean response interval
lbym = heelpred - tinv(1-alpha/2, n-p) * sym;
rbym = heelpred + tinv(1-alpha/2, n-p) * sym;
% prediction interval
lbyp = heelpred - tinv(1-alpha/2, n-p) * syp;
rbyp = heelpred + tinv(1-alpha/2, n-p) * syp;
%print the intervals
[lbym rbym]
  % 18.5534   21.3466

[lbyp rbyp]
  % 13.6905   26.2094

% PLOTS
close all
figure(1)
plot([min(cath), max(cath)], [min(cath), max(cath)], 'r:','LineWidth',lw)
hold on
% scatterplot of (x,y)'s and regression
plot(cath, heel, 'ro',...
'MarkerSize',msize, 'MarkerEdgeColor','k', 'MarkerFaceColor','g')
plot(cath, heelhat, 'k-','LineWidth',lw)


figure(2)
%residuals against x
stem(cath, res) 

figure(3)
%residuals against order
stem((1:n), res) 

%
% built in matlab's black-box m-file


vecones = ones(size(heel));
X = [vecones cath];
[bx,bintx,rx,rintx,statsx] = regress(heel, X);

bx
bintx
rx
rintx
statsx

% bx =
%    4.2101
%    0.7870
% 
% 
% bintx =
%    -1.4433    9.8635
%     0.5276    1.0464
% 
% 
% rx =
% 
%    -1.0423
%     8.5281
%    -0.4381
%     0.6837
%    -2.4683
%     1.5453
%    -1.0813
%     1.1536
%    -1.2056
%     2.8832
%     2.6725
%    -3.5015
%    -1.4298
%     1.4637
%     1.5453
%    -3.8381
%    -0.4430
%    -2.8643
%     0.5240
%    -2.6872
% 
% 
% rintx =
% 
%    -7.0422    4.9575
%     4.4184   12.6378
%    -6.4150    5.5388
%    -5.2087    6.5760
%    -8.4322    3.4955
%    -4.4685    7.5592
%    -7.1696    5.0069
%    -4.8457    7.1529
%    -7.2618    4.8507
%    -2.9205    8.6869
%    -3.2136    8.5585
%    -9.1576    2.1546
%    -7.3088    4.4492
%    -4.3766    7.3040
%    -4.4685    7.5592
%    -9.4904    1.8142
%    -5.9255    5.0395
%    -8.1674    2.4388
%    -5.5316    6.5797
%    -8.4647    3.0903
% 
% 
% statsx =
% 
%     0.6930   40.6276    0.0000    8.4347

figure(4)
errorbar((1:n),rx(1:n), rintx(1:n,1), rintx(1:n,2))
hold on
plot((1:n),rx(1:n),'ro','markersize',5,'MarkerFaceColor','g') 
plot((1:n),1.96 * sqrt(statsx(4))*ones(1,n), 'r:')
plot((1:n),-1.96 * sqrt(statsx(4))*ones(1,n), 'r:')
axis tight


