
clear all
close all


 disp('Multivariate Linear Regression: Silver-Zinc Batteries')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;
load   'C:\STAT1\Reg\Regdat\silverzinc.mat' 
%%
 x1=silverzinc.chr;
 x2=silverzinc.dchr.^2;
 x3=silverzinc.ddch;
 x4=silverzinc.temp.^0.5;
 x5=silverzinc.ecv; 
 y=log(silverzinc.ctf)
vecones = ones(size(y));
 X = [x1 x2 x3 x4 x5]
 %%
Xdes=[vecones X]
[betahat,Ibeta,res,Ires,stats] = regress(y,Xdes)
%%
stepwise(X, y)

qqplot(res)