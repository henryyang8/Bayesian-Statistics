
%Data from http://lib.stat.cmu.edu/DASL/Stories/CheddarCheeseTaste.html
%As cheddar cheese matures, a variety of chemical processes take place. 
%The taste of matured cheese is related to the concentration of several 
%chemicals in the final product. In a study of cheddar cheese from the 
%LaTrobe Valley of Victoria, Australia, 
%samples of cheese were analyzed for their chemical composition and were 
%subjected to taste tests. Overall taste scores were obtained by combining 
%the scores from several tasters.



clear all
close all

Taste  = [12.3 39.0  5.6 37.3 18.1 34.9  0.7 54.9 15.9 18.0 14.0 ...
          32.0 16.8 26.5 13.4 20.9 47.9 25.9 21.9 21.0 57.2 25.9 ...
          40.9  6.4 38.9 15.2 56.7 11.6  0.7  5.5];
Acetic = [ 4.54 5.37 4.66 5.89 4.90 5.74 4.48 6.15 4.79 5.25 ...
           4.56 5.46 5.37 6.46 5.80 5.16 5.76 5.70 6.08 5.24 ...
           6.45 5.24 6.37 5.41 5.44 5.30 5.86 6.04 5.33 6.18];
H2S   =  [ 3.13 5.44 3.81 8.73 3.85 6.14 3.00 6.75 3.91 6.17 ...
           4.95 9.24 3.66 6.92 6.69 5.04 7.59 7.60 7.97 4.17 ...
           7.91 4.94 9.59 4.70 9.06 5.22 10.2 3.22 3.91 4.79];
Lactic = [ 0.86 1.57 0.99 1.29 1.29 1.68 1.06 1.52 1.16 1.63 ...
           1.15 1.44 1.31 1.72 1.08 1.53 1.81 1.09 1.78 1.58 ...
           1.90 1.30 1.74 1.49 1.99 1.33 2.01 1.46 1.25 1.25];
y =Taste'; 
Z =[Acetic' H2S' Lactic'];
[index,res,stud_res,lev,DFFITS,Cooks_D,DFBETAS]=diagnostics(y,Z);

Y =Taste';
p = 4;
n = length(Taste);
X=[ones(n,1) Acetic' H2S' Lactic'];
b = (X' * X)^-1 * X'* Y;
H = X * (X' * X)^-1 * X';
max(max(H * H - H));
J=ones(n); I = eye(n);
SSR = Y' * (H - 1/n * J) * Y;
SSE = Y' * (I - H) * Y;
SST = Y' * (I - 1/n * J) * Y;
MSR = SSR/(p-1);
MSE = SSE/(n-p);
F = MSR/MSE;
sig2 = MSE * (X' * X)^-1;
sb=sqrt(diag(sig2));
tstats = b./sb;
Xh = [1  5 8 2]';
Yh = Xh' * b;
sig2h = MSE *  Xh' * (X' * X)^-1 * Xh; 
sig2hpre = MSE * (1 + Xh' * (X' * X)^-1 * Xh);
sigh = sqrt(sig2h);
sighpre = sqrt(sig2hpre);
%============ testing H_0: beta2 = 0 =======================
Xal = [ones(n, 1) Acetic' Lactic'];
Hal = Xal * (Xal' * Xal)^-1 * Xal';
SSEal = Y' * (I - Hal) * Y;
SSR2_13 = SSEal - SSE;
F2_13 = (SSR2_13/1)/MSE   
pval = 1-fcdf(F2_13, 1, n-p)   


% check that (b2/sb2)^2 = F2_13. b2/sb2 is tstats(3)
%============ testing H_0: beta2 = beta3 = 0 =======================
Xa = [ones(n, 1) Acetic'];
Ha = Xa * (Xa' * Xa)^-1 * Xa';
SSEa = Y' * (I - Ha) * Y;
SSR23_1 = SSEa - SSE;
F23_1 = (SSR23_1/2)/MSE 
pval = 1-fcdf(F23_1, 2, n-p)  


%===================================================================
ind = 1:n;
 Yhati = []; esternalstd =[];
 for i = 1:n
     indi = find(ind ~= i);
     Yi = Y(indi);
     Xi=X(indi,:);
     bi = inv(Xi' * Xi) * Xi'* Yi;  
          Hi = Xi * inv(Xi' * Xi)  * Xi';
         SSEi  = Yi' * (eye(n-1) - Hi) * Yi;
         MSEi  = SSEi./(n-p-1); 
         externalstd=[esternalstd MSEi]  
     Yhatii = X(i,:) * bi;
     Yhati =[Yhati; Yhatii];   
 end
 Yhati  %prediction of y_i with ith observationb removed
 %================== outliers============================
 %---------- studentized residuals-------------
 hii = diag(H); %leverage
 resid = (I - H)*Y;          %ordinary residuals
 sresid = sqrt(MSE .* (1-hii));
 stresid = resid./sresid  %studentized residuals
 %-----------studentized deleted residuals-----
 di = Y - Yhati;
 dii = resid./(1-hii) 
 sminusi = sqrt(((n-p-1)*MSE*ones(n,1) -...
     resid.^2./(1-hii))/(n-p)); 
 ti  = resid ./(sminusi .* sqrt(1-hii))
 % externally studentized residuals
 %-----------leverage--------------------------
 outli=hii/mean(hii);
 find(outli > 2)
 %-----------leverage----------------------------
 % outliers based on leverage = hii
 outli=hii/mean(hii);
 find(outli > 2)

 %===================== influential observations=====
 Dffits = ti .* sqrt( hii ./(1-hii))  % infl ith to ith
 find(Dffits > 2 * sqrt(p/n));
 CooksD = resid.^2 .* (hii./(1-hii).^2)/(p * MSE) 
 % influence if ith to all; 
 find(CooksD > 4/n) %find influential
 cii = diag(inv(X' * X));
 DFBetas =[];
 for i = 1:n
     indi = find(ind ~= i);
     Yi = Y(indi);
     Xi=X(indi,:);
     bi  = inv(Xi' * Xi) * Xi'* Yi;
     Hi = Xi * inv(Xi' * Xi)  * Xi';
     SSEi  = Yi' * (eye(n-1) - Hi) * Yi;
     MSEi  = SSEi./(n-p-1); 
 DFBetasi = (b - bi)./sqrt(MSEi .* cii) ;  
 DFBetas = [DFBetas; DFBetasi'];  
 end
 
  infl = find(abs(DFBetas) > 2/sqrt(n) ) 
 % check if any obs is influential? 
 
 PRESS =sum(di.^2)
 R2pred = 1-PRESS/SST
 
 %================== multicollinearity========================
 Xprime = [];
   for k=2:p
       Xkbar = mean(X(:,k));
       sk = std( X(:,k));
       Xprimek = 1/sqrt(n-1) .* (X(:,k)- Xkbar * ones(n,1) )/sk;
       Xprime = [Xprime Xprimek];
   end
   RXX = Xprime' * Xprime;
   VIF = diag (inv(RXX)) 
   
   
   
y =Taste'; 
Z =[Acetic' H2S' Lactic'];
[index,res,stud_res,lev,DFFITS1,Cooks_D,DFBETAS]=diagnostics(y,Z);

%MATLAB's   regstats

s = regstats(y,Z,'linear',{'all'});
 
% KEY for option {'all'}
%        'Q'           Q from the QR Decomposition of the design matrix
%        'R'           R from the QR Decomposition of the design matrix
%        'beta'        Regression coefficients
%        'covb'        Covariance of regression coefficients
%        'yhat'        Fitted values of the response data
%        'r'           Residuals
%        'mse'         Mean squared error
%        'rsquare'     R-square statistic
%        'adjrsquare'  Adjusted R-square statistic
%        'leverage'    Leverage
%        'hatmat'      Hat (projection) matrix
%        's2_i'        Delete-1 variance
%        'beta_i'      Delete-1 coefficients
%        'standres'    Standardized residuals
%        'studres'     Studentized residuals
%        'dfbetas'     Scaled change in regression coefficients
%        'dffit'       Change in fitted values
%        'dffits'      Scaled change in fitted values
%        'covratio'    Change in covariance
%        'cookd'       Cook's distance
%        'tstat'       t statistics and p-values for coefficients
%        'fstat'       F statistic and p-value
%        'dwstat'      Durbin-Watson statistic and p-value
%        'all'         Create all of the above statistics
