clear all
close all


 disp('Multiple Linear Regression Diagnostics: Body Fat Example')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;
load('C:\BESTAT\Reg\Regdat\fat.dat')
%load('C:\Courses\bmestatu\fat.dat')
%load  'fat.dat'
casen = fat(:,1); %case number
broz = fat(:,2); %dependent variable Y
siri = fat(:,3); %function of densi
densi = fat(:,4); %invasive (makes people wet!)
age = fat(:,5); %below are predictors
weight = fat(:,6);   
height = fat(:,7);   
adiposi = fat(:,8); 
    % ffwei is functionally dependent 
    % on density//brozek and siri indices.
neck = fat(:,10); 
chest = fat(:,11);
abdomen  = fat(:,12);   
hip = fat(:,13); 
thigh = fat(:,14);
knee = fat(:,15);   
ankle = fat(:,16);
biceps = fat(:,17);
forearm = fat(:,18);     
wrist = fat(:,19);
vecones = ones(size(broz)); % necessary for intercept
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 disp('========================================================')
 disp('        p = 15, 14 variables + intercept')
 disp('========================================================')
 
Z =[age  weight  height  adiposi   neck    chest    abdomen ...
    hip  thigh   knee    ankle     biceps  forearm  wrist];
X =[vecones Z];
Y = broz 
% X is design matrix, n x p where n is the number of subjects and p is
% the number of parameters, or number of predictors + 1.
% varnames = ['intercept=0' 'age=1' 'weight=2' 'height=3' 'adiposi=4' 
% 'neck=5' 'chest=6' 'abdomen=7' 'hip=8' 'thigh=9' 'knee=10' 'ankle=11' 
% 'biceps=12' 'forearm=13' 'wrist=14'];
[n, p] = size(X)
b = inv(X' * X) * X'* Y;
H = X * inv(X' * X) * X';
Yhat = H * Y;  %or Yhat  = X * b;

%----------------------------------------------------
%Sums of Squares
J=ones(n); I = eye(n);
SSR = Y' * (H - 1/n * J) * Y;
SSE = Y' * (I - H) * Y;
SST = Y' * (I - 1/n * J) * Y;
MSR = SSR/(p-1) 
MSE = SSE/(n-p) 
F = MSR/MSE 
 pval = 1-fcdf(F, p-1, n-p)
 Rsq = 1 - SSE/SST
 Rsqadj = 1 - (n-1)/(n-p) * SSE/SST
 s = sqrt(MSE)
 
%------------------------------------------------------

sig2 = MSE * inv(X' * X);% covariances among b's}
sb=sqrt(diag(sig2));
tstats = b./sb;
pvals = 2 * tcdf(-abs(tstats), n-p); 
 disp('--------------------------------------')
 disp('     var#       t        pval         ')
 disp('--------------------------------------')
 [ (0:p-1)'   tstats     pvals   ]
  %--------------------------------------------------- 
  % predicting mean responses and individual response
  % with 95% confidence intervals
  % Measurements for an individual starting with 1 for 
  % the intercept
Xh=[1  38  191  72  26  41   104  ...
   95  101.5   66  39   24   31   30  18.5];
%predicted brozek index
  Yh = Xh  * b   %19.5143
 %variances for mean response and for predicted response
  sig2h = MSE *  Xh  * inv(X' * X)  * Xh'; 
  sig2hpre = MSE * (1 + Xh  * inv(X' * X)  * Xh');
  sigh = sqrt(sig2h);
  sighpre = sqrt(sig2hpre);
  [Yh - tinv(0.975, n-p) * sigh,  Yh + tinv(0.975, n-p) * sigh]  
      %17.4347   21.5940
  [Yh - tinv(0.975, n-p) * sighpre,  Yh + tinv(0.975, n-p) * sighpre]
      %11.3721   27.6566
      
%

%============ testing H_0: beta1 = 0 =======================
Xno1 = [vecones weight  height adiposi  neck chest abdomen ...
    hip thigh knee ankle biceps forearm wrist]; %age excluded
Hno1 = Xno1 * inv(Xno1' * Xno1) * Xno1';
SSEno1 = Y' * (I - Hno1) * Y;
SSRno1  = SSEno1 - SSE;  %important -- increase in SSE when var 1 is omitted
Fno1 = (SSRno1/1)/MSE    %3.5881
pval = 1-fcdf(Fno1, 1, n-p)   %0.0594; beta1 not significantly  different
                                        % from zero; age may not be a good predictor
    
% check that (b2/sb2)^2 = F2_13. b2/sb2 is tstats(3)
%====== testing H_0: beta2 = beta3 = beta4 = 0 ==============

Xno234 = [vecones age neck chest abdomen ...
  hip thigh knee ankle biceps forearm wrist]; %weight, height, adiposi excluded
Hno234 = Xno234 * inv(Xno234' * Xno234) * Xno234';
SSEno234 = Y' * (I - Hno234) * Y;
SSRno234  = SSEno234 - SSE;
Fno234 = (SSRno234/3)/MSE     %1.7800
pval = 1-fcdf(Fno234, 3, n-p)      %0.1516

%======================== Influence diagnostics===================
 %prediction of y_i with ith observation removed 
 %hat y_i(-i). Check how close hat y_i(-i) and hat y_i are.
 ind = 1:n;
 Yhati = [];  
 for i = 1:n
     indi = find(ind ~= i);
     Yi = Y(indi);
     Xi=X(indi,:);
     bi = inv(Xi' * Xi) * Xi'* Yi;  
     Yhatii = X(i,:) * bi;
     Yhati =[Yhati; Yhatii];   
 end
 Yhati  %prediction of y_i without i-th observation 
 
 figure(10)
 plot(Yhat, 'bo')
 hold on
 plot(Yhati, 'ro')
 for i = 1:n
 plot([i i],[Yhat(i) Yhati(i)],'k-','linewidth',3)
 end
 xlabel('Index','Interpreter','LaTeX')
 ylabel('$\hat y$ (blue) and $\hat{y}(-i)$ (red)','Interpreter','LaTeX')
 axis tight
 %42nd and 39th have the maximum absolute difference: 
 %20.9876 and 6.5254. Median abs difference is 0.1265
%  print -depsc 'C:\Springer\Reg\Regeps\yhatminusi.eps'
%%
 %============ residual analysis================== 
 %---------- studentized residuals----------------
 hii = diag(H);                  %``leverages''
 resid = (I - H)*Y;            %ordinary residuals
 sresid = sqrt(MSE .* (1-hii));
 stresid = resid./sresid;  %studentized residuals
 %----------studentized deleted residuals---------
 di = Y - Yhati; %or  equivalently di  = resid./(1-hii) 
 % d's are also called PRESS residuals
 PRESS =sum(di.^2)  %PRESS=PRedictive Error Sum of Squares"
 R2pred = 1-PRESS/SST %R^2 predictive
 %---------externally studentized residuals -------------------
 sminusi = sqrt(((n-p)*MSE*ones(n,1) -...
         resid.^2./(1-hii))/(n-p-1));     %stdev(-i)
 ti  = resid ./(sminusi .* sqrt(1-hii));
 %------ outliers judged by leverages  hii
 outli=hii/mean(hii);
 find(outli > 3)
   % 31 36 39 41 42  86 159 175  206
 find(outli > 6)
   % 39 42
   
%==============influential observations===========
 Dffits = ti .* sqrt( hii ./(1-hii)); %influ ith to ith
 find(abs(Dffits) > 2 * sqrt(p/n))  
 %    31  39  42  82  86 128 140
 %   175 207 216 221 231 250

 
 figure(1)
 scatter(Yhat ,resid,100*abs(Dffits),'green','filled')
 hold on
 scatter(Yhat,resid,100*abs(Dffits),'black')
 xlabel('$\hat y$','Interpreter','Latex') 
 ylabel('Residuals scaled by $|\mbox{\it Dffits}|$','Interpreter','Latex')
 axis([0 46 -12 10])
%     39  2.5162
%     42  5.3767
%print -depsc 'C:\BESTAT\Reg\Regeps\fatdffits.eps'


 CooksD = resid.^2 .* (hii./(1-hii).^2)/(p * MSE);
 % influence if ith to all; 
 find(CooksD > 4/n) %find influential
 %31 39 42 82  86 128 175 207 216 221 250

 %
 figure(2)
 scatter(Yhat,resid, 1000*CooksD,'green','filled')
 hold on
 scatter(Yhat,resid, 1000*CooksD,'black')
  xlabel('$\hat y$','Interpreter','Latex') 
 ylabel('Residuals scaled by {\it CooksD}','Interpreter','Latex')  
 axis([0 46 -12 10])
%     39 0.4063
%     42 1.9101
 %print -depsc 'C:\BESTAT\Reg\Regeps\fatcooksd.eps'

 %%
 %DFBetas - influence if ith obs on jth coefficient
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
 %
  figure(3)
  imagesc( 1- (abs(DFBetas) > 2/sqrt(n))   ) 
 % check if any obs is influential? 
  colormap gray
  %print -depsc 'C:\BESTAT\Reg\Regeps\fatdfbetas.eps'
 %%
 
 %================== multicollinearity========================
 lambdas = eig(X' * X);
K = max(lambdas)/min(lambdas)  %3.2114e+008 >>> 10-100
Ki = max(lambdas)./lambdas
  Xprime = [];
   for k=2:p
       Xkbar = mean(X(:,k));
       sk = std( X(:,k));
       Xprimek = 1/sqrt(n-1) .* (X(:,k)- Xkbar * ones(n,1) )/sk;
       Xprime = [Xprime Xprimek];
   end
   RXX = Xprime' * Xprime;
   VIF = diag (inv(RXX)) 
 %   2.2509   33.7869    2.2566   16.1634    4.4307   10.6846   13.3467
 %  15.1583    7.9615    4.8288    1.9455    3.6745    2.1934    3.3796

   
   %%
   
[index,res,stud_res,lev,DFFITS1,Cooks_D,DFBETAS]=diagnostics(Y,Z); 
s = regstats(Y,Z,'linear',{'all'});
   
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
 
%    plot(CooksD+1,'k-')
%    hold on
%    plot(Cooks_D-1,'r-')
%    plot(s.cookd,'b-')
   
%    plot(Dffits,'k-')
%    hold on
%    plot(DFFITS1,'r-')
%    plot(s.dffits,'b-')
   
%    plot(stresid,'k-')
%    hold on
%    plot(stud_res,'r-')
%    plot(s.studres,'b-')

%    plot(hii,'k-')
%    hold on
%    plot(lev,'r-')
%    plot(s.leverage,'b-')

%    plot(DFBetas,'k-')
%    hold on
%    plot(DFBETAS,'ko')
%    plot(DFBETAS,'r-')
%    hold on
%    plot(s.dfbetas','g-')

   
%     plot(di,'k-')
%     hold on
%     plot(dii,'r-') 

 
 
 
 
 
 