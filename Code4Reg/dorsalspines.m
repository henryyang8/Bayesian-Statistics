
%BENKA LAKE  + GARDEN BAY LAKE + BIG LAKE
%               AS REGRESSION
clear all
close all force
  disp('Dorsal Spines')
  set(0, 'DefaultAxesFontSize', 16);
  lw =2.5;  fs = 16;  msize = 10;

stickleback =[...
4.2 4.4 4.9;  ...
4.1 4.6 4.6;  ...
4.2 4.5 4.3;  ...
4.3 4.2 4.9;  ...
4.5 4.4 4.7;  ...
4.4 4.2 4.4;  ...
4.5 4.5 4.5;  ...
4.3 4.7 4.4 ];
lakes = [1 1 1 1 1 1 1 1 ...
2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3];
[p table stats] = anova1(stickleback(:), lakes(:))

stats
% p =    0.0267
% 
% 
% table = 
% 
%     'Source'    'SS'        'df'    'MS'        'F'         'Prob>F'
%     'Groups'    [0.3033]    [ 2]    [0.1517]    [4.3260]    [0.0267]
%     'Error'     [0.7363]    [21]    [0.0351]          []          []
%     'Total'     [1.0396]    [23]          []          []          []
% 
% 
% stats = 
%     gnames: {3x1 cell}
%          n: [8 8 8]
%     source: 'anova1'
%      means: [4.3125 4.4375 4.5875]
%         df: 21
%          s: 0.1872
%[mu; alpha1; alpha2; alpha3]=
[1  1  0 0; 1 0 1 0; 1 0 0 1; 0 1 1 1]\[4.3125; 4.4375; 4.5875;  0] 
%     4.4458
%    -0.1333
%    -0.0083
%     0.1417
% 
% (a) Hypothesis $H_0$ stating that the mean lengths of dorsal spines are the same for the three lakes
% is rejected at the level $\alpha = 0.05$ since the $p$-value is 0.0267.
% 
% (b) If the significance level was $\alpha = 0.01$ the hypothesis $H_0$ would not be rejected.

multcompare(stats)
%     1.0000    2.0000   -0.3610   -0.1250    0.1110    0.3921
%     1.0000    3.0000   -0.5110   -0.2750   -0.0390    0.0206
%     2.0000    3.0000   -0.3860   -0.1500    0.0860    0.2668

%%
stickleback =[...
4.2 4.4 4.9;  ...
4.1 4.6 4.6;  ...
4.2 4.5 4.3;  ...
4.3 4.2 4.9;  ...
4.5 4.4 4.7;  ...
4.4 4.2 4.4;  ...
4.5 4.5 4.5;  ...
4.3 4.7 4.4 ];
y = stickleback(:); %column vector
X = [ 1 1 1 1 1 1 1 1    0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0    1 1 1 1 1 1 1 1    0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0    1 1 1 1 1 1 1 1];
X=X';
 [B,~]=regress(y, X)
%     4.3125
%     4.4375
%     4.5875

 n=3; I=eye(3); e=[1;1;1];

C=[1/n*e'; I - 1/n * e*e']
    

 C * B
%     4.4458
%    -0.1333
%    -0.0083
%     0.1417
 %%
X = [ 1 1 1 1 1 1 1 1    1 1 1 1 1 1 1 1    1 1 1 1 1 1 1 1;
      -1 -1 -1 -1 -1 -1 -1 -1    0 0 0 0 0 0 0 0         0 0 0 0 0 0 0 0; ...
       0 0 0 0 0 0 0 0       -1 -1 -1 -1 -1 -1 -1 -1     0 0 0 0 0 0 0 0];
X= X';
 [B,~]=regress(y, X);
 
% B =
% 
%     4.5875  
%     0.2750
%     0.1500
% third mean   4.5875 baseline, subtract  0.2750  and 0.1500
% to otain the firts and second mean. So here 0.2750=mu3-m1  
% and 0.1500=mu3 - mu2.




 
%% 
X = [ 1 1 1 1 1 1 1 1    1 1 1 1 1 1 1 1    1 1 1 1 1 1 1 1;
        0 0 0 0 0 0 0 0    1 1 1 1 1 1 1 1    0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0    1 1 1 1 1 1 1 1];
X=X';
 [B,~]=regress(y, X)

% B =
% 
%     4.3125
%     0.1250
%     0.2750
% first mean is    4.3125, add  0.1250 and 0.2750
% to obtain the second and thirs mean. That is,
% 0.1250= mu2-mu1 and 0.2750=mu3-mu1.
  
 