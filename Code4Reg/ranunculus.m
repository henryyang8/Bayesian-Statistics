%ranunculus.m
clear all
close all
% 
% Data by MacLeod, Julius (1901)

% Weldon, W. F. R. (1901). Chage in organic correlation
% of {\it Ficaria Rununculus } during the flowering season.
% {\it Biometrika}, {\bf 1}, 1, 125--128.

 disp('Ranunculus Ficaria: Nymber of Pistils & Stamens')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 17);
 fs = 16;
 msize = 10;


load 'C:\BESTAT\Reg\Regdat\ranunculus.mat';
[m n]=size(ranunculus);
rexp=[];  %ranunculus expanded...
for i = 1:m
    re=ranunculus(i,:);
 rexp = [rexp;...
        repmat( [re(1)  re(2) re(3)],[re(4), 1])];
end

% 0 in the first column: Early Flowers
figure(1)
%cumulative split for scaterplots
ranun0 =  ranunculus(ranunculus(:,1)==0,:);
scatter(ranun0(:,2), ranun0(:,3),  10*ranun0(:,4))
% expanded split for coefficient of correlation
rexp0 = rexp(rexp(:,1)==0,:);
corrcoef(rexp0(:,2), rexp0(:,3))

% 1 in the first column: Late Flowers
figure(2)
ranun1 =  ranunculus(ranunculus(:,1)==1,:);
rexp1 = rexp(rexp(:,1)==1,:);
scatter(ranun1(:,2), ranun1(:,3),  10*ranun1(:,4))
corrcoef(rexp1(:,2), rexp1(:,3))




