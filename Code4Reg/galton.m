close all
clear all
 disp('Galton: Mid-Parents - Sons')
 lw = 2; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 9;
 % 

figure(1)
load 'C:\STAT\Reg\Regdat\galton.dat'
midparents = galton(:,1);
  % midparent = (father + 1.08 * mother)/2
sons = galton(:,2);
scatterhist(midparents, sons) %,'ks','MarkerSize',12)
hold on
plot([55,80], [55,80],'k-','LineWidth',2)
a = regress(sons, [ones(size(midparents))  midparents] );
reglin = @(x) a(1) + a(2) .* x;
plot(55:0.1:80,  reglin(55:0.1:80),'r-','LineWidth',2)
xlabel('Height of Midparents (in)')
ylabel('Height of Sons (in)')



%%
% Pearson's Version of Galton's Data
clear all
%close all
  disp('Galton: Fathers - Sons')
  lw = 3; 
  set(0, 'DefaultAxesFontSize', 16);
  fs = 15;
  msize = 10;

figure(2)
load 'C:\BESTAT\Reg\Regdat\pearson.dat'
fathers = pearson(:,1);
sons = pearson(:,2);
%plot(fathers, sons,'o','MarkerSize',8)
scatterhist(fathers, sons) %,'ks','MarkerSize',12)
hold on
plot([55,80], [55,80],'k-','LineWidth',2)
a = regress(sons, [ones(size(fathers))  fathers] );
reglin = @(x) a(1) + a(2) .* x;
plot(55:0.1:80,  reglin(55:0.1:80),'r-','LineWidth',2)
plot(55:0.1:80,  1.0139 * (55:0.1:80),'g-','LineWidth',2)
xlabel('Height of Fathers (in)')
ylabel('Height of Sons (in)')
legend('(father, son)','y=x line','regr line','regr line int=0',2)
%print -depsc 'C:\Springer\Reg\Regeps\galtonfs.eps'


