


clear all
close all
 %%%%%%%%%%%%% figure setup %%%%%%%%%%%%%
 lw = 2;  
 set(0, 'DefaultAxesFontSize', 15);
 fs = 15;
 msize = 6;


xlsread C:\BESTAT\Reg\Regdat\vitalcapacity.xlsx;
twos = ans;
  x1 = twos( twos(:,3) > 0, 1);
  y1 = twos( twos(:,3) > 0, 2);
  x2 = twos( twos(:,3) ==0, 1);
  y2 = twos( twos(:,3) ==0, 2);
n1=length(x1); n2 = length(x2);
SXX1 = sum((x1 - mean(x1)).^2)  %4.3974e+003
SXX2 = sum((x2 - mean(x2)).^2)  %6.1972e+003
SYY1 = sum((y1 - mean(y1)).^2)  %26.5812
SYY2 = sum((y2 - mean(y2)).^2)  %20.6067
SXY1 = sum((x1 - mean(x1)).*(y1 - mean(y1))) %-236.3850
SXY2 = sum((x2 - mean(x2)).*(y2 - mean(y2))) %-189.7116
b1_1 = SXY1/SXX1   %-0.0538
b1_2 = SXY2/SXX2   %-0.0306

SSE1 = SYY1 - (SXY1)^2/SXX1  %13.8741
SSE2 = SYY2 - (SXY2)^2/SXX2  %14.7991

s2 = (SSE1 + SSE2)/(n1 + n2 - 4)  %0.3584
s = sqrt(s2)  %0.5987
seb1b2 = s * sqrt( 1/SXX1 + 1/SXX2 ) %0.0118
t = (b1_1 - b1_2)/seb1b2     %-1.9606
pval = tcdf(t, n1 + n2 - 4)   %0.0267


s22 = (SYY1 + SYY2 - ...
 (SXY1 + SXY2)^2/(SXX1 + SXX2))/(n1 + n2 - 3) %0.3710
s = sqrt(s22)  %0.6091
seb1b2 = s * sqrt( 1/SXX1 + 1/SXX2 ) %0.0120
t = (b1_1 - b1_2)/seb1b2     %-1.9270
pval = tcdf(t, n1 + n2 - 3)   %0.0287

%============== graphics =======================================  
  figure(1)
  plot( x1, y1, 'ro' )
  hold on
  dkg =[0 0.6 0];
  plot( x2, y2, 'o' , 'color', dkg)
  xx = 15:70;
  plot(xx, polyval(polyfit(x1, y1, 1), xx),'r-','linewidth',lw)
  plot(xx, polyval(polyfit(x2, y2, 1), xx),'-','color', dkg,'linewidth',lw)
  xlabel('age (years)'); ylabel('VC (liters)');
  axis([10 75 2.5 6])
  legend('exposed','unexposed')
  %print -depsc 'C:\Springer\Reg\Regeps\twoslopes.eps'
  