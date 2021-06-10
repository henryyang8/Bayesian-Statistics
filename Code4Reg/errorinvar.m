clear all
close all


 disp('Predicting Albumin from Plasma Volume: EIV Regression')
 lw = 2; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;


load 'C:\BESTAT\Reg\Regdat\circalbumin.dat'
  %First column:  
  %Second Column:  
  
x = circalbumin(:,1);
y = circalbumin(:,2);
n = length(x);
xbar = mean(x); ybar=mean(y);
Sxx = sum( (x - xbar).^2 );
Syy = sum( (y - ybar).^2 );
Sxy = sum( (x - xbar).* (y - ybar) );

figure(1)
plot(x, y,'o','Markersize',10,...
    'MarkerEdgeColor','k',...
    'Markerface','g')

%traditional regression
b1 = Sxy/Sxx; b0 = ybar - b1 * xbar;
xx = 1700:10:4300;
hold on
plot(xx, b0 + b1 * xx, 'k-','LineWidth',2)

%EIV regression
eta = 200 
bb1= (-(Sxx -  eta * Syy) +  sqrt( (Sxx  -  eta * Syy)^2 ...
+ 4 * eta * Sxy^2 ))/(2 *  eta * Sxy); 
bb0 = ybar - bb1 * xbar;
plot(xx, bb0 + bb1 * xx, 'r-','LineWidth',2)
xlabel('Plasma Volume')
ylabel('Circulating Albumin')
axis tight
yhat = bb0 + bb1 * x;

s2x = 1/(2 * n) * eta/(1 + eta * bb1^2) * sum( (y - yhat).^2 )
s2y = s2x/eta

%print -depsc 'C:\BESTAT\Reg\Regeps\errorinvar.eps'
b0   %-5.7871
b1   %0.0494
bb0  %-13.1619
bb1  %0.0521
s2x  %4852.8
s2y  %24.2641


