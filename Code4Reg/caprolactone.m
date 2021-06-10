   close all force
   clear all


data=[
    9.71  24.276  
    9.71  24.083  
    9.71  24.276
    8.52  20.206   
    8.52  20.199  
    8.52  20.223
    7.96 19.773   
    7.96 19.759  
    7.96 19.765
    6.82 16.743  
    6.82 16.587  
    6.82 16.744
    5.85 15.081  
    5.85 15.121   
    5.85 15.274
    4.95 12.636  
    4.95 12.641 
    4.95 12.682
     3.91 9.869 
     3.91 9.906 
     3.91 9.883
     2.98 7.624
     2.98 7.592 
     2.98 7.585
     2.07 4.638
     2.07 4.666 
     2.07 4.649
     1.02 2.860
     1.02 2.859  
     1.02 2.896
     0.00 0.000 
     0.00 0.000
     0.00 0.000];
 x = data(:,1); y=data(:,2); 
 n=length(x)  %33
SXX = sum( (x - mean(x)).^2 ) %305.8026
SXY = sum( (x - mean(x)).* (y - mean(y)) ) %749.4972
b1 = SXY/SXX %2.4509
b0 = mean(y) - b1 * mean(x)  %0.1694
s=sqrt(sum(  (y - (b0+b1*x)).^2 )/(n-2)) %0.4216
% 
ystars=[1.582  1.793 1.787]; m=length(ystars);
ystar = mean(ystars);
xbar=mean(x);
xstar = 1/b1 * (ystar - b0)  %0.6329
LB=xstar - tinv(0.975, n-2) * s/b1 * sqrt(1/m+1/n+(xstar-xbar)^2/SXX);
UB=xstar + tinv(0.975, n-2) * s/b1 * sqrt(1/m+1/n+(xstar-xbar)^2/SXX);
[LB xstar UB] %[0.4047 0.6329 0.8611]

% figure(1)
% plot(x,y,'*')
% hold on
% plot(x,b1*x + b0,'r-')
% plot(x, ystars(1) * ones(length(x)), 'g-')
% plot(x, ystars(2) * ones(length(x)), 'g-')
% plot(x, ystars(3) * ones(length(x)), 'g-')
% plot([xstar xstar],[0 25],'r-')
