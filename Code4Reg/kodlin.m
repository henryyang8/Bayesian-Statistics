close all
clear all
% Kodlin (1951) reported an experiment that compared two 
% substances for lowering blood pressure, denoted as 
% Substances A and B. Two groups of animals are randomized
% to the two substances and a decrease in pressure is recorded.
% The initial pressure is recorded.  Compare the two
% substances by accounting for possible effect of the initial
% pressure on the decrease in pressure. Discuss the results 
% of the ANCOVA analysis and compare them with those obtained
% by ignoring the potential effect of initial pressure.
% 
x = [135 125 125 130 105 130 140  93 110 100 ...
      90 135 130 115 110 140 130  95  90 105];
y = [45 45 20 50 25 37 50 20 25 15 ...
     34 55 50 45 30 45 45 23 40 35];
g = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2];

a = 2; n = 10;

x1=x(g==1);  x2=x(g==2); 
y1=y(g==1);  y2=y(g==2); 
x1b=mean(x1); x2b=mean(x2);  
y1b=mean(y1); y2b=mean(y2); 
xb = mean(x); yb = mean(y);

SXX = sum( (x - xb).^2 );
SXY = sum( (x - xb).* (y - yb) );
SYY = sum( (y - yb).^2 );

QXX = sum( (x1-x1b).^2 + (x2-x2b).^2  );
QXY = sum( (x1-x1b).* (y1 - y1b) + ...
    (x2-x2b).*(y2 - y2b)  );
QYY = sum( (y1-y1b).^2 + (y2-y2b).^2  );
% estimators in the model parameters mu, alpha_i, beta
mu = yb
b = QXY/QXX     % 0.5175
alpha = [y1b y2b] - [yb yb] - b*([x1b x2b]-[xb xb])     %-4.8714    4.8714
intercept = mu - b*xb   %  -23.6659
% two group intercepts  -23.6659-4.8714= -28.5373  
%                                    -23.6659+4.8714= -18.7946

TXX = SXX - QXX;
TXY = SXY - QXY;
TYY = SYY - QYY;

SSE = QYY - QXY^2/QXX;
MSE = SSE/(a * (n-1) - 1 );
SSEp = SYY - SXY^2/SXX;

% F test for testing H_0: alpha_i = 0 (all i)
F1 = ((SSEp - SSE)/(a-1)) /MSE   %   7.6321
pvalF1 = 1- fcdf(F1,a-1, a*(n-1) - 1)  % 0.0133

% F test for testing H_0: beta = 0
F2 = (QXY^2/QXX)/MSE    %  24.5668
pvalF2 = 1 - fcdf(F2, 1, a * (n-1)-1 )  % 1.2003e-04

  
plot(x1, y1, 'ro','markersize',7, 'MarkerFaceColor','r')
hold on
plot(x2, y2, 'ko','markersize',7, 'MarkerFaceColor','k')
legend('A','B',2)
plot(x, mu + alpha(1) + b*(x - xb ),'r-', 'linewidth',2)
plot(x, mu + alpha(2) + b*(x - xb ),'k-', 'linewidth',2)
hold off

%print -depsc 'C:\Springer\Reg\Regeps\ancovakodlin.eps'

%%
aoctool(x, y, g)    %intercepts p=0.0133; slope p=0.0001

%% 
% Ignore x ... no significant differences.
y = [45 45 20 50 25 37 50 20 25 15 ...
     34 55 50 45 30 45 45 23 40 35];
g = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2];

y1=y(g==1); y2=y(g==2);
n1=length(y1); n2=length(y2); %10  10
s1=std(y1); s2=std(y2); %13.6284    9.7160
sp = sqrt( ((n1-1)* s1^2 + (n2-1)*s2^2)/(n1+n2-2));  %11.8350
t = (mean(y1)-mean(y2))/(sp * sqrt( 1/n1 + 1/n2)) %-1.3226
%Pval for onesided alternative H1: mu1 < mu2
pval = tcdf(t, n1+n2-2)   %0.1013

