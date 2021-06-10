
% Breaking Strength of Medical Fibers
% 
% Three machines are producing monofilament sinthetic fibre for
% medical use (surgery, implants, devices, etc).
% 
% The measured response is strength y (in lb) and a covariate
% is the diameter x (in $10^{-3}$ inches)
% 
clear
close all force
%
y = [36 41 39 42 49   40 48 39 45 44     35 37 42 34 32];
x = [20 25 24 25 32   22 28 22 30 28     21 23 26 21 15];
g = [1  1   1  1  1    2  2  2  2  2      3  3  3  3  3];

plot(x(g==1),y(g==1),'s','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
hold on
plot(x(g==2),y(g==2),'s','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10)
plot(x(g==3),y(g==3),'s','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10)
            legend('1','2','3','Location','NorthWest')
hold off


[p,table,stats]=anova1(y, g)
multcompare(stats)
%%

%%
aoctool(x, y, g)

%% 

%%
%=======================
a =3; n=5;

x1=x(g==1);  x2=x(g==2); x3=x(g==3);
y1=y(g==1);  y2=y(g==2); y3=y(g==3);
x1b=mean(x1); x2b=mean(x2); x3b=mean(x3); 
y1b=mean(y1); y2b=mean(y2); y3b=mean(y3);
xb = mean(x); yb = mean(y);

SXX = sum( (x - xb).^2 );
SXY = sum( (x - xb).* (y - yb) );
SYY = sum( (y - yb).^2 );

QXX = sum( (x1-x1b).^2 + (x2-x2b).^2 + (x3-x3b).^2 );
QXY = sum( (x1-x1b).* (y1 - y1b) + ...
    (x2-x2b).*(y2 - y2b) + (x3-x3b).* (y3-y3b) );
QYY = sum( (y1-y1b).^2 + (y2-y2b).^2 + (y3-y3b).^2 );


TXX = SXX - QXX;
TXY = SXY - QXY;
TYY = SYY - QYY;

SSE = QYY - QXY^2/QXX;
MSE = SSE/(a * (n-1) -1 );
SSEp = SYY - SXY^2/SXX;
b = QXY/QXX    %0.9540
a1 = y1b - yb  - b * (x1b - xb);
a2 = y2b - yb  - b * (x2b - xb);
a3 = y3b - yb  - b * (x3b - xb);
[a1 a2 a3]   %0.1824    1.2192   -1.4016

F1 = ((SSEp - SSE)/(a-1)) /MSE   %2.6106

F2 = (QXY^2/QXX)/MSE    %69.9694

xx = 15:0.2:33;

plot(x1, y1, 'ro','markersize',7, 'MarkerFaceColor','g')
hold on
plot(xx, a1+yb  + b*(xx - xb), 'r-')
plot(x2, y2, 'ko','markersize',7, 'MarkerFaceColor','r')
plot(x3, y3, 'ko','markersize',7, 'MarkerFaceColor','b')
plot(xx, a2+yb + b*(xx-xb), 'g-')
plot(xx, a3+yb + b*(xx-xb), 'b-')
%%

close all
