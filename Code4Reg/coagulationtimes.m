%Box, Hunter, Hunter Coagulation Times Example
%Here 24 animals (Box, Hunter, Hunter; Statistics for Experimenters,
%p. 166) are randomly allocated to 4 different diets,
%but the numbers allocated to different diets are not the same.
%The coagulation time for blood is measured for each animal.
%Are the diet-based differences significant?
% DATA
 times = [62, 60, 63, 59, ...
          63, 67, 71, 64, 65, 66, ...
          68, 66, 71, 67, 68, 68, ...
          56, 62, 60, 61, 63, 64, 63, 59];
 diets = {'dietA', 'dietA','dietA','dietA', ...
'dietB','dietB','dietB','dietB','dietB','dietB',...
'dietC','dietC','dietC','dietC','dietC','dietC',...
'dietD','dietD','dietD','dietD','dietD','dietD','dietD','dietD'};
 [p,table,stats] = anova1(times, diets)
 c = multcompare(stats)
 
 
 %%
 
m = stats.means 
c = [ 1  1 -1 -1 ];
L = c(1)*m(1) + c(2)*m(2)+c(3)*m(3) + c(4)*m(4) %L=-2 , or
LL= m * c'  %LL=-2 
stdL = stats.s * sqrt(c(1)^2/4+c(2)^2/6+c(3)^2/6+c(4)^2/8)
 %stdL = 1.9916 
t = LL/stdL     %t =-1.0042 

 %test H_o: mu * c' = 0  H_1: mu * c' < 0 
 % p-value 
tcdf(t, 23)   %0.1629 

 %or 95% confidence interval for population contrast 
[LL -  tinv(0.975, 23)*stdL, LL +  tinv(0.975, 23)*stdL]
 %   -6.1200    2.1200 

%%
 

[pval stats]=vartestn(times', diets','on') %bartlet 

  %  pval = 0.6441
  %  chisqstat: 1.6680
  %  df: 3

[pval stats]=vartestn(times', diets','on','robust')  %levene 

  %  pval = 0.6237 
  %  fstat: 0.5980 
  %  df: [3 20] 


 
