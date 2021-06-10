%
% 
% 
% Duncan J. S., Schofield G.,  Duncan E. K., Hinckson E. A. (2007).
% Effects of Age, Walking Speed, and Body Composition on Pedometer Accuracy in Children
% {\it Research Quarterly for Exercise and Sport}, {\bf 78}, 5, 420--428.
% 
% %Duncan et al (2007) investigated the effects of age group, walking speed, 
% and body composition on the accuracy of pedometer-determined step counts in children.
% Eighty-five participants (43 boys, 42 girls), ages 5-7 and 9-11 years, walked on a treadmill 
% for two-minute bouts at speeds of 42, 66, and 90 m/min  while wearing a spring-levered
%     (Yamax SW-200) and a piezoelectric (New Lifestyles NL-2000) pedometers. The number 
%     of steps taken during each bout was also recorded using a hand counter. 
% Body mass index (BMI) was calculated from height and mass, and percentage of body fat (\%BF) was 
% determined using hand-to-foot bioelectrical impedance analysis. The tilt angle of the
% pedometer was assessed using a magnetic protractor. Both pedometers performed
% well at 66 and 90 m/min, but undercounted steps by approximately 20 \% 
% at 42 m/min. Although age group, BMI, waist circumference, and \%BF did 
% not affect pedometer accuracy, children with large pedometer tilt angles ($> 10^\circ$)
% showed significantly greater percent bias than those with small tilt angles ($< 10^\circ$).
% Findings:
% We suggest that the style of waistband on the child's clothing is a more important
% determinant of tilt angle and thus pedometer accuracy than body composition.
% Our results also indicate that the NL-2000 pedometer provides similar accuracy 
% and better precision than the SW-200 pedometer, especially in children with large
% tilt angles. We conclude that fastening pedometers to a firm elastic belt may 
% improve stability and reduce undercounting in young people.
% 
% 
% 


clear all
close all
%
load 'C:\BESTAT\Reg\Regdat\pmr1.mat'  %data set
diffs = pmr1(:,1); 
X = [ones(length(diffs),1)  pmr1(:, 2:end) ];

% 
[b,bint,res,resint,stat] = regress(diffs,X)
rcoplot(res,resint)


%  plot
x1=X(:,2); %tilt
x2=X(:,4); %BMI
scatter3(x1,x2,diffs,'filled')
hold on
  x1fit = min(x1):0.01:max(x1);
  x2fit = min(x2):0.01:max(x2);
  [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
  YFIT = b(1) + b(2)*X1FIT + b(4)*X2FIT ;
  mesh(X1FIT,X2FIT,YFIT);
  xlabel('Tilt')
  ylabel('BMI')
  zlabel('diff')
  view(40,30)
hold off
