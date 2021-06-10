
clear all
close all
% Von Willebrand Disease is a bleeding disorder caused by a defect or deficiency of a blood clotting protein,
% called von Willebrand Factor. This glue-like protein, produced by the cells that line the blood vessel walls,
% interacts with blood cells called platelets to form a plug which prevents bleeding.
% In order to understand the differential bonding mechanics
% underlying von Willebrand type bleeding disorders, researchers in Dr. Larry McIntire Lab at Georgia Tech
% studied the interactions between
% the wild-type platelet GPIba molecule (receptor) and wild-type von
% Willebrand Factor (ligand).
% 
%  The mean stop time rolling parameter
%  was calculated from frame-by-frame rolling velocity data collected at 250 frames per second.
%  Mean stop time indicates the amount of time a cell spends stopped, so it is analogous
%  to the bond lifetime.  This parameter,  being an indicator for how
%   long the bond is stopped before the platelet moves again,
%   can be used to assess bond lifetime and off-rate (Yago et al., 2008).
% 
% For the purpose of exploring interactions between the force and
% the mean  stop times, Ficoll 6\% is added to increase the viscosity.
% 
% Data are courtesy of Dr. Leslie Coburn.
% The mat-file {\tt coburn.mat} contains the structure {\tt coburn} with data vectors
% {\tt coburn.fxssy} where {\tt x = 0, 6} is a code for Ficoll absence/presence and
% {\tt y = 1, 2, 4, ..., 256} denotes the shear stress (in $dyn/cm^2$).
% For example {\tt coburn.f0ss16} is a $243 \times 1$ vector of mean stop times
% obtained with no Ficoll, under shear stress of 16 $dyn/cm^2.$
% 
% {\small
% \begin{center}
% \begin{tabular}{ccccccccc}
% \hline
% Shear Stress      & 2 & 4 & 8 & 16 & 32 & 64 & 128 & 256 \\
% Shear Number      & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 \\
% Mean Stop Time  &f0ss2& f0ss4&  f0ss8 &f0ss16& f0ss32& f0ss64 &f0ss128 &f0ss256 \\
% Sample Size         &26& 57 & 157 & 243 & 256 & 185 & 62 & 14 \\
% \hline
% \end{tabular}
% \end{center}
% }
% 
% 
% (a) Fit a regression on the logarithm (base 2) of mean stop time, with no Ficoll present, as a quadratic function of
% share stress number ($log2(dyn/cm^2)$).
% This regression is linear (in parameters) with two predictors, \ttmat{ shear} and the squared shear, \ttmat{ shear2}.
% Figure \ref{fig:coburnreg} is plotted by \matlabico \ttmat{coburnreg.m}.
% 


 disp('Mean Stop Time in Von Willebrand Factor')
 lw = 3; 
 set(0, 'DefaultAxesFontSize', 17);
 fs = 16;
 msize = 10;


load 'C:\BESTAT\Reg\Regdat\coburn.mat';

mst=[coburn.f0ss2;  coburn.f0ss4;  coburn.f0ss8;   coburn.f0ss16; ...
     coburn.f0ss32; coburn.f0ss64; coburn.f0ss128; coburn.f0ss256];

shearn = [1 * ones( 26,1); 2 * ones( 57,1); 3 * ones(157,1); ...
          4 * ones(243,1); 5 * ones(256,1); 6 * ones(185,1); ...
          7 * ones( 62,1); 8 * ones( 14,1)];

shearn2 = shearn.^2;  %quadratic term

% design matrix
X = [ones(length(shearn),1) shearn  shearn2];
[b,bint,res,resint,stats] = regress(log2(mst), X) ; 
   % b0=-6.2423, b1=0.8532, b2 =-0.0978
 
xx = 0:0.01:9;
yy = -6.2423 +  0.8532*xx  - 0.0978 * xx.^2;
    
close all;
figure(1)
plot(shearn, log2(mst), '*')
hold on
plot(xx, yy, 'r-','linewidth',2)
xlabel('$\log_2(dyn/cm^2)$','interpreter','latex')
ylabel('$\log_2(mst)$','interpreter','latex')
axis([0 9 -7 -1])
%print -depsc 'C:\Springer\Reg\Regeps\coburnreg1.eps'

figure(2)
histn(res, -2,0.25,2)
hold on
std(res) %  0.7645
plot( (-2:0.01:2), normpdf((-2:0.01:2), mean(res), std(res)), 'r-','linewidth',3)
axis([-3 3 0 0.6])
%print -depsc 'C:\Springer\Reg\Regeps\coburnreg2.eps'
