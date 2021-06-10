%clf; 
 clear all;
 close all; 
 
 
 disp('Hubble Regression')
  lw = 2; 
  set(0, 'DefaultAxesFontSize', 17);
  fs = 17;
  msize = 7;
  
distance = [.032  .034   .214   .263   .275 .275 ...
       .45   .5     .5     .63    .8   .9   ...
       .9    .9    .9     1.0    1.1  1.1   ...
      1.4   1.7   2.0     2.0   2.0   2.0 ]';

velocity = [170   290   -130   -70   -185   -220 ...
      200    290    270   200    300    -30 ...
      650    150    500   920    450    500 ...
      500    960    500   850    800   1090 ]';

%mod = x2fx(distance, [0; 1]);
      
       s = regstats(velocity,distance,[0; 1],{'yhat','r','beta'});
       figure(1)
       plot(s.yhat,s.r,'o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',8) 
            hold on
        plot([0 1000],0.*[0 1000], 'k--','linewidth',2)
       xlabel('Fitted Values'); ylabel('Residuals');
       axis tight
       %
        % print -depsc 'C:\Springer\Reg\Regeps\hubbleresid.eps'

       figure(2)
       plot(distance, velocity, 'o', 'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',8) 
       hold on
       b = regress(velocity,distance);
       di=[0  2.1];
       plot(di, b*di, 'r-','linewidth',2)
       plot(di, s.beta(1) + s.beta(2)*di, 'b-','linewidth',2)
       axis tight
       ylabel('Velocity (km/s)')
       xlabel('Distance (mps)') 
       %  print -depsc 'C:\Springer\Reg\Regeps\hubblescatter.eps'

       
%        'Q'           Q from the QR Decomposition of the design matrix
%        'R'           R from the QR Decomposition of the design matrix
%        'beta'        Regression coefficients
%        'covb'        Covariance of regression coefficients
%        'yhat'        Fitted values of the response data
%        'r'           Residuals
%        'mse'         Mean squared error
%        'rsquare'     R-square statistic
%        'adjrsquare'  Adjusted R-square statistic
%        'leverage'    Leverage
%        'hatmat'      Hat (projection) matrix
%        's2_i'        Delete-1 variance
%        'beta_i'      Delete-1 coefficients
%        'standres'    Standardized residuals
%        'studres'     Studentized residuals
%        'dfbetas'     Scaled change in regression coefficients
%        'dffit'       Change in fitted values
%        'dffits'      Scaled change in fitted values
%        'covratio'    Change in covariance
%        'cookd'       Cook's distance
%        'tstat'       t statistics and p-values for coefficients
%        'fstat'       F statistic and p-value
%        'dwstat'      Durbin-Watson statistic and p-value
%        'all'         Create all of the above statistics


% R has units of distance we discover that Ho has units of 1/time. Inverting the Hubble constant gives us an approximate age of the universe:
% 
% Ho = 65 km.s-1.Mpc-1
% 
% 1 Mpc = 3.26 million light years = 30841981340613408000km
% 
% Ho = 65 km.s-1.Mpc-1 / 30841981340613408000km/Mpc
% 
% 1/Ho = 30841981340613408000/65 seconds
% 
% This gives an age for the universe of about 15 billion years. More accurate determinations of the constant would give more accurate universal ages - and vice versa.


1/424 * 3.086 * 10^19 /(60*60*24*365) %2.3079e+009
  

% The slope of the fitted line is 464 km/sec/Mpc,
% and is now known as the Hubble constant, Ho. 
%     Since both kilometers and Megaparsecs (1 Mpc = 3.086E19 km) 
%     are units of distance, the simplified units of Ho are 1/time,
%     and the conversion is given by
% 1/Ho = (978 Gyr)/(Ho in km/sec/Mpc)
% Thus Hubble's value is equivalent to approximately 2 Gyr.
% Since this should be close to the age of the Universe, 
% and we know (and it was known in 1929) that the age of 
% the Earth is larger than 2 billion years, Hubble's value
% for Ho led to considerable skepticism about cosmological 
%     models, and motivated the Steady State model. However, 
%     later work found that Hubble had confused two different
%     kinds of Cepheid variable stars used for calibrating distances, 
%     and also that what Hubble thought were bright stars in distant
%     galaxies were actually H II regions. Correcting for these 
%     errors has led to a lowering of the value of the Hubble
%     constant: there are now primarily two groups using Cepheids:
%     the HST Distance Scale Key Project team 
%     (Freedman, Kennicutt, Mould etal) which gets 72+/-8 km/sec/Mpc, 
%     while the Sandage team, also using HST observations of
%         Cepheids to calibrate Type Ia supernovae, gets 57+/-4 km/sec/Mpc.
%         Other methods to determine the distance scale include the time
%         delay in gravitational lenses and the Sunyaev-Zeldovich effect i
%         n distant clusters: both are independent of the Cepheid calibration 
%         and give values consistent with the average of the two HST groups: 
%         65+/-8 km/sec/Mpc. These results are consistent with a combination 
%         of results from CMB anisotropy and the accelerating expansion of 
%         the Universe which give 71+/-3.5 km/sec/Mpc. With this value for 
%         Ho, the "age" 1/Ho is 14 Gyr while the actual age from the 
%         consistent model is 13.7+/-0.2 Gyr.
% 


% Intrinsic Redshifts and the Hubble Constant
% M.B. Bell1 and S.P. Comeau1
% ABSTRACT
% We show that the VCMB velocities of the Fundamental Plane (FP) clusters studied in the
% Hubble Key Project appear to contain the same discrete "velocities" found previously by us
% and by Ti to be present in normal galaxies. Although there is a particular Hubble constant
% associated with our findings we make no claim that its accuracy is better than that found by the
% Hubble Key Project. We do conclude, however, that if intrinsic redshifts are present and are not
% taken into account, the Hubble constant obtained will be too high.
% 
% 
% Cluster/Group D (Mpc) VCMB (km s?1)a Transit.(Disc.Vel.)(km s?1)b VH (km s?1)c
% Dorado 13.8 1131 ziG[1,7](145.2) 986
% Grm 15 47.4 4530 ziG[1,4](1157.9) 3372
% Hydra 49.1 4061 ziG[1,5](580.1) 3481
% Abell S753 49.7 4351 ziG[2,6](725.2) 3626
% Abell 3574 51.6 4749 ziG[1,4](1157.9) 3591
% Abell 194 55.9 5100 ziG[1,4](1157.9) 3942
% Abell S639 59.6 6533 ziG[1,3](2314.3) 4219
% Coma 85.8 7143 ziG[1,4](1157.9) 5985
% Abell 539 102.0 8792 ziG[2,7](1448.6) 7343
% DC 2345-28 102.1 8500 ziG[1,4](1157.9) 7342
% Abell 3381 129.8 11536 ziG[1,3](2314.3) 9222
D =[13.8 47.7 49.1 49.7 51.6 55.9 59.6 85.8 102.0 102.1 129.8]; %Mpc
VH=[986 3372 3481 3626 3591 3942 4219 5985 7343 7342 9222]; %km/sec
figure(3)
scatter(D, VH)
H = mean(D .* VH)/mean(VH) %82.6944
1/H * 3.086 * 10^19 /(60*60*24*365)   % 1.1833e+010 
% The universe is 11.8 billion years old