% Reference from Brani's Engineering Biostatistics textbook ch.14
% Data reference: https://www.kaggle.com/mirichoi0218/insurance
% Description of data: (col 1-5: pred, col 6: response)
% Column 1: age
% Column 2: sex -> male: 0 female: 1
% Column 3: bmi
% Column 4: number of children
% Column 5: smoker -> no: 0 yes: 1
% Column 6: insurance charges (in ten thousands)
clear;

data = readtable("insurance.csv");
head(data,3)
x1 = table2array(data(:,1));
x2 = table2array(data(:,2));
x3 = table2array(data(:,3));
x4 = table2array(data(:,4));
x5 = table2array(data(:,5));
intercept = ones(size(x1));
Z = [x1 x2 x3 x4 x5];
X = [intercept Z];
Y = table2array(data(:,6));
[n,p] = size(X) %n=200 p=6
b = inv(X' *X)* X'* Y; %intercept=-1.2041 b1=0.0282 b2=-0.1356 b3=0.0305 b4=0.0442 b5=2.3387
H = X * inv(X' * X) * X';
max(max(H * H - H));
Yhat = H * Y;
J = ones(n); 
I = eye(n);
SSR = Y' * (H - 1/n * J) * Y;
SSE = Y' * (I - H) * Y;
SST = Y' * (I - 1/n * J) * Y; 
MSR = SSR/(p-1) %47.3422
MSE = SSE/(n-p) %0.3459
F = MSR/MSE %136.8589
pval = 1-fcdf(F, p-1, n-p) %0
Rsq = 1 - SSE/SST %0.7791
Rsqadj = 1 - (n-1)/(n-p) * SSE/SST %0.7734
s = sqrt(MSE)  %0.5881

sig2 = MSE * inv(X' * X); 
sb=sqrt(diag(sig2)); 
tstats = b./sb; 
pvals = 2 * tcdf(-abs(tstats), n-p);

Xh=[1 30 0 30 1 1]; %assume 30 yr old male smoker with 30 bmi & 1 child
Yh = Xh * b %2.9410
sig2h = MSE * Xh * inv(X' * X) * Xh';
sig2hpre = MSE * (1 + Xh * inv(X' * X) * Xh');
sigh = sqrt(sig2h); 
sighpre = sqrt(sig2hpre); 
%95% credible sets on the mean and predicted responses 
[Yh-tinv(0.975,n-p)*sigh, Yh+tinv(0.975,n-p)*sigh] %[2.7522  3.1298]
[Yh-tinv(0.975,n-p)*sighpre,Yh+tinv(0.975,n-p)*sighpre] %[1.7657  4.1162]

%Diagnosis: prediction of yi with ith observation removed 
ind = 1:n; 
Yhati = []; 
for i = 1:n
    indi = find(ind ~= i); 
    Yi = Y(indi); 
    Xi=X(indi,:); 
    bi = inv(Xi' * Xi) * Xi'* Yi; 
    Yhatii = X(i,:) * bi; 
    Yhati =[Yhati; Yhatii];
end
Yhati %prediction of yi without ith observation

hii = diag(H); %leverages 
resid = (I - H)*Y; %ordinary residuals
sresid = sqrt(MSE .* (1-hii)); 
stresid = resid./sresid %studentized residuals
di = Y - Yhati;
PRESS = sum(di.^2) 
R2pred = 1-PRESS/SST %R^2 predictive
sminusi = sqrt(((n-p)*MSE*ones(n,1) - resid.^2./(1-hii))/(n-p-1));
ti = resid ./(sminusi .* sqrt(1-hii)) 
outli=hii/mean(hii); 
find(outli > 3) %None

%influential observations 
CooksD = resid.^2 .* (hii./(1-hii).^2)/(p * MSE) 
%influence if ith to all; 
find(CooksD > 4/n) %4 10 35 59 63 70 86 99 100 103 116 117 139 141 157 162 176 186



