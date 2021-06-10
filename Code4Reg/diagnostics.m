function [index,res,stud_res,lev,DFFITS,Cooks_D,DFBETAS]=diagnostics(y,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [index,res,stud_res,lev,DFFITS,Cooks_D,DFBETAS]=diagnostics(y,Z)
%
% y is an (n x 1) vector of observations on the dependent variable
% Z is an (n x r) design matrix of explanatory variables without intercept
% The output:
%
%    index    : the number of the observation (index from 1 to n)
%    res      : the residuals 
%    stud_res : the studentized residuals 
%    lev      : leverage-value 
%    DFFITS   : DFFITS 
%    Cook_D   : Cook´s distance 
%    DFBETAS  : DFBETAS (for each the explanatory variables in the model) 
%    Example:
%    y = sqrt([1:10]);
%    y = y' + 0.3 * randn(10,1)
%    Z = [y + 0.4 *randn(10,1)    y + 0.5 *randn(10,1)   y + 0.4 *randn(10,1) ]
%     [index,res,stud_res,lev,DFFITS,Cooks_D,DFBETAS]=island(y,Z)
% AUTHOR:
% Birgir Hrafnkelsson, Ph.D.
% Associate Research Professor
% Applied Mathematics Division
% Science Institute
% University of Iceland
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Pre-calculations

n=size(Z,1);
r=size(Z,2);
p=r+1;
df=n-p;

% Pre-cacluations for Regression 

syy=(y-mean(y))'*(y-mean(y));
Z0=[ones(n,1) Z];
A=Z0'*Z0;
Ainv=inv(A);
ckk=diag(Ainv);
b=Z0'*y;
beta0=A\b;
yhat=Z0*beta0;
res=y-yhat;
SSres0=res'*res;
s2full=SSres0/(n-r-1);

% Calculations start

index=(1:n)';
H=Z0*Ainv*Z0';
lev=diag(H);
stud_res=res./sqrt(s2full*(1-lev));
DFFITS=res.*sqrt((n-r)./(SSres0*(1-lev)-res.^2)).*sqrt(lev./(1-lev));
Cooks_D=res.^2/(p*s2full).*(lev./(1-lev).^2);

% Calculations for DFBETAS

DFBETAS=zeros(n,p);

for j=1:n
    indexZ=[(1:j-1) (j+1:n)];
    Z0=[ones(n-1,1) Z(indexZ,:)];
    A=Z0'*Z0;
    Ainv=inv(A);
    ckk=diag(Ainv);
    b=Z0'*y(indexZ);
    betaj=A\b;
    yhat=Z0*betaj;
    e=y(indexZ)-yhat;
    MSE=(e'*e)/(n-1-p);    
    DFBETAS(j,:)=((beta0-betaj)./sqrt(MSE.*ckk))';
end
