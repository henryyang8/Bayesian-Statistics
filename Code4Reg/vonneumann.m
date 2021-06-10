

%von Neumann Ratio Statistic
%Q = @(x)  sum(diff(x).^2)/(2 * sum( (x -mean(x)).^2 ))

Q = @(x)  sum(diff(x).^2)/(2 * sum( (x).^2 ))

n=100;
randn('state',123);  x = randn(1, (n+1));
y = x(1:n);
z = x(1:n) + 0.2*x(2:(n+1));  %positive  autocorr
w =  x(1:n) - 0.2*x(2:(n+1)); %negative autocorr

q1=Q(y)   % 1.0471
q2=Q(z)   % 0.8648
q3=Q(w)  % 1.2215

 pval = @(Q,n) normcdf( -(1-Q) * sqrt( (2*n + 1)/(2 - (1-Q)^2)) ) 
 
 2*min(pval(q1,100) , 1-pval(q1,n) )     %0.6819
 pval(q2, n)      %0.0866
 1-pval(q3, n)  %0.0123
 
 p=dwtest(y', ones(length(y),1), 'approximate','both')

 p=dwtest(z', ones(length(z),1), 'approximate','right')
 
  p=dwtest(w', ones(length(w),1), 'approximate','left')
 
  
  %ABBE TEST
n=100;
randn('state',123);  
x = randn(1, n);
y = x(1:n) + sin(3.14159215*(1:n)./n);
z = x(1:n) + 0.2 * repmat([1 -1], 1, n/2)

q1=Q(x)    
q2=Q(y)  
q3=Q(z)   


