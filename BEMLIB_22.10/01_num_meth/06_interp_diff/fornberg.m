%================================
% Fornberg's interpolation method
%================================

N=2;

x = [0.0, 0.5, 1.0];
xint = 0.5;
m = 2;

for k=1:m+1
 for j=1:N+1
  for i=1:N+1
   c(k,j,i)=0.0;
  end
 end
end

%---
c(1,1,1) = 1;
alpha = 1.0;
%---

for i=2:N+1

 beta = 1.0;
 p = min(i,m)+1;

 for j=1:i
  beta = beta*(x(i)-x(j));
  for k=1:p
   c(k,i,j) = ((x(i)-xint)*c(k,i-1,j)-k*c(k-1,i-1,j))/(x(i)-x(j));
  end
 end

 for k=1:p
   c(k,i,j) = alpha/beta *(k*c(k-1,i-1,j-1)-(x(i-1)-xint)*c(k,i-1,j-1));
 end

 alpha = beta;

end
%---
