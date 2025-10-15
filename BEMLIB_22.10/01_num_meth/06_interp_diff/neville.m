function [yint,A] = neville(N,x,f,xint)

%---
% neville interpolation
%---

for i=1:N+1
 u(i)=f(i);
 A(i,1)=u(i);
end

for m=2:N+1
 for i=1:N-m+2
  u(i) =((x(m+i-1)-xint)*u(i)+(xint-x(i))*u(i+1))/(x(m+i-1)-x(i));
%  u(i)=u(i) +(u(i+1)-u(i))/(1+(x(m+i-1)-xint)/(xint-x(i)));
  A(i,m)=u(i); 
 end
end

yint=u(1);

return

