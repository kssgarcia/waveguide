function [xdouble, edouble] = mg_prolongate(N,x,e)

%---
% multigrid method: double the grid
%
% prolongate the position (x), error(e)
%---

xdouble(1)=x(1);
edouble(1)=e(1);

Ic=1;
for i=2:N+1
 Ic=Ic+1;
 xdouble(Ic)=0.5*(x(i)+x(i-1));
 edouble(Ic)=0.5*(e(i)+e(i-1));
 Ic=Ic+1;
 xdouble(Ic)=x(i);
 edouble(Ic)=e(i);
end

%---
% done
%---

return
