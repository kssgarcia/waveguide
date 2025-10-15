function [xhalf, rhalf] = mg_restrict(N,x,r)

%---
% multigrid method: half the grid
%
% position (x), residual(r)
%---

xhalf(1)=x(1);
rhalf(1)=2*r(1)+2*r(2);

Ic=1;
for i=3:2:N-1
 Ic=Ic+1;
 xhalf(Ic)=x(i);
 rhalf(Ic)=r(i-1)+2*r(i)+r(i+1);
end

xhalf(Ic+1)=x(N+1);

%-----
% done
%-----

return
