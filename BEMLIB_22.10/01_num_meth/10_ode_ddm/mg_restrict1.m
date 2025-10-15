function rhalf = mg_restrict1(N,r)

%---
% half the grid
%---

rhalf(1)=2*r(1)+2*r(2);

Ic=1;
for i=3:2:N-1
 Ic=Ic+1;
 rhalf(Ic)=r(i-1)+2*r(i)+r(i+1);
end

return
