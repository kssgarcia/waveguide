function [xhalf, yhalf, rhalf] = mg_restrict(Nx,Ny,x,y,r)

%--------------------------------
% multigrid method: half the grid
%
% position (x,y), residual(r)
%----------------------------

%---
% grid lines
%---

Ic=0;
for i=1:2:Nx+1
 Ic=Ic+1;
 xhalf(Ic)=x(i);
end

Ic=0;
for j=1:2:Ny+1
 Ic=Ic+1;
 yhalf(Ic)=y(j);
end

%---
% full weighting
%---

rhalf=zeros(Ny/2+1,Ny/2+1);

Jc=1;
for j=3:2:Ny-1
 Jc=Jc+1;
 Ic=1;
 for i=3:2:Nx-1
  Ic=Ic+1;
  rhalf(Ic,Jc)= 0.25*(r(i-1,j+1)+2*r(i,j+1)+r(i+1,j+1) ...
              +2*r(i-1,j)+4*r(i,j)+2*r(i+1,j) ...
              +r(i-1,j-1)+2*r(i,j-1)+r(i+1,j-1));
 end
end

%---
% done
%---

return
