function [xdouble, ydouble, edouble] = mg_prolongate(Nx,Ny,x,y,e)

%---------------------------------------
% multigrid method: double the grid
%
% prolongate the position (x,y), error(e)
%---------------------------------------

%---
% interpolate the x grid lines
%---

Ic=1;
for i=1:Nx
 xdouble(Ic)=x(i);
 Ic=Ic+1;
 xdouble(Ic)=0.5*(x(i)+x(i+1));
 Ic=Ic+1;
end
xdouble(Ic)=x(Nx+1);

%---
% interpolate the y grid lines
%---

Ic=1;
for i=1:Ny
 ydouble(Ic)=y(i);
 Ic=Ic+1;
 ydouble(Ic)=0.5*(y(i)+y(i+1));
 Ic=Ic+1;
end
ydouble(Ic)=y(Nx+1);

%---
% interpolate the error
%---

edouble=zeros(2*Nx+1,2*Ny+1);

for j=2:Ny
 jfn=2*j-1;
 for i=2:Nx
 ifn=2*i-1;
 edouble(ifn-1,jfn+1)=edouble(ifn-1,jfn+1)+e(i,j)/4.0;
 edouble(ifn,jfn+1)=edouble(ifn,jfn+1)+e(i,j)/2.0;
 edouble(ifn+1,jfn+1)=edouble(ifn+1,jfn+1)+e(i,j)/4.0;
 edouble(ifn-1,jfn)=edouble(ifn-1,jfn)+e(i,j)/2.0;
 edouble(ifn,jfn)=edouble(ifn,jfn)+e(i,j);
 edouble(ifn+1,jfn)=edouble(ifn+1,jfn)+e(i,j)/2.0;
 edouble(ifn-1,jfn-1)=edouble(ifn-1,jfn-1)+e(i,j)/4.0;
 edouble(ifn,jfn-1)=edouble(ifn,jfn-1)+e(i,j)/2.0;
 edouble(ifn+1,jfn-1)=edouble(ifn+1,jfn-1)+e(i,j)/4.0;
 end
end

%---
% done
%---

return
