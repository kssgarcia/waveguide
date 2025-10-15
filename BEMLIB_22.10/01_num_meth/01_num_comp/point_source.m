%=========================
% Anmation of a diffusive
% and convective field
% due to a point source
%========================

clear all;
close all;

%--- parameters

U=1.0; k=0.0; D=0.1;
x0=0.0; y0=0.0; t0=0.0;
Nx=32; Ny=32; Dt=0.01;

%--- grid lines

for i=1:Nx+1
 x(i)=-2+(i-1.0)*10/Nx;
end

for j=1:Ny+1
 y(j)=-5+(j-1.0)*10/Ny;
end

%------------------------------
for step=1:1000  % loop over time
%------------------------------

t=0.01+(step-1.0)*Dt;

for i=1:Nx+1
 for j=1:Ny+1

 th=t-t0;
 xh=x(i)-x0-(U-k*y0)*th;
 yh=y(j)-y0;

 G(i,j) = 1.0/( (4*pi*D*th) * sqrt(1+k^2*th^2/12.0) ) ...
         *exp(-0.25*(xh-0.5*k*th*yh)^2/(1+k^2*th^2/12)/D/th ...
         -0.25*yh^2/D/th);

 end
end

if(step==1)
 Handle1 = mesh(G');
 axis([-2 8 -5 5 0 0.4])
 xlabel('x','fontsize',15)
 ylabel('y','fontsize',15)
 zlabel('G','fontsize',15)
 set(gca,'fontsize',15)
else
  set(Handle1,'XData',x,'YData',y,'ZData',G');
  drawnow
end

%---
end  % of time stepping
%--
