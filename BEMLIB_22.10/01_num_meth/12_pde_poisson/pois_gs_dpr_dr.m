close all
clear all

ax=0; bx=1.4;
ay=0; by=0.8;

Nx=64; Ny=64;
itermax=128*128;

tol=0.0000001;
relax=1.5;

fbot=0.0;
ftop=0.0;

Dx=(bx-ax)/Nx;
Dy=(by-ay)/Ny;

%---
% initialize
%---

for i=1:Nx+2
 for j=1:Ny+1
   g(i,j)=1.0;
   g(i,j)=sin(2*pi*(i-1)/Nx)*sin(2*pi*(j-1)/Ny);
 end
end

%---
% solve
%---

[f, iter, Iflag] = pois_gs_dpr ...
 ...
 (Nx,Ny,Dx,Dy,g,itermax,tol,relax,fbot,ftop);

%---
% plot
%---

for i=1:Nx+1
 x(i)=(i-1)*Dx;
end

for j=1:Ny+1
 y(j)=(j-1)*Dy;
end

mesh(x,y,f(1:Nx+1,1:Ny+1)');
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('f','fontsize',15)
axis([ax bx ay by])
