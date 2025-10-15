close all
clear all
hold on

%=======
% stability graphs of the convection equation
% in two dimensions
%=====

Nx=2*2*2*2*8;
Ny=2*2*2*2*8;
Nth=2*2*8;

ax=-1;bx=1;
ay=-1;by=1;

xx = linspace(ax, bx, Nx+1);
yy = linspace(ay, by, Ny+1);

pih=0.5*pi;
im = sqrt(-1);

opt=1;

%---
for i=1:Nx+1
 x=xx(i);
 for j=1:Ny+1
  y=yy(j);
%---

  Gmax = 0.0;
  for lx=1:Nth+1
   thetax=(lx-1)*pih/Nth;
   for ly=1:Nth+1
    thetay=(ly-1)*pih/Nth;

%---
    if(opt==1) % LAX
    G = 0.5*(cos(thetax)+cos(thetay)) + im*(x*sin(thetax)+y*sin(thetay));
%    G = cos(thetax) + im*x*sin(thetax);
    end
%---
    Gtest = abs(G);
    if(Gtest>Gmax)
      Gmax=Gtest;
     end
   end
  end

  if(Gmax<1.00001)
   plot(x,y,'k+');
  end

%---
 end
end
%---

rad=1/sqrt(2);
Nc=64
for icircle=1:Nc+1
 thc = 2*pi*(icircle-1)/Nc;
  xcr(icircle)=rad*cos(thc);
  ycr(icircle)=rad*sin(thc);
end
plot(xcr,ycr,'r')

axis([ax bx ay by])
box
set(gca,'fontsize',15)
xlabel('c_x','fontsize',15)
ylabel('c_y','fontsize',15)
axis square

