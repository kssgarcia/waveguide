close all
clear all
hold on

%=======
% stability graphs in 3D
%=====

im = sqrt(-1);

Nx=2*2*8;
Ny=2*2*8;
Nz=2*2*8;
Nth=2*8;

ax=0;bx=1;
ay=0;by=1;
az=0;bz=1;

xx = linspace(ax, bx, Nx+1);
yy = linspace(ay, by, Ny+1);
zz = linspace(az, bz, Nz+1);

opt=1;

%---
for i=1:Nx+1
 x=xx(i);
 for j=1:Ny+1
  y=yy(j);
   for k=1:Nz+1
    z=yy(k);
%---

  Gmax = 0.0;
  for lx=1:Nth+1
   thetax=(lx-1)*pi/Nth;
   for ly=1:Nth+1
    thetay=(ly-1)*pi/Nth;
     for lz=1:Nth+1
      thetaz=(lz-1)*pi/Nth;

      if(opt==1) % LAX
      G = (cos(thetax)+cos(thetay)+cos(thetaz))/3.0 ...
         + im*(x*sin(thetax)+y*sin(thetay)+z*sin(thetaz));
      end

    Gtest = abs(G);

    if(Gtest>Gmax)
      Gmax=Gtest;
    end
   end
  end
  end

  if(Gmax<1.0)
  Gmax
   plot3(x,y,z,'+');
  end

%---
  end
 end
end
%---
