close all
clear all

%=====================
% evolve a closed line (loop)
% with point redistribution
% in the xy plane
%
% interpolate properties P1 P2
%
% P1 is arc length like
% P2 is periodic
%=====================

N = 16;

%---
% initial shape
%---

Dtheta = 2*pi/N;

for i=1:N+1
 theta = (i-1)*Dtheta;
 x(i) = cos(theta);
 y(i) = sin(theta);
 P1(i) = theta;
 P2(i) = 0.5+cos(theta+pi/3)^2;
end

%---
% parameters
%---

Ich1 = 1; thmax = pi/4;
Ich2 = 1; spmax = 0.20;
Ich3 = 1; spmin = 0.01;

Isym = 0;
Dt = 0.025;
nstep = 38;
Italk = 0;

%---
% evolve the line
%---

for istep=1:nstep

  istep
  Irepeat = 1;

  while(Irepeat==1)

         [N,x,y ...
         ,vnx,vny ...
         ,crv,s,Xint ...
         ,Axint,Bxint,Cxint ...
         ,Ayint,Byint,Cyint ...
         ,area,centerx,centery ...
         ,P1,P2 ...
         ,Istop ...
         ,Irepeat ...
         ,Iaction] ...
...
   = prd_2d (N ...
            ,x,y ...
            ,P1,P2 ...
            ,Ich1,thmax ...
            ,Ich2,spmax ...
            ,Ich3,spmin ...
            ,Isym ...
            ,Italk);

     if(Istop==1)
      disp(' prd_2d_dr: something went wrong');
      disp('            code:');
      disp(Iaction)
      break
    end

  end

  if(Istop==1)
    break
  end

  figure(1)
  clf
  hold on
  axis([-5 5 -5 5 0 2])
  axis square
  plot(x,y,'k.-')
  patch(x,y,'c')
  plot3(x,y,0.2*P1,'b-')
  plot3(x,y,P2,'r-')
  box on
  view(6,32)
  pause(0.01)

   for i=1:N+1
    uxvel = x(i)+y(i)^3;
    uyvel =-y(i)+x(i)^2;
    x(i) = x(i)+ uxvel*Dt;
    y(i) = y(i)+ uyvel*Dt;
  end

%---
end  % of time
%---

