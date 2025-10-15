close all
clear all

%=====================
% evolve an open line with point redistribution
% interpolate properties P1 P2
%
% P1 is arc length like
% P2 is periodic
%=====================

N = 16;

%---
% initial shape
%---

Dx = 1.0/N;

for i=1:N+1
 xflat = (i-1)*Dx;
 x(i) = xflat;
 y(i) = 0.1*sin(pi*xflat);
 P1(i) = xflat;
 P2(i) = xflat^4;
end

%---
% parameters
%---

Ich1 = 1; thmax = 0.25*pi;
Ich2 = 1; spmax = 0.075;
Ich3 = 1; spmin = 0.010;

Isym = 0;
Dt = 0.025;
nstep = 100;
Italk = 0;
Italk = 1;

%---
% evolve the line
%---

for istep=1:nstep

  istep
  Irepeat = 1;

  while(Irepeat==1)

         [N,x,y ...
         ,vnx,vny ...
         ,crv ...
         ,s ...
         ,Xint ...
         ,Axint,Bxint,Cxint ...
         ,Ayint,Byint,Cyint ...
         ,area,centerx,centery ...
         ,P1,P2 ...
         ,Istop ...
         ,Irepeat ...
         ,Iaction] ...
...
   = prd_2d_open (N ...
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
  axis([0 2 -0.5 0.5 0 1])
  axis square
  plot(x,y,'k.-')
  plot3(x(1:N+1),y(1:N+1),P1(1:N+1),'b-')
  plot3(x(1:N+1),y(1:N+1),P2(1:N+1),'r-')
  box on
  view(6,32)
  pause(0.01)

  Ido = 0;
  Ido = 1;
  if(Ido==1)
   figure(2)
   clf
   hold on
%  axis([0 2 -0.5 0.5 0 1])
   axis square
   plot3(x(1:N+1),y(1:N+1),crv(1:N+1),'k-o')
   box on
   view(6,32)
   pause(0.01)
  end

   for i=1:N+1
    uxvel = x(i)+y(i)^3;
    uyvel =-y(i)+x(i)^2;
    uxvel = 10*y(i);
    uyvel = x(i)*(1-x(i));
    x(i) = x(i)+ uxvel*Dt;
    y(i) = y(i)+ uyvel*Dt;
  end

%---
end  % of time
%---
