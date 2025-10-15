close all;
clear all;

%===============================================
% Motion of two intercepting spherical particles
% in simple shear flow
%
% u_x = y, u_y=0, u_z=0
%
% SYMBOLS
%
% x,y,z: particle center position
% xsep,ysep,zsep: particle center separation
%===============================================

%----
% integrate particle center
%----

x(1) =-10.0; % for the shaded area:
y(1) = 0.01;
z(1) = 0.0;
x(1) =-10.0; % arbitrary
y(1) = 0.60;
z(1) = 0.0;
x(1) =-5.0; % for the shaded area:
y(1) = 0.0;
z(1) = 0.0;

% time step:

%Dt = 0.0025;
Dt = 0.0075;

%-----
% prepare
%----

xsep(1) = 2.0*x(1);
ysep(1) = 2.0*y(1);
zsep(1) = 2.0*z(1);

Dth = 0.5*Dt;

%-----
% launch
%----

nstep = 80000;

for i=1:nstep

 sep = [xsep(i),ysep(i),zsep(i)];
 f = two_spheres_vel(sep);

 xsepsave(i) = xsep(i);
 ysepsave(i) = ysep(i);
 zsepsave(i) = zsep(i);
 fsave = f;
 xseptmp = xsep(i)+Dt*f(1);
 yseptmp = ysep(i)+Dt*f(2);
 zseptmp = zsep(i)+Dt*f(3);

 sep = [xseptmp,yseptmp,zseptmp];

 f = two_spheres_vel(sep);

 xsep(i+1) = xsep(i) + Dth*(f(1)+fsave(1));
 ysep(i+1) = ysep(i) + Dth*(f(2)+fsave(2));
 zsep(i+1) = zsep(i) + Dth*(f(3)+fsave(3));

% if(xsep(i)>0) break; end
  if(xsep(i)>0) break; end

end

%------------
% end of time integration
%------------

%---
% number of steps
%---

idle = size(xsep);
m = idle(2);

%---
% particle position
%---

for i=1:m
 x(i) = 0.5*xsep(i);
 y(i) = 0.5*ysep(i);
 z(i) = 0.5*zsep(i);
 gap(i) = sqrt( x(i)^2+y(i)^2+z(i)^2 )-1.0;
end

%------------
% graphics
%------------

Ido = 0;
Ido = 1;

if(Ido==1)

 figure(44)
 hold on
 axis([-2 2 -2 2])
 axis equal
 %axis([-1.5 1.5 0.6 1.2])
 xlabel('x/a','fontsize',15)
 ylabel('y/a','fontsize',15)
 title('sphere center trajectories in the xy plane')
 set(gca,'fontsize',15)
 box on

% shade the area

 Ic=0;
 for i=1:m
  Ic=Ic+1;
  xp(Ic) = x(i);
  yp(Ic) = y(i);
 end
 for i=1:m-1
  Ic=Ic+1;
  xp(Ic) =-x(m-i);
  yp(Ic) = y(m-i);
 end
 fill(xp,yp,'y');

  Ic=0;
  for i=1:m
   Ic=Ic+1;
   xp(Ic) =-x(i);
   yp(Ic) =-y(i);
  end
  for i=1:m-1
   Ic=Ic+1;
   xp(Ic) = x(m-i);
   yp(Ic) =-y(m-i);
  end
  fill(xp,yp,'y');

  %plot(xx, yy,'.')
  plot( x, y)
  plot(-x, y)
  plot( x,-y)
  plot(-x,-y)

 % draw a circle

 ncrcl = 64;
 step = 2*pi/ncrcl;
 for i=1:ncrcl
  th = step*(i-1.0);
  xcrcl(i) = cos(th);
  ycrcl(i) = sin(th);
 end
 fill(xcrcl,ycrcl,'--w');
% plot(xcrcl,ycrcl,'--k');

%----
end
%----

%---
% plot the gap
%---

figure(45)
hold on
gap = log10(gap);
xlabel('x/a','fontsize',15)
ylabel('log_{10}(\epsilon/a)','fontsize',15)
set(gca,'fontsize',15)
box on

plot( x,gap,'k');
plot(-x,gap,'k');
