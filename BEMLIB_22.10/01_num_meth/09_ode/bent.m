function bent

%---
% solution of a boundary-value problem
% using pvp4c
%---

clear all
close all

global alpha L slope

%---
% parameters
%---

L = pi;
slope = 0.0;
%---

%---
% mesh
%---

Ndiv = 64;

dl = L/Ndiv;
for i=1:Ndiv+1
  l(i) = (i-1.0)*dl;
end

delta0 = -0.5*(pi/L)^2;
solinit = bvpinit(l, @initialize, delta0);
solinit0 = solinit;
options = bvpset('RelTol',1e-2);

%---
% plot
%---

figure(1)
axis equal
box
plot([-1.5 1.5],[0,0],'k')
hold on

%---
for repeat=1:16
%---

 if(repeat==1)  alpha = 1; end
 if(repeat==2)  alpha = 1.1; end
 if(repeat==3)  alpha = 1.2; end
 if(repeat==4)  alpha = 1.3; end
 if(repeat==5)  alpha = 1.4; end
 if(repeat==6)  alpha = 1.5; end
 if(repeat==7)  alpha = 1.5; end
 if(repeat==8)  solinit = solinit0;alpha = 0.90; end
 if(repeat==9)  alpha = 0.80; end
 if(repeat==10)  alpha = 0.70; end
 if(repeat==11) alpha = 0.60; end
 if(repeat==12) alpha = 0.50; end
 if(repeat==13) alpha = 0.40; end
 if(repeat==14) alpha = 0.30; end
 if(repeat==15) alpha = 0.20; end
 if(repeat==16) alpha = 0.10; end

 sol = bvp4c(@odefun,@bcfun,solinit,options);
 solinit =sol;

 arclength= sol.x(1,:);
 x = sol.y(1,:);
 y = sol.y(2,:);
 plot(x, y,'k')
 plot(x,-y,'k')
 plot(x,y,'r.')

end

return

%----
% initialize
%----

function v = initialize(l)

global L

  theta = pi*l/L;
  cs = cos(theta);
  sn = sin(theta);
  v = [ -L/pi *cs
         L/pi *sn
               sn
               cs
         -pi/L
         0.0 ];

return

%---
% odes
%---

function dydx = odefun(l,x,delta)

dydx = [x(3)
        x(4)
       -x(5)*x(4)
        x(5)*x(3)
        x(6)
       -0.5*x(5)^3-delta*x(5)
        ];

return

%---
% boundary conditions
%---

function res = bcfun(xa,xb,delta)

global alpha L slope

cf1 = slope;
cf2 = sqrt(1.0-slope^2);

res = [ xa(1) + alpha*L/pi
        xa(2) 
        xa(3)-cf1
        xa(4)-cf2
        xb(1) - alpha*L/pi
        xb(2)
        xb(3)-cf1
       ];
%        xb(4)+cf2

return
