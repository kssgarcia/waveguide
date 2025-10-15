function [xi,eta,x,y,h] = conf2_grid ...
...
      (alpha,L,Nxil,Nxic,Nxir,Neta)

%==============================================
% Orthogonal grid in the exterior of a periodic
% array of cylinders using conformal mapping
%==============================================

%---
% prepare
%---

k=2*pi/L;
xs=acos(1.0-2.0*alpha)/k;
xis=xs + sin(k*xs)/k;

%---
% xi divisions
%
% left, central, and right
%---

Ic=0;
step=(L/2-xis)/Nxil;
for i=1:Nxil+1
  Ic=Ic+1;
  xi(Ic)=-L/2+(i-1.0)*step;
end

step=2*xis/Nxic;
for i=2:Nxic
  Ic=Ic+1;
  xi(Ic)=-xis+(i-1.0)*step;
end

step=(L/2-xis)/Nxir;
for i=1:Nxir+1
  Ic=Ic+1;
  xi(Ic)=xis+(i-1.0)*step;
end
Nxi=Ic-1;

Nxi=Nxil+Nxic+Nxir;

%---
% eta divisions
%---

eta_top=0.5*L;
eta_top=2.0*L;
eta_top=L;
eta_top=1.5*L;

eta=linspace(0,eta_top,Neta+1);

%===============
% Cartesian grid
%===============

Iskip=1;
if(Iskip==0)

hold on

for j=1:Neta+1
 for i=1:Nxi
  plot([xi(i),xi(i+1)],[eta(j),eta(j)],'-r')
 end
end
for i=1:Nxi+1
 for j=1:Neta
  plot([xi(i),xi(i)],[eta(j),eta(j+1)],'-r')
 end
end

axis([-L/2, L/2, 0, L])
set(gca,'fontsize',15)
xlabel('\xi','fontsize',15)
ylabel('\eta','fontsize',15)
axis('equal')
box

end

%============
% mapped grid
%============

%---
% eta=0
%---

 for i=1:Nxil
   [x(i,1), y(i,1)] = conf2_root(alpha,k,xi(i),0.0,xi(i),0.0);
 end
 x(Nxil+1,1)=-xs;
 y(Nxil+1,1)=0.0;

 for i=Nxil+2:Nxil+Nxic
%   [x(i,1), y(i,1)] = conf2_root(alpha,k,xi(i),0.00001,xi(i),0.01);
   [x(i,1), y(i,1)] = conf2_root(alpha,k,xi(i),0.000001,x(i-1),y(i-1));
 end

 x(Nxil+Nxic+1,1) = xs;
 y(Nxil+Nxic+1,1) = 0.0;
 for i=Nxil+Nxic+2:Nxi+1
   [x(i,1), y(i,1)] = conf2_root(alpha,k,xi(i),0.0,xi(i),0.0);
 end

%---
% eta>0
%---

 for j=2:Neta+1
  i=1;
  [x(i,j), y(i,j)] = conf2_root(alpha,k,xi(i),eta(j),x(i,j-1),y(i,j-1));
  for i=2:Nxi+1
%  [x(i,j), y(i,j)] = conf2_root(alpha,k,xi(i),eta(j),xi(i),eta(j));
   [x(i,j), y(i,j)] = conf2_root(alpha,k,xi(i),eta(j),x(i-1,j),y(i-1,j));
  end
 end

%---
% plot
%---

Iskip=1;
if(Iskip==0)

figure
hold on

for j=1:Neta+1
 for i=1:Nxi
  plot([x(i,j),x(i+1,j)],[y(i,j),y(i+1,j)],'-')
 end
end
for i=1:Nxi+1
 for j=1:Neta
  plot([x(i,j),x(i,j+1)],[y(i,j),y(i,j+1)],'-')
 end
end

axis([-L/2, L/2, 0, L])
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
axis equal
box

end

%==============
% some solution
%==============

Iskip=1;
if(Iskip==0)

figure
hold on

for i=1:Nxi
 for j=1:Neta
 A=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)];
 B=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1)];
 C=[eta(j),eta(j),eta(j+1),eta(j+1)];
 patch(A,B,C,C);
 end
end
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('f','fontsize',15)
axis equal
box

end

%===================
% metric coefficient
%===================

cf = 2.0*alpha;

for i=1:Nxi+1
 for j=1:Neta+1
  kx=k*x(i,j);
  ky=k*y(i,j);
  den = (cosh(ky)-cos(kx))^2;
  dxidx = 1.0 + cf*(cos(kx)*cosh(ky)-1.0)/den;
  dxidy = -cf*sin(kx)*sinh(ky)/den;
  h(i,j) = 1.0/sqrt(dxidx^2+dxidy^2);
 end
end

%---
% plot the metric coefficient
%---

Iskip=1;
if(Iskip==0)

figure
hold on

for i=1:Nxi
 for j=1:Neta
 A=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)];
 B=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1)];
 C=[h(i,j),h(i+1,j),h(i+1,j+1),h(i,j+1)];
 D=[0.1, 0.1, 0.1, 0.1];
 patch(A,B,C,D);
 end
end

set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('h','fontsize',15)
axis([-0.5 0.5 0.0 1.0 0.5 1.5])
axis equal
box

end

%-----
% Done
%-----

return
