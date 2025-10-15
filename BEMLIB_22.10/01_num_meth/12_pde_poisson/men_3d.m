%-------
clear all
close all
%------

%=======================================
% meniscus in the exterior of an ellipse
%=======================================

b=0.5; % ellipse semi-axis
c=1.0; % ellipse semi-axis

capl = 0.75;   % capillary length
capl = 0.075;   % capillary length
capl = 0.50;   % capillary length

xcline=1.0; % height of the contact line

tol=0.000001;  % tolerance
Niter=1000;  % max number of iterations

%----------
% divisions
%----------

Nu = 16;
Nphi = 16;

Nu=32;
Nphi=32;

%---
% prepare
%---

Dphi = 2*pi/Nphi;
capls = capl*capl;
Dphi2 = 2.0*Dphi;

u0 = atanh(b/c);
snhu0 = sinh(u0);
cshu0 = cosh(u0);
A = b/snhu0;

umax = log(32.0*b/A);

Du = (umax-u0)/Nu;
Du2 = 2.0*Du;

%---
% grid
%---

for i=1:Nu+1

 u(i) = u0+(i-1.0)*Du;
 snhu(i) = sinh(u(i));
 cshu(i) = cosh(u(i));

 for j=1:Nphi+1

  phi(j)=(j-1.0)*Dphi;
  csphi(j)=cos(phi(j));
  snphi(j)=sin(phi(j));
  y(i,j) = A*snhu(i)*csphi(j);
  z(i,j) = A*cshu(i)*snphi(j);
%  h(i,j) = A*sqrt(cshu^2-snphi^2)/cosh(u);
  h(i,j) = A*sqrt(cshu(i)^2-snphi(j)^2);

 end

end

%---
% plot the grid
%---

Ido = 0;

if(Ido==1)

hold on;
for i=1:Nu+1
 plot(y(i,:),z(i,:),'r')
end
for j=1:Nphi
 plot(y(:,j),z(:,j),'r')
end

%for i=1:Nu+1
% plot3(y(i,:),z(i,:),h(i,:))
%end
%for j=1:Nphi
% plot3(y(:,j),z(:,j),h(:,j))
%end
%for j=1:Nphi
% plot3(y(:,j),z(:,j),h(:,j))
%end

xlabel('y','fontsize',15);
ylabel('z','fontsize',15);
zlabel('x','fontsize',15);
set(gca,'fontsize',15)

axis square
axis([-2 2 -2 2 -2 2])
view(140,30)

figure

end

%---
% initialize
%---

for i=1:Nu+1
 for j=1:Nphi+2
  x(i,j) = 0.0;
 end
end

%---
% boundary conditions
%---

for j=1:Nphi+1
  x(1,j) = xcline;  % constant height
  x(1,j) = xcline*y(1,j);  % rotated
  x(Nu+1,j) = 0.0;
end

x(1,   Nphi+2)=x(1,2);
x(Nu+1,Nphi+2) = 0.0;

%---
% iterations
%---

for iterations=1:Niter

%---
% compute the first derivatives
%---

for i=2:Nu

 for j=2:Nphi+1
  MAT(1,1) =  A*cshu(i)*csphi(j);
  MAT(1,2) =  A*snhu(i)*snphi(j);
  MAT(2,1) = -MAT(1,2);
  MAT(2,2) =  MAT(1,1);
  RHS(1)=(x(i+1,j)-x(i-1,j))/Du2;
  RHS(2)=(x(i,j+1)-x(i,j-1))/Dphi2;
  SOL=RHS/MAT';
  dxdy(i,j)=SOL(1);
  dxdz(i,j)=SOL(2);
 end

 dxdy(i,1)=dxdy(i,Nphi+1);
 dxdz(i,1)=dxdz(i,Nphi+1);
 dxdy(i,Nphi+2)=dxdy(i,2);
 dxdz(i,Nphi+2)=dxdz(i,2);

end

for j=1:Nphi+2
  dxdy(1,j)   =2.0*dxdy(2,j) -dxdy(3,j);
  dxdy(Nu+1,j)=2.0*dxdy(Nu,j)-dxdy(Nu-1,j);
  dxdz(1,j)   =2.0*dxdz(2,j) -dxdz(3,j);
  dxdz(Nu+1,j)=2.0*dxdz(Nu,j)-dxdz(Nu-1,j);
end

%---
% compute the second derivatives
%---

for i=2:Nu
 for j=2:Nphi+1
  MAT(1,1) =  A*cshu(i)*csphi(j);
  MAT(1,2) =  A*snhu(i)*snphi(j);
  MAT(2,1) = -MAT(1,2);
  MAT(2,2) =  MAT(1,1);
  RHS(1) = (dxdy(i+1,j)-dxdy(i-1,j))/Du2;
  RHS(2) = (dxdy(i,j+1)-dxdy(i,j-1))/Dphi2;
  SOL = RHS/MAT';
  dxdyy(i,j) = SOL(1);
  dxdyz(i,j) = SOL(2);
  RHS(1)=(dxdz(i+1,j)-dxdz(i-1,j))/Du2;
  RHS(2)=(dxdz(i,j+1)-dxdz(i,j-1))/Dphi2;
  SOL=RHS/MAT';
  dxdzy(i,j)=SOL(1);
  dxdzz(i,j)=SOL(2);
 end

end

%------
% scan the grid points
%------

errr=0.0;

for i=2:Nu
 for j=2:Nphi+1
  H=h(i,j);
  HS=H*H;
  tmp = 1.0+dxdy(i,j)^2+dxdz(i,j)^2;
  G = 2.0*(1.0/Du^2+1.0/Dphi^2)/HS + tmp^(3/2)/capls;
  xnew = (x(i+1,j)+x(i-1,j))/(HS*Du^2) ...
        +(x(i,j+1)+x(i,j-1))/(HS*Dphi^2);
  xnew = xnew + dxdz(i,j)^2 * dxdyy(i,j);
  xnew = xnew - 2.0*dxdy(i,j)*dxdz(i,j)*dxdyz(i,j);
  xnew = xnew + dxdy(i,j)^2 * dxdzz(i,j);
  xnew=xnew/G;
  corr = abs(xnew-x(i,j));
  x(i,j) = xnew;
  if(corr>errr)
    errr = corr;
  end
 end
end

 for i=2:Nu+1
  x(i,1)     =x(i,Nphi+1);
  x(i,Nphi+2)=x(i,2);
 end

 if(errr<tol)
  break
 end

%--
end
%--

errr

if(errr>tol)
 disp('the iterations did not converge')
 errr
 return
end

%---
% plotting
%---

hold on;

%for i=1:Nu+1
% plot3(y(i,:),z(i,:),x(i,1:Nphi+1))
%end
%for j=1:Nphi
% plot3(y(:,j),z(:,j),x(1:Nu+1,j))
%end

for i=1:Nu
 for j=1:Nphi
  patch([y(i,j), y(i,j+1), y(i+1,j+1),  y(i+1,j)], ...
        [z(i,j), z(i,j+1), z(i+1,j+1),  z(i+1,j)], ...
        [x(i,j), x(i,j+1), x(i+1,j+1),  x(i+1,j)], ...
        [x(i,j), x(i,j+1), x(i+1,j+1),  x(i+1,j)]);
 end
end

xlabel('y','fontsize',15);
ylabel('z','fontsize',15);
zlabel('x','fontsize',15);
set(gca,'fontsize',15)

axis square
axis([-1 1 -1 1 -1 1])
view(162, 12)

%---
% compute dxdu, dxdphi,
% the surface normal and tension vector
% around the contact line
%---

x0=0.0;
y0=0.0;
z0=0.0;

for j=1:Nphi+1

 dxdu(j) = (-x(3,j)+4.0*x(2,j)-3.0*x(1,j))/Du2;
 dydu(j) = A*cshu0*csphi(j);
 dzdu(j) = A*snhu0*snphi(j);
 dmdu(j) = sqrt(dxdu(j)^2+dydu(j)^2+dzdu(j)^2);

 if(j==1)
  dxdphi(j) = (x(1,2)-x(1,Nphi))/Dphi2;
 else
  dxdphi(j) = (x(1,j+1)-x(1,j-1))/Dphi2;
 end
 dydphi(j) =-A*snhu0*snphi(j);
 dzdphi(j) = A*cshu0*csphi(j);
 dmdphi(j) = sqrt(dxdphi(j)^2+dydphi(j)^2+dzdphi(j)^2);

 vnx(j) = dydu(j)*dzdphi(j)-dzdu(j)*dydphi(j);
 vny(j) = dzdu(j)*dxdphi(j)-dxdu(j)*dzdphi(j);
 vnz(j) = dxdu(j)*dydphi(j)-dydu(j)*dxdphi(j);
 vnm(j) = sqrt(vnx(j)^2+ vny(j)^2+vnz(j)^2);
 vnx(j) = vnx(j)/vnm(j);
 vny(j) = vny(j)/vnm(j);
 vnz(j) = vnz(j)/vnm(j);

 tngx(j) = dydphi(j)*vnz(j)-dzdphi(j)*vny(j);
 tngy(j) = dzdphi(j)*vnx(j)-dxdphi(j)*vnz(j);
 tngz(j) = dxdphi(j)*vny(j)-dydphi(j)*vnx(j);

 crsx(j) = (y(1,j)-y0)*tngz(j)-(z(1,j)-z0)*tngy(j);
 crsy(j) = (z(1,j)-z0)*tngx(j)-(x(1,j)-x0)*tngz(j);
 crsz(j) = (x(1,j)-x0)*tngy(j)-(y(1,j)-y0)*tngx(j);

end

%---
% force and torque
%---

forcex=0.0;
forcey=0.0;
forcez=0.0;
torqux=0.0;
torquy=0.0;
torquz=0.0;

for j=1:Nphi
  forcex = forcex+tngx(j);
  forcey = forcey+tngy(j);
  forcez = forcez+tngz(j);
  torqux = torqux+crsx(j);
  torquy = torquy+crsy(j);
  torquz = torquz+crsz(j);
end

forcex = forcex*Dphi;
forcey = forcey*Dphi;
forcez = forcez*Dphi;
torqux = torqux*Dphi;
torquy = torquy*Dphi;
torquz = torquz*Dphi;

[forcex, forcey, forcez]
[torqux, torquy, torquz]

