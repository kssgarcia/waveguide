%=======
% psi_t3
%
% Plot the element interpolation functions
% for a 3-node triangle
%=======

hold on

N=32;
M=32;

Dx=1.0/N;
Dy=1.0/M;

for i=1:N+1 
  x(i)=Dx*(i-1.0);
%   fprintf (1,'%6.2f',X(i))
end

for j=1:M+1
  y(j)=Dy*(j-1.0);
% fprintf (1,'%6.2f',Y(j))
end

%--------
% first set of grid lines
%--------

for i=1:N
 for j=1:M+2-i
  xp(j) = x(i);
  yp(j) = y(j);
  zeta= 1.0-xp(j)-yp(j);
  zp(j)= zeta;   % first node
 end
 plot3(xp,yp,zp);
end

%--------
% second set of grid lines
%--------

for j=1:M
 for i=1:N+2-j
  ypp(i) = y(j);
  xpp(i) = x(i);
  zeta= 1.0-xpp(i)-ypp(i);
  zpp(i)= zeta;   % first node
 end
plot3(xpp,ypp,zpp);
end

%------------------
% draw the triangle
%------------------

xx(1)=0.0;yy(1)=0.0;zz(1)=0.0;
xx(2)=1.0;yy(2)=0.0;zz(2)=0.0;
xx(3)=0.0;yy(3)=1.0;zz(3)=0.0;
xx(4)=0.0;yy(4)=0.0;zz(4)=0.0;

plot3(xx,yy,zz); hold on;
set(gca,'fontsize',15)
xlabel('\xi','fontsize',15)
ylabel('\eta','fontsize',15)

%-----
% done
%-----
