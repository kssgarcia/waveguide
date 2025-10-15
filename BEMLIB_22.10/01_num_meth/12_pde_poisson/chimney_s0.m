clear all
close all

%===============
% chimney_s
%
% Steady-state temperature distribution
% over the cross-section
% of a square chimney wall
% with Dirichlet BC
%===============

%--------------------------
% prameters and conditions:
%--------------------------

T1=300;         % temperature of inside wall in deg C
T2=10;          % temperature of outside wall in deg C
a=0.5;          % thickness of chimney wall in meters
k = 0.15;       % conductivity
tolerance = 0.00001;

N=32;   % number of nodes in x-direction
maxiter=3000;% maximum # of iterations

%--------
% prepare
%--------


M=N/2;  %  number of nodes in y-direction
b=2*a;

dx=b/N;
dy=a/M;

beta=(dx/dy)^2;
beta1=2.0*(1.0+beta);

%--------------
% initial guess
%--------------

for i=1:N+1
   for j=1:M+1
      T(i,j)=T2;
   end
end

%--------------------
% boundary conditions
%--------------------

for j=1:M+1
   T(1,j)=T2;          % Dirichlet on left edge (outer wall),
   T(N+2,j)=T(N,j);     % and right edge (x=b; also using symmetry)
end

for i=1:M+1
   T(i,1)=T(M+1,M+2-i); % lower edge (from 0 to a; using symmetry)
end

for i=1:N+1     % Dirichlet b.c. and on upper edge (outer wall)
   T(i,M+1)=T2;
end

for i=M+1:N+1   % ...and also on lower edge (bordering inner
   T(i,1)=T1;   % wall; from a to b)
end

%-----------
% iterations
%-----------

for n=1:maxiter         % iterate until convergence

   correction = 0.0;

   for i=2:N+1  % central finite-diff discretization
      for j=2:M         % del^2(T)=0 combined w/ point-Gauss-Siedel
        Told = T(i,j);
        T(i,j)=(T(i+1,j)+T(i-1,j)+beta*(T(i,j+1)+T(i,j-1)))/beta1;
        diff = abs(T(i,j)-Told);
          if (diff>correction)
           correction = diff;
         end
      end
   end

   for i=1:M+1          % reset the bottom (from 0 to a)
      T(i,1)=T(N/2+1,M+2-i);
   end
   for j=1:M+1                  % reset the right edge
      T(N+2,j)=T(N,j);
   end

  correction
  if(correction<tolerance) break; end;

end

%------------------
% end of iterations
%------------------

for i=1:N+1             % set up plotting vectors
  X(i)=dx*(i-1);
end
for j=1:M+1
  Y(j)=dy*(j-1);
end

for i=1:N+1
   X1(i)=4*a-X(i);
end

%---------
% plotting
%---------

hold on

mesh(X,Y,T(1:N+1,:)')   % plotting, labelling, and formatting
mesh(X1,Y,T(1:N+1,:)')
mesh(Y+3.0*a,X-3*a,T(1:N+1,:))
mesh(Y+3.0*a,X1-3*a,T(1:N+1,:))
mesh(X,-2*a-Y,T(1:N+1,:)')
mesh(X1,-2*a-Y,T(1:N+1,:)')
mesh(a-Y,X-3*a,T(1:N+1,:))
mesh(a-Y,X1-3*a,T(1:N+1,:))

xlabel('x (m)','fontsize',15)
ylabel('y (m)','fontsize',15)
zlabel('T (C)','fontsize',15)
set(gca,'fontsize',15)

%title('Temperature over the cross-section of chimney wall')
%axis([0 2.0 0 0.0 0 2.0])
view(-80,40)

%------------------------------------
% compute the flux on the outer wall
%------------------------------------

%---
% top segment:
%---

for i=1:N+1
 flux_top(i) = -k*(T(i,M+1)-T(i,M))/dy;
end

%------------------------------------
% build the flux around the outer wall
%------------------------------------

Ic = 1;  % counter

fluxo(1)=flux_top(1);
perimo(1)=0;
xo(1)=0;
yo(1)=a;

for i=2:N+1
 Ic=Ic+1;
 fluxo(Ic)=flux_top(i);
 perimo(Ic)=perimo(Ic-1)+dx;
 xo(Ic)=xo(Ic-1)+dx;
 yo(Ic)=yo(Ic-1);
end
for i=2:N+1
 Ic=Ic+1;
 fluxo(Ic)=flux_top(N+2-i);
 perimo(Ic)=perimo(Ic-1)+dx;
 xo(Ic)=xo(Ic-1)+dx;
 yo(Ic)=yo(Ic-1);
end
for j=2:N+1
 Ic=Ic+1;
 fluxo(Ic)=flux_top(j);
 perimo(Ic)=perimo(Ic-1)+dy;
 xo(Ic)=xo(Ic-1);
 yo(Ic)=yo(Ic-1)-dy;
end
for j=2:N+1
 Ic=Ic+1;
 fluxo(Ic)=flux_top(N+2-j);
 perimo(Ic)=perimo(Ic-1)+dy;
 xo(Ic)=xo(Ic-1);
 yo(Ic)=yo(Ic-1)-dy;
end
for i=2:N+1
 Ic=Ic+1;
 fluxo(Ic)=flux_top(i);
 perimo(Ic)=perimo(Ic-1)+dx;
 xo(Ic)=xo(Ic-1)-dx;
 yo(Ic)=yo(Ic-1);
end
for i=2:N+1
 Ic=Ic+1;
 fluxo(Ic)=flux_top(N+2-i);
 perimo(Ic)=perimo(Ic-1)+dx;
 xo(Ic)=xo(Ic-1)-dx;
 yo(Ic)=yo(Ic-1);
end
for j=2:N+1
 Ic=Ic+1;
 fluxo(Ic)=flux_top(j);
 perimo(Ic)=perimo(Ic-1)+dy;
 xo(Ic)=xo(Ic-1);
 yo(Ic)=yo(Ic-1)+dy;
end
for j=2:N+1
 Ic=Ic+1;
 fluxo(Ic)=flux_top(N+2-j);
 perimo(Ic)=perimo(Ic-1)+dy;
 xo(Ic)=xo(Ic-1);
 yo(Ic)=yo(Ic-1)+dy;
end

%-----------------------------------
% compute the flux on the inner wall
%-----------------------------------

%------------
% top segment 
%------------

for i=1:M+1
 flux_top(i) = -k*(T(M+i,2)-T(M+i,1))/dy;
end

%------------------------------------
% build the flux around the inner wall
%------------------------------------

Ic = 1;  % counter
fluxi(1)=flux_top(1);
perimi(1)=0;
xi(1)=a;
yi(1)=0;

for i=2:M+1
 Ic=Ic+1;
 perimi(Ic)=perimi(Ic-1)+dx;
 fluxi(Ic)=flux_top(i);
 xi(Ic)=xi(Ic-1)+dx;
 yi(Ic)=yi(Ic-1);
end
for i=2:M+1
 Ic=Ic+1;
 perimi(Ic)=perimi(Ic-1)+dx;
 fluxi(Ic)=flux_top(M+2-i);
 xi(Ic)=xi(Ic-1)+dx;
 yi(Ic)=yi(Ic-1);
end
for j=2:M+1
 Ic=Ic+1;
 perimi(Ic)=perimi(Ic-1)+dy;
 fluxi(Ic)=flux_top(j);
 xi(Ic)=xi(Ic-1);
 yi(Ic)=yi(Ic-1)-dy;
end
for j=2:M+1
 Ic=Ic+1;
 perimi(Ic)=perimi(Ic-1)+dy;
 fluxi(Ic)=flux_top(M+2-j);
 xi(Ic)=xi(Ic-1);
 yi(Ic)=yi(Ic-1)-dy;
end
for i=2:M+1
 Ic=Ic+1;
 fluxi(Ic)=flux_top(i);
 perimi(Ic)=perimi(Ic-1)+dx;
 xi(Ic)=xi(Ic-1)-dx;
 yi(Ic)=yi(Ic-1);
end
for i=2:M+1
 Ic=Ic+1;
 perimi(Ic)=perimi(Ic-1)+dx;
 fluxi(Ic)=flux_top(M+2-i);
 xi(Ic)=xi(Ic-1)-dx;
 yi(Ic)=yi(Ic-1);
end
for j=2:M+1
 Ic=Ic+1;
 fluxi(Ic)=flux_top(j);
 perimi(Ic)=perimi(Ic-1)+dy;
 xi(Ic)=xi(Ic-1);
 yi(Ic)=yi(Ic-1)+dy;
end
for j=2:M+1
 Ic=Ic+1;
 fluxi(Ic)=flux_top(M+2-j);
 perimi(Ic)=perimi(Ic-1)+dy;
 xi(Ic)=xi(Ic-1);
 yi(Ic)=yi(Ic-1)+dy;
end

figure
hold on
%title('Wall flux')
%plot(perimo,fluxo,perimi,fluxi,'--');
%legend('outer','inner')
%xlabel('Perimeter (m)','fontsize',15)
%ylabel('flux (Watt)','fontsize',15)
plot3(xo,yo,zeros(size(xo)))
plot3(xo,yo,fluxo)
plot3(xi,yi,zeros(size(xi)),'r')
plot3(xi,yi,fluxi,'r')
set(gca,'fontsize',15)
xlabel('x(m)','fontsize',15)
ylabel('y(m)','fontsize',15)
zlabel('flux (Watt)','fontsize',15)

