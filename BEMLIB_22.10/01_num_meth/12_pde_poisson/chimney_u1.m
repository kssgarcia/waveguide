clear all;
close all;

ho=10;		% W/(m^2*K), outer heat xfer coeff.
hi=25;		% W/(m^2*K), inner heat xfer coeff.
TA=10;	% K, ambient temperature
TS=300;	% K, smoke temperature
k=0.72;		% W/(m*K), thermal conductivity of bricks
rho=1920;	% kg/m^3, brick density
cp=835;		% J/(kg*K), brick heat capacity at constant pressure
a=0.4;		% width of chimney wall
N=32;		% number of divisions
Nstep=5000;   % number of steps

dx=a/N; 
dy=dx; 
alpha=0.2;		% diffusion number (k*dt)/(rho*cp*dx^2)
kappa=k/(rho*cp);	% diffusivity
dt=alpha*dx^2/kappa;	% time step

%---
% grid lines
%---

for i=1:2*N+1
  X(i)=dx*(i-1);
end
for j=1:N+1
  Y(j)=dy*(j-1);
end

%---
% fictitious lines
%---

for i=1:2*N+1
  Tj0(i)=TA;
end
for j=1:N+1
  Ti0(j)=TA;
end

%---
% initialize
%---

for i=1:2*N+2
 for j=1:N+2
   T(i,j)=TA;
 end
end

Tnew=T;

%==============
% time stepping
%==============

for step=1:Nstep

%---
% FTCS
%---

 for i=1:2*N+1
   for j=2:N+1
    if(i==1)
     Tnew(i,j)=T(i,j)+alpha*(Ti0(j)-4*T(i,j)+T(i+1,j)+T(i,j-1)+T(i,j+1));
    else
     Tnew(i,j)=T(i,j)+alpha*(T(i-1,j)-4*T(i,j)+T(i+1,j)+T(i,j-1)+T(i,j+1));
    end
   end
 end

 for i=N+1:2*N+1
  Tnew(i,1)=T(i,1)+alpha*(T(i-1,1)-4*T(i,1)+T(i+1,1)+Tj0(i)+T(i,j+1));
 end

 T=Tnew;  % update

%---
% Boundary and symmetry conditions
%---

for j=1:N+1 % left
    Ti0(j)=T(2,j)-(2*ho*dx/k)*(T(1,j)-TA);
end
for j=1:N+1 % right
    T(2*N+2,j)=T(2*N,j); 
end
for i=1:2*N+1 % top
     T(i,N+2)=T(i,N)-(2*ho*dy/k)*(T(i,N+1)-TA); % upper
end
for i=1:N+1   % first bottom
  T(i,1)=T(N+1,N+2-i);
end
for i=N+1:2*N+2 % second bottom
    Tj0(i)=T(i,2)-(2*hi*dy/k)*(T(i,1)-TS);
end

 time(step)=step*dt/3600;  % time vector for plotting
 temp(step)=T(1,1);

 Z = T(1:2*N+1,1:N+1)';

if(step==1)
  Handle1 = surf(X,Y,Z);
  set(Handle1, 'erasemode', 'xor');
  set(gca,'fontsize',15)
  axis([-0.0*a 2*a -0.0*a 2.0*a TA TS])
  xlabel('x','fontsize',15)
  ylabel('y','fontsize',15)
  zlabel('T(C)','fontsize',15)
  view(-140,40)
else
  set(Handle1,'XData',X,'YData',Y,'ZDATA',Z);
  drawnow
end

end

%---

figure(2)
plot(time,temp)
xlabel('Time (hours)')
ylabel('T(1,1) in C')
