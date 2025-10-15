% chimney_s1.m
%
% "Steady-state temperature distribution
%  over the cross-section
%  of a chimney wall."
%
%  With Neumann boundary conditions at the walls.
% 

ho=10;		% W/(m^2*K), outer heat xfer coeff.
hi=25;		% W/(m^2*K), inner heat xfer coeff.
To=10+273;	% K, ambient temperature
Ti=300+273;	% K, smoke temperature
k=0.72;		% W/(m*K), thermal conductivity of bricks
a=0.4;		% width of chimney wall
b=2*a;
c=a;
N=32;		% number of nodes in x-direction
M=N/2;		% ...and in y-direction
dx=b/N; dy=c/M;	% delta x and delta y step size
B=(dx/dy)^2;
maxiter=30;	% maximum # of iterations

for i=2:(N+2) 		% initialize all T values to To
   for j=2:(M+2)
      T(i,j)=To;
   end
end

for j=2:(M+2)
   T(1,j)=T(3,j)-(2*ho*dx/k)*(T(2,j)-To); % Neumann b.c. on left edge
					%   (outer wall),
   T(j,1)=T(N/2+2,M+3-j);	% lower edge (from 0 to a; using symmetry)
   T(N+3,j)=T(N+1,j);		% and right edge (x=b; also using symmetry)
end

for i=2:(N+2)		% Neumann b.c. and on upper edge (outer wall)
   T(i,M+3)=T(i,M+1)-(2*ho*dy/k)*(T(i,M+2)-To);
end

for i=(N/2+2):(N+2)	% ...and also on lower edge (bordering inner
   T(i,1)=T(i,3)-(2*hi*dy/k)*(T(i,2)-Ti);	% wall; from a to b)
end 

for n=1:maxiter		% iterate for enough times until steady state
   for i=2:(N+2)	% central finite-diff discretization
      for j=2:M+2	% del^2(T)=0 combined w/ point-Gauss-Siedel
	T(i,j)=(1/(2*(1+B)))*(T(i+1,j)+T(i-1,j)+B*(T(i,j+1)+T(i,j-1)));
      end
   end
end


for i=2:N+2		% set up plotting vectors
   X(i-1)=dx*(i-2);
end
for j=2:M+2
   Y(j-1)=(dy*(j-2));
end

mesh(X,Y,T(2:N+2,2:M+2)'-273)		% plotting, labelling, and formatting
xlabel('x (m)')
ylabel('y (m)')
zlabel('T (C)')
title('CHM2: Steady-state temp dist over x-section of chimney wall')
axis([0 0.9 0 0.45 0 400])
view(-80,40)
