clear all
close all
hold on

%===================================
% steepest descent for two equations
%===================================

A=[2.70 0.10; ...
   0.02 1.00];

b=[0.1 0.2]';

x=[-0.5 3.0]';
x=[3.0 -1.5]';

%--- precondition

b = A'*b;
A = A'*A;
eigenvalues = eig(A)

%--- minimization

xplot(1)=x(1); yplot(1)=x(2);

for i=1:10
 r = -A*x+b;
 alpha = r'*r/(r'*A*r);
 x = x+alpha*r
 xplot(2)=x(1); yplot(2)=x(2);
 plot(xplot,yplot,'--o')
 xplot(1)=x(1); yplot(1)=x(2);
end

%===
% contour plot
%===

Nx=63;
Dx = 6.0/Nx;
for i=1:Nx+1
 X(i)=-3.0+(i-1.0)*Dx;
end

Ny=64;
Dy = 6.0/Ny;
for j=1:Ny+1
 Y(j)=-3.0+(j-1.0)*Dy;
end

for i=1:Nx+1
 for j=1:Ny+1
  Xgrid=[X(i) Y(j)]';
  Z(j,i) = 0.5* Xgrid'*A*Xgrid - b'*Xgrid; 
 end
end

contour(X,Y,Z)
xlabel('x_1','fontsize',15)
ylabel('x_2','fontsize',15)
axis([-3 3 -3 3])
axis square
set(gca,'fontsize',15)
box
