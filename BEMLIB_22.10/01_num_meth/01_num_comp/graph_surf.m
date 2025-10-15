clear all
close all

%=======================
% 3D-graph of a function
%  z = f(x,y)
% and its Maclaurin series
%=======================

Nx = 12;
Ny = 16;

xmin = -1.0;
xmax =  1.0;

ymin = -1.0;
ymax =  1.0;

Dx = (xmax-xmin)/Nx;
Dy = (ymax-ymin)/Ny;

for i=1:Nx+1
 x(i) = xmin+(i-1)*Dx;
end

for i=1:Ny+1
 y(i) = ymin+(i-1)*Dy;
end

for i=1:Nx+1
 for j=1:Ny+1
   tmp = x(i)*y(j);
   z(j,i) = exp(tmp);
   zlin(j,i) = 1.0 + tmp;
   zquad(j,i) = 1.0 + tmp + 0.5*tmp^2;
 end
end

figure(1)
hold on
surf(x,y,z)
mesh(x,y,zlin)
mesh(x,y,zquad)
xlabel('x')
ylabel('y')
zlabel('z')
box on
view(60,30)
