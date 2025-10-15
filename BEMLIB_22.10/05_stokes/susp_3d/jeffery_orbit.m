e=4.0;
e=0.25;
e=100.0;
e=0.01;

phi=0.375;
phi=0.125;
phi=0.45;
phi=0.25;
phi=0.01;
phi=pi*phi;

% jeffery constant

c=e* cot(phi);

es=e^2;
N=64;

dchi = 2*pi/N;

for i=1:N+1
 chi=(i-1.0)*dchi;
  tanxi = c/sqrt(es*sin(chi)^2+cos(chi)^2);
 xi = atan(tanxi);
 rho = sin(xi);
 xj(i)=rho*cos(chi);
 yj(i)=rho*sin(chi);
 zj(i)=cos(xi);
end

plot3(xj,-zj,yj)
hold on
[z, y, x]= sphere;
mesh(0.99*x, 0.99*y, 0.99*z);

xlabel('x')
ylabel('z')
zlabel('y')
view(30,18)

