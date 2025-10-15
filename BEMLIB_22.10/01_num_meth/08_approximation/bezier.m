%---
% plot a bezier segment
%---

%---
clear all
close all
hold on
%---

Npl=64;
Dt = 1/Npl;

%---
%M=3;
%ydata=[0 1 1.6 0.8];
%xdata=[0 2 3 0.5];
%xdata=[0 0.5 2 3];
%xdata=[0 2 0.5 3];
%---

%---
M=2;
ydata=[0 0.75 1];
xdata=[0 1.5 2];
%---

%---
for k=1:Npl+1
 t=(k-1)*Dt;
 bern = bernstein(M,t);
 x(k) = 0.0;
 y(k) = 0.0;
 for j=1:M+1
  x(k)=x(k)+xdata(j)*bern(M,j);
  y(k)=y(k)+ydata(j)*bern(M,j);
 end

end
%---

plot(x,y)
plot(xdata,ydata,'--o')

%xlabel('x','fontsize',15)
%ylabel('y','fontsize',15)
set(gca,'fontsize',15)
%box
axis off
