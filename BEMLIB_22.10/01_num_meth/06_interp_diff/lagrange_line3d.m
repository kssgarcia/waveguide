%==========================
% Line reconstruction
% by Lagrange interpolation
%==========================

N=5;
fx=[-0.0 0.3 0.4 0.6 0.9 0.6];
fy=[0.1 0.5 -0.4 0.2 0.8 0.9];
fz=[0.3 0.9 -0.2 0.2 0.3 0.8];

%---
% polygonal arc length
%---

beta(1) = 0.0D0;

for i=2:N+1
 beta(i)=beta(i-1)+ sqrt((fx(i)-fx(i-1))^2+(fy(i)-fy(i-1))^2+(fz(i)-fz(i-1))^2);
end

%====
% plotting
%====

plot3(fx,fy,fz,'or')
hold on
kpl=16; % number of plotting intervals per segment

%---

for i=1:N    % loop over segments

 dbeta = (beta(i+1)-beta(i))/kpl;

 for j=1:kpl+1  % loop over plotting points
  xint = beta(i)+dbeta*(j-1.0);
  xx(j) = lagrange(N,beta,fx,xint);
  yy(j) = lagrange(N,beta,fy,xint);
  zz(j) = lagrange(N,beta,fz,xint);
 end

 plot3(xx,yy,zz)

end

%---
% done
%---

xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('z','fontsize',15)
set(gca,'fontsize',15)
grid
box
