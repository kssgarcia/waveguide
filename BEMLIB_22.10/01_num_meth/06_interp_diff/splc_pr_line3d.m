%=============================
% plot a closed line using
% periodic cubic splines
%=============================

N = 5;

fx = [0.0, 0.1, 0.9, 0.95 1.0];
fy = [0.1, 1.0, 1.1, 0.7 0.01];
fz = [0.2, 0.2, 0.3, 1.7 0.01];

%-----
% wrap
%-----

fx(N+1)=fx(1);
fy(N+1)=fy(1);
fz(N+1)=fz(1);

%---
% polygonal arc length
%---

beta(1) = 0.0D0;

for i=2:N+1
 beta(i)=beta(i-1)+ sqrt((fx(i)-fx(i-1))^2+(fy(i)-fy(i-1))^2+(fz(i)-fz(i-1))^2);
end

%----
% spline coefficients
%----

[ax,bx,cx] = splc_pr (N,beta,fx);
[ay,by,cy] = splc_pr (N,beta,fy);
[az,bz,cz] = splc_pr (N,beta,fz);

%====
% plotting
%====

plot3(fx,fy,fz,'or')
hold on
kpl=16; % number of plotting intervals per segment

%---
for i=1:N    % loop over segments
%---

 dbeta = (beta(i+1)-beta(i))/kpl;

 for j=1:kpl+1  % loop over plotting points
  xd = dbeta*(j-1.0);
  xx(j) = ( (ax(i)*xd+ bx(i) )*xd+ cx(i) )*xd + fx(i);
  yy(j) = ( (ay(i)*xd+ by(i) )*xd+ cy(i) )*xd + fy(i);
  zz(j) = ( (az(i)*xd+ bz(i) )*xd+ cz(i) )*xd + fz(i);
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

