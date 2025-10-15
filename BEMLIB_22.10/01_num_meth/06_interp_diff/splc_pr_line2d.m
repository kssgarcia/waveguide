close all
clear all

%=============================
% plot a periodic cubic spline
% line in the xy plane
%=============================

N = 4;

%---
% point coordinates
%---

fpx = [0.0, 0.1, 0.9, 1.0];
fpy = [0.1, 1.0, 1.1, 0.01];

%---
% plot the data
%---

plot(fpx,fpy,'or')
hold on

%-----
% wrap
%-----

fpx(N+1) = fpx(1);
fpy(N+1) = fpy(1);

%---
% polygonal arc length
% used as the interpolation variable
%---

xp(1) = 0.0D0;

for i=2:N+1
 xp(i)=xp(i-1)+ sqrt((fpx(i)-fpx(i-1))^2+(fpy(i)-fpy(i-1))^2);
end

%----
% spline coefficients
%----

[ax,bx,cx] = splc_pr(N,xp,fpx);
[ay,by,cy] = splc_pr(N,xp,fpy);

%====
% plotting
%====

kpl = 16; % number of plotting intervals per segment

%---
for i=1:N    % loop over segments
%---

dxint = (xp(i+1)-xp(i))/kpl;

 xint = xp(i); % initial point

 for j=1:kpl+1  % loop over plotting points
  xd = xint-xp(i);
  xx(j) = ( (ax(i)*xd+ bx(i) )*xd+ cx(i) )*xd + fpx(i);
  yy(j) = ( (ay(i)*xd+ by(i) )*xd+ cy(i) )*xd + fpy(i);
  xint = xint+dxint;
 end

 plot(xx,yy)

end

%===
% done
%===

xlabel('x')
ylabel('y')
