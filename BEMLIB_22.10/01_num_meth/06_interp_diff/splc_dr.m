close all
clear all

%====================
% plot a cubic spline
%====================

%------------------------------------------------
% This program accompanies the book:
%           C. Pozrikidis
% ``Numerical Computation in Science and Engineering''
%      Oxford University Press
%------------------------------------------------

%---
% sample data
%---


N = 5;
x = [0.0, 0.1, 0.6, 0.9, 1.0, 1.4];

y = [0.1, 1.0, 0.9, 1.1, 0.01, -0.3];
y = [0.1, 1.0, 0.9, 1.1, 0.01, 0.1];

% to test a true spline

% for i=1:N+1
%  y(i)=x(i)^3;
% end

kpl = 16; % number of plotting intervals per segment

%---
% prepare
%---

figure(1)
hold on
plot(x,y,'o','markersize',5,'color','red')
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
%axis([x(1) x(N+1) -1 2])
box

%----
%for Ichoose=1:6
 for Ichoose=5:5
%---

%----
% spline coefficients
%----

if(Ichoose==1)
 [a,b,c] = splc_pr (N,x,y);
 title('Periodic spline','fontsize',15)
elseif(Ichoose==2)
 [a,b,c] = splc_clm (N,x,y,0.0,0.0);
 title('Clamped spline','fontsize',15)
elseif(Ichoose==3)
 [a,b,c] = splc_nt (N,x,y);
 title('Natural spline','fontsize',15)
elseif(Ichoose==4)
 [a,b,c] = splc_true (N,x,y);
elseif(Ichoose==5)
 [a,b,c] = splc_nt_true (N,x,y);
 title('True spline','fontsize',15)
elseif(Ichoose==6)
 [a,b,c] = splc_clm_true (N,x,y,0.0);
end

%====
% plotting
%====

%---
for i=1:N    % loop over segments
%---

 dxint = (x(i+1)-x(i))/kpl;

 for j=1:kpl+1  % loop over plotting points
  xd = (j-1.0)*dxint;
  yplot(j) = ( (a(i)*xd+ b(i) )*xd+ c(i) )*xd + y(i);
  xplot(j) = x(i)+xd;
%  yexct(j) = xplot(j)^3;
 end

 plot(xplot,yplot,'k')

% plot(xplot,yexct,'r:')

end

%---
end % over Ichoose
%---

%---
% done
%---

