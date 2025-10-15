%---------
% quadratic B-spline reconstruction of a closed line
%-------

close all
clear all
axis square
hold on
axis([-0.1 1.1 -0.1 1.1])
box
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)

%---
% influence points
%---

N=4;
x=[-0.05 1.0 1.0 0.0];
y=[0.0 0.1 1.0 0.8];

x(N+1)=x(1); y(N+1)=y(1);   % wrap
x(N+2)=x(2); y(N+2)=y(2);
x(N+3)=x(3); y(N+3)=y(3);
plot(x,y,'o')

%---
% knots
%---

for i=1:N
  knotx(i)=0.5*(x(i)+x(i+1));
  knoty(i)=0.5*(y(i)+y(i+1));
end
plot(knotx,knoty,'s')

%---
% plotting
%---

npl=16; % points per segment

Dxi=1/npl;
Ic=0;

for i=2:5
 Ic=0;
for j=1:npl+1
 xi=(j-1.0)*Dxi;
 b2A=0.5*(1.0-xi)^2;
 b20=0.5*(1.0+2.0*xi-2.0*xi^2);
 b21=0.5*xi^2;
 Ic=Ic+1;
 xline(Ic) = x(i-1)*b2A+x(i)*b20+x(i+1)*b21;
 yline(Ic) = y(i-1)*b2A+y(i)*b20+y(i+1)*b21;
end
 if(i==2)
 plot(xline,yline)
 elseif(i==3)
 plot(xline,yline,'--')
 elseif(i==4)
 plot(xline,yline,':')
 elseif(i==5)
 plot(xline,yline,'-.')
 end

end
