%---------
% cubic B-spline reconstruction of a closed line
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
% knots are influence points
%---

N=4;
x=[-0.05 1.0 1.0 0.0];
y=[0.0 0.1 1.0 0.8];

x(N+1)=x(1); y(N+1)=y(1);   % wrap
x(N+2)=x(2); y(N+2)=y(2);
x(N+3)=x(3); y(N+3)=y(3);

plot(x,y,'o')

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
 b3A=(1.0-xi)^3/6.0;
 b30=(4.0-6.0*xi^2+3.0*xi^3)/6.0;
 b31=(1.0+3.0*xi+3.0*xi^2-3.0*xi^3)/6.0;
 b32=xi^3/6.0;
 Ic=Ic+1;
 xline(Ic) = x(i-1)*b3A+x(i)*b30+x(i+1)*b31+x(i+2)*b32;
 yline(Ic) = y(i-1)*b3A+y(i)*b30+y(i+1)*b31+y(i+2)*b32;
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
