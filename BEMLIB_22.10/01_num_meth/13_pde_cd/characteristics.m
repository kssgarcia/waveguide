clear all
close all

%------
% plot characteristics
%------

Ntime = 32;
Ntime = 2*32;

time1=0.0;
time2=0.5;
time2=2.0;

Dt=(time2-time1)/Ntime;

ximin=-5.0;
ximax=1.2;
ximax=5.0;

Nxi=32;

Dxi=(ximax-ximin)/Nxi;

figure(1)
hold on

xi=-0;

%---
for repeat=1:Nxi+1

 xi = ximin+(repeat-1.0)*Dxi;

for i=1:Ntime+1
 t(i) = time1 +(i-1)*Dt;
 X(i) = tan(atan(xi)+t(i))-t(i);
 X(i) = xi+t(i);
 Vel = exp(-xi*xi);
 X(i) = xi+Vel*t(i);
end

 plot(X,t);

end
%---

%axis([0, 1, 0, 10])
%axis([0, 3, 0, 1.2])
%axis([-1 1 0 1])
xlabel('X','fontsize',15)
ylabel('t','fontsize',15)
set(gca,'fontsize',15)
box
