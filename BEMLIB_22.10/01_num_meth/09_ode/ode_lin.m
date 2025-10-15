%close all
clear all
hold on

%========================
% integration of linear ODEs
% using different methods
%========================

%---
% definitions
%---

N=2;

A =[-1 2;-2 0];
A =[-1 0;0 -1];

b =[3 4]';
b =[0 0]';

x0=[0 0]';
x0=[1 1]';

eig(A)

tint = 12.0;  % time interval

Dt =     2*2*2*0.1;
Dt =       2*2*0.1;
Dt =         2*0.1;
Dt = 2*2*2*2*2*0.1;
Dt =   1.5;
Dt =   2.0;
Dt =   2.5;
Dt =   1.0;
Dt =   0.5;
Dt =   0.1;
Np=floor(tint/Dt)

x(:,1) = x0 + Dt*(A*x0-b);
time(1)=Dt;

%-----
% euler
%-----

%for k=1:Np
% x(:,k+1) = x(:,k) + Dt*(A*x(:,k)-b);
% time(k+1)=time(k)+Dt;
%end

%-----
% mid-point
%-----

x(:,2) = x0 + 2*Dt*(A*x(:,1)-b);
time(2)=2*Dt;
for k=2:Np
 x(:,k+1) = x(:,k-1) + 2*Dt*(A*x(:,k)-b);
 time(k+1)=time(k)+Dt;
end

%-----
% plotting
%-----

%plot(x0(1),x0(2),'s')
%plot(x(1,:),x(2,:),'-o')
%xlabel('x_1','fontsize',15)
%ylabel('x_2','fontsize',15)

plot([0 time(1)],[x0(1),x(1,1)],'-o');plot(time,x(1,:),'-o')
xlabel('t','fontsize',15)
ylabel('x','fontsize',15)
axis([0 10 -2 2])

set(gca,'fontsize',15)
