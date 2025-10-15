close all
clear all
hold on

%==================================
% solution of system of ODEs
% using the Pozrikidis/Prony method
%
% the starting series is produced
% by Taylor series of exp(tA)
%==================================

%---
% definitions
%---

N = 3;
A = [-1 2 1;-2 -5 -1; -3 4 -5];
b = [3 4 1]';
x0 = [0 0 0]'


N = 2;
A = [-1 2;-2 0];
b = [3 4]';
x0 =[0 0]'

N =2;
A =[-1  1;
     0 -1];
b =[3 4]';
x0 =[0 0]'

Dt=0.4;
Dt=0.2;
Dt=0.1;
Dt=0.01;
Dt=0.05;

tint = 12.0;  % time interval
Np=floor(tint/Dt)

%---
% generate the starting sequence
%---

X=b'/A'
I = eye(N,N);
A2 = A*A;
A3 = A2*A;
A4 = A3*A;

Dx0 = x0-X';

for i=1:N+1
 t = i*Dt;
 x(:,i) = X'+(I+t*A+t^2*A2/2.0+t^3*A3/6.0+t^4*A4/24.0)*Dx0;
end

%---
% solve for alphas
%---

mat(:,1)=x(:,1)-x0;
for i=2:N
  mat(:,i)=x(:,i)-x(:,i-1);
end
rhs = -x(:,N+1)+x(:,N)
sol = rhs'/mat'

al(1) = 1.0;
al0 = 1.0;
for i=1:N
 al(N+2-i) = sol(i);
 al0 = al0+sol(i);
end

%---
% continue the time series
%---

for j=N+2:Np
   x(:,j) = al0*X';
  for q=1:N
    x(:,j) = x(:,j) - al(q+1)* x(:,j-q);
  end
end

%---
% phase space portrait
%---

plot(x(1,1:N+1),x(2,1:N+1),'ro')
plot(x(1,:),x(2,:),'k.')
%plot(x(1,:),x(2,:),'--')
xlabel('x_1','fontsize',15)
ylabel('x_2','fontsize',15)
set(gca,'fontsize',15)
box
