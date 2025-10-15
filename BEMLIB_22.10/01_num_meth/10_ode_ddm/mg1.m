clear all
close all
hold on

%=====
% Richardson iterations for the 1D Poisson eqn
% on a single grid
%=====

a=0.0;
b=1.0;
qL=0.0;
fR=0.0;
ndiv=5;
ndiv=4;
Niter=32*64*64;
Niter=64*64;
Isee=2*64;
Isee=64*64;
Isee=8888*64;

L=b-a;
N=2^ndiv;
h=L/N;

omopt=0.5/(sin(0.5*pi/(N+1))^2+sin(0.5*N*pi/(N+1))^2);
omax=0.5/sin(0.5*N*pi/(N+1));
omega = 1.00*omax;
omega = omopt;

%---
% initialize
%---

f=zeros(N+1,1);

for i=1:N+1
 x(i) = a+(i-1)*h;
 f(i)=0.1*sin(2*pi*x(i)/L);
 f(i)=0.0;
 f(i)=0.1*rand-0.05;
end

f(N+1)=fR;

plot(x, f,'ro-');

%---
% right-hand side
%---

b=zeros(N,1);
for i=1:N
 xpos = (i-1)*h/L;
 g(i)=exp(-2*xpos);
 g(i)=sin(2*pi*xpos/L);
 g(i)=1.0;
 b(i)=h*h*g(i);
end
b(1)=b(1)-2*h*qL;
b(N)=b(N)+fR;
%-----------

%--- A ---
%A=2*eye(N,N);
%for i=2:N
% A(i-1,i)=-1;
% A(i,i-1)=-1;
%end
%A(1,2)=-2;
%P=eye(N,N)-omega*A;
%r=b-A*f;
%-----------

%---
% iterations
%---

fnew=f;

Ic=0;      % counter for plotting

for i=1:Niter

 fnew(1)= (1-2.0*omega)*f(1)+2.0*omega*f(2)+omega*b(1);
 for i=2:N
  fnew(i)= omega*f(i-1)+(1-2.0*omega)*f(i)+omega*f(i+1)+omega*b(i);
 end
 f=fnew;

 if(Ic==Isee)
  plot(x, f,'o-');
  Ic=0;
 end
 Ic=Ic+1;

end

%-----------
plot(x, f,'o-');

set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('f','fontsize',15)
axis([0 1 -0.1 1])
box
