clear all
close all

%-----------
% eigenvalues of the NxN projection matrix
% for the FTCS discretization of the convection equation
% with periodic boundary conditions
%-----------

N=3;

c =0.2;

P=eye(N,N);

for i=1:N-1
 P(i,i+1)=-0.5*c;
 P(i+1,i)= 0.5*c;
end

P(1,N)=0.5*c;
P(N,1)=-0.5*c;

%---
% numerical eigenvalues
%---

eigen=eig(P)
plot(real(eigen),imag(eigen),'o')
hold on

%---
% analytical eigenvalues
%---

for k=1:N
 lr(k)=1.0;
 li(k)=c*sin(k*2*pi/N);
end

plot(lr,li,'r+')

xlabel('\lambda_R','fontsize',15)
ylabel('\lambda_I','fontsize',15)
set(gca,'fontsize',15)

