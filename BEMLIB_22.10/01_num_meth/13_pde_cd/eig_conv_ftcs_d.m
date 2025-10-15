clear all
close all

%-----------
% eigenvalues of the NxN projection matrix
% for the FTCS discretization of the convection equation
% with Dirichlet boundary conditions
%-----------

N=9;
M=N-1;

c = 0.2;

P=eye(M,M);

for i=1:M-1
 P(i,i+1)=-0.5*c;
 P(i+1,i)= 0.5*c;
end

%---
% numerical eigenvalues
%---

eigen=eig(P)
plot(real(eigen),imag(eigen),'o')
hold on

%---
% analytical eigenvalues
%---

for k=1:M
 lr(k)=1.0;
 li(k)=c*cos(k*pi/N);
end

plot(lr,li,'r+')

xlabel('\lambda_R','fontsize',15)
ylabel('\lambda_I','fontsize',15)
set(gca,'fontsize',15)

