clear all
close all

%-----------
% eigenvalues of the NxN projection matrix
% for the FTBS discretization of the convection equation
% with periodic boundary conditions
%-----------

N=10;
c =0.5;

P=eye(N,N);

for i=1:N-1
 P(i,i)=1-c;
 P(i+1,i)= c;
end

P(N,N)=1-c; P(1,N)=c;

%---
% numerical eigenvalues
%---

eigen=eig(P)
plot(real(eigen),imag(eigen),'o')
hold on
axis equal

%---
% analytical eigenvalues
%---

hold on
for i=1:N
 rad = c+0.00;
 lr(i)=1-c+rad*cos(i*2*pi/N);
 li(i)=    rad*sin(i*2*pi/N);
end

plot(lr,li,'r+')

xlabel('\lambda_R','fontsize',15)
ylabel('\lambda_I','fontsize',15)
set(gca,'fontsize',15)

