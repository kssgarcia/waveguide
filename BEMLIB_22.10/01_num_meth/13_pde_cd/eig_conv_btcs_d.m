clear all
close all

%-----------
% eigenvalues of the NxN projection matrix
% for the FTBS discretization of the convection equation
% with Dirichlet boundary conditions
%-----------

N=10;
c=0.2;

M=N-1;

P=eye(M,M);

for i=1:M-1
 P(i,i)=1-c;
 P(i+1,i)= c;
end
P(M,M)=1-c;

[U,eigen]=eig(P)
plot(real(eigen),imag(eigen),'o')
hold on
axis equal

hold on
for i=1:N
 rad = c+0.00;
 lr(i)=1-c+rad*cos(i*2*pi/N);
 li(i)=rad*sin(i*2*pi/N);
end

plot(lr,li,'rs')
