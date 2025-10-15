clear all
close all

%---
% eigenvalues of the Gauss-Seidel projection matrix
% for the diffusion matrix
%---

N = 11;
A = 2*eye(N,N);
for i=2:N
 A(i,i-1) = -1;
end
B = inv(A)

U=zeros(N,N);
for i=1:N-1
 U(i,i+1) = -1;
end
D = B*U;
lam = eig(D)

for i=1:N
 xl(i) = i+1;
 l(i) = - cos(i*pi/(N+1))^2;
 xm(i) = i;
end
l'

figure(1)
hold on
plot(xl,l,'ko')
plot(xm,real(lam),'r+')
plot(xm,imag(lam),'c*')
