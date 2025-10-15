clear all

%================
%
% FDLIB
%
% eigenvalues and eigenvectors of a
% tridiagonal Toeplitz matrix A
%=============

N = 6;
a = 1.2;
b = 3.4;
c = 0.4;

%---
% construct the matrix
%---

A = a*eye(N,N);

for i=1:N-1
 A(i,i+1)=b;
end

for i=2:N
 A(i,i-1)=c;
end

[eivec,eival] = eig(A);

eival

%---
% exact solution
%---

d = sqrt(b*c);
f = sqrt(c/b);

for i=1:N

 l(i) = a - 2*d*cos(i*pi/(N+1));

 sum = 0;

 for j=1:N
  u(j,i)=f^j * sin(j*i*pi/(N+1));
  u(j,i)=f^j * sin(j*(N+1-i)*pi/(N+1));
  u(j,i)=f^j * (-1)^j* sin(j*i*pi/(N+1));
  sum = sum+u(j,i)^2;
 end

 sum = sqrt(sum);

 for j=1:N
  u(j,i) = u(j,i)/sum;
 end

end

l'
