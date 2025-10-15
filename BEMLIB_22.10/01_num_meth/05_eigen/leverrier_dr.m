clear all
close all

%---
% driver for coefficients
% of the characteristic polynomial
%---

N = 10;
N = 4;
N = 8;
N = 2;

for i=1:N
 for j=1:N
  A(i,j) = (-1)^(i+j)/(i+j-1);
  A(i,j) = 1/(i+j-1);
  A(i,j) = rand;
 end
end

c = leverrier(N,A)

%---
% roots of the characteristic polynomial
%---

cshifted(1) = 1.0;
for i=1:N
  cshifted(i+1) = c(i);
end

[eig(A) roots(cshifted)]

%---
% Cayley--Hamilton
%---

CH = c(N)*eye(N);
Z = A;

for i=N-1:-1:1
 CH = CH + c(i)*Z;
 Z=Z*A;
end

CH = (-1)^N*(CH+Z);

% [A^2 trace(A)*A-det(A)*eye(N)]

