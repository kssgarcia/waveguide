%-----
% eigenvalues of a hessenberg matrix
%---

a=1.2;
N=4;
A=eye(N,N);
for i=2:N
 A(i,i-1)=1.0;
end

%---
% matrix
%---

for i=1:N
 for j=i:N
  A(i,j)=a^(j-i+1);
 end
end
%eig(A)
[V,D]=eig(A)

%---
% exact
%---

m=floor(N/2);
for i=1:N-m
lam(i)=4*a*cos(i*pi/(N+2))^2;
end
lam
