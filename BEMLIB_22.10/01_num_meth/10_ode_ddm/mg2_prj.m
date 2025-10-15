nu1=3;
nu2=3;
omega=1.0/3.0
N=8;

%---
% fine grid matrix
%---

Nsys=N;
A=zeros(Nsys,Nsys);

for i=1:Nsys
 A(i,i)=2;
end
for i=1:Nsys-1
 A(i,i+1)=-1;
 A(i+1,i)=-1;
end
A(1,2)=-2;

%---
% coarse grid matrix
%---

Nsys=N/2;
Ap=zeros(Nsys,Nsys);

for i=1:Nsys
 Ap(i,i)=2;
end
for i=1:Nsys-1
 Ap(i,i+1)=-1;
 Ap(i+1,i)=-1;
end
Ap(1,2)=-2;

%---
% restriction matrix
%---

M=zeros(N/2,N);
M(1,1)=2;
M(1,2)=2;
for i=2:N/2
 M(i,2*i-2)=1;
 M(i,2*i-1)=2;
 M(i,2*i)=1;
end

%---
% prolongation matrix
%---

m=zeros(N,N/2);
m(1,1)=2;
m(2,1)=1;
for j=2:N/2
 m(2*j-2,j)=1;
 m(2*j-1,j)=2;
 m(2*j,j)=1;
end
m=0.5*m;

%---
% fine-grid Richardson projection matrix
%---

I=eye(N,N);
P=I-omega*A;

%---
% effective projection matrix
%---

Peff = P^nu1*(I - m*inv(Ap)*M*A) * P^nu2;
eig(Peff)
eig(P^(nu1+nu2))

