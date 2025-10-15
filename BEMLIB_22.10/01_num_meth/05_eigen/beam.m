clear all

%====================
% Euler beam buckling
%====================

N=32;
N=16;
N=8;

%---
% pentadiagonal matrix
%---

P = zeros(N-1,N-1);
P(1,1)= 5.0; P(1,2)=-4.0; P(1,3)= 1.0; 
P(2,1)=-4.0; P(2,2)= 6.0; P(2,3)=-4.0; P(2,4)= 1.0; 

for i=3:N-3
 P(i,i-2)=1; P(i,i-1)=-4; P(i,i)=6; P(i,i+1)=-4; P(i,i+2)=1;
end

P(N-2,N-4)= 1.0; P(N-2,N-3)=-4.0; P(N-2,N-2)= 6.0; P(N-2,N-1)=-4.0;
                 P(N-1,N-3)= 1.0; P(N-1,N-2)=-4.0; P(N-1,N-1)= 5.0;
%---
% tridiagonal matrix
%---

T = zeros(N-1,N-1);
T(1,1) = 2.0; T(1,2) =-1.0;

for i=2:N-2
 T(i,i-1)=-1.0; T(i,i)=2.0; T(i,i+1)=-1;
end
T(N-1,N-2)=-1.0; T(N-1,N-1)=2.0;

%---
% eigenvalues
%---

egv = eig(inv(T)*P)*N^2/pi^2;

sort(egv)
