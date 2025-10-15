clear all

%===============
% demonstrattion
% of the Sherman--Morrison formula
%==============

B = [1 2 3;
     2 3 4;
     1 4 5];

u = [3 4 9]';
v = [2 3 7]';

A = B+u*v';

invA = inv(A);
invB = inv(B);

s = v'*invB*u;

invA
invA1 = invB-invB*u*v'*invB/(1.0+s)

C = invA*u*v';
Omega = invB*u*v';
D = Omega-Omega^2/(1.0+s);

C
D
