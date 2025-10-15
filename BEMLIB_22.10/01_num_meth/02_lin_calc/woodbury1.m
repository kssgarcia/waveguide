clear all

%==============================
% demonstration of the generalized
% Woodbury formula
%==============================

N = 3;
p = 2;
K = 1;

A = [1 2 3;
     2 3 4;
     1 4 5];

u1 = [3 4 9]';
v1 = [2 3 7]';

u2 = [4 1 8]';
v2 = [2 2 8]';

U  =[u1  u2]
Vt =[v1'; v2']

B = A+u1*v1'+u2*v2';

invA = inv(A);
invB = inv(B);

G = Vt*invA*U; 

invB1 = invA-invA*U*inv(eye(p*K,p*K)+G)*Vt*invA;

invB
invB1
