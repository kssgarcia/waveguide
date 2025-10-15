clear all

%===============
% demonstration of the
% generalized Sherman--Morrison formula
%==============

B = [1 2 3;
     2 3 4;
     1 4 5];

z1 = 1.4;
u1 = [3 4 9]';
v1 = [2 3 7]';

z2 = 3.4;
u2 = [1 2 3]';
v2 = [6 5 4]';

U(:,1) = u1;
U(:,2) = u2;

V(:,1) = v1;
V(:,2) = v2;

w1 = u1'/B';
w2 = u2'/B';

G(1,1) = v1'*w1';
G(1,2) = v1'*w2';

G(2,1) = v2'*w1';
G(2,2) = v2'*w2';

Z = zeros(2,2);
Z(1,1) = z1;
Z(2,2) = z2;


A = B + c1*u1*v1' + c2*u2*v2' + c3*u3*v3' + c4*u4*v4';

invA = inv(A)
invB = inv(B);
invA1 = invB-invB*U*inv(inv(Z)+G)*V'*invB

%---
c1*u1*v1' + c2*u2*v2' 
U*Z*V'
%----

G
G = V'*invB*U
