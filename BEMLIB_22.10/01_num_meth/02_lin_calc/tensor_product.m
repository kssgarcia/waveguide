clear all

M = 4;
N = 5;
K = 7;

for i=1:M
 u1(i) = rand;
 u2(i) = rand;
 u3(i) = rand;
 u4(i) = rand;
 u5(i) = rand;
end

u1=u1';u2=u2';u3=u3';u4=u4';u5=u5';

A(:,1) = u1; 
A(:,2) = u2; 
A(:,3) = u3; 
A(:,4) = u4; 
A(:,5) = u5; 

for i=1:K
 v1(i) = rand;
 v2(i) = rand;
 v3(i) = rand;
 v4(i) = rand;
 v5(i) = rand;
end

B(1,:) = v1; 
B(2,:) = v2; 
B(3,:) = v3; 
B(4,:) = v4; 
B(5,:) = v5; 

A*B

u1*v1+u2*v2+u3*v3+u4*v4+u5*v5


