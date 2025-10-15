%----
% Node differentiation matrix
% for the second derivative
%---

%N=8;
%x=[1 2 3 4 5 6 7 8 9];
%j=5;

%N=5;
%x=[1 2 3 4 5 6];
%j=6;

N=4;
x=[1 2 3 4 5];
j=6;

%N=2;
%x=[4 5 6];

for j=1:N+1
 D = ndm(N,x,j);
 for i=1:N+1
  d(i,j)=D(i);
 end
end

j=3;
for k=1:N+1
 d2(k)=0.0;
 for i=1:N+1
  d2(k)=d2(k)+d(k,i)*d(i,j);
 end
end

d2

