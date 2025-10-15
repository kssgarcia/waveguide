clear all

%----
% Node differentiation matrix
% for the third derivative
%---

%N=8;
%x=[1 2 3 4 5 6 7 8 9];
%j=5;

%N=5;
%x=[1 2 3 4 5 6];
%j=6;

%N=3;
%x=[4 5 6 7];

N=6;
x=[1 2 3 4 5 6 7];

%-- compute d

for j=1:N+1
 D = ndm(N,x,j);
 for i=1:N+1
  d(i,j)=D(i);
 end
end

%-- compute d2

for k=1:N+1
 for j=1:N+1
  d2(k,j)=0.0;
  for i=1:N+1
   d2(k,j)=d2(k,j)+d(k,i)*d(i,j);
  end
 end
end

%-- compute d2

j=4;

for k=1:N+1
  d3(k)=0.0;
  for i=1:N+1
   d3(k)=d3(k)+d2(k,i)*d(i,j);
  end
end

d3

