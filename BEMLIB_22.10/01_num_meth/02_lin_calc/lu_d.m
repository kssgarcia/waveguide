%===========================
% Doolittle LU decomposition
%===========================


%---
n=4;
a=[ 3 2 3 7 ; 3 4 9 5 ; 4 5 6 5 ; -1 2 3 4 ];
%---

  for i=1:n
    l(i,i)=1.0;
    u(1,i)=a(1,i);
    l(i,1)=a(i,1)/u(1,1);
  end

%---
   for k=2:n
%---

   for j=k:n
     u(k,j)=a(k,j);
     for m=1:k-1
       u(k,j) = u(k,j) - l(k,m)*u(m,j);
     end
   end

   for i=k+1:n
      l(i,k) = a(i,k);
      for m=1:k-1
         l(i,k) = l(i,k) - l(i,m)*u(m,k);
      end
      l(i,k) = l(i,k)/u(k,k);
   end

%---
end
%---

a
l
u
l*u
