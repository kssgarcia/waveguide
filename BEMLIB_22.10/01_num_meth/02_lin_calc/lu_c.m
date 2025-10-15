%---
% Crout LU decomposition
%---
n=4;
a=[ 3 2 3 7 ; 3 4 9 5 ; 4 5 6 5 ; -1 2 3 4 ];
%---

  for i=1:n
    u(i,i)=1.0;
    l(i,1)=a(i,1);
    u(1,i)=a(1,i)/l(1,1);
  end

%---
   for k=2:n
%---

   for i=k:n
     l(i,k)=a(i,k);
     for m=1:k-1
       l(i,k) = l(i,k) - l(i,m)*u(m,k);
     end
   end

   for j=k+1:n
      u(k,j) = a(k,j);
      for m=1:k-1
         u(k,j) = u(k,j) - l(k,m)*u(m,j);
      end
      u(k,j) = u(k,j)/l(k,k);
   end

%---
end
%---

a
l
u
l*u

d=diag(diag(l));
lD = l*inv(d)
uD = d*u

