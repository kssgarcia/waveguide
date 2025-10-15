%---
% LU decomposition (compact)
%---

n =4;
a = [ 3 2 3 5 ; 3 4 5 5 ; 4 5 6 5 ; -1 2 3 4 ];
%---

  for i=2:n
    a(i,1) = a(i,1)/a(1,1);
  end

%---
   for k=2:n
%---

   for j=k:n
     a(k,j)=a(k,j);
     for m=1:k-1
       a(k,j) = a(k,j) - a(k,m)*a(m,j);
     end
   end

   for i=k+1:n
      a(i,k) = a(i,k);
      for m=1:k-1
         a(i,k) = a(i,k) - a(i,m)*a(m,k);
      end
      a(i,k) = a(i,k)/a(k,k);
   end

%---
end
%---

a

