function b=inv_l(N,a)
                                                                                
%-----------------------------------------
% Inverse of a lower triangular matrix "a"
%-----------------------------------------

%---
% initialize and compute the diagonals
%---
                                                                                
for i=1:N
 b(i,i) = 1.0/a(i,i);
 for j=i+1:N
   b(i,j) = 0.0;
 end
end

%---
% fill in the inverse
%---

for j=1:N-1
  for i=j+1:N
     sum = 0.0;
     for m=j:i-1
        sum = sum +a(i,m)*b(m,j);
      end
      b(i,j) = - sum/a(i,i);
  end
end
                                                                                
%-----
% done
%-----

return;
