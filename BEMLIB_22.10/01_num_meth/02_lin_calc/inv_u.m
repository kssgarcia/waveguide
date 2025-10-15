function b=inv_u(N,a)
                                                                                
%------------------------------------------
% inverse of an upper triangular matrix "a"
%
% Algorithm (2.5.8)
%------------------------------------------

%---
% initialize the inverse
% and compute the diagonals
%---
                                                                                
for i=1:N
 b(i,i) = 1.0/a(i,i);
 for j=1:i-1
   b(i,j) = 0.0;
 end
end

%---
% fill in
%---

for j=2:N
  for i=j-1:-1:1
     sum = 0.0;
     for m=i+1:j
        sum = sum +a(i,m)*b(m,j);
      end
      b(i,j) = -sum/a(i,i);
  end
end

%-----
% done
%-----

return;
