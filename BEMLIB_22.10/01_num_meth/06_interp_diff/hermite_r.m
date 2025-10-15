function r = r(xint,N,x,m,i,k)

%---
% evaluate the function r_{ik}(xint)
%---

r=1.0;

for j=1:N+1
  if(j~=i)
   r = r*(xint-x(j))^(1+m(j))/(x(i)-x(j))^(1+m(j));
  end
end

if(k>1)
 for q=1:k-1
  r = r*(xint-x(i))/q;
 end
end

%-----
% done
%-----

return
