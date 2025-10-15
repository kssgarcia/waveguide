function d = ndm(N,x,j)

%---
% node differentiation matrix
%---

 for i=1:N+1
  if(i~=j)

   d(i)=1.0;
   for k=1:N+1
     if(k~=j)
      d(i) = d(i)*(x(j)-x(k));
     end
     if(k~=i)
      d(i) = d(i)/(x(i)-x(k));
     end
   end
   d(i)=d(i)/(x(j)-x(i));

  end
 end

 d(j) = 0.0;
 for k=1:N+1
  if(k~=j)
   d(j) = d(j)+1.0/(x(j)-x(k));
  end
 end

%---
% Done
%---

return
