function edouble = mg_prolongate1(N,e)

%---
% double the grid
%
% prolongate the error(e)
%---

edouble(1)=e(1);

Ic=1;
for i=2:N+1
 Ic=Ic+1;
 edouble(Ic)=0.5*(e(i)+e(i-1));
 Ic=Ic+1;
 edouble(Ic)=e(i);
end

return

