function f = newton2_fun(menu,x)

%---
% function evaluation for Newton's method
%---

 if(menu==1)
   f(1) = x(1)^2+2.0*x(2)^2 -9.0;
   f(2) = x(1)*x(2) -2.0;
 else
   f(1) = x(1)+x(2) -3.0;
  f(2) = x(1)-x(2) -2.0;
 end

%---
% done
%---

return
