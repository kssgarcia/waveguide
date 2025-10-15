function f = newton3_fun(menu,x)

if(menu==1)
  f(1) = x(1)^2+2.0*x(2)^2-x(3)^2 -8.0;
  f(2) = x(1)*x(2)*x(3) -2.0;
  f(3) = x(1)*x(2)+x(1)*x(3) -3.0;
else
  disp('fun3: this option is not implemented')
end

%---
% done
%---

return

