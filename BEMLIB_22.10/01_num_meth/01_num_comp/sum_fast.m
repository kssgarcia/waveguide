function [Sum,a,b,kmax] = sum_fast (menu,s,N,p,kmax,tol)

%===============================
% Sum an infinite series whose
% terms decay like 1/i^s
% for i=1, 2, ...
% using Aitken extrapolation
%===============================

L = 0;
M = N;

Sum = 0.0;

%-------------
for k=1:kmax+2
%-------------

   for i=L+1:M

        if(menu==1)
          term = 1.0/i^s;
        elseif(menu==2)
          term = 1.0/(i*(i+2.0));
          s = 2;
        elseif(menu==3)
          is = i^2;
          term = 1.0/(1.0/is+is);
          s = 2;
        elseif(menu==4) 
          term = 1.0/(2.0*i-1.0)^2;
          s = 2;
        else
           disp('menu option is not available')
        end

    Sum = Sum + term;

   end

   b(k) = Sum - 0.5/M^s;

   if(k>2)
    a(k-1) = (b(k-2)*b(k)-b(k-1)*b(k-1))/(b(k)-2.0*b(k-1)+b(k-2));
   end

   if(k>3)
   %---
   if(abs(a(k-1)-a(k-2))<tol)
    Sum = a(k-1);
    kmax = k-2;
    break;
   end
   %---
   end

   L = M;
   M = M*p;

%---
end
%--

Sum = a(k-1);

%---
% done
%---

return
