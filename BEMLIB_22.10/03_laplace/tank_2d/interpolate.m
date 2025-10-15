function [Pint,Rint] = interpolate ...
...
   (X1,X2,X3,X4,Xint ...
   ,P1,P2,P3,P4 ...
   ,R1,R2,R3,R4 ...
   )

%---
% Lagrange interpolation
% of two sets of data (P & R)
%---

 N = 3;

  x(1)=X1; x(2)=X2; x(3)=X3; x(4)=X4;
  p(1)=P1; p(2)=P2; p(3)=P3; p(4)=P4;
  r(1)=R1; r(2)=R2; r(3)=R3; r(4)=R4;

  Pint = 0.0;
  Rint = 0.0;

  for i=1:N+1
     lpoly = 1.0;
     for j=1:N+1
       if (j~=i)
          lpoly = lpoly*(Xint-x(j))/(x(i)-x(j));
       end
     end
     Pint = Pint + lpoly*p(i);
     Rint = Rint + lpoly*r(i);
  end

%----
% done
%----

return
