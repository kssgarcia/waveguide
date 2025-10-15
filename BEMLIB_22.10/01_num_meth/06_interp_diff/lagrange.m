function yint = lagrange (N,x,y,xint)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% perform Lagrange interpolation
% on a specified set of N+1 points
%
% (x, y): data
% xint: interpolated x value
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 yint = 0.0;

 for i=1:N+1
    lpoly = 1.0;
    for j=1:N+1
      if (j~=i)
         lpoly = lpoly*(xint-x(j))/(x(i)-x(j));
      end
    end
    yint = yint + lpoly*y(i);
 end

%----
% done
%----

return
