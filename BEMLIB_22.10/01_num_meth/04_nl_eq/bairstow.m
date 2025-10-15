clear all

%======================
% Find the roots of a polynomial
% using Bairstow's method
%======================

% order of the polynomial

n = 5;

% polynomial coefficients

a=[1 1 -1 -1 1 1];

relax = 1.0;      % relaxation parameter
maximum = 100;  % maximum number of iterations

%---
% lopp
%---

k = 1;     % number of roots

while n>2     % loop over reduced-order polynomials

      flag = 0;
      j = 0;   % counter

      r = 2.0;    % initial guess
      s = 2.0;    % initial guess

% inner loop for Newton iteration on synthetic division

      while flag==0  % loop until convergence

        j = j+1;

        b(1) = a(1);
        b(2) = a(2)+r*b(1);

        for i=3:n+1
            b(i) = a(i)+r*b(i-1)+s*b(i-2);
        end

        c(1) = 0.0;
        c(2) = b(1);
        c(3) = b(2)+r*c(2);

        for i=4:n+1
            c(i) = b(i-1)+r*c(i-1)+s*c(i-2);
        end

%  check for convergence and increment (r,s)

        if(abs(b(n))<eps & abs(b(n+1))<eps)
            flag = 1;
        else
            fact =  c(n)^2-c(n-1)*c(n+1);
            dr   = -( b(n)*c(n)-c(n-1)*b(n+1) )/fact;
            ds   = -( c(n)*b(n+1)-c(n+1)*b(n) )/fact;
            r    = r + relax*dr;
            s    = s + relax*ds;
        end

       if j==maximum
            break
       end

      end

% given (r,s), find the next pair
% of roots of the polynomial

% roots of x^2-rx-s:

   x(k)=(r+sqrt(r^2+4*s))/2;
   k=k+1;
   x(k)=(r-sqrt(r^2+4*s))/2;
   k=k+1;

% reduce the size of the working polynomial
% and update the coefficient vector

   n = n-2;

   clear a
   for i=1:n+1
     a(i) = b(i);
   end

end   %  end of outer loop over polynomial order

if(n==1)
    x(k)=-b(2)/b(1);
else
    x(k)=(-b(2)+sqrt(b(2)^2+4*b(1)*b(3)))/(2*b(1));
    k=k+1;
    x(k)=(-b(2)-sqrt(b(2)^2+4*b(1)*b(3)))/(2*b(1));
    k=k+1;
end

%---
% done
%---

disp('roots:')
disp(x')
