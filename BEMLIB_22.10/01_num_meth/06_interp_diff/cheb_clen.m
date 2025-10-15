 function y=cheb_clen(N,c,t);

%---
% Clenshaw algorithm for Chebyshev interpolation
%---

  d(N+1) = c(N+1);
  d(N) = 2*t*c(N+1) + c(N);
  for i=N-1:-1:1
     d(i) = 2.0*t*d(i+1)-d(i+2)+c(i);
  end
  y = d(1) - t*d(2);

%---
% done
%---

  return
