function Cnm = binomial_coeff(n,m)

%---
% compute the binomial coefficient
%---

% kmax = min(m,n-m);
  kmax = m;
  if((n-m)<kmax) kmax = n-m; end

 Cnm = 1.0;

 for k=1:kmax
  Cnm = Cnm*(n-k+1)/k;
 end

%---
% done
%---

 return
