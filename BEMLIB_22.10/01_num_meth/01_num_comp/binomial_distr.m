function Bnm = binomial_distr(n,m,p)

%---
% compute the binomial distribution
%---

% kmax = min(m,n-m);
 kmax = m;
 if((n-m)==kmax) kmax = n-m; end

 combo = 1.0;

 for k=1:kmax
  combo = combo*(n-k+1)/k;
 end

 Bnm = combo* p^m *(1.0-p)^(n-m);

%---
% done
%---

 return
