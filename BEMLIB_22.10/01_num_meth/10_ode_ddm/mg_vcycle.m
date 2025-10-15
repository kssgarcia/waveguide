function e = mg_vcycle(N,e,r,nu1,nu2,omega)

%---------
 if(N>2)
%----------

 e = mg_rich(nu1,omega,N,e,r);   % presmoothing

 r(1)= r(1)-2.0*e(1)+2.0*e(2);   % new residual
 for i=2:N
   r(i)= r(i)+e(i-1)-2.0*e(i)+e(i+1);
 end

 rhalf = mg_restrict1(N,r);      % restrict
 ehalf=zeros(N/2+1,1);           % initial guess
 ehalf = mg_vcycle(N/2,ehalf,rhalf,nu1,nu2,omega);
 e = e + mg_prolongate1(N/2,ehalf);

%---------
 else
%----------

  e(1) = r(1)+r(2);
  e(2) = 0.5*(r(1)+2.0*r(2));
  e(3) = 0.0;

%----------
 end
%----------

%---
% post-smoothing
%---

if(nu2>0)
 e = mg_rich(nu2,omega,N,e,r);
end

return
