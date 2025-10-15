function enew = mg_rich(Niter,omega,N,e,r)

%---
% perform Niter Richardson iterations
% on the equation A e= r
%---

for i=1:Niter
 enew(1)= (1-2.0*omega)*e(1)+2.0*omega*e(2)+omega*r(1);
 for i=2:N
  enew(i)= omega*e(i-1)+(1-2.0*omega)*e(i)+omega*e(i+1)+omega*r(i);
 end
 enew(N+1)=0.0;
 e=enew;
end

%---
% done
%---

return

