function Bern = bernstein (M,t)

%---
% compute the Bernstein polynomials
%---

%----------------
Bern=zeros(M,M+1);
%----------------

Bern(1,1) = 1.0-t;
Bern(1,2) = t;

if(M>1)

Bern(2,1) = (1.0-t)^2;
Bern(2,2) = 2*(1.0-t)*t;
Bern(2,3) = t^2;

if(M>2)

for l=3:M
  Bern(l,1) = (1.0-t)*Bern(l-1,1);
  for i=2:M+1
      Bern(l,i) = (1.0-t)*Bern(l-1,i)+t*Bern(l-1,i-1);
  end
end

end
end

return
