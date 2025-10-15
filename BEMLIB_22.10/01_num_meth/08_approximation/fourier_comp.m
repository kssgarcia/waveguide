function [M,af,bf]=ftcomp(N,L,f)

%---
% prepare
%---

k  = 2*pi/L;
h  = L/N;
kh = 2.0*pi/N;
ff = 2.0/N;

%-------------------
% even/odd intervals
%-------------------

if(mod(N,2)==0) % even
  M  = N/2;
  fc = 0.50;
else    % odd
  M = (N-1)/2
  fc = 1.0;
end

%--------
% launch
%--------

for p=0:M         % loop over coefficients

   af(p+1) = 0.5D0*(f(1)+f(N+1));
   bf(p+1) = 0.0;

   for i=2:N
    arg   = p*(i-1.0)*kh;
    cs    = cos(arg);
    sn    = sin(arg);
    af(p+1) = af(p+1) + cs*f(i);
    bf(p+1) = bf(p+1) + sn*f(i);
   end

   af(p+1) = ff*af(p+1);
   bf(p+1) = ff*bf(p+1);

end

af(M+1) = fc*af(M+1);

%---
% done
%---

return
