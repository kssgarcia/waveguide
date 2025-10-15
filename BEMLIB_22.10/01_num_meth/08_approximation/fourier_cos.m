function ac=fourier_cos(N,L,f)

%---
% prepare
%---

k = 2*pi/L;
h = L/N;
kh = 2*pi/N;
ff = 2.0/N;

%--------
% launch
%--------

for p=0:N         % loop over coefficients

   ac(p+1) = 0.5D0*(f(1)+f(N+1)*cos(p*pi));

   for i=2:N
        arg   = 0.5*p*(i-1.0)*kh;
        cs    = cos(arg);
        sn    = sin(arg);
        ac(p+1) = ac(p+1) + cs*f(i);
    end

   ac(p+1) = ff*ac(p+1);

end

ac(N+1) = 0.5*ac(N+1);

return
