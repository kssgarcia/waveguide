function bs = fourier_sin(N,L,f)

%=======================================
% Evaluation of sine Fourier coefficients
%=======================================

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

for p=0:N-1         % loop over coefficients

   bs(p+1) = 0.0;

   for i=2:N
        arg = 0.5*p*(i-1.0)*kh;
        sn  = sin(arg);
        bs(p+1) = bs(p+1) + sn*f(i);
    end

   bs(p+1) = ff*bs(p+1);

end

bs(N) = 0.5*bs(N);

%---
% done
%---

return
