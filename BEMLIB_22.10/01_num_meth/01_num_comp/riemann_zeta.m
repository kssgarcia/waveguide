function ryemann = riemann_zeta(s)

%=======================================
% zeta = sum_{i=1}^infty 1/i^s
%
% Sum a series  that decays like 1/i^s
% using the Aitken extrapolation method
%
% Exact values for integer s are given by
% Abramowitz & Stegun (1972), p. 811.
%=======================================

format long

%---
% parameters
%---

kmax=15;
N=1;
p=2;

%---
% exact
%---

if(s==2)
 exact = pi^2/6.0;
elseif(s==3)
 exact = 1.20205690315959428540;
elseif(s==4)
 exact = pi^4/90.0;
else
 exact = 0; % not available
end

%---
% sum
%---

L=0;
M=N;

sum = 0.0;

for k=1:kmax+2
  for i=L+1:M
    sum = sum+1.0/(i^s);
   end
   x(k) = sum-0.5/(M^s);
   if(k>2)
    a(k-1) = (x(k-2)*x(k)-x(k-1)*x(k-1))/(x(k)-2.0*x(k-1)+x(k-2));
   end
   if(k>3)
   if(abs(a(k-1)-a(k-2))<0.00000001)
    q = k-1;
    ryemann = a(k-1);
    break;
   end
   end
   L = M;
   M = M*p;
end

%===
% verbose
%===

ispeak=1;
ispeak=0;

if(ispeak==1)

 dsp(1,1)=x(1);
 dsp(1,2)=x(1)-exact;

 for k=2:q
   dsp(k,1)=x(k);
   dsp(k,2)=x(k)-exact;
   dsp(k,3)=a(k);
   dsp(k,4)=a(k)-exact;
 end
 dsp(q+1,1)=x(q+1);
 dsp(q+1,2)=x(q+1)-exact;
 dsp

end
%===

return

