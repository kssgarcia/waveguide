x=[0.0 1.0 2.1 3.4];
f=[0.2 1.3 2.4 7.5];

%---
% Aitken interpolation
%---

N=3;

xint = 3.4

[yint, A] = aitken(N,x,f,xint)

