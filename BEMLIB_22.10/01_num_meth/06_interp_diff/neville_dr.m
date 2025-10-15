x=[0.0 1.0 2.1 3.4];
f=[0.2 1.3 2.4 7.5];
N=3;

%----
% Neville interpolation
%----

xint = 3.4

[yint, A] = neville(N,x,f,xint);

