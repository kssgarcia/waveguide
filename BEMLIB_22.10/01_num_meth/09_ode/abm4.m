%----
% testing the Adams-Bashforth-Moulton method
% against a known exponential solution
%----

h=0.1;

x0=1.0
f1=exp(h);
f0=1.0
fA=exp(-h);
fB=exp(-2*h);
fC=exp(-3*h);

% AB4

x1=x0+h/24*(-9*fC+37*fB-59*fA+55*f0);
er = x1-exp(h);
-er/h^5/(251/720)

% AM4

x1=x0+h/24*(fB-5*fA+19*f0+9*f1);
er = x1-exp(h);
er/h^5/(19/720)
