clear all
close all

%---
% graphs of the Stoer and Burlish functions
%---

x1=0.0;
x2=1.0;
h=x2-x1;

psi1_1=-4.0/h;
psi1_2=12.0/h^2;

nplt=64;
step=h/nplt;

for i=1:nplt+1
 x(i) = (i-1.0)*step;
 xx = x(i);
 q13(i) = (xx-x2)^4/h^4 * (xx-x1)^2/2.0;
 q12(i) = (xx-x2)^4/h^4 * (xx-x1 + 4/h * (xx-x1)^2/2.0 );
% q11(i) = (xx-x2)^4/h^4 * (1 + 4/h*(xx-x1) - 12/h^2 * (xx-x1)^2/2.0 );
 q11(i) = (xx-x2)^4/h^4 * (1 -psi1_1*(xx-x1)  ...
                              - (psi1_2-2.0*psi1_1^2) * (xx-x1)^2/2.0 );
end

hold on;
plot(x,q13)
plot(x,q12,'--')
plot(x,q11,':')

xplt(1)=x1;
yplt(1)=0.0;
xplt(2)=x2;
yplt(2)=x2;
plot(xplt,yplt)
