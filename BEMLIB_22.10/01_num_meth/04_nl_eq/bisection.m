%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: bisection
%
% To find the roots of a nonlinear
% equation using the bisection method
%
% Mark Blyth 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Choose initial interval [a,b]

a = 2.0;  
b = 3.0;

eps = 1.0;
x1  = a;  
x2  = b;

[f1] = bisection_fun(x1);
[f2] = bisection_fun(x2);

figure(1)
hold on

Ic = 0;

while eps > 0.001

 x3 = (x1+x2)/2.;
[f3] = bisection_fun(x3);

 if f3*f1<0.
  x2 = x3;
  f2 = f3;
 else
  x1 = x3;
  f1 = f3;
 end
 eps = abs(f3);

 Ic=Ic+1;
 plot(Ic,eps);

end

%---
% Confirm the final answer
%---

error = bisection_fun(x1)

%---
% display final answer
%---

x1

%---
% Show that location of the root
%---

x=[x1-0.5:0.01:x1+0.5];
plot(x,bisection_fun(x),'r.')
xlabel('x')
ylabel('f')
