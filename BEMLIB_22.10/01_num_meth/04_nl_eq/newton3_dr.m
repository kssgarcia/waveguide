clear all
close all

%---------
% driver for newton's method
% for 3 equations
%---------

menu  =1;
Niter = 8;
eps = 0.001;

%---
% initial guess
%---

x(1) = 0.5;
x(2) = 2.5;
x(3) = 0.5;

%---
% newton
%---

 [x,f,Iflag] = newton3 ...
     ...
     (menu ...
     ,Niter ...
     ,eps ...
     ,x ...
     )

%---
% done
%---
