 clear all
 close all

%===========================
% driver for Newton's method
% for two equations
%===========================

 Niter = 8;
 eps = 0.001;

 menu  =1;

%---
% initial guess
%---

 x(1) = 0.1;
 x(2) = 0.2;

%---
% newton' method
%---

 [x,f,Iflag] = newton2 ...
     ...
     (menu ...
     ,Niter ...
     ,eps ...
     ,x ...
     )

%---
% done
%---

