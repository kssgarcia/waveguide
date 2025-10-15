%=================================
% Matlab code: Gauss--Salamin
%
% Algorithm for rapidly computing
% the decimal expansion of pi
%
% M.G.Blyth 2007
%=================================

format long

%---
% initialise
%---

a = 1.0;
b = 1.0/sqrt(2);
c = 0.25;
x = 1.0;

%---
% iterate
%---

N = 4;     % number of iterations

for i=1:N
  y = a;
  a = 0.5*(a+b);
  b = sqrt(b*y);
  c = c - x*(a-y)^2; 
  xold = x;
  x    = 2.0*x;
  pisal(i) = (a+b)^2/(4*c);
end 

pisal'
