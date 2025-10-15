%===================================
% matlab code: Leibniz
%
% Algorithm for computing pi
% using the famous Leibniz series
%
% M.G.Blyth 2007
%===================================

format long

%----
% Initialise
%----

x = 0.0;
i = 1;

%----
% iterate
%----

N = 10^2;  % number of iterations

for m=1:N
  j = (i-1)/2;
  x = x + (-1)^j/i;
  i = i + 2;
end;
    
pii = 4.0*x;
[pii, pi]

