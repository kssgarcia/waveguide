function [F,E]= ell_int (k)

%-----------------------------------------
% FDLIB BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%--------------------------------------------
%  computation of complete elliptic integrals
%  of the first (F) and second (E) kind
%--------------------------------------------

tol=0.0000000000001;

%-----------
% iterations
%-----------

ks = k*k;

F = 0.5*pi;
E = 1.0;
G = 1.0;
B = k;
D = 10*tol;

%pause
%disp('one')

while abs(D)>tol,
  C = sqrt(1.0-B^2);
  B = (1.0-C)/(1.0+C);
  D = F*B;
  F = F+D;
  G = 0.50*G*B;
  E = E+G;
%pause
%disp('two')
end

E = F*(1.0-0.50*ks*E);

%-----
% done
%-----

return
