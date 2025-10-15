function [F,E] = ell_int(K)

%===========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%===========================================

%--------------------------------------------
% Complete elliptic integrals of the first
% and second kind
%--------------------------------------------

tol = 0.0000000000001; % tolerance

%--------
% prepare
%--------

F = 0.5*pi;
P = 1.0;
G = 1.0;
B = K;
D = 10.0*tol;

%-----------
% iterations
%-----------

while (abs(D)>tol)
  C = sqrt(1.0-B*B);
  B = (1.0-C)/(1.0+C);
  D = F*B;
  F = F+D;
  G = 0.50*G*B;
  P = P+G;
end

E = F*(1.0-0.5*K*K*P);

%-----
% done
%-----

return
