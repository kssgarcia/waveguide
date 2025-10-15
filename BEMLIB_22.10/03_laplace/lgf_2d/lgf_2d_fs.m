 function [G,Gx,Gy] = lgf_2d_fs ...
...
  (x,y ...
  ,x0,y0 ...
  ,Iopt)

%=========================================
% FDLIB, CFDLAB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%=========================================

%-------------------------------------------
% Free-space Green's function:
%
%  G = -(1/2*pi) * ln(r)
%
% Iopt = 1: compute only the Green's function
%     ne 1: compute the Green's function
%           and the gradient
%-------------------------------------------

%----------
% constants
%----------

pi2 = 2*pi;

%-----------
% initialize
%-----------

G = 0;
Gx = 0;
Gy = 0;

%-----------------
% Green's function
%-----------------

dx = x-x0;
dy = y-y0;
rs = dx*dx+dy*dy;
G  = -0.5*log(rs)/pi2;

%---------
% gradient
%---------

if(Iopt>1)

 den = rs*pi2;
 Gx = - dx/den;
 Gy = - dy/den;

end

%-----
% done
%-----

return
