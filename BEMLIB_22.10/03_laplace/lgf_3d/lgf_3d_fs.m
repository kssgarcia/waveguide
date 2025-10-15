 function [G,Gx,Gy,Gz] = lgf_3d_fs ...
...
  (x,y,z ...
  ,x0,y0,z0 ...
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
% Free-space Green's function of Laplace's
% equation:
%
%  G = 1/(4*pi*r)
%
% 
%  where: r = |x-x0|
%
% Iopt = 1: compute only the Green's function
%     ne 1: compute the Green's function
%           and the gradient
%-------------------------------------------

%----------
% constants
%----------

pi4 = 4*pi;

%-----------
% initialize
%-----------

G  = 0;
Gx = 0;
Gy = 0;
Gz = 0;

%-----------------
% Green's function
%-----------------

Dx = x-x0;
Dy = y-y0;
Dz = z-z0;

r = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
G  = 1.0/(pi4*r);

%---------
% gradient
%---------

if(Iopt>1)

 den = pi4*r^3;

 Gx = - Dx/den;
 Gy = - Dy/den;
 Gz = - Dz/den;

end

%-----
% done
%-----

return
