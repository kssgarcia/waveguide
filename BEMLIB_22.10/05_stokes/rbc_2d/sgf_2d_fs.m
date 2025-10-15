function [Gxx,Gxy ...
         ,Gyx,Gyy ...
         ,Px,Py ...
         ,Txxx,Txxy,Tyxx,Tyxy ...
         ,Txyx,Txyy,Tyyx,Tyyy ] ...
 ...
      = sgf_2d_fs  ...
       ...
        (x,y,x0,y0,Iopt)
         
%-----------------------------------------
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

Px=0;
Py=0;
Txxx=0;
Txxy=0;
Tyxx=0;
Tyxy=0;
Txyx=0;
Txyy=0;
Tyyx=0;
Tyyy=0;

%--------------------------------------------------
% Green's function for flow bounded by a plane wall
% located at y = wall
%
% see Pozrikidis (1992, p. 93)
%
% Iopt =  1 compute only G
%      ne 1 compute G, P, and T
%--------------------------------------------------

%====================
% primary point force
%====================

dx = x-x0;
dy = y-y0;
dxx = dx*dx;
dxy = dx*dy;
dyy = dy*dy;
r2  = dxx+dyy;
r   = sqrt(r2);
ri2 = 1.0D0/r2;
rl = log(r);

Gxx = -rl + dxx*ri2;
Gxy =       dxy*ri2;
Gyy = -rl + dyy*ri2;
Gyx = Gxy;

%--------------------
% pressure and stress:
%--------------------

if(Iopt~=1)

  cf = -4.0D0/(r2*r2);
  Txxx = dxx*dx * cf;
  Txxy = dxy*dx * cf;
  Tyxx = Txxy;
  Tyxy = dyy*dx * cf;

  Txyx = Txxy;
  Txyy = Tyxy;
  Tyyx = Txyy;
  Tyyy = dyy*dy * cf;

  cf = 2.0D0*ri2;
  Px = dx*cf;  % pressure vector
  Py = dy*cf;

end

%-----
% Done
%-----

return
