%     program sgf_3d_2p_w_dr

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%-------------------------------------
% Driver for the doubly periodic Green's function
% of 3D Stokes flow
% in a semi-infinite domain bounded by
% a plane wall located at x = wall
%-------------------------------------
 
%----------
% constants
%----------

    pi4 = 4.0D0 *pi;
    pi8 = 8.0D0 *pi;

%---------
% settings
%---------

    Iopt = 1; % will need velocity, pressure, and stress

    wall = 0.0D0;

    x0 = 1.750D0;
    y0 = 0.00D0;
    z0 = 0.00D0;

    a11 = 0.5; a12 = 0.0;
    a21 = 0.0; a22 = 0.8;

    Ns = 16;
    Np = 4;

%---------------
% one evaluation
%---------------

     x = 5.1D0;
     y = 0.1D0;
     z = 0.1D0;

     [Gxx,Gxy,Gxz ...
     ,Gyx,Gyy,Gyz ...
     ,Gzx,Gzy,Gzz ...
     ,px,py,pz ...
     ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz ...
     ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz ...
     ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz ...
     ] = ...
...
      sgf_3d_2p_w ...
...
       (Iopt ...
       ,x,y,z ...
       ,x0,y0,z0 ...
       ,wall ...
       ,a11,a12,a21,a22 ...
       ,Ns,Np ...
       );

      Tyxx = Txxy;
      Tzxx = Txxz;
      Tzxy = Tyxz;
      Tyyx = Txyy;
      Tzyx = Txyz;
      Tzyy = Tyyz;
      Tyzx = Txzy;
      Tzzx = Txzz;
      Tzzy = Tyzz;

%
% pring
%

  [Gxx,Gxy,Gxz;
   Gyx,Gyy,Gyz;
   Gzx,Gzy,Gzz]
