function [Gxx,Gxy,Gxz ...
         ,Gyx,Gyy,Gyz ...
         ,Gzx,Gzy,Gzz ...
         ,px,py,pz ...
         ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz ...
         ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz ...
         ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz ...
         ] = ...
...
      sgf_3d_w ...
...
       (Iopt ...
       ,x,y,z ...
       ,x0,y0,z0 ...
       ,wall ...
       )

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%----------------------------------------
% Green's function for semi-infinite flow
% bounded by a plane wall located at x = wall
%
% Ref: Pozrikidis (1992, p. 84)
%
% Iopt =  1 Compute only the Green's function
%      ne 1 Compute the Green's function,
%           pressure, and stress
%----------------------------------------

%-----------
% initialize
%-----------

 px = 0; py = 0; pz = 0;

 Txxx=0; Txxy=0; Txxz=0; Tyxy=0; Tyxz=0; Tzxz=0;
 Txyx=0; Txyy=0; Txyz=0; Tyyy=0; Tyyz=0; Tzyz=0;
 Txzx=0; Txzy=0; Txzz=0; Tyzy=0; Tyzz=0; Tzzz=0;

%-------------
% pimary part
%-------------

  dx = x-x0;
  dy = y-y0;
  dz = z-z0;

  dxx = dx*dx;
  dxy = dx*dy;
  dxz = dx*dz;
  dyy = dy*dy;
  dyz = dy*dz;
  dzz = dz*dz;

  r   = sqrt(dxx+dyy+dzz);
  r3  = r*r*r;
  r5  = r3*r*r;
  ri  = 1.0D0/r;
  ri3 = 1.0D0/r3;
  ri5 = 1.0D0/r5;

  Sxx = ri + dxx*ri3;
  Sxy =      dxy*ri3;
  Sxz =      dxz*ri3;
  Syy = ri + dyy*ri3;
  Syz =      dyz*ri3;
  Szz = ri + dzz*ri3;

%------------------------------
  if(Iopt==2) % compute the stress and pressure

      cf = -6.0D0*ri5;

      TSxxx = dxx*dx * cf;
      TSxxy = dxy*dx * cf;
      TSxxz = dxz*dx * cf;
      TSyxy = dyy*dx * cf;
      TSyxz = dyz*dx * cf;
      TSzxz = dzz*dx * cf;

      TSyyy = dyy*dy * cf;
      TSyyz = dyz*dy * cf;
      TSzyz = dzz*dy * cf;

      TSyzy = dyy*dz * cf;
      TSyzz = dyz*dz * cf;
      TSzzz = dzz*dz * cf;

      cf = 2.0D0*ri3;

      psx = dx * cf;
      psy = dy * cf;
      psz = dz * cf;

  end
      
%------------------------------

%-------------
% image system
%-------------

  x0im = 2.0*wall-x0;
  dx   = x-x0im;

  dxx = dx*dx;
  dxy = dx*dy;
  dxz = dx*dz;
      
  r   = sqrt(dxx+dyy+dzz);
  r3  = r*r*r;
  r5  = r3*r*r;
  ri  = 1.0D0/r;
  ri3 = 1.0D0/r3;
  ri5 = 1.0D0/r5;

%----------------
% image stokeslet 
%----------------

  Sxx = Sxx - ri - dxx*ri3;
  Sxy = Sxy      - dxy*ri3;
  Sxz = Sxz      - dxz*ri3;
  Syy = Syy - ri - dyy*ri3;
  Syz = Syz      - dyz*ri3;
  Szz = Szz - ri - dzz*ri3;

  Syx = Sxy;
  Szx = Sxz;
  Szy = Syz;

%-----------
   if(Iopt==2) % compute the stress and pressure

      cf = -6.0D0*ri5;

      TSxxx = TSxxx - dxx*dx * cf;
      TSxxy = TSxxy - dxy*dx * cf;
      TSxxz = TSxxz - dxz*dx * cf;
      TSyxy = TSyxy - dyy*dx * cf;
      TSyxz = TSyxz - dyz*dx * cf;
      TSzxz = TSzxz - dzz*dx * cf;

      TSxyx = TSxxy;
      TSxyy = TSyxy;
      TSxyz = TSyxz;
      TSyyy = TSyyy - dyy*dy * cf;
      TSyyz = TSyyz - dyz*dy * cf;
      TSzyz = TSzyz - dzz*dy * cf;

      TSxzx = TSxxz;
      TSxzy = TSyxz;
      TSxzz = TSzxz;
      TSyzy = TSyzy - dyy*dz * cf;
      TSyzz = TSyzz - dyz*dz * cf;
      TSzzz = TSzzz - dzz*dz * cf;

      cf = 2.0D0*ri3;

      psx = psx - dx * cf;
      psy = psy - dy * cf;
      psz = psz - dz * cf;

   end;
%-----------

%-----------------------
% image potential dipole
%-----------------------

  PDxx = - ri3 + 3.0D0*dxx*ri5;
  PDyx =         3.0D0*dxy*ri5;
  PDzx =         3.0D0*dxz*ri5;

  PDxy = - PDyx;
  PDyy =   ri3 - 3.0D0*dyy*ri5;
  PDzy =       - 3.0D0*dyz*ri5;

  PDxz = - PDzx;
  PDyz =   PDzy;
  PDzz =   ri3 - 3.0D0*dzz*ri5;

%-----------------------
      if(Iopt==2) % compute the stress and pressure

      r7  = r5*r*r;
      ri7 = 1.0D0/r7;

      cf  =  6.0D0*ri5;
      cf1 = 30.0D0*ri7;

      TPDxxx =  (dx+dx+dx) * cf - dxx*dx * cf1;
      TPDxxy =   dy        * cf - dxy*dx * cf1;
      TPDxxz =   dz        * cf - dxz*dx * cf1;
      TPDyxy =   dx        * cf - dyy*dx * cf1;
      TPDyxz =                  - dyz*dx * cf1;
      TPDzxz =   dx        * cf - dzz*dx * cf1;

      TPDxyx = - dy        * cf + dxx*dy * cf1;
      TPDxyy = - dx        * cf + dxy*dy * cf1;
      TPDxyz =                    dxz*dy * cf1;
      TPDyyy = -(dy+dy+dy) * cf + dyy*dy * cf1;
      TPDyyz = - dz        * cf + dyz*dy * cf1;
      TPDzyz = - dy        * cf + dzz*dy * cf1;

      TPDxzx = - dz        * cf + dxx*dz * cf1;
      TPDxzy =                    dxy*dz * cf1;
      TPDxzz = - dx        * cf + dxz*dz * cf1;
      TPDyzy = - dz        * cf + dyy*dz * cf1;
      TPDyzz = - dy        * cf + dyz*dz * cf1;
      TPDzzz = -(dz+dz+dz) * cf + dzz*dz * cf1;

      % no pressure

      end
%-----------------------

%-----------------------
% image stokeslet dipole
%-----------------------

      SDxx = dx * PDxx;
      SDyx = dx * PDyx - dy*ri3;
      SDzx = dx * PDzx - dz*ri3;

      SDxy = dx * PDxy - dy*ri3;
      SDyy = dx * PDyy;
      SDzy = dx * PDzy;

      SDxz = dx * PDxz - dz*ri3;
      SDyz = dx * PDyz;
      SDzz = dx * PDzz;

%-----------------------
   if(Iopt==2) % compute the stress and pressure

      cf = 6.0D0*ri5;

      TSDxxx = dx *TPDxxx;
      TSDxxy = dx *TPDxxy - (   -dxy) * cf;
      TSDxxz = dx *TPDxxz - (   -dxz) * cf;
      TSDyxy = dx *TPDyxy - (dxx-dyy) * cf;
      TSDyxz = dx *TPDyxz - (   -dyz) * cf;
      TSDzxz = dx *TPDzxz - (dxx-dzz) * cf;

      TSDxyx = dx *TPDxyx + dxy * cf;
      TSDxyy = dx *TPDxyy ;
      TSDxyz = dx *TPDxyz ;
      TSDyyy = dx *TPDyyy + dxy * cf;
      TSDyyz = dx *TPDyyz ;
      TSDzyz = dx *TPDzyz + dxy * cf;

      TSDxzx = dx *TPDxzx + dxz * cf;
      TSDxzy = dx *TPDxzy ;
      TSDxzz = dx *TPDxzz ;
      TSDyzy = dx *TPDyzy + dxz * cf;
      TSDyzz = dx *TPDyzz;
      TSDzzz = dx *TPDzzz + dxz * cf;

      PSDx =  - 2.0*ri3 + 6.0*dxx*ri5;
      PSDy =            - 6.0*dxy*ri5;
      PSDz =            - 6.0*dxz*ri5;

  end
%-----------------------

%---------
% assemble
%---------

      h0   = x0-wall ;
      h02  = 2.0*h0 ;
      h0s2 = 2.0*h0*h0 ;

      Gxx = Sxx + h0s2 * PDxx - h02 * SDxx ;
      Gxy = Sxy + h0s2 * PDxy - h02 * SDxy ;
      Gxz = Sxz + h0s2 * PDxz - h02 * SDxz ;

      Gyx = Syx + h0s2 * PDyx - h02 * SDyx ;
      Gyy = Syy + h0s2 * PDyy - h02 * SDyy ;
      Gyz = Syz + h0s2 * PDyz - h02 * SDyz ;

      Gzx = Szx + h0s2 * PDzx - h02 * SDzx ;
      Gzy = Szy + h0s2 * PDzy - h02 * SDzy ;
      Gzz = Szz + h0s2 * PDzz - h02 * SDzz ;

%-----------------------
      if(Iopt==2) % compute the stress and pressure

      Txxx = TSxxx + h0s2 * TPDxxx - h02 * TSDxxx ;
      Txxy = TSxxy + h0s2 * TPDxxy - h02 * TSDxxy ;
      Txxz = TSxxz + h0s2 * TPDxxz - h02 * TSDxxz ;
      Tyxy = TSyxy + h0s2 * TPDyxy - h02 * TSDyxy ;
      Tyxz = TSyxz + h0s2 * TPDyxz - h02 * TSDyxz ;
      Tzxz = TSzxz + h0s2 * TPDzxz - h02 * TSDzxz ;
      
      Txyx = TSxyx + h0s2 * TPDxyx - h02 * TSDxyx ;
      Txyy = TSxyy + h0s2 * TPDxyy - h02 * TSDxyy ;
      Txyz = TSxyz + h0s2 * TPDxyz - h02 * TSDxyz ;
      Tyyy = TSyyy + h0s2 * TPDyyy - h02 * TSDyyy ;
      Tyyz = TSyyz + h0s2 * TPDyyz - h02 * TSDyyz ;
      Tzyz = TSzyz + h0s2 * TPDzyz - h02 * TSDzyz ;

      Txzx = TSxzx + h0s2 * TPDxzx - h02 * TSDxzx ;
      Txzy = TSxzy + h0s2 * TPDxzy - h02 * TSDxzy ;
      Txzz = TSxzz + h0s2 * TPDxzz - h02 * TSDxzz ;
      Tyzy = TSyzy + h0s2 * TPDyzy - h02 * TSDyzy ;
      Tyzz = TSyzz + h0s2 * TPDyzz - h02 * TSDyzz ;
      Tzzz = TSzzz + h0s2 * TPDzzz - h02 * TSDzzz ;

      px = psx - h02 * PSDx ;
      py = psy - h02 * PSDy ;
      pz = psz - h02 * PSDz ;

      end
%----------------------

%-----
% done
%-----

return
