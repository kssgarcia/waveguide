function [x,s,volume] = drop_ax_ode ...
   ...
   (npts ...
   ,capls ...
   ,Isp ...
   ,dpsi ...
   ,shp ...
   )

%-----------------------------------------
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%--------------------------------------------
% Integrate ODEs by the modified Euler method
% with a uniform step size for the angle psi
%
% SYMBOLS:
% -------
%
% dpsi: increments in psi
% shp:  shooting parameter
%--------------------------------------------

%--------
% prepare
%--------

  dpsih = 0.5D0*dpsi;

%----------------
% top of the drop
%----------------

 psi  = 0.0D0;
 x(1) = 0.0D0;
 s(1) = 0.0D0;

 for i=1:npts

   if(i==1) 
     xp = 0.0D0;
     sp = 2.0D0/shp;
   else
     Q  = sin(psi)/s(i)+Isp*x(i)/capls-shp;
     xp = sin(psi)/Q;
     sp =-cos(psi)/Q;
   end
   xp1 = xp;
   sp1 = sp;

   psi    = psi    +dpsih;
   x(i+1) = x(i)+xp*dpsih;
   s(i+1) = s(i)+sp*dpsih;

   Q  = sin(psi)/s(i+1)+Isp*x(i+1)/capls-shp;
   xp = sin(psi)/Q;
   sp =-cos(psi)/Q;
   xp2 = xp;
   sp2 = sp;

   x(i+1) = x(i)+xp*dpsih;
   s(i+1) = s(i)+sp*dpsih;

   Q  = sin(psi)/s(i+1)+Isp*x(i+1)/capls-shp;
   xp = sin(psi)/Q;
   sp =-cos(psi)/Q;
   xp3 = xp;
   sp3 = sp;

   psi    = psi    +dpsih;
   x(i+1) = x(i)+xp*dpsi;
   s(i+1) = s(i)+sp*dpsi;

   Q  = sin(psi)/s(i+1)+Isp*x(i+1)/capls-shp;
   xp = sin(psi)/Q;
   sp =-cos(psi)/Q;
   xp4 = xp;
   sp4 = sp;

   x(i+1) = x(i) + (xp1+2*xp2+2*xp3+xp4)*dpsi/6.0;
   s(i+1) = s(i) + (sp1+2*sp2+2*sp3+sp4)*dpsi/6.0;

   end

%-------------------------------------------
% compute the volume of the integrated shape
% by the trapezoidal rule
%-------------------------------------------

 volume = 0.0D0;

 for i=1:npts
   volume = volume+(s(i+1)^2+s(i)^2)*abs(x(i+1)-x(i));
 end

 volume = 0.5D0*volume; % to account for trapezoidal weights
 volume = pi*volume;

%-----
% done
%-----

 return
