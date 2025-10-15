function [gxx,gxs  ...
         ,gsx,gss] ...
...
  = sgf_ax_1p_ct  ...
...
  (x,s     ...
  ,x0,s0   ...
  ,sc      ...
  ,RL      ...
  ,Nsum,Np ...
  )

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%----------------------------------------

%-----------------------------------------------------
%  Flow due to a periodic array of point-force rings
%  inside a cylinder of radius sc
%
%  The singularity is located at (x0, s0)
%  The Green's function is evaluated at (x, s)
%
%  SYMBOLS:
%  --------
%
%  RL: period along the x axis (cylinder axis)
%  sc: radius of cylinder
%
%  N0, N1, N2: 
%
%  Number of terms retained in the infinite sum
%  of sum in real space
%
%  Computations are performed for 
%
%  N0=Nsum, N1=N0*Np, N2=N1*Np
%
%  and the result is extapolated using Aitken's method
%--------------------------------------

  tol=0.0000001;

%----------
% constants
%----------

  pi2 = 2.0*pi;
  pi4 = 4.0*pi;

%----------------------------------
% if the points are on the axis,
% move them slightly off the axis
%----------------------------------

   if( s<tol) s  = tol; end
   if(s0<tol) s0 = tol; end

%--------
% prepare
%--------

  dx = x-x0;
  wn = pi2/RL; % wave number

%-------------------------------------
% first compute the complementary part
%-------------------------------------

 ssss = s+s0-2.0*sc; 

%---
% set truncation limit
%---

%     M = - int(25.0D0/ssss)
%     M = - int(15.0D0/ssss)
%     M = 20;
      M = 10;

 rmcxx = 0.0D0;
 rmcxs = 0.0D0;
 rmcsx = 0.0D0;
 rmcss = 0.0D0;

for i=0:M

   t = i*wn;

   if(i==0) t = 0.0001; end;

   oc  = sc*t;
   o0  = s0*t;
   o   = s *t;
   osn = -ssss *t;

   BI0c = besseli(0,oc);
   BI1c = besseli(1,oc);
   BK0c = besselk(0,oc);
   BK1c = besselk(1,oc);

   BI00 = besseli(0,o0);
   BI10 = besseli(1,o0);
   BK00 = besselk(0,o0);
   BK10 = besselk(1,o0);

   bxx =  (-2.0D0*BK0c+oc*BK1c)*BI00 - BK0c*o0*BI10;
   bxs =              -oc*BK0c *BI00 + BK1c*o0*BI10;
   bsx =              -oc*BK1c *BI10 + BK0c*o0*BI00;
   bss =  ( 2.0D0*BK1c+oc*BK0c)*BI10 - BK1c*o0*BI00;

   temp = oc*BI1c+2.0D0*BI0c;

   Det  = t*( oc*BI0c*BI0c-BI1c*temp );
   fc   = 4.0D0/Det;

   axx = fc*  ( oc*BI0c*bxx - bxs*temp );
   asx = fc*  ( oc*BI0c*bsx - bss*temp );
   axs = fc*t*(    BI0c*bxs - bxx*BI1c );
   ass = fc*t*(    BI0c*bss - bsx*BI1c );

   BI0 = besseli(0,o);
   BI1 = besseli(1,o);
   BK0 = besselk(0,o);
   BK1 = besselk(1,o);

   fxx = t*BI0*axx + (o*BI1+2.0*BI0)*axs;
   fxs = t*BI0*asx + (o*BI1+2.0*BI0)*ass;
   fsx = t*BI1*axx +  o*BI0         *axs;
   fss = t*BI1*asx +  o*BI0         *ass;

   BI0sn = besseli(0,osn);
   BI1sn = besseli(1,osn);
   BK0sn = besselk(0,osn);
   BK1sn = besselk(1,osn);

   fxx = fxx + 8.0D0*BK0sn;

   dxi = dx*t;
   dc  = cos(dxi);
   ds  = sin(dxi);

   fc = 1.0D0;
   if(i==0) fc = 0.50D0; end

   rmcxx = rmcxx + fxx *dc *fc;
   rmcxs = rmcxs + fxs *ds *fc;
   rmcsx = rmcsx + fsx *ds *fc;
   rmcss = rmcss - fss *dc *fc;

end

ccff  = s0*wn;

rmcxx = ccff*rmcxx;
rmcxs = ccff*rmcxs;
rmcsx = ccff*rmcsx;
rmcss = ccff*rmcss;

%---------------------------------
% sum in real space over the rings
%---------------------------------

      dx0 = dx;

      ss0   = s*s0;
      ss04  = 4.0D0*ss0;
      rss0  = sqrt(ss0);
      s2    = s*s;
      s02   = s0*s0;
      sc2   = sc*sc;
      rs0   = sqrt(s0);
      rs    = sqrt(s);
      sss02 = (s+s0)^2;
      ssd02 = (s-s0)^2;
      ssss2 = ssss^2;

%---
% prepare for fast summation
% with Aitken extrapolation
%---

  N0 = Nsum;
  N1 = N0*Np;
  N2 = N1*Np;

%---
% primary ring
%---

      dx  = dx0;
      dx2 = dx*dx;
      r2  = dx2+ssd02;
      rk2 = ss04/(dx2+sss02);
      rk  = sqrt(rk2);

      [f,e] = ellipke(rk2);

      hh  = sqrt(dx2+ssss2);

      rmrxx = 2.0D0*rk*(f+dx2*e/r2)*rs0/rs  - pi4*s0/hh;
      rmrxs = - rk*dx*(f-(s2-s02+dx2)*e/r2)/rss0;
      rmrsx = + rk*dx*(f+(s2-s02-dx2)*e/r2)*rs0/(s*rs);
      rmrss = + rk*((s02+s2+2.0*dx2)*f ...
              -( 2.0D0*dx2*dx2 +  3.0D0*dx2*(s02+s2) ...
              +  (s2-s02)*(s2-s02) )*e/r2 )/(rs0*s*rs);

%---
% sum over periodic rings
%---

      for Ir=1:N2

        dx = dx0+Ir*RL;    % ring to the left

        dx2 = dx*dx;
        r2  = dx2+ssd02;
        rk2 = ss04/(dx2+sss02);
        rk  = sqrt(rk2);

        [f,e] = ellipke(rk2);

        hh  = sqrt(dx2+ssss2);
        xx  = 2.0*rk*(f+dx2*e/r2)*rs0/rs  - pi4*s0/hh;
        xs  = - rk*dx*(f-(s2-s02+dx2)*e/r2)/rss0;
        sx  = + rk*dx*(f+(s2-s02-dx2)*e/r2)*rs0/(s*rs);
        ss  = + rk*((s02+s2+2.0*dx2)*f ...
                 -( 2.0*dx2*dx2 + 3.0*dx2*(s02+s2) ...
                 +  (s2-s02)*(s2-s02) )*e/r2 )/(rs0*s*rs);

        if(Ir==N0 | Ir==N1 | Ir==N2) % hold for extrapolation
         fxx = xx;
         fxs = xs;
         fsx = sx;
         fss = ss;
        end

        rmrxx = rmrxx + xx;
        rmrxs = rmrxs + xs;
        rmrsx = rmrsx + sx;
        rmrss = rmrss + ss;

        dx  = dx0-Ir*RL;    % ring to the right

        dx2 = dx*dx;
        r2  = dx2+ssd02;
        rk2 = ss04 / (dx2+sss02);
        rk  = sqrt(rk2);

        [f,e] = ellipke(rk2);

        hh = sqrt(dx2+ssss2);
        xx = 2.0D0 * rk * (f+dx2*e/r2)*rs0/rs  - pi4*s0/hh;
        xs = - rk*dx*(f-(s2-s02+dx2)*e/r2)/rss0;
        sx = + rk*dx*(f+(s2-s02-dx2)*e/r2)*rs0/(s*rs);
        ss = + rk*((s02+s2+2.0*dx2)*f ...
                 -( 2.0D0 * dx2*dx2 +  3.0D0 * dx2*(s02+s2) ...
                 +  (s2-s02)*(s2-s02) )*e/r2 )/(rs0*s*rs);

        if(Ir==N0 | Ir==N1 | Ir==N2) % hold for extrapolation
         fxx = fxx + xx;
         fxs = fxs + xs;
         fsx = fsx + sx;
         fss = fss + ss;
        end

        rmrxx = rmrxx + xx;
        rmrxs = rmrxs + xs;
        rmrsx = rmrsx + sx;
        rmrss = rmrss + ss;

      if(Ir==N0) 
        Gxx0 = rmrxx - 0.5D0*fxx;
        Gxs0 = rmrxs - 0.5D0*fxs;
        Gsx0 = rmrsx - 0.5D0*fsx;
        Gss0 = rmrss - 0.5D0*fss;
      elseif(Ir==N1)
        Gxx1 = rmrxx - 0.5D0*fxx;
        Gxs1 = rmrxs - 0.5D0*fxs;
        Gsx1 = rmrsx - 0.5D0*fsx;
        Gss1 = rmrss - 0.5D0*fss;
      elseif(Ir==N2)
        Gxx2 = rmrxx - 0.5D0*fxx;
        Gxs2 = rmrxs - 0.5D0*fxs;
        Gsx2 = rmrsx - 0.5D0*fsx;
        Gss2 = rmrss - 0.5D0*fss;
      end

     end

%---------------------
% Aitken extrapolation
%---------------------

 rmrxx = (Gxx0*Gxx2-Gxx1^2)/(Gxx2-2.0D0*Gxx1+Gxx0);
 rmrxs = (Gxs0*Gxs2-Gxs1^2)/(Gxs2-2.0D0*Gxs1+Gxs0);
 rmrsx = (Gsx0*Gsx2-Gsx1^2)/(Gsx2-2.0D0*Gsx1+Gsx0);
 rmrss = (Gss0*Gss2-Gss1^2)/(Gss2-2.0D0*Gss1+Gss0);

%---------------------------------------------
% add the primary and complementary components
%---------------------------------------------

 gxx = rmrxx + rmcxx;
 gxs = rmrxs + rmcxs;
 gsx = rmrsx + rmcsx;
 gss = rmrss + rmcss;

%-----
% done
%-----

return
