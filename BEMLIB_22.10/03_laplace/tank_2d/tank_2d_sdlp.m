function [SL, DL] = tank_2d_sdlp  ...
...
        (x0,y0 ...
        ,ip ...
        ,ie ...
        ,x1,y1 ...
        ,x2,y2 ...
        )

%=======================================
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%=======================================

%-------------------------------------------------------
% Compute the single-layer and double-layer potential
% of the Neumann function over straight segments
%
% LEGEND:
% ------
%
%  ip: node label where the potentials are computed
%  ie: element label
%
%  vnxm, vnym:  normal vector at element mid-point
%  elml: element length
%-------------------------------------------------------

 global vnxm vnym elml 
 global wall1 wall2 wall3
 global NGL ZZ WW

%----------
% constants
%----------

 pi2 = 2.0D0*pi;

%--------
% prepare
%--------

 Iopt = 2;  % for the Green's function

%------------------
% singular element?
%------------------

      Ising = 0;
      if(ip==ie) Ising = 1; end

%-------------------------
% Gauss-Legende quadrature
%-------------------------

      xm = 0.5*(x2+x1);
      xd = 0.5*(x2-x1);

      ym = 0.5*(y2+y1);
      yd = 0.5*(y2-y1);

      SL  = 0.0;
      DLx = 0.0;
      DLy = 0.0;

      %---
      for i=1:NGL
      %---

        x = xm + ZZ(i)*xd;
        y = ym + ZZ(i)*yd;

       [G,Gx,Gy] = lgf_2d_www ...
        ...
         (Iopt ...
         ,x,y ...
         ,x0,y0 ...
         ,wall1,wall2,wall3 ...
         );

%------------------------
% singular element: subtract out the singularity
%
% the free-space Green's function is: 
%
%  G = -1.0/(2*pi) * lnr
%------------------------

        if(Ising==1)
          dx = x-x0;
          dy = y-y0;
          r2 = dx*dx+dy*dy;
          G = G + 0.50D0*log(r2)/pi2;
        end 

        SL  = SL  + G *WW(i);
        DLx = DLx + Gx*WW(i);
        DLy = DLy + Gy*WW(i);

      %---
      end % over NGL
      %---

  elmlh = 0.5*elml(ie);

  SL = elmlh*SL;
  DL = elmlh*( DLx*vnxm(ie) + DLy*vnym(ie) );

%-------------------------------
% if the element is singular,
% add back singular contribution
% to the single layer
% (the contribution from the dlp is zero)
%-------------------------------

  if(Ising==1)
      SL = SL - elml(ie)*(log(elmlh)-1.0)/pi2;
  end

%-----
% done
%-----

  return
