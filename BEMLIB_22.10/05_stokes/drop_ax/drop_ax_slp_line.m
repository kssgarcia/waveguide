function [Qxx,Qxs ...
         ,Qsx,Qss] = pendant_slp_line  ...
...
  (X0,Y0 ...
  ,X1,Y1...
  ,X2,Y2...
  ,NGL ...
  ,Iflow ...
  ,wall ...
  ,sc      ...
  ,RL      ...
  ,Nsum,Np ...
  ,Ising)

%----------------------------------------
% integrate the Green's function over
% a straight segment
%----------------------------------------

Iopt=1;

[ZZ,WW]=gauss_leg(NGL);

%-----------
% initialize
%-----------

Qxx = 0.0D0;
Qxs = 0.0D0;
Qsx = 0.0D0;
Qss = 0.0D0;

%---
% prepare for the quadrature
%---

XM = 0.5D0*(X2+X1);
XD = 0.5D0*(X2-X1);
YM = 0.5D0*(Y2+Y1);
YD = 0.5D0*(Y2-Y1);
DR = sqrt(XD*XD+YD*YD);

%--------------------------
% loop over Gaussian points
%--------------------------

for i=1:NGL

  X = XM+XD*ZZ(i);
  Y = YM+YD*ZZ(i);

   if(Iflow==1)

   [Wxx,Wxs,Wsx,Wss ...
 ...
   ,QXXX,QXXY,QXYX,QXYY ...
   ,QYXX,QYXY,QYYX,QYYY ...
 ...
   ,PXX,PXY,PYX,PYY ...
   ,Iaxis] ...
 ...
   = sgf_ax_fs (Iopt,X0,Y0,X,Y);

   elseif(Iflow==2)

   [Wxx,Wxs,Wsx,Wss ...
  ...
   ,QXXX,QXXY,QXYX,QXYY ...
   ,QYXX,QYXY,QYYX,QYYY ...
  ...
   ,PXX,PXY,PYX,PYY ...
   ,Iaxis] ...
  ...
    = sgf_ax_w (Iopt,X0,Y0,X,Y,wall);

  elseif(Iflow==3)

   [Wxx,Wxs,Wsx,Wss]  ...
...
   = sgf_ax_1p_ct  ...
...
  (X0,Y0     ...
  ,X,Y   ...
  ,sc      ...
  ,RL      ...
  ,Nsum,Np ...
  );

   end

%---
% subtract out the singularity
%---

 if(Ising==1)
  Dist2 = (X-X0)^2+(Y-Y0)^2;
  DD = log(Dist2);
  Wxx = Wxx + DD;
  Wss = Wss + DD;
 end

%---
% carry out the quadrature
%---

WI = WW(i);

Qxx = Qxx + Wxx * WI;
Qxs = Qxs + Wxs * WI;
Qsx = Qsx + Wsx * WI;
Qss = Qss + Wss * WI;

end

%---
% finish up
%---

Qxx = Qxx*DR;
Qxs = Qxs*DR;

Qsx = Qsx*DR;
Qss = Qss*DR;

%---------------------
% add back the singularity
% the point X0 is a segment end-point
%---------------------

if(Ising==1)
 SEGL=2.0D0*DR;
 Qxx = Qxx - 2.0D0* SEGL*(log(SEGL)-1.0D0);
 Qss = Qss - 2.0D0* SEGL*(log(SEGL)-1.0D0);
end

%-----
% done
%-----

return
