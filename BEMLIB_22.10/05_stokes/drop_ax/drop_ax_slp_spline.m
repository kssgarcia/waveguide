function [Qx,Qy ...
         ,Wx,Wy] = pendant_slp_spline  ...
...
  (X0,Y0 ...
  ,X1,Y1...
  ,X2,Y2...
  ,NGL ...
  ,Xint1,Xint2 ...
  ,Axint,Bxint,Cxint ...
  ,Ayint,Byint,Cyint ...
  ,Dfn1,Dfn2 ...
  ,vnx1,vnx2 ...
  ,vny1,vny2 ...
  ,Dfn0 ...
  ,Iflow ...
  ,wall ...
  ,sc      ...
  ,RL      ...
  ,Nsum,Np ...
  ,Ising ...
  )

%----------------------------------------
% Integrate the Green's function over
% a spline element
%
% Only for constant surface tension
%----------------------------------------

Iopt=1;

[ZZ,WW]=gauss_leg(NGL);

%-----------
% initialize
%-----------

Qx = 0.0D0;
Qy = 0.0D0;

Wx = 0.0D0;
Wy = 0.0D0;

%---
% prepare for the quadrature
%---

XintM = 0.5D0*(Xint2+Xint1);
XintD = 0.5D0*(Xint2-Xint1);

DfnM = 0.5D0*(Dfn2+Dfn1);
DfnD = 0.5D0*(Dfn2-Dfn1);

vnxM = 0.5D0*(vnx2+vnx1);
vnxD = 0.5D0*(vnx2-vnx1);

vnyM = 0.5D0*(vny2+vny1);
vnyD = 0.5D0*(vny2-vny1);

%--------------------------
% loop over Gaussian points
%--------------------------

for i=1:NGL

  Xint = XintM + XintD*ZZ(i);
  DX = Xint-Xint1;
  X = ((Axint*DX + Bxint )*DX + Cxint)*DX+X1;
  Y = ((Ayint*DX + Byint )*DX + Cyint)*DX+Y1;

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
   ,QXXX,QXXY,QXYX,QXYY ...
   ,QYXX,QYXY,QYYX,QYYY ...
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
% carry out the quadrature
%---

  Dxp = (3.0D0*Axint*DX + 2.0D0*Bxint )*DX + Cxint;
  Dyp = (3.0D0*Ayint*DX + 2.0D0*Byint )*DX + Cyint;
  h = sqrt(Dxp*Dxp+Dyp*Dyp);

  WI = h*WW(i)*XintD;

  vnx  =  vnxM +  vnxD*ZZ(i);
  vny  =  vnyM +  vnyD*ZZ(i);
  norm = sqrt(vnx*vnx+vny*vny);
  vnx=vnx/norm;
  vny=vny/norm;

  Qx = Qx + (Wxx*vnx+Wxs*vny)*WI;
  Qy = Qy + (Wsx*vnx+Wss*vny)*WI;

  Dfn  =  DfnM +  DfnD*ZZ(i);

  if(Ising==1)
    Dfn=Dfn-Dfn0;
  end

  Wx = Wx + (Wxx*vnx+Wxs*vny)*Dfn*WI;
  Wy = Wy + (Wsx*vnx+Wss*vny)*Dfn*WI;

%---
end
%---

%-----
% done
%-----

return
