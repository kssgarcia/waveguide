function [SLPx, SLPy ...
         ,DLPxx1, DLPxy1 ...
         ,DLPyx1, DLPyy1 ...
         ,DLPxx2, DLPxy2 ...
         ,DLPyx2, DLPyy2 ...
...
] = rbc_sdlp_spline  ... 
...
  (X0,Y0 ...
  ,X1,Y1...
  ,X2,Y2...
  ,Xint1,Xint2 ...
  ,Axint,Bxint,Cxint ...
  ,Ayint,Byint,Cyint ...
  ,Dfx0,Dfy0 ...
  ,Dfx1,Dfx2 ...
  ,Dfy1,Dfy2 ...
  ,vnx1,vnx2 ...
  ,vny1,vny2 ...
  ,Ising)

%----------------------------------------
% Integrate the Green's function over
% a spline segment
% to compute the single- and double-layer
% potential
%----------------------------------------

%---
% prepare
%---

global NGL
global Iflow
global wall
global period

Iopt=2; % for the Green's function

%---
% Gauss quadrature
%---

[ZZ,WW] = gauss_leg(NGL);

%---
% singular element
%---

  if(Ising>0)

   if(Ising==1)     % singular point at the first point
    DX = 0.0;
   elseif(Ising==2) % singular point at the last point
    DX = Xint2-Xint1;
   end

   Dxp = (3.0*Axint*DX + 2.0*Bxint )*DX + Cxint;
   Dyp = (3.0*Ayint*DX + 2.0*Byint )*DX + Cyint;
   h0 = sqrt(Dxp*Dxp+Dyp*Dyp);

  end

%-----------
% initialize
%-----------

SLPx = 0.0;
SLPy = 0.0;

DLPxx1 = 0.0;
DLPxy1 = 0.0;
DLPyx1 = 0.0;
DLPyy1 = 0.0;

DLPxx2 = 0.0;
DLPxy2 = 0.0;
DLPyx2 = 0.0;
DLPyy2 = 0.0;

%---
% prepare for the quadrature
%---

XintM = 0.5D0*(Xint2+Xint1);
XintD = 0.5D0*(Xint2-Xint1);

DfxM = 0.5D0*(Dfx2+Dfx1);
DfxD = 0.5D0*(Dfx2-Dfx1);

DfyM = 0.5D0*(Dfy2+Dfy1);
DfyD = 0.5D0*(Dfy2-Dfy1);

%---
% normal vector
%---

vnxM = 0.5D0*(vnx2+vnx1);
vnxD = 0.5D0*(vnx2-vnx1);

vnyM = 0.5D0*(vny2+vny1);
vnyD = 0.5D0*(vny2-vny1);

%--------------------------
% loop over Gaussian points
%--------------------------

for i=1:NGL

  Xint = XintM + XintD*ZZ(i);

  Dfx = DfxM + DfxD*ZZ(i);
  Dfy = DfyM + DfyD*ZZ(i);

  vnx = vnxM + vnxD*ZZ(i);
  vny = vnyM + vnyD*ZZ(i);

  DX = Xint-Xint1;

  X = ((Axint*DX + Bxint )*DX + Cxint)*DX+X1;
  Y = ((Ayint*DX + Byint )*DX + Cyint)*DX+Y1;
  Dxp = (3.0*Axint*DX + 2.0*Bxint )*DX + Cxint;
  Dyp = (3.0*Ayint*DX + 2.0*Byint )*DX + Cyint;

  h = sqrt(Dxp*Dxp+Dyp*Dyp);

%---
  if(Iflow==0 | Iflow==1)
%---

  [Gxx,Gxy ...
  ,Gyx,Gyy ...
  ,Px,Py ...
  ,Txxx,Txxy,Tyxx,Tyxy ...
  ,Txyx,Txyy,Tyyx,Tyyy ] ...
 ...
   = sgf_2d_fs  ...
   ...
   (X,Y,X0,Y0,Iopt ...
   );

%---
   elseif(Iflow==2)
%---

   [Gxx,Gxy ...
   ,Gyx,Gyy ...
   ,Px,Py ...
   ,Txxx,Txxy,Tyxx,Tyxy ...
   ,Txyx,Txyy,Tyyx,Tyyy ] ...
 ...
    = sgf_2d_w  ...
     ...
   (X,Y,X0,Y0,wall,Iopt ...
   );

%---
   elseif(Iflow==3)
%---

    [Gxx,Gxy ...
    ,Gyx,Gyy ...
    ,Px,Py ...
    ,Txxx,Txxy,Tyxx,Tyxy ...
    ,Txyx,Txyy,Tyyx,Tyyy ] ...
 ...
     = sgf_2d_1p_w  ...
     ...
     (X,Y ...
     ,X0,Y0 ...
     ,wall ...
     ,period...
     ,Iopt ...
     );

%---
   end
%---

%---
% carry out the quadrature
%---

  WI = h*WW(i)*XintD;

  SLPx = SLPx + (Gxx*Dfx+Gyx*Dfy)*WI;
  SLPy = SLPy + (Gxy*Dfx+Gyy*Dfy)*WI;

  psi1 = 0.5*(1-ZZ(i)); % tent-1 function
  DLPxx1 = DLPxx1 + psi1*(Txxx*vnx + Txxy*vny)*WI;
  DLPxy1 = DLPxy1 + psi1*(Txyx*vnx + Txyy*vny)*WI;
  DLPyx1 = DLPyx1 + psi1*(Tyxx*vnx + Tyxy*vny)*WI;
  DLPyy1 = DLPyy1 + psi1*(Tyyx*vnx + Tyyy*vny)*WI;

  psi2 = 0.5*(1+ZZ(i)); % tent-2 function
  DLPxx2 = DLPxx2 + psi2*(Txxx*vnx + Txxy*vny)*WI;
  DLPxy2 = DLPxy2 + psi2*(Txyx*vnx + Txyy*vny)*WI;
  DLPyx2 = DLPyx2 + psi2*(Tyxx*vnx + Tyxy*vny)*WI;
  DLPyy2 = DLPyy2 + psi2*(Tyyx*vnx + Tyyy*vny)*WI;

%---
% singularity of the slp
%---

  if(Ising > 0)

    if(Ising == 1) DelX = Xint-Xint1; end
    if(Ising == 2) DelX = Xint2-Xint; end
    WI = log(DelX)*h0*WW(i)*XintD;
    SLPx = SLPx + Dfx0*WI;
    SLPy = SLPy + Dfy0*WI;

  end
%---

%---
end % over segments
%---

  if(Ising > 0)
    DelX = Xint2-Xint1;
    tmp = h0*DelX*(log(DelX)-1.0);
    SLPx = SLPx - tmp*Dfx0;
    SLPy = SLPy - tmp*Dfy0;
  end

%-----
% done
%-----

return
