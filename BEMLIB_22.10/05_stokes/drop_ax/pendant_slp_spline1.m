function [Qx,Qy ...
         ,Wx,Wy] = pendant_slp_spline1  ...
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
  ,wall ...
  ,Ising ...
  )

%----------------------------------------
% Integrate the Green's function over
% a spline element
% Desingularize by subtracting off the
% log term
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

%-------------------
% probe the singularity
%------------------

  if(Ising>0)

  if(Ising==1)
   Dfn0 = Dfn1;
   vnx0 = vnx1;
   vny0 = vny1;
   DX=0;
  end

  if(Ising==2)
   Dfn0 = Dfn2;
   vnx0 = vnx2;
   vny0 = vny2;
   DX=Xint2-Xint1;
  end

  Dxp = (3.0D0*Axint*DX + 2.0D0*Bxint )*DX + Cxint;
  Dyp = (3.0D0*Ayint*DX + 2.0D0*Byint )*DX + Cyint;
  h0 = sqrt(Dxp*Dxp+Dyp*Dyp);

%---
  end
%---

%--------------------------
% loop over Gaussian points
%--------------------------

for i=1:NGL

  Xint = XintM + XintD*ZZ(i);
  DX = Xint-Xint1;
  X = ((Axint*DX + Bxint )*DX + Cxint)*DX+X1;
  Y = ((Ayint*DX + Byint )*DX + Cyint)*DX+Y1;

  [Wxx,Wxs,Wsx,Wss ...
  ,QXXX,QXXY,QXYX,QXYY ...
  ,QYXX,QYXY,QYYX,QYYY ...
  ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
 ...
   = sgf_ax_w (Iopt,X,Y,X0,Y0,wall);

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

  Dfn  =  DfnM + DfnD*ZZ(i);

  Wx = Wx + (Wxx*vnx+Wxs*vny)*Dfn*WI;
  Wy = Wy + (Wsx*vnx+Wss*vny)*Dfn*WI;

  if(Ising>0)
   WI0 = h0*WW(i)*XintD;
   if(Ising==1)
    ldst = 2.0*log(abs(Xint-Xint1));
   end
   if(Ising==2)
    ldst = 2.0*log(abs(Xint-Xint2));
   end
   Qx = Qx + ldst*vnx0*WI0;
   Qy = Qy + ldst*vny0*WI0;
   Wx = Wx + ldst*vnx0*Dfn0*WI0;
   Wy = Wy + ldst*vny0*Dfn0*WI0;
  end

%---
end
%---

  if(Ising>0)
   DDD = abs(Xint2-Xint1);
   SING = 2.0*DDD*(log(DDD)-1.0);
   Qx = Qx - h0*vnx0*SING;
   Qy = Qy - h0*vny0*SING;
   Wx = Wx - h0*vnx0*Dfn0*SING;
   Wy = Wy - h0*vny0*Dfn0*SING;
  end

%-----
% done
%-----

return
