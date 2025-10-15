function [Qxx,Qxy  ...
         ,Qyx,Qyy] ...
...
   = prtcl_ax_slp ...
...
   (Iflow ...
   ,X0,Y0,T0 ...
   ,X1,Y1,T1 ...
   ,X2,Y2,T2 ...
   ,NGL ...
   ,Ising ...
   ,Itype ...
   ,xcntr,ycntr ...
   ,rad ...
   ,amaj,amin,tilt ...
   ,wall ...
   ,sc      ...
   ,RL      ...
   ,Nsum,Np ...
   )

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%----------------------------------------

%-----------------------------------------------------
% Computation of the single-layer potential at
% a collocation point located at the midpoint of
% a straight element or a native segment of an ellipse
%
% SYMBOLS:
% -------
%
% Qij:	components of the slp
%-----------------------------------------------------

%-----------
% quadrature
%-----------

[ZZ,WW] = gauss_leg(NGL);

%--------
% prepare
%--------

 Iopt = 1;

%-----------
% initialize
%-----------

  Qxx = 0.0;
  Qxy = 0.0;
  Qyx = 0.0;
  Qyy = 0.0;

%---
% prepare for the quadrture
%---

 if(Itype==1) % straight segments

   XM = 0.5*(X2+X1);
   XD = 0.5*(X2-X1);
   YM = 0.5*(Y2+Y1);
   YD = 0.5*(Y2-Y1);
   DR = sqrt(XD*XD+YD*YD);
   vnx =  YD/DR;    % unit normal vector
   vny = -XD/DR;    % points into the flow
   fcc = DR;        % factor for numerical integration

 elseif(Itype==2)  % circular elements

   TM = 0.5D0*(T2+T1);
   TD = 0.5D0*(T2-T1);
   ornt = 1.0;
   if(TD<0) ornt = -1.0; end
   DR   = rad*abs(TD);

 elseif(Itype==3)  % elliptical elements

   TM = 0.5D0*(T2+T1);
   TD = 0.5D0*(T2-T1);
   ornt = 1.0;
   if(TD<0) ornt = -1.0; end

   css  = cos(T0);     % for desingularization
   snn  = sin(T0);
   tmpx = amin*css;
   tmpy = amaj*snn;
   alm0 = sqrt(tmpx*tmpx+tmpy*tmpy); % arc length metric
   DR   = alm0*abs(TD);

   cs = cos(tilt);
   sn = sin(tilt);

  end

%---
% loop over Gaussian points
%---

  for i=1:NGL

   if(Itype==1)

     X = XM + XD*ZZ(i);
     Y = YM + YD*ZZ(i);

   elseif(Itype==2)

     T    = TM + TD*ZZ(i);
     css  = cos(T);
     snn  = sin(T);
     X    = xcntr + rad*css;
     Y    = ycntr + rad*snn;
     vnx  = css*ornt;
     vny  = snn*ornt;
     fcc  = rad*TD; 

   elseif(Itype==3)

     T    = TM + TD*ZZ(i);
     css  = cos(T);
     snn  = sin(T);
     tmpx = amaj*css;
     tmpy = amin*snn;
     X    = xcntr + tmpx*cs-tmpy*sn;
     Y    = ycntr + tmpx*sn+tmpy*cs;
     tmpx = amin*css;
     tmpy = amaj*snn;
     alm  = sqrt(tmpx*tmpx+tmpy*tmpy);   % arc length metric
     tmpx = tmpx/alm;
     tmpy = tmpy/alm;
     vnx  = tmpx*cs-tmpy*sn;
     vny  = tmpx*sn+tmpy*cs;
     vnx  = vnx * ornt;   % unit normal vector
     vny  = vny * ornt;   % points into the flow
     fcc  = alm*TD;       % factor for numerical integration

   end

%------------------
    if(Iflow==1)   % flow in free space
%-----------------

   [Gxx,Gxy ...
   ,Gyx,Gyy ...
   ,TXXX,TXXY,TXYX,TXYY ...
   ,TYXX,TYXY,TYYX,TYYY ...
   ,PXX,PXY,PYX,PYY ...
   ,Iaxis] ...
   ...
   = sgf_ax_fs (Iopt,X,Y,X0,Y0);

%------------------
    elseif(Iflow==2)   % flow bounded by a wall at x=wall
%-----------------

   [Gxx,Gxy ...
   ,Gyx,Gyy ...
   ,TXXX,TXXY,TXYX,TXYY ...
   ,TYXX,TYXY,TYYX,TYYY ...
   ,PXX,PXY,PYX,PYY ...
   ,Iaxis] ...
 ...
  = sgf_ax_w (Iopt,X,Y,X0,Y0,wall);

%------------------
    elseif(Iflow==3)   % periodic flow inside a tube
%-----------------

    [Gxx,Gxy  ...
    ,Gyx,Gyy] ...
 ...
 = sgf_ax_1p_ct  ...
...
  (X0,Y0    ...
  ,X,Y   ...
  ,sc      ...
  ,RL      ...
  ,Nsum,Np ...
  );

%------------
  end
%------------

%-------
% subtract off the singularity
%-------

   if(Ising==1) 

     if(Itype==1) 
       Dists = (X-X0)^2+(Y-Y0)^2;
       DD    = 0.5D0*log(Dists);
     elseif(Itype==2)
       Dists = ( rad*(T0-T) )^2;
       DD    = 0.5D0*log(Dists);
     elseif(Itype==3)
       Dists = ( alm0*(T0-T) )^2;
       DD    = 0.5D0*alm0/alm * log(Dists);
     end

     Gxx = Gxx + 2.0D0*DD;
     Gyy = Gyy + 2.0D0*DD;

  end

%----------
% end of subtract off the singularity
%----------

      WI = WW(i)*fcc;

      Qxx = Qxx + Gxx * WI;
      Qxy = Qxy + Gxy * WI;
      Qyx = Qyx + Gyx * WI;
      Qyy = Qyy + Gyy * WI;

%---
end % of Gaussian quadrature
%---

%--------------------------------
% add singularity back to the slp 
%--------------------------------

  if(Ising==1) 

    deduction = 4.0D0*DR*(log(DR)-1.0D0);
    Qxx = Qxx - deduction;
    Qyy = Qyy - deduction;

  end

%-----
% done
%-----

 return
