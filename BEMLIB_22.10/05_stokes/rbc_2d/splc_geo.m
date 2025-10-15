function [vnx,vny,crv,s ...
         ,Xint ...
         ,Axint,Bxint,Cxint ...
         ,Ayint,Byint,Cyint ...
         ,area,xcenter,ycenter ...
         ,aspect,angle1,angle2] = splc_geo (N,X,Y);

%--------------------------------
% Compute :
%
%  the normal vector (vnx, vny)
%  the curvature (crv)
%  the arc length (s)
%  the area (area)
%  the center 
%  the aspect ratio
%  the inclination angles
%
%  along a closed line
%  using cubic-spline interpolation
%
%  Normal vector points outward:
%
% curvature is positive for a circle
%--------------------------------

%--------
% parameters
%--------

Nstep=32;

%--------
% prepare
%--------

N1 = N+1;
N2 = N+2;
N3 = N+3;

%----------------------------
% interpolation variable
% is the polygonal arc length
%----------------------------

Xint(1) = 0.0D0;

for i=2:N1
 ia = i-1;
 Xint(i) = Xint(ia)+sqrt((X(i)-X(ia))^2 ...
                        +(Y(i)-Y(ia))^2);
end

%-------------
% interpolate X
%--------------

 for i=1:N1
  Yint(i) = X(i);
 end

[Axint,Bxint,Cxint] = splc_pr (N,Xint,Yint);

%--------------
% interpolate Y
%--------------

 for i=1:N1
   Yint(i) = Y(i);
 end

 [Ayint,Byint,Cyint] = splc_pr (N,Xint,Yint);

%-----------------------
% compute the arc length
% and the area
%-----------------------

   area = 0.0;
   xmom = 0.0;
   ymom = 0.0;
   xxmom = 0.0;
   xymom = 0.0;
   yymom = 0.0;

   s(1) = 0.0D0;
%   figure
%   hold on
%   plot(X,Y,'ro')

%---
 for i=1:N  % run over segments
%---

   XX1 = Xint(i);
   XX2 = Xint(i+1);

   step = (XX2-XX1)/Nstep;

%---
% integrate by the trapezoidal rule
%---

      DX = 0.0D0;  % first point

      xspl = ((Axint(i)*DX + Bxint(i) )*DX + Cxint(i))*DX+X(i);
      yspl = ((Ayint(i)*DX + Byint(i) )*DX + Cyint(i))*DX+Y(i);
%      plot(xspl,yspl)
      Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i);
      Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i);

      tmp = sqrt( Dxp*Dxp + Dyp*Dyp);
      arr   = 0.5*yspl;
      arl   = 0.5*tmp;
      sumx  = 0.5*xspl*tmp;
      sumy  = 0.5*yspl*tmp;
      sumxx = 0.5*xspl*xspl*tmp;
      sumxy = 0.5*xspl*yspl*tmp;
      sumyy = 0.5*yspl*yspl*tmp;

      for j=2:Nstep
        DX = (j-1.0D0)*step;
        xspl = ((Axint(i)*DX + Bxint(i) )*DX + Cxint(i))*DX+X(i);
        yspl = ((Ayint(i)*DX + Byint(i) )*DX + Cyint(i))*DX+Y(i);
%        plot(xspl,yspl)
        Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i);
        Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i);

        tmp = sqrt( Dxp*Dxp + Dyp*Dyp );
        arr   = arr   + yspl;
        arl   = arl   + tmp;
        sumx  = sumx  + xspl*tmp;
        sumy  = sumy  + yspl*tmp;
        sumxx = sumxx + xspl*xspl*tmp;
        sumxy = sumxy + xspl*yspl*tmp;
        sumyy = sumyy + yspl*yspl*tmp;
      end

      DX = XX2-XX1; % last point
      xspl = ((Axint(i)*DX + Bxint(i) )*DX + Cxint(i))*DX+X(i);
      yspl = ((Ayint(i)*DX + Byint(i) )*DX + Cyint(i))*DX+Y(i);
%      plot(xspl,yspl)
      Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i);
      Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i);

      tmp   = sqrt( Dxp*Dxp + Dyp*Dyp );
      arr   = arr   + 0.5*yspl;
      arl   = arl   + 0.5*tmp;
      sumx  = sumx  + 0.5*xspl*tmp;
      sumy  = sumy  + 0.5*yspl*tmp;
      sumxx = sumxx + 0.5*xspl*xspl*tmp;
      sumxy = sumxy + 0.5*xspl*yspl*tmp;
      sumyy = sumyy + 0.5*yspl*yspl*tmp;

      area = area - arr*(X(i+1)-X(i))/Nstep;

      s(i+1) = s(i)  + arl  *step;
      xmom   = xmom  + sumx *step;
      ymom   = ymom  + sumy *step;
      xxmom  = xxmom + sumxx*step;
      xymom  = xymom + sumxy*step;
      yymom  = yymom + sumyy*step;

 end

%---
% center
%----

 xcenter = xmom/s(N+1);
 ycenter = ymom/s(N+1);

 mom(1,1) = xxmom-2.0*xmom*xcenter+xcenter*xcenter*s(N+1);
 mom(1,2) = xymom-xmom*ycenter-ymom*xcenter+xcenter*ycenter*s(N+1);
 mom(2,1) = mom(1,2);
 mom(2,2) = yymom-2.0*ymom*ycenter+ycenter*ycenter*s(N+1);

 [V,D] = eig(mom);

 lambda1 = max(abs(D(1,1)),abs(D(2,2)));
 lambda2 = min(abs(D(1,1)),abs(D(2,2)));

 aspect = lambda1/lambda2;

 angle1 = atan(V(2,1)/V(1,1));
 angle2 = atan(V(2,2)/V(1,2));

% if(angle1<0)
%  angle1=angle1+pi;
% end
% if(angle2<0)
%  angle2=angle2+pi;
% end
  if(angle2>angle1)
   save=angle1;
   angle1=angle2;
   angle2=save;
  end

%-----------------------------------
% compute the outward normal vector
% and the curvature
% at the nodes
%-----------------------------------

  for i=1:N
    Den = sqrt(Cxint(i)^2+Cyint(i)^2);
    vnx(i) =  Cyint(i)/Den;
    vny(i) = -Cxint(i)/Den;
    crv(i) = -2.0D0*(Bxint(i)*Cyint(i) ...
                   -Byint(i)*Cxint(i))/Den^3;
  end

%------------
% wrap around
%------------

 s(N2) = s(N1)+s(2);
 s(N3) = s(N1)+s(3);

 vnx(N1) = vnx(1);
 vny(N1) = vny(1);
 crv(N1) = crv(1);

 Xint(N2) = Xint(N1)+Xint(2);
 vnx(N2)  = vnx(2);
 vny(N2)  = vny(2);
 crv(N2)  = crv(2);

 Xint(N3) = Xint(N1)+Xint(3);
 vnx(N3)  = vnx(3);
 vny(N3)  = vny(3);
 crv(N3)  = crv(3);

%-----
% done
%-----

return
