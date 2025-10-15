function [vnx,vny,crv,s ...
         ,Xint ...
         ,Axint,Bxint,Cxint ...
         ,Ayint,Byint,Cyint ...
         ,volume] = splc_geo (N,X,Y,alpha,Jsp);

%--------------------------------
% Compute :
%
%  the normal vector (vnx, vny)
%  the curvature (crv)
%  the arc length (s)
%
%  along a line
%  using cubic-spline interpolation
%
%
%  Normal vector points downward:
%
%            crv>0
%            ____
%   \ crv<0 /    \ crv<0 /
%    \_____/      \_____/
%       ||
%       ||
%       \/ 
%
% curvature is positive when the
% surface is downward parabolic
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

for i=2:N+1
 ia = i-1;
 Xint(i) = Xint(ia)+sqrt((X(i)-X(ia))^2 ...
                        +(Y(i)-Y(ia))^2);
end

%-------------
% interpolate X
%--------------

 for i=1:N+1
  Yint(i) = X(i);
 end

%slope2=-sin(alpha);
[Axint,Bxint,Cxint] = splc_clm (N,Xint,Yint,0.0,0.0);
%[Axint,Bxint,Cxint] = splc_clm_true (N,Xint,Yint,0.0);

%--------------
% interpolate Y
%--------------

 for i=1:N+1
   Yint(i) = Y(i);
 end

%slope2=cos(alpha);
[Ayint,Byint,Cyint] = splc_clm (N,Xint,Yint,1.0,-1.0);
%[Ayint,Byint,Cyint] = splc_clm (N,Xint,Yint,1.0,slope2);
%[Ayint,Byint,Cyint] = splc_nt_true (N,Xint,Yint);

%-----------------------
% compute the arc length
% and the volume
%-----------------------

   volume = 0.0;

   s(1) = 0.0D0;
%   figure
%   hold on
%   plot(X,Y,'ro')

    for i=1:N

      XX1 = Xint(i);
      XX2 = Xint(i+1);
      step = (XX2-XX1)/Nstep;

%---
% integrate by the trapezoidal rule
%---

      DX = 0.0D0;  % first point

      xspl = ((Axint(i)*DX + Bxint(i) )*DX + Cxint(i))*DX+X(i);
      yspl = ((Ayint(i)*DX + Byint(i) )*DX + Cyint(i))*DX+Y(i);
%       plot(yspl,Jsp*xspl,'r')
      Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i);
      Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i);

      vlm = 0.5D0*yspl*yspl;
      sum = 0.50D0*sqrt( Dxp^2 + Dyp^2 );

      for j=2:Nstep
        DX = (j-1.0D0)*step;
        xspl = ((Axint(i)*DX + Bxint(i) )*DX + Cxint(i))*DX+X(i);
        yspl = ((Ayint(i)*DX + Byint(i) )*DX + Cyint(i))*DX+Y(i);
%         plot(yspl,Jsp*xspl,'r')
        Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i);
        Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i);
        vlm = vlm+yspl*yspl;
        sum = sum + sqrt( Dxp^2 + Dyp^2 );
      end

      DX = XX2-XX1; % last point
      xspl = ((Axint(i)*DX + Bxint(i) )*DX + Cxint(i))*DX+X(i);
      yspl = ((Ayint(i)*DX + Byint(i) )*DX + Cyint(i))*DX+Y(i);
%       plot(yspl,Jsp*xspl,'r')
      Dxp = (3.0D0*Axint(i)*DX + 2.0D0*Bxint(i) )*DX + Cxint(i);
      Dyp = (3.0D0*Ayint(i)*DX + 2.0D0*Byint(i) )*DX + Cyint(i);
      vlm = vlm+0.5D0*yspl*yspl;
      sum = sum + 0.50D0*sqrt( Dxp^2 + Dyp^2 );

      volume = volume+vlm*(X(i+1)-X(i))/Nstep;

      s(i+1)=s(i)+sum*step;

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
    crv(i) = 2.0D0*(Bxint(i)*Cyint(i) ...
            -Byint(i)*Cxint(i))/Den^3;
  end

%  vnx(N1)= cos(alpha);
%  vny(N1)= sin(alpha);
  fc= (s(N1)-s(N-1))/(s(N)-s(N-1));
  vnx(N1)= vnx(N-1)+(vnx(N)-vnx(N-1))*fc;
  vny(N1)= vny(N-1)+(vny(N)-vny(N-1))*fc;
  crv(N1)= crv(N-1)+(crv(N)-crv(N-1))*fc;
  vnx(N1) = -1.0;
  vny(N1) =  0.0;

  crv(N+2) = crv(N);
  vnx(N+2) = vnx(N);
  vny(N+2) =-vny(N);
  s(N+2) = 2*s(N+1)-s(N);

%---
% reconstruct the shape
%---

  Ido=0;

  if(Ido==1)

  figure(88)
  hold on

  Nplt = 4;
  axis equal
  for i=1:N
   DXint = Xint(i+1)-Xint(i);
   for j=1:Nplt+1
    xxx =(j-1.0)*DXint/Nplt;
    xplt = X(i)+((Axint(i)*xxx)+Bxint(i)*xxx+Cxint(i))*xxx;
    yplt = Y(i)+((Ayint(i)*xxx)+Byint(i)*xxx+Cyint(i))*xxx;
    plot(xplt,yplt,'r.');
   end
  end
  plot(X,Y,'ob')

  end
 
%-----
% done
%-----

return
