function [Dfx,Dfy,elten,tsten] = rbc_df ...
...
    (NSG ...
    ,Xint ...
    ,srtn ...
    ,suns,crvuns ...
    ,s,crv,vnx,vny ...
    ,elst ...
    ,bmod ...
    ,Nsmooth)

%=========================
% compute the jump in traction across the membrane
%
% s:                 arc length
% suns:   unstressed arc length
% crv:               curvature
% crvuns: unstressed curvature
% srtn:   surface tension
%=========================

%---
% prepare
%---

NSG1 = NSG+1;
NSG2 = NSG+2;

%---
% wrap
%---

  suns(NSG2)=  suns(2)+suns(NSG1);
     s(NSG2)=     s(2)+   s(NSG1);
crvuns(NSG2)=crvuns(2);
   crv(NSG2)=   crv(2);
  srtn(NSG2)=  srtn(2);

%----------------------------
% compute elastic tensions
% by cubic spline interpolation
%----------------------------

 Ido=1;
 Ido=0;

%---
 if(Ido==1)
%---

 for i=1:NSG+1
   Yint(i) = suns(i);
 end
 [Aint,Bint,Cint] = splc_pr (NSG,Xint,Yint);
 for i=1:NSG
   dsunsdxint(i) = Cint(i);
 end

 for i=1:NSG+1
   Yint(i) = s(i);
 end
 [Aint,Bint,Cint] = splc_pr (NSG,Xint,Yint);
 for i=1:NSG
   dsdxint(i) = Cint(i);
 end

 for i=1:NSG
   dsdsuns(i) = dsdxint(i)/dsunsdxint(i);
   elten(i) = elst*(dsdsuns(i)-1.0);      % linear constit equation
 end
 dsdsuns(NSG+1)=dsdsuns(1);

 for i=1:NSG+1
   Yint(i) = crvuns(i)-crv(i);
 end
 [Aint,Bint,Cint] = splc_pr (NSG,Xint,Yint);
 for i=1:NSG
   tsten(i) = bmod*Cint(i)/dsdxint(i);
 end

 elten(NSG+1) = elten(1);
 elten(NSG+2) = elten(2); % wrap
 tsten(NSG+1) = tsten(1);
 tsten(NSG+2) = tsten(2); % wrap

% dsdsuns'
%  figure(76)
%  hold on
%  plot(s(1:NSG+1),dsdsuns(1:NSG+1),'kx-')

%---
  end % of Ido
%---

%-----------------------------
% compute the elastic tensions
% by parabolic interpolation
%-----------------------------

 Ido=0;
 Ido=1;

%---
 if(Ido==1)
%---

 crv = smooth(NSG,crv,Nsmooth);
 crv(NSG+2)=crv(2);

 for i=2:NSG+1

   ia = i-1;
   i1 = i+1;
   x0 = suns(ia)-suns(i);
   x1 = suns(i1)-suns(i);
   y0 = s(ia)-s(i);
   y1 = s(i1)-s(i);
   DsDs0(i) = (x0*y1/x1 - x1*y0/x0)/(x0-x1);

   elten(i) = elst*(DsDs0(i)-1.0D0);      % linear constit equation

   % transverse shear tension:

   x0 = s(ia)-s(i);
   x1 = s(i1)-s(i);
   crvia = crv(ia) - crvuns(ia);
   crvi  = crv(i)  - crvuns(i);
   crvi1 = crv(i1) - crvuns(i1);
   y0 = crvia - crvi;
   y1 = crvi1 - crvi;
   tsten(i) = bmod * (x0*y1/x1 - x1*y0/x0)/(x0-x1);

 end

   DsDs0(1) = DsDs0(NSG+1); % wrap

%   DsDs0'
%   figure(76)
%   plot(s(1:NSG+1),DsDs0(1:NSG+1),'r-o')

  elten(1)    = elten(NSG1); % wrap
  elten(NSG2) = elten(2); % wrap

  tsten(1)    = tsten(NSG1); % wrap
  tsten(NSG2) = tsten(2); % wrap

%----
  end % of Ido
%----

%-------------------------------------
%  Compute the derivative d(srtn+elten)/d(s)
%  at the nodes by parabolic interpolation
%
%  Compute Df
%------------------------------------

  elten = smooth(NSG,elten,Nsmooth);
  elten(NSG2) = elten(2);  % wrap

  tsten = smooth(NSG,tsten,Nsmooth);
  tsten(NSG2) = tsten(2);  % wrap

  for i=2:NSG+1

   ia = i-1;
   i1 = i+1;
   x0 = s(ia)-s(i);
   x1 = s(i1)-s(i);
   y0 = srtn(ia)- srtn(i) + elten(ia)-elten(i);
   y1 = srtn(i1)- srtn(i) + elten(i1)-elten(i);
   DtDs = (x0*y1/x1 - x1*y0/x0)/(x0-x1);

   y0 = tsten(ia)-tsten(i);
   y1 = tsten(i1)-tsten(i);
   DtstenDs = (x0*y1/x1 - x1*y0/x0)/(x0-x1);
   tnx(i) =-vny(i);
   tny(i) = vnx(i);

   tentot = srtn(i)+elten(i);

   FN = tentot*crv(i) - DtstenDs;

   Dfx(i) = FN*vnx(i) - tnx(i)*(DtDs + tsten(i) );
   Dfy(i) = FN*vny(i) - tny(i)*(DtDs + tsten(i) );

  end

%---
% wrap
%---

 Dfx(1) = Dfx(NSG+1);
 Dfy(1) = Dfy(NSG+1);

%---
% done
%---

 return
