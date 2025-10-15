function [N,X,Y ...
         ,vnx,vny ...
         ,crv,s,Xint ...
         ,Axint,Bxint,Cxint ...
         ,Ayint,Byint,Cyint ...
         ,area,centerx,centery ...
         ,P1,P2 ...
         ,Istop ...
         ,Irepeat ...
         ,Iaction] ...
...
   = prd_2d_open (N ...
                 ,X,Y ...
                 ,P1,P2 ...
                 ,Ich1,thmax ...
                 ,Ich2,spmax ...
                 ,Ich3,spmin ...
                 ,Isym ...
                 ,Italk)

%===============================================
%  FDLIB, BEMLIB
%
% Point redistribution on an OPEN planar line
% in the xy plane
%
% o-o-o-o-o-o-o-o
% 1           N N+1
%
% First point is labeled 1 
% Last point is labeled N+1 
%
% Check for MAXIMUM angle,
%           MAXIMUM point separation,
%           MINIMUM point separation
%
% N: number of interval defined by N+1 points
%    point 1 and N+1 are end-points
%
% Will ask for a redo if Irepeat=1
%
% P1: first interpolated property
% P2: second interpolated property
%
% Isym=1 for a symmetric line
% Isym=0 otherwise
%
% dependencies: splc_geo
%===============================================

%----------
% constants
%----------

   pih = 0.5*pi;

%--------
% prepare
%--------

  Irepeat = 0;
  Istop = 0;
  Iaction = 0;

%-----
% trap
%-----

  if(N>=512) 
    disp('prd_2d_open: too many points'); disp(N);
    Istop = 1;
    Iaction = 5;
    return
  end
 
%--------
% prepare
%--------

   Na = N-1;
   N1 = N+1;
   N2 = N+2;

%---------------------------
% cubic-spline interpolation
%---------------------------

  Itype = 2;     % natural splines

  [vnx,vny ...
  ,crv  ...  % curvature
  ,s ...
  ,Xint ...
  ,Axint,Bxint,Cxint ...
  ,Ayint,Byint,Cyint ...
  ,area,centerx,centery ...
  ,aspect,angle1,angle2] = splc_geo(N,X,Y,Itype);

%-------
% prepare to check
%-------

  loop1 = 1;
  loop2 = N;

  if(Isym==1)
   loop2 = N/2;
  end

%========================
for i=loop1:loop2   % cat
%========================

%----
% prepare
%----

   la = i-2;
   l  = i-1;
   k  = i+1;
   k1 = i+2;

    Ax = Axint(i);
    Bx = Bxint(i);
    Cx = Cxint(i);

    Ay = Ayint(i);
    By = Byint(i);
    Cy = Cyint(i);

   if(i<=N)
   sp = s(k)-s(i);
   end

  if(i==1)
      thtot = 2.0*s(2)*crv(1);
  elseif(i==N1)
      thtot = 2.0*(s(N1)-s(N))*crv(N1);
   else
      thtot = (s(k)-s(l))*crv(i);
  end

  Dth = abs(thtot);

%==========
if(Ich1==1) 
%==========

%-----------------------------------
%  Check for angles subtended by the arcs
%  for i=2, ..., N
%
%  If an arc is too large 
%  remove the middle point and add
%  two evenly-spaced points
%-----------------------------------

 if(Dth>pih)
   disp ('prd_2d: an arc is excessive'); disp(i);
   Istop = 1;
   Iaction = 4;
   return
 end

 %---
 if(Dth>thmax & i>=2 & i<=N)
 %---

 if(Italk==1)
   disp('prd_2d: an arc is too large'); disp(i);
 end

  Ax1 = Axint(l); Bx1 = Bxint(l); Cx1 = Cxint(l);
  Ay1 = Ayint(l); By1 = Byint(l); Cy1 = Cyint(l);

  Xinterp1 = (Xint(l)+2.0*Xint(i))/3.0;

  DX1 = Xinterp1-Xint(l);

  Xnew1 = ( ( Ax1*DX1 + Bx1) * DX1 + Cx1)*DX1 + X(l);
  Ynew1 = ( ( Ay1*DX1 + By1) * DX1 + Cy1)*DX1 + Y(l);

  if(i==2)
      [P1new1,P2new1] = interpolate ...
         ...
         (2*Xint(1)-Xint(2),Xint(l),Xint(i),Xint(k) ...
         ,Xinterp1 ...
         ,2*P1(1)-P1(2),P1(1),P1(2),P1(3) ...
         ,2*P2(1)-P2(2),P2(1),P2(2),P2(3) ...
         );
  else
      [P1new1,P2new1] = interpolate ...
         ...
         (Xint(la),Xint(l),Xint(i),Xint(k) ...
         ,Xinterp1 ...
         ,P1(la),P1(l),P1(i),P1(k) ...
         ,P2(la),P2(l),P2(i),P2(k) ...
         );
  end

   Ax2 = Axint(i); Bx2 = Bxint(i); Cx2 = Cxint(i);
   Ay2 = Ayint(i); By2 = Byint(i); Cy2 = Cyint(i);

   Xinterp2 = (2.0*Xint(i)+Xint(k))/3.0;

   DX2 = Xinterp2-Xint(i);

   Xnew2 = ( ( Ax2*DX2 + Bx2) * DX2 + Cx2 )*DX2 + X(i);
   Ynew2 = ( ( Ay2*DX2 + By2) * DX2 + Cy2 )*DX2 + Y(i);

   if(i==N)
       [P1new2,P2new2] = interpolate ...
    ...
          (Xint(N),Xint(N),Xint(N+1),2*Xint(N1)-Xint(N) ...
          ,Xinterp2 ...
          ,P1(Na),P1(N),P1(N1),2*P1(N1)-P1(N) ...
          ,P2(Na),P2(N),P2(N1),2*P2(N1)-P2(N) ...
          );
   else
       [P1new2,P2new2] = interpolate ...
    ...
          (Xint(l),Xint(i),Xint(k),Xint(k1) ...
          ,Xinterp2 ...
          ,P1(l),P1(i),P1(k),P1(k1) ...
          ,P2(l),P2(i),P2(k),P2(k1) ...
          );
   end

   %---
   % accommodate the new points
   %---

   for j=N1:-1:k
       j1 = j+1;
        X(j1) =  X(j);
        Y(j1) =  Y(j);
       P1(j1) = P1(j);
       P2(j1) = P2(j);
   end

       X(k) = Xnew2;
       Y(k) = Ynew2;
      P1(k) = P1new2;
      P2(k) = P2new2;

       X(i) = Xnew1;
       Y(i) = Ynew1;
      P1(i) = P1new1;
      P2(i) = P2new1;

       N = N1;

   if(Italk==1)
    disp ('one point inserted; segments :'); disp(N);
   end
   Irepeat = 1;
   Iaction = 1;
   return

 %---
   end
 %---

%==========
end
%==========

%=============
if(Ich2==1)
%=============

%-----------------------------------
%  check for maximum point separation
%
%  if a segment is too long,
%  add a point in the middle
%-----------------------------------

 %---
  if(sp>spmax)
 %---

  if(Italk==1)
   disp('prd_2d: a segment is too long'); disp(i);
  end

  Ax = Axint(i); Bx = Bxint(i); Cx = Cxint(i);
  Ay = Ayint(i); By = Byint(i); Cy = Cyint(i);

  Xinterp = 0.50*(Xint(i)+Xint(k));

  DX = Xinterp-Xint(i);

  Xnew = ( ( Ax*DX + Bx) * DX + Cx ) *DX + X(i);
  Ynew = ( ( Ay*DX + By) * DX + Cy ) *DX + Y(i);

  %---
  if(i==1)
  %---

   [P1new,P2new] = interpolate ...
...
    (2*Xint(1)-Xint(2),Xint(1),Xint(2),Xint(3) ...
    ,Xinterp ...
    ,2*P1(1)-P1(2),P1(1),P1(2),P1(3) ...
    ,2*P2(1)-P2(2),P2(1),P2(2),P2(3) ...
    );

  %---
  elseif(i==N)
  %---

   [P1new,P2new] = interpolate ...
...
    (Xint(Na),Xint(N),Xint(N1),2*Xint(N1)-Xint(N) ...
    ,Xinterp ...
    ,P1(Na),P1(N),P1(N1),2*P1(N1)-P1(N) ...
    ,P2(Na),P2(N),P2(N1),2*P2(N1)-P2(N) ...
    );

  %---
  else
  %---

   [P1new,P2new] = interpolate ...
...
    (Xint(l),Xint(i),Xint(k),Xint(k1) ...
    ,Xinterp ...
    ,P1(l),P1(i),P1(k),P1(k1) ...
    ,P2(l),P2(i),P2(k),P2(k1) ...
    );

  %---
  end
  %---

   %---
   % accommodate the new point
   %---

      for j=N1:-1:k
        j1     = j+1;
         X(j1) =  X(j);
         Y(j1) =  Y(j);
        P1(j1) = P1(j);
        P2(j1) = P2(j);
      end

       X(k) =  Xnew;
       Y(k) =  Ynew;
      P1(k) = P1new;
      P2(k) = P2new;

      N = N1;

      if(Italk==1)
        disp (' one point inserted: segments: '); disp(N);
      end
      Irepeat = 1;
      Iaction = 2;
      return

 %---
  end
 %---

%=============
end
%=============

%-----------------------------------
% check for minimum point separation
%-----------------------------------

%=============
if(Ich3==1)   
%=============

  if(sp<spmin & (i==1|i==2))  % remove the second point
      for j=2:N
          j1 = j+1;
        X(j) =  X(j1);
        Y(j) =  Y(j1);
       P1(j) = P1(j1);
       P2(j) = P2(j1);
      end
      N=N-1;
      if(Italk==1)
        disp('prd_2d: second point removed: segments: '); disp(N);
      end
   return
  end

  if(sp<spmin & i==N)  % remove the Nth point
        X(N) =  X(N1);
        Y(N) =  Y(N1);
       P1(N) = P1(N1);
       P2(N) = P2(N1);
       N=N-1;
      if(Italk==1)
       disp('prd_2d: Nth point removed: segments: '); disp(N);
      end
   return
  end

  if(sp<spmin & i>=3 & i<=N-1)   

%-------------------------------------
%  points i and i+1 will be removed
%  and replaced by a mid-point
%
%  But only if the resulting angles
%  and point separations do not exceed
%  the specified maxima
%-------------------------------------

    stemp   = 0.50*(  s(i)+  s(k));
    crvtemp = 0.50*(crv(i)+crv(k));

    sep1 = stemp - s(i-1);
    sep2 = s(i+2)- stemp;

       totang0 = crv(i-1)*(stemp-s(i-2));
       totang1 = crvtemp *(s(i+2)-s(i-1));
       totang2 = crv(i+2)*(s(i+3)-stemp);

       totang0 = abs(totang0);
       totang1 = abs(totang1);
       totang2 = abs(totang2);

%------------------------

   if(    totang0 < thmax ...
            & totang1 < thmax ...
            & totang2 < thmax ...
            & sep1 < spmax ...
            & sep2 < spmax) ...

      if(Italk==1)
         disp('prd_2d: two points are too close'); disp(i);
      end

%      figure
%      hold on
%      plot(X(1:N1),Y(1:N1),'ro')
%      plot(X(1),Y(1),'k+')
%      pause

      Xinterp = 0.50D0*(Xint(i)+Xint(k));
      DX = Xinterp-Xint(i);

      Xnew = ( ( Ax*DX + Bx) * DX + Cx )*DX + X(i);
      Ynew = ( ( Ay*DX + By) * DX + Cy )*DX + Y(i);

      if(i==1)
        [P1new,P2new] = interpolate ...
          ...
          (Xint(N)-Xint(N1),Xint(i),Xint(k),Xint(k1) ...
          ,Xinterp ...
          ,P1(N)-P1(N1),P1(i),P1(k),P1(k1) ...
          ,P2(N),P2(i),P2(k),P2(k1) ...
          );
       else
          [P1new,P2new] = interpolate ...
           ...
           (Xint(l),Xint(i),Xint(k),Xint(k1) ...
           ,Xinterp ...
          ,P1(l),P1(i),P1(k),P1(k1) ...
         ,P2(l),P2(i),P2(k),P2(k1) ...
         );
       end

      for j=k:N
          j1 = j+1;
        X(j) =  X(j1);
        Y(j) =  Y(j1);
       P1(j) = P1(j1);
       P2(j) = P2(j1);
      end

       X(i) = Xnew;
       Y(i) = Ynew;
      P1(i) = P1new;
      P2(i) = P2new;
	
      N  = N-1;
      N1 = N+1;

      if(Italk==1)
        disp('prd_2d: one point removed: segments: '); disp(N);
      end

      Irepeat = 1;
      Iaction = 3;
      return
   end
%------------------------

  end
end

%====================
     end % of cat
%====================

% for i=1:N1
%   Ni = N+i;
%   xmore(i) = X(i); xmore(Ni) = X(i); 
%   ymore(i) = Y(i); ymore(Ni) = Y(i); 
%   P1more(i) = P1(i); P1more(Ni) = P1(i)+P1save; 
%   P2more(i) = P2(i); P2more(Ni) = P2(i);
% end

%shift = floor(N/4)
%for i=1:N1
%   X(i)  = xmore(i+shift);
%   Y(i)  = ymore(i+shift);
%   P1(i) = P1more(i+shift)-P1save;
%   P2(i) = P2more(i+shift);
%end

%-----
% done
%-----

return
