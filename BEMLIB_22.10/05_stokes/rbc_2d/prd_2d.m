function [N,X,Y ...
         ,vnx,vny ...
         ,crv,s,Xint ...
         ,Axint,Bxint,Cxint ...
         ,Ayint,Byint,Cyint ...
         ,area,centerx,centery ...
         ,P1,P2 ...
         ,Istop,Irepeat,action] ...
...
   = prd_2d (N,X,Y,P1,P2 ...
            ,Ich1,thmax,Ich2,spmax,Ich3,spmin,Isym)

%------------------------------------------------
%  Point redistribution around a closed 2D line
%
%  Checks for MAXIMUM ANGLE, 
%             MAXIMUM POINT SEPARATION,
%             MINIMUM POINT SEPARATION,
%
% N: number of points
%
% will return to redo if Irepeat=1
%
% P1: arc length
% P2: another property
%------------------------------------------

%----------
% constants
%----------

   pih = 0.5*pi;

   P1save = P1(N+1)-P1(1);  % to be used for wrapping

%---
% wrap
%---

   P1(N+2) = P1save+P1(2);
   P2(N+2) =        P2(2);

%--------
% prepare
%--------

  Istop = 0;
  Irepeat = 0;
  action = 0;

   if(N>=512) 
      disp (' prd_2d: too many points'); disp(N);
      Istop = 1;
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

  [vnx,vny,crv,s ...
  ,Xint ...
  ,Axint,Bxint,Cxint ...
  ,Ayint,Byint,Cyint ...
  ,area,centerx,centery ...
  ,aspect,angle1,angle2] = splc_geo (N,X,Y);

   Xintsave = Xint(N+1)-Xint(1);   % to be used for wrapping

%-------
% Checks
%-------

loop1 = 1;
loop2 = N;

if(Isym==1)
 loop2=N/2;
end

%========================
for I=loop1:loop2   % cat
%========================

%----
% prepare
%----

   La = I-2;
   L  = I-1;
   K  = I+1;
   K1 = I+2;

   Ax = Axint(I);
   Bx = Bxint(I);
   Cx = Cxint(I);

   Ay = Ayint(I);
   By = Byint(I);
   Cy = Cyint(I);

   sp = s(K)-s(I);

   if(I==1)
      s0=s(N)-s(N+1);
      thtot = (s(2)-s0)*crv(1);
      Dth = abs(thtot);
    else
      thtot = (s(K)-s(L))*crv(I);
      Dth = abs(thtot);
   end

%-----------------------------------
%  check for angles subtended by the arcs
%
%  If I=1 add another point
%
%  If an arc is too large for I>1,
%  remove the middle point and add
%  two evenly-spaced points
%-----------------------------------


if(Ich1==1) 

 if(Dth>pih)
   disp (' prd_2d: an arc is excessive'); disp(I);
   Istop = 1;
   return
 end

 if(Dth>thmax) 

 disp (' prd_2d: an arc is too large'); disp(I);

  %~~~~~~
  if(I==1)
  %~~~~~~

      Ax = Axint(N); Bx = Bxint(N); Cx = Cxint(N);
      Ay = Ayint(N); By = Byint(N); Cy = Cyint(N);
      XINT1 = (Xint(N)+2.0*Xint(N1))/3.0;
      DX = XINT1-Xint(N);
      Xnew1 = ( ( Ax*DX + Bx) * DX + Cx)*DX + X(N);
      Ynew1 = ( ( Ay*DX + By) * DX + Cy)*DX + Y(N);

      [P1_new1,P2_new1] = interpolate ...
     ...
      (Xint(N-1)-Xintsave,Xint(N)-Xintsave,Xint(1),Xint(2) ...
      ,XINT1-Xintsave ...
      ,P1(N-1)-P1save,P1(N)-P1save,P1(1),P1(2) ...
      ,P2(N-1),P2(N),P2(1),P2(2) ...
      );

%     [Xint(N-1)-Xint(N1),Xint(N)-Xint(N1),Xint(1),Xint(2)]
%     [P2(N-1),P2(N),P2(1),P2(2)]
%     [P1(N-1)-P1(N1),P1(N)-P1(N1),P1(1),P1(2)]

      Ax = Axint(1);Bx = Bxint(1);Cx = Cxint(1);
      Ay = Ayint(1);By = Byint(1);Cy = Cyint(1);
      XINT2 = (2.0*Xint(1)+Xint(2))/3.0;
      DX = XINT2-Xint(1);
      Xnew2 = ( ( Ax*DX + Bx) * DX + Cx)*DX + X(1);
      Ynew2 = ( ( Ay*DX + By) * DX + Cy)*DX + Y(1);

      [P1_new2,P2_new2] = interpolate ...
      ...
      (Xint(N)-Xintsave,Xint(I),Xint(K),Xint(K1) ...
      ,XINT2 ...
      ,P1(N)-P1save,P1(1),P1(2),P1(3) ...
      ,P2(N),P2(1),P2(2),P2(3) ...
      );

      for j=N1:-1:2
       j1 = j+1;
        X(j1) =  X(j);
        Y(j1) =  Y(j);
        P1(j1) = P1(j);
        P2(j1) = P2(j);
      end

       X(K) = Xnew2;
       Y(K) = Ynew2;
       P1(K) = P1_new2;
       P2(K) = P2_new2;

       X(I) = Xnew1;
       Y(I) = Ynew1;
       P1(I) = P1_new1;
       P2(I) = P2_new1;

       X(N+2) = X(1);
       Y(N+2) = Y(1);
       P1(N+2) = P1(1)+P1save;
       P2(N+2) = P2(1);

  %~~~~~~
  else
  %~~~~~~

      Ax1 = Axint(L);Bx1 = Bxint(L);Cx1 = Cxint(L);
      Ay1 = Ayint(L);By1 = Byint(L);Cy1 = Cyint(L);
      XINT1 = (Xint(L)+2.0*Xint(I))/3.0;
      DX1 = XINT1-Xint(L);
      Xnew1 = ( ( Ax1*DX1 + Bx1) * DX1 + Cx1)*DX1 + X(L);
      Ynew1 = ( ( Ay1*DX1 + By1) * DX1 + Cy1)*DX1 + Y(L);

      if(I==2)
      [P1_new1,P2_new1] = interpolate ...
...
         (Xint(N)-Xintsave,Xint(L),Xint(I),Xint(K) ...
         ,XINT1 ...
         ,P1(N)-P1save,P1(L),P1(I),P1(K) ...
         ,P2(N),P2(L),P2(I),P2(K) ...
         );
      else
      [P1_new1,P2_new1] = interpolate ...
...
         (Xint(La),Xint(L),Xint(I),Xint(K) ...
         ,XINT1 ...
         ,P1(La),P1(L),P1(I),P1(K) ...
         ,P2(La),P2(L),P2(I),P2(K) ...
         );
      end

      Ax2 = Axint(I);Bx2 = Bxint(I);Cx2 = Cxint(I);
      Ay2 = Ayint(I);By2 = Byint(I);Cy2 = Cyint(I);
      XINT2 = (2.0*Xint(I)+Xint(K))/3.0;
      DX2 = XINT2-Xint(I);
      Xnew2 = ( ( Ax2*DX2 + Bx2) * DX2 + Cx2 )*DX2 + X(I);
      Ynew2 = ( ( Ay2*DX2 + By2) * DX2 + Cy2 )*DX2 + Y(I);

       [P1_new2,P2_new2] = interpolate ...
...
          (Xint(L),Xint(I),Xint(K),Xint(K1) ...
          ,XINT2 ...
          ,P1(L),P1(I),P1(K),P1(K1) ...
          ,P2(L),P2(I),P2(K),P2(K1) ...
          );

      for j=N1:-1:K
       j1 = j+1;
        X(j1) =  X(j);
        Y(j1) =  Y(j);
       P1(j1) = P1(j);
       P2(j1) = P2(j);
      end

       X(K) = Xnew2;
       Y(K) = Ynew2;
      P1(K) = P1_new2;
      P2(K) = P2_new2;

       X(I) = Xnew1;
       Y(I) = Ynew1;
      P1(I) = P1_new1;
      P2(I) = P2_new1;

  %~~~~~~
   end
  %~~~~~~

      N = N1;
      disp (' One point inserted; segments :'); disp(N);
      Irepeat=1;
      action=1;
      return

   end

end

%-----------------------------------
%  Check for maximum point separation
%
%  if a segment is too long,
%  add a point in the middle
%-----------------------------------

if(Ich2==1)

if(sp>spmax)

  disp (' prd_2d: segment too long'); disp(I);

  Ax = Axint(I);Bx = Bxint(I);Cx = Cxint(I);
  Ay = Ayint(I);By = Byint(I);Cy = Cyint(I);

  XINT = 0.50*(Xint(I)+Xint(K));

  DX = XINT-Xint(I);

  Xnew = ( ( Ax*DX + Bx)*DX + Cx )*DX + X(I);
  Ynew = ( ( Ay*DX + By)*DX + Cy )*DX + Y(I);

  if(I==1)

   [P1_new,P2_new] = interpolate ...
...
    (Xint(N)-Xintsave,Xint(1),Xint(2),Xint(3) ...
    ,XINT ...
    ,P1(N)-P1save,P1(1),P1(2),P1(3) ...
    ,P2(N),P2(1),P2(2),P2(3) ...
    );

  elseif(I==N)

   [P1_new,P2_new] = interpolate ...
...
    (Xint(L),Xint(I),Xint(K),Xint(N1)+Xint(2) ...
    ,XINT ...
    ,P1(L),P1(I),P1(K),P1(2)+P1(N1) ...
    ,P2(L),P2(I),P2(K),P2(2) ...
    );

  else

   [P1_new,P2_new] = interpolate ...
...
    (Xint(L),Xint(I),Xint(K),Xint(K1) ...
    ,XINT ...
    ,P1(L),P1(I),P1(K),P1(K1) ...
    ,P2(L),P2(I),P2(K),P2(K1) ...
    );

  end

      for j=N1:-1:K
        j1    = j+1;
         X(j1) =  X(j);
         Y(j1) =  Y(j);
        P1(j1) = P1(j);
        P2(j1) = P2(j);
      end

       X(K) =  Xnew;
       Y(K) =  Ynew;
      P1(K) = P1_new;
      P2(K) = P2_new;

      N = N1;
      N1 = N+1;

      disp (' one point inserted: segments: '); disp (N);

      action=2;
      Irepeat=1;
      return

end
end

%-----------------------------------
% Check for minimum point separation
%-----------------------------------

if(Ich3==1)   
if(sp<spmin)   

%-------------------------------------
%  POINTS I and I+1 WILL BE REMOVED
%  AND REPLACED BY A MIDPOINT
%
%  BUT ONLY IF THE RESULTING ANGLES
%  AND POINT SEPARATIONS DO NOT EXCEED
%  THE PRE-ESTABLISHED MAXIMA
%-------------------------------------

    stemp   = 0.50D0*(  s(I)+  s(K));
    crvtemp = 0.50D0*(crv(I)+crv(K));

    if(I==1)
     SEP1 = stemp-s(1)+s(N1)-s(N);
    else
     SEP1 = stemp - s(I-1);
    end

    if(I==N)
     SEP2 = s(2)+s(N1)- stemp;
    else
     SEP2 = s(I+2)- stemp;
    end

       if(I==1)
         TOTANG0 = crv(N)*(stemp+s(N1)-s(N-1));
       elseif(I==2)
         TOTANG0 = crv(I-1)*(stemp+s(N1)-s(N));
       else
         TOTANG0 = crv(I-1)*(stemp-s(I-2));
       end

       if(I==1)
         TOTANG1 = crvtemp *(s(3)+s(N1)-s(N));
       elseif(I==N)
         TOTANG1 = crvtemp *(s(2)+s(N1)-s(I-1));
       else
         TOTANG1 = crvtemp *(s(I+2)-s(I-1));
       end

       if(I==N)
         TOTANG2 = crv(2)*(s(I)+s(N1)-stemp);
       else
         TOTANG2 = crv(I+2)*(s(I+3)-stemp);
       end

       TOTANG0 = abs(TOTANG0);
       TOTANG1 = abs(TOTANG1);
       TOTANG2 = abs(TOTANG2);

%------------------------

   if(    TOTANG0 < thmax ...
            & TOTANG1 < thmax ...
            & TOTANG2 < thmax ...
            & SEP1 < spmax ...
            & SEP2 < spmax) ...

      disp (' prd_2d: two points are too close'); disp(I);

%      figure
%      hold on
%      plot(X(1:N1),Y(1:N1),'ro')
%      plot(X(1),Y(1),'k+')
%      pause

      XINT = 0.50D0*(Xint(I)+Xint(K));
      DX = XINT-Xint(I);

      Xnew = ( ( Ax*DX + Bx) * DX + Cx )*DX + X(I);
      Ynew = ( ( Ay*DX + By) * DX + Cy )*DX + Y(I);

   if(I==1)
   [P1_new,P2_new] = interpolate ...
   ...
    (Xint(N)-Xint(N1),Xint(I),Xint(K),Xint(K1) ...
    ,XINT ...
    ,P1(N)-P1(N1),P1(I),P1(K),P1(K1) ...
    ,P2(N),P2(I),P2(K),P2(K1) ...
    );
   else
   [P1_new,P2_new] = interpolate ...
   ...
    (Xint(L),Xint(I),Xint(K),Xint(K1) ...
    ,XINT ...
    ,P1(L),P1(I),P1(K),P1(K1) ...
    ,P2(L),P2(I),P2(K),P2(K1) ...
    );
   end

      for j=K:N
          j1 = j+1;
        X(j) =  X(j1);
        Y(j) =  Y(j1);
       P1(j) = P1(j1);
       P2(j) = P2(j1);
      end

       X(I) = Xnew;
       Y(I) = Ynew;
      P1(I) = P1_new;
      P2(I) = P2_new;
	
      N  = N-1;
      N1 = N+1;

      X(N1) = X(1);
      Y(N1) = Y(1);
      P1(N1) = P1(1)+P1save;
      P2(N1) = P2(1);

      disp (' one point removed: segments: '); disp (N);

      action=3;
      Irepeat = 1;
      return
   end
%------------------------

end
end

%====================
     end % of cat
%====================

%-----
% done
%-----

return
