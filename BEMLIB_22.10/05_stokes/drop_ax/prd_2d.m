function [N,X,Y,vnx,vny,crv,s ...
         ,Xint ...
         ,Axint,Bxint,Cxint ...
         ,Ayint,Byint,Cyint ...
         ,vlm ...
         ,Istop,Irepeat,action] ...
...
   = prd_2d (N,X,Y,Ich1,thmax,Ich2,spmax,Ich3,spmin)

%------------------------------------------------
%  Point redistribution on 2D line
%
%  Checks for MAXIMUM ANGLE, 
%             MAXIMUM POINT SEPARATION,
%             MINIMUM POINT SEPARATION,
%
% N: number of points
%
% will redo if Irepeat=1
%------------------------------------------

%----------
% constants
%----------

   pih = 0.5*pi;

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
  ,vlm] = splc_geo (N,X,Y);

%-------
% checks
%-------

loop1 = 1;
loop2 = N;

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
      s0 = -s(2);
      thtot = (s(2)-s0)*crv(1);
      Dth = abs(thtot);
    else
      thtot = (s(K)-s(L))*crv(I);
      Dth = abs(thtot);
   end

%-----------------------------------
%  Check for angles subtended by the arcs
%
%  If I=1 add another point
%
%  If an arc is too large for I>1,
%  remove the middle point and add
%  two evenly-spaced points
%-----------------------------------

if(Ich1==1) 

 if(Dth > pih)
   disp (' prd_2d: an arc is excessive:'); disp(I);
   Istop = 1;
   return
 end

 if(Dth>thmax) 
 if(I>loop1 & I<loop2)   

 disp (' prd_2d: an arc is too large'); disp(I);

  if(I==1)

      Ax = Axint(1);Bx = Bxint(1);Cx = Cxint(1);
      Ay = Ayint(1);By = Byint(1);Cy = Cyint(1);
      XINT = (Xint(1)+Xint(2))/2.0D0;
      DX = XINT-Xint(1);
      Xnew = ((Ax*DX + Bx) * DX + Cx)*DX + X(1);
      Ynew = ((Ay*DX + By) * DX + Cy)*DX + Y(1);

      for j=N1:-1:2
        j1 = j+1;
        X(j1) =  X(j);
        Y(j1) =  Y(j);
      end
      X(2) = Xnew;
      Y(2) = Ynew;

  else

      Ax1 = Axint(L);Bx1 = Bxint(L);Cx1 = Cxint(L);
      Ay1 = Ayint(L);By1 = Byint(L);Cy1 = Cyint(L);
      XINT1 = (Xint(L)+2.0D0*Xint(I))/3.0D0;
      DX1 = XINT1-Xint(L);
      Xnew1 = ( ( Ax1*DX1 + Bx1) * DX1 + Cx1)*DX1 + X(L);
      Ynew1 = ( ( Ay1*DX1 + By1) * DX1 + Cy1)*DX1 + Y(L);

%      call INT_CUB_dog
%     +
%     +    (Xint(La),Xint(L),Xint(I),Xint(K)
%     +    ,XINT1
%     +    ,P1(La),P1(L),P1(I),P1(K),P1_int1
%     +    ,P2(La),P2(L),P2(I),P2(K),P2_int1
%     +    ,P3(La),P3(L),P3(I),P3(K),P3_int1
%     +    ,P4(La),P4(L),P4(I),P4(K),P4_int1
%     +    ,P5(La),P5(L),P5(I),P5(K),P5_int1
%     +    )

      Ax2 = Axint(I);Bx2 = Bxint(I);Cx2 = Cxint(I);
      Ay2 = Ayint(I);By2 = Byint(I);Cy2 = Cyint(I);
      XINT2 = (2.0D0*Xint(I)+Xint(K))/3.0D0;
      DX2 = XINT2-Xint(I);
      Xnew2 = ( ( Ax2*DX2 + Bx2) * DX2 + Cx2 )*DX2 + X(I);
      Ynew2 = ( ( Ay2*DX2 + By2) * DX2 + Cy2 )*DX2 + Y(I);

%      call INT_CUB_dog
%     +
%     +    (Xint(L),Xint(I),Xint(K),Xint(K1)
%     +    ,XINT2
%     +    ,P1(L),P1(I),P1(K),P1(K1),P1_int2
%     +    ,P2(L),P2(I),P2(K),P2(K1),P2_int2
%     +    ,P3(L),P3(I),P3(K),P3(K1),P3_int2
%     +    ,P4(L),P4(I),P4(K),P4(K1),P4_int2
%     +    ,P5(L),P5(I),P5(K),P5(K1),P5_int2
%     +    )

      for j=N1:-1:K
       j1 = j+1;
        X(j1) =  X(j);
        Y(j1) =  Y(j);
%       P1(j1) = P1(j);
%       P2(j1) = P2(j);
%       P3(j1) = P3(j);
%       P4(j1) = P3(j);
%       P5(j1) = P3(j);
      end

       X(K) = Xnew2;
       Y(K) = Ynew2;
%      P1(K) = P1_int2;
%      P2(K) = P2_int2;
%      P3(K) = P3_int2;
%      P4(K) = P4_int2;
%      P5(K) = P5_int2;

       X(I) = Xnew1;
       Y(I) = Ynew1;
%      P1(I) = P1_int1;
%      P2(I) = P2_int1;
%      P3(I) = P3_int1;
%      P4(I) = P4_int1;
%      P5(I) = P5_int1;

   end

      N = N1;
      disp (' One point inserted; segments :'); disp(N);
      Irepeat=1;
      action=1;
      return

   end

end
end

%-----------------------------------
%  Check for maximum point separation
%  If a segment is too long,
%  add a point in the middle
%-----------------------------------

if(Ich2==1)
if(sp>spmax)

  disp (' prd_2d_pr: segment too long'); disp(I);

  Ax = Axint(I);Bx = Bxint(I);Cx = Cxint(I);
  Ay = Ayint(I);By = Byint(I);Cy = Cyint(I);

  XINT1 = 0.50D0*(Xint(I)+Xint(K));

  DX = XINT1-Xint(I);

  Xnew = ( ( Ax*DX + Bx) * DX + Cx ) *DX + X(I);
  Ynew = ( ( Ay*DX + By) * DX + Cy ) *DX + Y(I);

%        P1L = P1(L);
%        P2L = P2(L);
%        P3L = P3(L);
%        P4L = P4(L);
%        P5L = P5(L);

%      call INT_CUB_dog
%     +
%     +    (Xint1,Xint(I),Xint(K),Xint(K1)
%     +    ,XINT1
%     +    ,P1L,P1(I),P1(K),P1(K1),P1_int
%     +    ,P2L,P2(I),P2(K),P2(K1),P2_int
%     +    ,P3L,P3(I),P3(K),P3(K1),P3_int
%     +    ,P4L,P4(I),P4(K),P4(K1),P4_int
%     +    ,P5L,P5(I),P5(K),P5(K1),P5_int
%     +    )

      for j=N1:-1:K
        j1    = j+1;
         X(j1) =  X(j);
         Y(j1) =  Y(j);
%        P1(j1) = P1(j)
%        P2(j1) = P2(j)
%        P3(j1) = P3(j)
%        P4(j1) = P4(j)
%        P5(j1) = P5(j)
      end

       X(K) =  Xnew;
       Y(K) =  Ynew;
%      P1(K) = P1_INT
%      P2(K) = P2_INT
%      P3(K) = P3_INT
%      P4(K) = P4_INT
%      P5(K) = P5_INT

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
if(I>loop1 & I<loop2)   

%-------------------------------------
%  POINTS I and I+1 WILL BE REMOVED
%  AND REPLACED BY A MID-POINT
%
%  BUT ONLY IF THE RESULTING ANGLES
%  AND POINT SEPARATIONS DO NOT EXCEED
%  THE PRE-ESTABLISHED MAXIMA
%-------------------------------------

       stemp   = 0.50D0*(  s(I)+  s(K));
       crvtemp = 0.50D0*(crv(I)+crv(K));

       SEP1 = stemp - s(I-1);
       SEP2 = s(I+2)- stemp;

       if(I==2)
         TOTANG0 = crv(I-1)*(stemp+s(2));
       else
         TOTANG0 = crv(I-1)*(stemp-s(I-2));
       end

       TOTANG1 = crvtemp *(s(I+2)-s(I-1));
       TOTANG2 = crv(I+2)*(s(I+3)-stemp);

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

      for j=K:N
          j1 = j+1;
        X(j) =  X(j1);
        Y(j) =  Y(j1);
%       P1(j) = P1(j1)
%       P2(j) = P2(j1)
%       P3(j) = P3(j1)
%       P4(j) = P4(j1)
%       P5(j) = P5(j1)
      end

      Xint1 = 0.50D0*(Xint(I)+Xint(K));

      DX = Xint1-Xint(I);

      XNEW = ( ( Ax*DX + Bx) * DX + Cx )*DX + X(I);
      YNEW = ( ( Ay*DX + By) * DX + Cy )*DX + Y(I);

%      call INT_CUB_dog
%     +
%     +    (Xint(L),Xint(I),Xint(K),Xint(K1)
%     +    ,Xint1
%     +    ,P1(L),P1(I),P1(K),P1(K1),P1_int
%     +    ,P2(L),P2(I),P2(K),P2(K1),P2_int
%     +    ,P3(L),P3(I),P3(K),P3(K1),P3_int
%     +    ,P4(L),P4(I),P4(K),P4(K1),P4_int
%     +    ,P5(L),P5(I),P5(K),P5(K1),P5_int
%     +    )

       X(I) = XNEW;
       Y(I) = YNEW;
%      P1(I) = P1_int;
%      P2(I) = P2_int;
%      P3(I) = P3_int;
%      P4(I) = P4_int;
%      P5(I) = P5_int;
	
      N  = N-1;
      N1 = N+1;

      disp (' one point removed: segments: '); disp (N);

      action=3;
      Irepeat = 1;
      return
      end
%------------------------

end
end
end

%====================
     end % of cat
%====================

%-----
% done
%-----

 return
