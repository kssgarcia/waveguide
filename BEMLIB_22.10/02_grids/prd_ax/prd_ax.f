      subroutine prd_ax 
     +
     +   (NSG
     +   ,ICH1,THMAX
     +   ,ICH2,SPMAX
     +   ,ICH3,SPMIN
     +   ,ARCT,ARET,VOLT
     +   ,XCNT,YCNT
     +   ,P1,P2,P3,P4,P5
     +   ,P6,P7
     +   ,Istop
     +   )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c  Parametrizes a line in the xy plane beginning and
c  ending at the x axis.
c
c  The line is the trace an
c  axisymmetric surface in the xy plane
c
c  The line is defined by NSG segments
c
c  The end-points are numbered 1 and NSG+1
c  Both points lie on the x axis
c
c  CHECKS FOR MAXIMUM ANGLE SUBTENDED BY THE ARCS,
c             MAXIMUM POINT SEPARATION,
c             MINIMUM POINT SEPARATION,
c
c  AND GIVES AN IMPROVED DISTRIBUTION
c
c  IT ALSO INTERPOLATES FOR THE VARIABLES P1 - P7
c
c       Variable   Example
c
c            P1:  x-vel
c            P2:  s-vel
c            P3:  surface tension
c            P4:  reference arc length
c            P5:  reference y position
c            P6:  reference crvs
c            P7:  reference crvp
c
c  SYMBOLS:
c  -------
c
c  arclt: total arc length
c
c  tol:   tolerance for crossing the upper xy semi-plane
c 
c----------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension  X(0:900), Y(0:900)
      Dimension P1(0:900),P2(0:900)
      Dimension P3(0:900)
      Dimension P4(0:900),P5(0:900)
      Dimension P6(0:900),P7(0:900)

      Dimension XC   (900),YC   (900),R  (900),S   (900)
      Dimension TH1  (900),TH2  (900),TH3(900),ORNT(900)
      Dimension AREA1(900),AREA2(900),ARE(900)
      Dimension VOL1 (900),VOL2 (900),VOL(900)
      Dimension XCN1 (900),XCN2 (900),XCN(900)
      Dimension YCN1 (900),YCN2 (900),YCN(900)
 
      Parameter (Nmax=900,tol=-0.0001)

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT

c-----------
c  constants
c-----------

      pi  = 3.14159 265358 979323 84
      pih = 0.5*pi

c--------
c preliminaries
c--------

      Isafe = 100*NSG

      arclt = P4(NSG+1)-P4(1)

c-------
c launch
c-------

      Ipass = 0

  98  Continue

      Ipass = Ipass+1

      If(NSG.GE.Nmax) then
        write (6,*) 
        write (6,*) " prd_ax: more than ",Nmax," points "
        write (6,*) "         will abort"
        write (6,*) " ------------------"
        Istop = 1
        Return
      End If

      If(Ipass.gt.Isafe) then
        write (6,*) 
        write (6,*) " prd_ax: Caught in a loop redistributing points "
        write (6,*) "         will abort"
        write (6,*) " ------------------"
        Istop = 1
        Return
      End If

c--------
c prepare 
c--------

      NSGb = NSG-2
      NSGa = NSG-1

      NSG1 = NSG+1
      NSG2 = NSG+2
 
c------------
c wrap around
c------------

       X(0) =   X(2)
       Y(0) = - Y(2)
      P1(0) =  P1(2)
      P2(0) = -P2(2)
      P3(0) =  P3(2)
      P4(0) = -P4(2)
      P5(0) = -P5(2)
      P6(0) =  P6(2)
      P7(0) =  P7(2)

       X(NSG2) =   X(NSG)
       Y(NSG2) = - Y(NSG)
      P1(NSG2) =  P1(NSG)
      P2(NSG2) = -P2(NSG)
      P3(NSG2) =  P3(NSG)
      P4(NSG2) = -P4(NSG)+2.0*arclt
      P5(NSG2) = -P5(NSG)
      P6(NSG2) =  P6(NSG)
      P7(NSG2) =  P7(NSG)

c-----------
c initialize
c-----------
 
      Istop = 0
      Iloop = 1

c-----------------
c  DEFINE THE ARCS
c-----------------

      Do 3 I= Iloop,NSG1

        If(y(I).lt.tol) Go to 94    ! point crossed the zx plane

        M = I-2
        L = I-1
        K = I+1

        XX1 = X(L)
        YY1 = Y(L)
        XX2 = X(I)
        YY2 = Y(I)
        XX3 = X(K)
        YY3 = Y(K)

        call arc_ax 
     +
     +    (I
     +    ,XX1,XX2,XX3
     +    ,YY1,YY2,YY3
     +    ,XC(I),YC(I),R(I)
     +    ,TH1(I),TH2(I),TH3(I),ORNT(I)
     +    ,AREA1(I),AREA2(I),ARE(I)
     +    ,VOL1(I),VOL2(I),VOL(I)
     +    ,XCN1(I),XCN2(I),XCN(I)
     +    ,YCN1(I),YCN2(I),YCN(I)
     +    ,Istop
     +    )
 
c------------------------------------
c  CHECK ANGLES SUBTENDED BY THE ARCS
c------------------------------------

      If(Ich1.EQ.1) then

        THETOT = TH3(I)-TH1(I)
        Dth    = ABS(THETOT)

        If(Dth.GE.pih  ) Go to 96
        If(Dth.GE.THMAX) Go to 97

      End If

c--------------------------------
c  CHECK MAXIMUM POINT SEPARATION
c--------------------------------

      If(Ich2.EQ.1) then

       SP1 = R(I)*ABS(TH2(I)-TH1(I))
       If(SP1.gt.SPMAX.and.I.ne.1)    Go to 91

       SP2 = R(I)*ABS(TH3(I)-TH2(I))
       If(SP2.GT.SPMAX.and.I.ne.NSG1) Go to 92

      End If
 
c------------------------------------
c  CHECK FOR MINIMUM POINT SEPARATION
c------------------------------------
 
      If(Ich3.eq.0) GO TO 3
      If(I.eq.NSG1) GO TO 3

      SP2 = R(I)*ABS(TH3(I)-TH2(I))

      If(SP2.GT.SPMIN) GO TO 3

c---------------------------------------------------------
c A point will be removed only if the resulting
c angles and point separations do not violate
c the first two restrictions
c
c Will consider the cases I=1,
c                         I=NSG
c                         I.ne.1,NSG1 separately
c---------------------------------------------------------
 
c--------------------
      If(I.eq.1) then
c--------------------
c Point 1 will stay
c Point 2 may be removed

         XTMP =  X(3)
         YTMP = -Y(3)

         call arc_ax 
     +
     +      (I
     +      ,XTMP,X(1),X(3)
     +      ,YTMP,Y(1),Y(3)
     +      ,XCT,YCT,RT
     +      ,TH1T,TH2T,TH3T,ORNTT
     +      ,AREA1T,AREA2T,ARET
     +      ,VOL1T,VOL2T,VOLT
     +      ,XCN1T,XCN2T,XCNT
     +      ,YCN1T,YCN2T,YCNT
     +      ,Istop
     +      )

         TOTANG0 = ABS(TH3T-TH1T)

         call arc_ax 
     +
     +       (I
     +       ,X(1),X(3),X(4)
     +       ,Y(1),Y(3),Y(4)
     +       ,XCT,YCT,RT
     +       ,TH1T,TH2T,TH3T,ORNTT
     +       ,AREA1T,AREA2T,ARET
     +       ,VOL1T,VOL2T,VOLT
     +       ,XCN1T,XCN2T,XCNT
     +       ,YCN1T,YCN2T,YCNT
     +       ,Istop
     +       )

         TOTANG1 = ABS(TH3T-TH1T)
         SEP1    = RT * ABS(TH2T-TH1T)

         If(TOTANG0.LT.THMAX.AND.
     +      TOTANG1.LT.THMAX.AND.
     +         SEP1.LT.SPMAX) 
     +   Go to 93

c---------------------------
      Else If(I.eq.NSG) then
c---------------------------
c Point NSG1 will stay
c Point NSG  may be removed

         call arc_ax 
     +
     +       (I
     +       ,X(NSGb),X(NSGa),X(NSG1)
     +       ,Y(NSGb),Y(NSGa),Y(NSG1)
     +       ,XCT,YCT,RT
     +       ,TH1T,TH2T,TH3T,ORNTT
     +       ,AREA1T,AREA2T,ARET
     +       ,VOL1T,VOL2T,VOLT
     +       ,XCN1T,XCN2T,XCNT
     +       ,YCN1T,YCN2T,YCNT
     +       ,Istop
     +       )

         TOTANG0 = ABS(TH3T-TH1T)

         XTMP =  X(NSGa)
         YTMP = -Y(NSGa)

         call arc_ax 
     +               (I
     +               ,X(NSGa),X(NSG1),XTMP
     +               ,Y(NSGa),Y(NSG1),YTMP
     +               ,XCT,YCT,RT
     +               ,TH1T,TH2T,TH3T,ORNTT
     +               ,AREA1T,AREA2T,ARET
     +               ,VOL1T,VOL2T,VOLT
     +               ,XCN1T,XCN2T,XCNT
     +               ,YCN1T,YCN2T,YCNT
     +               ,Istop
     +               )

         TOTANG1 = ABS(TH3T-TH1T)
         SEP1    = RT * ABS(TH2T-TH1T)

         If( TOTANG0.LT.THMAX.AND.
     +       TOTANG1.LT.THMAX.AND.
     +          SEP1.LT.SPMAX)
     +   Go to 93

c---
      Else
c---
         ANGL = 0.5*(TH2(I)+TH3(I))
         XTMP = XC(I)+R(I)*COS(ANGL)
         YTMP = YC(I)+R(I)*SIN(ANGL)

         If(I.eq.2) then
           XXX =  XTMP
           YYY = -YTMP
         Else
           XXX = X(M)
           YYY = Y(M)
         End If

         call arc_ax 
     +
     +       (I
     +       ,XXX,X(L),XTMP
     +       ,YYY,Y(L),YTMP
     +       ,XCT,YCT,RT
     +       ,TH1T,TH2T,TH3T,ORNTT
     +       ,AREA1T,AREA2T,ARET
     +       ,VOL1T,VOL2T,VOLT
     +       ,XCN1T,XCN2T,XCNT
     +       ,YCN1T,YCN2T,YCNT
     +       ,Istop
     +       )

         TOTANG0 = ABS(TH3T-TH1T)

         call arc_ax 
     +
     +       (I
     +       ,X(L),XTMP,X(I+2)
     +       ,Y(L),YTMP,Y(I+2)
     +       ,XCT,YCT,RT
     +       ,TH1T,TH2T,TH3T,ORNTT
     +       ,AREA1T,AREA2T,ARET
     +       ,VOL1T,VOL2T,VOLT
     +       ,XCN1T,XCN2T,XCNT
     +       ,YCN1T,YCN2T,YCNT
     +       ,Istop
     +       )

         TOTANG1 = ABS(TH3T-TH1T)
         SEP1    = RT * ABS(TH2T-TH1T)
         SEP2    = RT * ABS(TH3T-TH2T)

         If(I.eq.NSGa) then
           XXX =  XTMP
           YYY = -YTMP
         Else
           XXX = X(I+3)
           YYY = Y(I+3)
         End If

         call arc_ax 
     +
     +       (I
     +       ,XTMP,X(I+2),XXX
     +       ,YTMP,Y(I+2),YYY
     +       ,XCT,YCT,RT
     +       ,TH1T,TH2T,TH3T,ORNTT
     +       ,AREA1T,AREA2T,ARET
     +       ,VOL1T,VOL2T,VOLT
     +       ,XCN1T,XCN2T,XCNT
     +       ,YCN1T,YCN2T,YCNT
     +       ,Istop
     +       )

         TOTANG2 = ABS(TH3T-TH1T)

         If(TOTANG0.LT.THMAX.AND.
     +      TOTANG1.LT.THMAX.AND.
     +      TOTANG2.LT.THMAX.AND.
     +         SEP1.LT.SPMAX.AND.
     +         SEP2.LT.SPMAX)
     +    Go to 93

c---
      End If
c---

  3   Continue
 
c----------------------------------------
c  We now have an acceptable distribution
c----------------------------------------
 
      Go to 29
 
c--------------------
c There is a problem
c execution will stop
c--------------------
 
  94  Continue

      write (6,207) I
      write (6,201)
      write (1,207) I
      write (1,201)

      Istop = 1
      Return
 
  96  Continue

      write (6,200) I,THETOT,THMAX
      write (6,201)
      write (1,200) I,THETOT,THMAX
      write (1,201)

      Istop = 1
      Return

c------

  97  Continue
 
c---------------------------------
c  IF AN ARC IS TOO LARGE, 
c  REMOVE MIDDLE POINT,
c  AND ADD TWO POINTS ALONG THE ARC
c
c  Except for points 1 and NSG1
c  that are treated in a special way
c---------------------------------

      write (6,*) 
      write (6,200) I,THETOT,THMAX

c---
      If(I.eq.1) then
c---
c will perform blended backward and forward interpolation
c---

      angl = 0.5D0*(TH3(1)+TH2(1))

      call INT_QUAD1
     +
     +    (TH1(1),TH2(1),TH3(1),angl
     +    ,P1(0),P1(1),P1(2),P1_int1
     +    ,P2(0),P2(1),P2(2),P2_int1
     +    ,P3(0),P3(1),P3(2),P3_int1
     +    ,P4(0),P4(1),P4(2),P4_int1
     +    ,P5(0),P5(1),P5(2),P5_int1
     +    ,P6(0),P6(1),P6(2),P6_int1
     +    ,P7(0),P7(1),P7(2),P7_int1
     +    )

      Xtmp1 = XC(1)+R(1)*Dcos(angl)
      Ytmp1 = YC(1)+R(1)*Dsin(angl)

      call arc_ax 
     +
     +     (I
     +     ,X(1),X(2),X(3)
     +     ,Y(1),Y(2),Y(3)
     +     ,XC(2),YC(2),R(2)
     +     ,TH1(2),TH2(2),TH3(2),ORNT(2)
     +     ,AREA1T,AREA2T,ARET
     +     ,VOL1T,VOL2T,VOLT
     +     ,XCN1T,XCN2T,XCNT
     +     ,YCN1T,YCN2T,YCNT
     +     ,Istop
     +     )

      angl = 0.5D0*(TH1(2)+TH2(2))

      call INT_QUAD1
     +
     +      (TH1(2),TH2(2),TH3(2),angl
     +      ,P1(1),P1(2),P1(3),P1_int2
     +      ,P2(1),P2(2),P2(3),P2_int2
     +      ,P3(1),P3(2),P3(3),P3_int2
     +      ,P4(1),P4(2),P4(3),P4_int2
     +      ,P5(1),P5(2),P5(3),P5_int2
     +      ,P6(1),P6(2),P6(3),P6_int2
     +      ,P7(1),P7(2),P7(3),P7_int2
     +      )

      Xtmp2 = XC(2)+R(2)*Dcos(angl)
      Ytmp2 = YC(2)+R(2)*Dsin(angl)

      DO J=NSG1,2,-1
        J1     = J+1
        X (J1) = X (J)
        Y (J1) = Y (J)
        P1(J1) = P1(J)
        P2(J1) = P2(J)
        P3(J1) = P3(J)
        P4(J1) = P4(J)
        P5(J1) = P5(J)
        P6(J1) = P6(J)
        P7(J1) = P7(J)
      END DO

      X(2) = 0.5D0*(Xtmp1+Xtmp2)
      Y(2) = 0.5D0*(Ytmp1+Ytmp2)

      P1(2) = 0.5D0*(P1_int1 + P1_int2)
      P2(2) = 0.5D0*(P2_int1 + P2_int2)
      P3(2) = 0.5D0*(P3_int1 + P3_int2)
      P4(2) = 0.5D0*(P4_int1 + P4_int2)
      P5(2) = 0.5D0*(P5_int1 + P5_int2)
      P6(2) = 0.5D0*(P6_int1 + P6_int2)
      P7(2) = 0.5D0*(P7_int1 + P7_int2)

c---
      Else If(I.eq.NSG1) then
c---
c  will perform blended backward and forward
c  interpolation

      angl = 0.5*(TH3(L)+TH2(L))

      La = L-1

      call INT_QUAD1
     +
     +    (TH1(L),TH2(L),TH3(L),angl
     +    ,P1(La),P1(L),P1(I),P1_int1
     +    ,P2(La),P2(L),P2(I),P2_int1
     +    ,P3(La),P3(L),P3(I),P3_int1
     +    ,P4(La),P4(L),P4(I),P4_int1
     +    ,P5(La),P5(L),P5(I),P5_int1
     +    ,P6(La),P6(L),P6(I),P6_int1
     +    ,P7(La),P7(L),P7(I),P7_int1
     +    )

      XTMP1 = XC(L)+R(L)*COS(angl)
      YTMP1 = YC(L)+R(L)*SIN(angl)

      angl = 0.5*(TH1(I)+TH2(I))

      call INT_QUAD1
     +
     +     (TH1(I),TH2(I),TH3(I),angl
     +     ,P1(L),P1(I),P1(K),P1_int2
     +     ,P2(L),P2(I),P2(K),P2_int2
     +     ,P3(L),P3(I),P3(K),P3_int2
     +     ,P4(L),P4(I),P4(K),P4_int2
     +     ,P5(L),P5(I),P5(K),P5_int2
     +     ,P6(L),P6(I),P6(K),P6_int2
     +     ,P7(L),P7(I),P7(K),P7_int2
     +     )

      XTMP2 = XC(I)+R(I)*Dcos(angl)
      YTMP2 = YC(I)+R(I)*Dsin(angl)

      X (K) = X (NSG1)
      Y (K) = Y (NSG1)
      P1(K) = P1(NSG1)
      P2(K) = P2(NSG1)
      P3(K) = P3(NSG1)
      P4(K) = P4(NSG1)
      P5(K) = P5(NSG1)
      P6(K) = P6(NSG1)
      P7(K) = P7(NSG1)

       X(I) = 0.5*(XTMP1+XTMP2)
       Y(I) = 0.5*(YTMP1+YTMP2)
      P1(I) = 0.5*(P1_int1+P1_int2)
      P2(I) = 0.5*(P2_int1+P2_int2)
      P3(I) = 0.5*(P3_int1+P3_int2)
      P4(I) = 0.5*(P4_int1+P4_int2)
      P5(I) = 0.5*(P5_int1+P5_int2)
      P6(I) = 0.5*(P6_int1+P6_int2)
      P7(I) = 0.5*(P7_int1+P7_int2)

c---
      Else
c---
      Dth  = THETOT/3.0D0
      TH13 = TH1(I)+    Dth
      TH23 = TH1(I)+2.0*Dth

      call INT_QUAD1
     +
     +    (TH1(I),TH2(I),TH3(I),TH13
     +    ,P1(L),P1(I),P1(K),P1_int1
     +    ,P2(L),P2(I),P2(K),P2_int1
     +    ,P3(L),P3(I),P3(K),P3_int1
     +    ,P4(L),P4(I),P4(K),P4_int1
     +    ,P5(L),P5(I),P5(K),P5_int1
     +    ,P6(L),P6(I),P6(K),P6_int1
     +    ,P7(L),P7(I),P7(K),P7_int1
     +    )

      call INT_QUAD1
     +
     +    (TH1(I),TH2(I),TH3(I),TH23
     +    ,P1(L),P1(I),P1(K),P1_int2
     +    ,P2(L),P2(I),P2(K),P2_int2
     +    ,P3(L),P3(I),P3(K),P3_int2
     +    ,P4(L),P4(I),P4(K),P4_int2
     +    ,P5(L),P5(I),P5(K),P5_int2
     +    ,P6(L),P6(I),P6(K),P6_int2
     +    ,P7(L),P7(I),P7(K),P7_int2
     +    )

      Do j=NSG1,K,-1
        j1 = j+1
        X (j1) = X(j)
        Y (j1) = Y(j)
        P1(j1) = P1(j)
        P2(j1) = P2(j)
        P3(j1) = P3(j)
        P4(j1) = P4(j)
        P5(j1) = P5(j)
        P6(j1) = P6(j)
        P7(j1) = P7(j)
      End Do

      X (K) = XC(I)+R(I)*DCOS(TH23)
      Y (K) = YC(I)+R(I)*DSIN(TH23)
      P1(K) = P1_int2
      P2(K) = P2_int2
      P3(K) = P3_int2
      P4(K) = P4_int2
      P5(K) = P5_int2
      P6(K) = P6_int2
      P7(K) = P7_int2

      X (I) = XC(I)+R(I)*DCOS(TH13)
      Y (I) = YC(I)+R(I)*DSIN(TH13)
      P1(I) = P1_int1
      P2(I) = P2_int1
      P3(I) = P3_int1
      P4(I) = P4_int1
      P5(I) = P5_int1
      P6(I) = P6_int1
      P7(I) = P7_int1

c---
      End If
c---

      NSG = NSG+1

      write (6,206) NSG

      Go to 98
 
c--------------------------
c  IF A SEGMENT IS TOO LONG,
c  ADD A POINT IN THE MIDDLE
C---------------------------
 
  91  Continue

      write (6,203) I

      ANGL1 = 0.5*(TH2(L)+TH3(L))
      Xint1 = XC(L)+R(L)*COS(ANGL1)
      Yint1 = YC(L)+R(L)*SIN(ANGL1)

      La = L-1

      call INT_QUAD1
     +
     +    (TH1(L),TH2(L),TH3(L),ANGL1
     +    ,P1(La),P1(L),P1(I),P1_int1
     +    ,P2(La),P2(L),P2(I),P2_int1
     +    ,P3(La),P3(L),P3(I),P3_int1
     +    ,P4(La),P4(L),P4(I),P4_int1
     +    ,P5(La),P5(L),P5(I),P5_int1
     +    ,P6(La),P6(L),P6(I),P6_int1
     +    ,P7(La),P7(L),P7(I),P7_int1
     +    )

      ANGL2 = 0.5*(TH1(I)+TH2(I))
      Xint2 = XC(I)+R(I)*COS(ANGL2)
      Yint2 = YC(I)+R(I)*SIN(ANGL2)

      call INT_QUAD1
     +
     +    (TH1(I),TH2(I),TH3(I),ANGL2
     +    ,P1(L),P1(I),P1(K),P1_int2
     +    ,P2(L),P2(I),P2(K),P2_int2
     +    ,P3(L),P3(I),P3(K),P3_int2
     +    ,P4(L),P4(I),P4(K),P4_int2
     +    ,P5(L),P5(I),P5(K),P5_int2
     +    ,P6(L),P6(I),P6(K),P6_int2
     +    ,P7(L),P7(I),P7(K),P7_int2
     +    )

        XINT = 0.5*(  Xint1 +   Xint2)
        YINT = 0.5*(  Yint1 +   Yint2)
      P1_int = 0.5*(P1_int1 + P1_int2)
      P2_int = 0.5*(P2_int1 + P2_int2)
      P3_int = 0.5*(P3_int1 + P3_int2)
      P4_int = 0.5*(P4_int1 + P4_int2)
      P5_int = 0.5*(P5_int1 + P5_int2)
      P6_int = 0.5*(P6_int1 + P6_int2)
      P7_int = 0.5*(P7_int1 + P7_int2)

      Do J=NSG1,I,-1
        J1 = J+1
        X(J1) = X(J)
        Y(J1) = Y(J)
        P1(j1) = P1(j)
        P2(j1) = P2(j)
        P3(j1) = P3(j)
        P4(j1) = P4(j)
        P5(j1) = P5(j)
        P6(j1) = P6(j)
        P7(j1) = P7(j)
      End Do

      X (I)  = XINT
      Y (I)  = YINT
      P1(I) = P1_INT
      P2(I) = P2_INT
      P3(I) = P3_INT
      P4(I) = P4_INT
      P5(I) = P5_INT
      P6(I) = P6_INT
      P7(I) = P7_INT

      NSG = NSG+1

      write (6,206) NSG

      Go to 98
 
c-------

  92  Continue

      write (6,204) I

      K1 = K+1

      call arc_ax 
     +
     +     (I
     +     ,X(I),X(K),X(K1)
     +     ,Y(I),Y(K),Y(K1)
     +     ,XC(K),YC(K),R(K)
     +     ,TH1(K),TH2(K),TH3(K),ORNT(K)
     +     ,AREA1(K),AREA2(K),ARE(K)
     +     ,VOL1(K),VOL2(K),VOL(K)
     +     ,XCN1(K),XCN2(K),XCN(K)
     +     ,YCN1(K),YCN2(K),YCN(K)
     +     ,Istop
     +     )

      ANGL1 = 0.5*(TH2(I)+TH3(I))

      call INT_QUAD1
     +
     +    (TH1(I),TH2(I),TH3(I),ANGL1
     +    ,P1(L),P1(I),P1(K),P1_int1
     +    ,P2(L),P2(I),P2(K),P2_int1
     +    ,P3(L),P3(I),P3(K),P3_int1
     +    ,P4(L),P4(I),P4(K),P4_int1
     +    ,P5(L),P5(I),P5(K),P5_int1
     +    ,P6(L),P6(I),P6(K),P6_int1
     +    ,P7(L),P7(I),P7(K),P7_int1
     +    )

      XINT1 = XC(I)+R(I)*COS(ANGL1)
      YINT1 = YC(I)+R(I)*SIN(ANGL1)

      ANGL2 = 0.5*(TH1(K)+TH2(K))

      call INT_QUAD1
     +
     +    (TH1(K),TH2(K),TH3(K),ANGL2
     +    ,P1(I),P1(K),P1(K1),P1_int2
     +    ,P2(I),P2(K),P2(K1),P2_int2
     +    ,P3(I),P3(K),P3(K1),P3_int2
     +    ,P4(I),P4(K),P4(K1),P4_int2
     +    ,P5(I),P5(K),P5(K1),P5_int2
     +    ,P6(I),P6(K),P6(K1),P6_int2
     +    ,P7(I),P7(K),P7(K1),P7_int2
     +    )

      XINT2 = XC(K)+R(K)*COS(ANGL2)
      YINT2 = YC(K)+R(K)*SIN(ANGL2)

        XINT = 0.5*(  XINT1 +   XINT2)
        YINT = 0.5*(  YINT1 +   YINT2)
      P1_INT = 0.5*(P1_INT1 + P1_INT2)
      P2_INT = 0.5*(P2_INT1 + P2_INT2)
      P3_INT = 0.5*(P3_INT1 + P3_INT2)
      P4_INT = 0.5*(P4_INT1 + P4_INT2)
      P5_INT = 0.5*(P5_INT1 + P5_INT2)
      P6_INT = 0.5*(P6_INT1 + P6_INT2)
      P7_INT = 0.5*(P7_INT1 + P7_INT2)

      Do J=NSG1,K,-1
        J1 = J+1
        X (J1) = X (J)
        Y (J1) = Y (J)
        P1(j1) = P1(j)
        P2(j1) = P2(j)
        P3(j1) = P3(j)
        P4(j1) = P4(j)
        P5(j1) = P5(j)
        P6(j1) = P6(j)
        P7(j1) = P7(j)
      End Do

      X (K) =   XINT
      Y (K) =   YINT
      P1(K) = P1_INT
      P2(K) = P2_INT
      P3(K) = P3_INT
      P4(K) = P4_INT
      P5(K) = P5_INT
      P6(K) = P6_INT
      P7(K) = P7_INT

      NSG = NSG+1

      write (6,206) NSG

      GO TO 98
 
c-------------------------------
c  IF A SEGMENT IS TOO SHORT,
c  REMOVE END-POINTS
c  AND ADD A POINT IN THE MIDDLE
c-------------------------------
 
  93  Continue

      write (6,209) I

c---
      If(I.ne.NSG) then
c---
      Do J=K,NSG1
       J1 = J+1
       X (J) = X(J1)
       Y (J) = Y(J1)
       P1(j) = P1(j1)
       P2(j) = P2(j1)
       P3(j) = P3(j1)
       P4(j) = P4(j1)
       P5(j) = P5(j1)
       P6(j) = P6(j1)
       P7(j) = P7(j1)
      End Do

      If(I.NE.1) then

        X(I) = XTMP
        Y(I) = YTMP
        ANGL = 0.5*(TH2(I)+TH3(I))

        call INT_QUAD1
     +
     +   (TH1(I),TH2(I),TH3(I),ANGL
     +   ,P1(L),P1(I),P1(K),P1_int
     +   ,P2(L),P2(I),P2(K),P2_int
     +   ,P3(L),P3(I),P3(K),P3_int
     +   ,P4(L),P4(I),P4(K),P4_int
     +   ,P5(L),P5(I),P5(K),P5_int
     +   ,P6(L),P6(I),P6(K),P6_int
     +   ,P7(L),P7(I),P7(K),P7_int
     +   )

        P1(I) = P1_int
        P2(I) = P2_int
        P3(I) = P3_int
        P4(I) = P4_int
        P5(I) = P5_int
        P6(I) = P6_int
        P7(I) = P7_int

      End If
c---
      Else
c---
        X (NSG) = X (NSG1)
        Y (NSG) = Y (NSG1)
        P1(NSG) = P1(NSG1)
        P2(NSG) = P2(NSG1)
        P3(NSG) = P3(NSG1)
        P4(NSG) = P4(NSG1)
        P5(NSG) = P5(NSG1)
        P6(NSG) = P6(NSG1)
        P7(NSG) = P7(NSG1)
c---
      End If
c---

      NSG = NSG-1

      write (6,208) NSG

      Go to 98
 
c---------------------------------------------
 
  29  Continue

c-----------------------------
c THE SUCCESSFUL DISTRIBUTION
c-----------------------------

c-----------------------
c Compute the arc length
c-----------------------

      S(1) = 0.0D0

      Do i=2,NSG1
       ia = i-1
       S(i) = S(ia)+0.5D0*(R(i)*ABS(TH2(i) -TH1(i) )
     +                   +R(ia)*ABS(TH3(ia)-TH2(ia)))
      End Do

      S(NSG2) = 2.0*S(NSG1)- S(NSG)

c--------------------------
c Compute:
c           Surface area,
c           Volume,
c           Volume-centroid
c--------------------------

      ARET = 0.0D0
      VOLT = 0.0D0
      XCNT = 0.0D0
      YCNT = 0.0D0

      Do i=2,NSG1
        ia = i-1
        ARET = ARET + 0.5D0*(AREA1(i)+AREA2(ia))
        VOLT = VOLT + 0.5D0*(VOL1 (i)+VOL2 (ia))
        XCNT = XCNT + 0.5D0*(XCN1 (i)+XCN2 (ia))
        YCNT = YCNT + 0.5D0*(YCN1 (i)+YCN2 (ia))
      End Do
 
      ARCT = S(NSG1)
      ARET = ARET
      VOLT = VOLT
      XCNT = XCNT/VOLT
      YCNT = YCNT/VOLT

c-----
c Done
c-----

 100  Format (1X,F10.5,2X,F10.5,10X,F10.5,2X,F10.5)
 200  Format (' ARC',I3,' IS TOO LARGE',2(F10.5))
 201  Format (' EXECUTION WILL BE INTERRUPTED')
 203  Format (' SEGMENT 1 OF ARC',I3,' IS TOO LARGE')
 204  Format (' SEGMENT 2 OF ARC',I3,' IS TOO LARGE')
 206  Format (' ONE POINT INSERTED',/,' TOTAL NSG ',I3)
 207  Format (' POINT',I3,' CROSSED THE AXIS')
 208  Format (' ONE POINT REMOVED',/,' TOTAL NSG ',I3)
 209  Format (' SEGMENT 2 OF ARC',I3,' IS TOO SMALL')

      Return
      End

c===========================================

      subroutine INT_QUAD1
     +
     +  (x1,x2,x3,x
     +  ,u1,u2,u3,u
     +  ,v1,v2,v3,v
     +  ,c1,c2,c3,c
     +  ,w1,w2,w3,w
     +  ,z1,z2,z3,z
     +  ,t1,t2,t3,t
     +  ,q1,q2,q3,q
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------
c Quadratic Lagrange Interpolation
c---------------------------------

      Implicit Double Precision (A-H,O-Z)

      p1 = (x-x2)*(x-x3)/((x1-x2)*(x1-x3))
      p2 = (x-x1)*(x-x3)/((x2-x1)*(x2-x3))
      p3 = (x-x1)*(x-x2)/((x3-x1)*(x3-x2))

      u = u1*p1 + u2*p2 + u3*p3
      v = v1*p1 + v2*p2 + v3*p3
      c = c1*p1 + c2*p2 + c3*p3
      w = w1*p1 + w2*p2 + w3*p3
      z = z1*p1 + z2*p2 + z3*p3
      t = t1*p1 + t2*p2 + t3*p3
      q = q1*p1 + q2*p2 + q3*p3

c-----
c Done
c-----

      Return
      End
