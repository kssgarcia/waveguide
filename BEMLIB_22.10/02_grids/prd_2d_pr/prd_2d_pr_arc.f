      subroutine prd_2d_pr_arc
     +
     +   (NSG
     +   ,RL
     +   ,ICH1,THMAX
     +   ,ICH2,SPMAX
     +   ,ICH3,SPMIN
     +   ,ICH4,thrs1,thrs2
     +   ,ARCT,ARET
     +   ,XCNT,YCNT
     +   ,P1,P2,P3
     +   ,P4,P5
     +   ,Istop
     +   )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c  Point redistribution on a periodic line with
c  period RL
c
c  Point NSG+1 is the periodic repetition of point 1
c
c  Checks for MAXIMUM ANGLE, 
c             MAXIMUM POINT SEPARATION,
c             MINIMUM POINT SEPARATION,
c
c  The algorithm
c  resets the periodic line so that the first point
c  lies inside the window (thrs1, thrs2)
c  
c  and gives an improved distribution
c
c  It also interpolates for the property variables P1-P5
c
c    Variable Example
c
c    P1:  x-vel
c    P2:  s-vel
c    P3:  surface tension
c    P4:  arc length
c    P5:  unused
c
c SYMBOLS:
c --------
c
c XC, YC: Arc center
c R     : Arc radius 
c
c----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension X (0:513),Y (0:513)
      Dimension P1(0:513),P2(0:513),P3(0:513)
      Dimension P4(0:513),P5(0:513)

      Dimension    XC(513),   YC(513),  R(513),   S(513)
      Dimension   TH1(513),  TH2(513),TH3(513),ORNT(513)
      Dimension AREA1(513),AREA2(513),ARE(513)
      Dimension  XCN1(513), XCN2(513),XCN(513)
      Dimension  YCN1(513), YCN2(513),YCN(513)

      Parameter (Nmax=513)

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT

c----------
c constants
c----------

      pi  = 3.14159 265358 979 32384  D0
      pih = 0.5D0*pi

c--------
c prepare
c--------

      Isafe = 100*NSG

      arclt = P4(NSG+1)-P4(1)

c--------------------------

      Ipass = 0
  98  Continue
      Ipass = Ipass+1

      If(NSG.GE.Nmax) then
        write (6,*) 
        write (6,*) " More than ",Nmax," points "
        write (6,*) " I will stop "
        write (6,*) 
        Istop = 1
        Return
      End If
 
      If(Ipass.gt.Isafe) then
        write (6,*) 
        write(6,*) " Caught in a loop redistributing points "
        write(6,*) " I will stop "
        write (6,*) 
        Istop = 1
        Return
      End If

      NSG1 = NSG+1
      NSG2 = NSG+2

c---
c wrap around
c---

       X(0) =  X(NSG)-RL
       Y(0) =  Y(NSG)
      P1(0) = P1(NSG)
      P2(0) = P2(NSG)
      P3(0) = P3(NSG)
      P4(0) = P4(NSG)-arclt
      P5(0) = P5(NSG)

       X(NSG1) =  X(1)+RL
       Y(NSG1) =  Y(1)
      P1(NSG1) = P1(1)
      P2(NSG1) = P2(1)
      P3(NSG1) = P3(1)
      P4(NSG1) = P4(1)+arclt
      P5(NSG1) = P5(1)

       X(NSG2) =  X(2)+RL
       Y(NSG2) =  Y(2)
      P1(NSG2) = P1(2)
      P2(NSG2) = P2(2)
      P3(NSG2) = P3(2)
      P4(NSG2) = P4(2)+arclt
      P5(NSG2) = P5(2)

      Istop = 0
      Iloop = 1

c----------------
c Define the arcs
c----------------

      Do 3 I=Iloop,NSG1

      M = I-2
      L = I-1
      K = I+1

      XX1 = X(L)
      YY1 = Y(L)
      XX2 = X(I)
      YY2 = Y(I)
      XX3 = X(K)
      YY3 = Y(K)

      call arc_2d 
     +
     +  (I
     +  ,XX1,XX2,XX3
     +  ,YY1,YY2,YY3
     +  ,XC(I),YC(I),R(I)
     +  ,TH1(I),TH2(I),TH3(I),ORNT(I)
     +  ,AREA1(I),AREA2(I),ARE(I)
     +  ,XCN1(I),XCN2(I),XCN(I)
     +  ,YCN1(I),YCN2(I),YCN(I)
     +  ,Istop
     +  )

      If(I.eq.NSG1) GO TO 3

c---
c prepare
c---

      SP1 = R(I)*Dabs(TH2(I)-TH1(I))
      SP2 = R(I)*Dabs(TH3(I)-TH2(I))

c---------------------------------------
c Check for angles subtended by the arcs
c---------------------------------------
 
      If(Ich1.EQ.1) then

        THETOT  = TH3(I)-TH1(I)
        Dth     = abs(THETOT)

        If(Dth.GE.pih  ) Go to 96
        If(Dth.GE.THMAX) Go to 97

      End If

c-----------------------------------
c Check for maximum point separation
c-----------------------------------

      If(Ich2.EQ.1) then

        If(SP1.GT.SPMAX.and.I.NE.1)  Go to 91
        If(SP2.GT.SPMAX)             Go to 92

      End If

c-----------------------------------
c Check for minimum point separation
c-----------------------------------

      If(Ich3.eq.0) Go to 3    ! skip
      If(I.eq.NSG1) Go to 3

      If(SP2.GT.SPMIN) Go to 3  ! skip

c-------------------------------------
c  TWO POINTS WILL BE REMOVED
c  AND REPLACED BY A MID_POINT
c  ONLY IF THE RESULTING ANGLES
c  AND POINT SEPARATIONS DO NOT EXCEED
c  THE PRE-ESTABLISHED THRESHOLD
c-------------------------------------

       ANGLIA = 0.5D0*(TH2(I)+TH3(I))

       XTMP = XC(I) + R(I)*Dcos(ANGLIA)  ! mid-point 
       YTMP = YC(I) + R(I)*Dsin(ANGLIA)  ! mid-point 
 
c---
c properties of the candidate distriution:
c---

       If(I.NE.1) then
         XXX = X(M)
         YYY = Y(M)
       Else
         XXX = X(NSG-1)-RL
         YYY = Y(NSG-1) 
       End If

       call arc_2d 
     +
     +   (I
     +   ,XXX,X(L),XTMP
     +   ,YYY,Y(L),YTMP
     +   ,XCT,YCT,RT
     +   ,TH1T,TH2T,TH3T,ORNTT
     +   ,AREA1T,AREA2T,ARET
     +   ,XCN1T,XCN2T,XCNT
     +   ,YCN1T,YCN2T,YCNT
     +   ,Istop
     +   )

       TOTANG0 = Dabs(TH3T-TH1T)

       call arc_2d 
     +
     +    (I
     +    ,X(L),XTMP,X(I+2)
     +    ,Y(L),YTMP,Y(I+2)
     +    ,XCT,YCT,RT
     +    ,TH1T,TH2T,TH3T,ORNTT
     +    ,AREA1T,AREA2T,ARET
     +    ,XCN1T,XCN2T,XCNT
     +    ,YCN1T,YCN2T,YCNT
     +    ,Istop
     +    )

       TOTANG1 = Dabs(TH3T-TH1T)

       SEP1 = RT * abs(TH2T-TH1T)
       SEP2 = RT * abs(TH3T-TH2T)

       If(I.NE.NSG) then
         XXX = X(I+3)
         YYY = Y(I+3)
       Else
         XXX = X(3)+RL
         YYY = Y(3) 
       End If

       call arc_2d 
     +
     +      (I
     +      ,XTMP,X(I+2),XXX
     +      ,YTMP,Y(I+2),YYY
     +      ,XCT,YCT,RT
     +      ,TH1T,TH2T,TH3T,ORNTT
     +      ,AREA1T,AREA2T,ARET
     +      ,XCN1T,XCN2T,XCNT
     +      ,YCN1T,YCN2T,YCNT
     +      ,Istop
     +      )

       TOTANG2 = abs(TH3T-TH1T)

       If(     TOTANG0.LT.THMAX
     +    .AND.TOTANG1.LT.THMAX
     +    .AND.TOTANG2.LT.THMAX
     +    .AND.   SEP1.LT.SPMAX
     +    .AND.   SEP2.LT.SPMAX) GO TO 93

  3   Continue

c----------------------------
c We now have an acceptable
c point distribution
c----------------------------

      Go to 29

c--------------------
c  There is a problem;
c  execution will stop
c--------------------

  96  Continue

      write (6,200) I,Dth,THMAX
      write (6,201)

      Istop = 1
      Return

c-----

  97  Continue

c-----------------------------------
c  If an arc is too large,
c  remove the middle point and add
c  two evenly-spaced points
c-----------------------------------

      write (6,*)
      write (6,200) I,THETOT,THMAX
      write (6,*)

      Dth  = THETOT/3.0
      TH13 = TH1(I)+    Dth
      TH23 = TH1(I)+2.0*Dth

      call prd_2d_pr_arc_interp
     +
     +    (TH1(I),TH2(I),TH3(I),TH13
     +    ,P1(L),P1(I),P1(K),P1_int1
     +    ,P2(L),P2(I),P2(K),P2_int1
     +    ,P3(L),P3(I),P3(K),P3_int1
     +    ,P4(L),P4(I),P4(K),P4_int1
     +    ,P5(L),P5(I),P5(K),P5_int1
     +    )

      call prd_2d_pr_arc_interp
     +
     +    (TH1(I),TH2(I),TH3(I),TH23
     +    ,P1(L),P1(I),P1(K),P1_int2
     +    ,P2(L),P2(I),P2(K),P2_int2
     +    ,P3(L),P3(I),P3(K),P3_int2
     +    ,P4(L),P4(I),P4(K),P4_int2
     +    ,P5(L),P5(I),P5(K),P5_int2
     +    )

      Do j=NSG2,K,-1
       j1     = j+1
        X(j1) =  X(j)
        Y(j1) =  Y(j)
       P1(j1) = P1(j)
       P2(j1) = P2(j)
       P3(j1) = P3(j)
       P4(j1) = P3(j)
       P5(j1) = P3(j)
      End Do

       X(K) = XC(I)+R(I)*COS(TH23)
       Y(K) = YC(I)+R(I)*SIN(TH23)
      P1(K) = P1_int2
      P2(K) = P2_int2
      P3(K) = P3_int2
      P4(K) = P4_int2
      P5(K) = P5_int2

       X(I) = XC(I)+R(I)*COS(TH13)
       Y(I) = YC(I)+R(I)*SIN(TH13) 
      P1(I) = P1_int1
      P2(I) = P2_int1
      P3(I) = P3_int1
      P4(I) = P4_int1
      P5(I) = P5_int1

      NSG = NSG1

      write (6,206) NSG

      Go to 98

c---------------------------
c  If a segment is too long,
c  add a point in the middle
c---------------------------
 
  91  Continue

      write (6,203) I

      ANGL1 = 0.5D0*(TH2(L)+TH3(L))

      X_int1 = XC(L)+R(L)*COS(ANGL1)
      Y_int1 = YC(L)+R(L)*SIN(ANGL1)

      La = L-1

      call prd_2d_pr_arc_interp
     +
     +     (TH1(L),TH2(L),TH3(L),ANGL1
     +     ,P1(La),P1(L),P1(I),P1_int1
     +     ,P2(La),P2(L),P2(I),P2_int1
     +     ,P3(La),P3(L),P3(I),P3_int1
     +     ,P4(La),P4(L),P4(I),P4_int1
     +     ,P5(La),P5(L),P5(I),P5_int1
     +     )

      ANGL2 = 0.5D0*(TH1(I)+TH2(I))

      X_int2 = XC(I)+R(I)*COS(ANGL2)
      Y_int2 = YC(I)+R(I)*SIN(ANGL2)

      call prd_2d_pr_arc_interp
     +
     +     (TH1(I),TH2(I),TH3(I),ANGL2
     +     ,P1(L),P1(I),P1(K),P1_int2
     +     ,P2(L),P2(I),P2(K),P2_int2
     +     ,P3(L),P3(I),P3(K),P3_int2
     +     ,P4(L),P4(I),P4(K),P4_int2
     +     ,P5(L),P5(I),P5(K),P5_int2
     +     )

       X_int = 0.5D0*( X_int1 +  X_int2)
       Y_int = 0.5D0*( Y_int1 +  Y_int2)
      P1_int = 0.5D0*(P1_int1 + P1_int2)
      P2_int = 0.5D0*(P2_int1 + P2_int2)
      P3_int = 0.5D0*(P3_int1 + P3_int2)
      P4_int = 0.5D0*(P4_int1 + P4_int2)
      P5_int = 0.5D0*(P5_int1 + P5_int2)

      Do j=NSG1,I,-1
        j1 = j+1
         X(j1) =  X(j)
         Y(j1) =  Y(j)
        P1(j1) = P1(j)
        P2(j1) = P2(j)
        P3(j1) = P3(j)
        P4(j1) = P4(j)
        P5(j1) = P5(j)
      End Do

       X(I) =  X_int
       Y(I) =  Y_int
      P1(I) = P1_int
      P2(I) = P2_int
      P3(I) = P3_int
      P4(I) = P4_int
      P5(I) = P5_int

      NSG = NSG1

      write (6,206) NSG

      Go to 98

c---------------------------
c  If a segment is too long,
c  add a point in the middle
c---------------------------

  92  Continue

      write (6,204) I

      K1 = K+1

      call arc_2d 
     +
     +    (I
     +    ,X(I),X(K),X(K1)
     +    ,Y(I),Y(K),Y(K1)
     +    ,XC(K),YC(K),R(K)
     +    ,TH1(K),TH2(K),TH3(K),ORNT(K)
     +    ,AREA1(K),AREA2(K),ARE(K)
     +    ,XCN1 (K),XCN2 (K),XCN(K)
     +    ,YCN1 (K),YCN2 (K),YCN(K)
     +    ,Istop
     +    )

      ANGL1  = 0.5D0*(TH2(I)+TH3(I))

      X_INT1 = XC(I)+R(I)*COS(ANGL1)
      Y_INT1 = YC(I)+R(I)*SIN(ANGL1)

      call prd_2d_pr_arc_interp
     +
     +  (TH1(I),TH2(I),TH3(I),ANGL1
     +  ,P1(L),P1(I),P1(K),P1_int1
     +  ,P2(L),P2(I),P2(K),P2_int1
     +  ,P3(L),P3(I),P3(K),P3_int1
     +  ,P4(L),P4(I),P4(K),P4_int1
     +  ,P5(L),P5(I),P5(K),P5_int1
     +  )

      ANGL2  = 0.5D0*(TH1(K)+TH2(K))

      X_INT2 = XC(K)+R(K)*COS(ANGL2)
      Y_INT2 = YC(K)+R(K)*SIN(ANGL2)

      call prd_2d_pr_arc_interp
     +
     +    (TH1(K),TH2(K),TH3(K),ANGL2
     +    ,P1(I),P1(K),P1(K1),P1_int2
     +    ,P2(I),P2(K),P2(K1),P2_int2
     +    ,P3(I),P3(K),P3(K1),P3_int2
     +    ,P4(I),P4(K),P4(K1),P4_int2
     +    ,P5(I),P5(K),P5(K1),P5_int2
     +    )

       X_INT = 0.5D0*( X_INT1 +  X_INT2)
       Y_INT = 0.5D0*( Y_INT1 +  Y_INT2)
      P1_INT = 0.5D0*(P1_INT1 + P1_INT2)
      P2_INT = 0.5D0*(P2_INT1 + P2_INT2)
      P3_INT = 0.5D0*(P3_INT1 + P3_INT2)
      P4_INT = 0.5D0*(P4_INT1 + P4_INT2)
      P5_INT = 0.5D0*(P5_INT1 + P5_INT2)

      Do j=NSG1,K,-1
        j1 = j+1
         X(j1) =  X(j)
         Y(j1) =  Y(j)
        P1(j1) = P1(j)
        P2(j1) = P2(j)
        P3(j1) = P3(j)
        P4(j1) = P4(j)
        P5(j1) = P5(j)
      End Do

       X(K) =  X_INT
       Y(K) =  Y_INT
      P1(K) = P1_INT
      P2(K) = P2_INT
      P3(K) = P3_INT
      P4(K) = P4_INT
      P5(K) = P5_INT

      NSG = NSG1

      write (6,206) NSG

      Go to 98
 
c-------------------------------
c  If a segment is too short,
c  remove the end-poitns and add
c  a point in the middle
c-------------------------------
 
  93  Continue

      write (6,209) I

      Do j = K,NSG
        j1  = j+1
        X(j) =  X(j1)
        Y(j) =  Y(j1)
       P1(j) = P1(j1)
       P2(j) = P2(j1)
       P3(j) = P3(j1)
       P4(j) = P4(j1)
       P5(j) = P5(j1)
      End Do

      ANGL = 0.5D0*(TH2(I)+TH3(I))

      call prd_2d_pr_arc_interp
     +
     +    (TH1(I),TH2(I),TH3(I),ANGL
     +    ,P1(L),P1(I),P1(K),P1_int
     +    ,P2(L),P2(I),P2(K),P2_int
     +    ,P3(L),P3(I),P3(K),P3_int
     +    ,P4(L),P4(I),P4(K),P4_int
     +    ,P5(L),P5(I),P5(K),P5_int
     +    )

       X(I) = XTMP
       Y(I) = YTMP
      P1(I) = P1_int
      P2(I) = P2_int
      P3(I) = P3_int
      P4(I) = P4_int
      P5(I) = P5_int

      NSG  = NSG-1
      NSG1 = NSG+1

c---
c to account for removal when i = NSG
c---

      If(I.eq.NSG1) then
       X(1) = X(NSG1)-RL
       Y(1) = Y(NSG1)
      End If

      write (6,208) NSG

      Go to 98
 
c=========================================

  29  Continue

      If(Ich4.eq.1) then

c--------------------
c Redefine the period
c--------------------

      Iloop = 1

  89  Continue

      Iloop = Iloop + 1

      IF(Iloop.eq.10) Go to 88   ! prevent infinite looping

c---
      If(X(1).lt.thrs1) then
c---
      write (6,*) " Period redefined, 1"

        Do I = 0,NSG1
          I1 = I+1
           X(I) =  X(I1) 
           Y(I) =  Y(I1) 
          P1(I) = P1(I1) 
          P2(I) = P2(I1) 
          P3(I) = P3(I1) 
          P4(I) = P4(I1) 
          P5(I) = P5(I1) 
        End Do

         X(NSG2) =  X(2) + RL 
         Y(NSG2) =  Y(2) 
        P1(NSG2) = P1(2) 
        P2(NSG2) = P2(2) 
        P3(NSG2) = P3(2) 
        P4(NSG2) = P4(2) 
        P5(NSG2) = P5(2) 

        Go to 89
c---
      Else If(X(1).gt.thrs2) then
c---

      write (6,*) " Period redefined, 2"

        Do I=1,NSG1
          I2 = NSG2-I
          I1 = NSG1-I
           X(I2) =  X(I1) 
           Y(I2) =  Y(I1) 
          P1(I2) = P1(I1) 
          P2(I2) = P2(I1) 
          P3(I2) = P3(I1) 
          P4(I2) = P4(I1) 
          P5(I2) = P5(I1) 
        End Do

        X (NSG2) = X (2) + RL 
        Y (NSG2) = Y (2) 
        P1(NSG2) = P1(2) 
        P2(NSG2) = P2(2) 
        P3(NSG2) = P3(2) 
        P4(NSG2) = P4(2)+arclt
        P5(NSG2) = P5(2) 

         X(0) = X (NSG) - RL
         Y(0) = Y (NSG) 
        P1(0) = P1(NSG) 
        P2(0) = P2(NSG) 
        P3(0) = P3(NSG) 
        P4(0) = P4(NSG)-arclt
        P5(0) = P5(NSG) 

        Go to 89
c---
      End If
c---

  88  Continue

      End If

c-----------------------------
c THE SUCCESSFUL DISTRIBUTION
c-----------------------------
 
c-----------------------
c Compute the arc length
c at the marker points starting
c from the first point
c-----------------------

      S(1) = 0.0D0

      Do i=2,NSG1
       ia = i-1
       S(i) = S(ia)+0.5 *(R(i) *abs(TH2(i )-TH1(i ))
     +                   +R(ia)*abs(TH3(ia)-TH2(ia)))
      End Do

      S(NSG2) = S(NSG1)+S(2)

      ARCT = S(NSG1)

c------------------
c Area and Centroid
c------------------

      ARET = 0.0D0
      XCNT = 0.0D0
      YCNT = 0.0D0

      Do I=2,NSG1
        ia = i-1
        ARET = ARET + 0.5D0*(AREA1(i)+AREA2(ia))
        XCNT = XCNT + 0.5D0*( XCN1(i)+ XCN2(ia))
        YCNT = YCNT + 0.5D0*( YCN1(i)+ YCN2(ia))
      End Do
 
      ARET = - ARET
      XCNT = XCNT/ARET
      YCNT = YCNT/ARET

c---
c Done
c---

 100  Format (1X,F10.5,2X,F10.5,10X,F10.5,2X,F10.5)
 101  Format (1X,I3,10(1x,F10.5))

 200  Format (' ARC',I3,' IS TOO BIG',2(F10.5))
 201  Format (' EXECUTION WILL BE INTERRUPTED')
 203  Format (' SEGMENT 1 OF ARC',I3,' IS TOO LARGE')
 204  Format (' SEGMENT 2 OF ARC',I3,' IS TOO LARGE')
 206  Format (' ONE POINT INSERTED, TOTAL Points ',I3,/)
 207  Format (' POINT',I3,' CROSSED THE AXIS')
 208  Format (' ONE POINT REMOVED , TOTAL NSG ',I3,/)
 209  Format (' SEGMENT 2 OF ARC',I3,' IS TOO SMALL')

      return
      end
