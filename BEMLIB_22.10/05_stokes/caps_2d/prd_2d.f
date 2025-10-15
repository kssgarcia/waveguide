      subroutine prd_2d
     +
     +  (NSG
     +  ,ICH1,THMAX
     +  ,ICH2,SPMAX
     +  ,ICH3,SPMIN
     +  ,ARCT,ARET
     +  ,XCNT,YCNT
     +  ,P1,P2,P3,P4,P5
     +  ,Istop
     +  )

c=====================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=====================================

c------------------------------------------------
c  Parametrize a CLOSED line defined by NSG+1 points
c  defining NSG segments
c
c  Point 1 coincides with with point NSG+1
c
c  Checks the primary point distribution for:
c
c  the maximum angle subtended by the arcs
c  the maximum point separation
c  the minimum point separation
c
c  and returns an improved distribution
c
c  It also interpolates for the variables P1-P5
c
c      Property   Example
c
c            P1:  x-vel
c            P2:  y-vel
c            P3:  surface tension
c            P4:  arc length 
c            P5:  unused
c
c  SYMBOLS:
c  --------
c
c  X,Y:   coordinates of the marker points
c
c  ICh1:  check 1 flag
c  ICh2:  check 2 flag
c  ICh3:  check 3 flag
c
c  arct:  total arc length
c
c  Istop:  stopping flag; set Iflag = 1 when stop
c----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   X(0:900), Y(0:900)
      Dimension  P1(0:900),P2(0:900),P3(0:900)
      Dimension  P4(0:900),P5(0:900)

      Dimension    XC(900),   YC(900),  R(900),   S(900)
      Dimension   TH1(900),  TH2(900),TH3(900),ORNT(900)
      Dimension AREA1(900),AREA2(900),ARE(900)
      Dimension  XCN1(900), XCN2(900),XCN(900)
      Dimension  YCN1(900), YCN2(900),YCN(900)

      Parameter (Nmax=900)

c--------------
c common blocks
c--------------

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT

c----------
c constants
c----------

      pi  = 3.14159 265358 979323 84 D0
      pih = 0.5D0*pi

c--------
c prepare
c--------

      Isafe = 10*NSG      ! safety valve

      thtol = 0.99*thmax  ! adjusted for tolerance

c-----------------
c total arc length
c-----------------

      arclt = P4(NSG+1)-P4(1)

c-------
c launch
c-------

      Ipass = 0

  98  Continue

      Ipass = Ipass+1

      if(NSG.GE.Nmax) then
        write (6,*) 
        write (6,*) " prd_2d_cl: more than ",Nmax," points "
        write (6,*) "            execution interrupted"
        Istop = 1
        return
      end if

      if(Ipass.gt.Isafe) then
        write (6,*) 
        write (6,*) " prd_2d_cl: caught in a loop redistributing"
        write (6,*) "            points; execution interrupted"
        Istop = 1
        return
      end if

c--------
c prepare
c--------

      NSG1 = NSG+1
      NSG2 = NSG+2
      NSG3 = NSG+3

c------------
c wrap around
c------------

       X(0) =  X(NSG)
       Y(0) =  Y(NSG)
      P1(0) = P1(NSG)
      P2(0) = P2(NSG)
      P3(0) = P3(NSG)
      P4(0) = P4(NSG)-arclt
      P5(0) = P5(NSG)

       X(NSG1) =  X(1)
       Y(NSG1) =  Y(1)
      P1(NSG1) = P1(1)
      P2(NSG1) = P2(1)
      P3(NSG1) = P3(1)
      P4(NSG1) = P4(1)+arclt
      P5(NSG1) = P5(1)

       X(NSG2) =  X(2)
       Y(NSG2) =  Y(2)
      P1(NSG2) = P1(2)
      P2(NSG2) = P2(2)
      P3(NSG2) = P3(2)
      P4(NSG2) = P4(2)+arclt
      P5(NSG2) = P5(2)

       X(NSG3) =  X(3)
       Y(NSG3) =  Y(3)
      P1(NSG3) = P1(3)
      P2(NSG3) = P2(3)
      P3(NSG3) = P3(3)
      P4(NSG3) = P4(3)+arclt
      P5(NSG3) = P5(3)

c----------
c set flags
c----------

      Istop = 0
      Iloop = 1

c----------------
c DEFINE THE ARCS
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
     +   (I
     +   ,XX1,XX2,XX3
     +   ,YY1,YY2,YY3
     +   ,XC(I),YC(I),R(I)
     +   ,TH1(I),TH2(I),TH3(I),ORNT(I)
     +   ,AREA1(I),AREA2(I),ARE(I)
     +   , XCN1(I),XCN2(I),XCN(I)
     +   , YCN1(I),YCN2(I),YCN(I)
     +   ,Istop
     +   )

      if(Istop.eq.1) return

      if(I.eq.NSG1) Go to 3

c-----------------------------------
c CHECK ANGLES SUBTENDED BY THE ARCS
c-----------------------------------
 
      if(Ich1.eq.1) then

        THETOT = TH3(i)-TH1(i)
        Dth    = ABS(THETOT)

        If(Dth.GE.pih  ) Go to 96
        If(Dth.GE.THMAX) Go to 97

      end if

c-------------------------------
c CHECK MAXIMUM POINT SEPARATION
c-------------------------------

      if(Ich2.eq.1) then

        SP1 = R(I)*ABS(TH2(I)-TH1(I))
        If(SP1.GT.SPMAX.and.I.ne.1) Go to 91

        SP2 = R(I)*ABS(TH3(I)-TH2(I))
        If(SP2.GT.SPMAX)            Go to 92

      end if

c-----------------------------------
c CHECK FOR MINIMUM POINT SEPARATION
c-----------------------------------

      if(Ich3.eq.0) GO TO 3

      SP2 = R(i)*ABS(TH3(i)-TH2(i))
      if(SP2.GT.SPMIN) GO TO 3

c---------------------------------
c  A point will be removed only of
c  the resulting angles and point
c  separations do not violate preset
c  thresholds
c---------------------------------
 
      ANGL = 0.5D0*(TH2(I)+TH3(I))
      XTMP = XC(I) + R(I) * COS(ANGL)
      YTMP = YC(I) + R(I) * SIN(ANGL)

      if(I.NE.1) then
        XXX = X(M)
        YYY = Y(M)
      else
        XXX = X(NSG-1)
        YYY = Y(NSG-1)
      end if

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

      if(Istop.eq.1) Return

      TOTANG0 = ABS(TH3T-TH1T)

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

      If(Istop.eq.1) Return

      TOTANG1 = ABS(TH3T-TH1T)
 
      SEP1 = RT * ABS(TH2T-TH1T)
      SEP2 = RT * ABS(TH3T-TH2T)

       if(I.ne.NSG) then
         XXX = X(I+3)
         YYY = Y(I+3)
       else
         XXX = X(3)
         YYY = Y(3)
       end if

       call arc_2d 
     +
     +    (I
     +    ,XTMP,X(I+2),XXX
     +    ,YTMP,Y(I+2),YYY
     +    ,XCT,YCT,RT
     +    ,TH1T,TH2T,TH3T,ORNTT
     +    ,AREA1T,AREA2T,ARET
     +    ,XCN1T,XCN2T,XCNT
     +    ,YCN1T,YCN2T,YCNT
     +    ,Istop
     +    )

       if(Istop.eq.1) return

       TOTANG2 = ABS(TH3T-TH1T)

      if( TOTANG0.LT.THtol
     +    .and.TOTANG1.LT.THtol
     +    .and.TOTANG2.LT.THtol
     +    .and.SEP1.LT.SPMAX
     +    .and.SEP2.LT.SPMAX)  Go to 93

  3   Continue

c----------------------------------
c  WE NOW HAVE 
c  AN ACCEPTABLE POINT DISTRIBUTION
c----------------------------------

      Go to 29      ! proceed to finish up

c--------------------
c  There is a problem
c  execution will stop
c--------------------

  96  Continue

      write (6,200) I,THETOT,THMAX
      write (6,201)

      write (1,200) I,THETOT,THMAX
      write (1,200) I
      write (1,201)

      Istop = 1
      Return

c----

  97  Continue

c-----------------------------------
c IF AN ARC IS TOO LARGE, 
c REMOVE THE MIDDLE POINT
c AND ADD TWO POINTS ALONG THE ARC
c-----------------------------------

      write (6,200) I,THETOT,THMAX

      Dth  = THETOT/3.0

      TH13 = TH1(i)+    Dth
      TH23 = TH1(i)+2.0*Dth

      call prd_2d_interp
     +
     +   (TH1(I),TH2(I),TH3(I),TH13
     +   ,P1(L),P1(I),P1(K),P1_int1
     +   ,P2(L),P2(I),P2(K),P2_int1
     +   ,P3(L),P3(I),P3(K),P3_int1
     +   ,P4(L),P4(I),P4(K),P4_int1
     +   ,P5(L),P5(I),P5(K),P5_int1
     +   )

      call prd_2d_interp
     +
     +   (TH1(I),TH2(I),TH3(I),TH23
     +   ,P1(L),P1(I),P1(K),P1_int2
     +   ,P2(L),P2(I),P2(K),P2_int2
     +   ,P3(L),P3(I),P3(K),P3_int2
     +   ,P4(L),P4(I),P4(K),P4_int2
     +   ,P5(L),P5(I),P5(K),P5_int2
     +   )

      Do j=NSG,K,-1
        j1 = j+1
         X(j1) =  X(j)
         Y(j1) =  Y(j)
        P1(j1) = P1(j)
        P2(j1) = P2(j)
        P3(j1) = P3(j)
        P4(j1) = P4(j)
        P5(j1) = P5(j)
      End Do

       X(K) = XC(I)+R(i)*DCOS(TH23)
       Y(K) = YC(I)+R(i)*DSIN(TH23)
      P1(K) = P1_int2
      P2(K) = P2_int2
      P3(K) = P3_int2
      P4(K) = P4_int2
      P5(K) = P5_int2

       X(I) = XC(i)+R(i)*DCOS(TH13)
       Y(I) = YC(i)+R(i)*DSIN(TH13) 
      P1(I) = P1_int1
      P2(I) = P2_int1
      P3(I) = P3_int1
      P4(I) = P4_int1
      P5(I) = P5_int1

      NSG = NSG+1

      write (6,206) NSG

      Go to 98

c---------------------------
c IF A SEGMENT IS TOO LONG,
c ADD A POINT IN THE MIDDLE
c---------------------------
 
  91  Continue

      write (6,203) I

      ANGL1 = 0.5D0*(TH2(L)+TH3(L))
      X_int1 = XC(L)+R(L)*DCOS(ANGL1)
      Y_int1 = YC(L)+R(L)*DSIN(ANGL1)

      La = L-1

      call prd_2d_interp
     +
     +    (TH1(L),TH2(L),TH3(L),ANGL1
     +    ,P1(La),P1(L),P1(I),P1_int1
     +    ,P2(La),P2(L),P2(I),P2_int1
     +    ,P3(La),P3(L),P3(I),P3_int1
     +    ,P4(La),P4(L),P4(I),P4_int1
     +    ,P5(La),P5(L),P5(I),P5_int1
     +    )

      ANGL2 = 0.5*(TH1(I)+TH2(I))
      X_int2 = XC(I)+R(I)*DCOS(ANGL2)
      Y_int2 = YC(I)+R(I)*DSIN(ANGL2)

      call prd_2d_interp
     +
     +   (TH1(I),TH2(I),TH3(I),ANGL2
     +   ,P1(L),P1(I),P1(K),P1_int2
     +   ,P2(L),P2(I),P2(K),P2_int2
     +   ,P3(L),P3(I),P3(K),P3_int2
     +   ,P4(L),P4(I),P4(K),P4_int2
     +   ,P5(L),P5(I),P5(K),P5_int2
     +   )

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

       X(I) =  X_INT
       Y(I) =  Y_INT
      P1(I) = P1_INT
      P2(I) = P2_INT
      P3(I) = P3_INT
      P4(I) = P4_INT
      P5(I) = P5_INT

      NSG = NSG+1

      write (6,206) NSG

      Go to 98

c------

  92  Continue

      write (6,204) I

      K1 = K+1

      call arc_2d 
     +
     +   (I
     +   ,X(I),X(K),X(K1)
     +   ,Y(I),Y(K),Y(K1)
     +   ,XC(K),YC(K),R(K)
     +   ,TH1(K),TH2(K),TH3(K),ORNT(K)
     +   ,AREA1(K),AREA2(K),ARE(K)
     +   , XCN1(K), XCN2(K),XCN(K)
     +   , YCN1(K), YCN2(K),YCN(K)
     +   ,Istop
     +   )

      if(Istop.eq.1) Return

      ANGL1  = 0.5D0*(TH2(I)+TH3(I))
      X_INT1 = XC(I)+R(I)*COS(ANGL1)
      Y_INT1 = YC(I)+R(I)*SIN(ANGL1)

      call prd_2d_interp
     +
     +   (TH1(I),TH2(I),TH3(I),ANGL1
     +   ,P1(L),P1(I),P1(K),P1_int1
     +   ,P2(L),P2(I),P2(K),P2_int1
     +   ,P3(L),P3(I),P3(K),P3_int1
     +   ,P4(L),P4(I),P4(K),P4_int1
     +   ,P5(L),P5(I),P5(K),P5_int1
     +   )

      ANGL2 = 0.5D0*(TH1(K)+TH2(K))
      X_INT2 = XC(K)+R(K)*DCOS(ANGL2)
      Y_INT2 = YC(K)+R(K)*DSIN(ANGL2)

      call prd_2d_interp
     +
     +   (TH1(K),TH2(K),TH3(K),ANGL2
     +   ,P1(I),P1(K),P1(K1),P1_int2
     +   ,P2(I),P2(K),P2(K1),P2_int2
     +   ,P3(I),P3(K),P3(K1),P3_int2
     +   ,P4(I),P4(K),P4(K1),P4_int2
     +   ,P5(I),P5(K),P5(K1),P5_int2
     +   )

       X_INT = 0.5D0*( X_INT1 +  X_INT2)
       Y_INT = 0.5D0*( Y_INT1 +  Y_INT2)
      P1_INT = 0.5D0*(P1_INT1 + P1_INT2)
      P2_INT = 0.5*(P2_INT1 + P2_INT2)
      P3_INT = 0.5*(P3_INT1 + P3_INT2)
      P4_INT = 0.5*(P4_INT1 + P4_INT2)
      P5_INT = 0.5*(P5_INT1 + P5_INT2)

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

      NSG = NSG+1

      write (6,206) NSG

      Go to 98
 
c-------------------------------
c IF A SEGMENT IS TOO SHORT,
c REMOVE END-POINTS
c AND ADD A POINT IN THE MIDDLE
c-------------------------------
 
  93  Continue

      write (6,209) I

      Do j=K,NSG
       j1   = j+1
        X(j) =  X(j1)
        Y(j) =  Y(j1)
       P1(j) = P1(j1)
       P2(j) = P2(j1)
       P3(j) = P3(j1)
       P4(j) = P4(j1)
       P5(j) = P5(j1)
      End Do

      X(I) = XTMP
      Y(I) = YTMP

      ANGL = 0.5D0*(TH2(I)+TH3(I))

      call prd_2d_interp
     +
     +   (TH1(I),TH2(I),TH3(I),ANGL
     +   ,P1(L),P1(I),P1(K),P1_int
     +   ,P2(L),P2(I),P2(K),P2_int
     +   ,P3(L),P3(I),P3(K),P3_int
     +   ,P4(L),P4(I),P4(K),P4_int
     +   ,P5(L),P5(I),P5(K),P5_int
     +   )

      P1(I) = P1_int
      P2(I) = P2_int
      P3(I) = P3_int
      P4(I) = P4_int
      P5(I) = P5_int

c---
c to account for removal when i = nsg
c---

      if(I.eq.NSG) then
       X(1) = X(NSG)
       Y(1) = Y(NSG)
      end if

      NSG = NSG-1

      write (6,208) NSG

      Go to 98

c-------------------------------------------

  29  Continue

c-----------------------------
c THE SUCCESSFUL DISTRIBUTION
c-----------------------------

c-----------------------
c Compute the arc length
c-----------------------

      S(1) = 0.0D0

      Do i=2,NSG
       ia=i-1
       S(i) = S(ia)+0.5*(R(i) *ABS(TH2(i )-TH1(i) )
     +                  +R(ia)*ABS(TH3(ia)-TH2(ia)))
      End Do

      i = NSG1

      S(i) = S(NSG)+0.5D0*(R(i)*ABS(TH2(i)-TH1(i))
     +                    +R(1)*ABS(TH3(1)-TH2(1)))
      S(NSG2) = S(NSG1)+S(2)

c--------------------------
c Compute Area and Centroid
c--------------------------

      ARET = 0.0D0
      XCNT = 0.0D0
      YCNT = 0.0D0

      Do i=2,NSG1
       ARET = ARET + AREA1(i)
       XCNT = XCNT + XCN1 (i)
       YCNT = YCNT + YCN1 (i)
      End Do

      ARCT = S(NSG1)
      XCNT = XCNT/ARET
      YCNT = YCNT/ARET

c-----
c Done
c-----

 100  Format (1X,F10.5,2X,F10.5,10X,F10.5,2X,F10.5)
 101  Format (1x,i3,10(1X,F10.5))
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
