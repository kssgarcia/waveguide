      subroutine prd_2d_pr_splc
     +
     +   (NSG
     +   ,X,Y
     +   ,RL
     +   ,ICH1,thmax
     +   ,ICH2,spmax
     +   ,ICH3,spmin
     +   ,ICH4,thrs1,thrs2
     +   ,P1
c    +   ,P2,P3
c    +   ,P4,P5
     +   ,Istop
     +   )

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

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
c----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  X(0:513), Y(0:513)

      Dimension vnx(0:513),vny(0:513)
      Dimension crv(0:513)
      Dimension   S(0:513)

      Dimension P1(0:513),P2(0:513),P3(0:513)
      Dimension P4(0:513),P5(0:513)

      Dimension Xint(513)

      Dimension Axint(513),Bxint(513),Cxint(513)
      Dimension Ayint(513),Byint(513),Cyint(513)

      Parameter (Nmax=513)

c-------
c common
c-------

      common/splc_xy0/Xint
      common/splc_xy1/Axint,Bxint,Cxint
      common/splc_xy2/Ayint,Byint,Cyint

c----------
c constants
c----------

      pi  = 3.14159 265358 979 32384D0
      pih = 0.5D0*pi

c--------
c prepare
c--------

      Ipass = 0

      Isafe = 100*NSG

      arclt = P4(NSG+1)-P4(1)

c--------------------------

  98  Continue

      Ipass = Ipass+1

      if(NSG.GE.Nmax) then
        write (6,*) 
        write (6,*) " prd_2d_pr_splc: more than ",Nmax," points "
        Istop = 1
        return
      end if
 
      if(Ipass.gt.Isafe) then
        write (6,*)
        write(6,*) " prd_2d_pr_splc: Caught in a loop"
        write(6,*) "                 redistributing points "
        Istop = 1
        return
      end if

c--------
c prepare
c--------

      NSGa = NSG-1
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

c---------------------------
c cubic-spline interpolation
c---------------------------

      call splc_pr_geo
     +
     +   (NSG
     +   ,X,Y
     +   ,vnx,vny
     +   ,crv
     +   ,S
     +   )

c-------
c checks
c-------

      loop1 = 1
      loop2 = NSG

      Do 3 i=loop1,loop2

c----
c prepare
c----

       La = L-1
       L  = I-1
       K  = I+1
       K1 = K+1

       Ax = Axint(I)
       Bx = Bxint(I)
       Cx = Cxint(I)

       Ay = Ayint(I)
       By = Byint(I)
       Cy = Cyint(I)

c---
c arc length of the Ith segment
c---

       sp2 = S(K)-S(I)

c---------------------------------------
c check for angles subtended by the arcs
c---------------------------------------
 
      if(Ich1.eq.1) then

        if(i.gt.1) then

        thtot = (S(K)-S(L))*crv(i)

        Dth = Dabs(thtot)

        if(Dth.GE.pih  ) Go to 96
        if(Dth.GE.thmax) Go to 97

        end if

      end if

c-----------------------------------
c Check for maximum point separation
c-----------------------------------

      if(Ich2.eq.1) then

        if(sp2.GT.spmax) Go to 92

      end if

c-----------------------------------
c Check for minimum point separation
c-----------------------------------

      if(Ich3.eq.0)    Go to 3 ! skip
      if(sp2.GT.spmin) Go to 3 ! skip
      if(i.eq.1)       Go to 3 ! skip

c-------------------------------------
c  POINTS I and I+1 WILL BE REMOVED
c  AND REPLACED BY A MID-POINT
c
c  BUT ONLY IF THE RESULTING ANGLES
c  AND POINT SEPARATIONS DO NOT EXCEED
c  THE PRE-ESTABLISHED MAXIMA
c-------------------------------------

       Stemp   = 0.50D0*(  S(I)+  S(K))
       crvtemp = 0.50D0*(crv(I)+crv(K))

       SEP1 = Stemp - S(I-1)
       SEP2 = S(I+2)- Stemp

       TOTANG0 = crv(I-1)*(Stemp -S(I-2))
       TOTANG1 = crvTEMP *(S(I+2)-S(I-1))
       TOTANG2 = crv(I+2)*(S(I+3)-Stemp)

       TOTANG0 = Dabs(TOTANG0)
       TOTANG1 = Dabs(TOTANG1)
       TOTANG2 = Dabs(TOTANG2)

       if(      TOTANG0.LT.thmax
     +     .AND.TOTANG1.LT.thmax
     +     .AND.TOTANG2.LT.thmax
     +     .AND.   SEP1.LT.spmax
     +     .AND.   SEP2.LT.spmax) Go to 93

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

      write (6,*) " prd_2d_pr_splc: arc ",I," too large",Dth,thmax
      write (6,*) "             execution aborted"

      Istop = 1

      return

c--------------------

  97  Continue

c-----------------------------------
c  if an arc is too large,
c  remove the middle point and add
c  two evenly-spaced points
c-----------------------------------

      write (6,*) " prd_2d_pr_splc: arc ",I," too large",Dth,thmax

      Ax1 = Axint(L)
      Bx1 = Bxint(L)
      Cx1 = Cxint(L)

      Ay1 = Ayint(L)
      By1 = Byint(L)
      Cy1 = Cyint(L)

      XINT1 = (Xint(L)+2.0D0*Xint(I))/3.0D0

      DX1 = XINT1-Xint(L)

      Xnew1 = ( ( Ax1*DX1 + Bx1) * DX1 + Cx1 )*DX1 + X(L)
      Ynew1 = ( ( Ay1*DX1 + By1) * DX1 + Cy1 )*DX1 + Y(L)

      call prd_2d_pr_splc_interp
     +
     +    (Xint(La),Xint(L),Xint(I),Xint(K)
     +    ,XINT1
     +    ,P1(La),P1(L),P1(I),P1(K),P1_int1
     +    ,P2(La),P2(L),P2(I),P2(K),P2_int1
     +    ,P3(La),P3(L),P3(I),P3(K),P3_int1
     +    ,P4(La),P4(L),P4(I),P4(K),P4_int1
     +    ,P5(La),P5(L),P5(I),P5(K),P5_int1
     +    )

      Ax2 = Axint(I)
      Bx2 = Bxint(I)
      Cx2 = Cxint(I)

      Ay2 = Ayint(I)
      By2 = Byint(I)
      Cy2 = Cyint(I)

      XINT2 = (2.0D0*Xint(I)+Xint(K))/3.0D0

      DX2 = XINT2-Xint(I)

      Xnew2 = ( ( Ax2*DX2 + Bx2) * DX2 + Cx2 )*DX2 + X(I)
      Ynew2 = ( ( Ay2*DX2 + By2) * DX2 + Cy2 )*DX2 + Y(I)

      call prd_2d_pr_splc_interp
     +
     +    (Xint(L),Xint(I),Xint(K),Xint(K1)
     +    ,XINT2
     +    ,P1(L),P1(I),P1(K),P1(K1),P1_int2
     +    ,P2(L),P2(I),P2(K),P2(K1),P2_int2
     +    ,P3(L),P3(I),P3(K),P3(K1),P3_int2
     +    ,P4(L),P4(I),P4(K),P4(K1),P4_int2
     +    ,P5(L),P5(I),P5(K),P5(K1),P5_int2
     +    )

      Do j=NSG2,K,-1
       j1 = j+1
        X(j1) =  X(j)
        Y(j1) =  Y(j)
       P1(j1) = P1(j)
       P2(j1) = P2(j)
       P3(j1) = P3(j)
       P4(j1) = P3(j)
       P5(j1) = P3(j)
      End Do

       X(K) = Xnew2
       Y(K) = Ynew2
      P1(K) = P1_int2
      P2(K) = P2_int2
      P3(K) = P3_int2
      P4(K) = P4_int2
      P5(K) = P5_int2

       X(I) = Xnew1
       Y(I) = Ynew1
      P1(I) = P1_int1
      P2(I) = P2_int1
      P3(I) = P3_int1
      P4(I) = P4_int1
      P5(I) = P5_int1

      NSG = NSG1

      write (6,*) "  One point inserted; segments : ",NSG

      Go to 98

c---------------------------
c  if a segment is too long,
c  add a point in the middle
c---------------------------

  92  Continue

      write (6,*) ' segment ',I,' too long'

      XINT1 = 0.50D0*(Xint(I)+Xint(K))

      DX = XINT1-Xint(I)

      Xnew = ( ( Ax*DX + Bx) * DX + Cx ) *DX + X(I)
      Ynew = ( ( Ay*DX + By) * DX + Cy ) *DX + Y(I)

      if(I.eq.1) then
        XintL = Xint(NSGa)-RL
        P1L = P1(NSGa)
        P2L = P2(NSGa)
        P3L = P3(NSGa)
        P4L = P4(NSGa)
        P5L = P5(NSGa)
      else
        Xint(L) = Xint(L)
        P1L = P1(L)
        P2L = P2(L)
        P3L = P3(L)
        P4L = P4(L)
        P5L = P5(L)
      end if

      call prd_2d_pr_splc_interp
     +
     +    (XintL,Xint(I),Xint(K),Xint(K1)
     +    ,XINT1
     +    ,P1L,P1(I),P1(K),P1(K1),P1_int
     +    ,P2L,P2(I),P2(K),P2(K1),P2_int
     +    ,P3L,P3(I),P3(K),P3(K1),P3_int
     +    ,P4L,P4(I),P4(K),P4(K1),P4_int
     +    ,P5L,P5(I),P5(K),P5(K1),P5_int
     +    )

      Do j=NSG1,K,-1
        j1    = j+1
         X(j1) =  X(j)
         Y(j1) =  Y(j)
        P1(j1) = P1(j)
        P2(j1) = P2(j)
        P3(j1) = P3(j)
        P4(j1) = P4(j)
        P5(j1) = P5(j)
      End Do

       X(K) =  Xnew
       Y(K) =  Ynew
      P1(K) = P1_INT
      P2(K) = P2_INT
      P3(K) = P3_INT
      P4(K) = P4_INT
      P5(K) = P5_INT

      NSG = NSG1

      NSG1 = NSG+1

      write (6,*) " one point inserted: segments: ",NSG

c     Do j=1,NSG1
c       write (6,101) j,X(j),Y(j),P1(j)
c     End Do
c     pause

      Go to 98
 
c-------------------------------
c  if a segment is too short,
c  remove the end-points and add
c  a point in the middle
c-------------------------------
 
  93  Continue

      write (6,*) ' prd_2d_pr1: segment ',I,' is too small', sp2,spmin

      Do j=K,NSG
          j1 = j+1
        X(j) =  X(j1)
        Y(j) =  Y(j1)
       P1(j) = P1(j1)
       P2(j) = P2(j1)
       P3(j) = P3(j1)
       P4(j) = P4(j1)
       P5(j) = P5(j1)
      End Do

      Xint1 = 0.50D0*(Xint(I)+Xint(K))

      DX = Xint1 - Xint(I)

      XNEW = ( ( Ax*DX + Bx) * DX + Cx )*DX + X(I)
      YNEW = ( ( Ay*DX + By) * DX + Cy )*DX + Y(I)

      call prd_2d_pr_splc_interp
     +
     +    (Xint(L),Xint(I),Xint(K),Xint(K1)
     +    ,Xint1
     +    ,P1(L),P1(I),P1(K),P1(K1),P1_int
     +    ,P2(L),P2(I),P2(K),P2(K1),P2_int
     +    ,P3(L),P3(I),P3(K),P3(K1),P3_int
     +    ,P4(L),P4(I),P4(K),P4(K1),P4_int
     +    ,P5(L),P5(I),P5(K),P5(K1),P5_int
     +    )

       X(I) = XNEW
       Y(I) = YNEW
      P1(I) = P1_int
      P2(I) = P2_int
      P3(I) = P3_int
      P4(I) = P4_int
      P5(I) = P5_int

      NSG  = NSG-1
      NSG1 = NSG+1

c---
c to account for removal when i = nsg
c---

      if(I.eq.NSG1) THEN
       X(1) = X(NSG1)-RL
       Y(1) = Y(NSG1)
      end if

      write (6,208) NSG

      Go to 98
 
c=========================================

  29  Continue

      if(Ich4.eq.1) then

c     write (6,*) " prd_2d_pr1: check 4:",thrs1,thrs2,X(1),X(NSG1)
      
c--------------------
c Redefine the period
c--------------------

      Iloop = 1

  89  Continue

      Iloop = Iloop + 1

      if(Iloop.eq.10) Go to 88   ! escape to prevent infinite looping

c---
      if(X(1).lt.thrs1) then
c---

      write (6,*) " prd_2d_pr1: period redefined from the left"

        Do I=0,NSG1
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
      else if(X(NSG1).gt.thrs2) then
c---

      write (6,*) " prd_2d_pr_splc: period redefined from the right"

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

         X(NSG2) =  X(2) + RL 
         Y(NSG2) =  Y(2) 
        P1(NSG2) = P1(2) 
        P2(NSG2) = P2(2) 
        P3(NSG2) = P3(2) 
        P4(NSG2) = P4(2)+arclt
        P5(NSG2) = P5(2) 

         X(0) =  X(NSG) - RL
         Y(0) =  Y(NSG) 
        P1(0) = P1(NSG) 
        P2(0) = P2(NSG) 
        P3(0) = P3(NSG) 
        P4(0) = P4(NSG)-arclt
        P5(0) = P5(NSG) 

        Go to 89
c---
      end if
c---

  88  Continue

c-----------
      end if
c-----------

c-----
c done
c-----

 100  Format (1X,F10.5,2X,F10.5,10X,F10.5,2X,F10.5)
 101  Format (1X,I3,10(1x,F10.5))

 200  Format (' ARC',I3,' IS TOO BIG',2(F10.5))
 201  Format (' EXECUTION WILL BE INTERRUPTED')
 203  Format (' SEGMENT 1 OF ARC',I3,' IS TOO LARGE')
 207  Format (' POINT',I3,' CROSSED THE AXIS')
 208  Format (' ONE POINT REMOVED , TOTAL NSG ',I3,/)

      Return
      End
