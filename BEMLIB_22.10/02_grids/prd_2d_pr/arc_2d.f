      subroutine arc_2d 
     +
     +   (Iarc
     +   ,X1,X2,X3
     +   ,Y1,Y2,Y3
     +   ,XC,YC,R
     +   ,TH1,TH2,TH3,ORNT
     +   ,AREA1,AREA2,AREA
     +   ,XCN1,XCN2,XCN
     +   ,YCN1,YCN2,YCN
     +   ,Istop
     +   )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c  Properties of the circular arc passing through
c  three points numbered 1,2,3 in the xy plane
c
c  SYMBOLS:
c  -------
c
c  (xc, yc):  coordinates of arc center
c
c  th1,th2,th3:	angles subtended by the arc center and
c		and each one of the three points
c
c  ornt:	arc orientation: 1 for counterclockwise
c                               -1 for clockwise
c
c  area1	area between the points 12 and the x-axis
c  area2	area between the points 23 and the x-axis
c  area         area1+area2
c
c  (xcn1,ycn1)	areal centroid integral of first part
c               of the arc
c  (xcn2,ycn2)	areal centroid integral of second part
c               of the arc
c  (xcn,ycn)	areal centroid integral of the whole arc
c
c  eps:		tolerance for identifying the three points
c		as colinear
c  eps1:	tolerance for the computation of the arc center
c------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (eps=0.000001D0,eps1=0.000001D0)
      Parameter (tol=0.000001D0)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi

      epsm = -eps

c---------------------------
c deduce the arc orientation
c---------------------------
 
      outp = (X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1) ! outer product
      D1   = Dsqrt((X2-X1)**2+(Y2-Y1)**2)
      D2   = Dsqrt((X3-X1)**2+(Y3-Y1)**2)
      ORNT = outp/(D1*D2)

      if(ORNT.gt.eps) then
        ORNT = 1.0D0
      else if(ORNT.lt.epsm) then
        ORNT =-1.0D0
      else
        write (6,*) 
        write (6,*) ' arc_2d: arc ',Iarc,' is not defined'
        write (6,*) 
        Istop = 1
        Go to 99
      end if

c-----------------------
c Compute the arc center
c-----------------------

      if(abs(X2-X1).lt.eps1) then
       AA1 = 0.0D0
      else 
       A1 = (Y2-Y1)/(X2-X1)
       If(abs(A1).le.eps1) A1 = eps1
       AA1 = -1.0D0/A1
      end if

      if(abs(X3-X2).lt.eps1) then
       AA2 = 0.0D0
      else
       A2 = (Y3-Y2)/(X3-X2)
       If(abs(A2).le.eps1) A2 = eps1
       AA2 = -1.0D0/A2
      end if

      XX1 = 0.5D0*(X1+X2)   ! mid-points
      YY1 = 0.5D0*(Y1+Y2)
      XX2 = 0.5D0*(X2+X3)
      YY2 = 0.5D0*(Y2+Y3)

      GG1 = YY1-AA1*XX1
      GG2 = YY2-AA2*XX2

      XC = (GG1-GG2)/(AA2-AA1)
      YC = AA1*XC+GG1

c-----------
c arc radius
c-----------

      R2 = (X1-XC)*(X1-XC)+(Y1-YC)*(Y1-YC)
      R  = Dsqrt(R2)

c-------------------------------------
c angles subtended by the three points
c from the arc center
c-------------------------------------

      ARG = (X1-XC)/R

      If(ARG.GT.1.0D0) then
        write (6,*)
        write (6,300) ARG
        write (6,*)
        ARG = 0.9999999999999D0
      Else If(ARG.LT.-1.0D0) then
        write (6,*)
        write (6,*) ARG
        write (6,*)
        ARG = -0.9999999999999D0
      End If

      TH1  = ACOS(ARG)
      SINA = Y1-YC

      If(SINA.LT.-0.000000001) TH1 = pi2-TH1

      PROD2 = (X1-XC)*(X2-XC)+(Y1-YC)*(Y2-YC)
      PROD3 = (X1-XC)*(X3-XC)+(Y1-YC)*(Y3-YC)

      ANG2 = PROD2/R2
      ANG3 = PROD3/R2

      if(ANG2.GT.1.0D0) then
        write (6,*)
        write (6,300) ANG2
        write (6,*)
        ANG2 = 0.9999999999999D0
      else if(ANG2.LT.-1.0D0) then
        write (6,*)
        write (6,*) ANG2
        write (6,*)
        ANG2 = -0.9999999999999D0
      end if

      if(ANG3.GT.1.0D0) then
        write (6,*)
        write (6,300) ANG3
        write (6,*)
        ANG3 = 0.9999999999999D0
      else if(ANG3.LT.-1.0D0) then
        write (6,*)
        write (6,*) ANG3
        write (6,*)
        ANG3 = -0.9999999999999D0
      end if

      DTH2 = ORNT*ACOS(ANG2)
      DTH3 = ORNT*ACOS(ANG3)

      TH2  = TH1+DTH2
      TH3  = TH1+DTH3

c----------------------
c Compute the integrals:
c
c Int(y*n_y dl) = Int( y dx) = Int(dA)
c
c Int(xy dx) = Int(x dA)
c
c Int(yy dx) = Int(y dA)
c
c----------------------

      COS1 = DCOS(TH1)
      COS2 = DCOS(TH2)
      COS3 = DCOS(TH3)
      SIN1 = DSIN(TH1)

      SIN2 = DSIN(TH2)
      SIN3 = DSIN(TH3)

      CS1C = COS1**3
      CS2C = COS2**3
      CS3C = COS3**3

      SN1C = SIN1**3
      SN2C = SIN2**3
      SN3C = SIN3**3

      TH12 = 2.0D0*TH1
      TH22 = 2.0D0*TH2
      TH32 = 2.0D0*TH3

c---
c first part of the arc
c---

      DTH   = TH2       -TH1
      DSIN2 = SIN(TH22) -SIN(TH12)
      DCS   = COS2      -COS1
      DCS3  = CS2C      -CS1C
      DSN   = SIN2      -SIN1
      DSN3  = SN2C      -SN1C

      AREA1 = -R*(YC*DCS - 0.5D0*R*(DTH-0.5D0*DSIN2))

      XCN1 = 0.5D0*R*( XC**2*DSN
     +                +R*XC*(DTH+0.5D0*DSIN2)
     +                +R**2*(DSN-DSN3/3.0D0) )
      YCN1 = 0.5D0*R*(-YC**2*DCS
     +                +R*YC*(DTH-0.5D0*DSIN2)
     +                -R**2*(DCS-DCS3/3.0D0) )

c---
c second part of the arc
c---

      DTH   = TH3       -TH2
      DSIN2 = SIN(TH32) -SIN(TH22)
      DCS   = COS3      -COS2
      DCS3  = CS3C      -CS2C
      DSN   = SIN3      -SIN2
      DSN3  = SN3C      -SN2C

      AREA2 = -R*( YC*DCS - 0.5D0*R*(DTH-0.5D0*DSIN2))

      XCN2 = 0.5D0*R*(XC**2 * DSN
     +             +2.0D0*R*XC*(0.5D0*DTH+0.25D0*DSIN2)
     +             +R**2    *(DSN-DSN3/3.0D0) )

      YCN2 = 0.5*R*(-YC**2 * DCS
     +             +R*YC*(DTH-0.5D0*DSIN2)
     +             -R**2*(DCS-DCS3/3.0D0) )

c---
c whole of the arc
c---

      AREA  = AREA1 + AREA2
      XCN   = XCN1  + XCN2
      YCN   = YCN1  + YCN2

c----------------------
c check for consistency
c----------------------

      XX1 = XC+R*COS(TH1)
      YY1 = YC+R*SIN(TH1)

      XX2 = XC+R*COS(TH2)
      YY2 = YC+R*SIN(TH2)

      XX3 = XC+R*COS(TH3)
      YY3 = YC+R*SIN(TH3)

      D1 = XX1-X1       ! error
      D2 = XX2-X2
      D3 = XX3-X3

      F1 = YY1-Y1       ! error
      F2 = YY2-Y2
      F3 = YY3-Y3

      If(abs(D1).GE.tol.or
     +  .abs(F1).GE.tol) then
        write (6,*)
        write (6,101) Iarc,X1,XX1,Y1,YY1
        write (6,*)
        Istop = 1
        Go to 99
      End If

      if(abs(D2).GE.tol.or
     +  .abs(F2).GE.tol) then
        write (6,*)
        write (6,102) Iarc,X2,XX2,Y2,YY2
        write (6,*)
        Istop = 1
        Go to 99
      end if

      if(abs(D3).GE.tol.or
     +  .abs(F3).GE.tol) then
        write (6,*)
        write (6,103) Iarc,X3,XX3,Y3,YY3
        write (6,*)
        Istop = 1
      End If

c-----
c Done
c-----

  99  Continue

 100  Format (10(1X,F15.10))
 101  Format (' Problem with point 1 of arc ',I3,5(1x,F13.9))
 102  Format (' Problem with point 2 of arc ',I3,5(1x,F13.9))
 103  Format (' Problem with point 3 of arc ',I3,5(1x,F13.9))

 300  Format (' WARNING IN ARC',/,
     +        ' POINT 1 IS ALIGNED WITH THE X AXIS',/,
     +        ' DIRECTION COSINE WAS FOUND TO BE ',F15.12,/,
     +        ' WILL BE SET EQUAL TO +1 OR -1')

      Return
      End
