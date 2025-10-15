       program prd_2d_dr

c===========================================
c FDLIB, BELMIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c-------------------------------------
c  Driver for subroutine: prd_2d
c
c  Line parametrization and
c  point redistribution around a closed line
c  defined by NSG+1 points
c
c  The redistribution is based on three criteria:
c
c  1: maximum angle subtended by arcs defined
c     by three consecutive points
c
c  2: maximum point separation
c
c  3: minimum pooint separation
c
c  SYMBOLS
c  -------
c
c  X, Y:	coordinates of the marker points
c  P1-P5:	various properties defined along
c		the line
c  XC, YC:	arc centers
c  R:		arc radii
c  ORNT:	arc orientation: 1 or -1
c
c  thmax:	maximum angle subtended by the arcs
c  spmax:	maximum point separation
c  spmin:	minimum point separation
c
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   X(0:900),Y(0:900)
      Dimension  P1(0:900),P2(0:900)
      Dimension  P3(0:900)
      Dimension  P4(0:900),P5(0:900)
      Dimension  XC(900), YC(900),  R(900),  S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

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
      pi2 = 2.0D0*pi

      Null = 0

c------------
c preferences
c------------

 98   Continue

      open (1,file="prd_2d.dat")
       read (1,*) thmax
       read (1,*) spmax
       read (1,*) spmin
      close (1)

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to generate the primary distribution '
      write (6,*) ' 2 to read the point distribution'
      write (6,*) '   from file: prd_2d.inp'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read (5,*) Iread

      if(Iread.eq.0) Go to 99

c--------
c prepare 
c--------

      thmax = thmax * pi

      open (2,file="prd_2d.out")

c----------------------------------
c generate the primary distribution
c----------------------------------

      if(Iread.eq.1) then

       write (6,*)
       write (6,*) ' Will generate an elliptical shape'
       write (6,*)
       write (6,*) ' Enter the coordinates of the center'
       write (6,*) ' of the ellipse'
       write (6,*) ' --------------'
       read  (5,*) cntrx,cntry

       write (6,*) ' Enter the x semi-axis'
       write (6,*) ' ---------------------'
       read  (5,*) a

       write (6,*)
       write (6,*) ' Enter the y semi-axis'
       write (6,*) ' ---------------------'
       read  (5,*) b

       write (6,*)
       write (6,*) ' Enter the number of segments'
       write (6,*) ' ----------------------------'
       read  (5,*) NSG

       NSG1 = NSG+1
       NSG2 = NSG+2

       Dth = pi2/(NSG1-1.0D0)

       Do i=0,NSG2
        thet = (i-1.0)*Dth
         X(i) = cntrx + a * cos(thet)
         Y(i) = cntry + b * sin(thet)
         S(i) = a*thet
        P1(i) = cos(thet)**3         ! example
        P2(i) = sin(thet)**3         ! example
        P3(i) = 0.5D0*sin(thet)**2   ! example
        P4(i) = S(i)                 ! arc length like
        P5(i) = 0.10*cos(thet)+0.30*cos(thet)    ! example
      End Do

c---------
      else     ! read the primary distribution
c---------

       open (8,file="prd_2d.inp")

        read (8,*) NSG1

        Do i=1,NSG1
         read (8,*) idle,X(i),Y(i),P1(i),P2(i),P3(i),P5(i)
        End Do

        NSG  = NSG1-1
        NSG2 = NSG +2

        X(0) = X(NSG)      ! wrap
        Y(0) = Y(NSG)

        X(NSG2) = X(2)     ! wrap
        Y(NSG2) = Y(2)

      close (8)

c-----------
      end if
c-----------

c-----------------
c printing session
c-----------------

      write (6,*)
      write (6,*) ' Primary distribution and properties'
      write (6,*)

      write (2,104) NSG1

      Do i=1,NSG1
        write (2,104) i,X(i),Y(i),P1(i),P2(i),P3(i)
     +                ,P4(i),P5(i)
        write (6,104) i,X(i),Y(i),P1(i),P2(i),P3(i)
     +                ,P4(i),P5(i)
      End Do

c---------------
c redistribution
c---------------

      ICH1 = 1
      ICH2 = 1
      ICH3 = 1

      call prd_2d
     +
     +   (NSG
     +   ,ICH1,THMAX
     +   ,ICH2,SPMAX
     +   ,ICH3,SPMIN
     +   ,arcl,area
     +   ,xcentr,ycentr
     +   ,P1,P2,P3,P4,P5
     +   ,Istop
     +   )

      NSG1 = NSG+1

c---------
c printing
c---------

      write (6,*)
      write (6,*) ' Revised distribution and properties'
      write (6,*)

      write (2,104) NSG1

      Do i=1,NSG1
       write (2,104) i,X(i),Y(i),P1(i),P2(i),P3(i)
     +                ,P4(i),P5(i)
       write (6,104) i,X(i),Y(i),P1(i),P2(i),P3(i)
     +                ,P4(i),P5(i)
      End Do

      write (2,104) Null

c--------------
c more printing
c--------------

      write (6,*)
      write (6,200) arcl
      write (6,201) area
      write (6,202) xcentr
      write (6,203) ycentr
      write (6,*)

c-------
c repeat
c-------

      Go to 98   ! repeat

c-----
c done
c-----

  99  Continue

      close (2)

 104  Format (1X,I3,20(1X,F8.4))

 200  Format (1X," Arc length =",f15.10)
 201  Format (1X," Area       =",f15.10)
 202  Format (1X," x-centroid =",f15.10)
 203  Format (1X," y-centroid =",f15.10)

      stop
      end
