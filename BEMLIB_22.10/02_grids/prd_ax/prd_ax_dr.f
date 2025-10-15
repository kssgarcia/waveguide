      program prd_ax_dr

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------
c  Point distribution on a line
c  representing the trace of an
c  axisymmetric surface in the xy plane
c
c LEGEND:
c ------
c
c P1 - P5:  properties of marker points
c
c-------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension  X (0:900), Y(0:900)
      Dimension  P1(0:900),P2(0:900),P3(0:900)
      Dimension  P4(0:900),P5(0:900),P6(0:900),P7(0:900)
      Dimension  XC(900), YC(900),  R(900),   S(900)
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

c-----------------
c input parameters
c-----------------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 0 to quit'
      write (6,*) ' 1 to generate shape '
      write (6,*) ' 2 to read it from file: prd_ax.inp'
      write (6,*) ' ----------------------------------'
      read (5,*) Iread

      open (1,file="prd_ax.dat",status="unknown")

       read (1,*) THmax
       read (1,*) SPmin
       read (1,*) SPmax

      close (1)

c---------------
c prepare to run
c---------------

      THmax = THmax *pi

      open (2,file="prd_ax.out")

c------------------------------
c generate the primary distribution
c------------------------------

      If(Iread.eq.1) then

       write (6,*)
       write (6,*) ' Will generate an elliptical shape'
       write (6,*)
       write (6,*) ' Please enter the x position of the center'
       write (6,*) ' ----------------------------------------'
       read  (5,*) xcntr

       write (6,*)
       write (6,*) ' Please enter the x semi-axis'
       write (6,*) ' ----------------------------'
       read  (5,*) a

       write (6,*)
       write (6,*) ' Please enter the y semi-axis'
       write (6,*) ' ----------------------------'
       read  (5,*) b

       write (6,*)
       write (6,*) ' Please enter the number of segments'
       write (6,*) ' -----------------------------------'
       read  (5,*) NSG

c--------
c distribution:
c--------

       NSG1 = NSG+1
       Dth  = pi/(NSG1-1.0D0)

       Do i=0,NSG1
        thet  = (i-1.0)*Dth
         X(i) = xcntr + a * Dcos(thet)
         Y(i) =         b * Dsin(thet)
         S(i) = a*thet
        P1(i) = Dcos(thet)**3
        P2(i) = Dsin(thet)**3
        P3(i) = 0.50*sin(thet)**2
        P4(i) = S(i)
        P5(i) = 0.10*cos(thet)+0.30*cos(thet)
        P6(i) = 0.0D0
        P7(i) = 0.0D0
       End Do

c---------
      Else
c---------

      open (8,file="prd_ax.inp")

       read (8,*) NSG1

       NSG  = NSG1-1
       NSG2 = NSG +2

       Do i=1,NSG1
        read (8,*) idle,X(i),Y(i),P1(i),P2(i)
        P3(i) = 0.0D0
        P4(i) = 0.0D0
        P5(i) = 0.0D0
       End Do

      Close (8)

c-----------
      End If
c-----------

c-----------------
c Printing session
c-----------------

      write (6,*)
      write (6,*) 'Primary Distribution' 
      write (6,*)

      write (2,104) NSG1

      Do i=1,NSG1
        write (2,104) i,X(i),Y(i),P1(i),P2(i),P3(i)
     +                ,P4(i),P5(i),P6(i),P7(i)
        write (6,104) i,X(i),Y(i),P1(i),P2(i),P3(i)
     +                ,P4(i),P5(i),P6(i),P7(i)
      End Do

c---------------
c redistribution
c---------------

      ICH1 = 1    ! check enabled
      ICH2 = 1    ! check enabled
      ICH3 = 1    ! check enabled

      call prd_ax
     +
     +  (NSG
     +  ,ICH1,THMAX
     +  ,ICH2,SPMAX
     +  ,ICH3,SPMIN
     +  ,Arct,Aret,Volt
     +  ,Xcnt,Ycnt
     +  ,P1,P2,P3,P4,P5,P6,P7
     +  ,Istop
     +  )

      NSG1 = NSG+1

c---------
c printing
c---------

      write (6,*)
      write (6,*) ' Revised Distribution'
      write (6,*)
      write (2,104) NSG1

      Do i=1,NSG1
       write (2,104) i,X(i),Y(i),P1(i),P2(i),P3(i),P4(i),P5(i)
       write (6,104) i,X(i),Y(i),P1(i),P2(i),P3(i),P4(i),P5(i)
      End Do
      write (2,104) Null

c---------------------------
c normalization and printing
c---------------------------

      Arct = Arct/(pi*a)               ! arc length
      Aret = Aret/(4.0D0*pi*a**2)        ! surface area
      Volt = Volt/(4.0D0*pi*a**3/3.0D0)    ! volume

      write (6,*)
      write (6,200) Arct
      write (6,201) Aret
      write (6,202) volt
      write (6,203) Xcnt,Ycnt
       
c-------------------
c Smooth property P1
c-------------------

      call sm_ax (NSG,P1)

      write (6,*)
      write (6,*)  " Smoothed property P1"
      write (6,*)
      write (2,104) NSG1

      Do i=1,NSG1
       write (2,104) i,X(i),Y(i),P1(i),P2(i),P3(i),P4(i),P5(i)
       write (6,104) i,X(i),Y(i),P1(i),P2(i),P3(i),P4(i),P5(i)
      End Do

c-----
c Done
c-----

      close (2)

 104  Format (1X,I3,20(1X,F8.4))
 200  Format (1x,"Arc length   =",f15.10)
 201  Format (1x,"Surface area =",f15.10)
 202  Format (1x,"Volume       =",f15.10)
 203  Format (1x,"Centroid     =",f15.10,f15.10)

      Stop
      End
