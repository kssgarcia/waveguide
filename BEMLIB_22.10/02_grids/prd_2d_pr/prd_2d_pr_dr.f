      program prd_2d_pr_dr

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------
c  Point distribution on a line
c  that is periodic along the x-axis
c--------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  X(0:513), Y(0:513)
      Dimension P1(0:513),P2(0:513)
      Dimension P3(0:513)
      Dimension P4(0:513),P5(0:513)

      Dimension  XC(513), YC(513),  R(513),   S(513)
      Dimension TH1(513),TH2(513),TH3(513),ORNT(513)

c--------------
c common blocks
c--------------

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/PPII/PI,PIH,PI2,PI4,PI8

c----------
c constants
c----------

      pi  = 3.14159 265358 979323 84D0
      pi2 = 2.0D0*pi

      Null = 0

c-----------------
c input parameters
c-----------------

      open (1,file="prd_2d_pr.dat",status="unknown")

       read (1,*) thmax
       read (1,*) spmin
       read (1,*) spmax
       read (1,*) thrs1,thrs2
       read (1,*) 
       read (1,*) method
       read (1,*) 

c---------------
c prepare to run
c---------------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to generate periodic shape '
      write (6,*) ' 2 to read from file: prd_2d_pr.inp '
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'

      read (5,*) Iread

      if(Iread.eq.0) Go to 99

c----
c open output files
c----

      open (2,file="prd_2d_pr.out",status="unknown")

c----------------------------------
c Generate the primary distribution
c----------------------------------

      if(Iread.eq.1) then

c      write (6,*)
c      write (6,*) ' Enter the period: L '
c      write (6,*) ' --------------------'
c      read  (5,*)   RL
       read (1,*) RL

c      write (6,*) 
c      write (6,*) ' Will generate a sinusoidal wave'
c      write (6,*) 
c      write (6,*) ' Enter the amplitude'
c      write (6,*) ' -------------------'
c      read  (5,*)  amp
       read (1,*) amp

c      write (6,*) 
c      write (6,*) ' Enter the number of uniquepoints'
c      write (6,*) ' maximum is 512'
c      write (6,*) ' --------------'
c      read  (5,*)   NSG
       read (1,*) NSG

       NSG1 = NSG+1
       NSG2 = NSG+2
       NSG3 = NSG+3

       step = RL/(NSG1-1.0)
       wn = pi2/RL

       ampx = 0.5*step

       Do i=0,NSG2

        xflat   = (i-1.0D0)*step
        tmp     = xflat*pi2
        X(i) = xflat + ampx*Dsin(xflat*wn)
        Y(i) =          amp*Dcos(X(i)*wn)

       End Do

      write (6,104) 

c---------
      else
c---------

      open (8,file="prd_2d_pr.inp",status="unknown")

       read (8,*) NSG1

       NSG  = NSG1-1
       NSG2 = NSG +2
       NSG3 = NSG +3

       Do i=1,NSG1
         j = NSG1+1-i
         read (8,*) idle,Y(j),X(j)
       End Do

      close (8)

c-----
c wrap
c-----

       RL = X(NSG1)-X(1)
 
       X(0) = X(NSG)-RL
       Y(0) = Y(NSG)

       X(NSG2) = X(2)+RL
       Y(NSG2) = Y(2)

c-----------
      end if
c-----------

c---------
c printing 
c---------

      write (6,*)
      write (6,*) " prd_2d_pr_dr: primary distribution:"
      write (6,*)
      write (2,104) NSG3

      Do i=0,NSG2
        write (2,104) i,X(i),Y(i)
        write (6,104) i,X(i),Y(i)
      End Do

c---
c prepare for redistribution
c---

      thmax = thmax * pi
      thrs1 = thrs1 * RL
      thrs2 = thrs2 * RL

      ICH1 = 1
      ICH2 = 1
      ICH3 = 1
      ICH4 = 1

c---------------
c redistribution
c---------------

c-------------------------
      if(method.eq.1) then
c-------------------------

      call prd_2d_pr_arc
     +
     +   (NSG,RL
     +   ,ICH1,thmax
     +   ,ICH2,SPMAX
     +   ,ICH3,SPMIN
     +   ,ICH4,thrs1,thrs2
     +   ,arc_length,area
     +   ,unused1,unused2
     +   ,P1,P2,P3
     +   ,P4,P5
     +   ,Istop
     +   )

c------------------------------
      else if(method.eq.2) then
c------------------------------

      call prd_2d_pr_splc
     +
     +   (NSG
     +   ,X,Y
     +   ,RL
     +   ,ICH1,THMAX
     +   ,ICH2,SPMAX
     +   ,ICH3,SPMIN
     +   ,ICH4,thrs1,thrs2
     +   ,P1,P2,P3
     +   ,P4,P5
     +   ,Istop
     +   )

c-----------
      end if
c-----------

      NSG1 = NSG+1
      NSG2 = NSG+2
      NSG3 = NSG+3

c-------------------------------
c print the revised distribution
c-------------------------------

      write (6,*)
      write (6,*) " prd_2d_pr: revised distribution"
      write (6,*)

      write (2,104) NSG3

      Do i=0,NSG2
        write (2,104) i,X(i),Y(i)
        write (6,104) i,X(i),Y(i)
      End Do

      write (2,104) Null

      close (2)

c----
c print area and arc length
c----
       
      if(method.eq.1) then
        write (6,200) arc_length
        write (6,201) area
      end if

c-----
c done
c-----
 
  99  Continue

 104  Format (1X,I3,1X,F12.8,1X,F12.8,1X,I1)

 200  Format (" Arc length                            = ",f15.10)
 201  Format (" Area between the curve and the x axis = ",f15.10)

      Stop
      End
