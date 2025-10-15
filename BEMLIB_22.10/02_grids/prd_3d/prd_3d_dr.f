      program prd_3d_dr

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c----------------------------------------
c  Point distribution on a closed 3D line
c
c  The line is approximated locally with
c  a circular arc passing through 3
c  successive points
c
c  LEGEND:
c  ------
c
c  (xc,yc,zc) arc center
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension X(0:501),Y(0:501),Z(0:501)

      Dimension P1(0:501),P2(0:501),P3(0:501)
      Dimension P4(0:501),P5(0:501)

      Dimension xc(500),yc(500),zc(500),r(500),s(500)
      Dimension chi1(500),chi3(500)

      common/xxyyzz/X,Y,Z
      common/ARCC/xc,yc,zc,r,s,chi1,chi3

c----------
c constants
c----------

      pi  = 3.14159 265358 979323 84 D0
      pi2 = 2.0D0*pi

      Null = 0

c-----------------
c input parameters
c-----------------

      write (6,*) 
      write (6,*) ' Enter: '
      write (6,*)
      write (6,*) ' 1 to generate a shape '
      write (6,*) ' 2 to read from file: prd_3d.inp'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read (5,*) Iread

      if(Iread.eq.0) Go to 99

      open (1,file="prd_3d.dat")

       read (1,*) phimax
       read (1,*) arcmin
       read (1,*) arcmax

      close (1)

c--------
c prepare
c--------

      phimax= phimax *pi

c--------
c recording file
c--------

      open (2,file="prd_3d.net")

c-----------------------

      if(Iread.eq.1) then

c----------------------------------
c Generate the primary distribution
c----------------------------------

      write (6,*) 
      write (6,*) ' Enter the number of points'
      write (6,*) ' Try 16 (maximum is 500)'
      write (6,*) ' --------------'
      read  (5,*)  NSG

      NSG1 = NSG+1
      step = pi2/(NSG1-1.0D0)

      Do i=0,NSG1
        arg   = (i-1.0D0)*step
        X (i) = Dcos(arg)
c       Y (i) = Dsin(arg)
c       Z (i) = 0.
        Y (i) = 2.0*sin(arg)
        Z (i) = cos(2.0*arg)
        P1(i) = cos(2.0*arg)**2
        P2(i) = cos(2.0*arg)**3
        P3(i) = cos(2.0*arg)**4
        P4(i) = cos(2.0*arg)**5
        P5(i) = cos(2.0*arg)**6
      End Do

c---------
      else  !  read the primary distribution
c---------

      open (8,file="prd_3d.inp")

       read (8,*) NSG1

       NSG = NSG1-1

       Do i=1,NSG1
         read (8,*) X(i),Y(i),Z(i)
       End Do

       close (8)

c-----------
      end if
c-----------

c---
c printing session
c---

      write (6,*)
      write (6,*) "Original distribution"
      write (6,*)

      write (2,104) NSG1

      Do i=1,NSG1
        write (2,104) i,X(i),Y(i),Z(i),P1(i),P2(i)
c    +                ,P3(i),P4(i),P5(i)
        write (6,104) i,X(i),Y(i),Z(i),P1(i)
c    +                ,P3(i),P4(i),P5(i)
      End Do

c---
c redistribution
c---

      write (6,*) " prd_3d_dr: entering prd_3d_cl"

      Ich1 = 1
      Ich2 = 1
      Ich3 = 1

      call prd_3d
     +
     +   (NSG
     +   ,Ich1,phimax
     +   ,Ich2,arcmax
     +   ,Ich3,arcmin
     +   ,P1,P2,P3,P4,P5
     +   ,Istop
     +   )

      NSG1 = NSG+1

c---
c printing
c---

      write (6,*)
      write (6,*) "Revised distribution"
      write (6,*)
      write (2,104) NSG1

      Do i=1,NSG1
        write (2,104) i,X(i),Y(i),Z(i),s(i),P1(i),P2(i)
c    +                ,P3(i),P4(i),P5(i)
        write (6,104) i,X(i),Y(i),Z(i),s(i),P1(i),P2(i)
c    +                ,P3(i),P4(i),P5(i)
      End Do

c-----
c Done
c-----

   99 Continue

      write (2,104) Null
      close (2)
       
 104  Format (1X,I3,10(1X,F8.4))

      Stop
      End
