      program prd_3d_pr_dr

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c----------------------------------------
c  Adaptive representation of a periodic 3D line
c
c  If (x,y,z) is a point on the line,
c  then (x+RLx, y+RLy, x+RLz) is also a point
c  on the line
c
c  (xc,yc,zc) arc centers
c
c  Pj: Properties defined along the line
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension X(0:501),Y(0:501),Z(0:501)
      Dimension P1(0:501),P2(0:501),P3(0:501),P4(0:501),P5(0:501)

      Dimension xc(500),yc(500),zc(500),r(500),s(500)
      Dimension chi1(500),chi3(500)

      common/xxyyzz/x,y,z
      common/ARCC/xc,yc,zc,r,s,chi1,chi3

c---
c constants
c---

      pi  = 3.14159 265358 979323 84 D0
      pi2 = 2.0D0*pi

      Null = 0

c-----------------
c input parameters
c-----------------

      open (1,file="prd_3d_pr.dat")

       read (1,*) RLX,RLY,RLZ
       read (1,*) phimax
       read (1,*) arcmax
       read (1,*) arcmin

      close (1)

      write (6,*) 
      write (6,*) ' Enter: '
      write (6,*) 
      write (6,*) ' 1 to generate a shape '
      write (6,*) ' 2 to read from file: prd_3d_pr.inp'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) Iread

      if(Iread.eq.0) Go to 99

c---
c prepare to run
c---

      phimax = phimax *pi

c---
c recording file
c---

      open (2,file="prd_3d_pr.net")

c---------------
c generate the primary distribution
c---------------

      if(Iread.eq.1) then

      write (6,*)
      write (6,*)  " Please enter the number of segments"
      write (6,*)
      write (6,*)  " Must be less than 500"
      write (6,*)  "----------------------"
      read  (5,*) NSG

      NSG1 = NSG+1
      step = pi2/(NSG1-1.0)

     
      if(RLX.gt.0) then
         wnx = pi2/RLX      ! wave number
      else
         wnx = 10000000.0   ! fake infinity
      end if

      if(RLY.gt.0) then
         wny = pi2/RLY      ! wave number
      else
         wny = 10000000.0   ! fake infinity
      end if

      if(RLZ.gt.0) then
         wnz = pi2/RLZ      ! wave number
      else
         wnz = 10000000.0   ! fake infinity
      end if

      Do i=1,NSG1
        arg   = (i-1.0)*step
        X (i) = arg/wnx + 0.10*cos(arg)
c       Y (i) = arg/wny + sin(arg)
        Y (i) = arg/wny + 0.10*sin(arg)
c       Z (i) = arg/wny
        Z (i) = arg/wnz + 0.02*cos(2.0*arg)
        P1(i) = cos(2.0*arg)**2
        P2(i) = cos(2.0*arg)**3
        P3(i) = cos(2.0*arg)**4
        P4(i) = cos(2.0*arg)**5
        P5(i) = cos(2.0*arg)**6
      End Do

c---------------
c read the primary distribution
c---------------

      else

      open (8,file="prd_3d_pr.inp")

        read (8,*) NSG1
        NSG  = NSG1-1
        Do i = 1,NSG1
          read (8,*) X(i),Y(i),Z(i),P1(i),P2(i),P3(i)
     +              ,P4(i),P5(i)
        End Do

      close (8)

c---
      end if
c---

c---
c printing session
c---

      write (6,*)
      write (6,*) "Original distribution"
      write (6,*)
      write (2,104) NSG1
      write (6,104) NSG1

      Do i=1,NSG1
        write  (2,104) i,X(i),Y(i),Z(i)
c                     ,P1(i),P2(i)
c    +                ,P3(i),P4(i),P5(i)
        write (6,104) i,X(i),Y(i),Z(i)
c                     ,P1(i),P2(i)
c    +                ,P3(i),P4(i),P5(i)
      End Do

c---------------
c redistribution
c---------------

      Ich1 = 1
      Ich2 = 1
      Ich3 = 1

      call prd_3d_pr
     +
     +   (NSG
     +   ,RLX,RLY,RLZ
     +   ,Ich1,phimax
     +   ,Ich2,arcmax
     +   ,Ich3,arcmin
     +   ,P1,P2,P3,P4,P5
     +   ,Istop
     +   )

c---
c printing
c---

      NSG1 = NSG+1

      write (6,*)
      write (6,*) "Revised distribution"
      write (6,*)
      write (2,104) NSG1

      Do i=1,NSG1
        write (2,104) i,X(i),Y(i),Z(i)
c                     ,s(i),P1(i),P2(i)
c    +                ,P3(i),P4(i),P5(i)
        write (6,104) i,X(i),Y(i),Z(i)
c                     ,s(i),P1(i),P2(i)
c    +                ,P3(i),P4(i),P5(i)
      End Do

c-------
c finish
c-------

   99 Continue

      write (2,104) Null
      Close (2)
       
 104  Format (1X,I3,10(1X,F8.4))

      stop
      end
