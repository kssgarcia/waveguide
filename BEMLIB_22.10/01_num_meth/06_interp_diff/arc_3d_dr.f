      program arc_3d_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press, 1998
c------------------------------------------------

c------------------------------------------------
c Driver for arc_3d to compute the arc
c passing through three points x1,x2,x3
c in 3D space
c
c Generates nprof+1 points along the arc for plotting purposes
c
c SYMBOLS
c -------
c
c  xc,yc,zc     coordinates of the arc center
c  a            radius of the arc
c  chi1, chi3   angles subtended by the second-first 
c               and second-third point
c               as described in text
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (nprof=16)
   
c----------
c constants
c----------

      null = 0

c-----------------------------------
c Read the coordinates of the points
c-----------------------------------

      open (1,file="arc_3d.dat")

        read (1,*) x1,y1,z1   ! coordinates of the first point
        read (1,*) x2,y2,z2   ! coordinates of the second point
        read (1,*) x3,y3,z3   ! coordinates of the third point

      close (1)

c---
c Display the points
c---

      write (6,*)
      write (6,*) " Points defining the arc"
      write (6,*)

      write (6,100) x1,y1,z1
      write (6,100) x2,y2,z2
      write (6,100) x3,y3,z3

c---
c compute the arc center radius
c and subtended angles
c---

      call arc_3d 
     +
     +   (x1,x2,x3
     +   ,y1,y2,y3
     +   ,z1,z2,z3
     +   ,xc,yc,zc
     +   ,a
     +   ,chi1,chi3
     +   )

      write (6,*)
      write (6,110) xc,yc,zc,a
      write (6,*)

      write (6,100)
      write (6,111) chi1,chi3
      write (6,100)

c---
c Display the arc
c
c nprof:  number of plotting markers
c---

      open (2,file="arc_3d.out")

      Dth    = (chi3-chi1)/nprof
      as     = a**2
      nprof1 = nprof+1

      a11 = x1-xc
      a12 = y1-yc
      a13 = z1-zc

      a21 = x2-xc
      a22 = y2-yc
      a23 = z2-zc

      a31 = (y3-y1)*(z2-z1)-(z3-z1)*(y2-y1)
      a32 = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1)
      a33 = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1)

      write (6,*)
      write (6,*) " Points along the arc"
      write (6,*)

      write (2,101) nprof1

      Do i=1,nprof1

       chi = chi1+(i-1.0)*Dth
       b1 = as*cos(chi-chi1)
       b2 = as*cos(chi)
       b3 = 0.0

c---
c solve for x-xc
c---

      call cramer_33
     +
     +   (A11,A12,A13
     +   ,A21,A22,A23
     +   ,A31,A32,A33
     +   ,B1,B2,B3
     +   ,x,y,z)

       x = x +xc
       y = y +yc
       z = z +zc

       write (6,100) chi,x,y,z
c      write (2,101) i,x,y,z
       write (2,100) x,y,z

      End Do

c-----
c Done
c-----

      write (2,101) null

      close (2)

 100  Format (10(1x,f15.10))
 101  Format (1x,i3,10(1x,f15.10))
 110  Format  ("Center coordinates: ",3(1x,f10.5),/,
     +         "Radius:             ",  1x,f10.5)
 111  Format  ("Subtended angles  : ",2(1x,f10.5))

      Stop 
      End
