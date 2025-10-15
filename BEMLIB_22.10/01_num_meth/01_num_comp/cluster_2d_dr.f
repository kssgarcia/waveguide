      program cluster_2d_dr

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c========================================

c----------------------------------------
c Arrange particles into clusters
c
c SYMBOLS:
c --------
c
c Nprtcl:  number of particles
c Ipl(i):  original label of the ith particle
c Ncls:    Number of clusters
c
c eps:   parameter for cluster identification
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension     x(25),  y(25)
      Dimension  Lump(25),Ipl(25)

      Dimension x_orig(25),y_orig(25)

      common/points/x,y

c----------
c constants
c----------

      Null = 0

c----
c read the particle positions
c----

      open (1,file="cluster_2d.dat")

       read (1,*) eps
       read (1,*) Nprtcl

       Do i=1,Nprtcl
        read (1,*) idle,x(i),y(i)
        Ipl(i) = i
       End Do

      close (1)

c---
c printing
c---

      write (6,*)
      write (6,*) " Initial particle labels and position"
      write (6,*)

      Do i=1,Nprtcl
        write (6,101) Ipl(i),x(i),y(i)
      End Do

c----------------------
c arrange into clusters
c----------------------

      call cluster_2d
     +
     +   (Nprtcl
     +   ,eps
     +   ,Ncls,Lump
     +   ,Ipl
     +   )

c-----------------------------------
c save the original position vectors
c-----------------------------------

      Do i=1,Nprtcl
       x_orig(i) = x(i)
       y_orig(i) = y(i)
      End Do

c---
c  redefine the position vector
c---

      Do i=1,Nprtcl
       j    = Ipl(i)
       x(i) = x_orig(j)
       y(i) = y_orig(j)
      End Do

c---
c  A property of the ith new particle
c  is also a property of the Ipl(i)
c  original particle;
c
c  for example,
c  if u(i) is the x velocity of the ith new particle,
c  then dx(x_orig(Ipl(i))) dt = u(i)
c---

c---
c printing
c---

      open (2,file="cluster_2d.out")

      write (6,*)
      write (6,*) "Number of clusters: ", Ncls
      write (6,*)
      write (6,*) "Cluster sizes:"
      write (6,*)

      Icount = 0
      Do i=1,Ncls
        write (6,100) i,Lump(i)
        Icount = Icount+Lump(i)
      End Do

      write (6,*)
      write (6,*) "total:",Icount

c----------------------------
c print in groups of clusters
c----------------------------

      write (6,*)
      write (6,*) " Particle labels and position"
      write (6,*)

      Il = 1  ! counts lumps
      Ic = 1  ! counter particles in a lump

      write (6,100) Lump(Il)

      write (2,100) Lump(Il)

      Do i=1,Nprtcl

        write (2,100) i,Ipl(i),x(i),y(i)
        write (6,100) i,Ipl(i),x(i),y(i)

        if(Ic.eq.Lump(Il)) then
         Ic = 1
         Il = Il+1
         if(Il.le.Ncls) then 
            write (6,100) Lump(Il)
            write (2,100) Lump(Il)
         end if
        else
         Ic = Ic+1
        end if

      End Do

      write (2,*) Null
      close (2)

c-----
c done
c-----

 100  Format (2(1x,i3),3(1x,f10.5))
 101  Format (1x,i3,3(1x,f10.5))

      stop
      end
