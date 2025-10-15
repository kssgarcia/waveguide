      program pois_fds_DNDD_dr

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c-----------------------------------------
c Driver for generating a linear system
c for the finite-difference solution
c of the Poisson equation lapl(f)+g = 0,
c as discussed in the text
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision MAT(512,512),rhs(512)

      Dimension g(0:128,0:128)
      Dimension w(512),q(512),z(512),v(512)

c---
c input
c---

      write (6,*)
      write (6,*) " Please enter the box dimensions ax, bx, ay, by"
      write (6,*) " ----------------------------------------------"
      read  (5,*) ax,bx,ay,by

      write (6,*)
      write (6,*) " Please enter the discretization levels Nx and Ny"
      write (6,*) " ----------------------------------------------"
      read  (5,*) Nx,Ny

c--------------------
c boundary conditions
c and source term
c--------------------

      Do j=1,Ny+1
       w(j) = 0.0D0   ! example
       q(j) = 1.0D0   ! example
      End Do

      Do i=1,Nx+1
       z(i) = 0.0D0   ! example
       v(i) = 0.0D0   ! example
      End Do

      Do i=1,Nx+1
       Do j=1,Ny+1
        g(i,j) = 0.0D0  ! example
       End Do
      End Do

c---------------------------
c generate the linear system
c---------------------------

      call pois_fds_DNDD 
     +
     +   (ax,bx
     +   ,ay,by
     +   ,Nx
     +   ,Ny
     +   ,g
     +   ,w,q,z,v
     +   ,MAT,Rhs
     +   ,mats
     +   )

c-----------------
c printing session
c-----------------

      write (6,*)
      write (6,*)  "Linear system:"
      write (6,*)

      Do i=1,mats
       write (6,100) (MAT(i,j),j=1,mats),rhs(i)
      End Do

c-----
c done
c-----

  100 Format (100(1x,f6.3))

      Stop
      End
