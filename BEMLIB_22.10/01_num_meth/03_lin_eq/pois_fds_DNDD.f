      subroutine pois_fds_DNDD 
     +
     +   (ax,bx
     +   ,ay,by
     +   ,Nx,Ny
     +   ,g
     +   ,w,q,z,v
     +   ,mat,rhs
     +   ,mats
     +   )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c              C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c------------------------------------------------

c------------------------------------------------
c Generate the finite-difference linear system 
c for the Poisson equation
c in a rectangle confined in ax<x<bx, ay<y<by
c
c Equation is:    Lalp(f) + g = 0
c
c System is :     mat * x = rhs
c
c Boundary conditions:   f = w   at x = ax
c                    df/dx = q   at x = bx
c                        f = z   at y = ay
c                        f = v   at y = by
c
c System is produced by the method of impulses
c described by Pozrikidis (1998)
c
c Phantom nodes are used for the implementation
c of the Neumann boundary condition
c
c SYMBOLS:
c -------
c
c Nx...	intervals in x direction
c Ny...	intervals in y direction
c mats .. System size
c mat..	  Finite-difference matrix
c rhs..   Right-hand side
c
c unknown vector is comprised of sequential values (i,j):
c (horizonal and then up)
c
c  2,2   3,2   4,2    ...  N-1,2   N,2   N+1,2
c  2,3   3,3   4,3    ...  N-1,3   N,2   N+1,3
c  ....
c
c  2,M   3,M   4,M    ...  N-1,M   N,M   N+1,M
c
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision mat(512,512),rhs(512)

      Dimension f(0:128,0:128),g(0:128,0:128)
      Dimension w(512),q(512),z(512),v(512)

      Integer p,s,t

c-------------
c preparations
c-------------

      Nx1 = Nx+1
      Ny1 = Ny+1
      Dx = (bx-ax)/(Nx1-1.0)
      Dy = (by-ay)/(Ny1-1.0)

c------------------
c more preparations
c------------------

      Dx2  = 2.0D0*Dx
      Dy2  = 2.0D0*Dy
      Dxs  = Dx**2
      Dys  = Dy**2
      beta = Dxs/Dys

c-------------------
c initialize to zero
c-------------------

      Do j=2,Ny
        Do i=2,Nx+1
         f(i,j)=0.0D0
        End Do
      End Do

c------------------------------
c Dirichlet boundary conditions
c------------------------------

      Do j=2,Ny
       f(1,j) = w(j)     ! on left side
      End Do

      Do i=2,Nx+1
       f(i,   1) = z(i)     ! down
       f(i,Ny+1) = v(i)     ! up
      End Do

c---
c Neumann boundary condition for zero impulse
c to compute the rhs
c---

      Do j=2,Ny
       f(Nx+2,j) = Dx2*q(j)
      End Do

      p = 0                              ! counter

      Do j=2,Ny
        Do i=2,Nx+1
         p = p+1
         R =                 f(i+1,j)
     +       -2.0D0*(1.0D0+beta)*f(i  ,j)
     +       +               f(i-1,j)
     +       +         beta *f(i,j-1)
     +       +         beta *f(i,j+1)
     +       + Dxs*g(i,j)
         rhs(p) = - R
        End Do
      End Do

      mats = p                           ! system size

c------------------------------------
c will scan row-by-row to compute mat
c------------------------------------

      t = 0    ! counter

      Do s=2,Ny
        Do l=2,Nx+1

          f(l,s) = 1.0D0  ! impulse
          t = t+1
c---
c phantom nodes
c---

          Do k=2,Ny
           f(Nx+2,k) = f(Nx,k) + Dx2*q(k)
          End Do

c---
c regular nodes
c---

          p = 0      ! counter

          Do j=2,Ny
            Do i=2,Nx+1
             p = p+1
             R = f(i+1,j)
     +           -2.0D0*(1.0D0+beta)*f(i  ,j)
     +           +                   f(i-1,j)
     +           +              beta*f(i,j+1)
     +           +              beta*f(i,j-1)
     +           + Dxs*g(i,j)
            mat(p,t) = R+rhs(p)
            End Do
          End Do

        f(l,s) = 0.0D0   ! reset

        End Do
      End Do

c-----
c Done
c-----

      Return
      End
