      subroutine pois_gs_dpr
     +
     +  (N,M
     +  ,Dx,Dy
     +  ,rhs
     +  ,Itermax
     +  ,tol
     +  ,relax
     +  ,Ispeak1
     +  ,Ispeak2
     +  ,T
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------
c Gauss-Siedel solution of Poisson's equation
c in a rectangular domain
c with the homogeneous Dirichlet boundary condition
c at the top and bottom boundaries
c and the periodic boundary condition
c at the left and right boundaries
c
c
c The solution is found by point GS iterations
c ...........................................
c
c SYMBOLS:
c -------
c
c N=Nx:   Grid size in the x direction
c M=Ny:   Grid size in the y direction
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension     T(129,129)
      Dimension Tsave(129,129)

      Dimension rhs(16666)

c--------
c prepare
c--------

      Na = N-1
      N1 = N+1
      N2 = N+2

      Ma = M-1
      M1 = M+1
      M2 = M+2

      Dxs = Dx**2
      Dys = Dy**2

      beta  = Dxs/Dys
      beta1 = 2.0D0*(beta+1.0D0)

      open (9,file="error.out")

c--------
c verbose
c--------

      If(Ispeak1.eq.1) then
       write (6,*) " pois_gs_pr3: Gauss-Siedel iterations:"
      End If

c---
c Implement the homogeneous Dirichlet
c boundary condition
c---

       Do i=1,N2
         T(i, 1) = 0.0D0  ! bottom
         T(i,M1) = 0.0D0  ! top
       End Do

c------------------------
c Gauss-Siedel iterations
c------------------------

      Do 97 Iter=1,Itermax

c------------------------
c update nodes row-by-row
c------------------------

c---
c save
c---
       Do J=1,M1
        Do I=2,N1
          Tsave(I,J) = T(I,J)
        End Do
       End Do

c---
c iterate
c---

       Ic = 0

       Do J=2,M

        Do I=2,N1
         Ic = Ic+1

         cor = (        T(I+1,J) + T(I-1,J)
     +          +beta*( T(I,J+1) + T(I,J-1))
     +          -rhs(Ic) 
     +          )/beta1 - T(I,J)

         T(I,J) = Tsave(I,J) + relax*cor

c        write (6,113) Ic,rhs(Ic)
        End Do
c       write (6,108)

       End Do
c      pause

c-------------------
c maximum correction
c-------------------

       cormax = 0.0D0

       Do J=1,M1
        Do I=2,N1
         testo = T(I,J)-Tsave(I,J)
         If(Dabs(testo).gt.cormax) cormax = Dabs(testo)
        End Do
       End Do

c      write (6,*) Iter,T(3,3)

c------------------------------------------------
c wrap: implement the periodic boundary condition
c------------------------------------------------

      Do j=1,M1
        T(1,j)  = T(N1,j)
        T(N2,j) = T(2,j)
      End Do

c---
c      write (6,*)
c      write (6,117) iter
c      Do i=1,N1
c        write (6,110) (T(i,j),j=1,M1)
c      End Do
c---

       If(Ispeak1.eq.1) then
         write (6,201) Iter,cormax,log10(abs(cormax))
       End If

       If(cormax.lt.tol) Go to 98

  97  Continue

c-------
c failed
c-------

      write (6,*) " pois_gs_pr3: Gauss-Siedel iterations"
      write (6,*) "              did not converge"

c     stop

c------------------
c solution computed
c------------------

  98  Continue

c---
c wrap
c---

      Do j=1,M1
       T(N1,j) = T(1,j)
       T(N2,j) = T(2,j)
      End Do

      If(Ispeak2.eq.1) then
        write (6,*) " pois_gs_pr3: Gauss-Siedel iterations: ",Iter
      End If

      close (9)

c-----
c Done
c-----

 100  Format (50(1x,f8.4))
 101  Format (2(1x,i3),2(1x,f5.3))
 102  Format (" Correction factor           = ",f10.5)
 103  Format (" Residuals = ",2(1x,f15.4))
 108  Format (2(1x,i3),6(1x,f15.10))

 110  Format (30(1x,f6.3))
 113  Format (1x,i4,f15.10,1x,f15.10)
 117  Format ("Iteration: ",i5,1x,f15.8)

 201  Format (1x,i5,3(1x,f15.10))

      Return
      End
