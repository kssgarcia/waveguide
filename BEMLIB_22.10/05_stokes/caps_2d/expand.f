      subroutine expand
     +
     +   (Imax
     +   ,Amat,rhs
     +   ,roexp
     +   ,method
     +   )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

      Implicit Double Precision (a-h,o-z)

      Dimension XC (900),YC (900),  R(900),   S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

      Dimension Amat(500,500),rhs(500)

      Dimension evt(500)

c-------
c common
c-------

      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT

c--------
c prepare
c--------

      Imax2 = 2*Imax

c-----------------------------------------
c first method: singular preconditioning
c                of the linear system
c-----------------------------------------

      If(method.eq.1) then

       write (6,*) " expand: singular preconditioning"

c---
c eigenvector of the adjoint
c---

      evtn = 0.0D0     ! norm

      Do i=1,Imax
       iImax = i+Imax
       DSS        = 0.5D0*R(i)*(th3(i)-th1(i))
       evt(i)     = DSS * Dcos(th2(i))
       evt(iImax) = DSS * Dsin(th2(i))
       evtn = evtn + evt(i)**2+evt(iImax)**2
      End Do

      evtn = Dsqrt(evtn)

      Do i=1,Imax2
       evt(i) = evt(i) /evtn
      End Do

c---
c precondition the matrix
c---

      Do j=1,Imax2
        sum = 0.0D0
        Do k=1,Imax2
         sum = sum+evt(k)*Amat(k,j)
        End Do
c       write (6,101) sum
        Do m=1,Imax2
         Amat(m,j) = Amat(m,j)-evt(m)*sum
        End Do
      End Do

c---
c precondition the right-hand side
c---

      sum = 0.0D0
      Do j=1,Imax2
        sum = sum+evt(j)*rhs(j)
      End Do
c     write (6,101) sum

      Do i=1,Imax2
       rhs(i) = rhs(i) - evt(i)*sum
      End Do

c---
c last equation
c is constraint on flow rate
c---

      Do j=1,Imax2
       Amat(Imax2,j) = evt(j)*evtn
      End Do

      rhs(Imax2) = roexp

c--------------------------------------------
c second method: deflate the integral equation
c--------------------------------------------

      Else If(method.eq.2) then

       write (6,*) " expand: deflation"

       Do j=1,Imax

        jImax = j+Imax
        DSS   = 0.50D0*R(j)*(th3(j)-th1(j))
        ttx = DSS * Dcos(th2(j))
        tty = DSS * Dsin(th2(j))

        Do i=1,Imax
         iImax = i+Imax
         fc = ornt(i)
         vnx = fc*Dcos(th2(i))
         vny = fc*Dsin(th2(i))
         Amat(i    ,j)     = Amat(i    ,j)    - vnx*ttx
         Amat(iImax,j)     = Amat(iImax,j)    - vny*ttx
         Amat(i    ,jImax) = Amat(i    ,jImax)- vnx*tty
         Amat(iImax,jImax) = Amat(iImax,jImax)- vny*tty
        End Do

       End Do

       Do i=1,Imax
         iImax = i+Imax
         fc = ornt(i)
         vnx = fc*Dcos(th2(i))
         vny = fc*Dsin(th2(i))
         rhs(i)     = rhs(i)     - vnx*roexp
         rhs(iImax) = rhs(iImax) - vny*roexp
       End Do

c-----------
      End If
c-----------

c-----
c Done
c-----

      return
      end
