      subroutine eigenvectors (ev, evt)

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------------
c Generate the right and right (traspose)
c global eigenvectors
c
c SYMBOLS:
c -------
c
c Ncl:	Number of collocation points
c
c ev:  Eigenvectors of the influence matrix
c
c evt: Eigenvectors of the transpose
c      of the influence matrix
c
c NOTES:
c -----
c
c For straight segments, the eigenvectors are exact,
c whereas for the native parametrization of an ellipse,
c the eigenvectors are approximate.
c
c
c Capacity:
c --------
c
c   72 particles
c   64 elements
c
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(72),Itp(72)
 
      Dimension vnx0(4608),vny0(4608),elml(4608)
      Dimension   ev(4608), evt(4608)

      Dimension hold(72),holdt(72)    ! for normalizing eigenvectors

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl
      common/colloc2/vnx0,vny0,elml

c--------
c prepare
c--------

      Ncl2 = 2*Ncl

c--------------------------------
c Compute eigenvectors,
c         transpose eigenvectors,
c and hold their norms
c         for normalization
c--------------------------------

      Ic = 0         ! collocation point counter
      Ib = 0         ! block counter

      Do i=1,Nprtcl    ! run over particles

         hold(i) = 0.0D0  ! for normalization
        holdt(i) = 0.0D0  ! for normalization

        Do j=1,NE(i)

          Ic = Ic+1

          Ibj  = Ib  + j
          Ibji = Ibj + NE(i)

          ev(Ibj)  = vnx0(Ic)
          ev(Ibji) = vny0(Ic)

          hold(i) = hold(i) + ev(Ibj)**2+ ev(Ibji)**2

          evt(Ibj)  = vnx0(Ic)*elml(Ic)
          evt(Ibji) = vny0(Ic)*elml(Ic)

          holdt(i) = holdt(i) + evt(Ibj)**2+evt(Ibji)**2

        End Do

        Ib = Ib+2*NE(i)

        hold (i) = Dsqrt( hold(i))
        holdt(i) = Dsqrt(holdt(i))

      End Do

c----------------------------
c  Normalize the eigenvectors
c----------------------------

      Ic = 0    ! collocation point counter
      Ib = 0    ! block counter

      Do i=1,Nprtcl
       Do j=1,NE(i)

        Ic = Ic+1

        Ibj  = Ib  + j
        Ibji = Ibj + NE(i)

        ev(Ibj)  = ev(Ibj) /hold(i)
        ev(Ibji) = ev(Ibji)/hold(i)

        evt(Ibj)  = evt(Ibj) /holdt(i)
        evt(Ibji) = evt(Ibji)/holdt(i)

       End Do

       Ib = Ib+2*NE(i)

      End Do

c-------------------------
c  sum1 and sum2 should be equal to Nprtcl
c-------------------------

      sum1 = 0.0
      sum2 = 0.0

      Do i=1,Ncl2
       sum1 = sum1 + ev(i)**2
       sum2 = sum2 + evt(i)**2
c      write (6,112) ev(i),evt(i)
      End Do

      write (6,*)
      write (6,*) "eigenvectors: both should be"
      write (6,*) "              equal to", Nprtcl," :",sum1,sum2
      write (6,*)
 
c---------------------------

c-----
c Done
c-----

 112  Format (80(1x,f10.7))

      Return
      End

c=============================================

      subroutine precondition ()

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------------
c Singular preconditioning of the 
c Master Linear System (MLS)
c
c SYMBOLS:
c -------
c
c Ncl:	Number of collocation points
c
c ev:  Eigenvectors of the influence matrix
c
c evt: Eigenvectors of the transpose
c      of the influence matrix
c
c ALP:	Preconditioned Matrix
c
c BLP:	Preconditioned rhs
c
c P indicates preconditioning
c
c
c NOTES:
c -----
c
c For straight segments, the eigenvectors are exact,
c whereas for the native parametrization of an ellipse,
c the eigenvectors are approximate.
c
c Capacity:
c -----
c
c   72 particles
c   64 elements
c
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(72),Itp(72)
 
      Dimension   ev(4608), evt(4608)

      Dimension  AL(2304,2304), BL(2304)
      Dimension ALP(2304,2304),BLP(2304)
c     Dimension  AL(4608,4608), BL(4608)
c     Dimension ALP(4608,4608),BLP(4608)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl

      common/matrix/AL,BL

c--------
c prepare
c--------

      Ncl2 = 2*Ncl

c----------------------------------------- 
c Initialize preconditioned matrix and rhs
c----------------------------------------- 

      Do i=1,Ncl2

       BLP(i) = BL(i)

       Do j=1,Ncl2
        ALP(i,j) = AL(i,j)
       End Do

      End Do

c---------------------------------
c generate the global eigenvectors
c---------------------------------

      call eigenvectors (ev, evt)

c---------------------------------
c Test the discrete eigenvector
c
c sum should be very close to zero
c---------------------------------
c
c     write (6,*)
c     write (6,*) " Matrix*eigenvector"
c     write (6,*)
c     Do i=1,Ncl2
c      sum = 0.0
c      Do j=1,Ncl2
c       sum = sum+AL(i,j)*ev(j)
c      End Do
c      write (6,112) sum
c     End Do
c--------------------------

c---------------------------------
c Test the adjoint eigenvector
c
c sum should be very close to zero
c---------------------------------
c
c     write (6,*)
c     write (6,*) " Adjoint eigenvector * Matrix"
c     write (6,*)
c     Do i=1,Ncl2
c     sum = 0.0
c     Do j=1,Ncl2
c      sum = sum+AL(j,i)*evt(j)
c     End Do
c      write (6,112) sum
c     End Do
c--------------------------

c--------------------------
c Precondition the linear system
c by projecting it onto the space
c that is normal to the particle blocks
c of the global transpose eigenvector
c
c P indicates preconditioning
c-----------------------------

      Ib = 0

      Do 19 i=1,Nprtcl  ! run over particles

       Min = Ib+1          ! lower limit of ith partition
       Max = Ib+2*NE(i)    ! upper limit of ith partition

c---------------------------------
c precondition the right-hand side
c
c sum is the projection of the adjoint eigenvector
c on the right-hand side
c---------------------------------

       sum = 0.0D0
       Do k=Min,Max
        sum = sum + evt(k)*BLP(k)
       End Do

c      write (6,*) "Projection of adj eigvctr on the rhs:",sum

       Do k=Min,Max
         BLP(k) = BLP(k)-sum*evt(k)
       End Do

c------------------------
c precondition the matrix
c------------------------

       Do j=1,Ncl2

        sum = 0.0D0

        Do k=Min,Max
         sum = sum+evt(k)*ALP(k,j)
        End Do

c       write (6,*) "Projection of adj eigvctr on jth column:",sum

        Do m=Min,Max
         ALP(m,j) = ALP(m,j)-evt(m)*sum
        End Do

       End Do

      Ib = Ib +2*NE(i)

  19  Continue            ! target for loop over particles

c-------------------------
c  confirm preconditioning
c  sums should zero
c-------------------------
c
c     sum = 0.0
c     Do k=1,Ncl2
c      sum = sum + evt(k)*BLP(k)
c     End Do
c     write (6,112) sum
c
c     Do i=1,Ncl2
c      sum = 0.0
c      Do j=1,Ncl2
c       sum = sum+ALP(j,i)*evt(j)
c      End Do
c      write (6,112) sum
c     End Do
c-------------------------

c----------
c Display
c
c    Do I=1,Ncl2
c      write (6,101) (ALP(I,J),J=1,Ncl2),BLP(i)
c    End Do
c---

c----------------
c Put ALP into AL
c     BLP into BL
c----------------

      Do i=1,Ncl2
       Do j=1,Ncl2
        AL(i,j) = ALP(i,j)
       End Do
       BL(i) = BLP(i)
      End Do

c-----
c Done
c-----

 112  Format (80(1x,f10.7))

      Return
      End

c===========================================

      subroutine reduce ()

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------------------
c Reduce the primary or preconditioned system
c
c Replace the last equation of each particle block
c with the constraint: Int f . n dl = 1.0
c
c SYMBOLS:
c -------
c
c AL:	Master Matrix
c BL:	Master rhs
c
c--------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(72),Itp(72)
      Dimension vnx0(4608),vny0(4608),elml(4608)

c     Dimension  AL(4608,4608),BL(4608)
      Dimension  AL(2304,2304),BL(2304)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl
      common/colloc2/vnx0,vny0,elml
      common/matrix/AL,BL

c--------
c prepare
c-------

      write (6,*) " reduce: Reducing the MLS"

      Ncl2 = 2*Ncl

c-----------
c initialize
c-----------

      Ir  = 0     ! row pointer
      Ic  = 0     ! column counter
      Icp = 0     ! collocation point counter

c-------
c reduce
c-------

      Do i=1,Nprtcl

       Ir = Ir+2*NE(i)

       Do j=1,Ncl2
        AL(Ir,j) = 0.0D0
       End Do

       BL(Ir) = 1.0D0

       Do j=1,NE(i)
        Icp = Icp+1
        Icj = Ic+j
        AL(Ir,Icj) = vnx0(Icp)*elml(Icp)
        Icj = Icj+NE(i)
        AL(Ir,Icj) = vny0(Icp)*elml(Icp)
       End Do

       Ic = Ic + 2*NE(i)

      End Do

c-----
c Done
c-----

      Return
      End
