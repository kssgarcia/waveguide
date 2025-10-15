      subroutine precondition ()

c-----------------------------------------
c FDLIB and  BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c-------------------------------------------------
c Preconditioning the Master Linear Matrix
c
c LEGEND:
c -------
c
c Ncl:  Number of collocation points
c
c ev:  Eigenvectors of the influence matrix
c evt: Eigenvectors of the transpose
c      of the influence matrix
c
c ALP:  Preconditioned Matrix
c BLP   Preconditioned rhs
c
c NOTES:
c -----
c
c For straight segments, the eigenvectors are exact
c
c For native parametrization of an ellipse,
c the eigenvectors are approximate
c
c Capacity:
c ---------
c
c   25 particles
c   64 elements
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(25),Itp(25)
 
      Dimension vnx0(3200),vny0(3200),elar(3200)
      Dimension   ev(3200), evt(3200)

      Dimension  AL(3200,3200), BL(3200)
      Dimension ALP(1600,1600),BLP(1600)

      Dimension hold(25),holdt(25)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NGL,NE,Itp,Ncl
      common/colloc2/vnx0,vny0,elar
      common/matrix/AL,BL

      write (6,*) " preconditioning: started"

c--------
c prepare
c--------

      Ncl2 = 2*Ncl

c----------------------------------------- 
c initialize preconditioned matrix and rhs
c----------------------------------------- 

      Do i=1,Ncl2

       BLP(i) = BL(i)

       Do j=1,Ncl2
        ALP(i,j) = AL(i,j)
       End Do

      End Do

c---------------------------------------------
c compute eigenvectors, transpose eigenvectors
c and hold their norms
c---------------------------------------------

      Ic = 0         ! collocation point counter
      Ib = 0         ! block counter

      Do i=1,Nprtcl

        hold (i) = 0.0D0     ! for normalization
        holdt(i) = 0.0D0     ! for normalization

        Do j=1,NE(i)

          Ic = Ic+1

          Ibj  = Ib  + j
          Ibji = Ibj + NE(i)

          ev(Ibj)  = vnx0(Ic)
          ev(Ibji) = vny0(Ic)

          hold(i) = hold(i) + ev(Ibj)**2+ ev(Ibji)**2

          evt(Ibj)  = vnx0(Ic)*elar(Ic)
          evt(Ibji) = vny0(Ic)*elar(Ic)

          holdt(i)  = holdt(i) + evt(Ibj)**2+evt(Ibji)**2

        End Do

        Ib = Ib+2*NE(i)

        hold (i) = dsqrt(hold (i))
        holdt(i) = dsqrt(holdt(i))

      End Do

c------------------------
c  normalize eigenvectors
c------------------------

      Ic = 0
      Ib = 0

      Do i=1,Nprtcl
       Do j=1,NE(i)

        Ic = Ic+1

        Ibj  = Ib  + j
        Ibji = Ibj + NE(i)

        ev(Ibj)  =  ev(Ibj) /hold(i)
        ev(Ibji) =  ev(Ibji)/hold(i)

        evt(Ibj)  = evt(Ibj) /holdt(i)
        evt(Ibji) = evt(Ibji)/holdt(i)
       End Do

       Ib = Ib+2*NE(i)

      End Do

c-------------------------
c  sum1 and sum2 should be equal to Nprtcl
c
c     sum1 = 0.0
c     sum2 = 0.0
c     Do i=1,Ncl2
c      write (6,112) ev(i),evt(i)
c      sum1 = sum1 + ev(i)**2
c      sum2 = sum2 + evt(i)**2
c     End Do
c     write (6,*)
c     write (6,*) "Check in preconditioning"
c     write (6,*) "Both should be equal to", Nprtcl," :",sum1,sum2
c     write (6,*)
c---------------------------

c---------------------------------
c Test the discrete eigenvector
c
c sum should be very close to zero
c---------------------------------
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
c-----------------------------------

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
c---------------------------------

c-------------------------------
c precondition the linear system
c by projecting it onto the space
c that is normal to the transpose eigenvectors
c
c P indicates preconditioning
c-------------------------------

      Ib = 0

      Do 19 i=1,Nprtcl

       Min = Ib+1          ! lower limit of ith partition
       Max = Ib+2*NE(i)    ! upper limit of ith partition

c---------------------------------
c precondition the right-hand side
c
c sum is the projection of the adjoint eigenvector
c on the right-hand side
c---------------------------------

       sum = 0.0
       Do k=Min,Max
        sum = sum + evt(k)*BLP(k)
       End Do

c      write (6,*) "precondition: projection of adj egv on the rhs:",sum

       Do k=Min,Max
         BLP(k) = BLP(k)-sum*evt(k)
       End Do

c------------------------
c precondition the matrix
c------------------------

       Do j=1,Ncl2

        sum = 0.0

        Do k=Min,Max
         sum = sum+evt(k)*ALP(k,j)
        End Do

c       write (6,*) "precondition: projection of adj egv on jth column:",sum

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
c----------------------------

c----------
c Display
c
c    Do I=1,Ncl2
c      write (6,101) (ALP(I,J),J=1,Ncl2),BLP(i)
c    End Do
c---

c------
c Put ALP into AL
c     BLP into BL
c------

      Do i=1,Ncl2
       Do j=1,Ncl2
        AL(i,j) = ALP(i,j)
       End Do
       BL(i) = BLP(i)
      End Do

c-----
c Done
c-----

      write (6,*) " precondition: exiting"

 112  Format (80(1x,f10.7))

      Return
      End
