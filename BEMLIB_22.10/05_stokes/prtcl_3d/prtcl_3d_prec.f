      subroutine precondition
     +
     +  (nelm
     +  ,Mdim,nrhs
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------
c Compute the eigenvector
c         and the transpose eigenvector
c         of the influence matrix
c
c and precondition the linear system
c
c SYMBOLS:
c -------
c
c  vnx0, vny0, vnz0:  
c
c       normal vector
c       at collocation points
c
c  arel(i)    surface area of ith element
c
c  eigen:     Approximate eigenvector of the discrete system
c  eigent:    Approximate eigenvector of the transpose
c             of the discrete system
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension arel(512)

      Dimension vnx0(512),vny0(512),vnz0(512)

      Dimension eigen (512),resdl (512)
      Dimension eigent(512),resdlt(512)

      Dimension RM(1600,1600),rhs(1600,10)

c--------------
c common blocks
c--------------

      common/geo1/arel
      common/geo2/vnx0,vny0,vnz0

      common/prec/eigen,eigent,resdl,resdlt
      common/sys1/rm,rhs

c--------------------------------------------
c The normal vector is an eigenvector of the
c single-layer operator corresponding to the
c zero eigenvalue.
c
c Compute and normalize the approximate eigenvector
c of the matrix and its transpose
c
c Compute the residual : resdl = Matrix*eigen
c--------------------------------------------

      rnorm  = 0.0D0    ! for normalization
      rnormt = 0.0D0    ! for normalization

      Do i=1,nelm

       inelm  = i+nelm
       inelm2 = i+nelm+nelm

       eigen(i)      = vnx0(i)
       eigen(inelm)  = vny0(i)
       eigen(inelm2) = vnz0(i)

       eigent(i)      = vnx0(i)*arel(i)
       eigent(inelm)  = vny0(i)*arel(i)
       eigent(inelm2) = vnz0(i)*arel(i)

       rnorm = rnorm +eigen (i)**2+eigen (inelm)**2
     +                            +eigen (inelm2)**2
       rnormt= rnormt+eigent(i)**2+eigent(inelm)**2
     +                            +eigent(inelm2)**2

      End Do

      rnorm  = sqrt(rnorm)
      rnormt = sqrt(rnormt)

c---------------------------
c normalize the eigenvectors
c---------------------------

      Do i=1,Mdim
       eigen (i) = eigen (i)/rnorm
       eigent(i) = eigent(i)/rnormt
      End Do

c----------------------
c compute the residuals
c----------------------

      Do i=1,Mdim

        resdl (i) = 0.0D0
        resdlt(i) = 0.0D0

        Do j=1,Mdim
         resdl (i) = resdl (i) + rm(i,j)*eigen (j)
         resdlt(i) = resdlt(i) + rm(j,i)*eigent(j)
        End Do

      End Do

c-------------------
c a room with a view
c-------------------

      write (6,*)
      write (6,*) " Display discrete eigenvectors and residual ?"
      write (6,*)
      write (6,*) " Enter 0 for NO, 1 for YES"
      write (6,*) " -------------------------"
      write (6,*)
      read  (5,*) Isee

      if(Isee.eq.1) then

        Do i=1,Mdim
          write (6,100) i,eigen(i),resdl(i),eigent(i),resdlt(i)
        End Do

      end if

c-----------------------------------------------
c Make the coefficient matrix singular
c by premultiplying the linear system
c by the projection matrix:
c
c        I-eigent*eigent
c-----------------------------------------------

      Do i=1,Mdim
        Do j=1,Mdim
          rm(i,j) = rm(i,j)-eigent(i)*resdlt(j)
        End Do
      End Do

      Do k=1,nrhs
       sum = 0.0
       Do j=1,Mdim
        sum = sum+rhs(j,k)*eigent(j)
       End Do
       Do j=1,Mdim
        rhs(j,k) = rhs(j,k)-eigent(j)*sum
       End Do
      End Do

c-----
c done
c-----

  100 Format(1x,i4,10(1x,f12.5))

      Return
      End
