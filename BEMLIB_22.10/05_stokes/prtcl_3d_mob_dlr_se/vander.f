      subroutine vander(Idist)

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c---------------------------------------
c Compute the inverse of the generalized
c Vandermonde matrix for:
c
c the lobatto triangle nodes
c the Proriol basis fuctions
c
c m: polynomial order
c
c Idist = 0  uniform distribution
c       = 1  spectral distribution
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension v(8),vmaster(8)

      Dimension Zlob(20),Wlob(20)

      Dimension    van(500,500)
      Dimension vaninv(500,500)

      common/spectral1/mpoly,npoly,vmaster,van,vaninv

c--------------
c triangle grid
c--------------

      m = mpoly

c------------------------
      if(Idist.eq.0) then  ! uniform grid
c------------------------

       Do i=1,m+1
        v(i) = (i-1.0D0)/m
       End Do

c----------------------------
      else if(Idist.eq.1) then   ! spectral grid
c----------------------------

      call lobatto (m+1,Zlob,Wlob)    ! Wlob will not be needed

      Do i=1,m+1
        v(i) = 0.5D0*(1.0D0+Zlob(i))
      End Do

c-----------
      end if
c-----------

      write (6,*) "vander: master grid:"

      Do i=1,m+1
       write (6,100) v(i)
      End Do

      Do i=1,m+1
       vmaster(i) = v(i)
      End Do

c-------------------------------------------
c compute the generalized vandermonde matrix
c-------------------------------------------

      Ic = 0

      Do i=1,m+1        ! run over the nodes
        Do j=1,m+2-i

         Ic = Ic+1

         k = m+3-i-j
         xi  = (1.0D0 + 2.0D0*v(i)  -v(j)-v(k))/3.0D0
         eta = (1.0D0 -   v(i)+2.0D0*v(j)-v(k))/3.0D0

c        write (6,*) xi,eta

         if(eta.gt.0.999999) then
          eta=0.999999
         end if

         xip = 2.0D0*xi/(1.0D0-eta)-1.0D0
         etap = 2.0D0*eta-1.0D0

         Jc = 0
         Do k=0,m
          Do l=0,m-k
            Jc = Jc+1
            call jacobi(0D0,0D0,k,xip,rjac1)
            call jacobi(2.0D0*k+1.0D0,0.0D0,l,etap,rjac2)
            van(Jc,Ic) = rjac1 *(1.0D0-eta)**k *rjac2
c           write (6,*) xi,eta,rjac1,rjac2
          End Do
        End Do

        End Do
      End Do

      npoly = Ic

c     Do i=1,npoly
c       write (6,*) (van(i,j),j=1,nlob)
c     End Do
c     pause

c----
c compute the inverse
c----

      Isym = 0
      Iwlpvt = 0

      call gel_inv
     +
     +  (npoly
     +  ,van
     +  ,vaninv
     +  ,Isym,Iwlpvt
c    +  ,l,u
c    +  ,det
     +  ,Istop
     +  )

c----
c printing session
c----

c      write (6,*) npoly
c      Do i=1,npoly
c        write (6,102) (van(i,j),j=1,npoly)
c      End Do
c      write (6,*) npoly
c      Do i=1,npoly
c        write (6,102) (vaninv(i,j),j=1,npoly)
c      End Do
c      pause

c-----
c done
c-----
 
 100  Format (1x,f10.5)

      return
      end
