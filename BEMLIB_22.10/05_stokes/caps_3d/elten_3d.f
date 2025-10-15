      subroutine elten_3d
     +
     +   (nelm
     +   ,npts
     +   ,elst
     +   ,elten
     +   )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c-----------------------------------------------
c Compute elastic tensions at the
c nodes using shell theory
c
c SYMBOLS:
c -------
c
c elten(i,j,k):   elastic tensions at the ith node
c Itally:  counter for averaging over elements
c elst:   in-plane modulus of elasticity
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)

      Dimension   pr(1026,3)
      Dimension vnar(1026,3)

      Dimension  elten(1026,3,3)
      Dimension    eltenten(3,3)

      Dimension prj(3,3),save(3,3)

      Dimension Itally(1026)     ! internal

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xxi(6), eet(6)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/elten1/pr,vnar

c-----------
c initialize
c-----------

      Do i=1,npts

       Do j=1,3
        Do k=1,3
         elten(i,j,k) = 0.0D0
        End Do
       End Do
 
       Itally(i) = 0 

      End Do

c-------------------
c loop over elements
c-------------------

      Do k=1,nelm

       i1 = n(k,1)     ! global node labels
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       al = alpha(k)
       be =  beta(k)
       ga = gamma(k)

c---------------------
c triangle coordinates
c at the nodes
c---------------------

       xxi(1) = 0.0D0
       eet(1) = 0.0D0

       xxi(2) = 1.0D0
       eet(2) = 0.0D0

       xxi(3) = 0.0D0
       eet(3) = 1.0D0

       xxi(4) = al
       eet(4) = 0.0D0

       xxi(5) = ga
       eet(5) = 1.0D0-ga

       xxi(6) = 0.0D0
       eet(6) = be

c---------------------------------------
c compute the vector tangential gradient
c---------------------------------------

       Do i=1,6    ! loop over element-nodes

        m = n(k,i)   ! global index of local node i on element k

        xi  = xxi(i)
        eta = eet(i)

        call elten_3d_interp
     +
     +    (elst
     +    ,p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +
     +    ,vna(i1,1),vna(i1,2),vna(i1,3)
     +    ,vna(i2,1),vna(i2,2),vna(i2,3)
     +    ,vna(i3,1),vna(i3,2),vna(i3,3)
     +    ,vna(i4,1),vna(i4,2),vna(i4,3)
     +    ,vna(i5,1),vna(i5,2),vna(i5,3)
     +    ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +    ,pr(i1,1),pr(i1,2),pr(i1,3)
     +    ,pr(i2,1),pr(i2,2),pr(i2,3)
     +    ,pr(i3,1),pr(i3,2),pr(i3,3)
     +    ,pr(i4,1),pr(i4,2),pr(i4,3)
     +    ,pr(i5,1),pr(i5,2),pr(i5,3)
     +    ,pr(i6,1),pr(i6,2),pr(i6,3)
     +
     +    ,vnar(i1,1),vnar(i1,2),vnar(i1,3)
     +    ,vnar(i2,1),vnar(i2,2),vnar(i2,3)
     +    ,vnar(i3,1),vnar(i3,2),vnar(i3,3)
     +    ,vnar(i4,1),vnar(i4,2),vnar(i4,3)
     +    ,vnar(i5,1),vnar(i5,2),vnar(i5,3)
     +    ,vnar(i6,1),vnar(i6,2),vnar(i6,3)
     +
     +    ,al,be,ga
     +    ,xi,eta
     +    ,eltenten
     +    )

        Do j=1,3
         Do l=1,3
           elten(m,j,l) = elten(m,j,l) + eltenten(j,l)
          End Do
        End Do

        Itally(m) = Itally(m)+1

      End Do

      End Do ! end of loop over elements

c---------------------------
c Average the elastic tensions
c at the nodes
c---------------------------

      Do i=1,npts

        par = float(Itally(i))

        Do j=1,3
         Do l=1,3
          elten(i,j,l) = elten(i,j,l)/par
         End Do
        End Do

      End Do

c--------
c project
c--------

      Do i=1,npts

      prj(1,1) = 1.0D0 - vna(i,1)*vna(i,1)
      prj(1,2) =       - vna(i,1)*vna(i,2)
      prj(1,3) =       - vna(i,1)*vna(i,3)
      prj(2,2) = 1.0D0 - vna(i,2)*vna(i,2)
      prj(2,3) =       - vna(i,2)*vna(i,3)
      prj(3,3) = 1.0D0 - vna(i,3)*vna(i,3)

      prj(2,1) = prj(1,2)
      prj(3,1) = prj(1,3)
      prj(3,2) = prj(2,3)

      Do j=1,3
       Do l=1,3
        save(j,l) = elten(i,j,l)
       End Do
      End Do

      Do j=1,3
       Do l=1,3
         elten(i,j,l) = 0.0
         Do k=1,3
          elten(i,j,l) = elten(i,j,l) + prj(j,k)*save(k,l)
         End Do
       End Do
      End Do

      Do j=1,3
       Do l=1,3
        save(j,l) = elten(i,j,l)
       End Do
      End Do

      Do j=1,3
       Do l=1,3
         elten(i,j,l) = 0.0
         Do k=1,3
          elten(i,j,l) = elten(i,j,l) + save(j,k)*prj(k,l)
         End Do
       End Do
      End Do

      End Do

c-----------------
c printing session
c-----------------

c     Do i=1,npts
c       testx = vna(i,1)*elten(i,1,1)+vna(i,2)*elten(i,2,1)+
c    +        + vna(i,3)*elten(i,3,1)
c       testy = vna(i,1)*elten(i,1,2)+vna(i,2)*elten(i,2,2)+
c    +        + vna(i,3)*elten(i,3,2)
c       testz = vna(i,1)*elten(i,1,3)+vna(i,2)*elten(i,2,3)+
c    +        + vna(i,3)*elten(i,3,3)
c       write (6,*) i
c       write (6,100) elten(i,1,1),elten(i,1,2),elten(i,1,3)
c       write (6,100) elten(i,2,1),elten(i,2,2),elten(i,2,3)
c       write (6,100) elten(i,3,1),elten(i,3,2),elten(i,3,3)
c       write (6,100) testx,testy,testz
c     End Do

c-----
c Done
c-----

 100  Format (10(1x,f10.5))

      return
      end
