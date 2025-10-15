      subroutine vel
     +
     +  (npts
     +  ,nelm
     +  ,mint
     +  ,NGL
     +  ,Idfl
     +  ,Iflow
     +  ,Isym_xy
     +  ,nter
     +  ,tol
     +  ,Istop
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

c-------------------------------------------
c Solve an integral equation of the
c second kind for the interfacial velocity
c
c SYMBOLS:
c -------
c
c Isym_xy = 0: velocity is computed at all nodes
c Isym_xy = 1: velocity is computed only at nodes that lie
c              on the right of the xy plane
c              The velocity at the other nodes is obtained
c              by reflection
c
c slp:   single-layer potential
c dlp:   double-layer potential
c
c vna:   unit normal vector at the nodes
c-----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension    u(1026,3)
      Dimension uinf(1026,3)
      Dimension upvt(1026,3)

      Dimension nvel(1026),lxy(1026,2)

      Dimension slp(1026,3)
      Dimension dlp(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/viscr/vs1,vs2,vsrt,vsrtp,vsrtm,vsf,vsk
      common/visci/Ivs

      common/var/shrt,wall

      common/veloc0/u
      common/veloc1/nvelt,nvel
      common/veloc2/nvelr,lxy

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c---------------------
c unperturbed velocity:
c
c simple shear flow: ux = k y
c---------------------

       Do i=1,npts
        uinf(i,1) = shrt*p(i,2)
        uinf(i,2) = 0.0D0
        uinf(i,3) = 0.0D0
c       write (6,100) i,uinf(i,1),uinf(i,2),uinf(i,3)
       End Do

c---------------------------------------------
c contribution from the single-layer potential
c---------------------------------------------

      write (6,*) " vel: entering slp"

      call sslp_3d
     +
     +  (nelm
     +  ,npts
     +  ,mint
     +  ,NGL
     +  ,Iflow
     +  ,slp
     +  )

      write (6,*) " vel: exited slp"

c--------------------
c reflect a symmetric interface
c--------------------

      if(Isym_xy.eq.1) then

        Do node=1,nvelr
         i = lxy(node,1)
         j = lxy(node,2)
         slp(i,1) =  slp(j,1)
         slp(i,2) =  slp(j,2)
         slp(i,3) = -slp(j,3)
        End Do

      end if

c--------
c inspect
c--------

c     Do i=1,npts
c      write (6,100) i,slp(i,1),slp(i,2),slp(i,3)
c     End Do

c-----------------------------
c no iterations for lamdba = 1
c-----------------------------

      if(Ivs.eq.0) then

       Do i=1,npts
        u(i,1) = uinf(i,1) + slp(i,1)/vs1
        u(i,2) = uinf(i,2) + slp(i,2)/vs1
        u(i,3) = uinf(i,3) + slp(i,3)/vs1
       End Do

       Go to 99

      end if

c-----------------------------------------
c viscosity ratio is different than 1
c continue with the double-layer potential
c-----------------------------------------

      Do i=1,npts

       upvt(i,1) = vsf*( uinf(i,1) + slp(i,1)/vs1 )
       upvt(i,2) = vsf*( uinf(i,2) + slp(i,2)/vs1 )
       upvt(i,3) = vsf*( uinf(i,3) + slp(i,3)/vs1 )

      End Do

c---
c begin the Neumann iterations
c---

      Istop = 0

      iter  = 1    ! iterations counter

 96   Continue

      CF1 = 0.0D0    ! deflation coefficient

c---
c deflate one eigenvalue
c---

      if(Idfl.ne.0) then

        call deflation  ! to compute CF1
     +
     +   (Idfl
     +   ,npts
     +   ,nelm
     +   ,mint
     +   ,vsk
     +   ,u
     +   ,CF1
     +   )

      end if

c---
c compute the double-layer potential
c---

      call sdlp_3d 
     +
     + (npts
     + ,nelm
     + ,mint
     + ,Iflow
     + ,u
     + ,dlp
     + )

c---
c reflect a symmetric interface
c---

      if(Isym_xy.eq.1) then

        Do node = 1,nvelr
         i = lxy(node,1)
         j = lxy(node,2)
         dlp(i,1) =  dlp(j,1)
         dlp(i,2) =  dlp(j,2)
         dlp(i,3) = -dlp(j,3)
        End Do

      end if

c---
c Gauss--Siedel updating
c---

      Diff = 0.0D0

      Do i=1,npts

       unew = upvt(i,1) + vsk*dlp(i,1) - CF1*vna(i,1)
       vnew = upvt(i,2) + vsk*dlp(i,2) - CF1*vna(i,2)
       wnew = upvt(i,3) + vsk*dlp(i,3) - CF1*vna(i,3)

       Dev  = Dsqrt((unew-u(i,1))**2     ! maximum correction
     +            + (vnew-u(i,2))**2
     +            + (wnew-u(i,3))**2)

       if(Dev.gt.Diff) Diff = Dev

       u(i,1) = unew
       u(i,2) = vnew
       u(i,3) = wnew

      End Do

      if(Idfl.eq.0) write (6,107) iter,Diff
      if(Idfl.eq.1) write (6,108) iter,Diff,CF1

c---------------------
c stop the iterations ?
c---------------------

      if(iter.gt.nter) then

       write (6,*) " vel: iterations failed to converge"

       Istop = 1
       Go to 99

      end if

      if(diff.gt.tol) then
       iter = iter+1
       Go to 96
      end if

c-----
c done
c-----

  99  Continue

 100  Format (1x,i3,10(f12.8))
 107  Format (" iter:",I3," Max corr =",f15.10)
 108  Format (" iter:",I3," Max corr =",f15.10," CF1 =",f15.10)

      return
      end
