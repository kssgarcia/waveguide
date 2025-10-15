      subroutine lslp3_3d_integrate
     +
     +    (x0,y0,z0
     +    ,ipoint,k
     +    ,mint
     +    ,zeta
     +    ,ptlx,ptly,ptlz
     +    ,area
     +    )

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=============================================

c--------------------------------------------------
c Compute the vector laplace single-layer potential
c over a non-singular triangle numbered k
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension zeta(1026,3)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),beta(512),gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/var1/Iflow,Ign,wall

      common/geo2/vna

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c-----------
c initialize
c-----------

      Iopt_int = 2   ! need the surface metric
      Iopt_lgf = 1   ! need only G  - not grad(G)

      area = 0.0D0
      ptlx = 0.0D0
      ptly = 0.0D0
      ptlz = 0.0D0

c---
c define triangle nodes
c and mapping parameters
c---

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be = beta (k)
      ga = gamma(k)

c---
c loop over integration points
c---

      Do 1 i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        call lslp3_3d_interp
     +  
     +      (Iopt_int
     +      ,p(i1,1),p(i1,2),p(i1,3)
     +      ,p(i2,1),p(i2,2),p(i2,3)
     +      ,p(i3,1),p(i3,2),p(i3,3)
     +      ,p(i4,1),p(i4,2),p(i4,3)
     +      ,p(i5,1),p(i5,2),p(i5,3)
     +      ,p(i6,1),p(i6,2),p(i6,3)
     +      ,zeta(i1,1),zeta(i1,2),zeta(i1,3)
     +      ,zeta(i2,1),zeta(i2,2),zeta(i2,3)
     +      ,zeta(i3,1),zeta(i3,2),zeta(i3,3)
     +      ,zeta(i4,1),zeta(i4,2),zeta(i4,3)
     +      ,zeta(i5,1),zeta(i5,2),zeta(i5,3)
     +      ,zeta(i6,1),zeta(i6,2),zeta(i6,3)
     +      ,al,be,ga
     +      ,xi,eta
     +      ,x,y,z
     +      ,hs
     +      ,zet1,zet2,zet3
     +      )

c---
c compute the Green's function
c---

c--------
       if(Iflow.eq.1) then
c--------

       call lgf_3d_fs 
     +
     +    (Iopt_lgf
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,G
     +    ,Gx,Gy,Gz
     +    )

c--------
       else if(Iflow.eq.2) then
c--------

        call lgf_3d_w
     +
     +   (Iopt_lgf
     +   ,Ign
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,wall
     +   ,G
     +   ,Gx,Gy,Gz
     +   )

c-------
       end if
c-------

c---
c apply triangle quadrature
c---

      area = area+ 0.5*hs*wq(i)

      cf = G*0.5*hs*wq(i)

      ptlx = ptlx + zet1*cf 
      ptly = ptly + zet2*cf 
      ptlz = ptlz + zet3*cf 

  1   Continue      ! end of loop over integration points

c-----
c Done
c-----

 100  Format (4(1x,f10.5))

      return
      end
