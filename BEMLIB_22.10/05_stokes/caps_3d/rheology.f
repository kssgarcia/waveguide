      subroutine rheology 
     +
     + (nelm
     + ,npts
     + ,mint
     + ,vlm
     + ,u
     + ,sxy,sd1,sd2
     + )

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c----------------------------------------
c Computation of the effective stress tensor,
c shear stress, first and second normal stress
c differences of a dilute suspension
c
c see:   Kennedy et al. (1994), eq (8)
c        or Pozrikidis (1992) p. 143
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)
      Dimension   u(1026,3)
      Dimension Df(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)
      Dimension  Dfel(512,4)

      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/dff/Df
      common/dffel/Dfel

      common/viscr/vs1,vs2,vsrt,vsrtp,vsrtm,vsf,vsk

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c-----------
c initialize
c-----------
      
      strxx = 0.0D0
      strxy = 0.0D0
      strxz = 0.0D0

      stryx = 0.0D0
      stryy = 0.0D0
      stryz = 0.0D0

      strzx = 0.0D0
      strzy = 0.0D0
      strzz = 0.0D0

c-----------------------------
c integrate over the interface
c-----------------------------

      Do 1 k=1,nelm

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be =  beta(k)
      ga = gamma(k)

      sxx = 0.0D0
      sxy = 0.0D0
      sxz = 0.0D0

      syx = 0.0D0
      syy = 0.0D0
      syz = 0.0D0

      szx = 0.0D0
      szy = 0.0D0
      szz = 0.0D0

      Do 2 i=1,mint

        xi  = xiq(i)
        eta = etq(i)

c------------------------------
c compute df and u
c at triangle quadrature points
c------------------------------

        call rheology_interp 
     +
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
     +
     +  ,df(i1,1),df(i1,2),df(i1,3)
     +  ,df(i2,1),df(i2,2),df(i2,3)
     +  ,df(i3,1),df(i3,2),df(i3,3)
     +  ,df(i4,1),df(i4,2),df(i4,3)
     +  ,df(i5,1),df(i5,2),df(i5,3)
     +  ,df(i6,1),df(i6,2),df(i6,3)
     +
     +  ,vna(i1,1),vna(i1,2),vna(i1,3)
     +  ,vna(i2,1),vna(i2,2),vna(i2,3)
     +  ,vna(i3,1),vna(i3,2),vna(i3,3)
     +  ,vna(i4,1),vna(i4,2),vna(i4,3)
     +  ,vna(i5,1),vna(i5,2),vna(i5,3)
     +  ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +  ,u(i1,1),u(i1,2),u(i1,3)
     +  ,u(i2,1),u(i2,2),u(i2,3)
     +  ,u(i3,1),u(i3,2),u(i3,3)
     +  ,u(i4,1),u(i4,2),u(i4,3)
     +  ,u(i5,1),u(i5,2),u(i5,3)
     +  ,u(i6,1),u(i6,2),u(i6,3)
     +
     +  ,al,be,ga
     +  ,xi,eta
     +  ,x,y,z
     +  ,dfx,dfy,dfz
     +  ,vnx,vny,vnz
     +  ,ux,uy,uz
     +  ,hs
     +  )

      cf = 0.5D0*hs*wq(i)

      sxx = sxx + (dfx*x - vsrtm*(ux*vnx+ux*vnx))*cf
      sxy = sxy + (dfx*y - vsrtm*(ux*vny+uy*vnx))*cf
      sxz = sxz + (dfx*z - vsrtm*(ux*vnz+uz*vnx))*cf

      syx = syx + (dfy*x - vsrtm*(uy*vnx+ux*vny))*cf
      syy = syy + (dfy*y - vsrtm*(uy*vny+uy*vny))*cf
      syz = syz + (dfy*z - vsrtm*(uy*vnz+uz*vny))*cf
 
      szx = szx + (dfz*x - vsrtm*(uz*vnx+ux*vnz))*cf
      szy = szy + (dfz*y - vsrtm*(uz*vny+uy*vnz))*cf
      szz = szz + (dfz*z - vsrtm*(uz*vnz+uz*vnz))*cf

  2   Continue

      strxx = strxx + sxx
      strxy = strxy + sxy
      strxz = strxz + sxz

      stryx = stryx + syx
      stryy = stryy + syy
      stryz = stryz + syz

      strzx = strzx + szx
      strzy = strzy + szy
      strzz = strzz + szz

  1   Continue

c---
c divide by drop volume
c---

c     strxx = strxx/vlm
c     strxy = strxy/vlm
c     strxz = strxz/vlm
c     stryx = stryx/vlm
c     stryy = stryy/vlm
c     stryz = stryz/vlm
c     strzx = strzx/vlm
c     strzy = strzy/vlm
c     strzz = strzz/vlm
c     write (6,100) strxx,strxy,strxz
c     write (6,100) stryx,stryy,stryz
c     write (6,100) strzx,strzy,strzz

c---
c compute normal stress differences
c---

      sxy = strxy
      sd1 = strxx-stryy
      sd2 = stryy-strzz

c-----
c Done
c-----

 100  Format (3(1x,f15.8))

      return
      end
