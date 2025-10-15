      subroutine susp_3d_dlp_integral_mutual
     +
     +     (x0,y0,z0
     +     ,mpoly
     +     ,pcl
     +     ,k
     +     ,intm
     +     ,mint
     +     ,necl,ngcl
     +     ,q
     +     ,ptlx,ptly,ptlz
     +     ,area
     +     ,tstxx,tstxy,tstxz
     +     ,tstyx,tstyy,tstyz
     +     ,tstzx,tstzy,tstzz
     +     )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c--------------------------------------
c Compute the double-layer laplace potential
c over a non-singular triangle numbered k
c at the point x0,y0,z0
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),beta(512),gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c collocation:

      Dimension vaninv(500,500)
      Dimension    xcl(512,100),ycl(512,100),zcl(512,100)

      Dimension   ncl(512,100)   ! connectivity
      Dimension   pcl(2306,3)
      Dimension     q(2306,3)

      Dimension aux(100),psicl(100)

      Dimension xir(2400),etr(2400),wr(2400)    ! integration weights

c periodic flow:

      Dimension vxxx(-9:9,-9:9,-9:9)
      Dimension vxxy(-9:9,-9:9,-9:9)
      Dimension vxxz(-9:9,-9:9,-9:9)
      Dimension vyxy(-9:9,-9:9,-9:9)
      Dimension vyxz(-9:9,-9:9,-9:9)
      Dimension vzxz(-9:9,-9:9,-9:9)

      Dimension vxyx(-9:9,-9:9,-9:9)
      Dimension vxyy(-9:9,-9:9,-9:9)
      Dimension vxyz(-9:9,-9:9,-9:9)
      Dimension vyyy(-9:9,-9:9,-9:9)
      Dimension vyyz(-9:9,-9:9,-9:9)
      Dimension vzyz(-9:9,-9:9,-9:9)

      Dimension vxzx(-9:9,-9:9,-9:9)
      Dimension vxzy(-9:9,-9:9,-9:9)
      Dimension vxzz(-9:9,-9:9,-9:9)
      Dimension vyzy(-9:9,-9:9,-9:9)
      Dimension vyzz(-9:9,-9:9,-9:9)
      Dimension vzzz(-9:9,-9:9,-9:9)

      Dimension ppx(-9:9,-9:9,-9:9)
      Dimension ppy(-9:9,-9:9,-9:9)
      Dimension ppz(-9:9,-9:9,-9:9)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/conncl/ncl

      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/var1/Iflow,Iwall
      common/var3/wall

      common/spectral2/vaninv

      common/trq/xiq,etq,wq
      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c: periodic flow

      common/aaaa/a11,a12,a13,a21,a22,a23,a31,a32,a33
      common/bbbb/b11,b12,b13,b21,b22,b23,b31,b32,b33
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

      common/vvvv_3d/vxxx,vxxy,vxxz,vyxy,vyxz,vzxz
     +              ,vxyx,vxyy,vxyz,vyyy,vyyz,vzyz
     +              ,vxzx,vxzy,vxzz,vyzy,vyzz,vzzz
     +              ,ppx,ppy,ppz

c----------
c set flags
c----------

      Iopt_int = 2   ! need the surface metric
      Iopt_sgf = 2   ! need T

c-----------
c initialize
c-----------

      ptlx = 0.0D0
      ptly = 0.0D0
      ptlz = 0.0D0

      area = 0.0D0

      tstxx = 0.0D0
      tstxy = 0.0D0
      tstxz = 0.0D0
      tstyx = 0.0D0
      tstyy = 0.0D0
      tstyz = 0.0D0
      tstzx = 0.0D0
      tstzy = 0.0D0
      tstzz = 0.0D0

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

c-----------------------------
c loop over integration points
c-----------------------------

      if(intm.eq.1) then    ! first integration method: quadrature

        mintt=mint

      else if(intm.eq.2) then ! second method: Newton-Cotes

       call trgl_rule (mint,Nbase,xir,etr,wr)

c      write (6,*) Nbase
c      write (6,*)
c      write (6,*) (xir(i),i=1,Nbase)
c      write (6,*)
c      write (6,*) (etr(i),i=1,Nbase)
c      write (6,*)
c      write (6,*) ( wr(i),i=1,Nbase)
c      stop

       mintt=Nbase

      End If

c---
c base points and weights
c---

      Do 1 i=1,mintt    ! quadrature

        if(intm.eq.1) then
          xi  = xiq(i)    ! quadrature
          eta = etq(i)     ! quadrature
        else if(intm.eq.2) then
          xi  = xir(i)    ! quadrature
          eta = etr(i)     ! quadrature
        end if

        if(xi.lt.0.000001)  xi=0.001D0
        if(xi.gt.0.999999)  xi=0.999D0
        if(eta.lt.0.000001) eta=0.001D0
        if(eta.gt.0.999999) eta=0.999D0

        call susp_3d_dlp_interp
     +
     +    (Iopt_int
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
     +    ,al,be,ga
     +    ,xi,eta
     +
     +    ,x,y,z
     +    ,vnx,vny,vnz
     +    ,hs
     +    )

c---
c compute the Green's function
c---

c------
      if(Iflow.lt.10) then
c------

       call sgf_3d_fs
     +
     +    (Iopt_sgf
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,Gxx,Gxy,Gxz
     +    ,Gyx,Gyy,Gyz
     +    ,Gzx,Gzy,Gzz
     +    ,px,py,pz
     +    ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +    ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +    ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +    )

      Tyxx = Txxy
      Tzxx = Txxz
      Tzxy = Tyxz

      Tyyx = Txyy
      Tzyx = Txyz
      Tzyy = Tyyz

      Tyzx = Txzy
      Tzzx = Txzz
      Tzzy = Tyzz

c------
      else if(Iflow.eq.10) then
c------

       call sgf_3d_w
     +
     +  (Iopt_sgf
     +  ,y,z,x
     +  ,y0,z0,x0
     +  ,wall
     +  ,gyy,gyz,gyx
     +  ,gzy,gzz,gzx
     +  ,gxy,gxz,gxx
     +  ,py,pz,px
     +  ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +  ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +  ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +  )

      Txyy = Tyyx
      Txyz = Tzyx
      Tzyy = Tyyz

      Txzy = Tyzx
      Txzz = Tzzx
      Tzzy = Tyzz

      Txxy = Tyxx
      Txxz = Tzxx
      Tzxy = Tyxz

c------
      else if(Iflow.eq.20) then
c------

      call sgf_3d_3p
     +
     +   (Iopt_sgf
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,a11,a12,a13
     +   ,a21,a22,a23
     +   ,a31,a32,a33
     +   ,b11,b12,b13
     +   ,b21,b22,b23
     +   ,b31,b32,b33
     +   ,ew,tau
     +   ,Max1,Max2
     +   ,Gxx,Gxy,Gxz
     +   ,Gyx,Gyy,Gyz
     +   ,Gzx,Gzy,Gzz
     +   ,px,py,pz
     +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +   )

      Tyxx = Txxy
      Tzxx = Txxz
      Tzxy = Tyxz

      Tyyx = Txyy
      Tzyx = Txyz
      Tzyy = Tyyz

      Tyzx = Txzy
      Tzzx = Txzz
      Tzzy = Tyzz

c------
      else
c------

        write (6,*)
        write (6,*) " susp_3d_dlp: This flow has not yet implemented"
        stop

c------
      end if
c------

c        If(abs(G).gt.100000.0) then
c        write (6,100) x,y,z,x0,y0,z0,G,xi,eta
c        End If

c-----
c interpolate for the dipole
c using the Vandermonde matrix:
c-----

         xip = 2.0D0*xi/(1.D0-eta)-1.0D0
         etap = 2.0D0*eta-1.0D0

         Jc=0
         Do kcl=0,mpoly
          Do lcl=0,mpoly-kcl
            Jc=Jc+1
            call jacobi (0D0,0D0,kcl,xip,rjac1)
            call jacobi (2.0D0*kcl+1.0D0,0.0D0,lcl,etap,rjac2)
            aux(Jc) = rjac1 *(1.0D0-eta)**kcl *rjac2
c           write (6,*) aux(Jc)
           End Do
         End Do

        Do icl=1,necl
         psicl(icl) = 0.0D0
          Do jcl=1,necl
            psicl(icl) = psicl(icl) + vaninv(icl,jcl)*aux(jcl)
c           write (6,*) vaninv(jcl,icl)
          End Do
        End Do

c   construct q by Lagrange intepolation:

        qintx = 0.0D0
        qinty = 0.0D0
        qintz = 0.0D0

        Do icl=1,necl
         qnodex = q(ncl(k,icl),1)
         qnodey = q(ncl(k,icl),2)
         qnodez = q(ncl(k,icl),3)
c         nodex = 1.0D0
c        qnodey = 1.0D0
c        qnodez = 1.0D0
         qintx = qintx + psicl(icl)* qnodex
         qinty = qinty + psicl(icl)* qnodey
         qintz = qintz + psicl(icl)* qnodez
c        write (6,*) k,icl,ncl(k,icl),qnodex
        End Do

c     write (6,120) k,qintx,qinty,qintz

c---
c apply the triangle quadrature
c---

      If(intm.eq.1) then        ! quadrature
       cf = 0.5D0*hs*wq(i)
      Else If(intm.eq.2) then
       cf = 0.5D0*hs*wr(i)     ! integration rule
      End If

      ptlx = ptlx +(qintx*(Txxx*vnx+Txxy*vny+Txxz*vnz)
     +            + qinty*(Tyxx*vnx+Tyxy*vny+Tyxz*vnz)
     +            + qintz*(Tzxx*vnx+Tzxy*vny+Tzxz*vnz))*cf

      ptly = ptly +(qintx*(Txyx*vnx+Txyy*vny+Txyz*vnz)
     +            + qinty*(Tyyx*vnx+Tyyy*vny+Tyyz*vnz)
     +            + qintz*(Tzyx*vnx+Tzyy*vny+Tzyz*vnz))*cf

      ptlz = ptlz +(qintx*(Txzx*vnx+Txzy*vny+Txzz*vnz)
     +            + qinty*(Tyzx*vnx+Tyzy*vny+Tyzz*vnz)
     +            + qintz*(Tzzx*vnx+Tzzy*vny+Tzzz*vnz))*cf

c     write (6,120) k,ptlx,ptly,ptlz

      tstxx = tstxx + (Txxx*vnx+Txxy*vny+Txxz*vnz)*cf
      tstxy = tstxy + (Txyx*vnx+Txyy*vny+Txyz*vnz)*cf
      tstxz = tstxz + (Txzx*vnx+Txzy*vny+Txzz*vnz)*cf
      tstyx = tstyx + (Tyxx*vnx+Tyxy*vny+Tyxz*vnz)*cf
      tstyy = tstyy + (Tyyx*vnx+Tyyy*vny+Tyyz*vnz)*cf
      tstyz = tstyz + (Tyzx*vnx+Tyzy*vny+Tyzz*vnz)*cf
      tstzx = tstzx + (Tzxx*vnx+Tzxy*vny+Tzxz*vnz)*cf
      tstzy = tstzy + (Tzyx*vnx+Tzyy*vny+Tzyz*vnz)*cf
      tstzz = tstzz + (Tzzx*vnx+Tzzy*vny+Tzzz*vnz)*cf

c     write (6,*) tstxx,tstxy,tstxz

      area = area + cf

  1   Continue    ! end of loop over integration points

c-----
c done
c-----

 100  Format (100(1x,f10.5))
 120  Format (i3,1x,10(1x,f15.10))

      return
      end
