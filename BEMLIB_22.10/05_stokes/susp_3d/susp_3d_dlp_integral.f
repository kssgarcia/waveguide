      subroutine susp_3d_dlp_integral
     +
     +     (x0,y0,z0
     +     ,mpoly
     +     ,pcl
     +     ,k
     +     ,intm
     +     ,mint
     +     ,cx,cy,cz
     +     ,necl,ngcl
     +     ,q
     +     ,qx0,qy0,qz0
     +     ,ptlx,ptly,ptlz
     +     ,area
     +     ,qsintx,qsinty,qsintz
     +     ,qsinmx,qsinmy,qsinmz
     +     ,qdn
     +     ,stresslet
     +     )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c--------------------------------------
c Compute the double-layer laplace potential
c over a non-singular triangle numbered k
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),beta(512),gamma(512)

      Dimension xiq(20),etq(20),wq(20)

      Dimension stresslet(3,3)

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

      qsintx = 0.0D0
      qsinty = 0.0D0
      qsintz = 0.0D0

      qsinmx = 0.0D0
      qsinmy = 0.0D0
      qsinmz = 0.0D0

      qdn = 0.0D0

      Do i=1,3
       Do j=1,3
         stresslet(i,j)=0.0D0
       End Do
      End Do

      ptlx = 0.0D0
      ptly = 0.0D0
      ptlz = 0.0D0

      area = 0.0D0

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

      end if

c---
c base points and weights
c---

      Do 1 i=1,mintt    ! quadrature

        If(intm.eq.1) then
          xi  = xiq(i)    ! quadrature
          eta = etq(i)     ! quadrature
        Else If(intm.eq.2) then
          xi  = xir(i)    ! quadrature
          eta = etr(i)     ! quadrature
        End If

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

c        call sgf_3d_w
c    +
c    +   (Iopt_sgf
c    +   ,x,y,z
c    +   ,x0,y0,z0
c    +   ,wall
c    +   ,Gxx,Gxy,Gxz
c    +   ,Gyx,Gyy,Gyz
c    +   ,Gzx,Gzy,Gzz
c    +   ,px,py,pz
c    +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
c    +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
c    +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
c    +   )

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

c     write (6,*) ew,tau
c     write (6,*) a11,a12,a13
c     write (6,*) a21,a22,a23
c     write (6,*) a31,a32,a33
c     write (6,*) b11,b12,b13
c     write (6,*) b21,b22,b23
c     write (6,*) b31,b32,b33
c     write (6,*) Max1,Max2
c     stop

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
        write (6,*) " susp_3d_dlp: flow not implemented"
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

c   constuct q by Lagrange intepolation:

        qintx = 0.0D0
        qinty = 0.0D0
        qintz = 0.0D0

        Do icl=1,necl
         qintx = qintx + psicl(icl)*q(ncl(k,icl),1)
         qinty = qinty + psicl(icl)*q(ncl(k,icl),2)
         qintz = qintz + psicl(icl)*q(ncl(k,icl),3)
c        write (6,*) psicl(icl),q(ncl(k,icl),2)
        End Do

c---
c apply the triangle quadrature
c---

      If(intm.eq.1) then        ! quadrature
       cf = 0.5D0*hs*wq(i)   
      Else If(intm.eq.2) then
       cf = 0.5D0*hs*wr(i)     ! integration rule
      End If

      area = area + cf

      ptlx = ptlx +((qintx-qx0)*(Txxx*vnx+Txxy*vny+Txxz*vnz)
     +            + (qinty-qy0)*(Tyxx*vnx+Tyxy*vny+Tyxz*vnz)
     +            + (qintz-qz0)*(Tzxx*vnx+Tzxy*vny+Tzxz*vnz))*cf

      ptly = ptly +((qintx-qx0)*(Txyx*vnx+Txyy*vny+Txyz*vnz)
     +            + (qinty-qy0)*(Tyyx*vnx+Tyyy*vny+Tyyz*vnz)
     +            + (qintz-qz0)*(Tzyx*vnx+Tzyy*vny+Tzyz*vnz))*cf

      ptlz = ptlz +((qintx-qx0)*(Txzx*vnx+Txzy*vny+Txzz*vnz)
     +            + (qinty-qy0)*(Tyzx*vnx+Tyzy*vny+Tyzz*vnz)
     +            + (qintz-qz0)*(Tzzx*vnx+Tzzy*vny+Tzzz*vnz))*cf

c     write (6,*) qintx,qinty,qintz,qx0,qy0,qz0
c     write (6,*) ptlx,ptly,ptlz

      qdn = qdn + (qintx*vnx + qinty*vny + qintz*vnz)*cf

      qsintx = qsintx + qintx*cf
      qsinty = qsinty + qinty*cf
      qsintz = qsintz + qintz*cf

      qsinmx = qsinmx + ((y-cy)*qintz-qinty*(z-cz))*cf
      qsinmy = qsinmy + ((z-cz)*qintx-qintz*(x-cx))*cf
      qsinmz = qsinmz + ((x-cx)*qinty-qintx*(y-cy))*cf

      stresslet(1,1) = stresslet(1,1) + qintx*vnx*cf
      stresslet(1,2) = stresslet(1,2) + qintx*vny*cf
      stresslet(1,3) = stresslet(1,3) + qintx*vnz*cf

      stresslet(2,1) = stresslet(2,1) + qinty*vnx*cf
      stresslet(2,2) = stresslet(2,2) + qinty*vny*cf
      stresslet(2,3) = stresslet(2,3) + qinty*vnz*cf

      stresslet(3,1) = stresslet(3,1) + qintz*vnx*cf
      stresslet(3,2) = stresslet(3,2) + qintz*vny*cf
      stresslet(3,3) = stresslet(3,3) + qintz*vnz*cf

  1   Continue    ! end of loop over integration points

c-----
c Done
c-----

 100  Format (100(1x,f10.5))

      return
      end
