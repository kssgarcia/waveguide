      subroutine caps_3d_slp_integrate
     +
     +   (x0,y0,z0
     +   ,k
     +   ,mint
     +   ,Iflow
     +   ,uxel,uyel,uzel
     +   )

c======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c======================================

c-----------------------------------
c Compute the single-layer potential
c over the non-singular triangle
c numbered k
c-----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension srtn(1026)
      Dimension   Df(1026,3)

      Dimension     n(512,6),nbe(512,3)
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

      common/tension/srtn

      common/dff/Df
      common/dffel/Dfel

      common/var/shrt,wall

      common/trq/xiq,etq,wq

c---
c for triply-periodic flow
c---

      common/aaaa/a11,a12,a13,a21,a22,a23,a31,a32,a33
      common/bbbb/b11,b12,b13,b21,b22,b23,b31,b32,b33
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

c-----------
c initialize
c-----------

      Iopt_int = 1   ! need the surface metric
      Iopt_sgf = 1   ! only G is needed

c----------------------
c launch the quadrature
c----------------------

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be =  beta(k)
      ga = gamma(k)

      dfx1 = df(i1,1) + Dfel(k,1)
      dfy1 = df(i1,2) + Dfel(k,2)
      dfz1 = df(i1,3) + Dfel(k,3)

      dfx2 = df(i2,1) + Dfel(k,1)
      dfy2 = df(i2,2) + Dfel(k,2)
      dfz2 = df(i2,3) + Dfel(k,3)

      dfx3 = df(i3,1) + Dfel(k,1)
      dfy3 = df(i3,2) + Dfel(k,2)
      dfz3 = df(i3,3) + Dfel(k,3)

      dfx4 = df(i4,1) + Dfel(k,1)
      dfy4 = df(i4,2) + Dfel(k,2)
      dfz4 = df(i4,3) + Dfel(k,3)

      dfx5 = df(i5,1) + Dfel(k,1)
      dfy5 = df(i5,2) + Dfel(k,2)
      dfz5 = df(i5,3) + Dfel(k,3)

      dfx6 = df(i6,1) + Dfel(k,1)
      dfy6 = df(i6,2) + Dfel(k,2)
      dfz6 = df(i6,3) + Dfel(k,3)

      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        call caps_3d_slp_interp
     +
     +     (Iopt_int
     +
     +     ,p(i1,1),p(i1,2),p(i1,3)
     +     ,p(i2,1),p(i2,2),p(i2,3)
     +     ,p(i3,1),p(i3,2),p(i3,3)
     +     ,p(i4,1),p(i4,2),p(i4,3)
     +     ,p(i5,1),p(i5,2),p(i5,3)
     +     ,p(i6,1),p(i6,2),p(i6,3)
     +
     +     ,dfx1,dfy1,dfz1
     +     ,dfx2,dfy2,dfz2
     +     ,dfx3,dfy3,dfz3
     +     ,dfx4,dfy4,dfz4
     +     ,dfx5,dfy5,dfz5
     +     ,dfx6,dfy6,dfz6
     +
     +     ,al,be,ga
     +     ,xi,eta
     +     ,x,y,z
     +     ,hs
     +     ,dfx,dfy,dfz
     +
     +     ,vna(i1,1),vna(i1,2),vna(i1,3)
     +     ,vna(i2,1),vna(i2,2),vna(i2,3)
     +     ,vna(i3,1),vna(i3,2),vna(i3,3)
     +     ,vna(i4,1),vna(i4,2),vna(i4,3)
     +     ,vna(i5,1),vna(i5,2),vna(i5,3)
     +     ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +     ,srtn(i1),srtn(i2),srtn(i3)
     +     ,srtn(i4),srtn(i5),srtn(i6)
     +
     +     ,vnx,vny,vnz
     +     ,crvmn
     +     ,tn
     +     ,gradtnx,gradtny,gradtnz
     +     )

c-----------------------------
c compute the Green's function
c-----------------------------

c-------------------------
       if(Iflow.eq.1) then
c-------------------------

       call sgf_3d_fs 
     +
     +    (Iopt_sgf
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,gxx,gxy,gxz
     +    ,gyx,gyy,gyz
     +    ,gzx,gzy,gzz
     +    ,presx,presy,presz
     +    ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +    ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +    ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +    )

c------------------------------
       else if(Iflow.eq.2) then
c------------------------------

       call sgf_3d_w 
     +
     +      (Iopt_sgf
     +      ,y,z,x
     +      ,y0,z0,x0
     +      ,wall
     +      ,gyy,gyz,gyx
     +      ,gzy,gzz,gzx
     +      ,gxy,gxz,gxx
     +      ,presy,presz,presx
     +      ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +      ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +      ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +      )


c-----------------------------
      else if(Iflow.eq.5) then
c-----------------------------

      call sgf_3d_3p
     +
     +      (Iopt_sgf
     +      ,x,y,z
     +      ,x0,y0,z0
     +      ,a11,a12,a13
     +      ,a21,a22,a23
     +      ,a31,a32,a33
     +      ,b11,b12,b13
     +      ,b21,b22,b23
     +      ,b31,b32,b33
     +      ,ew,tau
     +      ,Max1,Max2
     +      ,Gxx,Gxy,Gxz
     +      ,Gyx,Gyy,Gyz
     +      ,Gzx,Gzy,Gzz
     +      ,px,py,pz
     +      ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +      ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +      ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +      )

c------------
       end if
c------------

      cf = 0.5*hs*wq(i)

      uxel = uxel + (dfx*gxx + dfy*gyx + dfz*gzx)*cf 
      uyel = uyel + (dfx*gxy + dfy*gyy + dfz*gzy)*cf 
      uzel = uzel + (dfx*gxz + dfy*gyz + dfz*gzz)*cf 

      End Do

c-----
c Done
c-----

      return
      end
