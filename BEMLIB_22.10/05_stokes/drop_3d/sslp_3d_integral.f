      subroutine sslp_3d_integral
     +
     +   (x0,y0,z0
     +   ,k
     +   ,mint
     +   ,Iflow
     +   ,uxel,uyel,uzel
     +   )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c--------------------------------------
c Compute the single-layer potential
c over the non-singular triangle labeled k
c using the Gauss triangle quadrature
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension   Df(1026,3)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/dff/Df
      common/var/shrt,wall
      common/trq/xiq,etq,wq

c---
c for triply-periodic flow
c---

      common/aaaa/a11,a12,a13,a21,a22,a23,a31,a32,a33
      common/bbbb/b11,b12,b13,b21,b22,b23,b31,b32,b33
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

c----------------------
c launch the quadrature
c----------------------

      i1 = n(k,1)   ! global element node labels
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be = beta (k)
      ga = gamma(k)

      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        call sslp_3d_interp 
     +
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +
     +    ,df(i1,1),df(i1,2),df(i1,3)
     +    ,df(i2,1),df(i2,2),df(i2,3)
     +    ,df(i3,1),df(i3,2),df(i3,3)
     +    ,df(i4,1),df(i4,2),df(i4,3)
     +    ,df(i5,1),df(i5,2),df(i5,3)
     +    ,df(i6,1),df(i6,2),df(i6,3)
     +
     +    ,al,be,ga
     +    ,xi,eta
     +
     +    ,x,y,z
     +    ,dfx,dfy,dfz
     +    ,hs
     +    )

c-----------------------------
c compute the Green's function
c-----------------------------

      Iopt_sgf = 1   ! need only G 

c-------------------------
       if(Iflow.eq.1) then
c-------------------------

       call sgf_3d_fs 
     +
     +   (Iopt_sgf
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,gxx,gxy,gxz
     +   ,gyx,gyy,gyz
     +   ,gzx,gzy,gzz
     +   ,px,py,pz
     +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +   )

c------------------------------
       else if(Iflow.eq.2) then
c------------------------------

       call sgf_3d_w 
     +
     +   (Iopt_sgf
     +   ,y,z,x
     +   ,y0,z0,x0
     +   ,wall
     +   ,gyy,gyz,gyx
     +   ,gzy,gzz,gzx
     +   ,gxy,gxz,gxx
     +   ,py,pz,px
     +   ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +   ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +   ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +   )

c-----------------------------
      else if(Iflow.eq.3) then
c-----------------------------

      call sgf_3d_3p
     +
     +  (Iopt_sgf
     +  ,x,y,z
     +  ,x0,y0,z0
     +  ,a11,a12,a13
     +  ,a21,a22,a23
     +  ,a31,a32,a33
     +  ,b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,ew,tau
     +  ,Max1,Max2
     +  ,Gxx,Gxy,Gxz
     +  ,Gyx,Gyy,Gyz
     +  ,Gzx,Gzy,Gzz
     +  ,px,py,pz
     +  ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +  ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +  ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +  )

c------------
       end if
c------------

      cf = 0.5D0*hs*wq(i)

      uxel = uxel + (dfx*gxx + dfy*gyx + dfz*gzx)*cf 
      uyel = uyel + (dfx*gxy + dfy*gyy + dfz*gzy)*cf 
      uzel = uzel + (dfx*gxz + dfy*gyz + dfz*gzz)*cf 

      End Do

c-----
c done
c-----

      return
      end
