      subroutine velocity 
     +
     +  (X00,Y00,phi
     +  ,Ux,Uy,Uz
     +  )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------
c Evaluates the velocity
c at a point in the flow
c-----------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  Xop(129) ,Yop(129)
      Dimension fopx(129),fops(129),fopf(129)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension XW(10,129),YW(10,129),TW(10,129)
      Dimension fx(10,129),fs(10,129),ff(10,129)

      Parameter (optol=0.00001)

c--------------
c common blocks
c--------------

      common/VEL00/Iflow,NEop,NSG,NGL,NE,Itp
      common/VEL01/opening,Xop,Yop,fopx,fops,fopf
      common/VEL02/XW,YW,TW,fx,fs,ff
      common/VEL03/actis,xcntr,ycntr
      common/VEL04/visc,shrt,delta

c----------
c constants
c----------

      pi  = 3.14159 265358
      pi8 = 8.0*pi

      Ising = 0

c---
c initialize
c---

      U00 = 0.0
      V00 = 0.0
      W00 = 0.0

c------------------------------
c Contribution from the opening
c------------------------------

      If(opening.gt.optol) then

      Itype = 1

        Do K=1,NEop

        X1 = Xop(K)
        Y1 = Yop(K)
        X2 = Xop(K+1)
        Y2 = Yop(K+1)

         call flow_3x_slp 
     +
     +      (X00,Y00,unused
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itype
     +      ,Rad,xcnt,ycnt
     +      ,QRxx,QRxs,QRsx
     +      ,QRss,QRff,QIxf
     +      ,QIfx,QIsf,QIfs
     +      )
 
        U00 = U00+QRxx*fopx(K)
     +           +QRxs*fops(K)
     +           -QIxf*fopf(K)

        V00 = V00+QRsx*fopx(K)
     +           +QRss*fops(K)
     +           -QIsf*fopf(K)

        W00 = W00+QIfx*fopx(K)
     +           +QIfs*fops(K)
     +           +QRff*fopf(K)

        End Do

      End If

c-------------------------------
c Contribution from the segments
c-------------------------------

      Do K=1,NSG

       rad  = actis(k)
       xcnt = xcntr(k)
       ycnt = ycntr(k)

        Do L=1,NE(K)

        X1 = XW(K,L)
        Y1 = YW(K,L)
        T1 = TW(K,L)

        X2 = XW(K,L+1)
        Y2 = YW(K,L+1)
        T2 = TW(K,L+1)

         call flow_3x_slp 
     +
     +       (X00,Y00,unused
     +       ,X1,Y1,T1
     +       ,X2,Y2,T2
     +       ,NGL
     +       ,Ising
     +       ,Itp(k)
     +       ,Rad,xcnt,ycnt
     +       ,QRxx,QRxs,QRsx
     +       ,QRss,QRff,QIxf
     +       ,QIfx,QIsf,QIfs
     +       )

        U00 = U00-(QRxx*fx(K,L)
     +                         +QRxs*fs(K,L)-QIxf*ff(K,L))
        V00 = V00-(QRsx*fx(K,L)
     +                         +QRss*fs(K,L)-QIsf*ff(K,L))
        W00 = W00-(QIfx*fx(K,L)
     +                         +QIfs*fs(K,L)+QRsf*ff(K,L))

        End Do

      End Do

c------------------------
c normalize and transform
c------------------------

      U00 = U00/pi8
      V00 = V00/pi8
      W00 = W00/pi8

      cs  = cos(phi)
      sn  = sin(phi)

      Ux =   U00*cs
      Us =   V00*cs
      Uf = - W00*sn

c---------------------
c Cartesian components
c---------------------

      Uy = Us*cs - Uf*sn
      Uz = Us*sn + Uf*cs

c-------------------
c Add reference flow
c-------------------

      If(X00.gt.0.0) Uy = Uy+shrt*X00
     +                  +0.5*delta/visc*X00**2

c-----
c Done
c-----

      Return
      End
