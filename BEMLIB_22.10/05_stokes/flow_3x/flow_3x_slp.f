      subroutine flow_3x_slp 
     +
     +  (X0,Y0,T0
     +  ,X1,Y1,T1
     +  ,X2,Y2,T2
     +  ,NGL
     +  ,Ising
     +  ,Itype
     +  ,rad,xcnt,ycnt
     +  ,QRxx1,QRxs1,QRsx1
     +  ,QRss1,QRff1,QIxf1
     +  ,QIfx1,QIsf1,QIfs1
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c----------------------------------------
c Integrates the Green's function over
c a straight segment or a circular arc
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

      common/ZZWW/ZZ,WW

c---
c initialize
c---

      QRxx1 = 0.0D0
      QRxs1 = 0.0D0
      QRsx1 = 0.0D0
      QRss1 = 0.0D0
      QRff1 = 0.0D0

      QIxf1 = 0.0D0
      QIfx1 = 0.0D0
      QIsf1 = 0.0D0
      QIfs1 = 0.0D0

c---
c prepare for the quadrture
c---

      If(Itype.eq.1) then     ! linear element

        XM = 0.5D0*(X2+X1)
        XD = 0.5D0*(X2-X1)
        YM = 0.5D0*(Y2+Y1)
        YD = 0.5D0*(Y2-Y1)
        DR = Dsqrt(XD*XD+YD*YD)

      Else                    ! circular arc

        TM = 0.5D0*(T2+T1)
        TD = 0.5D0*(T2-T1)
        DR = rad*Dabs(TD)

      End If

c--------------------------
c loop over Gaussian points
c--------------------------

      Do 1 i=1,NGL

        If(Itype.eq.1) then        ! linear element
          X = XM + XD*ZZ(i)
          Y = YM + YD*ZZ(i)
        Else                       ! arc element
          T = TM + TD*ZZ(i)
          X = xcnt + rad*Dcos(T)
          Y = ycnt + rad*Dsin(T)
        End If

        call sgf_3x 
     +
     +   (X,Y
     +   ,X0,Y0
     +   ,WRxx1,WRxs1,WRsx1
     +   ,WRss1,WRff1
     +   ,WIxf1
     +   ,WIfx1,WIsf1,WIfs1
     +   )

c---
c subtract off the singularity
c---

      If(Ising.eq.1) then

         If(Itype.eq.1) Dist2 = (X-X0)**2+(Y-Y0)**2
         If(Itype.eq.2) Dist2 = ( rad*(T0-T) )**2
         DD = log(Dist2)
         WRxx1 = WRxx1 +       DD
         WRss1 = WRss1 +       DD
         WRff1 = WRff1 + 2.0D0 * DD

      End If

c---
c carry out the quadrature
c---

      WI = WW(I)

      QRxx1 = QRxx1 + WRxx1 * WI
      QRxs1 = QRxs1 + WRxs1 * WI
      QRsx1 = QRsx1 + WRsx1 * WI
      QRss1 = QRss1 + WRss1 * WI
      QRff1 = QRff1 + WRff1 * WI
      QIxf1 = QIxf1 + WIxf1 * WI
      QIfx1 = QIfx1 + WIfx1 * WI
      QIsf1 = QIsf1 + WIsf1 * WI
      QIfs1 = QIfs1 + WIfs1 * WI

   1  Continue

c---
c finish up
c---

      QRxx1 = QRxx1 * DR
      QRxs1 = QRxs1 * DR
      QRsx1 = QRsx1 * DR
      QRss1 = QRss1 * DR
      QRff1 = QRff1 * DR
      QIxf1 = QIxf1 * DR
      QIfx1 = QIfx1 * DR
      QIsf1 = QIsf1 * DR
      QIfs1 = QIfs1 * DR

c---------------------
c add singularity back
c---------------------

      If(Ising.eq.1) then
        QRxx1 = QRxx1 - 4.0D0* DR*(Dlog(DR)-1.0D0)
        QRss1 = QRss1 - 4.0D0* DR*(Dlog(DR)-1.0D0)
        QRff1 = QRff1 - 8.0D0* DR*(Dlog(DR)-1.0D0)
      End If

c-----
c Done
c-----

 100  Format (3(1x,f15.10))

      Return
      End
