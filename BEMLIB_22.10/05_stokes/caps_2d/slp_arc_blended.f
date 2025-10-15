      subroutine slp_arc_blended
     +
     +  (XC1,YC1,RD1,T11,T21,OR1
     +  ,XC2,YC2,RD2,T12,T22,OR2
     +  ,FX1,FY1
     +  ,FX2,FY2
     +  ,X0,Y0
     +  ,Iflow
     +  ,Uarc,Varc
     +  ,Istop
     +  )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c========================================

c-----------------------------------
c Compute the single-layer potential
c over a blended circular arc
c-----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ZZ(20),WW(20)

c---
c common blocks
c---

      common/ancR5/RL

      common/hhhh/h,hh,h2,h3,h4,hs
      common/teihos/wall

      common/ancI2/NGFww,IQPD
      common/ancI3/NGL

      common/aaaa/a11,a12,a21,a22
      common/bbbb/b11,b12,b21,b22
      common/ewew/ew,tau
      common/mmmm/Max1,Max2
      common/sgf2d2p/Iglut

      common/ZZWW/ZZ,WW
      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c---
c constants
c---

      Null = 0
      Zero = 0.0

c---
c prepare
c---

      Istop = 0
      Iopt  = 1     ! only G is needed

      Uarc = 0.0D0
      Varc = 0.0D0

c----------------------------
c prepare for Gauss--Legendre
c----------------------------

      ST1 = 0.5D0*(T21 + T11)
      DT1 = 0.5D0*(T21 - T11)

      ST2 = 0.5D0*(T22 + T12)
      DT2 = 0.5D0*(T22 - T12)

      SFX = 0.5D0*(FX2 + FX1)
      DFX = 0.5D0*(FX2 - FX1)

      SFY = 0.5D0*(FY2 + FY1)
      DFY = 0.5D0*(FY2 - FY1)

c---
      Do i=1,NGL

      AN1 = ST1 + DT1*ZZ(i)
      CS1 = Dcos(AN1)
      SN1 = Dsin(AN1)
      XX1 = XC1 + RD1 * CS1
      YY1 = YC1 + RD1 * SN1

      AN2 = ST2 + DT2*ZZ(i)
      CS2 = Dcos(AN2)
      SN2 = Dsin(AN2)
      XX2 = XC2 + RD2 * CS2
      YY2 = YC2 + RD2 * SN2

      XX = 0.5D0*(XX1+XX2)
      YY = 0.5D0*(YY1+YY2)

c---
        if(Iflow.eq.1.or.Iflow.eq.2) then
c---
        call sgf_2d_fs
     +
     +   (Iopt
     +   ,XX,YY
     +   ,X0,Y0
     +   ,Sxx,Sxy
     +   ,Syx,Syy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

c---
        else if(Iflow.eq.5.or.Iflow.eq.6) then
c---

        if(Iglut.eq.0) then

        call sgf_2d_2p
     +
     +   (Iopt
     +   ,XX,YY
     +   ,X0,Y0
     +   ,a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,tau
     +   ,Max1,Max2
     +   ,Sxx,Sxy
     +   ,Syx,Syy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

         else if(Iglut.eq.1) then

         call sgf_2d_2p_int
     +
     +   (XX,YY
     +   ,X0,Y0
     +   ,a21
     +   ,Sxx,Sxy
     +   ,Syx,Syy
     +   )

         else

          write (6,*) "susp_2d_slp: invalid option"
          stop

         end if

c---
       else if(Iflow.eq.10) then
c---
       call sgf_2d_w
     +
     +   (Iopt
     +   ,XX,YY
     +   ,X0,Y0
     +   ,Wall
     +   ,Sxx,Sxy
     +   ,Syx,Syy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

c---
       else if(Iflow.eq.11) then
c---

       call sgf_2d_1p_w
     +
     +   (Iopt
     +   ,XX,YY
     +   ,X0,Y0
     +   ,Wall,RL,Null
     +   ,Sxx,Sxy
     +   ,Syx,Syy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )
c---

c---
       else if(Iflow.eq.20) then

c---
        call sgf_2d_1p_ww
     +
     +   (Iopt
     +   ,IQPD
     +   ,XX,YY
     +   ,X0,Y0
     +   ,RL,NGFww,H
     +   ,Sxx,Sxy
     +   ,Syx,Syy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )
c-- -
       end if
c---

c-----------
c compute Df
c-----------

      FX = SFX + DFX*ZZ(i)
      FY = SFY + DFY*ZZ(i)

c-----------------
c compute velocity
c-----------------

      U = FX*Sxx + FY*Syx
      V = FX*Sxy + FY*Syy

      Uarc = Uarc + U*WW(i)
      Varc = Varc + V*WW(i)

      End Do

c---------------------------

      cf = 0.5D0*( DT1*OR1*RD1+ DT2*OR2*RD2 ) 

      Uarc = Uarc * cf
      Varc = Varc * cf

c-----
c Done
c-----

 100  Format (5(1X,F10.5))
 110  Format (' CURVATURE =', F13.7)

      return
      end
