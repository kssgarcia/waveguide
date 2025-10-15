      subroutine slp_arc_host
     +
     +     (XC,YC,RD,TH1,TH2,ORNT
     +     ,FX1,FY1
     +     ,FX2,FY2
     +     ,X0,Y0,TH0
     +     ,FX0,FY0
     +     ,Ising
     +     ,Iflow
     +     ,Uarc,Varc
     +     ,Istop
     +     )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c---------------------------------------
c  Compute the single-layer potential
c  over a circular arc
c
c  INT[ F(i)*G(i,j) ]
c
c  If ISING = 1 THE SINGULARITY IS SUBTRACTED OFF
c             0 perform regular Gauss-Legendre
c----------------------------------------------

      Implicit Double Precision (A-H,O-Z)

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

c----------
c constants
c----------

      Null = 0
      Zero = 0.0D0

c---
c initialize
c---

      Istop = 0
      Iopt  = 1    ! velocity Green's function only needed

      Uarc = 0.0D0
      Varc = 0.0D0

c---------------------------
c prepare for the Gauss--Legendre quadrature
c---------------------------
      
      STH = 0.5D0*(TH1+TH2)
      DTH = 0.5D0*(TH2-TH1)

      SFX = 0.5D0*(FX2+FX1)
      DFX = 0.5D0*(FX2-FX1)

      SFY = 0.5D0*(FY2+FY1)
      DFY = 0.5D0*(FY2-FY1)

c----

      Do I=1,NGL

        AN = STH + DTH*ZZ(I)
        CS = Dcos(AN)
        SN = Dsin(AN)
        XX = XC + RD * CS
        YY = YC + RD * SN

c----------------------------------------
        if(Iflow.eq.1.or.Iflow.eq.2) then
c----------------------------------------

        call sgf_2d_fs
     +
     +    (Iopt
     +    ,XX,YY
     +    ,X0,Y0
     +    ,Sxx,Sxy
     +    ,Syx,Syy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

c-------------------------------
        else if(Iflow.eq.5) then
c-------------------------------

        if(Iglut.eq.0) then

        call sgf_2d_2p 
     +
     +    (Iopt
     +    ,XX,YY
     +    ,X0,Y0
     +    ,a11,a12,a21,a22
     +    ,b11,b12,b21,b22
     +    ,ew,tau
     +    ,Max1,Max2
     +    ,Sxx,Sxy
     +    ,Syx,Syy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

         else if(Iglut.eq.1) then

         call sgf_2d_2p_int
     +
     +     (XX,YY
     +     ,X0,Y0
     +     ,a21
     +     ,Sxx,Sxy
     +     ,Syx,Syy
     +     )

         else

          write (6,*) "slp_arc_host: invalid option"
          stop

         end if

c---
        else if(Iflow.eq.10) then
c---

        call sgf_2d_w
     +
     +    (Iopt
     +    ,XX,YY
     +    ,X0,Y0
     +    ,wall
     +    ,Sxx,Sxy
     +    ,Syx,Syy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

c---
      else if(Iflow.eq.11) then
c---

       call sgf_2d_1p_w
     +
     +    (Iopt
     +    ,XX,YY
     +    ,X0,Y0
     +    ,wall,RL
     +    ,Null
     +    ,Sxx,Sxy
     +    ,Syx,Syy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

c---
      else if(Iflow.eq.20) then
c---
        call sgf_2d_1p_ww
     +
     +    (Iopt
     +    ,IQPD
     +    ,XX,YY
     +    ,X0,Y0
     +    ,RL,NGFww,H
     +    ,Sxx,Sxy
     +    ,Syx,Syy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

c---
       end if
c---

c-----------
c compute Df
c-----------

      FX = SFX + DFX*ZZ(i)
      FY = SFY + DFY*ZZ(i)

c---------------------
c compute the velocity
c---------------------

      U = FX*Sxx + FY*Syx
      V = FX*Sxy + FY*Syy

c---
c singularity behaves like: - Df * ln(l)
c subtract it off
c---

      if(Ising.eq.1) then
        DES = Dlog(Dabs(AN-TH0))
        U   = U + DES*FX0
        V   = V + DES*FY0
      end if
c---

      Uarc = Uarc + U*WW(I) 
      Varc = Varc + V*WW(I)

      End Do

c----

      cf = DTH*ORNT*RD

      Uarc = Uarc * cf
      Varc = Varc * cf

c-------------------------
c add back the singularity
c-------------------------

      if(Ising.eq.1)then

        ARG2 = TH2-TH0
        ARG1 = TH1-TH0

        RES = 0.0D0

        if(Dabs(ARG2).GT.0.00000001) 
     +  RES = RES + ARG2*(Dlog(Dabs(ARG2))-1.0D0)

        if(Dabs(ARG1).GT.0.00000001) 
     +  RES  = RES - ARG1*(Dlog(Dabs(ARG1))-1.0D0)

        RES  = RD*ORNT*RES

        Uarc = Uarc - RES*FX0
        Varc = Varc - RES*FY0

      end if

c-----
c done
c-----

 100  Format (5(1X,F10.5))
 110  Format (' CURVATURE =', F13.7)

      return
      end
