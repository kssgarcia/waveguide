      subroutine caps_2d_dlp
     + 
     +  (J
     +  ,XC,YC,RAD
     +  ,TH1,TH2,TH3,ORNT
     +  ,Ising
     +  ,Imax
     +  ,Iflow
     +  ,Tpxx,Tpyx
     +  ,Tpxy,Tpyy
     +  ,Istop
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------------
c Compute the influence matrix for the
c double-layer potential over the 
c ciruclar arc numbered J
c
c The influence matrix is evaluated
c at nodes labeled: 1,...,Imax
c
c Assumes tent-function variation in velocity
c over each arc
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension X(0:900),Y(0:900)

      Dimension Tpxx(900),Tpyx(900)
      Dimension Tpxy(900),Tpyy(900)

      Dimension ZZ(20),WW(20)

c---
c common blocks
c---

      common/XXYY/X,Y

      common/ancR5/RL

      common/hhhh/h,hh,h2,h3,h4,hs
      common/teihos/wall

      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ
      common/ancI2/NGFww,IQPD
      common/ancI3/NGL

c---
c for doubly-periodic flow
c---

      common/aaaa/a11,a12,a21,a22
      common/bbbb/b11,b12,b21,b22
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

c---
c various
c---

      common/ZZWW/ZZ,WW
      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c----------
c constants
c----------

      null = 0
      zero = 0.0D0

c---
c Prepare to run
c---

      Istop = 0
      Iopt  = 2   ! for the Green's function

c---
c     if(Ising.eq.1) then
c       CT2  = COS(TH2)
c       ST2  = SIN(TH2)
c       CT22 = CT2**2
c       ST22 = ST2**2
c       CTST = CT2*ST2
c     end if
c---

      Do i=1,Imax
        tpxx(i) = 0.0D0
        tpyx(i) = 0.0D0
        tpxy(i) = 0.0D0
        tpyy(i) = 0.0D0
      End Do

c----------------------
c dirst part of the arc
c----------------------

      STH = 0.5D0*(TH1+TH2)
      DTH = 0.5D0*(TH2-TH1)
      FAC = DTH*ORNT*RAD

      Do k=1,NGL

        TH  = STH + DTH*ZZ(k)
        CS  = Dcos(TH)
        SN  = Dsin(TH)
        XX  = XC + RAD*CS
        YY  = YC + RAD*SN
        vnx = CS*ORNT
        vny = SN*ORNT

        U  = 0.5D0*(1.0D0+ZZ(k))
        UW = U*WW(k)*FAC

c---
        Do i=1,Imax

        X0 = X(i)
        Y0 = Y(i)

c---
        if(Iflow.eq.1.or.Iflow.eq.2) then 
c---
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

c---
        else if(Iflow.eq.5.or.Iflow.eq.6) then
c---
c       write (6,*) a11,a12,a21,a22
c       write (6,*) b11,b12,b21,b22

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
     +    ,Wall,RL
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

        tpxx(i) = tpxx(i) + UW*(Txxx*vnx+Txxy*vny)
        tpyx(i) = tpyx(i) + UW*(Tyxx*vnx+Tyxy*vny)
        tpxy(i) = tpxy(i) + UW*(Txyx*vnx+Txyy*vny)
        tpyy(i) = tpyy(i) + UW*(Tyyx*vnx+Tyyy*vny)

      End Do

c----------------------------------
c Subtract off free-space Stokeslet
c----------------------------------

c      IF(Ising.eq.1) THEN
c        UW   = WW(K)*FAC
c        VW   = WW(K)*FAC
c        X0   = X(M)
c        Y0   = Y(M)
c        DX   = XX-X0
c        DY   = YY-Y0
c        DR2  = DX**2+DY**2
c        fc   = -4.0/DR2**2
c        Txxx = fc*DX*DX*DX
c        Txxy = fc*DX*DX*DY
c        Tyxy = fc*DY*DX*DY
c        Tyyy = fc*DY*DY*DY
c        Txyx = Txxy
c        Txyy = Tyxy
c        tpxx(M) = tpxx(M) - UW*(Txxx*vnx+Txxy*vny)
c        tpyx(M) = tpyx(M) - VW*(Txxy*vnx+Tyxy*vny)
c        tpxy(M) = tpxy(M) - UW*(Txyx*vnx+Txyy*vny)
c        tpyy(M) = tpyy(M) - VW*(Txyy*vnx+Tyyy*vny)
c        fc      = 2.0*ORNT/RAD
c        Tpxx(J) = Tpxx(J) + fc*UW*ST22
c        Tpyx(J) = Tpyx(J) - fc*VW*CTST
c        Tpxy(J) = Tpxy(J) - fc*UW*CTST
c        Tpyy(J) = Tpyy(J) + fc*VW*CT22
c      ENDIF


       End Do ! loop over Gaussian points

c----

c----------------------
c second part of the arc
c----------------------

      STH = 0.5D0*(TH3+TH2)
      DTH = 0.5D0*(TH3-TH2)
      FAC = DTH*ORNT*RAD

      Do k=1,NGL

      TH  = STH + DTH*ZZ(k)
      CS  = Dcos(TH)
      SN  = Dsin(TH)
      XX  = XC + RAD*CS
      YY  = YC + RAD*SN
      vnx = CS*ORNT
      vny = SN*ORNT

      U  = 0.5D0*(1.0-ZZ(k))
      UW = U*WW(k)*FAC

c---
      Do i=1,Imax

        X0 = X(i)
        Y0 = Y(i)

c---
        if(Iflow.eq.1.or.Iflow.eq.2) then
c---
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

c---
        else if(Iflow.eq.5.or.Iflow.eq.6) then
c---
c       write (6,*) a11,a12,a21,a22

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

c---
      else if(Iflow.eq.10) then
c---
       call sgf_2d_w
     +
     +    (Iopt
     +    ,XX,YY
     +    ,X0,Y0
     +    ,Wall
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
     +    ,Wall,RL
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
      tpxx(i) = tpxx(i) + UW*(Txxx*vnx+Txxy*vny)
      tpyx(i) = tpyx(i) + UW*(Tyxx*vnx+Tyxy*vny)
      tpxy(i) = tpxy(i) + UW*(Txyx*vnx+Txyy*vny)
      tpyy(i) = tpyy(i) + UW*(Tyyx*vnx+Tyyy*vny)

      End Do

c----------------------------------
c Subtract off free-space Stokeslet
c----------------------------------

c      IF(Ising.eq.1) THEN
c          UW   = WW(K)*FAC
c          VW   = WW(K)*FAC
c          X0   = X(M)
c          Y0   = Y(M)
c          DX   = XX-X0
c          DY   = YY-Y0
c          DR2  = DX**2+DY**2
c          fc   = -4.0/DR2**2
c          Txxx = fc*DX*DX*DX
c          Txxy = fc*DX*DX*DY
c          Tyxy = fc*DY*DX*DY
c          Tyyy = fc*DY*DY*DY
c          Txyx = Txxy
c          Txyy = Tyxy
c          tpxx(J) = tpxx(J) - UW*(Txxx*vnx+Txxy*vny)
c          tpyx(J) = tpyx(J) - VW*(Txxy*vnx+Tyxy*vny)
c          tpxy(J) = tpxy(J) - UW*(Txyx*vnx+Txyy*vny)
c          tpyy(J) = tpyy(J) - VW*(Txyy*vnx+Tyyy*vny)
c        fc      = 2.0*ORNT/RAD
c        Tpxx(J) = Tpxx(J) + fc*UW*ST22
c        Tpyx(J) = Tpyx(J) - fc*VW*CTST
c        Tpxy(J) = Tpxy(J) - fc*UW*CTST
c       Tpyy(J) = Tpyy(J) + fc*VW*CT22
c 
c      ENDIF

      End Do          ! End of loop over Gaussian points

c-------------------------
c ADD Subtracted Stokeslet
c from host arc # J
c----------------------------------

c      If(Ising.eq.1) then
c        RES  =  - 2.0*(TH3-TH1)
c        Tpxx(J) = Tpxx(J) + RES*ST22
c        Tpyx(J) = Tpyx(J) - RES*CTST
c        Tpxy(J) = Tpxy(J) - RES*CTST
c        Tpyy(J) = Tpyy(J) + RES*CT22
c      End If
c---

 101  Format (1x,I2,5(1X,F10.5))

      return
      end
