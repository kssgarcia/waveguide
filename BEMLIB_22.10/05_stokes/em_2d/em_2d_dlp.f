      subroutine em_2d_dlp
     + 
     +  (Intrfc
     +  ,Ipoint
     +  ,XC,YC,RAD
     +  ,TH1,TH2,TH3,ORNT
     +  ,Ising
     +  ,Iflow
     +  ,Tpxx,Tpyx
     +  ,Tpxy,Tpyy
     +  ,Istop
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c--------------------------------------------
c  Compute the influence matrix for the
c  double-layer potential over the circular
c  arc numbered: Ipoint
c  located at the drop interface numbered: Intrfc
c
c  Assume linear tent-function variation 
c  in velocity over the arc
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)
 
      Dimension xg(0:066,49),yg(0:066,49)

      Dimension NSG(49),NSG1(49),NSG2(49)

      Dimension tpxx(066,49),tpyx(066,49)
      DImension tpxy(066,49),tpyy(066,49)

      Dimension ZZ(20),WW(20)

c--------------
c common blocks
c--------------

      common/XXYYg/Xg,Yg

      common/ancR5/RL

      common/hhhh/h,hh,h2,h3,h4,hs
      common/teihos/wall

      common/ancI1/Ndrops,NSG,NSG1,NSG2
      common/ancI2/NGFww,IQPD
      common/ancI3/NGL

c---
c doubly periodic flow
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

      Do j=1,Ndrops
       Do i=1,NSG(j)
        tpxx(i,j) = 0.0
        tpyx(i,j) = 0.0
        tpxy(i,j) = 0.0
        tpyy(i,j) = 0.0
       End Do
      End Do

c----------------------
c first part of the arc
c----------------------

      STH = 0.5*(TH1+TH2)
      DTH = 0.5*(TH2-TH1)
      FAC = DTH*ORNT*RAD

      Do k=1,NGL

        TH  = STH + DTH*ZZ(k)
        CS  = Dcos(TH)
        SN  = Dsin(TH)
        XX  = XC + RAD * CS
        YY  = YC + RAD * SN
        vnx = CS*ORNT
        vny = SN*ORNT

        U  = 0.5D0*(1.0+ZZ(k))
        UW = U*WW(k)*FAC

c------
        Do j=1,Ndrops     ! Influence on ith point
        Do i=1,NSG(j)     !  at the jth interface

        X0 = xg(i,j)
        Y0 = yg(i,j)

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

c-------------------------------
        else if(Iflow.eq.10) then
c-------------------------------

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

c-------------------------------
      else if(Iflow.eq.11) then
c-------------------------------

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

c-------------------------------
      else if(Iflow.eq.20) then
c-------------------------------

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

c------------
       end if
c------------

        tpxx(i,j) = tpxx(i,j) + UW*(Txxx*vnx+Txxy*vny)
        tpyx(i,j) = tpyx(i,j) + UW*(Tyxx*vnx+Tyxy*vny)
        tpxy(i,j) = tpxy(i,j) + UW*(Txyx*vnx+Txyy*vny)
        tpyy(i,j) = tpyy(i,j) + UW*(Tyyx*vnx+Tyyy*vny)

      End Do     ! over influence points
      End Do

      End Do! over Gaussian points

c------

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
      XX  = XC + RAD * CS
      YY  = YC + RAD * SN
      vnx = CS*ORNT
      vny = SN*ORNT

      U  = 0.5*(1.0-ZZ(k))
      UW = U*WW(k)*FAC

c---
        Do j=1,Ndrops     ! Influence on ith point
        Do i=1,NSG(j)     !  at the jth interface

        X0 = xg(i,j)
        Y0 = yg(i,j)

c-------------------------------
        if(Iflow.eq.1.or.Iflow.eq.2) then
c-------------------------------

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

c-------------------------------
        else if(Iflow.eq.10) then
c-------------------------------

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

c------------------------------
      else if(Iflow.eq.11) then
c------------------------------

       call sgf_2d_1p_w
     +
     +    (Iopt
     +    ,XX,YY
     +    ,X0,Y0
     +    ,Wall,RL,Null
     +    ,Sxx,Sxy
     +    ,Syx,Syy
     +    ,px,py
     +    ,Txxx,Txxy,Tyxx,Tyxy
     +    ,Txyx,Txyy,Tyyx,Tyyy
     +    )

c----------------------------------------
      else if(Iflow.eq.20) then
c----------------------------------------

        call sgf_2d_1p_ww
     +
     +     (Iopt
     +     ,IQPD
     +     ,XX,YY
     +     ,X0,Y0
     +     ,RL,NGFww,H
     +     ,Sxx,Sxy
     +     ,Syx,Syy
     +     ,px,py
     +     ,Txxx,Txxy,Tyxx,Tyxy
     +     ,Txyx,Txyy,Tyyx,Tyyy
     +     )

c-----------
      end if
c-----------

      tpxx(i,j) = tpxx(i,j) + UW*(Txxx*vnx+Txxy*vny)
      tpyx(i,j) = tpyx(i,j) + UW*(Tyxx*vnx+Tyxy*vny)
      tpxy(i,j) = tpxy(i,j) + UW*(Txyx*vnx+Txyy*vny)
      tpyy(i,j) = tpyy(i,j) + UW*(Tyyx*vnx+Tyyy*vny)

      End Do
      End Do

      End Do ! over Gaussian points

c-----
c done
c-----

 101  Format (1x,I2,5(1X,F10.5))

      return
      end
