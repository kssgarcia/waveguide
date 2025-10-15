	subroutine sgf_2d_1p_ww
     +
     +             (Iselect
     +             ,IQPD
     +             ,X,Y
     +             ,X0,Y0
     +             ,XP
     +             ,NP
     +             ,h
     +             ,Gxx,Gxy
     +             ,Gyx,Gyy
     +             ,px,py
     +             ,Txxx,Txxy,Tyxx,Tyxy
     +             ,Txyx,Txyy,Tyyx,Tyyy
     +             )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------
c  Green's functions of two-dimensional Stokes flow
c  representing the flow due to an
c  infinite x-periodic array of point forces
c  separated by the distance XP  (X-period)
c
c  The flow occurs within a channel
c  confined between two plane
c  parallel walls, a distance of 2h apart.
c
c  The walls are located at y = +-h
c
c  SYMBOLS:
c  --------
c
c  X0,Y0	:  location of the point force
c  X,Y		:  field point
c  XP		:  Period
c  NP		:  Limit of summation in real space
c                  Set at about 5
c  epsilon      :  used for the numerical computation
c		   of the zero-wave number contribution
c
c  IQPD ne 1    :  Pressure Drop (PD) = 0, Flow rate Q is finite
c  IQPD =  1    :  Pressure Drop (PD) is finite, and Q = 0
c
c  Iselect = 1: computes only G
c  Iselect = 2: computes      G, p, and T
c
c-------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Parameter (MMAX=1000000,eps=0.00000001)
      Parameter (epsilon = 0.00001)

c---
c constants
c---

      pi   = 3.14159 265358
      piq  = 0.25*pi
      pih  = 0.50*pi
      pi2  = 2.00*pi
      pi4  = 4.00*pi
      pi6  = 6.00*pi
      pi8  = 8.00*pi

      srpi = sqrt(pi)

c-----------------------
c Prepare and Initialize
c-----------------------

      H4 = 4.0*H

      GPP11 = 0.0
      GPP12 = 0.0
      GPP21 = 0.0
      GPP22 = 0.0

      If(Iselect.ne.1) then
        PP1   = 0.0
        PP2   = 0.0
        DG111 = 0.0
        DG112 = 0.0
        DG121 = 0.0
        DG221 = 0.0
      End If

c-------------------------------
c  first compute the primary part
c  with the image singularities
c  running from -NP to NP
c-------------------------------

      H1 = H+Y0
      H2 = H-Y0
      YM = Y0
      YS = Y-YM
      YY = Y+YM+2.0*H

      Do 10 M=-NP,NP

	  XM = X0+M*XP
	  XS = X-XM

	  If(sqrt(xs**2+ys**2).le.eps) Go to 10

	  SC = pih/H

	  XSPI = XS*SC
	  YSPI = YS*SC
	  YYPI = YY*SC

	  COYY = COS (YYPI)
	  COYS = COS (YSPI)
	  SIYY = SIN (YYPI)
	  SIYS = SIN (YSPI)
	  SHXS = SINH(XSPI)
	  CHXS = COSH(XSPI)

          xsph = 0.5*XSPI
          dxys = CHXS-COYS
          dxyc = CHXS-COYY

	  DG1 = 0.5*LOG(dxyc/dxys)
	  DG2 = SHXS*xsph*( 1.0/dxys- 1.0/dxyc)
	  DG3 =      xsph*(SIYS/dxys-SIYY/dxyc)

	  GPP11 = GPP11 + DG1+DG2
	  GPP22 = GPP22 + DG1-DG2
	  GPP12 = GPP12 + DG3

c----
      If(Iselect.ne.1) then

          tmp12 = dxys**2
          tmp22 = dxyc**2
          tmp3  = xsph*SC
	  DG4   = tmp3*((1.0-CHXS*COYS)/tmp12-
     +                  (1.0-CHXS*COYY)/tmp22 )
	  DG5   = SHXS*tmp3*(SIYY/tmp22-SIYS/tmp12)
          tmp2  = SHXS*SC*( 1.0/dxys- 1.0/dxyc)
          tmp3  =  0.5*SC*(SIYS/dxys-SIYY/dxyc)
	  PP1   = PP1 +   tmp2
	  PP2   = PP2 + 2*tmp3
	  DG111 = DG111 + DG4
	  DG112 = DG112 + DG5 - tmp3
	  DG121 = DG121 + DG5 + tmp3
	  DG221 = DG221 - DG4 - tmp2

      End If
c----

10    Continue

      GPP21 = GPP12

c----
      If(Iselect.ne.1) then

        DG211  =   DG121
        DG222  = - DG121
        DG122  = - DG111
        DG212  =   DG122

        TPP111 = - PP1 + DG111 + DG111
        TPP112 =         DG211 + DG112
        TPP212 = - PP1 + DG212 + DG212

        TPP121 = - PP2 + DG121 + DG121
        TPP122 =         DG122 + DG221
        TPP222 = - PP2 + DG222 + DG222

      End If
c----

c-------------------------------
c Complementary part GCP and TCP
c-------------------------------

c------------
c Preparation
c------------

      CXP = pi2/XP

      BB11 =  0.5*H
      BB13 =  1.0
      BB22 =  BB11
      BB24 =  1.0
      BB33 =  1.0
      BB44 = -1.0

      A_1 = 0.0
      A_2 = 0.0
      B_3 = 0.0
      B_4 = 0.0

      GCP11 = 0.0
      GCP12 = 0.0
      GCP21 = 0.0
      GCP22 = 0.0

c------
      If(Iselect.ne.1) then

        PC1   = 0.0
        PC2   = 0.0
        DG111 = 0.0
        DG112 = 0.0
        DG211 = 0.0
        DG212 = 0.0
        DG121 = 0.0
        DG122 = 0.0
        DG221 = 0.0
        DG222 = 0.0

      End If
c-------

c---------------------------------------
c Will do the zero wave number 
c only when flow pressure drop is finite
c
c This term will then be multipled by 0.5
c----------------------------------------

      If(IQPD.eq.1) then

	  FM   = epsilon*CXP
          FMH  = FM*H
          tmp  = -2.0*FMH
	  thet = EXP(tmp)

	  BB12 = -0.5*H*thet
	  BB14 = thet
	  BB21 = BB12
	  BB23 = thet
	  BB31 = 0.5*     (H-1.0/FM)
	  BB32 = 0.5*thet*(H+1.0/FM)
	  BB34 = - thet
	  BB41 = - BB32
	  BB42 = - BB31
	  BB43 =  thet

          EFMH = EXP(FMH)

          call FGA (H,FM,H1,out1)
          call FGA (H,FM,H2,out2)

          A_3 = - out2/EFMH
          A_4 =   out1/EFMH
	  B_1 = - A_3
	  B_2 = - A_4

c-------------------------------
c SOLUTION BY GAUSS ELIMINATION
c-------------------------------

      fc   = BB21/BB11
      BA22 = BB22-fc*BB12
      BA23 = BB23-fc*BB13
      BA24 = BB24-fc*BB14

      A1_2 = A_2-fc*A_1
      B1_2 = B_2-fc*B_1

      fc   = BB31/BB11
      BA32 = BB32-fc*BB12
      BA33 = BB33-fc*BB13
      BA34 = BB34-fc*BB14

      A1_3 = A_3-fc*A_1
      B1_3 = B_3-fc*B_1

      fc    = BB41/BB11
      BA42  = BB42-fc*BB12
      BA43  = BB43-fc*BB13
      BA44  = BB44-fc*BB14

      A1_4 = A_4-fc*A_1
      B1_4 = B_4-fc*B_1

      fc    = BA32/BA22
      BC33  = BA33-fc*BA23
      BC34  = BA34-fc*BA24

      A2_3 = A1_3-fc*A1_2
      B2_3 = B1_3-fc*B1_2

      fc   = BA42/BA22
      BC43 = BA43-fc*BA23
      BC44 = BA44-fc*BA24
      A2_4 = A1_4-fc*A1_2
      B2_4 = B1_4-fc*B1_2

      fc   = BC43/BC33
      BD44 = BC44-fc*BC34
      A3_4 = A2_4-fc*A2_3
      B3_4 = B2_4-fc*B2_3

      AP4 =  A3_4                            /BD44
      AP3 = (A2_3-BC34*AP4)                  /BC33
      AP2 = (A1_2-BA24*AP4-BA23*AP3)         /BA22
      AP1 = (A _1-BB14*AP4-BB13*AP3-BB12*AP2)/BB11

c----
     
      BP4 =  B3_4                            /BD44
      BP3 = (B2_3-BC34*BP4)                  /BC33
      BP2 = (B1_2-BA24*BP4-BA23*BP3)         /BA22
      BP1 = (B _1-BB14*BP4-BB13*BP3-BB12*BP2)/BB11

      CE   = EXP(FM*Y)
      tmp  = FM*(X-X0)
      SINX = SIN(tmp)
      COSX = COS(tmp)

      tmp1  = CE*AP1-AP2/CE
      tmp2  = CE*BP1-BP2/CE
      tmp3  = CE*AP1+AP2/CE
      tmp4  = CE*BP1+BP2/CE

      tmp11 = CE*AP3-AP4/CE
      tmp21 = CE*BP3-BP4/CE
      tmp31 = CE*AP3+AP4/CE
      tmp41 = CE*BP3+BP4/CE

      tmp91 = Y*tmp1+2.0*tmp31
      tmp92 = Y*tmp2+2.0*tmp41
      tmp93 = Y*tmp3+2.0*tmp11
      tmp94 = Y*tmp4+2.0*tmp21

      DG1 =  tmp91*COSX
      DG2 =  tmp92*SINX
      DG3 = (tmp93+(AP2/CE-CE*AP1)/FM)*SINX
      DG4 = (tmp94+(BP2/CE-CE*BP1)/FM)*COSX

      GCP11 =  0.5*DG1
      GCP12 = -0.5*DG2
      GCP21 =  0.5*DG3
      GCP22 =  0.5*DG4

c-----------
      If(Iselect.ne.1) then

	  DG5  =  2.0 * tmp3 * SINX
	  DG6  =  2.0 * tmp4 * COSX
          COSXFM = COSX*FM
          SINXFM = SINX*FM
	  DG7  =  tmp91         *SINXFM
	  DG8  = (tmp93+tmp1/FM)*COSXFM
	  DG9  = (tmp93-tmp1/FM)*COSXFM
	  DG10 =  tmp92         *COSXFM
	  DG11 = (tmp94+tmp2/FM)*SINXFM
	  DG12 = (tmp94-tmp2/FM)*SINXFM

	  PC1   =   0.5*DG5
	  PC2   =   0.5*DG6
	  DG111 = - 0.5*DG7
	  DG112 =   0.5*DG8
	  DG211 =   0.5*DG9
	  DG121 = - 0.5*DG10
	  DG122 = - 0.5*DG11
	  DG221 = - 0.5*DG12

        End If
c-----------

	End If

c--------------------------------------
c Will do the wave numbers from 1 to NP
c--------------------------------------

	Do 20 M=1,NP

	  FM    = M*CXP
	  FM2   = FM*FM
          FMH   = FM*H
	  En2HW = EXP(-2*FMH)
	  En4HW = EXP(-4*FMH)
	  En6HW = EXP(-6*FMH)
	  En8HW = EXP(-8*FMH)

	  det  = -1/(4*FM2) + (1/(2*FM2) + 4*H*H)*En4HW 
     +          - En8HW/(4*FM2)

	  hn2  = H*En2HW
	  hn4  = H*En4HW
	  h2n2 = H*H*En2HW
	  h2n4 = H*H*En4HW
	  wn   = 1/FM
	  wn2  = En2HW/FM
	  wn4  = En4HW/FM
	  wn6  = En6HW/FM

	  A11 =  2*hn4 - wn/2  + wn4/2
	  A12 = -2*hn2 - wn6/2 + wn2/2 
          A13 =  2*hn4 + wn/2  - wn4/2 
          A14 = -2*hn2 + wn6/2 - wn2/2
	  A21 =  A12
	  A22 =  A11
	  A23 = -A14
	  A24 = -A13
          A31 =  h2n4 - wn*wn/4   + wn*wn4/4 + H*wn/4 + 3*H*wn4/4
	  A32 =  h2n2 - wn*wn6/4  + wn*wn2/4 - H*wn6/4 - 3*H*wn2/4 
          A33 =  h2n4 - H*wn/4  + H*wn4/4
          A34 =  h2n2 + H*wn6/4 - H*wn2/4 
	  A41 =  A32
          A42 =  A31
          A43 = -A34
          A44 = -A33

          call FGA (H,FM,H1,out1)
          call FGA (H,FM,H2,out2)

	  C3 = -out2
	  C4 =  out1
	  D1 =  out2
	  D2 = -out1

	  AP1 = (A13*C3 + A14*C4)/det
	  AP2 = (A23*C3 + A24*C4)/det
	  AP3 = (A33*C3 + A34*C4)/det
	  AP4 = (A43*C3 + A44*C4)/det

	  BP1 = (A11*D1 + A12*D2)/det
	  BP2 = (A21*D1 + A22*D2)/det
	  BP3 = (A31*D1 + A32*D2)/det
	  BP4 = (A41*D1 + A42*D2)/det

c---

      CE   = EXP(FM*(Y-H))
      CEn  = EXP(FM*(-Y-H))
      tmp  = FM*(X-X0)
      SINX = SIN(tmp)
      COSX = COS(tmp)

      tmp1  = CE*AP1-AP2*CEn
      tmp2  = CE*BP1-BP2*CEn
      tmp3  = CE*AP1+AP2*CEn
      tmp4  = CE*BP1+BP2*CEn

      tmp11 = CE*AP3-AP4*CEn
      tmp21 = CE*BP3-BP4*CEn
      tmp31 = CE*AP3+AP4*CEn
      tmp41 = CE*BP3+BP4*CEn

      tmp91 = Y*tmp1 + 2.0*tmp31
      tmp92 = Y*tmp2 + 2.0*tmp41
      tmp93 = Y*tmp3 + 2.0*tmp11
      tmp94 = Y*tmp4 + 2.0*tmp21

      DG1 =  tmp91*COSX
      DG2 =  tmp92*SINX
      DG3 = (tmp93+(AP2*CEn-CE*AP1)/FM)*SINX
      DG4 = (tmp94+(BP2*CEn-CE*BP1)/FM)*COSX

      GCP11 = GCP11 + DG1
      GCP12 = GCP12 - DG2
      GCP21 = GCP21 + DG3
      GCP22 = GCP22 + DG4

c-----

      If(Iselect.ne.1) then

	  DG5    =  2.0* tmp3 * SINX
	  DG6    =  2.0* tmp4 * COSX
          COSXFM = COSX*FM
          SINXFM = SINX*FM
	  DG7    =  tmp91         *SINXFM
	  DG8    = (tmp93+tmp1/FM)*COSXFM
	  DG9    = (tmp93-tmp1/FM)*COSXFM
	  DG10   =  tmp92         *COSXFM
	  DG11   = (tmp94+tmp2/FM)*SINXFM
	  DG12   = (tmp94-tmp2/FM)*SINXFM

	  PC1   = PC1   + DG5
	  PC2   = PC2   + DG6
	  DG111 = DG111 - DG7
	  DG112 = DG112 + DG8
	  DG211 = DG211 + DG9
	  DG121 = DG121 - DG10
	  DG122 = DG122 - DG11
          DG221 = DG221 - DG12

      End If

c----------

   20 Continue

c--------------------------------------------------
c Will do the wave numbers from NP+1 to NMAX
c will stop when adequate accuracy has been reached
c--------------------------------------------------

	DO 30 M=NP+1, MMAX

	  DGC  = 0.

	  FM   = M*CXP
	  FM2  = FM*FM
          FMH  = FM*H
	  En2HW  = EXP(-2*FMH)
	  En4HW  = EXP(-4*FMH)
	  En6HW  = EXP(-6*FMH)
	  En8HW  = EXP(-8*FMH)

	  det  = -1/(4*FM2) + (1/(2*FM2) + 4*H*H)*En4HW 
     +           - En8HW/(4*FM2)
	  hn2  = H*En2HW
	  hn4  = H*En4HW
	  h2n2 = H*H*En2HW
	  h2n4 = H*H*En4HW
	  wn   = 1/FM
	  wn2  = En2HW/FM
	  wn4  = En4HW/FM
	  wn6  = En6HW/FM

	  A11 =  2*hn4 - wn/2  + wn4/2
	  A12 = -2*hn2 - wn6/2 + wn2/2 
          A13 =  2*hn4 + wn/2  - wn4/2 
          A14 = -2*hn2 + wn6/2 - wn2/2
	  A21 =  A12
	  A22 =  A11
	  A23 = -A14
	  A24 = -A13
          A31 =  h2n4 - wn*wn/4   + wn*wn4/4 + H*wn/4 + 3*H*wn4/4
	  A32 =  h2n2 - wn*wn6/4  + wn*wn2/4 - H*wn6/4 - 3*H*wn2/4 
          A33 =  h2n4 - H*wn/4  + H*wn4/4
          A34 =  h2n2 + H*wn6/4 - H*wn2/4 
	  A41 =  A32
          A42 =  A31
          A43 = -A34
          A44 = -A33

          call FGA (H,FM,H1,out1)
          call FGA (H,FM,H2,out2)

	  C3 = -out2
	  C4 =  out1
	  D1 =  out2
	  D2 = -out1

	  AP1 = (A13*C3 + A14*C4)/det
	  AP2 = (A23*C3 + A24*C4)/det
	  AP3 = (A33*C3 + A34*C4)/det
	  AP4 = (A43*C3 + A44*C4)/det

	  BP1 = (A11*D1 + A12*D2)/det
	  BP2 = (A21*D1 + A22*D2)/det
	  BP3 = (A31*D1 + A32*D2)/det
	  BP4 = (A41*D1 + A42*D2)/det

      CE   = EXP(FM*(Y-H))
      CEn  = EXP(FM*(-Y-H))
      tmp  = FM*(X-X0)
      SINX = SIN(tmp)
      COSX = COS(tmp)

      tmp1  = CE*AP1-AP2*CEn
      tmp2  = CE*BP1-BP2*CEn
      tmp3  = CE*AP1+AP2*CEn
      tmp4  = CE*BP1+BP2*CEn

      tmp11 = CE*AP3-AP4*CEn
      tmp21 = CE*BP3-BP4*CEn
      tmp31 = CE*AP3+AP4*CEn
      tmp41 = CE*BP3+BP4*CEn

      tmp91 = Y*tmp1 + 2.0*tmp31
      tmp92 = Y*tmp2 + 2.0*tmp41
      tmp93 = Y*tmp3 + 2.0*tmp11
      tmp94 = Y*tmp4 + 2.0*tmp21

      DG1 =  tmp91*COSX
      DG2 =  tmp92*SINX
      DG3 = (tmp93+(AP2*CEn-CE*AP1)/FM)*SINX
      DG4 = (tmp94+(BP2*CEn-CE*BP1)/FM)*COSX

      GCP11 = GCP11 + DG1
      GCP12 = GCP12 - DG2
      GCP21 = GCP21 + DG3
      GCP22 = GCP22 + DG4

      IF(DGC.LT.ABS(DG1)) DGC = ABS(DG1)
      IF(DGC.LT.ABS(DG2)) DGC = ABS(DG2)
      IF(DGC.LT.ABS(DG3)) DGC = ABS(DG3)
      IF(DGC.LT.ABS(DG4)) DGC = ABS(DG4)
 
c-----

      If(Iselect.ne.1) then

	  DG5    =  2.0* tmp3 * SINX
	  DG6    =  2.0* tmp4 * COSX
          COSXFM = COSX*FM
          SINXFM = SINX*FM
	  DG7    =  tmp91         *SINXFM
	  DG8    = (tmp93+tmp1/FM)*COSXFM
	  DG9    = (tmp93-tmp1/FM)*COSXFM
	  DG10   =  tmp92         *COSXFM
	  DG11   = (tmp94+tmp2/FM)*SINXFM
	  DG12   = (tmp94-tmp2/FM)*SINXFM

	  PC1   = PC1   + DG5
	  PC2   = PC2   + DG6
	  DG111 = DG111 - DG7
	  DG112 = DG112 + DG8
	  DG211 = DG211 + DG9
	  DG121 = DG121 - DG10
	  DG122 = DG122 - DG11
          DG221 = DG221 - DG12

          IF(DGC.LT.ABS(DG5))  DGC=ABS(DG5)
          IF(DGC.LT.ABS(DG6))  DGC=ABS(DG6)
          IF(DGC.LT.ABS(DG7))  DGC=ABS(DG7)
          IF(DGC.LT.ABS(DG8))  DGC=ABS(DG8)
          IF(DGC.LT.ABS(DG9))  DGC=ABS(DG9)
          IF(DGC.LT.ABS(DG10)) DGC=ABS(DG10)
          IF(DGC.LT.ABS(DG11)) DGC=ABS(DG11)
          IF(DGC.LT.ABS(DG12)) DGC=ABS(DG12)
 

      End If

c---------

        If(DGC.LE.eps) Go to 40

  30  Continue

c------
c  Done
c------

  40  Continue

        Gxx = GPP11 + CXP * GCP11
        Gxy = GPP12 + CXP * GCP12
        Gyx = GPP21 + CXP * GCP21
        Gyy = GPP22 + CXP * GCP22

c-----------

      If(Iselect.ne.1) then

      DG212 = - DG111
      DG222 = - DG121

      TCP111 = - PC1 + DG111 + DG111
      TCP112 =         DG112 + DG211
      TCP212 = - PC1 + DG212 + DG212

      TCP121 = - PC2 + DG121 + DG121
      TCP122 =         DG122 + DG221
      TCP222 = - PC2 + DG222 + DG222

      px   = PP1 + CXP * PC1
      py   = PP2 + CXP * PC2

      Txxx = TPP111 
     +               + CXP * TCP111
      Txxy = TPP112 
     +               + CXP * TCP112
      Tyxx = Txxy
      Tyxy = TPP212 
     +               + CXP * TCP212
      Txyx = TPP121 
     +               + CXP * TCP121
      Txyy = TPP122 
     +               + CXP * TCP122
      Tyyx = Txyy
      Tyyy = TPP222 
     +               + CXP * TCP222

      End If

c------

   99 Continue

      Return
      End

c=======================================

      subroutine FGA (H,om,sh,out)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      tmp = 4*H-sh
      out = (sh*(exp(-sh*om) - exp((sh-8*H)*om))
     +       + (sh-4*H)*(exp((sh-4*H)*om) - exp(-(sh+4*H)*om)))
     +      /(1 - exp(-4*H*om))**2

c-----
c Done
c-----

      Return
      End
