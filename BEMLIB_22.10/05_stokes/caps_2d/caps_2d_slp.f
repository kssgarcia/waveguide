      subroutine caps_2d_slp
     +
     +  (IS_slp
     +  ,Isym
     +  ,Imax
     +  ,Iflow
     +  ,Uslp,Vslp
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

c----------------------------------------------
c  Compute the single-layer potential over
c  an interface
c
c  If IS_slp = 1, singularity is subtracted out
c
c  Potential computed at points i = 1, ...., Imax
c
c  If Isym = 1, rest of it is set by reflection
c
c  SYMBOLS:
c  -------
c  suns:  unstressed arc length
c  s:     current    arc length
c  srtn:  surface tension
c  elten: elastic tension
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension X(0:900),Y(0:900)

      Dimension  XC(900), YC(900),  R(900),  S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

      Dimension srtn(0:900)
      Dimension suns(0:900),elten(0:900)

      Dimension FX(0:900),FY(0:900)

      Dimension Uslp(900),Vslp(900)

c--------------
c common blocks
c--------------

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT

      common/TENSION/srtn
      common/ELAST/elst,suns

      common/DDFF/FX,FY

      common/ancR1/thet0,gx,gy
      common/ancR4/Drho

      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ
      common/ancI7/Ielst

      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c-------------------------------------
c compute Df at the nodes
c
c Constitutive equation incorporates:
c
c   (a) Body force
c   (b) Isotropic tension
c   (c) Elastic tension
c
c-------------------------------------

c-------------------------------------
c  Elastic tension at nodes.
c
c  Compute the derivative d(s)/d(suns)
c  by parabolic interpolation
c---------------------------------------

      Do i=0,NSG2
        elten(i) = 0.0D0
      End Do

      if(Ielst.eq.1) then

        Do i=2,NSG1
          ia = i-1
          i1 = i+1
          x0 = suns(ia)-suns(i)
          x1 = suns(i1)-suns(i)
          y0 =    s(ia)-   s(i)
          y1 =    s(i1)-   s(i)
          DsDs0    = (x0*y1/x1 - x1*y0/x0)/(x0-x1)
          elten(i) = elst*(DsDs0-1.0D0)      ! linear constit equation
        End Do

        elten(0)    = elten(NSG)  ! wrap
        elten(1)    = elten(NSG1) ! wrap
        elten(NSG2) = elten(2)

      end if

c-------------------------------------
c  Compute the derivative d(srtn+elten)/d(s)
c  at the nodes by parabolic interpolation
c
c  Compute Df
c------------------------------------

      Do i=2,NSG1

       ia = i-1
       i1 = i+1
       x0 =     s(ia)-    s(i)
       x1 =     s(i1)-    s(i)
       y0 =  srtn(ia)- srtn(i)
     +    + elten(ia)-elten(i)
       y1 =  srtn(i1)- srtn(i)
     +    + elten(i1)-elten(i)
       DtDs = (x0*y1/x1 - x1*y0/x0)/(x0-x1)

       crv  = ORNT(i)/R(i)
       vnx  = ORNT(i) * Dcos(TH2(i))
       vny  = ORNT(i) * Dsin(TH2(i))
       tnx  =-vny
       tny  = vnx

       tentot = srtn(i)+elten(i)

       FN = tentot*crv + Drho*(X(i)*gx+Y(i)*gy)

       FX(i) = FN * vnx - tnx * DtDs
       FY(i) = FN * vny - tny * DtDs

      End Do

      FX(0) = FX(NSG)
      FY(0) = FY(NSG)

      FX(1) = FX(NSG1)
      FY(1) = FY(NSG1)

c-------------
c Done with Df
c-------------

c-------------------------------------------
c  compute the slp at the nodes 1, ..., Imax
c-------------------------------------------

      Do j=1,Imax

      ja = j-1
      j1 = j+1

      X0  =   X(j)
      Y0  =   Y(j)
      TH0 = TH2(j)
      FX0 =  FX(j)
      FY0 =  FY(j)

      Utot = 0.0D0
      Vtot = 0.0D0

c--------------------------------------
c non-singular integrals
c will be computed over the blended arc
c--------------------------------------

      Do 2 i=1,NSG

      i1 = i+1

      if(j.eq.i.or.j.eq.i1)   Go to 2   ! arc is singular
      if(j.eq.1.and.i.eq.NSG) Go to 2   ! arc is singular

      XC1 =   XC(i)
      YC1 =   YC(i)
      RD1 =    R(i)
      T11 =  TH2(i)
      T21 =  TH3(i)
      OR1 = ORNT(i)
      FX1 =   FX(i)
      FY1 =   FY(i)

      XC2 =   XC(i1)
      YC2 =   YC(i1)
      RD2 =    R(i1)
      T12 =  TH1(i1)
      T22 =  TH2(i1)
      OR2 = ORNT(i1)
      FX2 =   FX(i1)
      FY2 =   FY(i1)

      call slp_arc_blended
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

      Utot = Utot + Uarc
      Vtot = Vtot + Varc

  2   Continue

c--------------------------------------------------
c integrate over the two sides of the host arc
c numbered j
c-------------------------------------------------

      Ising = 1

      if(IS_slp.eq.0) Ising = 0

      XCS =   XC(j)
      YCS =   YC(j)
      RDS =    R(j)
      T1S =  TH1(j)
      T2S =  TH2(j)
      T3S =  TH3(j)
      ORS = ORNT(j)

      FX1 = FX(ja)
      FX2 = FX(j)
      FX3 = FX(j1)

      FY1 = FY(ja)
      FY2 = FY(j)
      FY3 = FY(j1)

c---

      call slp_arc_host
     +
     +  (XCS,YCS,RDS,T1S,T2S,ORS
     +  ,FX1,FY1
     +  ,FX2,FY2
     +  ,X0,Y0,TH0
     +  ,FX0,FY0
     +  ,Ising
     +  ,Iflow
     +  ,Uarc,Varc
     +  ,Istop
     +  )

      Utot = Utot + Uarc
      Vtot = Vtot + Varc

c---

      call slp_arc_host
     +
     +  (XCS,YCS,RDS,T2S,T3S,ORS
     +  ,FX2,FY2
     +  ,FX3,FY3
     +  ,X0,Y0,TH0
     +  ,FX0,FY0
     +  ,Ising
     +  ,Iflow
     +  ,Uarc,Varc
     +  ,Istop
     +  )

      Utot = Utot + Uarc
      Vtot = Vtot + Varc

c---

      Uslp(j) = - Utot/pi4
      Vslp(j) = - Vtot/pi4

      End Do

c------------------------------------------------------

c------------------------------
c Quarterly symmetric interface
c------------------------------

      if(Isym.eq.2) then

        Do i=1,NSGQ
         Uslp(NSGM-i+1) = -Uslp(i)
         Vslp(NSGM-i+1) =  Vslp(i)
         Uslp(NSGM+i-1) = -Uslp(i)
         Vslp(NSGM+i-1) = -Vslp(i)
         Uslp(NSG -i+2) =  Uslp(i)
         Vslp(NSG -i+2) = -Vslp(i)
        End Do

      end if

c-----
c done
c-----

 100  Format (1x,i3,3(1x,f15.10))
 200  Format (' Evaluating the SL potential at point:',I2)
 201  Format (' Arc:',I2)

      return
      end
