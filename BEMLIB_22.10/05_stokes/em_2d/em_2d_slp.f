      subroutine em_2d_slp
     +
     + (IS_slp
     + ,Iflow
     + ,uslpg,vslpg
     + ,Istop
     + )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------------------
c Compute the single-layer potential over
c the interfaces of Ndrops
c
c Maximum capacity 25 drops
c
c If IS_slp = 1, singularity is subtracted out
c                otherwise use Gauss quadrature
c
c  suns:  unstressed arc length
c  s:     current    arc length
c  srtn:  surface tension
c  elten: elastic tension
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c global variables
c---
 
      Dimension xg(0:066,49),yg(0:066,49)

      Dimension  XCg(066,49), YCg(066,49),  Rg(066,49),   Sg(066,49)
      Dimension TH1g(066,49),TH2g(066,49),TH3g(066,49),orntg(066,49)

      Dimension srtng(0:066,49)
      Dimension sunsg(0:066,49)
      Dimension   fxg(0:066,49),fyg(0:066,49)

      Dimension Uslpg(066,49),Vslpg(066,49)
 
c---
c properties of drops and interfaces
c---
 
      Dimension  NSG(49),NSG1(49),NSG2(49)
      Dimension rhod(49)
      Dimension elst(49),Ielstg(49)
 
c---
c local variables
c---

      Dimension X(0:200),Y(0:200)
      Dimension  XC(200), YC(200),  R(200),   S(200)
      Dimension TH1(200),TH2(200),TH3(200),ORNT(200)

      Dimension srtn(0:200)
      Dimension suns(0:200),elten(0:200)

c--------------
c common blocks
c--------------

c---
c for global variables
c---
 
      common/XXYYg/Xg,Yg
      common/ARCCg/XCg,YCg,Rg,Sg,TH1g,TH2g,TH3g,ORNTg

      common/TENSIONg/srtng
      common/ELASTg/sunsg
      common/DDFFg/fxg,fyg

c---
c various
c---
 
      common/ancR1/gx,gy
      common/ancR4/rho1,rhod,elst
 
      common/ancI1/Ndrops,NSG,NSG1,NSG2
      common/ancI7/Ielst,Ielstg

      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c-------------------------------------
c Compute Df at the nodes of all interfaces
c
c Constitutive equation incorporates:
c
c   (a) Body force
c   (b) Isotropic tension
c   (c) Elastic tension
c-------------------------------------

c-------------------------------------
c     Do j=1,Ndrops
c      Do i=1,NSG2(j)
c       write (6,*) i,j,sg(i,j)
c      End Do
c     End Do
c-------------------------------------

c-------------------------------------

      Do j=1,Ndrops              ! Run over drops
 
c---
c transfer into local vectors
c---
 
       Do i=1,NSG2(j)
           x(i) =    xg(i,j)
           y(i) =    yg(i,j)
        srtn(i) = srtng(i,j)
        suns(i) = sunsg(i,j)
          xc(i) =   xcg(i,j)
          yc(i) =   ycg(i,j)
           r(i) =    rg(i,j)
           s(i) =    sg(i,j)
         th1(i) =  th1g(i,j)
         th2(i) =  th2g(i,j)
         th3(i) =  th3g(i,j)
        ornt(i) = orntg(i,j)
       End Do

c---
c  Elastic tensions
c
c  Compute the derivative d(s)/d(suns)
c  at the nodes by parabolic interpolation
c---

        Do i=0,NSG1(j)
          elten(i) = 0.0D0
        End Do

       if(Ielstg(j).eq.1) then

          Do i=2,NSG1(j)
            ia = i-1
            i1 = i+1
            x0    = suns(ia)-suns(i)
            x1    = suns(i1)-suns(i)
            y0    =    s(ia)-   s(i)
            y1    =    s(i1)-   s(i)
            DsDs0 = (x0*y1/x1 - x1*y0/x0)/(x0-x1)
            elten(i) = elst(j)*(DsDs0-1.0D0)  ! linear constit equation
          End Do

          elten(0)       = elten(NSG(j))  ! wrap
          elten(1)       = elten(NSG1(j))
          elten(NSG2(j)) = elten(2)

       end if

c---------------
c  Compute the derivative d(srtn+elten)/d(s)
c  at the nodes by parabolic interpolation
c
c  Compute Df
c--------------

       Do i=2,NSG1(j)
 
        ia = i-1
        i1 = i+1
        x0   =    s(ia)-   s(i)
        x1   =    s(i1)-   s(i)
        y0   = srtn(ia)-srtn(i)+elten(ia)-elten(i)
        y1   = srtn(i1)-srtn(i)+elten(i1)-elten(i)
        DtDs = (x0*y1/x1 - x1*y0/x0)/(x0-x1)

        crv = ornt(i)/r(i)
        vnx = ornt(i) * Dcos(th2(i))
        vny = ornt(i) * Dsin(th2(i))
        tnx =-vny
        tny = vnx

        tentot = srtn(i)+elten(i)

        FN = tentot*crv +(rho1-rhod(j))*(x(i)*gx+y(i)*gy) 

        fxg(i,j) = FN*vnx - tnx*DtDs
        fyg(i,j) = FN*vny - tny*DtDs

c        write (6,100) j,i,elten(i),fxg(i,j),fyg(i,j)

       End Do

       fxg(0,j) = fxg(NSG(j),j)
       fyg(0,j) = fyg(NSG(j),j)

       fxg(1,j) = fxg(NSG1(j),j)
       fyg(1,j) = fyg(NSG1(j),j)

      End Do              ! -------- end of Do over interfaces

c---
c Inspect Df
c---

c       Do j=1,Ndrops
c        write (6,*)
c        Do i=0,NSG1(j)
c          write (6,100) j,i,fxg(i,j),fyg(i,j)
c        End Do
c       End Do

c-----------------
c compute the slp 
c-----------------

      Do l=1,Ndrops          ! run over drops
      Do j=1,NSG(l)          ! run over nodes
 
       ja = j-1
       j1 = j+1

        X0 =   xg(j,l)
        Y0 =   yg(j,l)
       TH0 = th2g(j,l)
       FX0 =  fxg(j,l)
       FY0 =  fyg(j,l)

       Utot = 0.0D0
       Vtot = 0.0D0

c--------------------------------------
c non-singular integrals
c will be computed over the blended arc
c--------------------------------------

       Do 20 k=1,Ndrops      ! run over interfaces

c      accm1 = 0     ! for debugging
c      accm2 = 0

       Do 21 i=1,NSG(k)      ! run over elements

       i1 = i+1

       if(l.eq.k) then
        if(j.eq.i.or.j.eq.i1)      Go to 21
        if(j.eq.1.and.i.eq.NSG(k)) Go to 21
       end if

       XC1 =   XCg(i,k)
       YC1 =   YCg(i,k)
       RD1 =    Rg(i,k)
       T11 =  TH2g(i,k)
       T21 =  TH3g(i,k)
       OR1 = ORNTg(i,k)
       FX1 =   FXg(i,k)
       FY1 =   FYg(i,k)

       XC2 =   XCg(i1,k)
       YC2 =   YCg(i1,k)
       RD2 =    Rg(i1,k)
       T12 =  TH1g(i1,k)
       T22 =  TH2g(i1,k)
       OR2 = ORNTg(i1,k)
       FX2 =   FXg(i1,k)
       FY2 =   FYg(i1,k)

       call slp_arc_blended
     +
     +     (XC1,YC1,RD1,T11,T21,OR1
     +     ,XC2,YC2,RD2,T12,T22,OR2
     +     ,FX1,FY1
     +     ,FX2,FY2
     +     ,X0,Y0
     +     ,Iflow
     +     ,Uarc,Varc
     +     ,Istop
     +     )

c---
c      write (6,100) k,i,Uarc,Varc
c---

       Utot = Utot + Uarc
       Vtot = Vtot + Varc

c       accm1 = accm1 + Uarc
c       accm2 = accm2 + Varc

  21  Continue

c     write (6,100) k,i,accm1,accm2

  20  Continue

c---------------------------------------------
c integrate over the two sides of the host arc
c numbered j located at the lth interface
c---------------------------------------------

      Ising = 1

      if(IS_slp.eq.0) Ising = 0

      XCS =   XCg(j,l)
      YCS =   YCg(j,l)
      RDS =    Rg(j,l)
      T1S =  TH1g(j,l)
      T2S =  TH2g(j,l)
      T3S =  TH3g(j,l)
      ORS = ORNTg(j,l)

      FX1 = FXg(ja,l)
      FX2 = FXg( j,l)
      FX3 = FXg(j1,l)

      FY1 = FYg(ja,l)
      FY2 = FYg( j,l)
      FY3 = FYg(j1,l)

c---

      call slp_arc_host
     +
     +     (XCS,YCS,RDS,T1S,T2S,ORS
     +     ,FX1,FY1
     +     ,FX2,FY2
     +     ,X0,Y0,TH0
     +     ,FX0,FY0
     +     ,Ising
     +     ,Iflow
     +     ,Uarc,Varc
     +     ,Istop
     +     )

      Utot = Utot + Uarc
      Vtot = Vtot + Varc

c---
c     write (6,100) l,j,Uarc,Varc
c---

      call slp_arc_host
     +
     +     (XCS,YCS,RDS,T2S,T3S,ORS
     +     ,FX2,FY2
     +     ,FX3,FY3
     +     ,X0,Y0,TH0
     +     ,FX0,FY0
     +     ,Ising
     +     ,Iflow
     +     ,Uarc,Varc
     +     ,Istop
     +     )

      Utot = Utot + Uarc
      Vtot = Vtot + Varc

c---
c     write(6,100) l,j,Uarc,Varc
c---

      Uslpg(j,l) = - Utot/pi4
      Vslpg(j,l) = - Vtot/pi4

c      write (6,100) l,j,Uslpg(j,l),Vslpg(j,l)

      End Do
      End Do           ! end of loop over drops

c-----
c done
c-----

 100  Format (1x,i3,1x,i3,10(1x,f15.10))
 200  Format (' EVALUATING SL potential AT POINT:',I2)
 201  Format (' Arc:',I2)

      return
      end
