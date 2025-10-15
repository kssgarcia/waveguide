      program caps_2d 
 
c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c--------------------------------------------------
c Motion of a wo-dimensional 
c drop or capsule, a file of drops or capsules,
c or a two-dimensional array of drops or capsules,
c for a variety of flow configurations
c
c
c  SYMBOLS:
c  --------
c
c  NSG  :   Number of segments (nodes)
c  x,y  :   coordinates of nodes
c  u,v  :   velocity components
c
c  tinit :   initial constant surface tension
c  cinit :   initial constant surfactant concentration
c  srfin :   initial amount of surfactant
c  elst  :   modulus of elasticity
c
c  srtn :   surface tension at nodes
c
c  suns :   arc length around unstressed shape
c  s    :   arc length around the interface
c
c  c    :   surfactant concentration at nodes
c  cm   :   surfactant concentration at middle nodes
c  Dcds   : d(concentration)/d(arc length)
c
c  vnx,vny : normal vector at nodes
c  vtx,vty : tangential vector at nodes
c
c  c	: surfactant concentration at end-nodes
c  cm	: surfactant concentration at mid-nodes
c  crv	: curvature at end-nodes
c
c  sxy	:	effective shear stress
c  sd1	:	effective first normal-stress difference
c  efv	:	effective viscosity in channel flow
c
c  crvmax:  maximum curvature
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    X(0:900),Y(0:900)
      Dimension    U(0:900),V(0:900)
      Dimension srtn(0:900)
      Dimension    c(0:900),cm(0:900)
      Dimension suns(0:900),Unused(0:900)

      Dimension  XC(900), YC(900),  R(900),   S(900)
      Dimension TH1(900),TH2(900),TH3(900),ORNT(900)

      Dimension vnx(900),vny(900),vtx(900),vty(900)
      Dimension crv(900),dcds(0:900)

      Dimension Xsv(900),Ysv(900)
      Dimension Usv(900),Vsv(900)
      Dimension csv(900)

c---
c related to time stepping
c---

      Dimension time(0:2000),peri(0:2000),epif(0:2000)
      Dimension xcen(0:2000),ycen(0:2000),TAYL(0:2000)
      Dimension crvmax(0:2000)

      Dimension THmx(0:2000),THmn(0:2000)
      Dimension  sxy(0:2000),sd1(0:2000)
      Dimension  efv(0:2000)
 
c---
c Gauss Legendre
c---

      Dimension ZZ(20),WW(20)

c---
c for sgf_2d_2p
c---

      Dimension glut_xx(-64:64,0:64,-64:64)
      Dimension glut_xy(-64:64,0:64,-64:64)
      Dimension glut_yy(-64:64,0:64,-64:64)

c---
c common blocks
c---

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/UUVV/U,V
      common/TENSION/srtn
      common/ELAST/elst,suns

      common/ancR1/thet0,gx,gy
      common/ancR2/vs1,vs2,vsr,vsr1,vsr2,rkap,rkap2,eps
      common/ancR4/Drho
      common/ancR5/RL
      common/ancR6/shrt,yref,pg,ant_c1,ant_c2

      common/hhhh/h,hh,h2,h3,h4,hs
      common/teihos/wall

      common/ancI1/NSG,NSG1,NSG2,NSGM,NSGQ
      common/ancI2/NGFww,IQPD
      common/ancI3/NGL
      common/ancI4/Ivs,Nter
      common/ancI7/Ielst

c---
c for doubly-periodic flow
c---

      common/aaaa/a11,a12,a21,a22
      common/bbbb/b11,b12,b21,b22
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

      common/sgf2d2p/Iglut
      common/glut_r/glut_xx,glut_xy,glut_yy
      common/glut_i/Na21,Nxxx,Nyyy

c---
c various
c---

      common/ZZWW/ZZ,WW
      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi
 
c---
c constants
c---

      pi   = 3.14159 265358 979 32384 D0
      pih  = 0.50D0*pi 
      piq  = 0.25D0*pi
      pi2  = 2.00D0*pi
      pi4  = 4.00D0*pi
      pi6  = 6.00D0*pi
      pi8  = 8.00D0*pi
      srpi = sqrt(pi)
 
      Null = 0
      Ntwo = 2

      zero  = 00.0D0
      three = 3.0D0
      ten   = 10.0D0
      tenm  =-10.0D0
      oth   = 01.0D0/3.0D0
      srth  = sqrt(three)

c------
c input
c------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to read the date from file: caps_2d.dat'
      write (6,*) ' 2 to type data on the keyboard'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*)  Ienrd

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (3,file='caps_2d.dat')

        read (3,*) Iflow
        read (3,*)
        read (3,*) rho1
        read (3,*) rho2
        read (3,*)
        read (3,*) vs1
        read (3,*) vs2
        read (3,*)
        read (3,*) cinit
        read (3,*) tinit
        read (3,*)
        read (3,*) Ds
        read (3,*) betas
        read (3,*) 
        read (3,*) elst
        read (3,*) 
        read (3,*) gac
        read (3,*) thet0
        read (3,*)
        read (3,*) Isolve
        read (3,*)
        read (3,*) eps
        read (3,*) Nter
        read (3,*) Idfl
        read (3,*) JGS
        read (3,*)
        read (3,*) shrt
        read (3,*) U1
        read (3,*) U2
        read (3,*) pg
        read (3,*)
        read (3,*) ant_c1
        read (3,*) ant_c2
        read (3,*) 
        read (3,*) roexp
        read (3,*) mexp
        read (3,*) mexp_foam
        read (3,*)
        read (3,*) wall
        read (3,*)
        read (3,*) h
        read (3,*) Ngfww
        read (3,*) IQPD
        read (3,*)
        read (3,*) IS_slp,IS_dlp
        read (3,*) NGL
        read (3,*)  
        read (3,*) THmax
        read (3,*) SPmin
        read (3,*) SPmax
        read (3,*)  
        read (3,*) RL
        read (3,*) 
        read (3,*) Iread
        read (3,*)  
        read (3,*) Ishape
        read (3,*)  
        read (3,*) aob 
        read (3,*) req
        read (3,*) awob,bwob
        read (3,*)  
        read (3,*) Xdc
        read (3,*) Ydc
        read (3,*) 
        read (3,*) NSG
        read (3,*)
        read (3,*) a11,a12
        read (3,*) a21,a22
        read (3,*) Iglut
        read (3,*) Max1,Max2
        read (3,*)
        read (3,*) Isym
        read (3,*) Iadjust
        read (3,*) Eadj
        read (3,*)
        read (3,*) Norm
        read (3,*)  
        read (3,*) Nprint
        read (3,*)   
        read (3,*) IRK
        read (3,*) Dt
        read (3,*) move
        read (3,*) 

c-----------------------------
      else if(Ienrd.eq.2) then
c-----------------------------

      call verbal
     +
     + (Iflow
     + ,rho1,rho2
     + ,vs1,vs2
     + ,tinit
     + ,cinit
     + ,Ds
     + ,betas
     + ,elst
     + ,gac
     + ,shrt
     + ,ant_c1,ant_c2
     + ,a11,a12
     + ,a21,a22
     + ,Iglut
     + ,Max1,Max2
     + ,roexp,mexp,mexp_foam
     + ,wall
     + ,U1,U2
     + ,pg
     + ,RL
     + ,h,Ngfww,IQPD
     + ,thet0
     + ,IS_slp,IS_dlp
     + ,NGL
     + ,thmax
     + ,spmin,spmax
     + ,Iread
     + ,Ishape
     + ,aob,req
     + ,awob,bwob
     + ,xdc,udc
     + ,NSG
     + ,Isym
     + ,Isolve
     + ,eps,Nter,Idfl,JGS
     + ,Iadjust,Eadj,Norm
     + ,Nprint
     + ,IRK
     + ,Dt
     + ,Istop
     + )

c-----------
      end if     ! end of reading parameters
c-----------


c---
c  display
c
c  Uncomment selected items
c  to verify correctness of input
c---

c       write (i,102) Iflow
c       write (6,100) rho1
c       write (6,100) rho2
c       write (6,100) cinit
c       write (6,100) tinit
c       write (6,100) Ds
c       write (6,100) betas
c       write (6,100) Ielst
c       write (6,100) gac
c       write (6,100) thet0
c       write (6,100) vs1,vs2
c       write (6,112) Isolve
c       write (6,100) eps
c       write (6,112) Nter
c       write (6,112) Idfl
c       write (6,112) JGS
c       write (6,100) shrt
c       write (6,100) ant_c1
c       write (6,100) ant_c2
c       write (6,100) roexp
c       write (6,100) pg
c       write (6,100) H
c       write (6,102) Ngfww
c       write (6,102) IQPD
c       write (6,102) IS_slp,IS_dlp
c       write (6,102) NGL
c       write (6,100) Thmax
c       write (6,100) SPmin,SPmax
c       write (6,100) RL
c       write (6,102) Iread
c       write (6,102) Ishape
c       write (6,100) aob
c       write (6,100) req
c       write (6,100) awob,bwob
c       write (6,100) XDC,YDC
c       write (6,102) NSG
c       write (6,100) a11,a12
c       write (6,100) a21,a22
c       write (6,102) Max1, Max2
c       write (6,102) Isym
c       write (6,102) Iadjust
c       write (6,100) Eadj
c       write (6,102) Norm
c       write (6,102) Nprint
c       write (6,102) IRK
c       write (6,100) Dt
c       write (6,102) move

c---------------------------
c doubly-periodic shear flow
c---------------------------

      if(Iflow.eq.5.or.Iflow.eq.6) then  ! doubly periodic flow

        a11h  = 0.5*a11    ! used for reseting
        a11hm = - a11h     ! the second lattice vector

        box_area = a11*a22-a12*a21     ! area of a periodic box
        fc_rhe = vs1*shrt*box_area     ! for rheology

c---
c read the look-up tables
c---

        if(Iglut.eq.1) then

        open (4,file="sgf_2d_2p.glut")

        write (6,*)
        write (6,*) " caps_2d: Reading the look-up tables"
        write (6,*) "          This may take a while"
        write (6,*) "          Please twiddle your thumbs"
        write (6,*) "          in anticipation"

        read (4,*) Na21,Nxxx,Nyyy

        Do k=-Na21,Na21
          Do j=0,Nyyy
           Do i=-Nxxx,Nxxx
             read (4,*) Glut_xx(i,j,k)
     +                 ,Glut_xy(i,j,k)
     +                 ,Glut_yy(i,j,k)
           End Do
          End Do
        End Do

       close (4)

       end if

      end if

c---------------------
c singly-periodic flow
c---------------------

      if(Iflow.eq.11.or.Iflow.eq.20
     +              .or.Iflow.eq.21
     +              .or.Iflow.eq.22
     +  ) then

        RLM  = -RL
        RL2  =  2.0*RL
        RL2M = -RL2

      end if

c-----------------------
c integration quadrature
c-----------------------

      call gauss_leg (NGL,ZZ,WW)

c-------------------------
c prepare for channel flow
c-------------------------

      if(Iflow.eq.20.or.Iflow.eq.21
     +              .or.Iflow.eq.22) then  ! channel flow

        hh = 0.5*h
        h2 = h+h
        h3 = h2+h
        h4 = h2+h2
        hs = h**2
        hm = -h

        if(Iflow.eq.20) then     ! Couette flow
          DU   = U2-U1
          shrt = DU/(2.0*h)
          yref = -h*(U1+U2)/DU
        end if

      end if

c--------------------------- 
c initial perimeter and area
c--------------------------- 

      if(Ishape.eq.1) then
       peri0 = pi2 * req
       epif0 = pi  * req**2
      else
       req   = 1.0D0
       peri0 = pi2 
       epif0 = pi 
      end if

c------------------------
c additional preparations
c------------------------

      thmax = thmax*pi
      thet0 = thet0*pi        ! inclination angle

      gx =   gac*Dsin(thet0)   ! x-component of gravity
      gy = - gac*Dcos(thet0)   ! y-component of gravity

      NSG1 = NSG+1
      NSG2 = NSG+2

c--------------------
c  density difference
c--------------------

      Drho = rho1-rho2

c--------------------
c viscosity ratio etc
c--------------------

      vsr = vs2/vs1

      vsr1  = 2.0/(1.0D0+vsr)
      vsr2  = 1.0D0-vsr
      rkap  = vsr2/(1.0+vsr)
      rkap2 = vsr1*vsr2/pi4

c---
      if(vsr.gt.999999) then ! if vsr>999999,
                             ! set it equal to infinity
        vsr1  =  2.0D0
        rkap  = -1.0D0
        rkap2 = -1.0D0/pi2
      end if
c---

      Ivs   = 0                       ! Ivs = 0 if vs1 = vs2
      if(abs(vsr2).gt.0.0001) Ivs = 1

c-----------------
c elasticity index
c-----------------

      Ielst = 0                     ! Ielst = 0 for zero elasticity
      if(elst.gt.0.0000001) Ielst = 1

c-----------------------------
c generate the initial shape
c and assign initial velocity
c and surfactant concentration
c-----------------------------

      time(1) = 0.0
      step    = pi2/(NSG1-1.0)

c-------------------------
      if(Ishape.eq.1) then    ! elliptical shape
c-------------------------

      rrx = req * sqrt(aob)
      rry = rrx/aob

      Do i=0,NSG2
        tt   = (i-1.0)*step
        cs   = Dcos(tt)
        sn   = Dsin(tt)
        X(i) = XDC + rrx*cs
        Y(i) = YDC + rry*sn
        U(i) = 0.0D0
        V(i) = 0.0D0
        c(i) = cinit
      End Do

c---------
      else    ! wobbly shape
c---------

      Do i=0,NSG2
        TT   = (i-1.0)*step
        TH   = TT  + awob*Dsin(8.0*TT)
        RR   = 1.0 + bwob*Dcos(4.0*TT)
        X(i) = XDC + RR * Dcos(TH)
        Y(i) = YDC + RR * Dsin(TH)
        U(i) = 0.0D0
        V(i) = 0.0D0
        c(i) = cinit
      End Do

c-----------
      end if
c-----------

c-------------------------------------------
c Compute the arc length around unstressed shape
c to be used for computing elastic tensions
c and total amount of surfactant
c-------------------------------------------

       ICH1 = 0
       ICH2 = 0
       ICH3 = 0

       call prd_2d
     +
     +   (NSG
     +   ,ICH1,THMAX
     +   ,ICH2,SPMAX
     +   ,ICH3,SPMIN
     +   ,peri,epif
     +   ,xcen,ycen
     +   ,U,V,c
     +   ,Suns,Unused
     +   ,Istop
     +   )

       if(Istop.eq.1) Go to 99

       srfin = s(NSG1)*cinit  ! total amount of surfactant

       Do i=0,NSG2
        suns(i) = s(i)
       End Do 

c--------
c prepare
c--------

      if(Isym.eq.0) then           ! non-symmetric shape

        NSGQ = 0
        NSGM = 0
        NPRP = NSG1

      else if(Isym.eq.2) then      ! symmetric shape

        NSGQ = NSG/4+1
        NSGM = NSG/2+1
        NPRP = NSGQ

      end if

c-----------------------------------
c Read initial shape
c      surfactant concentration
c      velocity
c      arc length around the unstressed shape
c
c      from file: caps_2d.inp
c-----------------------------------

      if(Iread.eq.1) then

       open (2,file='caps_2d.inp')

       if(Iflow.eq.5) then
         read (2,*) NPRP,time(1),a21
       else if(Iflow.eq.6) then
         read (2,*) NPRP,time(1),a11,a12,a21,a22
       else
         read (2,*) NPRP,time(1)
       end if

       if(Ielst.eq.1) then
         Do i=1,NPRP
           read (2,*) Idle,X(i),Y(i),U(i),V(i),suns(i),c(i)
        End Do
       else
         Do i=1,NPRP
           read (2,*) Idle,X(i),Y(i),U(i),V(i),c(i)
         End Do
       end if

       if(Iflow.eq.6) then
         Do i=1,NPRP
           X(i) = X(i)*a11
           Y(i) = Y(i)*a11
         End Do
       end if

      close (2)

c---
c prepare and reflect
c---

      if(Isym.eq.0) then

         NSG1 = NPRP
         NSG  = NSG1-1
         NSG2 = NSG+2

      else if(Isym.eq.2) then

         NSGQ = NPRP
         NSG  = 4*(NSGQ-1)
         NSG1 = NSG+1
         NSG2 = NSG+2
         NSGM = NSG/2+1

         Y(1)    = 0.0D0    ! reset to avoid drifting
         X(NSGQ) = 0.0D0    ! reset to avoid drifting

         Do i=1,NSGQ
          i1    = NSGM-i+1
          X   (i1) = -X   (i)
          Y   (i1) =  Y   (i)
          U   (i1) = -U   (i)
          V   (i1) =  V   (i)
          c   (i1) =  c   (i)
          suns(i1) = 2.0D0*suns(NSGQ)-suns(i)
          i1    = NSGM+i-1
          X   (i1) = -X   (i)
          Y   (i1) = -Y   (i)
          U   (i1) = -U   (i)
          V   (i1) = -V   (i)
          c   (i1) =  c   (i)
          suns(i1) = 2.0D0*suns(NSGQ)+suns(i)
          i1    = NSG1-i+1
          X   (i1) =  X   (i)
          Y   (i1) = -Y   (i)
          U   (i1) =  U   (i)
          V   (i1) = -V   (i)
          c   (i1) =  c   (i)
          suns(i1) = 4.0D0*suns(NSGQ)-suns(i)
        End Do

       end if

      end if

c---
c wrap
c---

       X(0) = X(NSG)
       Y(0) = Y(NSG)
       U(0) = U(NSG)
       V(0) = V(NSG)
       c(0) = c(NSG)
       suns(0) = suns(NSG)-suns(NSG1)
  
       X(NSG1) = X(1)
       Y(NSG1) = Y(1)
       U(NSG1) = U(1)
       V(NSG1) = V(1)
       c(NSG1) = c(1)

       X(NSG2) = X(2)
       Y(NSG2) = Y(2)
       U(NSG2) = U(2)
       V(NSG2) = V(2)
       c(NSG2) = c(2)
       suns(NSG2) = suns(NSG1)+suns(2)

c----------------------
c Display initial shape
c on screen
c----------------------

       write (6,*)
       write (6,*)  " Starting shape:"
       write (6,*)

       Do i=0,NSG2
        if(Ielst.eq.1) then
         write (6,120) i,X(i),Y(i),U(i),V(i),suns(i)
        else
         write (6,120) i,X(i),Y(i),U(i),V(i),c(i)
        end if
       End Do

c---------------------
c prepare output files
c---------------------

      open (1,file='caps_2d.xy')    ! for profiles
      open (4,file='caps_2d.diag')  ! for deformation and orientation
      open (7,file="caps_2d.rhe")   ! for rheology

c---------------------------------
c plot walls or channel boundaries
c---------------------------------

      if(Iflow.eq.10.or.Iflow.eq.11) then   ! wall-bounded flow

        write (1,103) Ntwo
        write (1,102) None,tenm,wall
        write (1,102) Ntwo,ten ,wall

      else if(Iflow.eq.20.or.Iflow.eq.21    ! channel flow
     +                   .or.Iflow.eq.22
     +       ) then

        write (1,103) Ntwo
        write (1,102) None,RLm,hm
        write (1,102) Ntwo,RL ,hm
        write (1,103) Ntwo
        write (1,102) None,RLm,h
        write (1,102) Ntwo,RL ,h

      end if

c-.-.-.-.-.-.-.-.-.-.-
c BEGIN THE SIMULATION
c-.-.-.-.-.-.-.-.-.-.-

      Kstep  = 1        ! Global time step counter

      Iprint = Nprint   ! Printing counter

      Itry = 1          ! graphics

  90  Continue

      if(Ienrd.eq.0) then
        write (6,*)
        write (6,*) 'Enter number of steps before pausing'
        write (6,*)
        write (6,*) '0 to quit'
        write (6,*) '---------'
        read  (5,*)  Nstep
      else
        read  (3,*)   Nstep
      end if

      if(Nstep.eq.0) Go to 99

c-----------
c initialize
c-----------

      Istep = 1       ! Batch step counter

  97  Continue

      write (6,*)
      write (6,*) "--------------------"
      write (6,111) Istep,Nstep

c-------------------------------
c check the point distribution
c compute properties of arcs
c-------------------------------

      ICH1 = 1
      ICH2 = 1
      ICH3 = 1

      if(Isym.eq.0) then

        call prd_2d
     +
     +    (NSG
     +    ,ICH1,THMAX
     +    ,ICH2,SPMAX
     +    ,ICH3,SPMIN
     +    ,peri(Kstep),epif(Kstep)
     +    ,xcen(Kstep),ycen(Kstep)
     +    ,U,V,c
     +    ,suns,unused
     +    ,Istop
     +    )

        NSG1 = NSG+1
        NSG2 = NSG+2
        NPRP = NSG1

      else if(Isym.eq.2) then

        call prd_2d_qs
     +
     +     (NSG
     +     ,ICH1,THMAX
     +     ,ICH2,SPMAX
     +     ,ICH3,SPMIN
     +     ,peri(Kstep),epif(Kstep)
     +     ,xcen(Kstep),ycen(Kstep)
     +     ,U,V,srtn,suns,unused
     +     ,Istop
     +     )

        NSG1 = NSG+1
        NSG2 = NSG+2
        NSGQ = NSG/4+1
        NSGM = NSG/2+1
        NPRP = NSGQ

      end if

      if(Istop.eq.1) Go to 99

c---
c wrap
c---

      X(0) = X(NSG)
      Y(0) = Y(NSG)
      U(0) = U(NSG)
      V(0) = V(NSG)
      c(0) = c(NSG)
      suns(0) = suns(NSG)-suns(NSG1)

      X(NSG1) = X(1)
      Y(NSG1) = Y(1)
      U(NSG1) = U(1)
      V(NSG1) = V(1)
      c(NSG1) = c(1)

      X(NSG2) = X(2)
      Y(NSG2) = Y(2)
      U(NSG2) = U(2)
      V(NSG2) = V(2)
      c(NSG2) = c(2)
      suns(NSG2) = suns(NSG1)+suns(2)

c-----------------------------
c normalize the capsule area
c by expanding or contracting
c with respect to the centroid
c-----------------------------

      if(Norm.eq.1) then

      write (6,*)
      write (6,*) " Normalizing capsule area"
      write (6,*) " ------------------------"

      scale  = 1.0/sqrt(epif(Kstep)/epif0)

      Do i=0,NSG2
        X(i) = xcen(Kstep) + (X(i)-xcen(Kstep))*scale
        Y(i) = ycen(Kstep) + (Y(i)-ycen(Kstep))*scale
      End Do

      ICH1 = 0
      ICH2 = 0
      ICH3 = 0

      if(Isym.eq.0) then

        call prd_2d
     +
     +  (NSG
     +  ,ICH1,THMAX
     +  ,ICH2,SPMAX
     +  ,ICH3,SPMIN
     +  ,peri(Kstep),epif(Kstep)
     +  ,xcen(Kstep),ycen(KStep)
     +  ,U,V,c
     +  ,Suns,Unused
     +  ,Istop
     +  )

      else if(Isym.eq.2) then

        call prd_2d_qs
     +
     +  (NSG
     +  ,ICH1,THMAX
     +  ,ICH2,SPMAX
     +  ,ICH3,SPMIN
     +  ,peri(Kstep),epif(Kstep)
     +  ,xcen(Kstep),ycen(Kstep)
     +  ,U,V,c
     +  ,Suns,Unused
     +  ,Istop
     +  )

      end if

      if(Istop.eq.1) Go to 99

      end if

c---
c wrap
c---

      X(0) = X(NSG)
      Y(0) = Y(NSG)
      U(0) = U(NSG)
      V(0) = V(NSG)
      c(0) = c(NSG)
      suns(0) =-suns(2)

      X(NSG1) = X(1)
      Y(NSG1) = Y(1)
      U(NSG1) = U(1)
      V(NSG1) = V(1)
      c(NSG1) = c(1)
      suns(NSG1) = suns(NSG1)+suns(1)

      X(NSG2) = X(2)
      Y(NSG2) = Y(2)
      U(NSG2) = U(2)
      V(NSG2) = V(2)
      c(NSG2) = c(2)
      suns(NSG2) = suns(NSG1)+suns(2)

c---------------------------------
c Prepare for the
c doubly-periodic Green's function
c---------------------------------

      if((Iflow.eq.5..and.Iglut.eq.0)
     +                .or.Iflow.eq.6) then

        call sgf_2d_2p_ewald
     +
     +   (a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,tau
     +   )

        call sgf_2d_2p_qqq (b11,b12,b21,b22,Max2,ew)

        if(Ivs.eq.1) then
         call sgf_2d_2p_vvv (b11,b12,b21,b22,Max2,ew)
        end if

      end if

c----------------------------------
c scale arc length, area, centroid
c and print
c----------------------------------

      if(Iflow.eq.6) then
       if(a21.lt.0.000001) then
           peric = 4.0*a11        ! square cell
       else
           peric = 6.0*a11/srth   ! hexagonal cell
       end if
       peri(Kstep) = peri(Kstep)/peric
       epif(Kstep) = epif(Kstep)/tau
      else
       peri(Kstep) = peri(Kstep)/peri0
       epif(Kstep) = epif(Kstep)/epif0
      end if

      write (6,305) peri(Kstep) 
      write (6,306) epif(Kstep)
      write (6,308) xcen(Kstep),ycen(Kstep)
  
c-------------------------------------
c Compute Taylor deformation parameter
c         capsile inclination
c         and maximum curvature
c
c Anmax: inclination of major axis
c Anmin: inclination of minor axis
c------------------------------------

      if(Isym.eq.2) then

        Dxy   = (X(1)-Y(NSGQ))/(X(1)+Y(NSGQ))
        ANmax = 0.0
        ANmin = pih

c       crvmax(Kstep) = ORNT(1)/R(1)   ! maximum curvature

        call crv_max (crvmax(Kstep))

      else

        call Taylor 
     +
     +    (xcen,ycen
     +    ,Dxy
     +    ,ANmax,ANmin
     +    )

        call crv_max 
     +
     +    (crvmax(Kstep)
     +    )

      end if

      tayl(Kstep) = Dxy
      thmx(Kstep) = ANmax
      thmn(Kstep) = ANmin

c---------------------
c Aitken extrapolation
c for the deformation
c---------------------

      if(Kstep.ge.3) then

        q0  = TAYL(Kstep-2)
        q1  = TAYL(Kstep-1)
        q2  = TAYL(Kstep)
        AIT = q2-(q2-q1)**2/(q2-2.0*q1+q0)

      end if

c-----------------
c printing session
c-----------------

      if(Isym.eq.2) then

       if(Iflow.eq.6) then
         ftc = (a11-2.0*x(1))/a11    ! film thickness
         crvmax(Kstep) = crvmax(Kstep)*a11
         A_L = (1.0-epif(Kstep))*tau
         write (6,*) " Amount of liquid:", A_L
       end if

       write (4,145) Kstep,time(Kstep)
c    +              ,epif(Kstep),peri(Kstep)
c    +              ,ftc
c    +              ,A_L
c    +              ,a11
     +              ,TAYL(Kstep),AIT
     +              ,crvmax(Kstep)
c      write (6,145) Kstep,time(Kstep)
c    +              ,epif(Kstep),peri(Kstep)
c    +              ,ftc
c    +              ,crvmax(Kstep)
c    +              ,A_L
c    +              ,TAYL(Kstep),AIT

      else

       write (4,140) Kstep,time(Kstep)
c    +              ,epif(Kstep),peri(Kstep)
c    +              ,ftc
c    +              ,crvmax(Kstep)
c    +              ,A_L
c    +              ,a11
     +              ,TAYL(Kstep),AIT
c    +              ,THmx(Kstep),THmn(Kstep)
       write (6,140) Kstep,time(Kstep)
c    +              ,epif(Kstep),peri(Kstep)
c    +              ,ftc
c    +              ,crvmax(Kstep)
c    +              ,A_L
     +              ,TAYL(Kstep),AIT
c    +              ,THmx(Kstep),THmn(Kstep)

      end if

c---------------------------------------------
c adjust the time step according the curvature
c in two alternative ways
c---------------------------------------------

c-----------
      if(Iadjust.eq.1) then
c-----------

        if(Kstep.eq.1) then
         adjust = 1.0D0
        else
         adjust = abs(crvmax(Kstep)-crvmax(Kstep-1))
        end if

        if(adjust.gt.Eadj) Dt = Dt/adjust

        write (6,*) "Adjust = ",adjust

c-----------
      else if(Iadjust.eq.2) then
c-----------

        if(Kstep.eq.1) then
         adjust = 1.0
        else
         adjust = crvmax(Kstep-1)/crvmax(Kstep)
        end if

        if(adjust.le.Eadj) Dt = Dt*adjust

        write (6,*) "Adjust = ",adjust

c-----------
      end if
c-----------

c-----------------------------------------------  
c Compute the surface tension from concentration   
c at the end-nodes, from a constitutive equation
c
c Note that when betas = 0, srtn(i)=tinit
c-----------------------------------------------   

      Do i=0,NSG2
       srtn(i) = tinit*(1.0D0-betas*c(i)/cinit)
     +                /(1.0D0-betas)
      End Do

c---------------------
c Doubly-periodic shear flow:
c reset the second lattice vector
c---------------------

      if(Iflow.eq.5) then
        if(a21.gt.a11h)  a21 = a21-a11
        if(a21.lt.a11hm) a21 = a21+a11
      end if

c---------------------
c compute the velocity
c at the marker points
c---------------------

      call caps_2d_vel
     +
     +  (IS_slp,IS_dlp
     +  ,Isolve,JGS
     +  ,Iflow,Isym,Idfl
     +  ,roexp,mexp,mexp_foam
     +  ,Istop
     +  )

      if(Istop.eq.1) Go to 99

c----------------------------------
c RHEOLOGY
c
c Compute the effective particle stress
c tensor  in the dilute limit
c
c Drop in unbounded flow only
c----------------------------------

c-----
      if(Iflow.eq.1.or.Iflow.eq.2
     +             .or.Iflow.eq.5
     +             .or.Iflow.eq.10
     +             .or.Iflow.eq.20
     +             ) then
c-----

        call caps_2d_rhe (strxx,strxy,stryy)

c-----
        if(Iflow.eq.1.or.Iflow.eq.2
     +              .or.Iflow.eq.10
     +              .or.Iflow.eq.20
     +    ) then
c-----

         sxy(Kstep) = strxy/(shrt*vs1)
         sd1(Kstep) = (strxx-stryy)/(shrt*vs1)

        else if(Iflow.eq.5) then       ! doubly-periodic flow

         sxy(Kstep) =  strxy       /fc_rhe
         sd1(Kstep) = (strxx-stryy)/fc_rhe

        end if

        write (6,*)
        write (6,*) " time and rheological properties:"
        write (6,141) Kstep,time(Kstep),sxy(Kstep),sd1(Kstep)
        write (7,141) Kstep,time(Kstep),sxy(Kstep),sd1(Kstep)
        write (6,*) " --------------------------------"

c-----
      end if
c-----

c----------------------------------
c Pressure-driven flow in a channel
c with viscosity ratio equal to unity
c
c Compute the disturbance flow rate for
c vanishing pressure drop (IQPD = 0)
c
c or the disturbance pressure drop for
c vanishing disturbance flow rate (IQPD = 1)
c 
c Define the effective viscosity
c----------------------------------

      if(Ivs.eq.0) then

       if(Iflow.eq.21.or.Iflow.eq.22) then

        call channel (IQPD,RL,Prgr,Flrt)

        if(Iflow.eq.21) then
          flrt = 2.0/3.0*pg*h**3/vs1 + flrt/vs1 
          prgr = pg                  + prgr
        else if(Iflow.eq.22) then
          flrt = 1.0
          prgr = 1.0
        end if

        efv(Kstep)= 2.0/3.0*h**3*prgr/flrt

        write (6,141) Kstep,time(Kstep),efv(Kstep)
        write (7,141) Kstep,time(Kstep),efv(Kstep)

       end if

      end if

c-----------------
c printing session
c-----------------

      if(Iprint.eq.Nprint) then

        if(Iflow.eq.5) then
         write (1,103) NPRP,time(Kstep),a21,Dt
c        write (6,103) NPRP,time(Kstep),a21,Dt
        else if(Iflow.eq.6) then
         write (1,103) NPRP,time(Kstep),a11,a12,a21,a22,Dt
c        write (6,103) NPRP,time(Kstep),a11,a12,a21,a22,Dt
        else
         write (1,103) NPRP,time(Kstep),Dt
c        write (6,103) NPRP,time(Kstep),Dt
        end if

        if(Ielst.eq.1) then
          Do i=1,NPRP
            Xpr = X(i)
            Ypr = Y(i)
            If(Iflow.eq.6) Xpr = Xpr/a11
            If(Iflow.eq.6) Ypr = Ypr/a11
            write (1,120) i,Xpr,Ypr,U(i),V(i),suns(i)
c           write (6,120) i,Xpr,Ypr,U(i),V(i),suns(i)
          End Do
        else
          Do i=1,NPRP
            Xpr = X(i)
            Ypr = Y(i)
            If(Iflow.eq.6) Xpr = Xpr/a11
            If(Iflow.eq.6) Ypr = Ypr/a11
            write (1,120) i,Xpr,Ypr,U(i),V(i),c(i)
c           write (6,120) i,Xpr,Ypr,U(i),V(i),c(i)
          End Do
        end if

       Iprint = 0

      end if

c--------------------------------------------------
c The following block integrates the 
c convection-diffusion equation for the surfactant
c concentration using a finite-volume method
c
c Block is bypassed if betas = 0, 
c for then the surface tension remains constant
c---------------------------------------------------

      If(betas.gt.0.0000001) then

c------------------------------------------
c At the first step,
c interpolate concentration at middle-nodes
c from end-nodes
c------------------------------------------

      if(Kstep.eq.1) call interp_mn (NSG,c,s,cm)

c-------------------------------
c compute normal vector,
c         tangential vectors
c         curvature at end-nodes
c-------------------------------

      Do i=1,NSG
        vnx(i) = Dcos(TH2(i)) * ornt(i)
        vny(i) = Dsin(TH2(i)) * ornt(i)
        vtx(i) = - vny(i)
        vty(i) =   vnx(i)
        crv(i) = ornt(i)/R(i)
      End Do

c---
c wrap
c---

      vnx(NSG1) = vnx(1)
      vny(NSG1) = vny(1)
      vtx(NSG1) = vtx(1)
      vty(NSG1) = vty(1)
      crv(NSG1) = crv(1)

      vnx(NSG2) = vnx(2)
      vny(NSG2) = vny(2)
      vtx(NSG2) = vtx(2)
      vty(NSG2) = vty(2)
      crv(NSG2) = crv(2)

c---
c compute and print tangential velocity
c---

c     If(Iflow.eq.5) then
c      write (1,103) NPRP,time(Kstep),a21,Dt
c     Else If(Iflow.eq.6) then
c      write (1,103) NPRP,time(Kstep),a11,a12,a21,a22,Dt
c     Else
c      write (6,103) NPRP,time(Kstep),Dt
c     End If
c
c     Do i=1,NPRP
c      tanvel = u(i)*vtx(i) + v(i)*vty(i)
c      If(Isym.eq.0) then
c       spr = s(i)/s(NSG1)
c      Else if(Isym.eq.2) then
c       spr = 1.0-s(i)/s(NSGQ)
c      End If
c      write(1,120) i,spr,tanvel
c     End Do

c-----------------------------------------
c  compute the derivative d(c)/d(l)
c  at end-nodes, by parabolic interpolation
c  with respect to arc length
c-----------------------------------------

      Do i=2,NSG1
       x0      = s(i-1)-s(i)
       x2      = s(i+1)-s(i)
       y0      = c(i-1)-c(i)
       y2      = c(i+1)-c(i)
       dcds(i) = (x0/x2*y2-x2/x0*y0)/(x0-x2)
      End Do

c---
c extend
c---

      dcds(0)    = dcds(NSG)
      dcds(1)    = dcds(NSG1)
      dcds(NSG2) = dcds(2)

c---
c update the mean surfactant concentration
c at mid-nodes using
c an implicit finite-volume method
c---

      call caps_2d_cd
     +
     +   (Move
     +   ,Isym
     +   ,vnx,vny
     +   ,vtx,vty
     +   ,cm,c,crv
     +   ,Ds
     +   ,Dt
     +   ,srfam      ! total amount of surfactant
     +   ,Istop
     +   )

      if(Istop.eq.1) Go to 99

c-----------------------------
c normalize the concentration
c of a surfactant
c-----------------------------

      fc = srfin/srfam

      Do i=0,NSG2
        c(i) = fc* c(i)
       cm(i) = fc*cm(i)
      End Do

c---------------------------------------
c interpolate concentration at end-nodes
c from mid-nodes
c
c end-node concentration is used to compute
c the velocity
c---------------------------------------

      call interp_en (NSG,cm,s,c)

      end if                           ! End of surfactant module

c-------------------------------------------

c------------------------
c Runge--Kutta integration
c------------------------

c----------------------
      if(IRK.eq.1) then    ! first order
c----------------------

      if(Move.eq.0) then         ! points move with total velocity

        Do i=1,NPRP
          X(i) = X(i) + Dt*U(i)
          Y(i) = Y(i) + Dt*V(i)
        End Do

      else             ! points move with normal velocity

        Do i=1,NSG
          cs   = cos(th2(i))
          sn   = sin(th2(i))
          Utot = U(i)*cs+V(i)*sn
          Umv  = Utot*cs
          Vmv  = Utot*sn
          X(i) = X(i) + Dt*Umv
          Y(i) = Y(i) + Dt*Vmv
        End Do

      end if

      if(Iflow.eq.5) then    ! Doubly-periodic shear flow
        a21 = a21 + Dt*shrt*a22
      else if(Iflow.eq.6) then
        STRETCH = 1.0D0+Dt*0.5D0*roexp/tau
        a11 = a11 * STRETCH
        a12 = a12 * STRETCH
        a21 = a21 * STRETCH
        a22 = a22 * STRETCH
      end if

c----------------------
      else if(IRK.eq.2) then    ! second order
c----------------------

      if(Move.eq.0) then

        Do i=1,NPRP
             Xsv(i) = X(i)
             Ysv(i) = Y(i)
             csv(i) = c(i)
             Usv(i) = U(i)
             Vsv(i) = V(i)
c         srtnsv(i) = srtn(i)
               X(i) = X(i) + Dt*U(i)
               Y(i) = Y(i) + Dt*V(i)
        End Do

       else

        Do i=1,NSG
          Xsv(i) = X(i)
          Ysv(i) = Y(i)
          csv(i) = c(i)
          cs     = cos(th2(i))
          sn     = sin(th2(i))
          Utot   = U(i)*cs+V(i)*sn
          Umv    = Utot*cs
          Vmv    = Utot*sn
          Usv(i) = Umv
          Vsv(i) = Vmv
            X(i) = X(i) + Dt*Umv
            Y(i) = Y(i) + Dt*Vmv
        End Do

       end if

       if(Iflow.eq.5) then          ! Doubly-periodic shear flow
        a21 = a21 + Dt*shrt*a22
       else if(Iflow.eq.6) then     ! Expanding foam
        STRETCH = 1.0+Dt*0.5*roexp/tau
        a11 = a11 * STRETCH
        a12 = a12 * STRETCH
        a21 = a21 * STRETCH
        a22 = a22 * STRETCH
       end if

c---
c reflect a symmetric interface
c---

      if(Isym.eq.2) then

        Y(1)    = 0.0D0   ! reset to avoid drifting
        X(NSGQ) = 0.0D0   ! reset to avoid drifting

        Do i=1,NSGQ
          j = NSGM-i+1
             X(j) = -X(i)
             Y(j) =  Y(i)
             U(j) = -U(i)
             V(j) =  V(i)
             c(j) =  c(i)
          srtn(j) = srtn(i)
          suns(j) = 2.0*suns(NSGQ)-suns(i)
          j = NSGM+i-1
          X   (j) = -X(i)
          Y   (j) = -Y(i)
          U   (j) = -U(i)
          V   (j) = -V(i)
          c   (j) =  c(i)
          srtn(j) = srtn(i)
          suns(j) = 2.0*suns(NSGQ)+suns(i)
          j = NSG1-1+1
             X(j) =  X(i)
             Y(j) = -Y(i)
             U(j) =  U(i)
             V(j) = -V(i)
             c(j) =  c(i)
          srtn(j) = srtn(i)
          suns(j) = 4.0*suns(NSGQ)-suns(i)
        End Do

      end if

c---
c wrap
c---

      X(0) = X(NSG)
      Y(0) = Y(NSG)
      U(0) = U(NSG)
      V(0) = V(NSG)
      c(0) = c(NSG)

      srtn(0) = srtn(NSG)
      suns(0) = suns(NSG1)-suns(NSG)

      X(NSG1) = X(1)
      Y(NSG1) = Y(1)
      U(NSG1) = U(1)
      V(NSG1) = V(1)
      c(NSG1) = c(1)

      srtn(NSG1) = srtn(1)
      suns(NSG1) = suns(NSG1)+suns(1)

      X(NSG2) = X(2)
      Y(NSG2) = Y(2)
      U(NSG2) = U(2)
      V(NSG2) = V(2)
      c(NSG2) = c(2)

      srtn(NSG2) = srtn(2)
      suns(NSG2) = suns(NSG1)+suns(2)

c---
c Geometry at intermediate step
c---

      ICH1 = 0   ! no redistribution
      ICH2 = 0
      ICH3 = 0

      if(Isym.eq.0) then

        call prd_2d
     +
     +     (NSG
     +     ,ICH1,THMAX
     +     ,ICH2,SPMAX
     +     ,ICH3,SPMIN
     +     ,peri,epif,xcen,ycen
     +     ,U,V,c
     +     ,Suns,Unused
     +     ,Istop
     +     )

      else if(Isym.eq.2) then

        call prd_2d_qs
     +
     +     (NSG
     +     ,ICH1,THMAX
     +     ,ICH2,SPMAX
     +     ,ICH3,SPMIN
     +     ,peri,epif,xcen,ycen
     +     ,U,V,c
     +     ,Suns,Unused
     +     ,Istop
     +     )

      end if

      if(Istop.eq.1) Go to 99

c---
c  prepare for the doubly-periodic Green's function
c---

      if(Iflow.eq.5.or.Iflow.eq.6) then

        if(Iglut.eq.0) then

        call sgf_2d_2p_ewald
     +
     +     (a11,a12,a21,a22
     +     ,b11,b12,b21,b22
     +     ,ew,tau
     +     )

        call sgf_2d_2p_qqq (b11,b12,b21,b22,Max2,ew)

        if(Ivs.eq.1) then
         call sgf_2d_2p_vvv (b11,b12,b21,b22,Max2,ew)
        end if

        end if

      end if

c-----------------------------------
c Surface tension from concentration
c at end-nodes
c-----------------------------------

      Do i=0,NSG2
       srtn(i) = tinit*(1.0D0-betas*c(i)/cinit)
     +                /(1.0D0-betas)   ! SRF
      End Do

c---------------------
c Doubly-periodic shear flow:
c reset the second lattice vector
c---------------------

      if(Iflow.eq.5) then
        if(a21.gt.a11h)  a21 = a21-a11
        if(a21.lt.a11hm) a21 = a21+a11
      end if

c---
c Second velocity evaluation
c---

      call caps_2d_vel
     +
     +  (IS_slp,IS_dlp
     +  ,Isolve,JGS
     +  ,Iflow,Isym,Idfl
     +  ,roexp,mexp,mexp_foam
     +  ,Istop
     +  )

      if(Istop.eq.1) Go to 99

c---
c second step in RK2
c---

      Dth = 0.5*Dt

      if(Move.eq.0) then
        Do i=1,NSG
          X(i) = Xsv(i)+Dth*(Usv(i)+U(i))
          Y(i) = Ysv(i)+Dth*(Vsv(i)+V(i))
        End Do
      else
        Do i=1,NSG
          cs   = cos(th2(i))
          sn   = sin(th2(i))
          Utot = U(i)*cs+V(i)*sn
          Umv  = Utot*cs
          Vmv  = Utot*sn
          X(i) = Xsv(i)+Dth*(Usv(i)+Umv)
          Y(i) = Ysv(i)+Dth*(Vsv(i)+Vmv)
        End Do
      end if

c-----------
      end if          ! End of RK2
c-----------

c-------------------------------------------

c-------------------
c End of a time step
c-------------------

  96  Continue

c---
c reflect a symmetric interface
c---

      if(Isym.eq.2) then

        Y(1)    = 0.0D0   ! reset to avoid drifting
        X(NSGQ) = 0.0D0   ! reset to avoid drifting

        Do i=1,NSGQ
          j=NSGM-i+1
          X(j) = -X(i)
          Y(j) =  Y(i)
          U(j) = -U(i)
          V(j) =  V(i)
          c(j) =  c(i)
          srtn(j) = srtn(i)
          suns(j) = 2.0*suns(NSGQ)-suns(i)
          j = NSGM+i-1
          X(j) = -X(i)
          Y(j) = -Y(i)
          U(j) = -U(i)
          V(j) = -V(i)
          c(j) =  c(i)
          srtn(j) = srtn(i)
          suns(j) = 2.0*suns(NSGQ)+suns(i)
          j = NSG1-1+1
          X(j) =  X(i)
          Y(j) = -Y(i)
          U(j) =  U(i)
          V(j) = -V(i)
          c(j) =  c(i)
          srtn(j) = srtn(i)
          suns(j) = 4.0*suns(NSGQ)-suns(i)
        End Do

      end if

c---
c wrap
c---

      X(0) = X(NSG)
      Y(0) = Y(NSG)
      U(0) = U(NSG)
      V(0) = V(NSG)
      c(0) = c(NSG)
      srtn(0) = srtn(NSG)
      suns(0) = suns(NSG1)-suns(NSG)

      X(NSG1) = X(1)
      Y(NSG1) = Y(1)
      U(NSG1) = U(1)
      V(NSG1) = V(1)
      c(NSG1) = c(1)
      srtn(NSG1) = srtn(1)
      suns(NSG1) = suns(NSG1)+suns(1)

      X(NSG2) = X(2)
      Y(NSG2) = Y(2)
      U(NSG2) = U(2)
      V(NSG2) = V(2)
      c(NSG2) = c(2)
      srtn(NSG2) = srtn(2)
      suns(NSG2) = suns(NSG1)+suns(2)

c---
c reset counters and time
c---

      Istep  = Istep +1
      Kstep  = Kstep +1
      Iprint = Iprint+1

      time(Kstep) = time(Kstep-1) + Dt

      if(Istep.le.Nstep) Go to 97

      Go to 90

c-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
c  Simulation has ended           c
c-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

  99  Continue

c----------------------
c finish up by printing
c----------------------

      open (9,file="caps_2d.rst")     ! restart file

      if(Iflow.eq.5) then
       write (1,103) NPRP,time(Kstep),a21,Dt
       write (6,103) NPRP,time(Kstep),a21,Dt
       write (9,103) NPRP,time(Kstep),a21,Dt
      else if(Iflow.eq.6) then
       write (1,103) NPRP,time(Kstep),a11,a12,a21,a22,Dt
       write (6,103) NPRP,time(Kstep),a11,a12,a21,a22,Dt
       write (9,103) NPRP,time(Kstep),a11,a12,a21,a22,Dt
      else
       write (1,103) NPRP,time(Kstep),Dt
       write (6,103) NPRP,time(Kstep),Dt
       write (9,103) NPRP,time(Kstep),Dt
      end if

      if(Ielst.eq.1) then
        Do j=1,NPRP
         Xpr = X(j)
         Ypr = Y(j)
         If(Iflow.eq.6) Xpr = Xpr/a11
         If(Iflow.eq.6) Ypr = Ypr/a11
         write (1,120) j,Xpr,Ypr,U(j),V(j),suns(j),c(j)
c        write (6,120) j,Xpr,Ypr,U(j),V(j),suns(j),c(j)
         write (9,120) j,Xpr,Ypr,U(j),V(j),suns(j),c(j)
        End Do
      else
        Do j=1,NPRP
         Xpr = X(j)
         Ypr = Y(j)
         If(Iflow.eq.6) Xpr = Xpr/a11
         If(Iflow.eq.6) Ypr = Ypr/a11
         write (1,120) j,Xpr,Ypr,U(j),V(j),c(j)
c        write (6,120) j,Xpr,Ypr,U(j),V(j),c(j)
         write (9,120) j,Xpr,Ypr,U(j),V(j),c(j)
        End Do
      end if

      write (1,*) Null,Null,Null,Null,Null,Null,Null,Null
      write (9,*) Null,Null,Null,Null,Null,Null,Null,Null

      close (9)

c---
c compute final centroid and volume
c---

      ICH1 = 0
      ICH2 = 0
      ICH3 = 0

      if(Isym.eq.0) then

        call prd_2d
     +
     +    (NSG
     +    ,ICH1,THMAX
     +    ,ICH2,SPMAX
     +    ,ICH3,SPMIN
     +    ,peri(Kstep),epif(Kstep)
     +    ,xcen(Kstep),ycen(Kstep)
     +    ,U,V,SRTN,Suns,Unused
     +    ,Istop
     +    )

      else if(Isym.eq.2) then

        call prd_2d_qs
     +
     +    (NSG
     +    ,ICH1,THMAX
     +    ,ICH2,SPMAX
     +    ,ICH3,SPMIN
     +    ,peri(Kstep),epif(Kstep)
     +    ,xcen(Kstep),ycen(Kstep)
     +    ,U,V,SRTN,Suns,Unused
     +    ,Istop
     +    )

      end if

      peri(Kstep) = peri(Kstep)/peri0
      epif(Kstep) = epif(Kstep)/epif0

      write (6,305) peri(Kstep) 
      write (6,306) epif(Kstep)
      write (6,308) xcen(Kstep),ycen(Kstep)

c---
c record diagnostics
c---

      Do i=1,Kstep
        write (1,109) time(i),epif(i),peri(i),xcen(i),ycen(i)
        write (4,109) time(i),epif(i),peri(i),xcen(i),ycen(i)
        write (6,109) time(i),epif(i),peri(i),xcen(i),ycen(i)
      End Do

c----------------------
c record run parameters
c----------------------

      thet0 = thet0/pi   ! unscale
      thmax = thmax/pi   ! unscale

      write (1,105) 
     +
     +  Iflow,Isolve,Isym,RL,AoB,req,XDC,YDC
     +  ,shrt,pg,ant_c1,ant_c2
     +  ,cinit,tinit,Ds,betas
     +  ,elst
     +  ,rho1,rho2,g,vsr
     +  ,eps,Nter,alp,Norm
     +  ,Iadjust,Eadj
     +  ,roexp,mexp,mexp_foam

      write (1,107) THmax,SPmax,SPmin
      write (1,108) IS_slp,IS_dlp,NGL,Dt,IRK,MOVE,IDFL,JGS

      If(Iflow.le.9) write (1,106) h,Ngfww,IQPD

      write (4,105) 
     +
     +  Iflow,Isolve,Isym,RL,AoB,req,XDC,YDC
     +  ,shrt,pg,ant_c1,ant_c2
     +  ,cinit,tinit,Ds,betas
     +  ,elst
     +  ,rho1,rho2,g,vsr
     +  ,eps,Nter,alp,Norm
     +  ,Iadjust,Eadj
     +  ,roexp,mexp,mexp_foam

      write (4,107) THmax,SPmax,SPmin
      write (4,108) IS_slp,IS_dlp,NGL,Dt,IRK,MOVE,IDFL,JGS

      If(Iflow.le.9) write (4,106) h,Ngfww,IQPD

      write (7,105) 
     +
     +  Iflow,Isolve,Isym,RL,AoB,req,XDC,YDC
     +  ,shrt,pg,ant_c1,ant_c2
     +  ,cinit,tinit,Ds,betas
     +  ,elst
     +  ,rho1,rho2,g,vsr
     +  ,eps,Nter,alp,Norm
     +  ,Iadjust,Eadj
     +  ,roexp,mexp,mexp_foam
      write (7,107) THmax,SPmax,SPmin
      write (7,108) IS_slp,IS_dlp,NGL,Dt,IRK,MOVE,IDFL,JGS

      If(Iflow.le.9) write (7,106) h,Ngfww,IQPD

c--------
c wrap up
c--------

      write (1,110) 
      write (4,*) Null,Null,Null,Null,Null,Null,Null

c------------
c close files
c------------

      close (1)
      If(Ienrd.eq.1) close (3)
      close (4)
      close (7)
 
c-----
c Done
c-----

 100  Format (1X,F10.5,2X,F10.5,10X,F10.5,2X,F10.5)
 101  Format (1X,'U = ',F15.10,'  V = ',F15.10)
 102  Format (1X,I4,5(1X,F12.8))
 120  Format (1X,I4,2(1X,F12.8),2(1X,F8.4),2(1X,F12.8))
 103  Format (1X,I4,9(1X,F10.6))
 104  Format (1X,I4,5(1X,F12.8))
 105  Format (/,
     +        ' Iflow  = ',I2,/,/,
     +        ' Isolve = ',I2,/,/,
     +        ' Isym   = ',I1,/,/,
     +        ' L      = ',F7.4,/,
     +        ' A/B    = ',F7.4,/,
     +        ' Req    = ',F7.4,/,
     +        ' XDC    = ',F7.4,/,
     +        ' YDC    = ',F7.4,/,/,
     +        ' shrt   = ',F7.4,/,
     +        ' pg     = ',F7.4,/,/,
     +        ' c1_ant = ',F7.4,/,
     +        ' c2_ant = ',F7.4,/,/,
     +        ' cinit  = ',F7.4,/,
     +        ' tinit  = ',F7.4,/,
     +        ' Ds     = ',F7.4,/,
     +        ' betas  = ',F7.4,/,
     +        ' Elast  = ',F7.4,/,
     +        ' rho1   = ',F7.4,/,
     +        ' rho2   = ',F7.4,/,
     +        ' g      = ',F7.4,/,
     +        ' lamda  = ',F15.5,/,/,
     +        ' eps    = ',F15.13,/,
     +        ' Niter  = ',I3,/,
     +        ' thet0  = ',F7.4,/,
     +        ' Norm   = ',I1,/,
     +        ' Iadjust= ',I1,/,
     +        ' Eadj   = ',F7.4,/
     +        ' Roexp  = ',F7.4,/
     +        ' Mexp   = ',I3,
     +        ' Mexp_foam  = ',I3
     +        )
 106  Format  (/,
     +        ' h     = ',F7.4,/,
     +        ' Ngfww = ',I2,/,
     +        ' IQPD  = ',I2,/
     +        )
 107  Format (' THMAX = ',F6.3,/,
     +        ' SPMAX = ',F6.3,/,
     +        ' SPMIN = ',F6.3)
 108  Format  (
     +        ' IS_slp= ',I1,/,
     +        ' IS_dlp= ',I1,/,
     +        ' NGL   = ',I2,/,
     +        ' Dt    = ',F8.6,/,
     +        ' RUNGE = ',I1,/,
     +        ' MOVE  = ',I1,/,
     +        ' Defl  = ',I1,/,
     +        ' JGS   = ',I1,/
     +        )
 109  Format (' time=',F7.3,' area=',F10.7,' arc_l=',F10.7,
     +        ' Xc=',F8.5,' Yc=',F8.5)

 110  Format ( " PROGRAM CAPS_2D",/,/)
 111  Format ( " Executing step ",I4, " Total steps: ",I4)
 112  Format (1X,I4,1X,I4)
 113  Format (1X,F8.5,1X,I3)
 115  Format (10(1X,f10.5))

 140  Format (1X,I4,1x,f10.6,20(1X,F10.5))
 141  Format (1X,I4,1x,f10.6,20(1X,F11.5))
c145  Format (1X,I4,1x,f15.6,2(1X,F10.5),1X,f15.5)
 145  Format (1X,I4,1x,f10.7,10(1X,F10.7),1X,f15.5)

 305  Format (' ARC LENGTH   = ',F12.8)
 306  Format (' SURFACE AREA = ',F12.8)
 308  Format (' CENTER       = ',F12.8,'  ',F12.8)

      stop
      end
