      program em_2d 

c=============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=============================================

c--------------------------------------------------
c Dynamic simulation of the motion of an emulsion (suspension)
c of two-dimensional deformable liquid drops, bubbles, or capsules
c for a variety of configurations.
c
c The interfaces may be occupied by an insoluble surfactant
c
c The interfaces exhibit a combination of
c (a) isotropic surface tensions
c (b) elastic tensions
c-------------------------------------------
c
c  Default dimension capacity:
c
c     49 interfaces
c    64 points along each interface
c    2000 time steps
c---------
c
c SYMBOLS:
c -------
c
c NSG  :  Number of segments (nodes) along an interface
c x,y  :  coordinates of nodes
c u,v  :  velocity components
c
c tinit :  initial constant surface tension
c cinit :  initial constant surfactant concentration
c srfin :  initial amount of surfactant
c elst  :  modulus of elasticity
c
c srtn :  surface tension at nodes
c
c suns :  arc length around unstressed shape
c s    :  arc length around the interface
c
c c    :  surfactant concentration at nodes
c cm   :  surfactant concentration at middle nodes
c Dcds    d(concentration)/d(arc length)
c
c vnx,vny : normal vector at nodes
c vtx,vty : tangential vector at nodes
c
c c	: surfactant concentration at end-nodes
c cm	: surfactant concentration at mid-nodes
c crv	: curvature at end-nodes
c
c sxy	: effective shear stress
c sd1	: effective first normal-stress difference
c efv	: effective viscosity for Poiseuiile flow
c
c
c---
c Variables ending with g are global
c---
c
c For example:
c
c  srft(l):     surface tension of the lth interface
c
c  Xg(i,l):  x-position of ith point on lth interface
c  Yg(i,l):  y-position of ith point on lth interface
c  Ug(i,l):  u-velocity of ith point on lth interface
c  Vg(i,l):  v-velocity of ith point on lth interface
c
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c global variables
c---

      Dimension    xg(0:066,49),   yg(0:066,49)
      Dimension    ug(0:066,49),   vg(0:066,49)
      Dimension srtng(0:066,49),sunsg(0:066,49)
      Dimension    cg(0:066,49)
      Dimension               unusedg(0:066,49)

      Dimension  XCg(066,49), YCg(066,49),  Rg(066,49),   sg(066,49)
      Dimension TH1g(066,49),TH2g(066,49),TH3g(066,49),orntg(066,49)

      Dimension fxg(0:066,49),fyg(0:066,49)

c---
c properties of drops and interfaces
c---

      Dimension   NSG(49),NSG1(49),NSG2(49)
      Dimension  rhod(49), vsd(49),Ivsg(49)
      Dimension  elst(49),Ielstg(49)
      Dimension cinit(49), tinit(49),Ds(49),betas(49)
      Dimension srfin(49)

      Dimension xdc(49),ydc(49)
      Dimension aob(49),req(49),awob(49),bwob(49)
      Dimension peri0(49),epif0(49)

c---
c local variables
c---

      Dimension X(0:200),Y(0:200)
      Dimension U(0:200),V(0:200)
      Dimension suns(0:200),Unused(0:200)

      Dimension fx(0:200),fy(0:200)

      Dimension  XC(200), YC(200),  R(200),   S(200)
      Dimension TH1(200),TH2(200),TH3(200),ORNT(200)

      Dimension   c(0:200), cm(0:200),dcds(0:200)
      Dimension vnx(0:200),vny(0:200)
      Dimension vtx(0:200),vty(0:200)
      Dimension crv(0:200)

c---
c  save for the second-order Runge-Kutta (RK2)
c---

      Dimension Xsv(066,49),Ysv(066,49)
      Dimension Usv(066,49),Vsv(066,49)

c---
c  related to time stepping
c---

      Dimension time(0:2000)
      Dimension  sxy(0:2000),sd1(0:2000)
      Dimension  efv(0:2000)

      Dimension peri(0:2000,49),epif (0:2000,49)
      Dimension xcen(0:2000,49),ycen (0:2000,49)
      Dimension Dxy (0:2000,49),anmax(0:2000,49),anmin(0:2000,49)
 
c---
c  various
c---

      Dimension ZZ(20),WW(20)      ! Gauss-Legendre weights

c---
c for sgf_2d_2p
c---

      Dimension glut_xx(-64:64,0:64,-64:64)
      Dimension glut_xy(-64:64,0:64,-64:64)
      Dimension glut_yy(-64:64,0:64,-64:64)

c--------------
c common blocks
c--------------

c---
c common blocks for global variables
c---

      common/XXYYg/Xg,Yg
      common/ARCCg/XCg,YCg,Rg,Sg,TH1g,TH2g,TH3g,ORNTg
      common/UUVVg/Ug,Vg
      common/TENSIONg/srtng
      common/ELASTg/sunsg
      common/DDFFg/fxg,fyg

c---
c common blocks for local variables
c---

      common/XXYY/X,Y
      common/ARCC/XC,YC,R,S,TH1,TH2,TH3,ORNT
      common/UUVV/U,V
      common/DDFF/fx,fy

c---
c various common blocks
c---

      common/ancR1/gx,gy
      common/ancR2/vs1,vsd
      common/ancR4/rho1,rhod,elst
      common/ancR5/RL
      common/ancR6/shrt,U1,U2,pg

      common/hhhh/h,hh,h2,h3,h4,hs
      common/teihos/wall

      common/ancI1/Ndrops,NSG,NSG1,NSG2
      common/ancI2/NGFww,IQPD
      common/ancI3/NGL
      common/ancI4/Ivs,Ivsg
      common/ancI7/Ielst,Ielstg

c---
c doubly periodic flow
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
 
c----------
c constants
c----------

      pi   = 3.14159 265358 979 32384 D0
      pih  = 0.50D0*pi 
      piq  = 0.25D0*pi
      pi2  = 2.00D0*pi
      pi4  = 4.00D0*pi
      pi6  = 6.00D0*pi
      pi8  = 8.00D0*pi
      srpi = sqrt(pi)
 
      Null = 0
      None = 1
      Ntwo = 2

      zero = 0.0D0

c-----------
c parameters
c-----------

      write (6,*) 
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 1 to read data from input files'
      write (6,*) ' 2 to type in parameters'
      write (6,*) ' 0 to quit'
      write (6,*)
      write (6,*) ' Warning: if you type 2 will spend lots'
      write (6,*) ' of time typing numbers'
      write (6,*) ' ----------------------'
      read (5,*)   Ienrd

c------------------------
      if(Ienrd.eq.1) then
c------------------------

        open (3,file='em_2d.dat',status='old')

        read (3,*) Iflow
        read (3,*)
        read (3,*) rho1
        read (3,*) vs1
        read (3,*) 
        read (3,*) gac
        read (3,*) th0
        read (3,*)
        read (3,*) Isolve
        read (3,*)
        read (3,*) eps
        read (3,*) Nter
        read (3,*) Idfl
        read (3,*) JGS
        read (3,*)
        read (3,*) shrt
        read (3,*) U1,U2
        read (3,*) pg
        read (3,*)
        read (3,*) wall
        read (3,*)
        read (3,*) h
        read (3,*) NGFww
        read (3,*) IQPD
        read (3,*)
        read (3,*) IS_slp,IS_dlp
        read (3,*) NGL
        read (3,*)  
        read (3,*) THmax
        read (3,*) SPmin
        read (3,*) SPmax
        read (3,*) 
        read (3,*) a11,a12
        read (3,*) a21,a22
        read (3,*) Iglut
        read (3,*) Max1,Max2
        read (3,*)
        read (3,*) RL
        read (3,*)  
        read (3,*) Iread
        read (3,*)
        read (3,*) Ishape
        read (3,*)
        read (3,*) Norm
        read (3,*)  
        read (3,*) Nprint
        read (3,*)   
        read (3,*) IRK
        read (3,*) Dt
        read (3,*) move
        read (3,*) 

c------------------------------------------
c read initial capsule positions and shapes
c from file: em_2d_geo.dat
c
c If Index = 0, shapes are uniform
c------------------------------------------

c-----------------
        open (8,file='em_2d_geo.dat')
        read (8,*) Ndrops
c-----------------
        Do i=1,Ndrops
         read (8,*) idle,xdc(i),ydc(i)
        End Do 
c-----------------
        read (8,*) 
c-----------------
        read (8,*) Index,Nsgr   ! points around interfaces
        Do i=1,Ndrops
          Nsg(i) = Nsgr
          If(Index.eq.1) read (8,*) idle,Nsg(i)
        End Do
c-----------------
        read (8,*) Index,aobr     ! axes ratio
        Do i=1,Ndrops
          aob(i) = aobr
          If(Index.eq.1) read (8,*) idle,aob(i)
        End Do
c-----------------
        read (8,*) Index,reqr     ! equivalent radius
        Do i=1,Ndrops
          req(i) = reqr
          If(Index.eq.1) read (8,*) idle,req(i)
        End Do
c-----------------
        read (8,*) Index,awobr     ! first shape parameter 
        Do i=1,Ndrops
          awob(i) = awobr
          If(Index.eq.1) read (8,*) idle,awob(i)
        End Do
c-----------------
        read (8,*) Index,bwobr     ! second shape parameter
        Do i=1,Ndrops
          bwob(i) = bwobr
          If(Index.eq.1) read (8,*) idle,bwob(i)
        End Do
c-----------------
        close (8)
c-----------------

c-------------------------------------
c read drop properties from file: em_2d_prop.dat
c
c if Index = 0, properties are uniform
c-------------------------------------

        open (7,file='em_2d_prop.dat')

c---
        read (7,*) Index,vsdr      ! viscosity
        Do i=1,Ndrops
          vsd(i) = vsdr
          If(Index.eq.1) read (7,*) idle,vsd(i)
        End Do
c---
        read (7,*) Index,rhodr      ! density
        Do i=1,Ndrops
          rhod(i) = rhodr
          If(Index.eq.1) read (7,*) idle,rhod(i)
        End Do
c---
        read (7,*) Index,tinitr   ! initial surface tension
        Do i=1,Ndrops
          tinit(i) = tinitr
          If(Index.eq.1) read (7,*) idle,tinit(i)
        End Do
c---
        read (7,*) Index,cinitr   ! initial surf conc
        Do i=1,Ndrops
          cinit(i) = cinitr
          If(Index.eq.1) read (7,*) idle,cinit(i)
        End Do
c---
        read (7,*) Index,Dsr   ! surfactant diffusivity
        Do i=1,Ndrops
          Ds(i) = Dsr
          If(Index.eq.1) read (7,*) idle,Ds(i)
        End Do
c---
        read (7,*) Index,betasr   ! surfactant sensitivity beta
        Do i=1,Ndrops
          betas(i) = betasr
          If(Index.eq.1) read (7,*) idle,betas(i)
        End Do
c---
        read (7,*) Index,elstr   ! elasticity
        Do i =1,Ndrops
          elst(i) = elstr
          If(Index.eq.1) read (7,*) idle,elst(i)
        End Do
c---
        close(7)
c---------
      else
c---------

      call verbal
     +
     + (Iflow
     + ,cinit
     + ,tinit
     + ,Ds
     + ,betas
     + ,gac
     + ,th0
     + ,IS_slp,IS_dlp
     + ,thmax,SPmin,SPmax
     + ,Iread,Ishape
     + ,aob,req,awob,bwob
     + ,xdc,ydc
     + ,Isolve
     + ,eps,Nter,Idfl,JGS
     + ,Norm,Nprint,IRK,Dt,Move
     + )

c-----------
      end if     ! end of reading parameters
c-----------

c     write (6,105) 
c    +
c    +  Iflow
c    +  ,rho1,vs1
c    +  ,gac,th0
c    +  ,Isolve
c    +  ,eps,Nter,Idfl,JGS
c    +  ,shrt,U1,U2,pg
c    +  ,wall
c    +  ,h,NGFww,IQPD
c    +  ,IS_slp,IS_dlp
c    +  ,NGL
c    +  ,THmax,SPmax,SPmin
c    +  ,a11,a12
c    +  ,a21,a22
c    +  ,Max1,Max2
c    +  ,RL
c    +  ,Iread
c    +  ,Ishape
c    +  ,Norm
c    +  ,Nprint
c    +  ,IRK,DT,move

c-----------------
c printing session
c-----------------

      write (6,*) 
      write (6,*) Ndrops," capsules"
      write (6,*) 
      write (6,108) 
      write (6,*) 

      Do j=1,Ndrops
        write (6,107) j,xdc(j),ydc(j),req(j),aob(j),awob(j),bwob(j)
      End Do

      write (6,*) 
      write (6,106) 
      write (6,*) 

      Do j=1,Ndrops
        write (6,107) j,vsd(j),rhod(j),tinit(j),cinit(j)
     +                 ,Ds(j),betas(j),elst(j)
      End Do

c------------
c Flag for surfactant concentration
c Check-point
c------------

      Do j=1,Ndrops

        if(abs(betas(j)-1.0).lt.0.00000001) then
         write (6,*) 
         write (6,*) ' em_2d: beta for interface',j
         write (6,*) '        cannot be 1.'
         write (6,*) 
         stop
        else if(betas(j).gt.0.00000001) then
         Ibetas = 1
        else
         Ibetas = 0
         cinit(j) = 1.0D0
        end if

      End Do

c---------------------------------
c prepare for doubly periodic flow
c---------------------------------

      if(Iflow.eq.5) then

        a11h  = 0.5D0*a11
        a11hm = - a11h

        box_area = a11*a22-a12*a21     ! area of a periodic box
        fc_rhe = vs1*shrt*box_area     ! for rheology

c---
c read the look-up tables
c---

        if(Iglut.eq.1) then

        open (4,file="sgf_2d_2p.glut")

        write (6,*)
        write (6,*)  " Reading the look-up tables"
        write (6,*)  " This may take a while"
        write (6,*)  " Please twiddle your thumbs in anticipation"
        write (6,*)

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

c---------------------------------
c prepare for singly periodic flow
c---------------------------------

      if(Iflow.eq.11.or.Iflow.eq.20) then

        RLM  = -RL
        RL2  =  2.0*RL
        RL2M = -RL2

      end if

c-------------------------
c prepare for channel flow
c-------------------------

      if(Iflow.eq.20) then

        hh = 0.5*h
        h2 = h+h
        h3 = h2+h
        h4 = h2+h2
        hs = h**2
        hm = -h

      end if

c-----------
c quadrature
c-----------

      call gauss_leg (NGL,ZZ,WW)

c---------------------------------
c initial drop perimeter and area
c---------------------------------

      Do i=1,Ndrops

        if(Ishape.eq.1) then
         peri0(i) = pi2 * req(i)
         epif0(i) = pi  * req(i)**2
        else
         peri0(i) = pi2 
         epif0(i) = pi 
        end if

      End Do

c------------------------
c additional preparations
c------------------------

      thmax = thmax*pi

      gx =   gac*Dsin(th0*pi)
      gy = - gac*Dcos(th0*pi)

c----------------------
c viscosity ratio flags
c----------------------

      Ivs = 0        ! global index

      Do i=1,Ndrops

        Ivsg(i) = 0
        if(abs(vsd(i)-vs1).gt.0.000001) then
          Ivsg(i) = 1
          Ivs     = 1
        end if

      End Do

c-----------------
c elasticity flags
c-----------------

      Ielst = 0       ! global index

      Do i=1,Ndrops

        Ielstg(i) = 0
        if(elst(i).gt.0.000001) then
          Ielstg(i) = 1
          Ielst     = 1
        end if

      End Do

c------------------------
c Generate initial shapes
c
c Assign initial velocity
c surfactant concentration 
c------------------------

        time(1) = 0.0D0

        Do j=1,Ndrops      !--------- over drops

        write (6,*)
        write (6,*) " Generating interface:",j
        write (6,*)

c-------------------------
      if(Ishape.eq.1) then     ! ellipses
c-------------------------

        NSG1(j) = NSG(j)+1
        NSG2(j) = NSG(j)+2

        step = pi2/(NSG1(j)-1.0)
        rrx  = req(j) * sqrt(aob(j))
        rry  = rrx/aob(j)

        Do i=0,NSG2(j)
          tt = (i-1.0)*step
          cs = Dcos(tt)
          sn = Dsin(tt)
          xg(i,j) = xdc(j) + rrx * cs
          yg(i,j) = ydc(j) + rry * sn
          Ug(i,j) = 0.0
          Vg(i,j) = 0.0
          cg(i,j) = cinit(j)
          unusedg(i,j) = 1.0D0
        End Do

c------------------------------
      else if(Ishape.eq.2) then  ! wobbly shapes
c------------------------------

        NSG1(j) = NSG(j)+1
        NSG2(j) = NSG(j)+2

        step = pi2/(NSG1(j)-1.0)

        Do i=0,NSG2(j)
          TT = (i-1.0)*step
          TH = TT  + awob(j)*sin(8.0*TT)
          RR = 1.0 + bwob(j)*cos(4.0*TT)
          xg(i,j) = xdc(j) + RR*cos(TH)
          yg(i,j) = ydc(j) + RR*sin(TH)
          ug(i,j) = 0.0
          vg(i,j) = 0.0
          cg(i,j) = cinit(j)
          unusedg(i,j) = 1.0D0
        End Do

c-----------
      end if
c-----------

      End Do          !---------- over drops

c-------------------------------------------
c Compute arc length around the unstressed shape
c to be used for elastic tensions
c and for computing the total amount
c of a surfactant
c-------------------------------------------

        ICH1 = 0
        ICH2 = 0
        ICH3 = 0

        write (6,*)

        Do j=1,Ndrops    !--------------- over drops

        write (6,*) "Computing untressed arc length ",j

c---
c transfer into prd (point redistribution)
c---

        Do i=0,NSG2(j)
              x(i) =      xg(i,j)
              y(i) =      yg(i,j)
              u(i) =      ug(i,j)
              v(i) =      vg(i,j)
              c(i) =      cg(i,j)
         unused(i) = unusedg(i,j)
        End Do

        call prd_2d
     +
     +    (NSG(j)
     +    ,ICH1,THMAX
     +    ,ICH2,SPMAX
     +    ,ICH3,SPMIN
     +    ,dm1,dm2
     +    ,dm3,dm4
     +    ,U,V,c
     +    ,suns
     +    ,unused
     +    ,Istop
     +    )

        if(Istop.eq.1) Go to 99

        srfin(j) = s(NSG1(j))*cinit(j)  ! total amount of surfactant

        Do i=0,NSG2(j)
         sunsg(i,j) = s(i)
        End Do 

       End Do            ! over drops

c-----------------------------------
c read position of marker points
c      velocity
c      surfactant concentration
c      arc length around the unstressed shape
c
c      from file: em_2d.inp
c-----------------------------------

      if(Iread.eq.1) then

        open (2,file='em_2d.inp')

        Do j=1,Ndrops             !------ loop over drops

          write (6,*) "Reading Interface ",j

          if(Iflow.eq.5) then
           read (2,*) NSG1(j),idle,idle,time(1),a21
          else
           read (2,*) NSG1(j),idle,idle,time(1)
          end if

          if(Ielst.eq.1) then
            Do i=1,NSG1(j)
            read (2,*) idle,xg(i,j),yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),sunsg(i,j),cg(i,j)
            End Do
          else                      ! zero surface elasticity
            if(Ibetas.eq.0) then  ! surface tension is constant
             Do i=1,NSG1(j)
               read (2,*) idle,xg(i,j),yg(i,j),Ug(i,j),Vg(i,j)
               cg(i,j) = cinit(j)
             End Do
            else
             Do i=1,NSG1(j)
               read (2,*) idle,Xg(i,j),Yg(i,j),Ug(i,j),Vg(i,j),cg(i,j)
             End Do
            end if
          end if

          NSG (j) = NSG1(j)-1
          NSG2(j) = NSG (j)+2

        End Do                    ! end of loop over interfaces

        close (2)

c---
      end if   !------------ End of Iread
c--

c-----
c wrap
c-----

      Do j=1,Ndrops

           Xg(0,j) =    Xg(NSG(j),j)
           Yg(0,j) =    Yg(NSG(j),j)
           Ug(0,j) =    Ug(NSG(j),j)
           Vg(0,j) =    Vg(NSG(j),j)
           cg(0,j) =    cg(NSG(j),j)
        sunsg(0,j) = sunsg(NSG(j),j)-sunsg(NSG1(j),j)

        Xg(NSG1(j),j) = Xg(1,j)
        Yg(NSG1(j),j) = Yg(1,j)
        Ug(NSG1(j),j) = Ug(1,j)
        Vg(NSG1(j),j) = Vg(1,j)
        cg(NSG1(j),j) = cg(1,j)

           Xg(NSG2(j),j) =    Xg(2,j)
           Yg(NSG2(j),j) =    Yg(2,j)
           Ug(NSG2(j),j) =    Ug(2,j)
           Vg(NSG2(j),j) =    Vg(2,j)
           cg(NSG2(j),j) =    cg(2,j)
        sunsg(NSG2(j),j) = sunsg(2,j)+sunsg(NSG1(j),j)

      End Do

c----------------------
c Display initial shape
c on the screen
c----------------------

c      write (6,*)
c      write (6,*)  " Starting shapes"
c      write (6,*)

c      Do j = 1,Ndrops
c       write (6,*)
c       Do i = 0,NSG2(j)
c        If(Ielst.eq.1) then
c         write (6,100) j,i,Xg(i,j),Yg(i,j),cg(i,j),sunsg(i,j)
c        Else
c         write (6,100) j,i,xg(i,j),yg(i,j),cg(i,j)
c        End If
c       End Do
c      End Do

c---------------------
c Prepare output files
c---------------------

      open (1,file='em_2d.xy')      ! for profiles
      open (4,file='em_2d.diag')    ! for deformation and orientation
      open (7,file="em_2d.rhe")     ! for rheology

c---------------------------------
c plot walls or channel boundaries
c---------------------------------

      if(Iflow.eq.10.or.Iflow.eq.11) then

        write (1,103) Ntwo
        write (1,102) None,zero,wall
        write (1,102) Ntwo,RLM ,wall

      else if(Iflow.eq.20) then

        write (1,103) Ntwo
        write (1,102) None,RLM,hm
        write (1,102) Ntwo,RL ,hm
        write (1,103) Ntwo
        write (1,102) None,RLm,h
        write (1,102) Ntwo,RL ,h

      end if

c-.-.-.-.-.-.-.-.-.-.-
c BEGINNING OF THE SIMULATION
c-.-.-.-.-.-.-.-.-.-.-

      Kstep  = 1        ! Global time step counter
      Iprint = Nprint   ! Printing counter
      Itry = 1          ! graphics

  90  Continue

      if(Ienrd.eq.2) then
        write (6,*)
        write (6,*) 'Enter number of steps before pausing'
        write (6,*) '0 to quit'
        write (6,*) '---------'
        read  (5,*)  Nstep
      else
        read  (3,*)    Nstep
      end if

      if(Nstep.eq.0) Go to 99
      write (4,102) Nstep

c---
c batch session
c---

      Istep = 1     ! Batch step counter

  97  Continue

      write (6,*)
      write (6,*) "------------------"
      write (6,*) 
      write (6,111) Istep,Nstep

c-----------------------------
c Check the point distribution
c Compute properties of arcs
c-----------------------------

      Do j=1,Ndrops     !----- over drops

      write (6,*) "Checking interface :",j

c---
c transfer into local variables
c---

        Do i=0,NSG2(j)
              x(i) =      xg(i,j)
              y(i) =      yg(i,j)
              u(i) =      ug(i,j)
              v(i) =      vg(i,j)
              c(i) =      cg(i,j)
           suns(i) =   sunsg(i,j)
         unused(i) = unusedg(i,j)
        End Do

      ICH1 = 1
      ICH2 = 1
      ICH3 = 1

        call prd_2d
     +
     +   (NSG(j)
     +   ,ICH1,THMAX
     +   ,ICH2,SPMAX
     +   ,ICH3,SPMIN
     +   ,peri(Kstep,j),epif(Kstep,j)
     +   ,xcen(Kstep,j),ycen(Kstep,j)
     +   ,U,V,c
     +   ,suns,unused
     +   ,Istop
     +   )

        if(Istop.eq.1) Go to 99

        NSG1(j) = NSG(j)+1
        NSG2(j) = NSG(j)+2
   
c-------------------------
c normalize the drop areas
c-------------------------

      if(Norm.eq.1) then

        write (6,*) " em_2d: Normalizing area of drop : ",j

        epf = epif(Kstep,j)
        xcn = xcen(Kstep,j)
        ycn = ycen(Kstep,j)

        scale = 1.0D0/sqrt(epf/epif0(j))

        Do i=0,NSG2(j)
         X(i) = xcn+(x(i)-xcn)*scale
         Y(i) = ycn+(y(i)-ycn)*scale
        End Do

        ICH1 = 0  ! no redistribution
        ICH2 = 0
        ICH3 = 0

        call prd_2d
     +
     +    (NSG(j)
     +    ,ICH1,THMAX
     +    ,ICH2,SPMAX
     +    ,ICH3,SPMIN
     +    ,peri(Kstep,j),epif(Kstep,j)
     +    ,xcen(Kstep,j),ycen(Kstep,j)
     +    ,U,V,c
     +    ,suns,Unused
     +    ,Istop
     +    )

        if(Istop.eq.1) Go to 99

      end if

c----
c transfer into global variables
c---

      Do i=1,NSG1(j)
            xg(i,j) =      x(i)
            yg(i,j) =      y(i)
            ug(i,j) =      u(i)
            vg(i,j) =      v(i)
            cg(i,j) =      c(i)
         sunsg(i,j) =   suns(i)
       unusedg(i,j) = unused(i)
           xcg(i,j) =     xc(i)
           ycg(i,j) =     yc(i)
            rg(i,j) =      r(i)
            sg(i,j) =      s(i)
          th1g(i,j) =    th1(i)
          th2g(i,j) =    th2(i)
          th3g(i,j) =    th3(i)
         orntg(i,j) =   ornt(i)
      End Do

c---
c extend
c---

         xg(0,j) =    xg(NSG(j),j)
         yg(0,j) =    yg(NSG(j),j)
         ug(0,j) =    ug(NSG(j),j)
         vg(0,j) =    vg(NSG(j),j)
         cg(0,j) =    cg(NSG(j),j)
      sunsg(0,j) = sunsg(NSG(j),j)-sunsg(NSG1(j),j)

      it = NSG1(j)
      xg(it,j) = xg(1,j)
      yg(it,j) = yg(1,j)
      ug(it,j) = ug(1,j)
      vg(it,j) = vg(1,j)
      cg(it,j) = cg(1,j)

      it = NSG2(j)
      xg(it,j) = xg(2,j)
      yg(it,j) = yg(2,j)
      ug(it,j) = ug(2,j)
      vg(it,j) = vg(2,j)
      cg(it,j) = cg(2,j)
      sg(it,j) = sg(2,j)+sg(NSG1(j),j)
      sunsg(it,j) = sunsg(2,j)+sunsg(NSG1(j),j)

c---------------------------------------
c scale drop arc length, area, centroid
c---------------------------------------

      peri(Kstep,j) = peri(Kstep,j)/peri0(j)
      epif(Kstep,j) = epif(Kstep,j)/epif0(j)
  
c-------------------------------------
c compute the Taylor deformation parameter
c------------------------------------

        call taylor
     +
     +    (NSG(j)
     +    , xcen(Kstep,j)
     +    , ycen(Kstep,j)
     +    ,  Dxy(Kstep,j)
     +    ,anmax(Kstep,j)
     +    ,anmin(Kstep,j)
     +    ) 

      End Do                !------ loop over drops

c-----------------
c Printing session
c-----------------

      write (6,*)
      write (6,113) time(Kstep),Kstep
      write (6,*)

      write (6,305) (peri(Kstep,j),j=1,Ndrops)
      write (6,306) (epif(Kstep,j),j=1,Ndrops)
      write (6,307) (xcen(Kstep,j),j=1,Ndrops)
      write (6,308) (ycen(Kstep,j),j=1,Ndrops)
  
      write (4,113) time(Kstep),Kstep

      Do j=1,Ndrops
       write (4,104) j,xcen(Kstep,j), ycen(Kstep,j)
     +                ,Dxy(Kstep,j),anmax(Kstep,j),anmin(Kstep,j)
     +                ,peri(Kstep,j), epif(Kstep,j)
      End Do

      write (6,*)
      write (6,*) "xce-yce-Dxy-thmax-thmin-per-ar"
      write (6,*)

      Do j=1,Ndrops
       write (6,104) j,xcen(Kstep,j),ycen(Kstep,j)
     +                 ,Dxy(Kstep,j),anmax(Kstep,j),anmin(Kstep,j)
     +                ,peri(Kstep,j), epif(Kstep,j)
      End Do

c---------------------------------
c Prepare for the
c doubly-periodic Green's function
c---------------------------------

      if(Iflow.eq.5) then

        if(a21.gt.a11h)  a21 = a21-a11
        if(a21.lt.a11hm) a21 = a21+a11

        if(Iglut.eq.0) then

         call sgf_2d_2p_ewald
     +
     +    (a11,a12,a21,a22
     +    ,b11,b12,b21,b22
     +    ,ew,tau
     +    )

         call sgf_2d_2p_qqq (b11,b12,b21,b22,Max2,ew)

         if(Ivs.eq.1)
     +   call sgf_2d_2p_vvv (b11,b12,b21,b22,Max2,ew)

        end if

      end if

c-----------------------------------                   ! SRF
c Compute the surface tension from the concentration   
c at the end-nodes, using a constitutive equation
c
c Note that when betas(j) = 0: 
c srtn(i,j) = tinit(j)
c-----------------------------------   

      Do j=1,Ndrops
       Do i=0,NSG2(j)
       srtng(i,j) = tinit(j)*(1.0-betas(j)*cg(i,j)/cinit(j))
     +                      /(1.0-betas(j))
       End Do
      End Do

c----------------------
c Display shape on the screen
c----------------------

c      Do j=1,Ndrops
c       write (6,*)
c       Do i=0,NSG2(j)
c        If(Ielst.eq.1) then
c         write (6,100) j,i,Xg(i,j),Yg(i,j),cg(i,j)
c    +                     ,srtng(i,j),sunsg(i,j)
c        Else
c         write (6,100) j,i,xg(i,j),yg(i,j),cg(i,j),srtng(i,j)
c        End If
c       End Do
c      End Do

c-----------------
c compute velocity
c-----------------

      call em_2d_vel
     +
     +  (IS_slp,IS_dlp
     +  ,Isolve,JGS
     +  ,eps,Nter
     +  ,Iflow,Idfl
     +  ,Istop
     +  )

      write (6,*) ' em_2d: integral equation solved'

      if(Istop.eq.1) Go to 99

c-----------------------------------------
c RHEOLOGY
c
c Effective stress tensor will be computed
c only when Iswitch = 1 or Iflow = 5
c
c Its physical relevance is poorly
c founded otherwise
c-----------------------------------------

      Iswitch = 0

      if(Iflow.eq.20.and.Ivs.eq.0) Iswitch = 1

      Iskip = 0

      if(Iswitch.eq.1) Iskip = 1
      if(Iflow.eq.5)   Iskip = 1
      if(Iflow.eq.20)  Iskip = 1

      if(Iskip.eq.0) Go to 49

c-----------
c initialize
c-----------

      if(Iflow.eq.5) then
       strxx = 0.0
       strxy = 0.0
       stryy = 0.0
      end if

      if(Iswitch.eq.1) then
       prgr = 0.0
       flrt = 0.0
      end if

c------------------
      Do j=1,Ndrops
c------------------

c---
c transfer into local variables
c---

       Do i=1,NSG2(j)
           x(i) =    xg(i,j)
           y(i) =    yg(i,j)
           u(i) =    ug(i,j)
           v(i) =    vg(i,j)
          fx(i) =   fxg(i,j)
          fy(i) =   fyg(i,j)
          xc(i) =   xcg(i,j)
          yc(i) =   ycg(i,j)
           r(i) =    rg(i,j)
         th1(i) =  th1g(i,j)
         th2(i) =  th2g(i,j)
         th3(i) =  th3g(i,j)
        ornt(i) = orntg(i,j)
       End Do

c---------
c compute the particle stress tensor
c---------

        if(Iflow.eq.5.or.Iflow.eq.20) then

         call drop_2d_rhe 
     +
     +      (NSG(j)
     +      ,vs1,vsd(j)
     +      ,strxxd,strxyd,stryyd
     +      )

         strxx = strxx + strxxd
         strxy = strxy + strxyd
         stryy = stryy + stryyd

        end if

c------
        if(Iswitch.eq.1) then

c------------------------------
c Pressure-driven flow in a channel
c
c Compute disturbance pressure drop for
c vanishing disturbance flow rate (IQPD = 1)
c or the disturbance flow rate for
c vanishing pressure drop (IQPD = 0)
c------------------------------

         call em_2d_pdfl (NSG(j),IQPD,RL,prgrd,flrtd)

         prgr = prgr + prgrd
         flrt = flrt + flrtd

        end if

c-------------
        End Do      ! end of loop over drops
c-------------

      if(Iflow.eq.5) then

        sxy(Kstep) =  strxy       /fc_rhe
        sd1(Kstep) = (strxx-stryy)/fc_rhe
        write (6,140) Kstep,time(Kstep),sxy(Kstep),sd1(Kstep)
        write (7,140) Kstep,time(Kstep),sxy(Kstep),sd1(Kstep)

      else if(Iflow.eq.20) then

        sxy(Kstep) =  strxy 
        sd1(Kstep) =  strxx-stryy
        write (6,140) Kstep,time(Kstep),sxy(Kstep),sd1(Kstep)
        write (7,140) Kstep,time(Kstep),sxy(Kstep),sd1(Kstep)

      end if

      if(Iswitch.eq.1) then

        if(Iflow.eq.20) then
          flrt = 2.0/3.0*pg*h**3/vs1 + flrt/vs1
          prgr = pg                  + prgr
        end if

        efv(Kstep)=2.0/3.0*h**3*prgr/flrt 

        write (6,140) Kstep,time(Kstep),efv(Kstep)
        write (7,140) Kstep,time(Kstep),efv(Kstep)

      end if

c--------------------------------------- End of rheology module

  49  Continue

c-----------------
c printing session
c-----------------

      if(Iprint.eq.Nprint) then

       Do j=1,Ndrops

        if(Iflow.eq.5) then
         write (1,103) NSG1(j),j,Ndrops,time(Kstep),a21
         write (6,103) NSG1(j),j,Ndrops,time(Kstep),a21
        else
         write (1,103) NSG1(j),j,Ndrops,time(Kstep),Dt
         write (6,103) NSG1(j),j,Ndrops,time(Kstep),Dt
        end if

        if(Ielst.eq.1) then

          Do i=1,NSG1(j)
            write (1,120) i,Xg(i,j),Yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),sunsg(i,j),cg(i,j)
            write (6,120) i,Xg(i,j),Yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),sunsg(i,j),cg(i,j)
          End Do

        else

         if(Ibetas.eq.0) then
          Do i=1,NSG1(j)
            write (1,120) i,Xg(i,j),Yg(i,j),Ug(i,j),Vg(i,j)
            write (6,120) i,Xg(i,j),Yg(i,j),Ug(i,j),Vg(i,j)
          End Do
         else
          Do i=1,NSG1(j)
            write (1,120) i,Xg(i,j),Yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),cg(i,j)
            write (6,120) i,Xg(i,j),Yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),cg(i,j)
          End Do
         end if

        End If
       End Do

       Iprint = 0

      end if

c--------------------------------------------------
c The following block integrates the surfactant
c convection-diffusion equation for the surfactant
c concentration using a finite-volume method
c
c Block is bypassed if betas = 0,
c for then the surface tension remains constant
c---------------------------------------------------

      if(Ibetas.eq.0) Go to 88

      Do j=1,Ndrops                !----- start of loop over drops

      if(betas(j).lt.0.0000001) Go to 33

c---
c transfer into local variables
c---

        Do i=1,NSG2(j)
            x(i) =    xg(i,j)
            y(i) =    yg(i,j)
            u(i) =    ug(i,j)
            v(i) =    vg(i,j)
            c(i) =    cg(i,j)
           xc(i) =   xcg(i,j)
           yc(i) =   ycg(i,j)
            r(i) =    rg(i,j)
            s(i) =    sg(i,j)
          th1(i) =  th1g(i,j)
          th2(i) =  th2g(i,j)
          th3(i) =  th3g(i,j)
         ornt(i) = orntg(i,j)
        End Do

c------------------------------------------
c At the first step,
c interpolate concentration at middle-nodes
c from end-nodes
c------------------------------------------

      if(Kstep.eq.1) call interp_mn (NSG(j),c,s,cm)

c-------------------------------
c compute normal vector,
c         tangential vectors
c         curvature at end-nodes
c-------------------------------

      Do i=1,NSG(j)
        vnx(i) = cos(th2(i)) * ornt(i)
        vny(i) = sin(th2(i)) * ornt(i)
        vtx(i) = - vny(i)
        vty(i) =   vnx(i)
        crv(i) = ornt(i)/r(i)
      End Do

c---
c extend
c---

      vnx(0) = vnx(NSG(j))
      vny(0) = vny(NSG(j))
      vtx(0) = vtx(NSG(j))
      vty(0) = vty(NSG(j))
      crv(0) = crv(NSG(j))

      vnx(NSG1(j)) = vnx(1)
      vny(NSG1(j)) = vny(1)
      vtx(NSG1(j)) = vtx(1)
      vty(NSG1(j)) = vty(1)
      crv(NSG1(j)) = crv(1)

      vnx(NSG2(j)) = vnx(2)
      vny(NSG2(j)) = vny(2)
      vtx(NSG2(j)) = vtx(2)
      vty(NSG2(j)) = vty(2)
      crv(NSG2(j)) = crv(2)

c------------------------------------------
c  compute the derivative d(c)/d(l)
c  at end-nodes, by parabolic interpolation
c  with respect to arc length
c------------------------------------------

      Do i=2,NSG1(j)
       x0      = s(i-1)-s(i)
       x2      = s(i+1)-s(i)
       y0      = c(i-1)-c(i)
       y2      = c(i+1)-c(i)
       dcds(i) = (x0/x2*y2-x2/x0*y0)/(x0-x2)
      End Do

c---
c extend
c---

      dcds(0)       = dcds(NSG(j))
      dcds(1)       = dcds(NSG1(j))
      dcds(NSG2(j)) = dcds(2)

c---
c update mean surfactant concentration
c at mid-nodes using
c an implicit finite-volume method
c---

      call em_2d_cd
     +
     +  (j
     +  ,Move
     +  ,NSG(j)
     +  ,vnx,vny
     +  ,vtx,vty
     +  ,cm,c,crv
     +  ,Ds(j)
     +  ,Dt
     +  ,srfam    ! total amount of surfactant
     +  ,Istop
     +  )

      if(Istop.eq.1) Go to 99

c-----------------------------
c normalize the concentration
c of a surfactant
c-----------------------------

      fc = srfin(j)/srfam

      Do i=0,NSG2(j)
        c(i) = fc* c(i)
       cm(i) = fc*cm(i)
      End Do

c---------------------------------------
c Interpolate concentration at end-nodes
c from mid-nodes
c
c End-node concentration is used to compute
c the velocity
c---------------------------------------

      call interp_en (NSG(j),cm,s,c)

c---------------------------------------
c transfer concentration to global array
c---------------------------------------

      Do i=0,NSG2(j)
       cg(i,j) = c(i)
      End Do

   33 Continue  

      End Do                ! end of loop over drops

c------------------------------------------   ! End of surfactant module

 88   Continue

c------------------------
c Runge-Kutta integration
c------------------------

c----------------------
      if(IRK.eq.1) then    ! first order
c----------------------

        if(Move.eq.0) then

          Do j=1,Ndrops
            Do i=1,NSG1(j)
              xg(i,j) = xg(i,j) + Dt*Ug(i,j)
              yg(i,j) = yg(i,j) + Dt*Vg(i,j)
            End Do
          End Do

        else

          Do i=1,NSG1(j)
            cs   = cos(th2g(i,j))
            sn   = sin(th2g(i,j))
            Utot = Ug(i,j)*cs+Vg(i,j)*sn
            Umv  = Utot*cs
            Vmv  = Utot*sn
            xg(i,j) = xg(i,j) + Dt*Umv
            yg(i,j) = yg(i,j) + Dt*Vmv
          End Do

        End If

      if(Iflow.eq.5) a21 = a21+Dt*shrt*a22

c---------------------------
      else if(IRK.eq.2) then    ! second order
c---------------------------

      if(Move.eq.0) then

        Do j=1,Ndrops
          Do i=1,NSG1(j)
           Xsv(i,j) = Xg(i,j)  ! save
           Ysv(i,j) = Yg(i,j)
           Usv(i,j) = Ug(i,j)
           Vsv(i,j) = Vg(i,j)
            Xg(i,j) = Xg(i,j) + Dt*Ug(i,j)
            Yg(i,j) = Yg(i,j) + Dt*Vg(i,j)
          End Do
        End Do

       else

        Do j=1,Ndrops
          Do i=1,NSG1(j)
            Xsv(i,j) = Xg(i,j)  ! save
            Ysv(i,j) = Yg(i,j)
            cs       = cos(th2g(i,j))
            sn       = sin(th2g(i,j))
            Utot     = Ug(i,j)*cs+Vg(i,j)*sn
            Umv      = Utot*cs
            Vmv      = Utot*sn
            Usv(i,j) = Umv
            Vsv(i,j) = Vmv
             Xg(i,j) = Xg(i,j) + Dt*Umv
             Yg(i,j) = Yg(i,j) + Dt*Vmv
          End Do
        End Do

       end if

       if(Iflow.eq.5) a21 = a21+Dt*shrt*a22

c---
c Second step
c---

      ICH1 = 0  ! will not redistribute points
      ICH2 = 0
      ICH3 = 0

      Do j=1,Ndrops          ! loop over drops

c---
c wrap around
c---

         Xg(0,j) =    Xg(NSG(j),j)
         Yg(0,j) =    Yg(NSG(j),j)
         Ug(0,j) =    Ug(NSG(j),j)
         Vg(0,j) =    Vg(NSG(j),j)
         cg(0,j) =    cg(NSG(j),j)
      srtng(0,j) = srtng(NSG(j),j)
      sunsg(0,j) = sunsg(NSG1(j),j)-sunsg(NSG(j),j)

         Xg(NSG1(j),j) =    Xg(1,j)
         Yg(NSG1(j),j) =    Yg(1,j)
         Ug(NSG1(j),j) =    Ug(1,j)
         Vg(NSG1(j),j) =    Vg(1,j)
         cg(NSG1(j),j) =    cg(1,j)
      srtng(NSG1(j),j) = srtng(1,j)
      sunsg(NSG1(j),j) = sunsg(NSG1(j),j)+sunsg(1,j)

         Xg(NSG2(j),j) =    Xg(2,j)
         Yg(NSG2(j),j) =    Yg(2,j)
         Ug(NSG2(j),j) =    Ug(2,j)
         Vg(NSG2(j),j) =    Vg(2,j)
         cg(NSG2(j),j) =    cg(2,j)
      srtng(NSG2(j),j) = srtng(2,j)
      sunsg(NSG2(j),j) = sunsg(NSG1(j),j)+sunsg(2,j)

c---
c Transfer into prd
c---

      Do i=0,NSG2(j)
       x(i) = xg(i,j)
       y(i) = yg(i,j)
      End Do

c---------------
c Geometry at intermediate step
c--------------

      write (6,*) "Checking interface :",j

      call prd_2d
     +
     +   (NSG(j)
     +   ,ICH1,THMAX
     +   ,ICH2,SPMAX
     +   ,ICH3,SPMIN
     +   ,peri,epif,xcen,ycen
     +   ,U,V,c
     +   ,Suns,Unused
     +   ,Istop
     +   )

      if(Istop.eq.1) Go to 99

c---
c transfer from prd
c---

        Do i=1,NSG1(j)
            xg(i,j) =    x(i)
            yg(i,j) =    y(i)
           xcg(i,j) =   xc(i)
           ycg(i,j) =   yc(i)
            rg(i,j) =    r(i)
            sg(i,j) =    s(i)
          th1g(i,j) =  th1(i)
          th2g(i,j) =  th2(i)
          th3g(i,j) =  th3(i)
         orntg(i,j) = ornt(i)
      End Do

c---
      End Do                    ! End of loop over drops
c---

c------------------------------------------------ 
c prepare for the doubly-periodic Green's function
c------------------------------------------------ 

      if(Iflow.eq.5) then

        if(Iglut.eq.0) then

         call sgf_2d_2p_ewald
     +
     +    (a11,a12,a21,a22
     +    ,b11,b12,b21,b22
     +    ,ew,tau
     +    )

         call sgf_2d_2p_qqq (b11,b12,b21,b22,Max2,ew)

         if(Ivs.eq.1)
     +   call sgf_2d_2p_vvv (b11,b12,b21,b22,Max2,ew)

        end if

      end if
c---

c-----------------------------------              ! SRF
c Compute surface tension from concentration   
c at the end-nodes, using a constitutive equation
c
c Note that when betas(j) = 0,
c  srtn(i,j) = tinit(j)
c-----------------------------------   

      Do j=1,Ndrops

       Do i=0,NSG2(j)
       srtng(i,j) = tinit(j)*(1.0-betas(j)*cg(i,j)/cinit(j))
     +                      /(1.0-betas(j))
       End Do

      End Do

c------------------
c Second velocity evaluation
c------------------

      call em_2d_vel
     +
     +  (IS_slp,IS_dlp
     +  ,Isolve,JGS
     +  ,eps,Nter
     +  ,Iflow,Idfl
     +  ,Istop
     +  )

      if(Istop.eq.1) Go to 99

c---
c second step in RK2
c---

      Dth = 0.5*Dt

      if(Move.eq.0) then

        Do j=1,Ndrops
          Do i=1,NSG(j)
           Xg(i,j) = Xsv(i,j)+Dth*(Usv(i,j)+Ug(i,j))
           Yg(i,j) = Ysv(i,j)+Dth*(Vsv(i,j)+Vg(i,j))
          End Do
        End Do

      else

        Do j=1,Ndrops
          Do i=1,NSG(j)
            cs   = cos(th2g(i,j))
            sn   = sin(th2g(i,j))
            Utot = Ug(i,j)*cs+Vg(i,j)*sn
            Umv  = Utot*cs
            Vmv  = Utot*sn
            Xg(i,j) = Xsv(i,j)+Dth*(Usv(i,j)+Umv)
            Yg(i,j) = Ysv(i,j)+Dth*(Vsv(i,j)+Vmv)
          End Do
        End Do

      end if

c---
      End If                        ! End of RK2
c---

c------------------------------------------------

c----------
c End of a time step
c----------

  96  Continue

c---
c wrap around
c---

      Do j=1,Ndrops

         Xg(0,j) =    Xg(NSG(j),j)
         Yg(0,j) =    Yg(NSG(j),j)
         Ug(0,j) =    Ug(NSG(j),j)
         Vg(0,j) =    Vg(NSG(j),j)
         cg(0,j) =    cg(NSG(j),j)
      srtng(0,j) = srtng(NSG(j),j)
      sunsg(0,j) = sunsg(NSG1(j),j)-sunsg(NSG(j),j)

         Xg(NSG1(j),j) =    Xg(1,j)
         Yg(NSG1(j),j) =    Yg(1,j)
         Ug(NSG1(j),j) =    Ug(1,j)
         Vg(NSG1(j),j) =    Vg(1,j)
         cg(NSG1(j),j) =    cg(1,j)
      srtng(NSG1(j),j) = srtng(1,j)
      sunsg(NSG1(j),j) = sunsg(NSG1(j),j)+sunsg(1,j)

      Xg   (NSG2(j),j) =    Xg(2,j)
      Yg   (NSG2(j),j) =    Yg(2,j)
      Ug   (NSG2(j),j) =    Ug(2,j)
      Vg   (NSG2(j),j) =    Vg(2,j)
      cg   (NSG2(j),j) =    cg(2,j)
      srtng(NSG2(j),j) = srtng(2,j)
      sunsg(NSG2(j),j) = sunsg(NSG1(j),j)+sunsg(2,j)

      End Do

c------------------------
c Reset position of drops
c
c Doubly-periodic flow
c------------------------

      if(Iflow.eq.5) then         ! doubly-periodic array

      Do j=1,Ndrops

       if(ycen(Kstep,j).gt.a22) then  ! move along second base vector
        Do i=0,NSG2(j)
         xg(i,j) = xg(i,j)-a21
         yg(i,j) = yg(i,j)-a22
        End Do
        xcen(Kstep,j) = xcen(Kstep,j)-a21
        ycen(Kstep,j) = ycen(Kstep,j)-a22
       end if

       if(ycen(Kstep,j).lt.0) then  ! move along second base vector
        Do i=0,NSG2(j)
         xg(i,j) = xg(i,j)+a21
         yg(i,j) = yg(i,j)+a22
        End Do
        xcen(Kstep,j) = xcen(Kstep,j)+a21
        ycen(Kstep,j) = ycen(Kstep,j)+a22
       end if

       if(xcen(Kstep,j).gt.a11) then  ! move along first base vector
        Do i=0,NSG2(j)
         xg(i,j) = xg(i,j)-a11
         yg(i,j) = yg(i,j)-a12
        End Do
        xcen(Kstep,j) = xcen(Kstep,j)-a11
        ycen(Kstep,j) = ycen(Kstep,j)-a12
       end if

       if(xcen(Kstep,j).lt.0) then ! move along first base vector
        Do i=0,NSG2(j)
         xg(i,j) = xg(i,j)+a11
         yg(i,j) = yg(i,j)+a12
        End Do
        xcen(Kstep,j) = xcen(Kstep,j)+a11
        ycen(Kstep,j) = ycen(Kstep,j)+a12
       end if

      End Do

      end if

c------------------------
c Reset position of drops
c
c Singly-periodic flow
c------------------------

      if(Iflow.eq.11.or.Iflow.eq.20) then

      Do j=1,Ndrops

       if(xcen(Kstep,j).gt.RL) then
        Do i=0,NSG2(j)
         xg(i,j) = xg(i,j)-RL
        End Do
       end if

       if(xcen(Kstep,j).lt.0) then
        Do i=0,NSG2(j)
         xg(i,j) = xg(i,j)+RL
        End Do
       end if

      End Do

      end if

c------------------------
c reset counters and time
c------------------------

      Istep  = Istep+1
      Kstep  = Kstep+1
      Iprint = Iprint+1

      time(Kstep) = time(Kstep-1) + Dt

      if(Istep.le.Nstep) Go to 97

      Go to 90

c-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
c---------------------------------c
c  Simulation has ended           c
c---------------------------------c
c-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

  99  Continue

c----------------------
c finish up by printing
c----------------------

       open (9,file="em_2d.rst")     ! restart file

       Do j=1,Ndrops

        if(Iflow.eq.5) then
         write (1,103) NSG1(j),j,Ndrops,time(Kstep),a21
         write (9,103) NSG1(j),j,Ndrops,time(Kstep),a21
        else
         write (1,103) NSG1(j),j,Ndrops,time(Kstep)
         write (9,103) NSG1(j),j,Ndrops,time(Kstep)
        end if

        if(Ielst.eq.1) then
          Do i=1,NSG1(j)
            write (1,120) i,Xg(i,j),Yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),sunsg(i,j),cg(i,j)
            write (9,120) i,Xg(i,j),Yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),sunsg(i,j),cg(i,j)
          End Do
        else
          If(Ibetas.eq.0) then
           Do i=1,NSG1(j)
             write (1,120) i,Xg(i,j),Yg(i,j)
     +                      ,Ug(i,j),Vg(i,j),srtng(i,j)
             write (9,120) i,Xg(i,j),Yg(i,j)
     +                      ,Ug(i,j),Vg(i,j),srtng(i,j)
     +                      
           End Do
          else
           Do i=1,NSG1(j)
            write (1,120) i,Xg(i,j),Yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),cg(i,j)
            write (9,120) i,Xg(i,j),Yg(i,j)
     +                     ,Ug(i,j),Vg(i,j),cg(i,j)
           End Do
          End If
        end if

       End Do

      write (1,112) Null,Null,Null,Null,Null,Null,Null
      write (9,112) Null,Null,Null,Null,Null,Null,Null

      close (9)

c---
c Print run parameters
c---

      thmax = thmax/pi

      write (1,105) Iflow,Ndrops
     +              ,rho1,vs1
     +              ,gac,th0
     +              ,Isolve
     +              ,eps,Nter,Idfl,JGS
     +              ,shrt,U1,U2,pg
     +              ,wall
     +              ,h,NGFww,IQPD
     +              ,IS_slp,IS_dlp
     +              ,NGL
     +              ,THmax,SPmax,SPmin
     +              ,a11,a12
     +              ,a21,a22
     +              ,Max1,Max2
     +              ,RL
     +              ,Iread
     +              ,Ishape
     +              ,Norm
     +              ,Nprint
     +              ,IRK,DT,move

c---
c print parameters
c---

      write (1,*) 
      write (1,106) 
      write (1,*) 

      Do j=1,Ndrops
       write (1,107) j,vsd(j),rhod(j),tinit(j),cinit(j),Ds(j)
     +               ,betas(j),elst(j)
       write (4,107) j,vsd(j),rhod(j),tinit(j),cinit(j),Ds(j)
     +               ,betas(j),elst(j)
      End Do

      write (1,*) 
      write (1,108) 
      write (1,*) 

      Do j=1,Ndrops
       write (1,107) j,xdc(j),ydc(j),req(j),awob(j),bwob(j)
       write (4,107) j,xdc(j),ydc(j),req(j),awob(j),bwob(j)
      End Do

c---
c close down
c---

      write (1,110) 
      write (4,112) Null,Null,Null,Null,Null,Null,Null

      close (1)
      if(Ienrd.eq.1) close(3)
      close (4)
      close (7)
 
c-----
c Done
c-----

 100  Format (1X,I2,1X,I2,10(1X,F10.5))
 101  Format (1X,'U = ',F15.10,'  V = ',F15.10)
 102  Format (1X,I4,5(1X,F12.8))
 120  Format (1X,I4,2(1X,F12.8),2(1X,F8.4),2(1X,F12.8))
 103  Format (1X,I2,1X,I2,1X,I2,3(1X,F15.5))
 104  Format (1X,I4,20(1X,F9.5))
 105  Format (/,
     +        ' Iflow  = ',I2,/,/,
     +        ' Ndrops = ',I2,/,/,
     +        ' rho1   = ',F7.4,/,
     +        ' vs1    = ',F7.4,/,
     +        ' g      = ',F7.4,/,
     +        ' th0    = ',F7.4,/,
     +        ' Isolve = ',I2,/,/,
     +        ' eps    = ',F15.13,/,
     +        ' Nter   = ',I3,/,
     +        ' Idfl   = ',I3,/,
     +        ' JGS    = ',I3,/,
     +        ' shrt   = ',F7.4,/,
     +        ' U1     = ',F7.4,/,
     +        ' U2     = ',F7.4,/,
     +        ' pg     = ',F7.4,/,/,
     +        ' wall   = ',F7.4,/,
     +        ' h      = ',F7.4,/,
     +        ' NGFww  = ',I3,/,
     +        ' IQPD   = ',I1,/,
     +        ' IS_slp = ',I1,/,
     +        ' IS_dlp = ',I1,/,
     +        ' NGL    = ',I1,/,
     +        ' THMAX  = ',F6.3,/,
     +        ' SPMAX  = ',F6.3,/,
     +        ' SPMIN  = ',F6.3,/,
     +        ' a11,a12= ',F6.3,1x,F6.3,/,
     +        ' a21,a22= ',F6.3,1x,F6.3,/,
     +        ' Max1   = ',I3,/,
     +        ' Max2   = ',I3,/,
     +        ' RL     = ',F7.4,/,
     +        ' Iread  = ',I1,/,
     +        ' Ishape = ',I1,/,
     +        ' Norm   = ',I1,/,
     +        ' Nprint = ',I3,/,
     +        ' IRK    = ',I1,/,
     +        ' Dt     = ',F7.4,/,
     +        ' Move   = ',I1,/
     +        )

 106  Format ('vsd, rhod, tinit, cinit, Ds, betas, elst')
 107  Format (1x,i2,10(1x,f5.3))
 108  Format ('xdc, ydc, req, awob, bwob')

 109  Format (' T=',F7.3,' S=',F10.7,' A=',F10.7,
     +        ' X=',F8.5,' Y=',F8.5)
 110  Format ( " PROGRAM EM_2D",/,/)
 111  Format ( " Executing step ",I4, " Total steps: ",I4)
 112  Format (10(1X,I4))
 113  Format (1X,F8.5,1X,I6)
 115  Format (10(1X,f10.5))

 140  Format (1X,I4,1x,f10.6,20(1X,F8.5))
 145  Format (1X,I4,1x,f15.6,2(1X,F8.5),1X,f15.5)

 169  Format (1X,F8.4,100(1x,f7.5))

 305  Format (' ARC LENGTHS   ',10(/,10(1x,F6.3)))
 306  Format (' SURFACE AREAS ',10(/,10(1x,F6.3)))
 307  Format (' X-CENTERS     ',10(/,10(1x,F6.3)))
 308  Format (' Y-CENTERS     ',10(/,10(1x,F6.3)))

      Stop
      End
