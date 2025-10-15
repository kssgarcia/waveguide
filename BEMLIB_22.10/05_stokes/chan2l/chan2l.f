      program chan2l

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------
c Simulation of two-layer Stokes flow in a channel
c in the presence of an insoluble surfactant
c
c NOTATION:
c --------
c
c th0: inclination angle
c RL:  period
c
c NSG:   Number of segments along the free surface
c x, y:  Interfacial nodes
c c:     surfactant concentration
c srtn:  surface tension
c
c X_sv, Y_sv: saved for RK2
c-------------------------------------------------

      Implicit double precision (a-h,o-z)

      Dimension   X(0:513),  Y(0:513)

      Dimension vnx(0:513), vny(0:513)
      Dimension   s(0:513), crv(0:513)
      Dimension   c(0:513),srtn(0:513)
      Dimension  c0(0:513)

      Dimension  Ux(0:513), Uy(0:513)
      Dimension  Un(0:513), Ut(0:513)

      Dimension    X_sv(0:513), Y_sv(0:513)
      Dimension   Ux_sv(0:513),Uy_sv(0:513)

      Dimension time(0:10000)

      Dimension ZZ(20),WW(20)

c--------------
c common blocks
c--------------

      common/XXYY/X,Y

      common/ZZWW/ZZ,WW

      common/ppii/pi,pih,pi2,pi4,pi6,pi8

c----------
c constants
c----------

      pi  = 3.14159 265358 979323 84 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi
      pi8 = 8.0D0*pi

      None = 1
      Ntwo = 2

      zero = 0.0D0

c-----
c input
c------

      open (1,file="chan2l.dat")

      read (1,*) h
      read (1,*) th0
      read (1,*) 
      read (1,*) RL
      read (1,*) 
      read (1,*) Y0
      read (1,*) 
      read (1,*) a0
      read (1,*) ph0
      read (1,*) 
      read (1,*) a0c
      read (1,*) ph0c
      read (1,*) 
      read (1,*) rho1,rho2
      read (1,*) vs1,vs2
      read (1,*)
      read (1,*) tinit
      read (1,*) cinit
      read (1,*) 
      read (1,*) beta
      read (1,*) Ds
      read (1,*) 
      read (1,*) V1
      read (1,*) V2
      read (1,*) 
      read (1,*) chi
      read (1,*) 
      read (1,*) gac
      read (1,*) 
      read (1,*) NGL
      read (1,*) 
      read (1,*) NSG
      read (1,*) Dt
      read (1,*) Nstep
      read (1,*) tend
      read (1,*) 
      read (1,*) Nprint
      read (1,*) 
      read (1,*) Iread
      read (1,*) 
      read (1,*) Move
      read (1,*) 
      read (1,*) Ich1,THmax
      read (1,*) Ich2,SPmin
      read (1,*) Ich3,SPmax
      read (1,*) Ich4,thrs1,thrs2

      close(1)
       
      write (6,*) " chan2l: input data read"

c--------
c prepare - wave length
c--------

      if(RL.eq.95) then
           RL=pi2
           write (6,*) "Period L =",RL
           pause
      else if(RL.eq.96) then
           RL=pi4
           write (6,*) "Period L =",RL
           pause
      else if(RL.eq.97) then
           RL=pi8
           write (6,*) "Period L =",RL
           pause
      end if

c--------
c prepare - inclination
c--------

      cth0 = Dcos(th0*pi)
      sth0 = Dsin(th0*pi)

      gx =   gac*sth0
      gy = - gac*cth0

      RLH = 0.5D0*RL

c--------
c prepare: channel
c--------

      hm = -h
      hs = h**2

c--------
c prepare: viscosity and density
c--------

      vsr = vs2/vs1      ! viscosity ratio

      Drho  = rho1-rho2    ! density difference

      betac = 1.0D0-beta        ! surfactant
      alpha = tinit*RL/(vs1*Ds) ! surfactant

c--------
c prepare - more
c--------

c     shiftx = - RL*cth0
c     shifty =   RL*sth0 

      shiftx = 0.0D0
      shifty = 0.0D0

      RLM  = -RL

      thmax = thmax * pi

c     ths1 = thrs1 * RL + 0.5*RL
c     ths2 = thrs2 * RL + 0.5*RL
      ths1 = thrs1 * RL 
      ths2 = thrs2 * RL 

      srfam0 = cinit*RL

      call gauss_leg (NGL,ZZ,WW)

c-------
c recite
c-------

      write (6,*)
      write (6,*) "chan2l: chanel semi-width    = ",h
      write (6,*) "chan2l: property group alpha = ",alpha
      write (6,*) "chan2l: beta                 = ",beta
      write (6,*)

c-------
c open files
c-------

      open (2,file="chan2l.xy")    ! profiles
      open (7,file="chan2l.diag")  ! diagnostics

c-----------------------------
c print inclined channel walls
c-----------------------------

      tmp1 =  RLM*cth0 + hm*sth0
      tmp2 = -RLM*sth0 + hm*cth0
      tmp3 =  RL *cth0 + hm*sth0
      tmp4 = -RL *sth0  + hm*cth0
      tmp5 =  RLM*cth0 + h *sth0
      tmp6 = -RLM*sth0 + h *cth0
      tmp7 =  RL *cth0 + h *sth0
      tmp8 = -RL *sth0 + h *cth0

      write (2,103) Ntwo
      write (2,102) None,tmp1,tmp2,cinit
      write (2,102) Ntwo,tmp3,tmp4,cinit

      write (2,103) Ntwo
      write (2,102) None,tmp5,tmp6,cinit
      write (2,102) Ntwo,tmp7,tmp8,cinit

c----------------------------------
c generate the free surface profile 
c----------------------------------

c------------------------
      if(Iread.eq.1) then        ! generate
c------------------------

        time(1) = 0.0D0

        NSGd = NSG-4
        NSGc = NSG-3
        NSGb = NSG-2
        NSGa = NSG-1

        NSG1 = NSG+1
        NSG2 = NSG+2
        NSG3 = NSG+3
        NSG4 = NSG+4

        wn = pi2/RL               ! wave number

        Dx = RL/(NSG1-1.0D0)

        Do i=0,NSG2

          X(i) = -RLH+(i-1.0D0)*Dx

          arg  = wn*x(i) - ph0 * pi
          Y(i) = Y0 + a0*Dcos(arg)

          arg  = wn*x(i) - ph0c * pi
          c(i) = cinit*(1.0D0 + a0c*Dcos(arg))

c         write (6,100) i,x(i),y(i),c(i)

        End Do

c       pause
        
c---------
      else   ! read
c---------

        open (4,file='chan2l.inp')

        read (4,*) NSG1,time(1)

        Do j=1,NSG1

          read (4,*) idle,X(j),Y(j),c(j)

          X(j) = X(j)-shiftx         ! it was shifted and tilted
          Y(j) = Y(j)-shifty
          tmx  = X(j)*cth0 - Y(j)*sth0
          tmy  = X(j)*sth0 + Y(j)*cth0
          X(j) = tmx
          Y(j) = tmy

          write (6,100) J,X(j),Y(j),c(j)

        End Do

        close (4)

        NSG  = NSG1-1
        NSGd = NSG-4
        NSGc = NSG-3
        NSGb = NSG-2
        NSGa = NSG-1

        NSG1 = NSG+1
        NSG2 = NSG+2
        NSG3 = NSG+3
        NSG4 = NSG+4

c-----------
      end if                  ! end of defining the initial profile
c-----------

      write (6,*) " chan2l: surface nodes generated"
      write (6,*) " chan2l: number of surface segments: ",NSG

c-----
c wrap
c-----

       X(0) = X(NSG) -RL
       Y(0) = Y(NSG)
       c(0) = c(NSG)

       X(NSG1) = X(1)+RL
       Y(NSG1) = Y(1)
       c(NSG1) = c(1)

       X(NSG2) = X(2)+RL
       Y(NSG2) = Y(2)
       c(NSG2) = c(2)

       X(NSG3) = X(3)+RL
       Y(NSG3) = Y(3)
       c(NSG3) = c(3)

       X(NSG4) = X(4)+RL
       Y(NSG4) = Y(4)
       c(NSG4) = c(4)

c-----------
c initialize
c-----------

      Dth = 0.50D0*Dt

      Kstep  = 1       ! step counter
      Iprint = Nprint
      Itry   = 1       ! graphics

      write (6,*) " chan2l: begin time stepping"
      write (6,*) 
      write (6,*) "       time, amp,xmax,xmin, ampc,xmaxc,xminc"

c-.-.-.-.-.-.-.-.-.-.-.-.-.
c THE SIMULATION BEGINS HERE
c-.-.-.-.-.-.-.-.-.-.-.-.-.

  90  Continue

c---------------------
c point redistribution
c---------------------

c     write (6,*) " chan2l: point redistribution in progress"

      call prd_2d_pr_splc
     +
     +   (NSG
     +   ,X,Y
     +   ,RL
     +   ,Ich1,THMAX
     +   ,Ich2,SPMAX
     +   ,Ich3,SPMIN
     +   ,Ich4,ths1,ths2
     +   ,c
c    +   ,P2,P3
c    +   ,P4,P5
     +   ,Istop
     +   )

      if(Istop.eq.1) Go to 99

      NSG1 = NSG+1
      NSG2 = NSG+2
      NSG3 = NSG+3
      NSG4 = NSG+4

c     Do i=0,NSG2
c      write (6,300) i,XI(i),YI(i),c(i)
c     End Do
c     pause

c-------------------
c compute the area
c of the lower fluid
c-------------------

      area = 0.0D0

      Do i=1,NSG
       i1 = i+1
       Dx = X(i1)-X(i)
       area = area + 0.50D0*(Y(i)+Y(i1)+2.0D0*h)*Dx
      End Do

c-----------------------------
c compute the max and min of the interface
c by cubic spline interpolation
c-----------------------------

      call splc_mnmx
     +
     +   (NSG,X,Y
     +   ,Y
     +   ,Ymax,Xmax
     +   ,Ymin,Xmin
     +   )

       amp = 0.50D0*(Ymax-Ymin)

c      write (6,*) Ymax,Ymin
c      pause

       Pmax = Xmax
       if(Pmax.gt.RLH) Pmax = Pmax - RL

       Pmin = Xmin

c--------------------------------
c compute the max and min
c of the surfactant concentration
c--------------------------------

      call splc_mnmx
     +
     +   (NSG,X,Y
     +   ,c
     +   ,cmax,Xcmax
     +   ,cmin,Xcmin
     +   )

       ampc = 0.50D0*(cmax-cmin)

       Pcmax = Xcmax
       if(Pcmax.gt.RLH) Pcmax = Pcmax - RL

       Pcmin = Xcmin

c------------------
c printing diagnostics
c------------------

       if(a0.lt.0.000001) then
         a0red = h
       else
         a0red = a0
       end if

       if(amp.lt.0.0000001) then
         amp = 1.0D0
       end if

       if(ampc.lt.0.0000001) then
         ampc = 1.0D0
       end if

       if(a0c.lt.0.000001) then
         a0cred = 1.0D0
       else
         a0cred = a0c
       end if

       write (6,101) Kstep,time(Kstep),Dlog(amp/a0red),Pmax,Pmin
     +                  ,Dlog(ampc/a0cred),Pcmax,Pcmin
c    +                  ,area
c      write (6,101) Kstep,time(Kstep),Amp,Pmax,Pmin
c    +                  ,ampc,Pcmax,Pcmin
c    +                  ,area

       write (7,101) Kstep,time(Kstep),Dlog(amp/a0red),Pmax,Pmin
     +                  ,Dlog(ampc/a0cred),Pcmax,Pcmin
c    +                  ,area
c      write (7,101) Kstep,time(Kstep),amp,Pmax,Pmin
c    +                  ,Ampc,Pcmax,Pcmin
c    +                  ,area

c----------------------------------
c Surface tension equation of state
c----------------------------------

      Do i=0,NSG1
       srtn(i) = tinit/betac * (1.0D0-beta*c(i)/cinit)
      End Do

c----------------
c stopping checks
c----------------

      if(Kstep.eq.Nstep+1) Go to 99

      if(time(Kstep).ge.tend) Go to 99

c--------------------------------------
c compute the normal vector
c and curvature by spline interpolation
c--------------------------------------

      call splc_pr_geo
     +
     +   (NSG
     +   ,X,Y
     +   ,vnx,vny
     +   ,crv
     +   ,S
     +   )

c     Do i=1,NSG1
c       write (6,106) i,S(i),vnx(i),vny(i),crv(i)
c     End Do
c
c     pause

c-----------------------------
c solve the integral equations
c-----------------------------

      call chan2l_slv
     +
     +  (RL
     +  ,h
     +  ,NGL
     +  ,NSG
     +  ,rho1,vs1
     +  ,rho2,vs2
     +  ,V1,V2
     +  ,chi
     +  ,gx,gy
     +  ,vnx,vny
     +  ,crv
     +  ,srtn
     +  ,Ux,Uy
     +  ,Un
     +  ,Ut
     +  )

c----------------
c view velocities
c----------------

c     write (6,*) " chan2l: velocities:"
c
c     Do i=1,NSG1
c       write (6,113) i,X(i),Y(i),Ux(i),Uy(i)
c     End Do

c---------
c printing
c---------

      if(Iprint.eq.Nprint) then

c     write (6,*) NSG1,time(Kstep)
      write (2,*) NSG1,time(Kstep)

      Do i=1,NSG1
        xt = X(i)
        yt = Y(i)
        xpr =  cth0 * xt + sth0 * yt
        ypr = -sth0 * xt + cth0 * yt
        xpr = xpr+shiftx
        ypr = ypr+shifty
c       xpr = X(i)
c       ypr = Y(i)
        write (6,113) i,xpr,ypr,c(i)/cinit
c    +                      ,srtn(i)/tinit
     +                      ,Un(i),Ut(i)
        write (2,113) i,xpr,ypr,c(i)/cinit
c    +                      ,srtn(i)/tinit
     +                      ,Un(i),Ut(i)
      End Do

      Iprint = 0

      end if

c-------------------------------------
c update the  surfactant concentration 
c-------------------------------------

c---------------------------
c Compute the concentration at mid-nodes
c from end-nodes
c
c The mid-node distribution may have changed
c due to point redistribution
c----------------------------

       Do i=1,NSG
        c0(i) = 0.50D0*(c(i)+c(i+1))
       End Do

       c0(0)    = c0(NSG)
       c0(NSG1) = c0(1)
       c0(NSG2) = c0(2)
       c0(NSG3) = c0(3)
       c0(NSG4) = c0(4)

c-------------------------
c update the concentration
c-------------------------

      call chan2l_cd
     +
     + (NSG
     + ,c
     + ,c0
     + ,Ds
     + ,Dt
     + ,Move
     + ,X,Y
     + ,crv
     + ,Ux,Uy
     + ,Un,Ut
     + ,srfam0
     + ,srfam
     + )

c      write (6,101) Kstep,time(Kstep),srfam

c---------------------------
c update the surface tension 
c---------------------------

      Do i=0,NSG1
       srtn(i) = tinit/betac * (1.0D0-beta*c(i)/cinit)
      End Do

c---------------------------------
c update the interface node position
c by RK2
c---------------------------------

      Do i=1,NSG1

       if(Move.eq.0) then
         Uxmove = Ux(i)
         Uymove = Uy(i)
       else if(Move.eq.1) then
         Uxmove = Un(i) * vnx(i)
         Uymove = Un(i) * vny(i)
       end if

       X_sv(i) = X(i)
       Y_sv(i) = Y(i)

       Ux_sv(i) = Uxmove
       Uy_sv(i) = Uymove

       X(i) = X(i)+ Dt * Uxmove
       Y(i) = Y(i)+ Dt * Uymove

      End Do

      X(0) = X(NSG)-RL
      Y(0) = Y(NSG)

      X(NSG2) = X(2)+RL
      Y(NSG2) = Y(2)

      call splc_pr_geo
     +
     +   (NSG,X,Y
     +   ,vnx,vny
     +   ,crv
     +   ,s
     +   )

      call chan2l_slv
     +
     +  (RL
     +  ,h
     +  ,NGL
     +  ,NSG
     +  ,rho1,vs1
     +  ,rho2,vs2
     +  ,V1,V2
     +  ,chi
     +  ,gx,gy
     +  ,vnx,vny
     +  ,crv
     +  ,srtn
     +  ,Ux,Uy,Un,Ut
     +  )

      Do i=1,NSG1

       if(Move.eq.0) then
         Uxmove = Ux(i)
         Uymove = Uy(i)
       else if(Move.eq.1) then
         Uxmove = Un(i) * vnx(i)
         Uymove = Un(i) * vny(i)
       end if
       
       X(i) = X_sv(i)+ Dth * (Uxmove+Ux_sv(i))
       Y(i) = Y_sv(i)+ Dth * (Uymove+Uy_sv(i))

      End Do

      X(0) = X(NSG)-RL
      Y(0) = Y(NSG)

      X(NSG2) = X(2)+RL
      Y(NSG2) = Y(2)

c------------------------
c reset counters and time
c------------------------

      Kstep  = Kstep+1
      Iprint = Iprint+1

      time(Kstep) = time(Kstep-1)+Dt

c------------------------
c return for another step
c------------------------

      Go to 90

c-.-.-.-.-.-.-.-.-.-.-.
c END OF THE SIMULATION
c-.-.-.-.-.-.-.-.-.-.-.

 99   Continue

c-----
c Done
c-----

      write (2,*) Null
      write (7,*) Null

      close (2)
      close (7)

  100 Format (1X,I4,10(1x,F8.5))
  101 Format (1X,I4,1X,F12.5,20(1x,F10.6))
  102 Format (1X,I3,200(1x,F9.6))
  103 Format (1X,I4,10(1x,F10.5))
  104 Format (1X,I5,1X,F12.8,1X,F12.8,1X,I1)
  106 Format (1x,i4,20(1x,f15.5))
  109 Format (200(1x,F10.3))
  113 Format (1X,I4,10(1x,F10.6))

      stop
      end
