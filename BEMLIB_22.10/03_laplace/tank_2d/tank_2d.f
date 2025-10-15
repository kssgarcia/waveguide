      program tank_2d

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c--------------------------------------------------
c Dynamical simulation of the sloshing of an inviscid
c liquid inside a (two-dimensional) rectangular tank
c
c  LEGEND:
c  ------
c
c  N: 	Number of segments along the free surface
c
c  x ,y :  coordinates of points along the free surface
c  xs,ys:  saved coordinates
c  phi:	   potential at nodes
c  ps:	   saved potential at nodes
c  crv:	   curvature at nodes
c
c  wall1: One vertical wall located at x = wall1
c  wall2: One vertical wall located at x = wall2
c  wall3:   Horizontal wall located at y = wall3
c
c  NGL:   Number of Gauss-Legendre points for numerical integration
c         over the free-surface elements
c
c  xm,ym:    Coordinates of mid-points
c  vnx,vny:  Coordinates of normal vector at mid-points
c
c  rho:   fluid density
c  gamma: surface tension
c  rmu:   viscous damping coefficient
c
c  h:     unperturbed free surface located at y = h
c  a0:    initial amplitude of the free surface
c  a0phi: initial amplitude of the potential along the free surface
c
c  accx: x-component of the container acceleration
c  accy: y-component of the container acceleration
c
c  accxt: duration of the x-component of the container acceleration
c  accyt: duration of the y-component of the container acceleration

c  Unm:	 normal velocity at mid-points
c  phim: potential at mid-points
c
c  slp(i,j):  integral of the slp at ith point over jth segment
c  dlp(i,j):  integral of the dlp at ith point over jth segment
c
c  al:	arc-length of segments
c
c  AB * SLN = BM: Linear system for the normal velocity
c
c  Nprint:   Will print a profile after Nprint time steps
c  Nsmooth:  Will smooth profile and potential after Nsmooth steps
c
c  gac:   acceleration of gravity
c
c  Move = 0    marker points move with total velocity
c         1    marker points move with normal velocity
c
c  Iread = 0   Will generate the initial condition
c          1   Will read from file: tank_2d.inp
c
c  Iflp: index for the extrapolation of the contact line
c        1 for linear extrapolation
c        2 for quadratic extrapolation
c
c FILES:
c -----
c
c tank_2d.dat  parameters of computation 
c tank_2d.inp  initial condition
c tank_2d.xy   profiles of the free surface 
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(400),y(400),xs(400),ys(400) 
      Dimension phi(400),ps(400)
      Dimension xm(400),ym(400),elml(400),Unm(400),phim(400)
      Dimension vnx(400),vny(400),crv(400)

      Dimension Ut(400),Un(400)
      Dimension tnx(400),tny(400)
      Dimension vnxm(400),vnym(400)

      Dimension slp(400,400),dlp(400,400)

      Dimension AB(500,500),BM(500),SLN(500)

      Dimension ZZ(20),WW(20)

c---
c common blocks
c---

      common/vnxym/vnxm,vnym,elml
      common/GF/wall1,wall2,wall3
      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

      Null = 0
      None = 1
      Ntwo = 2
      Nthr = 3
      Nfou = 4

      zero = 0.0D0
      tot  = 3.0D0/2.0D0

c--------------------
c Read run parameters
c--------------------

      open (1,file="tank_2d.dat")

        read (1,*) wall1       ! first  wall at x = wall1
        read (1,*) wall2       ! second wall at x = wall2
        read (1,*) wall3       ! third  wall at y = wall3
        read (1,*) h
        read (1,*) a0
        read (1,*) a0phi
        read (1,*) gac
        read (1,*) 
        read (1,*) rho
        read (1,*) gamma
        read (1,*) rmu
        read (1,*) 
        read (1,*) N
        read (1,*) NGL
        read (1,*) Move
        read (1,*) Iflp
        read (1,*) 
        read (1,*) Dt          ! time step
        read (1,*) 
        read (1,*) accx,accy    ! protocol of acceleration
        read (1,*) accxt,accyt  ! protocol of acceleration
        read (1,*) 
        read (1,*) Nstep
        read (1,*) Nsmooth
        Read (1,*) Nprint
        read (1,*) 
        read (1,*) Iread

      close (1)

c--------
c prepare
c--------

      call gauss_leg (NGL,ZZ,WW)

      width = wall2-wall1

      top = wall3+h+h     ! top wall; for plotting purposes

c-----------------------------
c generate the initial
c position of the free surface
c distribution of potential
c over the free surface
c------------------------------

c-------------------------
      if(Iread.eq.0) then    ! Generate the initial distribution
c-------------------------

        N1 = N+1
        N7 = N+7
        wn = pi2/width     ! wave number

        Dx = width/(N1-1.0D0)

        Do i=1,N1
         tmp = (i-1.0D0)*Dx
           x(i) = wall1+tmp 
           y(i) = wall3+h+a0*cos(wn*tmp)
         phi(i) =      a0phi*cos(wn*tmp)
        End Do

      time = 0.0D0

c---------
      else
c---------

        open (4,file="tank_2d.inp",status="unknown")   ! read

          read (4,*) N1,time
          Do i=1,N1
            read (4,100) ic,x(i),y(i),phi(i)
          End Do

          N = N1-1
          N7= N+7

        close (4)

c-----------
      end If
c-----------

c-----------------
c open output file
c-----------------

      open (2,file="tank_2d.xy")

c-----------
c initialize
c-----------

      Istep   = 1      ! step counter
      Ismooth = 1      ! steps for smoothing

      Iprint  = Nprint

c-.-.-.-.-.-.-.-.-.-.-.-
c SIMULATION BEGINS HERE
c-.-.-.-.-.-.-.-.-.-.-.-

      Do Istep=1,Nstep

c---------------------------
c smooth the free surface
c and the potential
c by the five-point formula
c---------------------------

      if(Ismooth.eq.Nsmooth) then

c     write (6,*)
c     write (6,*) " tank_2d: smoothing the free surface and the potential"
c     write (6,*)

        Do i=1,N1        ! save old values
          xs(i) =   x(i)
          ys(i) =   y(i)
          ps(i) = phi(i)
        End Do

        Do i=3,N-1
          x(i) = (-xs(i-2)+4.0*xs(i-1)+10.0*xs(i)+4.0*xs(i+1)
     +            -xs(i+2))/16.0
          y(i) = (-ys(i-2)+4.0*ys(i-1)+10.0*ys(i)+4.0*ys(i+1)
     +            -ys(i+2))/16.0
        phi(i) = (-ps(i-2)+4.0*ps(i-1)+10.0*ps(i)+4.0*ps(i+1)
     +            -ps(i+2))/16.0
        End Do

        Ismooth = 0

      end if

c----------------------------------
c Printing the free-surface profile
c and the walls
c----------------------------------

      if(Iprint.eq.Nprint) then

        write (2,100) N7,time

        Do i=1,N1
         write (2,100) ic,x(i),y(i),phi(i)
        End Do

        ic = N1+1
        write (2,100) ic,wall2,y(N1)
        ic = ic+1
        write (2,100) ic,wall2,wall3
        ic = ic+1
        write (2,100) ic,wall1,wall3
        ic = ic+1
        write (2,100) ic,wall1,top
        ic = ic+1
        write (2,100) ic,wall2,top
        ic = ic+1
        write (2,100) ic,wall2,y(N1)

        Iprint = 0

      End If
   
      write (6,200) Istep,time

c------------------------
c Define segment mid-points
c
c Compute: phi at mid-points
c          element arc length
c          normal vector pointing into the fluid
c------------------------

      Do i=1,N

       xm(i) = 0.5D0*(x(i+1)+x(i))
       ym(i) = 0.5D0*(y(i+1)+y(i))
       elml(i) = Dsqrt((y(i+1)-y(i))**2+(x(i+1)-x(i))**2)

       phim(i) = 0.5D0*(phi(i+1)+phi(i))

       vnxm(i) =   (y(i+1)-y(i))/elml(i)   ! normal vector
       vnym(i) = - (x(i+1)-x(i))/elml(i)   ! normal vector

      End Do

c-----------------------------------------------
c Compute the single-layer and double-layer potential
c at the elements mid-points integrated over the segments
c
c Set up a linear system for the normal velocity
c-----------------------------------------------

      Do i=1,N

       x0 = xm(i)
       y0 = ym(i)

       Do j=1,N

        j1 = j+1
        x1 = x(j)
        y1 = y(j)
        x2 = x(j1)
        y2 = y(j1)

        call tank_2d_sdlp 
     +
     +   (x0,y0
     +   ,i,j
     +   ,x1,y1
     +   ,x2,y2
     +   ,NGL
     +   ,slp(i,j)
     +   ,dlp(i,j)
     +   )

        AB(i,j) = slp(i,j)       ! influence matrix

       End Do

c---
c compute the right-hand side (rhs)
c---

       accum = 0.0D0

       Do j=1,N
        accum = accum + phim(j)*dlp(i,j)
       End Do

       BM(i) = accum-0.5D0*phim(i)     ! right-hand side

      End Do

c---------------------------
c  display the linear system 
c---------------------------

c     Do  i=1,N
c       write (6,101) i,(AB(i,j),j=1,N),BM(i)
c     End Do

c-------------------------
c  Solve the linear system
c-------------------------

      Isym_gel = 0   ! system is not symmetric
      Iwlpvt   = 1   ! pivoting enabled

      call gel
     +
     +   (N,AB
     +   ,BM
     +   ,sln
     +   ,Isym_gel,Iwlpvt
c    +   ,l,u
     +   ,det
     +   ,Istop
     +   )

      if(Istop.eq.1) Go to 99

c------------------------
c Distribute the solution
c------------------------

      Do i=1,N
        Unm(i) = sln(i)
      End Do

c--------------------------------------------
c Compute variables at element end-nodes
c
c Compute: normal velocity by interpolation
c          tangential velocity by numerical differentiation
c          normal vector
c          tangential vector 
c          curvature
c--------------------------------------------

c     write (6,*)
c     write (6,*) "  un,phi,Ut,vnx,vny,tnx,tny,crv"
c     write (6,*)

      Do i=2,N

      dn = elml(i-1)+elml(i)
      Un(i) = (Unm(i)*elml(i-1)+Unm(i-1)*elml(i))/dn

      x1  = -elml(i-1)
      x2  =  elml(i)
      x21 =  x2-x1

c     Ut(i) = (phi(i)-phi(i-1))     ! linear differentiation
c    +       /(elml(i)+elml(i+1))

c---
c quadratic differentiation for phi
c with respect to arc length
c---

      y1  = phi(i-1)  
      y0  = phi(i)
      y2  = phi(i+1)
      aa  = ((y2-y0)/x2-(y1-y0)/x1)/x21
      bb  =  (y2-y0)/x2 - aa*x2
      Ut(i) = bb

c---
c quadratic differentiation for x
c with respect to arc length
c---

      y1 = x(i-1)                  ! quadratic differentiation
      y0 = x(i)
      y2 = x(i+1)
      aa = ((y2-y0)/x2-(y1-y0)/x1)/x21
      bb =  (y2-y0)/x2 - aa*x2
      xp = bb
      xpp= 2.0D0*aa

c---
c quadratic differentiation for y
c with respect to arc length
c---

      y1 = y(i-1)                  ! quadratic differentiation
      y0 = y(i)
      y2 = y(i+1)
      aa = ((y2-y0)/x2-(y1-y0)/x1)/x21
      bb =  (y2-y0)/x2 - aa*x2
      yp = bb
      ypp= 2.0D0*aa

      tnm    = dsqrt(xp*xp+yp*yp)
      tnx(i) = xp/tnm
      tny(i) = yp/tnm

      vnx(i) =  tny(i)    ! normal vector points into the fluid
      vny(i) = -tnx(i)

c---
c curvature
c---

      crv(i) = (xpp*yp-ypp*xp)/(xp*xp+yp*yp)**tot

c     write (6,101) i,Un(i),phi(i),Ut(i),vnx(i),vny(i)
c    +              ,tnx(i),tny(i),crv(i)

      End Do

c-------------------------------------------
c time stepping
c
c update phi using Bernoulli's equation
c-------------------------------------------

      Do i=2,N

       press  = gamma*crv(i)                ! capillary pressure
       Ums    = Un(i)**2 + Ut(i)**2
       phi(i) = phi(i) - Dt*rmu*phi(i)
       yref   = y(i)-(wall3+h)

c-----------------------
       if(Move.eq.0) then    ! total velocity
c-----------------------

        x(i)   = x(i) + Dt*( Un(i)*vnx(i)+Ut(i)*tnx(i) )
        y(i)   = y(i) + Dt*( Un(i)*vny(i)+Ut(i)*tny(i) )
        phi(i) = phi(i) + Dt*(0.5D0*Ums-press/rho-gac*yref) 

c----------
       else                  ! normal velocity
c----------

        x(i)   = x(i) + Dt*Un(i)*vnx(i)
        y(i)   = y(i) + Dt*Un(i)*vny(i)
        phi(i) = phi(i) + Dt* (Un(i)*Un(i) 
     +                        -0.5D0*Ums
     +                        -press/rho
     +                        -gac*yref) 

c------------
       end if
c------------

c---
c account for the tank acceleration
c---

        if(time.le.accxt) phi(i) = phi(i) - Dt*accx*x(i)
        if(time.le.accyt) phi(i) = phi(i) - Dt*accy*y(i)

      End Do

c------------------------------------------
c Compute y position and phi of first point
c------------------------------------------

      if(Iflp.eq.1) then   ! linear extrapolation

        x1 = elml(1)
        x2 = x1+elml(2)
        y1 = y(2)
        y2 = y(3)
        y(1) = y1*(0.0-x2)/(x1-x2)+y2*(0.0-x1)/(x2-x1)
        y1 = phi(2)
        y2 = phi(3)
        phi(1) = y1*(0.0-x2)/(x1-x2)
     +          +y2*(0.0-x1)/(x2-x1)

      else if(Iflp.eq.2) then   ! quadratic extrapolation

        x1 = elml(1)
        x2 = x1+elml(2)
        x3 = x2+elml(3)
        y1 = y(2)
        y2 = y(3)
        y3 = y(4)
        y(1) = y1*(0.0-x2)*(0.0-x3)/((x1-x2)*(x1-x3))
     +       + y2*(0.0-x1)*(0.0-x3)/((x2-x1)*(x2-x3))
     +       + y3*(0.0-x1)*(0.0-x2)/((x3-x1)*(x3-x2))
        y1 = phi(2)
        y2 = phi(3)
        y3 = phi(4)
        phi(1) = y1*(0.0-x2)*(0.0-x3)/((x1-x2)*(x1-x3))
     +         + y2*(0.0-x1)*(0.0-x3)/((x2-x1)*(x2-x3))
     +         + y3*(0.0-x1)*(0.0-x2)/((x3-x1)*(x3-x2))

       end if

c-----------------------------------------
c Compute y position and phi of last point
c numbered N+1
c-----------------------------------------

      if(Iflp.eq.1) then       ! linear extrapolation

        x1 = elml(N)
        x2 = x1+elml(N-1)
        y1 = y(N)
        y2 = y(N-1)
        y(N1) = y1*(0.0-x2)/(x1-x2)+y2*(0.0-x1)/(x2-x1)
        y1 = phi(N)
        y2 = phi(N-1)
        phi(N1) = y1*(0.0-x2)/(x1-x2)
     +           +y2*(0.0-x1)/(x2-x1)

      else if(Iflp.eq.2) then   ! quadratic extrapolation

        x1 = elml(N)
        x2 = x1+elml(N-1)
        x3 = x2+elml(N-2)
        y1 = y(N)
        y2 = y(N-1)
        y3 = y(N-2)
        y(N1) = y1*(0.0-x2)*(0.0-x3)/((x1-x2)*(x1-x3))
     +        + y2*(0.0-x1)*(0.0-x3)/((x2-x1)*(x2-x3))
     +        + y3*(0.0-x1)*(0.0-x2)/((x3-x1)*(x3-x2))
        y1 = phi(N)
        y2 = phi(N-1)
        y3 = phi(N-2)
        phi(N1) = y1*(0.0-x2)*(0.0-x3)/((x1-x2)*(x1-x3))
     +          + y2*(0.0-x1)*(0.0-x3)/((x2-x1)*(x2-x3))
     +          + y3*(0.0-x1)*(0.0-x2)/((x3-x1)*(x3-x2))

      end if

      time = time + Dt

      Iprint  = Iprint  +1
      Ismooth = Ismooth +1

c-------------
      End Do   ! over time steps
c-------------

c---------------------
c Simulation has ended
c---------------------

   99 Continue

c     write (6,*)                                  ! graphics
c     write (6,*) "Press any key to finish"        ! graphics
c     write (6,*)                                  ! graphics
c     call getkey                                  ! graphics
c     call vexit                                   ! graphics

c--------
c wrap up
c--------

      write (2,100) Null

      write (2,201) w,d,h,a0,a0phi
      write (2,202) gamma,rmu
      write (2,203) N
      write (2,204) NGL,Move,Iflp,Nsmooth
      write (2,205) Dt
      write (2,206) accx,accy
      write (2,207) accxt,accyt
      write (2,*)

      close (2)

c-------------
c restart file
c-------------

      open (2,file="tank_2d.rst")

        write (2,100) N7,time

        Do i=1,N1
         write (2,100) ic,x(i),y(i),phi(i)
        End Do

      close (2)

c-----
c Done
c-----

  100 Format (1x,i3,9(1x,f10.5))
  200 Format (1x,"Step: ",I3,"  time=",f10.5)
  101 Format (1x,i3,40(1x,f6.3))

  201 Format (1x," W=",f5.3," D=",f5.3," H=",f5.3,
     +           " a0=",f5.3," a0phi=",f5.3)
  202 Format (1x," Surf Tens=",f5.3," Damping=",f5.3)
  203 Format (1x," N=",I2," M=",I2," K=",I2," L=",I2)
  204 Format (1x," NGL=",I2," Move=",I1," Iflp=",I1,
     +           " Nsm=",I2)
  205 Format (1x," time step = ",f5.3)
  206 Format (1x," accx = ",f8.3," accy = ",f8.3)
  207 Format (1x," accxt= ",f8.3," accyt= ",f8.3)

      Stop
      End
