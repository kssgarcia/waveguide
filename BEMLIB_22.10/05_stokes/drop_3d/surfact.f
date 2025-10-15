      program surfact

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c-----------------------------------------------
c solve the convection -- diffusion 
c for an insoluble surfactant
c over a 3D interface.
c
c SYMBOLS:
c --------
c
c  npts	   total number of points
c  nelm	   total number of elements
c
c  p(i,j)   coordinates of nodes i (j=1,2,3)
c
c  ne(k,j)  ne(k,1) is the number of elements adjacent to point k
c           ne(k,2), ... are the elements numbers, j = 2, ..., 7
c           For this triangulation, ne(k,j) is up to six
c
c  n(k,i)    connectivity table: points for element k, i = 1,...,6
c
c  nbe(k,j)  the three neighboring elements of element k (j=1,2,3)
c
c  vna       average value of the normal vector at nodes
c  u         node velocity
c  srtn(i)   surface tension at node i   
c  c(i)      surfactant concentration at node i   
c
c  psv           save p
c  Usv,Vsv,Wsv   save U,V,W
c  csv           save c
c
c    arel(k): surface area of element k
c  crvmel(k): average value of the mean curvature of element k
c     cel(k): element surfactant concentration 
c
c    xmom(k)    x-moment of element k
c    ymom(k)    y-moment of element k
c    zmom(k)    z-moment of element k
c
c  jxy          number of sequential nodes in the xy plane
c
c  dxy          Taylor deformation parameter in the xy plane
c  thmax        inclination of maximum axis
c  thmin        inclination of minimum axis
c
c  vol         volume
c  ars         surface area 
c  crx         x position of centroid   
c  cry         y position of centroid   
c  crz         z position of centroid   
c
c  zz, ww         base points and weights of Gauss-Legenrde integr
c  xiq, etq, wq   base points and weights of Gauss-triangle integr
c
c  wall    wall located at y = wall
c
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension    u(1026,3)
      Dimension    c(1026)

      Dimension Umove(1026),Vmove(1026),Wmove(1026)

      Dimension lxy(1026,2)

      Dimension crvm(1026)
      Dimension dilt(1026)

      Dimension      n(512,6), nbe(512,3)
      Dimension  alpha(512),  beta(512), gamma(512)
      Dimension   arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension crvmel(512)
      Dimension    cel(512)

      Dimension jxy(100)

      Dimension time(1000),dxy(1000),thmax(1000),thmin(1000)
      Dimension  ars(1000),vol(1000),crx(1000),cry(1000),crz(1000)

      Dimension xiq(20),etq(20),wq(20)

c---
c axisymmetric solution
c---

      Dimension tha(-2:1026),cc(-2:1026)
      Dimension ATR(1026),BTR(1026),CTR(1026)
      Dimension ATM(1026),BTM(1026),CTM(1026)
      Dimension RHS(1026),SOL(1026)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna
      common/geo3/crvm
      common/geo6/crvmel
      common/geo9/xmom,ymom,zmom

      common/surfa/c

      common/viscr/vs1,vs2,vsrt,vsrtp,vsrtm,vsf,vsk
      common/visci/Ivs

      common/var/shrt,wall

      common/veloc0/u
      common/veloc1/nvelt,nvel
      common/veloc2/nvelr,lxy

c---
c various
c---

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------
c constants
c----------

      pi = 3.14159 265358 979 32384 D0

      piq = 0.25D0*pi
      pih = 0.50D0*pi
      pi2 = 2.00D0*pi
      pi4 = 4.00D0*pi
      pi6 = 6.00D0*pi
      pi8 = 8.00D0*pi

      Null = 0
      Nseven = 7
      Nfour  = 4

      oot = 1.0D0/3.0D0

      open (2,file="surfact.dat")

      read (2,*) Ivel
      read (2,*)
      read (2,*) Ioctaicos
      read (2,*) ndiv 
      read (2,*) 
      read (2,*) req
      read (2,*) boa,coa
      read (2,*) cx,cy,cz
      read (2,*) phi1,phi2,phi3
      read (2,*)
      read (2,*) mint
      read (2,*)
      read (2,*) cinit
      read (2,*) 
      read (2,*) Ds
      read (2,*) 
      read (2,*) Ismeth
      read (2,*)
      read (2,*) Isym_xy
      read (2,*)
      read (2,*) Dt
      read (2,*)
      read (2,*) Nprint_xy
      read (2,*) Nprint_xyz
      read (2,*)
      read (2,*) Move
      read (2,*)
      read (2,*) Vinf
      read (2,*) Nax
      read (2,*)

c---------------------------------------
c preparations, adjustments, definitions
c---------------------------------------

      betsc = 1.0D0-bets       ! surfactant

c----------------------
c Axisymmetric solution
c----------------------

      If(Ivel.eq.2) then

      write (6,*) " surfact: setting up the axisymmetric sln"

      Nax1 = Nax+1

      Do i=1,Nax1
       cc(i) = cinit
      End Do

c---
c set up the tridiagonal matrix
c---

      Dthe  = pi/Nax
      Dthes = Dthe**2

      tha(1) = 0.0D0
      uth =  0.0D0
      utp =  0.25D0*Vinf
      ATR(1) = - 4.0D0*Ds - 2.0D0*utp*Dthes
      BTR(1) =   4.0D0*Ds

      Do i=2,Nax
        tha(i) = (i-1.0D0)*Dthe
        cs = Dcos(tha(i))
        sn = Dsin(tha(i))
        cotn = cs/sn
        uth = 0.25D0*Vinf*sn
        utp = 0.25D0*Vinf*cs
        ATR(i) = - 2.0D0*Ds - Dthes * (utp+ cotn*uth)
        BTR(i) = - 0.50D0*uth*Dthe + Ds*(1.0D0+0.5D0*cotn*Dthe )
        CTR(i) =   0.50D0*uth*Dthe + Ds*(1.0D0-0.5D0*cotn*Dthe )
      End Do

      tha(Nax1) = pi
      uth =  0.0D0
      utp = -0.25D0*Vinf
      ATR(Nax1) = - 4.0D0*Ds - 2.0D0*utp*Dthes
      CTR(Nax1) =   4.0D0*Ds
 
      fc = 0.50D0*Dt/Dthes

      Do i=1,Nax1
        ATR(i) = fc*ATR(i)
        BTR(i) = fc*BTR(i)
        CTR(i) = fc*CTR(i)
c       write (6,100) i,ATR(i),BTR(i),CTR(i)
      End Do

      Do i=1,Nax1
        ATM(i) = 1.0D0-ATR(i)
        BTM(i) =      -BTR(i)
        CTM(i) =      -CTR(i)
      End Do

      write (6,*) " surfact: Done"

      End If

c------------------------------------------
c save the initial position of the centroid
c------------------------------------------

      cxs = cx
      cys = cy
      czs = cz

c------------------------------
c Define integration quadrature
c------------------------------

      call gauss_trgl 
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

      write (6,*) " surfact: quadrature read"

c----------------------------
c triangulate the unit sphere
c
c Run even at restart
c to generate the connectivity table
c-----------------------------------

      write (6,*) " surfact: triangulating"

      if(Ioctaicos.eq.1) then

       call trgl6_octa
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      else

      call trgl6_icos
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      end if

      write (6,*) " surfact: number of points: ",npts
      write (6,*) " surfact: number of elements: ",nelm

      nelmh = nelm/2
      nelm2 = 2*nelm

c--------------
c initial state   ! SPHEROID
c--------------

      time(1) = 0.0D0

c---------------------------
c expand to specified shape
c and equivalent radius
c
c rotate by the angles: phi1,phi2,phi3
c translate center to specified position
c---------------------------

      scale = req/(boa*coa)**oot

      Do i=1,npts
        p(i,1) = scale*p(i,1)
        p(i,2) = scale*p(i,2)*boa
        p(i,3) = scale*p(i,3)*coa
      End do

      phi1 = phi1*pi
      phi2 = phi2*pi
      phi3 = phi3*pi

      cs = Dcos(phi1)
      sn = Dsin(phi1)

      Do i=1,npts
       tmpx = p(i,1)
       tmpy = cs*p(i,2)+sn*p(i,3)
       tmpz =-sn*p(i,2)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = Dcos(phi2)
      sn = Dsin(phi2)

      Do i=1,npts
       tmpx = cs*p(i,1)-sn*p(i,3)
       tmpy = p(i,2)
       tmpz = sn*p(i,1)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = Dcos(phi3)
      sn = Dsin(phi3)

      Do i=1,npts
       tmpx = cs*p(i,1)+sn*p(i,2)
       tmpy =-sn*p(i,1)+cs*p(i,2)
       tmpz = p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      Do i=1,npts
        p(i,1) = p(i,1) + cx
        p(i,2) = p(i,2) + cy
        p(i,3) = p(i,3) + cz
      End Do

c---
c unscale to record
c---

      phi1 = phi1/pi
      phi2 = phi2/pi
      phi3 = phi3/pi

c--------------------------------
c Assign surfactant concentration
c--------------------------------

      Do i=1,npts
       c(i) = cinit
      End Do

c---------------------------------
c display half the initial surface
c---------------------------------

c     Do k=1,nelmh
c       call printel(k)
c     End Do

c------------------------------------------
c Find the points in the z=cz plane
c and order them
c
c nxy is number of points in the xy plane
c the vector jxy contains consecutive points
c-------------------------------------------

      call xy_slice 
     +
     +  (npts
     +  ,cx,cy,cz
     +  ,nxy
     +  ,jxy
     +  )

c---------------
c prepare to run
c---------------

      open (1,file="surfact.xyz")
      open (3,file="surfact.xy")
      open (4,file="surfact.comp")
      open (9,file="surfact.axis")

      Kstep = 1           ! time step counter

      Iprint_xy  = Nprint_xy       ! printing counter
      Iprint_xyz = Nprint_xyz      ! printing counter

c--------------------------------
c element from node concentration
c by averaging
c--------------------------------

      Do i=1,Nelm
       collect = 0.0D0
       Do k=1,6
        j = n(i,k)   ! global label
        collect = collect + c(j)
       End Do
       cel(i) = collect/6.0D0
      End Do

c-.-.-.-.-.-.-.-.-.-.-.-.-.-.
c  Time stepping begins here
c-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 90   Continue

      read (2,*)  Nstep
      If(Nstep.eq.0) Go to 99

c-----------
c initialize
c-----------

      Istep  = 1      ! cycle step counter

  97  Continue

      write (6,*)
      write (6,*) "--------------------"
      write (6,*)
      write (6,105) Istep,Nstep,time(Kstep)

  91  Continue

c----------------------------------------
c compute coefficients alpha, beta, gamma
c for quadratic xi-eta mapping
c over each element
c---------------------------------------

      Do k=1,nelm

        i1 = n(k,1)
        i2 = n(k,2)
        i3 = n(k,3)
        i4 = n(k,4)
        i5 = n(k,5)
        i6 = n(k,6)

        call abc 
     +
     +   (p(i1,1),p(i1,2),p(i1,3)
     +   ,p(i2,1),p(i2,2),p(i2,3)
     +   ,p(i3,1),p(i3,2),p(i3,3)
     +   ,p(i4,1),p(i4,2),p(i4,3)
     +   ,p(i5,1),p(i5,2),p(i5,3)
     +   ,p(i6,1),p(i6,2),p(i6,3)
     +   ,alpha(k),beta(k),gamma(k)
     +   )

      End Do

c------------------------------------------------
c compute:
c
c     surface area of the individual elements
c     x, y, and z surface moments over each element
c     total surface area and volume
c     mean curvature over each element
c     average normal vector at the nodes
c------------------------------------------------

      call elm_geom 
     +
     +  (nelm,npts,mint
     +  ,area,vlm
     +  ,cx,cy,cz
     +  )

c---------------
c     write (6,*) 
c     write (6,*) " Curvature of elements"
c     write (6,*) 
c     Do k = 1,nelm
c       write (6,100) k,crvmel(k)
c     End Do
c---------------

c------------------
c recording session
c------------------

      ars(Kstep) = area/(pi4*req**2)
      vol(Kstep) = vlm/ (pi4*req**3/3.0D0)

      crx(Kstep) = cx
      cry(Kstep) = cy
      crz(Kstep) = cz

      write (6,*)
      write (6,110) ars(Istep)
      write (6,111) vol(Istep)
      write (6,112) cx,cy,cz
      write (6,*)

c--------------------
c print cross-section
c--------------------

      If(Isym_xy.eq.1) then

c       write (6,*)
c       write (6,*) " surfact: x, y,  c at t=",time(Kstep)
c       write (6,*)

c       Do i=1,nxy
c         m = jxy(i)
c         write (6,100) i,(p(m,j),j=1,2),c(m)
c       End Do

c----
        If(Iprint_xy.eq.Nprint_xy) then

          Iprint_xy = 0

          If(Iflow.eq.3) then
           write (3,100) nxy,time(Kstep),a21
          Else
           write (3,100) nxy,time(Kstep)
          End If

          Do i=1,nxy
           m = jxy(i)
           radd = Dsqrt(p(m,1)**2+p(m,2)**2+p(m,3)**2)
           theta = acos(p(m,1)/radd)
c          write (3,100) i,p(m,1),p(m,2),c(m)
           write (3,100) i,theta/pi,c(m)
          End Do

c         write (6,*) Nax1
          write (9,100) Nax1,time(Kstep)
          Do i=1,Nax1
c           write (6,100) i,tha(i)/pi,cc(i)
            write (9,100) i,tha(i)/pi,cc(i)
           End Do


        End If
c----

      End If

c----------------
c print all nodes
c----------------

      If(Iprint_xyz.eq.Nprint_xyz) then

         Iprint_xyz = 0

         If(Iflow.eq.3) then
          write (1,100) npts,time(Kstep),a21
         Else
          write (1,100) npts,time(Kstep)
         End If

         Do i=1,npts
           write (1,108) i,p(i,1),p(i,2),p(i,3),c(i)
         End Do

      End If

c---------
c Velocity
c---------

      Do i=1,Npts

      If(Ivel.eq.1) then

c----------
c expansion
c----------

        u(i,1) = p(i,1)
        u(i,2) = p(i,2)
        u(i,3) = p(i,3)

      Else If(Ivel.eq.2) then

c---
c tangential motion
c---
        u(i,1) = - Vinf *  0.25D0 * (1.D0 - p(i,1)**2 )
        u(i,2) =   Vinf *  0.25D0 * p(i,1)*p(i,2)
        u(i,3) =   Vinf *  0.25D0 * p(i,1)*p(i,3)

       End If
      End Do

c----------------------------
c Mean curvature at the nodes
c----------------------------

      call crvm_3d (nelm,npts)

c-----------------------
c Advance the surfactant
c-----------------------

      If(Ismeth.eq.1) then ! finite volume method

      write (6,*) " surfact_dr: entering surfact_fvm"

      call surfact_fvm
     +
     +  (Nelm
     +  ,Npts
     +  ,mint
     +  ,Dt
     +  ,Ds
     +  ,cel
     +  ,Move
     +  )

      Else If(Ismeth.eq.2) then ! finite element method

      write (6,*) " surfact_fem: entering surfact_fem"

      call surfact_fem
     +
     +   (Nelm
     +   ,Npts
     +   ,mint
     +   ,Dt
     +   ,Ds
     +   ,Move
     +   )

      End If

      write (6,*) " surfact_dr: exited surfactant"

c----------------------
c Axisymmetric solution
c----------------------

      If(Ivel.eq.2) then

      RHS(1) = cc(1) + ATR(1)*cc(1) 
     +               + BTR(1)*cc(2)

      Do i=2,Nax
       RHS(i) = cc(i) + ATR(i)*cc(i) 
     +                + BTR(i)*cc(i+1)
     +                + CTR(i)*cc(i-1)
      End Do

      RHS(Nax1) = cc(Nax1) + ATR(Nax1)*cc(Nax1) 
     +                     + CTR(Nax1)*cc(Nax)

      call thomas
     +
     +  (Nax1 
     +  ,ATM
     +  ,BTM 
     +  ,CTM
     +  ,RHS 
     +  ,SOL   ! solution
     +  )

      Do i=1,Nax1
       cc(i) = SOL(i)
      End Do

       cc(0)  =   cc(2)
      tha(0)  = -tha(2)

       cc(-1) =   cc(3)
      tha(-1) = -tha(3)

       cc(Nax+2) =      cc(Nax)
      tha(Nax+2) = pi2-tha(Nax)

       cc(Nax+3) =      cc(Nax-1)
      tha(Nax+3) = pi2-tha(Nax-1)

c-----
c compute the ERROR
c-----

        error = 0.0D0

        Do i=1,Npts

          radd = Dsqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
          theta = acos(p(i,1)/radd)

          Do j=1,Nax
           Diff =theta-tha(j)
           If(Diff.le.0) then
            m=j
            Go to 334
           End If
          End Do

 334      Continue

          m1 = m-2
          m2 = m-1
          m3 = m
          m4 = m+1
          m5 = m+2
          q = theta
          cinterp 
     +    = cc(m1)* (q-tha(m2))*(q-tha(m3))*(q-tha(m4))*(q-tha(m5))
     +            /( (tha(m1)-tha(m2))*(tha(m1)-tha(m3))
     +            *(tha(m1)-tha(m4))*(tha(m1)-tha(m5)) )
     +     + cc(m2)* (q-tha(m1))*(q-tha(m3))*(q-tha(m4))*(q-tha(m5))
     +     /( (tha(m2)-tha(m1))*(tha(m2)-tha(m3))
     +     *(tha(m2)-tha(m4))*(tha(m2)-tha(m5)) )
     +     + cc(m3)* (q-tha(m1))*(q-tha(m2))*(q-tha(m4))*(q-tha(m5))
     +     /( (tha(m3)-tha(m1))*(tha(m3)-tha(m2))
     +     *(tha(m3)-tha(m4))*(tha(m3)-tha(m5)) )
     +     + cc(m4)* (q-tha(m1))*(q-tha(m2))*(q-tha(m3))*(q-tha(m5))
     +     /( (tha(m4)-tha(m1))*(tha(m4)-tha(m2))
     +     *(tha(m4)-tha(m3))*(tha(m4)-tha(m5)) )
     +     + cc(m5)* (q-tha(m1))*(q-tha(m2))*(q-tha(m3))*(q-tha(m4))
     +     /( (tha(m5)-tha(m1))*(tha(m5)-tha(m2))
     +     *(tha(m5)-tha(m3))*(tha(m5)-tha(m4)) )
          perror = c(i)-cinterp
          error = error + perror**2
c         write (6,100) i,c(i),cinterp,perror

        End Do

        error = Dsqrt(error)/Npts
        write (6,100) Kstep,time(Kstep)+Dt,error
        write (4,100) Kstep,time(Kstep)+Dt,error

      End If

c------------
c Move points
c------------

      Do i=1,npts

        If(Move.eq.0) then      ! points move with total velocity

          Umove(i) = u(i,1)
          Vmove(i) = u(i,2)
          Wmove(i) = u(i,3)

        Else If(Move.eq.1) then  ! points move with normal velocity

          Uvel  = u(i,1)*vna(i,1)   ! projection of velocity
     +          + u(i,2)*vna(i,2)   ! onto the normal vector
     +          + u(i,3)*vna(i,3)

          If(Ivel.eq.2) Uvel = 0.0D0

          Umove(i) = Uvel*vna(i,1)
          Vmove(i) = Uvel*vna(i,2)
          Wmove(i) = Uvel*vna(i,3)

        End If

        p(i,1) = p(i,1) + Dt*Umove(i)
        p(i,2) = p(i,2) + Dt*Vmove(i)
        p(i,3) = p(i,3) + Dt*Wmove(i)

      End Do

c-------------------
c End of a time step
c-------------------

c------------------------
c Reset counters and time
c------------------------

      Kstep  = Kstep+1
      Istep  = Istep+1

      Iprint_xy  = Iprint_xy+1
      Iprint_xyz = Iprint_xyz+1

      time(Kstep) = time(Kstep-1) + Dt

      If(Istep.le.Nstep) Go to 97

      Go to 90     ! return for another step

c-.-.-.-.-.-.-.-.-.-.-
c  Simulation has ended
c-.-.-.-.-.-.-.-.-.-.-

  99  Continue

c-------------------
c Record final shape
c-------------------

c     Do k=1,nelmh
c       call printel(k)
c     End Do
c     write (1,100) null

      write (1,*) Npts,time(Kstep)

      Do i=1,Npts
        write (1,108) i,p(i,1),p(i,2),p(i,3),c(i)
      End Do

      write (1,100) null
      write (3,100) null

c---
c Record run parameters
c---

      write (1,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,req,boa,coa,cxs,cys,czs
     +             ,vs1,vs2,tinit
     +             ,shrt

      write (3,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,req,boa,coa,cxs,cys,czs
     +             ,vs1,vs2,tinit
     +             ,shrt

      write (6,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,req,boa,coa,cxs,cys,czs
     +             ,vs1,vs2,tinit
     +             ,shrt

      write (1,208) mint,NGL,Dt,Irk,Move,Norm,Isym_xy
      write (3,208) mint,NGL,Dt,Irk,Move,Norm,Isym_xy
      write (6,208) mint,NGL,Dt,Irk,Move,Norm,Isym_xy

c------------
c close files
c------------

      close (1)
      close (2)
      close (3)

c-------------
c restart data
c-------------

      open (2,file="surfact.rst")

       write (2,*) Npts,time(Kstep)
       Do i=1,Npts
         write (2,108) i,p(i,1),p(i,2),p(i,3),c(i)
       End Do
       write (2,100) Null,Null

      close (2)

c---------------------
c Generate Matlab file
c---------------------

      open (2,file="drop_3d.net")

      Index = 2

      If(Index.eq.1) then          ! 6-node triangles
        write (2,*) Nseven
        write (2,*) 7*Nelm
        write (2,*) Nelm
      Else If(Index.eq.2) then     ! 3-node triangles
        write (2,*) Nfour
        write (2,*) 16*Nelm
        write (2,*) 4*Nelm
      End If

      Do k=1,Nelm
       call printel (k,Index)  ! print in file "drop_3d.net"
      End Do

      close (2)

c-----
c Done
c-----

  100 Format (1x,i4,10(1x,f12.5))
  101 Format (10(1x,f12.5))
  102 Format (10(1x,f10.6))
  103 Format (1x,i3,10(1x,f10.5))
  104 Format (1x,i3,10(1x,f8.5))

  105 Format (" Step ",i3," out of ",i4,"; time :",F15.10)
  106 Format (' T=',F7.3,' S=',F10.7,' V=',F10.7
     +       ,' X=',F8.5,' Y=',F 8.5,' Z=',F8.5
     +       )
  107 Format (' T=',F7.3,' D=',F10.7,' thmax=',F10.4,
     +                               ' thmin=',F10.4)
  108 Format (1x,i4,100(1x,f15.10))
  109 Format (' T=',F7.3,' sxy=',F10.6,' sd1=',F10.6,
     +                                 ' sd2=',F10.6)

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)
  112 Format(" Centroid     :",3(F15.10))
  113 Format("Time=",F12.5," Axes:",3(1x,F12.5))

  115 Format(" Deformation  :",F15.10)
  116 Format(" Max Incl     :",F10.4)
  117 Format(" Min Incl     :",F10.4)
  118 Format(" Eff shear st :",F10.4)
  119 Format(" Eff first nsd:",F10.4)
  120 Format(" Eff sec   nsd:",F10.4)

  200 Format(100(1x,f5.3))

  205 Format (/
     +       ,' ndiv   = ',I2,/
     +       ,' phi1   = ',F7.4,/
     +       ,' phi2   = ',F7.4,/
     +       ,' phi3   = ',F7.4,/,/
     +       ,' Iflow  = ',I2,/
     +       ,' wall   = ',F7.4,/,/
     +       ,' eq rad = ',F7.4,/
     +       ,' b/a    = ',F7.4,/
     +       ,' c/a    = ',F7.4,/
     +       ,' cx     = ',F7.4,/
     +       ,' cy     = ',F7.4,/
     +       ,' cz     = ',F7.4,/,/
     +       ,' vs1    = ',F10.5,/
     +       ,' vs2    = ',F10.5,/
     +       ,' tinit  = ',F7.4,/,/
     +       ,' shrt   = ',F7.4,/
     +       )

 208  Format
     +       (' mint   = ',I1,/
     +       ,' NGL     = ',I1,/
     +       ,' DT     = ',F8.6,/
     +       ,' RUNGE  = ',I1,/
     +       ,' MOVE   = ',I1,/
     +       ,' Norm   = ',I1,/
     +       ,' Isym_xy= ',I1,/
     +       )

      Stop
      End
