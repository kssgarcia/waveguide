      program drop_3d

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
c Dynamic simulation of the deformation of a liquid
c drop in an incident simple shear flow
c in the presence of an insoluble surfactant
c
c Three types of flow:
c
c (a) in an infinite domain
c (b) in a semi-infinite domain bounded by a plane wall
c     located at y=wall
c (c) triply periodic flow
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
c  vol         drop volume
c  ars         drop surface area 
c  crx         x position of centroid   
c  cry         y position of centroid   
c  crz         z position of centroid   
c
c  zz, ww         base points and weights of Gauss-Legenrde integr
c  xiq, etq, wq   base points and weights of Gauss-triangle integr
c
c  Ds: surfactant diffusivity
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension      p(1026,3)
      Dimension     ne(1026,7)
      Dimension    vna(1026,3)
      Dimension     u(1026,3)
      Dimension     c(1026),srtn(1026)
      Dimension Umove(1026),Vmove(1026),Wmove(1026)

      Dimension psv(1026,3),Usv(1026),Vsv(1026),Wsv(1026)
      Dimension csv(1026)

      Dimension nvel(1026),lxy(1026,2)

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
      Dimension  sxy(1000),sd1(1000),sd2(1000)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

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

      common/tension/srtn

c---
c for doubly-periodic flow
c---

      common/aaaa/a11,a12,a13,a21,a22,a23,a31,a32,a33
      common/bbbb/b11,b12,b13,b21,b22,b23,b31,b32,b33
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

c---
c various
c---

      common/zwl/zz,ww
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
      None = 1
      Nseven = 7
      Nfour  = 4

      oot = 1.0D0/3.0D0

c------
c input
c------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to read data from file: drop_3d.dat'
      write (6,*) ' 2 to enter the data'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*)   Ienrd

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (2,file="drop_3d.dat")

      read (2,*) Iflow
      read (2,*)
      read (2,*) Ioctaicos
      read (2,*) ndiv 
      read (2,*) 
      read (2,*) req
      read (2,*) boa,coa
      read (2,*) cxp,cyp,czp
      read (2,*) phi1,phi2,phi3
      read (2,*)
      read (2,*) mint
      read (2,*) NGL
      read (2,*)
      read (2,*) wall
      read (2,*)
      read (2,*) a11,a12,a13
      read (2,*) a21,a22,a23
      read (2,*) a31,a32,a33
      read (2,*) Max1,Max2
      read (2,*)
      read (2,*) shrt
      read (2,*)
      read (2,*) vs1
      read (2,*) vs2
      read (2,*)
      read (2,*) tinit
      read (2,*) cinit
      read (2,*) 
      read (2,*) Ds
      read (2,*) 
      read (2,*) Isurf
      read (2,*) betas
      read (2,*) psis
      read (2,*) 
      read (2,*) Ismeth
      read (2,*)
      read (2,*) Nter
      read (2,*) tol
      read (2,*) Idfl
      read (2,*)
      read (2,*) Iread
      read (2,*)
      read (2,*) Norm
      read (2,*) Isym_xy
      read (2,*)
      read (2,*) IRK
      read (2,*) Dt
      read (2,*)
      read (2,*) Move
      read (2,*)
      read (2,*) Nprint_xy
      read (2,*) Nprint_xyz
      read (2,*)

c     write (6,*) ndiv 
c     write (6,*) req
c     write (6,*) boa,coa
c     write (6,*) cxp,cyp,czp
c     write (6,*) phi1,phi2,phi3
c     write (6,*) mint
c     write (6,*) NGL
c     write (6,*) Iflow
c     write (6,*) wall
c     write (6,*) shrt
c     write (6,*) vs1
c     write (6,*) vs2
c     write (6,*) tinit
c     write (6,*) nter
c     write (6,*) tol
c     write (6,*) Idfl
c     write (6,*) Iread
c     write (6,*) Norm
c     write (6,*) Isym_xy
c     write (6,*) IRK
c     write (6,*) Dt
c     write (6,*) Nstep
c     write (6,*) Nprint_xy
c     write (6,*) Nprint_xyz
c     write (6,*) Move

c------------------------
      else if(Ienrd.eq.2) then
c------------------------

      call verbal
     +
     +  (Iflow
     +
     +  ,Ioctaicos
     +  ,ndiv
     +
     +  ,req
     +  ,boa,coa
     +  ,cxp,cyp,czp
     +  ,phi1,phi2,phi3
     +
     +  ,mint
     +  ,NGL
     +
     +  ,wall
     +
     +  ,a11,a12,a13
     +  ,a21,a22,a23
     +  ,a31,a32,a33
     +
     +  ,Max1,Max2
     +
     +  ,shrt
     +  ,vs1,vs2
     +
     +  ,tinit
     +  ,cinit
     +
     +  ,Ds
     +
     +  ,Isurf
     +  ,betas
     +  ,psis
     +
     +  ,Ismeth
     +
     +  ,Nter
     +  ,tol
     +  ,Idfl
     +
     +  ,Iread
     +
     +  ,Norm
     +  ,Isym_xy
     +
     +  ,IRK
     +  ,Dt
     +
     +  ,Nprint_xy
     +  ,Nprint_xyz
     +
     +  ,Move
     +  )

c-----------
      end if               ! End of reading parameters
c-----------

c------------------------
c prepare
c
c surfactant coefficients
c------------------------

      if(Isurf.eq.1) then
         betasc = 1.0D0-betas
      else if(Isurf.eq.2) then
         if(psis.ge.1.0) then
          write (6,*) " drop_3d: psi cannot be greater than 1"
          stop
         end if
         betasc = 1.0D0+betas/psis * Dlog(1.0D0-psis)
      end if

c---------------------------
c triply-periodic shear flow
c---------------------------

      if(Iflow.eq.3) then  ! triply periodic flow

        a11h  = 0.5*a11    ! used for reseting
        a11hm = - a11h     ! the second lattice vector

      end if

c-------------------------------
c set viscosity ratio index, etc
c-------------------------------

      vsrt = vs2/vs1   ! viscosity ratio

c---
c viscosity ratio not equal to unity
c---

      Ivs = 0   ! flag

      if(abs(vsrt-1.0).gt.0.0000001) then
       Ivs = 1
       vsrtm = 1.0D0-vsrt
       vsrtp = 1.0D0+vsrt
       vsf   = 2.0D0/vsrtp
       vsk   = vsrtm/vsrtp
      end if

c--------------------------------
c infinite simple shear flow:
c place the drop center at the origin
c--------------------------------

      if(Iflow.eq.1.or.Iflow.eq.3) then
        cxp = 0.0D0
        cyp = 0.0D0
        czp = 0.0D0
      end if

c-----------------------------
c Read integration quadratures
c-----------------------------

      call gauss_leg (NGL,zz,ww)

      call gauss_trgl 
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

c----------------------------
c triangulate the unit sphere
c
c run even at restart
c to generate the connectivity table
c-----------------------------------

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
      
      write (6,*) " drop_3d: number of points: ",npts
      write (6,*) " drop_3d: number of elements: ",nelm

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
        p(i,1) = p(i,1) + cxp
        p(i,2) = p(i,2) + cyp
        p(i,3) = p(i,3) + czp
      End Do

c---
c unscale to record
c---

      phi1 = phi1/pi
      phi2 = phi2/pi
      phi3 = phi3/pi

c--------------------------------
c assign surfactant concentration
c--------------------------------

      Do i=1,npts
       c(i) = cinit
      End Do

c------------------------
c initialize the velocity
c------------------------

      Do i=1,npts
       u(i,1) = 0.0D0
       u(i,2) = 0.0D0
       u(i,3) = 0.0D0
      End Do

c--------------------
      if(Iread.eq.1) then ! read data from file: drop_3d.inp
c--------------------

c-------------------------------------
c read coordinates of nodes and
c surfactant concentration from file drop_3d.inp
c Set surface tension at nodes
c-------------------------------------

       open (9,file="drop_3d.inp")

       if(Iflow.eq.1.or.Iflow.eq.2) then
         read (9,*) npts,time(1)
       else if(Iflow.eq.3) then
         read (9,*) npts,time(1),a21
       else
         write (6,*) " drop_3d: invalid choice for flow"
         stop
       end if

       Do i=1,npts
         read (9,*) idle,p(i,1),p(i,2),p(i,3),c(i)
       End Do

       close (9)

c-----------
      end if
c-----------

c--------------------------------
c element from node concentration
c by averaging over local nodes
c--------------------------------

      if(Ismeth.eq.1) then

       Do i=1,nelm
        collect = 0.0D0
        Do k=1,6
         j = n(i,k)   ! global label
         collect = collect + c(j)
        End Do
        cel(i) = collect/6.0D0
       End Do

      end if


c------------------------------
c compute surface tension
c from surfactant concentration
c------------------------------

      Do i=1,npts

c---
       if(Isurf.eq.1) then
c---
        srtn(i) = tinit/betasc * (1.0D0-betas*c(i)/cinit)
c---
       else if(Isurf.eq.2) then
c---
        arg = 1.0D0 - psis * c(i)/cinit

        if(arg.lt.0.00001) then
         write (6,*) " drop_3d: error in langmuir"
         Go to 99
        end if
        srtn(i) = tinit/betasc * (1.0D0 + betas/psis * Dlog(arg) )
c---
       else
c---
        srtn(i) = tinit
c---
       end if
c---

      End Do

c----
c compute the surface centroid
c---

      if(Iflow.eq.2) then

       Do k=1,nelm

        i1 = n(k,1)
        i2 = n(k,2)
        i3 = n(k,3)
        i4 = n(k,4)
        i5 = n(k,5)
        i6 = n(k,6)

        call abc 
     +
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
     +  ,alpha(k),beta(k),gamma(k)
     +  )

        End Do

        call elm_geom 
     +
     +   (nelm,npts,mint
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c--------------------
      end if 
c--------------------

c------------------------------------------
c find the points in the z=cz plane
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

c------------------------------------------
c Prepare to exploit symmetry with respect
c to the xy plane
c
c nvelt:  number of nodes where the velocity 
c                           will be computed
c nvelr:  number of nodes where the velocity 
c                          will be reflected 
c
c lxy(i,1): number of reflected node i
c lxy(i,2): its image
c------------------------------------------

c--------------------------
      if(Isym_xy.eq.0) then   ! will compute velocity at all nodes
c--------------------------

       nvelt = npts

       Do i=1,npts
        nvel(i) = i
       End Do

c---------
      else       ! will compute the velocity at roughly half nodes
c---------

      icount1 = 0
      icount2 = 0

      Do 66 i=1,npts

        if(p(i,3).gt.-0.0001) then
          icount1 = icount1+1
          nvel(icount1) = i
        else
          icount2 = icount2+1
          lxy(icount2,1) = i
          Do j=1,npts
           test1 = abs(p(i,1)-p(j,1))
           test2 = abs(p(i,2)-p(j,2))
           test3 = abs(p(i,3)+p(j,3))
           if(    test1.le.0.0001
     +       .and.test2.le.0.0001
     +       .and.test3.le.0.0001
     +       ) then
             lxy(icount2,2) = j
             Go to 66
           End If
          End Do
        end if

  66  Continue

      nvelt = icount1
      nvelr = icount2

c-----------
      end if
c-----------

c---------------
c prepare to run
c---------------

      open (3,file="drop_3d.xy")
      open (4,file="drop_3d.diag")
      open (8,file="drop_3d.xyz")
      open (9,file="drop_3d.rhe")

      Kstep = 1           ! time step counter

      Iprint_xy  = Nprint_xy       ! printing counter
      Iprint_xyz = Nprint_xyz      ! printing counter

c---------------------
c generate a Matlab file
c---------------------

      open (1,file="drop_3d.net")

      Index = 1  ! 6-node triangles
      Index = 2  ! 3-node triangles

      if(Index.eq.1) then          ! 6-node triangles
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      else if(Index.eq.2) then     ! 3-node triangles
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end if

      write (1,*) Iflow
      write (1,*) wall

      write (1,*) None 

      Do k=1,nelm
         call printel (k,Index,c)  ! print in file "drop_3d.net"
      End Do

c-.-.-.-.-.-.-.-.-.-.-.-.-.-.
c  time stepping begins here
c-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 90   Continue

      if(Ienrd.eq.1) then
        read  (2,*)  Nstep
      else if(Ienrd.eq.2) then       ! interactive
        write (6,*)
        write (6,*) 'Enter the number of steps before pausing'
        write (6,*) '0 to quit'
        write (6,*) '---------'
        read  (5,*)  Nstep
      end if

c-----------
c initialize
c-----------

      Istep  = 1      ! cycle step counter

  97  Continue

      write (6,*)
      write (6,*) "--------------------"
      write (6,*)
      write (6,105) Istep,Nstep,time(Kstep)

      Ipass = 0      ! Ipass will be used for volume normalization

  91  Continue

      Ipass = Ipass + 1 

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
c  (a)  surface area of the individual elements
c  (b)  x, y, and z surface moments over each element
c  (c)  total surface area and volume
c  (d)  mean curvature over each element
c  (e)  average normal vector at the nodes
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

c----------------------
c normalize drop volume
c----------------------

      if(Norm.eq.1) then
       if(Ipass.eq.1) then 

         cf = (pi4*req**3/(3.0D0*vlm))**oot

         Do i=1,npts
           p(i,1) = p(i,1)*cf
           p(i,2) = p(i,2)*cf
           p(i,3) = p(i,3)*cf
         End Do

         Go to 91

       end if
      end if

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

c-------------------------------------
c compute the Taylor deformation parameter
c         and the inclination                     
c-------------------------------------

      if(Kstep.eq.1) then ! special treatment of 
                          ! the spherical shape
                          ! commend out if desired
       rmax = 1.0D0
       rmin = 1.0D0
       zmax = 1.0D0

       dxy  (Istep) =  0.00D0
       thmax(Istep) =  0.25D0
       thmin(Istep) =  0.75D0

      else

       call taylor 
     +
     +  (npts,nxy,jxy,cx,cy,cz   ! deformation and
     +  ,rmax,rmin,zmax          ! inclination
     +  ,dxy(Kstep)
     +  ,thmax(Kstep),thmin(Kstep)
     +  )

       call inclination
     +
     +  (nxy,jxy,cx,cy   ! inclination by moment of
     +  ,thmax(Kstep)    ! inertia tensor eigenvalues
     +  ,thmin(Kstep)
     +  )

      end if

      write (6,*)
      write (4,104) Kstep,time(Kstep),rmax,rmin,Dxy(Kstep)
     +             ,thmax(Kstep),thmin(Kstep)
     +             ,zmax

      write (6,113)  time(Kstep),rmax,rmin,zmax
      write (6,115) dxy  (Kstep)
      write (6,116) thmax(Kstep)
      write (6,117) thmin(Kstep)
      write (6,*)

c---------------------------------
c Prepare for the
c triply-periodic Green's function
c---------------------------------

      if(Iflow.eq.3) then

        if(a21.gt.a11h)  a21 = a21-a11
        if(a21.lt.a11hm) a21 = a21+a11

         call sgf_3d_3p_ewald
     +
     +     (a11,a12,a13
     +     ,a21,a22,a23
     +     ,a31,a32,a33
     +     ,b11,b12,b13
     +     ,b21,b22,b23
     +     ,b31,b32,b33
     +     ,ew,tau
     +     )

          call sgf_3d_3p_qqq
     +
     +     (b11,b12,b13
     +     ,b21,b22,b23
     +     ,b31,b32,b33
     +     ,max2
     +     ,ew
     +     )

        if(Ivs.eq.1) then

         call sgf_3d_3p_vvv
     +
     +    (b11,b12,b13
     +    ,b21,b22,b23
     +    ,b31,b32,b33
     +    ,Max2
     +    ,ew
     +    )

        end if

      End If

c------------------------------
c Compute surface tension
c from surfactant concentration
c------------------------------

      Do i=1,npts
       if(Isurf.eq.1) then
        srtn(i) = tinit/betasc * (1.0D0-betas*c(i)/cinit)
       else if(Isurf.eq.2) then
        arg = 1.0D0 - psis * c(i)/cinit
        if(arg.lt.0.00001) then
         write (6,*) " drop_3d: error in langmuir"
         Go to 99
        end if
        srtn(i) = tinit/betasc * (1.0D0 + betas/psis*Dlog(arg) )
       else
        srtn(i) = tinit
       end if
      End Do

c----------------------------------
c compute the velocity at the nodes
c----------------------------------
c
c     write (6,*) 
c     write (6,*) " drop_3d: computing the velocity"
c     write (6,*) 

      call vel
     +
     +  (npts
     +  ,nelm
     +  ,mint
     +  ,NGL
     +  ,Idfl
     +  ,Iflow
     +  ,Isym_xy
     +  ,nter
     +  ,tol
     +  ,Istop
     +  )

      if(Istop.eq.1) Go to 99

c---------------------------------------
c SURFACE KINEMATICS
c
c Compute the rate of surface dilatation
c--------------------------------------

c     If(Iprint_xy.eq.Nprint_xy) then
c     call dilt_3d (nelm,npts,dilt)
c     sum = 0.0D0
c     write (9,*) npts,time(Kstep)
c     Do i=1,npts
c       qqqq = 2.0D0 * vs1 * dilt(i) + 2.0D0 * gamma(i) * crvm(i)
c       write (6,*) i,vs1,dilt(i),gamma(i),crvm(i),qqqq
c       write (9,*) i,qqqq
c       sum = sum+dilt(i)
c     End Do
c     sum = pi4*sum/npts
c     write (6,*) "surface dilatation: ",sum
c     pause
c     End If

c------------------------------------
c SUSPENSION RHEOLOGY:
c
c in the case of infinite shear flow,
c compute the effective stress tensor
c------------------------------------

      if(Iflow.eq.1.or.Iflow.eq.3) then

      call rheology 
     +
     +  (nelm
     +  ,npts
     +  ,mint
     +  ,vlm
     +  ,u
     +  ,sxy(Kstep)
     +  ,sd1(Kstep),sd2(Kstep)
     +  )

      if(Iflow.eq.1) then
       cf = 1.0D0/(vlm*shrt)
      else if(Iflow.eq.3) then
       cf = 1.0D0/shrt
      end if

      sxy(Kstep) = sxy(Kstep)*cf
      sd1(Kstep) = sd1(Kstep)*cf
      sd2(Kstep) = sd2(Kstep)*cf

      write (6,*)
      write (6,118) sxy(Kstep)
      write (6,119) sd1(Kstep)
      write (6,120) sd2(Kstep)

      write (6,*) "--------------------"

      end if

c--------------------
c print cross-section
c--------------------

      if(Isym_xy.eq.1) then

c       write (6,*)
c       write (6,*) " drop_3d: x, y, ux, uy, uz, c at t=",time(Kstep)
c       write (6,*)

c       Do i=1,nxy
c         m = jxy(i)
c         write (6,100) i,(p(m,j),j=1,2),(u(m,j),j=1,3),c(m)
c       End Do

c-----

        if(Iprint_xy.eq.Nprint_xy) then

          if(Iflow.eq.3) then
           write (3,100) nxy,time(Kstep),a21
          else
           write (3,100) nxy,time(Kstep)
          end if

          Do i=1,nxy
           m = jxy(i)
           write (3,100) i,p(m,1),p(m,2),c(m)
          End Do

          Iprint_xy = 0

        end if
c----

      end if

c----------------
c print all nodes
c----------------

      if(Iprint_xyz.eq.Nprint_xyz) then

         if(Iflow.eq.3) then
          write (8,100) npts,time(Kstep),a21
         else
          write (8,100) npts,time(Kstep)
         end If

         Do i=1,npts
           write (8,108) i,p(i,1),p(i,2),p(i,3),c(i)
         End Do

        ! print Matlab file (file 1)

        write (1,*) None 

        Do k=1,nelm
          call printel (k,Index,c)  ! print in file "drop_3d.net"
        End Do

        Iprint_xyz = 0

      end if

c-----------------
c time integration
c-----------------

      if(Nstep.eq.0) Go to 99

c-----------------------
c advance the surfactant
c-----------------------

      if(betas.gt.0.000001) then

      if(Ismeth.eq.1) then ! finite volume method

      write (6,*) " drop_3d: entering surfact_fvm"

      call surfact_fvm
     +
     +  (nelm
     +  ,npts
     +  ,mint
     +  ,Dt
     +  ,Ds
     +  ,cel
     +  ,Move
     +  )

      else if(Ismeth.eq.2) then ! finite element method

      write (6,*) " drop_3d: entering surfact_fem"

      call surfact_fem
     +
     +  (nelm
     +  ,npts
     +  ,mint
     +  ,Dt
     +  ,Ds
     +  ,Move
     +  )

      end if

      write (6,*) " drop_3d: exited"

      end if

c----------------------
      if(IRK.eq.2) then    ! RK2: save position

       Do i=1,npts
        psv(i,1) = p(i,1)  ! save position
        psv(i,2) = p(i,2)
        psv(i,3) = p(i,3)
       End Do

      end if
c----------------------

      Do i=1,npts

        if(Move.eq.0) then      ! points move with total velocity

          Umove(i) = u(i,1)
          Vmove(i) = u(i,2)
          Wmove(i) = u(i,3)

        else if(Move.eq.1) then  ! points move with normal velocity

          Uvel  = u(i,1)*vna(i,1)   ! projection of velocity
     +          + u(i,2)*vna(i,2)   ! onto the normal vector
     +          + u(i,3)*vna(i,3)

          Umove(i) = Uvel*vna(i,1)
          Vmove(i) = Uvel*vna(i,2)
          Wmove(i) = Uvel*vna(i,3)

        end if

c---
c  flow above a wall
c---

        if(Iflow.eq.2) then            ! add tangential component
          if(Move.eq.1) then           ! of centroid velocity

           Cntr_vel = shrt*(cy-wall)
           Umove(i) = Umove(i) + (1.0-vna(i,1)*vna(i,1))*Cntr_vel
           Vmove(i) = Vmove(i) + (   -vna(i,1)*vna(i,2))*Cntr_vel
           Wmove(i) = Wmove(i) + (   -vna(i,1)*vna(i,3))*Cntr_vel

          end if
        end if

c---
c advance in time
c---

        p(i,1) = p(i,1) + Dt*Umove(i)
        p(i,2) = p(i,2) + Dt*Vmove(i)
        p(i,3) = p(i,3) + Dt*Wmove(i)

c       write (6,100) i,(u(i,j),j=1,3) 

      End Do

c---
c enforce symmetry if required
c---

      if(Isym_xy.eq.1) then
        Do i=1,nxy
          p(jxy(i),3) = 0.0D0
        End Do
      End if

      if(Iflow.eq.3) a21 = a21 + Dt*shrt*a22    ! Triply-periodic flow

c-------------
c end of RK1
c-------------

      if(IRK.eq.1) Go to 87

c---------
c RK2 step
c---------

      Do i=1,npts          ! save velocity
        Usv(i) = Umove(i)
        Vsv(i) = Vmove(i)
        Wsv(i) = Wmove(i)
      End Do

c---
c Geometry at intermediate step
c---

      Do k=1,nelm

        i1 = n(k,1)
        i2 = n(k,2)
        i3 = n(k,3)
        i4 = n(k,4)
        i5 = n(k,5)
        i6 = n(k,6)

        call abc 
     +
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
     +  ,alpha(k),beta(k),gamma(k)
     +  )

      End Do

      call elm_geom 
     +
     +  (nelm,npts,mint
     +  ,area,vlm
     +  ,cx,cy,cz
     +  )

c----
c Prepare again for the
c triply-periodic Green's function
c-----

      if(Iflow.eq.3) then

        if(a21.gt.a11h)  a21 = a21-a11
        if(a21.lt.a11hm) a21 = a21+a11

         call sgf_3d_3p_ewald
     +
     +    (a11,a12,a13
     +    ,a21,a22,a23
     +    ,a31,a32,a33
     +    ,b11,b12,b13
     +    ,b21,b22,b23
     +    ,b31,b32,b33
     +    ,ew,tau
     +    )

          call sgf_3d_3p_qqq
     +
     +     (b11,b12,b13
     +     ,b21,b22,b23
     +     ,b31,b32,b33
     +     ,Max2
     +     ,ew
     +     )

        if(Ivs.eq.1) then

         call sgf_3d_3p_vvv
     +
     +    (b11,b12,b13
     +    ,b21,b22,b23
     +    ,b31,b32,b33
     +    ,Max2
     +    ,ew
     +    )

        end if

      end if

c------------------------------
c Compute the surface tension
c from surfactant concentration
c using a surface equation of state
c------------------------------

      Do i=1,npts

       if(Isurf.eq.1) then
        srtn(i) = tinit/betasc * (1.0D0-betas*c(i)/cinit)
       else if(Isurf.eq.2) then
        arg = 1.0D0 - psis * c(i)/cinit
        if(arg.lt.0.00001) then
         write (6,*) " drop_3d: error in langmuir"
         Go to 99
        end if
        srtn(i) = tinit/betasc * (1.0D0 + betas/psis * Dlog(arg) )
       else
        srtn(i) = tinit
       end if

      End Do
 
c---
c second velocity evaluation
c---

      call vel
     +
     +   (npts,nelm
     +   ,mint,NGL
     +   ,Idfl
     +   ,Iflow
     +   ,Isym_xy
     +   ,nter,tol
     +   ,Istop
     +   )

      if(Istop.eq.1) Go to 99

c---
c second step in RK2
c
c see previous comments for an explanation 
c of the individual steps
c---

      Dth = 0.5D0*Dt

      Do i=1,npts

        if(Move.eq.0) then

          Umove(i) = u(i,1)
          Vmove(i) = u(i,2)
          Wmove(i) = u(i,3)

        else if(Move.eq.1)  then

          Uvel  = u(i,1)*vna(i,1)
     +          + u(i,2)*vna(i,2)
     +          + u(i,3)*vna(i,3)
          Umove(i) = Uvel*vna(i,1)
          Vmove(i) = Uvel*vna(i,2)
          Wmove(i) = Uvel*vna(i,3)

        end if

        if(Iflow.eq.2) then

          if(Move.eq.1) then 

           Cntr_vel = shrt*(cy-wall)

           Umove(i) = Umove(i) + (1.0D0-vna(i,1)*vna(i,1))*Cntr_vel
           Vmove(i) = Vmove(i) + (     -vna(i,1)*vna(i,2))*Cntr_vel
           Wmove(i) = Wmove(i) + (     -vna(i,1)*vna(i,3))*Cntr_vel

          end if
        end if

        p(i,1) = psv(i,1) + Dth*(Usv(i)+Umove(i))
        p(i,2) = psv(i,2) + Dth*(Vsv(i)+Vmove(i))
        p(i,3) = psv(i,3) + Dth*(Wsv(i)+Wmove(i))

      End Do

c-----------
c end of RK2
c-----------

  87  Continue

c-------------------
c end of a time step
c-------------------

c------------------------------
c enforce symmetry if required
c------------------------------

      if(Isym_xy.eq.1) then
        Do i=1,nxy
         p(jxy(i),3) = 0.0D0
        End Do
      end if

c------------------------
c reset counters and time
c------------------------

      Kstep = Kstep+1
      Istep = Istep+1

      Iprint_xy  = Iprint_xy+1
      Iprint_xyz = Iprint_xyz+1

      time(Kstep) = time(Kstep-1) + Dt

      if(Istep.le.Nstep) Go to 97

      Go to 90     ! return for another step

c-.-.-.-.-.-.-.-.-.-.-
c  Simulation has ended
c-.-.-.-.-.-.-.-.-.-.-

  99  Continue

c-------------------
c record final shape
c-------------------

      write (8,*) npts,time(Kstep)

      Do i=1,npts
        write (8,108) i,p(i,1),p(i,2),p(i,3),c(i)
      End Do

      write (8,100) null
      write (3,100) null

c-------------------
c record diagnostics
c-------------------

      write (4,*) Kstep," time, surf_ar, vol, cx, cy, cz"

      Do i=1,Kstep
        write (3,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (4,104) i,time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (6,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (8,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
      End Do

      write (3,*)
      write (4,*) null
      write (4,*) Kstep," time, D, thmax, thmin"
      write (6,*)
      write (8,*)

      Do i=1,Kstep
        write (3,107)   time(i),dxy(i),thmax(i),thmin(i)
        write (4,104) i,time(i),dxy(i),thmax(i),thmin(i)
        write (6,107)   time(i),dxy(i),thmax(i),thmin(i)
        write (8,107)   time(i),dxy(i),thmax(i),thmin(i)
      End Do

      write (4,*) null

c---
c record rheology
c---

      if(Iflow.eq.1) then

        write (3,*)
        write (4,*) Kstep," time, sxy,sd1, sd2"
        write (6,*)
        write (8,*)

        Do i=1,Kstep
          write (3,109)   time(i),sxy(i),sd1(i),sd2(i)
          write (4,104) i,time(i),sxy(i),sd1(i),sd2(i)
          write (6,109)   time(i),sxy(i),sd1(i),sd2(i)
          write (8,109)   time(i),sxy(i),sd1(i),sd2(i)
        End Do
        write (4,*) null

      End If

      write (9,*) null,null,null

c---
c record run parameters
c---


      write (3,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,shrt
     +             ,req,boa,coa,cxp,cyp,czp
     +             ,vs1,vs2,tinit,cinit
     +             ,Ds,Isurf,betas,psis,Ismeth

      write (4,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,shrt
     +             ,req,boa,coa,cxp,cyp,czp
     +             ,vs1,vs2,tinit,cinit
     +             ,Ds,Isurf,betas,psis,Ismeth

      write (6,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,shrt
     +             ,req,boa,coa,cxp,cyp,czp
     +             ,vs1,vs2,tinit,cinit
     +             ,Ds,Isurf,betas,psis,Ismeth

      write (8,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,shrt
     +             ,req,boa,coa,cxp,cyp,czp
     +             ,vs1,vs2,tinit,cinit
     +             ,Ds,Isurf,betas,psis,Ismeth

      write (3,208) mint,NGL,Dt,Irk,Move,Norm,Isym_xy
      write (4,208) mint,NGL,Dt,Irk,Move,Norm,Isym_xy
      write (6,208) mint,NGL,Dt,Irk,Move,Norm,Isym_xy
      write (8,208) mint,NGL,Dt,Irk,Move,Norm,Isym_xy

c--------------------------
c print Matlab file (file 1)
c--------------------------

        write (1,*) None 

        Do k=1,nelm
          call printel (k,Index,c)  ! print in file "drop_3d.net"
        End Do

        write (1,*) Null

c------------
c close files
c------------

      close (1)
      close (2)
      close (3)
      close (4)
      close (9)

c-------------
c restart data
c-------------

      open (11,file="drop_3d.rst")

       write (11,*) npts,time(Kstep)

       Do i=1,npts
         write (11,108) i,p(i,1),p(i,2),p(i,3),c(i)
       End Do
       write (11,150) Null,Null

      close (11)

c-----
c done
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

  150 Format (1x,i4,1x,i4)

  200 Format(100(1x,f5.3))

  205 Format (/
     +       ,' ndiv   = ',I2,/
     +       ,' phi1   = ',F7.4,/
     +       ,' phi2   = ',F7.4,/
     +       ,' phi3   = ',F7.4,/,/
     +       ,' Iflow  = ',I2,/
     +       ,' wall   = ',F7.4,/,/
     +       ,' shrt   = ',F7.4,/
     +       ,' eq rad = ',F7.4,/
     +       ,' b/a    = ',F7.4,/
     +       ,' c/a    = ',F7.4,/
     +       ,' cx     = ',F7.4,/
     +       ,' cy     = ',F7.4,/
     +       ,' cz     = ',F7.4,/,/
     +       ,' vs1    = ',F10.5,/
     +       ,' vs2    = ',F10.5,/
     +       ,' tinit  = ',F7.4,/,/
     +       ,' cinit  = ',F7.4,/
     +       ,' Ds     = ',F7.4,/
     +       ,' Isurf  = ',I1,/
     +       ,' betas  = ',F7.4,/
     +       ,' psis   = ',F7.4,/,/
     +       ,' Ismeth = ',I1,/,/
     +       )

 208  Format
     +       (' mint   = ',I1,/
     +       ,' NGL    = ',I1,/
     +       ,' Dt     = ',F8.6,/
     +       ,' IRK    = ',I1,/
     +       ,' Move   = ',I1,/
     +       ,' Norm   = ',I1,/
     +       ,' Isym_xy= ',I1,/
     +       )

      Stop
      End
