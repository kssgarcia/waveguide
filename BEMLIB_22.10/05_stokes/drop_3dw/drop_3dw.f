      program drop_3dw

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c----------------------------------------------------
c Dynamic simulation of the deformation of a liquid
c drop adhering to a plane wall in simple shear flow
c
c The wall is located at y = 0
c
c In this implementation, the surface tension is constant.
c
c Surfactant concentration is assigned to the nodes, but is
c dynamically inactive; that is, it remains constant in time
c
c SYMBOLS:
c -------
c
c  npts	    total number of points
c  nelm	    total number of elements
c
c  p(i,j)   coordinates of nodes i (j=1,2,3)
c
c  ne(k,j)  ne(k,1) is the number of elements adjacent to point k
c           ne(k,2), ... are the elements numbers, j = 2, ..., 7
c           for this triangulation, up to six
c
c  n(k,i)   connectivity table: points for element k, i = 1,...,6
c
c  nbe(k,j) the three neighboring elements of element k (j=1,2,3)
c
c  vna      average value of the normal vector at nodes
c  u        node velocity
c  crvm(i)  mean curvature at node i   
c  srtn(i)  isotropic surface tension at node i   
c  c(i)     surfactant concentration at node i   
c
c  psv             save p
c  Usv, Vsv, Wsv   save U,V,W
c  csv             save c
c
c    arel(k)      surface area of element k
c  crvmel(k)    average value of the mean curvature of element k
c
c    xmom(k)      x-moment of element k
c    xmom(k)      y-moment of element k
c    zmom(k)      y-moment of element k
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
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension     p(1026,3)
      Dimension    ne(1026,7)
      Dimension   vna(1026,3)
      Dimension     u(1026,3)
      Dimension     c(1026),srtn(1026)
      Dimension Umove(1026),Vmove(1026),Wmove(1026)

      Dimension psv(1026,3),Usv(1026),Vsv(1026),Wsv(1026)

      Dimension nvel(1026),lxy(1026,2)

      Dimension crvm(1026)

      Dimension       n(512,6),nbe(512,3)
      Dimension   alpha(512), beta(512), gamma(512)
      Dimension    arel(512), xmom(512),  ymom(512),zmom(512)
      Dimension  crvmel(512)

      Dimension jxy(100)

      Dimension time(1000),dxy(1000),thmax(1000),thmin(1000)
      Dimension  ars(1000),vol(1000),crx(1000),cry(1000),crz(1000)


      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)
      Dimension xid(30),etd(30),wd(30)

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

      common/viscr/vs1,vs2,vsrt,vsrtp,vsrtm,vsf,vsk
      common/visci/Ivs

      common/var/shrt,wall

      common/veloc0/u
      common/veloc1/nvelt,nvel
      common/veloc2/nvelr,lxy

      common/tension/srtn

c---
c other
c---

      common/base/rbase
      common/dskr/xid,etd,wd
      common/dski/mdsk

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

      Null  = 0
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

      open (2,file="drop_3dw.dat")

      read (2,*) ndiv 
      read (2,*)
      read (2,*) req
      read (2,*) cont_angle
      read (2,*)
      read (2,*) mint
      read (2,*) NGL
      read (2,*) mdsk
      read (2,*)
      read (2,*) shrt
      read (2,*)
      read (2,*) vs1
      read (2,*) vs2
      read (2,*)
      read (2,*) tinit
      read (2,*)
      read (2,*) nter
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
      read (2,*) Nprint_xy
      read (2,*) Nprint_xyz
      read (2,*)
      read (2,*) Move
      read (2,*)

c------------------------
      elseif(Ienrd.eq.2) then
c------------------------

 88   write (6,*)
      write (6,*)
      write (6,*) " Please enter the level of triangulation"
      write (6,*) " Choose from: 0, 1, 2, 3"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) ndiv 

      if(ndiv.eq.99) Go to 99

      if(ndiv.gt.3) then
        write (6,*) 'Level is too high, should be less than 4'
        Go to 88
      end if

      write (6,*)
      write (6,*) " The volume of the drop will be equal to"
      write (6,*) " VD =  2*pi/3 * a**3,"
      write (6,*) " where a is the equivalent radius"
      write (6,*)
      write (6,*) " Please enter the drop equivalent radius"
      write (6,*) " ---------------------------------------"
      read  (5,*) req

      write (6,*)
      write (6,*) " The drop has a spherical initial shape"
      write (6,*)
      write (6,*) " Please enter the contact angle"
      write (6,*) "             in multiples of pi"
      write (6,*) " ------------------------------"
      read  (5,*) cont_angle

      write (6,*)
      write (6,*) " Non-singular integration over each triangle"
      write (6,*)
      write (6,*) " Will use the m-point rule."
      write (6,*)
      write (6,*) " Please enter m"
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*) "  0 to quit"
      write (6,*) " ---------------------------------"
      read  (5,*) mint

      if(mint.eq.0) Go to 99

      write (6,*)
      write (6,*) " Singular integration over each triangle"
      write (6,*) " Will use the m-point Gauss-Legendre quadrature"
      write (6,*) " Please enter m"
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) NGL

      if(NGL.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Integration over the diskoidal drop base"
      write (6,*) " for the desingularization of the"
      write (6,*) " double-layer potential"
      write (6,*) " Will use the m-point rule."
      write (6,*) " Please enter m; choose from 4, 5, 21"
      write (6,*) "  0 to quit"
      write (6,*) " ----------"
      read  (5,*) mdsk

      if(mdsk.eq.0) Go to 99

      write (6,*)
      write (6,*) " Simple shear flow"
      write (6,*)
      write (6,*) " Enter the shear rate" 
      write (6,*) " --------------------"
      read  (5,*) shrt

      write (6,*)
      write (6,*) " Enter the viscosity of the ambient fluid"
      write (6,*) " ----------------------------------------"
      read  (5,*) vs1

      write (6,*)
      write (6,*) " Enter the viscosity of the drop"
      write (6,*) " -------------------------------"
      read  (5,*) vs2

      write (6,*)
      write (6,*) " Enter the initial surface tension"
      write (6,*) " ---------------------------------"
      read  (5,*) tinit

      write (6,*)
      write (6,*) " Enter the number of max iterations "
      write (6,*) "       in solving the integral equation"
      write (6,*) " --------------------------------------"
      read  (5,*) nter

      write (6,*)
      write (6,*) " Enter the error tolerance"
      write (6,*) "       in solving the integral equation"
      write (6,*) " --------------------------------------"
      read  (5,*) tol

      write (6,*) " Enable deflation ? "
      write (6,*)
      write (6,*) " Enter 0 for no? "
      write (6,*) "       1 for deflation of one eigenvalue "
      write (6,*) " ----------------------------------------"
      read  (5,*) Idfl

      write (6,*) " Enter 0 to generate initial data"
      write (6,*) "       1 to read from file: drop_3dw.inp "
      write (6,*) "         (restart) "
      write (6,*) " ---------------------------------------"
      read  (5,*) Iread

      write (6,*)
      write (6,*) " Normalize drop volume after each step ?"
      write (6,*) " 0 for no, 1 for yes"
      write (6,*) " -------------------"
      read  (5,*) Norm

      write (6,*)
      write (6,*) " Exploit symmetry wr to the xy plane ?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      read  (5,*) Isym_xy

      write (6,*)
      write (6,*) " Time integration"
      write (6,*)
      write (6,*) " Enter 1 for the Euler explicit method"
      write (6,*) "       2 for the RK2 method"
      write (6,*) " -------------------------------------"
      read  (5,*) IRK

      write (6,*)
      write (6,*) " Choose the time step Dt"
      write (6,*) " -----------------------"
      read  (5,*) Dt

      write (6,*) " Will record a profile in the xy plane"
      write (6,*) "      after N steps; please enter N"
      write (6,*) " -------------------------------------"
      read  (5,*) Nprint_xy

      write (6,*)
      write (6,*) " Will record a 3D profile"
      write (6,*) "      after N steps; please enter N"
      write (6,*) " ----------------------------------"
      read  (5,*) Nprint_xyz

      write (6,*)
      write (6,*) ' Select point advancement policy '
      write (6,*)
      write (6,*) ' Enter 0 to move with the total velocity '
      write (6,*) '       1 to move with the normal velocity '
      write (6,*) ' -----------------------------------------'
      read  (5,*) Move

c---------
      end if
c---------

c---------
c uncomment to confirm
c---
c
c     write (6,*) ndiv 
c     write (6,*) req
c     write (6,*) cont_angle
c     write (6,*) mint
c     write (6,*) NGL
c     write (6,*) mdsk
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
c     write (6,*) Nprint_xy
c     write (6,*) Nprint_xyz
c     write (6,*) Move

c---------------------------------------
c preparations, adjustments, definitions
c---------------------------------------

c---
c set viscosity ratio index, etc
c---

      vsrt = vs2/vs1   ! viscosity ratio

c---
c viscosity ratio not equal to unity
c---

      Ivs = 0

      if(abs(vsrt-1.0).gt.0.0000001) then
       Ivs   = 1
       vsrtm = 1.0D0-vsrt
       vsrtp = 1.0D0+vsrt
       vsf   = 2.0D0/vsrtp
       vsk   = vsrtm/vsrtp
      end if

c------------------------
c integration quadratures
c------------------------

      call gauss_leg (NGL,zz,ww)      ! for singular slp

      call gauss_trgl
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

      call gauss_disk (mdsk,xid,etd,wd)     ! for dlp

c----------------------------------
c triangulate a hemisphere
c
c run even at restart
c to generate the connectivity table
c-----------------------------------

      call trgl6_hsph_octa 
     +
     +  (ndiv
     +  ,npts,nelm
     +  )
      
      write (6,*)
      write (6,*) "Number of points:   ",npts
      write (6,*)
      write (6,*) "Number of elements: ",nelm
      write (6,*)

      nelmh = nelm/2

c--------------
c initial state   ! SECTION OF A SPHERE
c--------------

      time(1) = 0.0D0

c-----------------------------------------------
c slide the hemisphere to conform with the chosen
c contact angle
c
c shift the shape so that the contact line is at y = 0
c
c Scale so that
c the volume of the spherical section
c is identical to that of the unit hemishpere
c-----------------------------------------------

      phiu  = cont_angle*pi
      cs    = cos(phiu)
      scale = 2.0D0/3.0D0 * 1.0D0/(1.0D0-cs-(1.0-cs**3)/3.0D0)
      scale = req*scale**oot

      shift = sin(pih-phiu)

      Do i=1,npts

       rr = sqrt(p(i,1)**2+p(i,3)**2)

       thet0 = atan2(rr,    p(i,2))
       thetz = atan2(p(i,1),p(i,3))

       thet = thet0*phiu/pih

       rrr =  sin(thet)

       p(i,1) = rrr*sin(thetz)
       p(i,2) =     cos(thet) 
       p(i,3) = rrr*cos(thetz)

       p(i,2) = p(i,2) - shift

       p(i,1) = p(i,1)*scale
       p(i,2) = p(i,2)*scale
       p(i,3) = p(i,3)*scale

      End Do

c-----------------------------------------
c assign surfactant concentration to nodes
c
c compute surface tension at nodes
c
c in this implementation, the surfactant
c concentration remains constant
c-----------------------------------------

      Do i=1,npts
          c(i) = 1.0D0
       srtn(i) = tinit
      End Do

c------------------------
c Initialize the velocity
c------------------------

      Do i=1,npts
       u(i,1) = 0.0D0
       u(i,2) = 0.0D0
       u(i,3) = 0.0D0
      End Do

c------------------------
      if(Iread.eq.1) then
c------------------------

c-------------------------------------
c read coordinates of nodes and
c surfactant concentration from file: drop_3dw.inp
c
c set surface tension at nodes
c-------------------------------------

      open (9,file="drop_3dw.inp")

      read (9,*) npts,time(1)

      Do i=1,npts
         write (1,100) idle,p(i,1),p(i,2),p(i,3),c(i)
         srtn(i) = tinit
      End Do

      close (9)

c-----------
      end if         ! END OF READING
c-----------

c---
c set the drop base center
c---

      cxbase = 0.0D0
      cybase = 0.0D0
      czbase = 0.0D0

c------------------------------------
c find the radius of the contact line
c------------------------------------

      Do i=1,npts
       if(abs(p(i,2)).lt.0.0000001) then
        rbase = sqrt(p(i,1)**2+p(i,3)**2)
        Go to 93
       end if
      End Do

 93   Continue

      write (6,*)
      write (6,*) " Base radius = ",rbase
      write (6,*)

c--------------------------------------------
c find the points in the xy slice,
c and order them
c
c nxy is number of points in the xy plane
c the vector jxy contains consecutive points
c-------------------------------------------

      call xy_slice
     +
     +  (npts
     +  ,cxbase,cybase,czbase
     +  ,nxy
     +  ,jxy
     +  )

      Do i=1,nxy
        nid  = jxy(i)
        write(6,102) p(nid,1),p(nid,2),p(nid,3)
      End Do

c----------------------------------------
c prepare to exploit symmetry with respect
c to the xy plane
c
c nvelt:  number of nodes where velocity 
c                       will be computed
c nvelr:  number of nodes where velocity 
c                      will be reflected 
c
c lxy(i,1): number of reflected node i
c lxy(i,2): its image
c------------------------------------------

c---------
      if(Isym_xy.eq.0) then
c---------

c----
c  will compute velocity at all nodes
c-----

       nvelt = npts

       Do i=1,npts
        nvel(i) = i
       End Do

c---------
      else
c---------

      icount1 = 0
      icount2 = 0

      Do 66 i=1,npts

        if(p(i,3).gt.-0.00001) then
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
           end if
          end Do
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

      open (3,file="drop_3dw.xy")
      open (4,file="drop_3dw.diag")
      open (8,file="drop_3dw.xyz")

      Kstep = 1               ! global time step counter

      Iprint_xy  = Nprint_xy
      Iprint_xyz = Nprint_xyz

c-----------------------
c generate a Matlab file
c-----------------------

      open (1,file="drop_3dw.net")

      nelmprint = nelmh
      nelmprint = nelm

      Index = 1  ! 6-node triangles
      Index = 2  ! 3-node triangles

      if(Index.eq.1) then          ! 6-node triangles
        write (1,*) Nseven
        write (1,*) 7*nelmprint
        write (1,*) nelmprint
      else if(Index.eq.2) then     ! 3-node triangles
        write (1,*) Nfour
        write (1,*) 16*nelmprint
        write (1,*) 4*nelmprint
      end if

      write (1,*) None

      Do k=1,nelmprint
        call printel(k,Index,c)
      End Do

      ! points around the mid-section

      write (1,*) nxy

      Do i=1,nxy
        nid  = jxy(i)
        write (1,102) p(nid,1),p(nid,2),p(nid,3)
      End Do

c-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
c=================================c
c   Time stepping begins here     c
c=================================c
c-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 90   Continue

      if(Ienrd.eq.1) then
        read (2,*)  Nstep
      else
        write (6,*)
        write (6,*) 'Enter the number of steps before pausing'
        write (6,*) '0 to quit'
        write (6,*) '---------'
        read  (5,*)  Nstep
      end if

c-----------
c initialize
c-----------

      Istep  = 1      ! batch step counter

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
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +    ,alpha(k),beta(k),gamma(k)
     +    )

      End Do

c------------------------------------------------
c compute:
c         surface area of the individual elements
c         x, y, and z moments over each element
c         total particle surface area and volume
c         mean curvature of each element
c         average normal vector
c------------------------------------------------

      call elm_geom
     +
     +   (nelm,npts,mint
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c---------------
c     write (6,*)
c     write (6,*) " Curvature of elements"
c     write (6,*)
c     Do k = 1,nelm
c       write (6,100) k,crvmel(k)
c     End Do
c---------------


c------------------------
c normalize drop volume ?
c------------------------

      if(Norm.eq.1) then
       if(Ipass.eq.1) then

         cf = (pi2*req**3/(3.0D0*vlm))**oot

         Do i=1,npts
           p(i,1) = p(i,1)*cf
           p(i,2) = p(i,2)*cf
           p(i,3) = p(i,3)*cf
         End Do

         Go to 91

       end if
      end if

c----------
c recording
c----------

      ars(Kstep) = area/(pi2*req**2)
      vol(Kstep) = vlm/ (pi2*req**3/3.0D0)

      crx(Kstep) = cx   ! drop volume centroid
      cry(Kstep) = cy
      crz(Kstep) = cz

      write (6,*)
      write (6,110) ars(Kstep)
      write (6,111) vol(Kstep)
      write (6,112) cx,cy,cz
      write (6,*)

c---------------
c     write (6,*) 
c     write (6,*) " Curvature of elements"
c     write (6,*) 
c     Do k = 1,nelm
c       write (6,100) k,crvm(k)
c     End Do
c---------------

c-------------------------------------
c Compute:
c         Taylor deformation parameter
c         drop inclination
c-------------------------------------

      if(Istep.eq.1) then      ! special treatment of
                               ! the spherical shape
                               ! commend out if desired
       rmax = 1.0
       rmin = 1.0
       zmax = 1.0
         dxy(Kstep) =   0.0
       thmax(Kstep) =  45.0
       thmin(Kstep) = 135.0

      else

       call taylor 
     +
     +    (npts,nxy,jxy,cx,cy,cz
     +    ,rmax,rmin,zmax
     +    ,dxy(Kstep)
     +    ,thmax(Kstep),thmin(Kstep)
     +    )

      end if

      write (6,*)
      write (6,115)   dxy(Kstep)
      write (6,116) thmax(Kstep)
      write (6,117) thmin(Kstep)
      write (6,*)

c----------------------------------
c compute the velocity at the nodes
c----------------------------------
c
c     write (6,*) 
c     write (6,*) " drop_3dw: computing the velocity"
c     write (6,*) 

      call vel
     +
     +   (npts
     +   ,nelm
     +   ,mint
     +   ,NGL
     +   ,Idfl
     +   ,Isym_xy
     +   ,nter
     +   ,tol
     +   ,Istop
     +   )

c      Do i=1,npts
c        write (6,100) i,(p(i,j),j=1,3),(u(i,j),j=1,3)
c      End Do

      if(Istop.eq.1) Go to 99

c--------------------
c print cross-section
c--------------------

      if(Isym_xy.eq.1) then

c       write (6,*)
c       write (6,*)  "position and velocity at time:",time(Istep)
c       write (6,*)

c       Do i=1,nxy
c         m = jxy(i)
c         write (6,100) i,(p(m,j),j=1,2),(u(m,j),j=1,3)
c       End Do

c-----

        if(Iprint_xy.eq.Nprint_xy) then

          Iprint_xy = 0

          write (3,100) nxy,time(Kstep)

          Do i=1,nxy
            m = jxy(i)
            write (3,100) i,(p(m,j),j=1,2)
          End Do

        end if

      end if

c----------------
c print all nodes
c----------------

      if(Iprint_xyz.eq.Nprint_xyz) then

         write (8,100) npts,time(Istep)

         Do i=1,npts
           write (8,108) i,p(i,1),p(i,2),p(i,3),c(i)
         End Do

         write (1,*) None

         Do k=1,nelmprint
           call printel(k,Index,c)
         End Do

      ! points around the mid-section

         write (1,*) nxy

         Do i=1,nxy
           nid  = jxy(i)
           write (1,102) p(nid,1),p(nid,2),p(nid,3)
         End Do

         Iprint_xyz = 0

      end if

c------------------------
c Runge-Kutta integration
c------------------------

      if(Nstep.eq.0) Go to 99

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
        Do i = 1,nxy
          p(jxy(i),3) = 0.0D0
        End Do
      end if

c-------------
c end of RK1
c-------------

      if(IRK.eq.1) Go to 87

c----------------------
c RK2 step
c----------------------

      Do i=1,npts
        Usv(i) = Umove(i)
        Vsv(i) = Vmove(i)
        Wsv(i) = Wmove(i)
      End Do

c---
c geometry at intermediate step
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
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +    ,alpha(k),beta(k),gamma(k)
     +    )

      End Do

      call elm_geom 
     +
     +    (nelm,npts,mint
     +    ,area,vlm
     +    ,cx,cy,cz
     +    )

c---
c second velocity evaluation
c---

      call vel
     +
     +   (npts
     +   ,nelm
     +   ,mint,NGL
     +   ,Idfl
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

      Dth = 0.5*Dt

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

        p(i,1) = psv(i,1) + Dth*(Usv(i)+Umove(i))
        p(i,2) = psv(i,2) + Dth*(Vsv(i)+Vmove(i))
        p(i,3) = psv(i,3) + Dth*(Wsv(i)+Wmove(i))

      End Do

c-----------
c end of RK2
c-----------

  87  Continue

c--------
c end of a time step
c--------

c-----------------------------
c enforce symmetry if required
c-----------------------------

      if(Isym_xy.eq.1) then
        Do i=1,nxy
         p(jxy(i),3) = 0.0D0
        End Do
      end if

c------------------------
c reset counters and time
c------------------------

      Kstep = Kstep + 1
      Istep = Istep + 1
      Iprint = Iprint + 1

      Iprint_xy  = Iprint_xy+1
      Iprint_xyz = Iprint_xyz+1

      time(Kstep) = time(Kstep-1) + Dt

      If(Istep.le.Nstep) Go to 97

      Go to 90       ! Return for another step

c-.-.-.-.-.-.-.-.--.-.-.-.-.-.-.-.-.
c----------------------------------c
c  simulation has ended            c
c----------------------------------c
c-.-.-.-.-.-.-.-.--.-.-.-.-.-.-.-.-.

  99  Continue

c------------------
c record final shape
c------------------

      write (8,*) npts,time(Kstep)," restart data"

      Do i=1,npts
        write (8,100) i,p(i,1),p(i,2),p(i,3),c(i)
      End Do

      write (3,100) null
      write (8,100) null

c------------------
c record diagnostics
c------------------

      write (4,*) Kstep," time, surf_ar, vol, cx, cy, cz"

      Do i=1,Kstep
        write (3,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (4,103) i,time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (6,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (8,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
      End Do

      write (3,*)
      write (4,*) null
      write (4,*) nstep," time, D, thmax, thmin"
      write (6,*)
      write (8,*)

      Do i=1,Kstep
        write (3,107)   time(i),dxy(i),thmax(i),thmin(i)
        write (4,103) i,time(i),dxy(i),thmax(i),thmin(i)
        write (6,107)   time(i),dxy(i),thmax(i),thmin(i)
        write (8,107)   time(i),dxy(i),thmax(i),thmin(i)
      End Do
      write (4,*) null

c---
c record run parameters
c---

      write (3,205) ndiv
     +             ,cont_angle
     +             ,vs1,vs2
     +             ,tinit
     +             ,shrt

      write (4,205) ndiv
     +             ,cont_angle
     +             ,vs1,vs2
     +             ,tinit
     +             ,shrt

      write (6,205) ndiv
     +             ,cont_angle
     +             ,vs1,vs2
     +             ,tinit
     +             ,shrt

      write (8,205) ndiv
     +             ,cont_angle
     +             ,vs1,vs2
     +             ,tinit
     +             ,shrt

      write (3,208) mint,NGL,dt,Irk,move,Norm,Isym_xy
      write (4,208) mint,NGL,dt,Irk,move,Norm,Isym_xy
      write (6,208) mint,NGL,dt,Irk,move,norm,Isym_xy
      write (8,208) mint,NGL,dt,Irk,move,Norm,Isym_xy

c--------------------------
c print Matlab file (file 1)
c--------------------------

       write (1,*) None

       Do k=1,nelmprint
         call printel (k,Index,c)  ! print in file "drop_3d.net"
       End Do

      ! points around the mid-section

      write (1,*) nxy

      Do i=1,nxy
        nid  = jxy(i)
        write (1,102) p(nid,1),p(nid,2),p(nid,3)
      End Do

      write (1,*) Null

c-----
c close files
c-----

      close (1)
      close (2)
      close (3)
      close (4)

c-----
c done
c-----

  100 Format(1x,i4,10(1x,f12.5))
  101 Format(10(1x,f12.5))
  102 Format(10(1x,f10.6))
  103 Format(1x,i3,10(1x,f10.5))

  105 Format (" Step ",i3," out of ",i4,"; time :",F15.10)
  106 Format (' T=',F7.3,' S=',F10.7,' V=',F10.7
     +       ,' X=',F8.5,' Y=',F 8.5,' Z=',F8.5
     +       )
  107 Format (' T=',F7.3,' D=',F10.7,' thmax=',F10.4,
     +                               ' thmin=',F10.4)
  108 Format (1x,i4,100(1x,f15.10))

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)
  112 Format(" Centroid     :",3(F15.10))

  115 Format(" Deformation  :",F15.10)
  116 Format(" Max Incl     :",F10.4)
  117 Format(" Min Incl     :",F10.4)
  118 Format(" Eff shear st :",F10.4)
  119 Format(" Eff first nsd:",F10.4)
  120 Format(" Eff sec   nsd:",F10.4)

  200 Format(100(1x,f5.3))

 205  Format (/
     +       ,' ndiv    = ',I2,/
     +       ,' cont_an = ',F7.4,/
     +       ,' visc1   = ',F10.5,/
     +       ,' visc2   = ',F10.5,/
     +       ,' tinit   = ',F7.4,/
     +       ,' shrt    = ',F7.4,/
     +       )

 208  Format (/
     +       ,' mint    = ',I3,/
     +       ,' NGL     = ',I3,/
     +       ,' DT      = ',F8.6,/
     +       ,' RUNGE   = ',I1,/
     +       ,' Move    = ',I1,/
     +       ,' Norm    = ',I1,/
     +       ,' Isym_xy = ',I1,/
     +                )

      stop
      end
