      program prtcl_3d_mob 

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c-----------------------------------------------
c Stokes flow past a three-dimensional rigid particle
c with a specified force and torque
c
c The flow is computed by solving the standard boundary integral
c equation involving the single-layer potential alone
c
c SYMBOLS:
c --------
c
c a:  x particle semi-axis
c b:  y particle semi-axis
c c:  z particle semi-axis 
c
c  p(i,j)	coordinates of the ith node (j=1,2,3)
c
c  ne (k,j)     ne(k,1) is the number of elements adjacent to point k
c               ne(k,2), ... are the elements numbers, j = 2, ..., 7
c               for this triangulation, up to six
c
c  n(k,i)       connectivity table: points for element k, i = 1,...,6
c
c  nbe(k,j)     the three neighboring elements of element k (j=1,2,3)
c
c  arel(i)	surface area of ith element
c
c  npts		total number of points
c
c  nelm		total number of elements
c
c  x0, y0, z0         coordinates of collocation points
c  vnx0, vny0, vnz0   coordinates of normal vector at collocation points
c
c  alpha, beta, gamma: parameters for quadratic 
c                      xi-eta isoparametric mapping
c
c  eigen:	Approximate eigenvector of the discrete system
c  eigent:	Approximate eigenvector of the transpose
c               of the discrete system
c
c  Iflow:
c
c  1: infinite flow
c  2: flow bounded by a walllocated at y=wall
c  3: triply periodic flow
c  7: Shear flow past a doubly periodic array of particles 
c     above a wall located at y=wall
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension p (1026,3)
      Dimension ne(1026,7)

      Dimension      n(512,6),nbe(512,3)
      Dimension  alpha(512),beta(512),gamma(512)
      Dimension   arel(512),xmom(512), ymom(512),zmom(512)
      Dimension crvmel(512)

      Dimension   x0(512),  y0(512),  z0(512)
      Dimension vnx0(512),vny0(512),vnz0(512)

      Dimension color(1026)

      Dimension RM(1542,1542),rhs(1542),sln(1542)
      Dimension RMinv(1542,1542)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vnx0,vny0,vnz0
      common/geo6/crvmel

      common/var/wall

      common/zwl/zz,ww

      common/trq/xiq,etq,wq

c---
c for the triply periodic Green's functions
c---

      common/sgfr/a11,a12,a13,a21,a22,a23,a31,a32,a33
     +           ,b11,b12,b13,b21,b22,b23,b31,b32,b33
     +           ,ew,cell_vlm

      common/sgfi/Max1,Max2,Method
      common/NsNp/Ns,Np

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi
      pi6 = 6.0D0*pi
      pi8 = 8.0D0*pi

      null = 0
      oot  = 1.0D0/3.0D0

      Nfour = 4
      Nseven = 7

c------------
c preferences
c------------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for flow past a single particle"
      write (6,*) "   in an infinite domain"
      write (6,*) " 2 for flow past a single particle"
      write (6,*) "   in a semi-infinite domain"
      write (6,*) "   above a plane wall"
      write (6,*) " 3 for flow past a triply periodic array"
      write (6,*) " 7 for flow past a doubly periodic array"
      write (6,*) "   above a plane wall"
      write (6,*) " 0 to quit"
      write (6,*) " ----------"

      read  (5,*) Iflow

      if(Iflow.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 to read parameters from file prtcl_3d_mob.dat" 
      write (6,*) " 2 to enter flow parameters"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"

      read  (5,*) Ienrd

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then   ! read parameters
c------------------------

      open (2,file="prtcl_3d_mob.dat")

      read (2,*) Ioctaicos
      read (2,*) ndiv 
      read (2,*) 
      read (2,*) boa,coa
      read (2,*) req
      read (2,*) cxp,cyp,czp
      read (2,*) phi1,phi2,phi3
      read (2,*) 
      read (2,*) Forcex,Forcey,Forcez
      read (2,*) Torqux,Torquy,Torquz
      read (2,*) 
      read (2,*) mint
      read (2,*) NGL
      read (2,*) 
      read (2,*) visc
      read (2,*) 
      read (2,*) shrt
      read (2,*) 
      read (2,*) wall
      read (2,*) 
      read (2,*) a11,a12,a13
      read (2,*) a21,a22,a23
      read (2,*) a31,a32,a33
      read (2,*) 
      read (2,*) Max1
      read (2,*) Max2
      read (2,*) 
      read (2,*) a11,a12
      read (2,*) a21,a22
      read (2,*) 
      read (2,*) Ns,Np

      close (2)

c---------
      else     ! type the parameters (this list is icomplete)
c---------

      call verbal
     +
     + (Iflow
     + ,ndiv
     + ,boa,coa
     + ,req
     + ,cxp,cyp,czp
     + ,phi1,phi2,phi3
     + ,mint
     + ,NGL
     + ,visc
     + ,Uinf
     + ,shrt
     + ,Iprec
     + ,Ireg
     + )

c-----------
      end if   ! end of reading parameters
c-----------

c---
c uncomment to verify successful input
c---

c     write (6,*) Iflow
c     write (6,*) ndiv 
c     write (6,*) boa,coa
c     write (6,*) req
c     write (6,*) cxp,cyp,czp
c     write (6,*) phi1,phi2,phi3
c     write (6,*) mint
c     write (6,*) NGL
c     write (6,*) visc
c     write (6,*) wall
c     write (6,*) a11,a12,a13
c     write (6,*) a21,a22,a23
c     write (6,*) a31,a32,a33
c     write (6,*) Max1,Max2

c---------------
c prepare to run
c---------------

      open (1,file="prtcl_3d_mob.net")

      phi1 = phi1*pi   ! scale
      phi2 = phi2*pi   ! scale
      phi3 = phi3*pi   ! scale

c---------------------
c triply periodic flow
c---------------------

      if(Iflow.eq.3) then

c----
c function sgf_3d_3p will generatethe reciprocal vectors
c and the optimal value of the splitting parameter xi
c called ew
c 
c function sgf_3d_3p_qqq will generate an array used to compute
c the sum of the velocity Green's function
c in reciprocal space
c----

      call sgf_3d_3p_ewald
     +
     +   (a11,a12,a13
     +   ,a21,a22,a23
     +   ,a31,a32,a33
     +   ,b11,b12,b13
     +   ,b21,b22,b23
     +   ,b31,b32,b33
     +   ,ew,cell_vlm
     +   )

      call sgf_3d_3p_ewald
     +
     +   (b11,b12,b13
     +   ,b21,b22,b23
     +   ,b31,b32,b33
     +   ,Max2,ew
     +   )

      end if

c----------------------------
c triangulate the unit sphere
c----------------------------

      !---
      if(Ioctaicos==1) then
      !---

      call trgl6_octa 
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      !---
      else
      !---

      call trgl6_icos
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      !---
      end if
      !---
      
      write (6,*)
      write (6,*) " prtcl_3d_mob: number of nodes   : ",npts
      write (6,*) " prtcl_3d_mob: number of elements: ",nelm
      write (6,*)

c---------------------------
c expand to specified shape
c and equivalent radius 
c--------------------------

      scale = req/(boa*coa)**oot

      x_axis = scale
      y_axis = scale*boa
      z_axis = scale*coa

      Do i=1,npts                  ! scale
        p(i,1) = x_axis*p(i,1)
        p(i,2) = y_axis*p(i,2)
        p(i,3) = z_axis*p(i,3)
      End Do

      write (6,*) " ellipsoid x semi-axis = ",x_axis
      write (6,*) " ellipsoid y semi-axis = ",y_axis
      write (6,*) " ellipsoid z semi-axis = ",z_axis

c--------------------------
c  rotate by phi1,phi2,phi3
c  around the x,y,z axes
c--------------------------

      cs = cos(phi1)
      sn = sin(phi1)

      Do i=1,npts                  ! rotate about the x axis
       tmpx = p(i,1)
       tmpy = cs*p(i,2)+sn*p(i,3)
       tmpz =-sn*p(i,2)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = cos(phi2)
      sn = sin(phi2)

      Do i=1,npts                   ! rotate about the y axis
       tmpx = cs*p(i,1)-sn*p(i,3)
       tmpy = p(i,2)
       tmpz = sn*p(i,1)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = cos(phi3)
      sn = sin(phi3)

      Do i=1,npts                  ! rotate about the z axis
       tmpx = cs*p(i,1)+sn*p(i,2)
       tmpy =-sn*p(i,1)+cs*p(i,2)
       tmpz = p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      phi1 = phi1/pi   ! unscale
      phi2 = phi2/pi
      phi3 = phi3/pi

c---------------------
c  translate center to
c  specified position
c---------------------

      Do i=1,npts                  ! translate
        p(i,1) = p(i,1) + cxp
        p(i,2) = p(i,2) + cyp
        p(i,3) = p(i,3) + czp
      End Do

c----------------------
c display the triangles
c----------------------

      index = 1  ! 6-node triangles
      index = 2  ! 3-node triangles

      if(index.eq.1) then
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      else if(index.eq.2) then
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end If

      Do i=1,npts    ! color
       color(i) = 0.4;
      End Do

      Do k=1,nelm
        call printel(k,index,color)
      End Do

      write (1,*) null

c---------------
c prepare to run
c---------------

      nelm2 = 2*nelm

c---------------------
c read the quadratures
c---------------------

      call gauss_legendre (NGL,zz,ww)

      call trgl_quad (mint,xiq,etq,wq)

c---------------------------------------------
c compute the coefficients alpha, beta, gamma,
c for quadratic xi-eta mapping
c of each element
c---------------------------------------------

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

c-----------------------------------------------
c compute:  
c
c  (a) The surface area of the individual elements
c        x, y, and z moments over each element
c  (b) The particle surface area and volume
c  (c) The mean curvature of each element: crvmel
c  (d) The averaged normal vector
c-----------------------------------------------

      call elm_geom 
     +
     +  (nelm,npts,mint
     +  ,xmom,ymom,zmom
     +  ,srf_area,prt_vlm
     +  ,cx,cy,cz
     +  )

c---
c normalize surface area and volume
c---

      srf_area_n = srf_area/(pi4*req*req)
      prt_vlm_n  = prt_vlm /(pi4*req*req*req/3.0D0)

      write (6,*)
      write (6,110) srf_area_n
      write (6,111) prt_vlm_n
      write (6,112) cx,cy,cz
      write (6,*)

c----------------------------------
c     write (6,*) 
c     write (6,*) " See the element mean curvature ?"
c     write (6,*) 
c     read  (5,*) Isee
c     If(Isee.eq.1) then
c      Do k = 1,nelm
c        write (6,100) k,crvmel(k)
c      End Do
c     End If
c----------------------------------

c---------------------------
c build the influence matrix
c---------------------------

c----------------------------------
c collocation points will be placed
c at the element centers
c
c Compute:
c
c coordinates of collocation points
c normal vector at collocation points
c----------------------------------

      xi  = 1.0D0/3.0D0
      eta = 1.0D0/3.0D0

      Do i=1,nelm

        Ichoose = 2   ! will interpolate for the normal vector

        i1 = n(i,1)
        i2 = n(i,2)
        i3 = n(i,3)
        i4 = n(i,4)
        i5 = n(i,5)
        i6 = n(i,6)

        call interp_p 
     +
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +    ,alpha(i),beta(i),gamma(i)
     +    ,xi,eta
     +    ,x0(i),y0(i),z0(i)
     +    ,DxDxi,DyDxi,DzDxi
     +    ,DxDet,DyDet,DzDet
     +    ,vnx0(i),vny0(i),vnz0(i)
     +    ,hxi,het,hs
     +    ,Ichoose
     +    )

      End Do

      write (6,*)
      write (6,*) " prtcl_3d_mob: properties of collocation points"
      write (6,*) "               computed"
      write (6,*)

c--------------------
c record the elements
c--------------------

      open (8,file="prtcl_3d_mob.out")

        write (8,*) nelm

        Do i=1,nelm
          i1 = n(i,1)
          i2 = n(i,2)
          i3 = n(i,3)
          i4 = n(i,4)
          i5 = n(i,5)
          i6 = n(i,6)
          write (8,103) p(i1,1),p(i1,2),p(i1,3)
          write (8,103) p(i2,1),p(i2,2),p(i2,3)
          write (8,103) p(i3,1),p(i3,2),p(i3,3)
          write (8,103) p(i4,1),p(i4,2),p(i4,3)
          write (8,103) p(i5,1),p(i5,2),p(i5,3)
          write (8,103) p(i6,1),p(i6,2),p(i6,3)
          write (8,103) arel(i)
        End Do

        write (8,*)
        Do i=1,nelm
         write (8,103) x0(i),y0(i),z0(i)
        End Do

      close (8)

c-------------------------------
c Generate the influence matrix
c three rows at a time
c corresponding to the x, y, z 
c components of the integral equation
c------------------------------------

      write (6,*) "Now generating the influence matrix"

      cf = 1.0D0/(pi8*visc)

      Do i=1,nelm

        write (6,*) " prtcl_3d_mob: collocating at point: ",i

        inelm  = i+nelm
        inelm2 = i+nelm+nelm

        Do j=1,nelm

          jnelm  = j+nelm
          jnelm2 = j+nelm+nelm

c-----------------------
          if(i.ne.j) then     ! regular element
c-----------------------
      
          call slp_trgl6
     +
     +     (x0(i),y0(i),z0(i)
     +     ,j
     +     ,GExx,GExy,GExz
     +     ,GEyx,GEyy,GEyz
     +     ,GEzx,GEzy,GEzz
     +     ,mint
     +     ,Iflow
     +     )

c-------------
          else     ! singular element
c-------------

          call slp_trgl6_sing 
     +
     +      (x0(i),y0(i),z0(i)
     +      ,j
     +      ,GExx,GExy,GExz
     +      ,GEyx,GEyy,GEyz
     +      ,GEzx,GEzy,GEzz
     +      ,NGL
     +      ,Iflow
     +      )

c---------------
          end if
c---------------

          RM(i,j)      = cf*GExx
          RM(i,jnelm)  = cf*GEyx
          RM(i,jnelm2) = cf*GEzx

          RM(inelm,j)      = cf*GExy
          RM(inelm,jnelm)  = cf*GEyy
          RM(inelm,jnelm2) = cf*GEzy

          RM(inelm2,j)      = cf*GExz
          RM(inelm2,jnelm)  = cf*GEyz
          RM(inelm2,jnelm2) = cf*GEzz

c------------------------------
c deflate the integral equation
c by adding the term: n int( Df . n dS)
c in discretized form
c------------------------------

         Ideflate = 1;

         if(Ideflate.eq.1) then
            aj = arel(j)

            RM(i,j)      = RM(i,j)      
     +                   + vnx0(i)*vnx0(j)*aj
            RM(i,jnelm)  = RM(i,jnelm)
     +                   + vnx0(i)*vny0(j)*aj
            RM(i,jnelm2) = RM(i,jnelm2)
     +                   + vnx0(i)*vnz0(j)*aj

            RM(inelm,j)      = RM(inelm,j)
     +                       + vny0(i)*vnx0(j)*aj
            RM(inelm,jnelm)  = RM(inelm,jnelm)
     +                       + vny0(i)*vny0(j)*aj
            RM(inelm,jnelm2) = RM(inelm,jnelm2)
     +                       + vny0(i)*vnz0(j)*aj

            RM(inelm2,j)      = RM(inelm2,j)
     +                        + vnz0(i)*vnx0(j)*aj
            RM(inelm2,jnelm)  = RM(inelm2,jnelm)
     +                        + vnz0(i)*vny0(j)*aj
            RM(inelm2,jnelm2) = RM(inelm2,jnelm2)
     +                        + vnz0(i)*vnz0(j)*aj
           end if

        End Do

      End Do
      
c--------------------------------------
c system size emerging from collocation
c--------------------------------------

      Mdim = 3*nelm

      Mls = Mdim  ! size of the linear system emerging from collocation
      Mls6 = Mls+6

c--------------------------------
c introduce six bordering columns
c to accomodate the linear
c and angular velocities
c--------------------------------

      Do i=1,nelm

         xhat = x0(i)-cxp
         yhat = y0(i)-cyp
         zhat = z0(i)-czp

         RM(i,Mls+1) = 1.0D0
         RM(i,Mls+2) = 0.0D0
         RM(i,Mls+3) = 0.0D0
         RM(i,Mls+4) = 0.0D0
         RM(i,Mls+5) = zhat
         RM(i,Mls+6) =-yhat 

         ii = nelm+i
         RM(ii,Mls+1) = 0.0D0
         RM(ii,Mls+2) = 1.0D0
         RM(ii,Mls+3) = 0.0D0
         RM(ii,Mls+4) = -zhat
         RM(ii,Mls+5) = 0.0D0
         RM(ii,Mls+6) = xhat

         ii = nelm+nelm+i
         RM(ii,Mls+1) = 0.0D0
         RM(ii,Mls+2) = 0.0D0
         RM(ii,Mls+3) = 1.0D0
         RM(ii,Mls+4) = yhat
         RM(ii,Mls+5) =-xhat
         RM(ii,Mls+6) = 0.0D0

      End Do

c--------------------------------
c introduce six bordering bottom rows
c to specify the force and torque
c--------------------------------

      Do j=1,Mls+6
       RM(Mls+1,j) = 0.0D0
       RM(Mls+2,j) = 0.0D0
       RM(Mls+3,j) = 0.0D0
       RM(Mls+4,j) = 0.0D0
       RM(Mls+5,j) = 0.0D0
       RM(Mls+6,j) = 0.0D0
      End Do

      Do j=1,nelm

       xhat = x0(j)-cxp
       yhat = y0(j)-cyp
       zhat = z0(j)-czp

       RM(Mls+1,j)           =  arel(j)   ! force
       RM(Mls+2,j+nelm)      =  arel(j)   ! force
       RM(Mls+3,j+nelm+nelm) =  arel(j)   ! force

       RM(Mls+4,j+nelm)      = -arel(j)*zhat ! torque
       RM(Mls+4,j+nelm+nelm) =  arel(j)*yhat

       RM(Mls+5,j)           =  arel(j)*zhat
       RM(Mls+5,j+nelm+nelm) = -arel(j)*xhat

       RM(Mls+6,j)           = -arel(j)*yhat
       RM(Mls+6,j+nelm)      =  arel(j)*xhat

      End Do

c     Do i=1,Mls+6
c         write (6,*) (RM(i,j),j=1,Mls+6)
c     pause
c     End Do

c------------------------
c set the right-hand side
c------------------------

      write (6,*) 
      write (6,*) " Compute and print the inverse?"
      write (6,*) 
      read (5,*) Icompute

      if(Icompute.eq.1) then

      call gel_inv
     +
     +  (Mls+6
     +  ,RM
     +  ,RMinv
     +  ,Isym,Iwlpvt
c    +  ,l,u
     +  ,det
     +  ,Istop
     +  )

       open (3,file="matrix_inverse.out")

        Do i=1,Mls+6
         Do j=1,Mls+6
          write (3,103) RMinv(i,j)
         End Do
        End Do

       close (3)

      end if

c------------------------
c set the right-hand side
c------------------------

        Do i=1,Mdim        ! initialize
          rhs(i) = 0.0D0
        End Do

        Do i=1,nelm

         inelm  = i+nelm
         inelm2 = i+nelm+nelm

         rhs(i)      = 0.0D0
         rhs(i)      = shrt*y0(i)   ! simple shear flow
         rhs(inelm)  = 0.0D0 
         rhs(inelm2) = shrt*y0(i)   ! simple shear flow
         rhs(inelm2) = 0.0D0 

        End Do

      rhs(Mdim+1) = -Forcex   ! hydrodynamic force is negative
      rhs(Mdim+2) = -Forcey   ! of external force
      rhs(Mdim+3) = -Forcez

      rhs(Mdim+4) = -Torqux
      rhs(Mdim+5) = -Torquy
      rhs(Mdim+6) = -Torquz

c------------------------
c solve the linear system
c------------------------

c     pause " Print the linear system ?"
c     Do i=1,Mdim
c      write (6,200) (RM(i,j),j=1,Mdim),(rhs(i,j),j=1,nrhs)
c     End Do

      Isymg  = 0    ! system is not symmetric
      Iwlpvt = 1    ! pivoting enabled

      write (6,*) " prtcl_3d_mob: system size: ",Mls+6
c     write (6,*) " prtcl_3d_mob: solving the linear system"

      call gel
     +
     +   (Mls+6
     +   ,RM,rhs,sln
     +   ,Isymg,Iwlpvt
     +   ,det
     +   ,Istop
     +   )

c     write (6,*) " prtcl_3d_mob: linear system solved"
c     write (6,*) " prtcl_3d_mob: determinant = ",det

c---------------------
c display the solution
c---------------------

c     write (6,*)
c     write (6,*) " prtcl_3d_mob: tractions for translation modes"
c     write (6,*)

c     Do i=1,nelm
c      inelm  = i+nelm
c      inelm2 = i+nelm2
c      write (6,210) i,sln(i),sln(inelm),sln(inelm2)
c     End Do

      Ux = sln(Mls+1)
      Uy = sln(Mls+2)
      Uz = sln(Mls+3)

      Ox = sln(Mls+4)
      Oy = sln(Mls+5)
      Oz = sln(Mls+6)

      write (6,*)
      write (6,*) " prtcl_3d_mob: Linear and Angular Velocity"
      write (6,*)

      write (6,205) Ux,Uy,Uz
      write (6,205) Ox,Oy,Oz

c-----
c Done
c-----

  99  continue

c---
c de-scale
c---

      phi1 = phi1/pi
      phi2 = phi2/pi
      phi3 = phi3/pi

      close (1)

  100 Format(1x,i4,10(1x,f12.5))
  101 Format(10(1x,f12.5))
  102 Format(10(1x,f10.6))
  103 Format(10(1x,f15.10))
  110 Format("Surface Area :",F15.10)
  111 Format("Volume       :",F15.10)
  112 Format("Centroid     :",3(F15.10))
  130 Format("Determinant =",f20.10)
  200 Format(100(1x,f5.3))
  201 Format(10(1x,f7.5))
  205 Format(10(1x,f15.10))

  210 Format(1x,i4,10(1x,f7.3))

      stop
      end
