      program prtcl_3d 

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
c  This program computes the following Stokes flows:
c
c  1) Flow due to the motion of
c     a three-dimensional rigid particle in free space
c     (Iflow=1)
c
c  2) Flow due to the motion of
c     a three-dimensional rigid particle in a semi-infinite
c     domain bounded by a plane wall located at y=wall
c     (Iflow=2)
c
c  3) Uniform flow past
c     a triply periodic lattice of rigid particles
c     (Iflow=3)
c
c  4) Uniform flow normal to a doubly periodic array of particles 
c     representing a screen
c     (Iflow=4)
c
c  5) Shear flow over a doubly periodic array of particles 
c     representing a screen.
c     (Iflow=5)
c
c  6) Shear flow past a spherical bump on a plane wall
c     located at y=wall
c     (Iflow=6)
c
c  7) Shear flow past a doubly periodic array of particles 
c     above a wall located at y=wall
c     (Iflow=7)
c
c  In cases 1--5, 7 the particle is an ellipsoid with:
c
c     x semi-axis equal to a
c     y semi-axis equal to b
c     z semi-axis equal to c
c
c SYMBOLS:
c --------
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
c  vnx0, vny0, vnz0   unit normal vector at collocation points
c
c  alpha, beta, gamma: parameters for quadratic 
c                      xi-eta isoparametric mapping
c
c  Resist(i,j):   Grand resistance matrix
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

      Dimension fx(512),fy(512),fz(512)

      Dimension vna(1026,3)
      Dimension fx_node(1026),fy_node(1026),fz_node(1026)
      Dimension fshear(1026)

      Dimension     RM(1600,1600)            ! Influence matrix
      Dimension  RMinv(1600,1600)            ! Inverse 
      Dimension rhs(1600,10),sln(1600,10)    ! RHSs and solution vectors

      Dimension dfl_rm(512),dfl_rhs(10)
      Dimension residual(7)

      Dimension Resist(6,6)

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
      common/geo3/vna
      common/geo6/crvmel

      common/var/wall

      common/zwl/zz,ww

      common/trq/xiq,etq,wq

      common/sys1/RM,rhs

c---
c for the Green's functions
c---

      common/sgfr/a11,a12,a13,a21,a22,a23,a31,a32,a33
     +           ,b11,b12,b13,b21,b22,b23,b31,b32,b33
     +           ,ew,cell_vlm,cell_area
     +           ,eps

      common/sgfi/Max1,Max2,Method
      common/NsNp/Ns,Np

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.500*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi
      pi6 = 6.0D0*pi
      pi8 = 8.0D0*pi

      null = 0
      Nfour = 4
      Nseven = 7

      oot  = 1.0D0/3.0D0

c------------
c preferences
c------------

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for flow due to the motion of particle"
      write (6,*) "   in free space"
      write (6,*) " 2 for flow due to the motion of particle"
      write (6,*) "   in a semi-infinite domain bounded"
      write (6,*) "   by a plane located at y=wall"
      write (6,*) " 3 for uniform flow past"
      write (6,*) "    a triply periodic lattice"
      write (6,*) " 4 for uniform flow normal to a doubly"
      write (6,*) "    periodic lattice representing latice"
      write (6,*) " 5 for shear flow over a doubly"
      write (6,*) "    periodic lattice representing a screen"
      write (6,*) " 6 for shear flow past a spherical bump"
      write (6,*) "   on a plane wall located at y=wall"
      write (6,*) " 7 for flow due to the motion "
      write (6,*) "   of a doubly periodic array of particles"
      write (6,*) "   in a semi-infinite domain bounded"
      write (6,*) "   by a plane located at y=wall"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"

      read (5,*) Iflow

      if(Iflow.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 to read parameters from file prtcl_3d.dat" 
      write (6,*) " 2 to enter flow parameters"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"

      read  (5,*) Ienrd

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then   ! read parameters
c------------------------

      open (2,file="prtcl_3d.dat")

      read (2,*) Ioctaicos
      read (2,*) ndiv 
      read (2,*) 
      read (2,*) boa,coa
      read (2,*) req
      read (2,*) cxp
      read (2,*) cyp
      read (2,*) czp
      read (2,*) phi1,phi2,phi3
      read (2,*) cont_angle
      read (2,*) 
      read (2,*) mint
      read (2,*) NGL
      read (2,*) 
      read (2,*) visc
      read (2,*) 
      read (2,*) Uinf
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
      read (2,*) Ns,Np
      read (2,*) 
      read (2,*) Method
      read (2,*) eps
      read (2,*) 
      read (2,*) Iprec
      read (2,*) Ireg

      close (2)

c-----------
      else     ! type in the parameters (this list is incomplete)
c-----------

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
c     write (6,*) Iprec
c     write (6,*) Ireg

c--------
c prepare
c--------

      open (1,file="prtcl_3d.net")  ! output for matlab visualization

c---------------------
c triply periodic flow
c---------------------

      if(Iflow.eq.3) then

c----
c function sgf_3d_3p_ewald will generate the reciprocal vectors
c and the optimal value of the splitting parameter xi
c called ew
c 
c function sgf_3d_2p_qqq will generate an array used to compute
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

      call sgf_3d_3p_qqq
     +
     +   (b11,b12,b13
     +   ,b21,b22,b23
     +   ,b31,b32,b33
     +   ,Max2,ew
     +   )

      end if

c---------------------
c doubly periodic flow
c---------------------

      if(Iflow.eq.4.or.Iflow.eq.5) then

c----
c ewald_3d_2p will produce the reciprocal vectors
c and the optimal value of xi
c called ew
c---

      call sgf_3d_2p_ewald
     +
     +  (a11,a12
     +  ,a21,a22
     +  ,b11,b12
     +  ,b21,b22
     +  ,ew,cell_area
     +  )

      write (6,*) " prtcl_3d: cell area = ",cell_area

      end if

c----------------------------
c triangulate the unit sphere
c
c Discretize the shape and generate
c the connectivity table
c----------------------------

c------------------------
      if(Iflow.ne.6) then
c------------------------

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
      
c--------------------------
c expand to specified shape
c and equivalent radius 
c--------------------------

      scale = req/(boa*coa)**oot
c     scale = 1.0D0/boa                 ! reset

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

c-------------------------
c rotate by phi1,phi2,phi3
c around the x,y,z axes
c-------------------------

      phi1 = phi1*pi   ! scale
      phi2 = phi2*pi   ! scale
      phi3 = phi3*pi   ! scale

      cs = Dcos(phi1)
      sn = Dsin(phi1)

      Do i=1,npts                  ! rotate about the x axis
       tmpx = p(i,1)
       tmpy = cs*p(i,2)+sn*p(i,3)
       tmpz =-sn*p(i,2)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = Dcos(phi2)
      sn = Dsin(phi2)

      Do i=1,npts                   ! rotate about the y axis
       tmpx = cs*p(i,1)-sn*p(i,3)
       tmpy = p(i,2)
       tmpz = sn*p(i,1)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = Dcos(phi3)
      sn = Dsin(phi3)

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
c translate center to
c specified position
c---------------------

      Do i=1,npts 
        p(i,1) = p(i,1) + cxp
        p(i,2) = p(i,2) + cyp
        p(i,3) = p(i,3) + czp
      End Do

c-----------
      end if
c-----------

c----------------------------
      if(Iflow.eq.6) then    !  sphere attached to a wall
c----------------------------

      call trgl6_octa_ss 
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

c---------------------
c expand to specified radius (req)
c slide to specified contact angle
c translate center to specified position
c---------------------

      Do i=1,npts
       p(i,1) = req*p(i,1)
       p(i,2) = req*p(i,2)
       p(i,3) = req*p(i,3)
      End Do

      phi   = pi*cont_angle
      cs    = Dcos(phi)
      shift = req*Dsin(pih-phi)

      Do i=1,npts
       rr    = Dsqrt(p(i,1)**2+p(i,3)**2)
       thet0 = Datan2(rr,p(i,2))
       thetz = Datan2(p(i,1),p(i,3))
       thet  = thet0*phi/pih
       rrr   = req*Dsin(thet)
       p(i,1) = rrr*Dsin(thetz)
       p(i,2) = req*Dcos(thet)
       p(i,3) = rrr*Dcos(thetz)
       p(i,2) = p(i,2)-shift
       p(i,2) = p(i,2) + wall
      End Do

c-----------
      end if
c-----------

      write (6,*)
      write (6,*) " prtcl_3d: number of nodes   : ",npts
      write (6,*) " prtcl_3d: number of elements: ",nelm
      write (6,*)

c---------------
c prepare to run
c---------------

      nelm2 = 2*nelm

c---------------------
c read the quadratures
c---------------------

      call gauss_leg (NGL,zz,ww)

      call gauss_trgl
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

c---------------------------------------------
c compute the coefficients alpha, beta, gamma,
c for the quadratic xi-eta mapping
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

c---------------------------------------------
c compute:  
c
c     vlm:  surface area of the individual elements
c     xmom,yomo,zmom:  x, y, and z surface moments
c                      over each element
c     area:   total surface area and volume
c     crvmel: mean curvature over each element
c     vna:    average normal vector at the nodes
c---------------------------------------------

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

c------------------------------
c  write (6,*) 
c  write (6,*) "  Proceed to compute the influence matrix ?"
c  write (6,*) 
c  read  (5,*) Igo
c  If(Igo.eq.0) Go to 99
c------------------------------

c----------------------------------
c Collocation points will be placed
c at the element centroids
c
c Compute:
c
c  1) coordinates of collocation points
c  2) normal vector at collocation points
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
     +
     +    ,alpha(i),beta(i),gamma(i)
     +
     +    ,xi,eta
     +
     +    ,x0(i),y0(i),z0(i)
     +
     +    ,DxDxi,DyDxi,DzDxi
     +    ,DxDet,DyDet,DzDet
     +    ,vnx0(i),vny0(i),vnz0(i)
     +    ,hxi,het,hs
     +
     +    ,Ichoose
     +    )

      End Do

      write (6,*)
      write (6,*) " prtcl_3d: collocation points computed"
      write (6,*)

c------------------------------
c record the collocation points
c------------------------------

      open (3,file="particle_points.out")

        write (3,*) nelm

        Do i=1,nelm
         write (3,103) x0(i),y0(i),z0(i)
        End Do

      close (3)

c--------------------
c record the elements
c--------------------

      open (8,file="particle_elements.out")

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
c Generate the influence matrix,
c three rows at a time,
c corresponding to the x, y, z 
c components of the integral equation
c------------------------------------

      cf = -1.0D0/(pi8*visc)

      Do i=1,nelm            ! run over the collocation points

       write (6,*) " prtcl_3d: collocating at point: ",i

       inelm  = i+nelm
       inelm2 = i+nelm+nelm

       Do j=1,nelm            ! run over the elements

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
          else               ! singular element
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

c--------
c Deflate the integral equation
c-------

         if(Ireg.eq.3) then

          aj = arel(j)

          RM(i,j)      = RM(i,j)      + vnx0(i)*vnx0(j)*aj
          RM(i,jnelm)  = RM(i,jnelm)  + vnx0(i)*vny0(j)*aj
          RM(i,jnelm2) = RM(i,jnelm2) + vnx0(i)*vnz0(j)*aj

          RM(inelm,j)      = RM(inelm,j)     + vny0(i)*vnx0(j)*aj
          RM(inelm,jnelm)  = RM(inelm,jnelm) + vny0(i)*vny0(j)*aj
          RM(inelm,jnelm2) = RM(inelm,jnelm2)+ vny0(i)*vnz0(j)*aj

          RM(inelm2,j)      = RM(inelm2,j)     + vnz0(i)*vnx0(j)*aj
          RM(inelm2,jnelm)  = RM(inelm2,jnelm) + vnz0(i)*vny0(j)*aj
          RM(inelm2,jnelm2) = RM(inelm2,jnelm2)+ vnz0(i)*vnz0(j)*aj

        end if
c--------

        End Do

      End Do
      
c------------
c system size
c------------

      Mdim = 3*nelm

      Mls = Mdim  ! size of the linear system emerging from collocation

c     Do i=1,Mls
c      write (6,*) (RM(i,j),j=1,Mls)
c      pause
c     End Do

c-------------------------------------------------
c test for uniform flow past a solitary sphere
c in free space.
c
c In this case, the exact solution of the integral
c equation for translation
c is known to be: 
c
c traction = 1.5 visc*U/req
c------------------------------------------------

      if (Iflow.eq.1        ! solitary spherical particle
     +   .and.boa.eq.1.0.and.coa.eq.1.0
     +   ) then

      Do i=1,nelm

       sum1 = 0.0D0
       sum2 = 0.0D0
       sum3 = 0.0D0
       sum4 = 0.0D0
       sum5 = 0.0D0
       sum6 = 0.0D0
       sum7 = 0.0D0
       sum8 = 0.0D0
       sum9 = 0.0D0

       inelm = i+nelm
       inelm2= inelm+nelm

       Do j=1,nelm

        jnelm  = j+nelm
        jnelm2 = jnelm+nelm
 
        sum1 = sum1 + RM(i,j)
        sum2 = sum2 + RM(inelm,j)
        sum3 = sum3 + RM(inelm2,j)

        sum4 = sum4 + RM(i,jnelm)
        sum5 = sum5 + RM(inelm,jnelm)
        sum6 = sum6+rm(inelm2,jnelm)

        sum7 = sum7 + RM(i,jnelm2)
        sum8 = sum8 + RM(inelm,jnelm2)
        sum9 = sum9 + RM(inelm2,jnelm2)

       End Do

       fc = 1.5D0*visc/(req)

       sum1 = fc*sum1
       sum2 = fc*sum2
       sum3 = fc*sum3
       sum4 = fc*sum4
       sum5 = fc*sum5
       sum6 = fc*sum6
       sum7 = fc*sum7
       sum8 = fc*sum8
       sum9 = fc*sum9

c     write (6,201) sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9

      End Do

      End If

c---------------------------
c set the right-hand side(s)
c---------------------------

c----------------------------------------------
c Flow due to particle rigid-body motion
c
c The six right-hand sides
c express the six modes of translation and rotation
c----------------------------------------------

      if(Iflow.eq.1.or.Iflow.eq.2
     +             .or.Iflow.eq.3
     +             .or.Iflow.eq.7
     +                           ) then 


        nrhs = 6           ! number of right-hand sides

        Do i=1,Mdim        ! initialize
         Do k=1,nrhs
           rhs(i,k) = 0.0D0
         End Do
        End Do

        Do i=1,nelm

         inelm  = i+nelm
         inelm2 = i+nelm+nelm

         rhs(i,     1) = 1.0D0  ! translation along the x axis
         rhs(inelm, 2) = 1.0D0  ! translation along the y axis
         rhs(inelm2,3) = 1.0D0  ! translation along the z axis

         xhat = x0(i)-cxp
         yhat = y0(i)-cyp
         zhat = z0(i)-czp
c        xhat = x0(i) ! reset
c        yhat = y0(i)
c        zhat = z0(i)

         rhs(inelm, 4) =-zhat      ! rigid-body rotation
         rhs(inelm2,4) = yhat
         rhs(i,     5) = zhat
         rhs(inelm2,5) =-xhat
         rhs(i,     6) =-yhat
         rhs(inelm, 6) = xhat

        End Do

      End If

c-------------------------------------------------
c Uniform flow normal to the doubly periodic array
c along the z axis
c
c Solve the integral equation for the
c disturbance flow
c--------------------------------------------------

      If(Iflow.eq.4) then

        nrhs = 1             ! number of right-hand sides

        Do i=1,Mdim          ! initialize
         Do k = 1,nrhs
           rhs(i,k) = 0.0D0
         End Do
        End Do

        Do i = 1,nelm
         inelm2 = i+2*nelm
         rhs(inelm2,1) = -Uinf    ! z-component
        End Do

      End If

c--------------------------------------------
c Shear flow over a planar array of particles
c along the x or y axis
c--------------------------------------------

      If(Iflow.eq.5) then

       nrhs = 2             ! number of right-hand sides

       Do i=1,Mdim          ! initialize
         Do k = 1,nrhs
           rhs(i,k) = 0.0
         End Do
       End Do

       Do i=1,nelm
        inelm = i+nelm
        rhs(i,1)     = -shrt*z0(i)
        rhs(inelm,2) = -shrt*z0(i)
       End Do

      End If

c---------------------------------
c Shear flow over a spherical bump
c along the x axis
c---------------------------------

      If(Iflow.eq.6) then

       nrhs = 1             ! number of right-hand sides

       Do i=1,Mdim          ! initialize
         Do k = 1,nrhs
           rhs(i,k) = 0.0D0
         End Do
       End Do

       Do i=1,nelm
        inelm = i+nelm
        rhs(i,1) = -shrt*y0(i)
       End Do

      End If

c-----------------------------------------------
c PRECONDITIONING
c
c Render the coefficient matrix singular
c by premultiplying the linear system by the matrix:
c
c  I-eigent*eigent
c
c eigent is produced in subroutine "precondition"
c------------------------------------------------

      if(Iprec.eq.1) then

       write (6,*) " prtcl_3d: preconditioning"

       call precondition (nelm,Mdim,nrhs)

      end if

c--------------------------
c CONDITION NUMBER
c
c compute the condition number 
c of the coefficient matrix
c--------------------------

      write (6,*) 
      write (6,*) " Compute the condition number"
      write (6,*) " of the influence matrix ?"
      write (6,*) 
      write (6,*) " Enter 0 for NO, 1 for YES"
      write (6,*) " -------------------------"
      write (6,*) 
      Read  (5,*) Icond

      If(Icond.eq.1) then

        maxit = 100        ! maximum number of power iterations
        epsln = 0.000001   ! error in power iterations

        call cond_num (Mls,RM,maxit,epsln,cond)

        write (6,*) 
        write (6,*) "Condition number = ",cond
        write (6,*) 

      End If

c----------------------------
c REGULARIZATION BY DEFLATION
c----------------------------

c------------------------------------
      if(Ireg.eq.1.or.Ireg.eq.2) then
c------------------------------------

       write (6,*) " prtcl_3d: regularizing"

c---
c save the last equation
c---

        Do j=1,Mdim
         dfl_rm(j) = rm(Mdim,j)
        End Do

        Do j=1,nrhs
         dfl_rhs(j) = rhs(Mdim,j)
        End Do

c---------------
c regularization:
c---------------

c-------
       if(Ireg.eq.1) then
c-------

c---
c Replace the last equation 
c with the constraint: Int{ f.n dS } = 0
c--

       Do j=1,nelm
         RM(Mdim,j)           = vnx0(j)*arel(j)
         RM(Mdim,j+nelm )     = vny0(j)*arel(j)
         RM(Mdim,j+nelm+nelm) = vnz0(j)*arel(j)
       End Do

       Do j=1,nrhs
         rhs(Mdim,j) = 0.0D0
       End Do

c-------
       else if(Ireg.eq.2) then
c-------

c-------------------
c Naive method:
c
c Remove the last equation 
c and set the last unknown equal to zero
c by reducing the size of the linear system
c by one unit
c-------------------

        Mls = Mdim-1      ! will solve the top Mdim-1 equations

      end if

c-----------
      End If
c-----------

c--------------------------------------------
c Shear flow past a planar array of particles:
c
c Append two columns and add two more equations.
c
c The two additional unknowns correspond to the
c added columns are the x and y velocity
c components above the array (slip velocity)
c
c The two additional equations specify the force
c exerted on each particle to compensate the
c simple shear flow
c--------------------------------------------

      If(Iflow.eq.5) then

       Mls1 = Mls+1
       Mls2 = Mls+2

c---
c append columns Mls+1 and Mls+2
c coefficents multiply the slip velocity
c components (above the array)
c---

       Do i=1,Mls2
        RM(i,Mls1) = 0.0D0
        RM(i,Mls2) = 0.0D0
       End Do

       Do i=1,nelm      
        inelm = i+nelm
        RM(i    ,Mls1) = 1.0D0
        RM(inelm,Mls2) = 1.0D0
       End Do

c---
c append rows Mls+1 and Mls+2
c---

       Do j=1,Mls2
        RM(Mls1,j) = 0.0D0
        RM(Mls2,j) = 0.0D0
       End Do

       Do j=1,nelm
        jnelm = j+nelm
        RM(Mls1,j    ) = arel(j)
        RM(Mls2,jnelm) = arel(j)
       End Do

       rhs(Mls1,1) = visc*shrt*cell_area    ! x-force constraint
       rhs(Mls2,1) = 0.0D0
       rhs(Mls1,2) = 0.0D0
       rhs(Mls2,2) = visc*shrt*cell_area    ! y-force constraint

       Mls = Mls+2       ! augment size of linear system by 2

      End If

c------------------------
c solve the linear system
c------------------------

c------------
c     pause " Print the linear system ?"
c     Do i = 1,Mdim
c      write (6,200) (RM(i,j),j=1,Mdim),(rhs(i,j),j=1,nrhs)
c     End Do
c------------

      Isymg  = 0    ! system is not symmetric
      Iwlpvt = 1    ! pivoting enabled

      write (6,*) " prtcl_3d: system size: ",Mls

c---------------
      write (6,*) 
      write (6,*) " Compute and print the inverse"
      write (6,*) " of the influence matrix?"
      write (6,*) 
      write (6,*) " Enter 0 for NO, 1 for YES"
      write (6,*) " -------------------------"
      write (6,*) 
      read  (5,*) Invert

      if(Invert.eq.1) then

      write (6,*) " prtcl_3d: computing the matrix inverse"

      call gel_inv
     +
     +  (Mls
     +  ,RM
     +  ,RMinv
     +  ,Isym,Iwlpvt
c    +  ,l,u
     +  ,det
     +  ,Istop
     +  )

       open (3,file="matrix_inverse.out")
        Do i=1,Mls
         Do j=1,Mls
          write (3,103) RMinv(i,j)
         End Do
        End Do
       close (3)

      end if
c-----------

      write (6,*) " prtcl_3d: solving the linear system"

      call gel_mrhs 
     +
     +   (Mls
     +   ,RM
     +   ,nrhs
     +   ,rhs
     +   ,sln
     +   ,Isymg
     +   ,Iwlpvt
     +   ,det
     +   ,Istop
     +   )

      write (6,*) "Determinant = ",det

c----------------------------------------
c Shear flow past a doubly periodic array
c
c Extract the x and y slip velocity 
c as the last two entries of the vector
c of unknowns
c----------------------------------------

      if(Iflow.eq.5) then

       Uslipx = sln(Mls1,1)
       Vslipx = sln(Mls2,1)

       Uslipy = sln(Mls1,2)
       Vslipy = sln(Mls2,2)

      end if

c------------------------------------
c If type 1 or 2 deflation is enabled,
c compute the last-equation residuals 
c-------------------------------------

      if(Ireg.eq.1.or.Ireg.eq.2) then         !    DEFLATION MODULE

        write (6,*)
        write (6,*) "last equation residuals (should be zero)"
        write (6,*)

        Do k=1,nrhs    ! run over right-hand sides

c         sln(Mdim,k) = 0.0D0

          residual(k) = 0.0D0

c         Do j=1,Mdim-1
          Do j=1,Mdim
           residual(k) = residual(k)+dfl_rm(j)*sln(j,k)
          End Do

          residual(k) = residual(k)-dfl_rhs(k)

          write (6,205) residual(k)

        End Do

      end if                      ! END OF DEFLATION MODULE

c--------------------------------------------------------

c---------------------
c Display the solution
c---------------------

c--------------------------------------
      if(Iflow.eq.1.or.Iflow.eq.2
     +             .or.Iflow.eq.3
     +             .or.Iflow.eq.7
     +             ) then
c--------------------------------------

      write (6,*)
      write (6,*) " prtcl_3d: tractions for translation modes:"
      write (6,*)

      Do i=1,nelm
       inelm  = i+nelm
       inelm2 = i+nelm2
       write (6,210) i,sln(i,1),sln(inelm,1),sln(inelm2,1)
     +                ,sln(i,2),sln(inelm,2),sln(inelm2,2)
     +                ,sln(i,3),sln(inelm,3),sln(inelm2,3)
      End Do

      write (6,*)
      write (6,*) " prtcl_3d: tractions for rotation modes:"
      write (6,*)

      Do i=1,nelm
       inelm  = i+nelm
       inelm2 = i+nelm2
       write (6,210) i,sln(i,4),sln(inelm,4),sln(inelm2,4)
     +                ,sln(i,5),sln(inelm,5),sln(inelm2,5)
     +                ,sln(i,6),sln(inelm,6),sln(inelm2,6)
      End Do

c-----------------------------
      else if(Iflow.eq.4) then
c-----------------------------

      write (6,*)
      write (6,*) " prtcl_3d: tractions"
      write (6,*)

      Do i=1,nelm
       inelm  = i+nelm
       inelm2 = i+nelm2
       write (6,210) i,sln(i,1),sln(inelm,1),sln(inelm2,1)
      End Do

c-----------------------------
      else if(Iflow.eq.5) then
c-----------------------------

      write (6,*)
      write (6,*) " prtcl_3d: tractions for shear flow"
      write (6,*) "           along the x or y axis"
      write (6,*)

      Do i=1,nelm
       inelm  = i+nelm
       inelm2 = i+nelm2
       write (6,210) i,sln(i,1),sln(inelm,1),sln(inelm2,1)
     +                ,sln(i,2),sln(inelm,2),sln(inelm2,2)
      End Do

c---
c compute the drift velocity below the array
c---

       sumxx = 0.0D0
       sumxy = 0.0D0
       sumyx = 0.0D0
       sumyy = 0.0D0

       Do i=1,nelm

        sumxx = sumxx + z0(i)*sln(i,1)*arel(i) 
        sumxy = sumxy + z0(i)*sln(i,2)*arel(i) 

        inelm  = i+nelm

        sumyx = sumyx + z0(i)*sln(inelm,1)*arel(i) 
        sumyy = sumyy + z0(i)*sln(inelm,2)*arel(i) 

       End Do

       fc = 1.0D0/(visc*cell_area)

       sumxx = fc*sumxx
       sumxy = fc*sumxy
       sumyx = fc*sumyx
       sumyy = fc*sumyy

       Udriftx = Uslipx + sumxx
       Udrifty = Uslipy + sumxy

       Vdriftx = Vslipx + sumyx
       Vdrifty = Vslipy + sumyy

       write (6,*)
       write (6,*) "Shear flow along the x axis"
       write (6,*)
       write (6,*) "x and y slip and drift velocities:"
       write (6,*) "----------------------------------"
       write (6,205) Uslipx,Vslipx
       write (6,205) Udriftx,Vdriftx

       write (6,*)
       write (6,*) "Shear flow along the y axis"
       write (6,*)
       write (6,*) "x and y slip and drift velocities:"
       write (6,*) "----------------------------------"
       write (6,205) Uslipy,Vslipy
       write (6,205) Udrifty,Vdrifty

c-----------
      end if
c-----------

c------------------------------------
c Assign traction to elements
c
c Compute force, torque, 
c and the grand resistance matrix
c------------------------------------

      cf = -1.0D0/(pi6*req*visc)

      Do k=1,nrhs

       Frcx = 0.0D0
       Frcy = 0.0D0
       Frcz = 0.0D0

       Trqx = 0.0D0
       Trqy = 0.0D0
       Trqz = 0.0D0

       Do i=1,nelm

        inelm  = i+nelm
        inelm2 = inelm+nelm

        fx(i) = sln(i,k)
        fy(i) = sln(inelm,k)
        fz(i) = sln(inelm2,k)

        Frcx = Frcx + fx(i)*arel(i)
        Frcy = Frcy + fy(i)*arel(i)
        Frcz = Frcz + fz(i)*arel(i)

        Trqx = Trqx + fz(i)*ymom(i)-fy(i)*zmom(i)
        Trqy = Trqy + fx(i)*zmom(i)-fz(i)*xmom(i)
        Trqz = Trqz + fy(i)*xmom(i)-fx(i)*ymom(i)

       End Do

c---
c compute the torque with respect to the particle center
c---

       Trcx = Trqx - Frcz*cy + Frcy*cz
       Trcy = Trqy - Frcx*cz + Frcz*cx
       Trcz = Trqz - Frcy*cx + Frcx*cy

       Resist(1,k) = cf*Frcx
       Resist(2,k) = cf*Frcy
       Resist(3,k) = cf*Frcz

       Resist(4,k) = cf*Trcx 
       Resist(5,k) = cf*Trcy
       Resist(6,k) = cf*Trcz

c      Resist(4,k) = cf*Trqx ! reset
c      Resist(5,k) = cf*Trqy
c      Resist(6,k) = cf*Trqz

      End Do

c---
c Display the grand resistance matrix
c---

      write (6,*)
      write (6,*) " prtcl_3d: Grand Resistance Matrix: "
      write (6,*) " ----------------------------------"
      write (6,*)

      Do k=1,6
        write (6,102) (Resist(k,md),md=1,nrhs)
      End Do

c---------------------
c compute nodal values
c of the traction
c---------------------

       mode = 6 ! pick a mode

       Do i=1,nelm
        inelm  = i+nelm
        inelm2 = inelm+nelm
        fx(i) = sln(i,mode)
        fy(i) = sln(inelm,mode)
        fz(i) = sln(inelm2,mode)
       End Do

       Do i=1,npts

        fx_node(i) = 0.0D0
        fy_node(i) = 0.0D0
        fz_node(i) = 0.0D0

        neigh_elm = ne(i,1)

        Do j=2,neigh_elm+1
          lelm = ne(i,j)
          fx_node(i) = fx_node(i)+fx(lelm)
          fy_node(i) = fy_node(i)+fy(lelm)
          fz_node(i) = fz_node(i)+fz(lelm)
        End Do

        fx_node(i) = fx_node(i)/neigh_elm
        fy_node(i) = fy_node(i)/neigh_elm
        fz_node(i) = fz_node(i)/neigh_elm

        fnormal = fx_node(i)*vna(i,1)
     +          + fy_node(i)*vna(i,2)
     +          + fz_node(i)*vna(i,3)

       ! shear stress:

       fx_node(i) = fx_node(i)-fnormal*vna(i,1)
       fy_node(i) = fy_node(i)-fnormal*vna(i,2)
       fz_node(i) = fz_node(i)-fnormal*vna(i,3)

       fshear(i) = Dsqrt(fx_node(i)**2+fy_node(i)**2
     +                  +fz_node(i)**2)
       End Do

c--------------------
c print the triangles
c and visualize the shear stress
c in a matlab file
c--------------------

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
      end if

      write (1,*) Iflow
      write (1,*) wall

      Do k=1,nelm
        call printel (k,index,fshear)
      End Do

c-----
c done
c-----

  99  Continue

      close (1)

  100 Format(1x,i4,10(1x,f12.5))
  101 Format(10(1x,f12.5))
  102 Format(10(1x,f10.6))
  103 Format(10(1x,f15.10))
  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)
  112 Format(" Centroid     :",3(F15.10))
  130 Format(" Determinant =",f20.10)
  200 Format(100(1x,f5.3))
  201 Format(10(1x,f7.5))
  205 Format(10(1x,f15.10))

  210 Format(1x,i4,10(1x,f7.3))

      stop
      end
