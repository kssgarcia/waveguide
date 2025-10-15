      program bump_3d 

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c-----------------------------------------------
c Shear flow past a spherical bump on a plane wall
c located at y=wall
c
c SYMBOLS:
c --------
c
c  p(i,j)	coordinates of the ith node (j=1,2,3)
c
c  ne(k,j)  ne(k,1) is the number of elements adjacent to point k
c           ne(k,2), ... are the elements label, j = 2, ..., 7
c           for this triangulation, up to six
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
      Dimension rhs(1600),sln(1600)    ! RHS and solution vectors

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

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.500*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi
      pi6 = 6.0D0*pi
      pi8 = 8.0D0*pi

      Null = 0
      None = 1
      Nfour = 4
      Nseven = 7

      oot  = 1.0D0/3.0D0

c------------
c preferences
c------------

      open (2,file="bump_3d.dat")

      read (2,*) ndiv 
      read (2,*) radius
      read (2,*) cont_angle
      read (2,*) mint
      read (2,*) NGL
      read (2,*) visc
      read (2,*) shrt
      read (2,*) wall

c--------------------------------
c boundary element discretization
c of the hemi-sphere
c--------------------------------

      call trgl6_octa_ss 
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      write (6,*)
      write (6,*) " bump_3d: number of elements: ",nelm
      write (6,*) " bump_3d: number of nodes   : ",npts
      write (6,*)

c     Do i=1,npts
c      neigh_elm = ne(i,1)
c      write (6,209) i,neigh_elm,(ne(i,j),j=2,neigh_elm+1)
c     End Do

c----------------------------------------------
c Expand to spedified radius
c Slide the surface grid to specified contact angle
c Translate center to specified position: y=wall
c----------------------------------------------

      Do i=1,npts
       p(i,1) = radius*p(i,1)
       p(i,2) = radius*p(i,2)
       p(i,3) = radius*p(i,3)
      End Do

      phi   = pi*cont_angle
      cs    = Dcos(phi)
      shift = radius*Dsin(pih-phi)

      Do i=1,npts
       rr    = Dsqrt(p(i,1)**2+p(i,3)**2)
       thet0 = Datan2(rr,p(i,2))
       thetz = Datan2(p(i,1),p(i,3))
       thet  = thet0*phi/pih
       rrr   = radius*Dsin(thet)
       p(i,1) = rrr   *Dsin(thetz)
       p(i,2) = radius*Dcos(thet)
       p(i,3) = rrr   *Dcos(thetz)
       p(i,2) = p(i,2)-shift
      End Do

      defx = 0.2D0
      defz = 0.0D0

      Do i=1,npts
       p(i,1) = p(i,1) * (1.0D0+defx*p(i,1))
       p(i,3) = p(i,3)
       write (6,*) p(i,1)
      End Do

      Do i=1,npts
       p(i,2) = p(i,2) + wall
      End Do

c---------------
c prepare to run
c---------------

      nelm2 = 2*nelm

c---------------------
c read the quadratures
c---------------------

      call gauss_leg (NGL,zz,ww)

      call gauss_trgl (mint,xiq,etq,wq)

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
c   The surface area of the individual elements
c   x, y, and z moments over each element
c
c   The particle surface area and volume
c
c   The mean curvature of each element: crvmel
c
c   The averaged element normal vector
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

      srf_area_n = srf_area/(pi4*radius**2)
      prt_vlm_n  = prt_vlm /(pi4*radius**3/3.0D0)

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
c  read  (5,*) Iproc
c  If(Iproc.eq.0) Go to 99
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
      write (6,*) " bump_3d: properties of collocation points computed"
      write (6,*)

c------------------------------
c record the collocation points
c------------------------------

      open (3,file="particle.points")

        write (3,*) nelm

        Do i=1,nelm
         write (3,103) x0(i),y0(i),z0(i)
        End Do

      close (3)

c-------------------------------
c Generate the influence matrix,
c three rows at a time,
c corresponding to the x, y, z 
c components of the integral equation
c------------------------------------

      cf = -1.0D0/(pi8*visc)

      Do i=1,nelm            ! run over the collocation points

       write (6,*) " bump_3d: collocating at point: ",i

       inelm  = i+nelm
       inelm2 = i+nelm+nelm

       Do j=1,nelm            ! run over the elements

         jnelm  = j+nelm
         jnelm2 = j+nelm+nelm

c-----------------------
          if(i.ne.j) then     ! regular element
c-----------------------
      
          call bump_3d_slp
     +
     +     (x0(i),y0(i),z0(i)
     +     ,j
     +     ,GExx,GExy,GExz
     +     ,GEyx,GEyy,GEyz
     +     ,GEzx,GEzy,GEzz
     +     ,mint
     +     )

c-------------
          else               ! singular element
c-------------

          call bump_3d_slp_sing
     +
     +      (x0(i),y0(i),z0(i)
     +      ,j
     +      ,GExx,GExy,GExz
     +      ,GEyx,GEyy,GEyz
     +      ,GEzx,GEzy,GEzz
     +      ,NGL
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

c--------

        End Do

      End Do
      
c------------
c system size
c------------

      mls = 3*nelm  ! size of the linear system emerging from collocation

c     Do i=1,mls
c      write (6,*) (RM(i,j),j=1,mls)
c     End Do

c------------------------
c Set the right-hand side
c------------------------

       Do i=1,mls
         rhs(i) = 0.0D0
       End Do

       Do i=1,nelm
        inelm = i+nelm
        rhs(i) = -shrt*(y0(i)-wall)
       End Do

c------------------------
c solve the linear system
c------------------------

c------------
c display
c
c     Do i=1,mls
c      write (6,200) (RM(i,j),j=1,Mdim),rhs(i)
c     End Do
c     stop
c------------

      Isymg  = 0    ! system is not symmetric
      Iwlpvt = 1    ! pivoting enabled

      write (6,*) " bump_3d: system size: ",mls

c---------------
      write (6,*) 
      write (6,*) " Compute the inverse"
      write (6,*) " of the influence matrix ?"
      write (6,*) 
      write (6,*) " Enter 0 for NO, 1 for YES"
      write (6,*) " -------------------------"
      write (6,*) 
      read  (5,*) Invert

      If(Invert.eq.1) then

      write (6,*) " bump_3d: computing the matrix inverse"

      call gel_inv
     +
     +  (mls
     +  ,RM
     +  ,RMinv
     +  ,Isym,Iwlpvt
c    +  ,l,u
     +  ,det
     +  ,Istop
     +  )

      write (6,*) " bump_3d: done"

      End If
c-----------

      write (6,*) " bump_3d: solving the linear system"

      call gel
     +
     +   (mls
     +   ,RM
     +   ,rhs
     +   ,sln
     +   ,Isymg
     +   ,Iwlpvt
     +   ,det
     +   ,Istop
     +   )

c     write (6,*) "Determinant = ",det

c--------------------------------
c assign the element traction
c
c compute force, torque, 
c and the grand resistance matrix
c--------------------------------

       Frcx = 0.0D0
       Frcy = 0.0D0
       Frcz = 0.0D0

       Trqx = 0.0D0
       Trqy = 0.0D0
       Trqz = 0.0D0

       write (6,*)
       write (6,*) " bump_3d: traction:"
       write (6,*)

       Do i=1,nelm

        inelm  = i+nelm
        inelm2 = inelm+nelm

        fx(i) = sln(i)
        fy(i) = sln(inelm)
        fz(i) = sln(inelm2)

        write (6,210) i,fx(i),fy(i),fz(i)

        Frcx = Frcx + fx(i)*arel(i)
        Frcy = Frcy + fy(i)*arel(i)
        Frcz = Frcz + fz(i)*arel(i)

        Trqx = Trqx + fz(i)*ymom(i)-fy(i)*zmom(i)
        Trqy = Trqy + fx(i)*zmom(i)-fz(i)*xmom(i)
        Trqz = Trqz + fy(i)*xmom(i)-fx(i)*ymom(i)

       End Do

       write (6,*)
       write (6,*) "Force:"
       write (6,*)

       fc = pi*radius**2
       write (6,101) Frcx/fc,Frcy/fc,Frcz/fc

       write (6,*)
       write (6,*) "Torque:"
       write (6,*)

       fc = pi*radius**3
       write (6,101) Trqx/fc,Trqy/fc,Trqz/fc

c---------------------
c compute nodal values of the traction
c by interpolation from element values
c---------------------

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

c-------------------------------
c display the triangles
c and visualize the shear stress
c-------------------------------

      open (1,file="bump_3d.net")

      index = 1  ! 6-node triangles
      index = 2  ! 3-node triangles arising by subdivision

      if(index.eq.1) then
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      elseif(index.eq.2) then
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end if

      Do k=1,nelm
        call printel(k,index,fshear)
      End Do

      write (1,*) wall

c-------------------------
c print collocation points
c in file: bump_3d.net
c-------------------------

      write (1,*) nelm

      Do i=1,nelm
        write (1,103) x0(i),y0(i),z0(i)
      End Do

c-----------------
c draw streamlines
c using RK2
c-----------------

      cf = -1.0D0/(pi8*visc)
      Nstep = 128
      Dt = 0.10D0
      Dth = 0.5D0*Dt

      read (2,*)
      write (6,*)
      Icount = 0

  91  Continue

      Icount = Icount+1

      read (2,*) xev,yev,zev  ! starting point
      If(xev.eq.99) Go to 99  ! stopping check

      write (6,*) " bump_3d: streamline: ",Icount

      write (1,*) Nstep+1
      write (1,101) xev,yev,zev
c     write (6,101) xev,yev,zev

c---
c first velocity evaluation
c---

      Do Istep=1,Nstep    ! step along a streamline

      velx = 0.0D0
      vely = 0.0D0
      velz = 0.0D0

      Do i=1,nelm            ! run over the elements

          call bump_3d_slp
     +
     +     (xev,yev,zev
     +     ,i
     +     ,GExx,GExy,GExz
     +     ,GEyx,GEyy,GEyz
     +     ,GEzx,GEzy,GEzz
     +     ,mint
     +     )

         velx = velx +GExx*fx(i) +GEyx*fy(i) +GEzx*fz(i)
         vely = vely +GExy*fx(i) +GEyy*fy(i) +GEzy*fz(i)
         velz = velz +GExz*fx(i) +GEyz*fy(i) +GEzz*fz(i)

c        write (6,*) velx,vely,velz

       End Do

       velx = velx*cf+shrt*(yev-wall)
       vely = vely*cf
       velz = velz*cf

c---
c save and advance
c---

        xev_save = xev
        yev_save = yev
        zev_save = zev

        velx_save = velx
        vely_save = vely
        velz_save = velz

        xev = xev+velx*Dt
        yev = yev+vely*Dt
        zev = zev+velz*Dt

c---
c second velocity evaluation
c---

      velx = 0.0D0
      vely = 0.0D0
      velz = 0.0D0

      Do i=1,nelm            ! run over the elements

         call bump_3d_slp
     +
     +   (xev,yev,zev
     +   ,i
     +   ,GExx,GExy,GExz
     +   ,GEyx,GEyy,GEyz
     +   ,GEzx,GEzy,GEzz
     +   ,mint
     +   )

         velx = velx + GExx*fx(i) + GEyx*fy(i) + GEzx*fz(i)
         vely = vely + GExy*fx(i) + GEyy*fy(i) + GEzy*fz(i)
         velz = velz + GExz*fx(i) + GEyz*fy(i) + GEzz*fz(i)

c        write (6,*) velx,vely,velz

        End Do

        velx = velx*cf+shrt*(yev-wall)
        vely = vely*cf
        velz = velz*cf

c----
c second step
c----
        xev = xev_save+(velx+velx_save)*Dth
        yev = yev_save+(vely+vely_save)*Dth
        zev = zev_save+(velz+velz_save)*Dth

        write (1,101) xev,yev,zev
c       write (6,101) xev,yev,zev

      End Do     ! step along a streamline

      Go to 91   ! another streamline

c-----
c done
c-----

  99  Continue

      write (1,*) Null

      close (1)
      close (2)

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
  209 Format(I4,100(1x,I3))

  210 Format(1x,i4,10(1x,f7.3))

      Stop
      End
