      program ldr_3d_ext_se

c===========================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c------------------------------------------
c Solution of the Dirichlet problem for 
c Laplace's equation in the exterior of
c a closed three-dimensional region
c
c This program solves the equation:
c
c Laplacian(f) = 0.0
c
c subject to boundary conditions for the function f
c at the (interior) closed boundary of the solution domain
c
c The method is based on the completed double-layer
c representation discussed by Pozrikidis (1997, p.516).
c
c SYMBOLS:
c -------
c
c f(i):         boundary value of f at the ith node
c
c gradn_f(i):   normal derivative of f at the ith node
c               gradn_f = n.grad(f)
c               where n is the outward normal vector
c
c dpl(i): dipole weight at the ith node
c
c dlp_pv(i): principal value of the double-layer potential
c---------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)

      Dimension      n(512,6), nbe(512,3)
      Dimension   alpha(512),  beta(512),gamma(512)
      Dimension    arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension  crvmel(512)

      Dimension xiq(20),etq(20),wq(20)

c spectral:

      Dimension vmaster(8),van(500,500),vaninv(500,500)
      Dimension  xsp(512,100),ysp(512,100),zsp(512,100)

      Dimension   nsp(512,100)   ! connectivity
      Dimension      psp(20000,3)
      Dimension        f(20000)
      Dimension      dpl(20000)
      Dimension   dlp_pv(20000)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/var1/Iflow,Ign,wall

      common/geo1/arel
      common/geo2/vna
      common/geo6/crvmel

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c spectral:

      common/spectral1/mpoly,Npoly,vmaster,van,vaninv
      common/spectral2/nesp,ngsp

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      piq = 0.25D0 *pi
      pih = 0.50D0 *pi
      pi2 = 2.00D0 *pi
      pi4 = 4.00D0 *pi
      pi6 = 6.00D0 *pi
      pi8 = 8.00D0 *pi
      spi4 = Dsqrt(pi4)

      null  = 0
      Nfour  = 4
      Nseven = 7
      oot  = 1.0D0/3.0D0

c-----------
c Input data
c-----------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to read data from file: ldr_3d_ext_se.dat'
      write (6,*) ' 2 to type data into the keyboard'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*)   Ienrd

      Ienrd = 1

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (2,file="ldr_3d_ext_se.dat")

       read (2,*) Ioctaicos
       read (2,*) ndiv
       read (2,*)
       read (2,*) boa,coa
       read (2,*) req
       read (2,*) cxp,cyp,czp
       read (2,*)
       read (2,*) intm
       read (2,*) mint
       read (2,*)
       read (2,*) beta_ps
       read (2,*)
       read (2,*) psx,psy,psz
       read (2,*)
       read (2,*) Nter
       read (2,*) tol
       read (2,*) 
       read (2,*) Iflow
       read (2,*) 
       read (2,*) wall
       read (2,*) Ign
       read (2,*) 
       read (2,*) mpoly
       read (2,*) Ipd

      close (2)

c---------
      else
c---------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to discretize an octahedron"
      write (6,*) " 2 to discretize an icosahedron"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ioctaicos

      if(Ioctaicos.eq.0) Go to 99

  91  Continue

      write (6,*)
      write (6,*) " Select the level of triangulation"
      write (6,*) " Choose from 0, 1, 2, 3; 99 to quit"
      write (6,*) " ----------------------------------"
      read  (5,*) ndiv

      if(ndiv.eq.99) Go to 99

      if((ndiv.lt.0).or.(lt.gt.3))  then
        write (6,*) 'out of range; please try again'
        Go to 91
      end if

      write (6,*)
      write (6,*) " The x,y,z semi-axes of the ellipsoid are: a,b,c"
      write (6,*) " Please enter the axes ratios b/a and c/a"
      write (6,*) " ----------------------------------------"
      read  (5,*) boa,coa

      write (6,*)
      write (6,*) " Please enter the eelipsoid equivalent radius"
      write (6,*) " --------------------------------------------"
      read  (5,*) req

      write (6,*)
      write (6,*) " Please enter the coordinates of the"
      write (6,*) " center of the ellipsoid"
      write (6,*) " -----------------------"
      read  (5,*) cxp,cyp,czp

      write (6,*)
      write (6,*) " Will use the m-point gauss-triangle rule"
      write (6,*) "      for integration over each triangle"
      write (6,*)
      write (6,*) " Please enter m"
      write (6,*)
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read  (5,*) mint

      write (6,*)
      write (6,*) " Please enter the order of the quadrature"
      write (6,*) " for integration of the singular triangles"
      write (6,*) " -----------------------------------------"
      read  (5,*) NGL

      write (6,*)
      write (6,*) " Please enter the point-source coefficient beta"
      write (6,*) " ----------------------------------------------"
      read  (5,*) beta_ps

      write (6,*)
      write (6,*) " Please enter the coordinates of the point source"
      write (6,*) " ------------------------------------------------"
      read  (5,*) psx,psy,psz

      write (6,*)
      write (6,*) " How many batch iterations in solving"
      write (6,*) " the integral equation ?"
      write (6,*) "------------------------"
      read  (5,*) Nter

      write (6,*)
      write (6,*) " Please enter the tolerance for stopping"
      write (6,*) " the iterations"
      write (6,*) " --------------"
      read  (5,*) tol

      write (6,*)
      write (6,*) " Choose the Green's function"
      write (6,*)
      write (6,*) " Enter 1 for free space"
      write (6,*) "       2 for domain bounded by a plane wall"
      write (6,*) " ------------------------------------------"
      read  (5,*) Iflow

      if(Iflow.eq.2) then

        write (6,*) " The wall is located at x=wall"
        write (6,*) " Please enter: wall"
        write (6,*) " ------------------"
        read  (5,*) wall

        write (6,*) " Enter 1 for the Green function"
        write (6,*) "       2 for the Neumann function"
        write (6,*) " --------------------------------"
        read  (5,*) Ign

      End If

c-----------
      end if
c-----------

c------------------
c open output files
c------------------

      open (4,file="ldr_3d_ext_se.out")

c----------------------------
c triangulate the unit sphere
c----------------------------

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

      write (6,*)
      write (6,*) "Number of points:   ",npts
      write (6,*)
      write (6,*) "Number of elements: ",nelm
      write (6,*)

c---
c Expand to specified radius;
c Deform to an ellipsoid;
c Translate center to specified position.
c---

      scale = req/(boa*coa)**oot

      Do i=1,npts
        p(i,1) = scale*    p(i,1) + cxp
        p(i,2) = scale*boa*p(i,2) + cyp
        p(i,3) = scale*coa*p(i,3) + czp
      End Do

c---------------------------------
c specify the boundary values of f
c---------------------------------

      write (6,*)
      write (6,*) 'CHOOSE THE BOUNDARY DATA'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for f = 1.0'
      write (6,*) ' 2 for f = x'
      write (6,*) ' 3 for f = x**2'
      write (6,*) ' 4 for f = sin(x)'
      write (6,*) ' 5 for f = x * cos(x)'
      write (6,*) ' 6 for f = exp(2x)'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) Ibc

      if(Ibc.eq.0) Go to 99

c-----------------------------
c read the triangle quadrature
c-----------------------------

      call gauss_trgl
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

c----------------------------
c compute alpha, beta, gamma,
c for the xi-eta mapping
c over each triangle
c---------------------------

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
c     vlm:  surface area of the individual elements
c     xmon, ymom,zmom: x, y, and z surface moments
c                      over each element
c     area:   total surface area and volume
c     crvmel: mean curvature over each element
c     vna:    average normal vector at the nodes
c------------------------------------------------

      call elm_geom
     +
     +   (nelm,npts
     +   ,mint
     +   ,xmom,ymom,zmom
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c---
c normalize surface area and volume
c---

      area_red = area/(pi4*req**2)        ! normalize
      vlm_red  = 3.0D0*vlm /(pi4*req**3)    ! normalize

      write (6,*)
      write (6,110) area_red
      write (6,111) vlm_red
      write (6,*)
      write (6,800) cx,cy,cz

c----------------------------
c Define the collocation nodes
c over the standard triangle
c and compute the inverse of
c the vandermonde matrix
c
c Ipd = 0 uniform distribution
c       1 spectral distribution
c---------------------------

       call vander (Ipd)

c      write (6,*) Npoly
c      Do i=1,Npoly
c        write (6,102) (van(i,j),j=1,Npoly)
c      End Do
c      write (6,*) Npoly
c      Do i=1,Npoly
c        write (6,102) (vaninv(i,j),j=1,Npoly)
c      End Do
c      pause

c-------------------------------------------------
c compute the collocation nodes over each triangle
c-------------------------------------------------

      Iopt_sp = 1   ! interpolation option

      Jc = 0

      Do k=1,nelm    ! loop over elements

       i1 = n(k,1)     ! global node labels
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       al = alpha(k)
       be = beta (k)
       ga = gamma(k)

       Ic = 0

       Do i=1,mpoly+1
         Do j=1,mpoly+2-i

         Jc = Jc +1

         l = mpoly+3-i-j

         xi  = (1.0D0+ 2.0D0*vmaster(i)  -vmaster(j)-vmaster(l))/3.0D0
         eta = (1.0D0-   vmaster(i)+2.0D0*vmaster(j)-vmaster(l))/3.0D0

         call interp_p
     +
     +     (Iopt_sp
     +
     +     ,p(i1,1),p(i1,2),p(i1,3)
     +     ,p(i2,1),p(i2,2),p(i2,3)
     +     ,p(i3,1),p(i3,2),p(i3,3)
     +     ,p(i4,1),p(i4,2),p(i4,3)
     +     ,p(i5,1),p(i5,2),p(i5,3)
     +     ,p(i6,1),p(i6,2),p(i6,3)
     +
     +     ,al,be,ga
     +     ,xi,eta
     +     ,x,y,z
     +
     +     ,DxDxiU,DyDxiU,DzDxiU
     +     ,DxDetU,DyDet,DzDetU
     +     ,vnxU,vnyU,vnzU
     +     ,hxiU,hetU,hsU
     +     )

         Ic = Ic+1

         xsp(k,Ic) = x
         ysp(k,Ic) = y
         zsp(k,Ic) = z

         End Do
       End Do

      End Do

      nesp = Ic
      nespt = Jc

      write (6,*) "number of element spectral nodes:",nesp
      write (6,*) "total of element spectral nodes:", nespt

c--------------------------------------------------
c Generate a list of unique global nodes by looping
c over all elements
c and adding nodes not found in the list.
c
c Fill in the connectivity table nsp(i,j)
c containing the global labels of element points 1-Nesp
c-------------------------------------------------

      Do i=1,nesp
       psp(i,1) = xsp(1,i)
       psp(i,2) = ysp(1,i)
       psp(i,3) = zsp(1,i)
       nsp(1,i) = i
      End Do

      ngsp = nesp

      eps = 0.0000000001

      Do i=2,nelm        ! loop over elements
        Do j=1,nesp         ! loop over element nodes

        Iflag=0

         Do k=1,ngsp

          if(abs(xsp(i,j)-psp(k,1)).le.eps) then
           if(abs(ysp(i,j)-psp(k,2)).le.eps) then
            if(abs(zsp(i,j)-psp(k,3)).le.eps) then

             Iflag = 1    ! the node has been recorded previously
             nsp(i,j) = k   ! the jth local node of element i
                          ! is the kth global node

            End If
           End If
          End If
         End Do

         if(Iflag.eq.0) then     ! record the node

          ngsp = ngsp+1          ! one more global node

          psp(ngsp,1) = xsp(i,j)
          psp(ngsp,2) = ysp(i,j)
          psp(ngsp,3) = zsp(i,j)

          nsp(i,j) = ngsp   ! the jth local node of element i
                             ! is the new global node
         end if

       End Do
      End Do                      !  end of loop over elements

      write (6,*) "number of global spectral nodes:",ngsp

c---------------------------------
c specify the boundary conditions
c---------------------------------

      Do i=1,ngsp
       if(Ibc.eq.1) f(i) = 1.0D0
       if(Ibc.eq.2) f(i) = psp(i,1)
       if(Ibc.eq.3) f(i) = psp(i,1)**2
       if(Ibc.eq.4) f(i) = Dsin( psp(i,1)*pi )
       if(Ibc.eq.5) f(i) = psp(i,1)*Dcos( psp(i,1)*pi )
       if(Ibc.eq.6) f(i) = Dexp(2.0D0*psp(i,1)) !  * psp(i,2)
      End Do

c----------------------
c initialize the dipole
c to an arbitrary distribution
c----------------------

      Do i=1,ngsp
       dpl(i) = 1.0D0
c      dpl(i) = 0.1d0
      End Do

c--------------------------------------------
c iterative solution of the integral equation
c for the dipole strength
c--------------------------------------------

      Iter = 0

 71   Continue

      Iter = Iter+1

c-----------------------------------
c compute the double-layer potential
c-----------------------------------

      call ldlp_3d_se
     +
     +   (npts
     +   ,nelm
     +   ,intm
     +   ,mint
     +   ,nsp,psp
     +   ,dpl
     +   ,dpl_sint
     +   ,dlp_pv
     +   )

c---
c strength of the point source
c---

      pss = -beta_ps * spi4 * dpl_sint/dsqrt(area)

c--------------------------------------------
c update the dipole and compute the rms error
c--------------------------------------------

      error = 0.0D0
     
      Do i=1,ngsp

c---
c Compute the Green's function
c for the point source
c---

       Iopt_gf = 1                 ! need only G

c---------
       if(Iflow.eq.1) then
c---------

       call lgf_3d_fs
     +
     +    (Iopt_gf
     +    ,psp(i,1),psp(i,2),psp(i,3)
     +    ,psx,psy,psz
     +    ,Gps
     +    ,Gpsx,Gpsy,Gpsz
     +    )

c-----
      else if(Iflow.eq.2) then   !  wall at x=wall
c-----

       call lgf_3d_w
     +
     +   (Iopt_gf
     +   ,Ign
     +   ,psp(i,1),psp(i,2),psp(i,3)
     +   ,psx,psy,psz
     +   ,wall
     +   ,Gps
     +   ,Gpsx,Gpsy,Gpsz
     +   )

c---------
       end if
c---------

       save = dpl(i)

       dpl(i) = 2.0D0*(                  ! update
     +                  f(i)
     +                 -dlp_pv(i)
     +                 +pss*Gps
     +                 )

       error = error + (save-dpl(i))**2

c      write (6,100) i,f(i),dpl(i),dlp_pv(i)

      End Do

c----------------
c end of updating
c----------------

      error = Dsqrt(error)/ngsp

      write (6,456) Iter,dpl_sint,error,pss/pi4

      if(error.lt.tol) Go to 711   ! solution found

      if(Iter.lt.Nter) Go to 71


      write (6,*)
      write (6,*) " One more iteration batch ?"
      write (6,*)
      write (6,*) " Enter 1 for Yes, 0 for No"
      write (6,*) "--------------------------"
      read  (5,*) more

      if(more.eq.1) then
        Go to 71
        Iter = 0
      end if

c-------------------------
c Done with the iterations
c-------------------------

 711   Continue

c-----
c strength of the point source
c-----

      pss = -beta_ps * spi4 * dpl_sint/dsqrt(area)

c--------
c recodring session
c--------

c     write (6,*)
c     write (6,*) " Converged solution: x-position, dipole"
c     write (6,*)

      write (4,*) ngsp
c     write (6,*) ngsp

      Do i=1,ngsp
       write (4,100) i,psp(i,1),dpl(i)
c      write (6,100) i,psp(i,1),dpl(i)
      End Do

      write (4,*) null

      write (4,*) ngsp
c     write (6,*) ngsp

      Do i=1,ngsp
        write (4,102) psp(i,1),psp(i,2),psp(i,3)
c       write (6,102) psp(i,1),psp(i,2),psp(i,3)
      End Do

c     write (1,*) ntespn
c     Do k=1,nelm
c      Do i=1,nesp
c       write (1,102) xsp(k,i),ysp(k,i),zsp(k,i)
c      End Do
c     End Do
c     write (1,*) 2*nesp
c     Do i=1,nesp
c      write (1,102) xsp(1,i),ysp(1,i),zsp(1,i)
c     End Do
c     Do i=1,nesp
c      write (1,102) xsp(2,i),ysp(2,i),zsp(2,i)
c     End Do

c--------------------
c print all triangles
c--------------------

      open (1,file="ldr_3d_ext_se.net")

      Index = 1  ! 6-node triangles
      Index = 2  ! 3-node triangles

      if(Index.eq.1) then          ! 6-node triangles
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*)   nelm
      else if(Index.eq.2) then     ! 3-node triangles
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*)  4*nelm
      end if

      write (1,*) Iflow
      write (1,*) wall

      Do k=1,nelm
       call printel (k,Index,p(:,1)) 
      End Do

      close (1)

c-----
c Done
c-----

 99   Continue

      close (4)

c------------
c Really Done
c------------

  100 Format (1x,i5,10(1x,f10.5))
  101 Format (f15.10)
  102 Format (10(1x,f12.8))

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)

  456 Format ("Iter: ",I3," Proj:",f12.8," Error: ",f12.10
     +       ," Point source: ",f12.8)
  459 Format ("Flow across the boundary:",f15.10)

  800 Format(" centroid:",3(1x,f10.5))

      stop
      end
