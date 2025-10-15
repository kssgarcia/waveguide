      subroutine prtcl_3d_mob_dlr_se_slv
     +
     +    (Ioctaicos
     +    ,ndiv
     +    ,boa,coa
     +    ,req
     +    ,cxp,cyp,czp
     +    ,phi1,phi2,phi3
     +    ,intm
     +    ,mint
     +    ,tol
     +    ,Nter
     +
     +    ,npts,nelm
     +
     +    ,dpl_dn
     +    ,Vx,Vy,Vz
     +    ,Ox,Oy,Oz
     +    )

c================================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c================================================

c------------------------------------------
c SYMBOLS:
c -------
c
c u(i,3):      incident velocity at the ith node
c dpl(i,3):    dipole weight at the ith node
c dlp_pv(i,3): principal value of the double-layer potential
c              at the ith node
c xsp(i,j):    x position of the jth node on the ith element
c---------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)

      Dimension       n(512,6), nbe(512,3)
      Dimension   alpha(512),  beta(512),gamma(512)
      Dimension    arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension  crvmel(512)

      Dimension xiq(20),etq(20),wq(20)

c spectral:

      Dimension vmaster(8),van(500,500),vaninv(500,500)

      Dimension    xsp(512,100),  ysp(512,100),  zsp(512,100)
      Dimension  xvnsp(512,100),yvnsp(512,100),zvnsp(512,100)
      Dimension    nsp(512,100)   ! connectivity

      Dimension    vnxsp(20000),vnysp(20000),vnzsp(20000)
      Dimension      psp(20000,3)   ! spectral node
      Dimension        u(20000,3)
      Dimension      dpl(20000,3)
      Dimension   dlp_pv(20000,3)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/var1/Iflow
      common/var2/Uinf,shrt
      common/var3/wall
      common/var4/Iflowinf

      common/geo1/arel
      common/geo2/vna
      common/geo6/crvmel

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c spectral:

      common/spectral1/mpoly,npoly,vmaster,van,vaninv
      common/spectral2/nesp,ngsp
      common/spectral5/psp,dpl

c----------
c constants
c----------

      null = 0
      oot  = 1.0D0/3.0D0

c----------------------------
c triangulate the unit sphere
c----------------------------

      if(Ioctaicos.eq.1) then

       call trgl6_octa
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      else if(Ioctaicos.eq.2) then

      call trgl6_icos
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      end if

      write (6,*)
      write (6,*) "number of surface nodes: ",npts
      write (6,*)
      write (6,*) "number of elements: ",nelm
      write (6,*)

c----------------------------
c expand to specified radius
c deform to an ellipsoid
c----------------------------

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

c-------------------------
c rotate by phi1,phi2,phi3
c around the x,y,z axes
c-------------------------

      cs = Dcos(phi1*pi)
      sn = Dsin(phi1*pi)

      Do i=1,npts                  ! rotate about the x axis
       tmpx = p(i,1)
       tmpy = cs*p(i,2)+sn*p(i,3)
       tmpz =-sn*p(i,2)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = Dcos(phi2*pi)
      sn = Dsin(phi2*pi)

      Do i=1,npts                   ! rotate about the y axis
       tmpx = cs*p(i,1)-sn*p(i,3)
       tmpy = p(i,2)
       tmpz = sn*p(i,1)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = Dcos(phi3*pi)
      sn = Dsin(phi3*pi)

      Do i=1,npts                  ! rotate about the z axis
       tmpx = cs*p(i,1)+sn*p(i,2)
       tmpy =-sn*p(i,1)+cs*p(i,2)
       tmpz = p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

c--------------------
c translate center to
c specified position
c--------------------

      Do i=1,npts 
        p(i,1) = p(i,1) + cxp
        p(i,2) = p(i,2) + cyp
        p(i,3) = p(i,3) + czp
      End Do

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
c vlm:  surface area of the individual elements
c xmon, ymom,zmom: x, y, and z surface moments
c                  over each element
c area:   total surface area and volume
c crvmel: mean curvature over each element
c vna:    average normal vector at the nodes
c------------------------------------------------

      call elm_geom
     +
     +   (nelm,npts,mint
     +   ,xmom,ymom,zmom
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

      area_red = area/(pi4*req**2)        ! scale
      vlm_red  = 3.0D0*vlm/(pi4*req**3)   ! scale

      write (6,*)
      write (6,110) area_red
      write (6,111) vlm_red
      write (6,*)
      write (6,800) cx,cy,cz

c-------------------------------------------------
c compute the collocation nodes over each triangle
c-------------------------------------------------

      Iopt = 2   ! interpolation option (need the normal vector)

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
         Ic = Ic+1

         l = mpoly+3-i-j

         xi = (1.0D0+2.0D0*vmaster(i)-vmaster(j)-vmaster(l))/3.0D0
         eta = (1.0D0-vmaster(i)+2.0D0*vmaster(j)-vmaster(l))/3.0D0

         call interp_p
     +
     +     (p(i1,1),p(i1,2),p(i1,3)
     +     ,p(i2,1),p(i2,2),p(i2,3)
     +     ,p(i3,1),p(i3,2),p(i3,3)
     +     ,p(i4,1),p(i4,2),p(i4,3)
     +     ,p(i5,1),p(i5,2),p(i5,3)
     +     ,p(i6,1),p(i6,2),p(i6,3)
     +
     +     ,al,be,ga
     +     ,xi,eta
     +
     +     ,xk,yk,zk
     +     ,DxDxik,DyDxik,DzDxik
     +     ,DxDetk,DyDetk,DzDetk
     +     ,vnxk,vnyk,vnzk
     +     ,hxik,hetk,hsk
     +     ,Iopt
     +     )

         xsp(k,Ic) = xk
         ysp(k,Ic) = yk
         zsp(k,Ic) = zk

         xvnsp(k,Ic) = vnxk
         yvnsp(k,Ic) = vnyk
         zvnsp(k,Ic) = vnzk

c        write (6,*) vnxU,vnyU,vnzU

         End Do
       End Do

      End Do

      nesp  = Ic     ! number of element spectral nodes
      nespt = Jc

      write (6,*) "number of element spectral nodes:",nesp
      write (6,*) "total of element spectral nodes:", nespt

c--------------------------------------------------
c Generate a list of unique global nodes by looping
c over all elements and adding nodes not found in the list
c
c Fill in the connectivity table nsp(i,j)
c containing the global labels of element points 1-Nesp
c
c Compute the normal vector at the unique global nodes
c-------------------------------------------------

c   first-element nodes are entered manually:

      Do j=1,nesp
         psp(j,1) =   xsp(1,j)
         psp(j,2) =   ysp(1,j)
         psp(j,3) =   zsp(1,j)
       vnxsp(j)   = xvnsp(1,j)
       vnysp(j)   = yvnsp(1,j)
       vnzsp(j)   = zvnsp(1,j)
       nsp(1,j)   = j
      End Do

      ngsp = nesp

      eps = 0.0000000001

      Do i=2,nelm        ! loop over element 2 to nelm
        Do j=1,nesp         ! loop over element nodes

        Iflag=0

         Do k=1,ngsp

          if(abs(xsp(i,j)-psp(k,1)).le.eps) then
           if(abs(ysp(i,j)-psp(k,2)).le.eps) then
            if(abs(zsp(i,j)-psp(k,3)).le.eps) then

             Iflag = 1    ! the node has been recorded previously
             nsp(i,j) = k   ! the jth local node of element i
                          ! is the kth global node

            end if
           end if
          end if

         End Do

         if(Iflag.eq.0) then     ! record the node

          ngsp = ngsp+1          ! one more global node

          psp(ngsp,1) = xsp(i,j)
          psp(ngsp,2) = ysp(i,j)
          psp(ngsp,3) = zsp(i,j)

          vnxsp(ngsp) = xvnsp(i,j)
          vnysp(ngsp) = yvnsp(i,j)
          vnzsp(ngsp) = zvnsp(i,j)

c         write (6,*) vnxsp(ngsp),vnysp(ngsp),vnzsp(ngsp)

          nsp(i,j) = ngsp   ! the jth local node of element i
                             ! is the new global node
         end if

       End Do
      End Do                      !  end of loop over elements

      write (6,*) "number of global spectral nodes:",ngsp

c--------------------------------
c specify the boundary conditions
c at the global spectral nodes
c--------------------------------

      Do i=1,ngsp
       u(i,1) = 0.0D0
       u(i,2) = 0.0D0
       u(i,3) = 0.0D0
       if(Iflowinf.eq.1) u(i,1) = Uinf
       if(Iflowinf.eq.2) u(i,2) = Uinf
       if(Iflowinf.eq.3) u(i,3) = Uinf
       if(Iflowinf.eq.4) u(i,2) = shrt*psp(i,1)
       if(Iflowinf.eq.5) u(i,3) = shrt*psp(i,1)
      End Do

c----------------------
c initialize the dipole
c to an arbitrary value
c----------------------

      Do i=1,ngsp
       dpl(i,1) = 1.30D0
       dpl(i,2) = 1.20D0
       dpl(i,3) = 1.40D0
       dpl(i,1) = 1.00D0
       dpl(i,2) = 1.00D0
       dpl(i,3) = 1.00D0
       dpl(i,1) = psp(i,1)
       dpl(i,2) = psp(i,2)
       dpl(i,3) = psp(i,3)
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

      call sdlp_3d
     +
     +   (npts
     +   ,nelm
     +   ,intm
     +   ,mint
     +   ,nsp,psp
     +
     +   ,cx,cy,cz
     +
     +   ,dpl
     +
     +   ,dpl_sintx
     +   ,dpl_sinty
     +   ,dpl_sintz
     +
     +   ,dpl_sinmx
     +   ,dpl_sinmy
     +   ,dpl_sinmz
     +
     +   ,dpl_dn
     +
     +   ,dlp_pv
     +   )

c----------------------------------------
c translational and rotational velocities
c----------------------------------------

      fc = - pi4/area

c     write (6,*) "hello1"
c     write (6,*) fc,dpl_sintx

      Vx = fc * dpl_sintx
      Vy = fc * dpl_sinty
      Vz = fc * dpl_sintz

      fc = - 1.5D0*(pi4/area)**2
      fc = - (pi4/area)**2

      Ox = fc * dpl_sinmx
      Oy = fc * dpl_sinmy
      Oz = fc * dpl_sinmz

c     write (6,*) "hello2"
c     write (6,*) Vx,Vy,Vz,qdn,Ox,Oy,Oz
c     pause

c--------------------------------------------
c update the dipole and compute the rms error
c--------------------------------------------

      error = 0.0D0
     
      Do i=1,ngsp

       savex = dpl(i,1)
       savey = dpl(i,2)
       savez = dpl(i,3)

       xx = psp(i,1)-cx
       yy = psp(i,2)-cy
       zz = psp(i,3)-cz

c      write (6,*) i, vnxsp(i)**2+vnysp(i)**2+vnzsp(i)**2

       dpl(i,1) = -dlp_pv(i,1)/pi4
     +                + vnxsp(i)*dpl_dn/area
     +                + Vx/pi4
     +                + (Oy*zz-Oz*yy)/pi4
     +                - u(i,1)/pi4

       dpl(i,2) = -dlp_pv(i,2)/pi4
     +                + vnysp(i)*dpl_dn/area
     +                + Vy/pi4
     +                + (Oz*xx-Ox*zz)/pi4
     +                - u(i,2)/pi4

       dpl(i,3) = -dlp_pv(i,3)/pi4
     +                + vnzsp(i)*dpl_dn/area
     +                + Vz/pi4 
     +                + (Ox*yy-Oy*xx)/pi4
     +                - u(i,3)/pi4

       error = error+(savex-dpl(i,1))**2
     +              +(savey-dpl(i,2))**2
     +              +(savez-dpl(i,3))**2

c      write (6,100) i,dlp_pv(i,1),dlp_pv(i,2),dlp_pv(i,3)
c      write (6,100) i,dpl(i,1),dpl(i,2),dpl(i,3)

      End Do

c     pause

c----------------
c end of updating
c----------------

      error = Dsqrt(error)/(3.0D0*ngsp)

      write (6,*)
c     write (6,457) Iter,error

      if(Iflow.eq.1) then
        write (6,456) Iter,dpl_dn,error,Vx,Vy,Vz,Ox,Oy,Oz
      else
        write (6,456) Iter,dpl_dn,error,Vx/cx,Vy/cx,Vz/cx
     +                ,2.0D0*Ox,2.0D0*Oy,2.0D0*Oz
      end if

 99   Continue

      if(error.lt.tol) Go to 711   ! solution found

      if(Iter.lt.Nter) Go to 71

      Iter = 0

      write (6,*)
      write (6,*) " One more iteration batch ?"
      write (6,*)
      write (6,*) " Enter 1 for Yes, 0 for No"
      write (6,*) "--------------------------"
      read  (5,*) more

      if(more.eq.1) Go to 71

c-------------------------
c done with the iterations
c-------------------------

 711   Continue

c--------
c display
c--------

c     write (6,*)
c     write (6,*) " Converged solution: x-position, dipole"
c     write (6,*)

      write (4,*) ngsp
c     write (6,*) ngsp

      Do i=1,ngsp
       write (4,100) i,psp(i,1),dpl(i,1),dpl(i,2),dpl(i,3)
c      write (6,100) i,psp(i,1),dpl(i,1),dpl(i,2),dpl(i,3)
      End Do

      write (4,*) null

c-----
c done
c-----

  100 Format (1x,i5,10(1x,f10.5))
  101 Format (f15.10)
  102 Format (10(1x,f12.8))

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)

  456 Format (" Iter: ",I3," Proj:",f12.8," Error: ",f12.10,/
     +       ," Vx: ",f12.8," Vy: ",f12.8," Vz: ",f12.8,/
     +       ," Ox: ",f12.8," Oy: ",f12.8," Oz: ",f12.8)
  457 Format (" prtcl_3ds_slv: Iter: ",I3," Error: ",f12.10)
  459 Format (" Flow across the boundary:",f15.10)

  800 Format(" centroid:",3(1x,f10.5))

      return
      end
