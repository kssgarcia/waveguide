      subroutine susp_3d_slv
     +
     +   (Ioctaicos
     +   ,ndiv
     +   ,boa,coa
     +   ,req
     +   ,Isize
     +   ,cx,cy,cz
     +   ,theta,phi
     +
     +   ,mpoly
     +   ,intm   ! integration method
     +   ,mint   ! quadrature order
     +   ,tol
     +   ,Nter
     +
     +   ,uclp
     +
     +   ,Iadhere
     +   ,torque,adhesion
     +
     +   ,npts,nelm
     +   ,points,nn
     +   ,area,vlm
     +   ,necl,ngcl
     +   ,pcl
     +   ,dpl
     +   ,dpl_dn
     +   ,stresslet
     +   ,Vx,Vy,Vz
     +   ,Ox,Oy,Oz
     +
     +   ,xforce_adh,yforce_adh,zforce_adh
     +   )

c==========================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------
c
c SYMBOLS:
c -------
c
c    ucl(i,3):  collocation point velocity
c   uclp(i,3):  particle-interaction collocation point velocity
c  dplcl(i,3): dipole weight at the ith node
c dlp_pv(i,3): principal value of the double-layer potential
c---------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3),points(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)

      Dimension      n(512,6), nn(512,6), nbe(512,3)
      Dimension   alpha(512),  beta(512),gamma(512)
      Dimension    arel(512),  xmom(512),  ymom(512),zmom(512)

      Dimension stresslet(3,3)

      Dimension xiq(20),etq(20),wq(20)

c collocation:

      Dimension vmaster(8)
      Dimension  xcl(512,100),ycl(512,100),zcl(512,100)
      Dimension  xvncl(512,100),yvncl(512,100),zvncl(512,100)

      Dimension   ncl(512,100)   ! connectivity

      Dimension     pcl(2306,3)   ! nodes
      Dimension     ucl(2306,3)
      Dimension uclsave(2306,3)
      Dimension    uclp(2306,3)
      Dimension     dpl(2306,3)
      Dimension  dlp_pv(2306,3)
      Dimension   vnxcl(2306),vnycl(2306),vnzcl(2306)

      Dimension     dpl1(2306,3)
      Dimension     dpl2(2306,3)
      Dimension     dpl3(2306,3)
      Dimension     dpl4(2306,3)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/conncl/ncl

      common/albega/alpha,beta,gamma

      common/var1/Iflow,Iwall
      common/var2/Uinf,shrt
      common/var3/wall
      common/var4/visc

      common/var9/dpl1,dpl2,dpl3,dpl4

      common/geo1/arel
      common/geo2/vna

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c collocation:

      common/spectral1/vmaster

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

c     write (6,*) " susp_3d_slv: number of points:   ",npts
c     write (6,*) " susp_3d_slv: number of elements: ",nelm

      Do i=1,nelm
        Do j=1,6
          nn(i,j) = n(i,j)
        End Do
      End Do

c     Do i=1,npts
c       write (6,100) i,p(i,1),p(i,2),p(i,3)
c     End Do

c----------------------------
c expand to specified radius.
c deform to an ellipsoid
c----------------------------

      scale = req/(boa*coa)**oot
                                                       
      x_axis = scale
      y_axis = scale*boa
      z_axis = scale*coa

      if(Isize.eq.1) then
        axis_max = x_axis
        if(y_axis.gt.axis_max) axis_max = y_axis
        if(z_axis.gt.axis_max) axis_max = z_axis
        x_axis = x_axis/axis_max
        y_axis = y_axis/axis_max
        z_axis = z_axis/axis_max
      end if

      Do i=1,npts                  ! scale
        p(i,1) = x_axis*p(i,1)
        p(i,2) = y_axis*p(i,2)
        p(i,3) = z_axis*p(i,3)
      End Do

c     write (6,*) " susp_3d_slv: x semi-axis = ",x_axis
c     write (6,*) " susp_3d_slv: y semi-axis = ",y_axis
c     write (6,*) " susp_3d_slv: z semi-axis = ",z_axis

c----------------------------------
c rotate by theta around the z axis
c----------------------------------

      cs = Dcos(theta)
      sn = Dsin(theta)

      Do i=1,npts                
       tmpx = cs*p(i,1)-sn*p(i,2) 
       tmpy = sn*p(i,1)+cs*p(i,2)
       tmpz = p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

c--------------------------------
c rotate by phi around the x axis
c--------------------------------

      cs = Dcos(phi)
      sn = Dsin(phi)

      Do i=1,npts        
       tmpx = p(i,1)
       tmpy = cs*p(i,2)-sn*p(i,3)
       tmpz = sn*p(i,2)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
c       write (6,100) i,p(i,1),p(i,2),p(i,3),phi
      End Do

c--------------------
c translate to
c specified position
c--------------------

      Do i=1,npts 
        p(i,1) = p(i,1) + cx
        p(i,2) = p(i,2) + cy
        p(i,3) = p(i,3) + cz
        points(i,1) = p(i,1)
        points(i,2) = p(i,2)
        points(i,3) = p(i,3)
      End Do

c----------------------------
c Compute alpha, beta, gamma,
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

c------------------------------
c compute:
c
c     vlm:  the particle volume
c
c------------------------------

      call elm_geom
     +
     +   (nelm,npts,mint
     +   ,xmom,ymom,zmom
     +   ,area,vlm
     +   ,cx,cy,cz
     +   ,cx1,cy1,cz1
     +   ,smomin_xx,smomin_xy,smomin_xz
     +   ,smomin_yx,smomin_yy,smomin_yz
     +   ,smomin_zx,smomin_zy,smomin_zz
     +   )

c     area_red = area/(pi4*req**2)        ! normalize
c     vlm_red  = 3.0D0*vlm /(pi4*req**3)    ! normalize
c     write (6,110) area_red
c     write (6,111) vlm_red
c     write (6,800) cx1,cy1,cz1

c     write (6,102) smomin_xx,smomin_xy,smomin_xz
c     write (6,102) smomin_yx,smomin_yy,smomin_yz
c     write (6,102) smomin_zx,smomin_zy,smomin_zz

c------------------------------
c compute the collocation nodes
c over each triangle
c------------------------------

      Iopt_col = 2   ! interpolation option (Need the normal vector)

      Jc = 0        ! total element node counter

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

       Ic = 0         ! element node counter

       Do i=1,mpoly+1
         Do j=1,mpoly+2-i

         Jc = Jc +1

         l = mpoly+3-i-j

         xi  = (1.0D0+ 2.0D0*vmaster(i)  -vmaster(j)-vmaster(l))/3.0D0
         eta = (1.0D0-   vmaster(i)+2.0D0*vmaster(j)-vmaster(l))/3.0D0

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
     +     ,x,y,z
     +
     +     ,DxDxiU,DyDxiU,DzDxiU
     +     ,DxDetU,DyDet,DzDetU
     +     ,vnxU,vnyU,vnzU
     +     ,hxiU,hetU,hsU
     +     ,Iopt_col
     +     )

         Ic = Ic+1

         xcl(k,Ic) = x
         ycl(k,Ic) = y
         zcl(k,Ic) = z

         xvncl(k,Ic) = vnxU
         yvncl(k,Ic) = vnyU
         zvncl(k,Ic) = vnzU

c        write (6,*) vnxU,vnyU,vnzU
c        write (6,100) k,xcl(k,Ic),ycl(k,Ic),zcl(k,Ic)

         End Do
       End Do

      End Do

      necl = Ic
      neclt = Jc

      write (6,*) " susp_3d_slv: No of elm col nodes:",necl
      write (6,*) " susp_3d_slv: Total of elm cl nodes:", neclt

c--------------------------------------------------
c Generate a list of unique global nodes by looping
c over all elements
c and adding nodes not found in the list.
c
c Fill in the connectivity table ncl(i,j)
c containing the global labels of element points 1-Necl
c
c Compute the normal vector at the unique global nodes
c-------------------------------------------------

      Do j=1,necl
            pcl(j,1) =   xcl(1,j)
            pcl(j,2) =   ycl(1,j)
            pcl(j,3) =   zcl(1,j)
          vnxcl(j)   = xvncl(1,j)
          vnycl(j)   = yvncl(1,j)
          vnzcl(j)   = zvncl(1,j)
       ncl(1,j)      = j
      End Do

      ngcl = necl

      eps = 0.00000001

      Do i=2,nelm        ! loop over elements
        Do j=1,necl         ! loop over element nodes

        Iflag=0

         Do k=1,ngcl
          if(abs(xcl(i,j)-pcl(k,1)).le.eps) then
           if(abs(ycl(i,j)-pcl(k,2)).le.eps) then
            if(abs(zcl(i,j)-pcl(k,3)).le.eps) then

             Iflag = 1    ! the node has been recorded previously
             ncl(i,j) = k   ! the jth local node of element i
                          ! is the kth global node

            end if
           end if
          end if
         End Do

         if(Iflag.eq.0) then     ! record the node

          ngcl = ngcl+1          ! one more global node

          pcl(ngcl,1) = xcl(i,j)
          pcl(ngcl,2) = ycl(i,j)
          pcl(ngcl,3) = zcl(i,j)

          vnxcl(ngcl) = xvncl(i,j)
          vnycl(ngcl) = yvncl(i,j)
          vnzcl(ngcl) = zvncl(i,j)

c         write (6,*) vnxcl(ngcl),vnycl(ngcl),vnzcl(ngcl)

          ncl(i,j) = ngcl   ! the jth local node of element i
                            ! is the new global node
         End If
       End Do
      End Do                      !  end of loop over elements

      write (6,*) " susp_3d_slv: num of global colloc nodes:",ngcl

c-----------------
c printing session
c-----------------

c     Do i=1,nelm
c     write (6,*) (ncl(i,j),j=1,necl)
c     End Do

c     Do i=1,ngcl
c     write (6,*) pcl(i,1),pcl(i,2),pcl(i,3)
c     End Do

c--------------------------------
c specify the boundary conditions
c at the global spectral nodes
c--------------------------------

      Do i=1,ngcl
       ucl(i,1) = uclp(i,1)      ! mutual particle interaction
       ucl(i,2) = uclp(i,2)
       ucl(i,3) = uclp(i,3)
       if(Iflow.eq.1)  ucl(i,1) = ucl(i,1)+Uinf
       if(Iflow.eq.2)  ucl(i,2) = ucl(i,2)+Uinf
       if(Iflow.eq.3)  ucl(i,3) = ucl(i,3)+Uinf
       if(Iflow.eq.4)  ucl(i,1) = ucl(i,1)+shrt*pcl(i,2)
       if(Iflow.eq.5)  ucl(i,3) = ucl(i,3)+shrt*pcl(i,2)
       if(Iflow.eq.6)  ucl(i,2) = ucl(i,2)+shrt*pcl(i,1)
       if(Iflow.eq.10) ucl(i,1) = ucl(i,1)+shrt*(pcl(i,2)-wall)
       if(Iflow.eq.20) ucl(i,1) = ucl(i,1)+shrt*pcl(i,2)
      End Do

c---------------------------------
c flow due to wall adhesion force
c applied at the particle centroid
c---------------------------------

      Do i=1,ngcl
       uclsave(i,1) = ucl(i,1)
       uclsave(i,2) = ucl(i,2)
       uclsave(i,3) = ucl(i,3)
      End Do

      Ipass = 1

  98  Continue

c---
      if(Iflow.eq.10) then  ! shear flow over a wall

      if(Iadhere.eq.1) then

       write (6,*) " susp_3d_slv: Ipass: ",Ipass

       xtip = cx + Dcos(theta-pih)
       ytip = cy + Dsin(theta-pih)
       ztip = 0.0D0

       ctx1 = 0.50D0*(cx+xtip)  ! mid-way from center to tip
       cty1 = 0.50D0*(cy+ytip)  ! mid-way
       ctz1 = 0.50D0*(cz+ztip)  ! mid-way

       pftx = -torque/(0.5D0*(cy-ytip)) ! force doublet producing torque

c      write (6,*) theta/pi
 
       cf = -1.0D0/(pi8*visc)

       Iopt_sgf = 1

       if(Ipass.eq.1) then    ! no point force

         pfx = 0.0D0
         pfy = 0.0D0
         pfz = 0.0D0

         Do i=1,ngcl      ! for the iterations
          dpl(i,1) = dpl1(i,1)
          dpl(i,2) = dpl1(i,2)
          dpl(i,3) = dpl1(i,3)
         End Do

       else if(Ipass.eq.2) then  ! point force in the x direction

         pfx = 1.0D0
         pfy = 0.0D0
         pfz = 0.0D0
         Do i=1,ngcl    ! for the iterations
          dpl(i,1) = dpl2(i,1)
          dpl(i,2) = dpl2(i,2)
          dpl(i,3) = dpl2(i,3)
         End Do

       else if(Ipass.eq.3) then    ! point force in the y direction

         pfx = 0.0D0
         pfy = 1.0D0
         pfz = 0.0D0

         Do i=1,ngcl    ! for the iterations
          dpl(i,1) = dpl3(i,1)
          dpl(i,2) = dpl3(i,2)
          dpl(i,3) = dpl3(i,3)
         End Do

       else if(Ipass.eq.4) then     ! correct the point force so that Utip=0

         Det =  A11*A22-A12*A21
         pfx = -(B1*A22-B2*A12)/Det
         pfy = -(B2*A11-B1*A21)/Det
         pfz = 0.0D0

         Do i=1,ngcl    ! for the iterations
          dpl(i,1) = dpl4(i,1)
          dpl(i,2) = dpl4(i,2)
          dpl(i,3) = dpl4(i,3)
         End Do

       end if

c      fil = -adhesion*Dsin(theta)
c      pfx = -fil*Dsin(theta)
c      pfy =  fil*Dcos(theta)
c      pfz = 0.0D0

       pfx1 = -2.0D0*pfx     ! cancel the torque due to the point force
       pfy1 = -2.0D0*pfy
       pfz1 = -2.0D0*pfz

       pfx  = pfx -pftx    ! force couplet
       pfx1 = pfx1+pftx    ! force couplet

       Do i=1,ngcl

        call sgf_3d_w
     +
     +  (Iopt_sgf
     +  ,pcl(i,2),pcl(i,3),pcl(i,1)
     +  ,cy,cz,cx
     +  ,wall
     +  ,gyy,gyz,gyx
     +  ,gzy,gzz,gzx
     +  ,gxy,gxz,gxx
     +  ,py,pz,px
     +  ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +  ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +  ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +  )

      ucl(i,1) = uclsave(i,1)+cf*(Gxx*pfx+Gxy*pfy+Gxz*pfz)
      ucl(i,2) = uclsave(i,2)+cf*(Gyx*pfx+Gyy*pfy+Gyz*pfz)
      ucl(i,3) = uclsave(i,3)+cf*(Gzx*pfx+Gzy*pfy+Gzz*pfz)

c---
        call sgf_3d_w
     +
     +  (Iopt_sgf
     +  ,pcl(i,2),pcl(i,3),pcl(i,1)
     +  ,cty1,ctz1,ctx1
     +  ,wall
     +  ,gyy,gyz,gyx
     +  ,gzy,gzz,gzx
     +  ,gxy,gxz,gxx
     +  ,py,pz,px
     +  ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +  ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +  ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +  )

      ucl(i,1) = ucl(i,1)+cf*(Gxx*pfx1+Gxy*pfy1+Gxz*pfz1)
      ucl(i,2) = ucl(i,2)+cf*(Gyx*pfx1+Gyy*pfy1+Gyz*pfz1)
      ucl(i,3) = ucl(i,3)+cf*(Gzx*pfx1+Gzy*pfy1+Gzz*pfz1)
c     write (6,*) ucl(i,1),ucl(i,2),ucl(i,3)

c---

      End Do

      End If
      End If

c--------------------------------------------
c iterative solution of the integral equation
c for the dipole strength
c--------------------------------------------

      write (6,*) " susp_3d_slv: inner iterations started"

      Iter = 0

 71   Continue

      Iter = Iter+1

c-----------------------------------
c compute the double-layer potential
c-----------------------------------

c     write (6,*) npts,nelm,intm,mint,cx,cy,cz

      call susp_3d_dlp
     +
     +   (npts
     +   ,nelm
     +   ,intm
     +   ,mint
     +   ,mpoly
     +   ,necl,ngcl
     +   ,pcl
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
     +   ,stresslet
     +
     +   ,dlp_pv
     +   )

c     Do i=1,ngcl
c       write (6,*) dlp_pv(i,1),dlp_pv(i,2),dlp_pv(i,3)
c       write (6,*) dpl(i,1),dpl(i,2),dpl(i,3)
c     End Do
c     write (6,*) dpl_sintx,dpl_sinty,dpl_sintz

c----------------------------------------
c translational and rotational velocities
c----------------------------------------

      fc = - pi4/area

      Vx = fc * dpl_sintx
      Vy = fc * dpl_sinty
      Vz = fc * dpl_sintz

c----
c     fc = -(pi4/area)**2
c     fc = -(4.0D0/3.0D0) * (pi4/area)**2
      fc = -1.5D0 * (pi4/area)**2
      Ox = fc * dpl_sinmx
      Oy = fc * dpl_sinmy
      Oz = fc * dpl_sinmz
c----

c----
      rhsx = -pi4*dpl_sinmx
      rhsy = -pi4*dpl_sinmy
      rhsz = -pi4*dpl_sinmz

      call cramer_33
     +
     +   (smomin_xx,smomin_xy,smomin_xz
     +   ,smomin_yx,smomin_yy,smomin_yz
     +   ,smomin_zx,smomin_zy,smomin_zz
     +   ,rhsx,rhsy,rhsz
     +   ,Ox,Oy,Oz)

c----

c     write (6,*) area
c     write (6,*) Vx,Vy,Vz,qdn,Ox,Oy,Oz

c--------------------------------------------
c update the dipole and compute the rms error
c--------------------------------------------

      error = 0.0D0
     
      Do i=1,ngcl

       savex = dpl(i,1)
       savey = dpl(i,2)
       savez = dpl(i,3)

       xx = pcl(i,1)-cx
       yy = pcl(i,2)-cy
       zz = pcl(i,3)-cz

c      write (6,*) vnxcl(i)**2+vnycl(i)**2+vnzcl(i)**2

       dpl(i,1) = -dlp_pv(i,1)/pi4
     +                + vnxcl(i)*dpl_dn/area
     +                + Vx/pi4
     +                + (Oy*zz-Oz*yy)/pi4
     +                - ucl(i,1)/pi4
       dpl(i,2) = -dlp_pv(i,2)/pi4
     +                + vnycl(i)*dpl_dn/area
     +                + Vy/pi4
     +                + (Oz*xx-Ox*zz)/pi4
     +                - ucl(i,2)/pi4
       dpl(i,3) = -dlp_pv(i,3)/pi4
     +                + vnzcl(i)*dpl_dn/area
     +                + Vz/pi4
     +                + (Ox*yy-Oy*xx)/pi4
     +                - ucl(i,3)/pi4

       error = error+(savex-dpl(i,1))**2
     +              +(savey-dpl(i,2))**2
     +              +(savez-dpl(i,3))**2

c      write (6,100) i,dlp_pv(i,1),dlp_pv(i,2),dlp_pv(i,3)
c      write (6,100) i,dpl(i,1),dpl(i,2),dpl(i,3)

      End Do

c----------------
c end of updating
c----------------

      error = Dsqrt(error)/(3.0D0*ngcl)

      write (6,457) Iter,error

      If(error.lt.tol) Go to 99   ! solution found
      If(Iter.lt.Nter) Go to 71   ! one more iteration

      ! one more iteration batch:

      Iter = 0

      write (6,*)
      write (6,*) " Not yet converged"
      write (6,*) " One more iteration batch ?"
      write (6,*)
      write (6,*) " Enter 1 for Yes, 0 for No"
      write (6,*) "--------------------------"
      read  (5,*) more

      if(more.eq.1) Go to 71

c-------------------------
c done with the iterations
c-------------------------

  99  Continue

      write (6,*) " susp_3d_slv: iterations :",Iter

c-----
c particle stress tensor
c-----

      fc = -pi8/vlm

      stresslet(1,2) = fc*(stresslet(1,2)+stresslet(2,1) )
      stresslet(1,3) = fc*(stresslet(1,3)+stresslet(3,1) )
      stresslet(2,3) = fc*(stresslet(2,3)+stresslet(3,2) )
      stresslet(2,1) = stresslet(1,2)
      stresslet(3,1) = stresslet(1,3)
      stresslet(3,2) = stresslet(2,3)

c----------------
c adhesion forces
c----------------

      If(Iflow.eq.10) then

      If(Iadhere.eq.1) then

c-----
c particle tip velocity
c-----

      Uxtip = Vx - Oz*sin(theta-pih)
      Uytip = Vy + Oz*cos(theta-pih)

      if(Ipass.eq.1) then
        B1 = Uxtip
        B2 = Uytip
         Do i=1,ngcl
          dpl1(i,1) = dpl(i,1)
          dpl1(i,2) = dpl(i,2)
          dpl1(i,3) = dpl(i,3)
         End Do
        Ipass = 2
        Go to 98
      else if(Ipass.eq.2) then
        A11 = Uxtip - B1
        A21 = Uytip - B2
         Do i=1,ngcl
          dpl2(i,1) = dpl(i,1)
          dpl2(i,2) = dpl(i,2)
          dpl2(i,3) = dpl(i,3)
         End Do
        Ipass = 3
        Go to 98
      else if(Ipass.eq.3) then
        A12 = Uxtip - B1
        A22 = Uytip - B2
         Do i=1,ngcl
          dpl3(i,1) = dpl(i,1)
          dpl3(i,2) = dpl(i,2)
          dpl3(i,3) = dpl(i,3)
         End Do
        Ipass = 4
        Go to 98
      else if(Ipass.eq.4) then
         Do i=1,ngcl
          dpl4(i,1) = dpl(i,1)
          dpl4(i,2) = dpl(i,2)
          dpl4(i,3) = dpl(i,3)
         End Do
      end if

c     write (6,102) Vx,Vy,Oz
c     write (6,102) Uxtip,Uytip
      write (6,*) " susp_3d_slv: particle tip vel: ",Uxtip,Uytip
      write (6,*) " susp_3d_slv: particle tip: ",xtip,ytip
      xforce_adh = pfx
      yforce_adh = pfy
      zforce_adh = pfz
      write (6,102) xforce_adh,yforce_adh,zforce_adh
c     pause

      end if
      end if

c-----
c Done
c-----

  100 Format (1x,i5,10(1x,f10.5))
  101 Format (f15.10)
  102 Format (10(1x,f15.8))

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)

  456 Format (" Iter: ",I3," Proj:",f12.8," Error: ",f12.10,/
     +       ," Vx: ",f12.8," Vy: ",f12.8," Vz: ",f12.8,/
     +       ," Ox: ",f12.8," Oy: ",f12.8," Oz: ",f12.8)
  457 Format (" susp_3d_slv: Inner iter: ",I3," Error: ",f12.10)
  459 Format (" Flow across the boundary:",f15.10)

  800 Format(" centroid:",3(1x,f10.5))

      Return
      End
