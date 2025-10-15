      subroutine susp_3d_col
     +
     +   (Ioctaicos
     +   ,ndiv
     +   ,boa,coa
     +   ,req
     +   ,Isize
     +   ,cx,cy,cz
     +   ,theta,phi
     +   ,mpoly
     +   ,npts,nelm
     +   ,necl
     +   ,ngcl
     +   ,pcl
     +   )

c=======================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c----------------------------
c
c SYMBOLS:
c -------
c
c Generate collocation points
c over a particle
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)

      Dimension      n(512,6), nbe(512,3)
      Dimension   alpha(512),  beta(512),gamma(512)
      Dimension    arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension  crvmel(512)

c collocation:

      Dimension vmaster(8)
      Dimension  xcl(512,100),ycl(512,100),zcl(512,100)

      Dimension   ncl(512,100)   ! connectivity
      Dimension    pcl(2306,3)   ! nodes

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna

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

c     write (6,*) " prtcl_3ds_slv: number of points:   ",npts
c     write (6,*) " prtcl_3ds_slv: number of elements: ",nelm

c----------------------------
c Expand to specified radius
c Deform to an ellipsoid
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

c     write (6,*) " prtcl_3ds_slv: x semi-axis = ",x_axis
c     write (6,*) " prtcl_3ds_slv: y semi-axis = ",y_axis
c     write (6,*) " prtcl_3ds_slv: z semi-axis = ",z_axis

c------------------------
c rotate by theta and phi
c------------------------

      cs = Dcos(theta)  ! rotate around the z axis
      sn = Dsin(theta)  ! by theta

      Do i=1,npts                
       tmpx = cs*p(i,1)-sn*p(i,2) 
       tmpy = sn*p(i,1)+cs*p(i,2)
       tmpz = p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = Dcos(phi)   ! rotate around the x axis
      sn = Dsin(phi)   ! by phi

      Do i=1,npts        
       tmpx = p(i,1)
       tmpy = cs*p(i,2)-sn*p(i,3)
       tmpz = sn*p(i,2)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

c--------------------
c translate center to
c specified position
c--------------------

      Do i=1,npts 
        p(i,1) = p(i,1) + cx
        p(i,2) = p(i,2) + cy
        p(i,3) = p(i,3) + cz
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
c compute the collocation nodes
c over each triangle
c------------------------------

      Iopt = 1   ! interpolation option (need only the position)

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
     +     ,Iopt
     +     )

         Ic = Ic+1

         xcl(k,Ic) = x
         ycl(k,Ic) = y
         zcl(k,Ic) = z

         End Do
       End Do

      End Do

      necl = Ic
      neclt = Jc

c     write (6,*) " prtcl_3ds_slv: No of elm col nodes:",necl
c     write (6,*) " prtcl_3ds_slv: Total of elm cl nodes:", neclt

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

      Do i=1,necl
       pcl(i,1) = xcl(1,i)
       pcl(i,2) = ycl(1,i)
       pcl(i,3) = zcl(1,i)
       ncl(1,i) = i
      End Do

      ngcl = necl

      eps = 0.0000000001

      Do i=2,nelm        ! loop over elements
        Do j=1,necl         ! loop over element nodes

        Iflag=0

         Do k=1,ngcl

          if(abs(xcl(i,j)-pcl(k,1)).le.eps) then
           if(abs(ycl(i,j)-pcl(k,2)).le.eps) then
            if(abs(zcl(i,j)-pcl(k,3)).le.eps) then

             Iflag = 1      ! the node has been recorded previously
             ncl(i,j) = k   ! the jth local node of element i
                            ! is the kth global node

            end if
           end if
          end if

         End Do

         if(Iflag.eq.0) then     ! accept the node

          ngcl = ngcl+1          ! one more global node

          pcl(ngcl,1) = xcl(i,j)
          pcl(ngcl,2) = ycl(i,j)
          pcl(ngcl,3) = zcl(i,j)

          ncl(i,j) = ngcl   ! the jth local node of element i
                            ! is the new global node
         end if

       End Do
      End Do                      !  end of loop over elements


c-----
c done
C-----

c     write (6,*) "prtcl_3ds_slv: No of global col nodes:",ngcl

c     Do i=1,nelm
c     write (6,*) (ncl(i,j),j=1,necl)
c     End Do

c     Do i=1,ngcl
c     write (6,*) pcl(i,1),pcl(i,2),pcl(i,3)
c     End Do


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

      Return
      End
