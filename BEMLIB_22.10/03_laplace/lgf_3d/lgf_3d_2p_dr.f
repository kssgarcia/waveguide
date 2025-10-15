      program lgf_3d_2p_dr

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------
c Driver of a subroutine for computing the 
c 3D doubly periodic Green's function
c of Laplace's equation
c
c The point sources are arranged in a plane
c that is perpendicular to the z axis
c (parallel to the xy plane)
c
c SYMBOLS:
c --------
c
c method = 1 for the Fourier series method
c method = 2 for the fast summation method of Hautman & Klein
c--------------------------------------

      Implicit double precision (a-h,o-z)

      Dimension xx(50),yy(50),zz(50),ww(50)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi

      Null = 0

c--------
c prepare
c--------

      open (1,file="lgf_3d_2p.out",status="unknown")

 96   Continue

      write (6,*)
      write (6,*) "  Enter:"
      write (6,*)
      write (6,*) "  1 for one evaluation"
      write (6,*) "  2 for a profile along the z axis"
      write (6,*) "  3 to test the integral identity on a sphere"
      write (6,*) "  4 to test the integral identity "
      write (6,*) "    over a periodic patch"
      write (6,*) "  0 to quit"
      write (6,*) "  ---------"
      write (6,*)
      read  (5,*) menu

      if(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) "  Please enter the components"
      write (6,*) "  of the first base vector: a11, a12"
      write (6,*) "  ----------------------------------"
      read  (5,*) a11,a12

      write (6,*)
      write (6,*) "  Please enter the components"
      write (6,*) "  of the second base vector: a21, a22"
      write (6,*) "  -----------------------------------"
      read  (5,*) a21,a22

c--------------------------
c ewald_3d_2p will produce:
c
c     the reciprocal vectors
c     the optimal value of xi (called ew)
c     and the area of the unit cell
c--------------------------

      call lgf_3d_2p_ewald
     +
     +   (a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,area
     +   )

      ewsave = ew    ! save the optimal value

c-----------
c displaying
c-----------

      write (6,*)
      write (6,*) " coordinates of the reciprocal lattice vectors:"
      write (6,*) " ----------------------------------------------"
      write (6,100) b11,b12
      write (6,100) b21,b22

      write (6,*)
      write (6,*) " enter the coordinates (x0,y0,z0)"
      write (6,*) " of a singularity"
      write (6,*) " ----------------"
      read  (5,*) x0,y0,z0

c---

 98   Continue

      write (6,*) "  Choose the method"
      write (6,*)
      write (6,*) "  Enter:"
      write (6,*)
      write (6,*) "  1 for the Fourier-series method"
      write (6,*) "  2 for the fast-summation method"
      write (6,*) "  0 to quit"
      write (6,*) " ----------"
      read  (5,*) method

      if(method.eq.0) Go to 99

c-------------------------
      if(method.eq.1) then
c-------------------------

       write (6,*)
       write (6,*) " Enter the truncation limit for summation"
       write (6,*) " over the reciprocal lattice"
       write (6,*) " 0 to quit"
       write (6,*) " ---------"
       read  (5,*) Max2

       if(Max2.eq.0) Go to 97

c------------------------------
      else if(method.eq.2) then
c------------------------------

       write (6,*)
       write (6,*) " Enter the truncation limit for summation"
       write (6,*) " over the physical lattice in the xy plane"
       write (6,*) " 0 to quit"
       write (6,*) " ---------"
       read  (5,*) Max1

       if(Max1.eq.0) Go to 97

       write (6,*)
       write (6,*) " Enter the truncation limit for summation"
       write (6,*) " over the reciprocal lattice"
       write (6,*) " 0 to quit"
       write (6,*) " ---------"
       read  (5,*) Max2

       if(Max2.eq.0) Go to 97

       write (6,*)
       write (6,*) " Enter the truncation limit for summation"
       write (6,*) " over the 3D physical lattice "
       write (6,*) " 0 to quit"
       write (6,*) " ---------"
       read  (5,*) Max3

       if(Max3.eq.0) Go to 97

       write (6,*)
       write (6,*) " Enter the value of ewald parameter"
       write (6,*) " 99 for the optimal value ",ewsave
       write (6,*) " 0 to quit"
       write (6,*) " ---------"
       read  (5,*) ew

       if(ew.eq.0) Go to 97

       if(ew.eq.99) ew = ewsave

c-----------
      end if
c-----------

c---------------
c One evaluation
c---------------

      if(menu.eq.1) then

      write (6,*)
      write (6,*) " Please enter field point coordinates: x,y,z"
      write (6,*) " -------------------------------------------"
      read  (5,*) x,y,z

      Iselect = 2

      call lgf_3d_2p
     +
     +  (Iselect
     +  ,method
     +  ,x,y,z
     +  ,x0,y0,z0
     +  ,a11,a12,a21,a22
     +  ,b11,b12,b21,b22
     +  ,ew,area
     +  ,Max1,Max2,Max3
     +  ,G
     +  ,Gx,Gy,Gz
     +  )

      write (6,*) ' ----------------------------'
      write (6,*) ' Green function and gradient '
      write (6,100)
      write (6,100) G
      write (6,100) Gx,Gy,Gz
      write (6,*) ' ----------------------------'

      End If

c--------------------------------
c Draw a profile along the z axis
c--------------------------------

      if(menu.eq.2) then

      write (6,*)
      write (6,*) " Please enter the x,y coordinates of the line"
      write (6,*) " --------------------------------------------"
      read  (5,*) x,y

      write (6,*)
      write (6,*) " Please z_min and z_max"
      write (6,*) " ----------------------"
      read  (5,*) zmin,zmax

      write (6,*)
      write (6,*) " Please the number of profile points"
      write (6,*) " -----------------------------------"
      read  (5,*) nprof

      Iselect = 2

      nprof1 = nprof+1
      zstp = (zmax-zmin)/nprof

      write (1,*) nprof1

      Do i=1,m1

        z = zmin+(i-1.0)*zstp

        call lgf_3d_2p
     +
     +   (Iselect
     +   ,method
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,area
     +   ,Max1,Max2,Max3
     +   ,G
     +   ,Gx,Gy,Gz
     +   )

        write (6,101) i,z,G,Gx,Gy,Gz
        write (1,101) i,z,G,Gx,Gy,Gz

        End Do

      End If

c--------------------------
c test of integral identity
c on a sphere
c--------------------------

      if(menu.eq.3) then

      write (6,*)
      write (6,*) ' Enter the center of the test sphere'
      write (6,*) ' 99 to quit'
      write (6,*) ' ----------'
      read  (5,*) xcnt,ycnt,zcnt

      if(xcnt.eq.99) Go to 99
      if(ycnt.eq.99) Go to 99
      if(zcnt.eq.99) Go to 99

      write (6,*)
      write (6,*) ' Enter the radius of the test sphere'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) rad

      if(rad.eq.0) Go to 99

      write (6,*)
      write (6,*) ' Choose the number of quadrature points '
      write (6,*)
      write (6,*) ' Select from 6, 18'
      write (6,*) '                      0 to quit'
      write (6,*) ' ------------------------------'
      read  (5,*) mint

      if(mint.eq.0) Go to 99

      call gauss_sph (mint,xx,yy,zz,ww)

c---
c integrate over the sphere
c---

      sm1 = 0.0D0

      Do i=1,mint

       vnx = xx(i)
       vny = yy(i)
       vnz = zz(i)

       x = xcnt+rad*vnx
       y = ycnt+rad*vny
       z = zcnt+rad*vnz

        call lgf_3d_2p
     +
     +    (Iselect
     +    ,method
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,a11,a12,a21,a22
     +    ,b11,b12,b21,b22
     +    ,ew,area
     +    ,Max1,Max2,Max3
     +    ,G
     +    ,Gx,Gy,Gz
     +    )

       sm1 = sm1 + (vnx*Gx + vny*Gy + vnz*Gz)*ww(i)

       End Do

      cf = pi4*rad**2

      sm1 = sm1 * cf

      write (6,*)
      write (6,*)  "Should be -1.0 or 0"
      write (6,*)
      write (6,100) sm1

c--------------------------
c Test of integral identity
c on a periodic patch
c
c Integrals computed by the trapezoidal rule
c--------------------------

      else if(menu.eq.4) then

       write (6,*)
       write (6,*) ' Enter the z elevation of the patch'
       write (6,*) ' 99 to quit'
       write (6,*) ' ----------'
       read  (5,*) zp

       if(zp.eq.99) Go to 99

       write (6,*)
       write (6,*) ' Enter the discretization level of the patch'
       write (6,*) ' for the x and y directions'
       write (6,*) ' --------------------------'
       read  (5,*) Ndivx,Ndivy

c---
c will integrate over a quadrilateral cell
c using the trapezoidal rule
c---

       Ndivx1 = Ndivx+1
       Ndivy1 = Ndivy+1

       xstep1 = a11/Ndivx
       xstep2 = a12/Ndivx
       ystep1 = a21/Ndivy
       ystep2 = a22/Ndivy

       vnx = 0.0D0   ! normal vector pointing upward
       vny = 0.0D0
       vnz = 1.0D0

       z = zp

       sumx = 0.0D0

       Do i=1,Ndivx1

        sumy = 0.0D0

        Do j=1,Ndivy1

        x = (i-1.0)*xstep1+(j-1.0)*xstep2
        y = (i-1.0)*ystep1+(j-1.0)*ystep2

        call lgf_3d_2p
     +
     +    (Iselect
     +    ,method
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,a11,a12,a21,a22
     +    ,b11,b12,b21,b22
     +    ,ew,area
     +    ,Max1,Max2,Max3
     +    ,G
     +    ,Gx,Gy,Gz
     +    )

        fc = 1.0D0                             ! trapezoidal rule weight
        If(j.eq.1.or.j.eq.Ndivy1) fc = 0.5D0   ! trapezoidal rule weight

        sumy = sumy + (vnx*Gx + vny*Gy + vnz*Gz)*fc

       End Do

       fc = 1.0D0                             ! trapezoidal rule weight
       If(i.eq.1.or.i.eq.Ndivx1) fc = 0.5D0   ! trapezoidal rule weight

       sumx = sumx + fc*sumy

      End Do

      sumx = sumx*area/(Ndivx*Ndivy)

      write (6,*)
      write (6,*)  " should be -0.5, 0, or 0.5"
      write (6,*)
      write (6,100) sumx

c-----------
      End If
c-----------

      Go to 98   ! return to the menu

 97   Continue

      Go to 96   ! return to the main menu

c-----
c Done
c-----

 99   Continue

      write (1,*) Null
      close (1) 

  100 Format (3(1x,f20.10))
  101 Format (1x,i3,10(1x,f9.5))

      Stop
      End
