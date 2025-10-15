      program sgf_3d_3p_dr

c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c Compute the triply periodic Green's function of
c three-dimensional Stokes flow
c with vanishing flow rate across the faces of the
c periodic cell
c
c Pozrikidis (1996, J. Eng. Math. vol 30, 79-96)
c
c Important note: 
c
c subroutines  sgf_3d_3p_ewald
c              sgf_3d_3p_qqq
c              sgf_3d_3p_vvv
c 
c must be called before sgf_3d_3d
c
c SYMBOLS:
c -------
c
c xx, yy, zz:	Gauss integration points over a sphere
c
c----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xx(50),yy(50),zz(50),ww(50)

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      pi4 = 4.0D0*pi
      pi8 = 8.0D0*pi
      
c---------------
c prepare to run
c---------------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for one value "
      write (6,*) " 2 to test integral identities on a sphere"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c--------
c prepare
c--------

      Iopt = 2   ! will compute G, p, and T

c------
c input
c------

      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to read the base vectors from file: sgf_3d_3p.dat"
      write (6,*) " 2 to type in input"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ienrd

c------------------------
      if(Ienrd.eq.2) then
c------------------------

       write (6,*) 
       write (6,*) " Enter the coordinates of the "
       write (6,*) " first lattice vector"
       write (6,*) " --------------------"
       read  (5,*) a11,a12,a13

       write (6,*) " Enter the coordinates of the "
       write (6,*) " second lattice vector"
       write (6,*) " ---------------------"
       read  (5,*) a21,a22,a23

       write (6,*) " Enter the coordinates of the "
       write (6,*) " third lattice vector"
       write (6,*) " --------------------"
       read  (5,*) a31,a32,a33

c---------
      else
c---------

        open (3,file="sgf_3d_3p.dat")

        read (3,*) a11,a12,a13
        read (3,*) a21,a22,a23
        read (3,*) a31,a32,a33

        close (3)

c-----------
      end if
c-----------

c----
c  sgf_3d_3p_ewald will produce the reciprocal vectors
c  and the optimal value of xi
c  according to Beenakker
c----

      call sgf_3d_3p_ewald
     +
     +  (a11,a12,a13
     +  ,a21,a22,a23
     +  ,a31,a32,a33
     +  ,b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,ewald,tau
     +  )

      write (6,*)
      write (6,*) " Lattice base vectors:"
      write (6,*) " --------------------"
      write (6,100) a11,a12,a13
      write (6,100) a21,a22,a23
      write (6,100) a31,a32,a33

      write (6,*)
      write (6,*) " Reciprocal lattice vectors:"
      write (6,*) " ---------------------------"
      write (6,100) b11,b12,b13
      write (6,100) b21,b22,b23
      write (6,100) b31,b32,b33

c---

 98   Continue

      write (6,*) 
      write (6,*) " Enter the value of xi"
      write (6,*) " 99 for the Beenakker value: ",ewald
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) ew
c     ew = 99

      if(ew.eq.0 ) Go to 99
      if(ew.eq.99) ew = ewald

      write (6,*) 
      write (6,*) " Enter max1, max2 for summation"
      write (6,*) " in real and reciprocal space. "
      write (6,*) " must be less than 9"
      write (6,*) " -------------------"
c     read  (5,*) max1,max2
      max1 = 4
      max2 = 4

      write (6,*) 
      write (6,*) " Enter (x0,y0,z0) of a point force"
      write (6,*) " ---------------------------------"
c     read  (5,*) x0,y0,z0
      x0 = 0.0D0
      y0 = 0.0D0
      z0 = 0.0D0

c--------------------------------------------
c  qqq_3d will produce an array used to compute
c  the sum of the velocity Green's function
c  in reciprocal space
c--------------------------------------------

      call sgf_3d_3p_qqq
     +
     +  (b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,max2,ew
     +  )

c-----------------------------------------
c  vvv_3d will produce an array used to compute
c  the sum of the pressure and stress Green's function
c  in reciprocal space
c-----------------------------------------

      if(Iopt.eq.2) then

      call sgf_3d_3p_vvv
     +
     +  (b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,max2,ew
     +  )

      end if

c---------------
c one evaluation
c---------------

      if(menu.eq.1) then

      write (6,*)
      write (6,*) ' Enter the (x,y,z) of the velocity point'
      write (6,*) ' ---------------------------------------'
      read  (5,*) x,y,z

c     x =0.50D0
c     y =0.50D0
c     z =0.50D0

      call sgf_3d_3p
     +
     +  (Iopt
     +  ,x,y,z
     +  ,x0,y0,z0
     +  ,a11,a12,a13
     +  ,a21,a22,a23
     +  ,a31,a32,a33
     +  ,b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,ew,tau
     +  ,max1,max2
     +  ,Gxx,Gxy,Gxz
     +  ,Gyx,Gyy,Gyz
     +  ,Gzx,Gzy,Gzz
     +  ,px,py,pz
     +  ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +  ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +  ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +  )

      Tyxx = Txxy
      Tzxx = Txxz
      Tzxy = Tyxz

      Tyyx = Txyy
      Tzyx = Txyz
      Tzyy = Tyyz

      Tyzx = Txzy
      Tzzx = Txzz
      Tzzy = Tyzz

      write (6,*) ' -------------------------------'
      write (6,*) ' Green function for the velocity'
      write (6,*) ' -------------------------------'

      write (6,100) Gxx,Gxy,Gxz
      write (6,100) Gyx,Gyy,Gyz
      write (6,100) Gzx,Gzy,Gzz

      write (6,*) ' -------------------------------'
      write (6,*) ' Green function for the pressure'
      write (6,*) ' -------------------------------'

      write (6,100) px,py,pz

      write (6,*) ' -----------------------------'
      write (6,*) ' Green function for the stress'
      write (6,*) ' -----------------------------'

      write (6,100) Txxx,Txxy,Txxz
      write (6,100) Tyxx,Tyxy,Tyxz
      write (6,100) Tzxx,Tzxy,Tzxz
      write (6,100)
      write (6,100) Txyx,Txyy,Txyz
      write (6,100) Tyyx,Tyyy,Tyyz
      write (6,100) Tzyx,Tzyy,Tzyz
      write (6,100)
      write (6,100) Txzx,Txzy,Txzz
      write (6,100) Tyzx,Tyzy,Tyzz
      write (6,100) Tzzx,Tzzy,Tzzz

      end if

c-------------------------
c test integral identities 
c-------------------------

      if(menu.eq.2) then

      write (6,*)
      write (6,*) ' Enter the coordinates of the center'
      write (6,*) ' of the test sphere'
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
      write (6,*) ' Will integrate using the sphere quadrature.'
      write (6,*) ' Choose the number of integration points.'
      write (6,*) ' Select from: 6, 18'
      write (6,*) ' Enter 0 to quit'
      write (6,*) ' ---------------'
      read  (5,*) mint

      if(mint.eq.0) Go to 99

c---
c read the gauss-sphere base points and weights
c---

      call gauss_sph 
     +
     +  (mint
     +  ,xx,yy,zz
     +  ,ww
     +  )

c---
c initialize
c---

      sm1 = 0.0D0
      sm2 = 0.0D0
      sm3 = 0.0D0

      sxx = 0.0D0
      sxy = 0.0D0
      sxz = 0.0D0

      syx = 0.0D0
      syy = 0.0D0
      syz = 0.0D0

      szx = 0.0D0
      szy = 0.0D0
      szz = 0.0D0

c---
c perform the quadrature
c---

      Do i=1,mint

       x = xcnt + rad*xx(i)
       y = ycnt + rad*yy(i)
       z = zcnt + rad*zz(i)

        call sgf_3d_3p
     +
     +   (Iopt
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,a11,a12,a13,a21,a22,a23,a31,a32,a33
     +   ,b11,b12,b13,b21,b22,b23,b31,b32,b33
     +   ,ew,tau
     +   ,max1,max2
     +   ,Gxx,Gxy,Gxz
     +   ,Gyx,Gyy,Gyz
     +   ,Gzx,Gzy,Gzz
     +   ,px,py,pz
     +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +   )

       Tyxx = Txxy
       Tzxx = Txxz
       Tzxy = Tyxz

       Tyyx = Txyy
       Tzyx = Txyz
       Tzyy = Tyyz

       Tyzx = Txzy
       Tzzx = Txzz
       Tzzy = Tyzz

       vnx = xx(i)
       vny = yy(i)
       vnz = zz(i)

       sm1 = sm1 + (vnx*Gxx + vny*Gyx + vnz*Gzx)*ww(i)
       sm2 = sm2 + (vnx*Gxy + vny*Gyy + vnz*Gzy)*ww(i)
       sm3 = sm3 + (vnx*Gxz + vny*Gyz + vnz*Gzz)*ww(i)

       sxx = sxx + (vnx*Txxx + vny*Tyxx + vnz*Tzxx)*ww(i)
       sxy = sxy + (vnx*Txxy + vny*Tyxy + vnz*Tzxy)*ww(i)
       sxz = sxz + (vnx*Txxz + vny*Tyxz + vnz*Tzxz)*ww(i)

       syx = syx + (vnx*Txyx + vny*Tyyx + vnz*Tzyx)*ww(i)
       syy = syy + (vnx*Txyy + vny*Tyyy + vnz*Tzyy)*ww(i)
       syz = syz + (vnx*Txyz + vny*Tyyz + vnz*Tzyz)*ww(i)

       szx = szx + (vnx*Txzx + vny*Tyzx + vnz*Tzzx)*ww(i)
       szy = szy + (vnx*Txzy + vny*Tyzy + vnz*Tzzy)*ww(i)
       szz = szz + (vnx*Txzz + vny*Tyzz + vnz*Tzzz)*ww(i)

      End Do

      cf = pi4*rad**2

      sm1 = sm1 * cf
      sm2 = sm2 * cf
      sm3 = sm3 * cf

      write (6,*)
      write (6,*) "Shoud be zero:"
      write (6,*)
      write (6,100) sm1,sm2,sm3

      cf = pi4*rad**2/pi8

      sxx = sxx * cf
      sxy = sxy * cf
      sxz = sxz * cf
      syx = syx * cf
      syy = syy * cf
      syz = syz * cf
      szx = szx * cf
      szy = szy * cf
      szz = szz * cf

      write (6,*)
      write (6,*) " Should be zero or negative of the identity matrix"
      write (6,*)
      write (6,100) sxx,sxy,sxz
      write (6,100) syx,syy,syz
      write (6,100) szx,szy,szz

c-------------
      end if
c-------------

      Go to 98    ! return to the main menu

 99   Continue

c-----
c done
c-----

 100  Format (3(1x,f15.10))

      Stop
      End
