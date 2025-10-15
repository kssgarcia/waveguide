      program sgf_3d_w_dr

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------
c Driver for the Green's function in a
c semi-infinite domain bounded by
c a plane wall located at x=wall
c-------------------------------------
 
      Implicit Double Precision (a-h,o-z)

      Dimension xx(50),yy(50),zz(50),ww(50)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0 *pi
      pi8 = 8.0D0 *pi

c----------
c constants
c----------

      Iopt = 2      ! will need velocity, pressure, and stress

c-----
c menu
c-----

  98  Continue

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for one evaluation'
      write (6,*) ' 2 for a test of integral identities'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) menu

      If(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) ' The wall is located at x = wall'
      write (6,*)
      write (6,*) ' Please enter: wall'
      write (6,*) ' ------------------'
      read  (5,*) wall

      write (6,*) 
      write (6,*) ' Enter (x0,y0,z0) of the point force'
      write (6,*) 
      write (6,*) ' 99 for anyone to quit'
      write (6,*) ' ---------------------'
      write (6,*)
      read  (5,*) x0,y0,z0

      If(x0.eq.99) Go to 99
      If(y0.eq.99) Go to 99
      If(z0.eq.99) Go to 99

c---------------
c one evaluation
c---------------

      if(menu.eq.1) then

      write (6,*)
      write (6,*) ' Enter the (x,y,z) of the velocity point'
      write (6,*) ' ---------------------------------------'
      read  (5,*) x,y,z

      call sgf_3d_w
     +
     +  (Iopt
     +  ,x,y,z
     +  ,x0,y0,z0
     +  ,wall
     +  ,gxx,gxy,gxz
     +  ,gyx,gyy,gyz
     +  ,gzx,gzy,gzz
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

      write (6,*) ' ---------------------------------'
      write (6,*) ' Green function '
      write (6,*) ' pressure and stress tensor:'
      write (6,*) 

      write (6,100)
      write (6,100) gxx,gxy,gxz
      write (6,100) gyx,gyy,gyz
      write (6,100) gzx,gzy,gzz
      write (6,100)

      write (6,100)
      write (6,100) px,py,pz
      write (6,100)

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

c----------------------------
c test of integral properties
c----------------------------

      if(menu.eq.2) then

      write (6,*) ' ---------'
      write (6,*) ' Enter the center of the test sphere'
      write (6,*) ' 99 to quit'
      write (6,*) ' ----------'
      read  (5,*) xcnt,ycnt,zcnt

      If(xcnt.eq.99) Go to 99
      If(ycnt.eq.99) Go to 99
      If(zcnt.eq.99) Go to 99

      write (6,*) ' ---------'
      write (6,*) ' Enter the radius of the test sphere'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) rad

      If(rad.eq.0) Go to 99

      write (6,*) ' ---------'
      write (6,*) ' Choose the integration formula'
      write (6,*) ' Select from 6, 18'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) mint

      if(mint.eq.0) Go to 99

      call gauss_sph (mint,xx,yy,zz,ww)

      sm1 = 0.D0
      sm2 = 0.D0
      sm3 = 0.D0

      sxx = 0.D0
      sxy = 0.D0
      sxz = 0.D0

      syx = 0.D0
      syy = 0.D0
      syz = 0.D0

      szx = 0.D0
      szy = 0.D0
      szz = 0.D0

      Do i=1,mint

      x = xcnt+rad*xx(i)
      y = ycnt+rad*yy(i)
      z = zcnt+rad*zz(i)

      call sgf_3d_w 
     +
     +    (Iopt
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,wall
     +    ,gxx,gxy,gxz
     +    ,gyx,gyy,gyz
     +    ,gzx,gzy,gzz
     +    ,px,py,pz
     +    ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +    ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +    ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +    )

      Tyxx = Txxy
      Tzxx = Txxz
      Tzxy = Tyxz

      Tyyx = Txyy
      Tzyx = Txyz
      Tzyy = Tyyz

      Tyzx = Txzy
      Tzzx = Txzz
      Tzzy = Tyzz

      vnx = xx(i)   ! coordinates on the unit sphere
      vny = yy(i)
      vnz = zz(i)

      sm1 = sm1 + (vnx*gxx + vny*gyx + vnz*gzx)*ww(i)
      sm2 = sm2 + (vnx*gxy + vny*gyy + vnz*gzy)*ww(i)
      sm3 = sm3 + (vnx*gxz + vny*gyz + vnz*gzz)*ww(i)

c     sm1 = sm1 + (gxx*vnx+gxy*vny+gxz*vnz)*ww(i)
c     sm2 = sm2 + (gyx*vnx+gyy*vny+gyz*vny)*ww(i)
c     sm3 = sm3 + (gzx*vnx+gzy*vny+gzz*vny)*ww(i)

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

      cf = pi4*rad*rad

c: continuity equation
 
      sm1 = sm1 * cf
      sm2 = sm2 * cf
      sm3 = sm3 * cf

      write (6,*) 
      write (6,*) " Should be 0:"
      write (6,*) 

      write (6,100) sm1,sm2,sm3

c: point force
 
      cf = pi4*rad*rad/pi8

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
      write (6,*) " Should be 0 or negative of I:"
      write (6,*) 

      write (6,*) 
      write (6,100) sxx,sxy,sxz
      write (6,100) syx,syy,syz
      write (6,100) szx,szy,szz
      write (6,*) 

c---
      end if
c---

      Go to 98    ! return to repeat

c------
c  done
c------

  99  Continue

 100  Format (3(2x,f15.10))

      stop
      end
