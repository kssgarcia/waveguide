      program integral_quadri

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------------
c integral a function f(x y,z)
c inside a quadrilateral volume
c
c ax<x<bx, ay<y<by, az<z<bz
c
c using the trapezoidal rule
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension weight(129,129,129),f(129,129,129)
      Dimension wx(129),wy(129),wz(129)

      common/lgf/m,m1,m2,m3

c---------
c constants
c---------

      pi = 3.14159 265358 D0

c------------
c parameters
c------------
 
   98 Continue

      write (6,*)
      write (6,*) ' Menu of functions'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) '  1 for f(x,y) = ...'
      write (6,*) '                 ...'
      write (6,*) '  2 for f(x,y) = ...'
      write (6,*) '                 ...'
      write (6,*) '  3 for f(x,y) = ...'
      write (6,*) '                 ...'
      write (6,*) '  4 for f(x,y) = ...'
      write (6,*) '                 ...'
      write (6,*) '  0 to quit'
      write (6,*) ' ----------'
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c----------
      if(menu.eq.1) then
       ax = 0.0D0
       bx = 2.0D0*pi
       ay = 0.0D0
       by = 2.0D0*pi
       az = 0.0D0
       bz = 2.0D0*pi
       m1   = 0
       m2   = 0
       m3   = 1
       Nx = 8*8
       Ny = 8*8
       Nz = 8*8
      end if
c----------

c----------
      if(menu.eq.2) then
       ax = 0.0D0
       bx = 1.0D0*pi
       ay = 0.0D0
       by = 1.0D0*pi
       az = 0.0D0
       bz = 1.0D0*pi
       m1   = 1
       m2   = 0
       m3   = 0
       Nx = 8*8
       Ny = 8*8
       Nz = 8*8
      end if
c----------

c----------
      if(menu.eq.3) then
       ax = 0.0D0
       bx = 2.0D0*pi
       ay = 0.0D0
       by = 2.0D0*pi
       az = 0.0D0
       bz = 2.0D0*pi
       m1   = 1
       m2   = 1
       m3   = 1
       Nx = 2*8*8
       Ny = 2*8*8
       Nz = 2*8*8
      end if
c----------

c----------
      if(menu.eq.4) then
       ax = 0.0D0
       bx = 1.0D0*pi
       ay = 0.0D0
       by = 1.0D0*pi
       az = 0.0D0
       bz = 1.0D0*pi
       m1   = 1
       m2   = 1
       m3   = 1
       Nx = 8*8
       Ny = 8*8
       Nz = 8*8
      end if
c----------

c---
c prepare
c---

      Dx = (bx-ax)/Nx
      Dy = (by-ay)/Ny
      Dz = (bz-az)/Nz

c----
c weight matrix
c----

      wx(1) = 0.5D0;

      Do i=2,Nx
        wx(i) = 1.0D0;
      End Do
      wx(Nx+1) = 0.5D0;

      wy(1) = 0.5D0;

      Do i=2,Ny
        wy(i) = 1.0D0;
      End Do
      wy(Ny+1) = 0.5D0;

      wz(1) = 0.5D0;

      Do i=2,Nz
        wz(i) = 1.0D0;
      End Do
      wz(Nz+1) = 0.5D0;

      Do i=1,Nx+1
       Do j=1,Ny+1
        Do k=1,Nz+1
         weight(i,j,k) = wx(i)*wy(j)*wz(k);
        End Do
       End Do
      End Do

c----
c print and test the weights
c----

c     sum = 0.0D0
c     Do i=1,Nx+1
c      Do j=1,Ny+1
c       Do k=1,Nz+1
c        sum = sum+weight(i,j)
c       End Do
c      End Do
c     End Do

c     write (6,*) sum,Nx*Ny*Nz

c----
c function values 
c----

      Do i=1,Nx+1
       Do j=1,Ny+1
        Do k=1,Nz+1
         x = ax + (i-1)*Dx
         y = ay + (j-1)*Dy
         z = az + (k-1)*Dz
         call integral_quadri_fun (x,y,z,f(i,j,k),menu)
        End Do
       End Do
      End Do

c----
c integral
c----
     
      res = 0.0;
      Do i=1,Nx+1
       Do j=1,Ny+1
        Do k=1,Nz+1
         res = res + weight(i,j,k)*f(i,j,k)
        End Do
       End Do
      End Do

      res = res*Dx*Dy*Dz

c----
c etc
c----

      if(menu.eq.1) then
       res = -res/(32.0D0*pi**3)
      elseif(menu.eq.2) then
       res = -res/(4.0D0*pi**3)
      elseif(menu.eq.3) then
       res = -res/(64.0D0*pi**3)
      elseif(menu.eq.4) then
       res = -res/(8.0D0*pi**3)
      end if

c----
c final
c----

      write (6,*) 'numer value: ',res

      Go to 98

c---
c done
c---

  99  Continue

 100  Format(1x,f15.10)
 102  Format(64(1x,f5.2))

      Stop
      End
