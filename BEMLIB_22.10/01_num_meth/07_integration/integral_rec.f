      program integral_rec

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------------
c Integrate a function f(x,y)
c over a rectangular area
c      ax<x<bx, ay<y<by
c using the trapezoidal rule
c-------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension weight(1025,1025),f(1025,1025)

      common/lgf/m,m1,m2

c---------
c constants
c---------

      pi = 3.14159265358D0

c------------
c parameters
c------------
 
   98 Continue

      write (6,*)
      write (6,*) ' Menu of functions'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) '  1  for f(x,y) = (1-cos(m1*x+m2*y))'
      write (6,*) '               / (sin(0.5*x)^2 + Dsin(0.5*y)^2 ) '
      write (6,*) '  2  for f(x,y) = (1-cos(m1*x)*cos(m2*y))'
      write (6,*) '                / (sin(0.5*x)^2 + sin(0.5*y)^2 ) '
      write (6,*) '  3  for f(x,y) = (1-cos(m1*x)*cos(m2*y))'
      write (6,*) '                / (1-cos(x)*cos(y) ) '
      write (6,*) '  4  for f(x,y) = (1-cos(2*m*x)))'
      write (6,*) '                / (1-cos(x)*cos(y))'
      write (6,*) '  5  for f(x,y) = (1-cos((m1+m2)x+(m1-m2)y))'
      write (6,*) '                / (1-cos(x)*cos(y) ) '
      write (6,*) '  10 for f(x,y) = honeycomb lattice'
      write (6,*) '  20 for f(x,y) = ln(4-2*cosx-2*cosy)'
      write (6,*) '  30 for f(x,y) = ln(6-2*cosx-2*cosy-2cos(x-y))'
      write (6,*) '  40 for f(x,y) = ln(sqr(sin^2x+sin^2y+1)+sqrt(...'
      write (6,*) '  41 for a shed'
      write (6,*)
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
       m1   = 1
       m2   = 2
       Nx = 4*8*8
       Ny = 4*8*8
      end if
c----------

c----------
      if(menu.eq.2) then
       ax =-0.0D0*pi
       bx = 1.0D0*pi
       ay =-0.0D0*pi
       by = 1.0D0*pi
       m1   = 1
       m2   = 2
       Nx = 8*8
       Ny = 8*8
      end if
c----------

c----------
      if(menu.eq.3) then
       ax =-1.0D0*pi
       bx = 1.0D0*pi
       ay =-1.0D0*pi
       by = 1.0D0*pi
       m1   = 1
       m2   = 2
       Nx = 8*8
       Ny = 8*8
      end if
c----------

c----------
      if(menu.eq.4) then
       ax = 0.0D0*pi
       bx = 2.0D0*pi
       ay = 0.0D0*pi
       by = 2.0D0*pi
       m = 1
        Nx = 4*8
        Ny = 4*8
      end if
c----------

c----------
      if(menu.eq.5) then
       ax =-0.0D0*pi
       bx = 2.0D0*pi
       ay =-0.0D0*pi
       by = 2.0D0*pi
       m1   = 1
       m2   = 2
       Nx = 8*8
       Ny = 8*8
      end if
c----------

c----------
      if(menu.eq.10) then
       ax = 0.0D0*pi
       bx = 2.0D0*pi
       ay = 0.0D0*pi
       by = 2.0D0*pi
       m1   = 2
       m2   = 0
       Nx = 8*8
       Ny = 8*8
      end if
c----------

c----------
      if(menu.eq.20) then
       ax = 0.0D0*pi
       bx = 1.0D0*pi
       ay = 0.0D0*pi
       by = 1.0D0*pi
       Nx = 1*16
       Ny = 1*16
      end if
c----------

c----------
      if(menu.eq.30) then
       ax = 0.0D0*pi
       bx = 2.0D0*pi
       ay = 0.0D0*pi
       by = 2.0D0*pi
       Nx = 8*16
       Ny = 8*16
      end if
c----------

c----------
      if(menu.eq.40) then
       ax = 0.0D0*pi
       bx = 0.5D0*pi
       ay = 0.0D0*pi
       by = 0.5D0*pi
       Nx = 2*16
       Ny = 2*16
      end if
c----------

c----------
      if(menu.eq.41) then
       ax =  0.0D0
       bx = 32.0D0
       ay =  0.0D0
       by = 28.0D0
       Nx = 4*16
       Ny = 4*16
      end if
c----------

c==================================
c proceed with the integration rule
c==================================

c---
c prepare
c---

      Dx = (bx-ax)/Nx
      Dy = (by-ay)/Ny

c----
c weight matrix
c----

      Do i=1,Nx+1
       Do j=1,Ny+1
        weight(i,j) = 0.5D0;
       End Do
      End Do

      weight(1,1)       = 0.25D0;
      weight(1,Nx+1)    = 0.25D0;
      weight(Ny+1,1)    = 0.25D0;
      weight(Ny+1,Nx+1) = 0.25D0;

      Do i=2,Nx
       Do j=2,Ny
        weight(i,j) = 1.0D0;
       End Do
      End Do

c----
c print and test the weights
c----

c     Do i=1,Nx+1
c       write (6,102) (weight(i,j),j=1,Ny+1)
c     End Do

c     sum = 0.0D0
c     Do i=1,Nx+1
c      Do j=1,Ny+1
c       sum = sum+weight(i,j)
c      End Do
c     End Do

c     write (6,*) sum,Nx*Ny

c----
c function values 
c----

      Do i=1,Nx+1
       Do j=1,Ny+1
        x = ax + (i-1)*Dx
        y = ay + (j-1)*Dy
        call integral_rec_fun (x,y,f(i,j),menu)
       End Do
      End Do

c----
c print the function values 
c----

c     Do i=1,Nx+1
c       write (6,102) (f(i,j),j=1,Ny+1)
c     End Do

c----
c integral
c----
     
      res = 0.0;
      Do i=1,Nx+1
       Do j=1,Ny+1
        res = res + weight(i,j)*f(i,j)
       End Do
      End Do
      res = res*Dx*Dy

c----
c etc
c----

      if(menu.eq.2) then
       res = 4.0D0*res
      end if

c----
c etc
c----

      if(menu.eq.1.or.menu.eq.2) then
         res = -res/16.0D0/pi**2
      end if

      if(menu.eq.4) then
         res = -res/16.0D0/pi**2
      end if

      if(menu.eq.5) then
         res = -res/16.0D0/pi**2
      end if

      if(menu.eq.10) then
         res = -res/16.0D0/pi**2
      end if

      if(menu.eq.20) then
       res = res/pi**2
      end if

      if(menu.eq.30) then
       res = res/(4.0*pi**2)
      end if

      if(menu.eq.40) then
       res = 8*res/pi**2
      end if

c----
c final
c----

      write (6,*) 'Integral: ',res

      Go to 98

c---
c done
c---

  99  Continue

 100  Format(1x,f15.10)
 102  Format(64(1x,f5.2))

      stop
      end
