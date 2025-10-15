      program ft_2d

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c        C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c------------------------------------------------

c-----------------------------------------
c Complex Fourier transform
c of a function of two variables
c in a rectangular domain,
c computed by the conventional method
c
c SYMBOLS:
c --------
c
c Lx:  length of the x domain of interpolation
c Ly:  length of the y domain of interpolation
c---------------------------------------------

      Implicit Double Precision (a-h, o-z)

      Double Precision Lx,kx,kxh
      Double Precision Ly,ky,kyh

      Integer px,py

      Dimension f(200,200)
      Dimension fv(200,200)
      Dimension cr(-200:200,-200:200)
      Dimension ci(-200:200,-200:200)
      Dimension qr(-200:200,-200:200)
      Dimension qi(-200:200,-200:200)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

      null = 0

c---
c input
c uncommend to enter mannualy
c---

      open (3,file="ft_2d.dat")
       read (3,*) ax,bx
       read (3,*) ay,by
       read (3,*) Nx,Ny
      close(3)

c     write (6,*)
c     write (6,*) " Enter limits of interpolation ax, bx"
c     write (6,*) " ------------------------------------"
c     read  (5,*) ax,bx

c     write (6,*) " Enter limits of interpolation ay, by"
c     write (6,*) " ------------------------------------"
c     read  (5,*) ay,by

c     write (6,*) " Enter Number of intervals  Nx and Ny"
c     write (6,*) " ------------------------------------"
c     read  (5,*) Nx,Ny

  98  Continue

      write (6,*)
      write (6,*) " MENU OF FUNCTIONS"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 to quit  "
      write (6,*) " 1 for f(x,y) = 1.0  "
      write (6,*) " 2 for f(x,y) = x**2+y**2 "
      write (6,*) " 3 for f(x,y) = sin(x**2+x*exp(y)) "
      write (6,*) " 4 for f(x,y) = x**2"
      write (6,*) " 5 for f(x,y) = y**2"
      write (6,*) " 6 for f(x,y) = exp(sin(2*pi*x)*exp(sin(2*pi*y) "
      write (6,*) " -------------------"
      read (5,*) menu

      If(menu.eq.0) Go to 99

c--------
c Prepare
c--------

      open (7,file="ft_2d.out")

      Nx1 = Nx+1
      Ny1 = Ny+1

      Lx = bx-ax
      Ly = by-ay

      hx  = Lx/(Nx1-1.0D0)
      hy  = Ly/(Ny1-1.0D0)

      kx  = pi2/Lx
      ky  = pi2/Ly

      kxh = kx*hx
      kyh = ky*hy

      ff = 1.0D0/(Nx*Ny)
 
c-------------------------
c generate the data points
c-------------------------

      Do i=1,Nx1
       x = ax +(i-1.0D0)*hx
       Do j=1,Ny1
         y = ay +(j-1.0D0)*hy
         if(menu.eq.1) f(i,j) = 1.0D0
         if(menu.eq.2) f(i,j) = x**2+y**2
         if(menu.eq.3) f(i,j) = sin(x**2+x*exp(y))
         if(menu.eq.4) f(i,j) = x**2
         if(menu.eq.5) f(i,j) = y**2
         if(menu.eq.6) f(i,j) = exp(sin(pi2*x))*exp(sin(pi2*y))
       End Do
      End Do
   
c---------------------
c Even/Odd Intervals ?
c---------------------

      If(mod(Nx,2).eq.0) then
         write (6,*)
         write (6,*) Nx," Nx is Even"
         Mx = Nx/2
      Else
         write (6,*)
         write (6,*) Nx," Nx is Odd"
         Mx = (Nx-1)/2
      End If

      If(mod(Ny,2).eq.0) then
         write (6,*)
         write (6,*) Ny," Ny is Even"
         My  = Ny/2
      Else
         write (6,*)
         write (6,*) Ny," Ny is Odd"
         My = (Ny-1)/2
      End If

c---------------------------------
c compute the Fourier coefficients
c---------------------------------

      Do py = -My,My
         Do mm = 1,Nx1

          qr(mm,py) = 0.5D0*(f(mm,1)+f(mm,Ny1))
          qi(mm,py) = 0.0D0

           Do j = 2,Ny
            arg   = py*(j-1.0D0)*kyh
            cs    = Dcos(arg)
            sn    = Dsin(arg)
            qr(mm,py) = qr(mm,py) + cs*f(mm,j)
            qi(mm,py) = qi(mm,py) + sn*f(mm,j)
           End Do
         End Do

      End Do

c---
c when Ny is even, correct the border elements
c---

      If(mod(Ny,2).eq.0) then
       Do i=1,Nx1
        qr(i,-My) = 0.5D0*qr(i,-My)
        qi(i,-My) = 0.5D0*qi(i,-My)
        qr(i, My) = 0.5D0*qr(i, My)
        qi(i, My) = 0.5D0*qi(i, My)
       End Do
      End If

      Do px=-Mx,Mx
        Do py=-My,My
          cr(px,py) = 0.5D0*(qr(1,py)+qr(Nx1,py))
          ci(px,py) = 0.5D0*(qi(1,py)+qi(Nx1,py))
          Do i=2,Nx
           arg = px*(i-1.0D0)*kxh
           cs  = Dcos(arg)
           sn  = Dsin(arg)
           cr(px,py) = cr(px,py) + cs*qr(i,py)-sn*qi(i,py)
           ci(px,py) = ci(px,py) + cs*qi(i,py)+sn*qr(i,py)
          End Do
          cr(px,py) = ff*cr(px,py)
          ci(px,py) = ff*ci(px,py)
        End Do
      End Do

c---
c when Nx is even,
c correct the border elements
c---

      If(mod(Nx,2).eq.0) then

        Do py=-My,My
         cr(-Mx,py) = 0.5*cr(-Mx,py)
         ci(-Mx,py) = 0.5*ci(-Mx,py)
         cr( Mx,py) = 0.5*cr( Mx,py)
         ci( Mx,py) = 0.5*ci( Mx,py)
        End Do

      End If

c-------------------------------
c print the Fourier coefficients
c-------------------------------

      write (6,*)
      write (6,*) " Fourier Coefficients"
      write (6,*) " --------------------"
      write (6,*)

      Do px=-Mx,Mx
        write (6,101) (cr(px,py),py=-My,My)
      End Do

      write (6,*)

      Do px=-Mx,Mx
        write (6,101) (ci(px,py),py=-My,My)
      End Do

c--------------------------
c Confirm the interpolation
c--------------------------

      write (6,*)
      write (6,*) " Original values"
      write (6,*) " ---------------"

      Do i=1,Nx1
         write (6,101) (f(i,j),j=1,Ny1)
      End Do

c---
c Reconstructed values
c---

      Do i=1,Nx1
       xh = (i-1.0D0)*hx
       Do j=1,Ny1
         yh = (j-1.0D0)*hy

         sumr = 0.0D0
         sumi = 0.0D0

         Do px=-Mx,Mx
           Do py=-My,My
             arg = -px*kx*xh-py*ky*yh
             cs  = Dcos(arg)
             sn  = Dsin(arg)
             sumr= sumr + cr(px,py)*cs-ci(px,py)*sn
             sumi= sumi + cr(px,py)*sn+ci(px,py)*cs
           End Do
         End Do

       fv(i,j) = sumr

       End Do
      End Do

c--------
c Display
c--------

      write (6,*)
      write (6,*) " Verified values:"
      write (6,*) " ----------------"

      Do i=1,Nx1
         write (6,101) (fv(i,j),j=1,Ny1)
         write (7,101) (fv(i,j),j=1,Ny1)
         write (7,102)
      End Do

      write (6,*)
      write (6,*) " The border values are different"
      write (6,*) " because of the implied periodicity"
      write (6,*) " of the Fourier representation"
      write (6,*) " -----------------------------"

c-----------------
c return to repeat
c-----------------

      Go to 98

c-----
c Done
c-----

 99   Continue

      write (7,101) Null
      close (7)

 100  Format (1x,f10.5,1x,f10.5,1x,f10.5)
 101  Format (52(1x,f8.5))
 102  Format (';')

      Stop
      End
