      program drop_ax

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------------
c Hydrostatic shape of an axisymmetric
c sessile drop resting on a horizontal plate,
c or pendant drop hanging underneath a horizontal plate,
c for a specified volume and contact angle.
c
c The numerical procedure is described in Section 4.4 
c of Pozrikidis (1997, pp. 173-177)
c
c The shooting parameter "shp" is the mean curvature
c of the interface at the top of the drop.
c
c The updates are done using secant iterations.
c
c Legend:
c ------
c
c npts: number of points along the interface.
c x,s:  axial and radial coordinates of interfacial marker points.
c shp:  shooting parameter.
c Isp:  orientation index.
c gac:  magnitude of the acceleration of gravity.
c psi:  interface slope angle.
c error: residual of object function defined in eq.(4.4.19).
c tol:  toleration for shooting iterations.
c Itmax: maximum number of shootings.
c alpha: contact angle
c
c Gravity is directed toward the negative x-axis
c------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(256),s(256)
      Dimension shp(300),error(300)
      Parameter (tol=0.00001,Itmax=1000)

c----------
c constants
c----------

      pi = 3.14159 265358D0

      null = 0
      none = 1
      ntwo = 2

      zero = 0.0D0
      one  = 1.0D0
      onem =-1.0D0
      two  = 2.0D0
      twom =-2.0D0

c---
c preferences
c---

      write (6,*)
      write (6,*) " Choose the data input mode"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 to read data from file: drop_ax.dat"
      write (6,*) "       2 to type in data"
      write (6,*) " --------------------------------------------"
      read  (5,*) Ienrd

      If(Ienrd.eq.0) Go to 999

c------------------------
      If(Ienrd.eq.2) then
c------------------------

  96  write (6,*) 
      write (6,*) 'Please enter 1 for a sessile drop'
      write (6,*) '            -1 for a pendant drop'
      write (6,*) '                       0  to quit'
      write (6,*) '----------------------------------'
      read  (5,*) Jsp

      If(Jsp.eq.0) Go to 999

      If(Jsp.ne.1.and.Jsp.ne.-1.and.Jsp.ne.0) then
       write (6,*)
       write (6,*) "Unfortunate selection; please try again"
       Go to 96
      End If

      write (6,*)
      write (6,*) 'Enter the magnitude of the acceleration of gravity'
      write (6,*) '--------------------------------------------------'
      read  (5,*) gac

      write (6,*)
      write (6,*) 'Please enter the interfacial tension'
      write (6,*) '------------------------------------'
      read  (5,*) gamma

      write (6,*)
      write (6,*) 'Please enter the density of the drop'
      write (6,*) '------------------------------------'
      read  (5,*) rhod

      write (6,*)
      write (6,*) 'Please enter the density of the ambient fluid'
      write (6,*) '---------------------------------------------'
      read  (5,*) rhoa

      write (6,*)
      write (6,*) 'Please enter the drop volume'
      write (6,*) '----------------------------'
      read  (5,*) volume

      write (6,*)
      write (6,*) 'Please enter contact angle'
      write (6,*) '        in multiples of pi'
      write (6,*) '--------------------------'
      read  (5,*) alpha

      write (6,*)
      write (6,*) 'Enter the number of interfacial marker points'
      write (6,*) '---------------------------------------------'
      read  (5,*) npts

      write (6,*)
      write (6,*) 'Please enter a small value for the'
      write (6,*) '        the shooting parameter'
      write (6,*) '----------------------------------'
      read  (5,*) epsilon

c----------
       Else
c----------

        open (1,file="drop_ax.dat")

         read (1,*) Jsp          ! sessile or pendant
         read (1,*) gac          ! acceleration of gravity
         read (1,*) gamma        ! surface tension
         read (1,*) rhod         ! density of the drop
         read (1,*) rhoa         ! density of the ambient fluid
         read (1,*) volume       ! drop volume
         read (1,*) alpha        ! contact angle
         read (1,*) npts         ! number of interfacial markers
         read (1,*) epsilon      ! for the shooting method

        close (1)

c----------
      End If
c----------

c--------
c prepare
c--------

      npts1 = npts+1
      drho  = rhod-rhoa              ! density difference
      capls = gamma/(gac*abs(drho))  ! square of the capillary number
      Isp = 1                        ! Isp is an orientation index
      If(drho.lt.0) Isp = -Isp
      If(Jsp.eq.-1) Isp = -Isp
      alpha = alpha*pi

c----------------------------------
c To start, assume drop shape
c is a truncated sphere
c and compute the sphere radius "a"
c as a function of 
c volume and contact angle.
c----------------------------------

      cosa = Dcos(alpha)

      a = (3.0D0*volume/pi)/(2.0D0+cosa**3-3.0D0*cosa)
      a = a**(1.0D0/3.0D0)

      shp(1) = 2.0D0/a     ! twice the mean curvature

      write (6,*)
      write (6,*) " Shooting parameter is "
      write (6,*) " the drop mean curvature at top"
      write (6,*)
      write (6,*) " Initial value set to: ",shp(1)
      write (6,*)
c     write (6,*) " Please enter 1 to proceed or 0 to change it"
c     write (6,*) " -------------------------------------------"
c     read  (5,*) Ichange
      Ichange = 1

      If(Ichange.eq.0) then  
       write (6,*)
       write (6,*) " Please enter the new value"
       read  (5,*) shp(1)
      End If

c----------------------------
c print data to plot the wall
c located at x=0
c----------------------------

      open (4,file="drop_ax.out",status="unknown")

      write (4,100) Ntwo
      write (4,100) none,twom,zero
      write (4,100) ntwo,two ,zero

c=================
c return to repeat
c=================

  94  Continue

      write (6,*)
      write (6,120) alpha
      write (6,*)
      write (6,*) " Volume, shooting parameter, error"
      write (6,*)

c---
c Compute initial solution of the odes
c to start-up the secant method
c---

      dpsi = alpha/(npts1-1.0)

      Ic=1             ! counter

      call drop_ax_ode
     +
     +   (npts
     +   ,capls,Isp
     +   ,dpsi
     +   ,shp(Ic)
     +   ,x,s
     +   ,volume_sh
     +   )

      error(Ic) = volume_sh - volume
      err       = abs(error(Ic))

      write (6,100) Ic,volume_sh,shp(Ic),err

c-------------------------
c second start-up solution 
c-------------------------

      Ic = 2

      shp(2) = shp(1)+epsilon

      call drop_ax_ode
     +
     +   (npts
     +   ,capls,Isp
     +   ,dpsi
     +   ,shp(Ic)
     +   ,x,s
     +   ,volume_sh
     +   )

      error(Ic) = volume_sh - volume
      err       = abs(error(Ic))

      write (6,100) Ic,volume_sh,shp(Ic),err

c---------------------------------------
c iterate on shp using the secant method
c until convergence
c---------------------------------------

      Do while(err.gt.tol)

       Ic = Ic+1

       If(Ic.gt.Itmax) then
        write (6,*) 
        write (6,*) 'Sorry: Iterations did not converge'
        write (6,*) 
        Go to 97
       End If

c---
c secant updating
c---
     
      Icb = Ic-2
      Ica = Ic-1

      dedc = (error(Ica)-error(Icb))
     +      /(  shp(Ica)-  shp(Icb))

      shp(Ic) = shp(Ica)-error(Ica)/dedc

       call drop_ax_ode
     +
     +   (npts
     +   ,capls,Isp
     +   ,dpsi
     +   ,shp(Ic)
     +   ,x,s
     +   ,volume_sh
     +   )

      error(Ic) = volume_sh - volume
      err       = abs(error(Ic))

      write (6,100) Ic,volume_sh,shp(Ic),err

      End Do                           ! iteration on shp

c------------------
c printing session
c
c draw drop interface
c and its reflection
c-------------------

      alpha = alpha/pi   ! unscale for printing
      shift = x(npts1)   ! shift to reset the origin at the wall

      write (4,*) 2*npts+1,alpha,a,Jsp

      Ic=0

      Do k=0,npts-1
       i=npts1-k
       xpr =Jsp*(x(i)-shift)
       spr =-s(i)
       Ic=Ic+1
       write (4,100) Ic,spr,xpr
      End Do

      Do i=1,npts1
       xpr = Jsp*(x(i)-shift)
       spr = s(i)
       Ic=Ic+1
       write (4,100) Ic,spr,xpr
      End Do

      alpha = alpha*pi   ! scale

c-------------------------
c return for another shape
c-------------------------

 97   Continue

      write (6,*)
      write (6,*)  " Current contact angle/pi is:",alpha/pi
      write (6,*)
      write (6,*)  " Please enter the new contact angle"
      write (6,*)  "              in multiples of pi"
      write (6,*)  "              0 to quit"
      write (6,*)  " ----------------------------------"
      read  (5,*) alpha
      alpha=alpha*pi
 
      If(alpha.eq.0) Go to 99
 
      Go to 94

c-----
c Done
c-----

  99  Continue

 999  Continue

      write (4,*) null
      close (4)

 100  Format (1x,i3,5(2x,f15.10))
 120  Format ("Contact angle: ",f15.10)

      Stop
      End
