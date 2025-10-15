      program rod

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c----------------------------------------
c Unsteady heat condition through a rod
c computed by the FTCS method
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(1000),f(1000),fnew(1000)

c---
c constants
c---

      null = 0

c---
c prepare
c---

      open (2,file="rod.out")

c---
c input
c---

      write (6,*)
      write (6,*) " Please enter the number of intervals"
      write (6,*) "-------------------------------------" 
      read  (5,*) N

      write (6,*)
      write (6,*) " Please enter time step Dt"
      write (6,*) "--------------------------" 
      read  (5,*) Dt

      write (6,*)
      write (6,*) " After how many steps a printout ?"
      write (6,*) "----------------------------------"
      read  (5,*) Nprint

c---
c prepare to run
c---

      arg = 30.0/Dt
      Ipause = int(arg)   ! will pause after Ipause steps

c---
c parameters
c---

      RL = 30
      Rk = 0.86

c---
c preparations
c---

      Dx  = RL/N
      al  = Dt*Rk/Dx**2    ! alpha
      al1 = (1.0D0-2.0D0*al)

      fnew(N+1) = 25.0D0

      write (6,140) al

      Do i=1,N+1
       x(i) = (i-1.0)*Dx
       f(i) = 25.0D0       ! initial condition
      End Do

      time   = 0
      Icount = 1    ! step conter
      Iprint = 0    ! printing counter

c---
c time stepping begins
c---

 98   Continue

      If(Icount.eq.Ipause) then

       write (6,*) 
       write (6,*) "Continue ? Enter 0 for NO"
       write (6,*) 
       read  (5,*) more

       If(more.eq.0) Go to 99
       If(more.eq.1) icount = 0

      End If

c---
c left-end bc
c---

      fnew(1) = 25.0*exp(-0.003*time)+80.0*(1.0-exp(-0.002*time))

      Do i=2,N
       fnew(i) = al*f(i-1) +al1*f(i)+al*f(i+1)
      End Do

      Do i=1,N+1
       f(i) = fnew(i)
      End Do

c---------
c printing
c---------

      If(iprint.eq.Nprint) then

       write (6,100) N+1,time
       Do i=1,N+1
        f(i) = fnew(i)
        write (6,100) i,x(i),f(i)
        write (2,100) i,x(i),f(i)
       End Do
       Iprint = 0

      End If

c---------------
c reset counters
c---------------

      icount = icount + 1
      iprint = iprint + 1
      time   = time+Dt

      Go to 98

c-----
c Done
c-----

  99  Continue

      write (2,100) null

 100  Format (1x,i3,3(1x,f10.5))
 140  Format ("alpha=",f10.5)

      Stop
      End
