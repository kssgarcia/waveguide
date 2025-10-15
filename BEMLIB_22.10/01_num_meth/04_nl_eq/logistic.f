      program logistic

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------
c Ligistic mapping
c-----------------

      Implicit Double precision (a-h,o-z)

c----------
c constants
c----------

      null = 0

c------------
c preferences
c------------

      write (6,*)
      write (6,*) " Please enter the minimum value of lambda"
      write (6,*)
      write (6,*) " Recommended: 3.5"
      write (6,*) " ----------------"
      read (5,*) rlam_min

      write (6,*)
      write (6,*) " Please enter the maximum value of lambda"
      write (6,*)
      write (6,*) " Recommended: 3.99"
      write (6,*) " -----------------"
      read (5,*) rlam_max

      write (6,*)
      write (6,*) " Please enter the step size in lambda"
      write (6,*)
      write (6,*) " Recommended: 0.002"
      write (6,*) " ------------------"
      read (5,*) rlam_step

      write (6,*)
      write (6,*) " How many iterations ?"
      write (6,*)
      write (6,*) " Recommended: 300"
      write (6,*) " ----------------"
      read (5,*) Niter

      write (6,*)
      write (6,*) " How many initial iterants to skip ?"
      write (6,*)
      write (6,*) " Recommended: 100"
      write (6,*) " ----------------"
      read (5,*) Nskip

c--------
c prepare
c--------

      open (1,file="PLOTDAT")

c-----------------
c loop over lambda
c-----------------

      rlam = rlam_min

  98  Continue
 
      x = 0.4         ! initial value
 
      write (1,*) Niter-Nskip

      Do i=1,Niter
       x = rlam*x*(1.0-x)
       If(i.gt.Nskip) write (1,100) rlam,x
      End Do

      rlam = rlam + rlam_step

      If(rlam.gt.rlam_max) Go to 99

      Go to 98

c-----
c Done
c-----

   99 Continue

      write (1,*) null

      write (6,*)
      write (6,*) "Done"
      write (6,*)

  100 Format (1x,f5.3,1x,f6.4)
 
      Stop
      End
