      program von_koch_dr

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c===============================================
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c===============================================

c-----------------------------------------------
c Generate a periodic von Koch line
c
c  SYMBOLS:
c  -------
c
c  m       :  order of the fractal line
c  RL      :  period 	
c  alpha   :  aspect ratio of spikes
c  mp	   :  number of vertices over one period
c  xfr,yfr :  coordinates of vertices
c  ifr	   :  vertex index  1 for a sharp corner
c  	                    2 for a blunt corner
c----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xfr(4001),yfr(4001),ifr(4001)

      common/Intg/ifr
      common/Real/xfr,yfr

c----------
c constants
c----------

      Null = 0

c------
c input
c------

      write (6,*) 
      write (6,*) " Enter the period L"
      write (6,*) " ------------------"
      read  (5,*) RL

      write (6,*) " Enter the generation index m;"
      write (6,*) " should be an integer"
      write (6,*) " --------------------"
      read  (5,*) m

  98  Continue

      write (6,*)
      write (6,*) " Enter alpha; 1/3 < alpha < 1.0 "
      write (6,*) " --------------------------------"
      read  (5,*) alpha

      If(alpha.le.(1.0D0/3.0D0).or.alpha.ge.1.0D0) then
        write (6,*) 
        write (6,*) " alpha must lie in (1/3, 1); please try again "
        Go to 98
      End If

c------------------
c generate the line
c------------------

      call von_koch (RL,m,alpha,mp)

c-------------------------------
c Print vertices in file PLOTDAT
c-------------------------------

      open (1,file="PLOTDAT")

      write (6,*)
      write (6,*) " Vertex coordinates and angle indices:"
      write (6,*)

      write (1,100) mp

      Do i=1,mp
        write (1,100) i,xfr(i),yfr(i),ifr(i)
        write (6,100) i,xfr(i),yfr(i),ifr(i)
      End Do

c-----
c Done
c-----

      write (1,100) Null
      close (1)

 100  Format (1x,i3,2(1x,f10.5),1x,i3)

      Stop
      End
