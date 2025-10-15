      program sgf_3d_2p_dr

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c--------------------------------------
c Three-dimensional doubly periodic
c Green's function of Stokes flow
c
c The point forces are located in a plane
c that is perpendicular to the z axis
c
c LEGEND:
c ------
c
c method = 1 implements the Fourier series method
c method = 2 implements the fast summation method
c            of Pozrikidis (1996)
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (eps=0.001)

c----------
c constants
c----------

      Null = 0

c--------
c prepare
c--------

      open (1,file="sgf.out")

      write (6,*)
      write (6,*) "  Enter:"
      write (6,*)
      write (6,*) " 1 for one evaluation"
      write (6,*) " 2 for a profile along the z axis"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c------
c input
c------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to read the base vectors from file sgf_3d_2p.dat"
      write (6,*) " 2 to enter input"
      write (6,*) " 0 to quit"
      write (6,*) " ----------------"
c     read  (5,*) Ienrd
      Ienrd = 1

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (3,file="sgf_3d_2p.dat")
        read (3,*) a11,a12
        read (3,*) a21,a22
      close (3)
c---
      else
c---
       write (6,*)
       write (6,*) " Enter the coordinates "
       write (6,*) " of the first lattice vector"
       write (6,*) " ---------------------------"
       read  (5,*) a11,a12

       write (6,*) " Enter the coordinates "
       write (6,*) " of the second lattice vector"
       write (6,*) " ----------------------------"
       read  (5,*) a21,a22

c-----------
      end if
c-----------

c------------------------------------------------
c sgf_3d_2p ewald will produce the reciprocal vectors
c and the optimal value of xi
c according to Benakker
c------------------------------------------------

      call sgf_3d_2p_ewald
     +
     +   (a11,a12
     +   ,a21,a22
     +   ,b11,b12
     +   ,b21,b22
     +   ,ewald,area
     +   )

c---
c display
c---

      write (6,*)
      write (6,*) " Coordinates of the lattice vectors"
      write (6,*) " ----------------------------------"
      write (6,100) a11,a12
      write (6,100) a21,a22

      write (6,*)
      write (6,*) " Reciprocal lattice vectors"
      write (6,*) " --------------------------"
      write (6,100) b11,b12
      write (6,100) b21,b22

c---

      write (6,*)
      write (6,*) " Enter (x0,y0,z0) of a point force"
      write (6,*) " ---------------------------------"
c     read  (5,*) x0,y0,z0

      x0 = 0.0D0
      y0 = 0.0D0
      z0 = 0.0D0

c---
c continue with input
c---

 98   Continue

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for the straight Fourier series method"
      write (6,*) " 2 for the Pozrikidis fast method"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
c     read  (5,*) method
      method = 2

      if(method.eq.0) Go to 99

      if(method.ne.1.and.method.ne.2) then
       write (6,*)
       write (6,*) " Invalid selection; please try again"
       Go to 98
      end if

      write (6,*)
      write (6,*) " Enter the truncation level for summation"
      write (6,*) " in wave number space"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
c     read  (5,*) max2
      max2 = 5

      if(max2.eq.0) Go to 99

c-------------------------

      if(method.eq.2) then

         write (6,*)
         write (6,*) " Enter truncation level for summation"
         write (6,*) " in real space"
         write (6,*) " 0 to quit"
         write (6,*) " ---------"
c        read  (5,*) max1
         max1 = 5

         if(max1.eq.0) Go to 99

         write (6,*)
         write (6,*) " Enter the value of xi"
         write (6,*) " 99 for the recommended value ",ewald
         write (6,*) " 0 to quit"
         write (6,*) " --------"
         read  (5,*) ew

         if(ew.eq.0 ) Go to 99
         if(ew.eq.99) ew = ewald

      End If

c-----------

c---------------
c one evaluation
c---------------

      if(menu.eq.1) then

      write (6,*)
      write (6,*) " Enter the coordinates of the field point x,y,z"
      write (6,*) " ----------------------------------------------"
      read  (5,*) x,y,z

      call sgf_3d_2p 
     +
     +  (method
     +  ,eps
     +  ,x,y,z
     +  ,x0,y0,z0
     +  ,a11,a12,a21,a22
     +  ,b11,b12,b21,b22
     +  ,ew,area
     +  ,max1,max2
     +  ,Gxx,Gxy,Gxz
     +  ,Gyx,Gyy,Gyz
     +  ,Gzx,Gzy,Gzz
     +  )

      write (6,*) ' --------------'
      write (6,*) ' Green function '
      write (6,*) ' --------------'

      write (6,100)
      write (6,100) Gxx,Gxy,Gxz
      write (6,100) Gyx,Gyy,Gyz
      write (6,100) Gzx,Gzy,Gzz

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
      write (6,*) " Enter zmin and zmax"
      write (6,*) " -------------------"
      read  (5,*) zmin,zmax

      write (6,*)
      write (6,*) " Enter the number of profile points"
      write (6,*) " ----------------------------------"
      read  (5,*) Nprf

      Nprf1 = Nprf+1
      zstp = (zmax-zmin)/Nprf

      write (1,*) Nprf1

      Do i=1,Nprf1

        z = zmin+(i-1.0)*zstp

        call sgf_3d_2p 
     +
     +     (Method
     +     ,eps
     +     ,x,y,z
     +     ,x0,y0,z0
     +     ,a11,a12,a21,a22
     +     ,b11,b12,b21,b22
     +     ,ew,area
     +     ,max1,max2
     +     ,gxx,gxy,gxz
     +     ,gyx,gyy,gyz
     +     ,gzx,gzy,gzz
     +     )

        write (6,101) i,z,gxx,gxy,gxz,gyy,gyz,gzz
        write (1,101) i,z,gxx,gxy,gxz,gyy,gyz,gzz

        End Do

      End If

c-------------
      Go to 98    ! return to repeat
c-------------

c-----
c done
c-----

 99   Continue

      write (1,*) null
      close (1) 

  100 Format (3(1x,f20.10))
  101 Format (1x,i3,10(1x,f9.5))

      stop
      end
