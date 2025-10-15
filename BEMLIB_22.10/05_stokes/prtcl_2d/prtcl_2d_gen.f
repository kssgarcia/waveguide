      program prtcl_2d_gen

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------------
c This program generates a random distribution
c of elliptical particles in a square box
c by perturbing particle centers placed 
c at the vertices of a regular lattice
c
c To be used as input for prtcl_2d.f
c---------------------------------------------

      open (1, file="prtcl_2d_gen.dat")
       read (1,*) RL
       read (1,*) Ncol
       read (1,*) Nrows
       read (1,*) NE
       read (1,*) Itp
       read (1,*) axis1
       read (1,*) axis2
      close (1)

c----
      Go to 888
c----

      write (6,*) 
      write (6,*) " Please enter the box size L"
      write (6,*) " ---------------------------"
      read  (5,*) RL

      write (6,*) 
      write (6,*) " Please enter number of particles in a row"
      write (6,*) " -----------------------------------------"
      read  (5,*) Ncol

      write (6,*) 
      write (6,*) " Please enter number of rows"
      write (6,*) " ---------------------------"
      read  (5,*) Nrows

      write (6,*) 
      write (6,*) " Please enter the number of elements"
      write (6,*) " per particle"
      write (6,*) " -----------------------------------"
      read  (5,*) NE


      write (6,*)
      write (6,*) " Please enter the element shape index"
      write (6,*)
      write (6,*) "      1 for a polygonal contour"
      write (6,*) "      2 for a elliptical contour"
      write (6,*) " -------------------------------"
      read  (5,*) Itp

 97   Continue

      write (6,*) 
      write (6,*) " Please enter the first particle axis"
      write (6,*) " ------------------------------------"
      read  (5,*) axis1

      write (6,*) 
      write (6,*) " Please enter the second particle axis"
      write (6,*) " -------------------------------------"
      read  (5,*) axis2

 888  Continue

c-----------------
c open output file
c-----------------

      open (1,file="prtcl_2d_gen.out")

c--------------------------------------
c compute the nominal particle size "a"
c--------------------------------------

      a = axis1
      If(axis2.gt.a) a = axis2

      b = RL/(2.0D0*Ncol)      ! grid size

      d_max = b-a             ! maximum displacement

      If(d_max.lt.0) then
       write (6,*)
       write (6,*) ' Particle size ',b,' is too large'
       Go to 97
      End If

c--------------------------
c Generate the distribution
c--------------------------

      b2 = 2.0D0*b
      Icount = 1

      Idum = 1  ! for the random number generator

      y = b

      Do i=1,Ncol

        x = b

        Do j=1,Nrows

         pert = 2.0*(ran2(idum)-0.5) * d_max
         xr = x+pert

         pert = 2.0*(ran2(idum)-0.5) * d_max
         yr = y+pert

         tilt = 2.0*(ran2(idum)-0.5)

         write (6,100) Icount,axis1,axis2,xr,yr,tilt,NE,Itp
         write (1,100) Icount,axis1,axis2,xr,yr,tilt,NE,Itp

         x = x+b2
         Icount = Icount+1

        End Do

        y = y+b2

      End Do

c-----
c Done
c-----

  100 Format (1x,i2,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5
     +       ,1x,f10.5,1x,i3,1x,i3)

      Stop
      End
