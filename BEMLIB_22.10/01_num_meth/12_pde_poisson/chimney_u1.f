      program chimney_u1

c-------------------------------------------
c Evolution of temperature distribution
c in a brick chimney wall
c
c We will assume: b = 2a, c = a
c
c Discretization: Dx = Dy = a/N
c
c                 (beta = 1)
c
c Time stepping by the explicit FTCS method
c
c SYMBOLS:
c -------
c
c k: thermal conductivity
c
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(0:500,0:500),y(0:500,0:500)
      Dimension T(0:500,0:500)

      Dimension H(0:500,0:500)  ! holder

      Double Precision k

c---------------
c default values
c---------------

      a      = 0.35
      ho     = 10.0
      hi     = 25.0
      k      = 0.72D0
      rho    = 1920.00
      cp     = 835.00
      Tinf_o = 20.0
      Tinf_i = 300.0
      N      = 4
      alpha  = 0.40

c---------------
c inquire input
c---------------
c
c     write (6,*)
c     write (6,*) " Please enter a"
c     write (6,*) " --------------"
c     read  (5,*) a

c     write (6,*)
c     write (6,*) " Please enter ho"
c     write (6,*) " ---------------"
c     read  (5,*) ho

c     write (6,*)
c     write (6,*) " Please enter hi"
c     write (6,*) " ---------------"
c     read  (5,*) hi

c     write (6,*)
c     write (6,*) " Please enter k"
c     write (6,*) " --------------"
c     read  (5,*) k

c     write (6,*)
c     write (6,*) " Please enter Tinf_o"
c     write (6,*) " -------------------"
c     read  (5,*) Tinf_o

c     write (6,*)
c     write (6,*) " Please enter Tinf_i"
c     write (6,*) " -------------------"
c     read  (5,*) Tinf_i

      write (6,*)
      write (6,*) " Please enter N"
      write (6,*) " --------------"
      read  (5,*) N

      write (6,*)
      write (6,*) " Please enter alpha"
      write (6,*) " ------------------"
      read  (5,*) alpha

      write (6,*)
      write (6,*) " Please enter the number of steps"
      write (6,*) " --------------------------------"
      read  (5,*) Niter

c----------
c constants
c----------

      Null = 0
      Ntwo = 2

c--------
c prepare
c--------

      Dx = a/N
      Dy = a/N

      N1  = N+1
      N2  = N+2

      NN  = 2*N
      NN1 = NN+1
      NN2 = NN+2

c------------
c grid points
c------------

      Do i=0,NN1
       xi = (i-1.0D0)*Dx
       Do j=0,N2
        x(i,j) = xi
        y(i,j) = (j-1.0D0)*Dy
       End Do
      End Do

c-----------
c initialize
c-----------
   
      Do i=1,NN1
       Do j=1,N1
        T(i,j) = 10.0D0
       End Do
      End Do

      open (4,file="chimney_u1.out1")

c------------------
c set the time step
c------------------

      Dt = alpha*rho*cp*Dx**2 /k

      time = 0.0

      Icount = 0
      Jcount = 0

      write (4,101) Icount,time,T(1,1)
      write (6,101) Icount,time,T(1,1)

c--------------
c Time stepping
c--------------

       Icount = 1
       Jcount = 1

  97   Continue

       Do iter=1,Niter

c------
c left wall
c------

       Do j=1,N1
        T(0,j) = T(2,j) - ho*2.0*Dx/k * (T(1,j)-Tinf_o)
       End Do

c------
c upper wall
c------

       Do i=1,NN+1
        T(i,N2) = T(i,N) - ho*2.0*Dy/k * (T(i,N+1)-Tinf_o)
       End Do

c------
c lower wall
c------

       Do i=1,N1
        T(i,0) = T(N2,N2-i)
       End Do

       Do i=N1+1,NN+1
        T(i,0) = T(i,2) - hi*2.0*Dy/k * (T(i,1)-Tinf_i)
       End Do

c------
c right wall: symmetry
c------

       Do j=1,N1
        T(NN2,j) = T(NN,j)
       End Do

c----
c Differential equation
c-----

        Do i=1,NN+1
         Do j=1,N1
          H(i,j) = T(i,j)
     +           + alpha*(T(i+1,j)-2.0*T(i,j)+T(i-1,j))
     +           + alpha*(T(i,j+1)-2.0*T(i,j)+T(i,j-1))
         End Do
        End Do

c-----
c Update
c-----

        Do i=1,NN+1
         Do j=1,N1
          T(i,j) = H(i,j)
         End Do
        End Do

c-----
c reset time and counters
c-----

       time = Icount*Dt/3600.0

       write (4,101) Jcount,time,T(1,1)
       write (6,101) Jcount,time,T(1,1)

       Icount = Icount + 1
       Jcount = Jcount + 1

       End Do

       write (6,*)
       write (6,*) " Continue the time stepping ?"
       write (6,*) "                    0 for NO"
       read  (5,*) Icon

       If(Icon.ne.0) Go to 97
       
c----------
c print out
c----------

       open (9,file="chimney_u1.out2")

       Do i=1,NN+1
        write (9,*) Ntwo
        write (9,100) x(i,1),y(i,1),Null-1.0
        write (9,100) x(i,N1),y(i,N1),Null-1.0
c       write (6,*) Ntwo
c       write (6,100) x(i,1),y(i,1),Null-1.0
c       write (6,100) x(i,N1),y(i,N1),Null-1.0
       End Do

       Do j=1,N1
        write (9,*) Ntwo
        write (9,100) x(1,j),y(1,j),Null-1.0
        write (9,100) x(NN1,j),y(NN1,j),Null-1.0
c       write (6,*) Ntwo
c       write (6,100) x(1,j),y(1,j),Null
c       write (6,100) x(NN1,j),y(NN1,j),Null
       End Do

       Do i=1,NN1
        write (9,*) N1
        Do j=1,N1
         write (9,100) x(i,j),y(i,j),T(i,j)
        End Do
       End Do

       Do j=1,N1
        write (9,*) NN1
        Do i=1,NN1
         write (9,100) x(i,j),y(i,j),T(i,j)
        End Do
       End Do

       write (9,*) Null
       close (9)

c-----
c Done
c-----

  100  Format (3(1x,f10.5))
  101  Format (1x,i4,3(1x,f10.5))

       Stop
       End
