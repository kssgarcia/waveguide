      program chimney_s1

      Implicit Double Precision (a-h,o-z)

      Dimension x(0:500,0:500),y(0:500,0:500)
      Dimension T(0:500,0:500)

      Double Precision k

c-------------------------------------------
c Steady-state temperature distribution
c in a chimney wall
c
c Will assume: b = 2a, c = a
c
c Discretization: Dx = Dy = a/N
c                 (beta = 1)
c
c Solution by the point Gauss-Siedel method
c
c SYMBOLS:
c -------
c
c k: thermal conductivity
c
c-------------------------------------------

c---------------
c Default values
c---------------

      a      = 0.35
      ho     = 10.0
      hi     = 25.0
      k      = 0.72D0
      Tinf_o = 20.0
      Tinf_i = 300.0
      N      = 16

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

c     write (6,*)
c     write (6,*) " Please enter N"
c     write (6,*) " --------------"
c     read  (5,*) N

      write (6,*)
      write (6,*) " Please enter the number of iterations"
      write (6,*) " -------------------------------------"
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
       xi = (i-1.0)*Dx
       Do j=0,N+2
        x(i,j) = xi
        y(i,j) = (j-1.0)*Dy
       End Do
      End Do

c-----------
c initialize
c-----------
   
      Do i=1,NN1
       Do j=1,N1
        T(i,j) = 20.0
       End Do
      End Do

c-----------
c Iterations
c-----------

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

       Do i=1,NN1
        T(i,N2) = T(i,N) - ho*2.0*Dy/k * (T(i,N+1)-Tinf_o)
       End Do

c------
c lower wall
c------

       Do i=1,N1
        T(i,0) = T(N2,N2-i)
       End Do

       Do i=N2,NN1
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

        Do i=1,NN1
         Do j=1,N1
         T(i,j) = 0.25*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1))
         End Do
        End Do

       write (6,102) Jcount,T(2,2)

       Icount = Icount + 1
       Jcount = Jcount + 1

       End Do

c------------------------
c end of batch iterations
c------------------------

       write (6,*)
       write (6,*) " Continue the iterations ?"
       write (6,*) "                  0 for NO"
       read  (5,*) Icon

       If(Icon.ne.0) Go to 97
       
c----------
c print out
c----------

       open (9,file="chimney_s1.out")

       Do i=1,NN1
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
  102  Format (1x,i4,3(1x,f12.8))

       Stop
       End
