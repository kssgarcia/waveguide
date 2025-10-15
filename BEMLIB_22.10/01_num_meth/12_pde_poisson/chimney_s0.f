      program chimney0

      Implicit Double Precision (a-h,o-z)

      Dimension x(500,500),y(500,500)
      Dimension T(500,500)

c--------------------------------------------
c Steady-state temperature distribution in a
c brick chimney wall
c
c Will assume: b = 2a, c = a
c
c Discretization: Dx = Dy = a/N
c                 (beta = 1)
c
c Solution by the point Gauss-Siedel method
c
c-------------------------------------------

c---------------
c Default values
c---------------

      a     = 0.35
      To    = 10.0
      Ti    = 300.0
      Niter = 300
      N     = 16

c     write (6,*)
c     write (6,*) " Please enter a"
c     write (6,*) " --------------"
c     read  (5,*) a

c     write (6,*)
c     write (6,*) " Please enter Touter"
c     write (6,*) " -------------------"
c     read  (5,*) To

c     write (6,*)
c     write (6,*) " Please enter Tinner"
c     write (6,*) " ---------------"
c     read  (5,*) Ti

c     write (6,*)
c     write (6,*) " Please enter the discretization level N"
c     write (6,*) " ---------------------------------------"
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

      N1 = N+1

      NN  = 2*N
      NN1 = NN+1

c------------
c grid points
c------------

      Do i=1,NN1
       xi = (i-1.0)*Dx
       Do j=1,N1
        x(i,j) = xi
        y(i,j) = (j-1.0)*Dy
       End Do
      End Do

c-----------
c initialize
c-----------
   
      Do i=1,NN1
       Do j=1,N1
        T(i,j) = 0.0D0
       End Do
      End Do

c----------
c left wall
c----------

       Do j=1,N1
        T(1,j) = To
       End Do

c-----------
c upper wall
c-----------

       Do i=1,NN1
        T(i,N1) = To
       End Do

c-----------
c lower wall: use the symmetry condition
c-----------

       Do i=1,N1
        T(i,1) = T(N1,N+2-i)
       End Do

       Do i=N1,NN1
        T(i,1) = Ti
       End Do

c-----------
c right wall
c-----------

       Do j=1,N1
        T(NN+2,j) = T(NN,j)
       End Do

c-----------
c Iterations
c-----------

       Icount = 1
       Jcount = 1

  97   Continue

       Do iter=1,Niter

        Do i=2,NN1
         Do j=2,N
          T(i,j) = 0.25*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1))
         End Do
        End Do

c---
c reset the bottom
c---

       Do i=1,N1
        T(i,1) = T(N1,N+2-i)
       End Do

c---
c reset the right 
c---

       Do j=1,N1
        T(NN+2,j) = T(NN,j)
       End Do

       write (6,102) Jcount,T(2,2)

       Icount = Icount+1
       Jcount = Jcount+1

       End Do

c------------------------
c end of batch iterations
c------------------------

       write (6,*)
       write (6,*) " Continue the iterations ?"
       write (6,*) " 1 for YES, 0 for NO"
       read  (5,*) Icon

       If(Icon.ne.0) Go to 97

c----------
c print out
c----------

       open (9,file="chimney_s0.out")

       Do i=1,NN1
        write (9,*) Ntwo
        write (9,100) x(i,1),y(i,1),Null
        write (9,100) x(i,N1),y(i,N1),Null
c       write (6,*) Ntwo
c       write (6,100) x(i,1),y(i,1),Null
c       write (6,100) x(i,N1),y(i,N1),Null
       End Do

       Do j=1,N1
        write (9,*) Ntwo
        write (9,100) x(1,j),y(1,j),Null
        write (9,100) x(NN1,j),y(NN1,j),Null
c       write (6,*) Ntwo
c       write (6,100) x(1,j),y(1,j),Null
c       write (6,100) x(NN1,j),y(NN1,j),Null
       End Do

       Do i=1,NN1
        write (9,*) N1
        Do j=1,N1
c        write (9,100) x(i,j),y(i,j),T(i,j)/300.0
         write (9,100) x(i,j),y(i,j),T(i,j)
        End Do
       End Do

       Do j=1,N1
        write (9,*) NN1
        Do i=1,NN1
c        write (9,100) x(i,j),y(i,j),T(i,j)/300.0
         write (9,100) x(i,j),y(i,j),T(i,j)
        End Do
       End Do

       write (9,*) Null
       close (9)

c-----
c Done
c-----

  100  Format (3(1x,f10.5))
  102  Format (1x,i4,3(1x,f12.8))

       Stop
       End
