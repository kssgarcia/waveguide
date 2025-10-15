      program richardson

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
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c------------------------------------------------
c Perform the Richardson extrpolation from
c a geometric series of data
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension p(20),b(100),A(100,100)

      write (6,*)
      write (6,*) " Enter epsilon and q"
      write (6,*) " -------------------"
      read  (5,*) eps,q

      write (6,*) " Enter m"
      write (6,*) " In this example, it should be less than 6"
      write (6,*) " -----------------------------------------"
      read  (5,*) m
      
c-----
c define sample exponents
c-----

      p(1) = 1.3 D0
      p(2) = 2.2 D0
      p(3) = 3.1 D0
      p(4) = 3.4 D0

c-----
c generate data 
c-----

      h = eps

      Do i=1,m
       b(i) = 1.234 +1.2*h**p(1)
     +              -1.4*h**p(2)
     +              +1.5*h**p(3)
     +              -1.9*h**p(4)
       h = h/q
      End Do

c----
c perform the extrapolation
c----

c----
c first column of A
c-----

      Do i=1,m
       A(i,1) = b(i)
      End Do
     
c-----
c extrapolation 
c-----

      Do k=2,m
       omega = q**p(k-1)
       Do i = k,m
        A(i,k) = (omega*A(i,k-1)-A(i-1,k-1))/(omega-1.0)
       End Do
      End Do

      Do i=1,m
        write (6,100) (A(i,j),j=1,m)
      End Do

c-----
c Done
c-----

 100  Format (100(1x,f10.5))

      Stop
      End
