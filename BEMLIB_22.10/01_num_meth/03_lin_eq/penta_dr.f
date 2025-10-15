      program penta_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c            C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press, 1998
c------------------------------------------------

c----------------------------------------
c Driver for penta-diagonal matrix solver
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension e(200),d(200),a(200),b(200),c(200),s(200)
      Dimension x(200)

c---------------------------
c read the diagonal vectors
c and the right-hand side
c---------------------------

      open (unit=8,file="penta.dat")

       read (8,*) n
       read (8,*) (a(i),i=1,n)
       read (8,*) (b(i),i=1,n-1)
       read (8,*) (c(i),i=1,n-2)
       read (8,*) (d(i),i=2,n)
       read (8,*) (e(i),i=3,n)
       read (8,*)
       read (8,*) (s(i),i=1,n)

      close (8)

c-----------------
c solve the system
c-----------------

      call penta (n,a,b,c,d,e,s,x)

      write (6,*)
      write (6,*) " ---------------"
      write (6,*) " Solution vector:"
      write (6,*)

      Do i=1,n
        write (6,100) i,x(i)
      End Do
 
c-------
c Verify
c-------

      write (6,*)
      write (6,*) " -----------------------"
      write (6,*) " Residuals:"

      i = 1
      res = s(1)-a(1)*x(1)-b(1)*x(2)-c(1)*x(3)
      write (6,101) i,res

      i = 2
      res = s(2)-d(2)*x(1)-a(2)*x(2)-b(2)*x(3)-c(2)*x(4)
      write (6,101) i,res

      Do i=3,n-2
        res = s(i)-e(i)*x(i-2)-d(i)*x(i-1)-a(i)*x(i)-b(i)*x(i+1)
     +            -c(i)*x(i+2)
        write (6,101) i,res
      End Do

      i = n-1
      res = s(n-1) - e(n-1)*x(n-3)-d(n-1)*x(n-2)-a(n-1)*x(n-1)
     +             - b(n-1)*x(n)
      write (6,101) i,res

      i = n
      res = s(n)   - e(n)  *x(n-2)-d(n)*x(n-1)-a(n)*x(n)
      write (6,101) i,res

c-----
c Done
c-----

 99   Continue

      write (6,*) " -----------------------"
      write (6,*) 
      write (6,*) "Thank you"

 100  Format (1x,i3,1x,f10.6)
 101  Format (1x,i3,2x,f15.10)

      Stop
      End
