        program qr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c               C. Pozrikidis
c ''Numerical Computation in Science and Engineering''
c       Oxford University Press
c------------------------------------------------

c----------------------------------
c QR decomposition of a real matrix
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),q(128,128),r(128,128),v(128,128)

  97  Continue

c----------------
c read the matrix
c----------------

      open (unit=8,file='matrix.dat')

       read (8,*) n
       Do i=1,n
         read(8,*) (a(i,j),j=1,n)
       End Do

      close(8)

      write (6,*)
      write (6,*) " original matrix:"
      write (6,*) " ----------------"

      Do i=1,n
        write (6,100) (a(i,j),j=1,n)
      End Do 

      nprint = 10
      Icount = 1

c---------------------
c Do the decomposition
c---------------------

      Do 87 irepeat=1,10*nprint

      call qr_rot (a,n,q,r)

      Do i=1,n
      Do j=1,n
       a(i,j)=0.0D0
       Do k=1,n
        a(i,j) = a(i,j)+r(i,k)*q(k,j)
       End Do
      End Do
      End Do

c------
c print
c------

      If(icount.eq.nprint) then

      write (6,*)
      write (6,*) " transformed matrix:"
      write (6,*) " -------------------"

      Do i=1,n
        write (6,100) (a(i,j),j=1,n)
      End Do 

      icount = 0

      End If

      icount=icount+1

  87  Continue

c-----
c Done
c-----

 99   Continue

 100  Format(20(3x,f10.5))

      Stop
      End
