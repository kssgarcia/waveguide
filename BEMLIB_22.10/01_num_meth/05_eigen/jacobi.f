      subroutine jacobi 
     +
     +    (n,a,b,u
     +    ,max
     +    ,Iflag
     +    ,icount
     +    )

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
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c---------------------------------------------------
c  Diagonalize a real symmetric matrix by 
c  Jacobi iterations, amounting to 
c  repeated plane rotations
c
c  Algorithms (5.7.5), (5.7.6), (5.7.7)
c
c  SYMBOLS:
c  -------
c
c  a .... square symmetric matrix
c  n .... size (rows/columns) of matrix a
c  b .... iterated matrix (tends to diagonal)
c  t .... temporary matrix
c  u .... matrix of eigenvectors
c
c  max .. maximum number of iterations
c  iflag. flag for convergence to diagonal form
c         (1 = converged)
c
c  tol... tolerance for skipping small elements
c         and stoping the rotations
c
c  tan_th tangent of rotation angle theta
c  sin_th sin of rotation angle theta
c  cos_th cos of rotation angle theta
c
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),t(10,10),b(10,10),u(10,10)

      Parameter (tol=0.00000001)

c-----------
c initialize
c-----------

      Iflag  = 0
      Icount = 0

c---------------------------------------
c copy matrix a onto b, that will become
c a similar diagonal matrix
c---------------------------------------

      Do i=1,n
        Do j=1,n
          b(i,j)=a(i,j)
        End Do
      End Do

c-------------------------------------
c introduce the matrix of eigenvectors
c and initialize to the unit matrix
c-------------------------------------

      Do i=1,n
        Do j=1,n
          u(i,j) = 0.0D0
        End Do
        u(i,i) = 1.0D0
      End Do

c----------------------------------------
c loop over the upper triangular elements
c until convergence or until a specified
c max number of iterations
c----------------------------------------

      Do 95 l=1,max           !  begin iteration loop

       icount = icount + 1

       Do i=1,n-1             !  loop over matrix rows
        Do j=i+1,n            !  loop over matrix columns

         If(abs(b(i,j)).gt.tol) then  !  skip small elements 

         w = 0.5D0*(b(i,i)-b(j,j))/b(i,j)

         If(w.ge.0) then
          sign= 1.0D0
         Else
          sign=-1.0D0
         End If

         tan_th = sign/(abs(w)+sqrt(w**2+1.0))
         cos_th = 1.0D0/(Dsqrt(tan_th**2+1.0))
         sin_th = cos_th*tan_th

         r = sin_th/(1.0D0+cos_th)

c-------------
c modify the elements of the
c temporary matrix t
c in the i,j rows and columns
c-------------

      Do k=1,n
c---
c rows
c---
        If(k.eq.i) then
            t(i,k)=b(i,i)+b(i,j)*tan_th 
        Else If(k.eq.j) then
            t(i,k)=0.0D0
        Else
            t(i,k)=b(i,k)+(b(j,k)-r*b(i,k))*sin_th
        End If

      End Do
 
c---
c columns
c---

      Do k=1,n

         If(k.eq.j) then
            t(k,j)=b(j,j)-b(i,j)*tan_th
         Else If(k.eq.i) then
            t(k,j)=0.0
         Else
            t(k,j)=b(k,j)-(b(i,k)+r*b(j,k))*sin_th
         End If

      End Do
 
c-------
c compute norm of     diagonal elements of b: d1
c compute norm of off-diagonal elements of b: s1
c-------

         d1 = 0.0D0
         s1 = 0.0D0

         Do k=1,n
          d1 = d1+b(k,k)**2
          Do m=k+1,n
           s1 = s1+b(k,m)**2
          End Do
         End Do

c---
c update the modified rows and columns
c---

         Do k=1,n
           b(i,k) = t(i,k)
           b(k,i) = b(i,k)
         End Do

         Do k=1,n
           b(k,j) = t(k,j)
           b(j,k) = b(k,j)
         End Do

c---
c update the eigenvector matrix
c---

         Do k=1,n
          tmp1 = u(k,i)
          tmp2 = u(k,j)
          u(k,i) =   tmp1*cos_th + tmp2*sin_th
          u(k,j) = - tmp1*sin_th + tmp2*cos_th
         End Do

        End If               !  end of small element check

c      write (6,*)
c      Do iprint=1,n
c         write (6,102) (b(iprint,jprint),jprint=1,n)
c      End do
c      pause

       End Do                !  end of column loop
      End Do                !  end of row loop

c-----
c End of a pass
c
c compute norm of     diagonal elements of b: d2
c         norm of off-diagonal elements of b: s2
c
c compute and print convergence diagnostics
c-----

      d2 = 0.0D0
      s2 = 0.0D0

      Do k=1,n
        d2 = d2+b(k,k)**2
        Do m=k+1,n
         s2 = s2+b(k,m)**2
        End Do
      End do

      diff1 = - 2.0D0*(s2-s1)
      diff2 =  d2-d1

      write (6,*)
      write (6,101) icount,s2
      write (6,105) diff1,diff2

c---
c check for stopping
c---

      If(s2.lt.tol) then
       Iflag=1
       Return
      End If

  95  Continue                  !  end of iteration loop

c-----
c Done
c-----

  101 Format (1x,"pass: ",i3,1x,"error: ",f15.10)
  102  Format (20(1x,f8.4))
  105 Format (1x,"Diagnostics :",2(1x,f8.5))

      Return
      End
