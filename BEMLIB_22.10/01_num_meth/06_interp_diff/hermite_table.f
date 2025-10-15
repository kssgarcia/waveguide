      subroutine hermite_coeff 
     +                         (n,m,L,xp
     +                         ,f,x,t,c
     +                         )

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
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press
c
c 1998
c------------------------------------------------

c----------------------------------------
c
c  HERMITE INTERPOLATION
c
c  Coefficients of the generalized interpolating
c  polynomial for prescribed data xp
c  by the method of modified Newton 
c  differences described in section 6.7
c
c  SYMBOLS:
c  -------
c
c  xp ... x coordinate of prescribed data 
c  f .... ordinate and derivatives of prescribed data
c  n .... n+1 equals number of prescribed data
c  m .... number of specified derivatives each point
c  c .... computed polynomial coefficients
c  eps... the maximum difference
c         for which two values are considered equal.
c
c  u .... ordinates and derivatives at prescribed data
c  t .... x coord.
c  G .... difference table element (k=0)
c  GG ... difference table element (k>0)
c
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension m(10)
      Dimension xp(50),f(50,0:4),u(80),c(0:50),t(80)
      Dimension G(80),GG(80,80)
      Integer   r

      Parameter (eps=0.00001)

c--------
c prepare
c--------

      n1 = n+1

c---
c compute L: order of the interpolating polynomial
c
c Table has L+1 entries
c---
	
      L = -1+m(1)
      Do i=2,n1
        L = L+m(i)
      End Do

c---
c first two columns of the table 
c ... t(i) and u(t(i))
c---

      icount = 0

      Do i=1,n1
        Do j = 0,m(i)-1
          icount = icount+1
          t(icount) = xp(i)
          u(icount) = f(i,j)
        End Do
      End Do

c---
c form the k=0 column of the table
c ... G(t(i))
c---

      k = 0
      icount=0
	
      Do i=1,n1
        Do j=0,m(i)-1  
          icount    = icount+1
          G(icount) = f(i,0)
        End Do
      End Do

c---------------------------------
c form the k=1 column of the table
c ... G(t(j),t(j+1))
c---------------------------------

      k      = 1
      fact   = 1
      icount = 1
	
      Do 1 i=1,L+1-k         !  begin loop over column k=1
	   
	icount=icount+1
	 
c---
c compute GG(i,icount)
c---
	 
        If(abs(t(icount)-t(i)).lt.eps) then 

c---
c find r, the minimum index for which t(r)=t(i)
c---
	    
	  jcount=i
	  jflag=0
	  
          Do while(jflag.eq.0)
	     
            jcount=jcount-1
	   
            If(jcount.eq.0) then
              r=1  
              jflag=1
            End If

            If(abs(t(jcount)-t(i)).gt.eps) then
              r = jcount+1
              jflag = 1
	     End If

           End Do
	
           GG(i,icount) = u(r+k)/fact

	 Else
           GG(i,icount)=(G(icount)-G(i))/(t(icount)-t(i))  
        End If
	
  1   Continue             !  end of loop over column k=1           
c----------------------------------------------
c  form the remaining columns k=2 through k = L
c----------------------------------------------
	
      Do 2 k=2,L
	   
       fact=1
       Do i=1,k
         fact=fact*real(i)
       End Do

       Do i=1,L+1-k

       If(abs(t(i+k)-t(i)).lt.eps) then
    
c---
c  find r, the minimum index for which t(r)=t(i)
c---

       jcount=i

       jflag=0
	  
         Do while(jflag.eq.0)
	     
           jcount=jcount-1
	   
           If(jcount.eq.0) then
             r=1  
             jflag=1
            End If

            If(abs(t(jcount)-t(i)).gt.eps) then
             r=jcount+1
             jflag=1
            End If

           End Do
     
           GG(i,i+k) = u(r+k)/fact

         Else  
           GG(i,i+k)=(GG(i+1,i+k)-GG(i,i+k-1))/(t(i+k)-t(i))

         End If

       End Do

  2   Continue

c---
c read the top line of the difference table for the
c coefficients of the interpolation polynomial.
c---
	
      c(0) = G(1)
      Do k = 1,L
        c(k) = GG(1,k+1)
      End Do

c-----
c Done
c-----

      Return
      End
