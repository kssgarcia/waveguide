      program msum

c--------------------
c compute the sum of the elements of a matrix
c in two different ways for timing purposes
c------------------

      Dimension A(2000,2000)

c---
c define the matrix
c---

      Do i=1,2000
       Do j=1,2000
        A(i,j)=0.10;
       End Do
      End Do

c---
c recall by row or column
c---

      sum = 0.0
      Do j=1,2000
       Do i=1,2000
c       sum = sum+A(i,j)
        sum = sum+A(j,i)
       End  Do
      End Do

c---
c done
c---
 
      Stop
      End
