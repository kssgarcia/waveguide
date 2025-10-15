      subroutine trgl6_octa
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------------
c Triangulation of the unit sphere
c by subdividing a regular octahedron
c into six-node quadratic triangles.
c
c-----------  variables   ---------
c
c  nelm .... number of surface elements
c  npts .... number of nodal points
c
c  x(i,j), y(i,j), z(i,j) .... Cartesian coords of point j on element i
c
c  p(i,j) .... (x,y,z) coords. of surface node labeled i  (j=1,2,3)
c              x = p(i,1)
c              y = p(i,2)
c              z = p(i,3)
c
c  n(i,j) .... node number of point j on element i, where j=1,...,6
c
c  ne(i,j) ... ne(i,1) is the number of elements touching node i.
c              ne(i,2:ne(i,1)) are the corresponding element labels 
c
c  nbe(i,j) .. label of element sharing side j of element i
c              where j = 1, 2, 3
c
c  ndiv ....  level of discretization of starting octahedron (0-ndiv)
c
c--------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  x(512,6), y(512,6), z(512,6)
      Dimension xn(512,6),yn(512,6),zn(512,6)
      Dimension p(1026,3)

      Dimension n(512,6),ne(1026,7),nbe(512,3)

      Parameter (eps=0.0000001)

      common/points/p,ne
      common/elmnts/n,nbe

c----------------------------------------
c  0'th level discretization (8 elements)
c  nodes are set manually
c----------------------------------------

      nelm = 8

c---
c  corner points .... upper half of xz plane
c---

      x(1,1)= 0.0
      y(1,1)= 0.0
      z(1,1)= 1.0

      x(1,2)= 1.0
      y(1,2)= 0.0
      z(1,2)= 0.0

      x(1,3)= 0.0
      y(1,3)= 1.0
      z(1,3)= 0.0
c---
      x(5,1)= 1.0
      y(5,1)= 0.0
      z(5,1)= 0.0

      x(5,2)= 0.0
      y(5,2)= 0.0
      z(5,2)=-1.0

      x(5,3)= 0.0
      y(5,3)= 1.0
      z(5,3)= 0.0
c---
      x(6,1)= 0.0
      y(6,1)= 0.0
      z(6,1)=-1.0 

      x(6,2)=-1.0
      y(6,2)= 0.0
      z(6,2)= 0.0

      x(6,3)= 0.0
      y(6,3)= 1.0
      z(6,3)= 0.0
c---
      x(2,1)=-1.0
      y(2,1)= 0.0
      z(2,1)= 0.0

      x(2,2)= 0.0
      y(2,2)= 0.0
      z(2,2)= 1.0

      x(2,3)= 0.0
      y(2,3)= 1.0
      z(2,3)= 0.0

c---
c  corner points .... lower half xz plane
c---

      x(4,1)= 0.0
      y(4,1)= 0.0
      z(4,1)= 1.0

      x(4,2)= 0.0
      y(4,2)=-1.0
      z(4,2)= 0.0

      x(4,3)= 1.0
      y(4,3)= 0.0
      z(4,3)= 0.0
c---
      x(8,1)= 1.0
      y(8,1)= 0.0
      z(8,1)= 0.0

      x(8,2)= 0.0
      y(8,2)=-1.0
      z(8,2)= 0.0

      x(8,3)= 0.0
      y(8,3)= 0.0
      z(8,3)=-1.0
c---
      x(7,1)= 0.0
      y(7,1)= 0.0
      z(7,1)=-1.0 

      x(7,2)= 0.0
      y(7,2)=-1.0
      z(7,2)= 0.0

      x(7,3)=-1.0
      y(7,3)= 0.0
      z(7,3)= 0.0

c---
      x(3,1)=-1.0
      y(3,1)= 0.0
      z(3,1)= 0.0

      x(3,2)= 0.0
      y(3,2)=-1.0
      z(3,2)= 0.0

      x(3,3)= 0.0
      y(3,3)= 0.0
      z(3,3)= 1.0

c---
c compute the mid-points of the sides
c numbered 4, 5, 6
c---

      Do i=1,nelm

       x(i,4)= 0.5*(x(i,1)+x(i,2))
       y(i,4)= 0.5*(y(i,1)+y(i,2))
       z(i,4)= 0.5*(z(i,1)+z(i,2))

       x(i,5)= 0.5*(x(i,2)+x(i,3))
       y(i,5)= 0.5*(y(i,2)+y(i,3))
       z(i,5)= 0.5*(z(i,2)+z(i,3))

       x(i,6)= 0.5*(x(i,3)+x(i,1))
       y(i,6)= 0.5*(y(i,3)+y(i,1))
       z(i,6)= 0.5*(z(i,3)+z(i,1))

      End Do

c----------------------------------------
c compute node coordinates on each element 
c for discretization levels
c 1 through ndiv
c----------------------------------------

      Do i=1,ndiv

       num = 1

       Do j=1,nelm                         !  loop over old elements

c---
c  assign corner points to sub-elements
c---

        xn(num,1)= x(j,1)                  !  first sub-element
        yn(num,1)= y(j,1)
        zn(num,1)= z(j,1)

        xn(num,2)= x(j,4)
        yn(num,2)= y(j,4) 
        zn(num,2)= z(j,4)

        xn(num,3)= x(j,6)
        yn(num,3)= y(j,6)
        zn(num,3)= z(j,6)

        xn(num,4)= 0.5*(xn(num,1)+xn(num,2))
        yn(num,4)= 0.5*(yn(num,1)+yn(num,2))
        zn(num,4)= 0.5*(zn(num,1)+zn(num,2))

        xn(num,5)= 0.5*(xn(num,2)+xn(num,3))
        yn(num,5)= 0.5*(yn(num,2)+yn(num,3))
        zn(num,5)= 0.5*(zn(num,2)+zn(num,3))

        xn(num,6)= 0.5*(xn(num,3)+xn(num,1))
        yn(num,6)= 0.5*(yn(num,3)+yn(num,1))
        zn(num,6)= 0.5*(zn(num,3)+zn(num,1))

        xn(num+1,1)= x(j,4)                !  second sub-element
        yn(num+1,1)= y(j,4)
        zn(num+1,1)= z(j,4)

        xn(num+1,2)= x(j,2)
        yn(num+1,2)= y(j,2)
        zn(num+1,2)= z(j,2)

        xn(num+1,3)= x(j,5)
        yn(num+1,3)= y(j,5)
        zn(num+1,3)= z(j,5)

        xn(num+1,4)= 0.5*(xn(num+1,1)+xn(num+1,2))
        yn(num+1,4)= 0.5*(yn(num+1,1)+yn(num+1,2))
        zn(num+1,4)= 0.5*(zn(num+1,1)+zn(num+1,2))

        xn(num+1,5)= 0.5*(xn(num+1,2)+xn(num+1,3))
        yn(num+1,5)= 0.5*(yn(num+1,2)+yn(num+1,3))
        zn(num+1,5)= 0.5*(zn(num+1,2)+zn(num+1,3))

        xn(num+1,6)= 0.5*(xn(num+1,3)+xn(num+1,1))
        yn(num+1,6)= 0.5*(yn(num+1,3)+yn(num+1,1))
        zn(num+1,6)= 0.5*(zn(num+1,3)+zn(num+1,1))

        xn(num+2,1)= x(j,6)                !  third sub-element
        yn(num+2,1)= y(j,6)
        zn(num+2,1)= z(j,6)

        xn(num+2,2)= x(j,5)
        yn(num+2,2)= y(j,5)
        zn(num+2,2)= z(j,5)

        xn(num+2,3)= x(j,3)
        yn(num+2,3)= y(j,3)
        zn(num+2,3)= z(j,3)

        xn(num+2,4)= 0.5*(xn(num+2,1)+xn(num+2,2))
        yn(num+2,4)= 0.5*(yn(num+2,1)+yn(num+2,2))
        zn(num+2,4)= 0.5*(zn(num+2,1)+zn(num+2,2))

        xn(num+2,5)= 0.5*(xn(num+2,2)+xn(num+2,3))
        yn(num+2,5)= 0.5*(yn(num+2,2)+yn(num+2,3))
        zn(num+2,5)= 0.5*(zn(num+2,2)+zn(num+2,3))

        xn(num+2,6)= 0.5*(xn(num+2,3)+xn(num+2,1))
        yn(num+2,6)= 0.5*(yn(num+2,3)+yn(num+2,1))
        zn(num+2,6)= 0.5*(zn(num+2,3)+zn(num+2,1))

        xn(num+3,1)= x(j,4)                !  fourth sub-element
        yn(num+3,1)= y(j,4)
        zn(num+3,1)= z(j,4)

        xn(num+3,2)= x(j,5)
        yn(num+3,2)= y(j,5)
        zn(num+3,2)= z(j,5)

        xn(num+3,3)= x(j,6)
        yn(num+3,3)= y(j,6)
        zn(num+3,3)= z(j,6)

        xn(num+3,4)= 0.5*(xn(num+3,1)+xn(num+3,2))    ! mid points
        yn(num+3,4)= 0.5*(yn(num+3,1)+yn(num+3,2))
        zn(num+3,4)= 0.5*(zn(num+3,1)+zn(num+3,2))

        xn(num+3,5)= 0.5*(xn(num+3,2)+xn(num+3,3))
        yn(num+3,5)= 0.5*(yn(num+3,2)+yn(num+3,3))
        zn(num+3,5)= 0.5*(zn(num+3,2)+zn(num+3,3))

        xn(num+3,6)= 0.5*(xn(num+3,3)+xn(num+3,1))
        yn(num+3,6)= 0.5*(yn(num+3,3)+yn(num+3,1))
        zn(num+3,6)= 0.5*(zn(num+3,3)+zn(num+3,1))

        num = num+4                ! four elements were generated

       End Do                      !  end of old element loop

       nelm=nelm*4

c---
c rename the new points
c and place them in the master list
c---

       Do k=1,nelm

        Do l=1,6

         x(k,l) = xn(k,l)
         y(k,l) = yn(k,l)
         z(k,l) = zn(k,l)

         xn(k,l) = 0.0   ! zero just in case
         yn(k,l) = 0.0
         zn(k,l) = 0.0

        End Do

       End Do

c----------

      End Do            !  end of level loop of i

c-----------------------------------
c Create a list of surface nodes by looping 
c over all elements
c and adding nodes not already found in the list.
c
c Fill the connectivity table n(i,j) 
c node numbers of element points 1-6
c-----------------------------------

c---
c first element is set mannualy
c---

      p(1,1) = x(1,1)
      p(1,2) = y(1,1)
      p(1,3) = z(1,1)

      p(2,1) = x(1,2)
      p(2,2) = y(1,2)
      p(2,3) = z(1,2)

      p(3,1) = x(1,3)
      p(3,2) = y(1,3)
      p(3,3) = z(1,3)

      p(4,1) = x(1,4)
      p(4,2) = y(1,4)
      p(4,3) = z(1,4)

      p(5,1) = x(1,5)
      p(5,2) = y(1,5)
      p(5,3) = z(1,5)

      p(6,1) = x(1,6)
      p(6,2) = y(1,6)
      p(6,3) = z(1,6)

      n(1,1) = 1
      n(1,2) = 2
      n(1,3) = 3
      n(1,4) = 4
      n(1,5) = 5
      n(1,6) = 6

      npts=6

      Do i=2,nelm                  !  loop over elements

       Do j=1,6                    ! loop over element nodes

        Iflag=0

         Do k=1,npts

          If(abs(x(i,j)-p(k,1)).le.eps) then
           If(abs(y(i,j)-p(k,2)).le.eps) then
            If(abs(z(i,j)-p(k,3)).le.eps) then

             Iflag = 1      ! the node has been previously recorded 
             n(i,j) = k     ! the jth local node of element i
                            ! is the kth global node 
            End If
           End If
          End If

         End Do
        
         If(Iflag.eq.0) then          ! record the node

          npts = npts+1         

          p(npts,1) = x(i,j)
          p(npts,2) = y(i,j)
          p(npts,3) = z(i,j)

          n(i,j) = npts

         End If

       End Do
      End Do                      !  end of loop over elements

c----------------------------------
c Generate connectivity table ne(i,j)
c for elements touching node i
c----------------------------------

      Do i=1,npts
       Do j=1,7
        ne(i,j) = 0
       End Do
      End Do 

      Do i=1,npts                 !  loop over global nodes

       ne(i,1) = 0

       icount = 1

       Do j=1,nelm                ! loop over local nodes
        Do k=1,6

         If(abs(p(i,1)-x(j,k)).le.eps) then
          If(abs(p(i,2)-y(j,k)).le.eps) then
           If(abs(p(i,3)-z(j,k)).le.eps) then

            icount=icount+1
            ne(i,1)=ne(i,1)+1
            ne(i,icount)=j

           End If
          End If
         End If

        End Do
       End Do
     
      End Do                       !  end of loop over surface points

c------------------------------------------
c  Create connectivity table nbe(i,j) for 
c  neighboring elements j of element i
c
c  Testing is done with respect to the mid-points
c
c (for boundary elements with only 2 neighbors,
c  the array entry will be zero)
c------------------------------------------

      Do i=1,nelm
       Do j=1,3
        nbe(i,j)=0
       End Do
      End Do

      Do i=1,nelm             !  loop over elements

       jcount=1

       Do j=4,6               !  loop over mid-points

        Do k=1,nelm           !  test element

         If(k.eq.i) Go to 91  !  loop over mid-points

         Do l=4,6

          If(abs(x(i,j)-x(k,l)).le.eps) then
           If(abs(y(i,j)-y(k,l)).le.eps) then
            If(abs(z(i,j)-z(k,l)).le.eps) then

             nbe(i,jcount)=k

            End If
           End If
          End If

         End Do

 91      Continue

        End Do                  !  end of test element

        If(nbe(i,jcount).ne.0) then
         jcount=jcount+1
        End If

       End Do
      End Do                    !  end of loop over elements

c---
c project points p(i,j) onto the unit sphere
c---

      Do i=1,npts

       r=sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)

       p(i,1) = p(i,1)/r
       p(i,2) = p(i,2)/r
       p(i,3) = p(i,3)/r

      End Do

c---------
c printing
c---------

c     write (6,*)
c     write (6,*) nelm,' grid elements'
c     write (6,*) npts,' grid points'
c     write (6,*)

c-----
c Done
c-----

      Return
      End
