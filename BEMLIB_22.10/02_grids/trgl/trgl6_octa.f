      subroutine trgl6_octa 
     +
     +  (ndiv  ! discretization level
     +  ,npts  ! number of points
     +  ,nelm  ! number of elements
     +  )

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------------
c Triangulation of the unit sphere
c by the successive subdivision
c of a regular octahedron
c into six-node quadratic triangles.
c
c SYMBOLS:
c -------
c
c  ndiv .... level of discretization of an octahedron
c            nvid = 0 gives 8 elements
c
c  npts .... number of nodes
c  nelm .... number of surface elements
c
c  x(i,j), y(i,j), z(i,j) .... 
c
c    Cartesian coordinates of local node j
c    on element i, j = 1,...,6 and i = 1,...,nelm
c
c  p(i,j) .... 
c
c    Cartesian coordinates of global node i
c    where j=1,2,3, with: x = p(i,1), y = p(i,2), z = p(i,3)
c                                 
c
c  n(i,j) .... global node number of local node number j on element i,
c              where j=1,...,6
c
c  ne(i,j) ... ne(i,1) is the number of elements touching global node i.
c              ne(i,2:ne(i,1)) are the corresponding element labels 
c
c  nbe(i,j) .. label of element sharing side j of element i
c              where j = 1, 2, 3
c
c--------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  x(512,6), y(512,6), z(512,6)
      Dimension xn(512,6),yn(512,6),zn(512,6)
      Dimension p(1026,3)

      Dimension n(512,6),ne(1026,7),nbe(512,3)

      Parameter (eps=0.00000001)

      common/points/p,ne
      common/elmnts/n,nbe

c----------------------------------------
c Begin with the zeroth-level discretization (8 elements)
c
c Nodes are set manually on the unit sphere
c----------------------------------------

      nelm = 8

c---
c corner nodes.... upper half of xz plane
c---

      x(1,1)= 0.0D0   ! first element
      y(1,1)= 0.0D0
      z(1,1)= 1.0D0

      x(1,2)= 1.0D0
      y(1,2)= 0.0D0
      z(1,2)= 0.0D0

      x(1,3)= 0.0D0
      y(1,3)= 1.0D0
      z(1,3)= 0.0D0
c---
      x(5,1)= 1.0D0  !  fifth element
      y(5,1)= 0.0D0
      z(5,1)= 0.0D0

      x(5,2)= 0.0D0
      y(5,2)= 0.0D0
      z(5,2)=-1.0D0

      x(5,3)= 0.0D0
      y(5,3)= 1.0D0
      z(5,3)= 0.0D0
c---
      x(6,1)= 0.0D0   ! sixth element
      y(6,1)= 0.0D0
      z(6,1)=-1.0D0

      x(6,2)=-1.0D0
      y(6,2)= 0.0D0
      z(6,2)= 0.0D0

      x(6,3)= 0.0D0
      y(6,3)= 1.0D0
      z(6,3)= 0.0D0
c---
      x(2,1)=-1.0D0    ! second element
      y(2,1)= 0.0D0
      z(2,1)= 0.0D0

      x(2,2)= 0.0D0
      y(2,2)= 0.0D0
      z(2,2)= 1.0D0

      x(2,3)= 0.0D0
      y(2,3)= 1.0D0
      z(2,3)= 0.0D0

c---
c  corner points .... lower half xz plane
c---

      x(4,1)= 0.0D0    ! fourth element
      y(4,1)= 0.0D0
      z(4,1)= 1.0D0

      x(4,2)= 0.0D0
      y(4,2)=-1.0D0
      z(4,2)= 0.0D0

      x(4,3)= 1.0D0
      y(4,3)= 0.0D0
      z(4,3)= 0.0D0
c---
      x(8,1)= 1.0D0    ! eighth element
      y(8,1)= 0.0D0
      z(8,1)= 0.0D0

      x(8,2)= 0.0D0
      y(8,2)=-1.0D0
      z(8,2)= 0.0D0

      x(8,3)= 0.0D0
      y(8,3)= 0.0D0
      z(8,3)=-1.0D0
c---
      x(7,1)= 0.0D0    ! seventh element
      y(7,1)= 0.0D0
      z(7,1)=-1.0D0

      x(7,2)= 0.0D0
      y(7,2)=-1.0D0
      z(7,2)= 0.0D0

      x(7,3)=-1.0D0
      y(7,3)= 0.0D0
      z(7,3)= 0.0D0
c---
      x(3,1)=-1.0D0    ! third element
      y(3,1)= 0.0D0
      z(3,1)= 0.0D0

      x(3,2)= 0.0D0
      y(3,2)=-1.0D0
      z(3,2)= 0.0D0

      x(3,3)= 0.0D0
      y(3,3)= 0.0D0
      z(3,3)= 1.0D0

c------------------------------------------
c Compute the mid-points of the three sides
c of the 8 first-generation elements
c
c mid-points are numbered 4, 5, 6
c------------------------------------------

      Do i=1,nelm

       x(i,4) = 0.5D0*(x(i,1)+x(i,2))
       y(i,4) = 0.5D0*(y(i,1)+y(i,2))
       z(i,4) = 0.5D0*(z(i,1)+z(i,2))

       x(i,5) = 0.5D0*(x(i,2)+x(i,3))
       y(i,5) = 0.5D0*(y(i,2)+y(i,3))
       z(i,5) = 0.5D0*(z(i,2)+z(i,3))

       x(i,6) = 0.5D0*(x(i,3)+x(i,1))
       y(i,6) = 0.5D0*(y(i,3)+y(i,1))
       z(i,6) = 0.5D0*(z(i,3)+z(i,1))

      End Do

c---
c project onto the unit sphere
c---

       Do k=1,nelm
        Do l=1,6
         rad = Dsqrt(x(k,l)**2+y(k,l)**2+z(k,l)**2) 
         x(k,l) = x(k,l)/rad
         y(k,l) = y(k,l)/rad
         z(k,l) = z(k,l)/rad
       End Do
      End Do

      If(ndiv.eq.0) Go to 98    ! octahedron done

c-------------------------------------------
c Compute the local element node coordinates
c for discretization levels 1 through ndiv
c-------------------------------------------

      Do i=1,ndiv     ! loop over refinement levels

       nm = 0        ! counts the new elements arising in each refinement loop
                      ! 4 elements will be generated in each pass
       Do j=1,nelm    ! loop over old elements

c---
c assign corner points to sub-elements
c these will become the "new" elements
c---

        nm = nm+1

        xn(nm,1) = x(j,1)      !  first sub-element
        yn(nm,1) = y(j,1)
        zn(nm,1) = z(j,1)

        xn(nm,2) = x(j,4)
        yn(nm,2) = y(j,4) 
        zn(nm,2) = z(j,4)

        xn(nm,3) = x(j,6)
        yn(nm,3) = y(j,6)
        zn(nm,3) = z(j,6)

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2))
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2))
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2))

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3))
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3))
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3))

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1))
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1))
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1))

        nm = nm+1

        xn(nm,1) = x(j,4)      !  second sub-element
        yn(nm,1) = y(j,4)
        zn(nm,1) = z(j,4)

        xn(nm,2) = x(j,2)
        yn(nm,2) = y(j,2)
        zn(nm,2) = z(j,2)

        xn(nm,3) = x(j,5)
        yn(nm,3) = y(j,5)
        zn(nm,3) = z(j,5)

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2))
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2))
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2))

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3))
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3))
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3))

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1))
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1))
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1))

        nm = nm+1

        xn(nm,1) = x(j,6)     !  third sub-element
        yn(nm,1) = y(j,6)
        zn(nm,1) = z(j,6)
 
        xn(nm,2) = x(j,5)
        yn(nm,2) = y(j,5)
        zn(nm,2) = z(j,5)

        xn(nm,3) = x(j,3)
        yn(nm,3) = y(j,3)
        zn(nm,3) = z(j,3)

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2))
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2))
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2))

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3))
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3))
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3))

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1))
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1))
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1))

        nm = nm+1

        xn(nm,1) = x(j,4)     !  fourth sub-element
        yn(nm,1) = y(j,4)
        zn(nm,1) = z(j,4)

        xn(nm,2) = x(j,5)
        yn(nm,2) = y(j,5)
        zn(nm,2) = z(j,5)

        xn(nm,3) = x(j,6)
        yn(nm,3) = y(j,6)
        zn(nm,3) = z(j,6)

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2))    ! mid points
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2))
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2))

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3))
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3))
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3))

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1))
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1))
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1))

       End Do                      !  end of loop over old elements

c--------------------------------------
c number of elements has been increased
c by a factor of 4
c--------------------------------------

       nelm = 4*nelm

c-----------------------------------
c relabel the new points
c and place them in the master list
c----------------------------------

       Do k=1,nelm
        Do l=1,6

         x(k,l) = xn(k,l)
         y(k,l) = yn(k,l)
         z(k,l) = zn(k,l)

c--- project onto the sphere

         rad = Dsqrt(x(k,l)**2+y(k,l)**2+z(k,l)**2)
         x(k,l) = x(k,l)/rad
         y(k,l) = y(k,l)/rad
         z(k,l) = z(k,l)/rad

         xn(k,l) = 0.0D0   ! zero just in case
         yn(k,l) = 0.0D0
         zn(k,l) = 0.0D0

        End Do
       End Do

c----------

      End Do   !  end of refinement loop

c-----------------------------------------

 98   Continue

c-----------------------------------
c Generate a list of global nodes by looping 
c over all elements
c and adding nodes not found in the list.
c
c Fill in the connectivity table n(i,j) 
c containing node numbers of element points 1-6
c-----------------------------------

c---
c six nodes of the first element are
c entered mannualy
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

      n(1,1) = 1  ! first  node of first element is global node 1
      n(1,2) = 2  ! second node of first element is global node 2
      n(1,3) = 3
      n(1,4) = 4
      n(1,5) = 5
      n(1,6) = 6  ! sixth node of first element is global node 6

      npts = 6

c---
c loop over further elements
c
c Iflag=0 will signal a new global node
c
c n(i,j): global label on jth node on ith element
c---

      Do i=2,nelm        ! loop over elements
       Do j=1,6          ! loop over element nodes

        Iflag=0

         Do k=1,npts
          If(Dabs(x(i,j)-p(k,1)).le.eps) then
           If(Dabs(y(i,j)-p(k,2)).le.eps) then
            If(Dabs(z(i,j)-p(k,3)).le.eps) then

             Iflag = 1    ! the node has been recorded previously
             n(i,j) = k   ! the jth local node of element i
                          ! is the kth global node 
            End If
           End If
          End If
         End Do
        
         If(Iflag.eq.0) then     ! record the node

          npts = npts+1          ! one more global node

          p(npts,1) = x(i,j)
          p(npts,2) = y(i,j)
          p(npts,3) = z(i,j)

          n(i,j) = npts   ! the jth local node of element i
                          ! is the new global node 

         End If

       End Do
      End Do                      !  end of loop over elements

c----------------------------------
c Generate connectivity table: ne(i,j)
c for elements touching global node i
c
c ne(i,1) is the number of elements touching  global node i
c ne(i,j) for j=2, ..., ne(i,1)+1
c         are the corresponding element labels
c----------------------------------

c---
c initialize
c---

      Do i=1,npts
       Do j=1,7
        ne(i,j) = 0
       End Do
      End Do 

c---
c loop over global nodes
c---

      Do i=1,npts 

       ne(i,1) = 0
       Icount = 1

       Do j=1,nelm     ! loop over elements
        Do k=1,6       ! loop over element nodes

           If(n(j,k).eq.i) then
             ne(i,1) = ne(i,1)+1
             Icount =Icount+1
             ne(i,Icount) = j
           End If

        End Do
       End Do
     
      End Do        !  end of loop over global nodes

c------------------------------------------
c Generate connectivity table nbe(i,j) for 
c
c  nbe(i,j) .. label of element sharing side j of element i
c              where j = 1, 2, 3
c
c For boundary elements with only 2 neighbors,
c the array entry will be zero.
c------------------------------------------

c---
c initialize
c---

      Do i=1,nelm
       Do j=1,3
        nbe(i,j) = 0
       End Do
      End Do

c---
c loop over elements
c---

      Do i=1,nelm             !  loop over elements
       jcount=1
       Do j=4,6               !  loop over mid-points

        Do k=1,nelm           !  test element

         If(k.eq.i) Go to 91  !  skip the self-element

         Do l=4,6             !  loop over mid-points

          If(abs(x(i,j)-x(k,l)).le.eps) then
           If(abs(y(i,j)-y(k,l)).le.eps) then
            If(abs(z(i,j)-z(k,l)).le.eps) then
             nbe(i,jcount) = k
            End If
           End If
          End If

         End Do

 91      Continue

        End Do                  !  end of test element

        If(nbe(i,jcount).ne.0) then
         jcount = jcount+1
        End If

       End Do
      End Do                    !  end of loop over elements

c-------------------------------------------
c project points p(i,j) onto the unit sphere
c-------------------------------------------

      Do i=1,npts

       rad = Dsqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)

       p(i,1) = p(i,1)/rad
       p(i,2) = p(i,2)/rad
       p(i,3) = p(i,3)/rad

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

      return
      end
