      subroutine trgl6_sqr
     +
     +  (ndiv,npts,nelm
     +  )

c==========================================
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement
c==========================================

c--------------------------------------------
c Triangulation of a square with unit sides
c in the xy plane
c
c The patch is confined inside: -side < x,y < side
c
c Will subdivide an eight-element pattern
c into six-node triangles.
c
c SYMBOLS:
c -------
c
c  ndiv .... level of discretization
c            nvid = 0 gives 8 elements
c
c  npts .... number of nodes
c  nelm .... number of surface elements
c
c  x(i,j), y(i,j) .... Cartesian coordinates of local node j
c                      on element i
c                      j = 1,...,6
c                      i = 1,...,nelm
c
c  p(i,j) .... Cartesian coordinates of global node i
c              where j=1,2 with
c                               x = p(i,1)
c                               y = p(i,2)
c
c  n(i,j) .... global node number of local node number j on element i,
c              where j=1,...,6
c
c  ne(i,j) ... ne(i,1) is the number of elements touching global node i.
c              ne(i,2:ne(i,1)) are the corresponding element labels 
c
c  nbe(i,j) .. label of element sharing side j of element i
c              where j = 1, 2
c
c Iedge: interior-edge-corner node index:
c -----
c
c   Iedge(i,1) = 0 if global node i is an interior node
c
c   Iedge(i,1) = 1 if global node i is an edge node
c                with one periodic image
c   Iedge(i,2): global index of image node
c
c   Iedge(i,1) = 3 if global node i is a corner node
c                with three periodic images
c   Iedge(i,2): global index of first image node
c   Iedge(i,3): global index of first image node
c   Iedge(i,4): global index of first image node
c--------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  x(512,6), y(512,6)
      Dimension xn(512,6),yn(512,6)

      Dimension   n(512,6)
      Dimension nbe(512,3)

      Dimension     p(1090,3)
      Dimension    ne(1090,9)
      Dimension Iedge(1090,4)

      Parameter (eps=0.000001)

      common/points/p,ne
      common/elmnts/n,nbe
      common/edgepoints/Iedge

c----------------------------------------
c  zeroth level discretization (8 elements)
c  nodes are set manually
c----------------------------------------

      side = 0.5D0
      nelm = 8

c---
c  corner nodes.... upper half of xz plane
c---

      x(1,1) = side   ! first element
      y(1,1) = 0.0

      x(1,2) = side
      y(1,2) = side

      x(1,3) = 0.0
      y(1,3) = 0.0

c---

      x(5,1) =-side  !  fifth element
      y(5,1) = 0.0

      x(5,2) =-side
      y(5,2) =-side

      x(5,3) = 0.0
      y(5,3) = 0.0

c---

      x(6,1) = 0.0   ! sixth element
      y(6,1) =-side

      x(6,2) = 0.0
      y(6,2) = 0.0

      x(6,3) = -side
      y(6,3) = -side

c---

      x(2,1) = 0.0    ! second element
      y(2,1) = side

      x(2,2) = 0.0
      y(2,2) = 0.0

      x(2,3)= side
      y(2,3)= side

c---

      x(4,1) =-side    ! fourth element
      y(4,1) = 0.0

      x(4,2) = 0.0
      y(4,2) = 0.0

      x(4,3) =-side
      y(4,3) = side

c---

      x(8,1) = side    ! eighth element
      y(8,1) = 0.0

      x(8,2) = 0.0
      y(8,2) = 0.0

      x(8,3) = side
      y(8,3) =-side

c---

      x(7,1) = 0.0    ! seventh element
      y(7,1) =-0.5

      x(7,2) = side
      y(7,2) =-side

      x(7,3) = 0.0
      y(7,3) = 0.0

c---

      x(3,1) = 0.0    ! third element
      y(3,1) = side

      x(3,2) =-side
      y(3,2) = side

      x(3,3) = 0.0
      y(3,3) = 0.0

c------------------------------------------
c compute the midpoints of the three sides
c of the eight parental elements
c
c mid-points are numbered 4, 5, 6
c------------------------------------------

      Do i=1,nelm

       x(i,4) = 0.5*(x(i,1)+x(i,2))
       y(i,4) = 0.5*(y(i,1)+y(i,2))

       x(i,5) = 0.5*(x(i,2)+x(i,3))
       y(i,5) = 0.5*(y(i,2)+y(i,3))

       x(i,6) = 0.5*(x(i,3)+x(i,1))
       y(i,6) = 0.5*(y(i,3)+y(i,1))

      End Do

c----------------------------------------
c compute node coordinates on each element 
c for discretization levels 1 through ndiv
c----------------------------------------

      if(ndiv.eq.0) Go to 98   ! skip if ndiv = 0

      Do i=1,ndiv

       num = 0       ! counts the new elements arising by sub-division

       Do j=1,nelm                         !  loop over old elements

c---
c assign corner points to sub-elements
c these are to become "new" elements
c---
        num = num+1                 !  first sub-element

        xn(num,1) = x(j,1) 
        yn(num,1)= y(j,1)

        xn(num,2) = x(j,4)
        yn(num,2) = y(j,4) 

        xn(num,3) = x(j,6)
        yn(num,3) = y(j,6)

        xn(num,4) = 0.5*(xn(num,1)+xn(num,2))
        yn(num,4) = 0.5*(yn(num,1)+yn(num,2))

        xn(num,5) = 0.5*(xn(num,2)+xn(num,3))
        yn(num,5) = 0.5*(yn(num,2)+yn(num,3))

        xn(num,6) = 0.5*(xn(num,3)+xn(num,1))
        yn(num,6) = 0.5*(yn(num,3)+yn(num,1))

        num = num+1                !  second sub-element

        xn(num,1) = x(j,4)
        yn(num,1) = y(j,4)

        xn(num,2) = x(j,2)
        yn(num,2) = y(j,2)

        xn(num,3) = x(j,5)
        yn(num,3) = y(j,5)

        xn(num,4) = 0.5*(xn(num,1)+xn(num,2))
        yn(num,4) = 0.5*(yn(num,1)+yn(num,2))

        xn(num,5) = 0.5*(xn(num,2)+xn(num,3))
        yn(num,5) = 0.5*(yn(num,2)+yn(num,3))

        xn(num,6) = 0.5*(xn(num,3)+xn(num,1))
        yn(num,6) = 0.5*(yn(num,3)+yn(num,1))

        num = num+1

        xn(num,1) = x(j,6)                !  third sub-element
        yn(num,1) = y(j,6)

        xn(num,2) = x(j,5)
        yn(num,2) = y(j,5)

        xn(num,3) = x(j,3)
        yn(num,3) = y(j,3)

        xn(num,4) = 0.5*(xn(num,1)+xn(num,2))
        yn(num,4) = 0.5*(yn(num,1)+yn(num,2))

        xn(num,5) = 0.5*(xn(num,2)+xn(num,3))
        yn(num,5) = 0.5*(yn(num,2)+yn(num,3))

        xn(num,6) = 0.5*(xn(num,3)+xn(num,1))
        yn(num,6) = 0.5*(yn(num,3)+yn(num,1))

        num = num+1

        xn(num,1) = x(j,4)                !  fourth sub-element
        yn(num,1) = y(j,4)

        xn(num,2) = x(j,5)
        yn(num,2) = y(j,5)

        xn(num,3) = x(j,6)
        yn(num,3) = y(j,6)

        xn(num,4) = 0.5*(xn(num,1)+xn(num,2))    ! mid points
        yn(num,4) = 0.5*(yn(num,1)+yn(num,2))

        xn(num,5) = 0.5*(xn(num,2)+xn(num,3))
        yn(num,5) = 0.5*(yn(num,2)+yn(num,3))

        xn(num,6) = 0.5*(xn(num,3)+xn(num,1))
        yn(num,6) = 0.5*(yn(num,3)+yn(num,1))

       End Do                      !  end of old element loop

c---------------------------------
c number of elements was increased
c by a factor of four
c---------------------------------

       nelm = 4*nelm

c---
c relabel the new points
c and place them in the master list
c---

       Do k=1,nelm
        Do l=1,6

         x(k,l) = xn(k,l)
         y(k,l) = yn(k,l)

         xn(k,l) = 0.0   ! zero just in case
         yn(k,l) = 0.0

        End Do
       End Do

c-----------
      End Do               !  end of discretization level loop
c-----------

c-----------------------------------

  98  Continue

c-----------------------------------
c Create a list of global nodes by looping 
c over all elements
c and adding nodes not already found in the list.
c
c Fill the connectivity table n(i,j) 
c containing node numbers of element points 1-6
c-----------------------------------

c---
c six nodes of the first element are set mannualy
c---

      p(1,1) = x(1,1)
      p(1,2) = y(1,1)

      p(2,1) = x(1,2)
      p(2,2) = y(1,2)

      p(3,1) = x(1,3)
      p(3,2) = y(1,3)

      p(4,1) = x(1,4)
      p(4,2) = y(1,4)

      p(5,1) = x(1,5)
      p(5,2) = y(1,5)

      p(6,1) = x(1,6)
      p(6,2) = y(1,6)

      n(1,1) = 1  ! first  node of first element is global node 1
      n(1,2) = 2  ! second node of first element is global node 2
      n(1,3) = 3
      n(1,4) = 4
      n(1,5) = 5
      n(1,6) = 6  ! sixth  node of first element is global node 6

      npts = 6

c---
c loop over further elements
c
c Iflag=0 will signal new global node
c---

      Do i=2,nelm          ! loop over element nodes
       Do j=1,6            ! loop over element nodes

        Iflag = 0          ! new node

         Do k=1,npts
          if(Dabs(x(i,j)-p(k,1)).le.eps) then
           if(Dabs(y(i,j)-p(k,2)).le.eps) then

             Iflag = 1    ! the node has been recorded previously
             n(i,j) = k   ! the jth local node of element i
                          ! is the kth global node 
           end if
          end if
         End Do
        
         if(Iflag.eq.0) then     ! record the node

          npts = npts+1          ! one more global node

          p(npts,1) = x(i,j)
          p(npts,2) = y(i,j)

          n(i,j) = npts  ! the jth local node of element i
                         ! is the new global mode
         end if

       End Do
      End Do                      !  end of loop over elements

c----------------------------------
c Generate the connectivity table
c    ne(i,j)
c for elements touching global node i
c
c  ne(i,j) ... ne(i,1) is the number of elements touching 
c                      global node i
c              ne(i,2:ne(i,1)) are the corresponding element labels
c----------------------------------

c---
c initialize
c---

      Do i=1,npts
       Do j=1,9
        ne(i,j) = 0
       End Do
      End Do 

c---
c loop over global nodes
c---

      Do i=1,npts 

       ne(i,1) = 0
       Icount = 1

       Do j=1,nelm         ! over elements
        Do k=1,6           ! over element nodes

         if(abs(p(i,1)-x(j,k)).le.eps) then
          if(abs(p(i,2)-y(j,k)).le.eps) then
            Icount = Icount+1
            ne(i,1) = ne(i,1)+1  ! one more element
            ne(i,Icount) = j
          end if
         end if

        End Do
       End Do

      End Do  !  over surface points

c------------------------------------------
c Create the connectivity table nbe(i,j) for 
c
c  nbe(i,j) .. label of element sharing side j of element i
c              where j = 1, 2, 3
c
c Testing is done with respect to the mid-points
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
         if(k.eq.i) Go to 91  ! skip self-element
         Do l=4,6             !  loop over mid-points
          If(abs(x(i,j)-x(k,l)).le.eps) then
           If(abs(y(i,j)-y(k,l)).le.eps) then
             nbe(i,jcount) = k
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

c-------------------------
c Generate the edge index
c
c loop over all nodes and set:
c
c Iedge(i,1) = 0 for interior nodes
c Iedge(i,1) = 1 for edge nodes
c--------------------------------

      Do i=1,npts
       Iedge(i,1) = 0
       if(abs(p(i,1)-0.5).lt.eps) Iedge(i,1) = 1
       if(abs(p(i,2)-0.5).lt.eps) Iedge(i,1) = 1
       if(abs(p(i,1)+0.5).lt.eps) Iedge(i,1) = 1
       if(abs(p(i,2)+0.5).lt.eps) Iedge(i,1) = 1
      End Do

c--------------------------------
c loop over the edge nodes and set:
c
c Iedge(i,1) = 1 for edge nodes
c Iedge(i,2): global label of image node
c
c Iedge(i,1) = 3 for corner nodes
c Iedge(i,2): global label of first image node
c Iedge(i,3): global label of first image node
c Iedge(i,4): global label of first image node
c--------------------------------

      Do 77 i=1,npts

       if(Iedge(i,1).eq.0) Go to 77

       Ic = 1      ! counter of image nodes

       Do 78 j=1,npts

       if(i.eq.j) Go to 78
       if(Iedge(j,1).eq.0) Go to 78

       test = abs(p(i,1)-p(j,1)+1.0)+abs(p(i,2)-p(j,2))
       if(test.lt.eps) then
        Ic = Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,1)-p(j,1)-1.0)+abs(p(i,2)-p(j,2))
       if(test.lt.eps) then
        Ic = Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)+1.0)+abs(p(i,1)-p(j,1))
       if(test.lt.eps) then
        Ic = Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)-1.0)+abs(p(i,1)-p(j,1))
       if(test.lt.eps) then
        Ic = Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)-1.0)+abs(p(i,1)-p(j,1)-1.0)
       if(test.lt.eps) then
        Ic = Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)-1.0)+abs(p(i,1)-p(j,1)+1.0)
       if(test.lt.eps) then
        Ic = Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)+1.0)+abs(p(i,1)-p(j,1)-1.0)
       if(test.lt.eps) then
        Ic = Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)+1.0)+abs(p(i,1)-p(j,1)+1.0)
       if(test.lt.eps) then
        Ic = Ic+1
        Iedge(i,Ic) = j
       end if

  78   Continue

       Iedge(i,1) = Ic-1

  77  Continue

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
