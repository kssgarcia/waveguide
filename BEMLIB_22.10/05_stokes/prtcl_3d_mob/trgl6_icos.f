      subroutine trgl6_icos
     +
     + (Ndiv
     + ,Npts,Nelm
     + )

c=============================================
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=============================================

c--------------------------------------------
c Triangulation of the unit sphere
c by subdividing a regular icosahedron
c into six-node quadratic triangles
c
c SYMBOLS:
c -------
c
c  Ndiv .... level of discretization of icosahedron
c            Nvid = 0 gives 20 elements
c
c  Npts .... number of nodes
c  Nelm .... number of surface elements
c
c  x(i,j), y(i,j), z(i,j) .... Cartesian coordinates of local node j
c                              on element i
c                              j = 1,...,6
c                              i = 1,...,Nelm
c
c  p(i,j) .... Cartesian coordinates of global node i
c              where j=1,2,3, with
c                                 x = p(i,1)
c                                 y = p(i,2)
c                                 z = p(i,3)
c
c  n(i,j) .... global node number of local node number j on element i,
c              where j=1,...,6
c
c  ne(i,j) ... ne(i,1) is the number of elements touching global node i.
c              ne(i, 2:ne(i,1)+1) are the corresponding element labels 
c
c  nbe(i,j) .. label of element sharing side j of element i
c              where j = 1, 2, 3
c
c--------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension VX(12), VY(12), VZ(12)

      Dimension  x(512,6), y(512,6), z(512,6)
      Dimension xn(512,6),yn(512,6),zn(512,6)
      Dimension p(1026,3)

      Dimension n(512,6),ne(1026,7),nbe(512,3)

      Parameter (eps=0.00000001)

      common/points/p,ne
      common/elmnts/n,nbe

c----------------------------------------
c Begin with the zeroth-level
c discretization (20 elements)
c
c Nodes are set manually on the unit sphere
c----------------------------------------

c---
c the icosahedron has 12 vertices
c---

      ru = 0.25D0*dsqrt(10.0D0+2.0D0*dsqrt(5.0D0))
      rm = 0.25D0*(1.0D0+dsqrt(5.0D0))   ! twice the golden ratio

      c0 = 0.0    D0
c     c1 = 0.9512 D0
      c1 = ru
c     c2 = 0.8507 D0
      c2 = 2.0*ru/dsqrt(5.0D0)
c     c3 = 0.8090 D0
      c3 = rm
c     c4 = 0.4253 D0
      c4 = ru/dsqrt(5.0D0)
c     c5 = 0.2629 D0
      c5 = dsqrt( ru**2-c3**2-c4**2)
      c6 = 0.5D0
c     c7 = 0.6882 D0
      c7 = dsqrt( ru**2-c4**2-c6**2)

c     write (6,*) c5,c7

      VX(1) =  c0
      VY(1) =  c0
      VZ(1) =  c1

      VX(2) =  c0
      VY(2) =  c2
      VZ(2) =  c4

      VX(3) =  c3
      VY(3) =  c5
      VZ(3) =  c4

      VX(4) =  c6
      VY(4) = -c7
      VZ(4) =  c4

      VX(5) = -c6
      VY(5) = -c7
      VZ(5) =  c4

      VX(6) = -c3
      VY(6) =  c5
      VZ(6) =  c4

      VX(7) = -c6
      VY(7) =  c7
      VZ(7) = -c4

      VX(8) =  c6
      VY(8) =  c7
      VZ(8) = -c4

      VX(9) =  c3
      VY(9) = -c5
      VZ(9) = -c4

      VX(10) =  c0
      VY(10) = -c2
      VZ(10) = -c4

      VX(11) = -c3
      VY(11) = -c5
      VZ(11) = -c4

      VX(12) =  c0
      VY(12) =  c0
      VZ(12) = -c1

c     Do i=1,12
c      VX(i) = VX(1)/c1
c      VY(i) = VX(1)/c1
c      VZ(i) = VX(1)/c1
c     End Do

c------------------------
c define the corner nodes
c------------------------

      x(1,1) = VX(1)   ! first element
      y(1,1) = VY(1)
      z(1,1) = VZ(1)
      x(1,2) = VX(3)
      y(1,2) = VY(3)
      z(1,2) = VZ(3)
      x(1,3) = VX(2)
      y(1,3) = VY(2)
      z(1,3) = VZ(2)
c---
      x(2,1) = VX(1)
      y(2,1) = VY(1)
      z(2,1) = VZ(1)
      x(2,2) = VX(4)
      y(2,2) = VY(4)
      z(2,2) = VZ(4)
      x(2,3) = VX(3)
      y(2,3) = VY(3)
      z(2,3) = VZ(3)
c---
      x(3,1) = VX(1)
      y(3,1) = VY(1)
      z(3,1) = VZ(1)
      x(3,2) = VX(5)
      y(3,2) = VY(5)
      z(3,2) = VZ(5)
      x(3,3) = VX(4)
      y(3,3) = VY(4)
      z(3,3) = VZ(4)
c---
      x(4,1) = VX(1)
      y(4,1) = VY(1)
      z(4,1) = VZ(1)
      x(4,2) = VX(6)
      y(4,2) = VY(6)
      z(4,2) = VZ(6)
      x(4,3) = VX(5)
      y(4,3) = VY(5)
      z(4,3) = VZ(5)
c---
      x(5,1) = VX(1)
      y(5,1) = VY(1)
      z(5,1) = VZ(1)
      x(5,2) = VX(2)
      y(5,2) = VY(2)
      z(5,2) = VZ(2)
      x(5,3) = VX(6)
      y(5,3) = VY(6)
      z(5,3) = VZ(6)
c---
      x(6,1) = VX(2)
      y(6,1) = VY(2)
      z(6,1) = VZ(2)
      x(6,2) = VX(3)
      y(6,2) = VY(3)
      z(6,2) = VZ(3)
      x(6,3) = VX(8)
      y(6,3) = VY(8)
      z(6,3) = VZ(8)
c---
      x(7,1) = VX(3)
      y(7,1) = VY(3)
      z(7,1) = VZ(3)
      x(7,2) = VX(4)
      y(7,2) = VY(4)
      z(7,2) = VZ(4)
      x(7,3) = VX(9)
      y(7,3) = VY(9)
      z(7,3) = VZ(9)
c---
      x(8,1) = VX(4)
      y(8,1) = VY(4)
      z(8,1) = VZ(4)
      x(8,2) = VX(5)
      y(8,2) = VY(5)
      z(8,2) = VZ(5)
      x(8,3) = VX(10)
      y(8,3) = VY(10)
      z(8,3) = VZ(10)
c---
      x(9,1) = VX(5)
      y(9,1) = VY(5)
      z(9,1) = VZ(5)
      x(9,2) = VX(6)
      y(9,2) = VY(6)
      z(9,2) = VZ(6)
      x(9,3) = VX(11)
      y(9,3) = VY(11)
      z(9,3) = VZ(11)
c---
      x(10,1) = VX(6)
      y(10,1) = VY(6)
      z(10,1) = VZ(6)
      x(10,2) = VX(2)
      y(10,2) = VY(2)
      z(10,2) = VZ(2)
      x(10,3) = VX(7)
      y(10,3) = VY(7)
      z(10,3) = VZ(7)
c---
      x(11,1) = VX(2)
      y(11,1) = VY(2)
      z(11,1) = VZ(2)
      x(11,2) = VX(8)
      y(11,2) = VY(8)
      z(11,2) = VZ(8)
      x(11,3) = VX(7)
      y(11,3) = VY(7)
      z(11,3) = VZ(7)
c---
      x(12,1) = VX(3)
      y(12,1) = VY(3)
      z(12,1) = VZ(3)
      x(12,2) = VX(9)
      y(12,2) = VY(9)
      z(12,2) = VZ(9)
      x(12,3) = VX(8)
      y(12,3) = VY(8)
      z(12,3) = VZ(8)
c---
      x(13,1) = VX(4)
      y(13,1) = VY(4)
      z(13,1) = VZ(4)
      x(13,2) = VX(10)
      y(13,2) = VY(10)
      z(13,2) = VZ(10)
      x(13,3) = VX(9)
      y(13,3) = VY(9)
      z(13,3) = VZ(9)
c---
      x(14,1) = VX(5)
      y(14,1) = VY(5)
      z(14,1) = VZ(5)
      x(14,2) = VX(11)
      y(14,2) = VY(11)
      z(14,2) = VZ(11)
      x(14,3) = VX(10)
      y(14,3) = VY(10)
      z(14,3) = VZ(10)
c---
      x(15,1) = VX(6)
      y(15,1) = VY(6)
      z(15,1) = VZ(6)
      x(15,2) = VX(7)
      y(15,2) = VY(7)
      z(15,2) = VZ(7)
      x(15,3) = VX(11)
      y(15,3) = VY(11)
      z(15,3) = VZ(11)
c---
      x(16,1) = VX(7)
      y(16,1) = VY(7)
      z(16,1) = VZ(7)
      x(16,2) = VX(8)
      y(16,2) = VY(8)
      z(16,2) = VZ(8)
      x(16,3) = VX(12)
      y(16,3) = VY(12)
      z(16,3) = VZ(12)
c---
      x(17,1) = VX(8)
      y(17,1) = VY(8)
      z(17,1) = VZ(8)
      x(17,2) = VX(9)
      y(17,2) = VY(9)
      z(17,2) = VZ(9)
      x(17,3) = VX(12)
      y(17,3) = VY(12)
      z(17,3) = VZ(12)
c--- 
      x(18,1) = VX(9)
      y(18,1) = VY(9)
      z(18,1) = VZ(9)
      x(18,2) = VX(10)
      y(18,2) = VY(10)
      z(18,2) = VZ(10)
      x(18,3) = VX(12)
      y(18,3) = VY(12)
      z(18,3) = VZ(12)
c---
      x(19,1) = VX(10)
      y(19,1) = VY(10)
      z(19,1) = VZ(10)
      x(19,2) = VX(11)
      y(19,2) = VY(11)
      z(19,2) = VZ(11)
      x(19,3) = VX(12)
      y(19,3) = VY(12)
      z(19,3) = VZ(12)
c---
      x(20,1) = VX(11)
      y(20,1) = VY(11)
      z(20,1) = VZ(11)
      x(20,2) = VX(7)
      y(20,2) = VY(7)
      z(20,2) = VZ(7)
      x(20,3) = VX(12)
      y(20,3) = VY(12)
      z(20,3) = VZ(12)

c------------------------------------------
c compute the mid-points of the three sides
c of the 20 first-generation elements
c
c midpoints are numbered 4, 5, 6
c------------------------------------------

      Nelm = 20

      Do i=1,Nelm

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
c project the nodes onto the unit sphere
c---

       Do k=1,Nelm
        Do l=1,6
         rad = dsqrt(x(k,l)**2+y(k,l)**2+z(k,l)**2) 
         x(k,l) = x(k,l)/rad
         y(k,l) = y(k,l)/rad
         z(k,l) = z(k,l)/rad
       End Do
      End Do

      if(Ndiv.eq.0) Go to 98    ! icosahedron done

c-------------------------------------------
c compute the local element node coordinates
c for discretization levels 1 through Ndiv
c-------------------------------------------

      Do i=1,Ndiv

       nm = 0      ! count the new elements arising by sub-division
                   ! four element will be generated during each pass
       Do j=1,Nelm ! over old elements

c---
c assign corner points to sub-elements
c these will become the "new" elements
c---

        nm = nm+1

        xn(nm,1) = x(j,1)                  !  first sub-element
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

        xn(nm,1) = x(j,4)                !  second sub-element
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

        xn(nm,1) = x(j,6)                !  third sub-element
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

        xn(nm,1) = x(j,4)                !  fourth sub-element
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

       End Do                      !  end of old-element loop

c--------------------------------------
c number of elements has been increased
c by a factor of four
c--------------------------------------

       Nelm = 4*Nelm

c---
c relabel the new points
c and place them in the master list
c---

       Do k=1,Nelm
        Do l=1,6

         x(k,l) = xn(k,l)
         y(k,l) = yn(k,l)
         z(k,l) = zn(k,l)

c--- project onto the unit sphere

         rad = dsqrt(x(k,l)**2+y(k,l)**2+z(k,l)**2)
         x(k,l) = x(k,l)/rad
         y(k,l) = y(k,l)/rad
         z(k,l) = z(k,l)/rad

         xn(k,l) = 0.0D0   ! zero just in case
         yn(k,l) = 0.0D0
         zn(k,l) = 0.0D0

        End Do
       End Do

c----------

      End Do            !  end of discretization- level loop

c-----------------------------------------

 98   Continue

c-----------------------------------
c Generate a list of global nodes by looping 
c over all elementsand adding nodes not found
c in the current list.
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
      n(1,3) = 3  ! third  node of first element is global node 3
      n(1,4) = 4
      n(1,5) = 5
      n(1,6) = 6  ! sixth  node of first element is global node 6

      Npts = 6

c---
c loop over further elements
c
c Iflag=0 will signal a new global node
c---

      Do i=2,Nelm        ! loop over elements
       Do j=1,6          ! loop over element nodes

        Iflag=0

         Do k=1,Npts
          If(abs(x(i,j)-p(k,1)).le.eps) then
           If(abs(y(i,j)-p(k,2)).le.eps) then
            If(abs(z(i,j)-p(k,3)).le.eps) then

             Iflag = 1    ! the node has been recorded previously
             n(i,j) = k   ! the jth local node of element i
                          ! is the kth global node 
            End If
           End If
          End If
         End Do
        
         If(Iflag.eq.0) then     ! record the node

          Npts = Npts+1          ! one more global node

          p(Npts,1) = x(i,j)
          p(Npts,2) = y(i,j)
          p(Npts,3) = z(i,j)

          n(i,j) = Npts   ! the jth local node of element i
                          ! is the new global node 

         End If

       End Do
      End Do                      !  end of loop over elements

c----------------------------------
c Generate connectivity table: ne(i,j)
c for elements touching global node i
c
c  ne(i,j) ... ne(i,1) is the number of elements touching 
c                      the ith global node
c
c  ne(i,2:ne(i,1)+1) are the corresponding element labels
c----------------------------------

c---
c initialize
c---

      Do i=1,Npts
       Do j=1,7
        ne(i,j) = 0
       End Do
      End Do 

c---
c loop over global nodes
c---

      Do i=1,Npts 

       ne(i,1) = 0

       Icount = 1

       Do j=1,Nelm     ! loop over elements
        Do k=1,6       ! loop over element nodes

         If(abs(p(i,1)-x(j,k)).le.eps) then
          If(abs(p(i,2)-y(j,k)).le.eps) then
           If(abs(p(i,3)-z(j,k)).le.eps) then

            Icount = Icount+1
            ne(i,1) = ne(i,1)+1
            ne(i,Icount) = j

           End If
          End If
         End If

        End Do
       End Do
     
      End Do        !  end of loop over global nodes

c------------------------------------------
c Generate connectivity table nbe(i,j)
c
c  nbe(i,j) .. label of element sharing side j of element i
c              where j = 1, 2, 3
c------------------------------------------

c---
c initialize
c---

      Do i=1,Nelm
       Do j=1,3
        nbe(i,j) = 0
       End Do
      End Do

c---
c loop over elements
c---

      Do i=1,Nelm             !  loop over elements
       jcount=1
       Do j=4,6               !  loop over mid-points

        Do k=1,Nelm           !  test element
         If(k.eq.i) Go to 91  !  self-element
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

      Do i=1,Npts

       rad = dsqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)

       p(i,1) = p(i,1)/rad
       p(i,2) = p(i,2)/rad
       p(i,3) = p(i,3)/rad

      End Do

c---------
c printing
c---------

c     write (6,*)
c     write (6,*) Nelm,' grid elements'
c     write (6,*) Npts,' grid points'
c     write (6,*)

c-----
c done
c-----

      return
      end
