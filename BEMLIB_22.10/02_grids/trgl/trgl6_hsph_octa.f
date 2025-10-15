      subroutine trgl6_hsph_octa
     +
     +   (ndiv
     +   ,npts,nelm
     +   )

c==========================================
c BEMLIB, FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------------
c Triangulation of the hemisphere
c resting on the xz plane
c
c The triangulation is effected
c by subdividing half a regular octahedron
c into six-node quadratic triangles.
c
c SYMBOLS:
c -------
c
c  ndiv .... level of discretization of octahedron
c            nvid = 0 gives 8 elements
c
c  npts .... number of nodes
c  nelm .... number of surface elements
c
c  x(i,j), y(i,j), z(i,j) .... Cartesian coordinates of local node j
c                              on element i
c                              j = 1,...,6
c                              i = 1,...,nelm
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
c              ne(i,2:ne(i,1)) are the corresponding element labels
c
c  nbe(i,j) .. label of element sharing side j of element i
c              where j = 1, 2, 3
c--------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  x(512,6), y(512,6), z(512,6)
      Dimension xn(512,6),yn(512,6),zn(512,6)
      Dimension p(1026,3)

      Dimension n(512,6),ne(1026,7),nbe(512,3)

      Parameter (eps=0.00000001)   ! used for point checks

      common/points/p,ne
      common/elmnts/n,nbe

c----------------------------------------
c Begin with the
c zeroth level discretization (4 elements)
c
c Nodes are set manually on the unit hemisphere
c in the upper half space
c----------------------------------------

      nelm = 4

c---
c vertex points
c---

      x(1,1) = 0.0D0  ! first element
      y(1,1) = 0.0D0
      z(1,1) = 1.0D0

      x(1,2) = 1.0D0
      y(1,2) = 0.0D0
      z(1,2) = 0.0D0
 
      x(1,3) = 0.0D0
      y(1,3) = 1.0D0
      z(1,3) = 0.0D0

c---

      x(2,3) =-1.0D0   ! second element
      y(2,3) = 0.0D0
      z(2,3) = 0.0D0

      x(2,1) = 0.0D0
      y(2,1) = 0.0D0
      z(2,1) = 1.0D0

      x(2,2) = 0.0D0
      y(2,2) = 1.0D0
      z(2,2) = 0.0D0

c---

      x(3,3) = 1.0D0   ! third element
      y(3,3) = 0.0D0
      z(3,3) = 0.0D0

      x(3,1) = 0.0D0
      y(3,1) = 0.0D0
      z(3,1) =-1.0D0

      x(3,2) = 0.0D0
      y(3,2) = 1.0D0
      z(3,2) = 0.0D0

c---

      x(4,1) = 0.0D0  ! fourth element
      y(4,1) = 0.0D0
      z(4,1) =-1.0D0

      x(4,2) =-1.0D0
      y(4,2) = 0.0D0
      z(4,2) = 0.0D0

      x(4,3) = 0.0D0
      y(4,3) = 1.0D0
      z(4,3) = 0.0D0

c------------------------------------------
c Compute the mid-points of the three edges
c of the four first-generation elements
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

c-----------------------------
c project onto the unit sphere
c-----------------------------

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

      Do i=1,ndiv   ! run over discretization levels

       num = 0        ! Counts the new elements arising by sub-dividing
                      ! Four element will be generated during each pass
       Do j=1,nelm    ! over old elements

c---
c assign corner points to sub-elements
c these will become the "new" elements
c---

        num = num+1

        xn(num,1) = x(j,1)                  !  first sub-element
        yn(num,1) = y(j,1)
        zn(num,1) = z(j,1)

        xn(num,2) = x(j,4)
        yn(num,2) = y(j,4) 
        zn(num,2) = z(j,4)

        xn(num,3) = x(j,6)
        yn(num,3) = y(j,6)
        zn(num,3) = z(j,6)

        xn(num,4) = 0.5D0*(xn(num,1)+xn(num,2))
        yn(num,4) = 0.5D0*(yn(num,1)+yn(num,2))
        zn(num,4) = 0.5D0*(zn(num,1)+zn(num,2))

        xn(num,5) = 0.5D0*(xn(num,2)+xn(num,3))
        yn(num,5) = 0.5D0*(yn(num,2)+yn(num,3))
        zn(num,5) = 0.5D0*(zn(num,2)+zn(num,3))

        xn(num,6) = 0.5D0*(xn(num,3)+xn(num,1))
        yn(num,6) = 0.5D0*(yn(num,3)+yn(num,1))
        zn(num,6) = 0.5D0*(zn(num,3)+zn(num,1))

        num = num+1

        xn(num,1) = x(j,4)                !  second sub-element
        yn(num,1) = y(j,4)
        zn(num,1) = z(j,4)

        xn(num,2) = x(j,2)
        yn(num,2) = y(j,2)
        zn(num,2) = z(j,2)

        xn(num,3) = x(j,5)
        yn(num,3) = y(j,5)
        zn(num,3) = z(j,5)

        xn(num,4) = 0.5D0*(xn(num,1)+xn(num,2))
        yn(num,4) = 0.5D0*(yn(num,1)+yn(num,2))
        zn(num,4) = 0.5D0*(zn(num,1)+zn(num,2))

        xn(num,5) = 0.5D0*(xn(num,2)+xn(num,3))
        yn(num,5) = 0.5D0*(yn(num,2)+yn(num,3))
        zn(num,5) = 0.5D0*(zn(num,2)+zn(num,3))

        xn(num,6) = 0.5D0*(xn(num,3)+xn(num,1))
        yn(num,6) = 0.5D0*(yn(num,3)+yn(num,1))
        zn(num,6) = 0.5D0*(zn(num,3)+zn(num,1))

        num = num+1

        xn(num,1) = x(j,6)                !  third sub-element
        yn(num,1) = y(j,6)
        zn(num,1) = z(j,6)

        xn(num,2)= x(j,5)
        yn(num,2)= y(j,5)
        zn(num,2)= z(j,5)

        xn(num,3) = x(j,3)
        yn(num,3) = y(j,3)
        zn(num,3) = z(j,3)

        xn(num,4) = 0.5D0*(xn(num,1)+xn(num,2))
        yn(num,4) = 0.5D0*(yn(num,1)+yn(num,2))
        zn(num,4) = 0.5D0*(zn(num,1)+zn(num,2))

        xn(num,5) = 0.5D0*(xn(num,2)+xn(num,3))
        yn(num,5) = 0.5D0*(yn(num,2)+yn(num,3))
        zn(num,5) = 0.5D0*(zn(num,2)+zn(num,3))

        xn(num,6) = 0.5D0*(xn(num,3)+xn(num,1))
        yn(num,6) = 0.5D0*(yn(num,3)+yn(num,1))
        zn(num,6) = 0.5D0*(zn(num,3)+zn(num,1))

        num = num+1

        xn(num,1) = x(j,4)                !  fourth sub-element
        yn(num,1) = y(j,4)
        zn(num,1) = z(j,4)

        xn(num,2) = x(j,5)
        yn(num,2) = y(j,5)
        zn(num,2) = z(j,5)

        xn(num,3) = x(j,6)
        yn(num,3) = y(j,6)
        zn(num,3) = z(j,6)

        xn(num,4) = 0.5D0*(xn(num,1)+xn(num,2))
        yn(num,4) = 0.5D0*(yn(num,1)+yn(num,2))
        zn(num,4) = 0.5D0*(zn(num,1)+zn(num,2))

        xn(num,5) = 0.5D0*(xn(num,2)+xn(num,3))
        yn(num,5) = 0.5D0*(yn(num,2)+yn(num,3))
        zn(num,5) = 0.5D0*(zn(num,2)+zn(num,3))

        xn(num,6) = 0.5D0*(xn(num,3)+xn(num,1))
        yn(num,6) = 0.5D0*(yn(num,3)+yn(num,1))
        zn(num,6) = 0.5D0*(zn(num,3)+zn(num,1))

       End Do                      !  end of old-element loop

c--------------------------------------
c number of elements has been increased
c by a factor of four
c--------------------------------------

       nelm = 4*nelm

c---
c relabel the new points
c and place them in the master list
c---

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

         xn(k,l) = 0.0D0     ! zero just in case
         yn(k,l) = 0.0D0
         zn(k,l) = 0.0D0

        End Do
       End Do

c----------

      End Do           !  end of discretization-level loop

  98  Continue

c------------------------------
c Generate a list of global nodes 
c by looping over all elements
c and adding nodes not found in the list. 
c
c Fill in the connectivity table n(i,j)
c containing node numbers of element points 1-6
c-----------------------------------

c---
c six nodes of the first element
c are entered mannualy
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
      n(1,2) = 2  ! second node of first element is global node 1
      n(1,3) = 3
      n(1,4) = 4
      n(1,5) = 5
      n(1,6) = 6  ! sixth node of first element is global node 6

      npts = 6

c---
c loop over further elements
c
c Iflag=0 will signal a new global node
c---

      Do i=2,nelm          !  loop over elements
       Do j=1,6            ! loop over element nodes

        iflag=0

         Do k=1,npts
          If(abs(x(i,j)-p(k,1)).lt.eps) then
           If(abs(y(i,j)-p(k,2)).lt.eps) then
            If(abs(z(i,j)-p(k,3)).lt.eps) then

             Iflag = 1  ! this node has been previously recorded            
             n(i,j)= k  ! the jth local node of element i
                        ! it the kth global node
            End If
           End If
          End If
         End Do
        
        If(iflag.eq.0) then       ! record the node

         npts = npts+1            ! one more global node

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
c  ne(i,1) is the number of elements touching
c          the global node i
c  ne(i,2:ne(i,1)) are the corresponding element labels
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

      Do i=1,npts                 !  loop over global nodes

       ne(i,1) = 0

       icount = 1

       Do j=1,nelm        ! loop over elements
        Do k=1,6          ! loop over element nodes

         If(abs(p(i,1)-x(j,k)).lt.eps) then
          If(abs(p(i,2)-y(j,k)).lt.eps) then
           If(abs(p(i,3)-z(j,k)).lt.eps) then

            Icount = Icount+1
            ne(i,1) = ne(i,1)+1
            ne(i,icount) = j

           End If
          End If
         End If

        End Do
       End Do
     
      End Do      !   end of loop over global nodes

c------------------------------------------
c Generate connectivity table nbe(i,j) 
c
c  nbe(i,j) is the label of the element
c           sharing side j of element i,
c           where j = 1, 2, 3
c
c For boundary elements with only 2 neighbors,
c this array entry will be zero.
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

      Do i=1,nelm        ! loop over elements
       jcount=1
       Do j=4,6          ! loop over mid-points

        Do k=1,nelm           ! test element
         If(k.eq.i) Go to 91  ! self-element
         Do l=4,6             ! loop over mid-points

          If(abs(x(i,j)-x(k,l)).lt.eps) then
           If(abs(y(i,j)-y(k,l)).lt.eps) then
            If(abs(z(i,j)-z(k,l)).lt.eps) then
             nbe(i,jcount)=k
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

c---------------------------------
c Project global nodes p(i,j) onto
c the unit hemi-sphere 
c---------------------------------

      Do i=1,npts

       rr=1.0D0/sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)

       p(i,1) = p(i,1)*rr
       p(i,2) = p(i,2)*rr
       p(i,3) = p(i,3)*rr

      End Do

c---------
c printing
c---------

c     write (6,*)
c     write (6,*) nelm,' grid elements'
c     write (6,*) npts,' grid points'
c     write (6,*)

c-----
c done
c-----

      return
      end
