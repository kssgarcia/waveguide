      subroutine trgl6_sc
     +
     +  (a,ndiv,npts,nelm
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c==========================================

c--------------------------------------------
c Triangulation of a square with unit sides
c containing a circular hole of radius "a"
c in the xy plane
c
c The patch is confined in: -0.5 < x,y <0.5
c
c Will subdivide a 12-element pattern
c into six-node triangles
c
c SYMBOLS:
c -------
c
c  ndiv .... level of discretization
c            nvid = 0 gives 12 elements
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
c              x = p(i,1)
c              y = p(i,2)
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
c
c     Iedge(i,1) = 0 if global node i is an interior node
c
c     Iedge(i,1) = 1 if global node i is an edge node
c     Iedge(i,2): global index of image node
c
c     Iedge(i,3) = 3 if global node i is a corner node
c     Iedge(i,2): global index of first image node
c     Iedge(i,3): global index of first image node
c     Iedge(i,4): global index of first image node
c--------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  x(512,6), y(512,6)
      Dimension xn(512,6),yn(512,6)

      Dimension   n(512,6),efl(512,6),efln(512,6)
      Dimension nbe(512,3)

      Dimension  p(1090,3), gfl(1090,1)
      Dimension ne(1090,9)
   
      Dimension Iedge(1090,4)

      Parameter (eps=0.000001)

      common/points/p,ne
      common/elmnts/n,nbe
      common/elmntf/efl
      common/edgepoints/Iedge

c----------------------------------------
c  zeroth level discretization (8 elements)
c  nodes are set manually
c----------------------------------------

      side = 0.5
      side2 = 2.0*side

      nelm = 12

        x(1,1) = a;  ! first element
        y(1,1) =-a;
        x(1,2) = side;
        y(1,2) =-side;
        x(1,3) = side;
        y(1,3) = 0.0;
      efl(1,1) = 2;
      efl(1,2) = 1;
      efl(1,3) = 1;
      efl(1,4) = 0;
      efl(1,5) = 1;
      efl(1,6) = 0;

        x(2,1) = a;  ! second element
        y(2,1) =-a;
        x(2,2) = side;
        y(2,2) = 0.0;
        x(2,3) = a;
        y(2,3) = a;
      efl(2,1) = 2;
      efl(2,2) = 1;
      efl(2,3) = 2;
      efl(2,4) = 0;
      efl(2,5) = 0;
      efl(2,6) = 2;

        x(3,1) =   a;  ! third element
        y(3,1) =   a;
        x(3,2) = side;
        y(3,2) = 0.0;
        x(3,3) = side;
        y(3,3) = side;
      efl(3,1) = 2;
      efl(3,2) = 1;
      efl(3,3) = 1;
      efl(3,4) = 0;
      efl(3,5) = 1;
      efl(3,6) = 0;

        x(4,1) =   a; ! fourth element
        y(4,1) =   a;
        x(4,2) = side;
        y(4,2) = side;
        x(4,3) = 0.0;
        y(4,3) = side;
      efl(4,1) = 2;
      efl(4,2) = 1;
      efl(4,3) = 1;
      efl(4,4) = 0;
      efl(4,5) = 1;
      efl(4,6) = 0;

        x(5,1) = 0.0; ! fifth element
        y(5,1) = side;
        x(5,2) =  -a;
        y(5,2) =   a;
        x(5,3) =   a;
        y(5,3) =   a;
      efl(5,1)=1;  
      efl(5,2)=2;
      efl(5,3)=2;
      efl(5,4)=0;
      efl(5,5)=2;
      efl(5,6)=0;

        x(6,1) =  -a;  ! sixth element
        y(6,1) =   a;
        x(6,2) = 0.0;
        y(6,2) = side;
        x(6,3) =-side;
        y(6,3) = side;
      efl(6,1)=2;  
      efl(6,2)=1;
      efl(6,3)=1;
      efl(6,4)=0;
      efl(6,5)=1;
      efl(6,6)=0;

c-----
c mid-nodes numbered 4,5,6
c-----

      Do i=1,6
        x(i,4) = 0.5*(x(i,1)+x(i,2))
        y(i,4) = 0.5*(y(i,1)+y(i,2))
        x(i,5) = 0.5*(x(i,2)+x(i,3))
        y(i,5) = 0.5*(y(i,2)+y(i,3))
        x(i,6) = 0.5*(x(i,3)+x(i,1))
        y(i,6) = 0.5*(y(i,3)+y(i,1))
      End Do

c-----
c rest of elements by reflection
c-----

      Do i=1,6
        Do j=1,6
          x(6+i,j)= -x(i,j)
          y(6+i,j)= -y(i,j)
        efl(6+i,j)=efl(i,j)
        End Do
      End Do

c-------
c project innermost nodes
c on a circle of radius a
c-------

       Do i=1,12
        Do j=1,6
         if(efl(i,j).eq.2) then
          rad = a/sqrt(x(i,j)**2+y(i,j)**2)
          x(i,j) = x(i,j)*rad
          y(i,j) = y(i,j)*rad
         end if
        End Do
       End Do

c----------------------------------------
c compute node coordinates on each element 
c for discretization levels 1 through ndiv
c----------------------------------------

      if(ndiv.eq.0) Go to 98   ! skip if ndiv = 0

      Do i=1,ndiv

       num = 0       ! counts the new elements arising by sub-division

       Do j=1,nelm   !  loop over current elements

c---
c assign corner points to sub-elements
c these are to become "new" elements
c---

        num = num+1                 !  first sub-element

        xn(num,1) = x(j,1) 
        yn(num,1) = y(j,1)

        xn(num,2) = x(j,4)
        yn(num,2) = y(j,4) 

        xn(num,3) = x(j,6)
        yn(num,3) = y(j,6)

        efln(num,1) = efl(j,1)
        efln(num,2) = efl(j,4)
        efln(num,3) = efl(j,6)

        xn(num,4)= 0.5*(xn(num,1)+xn(num,2))
        yn(num,4)= 0.5*(yn(num,1)+yn(num,2))

        xn(num,5)= 0.5*(xn(num,2)+xn(num,3))
        yn(num,5)= 0.5*(yn(num,2)+yn(num,3))

        xn(num,6)= 0.5*(xn(num,3)+xn(num,1))
        yn(num,6)= 0.5*(yn(num,3)+yn(num,1))

        efln(num,4) = 0
        if(efln(num,1).eq.1.and.efln(num,2).eq.1) efln(num,4) = 1
        if(efln(num,1).eq.2.and.efln(num,2).eq.2) efln(num,4) = 2
        efln(num,5) = 0
        if(efln(num,2).eq.1.and.efln(num,3).eq.1) efln(num,5) = 1
        if(efln(num,2).eq.2.and.efln(num,3).eq.2) efln(num,5) = 2
        efln(num,6) = 0
        if(efln(num,3).eq.1.and.efln(num,1).eq.1) efln(num,6) = 1
        if(efln(num,3).eq.2.and.efln(num,1).eq.2) efln(num,6) = 2

        num = num+1                !  second sub-element

        xn(num,1)= x(j,4)
        yn(num,1)= y(j,4)

        xn(num,2)= x(j,2)
        yn(num,2)= y(j,2)

        xn(num,3)= x(j,5)
        yn(num,3)= y(j,5)

        efln(num,1)=efl(j,4)
        efln(num,2)=efl(j,2)
        efln(num,3)=efl(j,5)

        xn(num,4)= 0.5*(xn(num,1)+xn(num,2))
        yn(num,4)= 0.5*(yn(num,1)+yn(num,2))

        xn(num,5)= 0.5*(xn(num,2)+xn(num,3))
        yn(num,5)= 0.5*(yn(num,2)+yn(num,3))

        xn(num,6)= 0.5*(xn(num,3)+xn(num,1))
        yn(num,6)= 0.5*(yn(num,3)+yn(num,1))

        efln(num,4) = 0
        if(efln(num,1).eq.1.and.efln(num,2).eq.1) efln(num,4) = 1
        if(efln(num,1).eq.2.and.efln(num,2).eq.2) efln(num,4) = 2
        efln(num,5) = 0
        if(efln(num,2).eq.1.and.efln(num,3).eq.1) efln(num,5) = 1
        if(efln(num,2).eq.2.and.efln(num,3).eq.2) efln(num,5) = 2
        efln(num,6) = 0
        if(efln(num,3).eq.1.and.efln(num,1).eq.1) efln(num,6) = 1
        if(efln(num,3).eq.2.and.efln(num,1).eq.2) efln(num,6) = 2

        num = num+1                !  third sub-element

        xn(num,1)= x(j,6)
        yn(num,1)= y(j,6)

        xn(num,2)= x(j,5)
        yn(num,2)= y(j,5)

        xn(num,3)= x(j,3)
        yn(num,3)= y(j,3)

        efln(num,1)=efl(j,6)
        efln(num,2)=efl(j,5)
        efln(num,3)=efl(j,3)

        xn(num,4)= 0.5*(xn(num,1)+xn(num,2))
        yn(num,4)= 0.5*(yn(num,1)+yn(num,2))

        xn(num,5)= 0.5*(xn(num,2)+xn(num,3))
        yn(num,5)= 0.5*(yn(num,2)+yn(num,3))

        xn(num,6)= 0.5*(xn(num,3)+xn(num,1))
        yn(num,6)= 0.5*(yn(num,3)+yn(num,1))

        efln(num,4) = 0
        if(efln(num,1).eq.1.and.efln(num,2).eq.1) efln(num,4) = 1
        if(efln(num,1).eq.2.and.efln(num,2).eq.2) efln(num,4) = 2
        efln(num,5) = 0
        if(efln(num,2).eq.1.and.efln(num,3).eq.1) efln(num,5) = 1
        if(efln(num,2).eq.2.and.efln(num,3).eq.2) efln(num,5) = 2
        efln(num,6) = 0
        if(efln(num,3).eq.1.and.efln(num,1).eq.1) efln(num,6) = 1
        if(efln(num,3).eq.2.and.efln(num,1).eq.2) efln(num,6) = 2

        num = num+1               !  fourth sub-element

        xn(num,1)= x(j,4) 
        yn(num,1)= y(j,4)

        xn(num,2)= x(j,5)
        yn(num,2)= y(j,5)

        xn(num,3)= x(j,6)
        yn(num,3)= y(j,6)

        efln(num,1)=efl(j,4)
        efln(num,2)=efl(j,5)
        efln(num,3)=efl(j,6)

        xn(num,4)= 0.5*(xn(num,1)+xn(num,2))    ! mid points
        yn(num,4)= 0.5*(yn(num,1)+yn(num,2))

        xn(num,5)= 0.5*(xn(num,2)+xn(num,3))
        yn(num,5)= 0.5*(yn(num,2)+yn(num,3))

        xn(num,6)= 0.5*(xn(num,3)+xn(num,1))
        yn(num,6)= 0.5*(yn(num,3)+yn(num,1))

        efln(num,4) = 0
        if(efln(num,1).eq.1.and.efln(num,2).eq.1) efln(num,4) = 1
        if(efln(num,1).eq.2.and.efln(num,2).eq.2) efln(num,4) = 2
        efln(num,5) = 0
        if(efln(num,2).eq.1.and.efln(num,3).eq.1) efln(num,5) = 1
        if(efln(num,2).eq.2.and.efln(num,3).eq.2) efln(num,5) = 2
        efln(num,6) = 0
        if(efln(num,3).eq.1.and.efln(num,1).eq.1) efln(num,6) = 1
        if(efln(num,3).eq.2.and.efln(num,1).eq.2) efln(num,6) = 2

       End Do                      !  end of old element loop

c---------------------------------
c number of elements was increased
c by a factor of four
c---------------------------------

       nelm=4*nelm

c---
c relabel the new points
c and place them in the master list
c---

       Do k=1,nelm
        Do l=1,6

           x(k,l) =  xn(k,l)
           y(k,l) =  yn(k,l)
         efl(k,l) = efln(k,l)

         !  project the innermost nodes

         if(efl(k,l).eq.2) then
          rad = a/sqrt(x(k,l)**2+y(k,l)**2)
          x(k,l) = x(k,l)*rad
          y(k,l) = y(k,l)*rad
         end if

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
      gfl(1,1) = efl(1,1)

      p(2,1) = x(1,2)
      p(2,2) = y(1,2)
      gfl(2,1) = efl(1,2)

      p(3,1) = x(1,3)
      p(3,2) = y(1,3)
      gfl(3,1) = efl(1,3)

      p(4,1) = x(1,4)
      p(4,2) = y(1,4)
      gfl(4,1) = efl(1,4)

      p(5,1) = x(1,5)
      p(5,2) = y(1,5)
      gfl(5,1) = efl(1,5)

      p(6,1) = x(1,6)
      p(6,2) = y(1,6)
      gfl(6,1) = efl(1,6)

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
          If(Dabs(x(i,j)-p(k,1)).le.eps) then
           If(Dabs(y(i,j)-p(k,2)).le.eps) then

             Iflag = 1    ! the node has been recorded previously
             n(i,j) = k   ! the jth local node of element i
                          ! is the kth global node 
           End If
          End If
         End Do
        
         If(Iflag.eq.0) then     ! record the node

          npts = npts+1          ! one more global node

            p(npts,1) =   x(i,j)
            p(npts,2) =   y(i,j)
          gfl(npts,1) = efl(i,j);

          n(i,j) = npts  ! the jth local node of element i
                         ! is the new global mode
         End If

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
         If(k.eq.i) Go to 91  ! skip self-element
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
c Iedge(i,1) = 2 for inner circle nodes
c--------------------------------

      Do i=1,npts
       Iedge(i,1) = 0
       if(abs(p(i,1)-side).lt.eps) Iedge(i,1) = 1
       if(abs(p(i,2)-side).lt.eps) Iedge(i,1) = 1
       if(abs(p(i,1)+side).lt.eps) Iedge(i,1) = 1
       if(abs(p(i,2)+side).lt.eps) Iedge(i,1) = 1
       dist = sqrt(p(i,1)**2+p(i,2)**2)
       if(abs(dist-a).lt.eps) Iedge(i,1) = 2  ! inner circle
      End Do

c--------------------------------
c loop over the edge nodes and set:
c
c Iedge(i,1) = 1 for edge nodes
c Iedge(i,2): global index of image node
c
c Iedge(i,1) = 3 for corner nodes
c Iedge(i,2): global index of first image node
c Iedge(i,3): global index of first image node
c Iedge(i,4): global index of first image node
c--------------------------------

      Do 77 i=1,npts

       if(Iedge(i,1).eq.1) then

       Ic = 1      ! counter of image nodes

       Do 78 j=1,npts

       if(i.eq.j) Go to 78
       if(Iedge(j,1).eq.0) Go to 78

       test = abs(p(i,1)-p(j,1)+side2)+abs(p(i,2)-p(j,2))

       if(test.lt.eps) then
        Ic=Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,1)-p(j,1)-side2)+abs(p(i,2)-p(j,2))

       if(test.lt.eps) then
        Ic=Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)+side2)+abs(p(i,1)-p(j,1))

       if(test.lt.eps) then
        Ic=Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)-side2)+abs(p(i,1)-p(j,1))

       if(test.lt.eps) then
        Ic=Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)-side2)+abs(p(i,1)-p(j,1)-side2)

       if(test.lt.eps) then
        Ic=Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)-side2)+abs(p(i,1)-p(j,1)+side2)

       if(test.lt.eps) then
        Ic=Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)+side2)+abs(p(i,1)-p(j,1)-side2)

       if(test.lt.eps) then
        Ic=Ic+1
        Iedge(i,Ic) = j
       end if

       test = abs(p(i,2)-p(j,2)+side2)+abs(p(i,1)-p(j,1)+side2)

       if(test.lt.eps) then
        Ic=Ic+1
        Iedge(i,Ic) = j
       end if

  78   Continue

       Iedge(i,1) = Ic-1

       end if

  77  Continue

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
