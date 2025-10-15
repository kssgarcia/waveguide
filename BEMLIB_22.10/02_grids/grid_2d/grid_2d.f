      program grid_2d

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------
c Division of a planar contour consisting of 
c line segments and circular arcs in boundary
c elements with corresponding shapes
c
c SYMBOLS:
c --------
c
c Xeg(i,j) is the x coordinate of the ith point
c          of the jth segment
c Yeg(i,j) is the y coordinate of the ith point
c          of the jth segment
c
c Isym = 0 if the element distribution
c          on a segment is not symmetric
c Isym = 1 if the element distribution on a segment 
c          is symmetric with respect to the mid-point
c
c If the jth segment is an arc, then
c        Teg(i,j) is the polar angle of the ith point subtended
c        from the arc center
c
c seg(i,j) is the cummulative arc length  at the ith point  of the
c          jth segment
c
c CAPACITY:
c --------
c
c  20 segments
c 128 elements per segment
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xeg(129,20),Yeg(129,20),Teg(129,20),seg(129,20)
      Dimension Xmg(129,20),Ymg(129,20),Tmg(129,20),smg(129,200)

      Dimension Xe(129),Ye(129),Te(129),se(129)
      Dimension Xm(128),Ym(128),Tm(128),sm(128)

c----------
c constants
c----------

      pi = 3.14159 265358 D0

c--------
c prepare
c--------
 
      factor = 1.00D0

      Itry = 1

  93  Continue

      open (4,file="grid_2d.dat")
      read (4,*) xstart,ystart     ! starting point
      read (4,*) 

      open (1,file="PLOTDAT")

      sinit = 0.0D0      ! origin of arc length
      Iseg = 1           ! counts the number of segments

c-------------------
c begin the gridding
c-------------------

  97  Continue           ! target for another segment

      read (4,*) Itype   ! 0 to quit, 1 for lines, 2 for arcs

      If(Itype.eq.0) Go to 98

c------------------------
      If(Itype.eq.1) then        ! line segment
c------------------------

        read (4,*) xend,yend,N,ratio,Isym

        ratio = ratio*factor

        x1 = xstart
        y1 = ystart
        x2 = xend
        y2 = yend

        call elm_line
     +
     +   (N,ratio
     +   ,X1,Y1
     +   ,X2,Y2
     +   ,sinit
     +   ,Isym
     +   ,Xe,Ye,se
     +   ,Xm,Ym,sm
     +   )

c------------------------------
       Else If(Itype.eq.2) then    ! circular arc
c------------------------------

        read (4,*) radius,angle1,angle2,N,ratio,Isym

        ratio = ratio*factor

        angle1 = angle1*pi
        angle2 = angle2*pi

        xcnt = xstart-radius*Dcos(angle1)
        ycnt = ystart-radius*Dsin(angle1)

        call elm_arc
     +
     +    (N,ratio
     +    ,Xcnt,Ycnt
     +    ,radius
     +    ,angle1,angle2
     +    ,sinit
     +    ,Isym
     +    ,Xe,Ye,Te,se
     +    ,Xm,Ym,Tm,sm
     +    )

c------------
       End If
c------------

c------------------------
c Save into global arrays
c------------------------

      Do i=1,N+1
        xeg(i,Iseg) = xe(i)
        yeg(i,Iseg) = ye(i)
        seg(i,Iseg) = se(i)
        xmg(i,Iseg) = xm(i)
        ymg(i,Iseg) = ym(i)
        smg(i,Iseg) = sm(i)
        tmg(i,Iseg) = tm(i)
        teg(i,Iseg) = te(i)
      End Do

      xstart = xe(N+1)    ! reset starting point
      ystart = ye(N+1)
      sinit  = se(N+1)

c---------
c printing
c---------

      write (1,101) N+1,Itype,Isym
c     write (6,101) N+1,Itype,Isym

      Do i=1,N+1
       write (1,100) i,xeg(i,Iseg),yeg(i,Iseg)
     +                ,teg(i,iseg),seg(i,Iseg)
c      write (6,100) i,xeg(i,Iseg),yeg(i,Iseg)
c    +                ,teg(i,iseg),seg(i,Iseg)
      End Do

      Iseg = Iseg+1

c-------------
      Go to 97      ! return for another segment
c-------------

  98  Continue

      write (1,101) Itype
      close (1)
      close (4)

  99  Continue

c-----
c done
c-----

 100  Format (1x,i3,10(1x,f10.5))
 101  Format (10(1x,i3))

      Stop
      End
