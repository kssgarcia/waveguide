      subroutine prd_3d_pr
     +
     +  (NSG
     +  ,RLX,RLY,RLZ
     +  ,ICH1,phimax
     +  ,ICH2,dstmax
     +  ,ICH3,dstmin
     +  ,P1,P2,P3,P4,P5
     +  ,Istop
     +  )

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c============================================

c--------------------------------------------------
c Adaptive represenation of a periodic 3D line.  
c
c Quadratic interpolation for properties P1, ..., P5
c
c Criteria for re-distribution include:
c
c TEST #1:  If the angle between three consecutive points is too
c           large, then an additional interior point is added.
c
c TEST #2:  If the arc length spanning two adjacent points is
c           too large, then an additional interior point is added.
c
c TEST #3:  If the arc length spanning two adjacent points is
c           too small, then the two points are
c           replaced by a single
c           point, unless the resulting distribution violates either
c           of the first two tests.
c
c SYMBOLS:
c -------
c
c xc, yc,zc :  centers of the arcs
c rad       :  radii of the arcs 
c chi1, chi3:  signed angles, as defined in the text
c              by Pozrikidis (1998)
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(0:501),y(0:501),z(0:501)
      Dimension P1(0:501),P2(0:501),P3(0:501),P4(0:501),P5(0:501)

      Dimension xc(500),yc(500),zc(500),r(500),s(500)
      Dimension chi1(500),chi3(500)

      Parameter (Nmax = 500)

      common/xxyyzz/x,y,z
      common/ARCC/xc,yc,zc,r,s,chi1,chi3

      common/SV33c/a11,a12,a13,a21,a22,a23,a31,a32,a33

c---
c constants
c---

      pi = 3.14159 265358 97932 384
      pih = 0.5*pi

      Nsafe = 100*NSG

c---
c prepare to run
c---

      Istop = 0
      Ipass = 0         !  counts passes

  98  Continue

      Ipass = Ipass+1

c---
c safety valves
c---

      If(NSG.gt.Nmax) then

        write (6,*) 
        write (6,*) 'More than',Nmax,'points'
        write (6,*) 'I will stop'
        write (6,*) '-----------------------'
        Istop = 1
        Return

      End if

      If(Ipass.gt.Nsafe) then

        write (6,*) 
        write(6,*) 'stuck in regridding loop'
        write(6,*) '.....exiting!'
        Istop = 1
        Return

      End If

      NSG1 = NSG+1
      NSG2 = NSG+2
      NSG3 = NSG+3

c---
c extend
c---

       x(0) = x (NSG)-RLX
       y(0) = y (NSG)-RLY
       z(0) = z (NSG)-RLZ
      p1(0) = p1(NSG)
      p2(0) = p2(NSG)
      p3(0) = p3(NSG)
      p4(0) = p4(NSG)
      p5(0) = p5(NSG)

       x(NSG1) = x (1)+RLX
       y(NSG1) = y (1)+RLY
       z(NSG1) = z (1)+RLZ
      p1(NSG1) = p1(1)
      p2(NSG1) = p2(1)
      p3(NSG1) = p3(1)
      p4(NSG1) = p4(1)
      p5(NSG1) = p5(1)

       x(NSG2) = x (2)+RLX
       y(NSG2) = y (2)+RLY
       z(NSG2) = z (2)+RLZ
      p1(NSG2) = p1(2)
      p2(NSG2) = p2(2)
      p3(NSG2) = p3(2)
      p4(NSG2) = p4(2)
      p5(NSG2) = p5(2)

       x(NSG3) = x (3)+RLX
       y(NSG3) = y (3)+RLY
       z(NSG3) = z (3)+RLZ
      p1(NSG3) = p1(3)
      p2(NSG3) = p2(3)
      p3(NSG3) = p3(3)
      p4(NSG3) = p4(3)
      p5(NSG3) = p5(3)

c-----------------
c Define the arcs
c-----------------

      Iloop = 1

      Do 3 i=Iloop,NSG1    

      ib = i-2
      ia = i-1
      i1 = i+1
      i2 = i+2
      i3 = i+3
c---
c  examine the arc passing through points i-1, i, i+1
c---

      call arc_3d
     +            (x(ia),x(i),x(i1)
     +            ,y(ia),y(i),y(i1)
     +            ,z(ia),z(i),z(i1)
     +            ,xc(i),yc(i),zc(i)
     +            ,r(i)
     +            ,chi1(i),chi3(i)
     +            )

      If(I.eq.NSG1) Go to 3

c---
c compute magnitude of subtended angles
c---

      phi12 =-chi1(i)
      phi23 = chi3(i)
      phi13 = phi12+phi23

c---
c prepare to test
c---

      x0 = xc(i)
      y0 = yc(i)
      z0 = zc(i)
      r0 = r (i)
      rs = r0**2

c---
c  compute arc lengths subtended by points 1-2 and 2-3
c---

      dst12 = r0*phi12
      dst23 = r0*phi23

      dst13 = dst12+dst23

c------
c  compute constants for quadratic interpolation of 
c  properties Pj
c
c  formula is: Pj = aj * phi**2 + bj *phi + cj
c
c  where phi is the angle measured from point i-1
c---
c
c      par   = phi13*phi12**2 - phi12*phi13**2
c 
c      rhs12 = P1(i) -P1(ia)
c      rhs13 = P1(i1)-P1(ia)
c      aa1   = (rhs12*phi13   -rhs13*phi12   )/par
c      bb1   = (rhs13*phi12**2-rhs12*phi13**2)/par
c      cc1   = P1(ia)
c
c      rhs12 = P2(i) -P2(ia)
c      rhs13 = P2(i1)-P2(ia)
c      aa2   = (rhs12*phi13   -rhs13*phi12   )/par
c      bb2   = (rhs13*phi12**2-rhs12*phi13**2)/par
c      cc2   = P2(ia)
c
c      rhs12 = P3(i) -P3(ia)
c      rhs13 = P3(i1)-P3(ia)
c      aa3   = (rhs12*phi13   -rhs13*phi12   )/par
c      bb3   = (rhs13*phi12**2-rhs12*phi13**2)/par
c      cc3   = P3(ia)
c
c      rhs12 = P4(i) -P4(ia)
c      rhs13 = P4(i1)-P4(ia)
c      aa4   = (rhs12*phi13   -rhs13*phi12   )/par
c      bb4   = (rhs13*phi12**2-rhs12*phi13**2)/par
c      cc4   = P4(ia)
c
c      rhs12 = P5(i) -P5(ia)
c      rhs13 = P5(i1)-P5(ia)
c      aa5   = (rhs12*phi13   -rhs13*phi12   )/par
c      bb5   = (rhs13*phi12**2-rhs12*phi13**2)/par
c      cc5   = P5(ia)
c
c------------------------------------------
c TEST # 1
c If angle is too large, add one point
c------------------------------------------

      If(Ich1.eq.1) then

        If(phi13.gt.pih   ) Go to 96
        If(phi13.gt.phimax) Go to 97

      End If

c----------------------------------------------------
c TEST #2 .... distance between points i and i+1
c---------------------------------------------------

      If(ich2.eq.1) then         

       If(dst23.gt.dstmax) Go to 92

      End If

c-------------------------------------------------------------------	
c   TEST #3 .... size of arc between points i and i+1
c
c  if spacing is too small, will remove points i and i+1 and
c  replace them with a single point
c  provided that check1 and check 2 pass
c
c-------------------------------------------------------------------	

      If(Ich3.eq.0) Go to 3

      If(dst23.gt.dstmin) Go to 3

c---
c compute the location of the point midway between i and i+1 
c---

      phi = 0.5*phi23
      csp = cos(phi)

      tx = (y(i1)-y0)*(z(i )-z0)-(y(i )-y0)*(z(i1)-z0)
      ty = (x(i1)-x0)*(z(i1)-z0)-(x(i1)-x0)*(z(i )-z0)
      tz = (x(i )-x0)*(y(i )-y0)-(x(i )-x0)*(y(i1)-y0)

      a11 = x(i1)-x0
      a12 = y(i1)-y0
      a13 = z(i1)-z0

      a21 = x(i)-x0
      a22 = y(i)-y0
      a23 = z(i)-z0

      a31 = tx
      a32 = ty
      a33 = tz

      b1 = rs*csp+x0*(x(i1)-x0)+y0*(y(i1)-y0)+z0*(z(i1)-z0)
      b2 = rs*csp+x0*(x(i )-x0)+y0*(y(i )-y0)+z0*(z(i )-z0)
      b3 = x0*tx+y0*ty+z0*tz

      call cramer_33 
     +
     +    (a11,a12,a13
     +    ,a21,a22,a23
     +    ,a31,a32,a33
     +    ,b1,b2,b3
     +    ,xt,yt,zt
     +    )

c---
c  arc passing through
c  points i-2, i-1, tmp
c---

      If(i.ne.1) then
        xxx = x(ib)
        yyy = y(ib)
        zzz = z(ib)
      Else
        xxx = x(NSG-1)-RLX
        yyy = y(NSG-1)-RLY
        zzz = z(NSG-1)-RLZ
      End If

      call arc_3d
     +            (xxx,x(ia),xt
     +            ,yyy,y(ia),yt
     +            ,zzz,z(ia),zt
     +            ,xx,yy,zz
     +            ,actis
     +            ,gn1,gn3
     +            )

      totang0 = gn3-gn1

c---
c  arc passing through
c  points i-1, xt, i+2
c---

      call arc_3d
     +            (x(ia),xt,x(i2)
     +            ,y(ia),yt,y(i2)
     +            ,z(ia),zt,z(i2)
     +            ,xx,yy,zz
     +            ,actis
     +            ,gn1,gn3
     +            )

      pt12 =-gn1
      pt23 = gn3

      totang1 = gn3-gn1

      at12 = actis*pt12
      at23 = actis*pt23

c---
c  arc passing through
c  points i-1, xt, i+2
c---

      If(i.ne.nsg) then
        xxx = x(i3)
        yyy = y(i3)
        zzz = z(i3)
      Else
        xxx = x(3)+RLX
        yyy = y(3)+RLY
        zzz = z(3)+RLZ
      End If

      call arc_3d
     +            (xt,x(i2),xxx
     +            ,yt,y(i2),yyy
     +            ,zt,z(i2),zzz
     +            ,xx,yy,zz
     +            ,actis
     +            ,gn1,gn3
     +            )

      totang2 = gn3-gn1

c------
c If the resulting angle/arc distribution satisfies angle 
c and arc length requirements, 
c then replace points 'i' and 'i+1' with the
c temporary point.
c------

      if 
     +      (totang0.lt.phimax
     +  .and.totang1.lt.phimax
     +  .and.totang2.lt.phimax
     +     .and.at12.lt.dstmax
     +     .and.at23.lt.dstmax
     +      ) 
     +     Go to 93

  3    Continue

c-------------------------------
c    We now have an acceptable
c    point distribution
c------------------------------

      Go to 29

c--------------------
c  There is a problem:
c  execution will stop
c--------------------

  96  Continue

      write (6,200) I,phi13,phimax
      write (6,201)

      Istop = 1
      Return
c----

  97  Continue

c-----------------------------------
c  If an arc is too large,
c  remove the mid-point and add
c  two points along the arc
c-----------------------------------

      write (6,*) " Arc ",i," is too large; will insert one point"

c---
c will remove point 'i' and add two points evenly-spaced
c with respect to angle along the arc 
c---

      phit1 = phi13/3.0
      phit2 = 2.0*phit1

      csp1 = cos(phit1)
      csp2 = cos(phit2)

      tx=(y(ia)-y0)*(z(i1)-z0)-(y(i1)-y0)*(z(ia)-z0)
      ty=(x(i1)-x0)*(z(ia)-z0)-(x(ia)-x0)*(z(i1)-z0)
      tz=(x(ia)-x0)*(y(i1)-y0)-(x(i1)-x0)*(y(ia)-y0)

      a11 = x(ia)-x0
      a12 = y(ia)-y0
      a13 = z(ia)-z0

      a21 = x(i1)-x0
      a22 = y(i1)-y0
      a23 = z(i1)-z0

      a31=tx
      a32=ty
      a33=tz

      b1=csp1*rs+x0*(x(ia)-x0)+y0*(y(ia)-y0)+z0*(z(ia)-z0)
      b2=csp2*rs+x0*(x(i1)-x0)+y0*(y(i1)-y0)+z0*(z(i1)-z0)
      b3=x0*tx+y0*ty+z0*tz

      call cramer_33                   ! first new point
     +
     +    (a11,a12,a13
     +    ,a21,a22,a23
     +    ,a31,a32,a33
     +    ,b1,b2,b3
     +    ,xt1,yt1,zt1
     +    )

c--- Interpolation by algebra
c
c     P1t1 = aa1*phit1**2 + bb1*phit1 + cc1
c     P2t1 = aa2*phit1**2 + bb2*phit1 + cc2
c     P3t1 = aa3*phit1**2 + bb3*phit1 + cc3
c     P4t1 = aa4*phit1**2 + bb4*phit1 + cc4
c     P5t1 = aa5*phit1**2 + bb5*phit1 + cc5
c
c     write (6,101) P1t1,P2t1,P3t1,P4t1,P5t1
c
c--- Lagrange Interpolation

      fc1 = (phit1-phi12)*(phit1-phi13)/(0.000-phi12)/(0.000-phi13)
      fc2 = (phit1-0.000)*(phit1-phi13)/(phi12-0.000)/(phi12-phi13)
      fc3 = (phit1-0.000)*(phit1-phi12)/(phi13-0.000)/(phi13-phi12)

      P1t1 = fc1*P1(ia) + fc2*P1(i) + fc3*P1(i1)  
      P2t1 = fc1*P2(ia) + fc2*P2(i) + fc3*P2(i1)  
      P3t1 = fc1*P3(ia) + fc2*P3(i) + fc3*P3(i1)  
      P4t1 = fc1*P4(ia) + fc2*P4(i) + fc3*P4(i1)  
      P5t1 = fc1*P5(ia) + fc2*P5(i) + fc3*P5(i1)  

c     write (6,101) P1t1,P2t1,P3t1,P4t1,P5t1
c---
c  second new point
c---

      b1=csp2*rs+x0*(x(ia)-x0)+y0*(y(ia)-y0)+z0*(z(ia)-z0)
      b2=csp1*rs+x0*(x(i1)-x0)+y0*(y(i1)-y0)+z0*(z(i1)-z0)

      call cramer_33                  ! second new point
     +
     +   (a11,a12,a13
     +   ,a21,a22,a23
     +   ,a31,a32,a33
     +   ,b1,b2,b3
     +   ,xt2,yt2,zt2
     +   )

c--- Interpolation by algebra
c
c     P1t2 = aa1*phit2**2 + bb1*phit2 + cc1
c     P2t2 = aa2*phit2**2 + bb2*phit2 + cc2
c     P3t2 = aa3*phit2**2 + bb3*phit2 + cc3
c     P4t2 = aa4*phit2**2 + bb4*phit2 + cc4
c     P5t2 = aa5*phit2**2 + bb5*phit2 + cc5
c
c     write (6,101) P1t2,P2t2,P3t2,P4t2,P5t2
c
c--- Lagrange Interpolation

      fc1 = (phit2-phi12)*(phit2-phi13)/(0.000-phi12)/(0.000-phi13)
      fc2 = (phit2-0.000)*(phit2-phi13)/(phi12-0.000)/(phi12-phi13)
      fc3 = (phit2-0.000)*(phit2-phi12)/(phi13-0.000)/(phi13-phi12)

      P1t2 = fc1*P1(ia) + fc2*P1(i) + fc3*P1(i1)  
      P2t2 = fc1*P2(ia) + fc2*P2(i) + fc3*P2(i1)  
      P3t2 = fc1*P3(ia) + fc2*P3(i) + fc3*P3(i1)  
      P4t2 = fc1*P4(ia) + fc2*P4(i) + fc3*P4(i1)  
      P5t2 = fc1*P5(ia) + fc2*P5(i) + fc3*P5(i1)  

c---
c  redistribution
c---

      NSG = NSG+1

      Do j=NSG,i2,-1           !  renumber grid points
        ja   =j-1
        x (j)=x (ja)
        y (j)=y (ja)
        z (j)=z (ja)
	P1(j)=P1(ja)         
	P2(j)=P2(ja)         
	P3(j)=P3(ja)         
	P4(j)=P4(ja)         
	P5(j)=P5(ja)         
      End Do
	  
c---
c  set location and property P1 of points i and i+1
c---

      x (i) = xt1
      y (i) = yt1
      z (i) = zt1
      p1(i) = p1t1
      p2(i) = p2t1
      p3(i) = p3t1
      p4(i) = p4t1
      p5(i) = p5t1
	  
      x (i1) = xt2
      y (i1) = yt2
      z (i1) = zt2
      p1(i1) = p1t2
      p2(i1) = p2t2
      p3(i1) = p3t2
      p4(i1) = p4t2
      p5(i1) = p5t2
	  
      Go to 98

c---------------------------
c If a segment is too long,
c add a point in the middle
c---------------------------

  92  Continue

      write (6,*) " Distance",i," is too large; will insert one point"

c---
c  compute location of new point
c  by blending
c---

c---
c forward arc
c---

      phi = 0.5*phi23
      csp = cos(phi)

      tx=(y(i1)-y0)*(z(i )-z0)-(y(i )-y0)*(z(i1)-z0)
      ty=(x(i )-x0)*(z(i1)-z0)-(x(i1)-x0)*(z(i )-z0)
      tz=(x(i1)-x0)*(y(i )-y0)-(x(i )-x0)*(y(i1)-y0)

      a11=x(i1)-x0
      a12=y(i1)-y0
      a13=z(i1)-z0

      a21 = x(i)-x0
      a22 = y(i)-y0
      a23 = z(i)-z0

      a31 = tx
      a32 =ty
      a33 = tz

      b1 = rs*csp+x0*(x(i1)-x0)+y0*(y(i1)-y0)+z0*(z(i1)-z0)
      b2 = rs*csp+x0*(x(i )-x0)+y0*(y(i )-y0)+z0*(z(i )-z0)
      b3 = x0*tx+y0*ty+z0*tz

      call cramer_33
     +
     +    (a11,a12,a13
     +    ,a21,a22,a23
     +    ,a31,a32,a33
     +    ,b1,b2,b3
     +    ,xt1,yt1,zt1
     +    )


      phit = phi12+phi

c--- Interpolation by algebra
c
c     p1t1 = aa1 * phit**2 + bb1 * phit + cc1
c     p2t1 = aa2 * phit**2 + bb2 * phit + cc2
c     p3t1 = aa3 * phit**2 + bb3 * phit + cc3
c     p4t1 = aa4 * phit**2 + bb4 * phit + cc4
c     p5t1 = aa5 * phit**2 + bb5 * phit + cc5
c
c     write (6,101) P1t1,P2t1,P3t1,P4t1,P5t1
c
c--- Lagrange interpolation

      fc1 = (phit-phi12)*(phit-phi13)/(0.000-phi12)/(0.000-phi13)
      fc2 = (phit-0.000)*(phit-phi13)/(phi12-0.000)/(phi12-phi13)
      fc3 = (phit-0.000)*(phit-phi12)/(phi13-0.000)/(phi13-phi12)

      P1t1 = fc1*P1(ia) + fc2*P1(i) + fc3*P1(i1)  
      P2t1 = fc1*P2(ia) + fc2*P2(i) + fc3*P2(i1)  
      P3t1 = fc1*P3(ia) + fc2*P3(i) + fc3*P3(i1)  
      P4t1 = fc1*P4(ia) + fc2*P4(i) + fc3*P4(i1)  
      P5t1 = fc1*P5(ia) + fc2*P5(i) + fc3*P5(i1)  

c     write (6,101) P1t1,P2t1,P3t1,P4t1,P5t1

c---
c backward arc
c---
      call arc_3d
     +            (x(i),x(i1),x(i2)
     +            ,y(i),y(i1),y(i2)
     +            ,z(i),z(i1),z(i2)
     +            ,xc(i1),yc(i1),zc(i1)
     +            ,r(i1)
     +            ,chi1(i1),chi3(i1)
     +            )

      phi23 = chi3(i1)
      phi12 =-chi1(i1)
      phi13 = phi12+phi23

      x0 = xc(i1)
      y0 = yc(i1)
      z0 = zc(i1)
      r0 = r (i1)
      rs = r0**2

      phi = 0.5*phi12
      csp = cos(phi)

      tx=(y(i2)-y0)*(z(i1)-z0)-(y(i1)-y0)*(z(i2)-z0)
      ty=(x(i1)-x0)*(z(i2)-z0)-(x(i2)-x0)*(z(i1)-z0)
      tz=(x(i2)-x0)*(y(i1)-y0)-(x(i1)-x0)*(y(i2)-y0)

      a11=x(i )-x0
      a12=y(i )-y0
      a13=z(i )-z0

      a21=x(i1)-x0
      a22=y(i1)-y0
      a23=z(i1)-z0

      a31=tx
      a32=ty
      a33=tz

      b1 =rs*csp+x0*(x(i )-x0)+y0*(y(i )-y0)+z0*(z(i )-z0)
      b2 =rs*csp+x0*(x(i1)-x0)+y0*(y(i1)-y0)+z0*(z(i1)-z0)
      b3 =x0*tx+y0*ty+z0*tz

      call cramer_33
     +
     +   (a11,a12,a13
     +   ,a21,a22,a23
     +   ,a31,a32,a33
     +   ,b1,b2,b3
     +   ,xt2,yt2,zt2
     +   )

c--- Lagrange interpolation

      fc1 = (phi-phi12)*(phi-phi13)/(0.000-phi12)/(0.000-phi13)
      fc2 = (phi-0.000)*(phi-phi13)/(phi12-0.000)/(phi12-phi13)
      fc3 = (phi-0.000)*(phi-phi12)/(phi13-0.000)/(phi13-phi12)

      P1t2 = fc1*P1(i) + fc2*P1(i1) + fc3*P1(i2)  
      P2t2 = fc1*P2(i) + fc2*P2(i1) + fc3*P2(i2)  
      P3t2 = fc1*P3(i) + fc2*P3(i1) + fc3*P3(i2)  
      P4t2 = fc1*P4(i) + fc2*P4(i1) + fc3*P4(i2)  
      P5t2 = fc1*P5(i) + fc2*P5(i1) + fc3*P5(i2)  

c---
c blending
c---

      xt = 0.5*(xt1+xt2)
      yt = 0.5*(yt1+yt2)
      zt = 0.5*(zt1+zt2)

      P1t = 0.5*(P1t1+P1t2)
      P2t = 0.5*(P2t1+P2t2)
      P3t = 0.5*(P3t1+P3t2)
      P4t = 0.5*(P4t1+P4t2)
      P5t = 0.5*(P5t1+P5t2)

c---
c redistribute
c---

      NSG=NSG+1

      Do j=NSG,i1,-1
        ja    = j-1
        x (j) = x (ja)
        y (j) = y (ja)
        z (j) = z (ja)
        p1(j) = p1(ja)
        p2(j) = p2(ja)
        p3(j) = p3(ja)
        p4(j) = p4(ja)
        p5(j) = p5(ja)
      End do

      x (i1)=xt
      y (i1)=yt
      z (i1)=zt
      p1(i1)=p1t
      p2(i1)=p2t
      p3(i1)=p3t
      p4(i1)=p4t
      p5(i1)=p5t

      Go to 98

c-------------------------------
c  If a segment is too short,
c  remove the end-points and
c  add a point in the middle
c-------------------------------

  93  Continue

      write (6,*) " Distance ",i," is too small; will remove one point"

      phit = phi12+phi

c--- Interpolation by algebra
c
c     p1t = aa1 * phit**2 + bb1 * phit + cc1
c     p2t = aa2 * phit**2 + bb2 * phit + cc2
c     p3t = aa3 * phit**2 + bb3 * phit + cc3
c     p4t = aa4 * phit**2 + bb4 * phit + cc4
c     p5t = aa5 * phit**2 + bb5 * phit + cc5
c
c     write (6,101) P1t,P2t,P3t,P4t,P5t
c
c--- Lagrange interpolation

      fc1 = (phit-phi12)*(phit-phi13)/(0.000-phi12)/(0.000-phi13)
      fc2 = (phit-0.000)*(phit-phi13)/(phi12-0.000)/(phi12-phi13)
      fc3 = (phit-0.000)*(phit-phi12)/(phi13-0.000)/(phi13-phi12)

      P1t = fc1*P1(ia) + fc2*P1(i) + fc3*P1(i1)  
      P2t = fc1*P2(ia) + fc2*P2(i) + fc3*P2(i1)  
      P3t = fc1*P3(ia) + fc2*P3(i) + fc3*P3(i1)  
      P4t = fc1*P4(ia) + fc2*P4(i) + fc3*P4(i1)  
      P5t = fc1*P5(ia) + fc2*P5(i) + fc3*P5(i1)  

      If(i.eq.NSG) then   !  special case 
        x (1) = xt
        y (1) = yt
        z (1) = zt
        p1(1) = p1t
        p2(1) = p2t
        p3(1) = p3t
        p4(1) = p4t
        p5(1) = p5t
      Else                 !  all other points
        x (i)=xt
        y (i)=yt
        z (i)=zt
        p1(i)=p1t
        p2(i)=p2t
        p3(i)=p3t
        p4(i)=p4t
        p5(i)=p5t

        Do j=i1,NSG
         j1 = j+1
         x (j)=x (j1)
         y (j)=y (j1)
         z (j)=z (j1)
         p1(j)=p1(j1)
         p2(j)=p2(j1)
         p3(j)=p3(j1)
         p4(j)=p4(j1)
         p5(j)=p5(j1)
        End Do
      End if

      NSG  = NSG-1
      NSG1 = NSG+1

c---
c to account for removal when i = NSG(old)
c---

      If(i.eq.NSG1) then
       x(1) = x(nsg1)-RLX
       y(1) = y(nsg1)-RLY
       z(1) = z(nsg1)-RLZ
      End If

      Go to 98

c-------------------------------------------

  29  Continue

c-----------------------------
c THE SUCCESSFUL DISTRIBUTION
c-----------------------------

c-----------------------
c Compute the arc length
c-----------------------

      S(1) = 0.
      Do i = 2,NSG
       ia = i-1
       S(i) = S(ia)+0.5*(R(i) *abs(chi1(i))
     +                  +R(ia)*abs(chi3(i)))
      End Do

      i = NSG1
      S(i) = S(NSG)+0.5*(R(i)*abs(chi1(i))
     +                  +R(1)*abs(chi3(1)))
      S(NSG2) = S(NSG1)+S(2)

c-----
c Done
c-----

  101 Format (10(1x,f10.5))
  200 Format (' ARC',I3,' IS TOO BIG',2(F10.5))
  201 Format (' EXECUTION WILL BE INTERRUPTED')

      return
      end
