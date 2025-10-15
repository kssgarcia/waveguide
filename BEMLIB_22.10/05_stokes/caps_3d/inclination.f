      subroutine inclination
     +
     +  (nxy,jxy
     +  ,cx,cy
     +  ,thmax,thmin
     +  )
 
c======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c======================================

c----------------------------------------
c  Compute the drop inclination angle
c
c  Diagonalizes the inertia tensor
c  of the cross-section of the interface
c  in the xy plane to find the directions of
c  max and min elongation
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension p(1026,3),ne(1026,7)
      Dimension x(100),y(100)
      Dimension jxy(100)

      common/points/ p,ne

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c---
c read points in the xy plane
c---

      Do i=1,nxy
        x(i) = p(jxy(i),1)-cx
        y(i) = p(jxy(i),2)-cy
      End Do

c---
c  compute line integrals around perimeter
c  of the trace on the xy plane.
c---

      axx = 0.0D0
      ayy = 0.0D0
      axy = 0.0D0

      Do i=1,nxy-1

        dl = sqrt((x(i+1)-x(i))**2+(y(i+1)-y(i))**2)  

        ayy = ayy+ 0.5*(x(i)**2+x(i+1)**2) * dl
        axx = axx+ 0.5*(y(i)**2+y(i+1)**2) * dl
        axy = axy+ 0.5*(x(i)*y(i)+x(i+1)*y(i+1)) * dl

      End Do

c---
c  find eigenvalues of inertia tensor
c---

      t1 = axx+ayy
      t2 = axx*ayy-axy*axy

      amax = 0.5D0*(t1+sqrt(t1**2-4.0D0*t2))
      amin = 0.5D0*(t1-sqrt(t1**2-4.0D0*t2))

c---------------------------------------------
c  find eigenvectors of inertia tensor 
c  and orientation of the
c  principal directions relative to the x axis
c---------------------------------------------

      v1 = 1.0D0
      v2 = -(axx-amin)/axy

      vnorm = sqrt(v1**2+v2**2)

      v1 = v1/vnorm
      v2 = v2/vnorm

c     thmax = acos(v1)*180.0D0/pi
      thmax = acos(v1)/pi

      v1=-1.0D0
      v2=-(axx-amax)/axy

      vnorm = sqrt(v1**2+v2**2)

      v1 = v1/vnorm
      v2 = v2/vnorm

c     thmin = acos(v1)*180.0/pi
      thmin = acos(v1)/pi

c-----
c done
c-----

      return
      end
