      subroutine taylor 
     +
     +  (npts,nxy,jxy
     +  ,cx,cy,cz
     +  ,rmax,rmin,zmax
     +  ,dxy,thmax,thmin
     +  )

c===========================================
c FDLIB, BEBLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c---------------------------------------------
c compute the Taylor deformation parameter 
c and inclination of the major and minor axes
c
c cx,cy,cz are the coordinates of the centroid
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  p(1026,3)
      Dimension ne(1026,7)
      Dimension  n(512,6),nbe(512,3)
      Dimension jxy(100)

      Dimension x(0:100),y(0:100),r(0:100)

c-------
c common
c-------

      common/points/p,ne
      common/elmnts/n,nbe

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c-----------------------------------------------
c put the points around the xy trace in a vector
c-----------------------------------------------

      Do i=1,nxy
        nid  = jxy(i)
        x(i) = p(nid,1)-cx
        y(i) = p(nid,2)-cy
        r(i) = sqrt(x(i)**2+y(i)**2)
      End Do

      x(0) = x(nxy-1)
      y(0) = y(nxy-1)
      r(0) = r(nxy-1)

c     Do i=0,nxy
c      write (6,*) i,x(i),y(i),r(i)
c     End Do

c---
c find points with
c maximum and minimum radial distance
c from the centroid
c---

      rmax  =    0.0  
      rmin  = 1000.0 
      imax  = 0
      imin  = 0

      Do i=1,nxy-1

         if(r(i).gt.rmax) then
            rmax  = r(i)
            imax  = i
         end if

         if(r(i).lt.rmin) then
            rmin = r(i)
            imin  = i
         end if

      End Do

c---
c find the coordinates of the point with the
c maximum distance by quadratic interpolation
c with respect to the poly-line
c---

      j0 = imax-1
      j1 = imax
      j2 = imax+1

      xx0 = -sqrt((x(j0)-x(j1))**2+(y(j0)-y(j1))**2)
      xx1 = 0.0
      xx2 =  sqrt((x(j2)-x(j1))**2+(y(j2)-y(j1))**2)

      ff0 = r(j0)
      ff1 = r(j1)
      ff2 = r(j2)
      tm2 = (ff2-ff1)/(xx2-xx1)
      tm1 = (ff0-ff1)/(xx0-xx1)
      bbb = (tm2-tm1)/(xx2-xx0)
      ccc =  tm2-bbb*(xx2-xx1)

      xxmx = -0.50*ccc   /bbb + xx1
      rmax = -0.25*ccc**2/bbb + ff1
     
      xxd  = xxmx-xx1

      ff0  = x(j0)
      ff1  = x(j1)
      ff2  = x(j2)
      tm2  = (ff2-ff1)/(xx2-xx1)
      tm1  = (ff0-ff1)/(xx0-xx1)
      bbb  = (tm2-tm1)/(xx2-xx0)
      ccc  =  tm2-bbb*(xx2-xx1)
      xmax = (bbb*xxd+ccc)*xxd+ff1

      ff0  = y(j0)
      ff1  = y(j1)
      ff2  = y(j2)
      tm2  = (ff2-ff1)/(xx2-xx1)
      tm1  = (ff0-ff1)/(xx0-xx1)
      bbb  = (tm2-tm1)/(xx2-xx0)
      ccc  =  tm2-bbb*(xx2-xx1)
      ymax = (bbb*xxd+ccc)*xxd+ff1

      rmax  = sqrt(xmax**2+ymax**2)

      thmax = 360.0*acos(xmax/rmax)/pi2

      if(ymax.lt.0)  thmax = -thmax
      if(thmax.lt.0) thmax = thmax + 180

c---
c find the coordinates of the point with the
c minimum distance by quadratic interpolation
c with respect to the poly-line
c---

      j0 = imin-1
      j1 = imin
      j2 = imin+1

      xx0 = -sqrt((x(j0)-x(j1))**2+(y(j0)-y(j1))**2)
      xx1 = 0.0
      xx2 =  sqrt((x(j2)-x(j1))**2+(y(j2)-y(j1))**2)

      ff0 = r(j0)
      ff1 = r(j1)
      ff2 = r(j2)
      tm2 = (ff2-ff1)/(xx2-xx1)
      tm1 = (ff0-ff1)/(xx0-xx1)
      bbb = (tm2-tm1)/(xx2-xx0)
      ccc =  tm2-bbb*(xx2-xx1)

      xxmn = -0.50*ccc   /bbb + xx1
      rmin = -0.25*ccc**2/bbb + ff1

      xxd  = xxmn-xx1

      ff0  = x(j0)
      ff1  = x(j1)
      ff2  = x(j2)
      tm2  = (ff2-ff1)/(xx2-xx1)
      tm1  = (ff0-ff1)/(xx0-xx1)
      bbb  = (tm2-tm1)/(xx2-xx0)
      ccc  =  tm2-bbb*(xx2-xx1)
      xmin = (bbb*xxd+ccc)*xxd+ff1

      ff0  = y(j0)
      ff1  = y(j1)
      ff2  = y(j2)
      tm2  = (ff2-ff1)/(xx2-xx1)
      tm1  = (ff0-ff1)/(xx0-xx1)
      bbb  = (tm2-tm1)/(xx2-xx0)
      ccc  =  tm2-bbb*(xx2-xx1)
      ymin = (bbb*xxd+ccc)*xxd+ff1

      rmin  = sqrt(xmin**2+ymin**2)
      thmin = 360.0*acos(xmin/rmin)/pi2

      if(ymin.lt.0)  thmin = -thmin
      if(thmin.lt.0) thmin = thmin + 180

c--
c final evaluation
c---

      dxy = (rmax-rmin)/(rmax+rmin)

c---
c find maximum distance along the z axis
c---

      zmax  = 0.0D0
      izmax = 0

      Do i=1,npts

        test = abs(p(i,3))

        if(test.gt.zmax) then
         zmax  = test
         izmax = i
        end if

      End Do

c-----
c done
c-----

      return
      end
