      subroutine taylor 
     +
     +  (NSG
     +  ,xc,yc
     +  ,Dxy
     +  ,thmax,thmin
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------------------------------
c computes the Taylor deformation parameter
c and inclination
c of the major and minor axes of a drop
c by parabolic interpolation
c
c SYMBOLS:
c -------
c
c xc,yc:   coordinates of the centroid (must be supplied)
c
c rr:	radial position from the centroid
c Dxy:	Deformation parameter
c----------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    X(0:200),  Y(0:200),rr(0:200)
      Dimension Xhat(0:200),Yhat(0:200)

      Parameter (tol=0.000001)

      common/XXYY/X,Y
      common/piii/pi,pih,piq,pi2,pi4,pi6,pi8,srpi

c---
c prepare
c---

      NSG1 = NSG+1
      NSG2 = NSG+2

      Do i=1,NSG1
       rr(i) = sqrt((x(i)-xc)**2+(y(i)-yc)**2 )
      End Do

       x(0) =  x(nsg)
       y(0) =  y(nsg)
      rr(0) = rr(nsg)

       x(nsg2) =  x(2)
       y(nsg2) =  y(2)
      rr(nsg2) = rr(2)

c---
c Find points with
c maximum and minimum radial distance
c from the centroid
c---

      rmax =      0.0  
      rmin = 100000.0 
      imax = 0
      imin = 0

      Do i=1,NSG-1
         if(rr(i).gt.rmax) then
           rmax  = rr(i)
           imax  = i
         end if
         if(rr(i).lt.rmin) then
           rmin = rr(i)
           imin  = i
         end if
      End Do

c---
c find the coordinates of the point with the
c maximum distance by quadratic interpolation
c with respect to the poly-line
c---

      ja = imax-1
      j1 = imax
      j2 = imax+1

      xx0 = - sqrt((x(ja)-x(j1))**2+(y(ja)-y(j1))**2)
      xx1 = 0.0
      xx2 =   sqrt((x(j2)-x(j1))**2+(y(j2)-y(j1))**2)

      ff0 = rr(ja)
      ff1 = rr(j1)
      ff2 = rr(j2)
      tm2 = (ff2-ff1)/(xx2-xx1)
      tm1 = (ff0-ff1)/(xx0-xx1)
      bbb = (tm2-tm1)/(xx2-xx0)
      ccc =  tm2-bbb*(xx2-xx1)

      if(abs(bbb).le.tol.and.abs(ccc).le.tol) then
       xxmx = xx1
       rmax = ff1
      else
       xxmx = -0.50*ccc   /bbb + xx1
       rmax = -0.25*ccc**2/bbb + ff1
      end if
     
      xxd  = xxmx-xx1

c---
c Interpolate for the x and y
c position with respect to polygonal
c arc length
c---

      ff0  = x(ja)
      ff1  = x(j1)
      ff2  = x(j2)

      tm2  = (ff2-ff1)/(xx2-xx1)
      tm1  = (ff0-ff1)/(xx0-xx1)

      bbb  = (tm2-tm1)/(xx2-xx0)
      ccc  =  tm2-bbb*(xx2-xx1)
      xmax = (bbb*xxd+ccc)*xxd+ff1

      ff0  = y(ja)
      ff1  = y(j1)
      ff2  = y(j2)

      tm2  = (ff2-ff1)/(xx2-xx1)
      tm1  = (ff0-ff1)/(xx0-xx1)

      bbb  = (tm2-tm1)/(xx2-xx0)
      ccc  =  tm2-bbb*(xx2-xx1)
      ymax = (bbb*xxd+ccc)*xxd+ff1

      thmax = xmax/rmax

      if(thmax.ge. 1.0) thmax = 0.9999999    ! safety valve
      if(thmax.le.-1.0) thmax =-0.9999999    ! safety valve

      thmax = 360.0*acos(thmax)/pi2

      if(ymax.lt.0)  thmax = -thmax
      if(thmax.lt.0) thmax = thmax + 180.0

c---
c compute the coordinates of the point with the
c minimum distance by quadratic interpolation
c---

      ja = imin-1
      j1 = imin
      j2 = imin+1

      xx0 = -sqrt((x(ja)-x(j1))**2+(y(ja)-y(j1))**2)
      xx1 = 0.0
      xx2 =  sqrt((x(j2)-x(j1))**2+(y(j2)-y(j1))**2)

      ff0 = rr(ja)
      ff1 = rr(j1)
      ff2 = rr(j2)

      tm2 = (ff2-ff1)/(xx2-xx1)
      tm1 = (ff0-ff1)/(xx0-xx1)

      bbb = (tm2-tm1)/(xx2-xx0)
      ccc =  tm2-bbb*(xx2-xx1)

      if(abs(bbb).le.tol.and.abs(ccc).le.tol) then
       xxmn = xx1
       rmin = ff1
      else
       xxmn = -0.50*ccc   /bbb + xx1
       rmin = -0.25*ccc**2/bbb + ff1
      end if

      xxd  = xxmn-xx1

c---
c Interpolate for the x and y
c position with respect to polygonal
c arc length
c---

      ff0  = x(ja)
      ff1  = x(j1)
      ff2  = x(j2)

      tm2  = (ff2-ff1)/(xx2-xx1)
      tm1  = (ff0-ff1)/(xx0-xx1)

      bbb  = (tm2-tm1)/(xx2-xx0)
      ccc  =  tm2-bbb*(xx2-xx1)
      xmin = (bbb*xxd+ccc)*xxd+ff1

      ff0  = y(ja)
      ff1  = y(j1)
      ff2  = y(j2)

      tm2  = (ff2-ff1)/(xx2-xx1)
      tm1  = (ff0-ff1)/(xx0-xx1)
      bbb  = (tm2-tm1)/(xx2-xx0)

      ccc  =  tm2-bbb*(xx2-xx1)
      ymin = (bbb*xxd+ccc)*xxd+ff1

      thmin = xmin/rmin

      if(thmin.ge. 1.0) thmin = 0.9999999
      if(thmin.le.-1.0) thmin =-0.9999999

      thmin = 360.0*acos(thmin)/pi2

      if(ymin.lt.0)  thmin = -thmin
      if(thmin.lt.0) thmin = thmin + 180.0

c--
c final evaluation
c---

      Dxy  = (rmax-rmin)/(rmax+rmin)

c--------------------------------------------------------------
c Alternative approach
c       diagonalize inertia tensor (2-D) to find directions of
c       max. and min. elongation.
c--------------------------------------------------------------

      Do i=1,NSG1
       xhat(i) = x(i)-xc
       yhat(i) = y(i)-yc
      End Do

      axx = 0.0
      ayy = 0.0
      axy = 0.0

      Do i=1,NSG
       i1 = i+1
       dl=sqrt((x(i1)-x(i))**2+(y(i1)-y(i))**2)
       axx=axx +(xhat(i)**2 +    xhat(i1)**2)      * dl
       ayy=ayy +(yhat(i)**2     + yhat(i1)**2)     * dl
       axy=axy +(xhat(i)*yhat(i)+xhat(i1)*yhat(i1))* dl
      End Do

      axx = 0.5*axx
      ayy = 0.5*ayy
      axy = 0.5*axy

c-----------------------------------
c find eigenvalues of inertia tensor
c-----------------------------------

        if(abs(axy).gt.0.0000000001) then

         t1    = axx+ayy
         t2    = axx*ayy-axy**2
         discr = t1**2-4.0*t2
         tmp   = sqrt(discr)
         amax  = 0.5D0*(t1+tmp)
         amin  = 0.5D0*(t1-tmp)

         aspra = amax/amin
c---
c find eigenvectors of inertia tensor and orientation of the
c principal directions relative to the x axis
c---

         v1    = 1.0D0
         v2    = -(axx-amax)/axy
         vnorm = sqrt(v1**2+v2**2)
         v1    = v1/vnorm
         v2    = v2/vnorm
         thmax = acos(v1)*180.0/pi
c---

         v1    = 1.0D0
         v2    = -(axx-amin)/axy
         vnorm = sqrt(v1**2+v2**2)
         v1    = v1/vnorm
         v2    = v2/vnorm
         thmin = acos(v1)*180.0/pi
         thmin = 180.0-thmin

        else

         aspra =   1.0
         thmax =  45.0
         thmin = 135.0

        end if

c-----
c done
c-----

 101  Format (1x,i3,3(1x,f10.5))

      return
      end
