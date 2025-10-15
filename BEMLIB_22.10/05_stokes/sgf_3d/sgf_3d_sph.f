      subroutine sgf_3d_sph
     +
     +   (x,y,z
     +   ,x0,y0,z0
     +   ,xc,yc,zc
     +   ,a
     +   ,gxx,gxy,gxz
     +   ,gyx,gyy,gyz
     +   ,gzx,gzy,gzz
     +   )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the lisensing agreement
c==========================================

c--------------------------------------------
c Green's function for flow in the exterior
c of a sphere centered
c at (xc,yc,zc)
c
c SYMBOLS:
c -------
c
c xc, yc, zc:  sphere center
c a:           sphere radius 
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

      pi = 3.14159 265358 D0

c------------------------------
c shift to the center and scale
c------------------------------

      xv = (x-xc)/a
      yv = (y-yc)/a
      zv = (z-zc)/a

      xs = (x0-xc)/a
      ys = (y0-yc)/a
      zs = (z0-zc)/a

c---------------
c velocity point
c---------------

      rv2 = xv*xv + yv*yv + zv*zv
      rv  = sqrt(rv2)

c---------------
c singular point
c---------------

      rs2 = xs*xs + ys*ys + zs*zs
      rs  = sqrt(rs2)
      rs3 = rs*rs2

c---------------
c distances etc
c---------------

      dx  = xv - xs
      dy  = yv - ys
      dz  = zv - zs
      dr2 = dx*dx + dy*dy + dz*dz
      dr  = Dsqrt(dr2)
      dr3 = dr*dr2

c---------------------
c image of singularity 
c---------------------

      xi  = xs/rs2
      yi  = ys/rs2
      zi  = zs/rs2
      ri2 = xi**2 + yi**2 + zi**2
      ri  = sqrt(ri2)

      dxi   = xv - xi
      dyi   = yv - yi
      dzi   = zv - zi
      dri2  = dxi**2 + dyi**2 + dzi**2
      dri   = sqrt(dri2)
      dri3  = dri*dri2
      dri3r = 1.0D0/dri3
      dri5  = dri3*dri2
      dri5r = 1.0D0/dri5

      f1 = (rs2-1.0)/rs
      f2 = 2.0*dri3r*(xi*dxi+yi*dyi+zi*dzi)
      f3 = (rv2-1.0)*(rs2-1.0)/(2.0*rs3) 
      f4 = dri3r/rs2
      f5 = f4/rs
      f6 = 3.0*f2/ dri2
      h7 = ri*dri+xv*xi+yv*yi+zv*zi-ri2
      f7 = 3.0*dri3r/(ri2*h7)
      f8 = f7*dri/h7
      h9 = rv*ri+xv*xi+yv*yi+zv*zi
      f9 = 3./(ri2*rv*h9)
      f0 = f9/h9

      rivx = ri*xv+rv*xi
      rivy = ri*yv+rv*yi
      rivz = ri*zv+rv*zi

      d53 = 3.0D0 * dri5r

      dxxi = d53*dxi*dxi 
      dxyi = d53*dxi*dyi 
      dxzi = d53*dxi*dzi 
      dyyi = d53*dyi*dyi 
      dyzi = d53*dyi*dzi 
      dzzi = d53*dzi*dzi 

      xsi3r = xs*dri3r
      ysi3r = ys*dri3r
      zsi3r = zs*dri3r
      driri  = dri*ri

      xidri2 = xi*dri2 - dxi*ri2 + driri*(dxi-xi)
      yidri2 = yi*dri2 - dyi*ri2 + driri*(dyi-yi)
      zidri2 = zi*dri2 - dzi*ri2 + driri*(dzi-zi)
      ridxi  = ri*dxi + dri*xi
      ridyi  = ri*dyi + dri*yi
      ridzi  = ri*dzi + dri*zi

      ooox = xi*dri2+dxi*ri2
      oooy = yi*dri2+dyi*ri2
      oooz = zi*dri2+dzi*ri2

      pppp = (dri-ri)*dri2*ri

      pxx = -3.*xsi3r*dxi + dri3r - dxxi 
     +      -2.*xsi3r* xi + xs*dxi*f6  
     +      +f7 *(dxi*ooox+pppp)
     +      -f8 * ridxi * xidri2
     +      -f9 *(xv*xi+rv*ri)
     +      +f0 * rivx  * rivx
      pxy = -3.*ysi3r*dxi         - dxyi
     +      -2.*ysi3r* xi + ys*dxi*f6  
     +      +f7 * dxi   * oooy
     +      -f8 * ridxi * yidri2
     +      -f9 * xv    * yi
     +      +f0 * rivx  * rivy
      pxz = -3.*zsi3r*dxi         - dxzi 
     +      -2.*zsi3r* xi + zs*dxi*f6  
     +      +f7 * dxi   * oooz
     +      -f8 * ridxi * zidri2
     +      -f9 * xv    * zi
     +      +f0 * rivx  * rivz
      pyx = -3.*xsi3r*dyi         - dxyi 
     +      -2.*xsi3r* yi + xs*dyi*f6  
     +      +f7 * dyi   * ooox
     +      -f8 * ridyi * xidri2
     +      -f9 * yv    * xi
     +      +f0 * rivy  * rivx
      pyy = -3.*ysi3r*dyi + dri3r - dyyi 
     +      -2.*ysi3r* yi + ys*dyi*f6  
     +      +f7 *(dyi*oooy+pppp)
     +      -f8 * ridyi * yidri2
     +      -f9 *(yv*yi+rv*ri)
     +      +f0 * rivy  * rivy
      pyz = -3.*zsi3r*dyi         - dyzi 
     +      -2.*zsi3r* yi + zs*dyi*f6  
     +      +f7 * dyi   * oooz
     +      -f8 * ridyi * zidri2
     +      -f9 * yv    * zi
     +      +f0 * rivy  * rivz
      pzx = -3.*xsi3r*dzi         - dxzi 
     +      -2.*xsi3r* zi + xs*dzi*f6  
     +      +f7 * dzi   * ooox
     +      -f8 * ridzi * xidri2
     +      -f9 * zv    * xi
     +      +f0 * rivz  * rivx
      pzy = -3.*ysi3r*dzi         - dyzi 
     +      -2.*ysi3r* zi + ys*dzi*f6  
     +      +f7 * dzi   * oooy
     +      -f8 * ridzi * yidri2
     +      -f9 * zv    * yi
     +      +f0 * rivz  * rivy
      pzz = -3.*zsi3r*dzi + dri3r - dzzi 
     +      -2.*zsi3r* zi + zs*dzi*f6  
     +      +f7 *(dzi*oooz+pppp)
     +      -f8 * ridzi * zidri2
     +      -f9 *(zv*zi+rv*ri)
     +      +f0 * rivz  * rivz

      pxx = pxx*f3
      pxy = pxy*f3
      pxz = pxz*f3
      pyx = pyx*f3 
      pyy = pyy*f3
      pyz = pyz*f3 
      pzx = pzx*f3 
      pzy = pzy*f3 
      pzz = pzz*f3 
c
c       write (6,*) ' ----------------'
c       write (6,*)   ' MATRIX P'
c       write (6,*) ' ----------------'
c       write (6,100) pxx,pxy,pxz
c       write (6,100) pyx,pyy,pyz
c       write (6,100) pzx,pzy,pzz
c
      drr    = 1.0D0/dr
      dr3r   = 1.0D0/dr3
      drir   = 1.0D0/dri
      rsdrir = drir/rs
      pppp   = drir+f2

      gxx = drr    + dx * dx * dr3r
     +     -rsdrir - dxi* dxi* f5
     + - f1* (xi*xi*pppp - f4*( xi*dxi + xi*dxi))
      gxy =          dx * dy * dr3r
     +             - dxi* dyi* f5
     + - f1* (xi*yi*pppp - f4*( xi*dyi + yi*dxi))
      gxz =          dx * dz * dr3r
     +             - dxi* dzi* f5
     + - f1* (xi*zi*pppp - f4*( xi*dzi + zi*dxi))
      gyy = drr    + dy * dy * dr3r
     +     -rsdrir - dyi* dyi* f5
     + - f1* (yi*yi*pppp - f4*( yi*dyi + yi*dyi))
      gyz =          dy * dz * dr3r
     +             - dyi* dzi* f5
     + - f1* (yi*zi*pppp - f4*( yi*dzi + zi*dyi))
      gzz = drr    + dz * dz * dr3r
     +     -rsdrir - dzi* dzi* f5
     + - f1* (zi*zi*pppp - f4*( zi*dzi + zi*dzi))

      gyx = gxy 
      gzx = gxz 
      gzy = gyz 

c-----------------------------------
c       write (6,*) ' MATRIX G:'
c       write (6,*)
c       write (6,100) gxx,gxy,gxz
c       write (6,100) gyx,gyy,gyz
c       write (6,100) gzx,gzy,gzz
c-----------------------------------

c---
c putting it together
c---

      gxx = gxx - pxx
      gxy = gxy - pxy
      gxz = gxz - pxz
      gyx = gyx - pyx
      gyy = gyy - pyy
      gyz = gyz - pyz
      gzx = gzx - pzx
      gzy = gzy - pzy
      gzz = gzz - pzz
 
      gxx = a * gxx
      gxy = a * gxy
      gxz = a * gxz
      gyx = a * gyx
      gyy = a * gyy
      gyz = a * gyz
      gzx = a * gzx
      gzy = a * gzy
      gzz = a * gzz

c-----
c done
c-----

      return
      end
