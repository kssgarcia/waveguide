      subroutine tank_2d_sdlp 
     +
     +   (x0,y0
     +   ,ip,ie
     +   ,x1,y1
     +   ,x2,y2
     +   ,NGL
     +   ,SL    ! single-layer potential
     +   ,DL    ! double-layer potential
     +   )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c----------------------------------------------------
c Compute the single-layer and double-layer potential
c of the Neumann function over straight segments
c
c LEGEND:
c ------
c
c  ip: node label where the potentials are computed
c  ie: element label
c
c  vnxm, vnym:  normal vector at element mid-point
c  elml: element length
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension vnxm(400),vnym(400),elml(400)
      Dimension ZZ(20),WW(20)

c--------------
c common blocks
c--------------

      common/vnxym/vnxm,vnym,elml
      common/GF/wall1,wall2,wall3
      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

c--------
c prepare
c--------

      Iopt = 2 !   for the Green's function

c------------------
c singular element?
c------------------

      Ising = 0
      if(ip.eq.ie) Ising = 1

c-------------------------
c Gauss-Legende quadrature
c-------------------------

      xm = 0.5D0*(x2+x1)
      xd = 0.5D0*(x2-x1)

      ym = 0.5D0*(y2+y1)
      yd = 0.5D0*(y2-y1)

      SL  = 0.0D0
      DLx = 0.0D0
      DLy = 0.0D0

      !---
      Do i=1,NGL
      !---

        x = xm + ZZ(i)*xd
        y = ym + ZZ(i)*yd

        call lgf_2d_www
     +
     +     (Iopt
     +     ,x,y
     +     ,x0,y0
     +     ,wall1,wall2,wall3
     +     ,G
     +     ,Gx,Gy
     +     )

c------------------------
c singular element: subtract off the singularity
c
c the free-space Green's function is: 
c
c   G = -1.0/(2*pi) * lnr
c------------------------

        if(Ising.eq.1) then
          dx = x-x0
          dy = y-y0
          r2 = dx*dx+dy*dy
          G = G + 0.50D0*log(r2)/pi2
        end if

        SL  = SL  + G *WW(i)
        DLx = DLx + Gx*WW(i)
        DLy = DLy + Gy*WW(i)

c       write (6,*)
c       write (6,*) wall1,wall2,wall3
c       write (6,*) x,y,x0,y0
c       write (6,*) Gx,Gy
c       write (6,*)

      !---
      End Do
      !---

      elmlh = 0.5D0*elml(ie)

      SL = elmlh*SL
      DL = elmlh*( DLx*vnxm(ie) + DLy*vnym(ie) )

c-------------------------------
c if the element is singular,
c add back singular contribution
c to the single layer
c (the contribution from the dlp is zero)
c-------------------------------

      if(Ising.eq.1) then
        SL = SL - elml(ie)*(log(elmlh)-1.0D0)/pi2
      end if

c-----
c done
c-----

      Return
      End
