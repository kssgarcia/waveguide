      subroutine printel_stl(k)

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c-------------------------------------------
c print element k in file unit 8 (trgl6.stl)
c The element divided into four subelements
c------------------------------------------

      Implicit double precision (a-h,o-z)

      Dimension  p(1026,3)
      Dimension ne(1026,7)

      Dimension  n(512,6),nbe(512,3)

      common/points/p,ne
      common/elmnts/n,nbe

c------------------------------------
c draw four three-node flat triangles
c------------------------------------

c--- first

      apos1x = p(n(k,4),1)-p(n(k,1),1)
      apos1y = p(n(k,4),2)-p(n(k,1),2)
      apos1z = p(n(k,4),3)-p(n(k,1),3)
      apos2x = p(n(k,6),1)-p(n(k,1),1)
      apos2y = p(n(k,6),2)-p(n(k,1),2)
      apos2z = p(n(k,6),3)-p(n(k,1),3)
      vnormx = apos1y*apos2z-apos1z*apos2y
      vnormy = apos1z*apos2x-apos1x*apos2z
      vnormz = apos1x*apos2y-apos1y*apos2x
      vnormm = Dsqrt(vnormx*vnormx + vnormy*vnormy + vnormz*vnormz)
      vnormx = vnormx/vnormm
      vnormy = vnormy/vnormm
      vnormz = vnormz/vnormm

      write (8,91) vnormx,vnormy,vnormz
      write (8,92)
      i = n(k,1)
      write (8,93) p(i,1),p(i,2),p(i,3)
      i = n(k,4)
      write (8,94) p(i,1),p(i,2),p(i,3)
      i = n(k,6)
      write (8,95) p(i,1),p(i,2),p(i,3)
      write (8,96)
      write (8,97)

c--- second

      apos1x = p(n(k,2),1)-p(n(k,4),1)
      apos1y = p(n(k,2),2)-p(n(k,4),2)
      apos1z = p(n(k,2),3)-p(n(k,4),3)
      apos2x = p(n(k,5),1)-p(n(k,4),1)
      apos2y = p(n(k,5),2)-p(n(k,4),2)
      apos2z = p(n(k,5),3)-p(n(k,4),3)
      vnormx = apos1y*apos2z-apos1z*apos2y
      vnormy = apos1z*apos2x-apos1x*apos2z
      vnormz = apos1x*apos2y-apos1y*apos2x
      vnormm = Dsqrt(vnormx*vnormx + vnormy*vnormy + vnormz*vnormz)
      vnormx = vnormx/vnormm
      vnormy = vnormy/vnormm
      vnormz = vnormz/vnormm

      write (8,91) vnormx,vnormy,vnormz
      write (8,92)
      i = n(k,4)
      write (8,93) p(i,1),p(i,2),p(i,3)
      i = n(k,2)
      write (8,94) p(i,1),p(i,2),p(i,3)
      i = n(k,5)
      write (8,95) p(i,1),p(i,2),p(i,3)
      write (8,96)
      write (8,97)

c--- third

      apos1x = p(n(k,5),1)-p(n(k,4),1)
      apos1y = p(n(k,5),2)-p(n(k,4),2)
      apos1z = p(n(k,5),3)-p(n(k,4),3)
      apos2x = p(n(k,6),1)-p(n(k,4),1)
      apos2y = p(n(k,6),2)-p(n(k,4),2)
      apos2z = p(n(k,6),3)-p(n(k,4),3)
      vnormx = apos1y*apos2z-apos1z*apos2y
      vnormy = apos1z*apos2x-apos1x*apos2z
      vnormz = apos1x*apos2y-apos1y*apos2x
      vnormm = Dsqrt(vnormx*vnormx + vnormy*vnormy + vnormz*vnormz)
      vnormx = vnormx/vnormm
      vnormy = vnormy/vnormm
      vnormz = vnormz/vnormm

      write (8,91) vnormx,vnormy,vnormz
      write (8,92)
      i = n(k,4)
      write (8,93) p(i,1),p(i,2),p(i,3)
      i = n(k,5)
      write (8,94) p(i,1),p(i,2),p(i,3)
      i = n(k,6)
      write (8,95) p(i,1),p(i,2),p(i,3)
      write (8,96)
      write (8,97)

c--- fourth

      apos1x = p(n(k,5),1)-p(n(k,6),1)
      apos1y = p(n(k,5),2)-p(n(k,6),2)
      apos1z = p(n(k,5),3)-p(n(k,6),3)
      apos2x = p(n(k,3),1)-p(n(k,6),1)
      apos2y = p(n(k,3),2)-p(n(k,6),2)
      apos2z = p(n(k,3),3)-p(n(k,6),3)
      vnormx = apos1y*apos2z-apos1z*apos2y
      vnormy = apos1z*apos2x-apos1x*apos2z
      vnormz = apos1x*apos2y-apos1y*apos2x
      vnormm = Dsqrt(vnormx*vnormx + vnormy*vnormy + vnormz*vnormz)
      vnormx = vnormx/vnormm
      vnormy = vnormy/vnormm

      write (8,91) vnormx,vnormy,vnormz
      write (8,92)
      i = n(k,6)
      write (8,93) p(i,1),p(i,2),p(i,3)
      i = n(k,5)
      write (8,94) p(i,1),p(i,2),p(i,3)
      i = n(k,3)
      write (8,95) p(i,1),p(i,2),p(i,3)
      write (8,96)
      write (8,97)

c-----
c done
c-----

   91 Format(" facet normal",10(1x,f12.5))
   92 Format("   outer loop")
   93 Format("    vertex",10(1x,f12.5))
   94 Format("    vertex",10(1x,f12.5))
   95 Format("    vertex",10(1x,f12.5))
   96 Format("   endloop")
   97 Format(" endfacet")

  100 Format(1x,i4,10(1x,f12.5))

      return
      end
