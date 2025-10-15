      subroutine printel (k,Index,c)

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c-------------------------------------
c Print successive nodes of element k
c in file unit 1
c
c Index = 1: print the whole element
c Index = 2: print the 4 subelements
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  p(1026,3)
      Dimension ne(1026,7)
      Dimension  c(1026)

      Dimension n(512,6),nbe(512,3)

      common/points/p,ne
      common/elmnts/n,nbe

c----------
c constants
c----------

      nfour  = 4
      nseven = 7

c---------------------------------
c Draw one 6-point curved triangle
c---------------------------------

      if(Index.eq.1) then

c     write (6,100) nseven
      write (1,100) nseven
        i = n(k,1)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,4)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,2)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,5)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,3)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,6)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)
        i = n(k,1)
c       write (6,101) p(i,1),p(i,2),p(i,3)
        write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c---------------------------------
c Draw four 3-point flat triangles
c---------------------------------

      else if(Index.eq.2) then

c--- first

c     write (1,100) nfour
      i = n(k,1)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,6)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,1)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c--- second

c     write (1,100) nfour
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,2)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,5)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c--- third

c     write (1,100) nfour
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,5)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,6)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,4)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)

c--- fourth

c     write (1,100) nfour
      i = n(k,6)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,5)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,3)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)
      i = n(k,6)
      write (1,101) p(i,1),p(i,2),p(i,3),c(i)

      End If

c-----
c done
c-----

  100 Format(1x,i4,10(1x,f12.5))
  101 Format(10(1x,f12.5))

      return
      end
