      subroutine cluster_2d
     +
     +  (Nprtcl
     +  ,eps
     +  ,Ncls,lump
     +  ,Ipl
     +  )

c======================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c======================================

c----------------------------------------
c Arrange a collection particles into clusters
c
c SYMBOLS:
c --------
c
c Ipl:        initial particle label
c Ncls:       number of clusters
c link(i,j):  association matrix
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension lump(25),Ipl(25)
      Dimension link(25,25)

c------------------------------
c Compute the associaton matrix
c
c link(i,j) = 1 if particles i and j
c             are associated
c------------------------------

       call Ilink 
     +
     +  (Nprtcl
     +  ,eps
     +  ,link
     +  )

      write (6,*)
      write (6,*) "Association matrix computed"
      write (6,*)

c-------------------------------
c begin by assuming one cluster
c consisting of the first particle
c-------------------------------

      Ncls = 1
      lump(Ncls) = 1

      l = 0
      m = 1          ! lower limit for search

  1   Continue

      if(l.eq.Nprtcl) Go to 99

      l  = l+1
      l1 = l+1

c---
c compare particles l+1, ..., Nprtcl
c with particles m,...l
c---

      Do j=l1,Nprtcl
      Do i=m,l

c---
c if an association is found,
c place the jth particle at the l+1 position
c---

       if(link(i,j).eq.1) then
        lump(Ncls) = lump(Ncls) +1
        Isave   = Ipl(l1)
        Ipl(l1) = Ipl(j)
        Ipl(j)  = Isave

        Do k=1,Nprtcl
         Isave = link(l1,k)
         link(l1,k) = link(j,k)
         link(j,k) = Isave
        End Do

        Do k=1,Nprtcl
         Isave = link(k,l1)
         link(k,l1) = link(k,j)
         link(k,j) = Isave
        End Do

        Go to 1

       end if

      End Do
      End Do
       
c---
c If no association is found,
c introduce one more cluster
c---

      m    = l+1       ! redifine the lower search limit
      Ncls = Ncls+1
      lump(Ncls) = 1

      Go to 1

  99  Continue

c-----
c done
c-----

      return
      end
