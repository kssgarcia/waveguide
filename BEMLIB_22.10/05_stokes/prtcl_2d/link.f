      subroutine Ilink 
     +
     +   (Nprtcl
     +   ,Iflow
     +   ,expn
     +   ,link
     +   )

c=============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=============================================

c---------------------------------------------
c Compute the pairwise associaton matrix: link
c
c link(i,j) = 1 if particles i and j
c               are associated
c             0 otherwise
c             0 for self-association
c
c Matrix ``link'' is symmetric
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xcntr(72),ycntr(72)
      Dimension axis1(72),axis2(72),tilt(72)

      Dimension link(72,72)

c--------------
c common blocks
c--------------

      common/CHANNELR/U1,U2,RL,h

      common/pax/axis1,axis2
      common/pap/xcntr,ycntr,tilt

      common/aaaa/a11,a12,a21,a22

c-----------
c initialize
c-----------

      Do i=1,Nprtcl
       link(i,i) = 2
      End Do

c--------------------------
c loop over particle pairs
c
c compare ith particle with all subsequent
c particles and their images
c--------------------------

      Do i=1,Nprtcl-1

         ax1_i = axis1(i)*expn
         ax2_i = axis2(i)*expn
       xcntr_i = xcntr(i)
       ycntr_i = ycntr(i)
        tilt_i =  tilt(i)

       Do j=i+1,Nprtcl

        ax1_j = axis1(j)*expn
        ax2_j = axis2(j)*expn
       tilt_j =  tilt(j)

       link(i,j) = 0

c---------------------
c singly periodic flow
c---------------------

       if(Iflow.eq.2.or.Iflow.eq.3) then

        ycntr_j = ycntr(j)

        Do ip=-1,1        ! periodic images

         xcntr_j = xcntr(j)+ip*RL

         call contact
     +
     +      (ax1_i,ax2_i,xcntr_i,ycntr_i,tilt_i
     +      ,ax1_j,ax2_j,xcntr_j,ycntr_j,tilt_j
     +      ,Tmin2,xmin2,ymin2
     +      ,Iloc
     +      )

         if(Iloc.eq.1) link(i,j) = 1

        End Do  ! periodic images

c---------------------
c doubly periodic flow
c---------------------

       else if(Iflow.eq.10) then

        Do ip=-1,1        ! periodic images
         Do jp=-1,1       ! periodic images

         xcntr_j = xcntr(j) + ip*a11+jp*a21
         ycntr_j = ycntr(j) + ip*a12+jp*a22

         call contact
     +
     +       (ax1_i,ax2_i,xcntr_i,ycntr_i,tilt_i
     +       ,ax1_j,ax2_j,xcntr_j,ycntr_j,tilt_j
     +       ,Tmin2,xmin2,ymin2
     +       ,Iloc
     +       )

         If(Iloc.eq.1) link(i,j) = 1

         End Do  ! periodic images
        End Do  ! periodic images

c------------------
c non-periodic flow
c------------------

       else

         xcntr_j = xcntr(j)
         ycntr_j = ycntr(j)

         call contact
     +
     +     (ax1_i,ax2_i,xcntr_i,ycntr_i,tilt_i
     +     ,ax1_j,ax2_j,xcntr_j,ycntr_j,tilt_j
     +     ,Tmin2,xmin2,ymin2
     +     ,Iloc
     +     )

         If(Iloc.eq.1) link(i,j) = 1

c-------------
       end if        ! end of flow options
c------------

       link(j,i) = link(i,j)

       End Do           ! loop over j
      End Do            ! loop over i

c-----
c Done
c-----

      return
      end
