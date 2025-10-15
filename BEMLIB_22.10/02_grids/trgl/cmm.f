      subroutine cmm
     +
     +   (Npts,Nelm
     +   ,mint
     +   ,Gcm,Gmm
     +   )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c---------------------------------
c Computatie the:
c
c Global conductivity matrix
c Global mass matrix
c
c LEGEND:
c ------
c
c ecm (512,6,6): element conductivity matrix
c emm (512,6,6): element mass matrix
c ph(i):         element basis function psi, i = 1, ..., 6
c gph(i,3):      surface gradient of psi(i)
c---------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)

      Dimension     n(512,6),nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xiq(20),etq(20),wq(20)

      Dimension ph(6),gph(6,3)
      Dimension ecm(6,6),emm(6,6)

      Dimension Gcm(1026,1026),Gmm(1026,1026)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/trq/xiq,etq,wq

c-----------
c initialize
c-----------

      Do k=1,Npts
       Do l=1,Npts
         Gcm(k,l) = 0.0D0
         Gmm(k,l) = 0.0D0
       End Do
      End Do

c-----------------------
c loop over the elements
c-----------------------

      Do k=1,Nelm

      i1 = n(k,1)   ! global element node labels
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be = beta (k)
      ga = gamma(k)

c---
c launch the quadrature
c---

       Do j1=1,6
        Do j2=1,6
          ecm(j1,j2) = 0.0D0
          emm(j1,j2) = 0.0D0
        End Do
       End Do

c---
c loop
c---

      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        call cmm_interp_elm
     +
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +
     +    ,al,be,ga
     +    ,xi,eta
     +
     +    ,ph,gph
     +    ,hs
     +    )

       cf = 0.5D0*hs*wq(i)

       Do j1=1,6
        Do j2=j1,6
          ecm(j1,j2) = ecm(j1,j2) +  ph(j1)*ph(j2)* cf
          emm(j1,j2) = emm(j1,j2) + (  gph(j1,1)*gph(j2,1)
     +                               + gph(j1,2)*gph(j2,2)
     +                               + gph(j1,3)*gph(j2,3) ) * cf
        End Do
       End Do

      End Do

c---

      Do j1=2,6
       Do j2=1,j1-1
         ecm(j1,j2) = ecm(j2,j1)
         emm(j1,j2) = emm(j2,j1)
       End Do
      End Do


c---
c quadrature finished
c---

      Gcm(i1,i1) = Gcm(i1,i1) + ecm(1,1)
      Gcm(i1,i2) = Gcm(i1,i2) + ecm(1,2)
      Gcm(i1,i3) = Gcm(i1,i2) + ecm(1,3)
      Gcm(i1,i4) = Gcm(i1,i4) + ecm(1,4)
      Gcm(i1,i5) = Gcm(i1,i5) + ecm(1,5)
      Gcm(i1,i6) = Gcm(i1,i6) + ecm(1,6)

      Gcm(i2,i1) = Gcm(i2,i1) + ecm(2,1)
      Gcm(i2,i2) = Gcm(i2,i2) + ecm(2,2)
      Gcm(i2,i3) = Gcm(i2,i2) + ecm(2,3)
      Gcm(i2,i4) = Gcm(i2,i4) + ecm(2,4)
      Gcm(i2,i5) = Gcm(i2,i5) + ecm(2,5)
      Gcm(i2,i6) = Gcm(i2,i6) + ecm(2,6)

      Gcm(i3,i1) = Gcm(i3,i1) + ecm(3,1)
      Gcm(i3,i2) = Gcm(i3,i2) + ecm(3,2)
      Gcm(i3,i3) = Gcm(i3,i2) + ecm(3,3)
      Gcm(i3,i4) = Gcm(i3,i4) + ecm(3,4)
      Gcm(i3,i5) = Gcm(i3,i5) + ecm(3,5)
      Gcm(i3,i6) = Gcm(i3,i6) + ecm(3,6)

      Gcm(i4,i1) = Gcm(i4,i1) + ecm(4,1)
      Gcm(i4,i2) = Gcm(i4,i2) + ecm(4,2)
      Gcm(i4,i3) = Gcm(i4,i2) + ecm(4,3)
      Gcm(i4,i4) = Gcm(i4,i4) + ecm(4,4)
      Gcm(i4,i5) = Gcm(i4,i5) + ecm(4,5)
      Gcm(i4,i6) = Gcm(i4,i6) + ecm(4,6)

      Gcm(i5,i1) = Gcm(i5,i1) + ecm(5,1)
      Gcm(i5,i2) = Gcm(i5,i2) + ecm(5,2)
      Gcm(i5,i3) = Gcm(i5,i2) + ecm(5,3)
      Gcm(i5,i4) = Gcm(i5,i4) + ecm(5,4)
      Gcm(i5,i5) = Gcm(i5,i5) + ecm(5,5)
      Gcm(i5,i6) = Gcm(i5,i6) + ecm(5,6)

      Gcm(i6,i1) = Gcm(i6,i1) + ecm(6,1)
      Gcm(i6,i2) = Gcm(i6,i2) + ecm(6,2)
      Gcm(i6,i3) = Gcm(i6,i2) + ecm(6,3)
      Gcm(i6,i4) = Gcm(i6,i4) + ecm(6,4)
      Gcm(i6,i5) = Gcm(i6,i5) + ecm(6,5)
      Gcm(i6,i6) = Gcm(i6,i6) + ecm(6,6)

c---

      Gmm(i1,i1) = Gmm(i1,i1) + emm(1,1)
      Gmm(i1,i2) = Gmm(i1,i2) + emm(1,2)
      Gmm(i1,i3) = Gmm(i1,i2) + emm(1,3)
      Gmm(i1,i4) = Gmm(i1,i4) + emm(1,4)
      Gmm(i1,i5) = Gmm(i1,i5) + emm(1,5)
      Gmm(i1,i6) = Gmm(i1,i6) + emm(1,6)

      Gmm(i2,i1) = Gmm(i2,i1) + emm(2,1)
      Gmm(i2,i2) = Gmm(i2,i2) + emm(2,2)
      Gmm(i2,i3) = Gmm(i2,i2) + emm(2,3)
      Gmm(i2,i4) = Gmm(i2,i4) + emm(2,4)
      Gmm(i2,i5) = Gmm(i2,i5) + emm(2,5)
      Gmm(i2,i6) = Gmm(i2,i6) + emm(2,6)

      Gmm(i3,i1) = Gmm(i3,i1) + emm(3,1)
      Gmm(i3,i2) = Gmm(i3,i2) + emm(3,2)
      Gmm(i3,i3) = Gmm(i3,i2) + emm(3,3)
      Gmm(i3,i4) = Gmm(i3,i4) + emm(3,4)
      Gmm(i3,i5) = Gmm(i3,i5) + emm(3,5)
      Gmm(i3,i6) = Gmm(i3,i6) + emm(3,6)

      Gmm(i4,i1) = Gmm(i4,i1) + emm(4,1)
      Gmm(i4,i2) = Gmm(i4,i2) + emm(4,2)
      Gmm(i4,i3) = Gmm(i4,i2) + emm(4,3)
      Gmm(i4,i4) = Gmm(i4,i4) + emm(4,4)
      Gmm(i4,i5) = Gmm(i4,i5) + emm(4,5)
      Gmm(i4,i6) = Gmm(i4,i6) + emm(4,6)

      Gmm(i5,i1) = Gmm(i5,i1) + emm(5,1)
      Gmm(i5,i2) = Gmm(i5,i2) + emm(5,2)
      Gmm(i5,i3) = Gmm(i5,i2) + emm(5,3)
      Gmm(i5,i4) = Gmm(i5,i4) + emm(5,4)
      Gmm(i5,i5) = Gmm(i5,i5) + emm(5,5)
      Gmm(i5,i6) = Gmm(i5,i6) + emm(5,6)

      Gmm(i6,i1) = Gmm(i6,i1) + emm(6,1)
      Gmm(i6,i2) = Gmm(i6,i2) + emm(6,2)
      Gmm(i6,i3) = Gmm(i6,i2) + emm(6,3)
      Gmm(i6,i4) = Gmm(i6,i4) + emm(6,4)
      Gmm(i6,i5) = Gmm(i6,i5) + emm(6,5)
      Gmm(i6,i6) = Gmm(i6,i6) + emm(6,6)

c------------------------------
c End of loop over the elements
c------------------------------

      End Do

c-----
c done
c-----

      return
      end
