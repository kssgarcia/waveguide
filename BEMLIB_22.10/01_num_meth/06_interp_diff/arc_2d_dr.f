      program arc_2d_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------------------------------------
c Properties of a circular arc passing through 
c three points
c----------------------------------------

      write (6,*) 
      write (6,*) " Please enter:"
      write (6,*) 
      write (6,*) "  X1, Y1"
      write (6,*) "  X2, Y2"
      write (6,*) "  X3, Y3"
      write (6,*) "-------"

      read  (5,*) X1,Y1
      read  (5,*) X2,Y2
      read  (5,*) X3,Y3

c-----------
c Initialize
c-----------

      Iarc = 1   ! arc label

c-----------------------
c compute the properties
c of the arc
c-----------------------

      call arc_2d
     +
     +   (Iarc
     +   ,X1,X2,X3
     +   ,Y1,Y2,Y3
     +   ,XC,YC,R
     +   ,TH1,TH2,TH3,ORNT
     +   ,AREA1,AREA2,AREA
     +   ,XCN1,XCN2,XCN
     +   ,YCN1,YCN2,YCN
     +   ,Istop
     +   )

      If(Istop.eq.1) then
       write (6,*) 
       write (6,*) ' arc_2d_dr: something went wrong'
      End If

c-----------------
c printing session
c-----------------

      write (6,*) 
      write (6,140) xc,yc,R
      write (6,141) th1,th2,th3
      write (6,142) ornt
      write (6,143) area1,area2,area
      write (6,144) xcn1,xcn2,xcn
      write (6,145) ycn1,ycn2,ycn
      write (6,*) 

c-----
c Done
c-----

 140  Format (1x,'xc  = ',f10.5,' yc  =',f10.5,' R  =',f10.5)
 141  Format (1x,'th1 = ',f10.5,' th2 =',f10.5,' th3=',f10.5)
 142  Format (1x,'ornt= ',f10.5)
 143  Format (1x,'ar1 = ',f10.5,' ar2 =',f10.5,' ar =',f10.5)
 144  Format (1x,'xcn1= ',f10.5,' xcn2=',f10.5,' xcn=',f10.5)
 145  Format (1x,'ycn1= ',f10.5,' ycn2=',f10.5,' ycn=',f10.5)

      Stop
      End
