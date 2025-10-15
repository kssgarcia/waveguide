      program lgf_3d_2p_glut

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c============================================

c--------------------------------------
c Generate a look up table for the non-singular
c part of the doubly-periodic green's function
c of three-dimensional Laplace's equation.
c
c Symbols:
c --------
c
c method = 1 for the straight fourier series method
c method = 2 for the fast summation method
c            of Hautman & Klein
c--------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension glut  (-64:64,-64:64,0:64)
      Dimension glut_x(-64:64,-64:64,0:64)
      Dimension glut_y(-64:64,-64:64,0:64)
      Dimension glut_z(-64:64,-64:64,0:64)

      common/glut_r/glut,glut_x,glut_y,glut_z
      common/glut_i/nxxx,nyyy,nzzz

c--------
c prepare
c--------

 96   continue

      write (6,*)
      write (6,*) "  Please enter the components of the"
      write (6,*) "         first base vector: a11, a12"
      write (6,*) "  ----------------------------------"
      read  (5,*) a11,a12

      write (6,*)
      write (6,*) "  Please enter the components of the"
      write (6,*) "         second base vector: a21, a22"
      write (6,*) "  -----------------------------------"
      read  (5,*) a21,a22

c--------------------------
c ewald_3d_2p will generate:
c
c     the reciprocal vectors
c     the optimal value of xi (called ew)
c     and the area of the unit cell
c--------------------------

      call lgf_3d_2p_ewald
     +
     +   (a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,area
     +   )

c---
c displaying
c---

      write (6,*)
      write (6,*) " Coordinates of the reciprocal lattice vectors"
      write (6,*) " ---------------------------------------------"
      write (6,100) b11,b12
      write (6,100) b21,b22

c---

 97   continue

      write (6,*) "  Choose the method"
      write (6,*)
      write (6,*) "  Enter:"
      write (6,*)
      write (6,*) "  1 for straight fourier method"
      write (6,*) "  2 for the fast-summation method"
      write (6,*) "  0 to quit"
      write (6,*) " ----------"
      read  (5,*) method

      if(method.eq.0) Go to 99

c------------------------------
      if(method.eq.1) then
c------------------------------

         write (6,*)
         write (6,*) " Enter the truncation limit for summation"
         write (6,*) "       over the reciprocal lattice"
         write (6,*) " 0 to quit"
         write (6,*) " ---------"
         read  (5,*) max2

         if(max2.eq.0) Go to 99

c------------------------------
      else if(method.eq.2) then
c------------------------------

         write (6,*)
         write (6,*) " Enter the truncation limit for summation"
         write (6,*) "       over the physical lattice in the xy plane"
         write (6,*) "  0 to quit"
         write (6,*) " ----------"
         read  (5,*) max1

         if(max1.eq.0) Go to 99

         write (6,*)
         write (6,*) " Enter the truncation limit for summation"
         write (6,*) "       over the reciprocal lattice"
         write (6,*) "  0 to quit"
         write (6,*) " ----------"
         read  (5,*) max2

         if(max2.eq.0) Go to 99

         write (6,*)
         write (6,*) " Enter the truncation limit for summation"
         write (6,*) "       over the 3d physical lattice "
         write (6,*) "  0 to quit"
         write (6,*) " ----------"
         read  (5,*) max3

         If(max3.eq.0) Go to 99

c-----------
      end if
c-----------

 98   Continue

      write (6,*)
      write (6,*) "will produce a look up table with respect to:"
      write (6,*) 
      write (6,*) "  1) x-x0 in the range of [-0.5, 0.5]"
      write (6,*) "  1) y-y0 in the range of [-0.5, 0.5]"
      write (6,*) "  1) z-z0 in the range of [ 0.0, 0.5]"
      write (6,*) 
      write (6,*) "Please enter half the number of divisions"
      write (6,*) "       for each tabulation variable"
      write (6,*) " ---------------------------------------"
      read  (5,*) nxxx,nyyy,nzzz 

c---------------
c initiate tabulation
c---------------

c---
c prepare
c---

      Iselect = 2

c---
c one point source at the origin
c---

      x0 = 0.0D0
      y0 = 0.0D0
      z0 = 0.0D0
  
c---
c scan the interpolation nodes
c---

      Dxxx = 0.5/nxxx
      Dyyy = 0.5/nyyy
      Dzzz = 0.5/nzzz

      Do k=0,Nzzz

         z = z0 + k*Dzzz

         Do j=-Nyyy,Nyyy

            y = y0 + j*Dyyy

            Do i=-Nxxx,Nxxx

               x = x0 + i*Dxxx

              if(i.eq.0.and.j.eq.0.and.k.eq.0)then
                x = x0 + 0.0000001
                y = y0 + 0.0000001
                z = z0 + 0.0000001
              end if 

              call lgf_3d_2p
     +
     +            (Iselect
     +            ,method
     +            ,x,y,z
     +            ,x0,y0,z0
     +            ,a11,a12,a21,a22
     +            ,b11,b12,b21,b22
     +            ,ew,area
     +            ,max1,max2,max3
     +            ,g
     +            ,gx,gy,gz
     +            )

              glut  (i,j,k) = G
              glut_x(i,j,k) = Gx
              glut_y(i,j,k) = Gy
              glut_z(i,j,k) = Gz

c      write (6,*) " -----------------------------"
c      write (6,*) " green's function and gradient "
c      write (3,*) i,j,k
c      write (3,100)
c      write (3,100) g
c      write (3,100) gx,gy,gz
c      write (6,*) " -----------------------------"

         End Do 
       End Do
      End Do 

c--------------------------------
c End of look-up table generation
c--------------------------------


c--------------------------------------------
c Record the data in file: lgf_3d_2p_glut.out
c--------------------------------------------

      write (6,*)
      write (6,*) " Record the interpolation tables in file:"
      write (6,*) "        lgf_3d_2p_glut.out ?"
      write (6,*)
      write (6,*) " Enter 0 for NO, 1 for YES"
      write (6,*) "--------------------------"
      read  (5,*) Irecord

      if(Irecord.eq.1) then

       open (3,file="lgf_3d_2p.glut")

       write (3,*) Na21,Nxxx,Nyyy

       Do k=0,Nzzz
         Do j=-Nyyy,Nyyy
          Do i=-Nxxx,Nxxx

            write (3,104) Glut(i,j,k)
     +                   ,Glut_x(i,j,k),Glut_y(i,j,k),Glut_z(i,j,k)

          End Do
         End Do
       End Do

       close (3)

      End If

c-----------------------------------------
c Evaluate the Green's function at a point
c by trilinear interpolation
c-----------------------------------------

      write (6,*)
      write (6,*) ' Will interpolate'
      write (6,*)

 95   Continue

      write (6,*)
      write (6,*) " Enter (x,y,z) of the field point" 
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) x,y,z 

      call lgf_3d_2p_int
     +
     +   (x,y,z
     +   ,x0,y0,z0
     +   ,a11,a12,a21,a22
     +   ,g
     +   ,Gx,Gy,Gz
     +   )
            
      write (6,*) ' --------------------------------'
      write (6,*) ' Green function and gradient'
      write (6,*)
      write (6,*) ' Interpolated values:'
      write (6,*)
      write (6,100) G
      write (6,100) Gx,Gy,Gz

c---
c Direct evaluation for verification
c---

       call lgf_3d_2p
     +
     +    (Iselect
     +    ,method
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,a11,a12,a21,a22
     +    ,b11,b12,b21,b22
     +    ,ew,area
     +    ,max1,max2,max3
     +    ,g
     +    ,gx,gy,gz
     +    )

      write (6,*)
      write (6,*) ' Directly computed values:'
      write (6,*)
      write (6,100) G
      write (6,100) Gx,Gy,Gz
      write (6,*)
      write (6,*) ' --------------------------------'


      Go to 95

c-----
c Done
c-----

 99   Continue

  100 Format (3(1x,f20.10))
  101 Format (1x,i3,10(1x,f9.5))
  102 Format (t1,3(3x,e16.10))
  104 Format (8(1x,f15.10))

      stop
      end
