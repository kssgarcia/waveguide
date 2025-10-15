      function ran2(idum)

c--------------------------------
c Returns a uniform random deviate
c between 0.0 and 1.0.
c
c Set IDUMto any negative value to initialize
c or reinitialize the sequence.
c--------------------------------

      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NRTAB,NDIV

      REAL ran2,AM,EPS,RNMX

      PARAMETER (IM1=2147483563
     +          ,IM2=2147483399
     +          ,AM=1./IM1
     +          ,IMM1=IM1-1
     +          ,IA1=40014
     +          ,IA2=40692
     +          ,IQ1=53668
     +          ,IQ2=52774
     +          ,IR1=12211
     +          ,IR2=3791
     +          ,NRTAB=32
     +          ,NDIV=1+IMM1/NRTAB
     +          ,EPS=1.2E-07
     +          ,RNMX=1.0-EPS
     +          )

      INTEGER idum2,j,k,iv(NRTAB),iy

      SAVE iv,iy,idum2

      DATA idum2/123456789/, iv/NRTAB*0/, iy/0/

c-------
c launch
c-------

      If (idum.le.0) then
          idum=max(-idum,1)
          idum2=idum
          Do j=NRTAB+8,1,-1
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           if (idum.lt.0) idum=idum+IM1
           if (j.le.NRTAB) iv(j)=idum
          End Do
          iy=iv(1)
      End If

      k=idum/IQ1

      idum=IA1*(idum-k*IQ1)-k*IR1

      If (idum.lt.0) idum=idum+IM1

      k=idum2/IQ2

      Idum2=IA2*(idum2-k*IQ2)-k*IR2

      If (idum2.lt.0) idum2=idum2+IM2

      j=1+iy/NDIV

      iy=iv(j)-idum2
      iv(j)=idum

      If(iy.lt.1)iy=iy+IMM1

      ran2=min(AM*iy,RNMX)

c-----
c Done
c-----

      Return
      End
