%---
close all
clear all
%---

figure(1)
hold on

figure(2)
hold on

file1 = fopen('flow_3x.out2');

Iflow=fscanf(file1,'%d ',[1,1]);

%----
NE1=fscanf(file1,'%d ',[1,1]);

Ic = 0;

while(NE1>0)

 Ic = Ic+1;
 clear X Y

 for i=1:NE1
  idle=fscanf(file1,'%d ',[1,1]);
  X(i)=fscanf(file1,'%f ',[1,1]);
  Y(i)=fscanf(file1,'%f ',[1,1]);
 end

 figure(1)
 plot( X,Y,'.-')
 plot(-X,Y,'.-')
 axis equal
 axis([-2 2 -2 2])

 for i=1:NE1-1
  slope(i)=(Y(i+1)-Y(i))/(X(i+1)-X(i));
  Ycl(i) = 0.5*(Y(i+1)+Y(i));
  Xcl(i) = 0.5*(X(i+1)+X(i));
 end

 figure(2)
 plot( Xcl,-2.0*slope)
 plot(-Xcl, 2.0*slope)
 axis([-1 1 -1 1])

 NE1=fscanf(file1,'%d ',[1,1]);

end
%----

%----
maxfny = 0.0;


NE=fscanf(file1,'%d ',[1,1]);

while(NE>0)

 seg = fscanf(file1,'%d ',[1,1]);

 clear Y0 fny

 for i=1:NE
  idle=fscanf(file1,'%d ',[1,1]);
  X0(i)=fscanf(file1,'%f ',[1,1]);
  Y0(i)=fscanf(file1,'%f ',[1,1]);
  S0(i)=fscanf(file1,'%f ',[1,1]);
  fx(i)=fscanf(file1,'%f ',[1,1]);
  fs(i)=fscanf(file1,'%f ',[1,1]);
  ff(i)=fscanf(file1,'%f ',[1,1]);
  fty(i)=fscanf(file1,'%f ',[1,1]);
  fny(i)=fscanf(file1,'%f ',[1,1]);
  ftz(i)=fscanf(file1,'%f ',[1,1]);
  if(abs(fny(i))>maxfny)
   maxfny = abs(fny(i));
  end
 end

 if(Iflow==41 & seg==1)
  win = 2.0*Y0(NE);
 elseif(Iflow==42 & seg==1)
  win = 2.0*Y0(NE);
 end

 figure(1)
 plot( Y0,X0,'o-')
 plot(-Y0,X0,'o-')
 axis equal
 axis([-win win -win win])
 xlabel('x')
 ylabel('y')
 box on

 figure(2)
 plot( Y0,-fny)
 plot(-Y0, fny)
 axis([-win win -1.2*maxfny 1.2*maxfny])
 xlabel('x')
 ylabel('f_n')
 box on

 NE=fscanf(file1,'%d ',[1,1]);

end
%----
