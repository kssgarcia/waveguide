close all
clear all

%===================
% driver for the free-space axisymmetric Green's function
% of two-dimensional Stokes flow
%===================

hold on

%---
X0=0.0;
Y0=1.2;
X =0.10;
Y =2.90;
%---

Iopt=1;

%==============
% confirm a symmetry property
%==============

[SXX,SXY,SYX,SYY ...
...
  ,QXXX,QXXY,QXYX,QXYY ...
  ,QYXX,QYXY,QYYX,QYYY ...
...
  ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
...
  = sgf_ax_fs (Iopt,X,Y,X0,Y0);

%---
% switch the singular point and the evaluation point
%---

[SXX1,SXY1,SYX1,SYY1 ...
...
  ,QXXX,QXXY,QXYX,QXYY ...
  ,QYXX,QYXY,QYYX,QYYY ...
...
  ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
...
  = sgf_ax_fs (Iopt,X0,Y0,X,Y);

[SXX,SXY;SYX,SYY]
Y/Y0*[SXX1,SXY1;SYX1,SYY1]

%==============
% assess the singular behavior
% by approaching the singular point
%==============

eps=0.50;

%-------

for i=1:22

 X = X0+eps;
 Y = Y0+eps;

 [SXX,SXY,SYX,SYY ...
 ...
   ,QXXX,QXXY,QXYX,QXYY ...
   ,QYXX,QYXY,QYYX,QYYY ...
 ...
   ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
 ...
  = sgf_ax_fs (Iopt,X,Y,X0,Y0);

%plot(log(eps),SXX,'o')
plot(log(eps),SXX/log(eps),'o')
plot(log(eps),SXY,'s')
plot(log(eps),SYX,'d')
plot(log(eps),SYY/log(eps),'+')
eps=0.50*eps;

end
%-------

%==============
% extract the singular behavior
%==============

 eps1=0.000001;
 X=X0+eps1;
 Y=Y0+eps1;

 [SXX1,SXY1,SYX1,SYY1 ...
 ...
   ,QXXX,QXXY,QXYX,QXYY ...
   ,QYXX,QYXY,QYYX,QYYY ...
 ...
   ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
 ...
  = sgf_ax_fs (Iopt,X,Y,X0,Y0);

 eps2=0.50*eps1;
 X=X0+eps2;
 Y=Y0+eps2;

 [SXX2,SXY2,SYX2,SYY2 ...
 ...
   ,QXXX,QXXY,QXYX,QXYY ...
   ,QYXX,QYXY,QYYX,QYYY ...
 ...
   ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
 ...
  = sgf_ax_fs (Iopt,X,Y,X0,Y0);

 A(1,1)=log(eps1); A(1,2)=1.0;
 A(2,1)=log(eps2); A(2,2)=1.0;
 rhs(1)=SXX1;
 rhs(2)=SXX2;
 cxx=rhs/A';
 rhs(1)=SYY1;
 rhs(2)=SYY2;
 cyy=rhs/A';

%---
% plot the asymptotics
%---

% SXXtest1 = cxx(1)*log(eps1)+cxx(2);
% SXXtest2 = cxx(1)*log(eps2)+cxx(2);
% plot(log(eps1),SXXtest1/log(eps1),'ro')
% plot(log(eps2),SXXtest2/log(eps2),'ro')

 for i=1:22
  epsl(i)=0.5/2^(i-1);
  SXXtest(i) = cxx(1)*log(epsl(i))+cxx(2);
  SYYtest(i) = cyy(1)*log(epsl(i))+cyy(2);
  plotx(i) = log(epsl(i));
  ploty1(i) = SXXtest(i);
  ploty1(i) = SXXtest(i)/log(epsl(i));
  ploty2(i) = SYYtest(i)/log(epsl(i));
 end

 plot(plotx,ploty1,'--')
 plot(plotx,ploty2,'--')
 xlabel('log(\epsilon)')
