%======================================
% driver for the axisymmetric
% Green's function of Stokes flow
% in semi-infinite space above an
% infinite plane wall located at x=wall
%======================================

wall=-2.0;
Iopt=2;

%---
% evaluation point
%---
X0=0.2;
Y0=0.5;

%---
% singular point
%---
X=0.0;
Y=1.0;

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
  = sgf_ax_w (Iopt,X0,Y0,X,Y,wall);

Y0*[SXX,SXY,SYX,SYY]

%---
% switch the singular point and the evaluation point
%---

[SXX,SXY,SYX,SYY ...
...
    ,QXXX,QXXY,QXYX,QXYY ...
    ,QYXX,QYXY,QYYX,QYYY ...
...
    ,PXX,PXY,PYX,PYY ...
    ,Iaxis] ...
...
  = sgf_ax_w (Iopt,X,Y,X0,Y0,wall);

Y*[SXX,SXY,SYX,SYY]

