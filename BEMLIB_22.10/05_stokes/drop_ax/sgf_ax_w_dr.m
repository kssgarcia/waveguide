close all
clear all

hold on

wall=0.0;
Iopt=1;
X=0.101;
Y=1.001;
X0=0.1;
Y0=1.0;

X=0.919579069589760; Y= 0.485528934340945;
X0=0.919581277763687; Y0= 0.485522319153874;

[SXX,SXY,SYX,SYY ...
...
        ,QXXX,QXXY,QXYX,QXYY ...
        ,QYXX,QYXY,QYYX,QYYY ...
...
        ,PXX,PXY,PYX,PYY ...
        ,Iaxis] ...
...
  = sgf_ax_w (Iopt,X0,Y0,X,Y,wall);

[SXX,SXY,SYX,SYY]

[QXXX,QXXY,QXYX,QXYY; ...
 QYXX,QYXY,QYYX,QYYY]

%==============
% assess the singular behavior
%==============

wall = 0.0;
X0=1.0;
Y0=1.0;
X0=0.919581277763687; Y0= 0.485522319153874;

eps=0.50;

%-------

for i=1:12

 X=X0+eps;
 Y=Y0+eps;

[SXX,SXY,SYX,SYY ...
...
        ,QXXX,QXXY,QXYX,QXYY ...
        ,QYXX,QYXY,QYYX,QYYY ...
...
        ,PXX,PXY,PYX,PYY ...
        ,Iaxis] ...
...
  = sgf_ax_w (Iopt,X,Y,X0,Y0,wall);

X
Y

%plot(log(eps),SXX,'o')
plot(log(eps),SXX/log(eps),'o')
plot(log(eps),SXY,'s')
plot(log(eps),SYX,'d')
plot(log(eps),SYY/log(eps),'+')
eps=0.50*eps;

end
%-------

