close all
clear all

%---
wall=-1.0;
Xring=0.0;
Yring=1.0;
%---
Iopt=1;

dir=1;  % axial
dir=2;  % radial

%---
hold on
plot(Xring, Yring,'ro')
plot(Xring,-Yring,'ro')
%---


if(dir==1)
 Xstream =[-0.99 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 ...
0.01 0.1 0.2 0.3 0.4 ...
-0.80 -0.5]';
 Ystream =[ 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 ...
1.99 1.99 1.99 1.99 1.99 ...
0.10 0.20]';
 Nstr=size(Xstream);
 Nstep=2*128;
 Ds=0.02;
else
 Xstream =[-0.99 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 ...
0.01 0.1 0.2 0.3 0.4 ...
0.5 0.6 0.7 0.8 0.9 ...
-0.80 -0.5 -0.5 -0.5]';
 Ystream =[ 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 ...
1.99 1.99 1.99 1.99 1.99 ...
1.99 1.99 1.99 1.99 1.99 ...
0.10 0.20 0.5 0.7]';
 Nstr=size(Xstream);
 Nstep=2*128;
 Ds=-0.02;
end


%---
for istr=1:Nstr

clear Xstr Ystr

if(dir==1)
 Xstr(1)=Xstream(istr);
 Ystr(1)=Ystream(istr);
else
 Xstr(1)=Xstream(istr);
 Ystr(1)=Ystream(istr);
end

%------
for step=1:Nstep

Xvel=Xstr(step);
Yvel=Ystr(step);

[SXX,SXY,SYX,SYY ...
...
  ,QXXX,QXXY,QXYX,QXYY ...
  ,QYXX,QYXY,QYYX,QYYY ...
...
  ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
...
  = sgf_ax_w (Iopt,Xring,Yring,Xvel,Yvel,wall);

if(dir==1)
 velx=SXX;
 vely=SYX;
else
 velx=SXY;
 vely=SYY;
end

Dt = Ds/sqrt(velx*velx+vely*vely);

Xvel=Xstr(step)+velx*Dt;
Yvel=Ystr(step)+vely*Dt;

[SXX,SXY,SYX,SYY ...
...
  ,QXXX,QXXY,QXYX,QXYY ...
  ,QYXX,QYXY,QYYX,QYYY ...
...
  ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
...
  = sgf_ax_w (Iopt,Xring,Yring,Xvel,Yvel,wall);

if(dir==1)
 velx1=SXX;
 vely1=SYX;
else
 velx1=SXY;
 vely1=SYY;
end

Xstr(step+1)=Xstr(step)+0.5*(velx+velx1)*Dt;
Ystr(step+1)=Ystr(step)+0.5*(vely+vely1)*Dt;

if(Xstr(step+1)>2)
 break
end
if(Xstr(step+1)<-2)
 break
end
if(Ystr(step+1)>2)
 break
end
if(Ystr(step+1)<-2)
 break
end

end

plot(Xstr,Ystr)
plot(Xstr,-Ystr)

end
%---

plot([wall, wall],[-2 2])
plot([-2, 2],[0 0])
axis([-2 2 -2 2])
axis square
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
box



