close all
clear all

%------------
% FDLIB BEMLIB
%
% draw streamlines due to a ring
% of point forces in front of a wall
% located at x=wall
%------------


Iopt=1; % compute only the Green's function

%---
wall=-1.0;
Xring=0.0;
Yring=1.0;
%---

%---
% preferences
%---

dir=1;  % axial
dir=2;  % radial

%---
figure(1)
hold on
plot([wall, wall],[-2 2],'k')
plot([-2, 2],[0 0],'k')
axis([-2 2 -2 2])
axis square
set(gca,'fontsize',14)
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
box
plot(Xring, Yring,'ko')
plot(Xring,-Yring,'ko')
%---


%---
if(dir==1)
%---
 Xstream =[-0.99 -0.9 -0.8 -0.7 ...
           -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 ...
            0.01 0.1 0.2 0.3 0.4 ...
           -0.80 -0.5]';
 Ystream =[ 1.99 1.99 1.99 1.99 ...
           1.99 1.99 1.99 1.99 1.99 1.99 ...
           1.99 1.99 1.99 1.99 1.99 ...
           0.10 0.20]';
 Nstr=size(Xstream);
 Nstep=2*128;
 Ds=0.02;
%---
else
%---
 Xstream =[-0.99 -0.9 -0.8 -0.7  ...
           -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 ...
            0.01 0.1 0.2 0.3 0.4 ...
            0.5 0.6 0.7 0.8 0.9 ...
            -0.80 -0.5 -0.5 -0.5]';
 Ystream =[ 1.99 1.99 1.99 1.99 ...
            1.99 1.99 1.99 1.99 1.99 1.99 ...
            1.99 1.99 1.99 1.99 1.99 ...
            1.99 1.99 1.99 1.99 1.99 ...
            0.10 0.20 0.5 0.7]';
 Nstr=size(Xstream);
 Nstep=2*128;
 Ds=-0.02;
%---
end
%---


%---
for istr=1:Nstr

 clear Xstr Ystr

 Xstr(1)=Xstream(istr);
 Ystr(1)=Ystream(istr);

%------
for step=1:Nstep

 Xvel=Xstr(step); % starting point
 Yvel=Ystr(step);

[SXX,SXY,SYX,SYY ...
...
  ,QXXX,QXXY,QXYX,QXYY ...
  ,QYXX,QYXY,QYYX,QYYY ...
...
  ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
...
  = sgf_ax_w (Iopt ...
             ,Xvel,Yvel ...
             ,Xring,Yring ...
             ,wall);

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
  = sgf_ax_w (Iopt ...
             ,Xvel,Yvel ...
             ,Xring,Yring ...
             ,wall);

 if(dir==1)
  velx1=SXX;
  vely1=SYX;
 else
  velx1=SXY;
  vely1=SYY;
 end

 Xstr(step+1) = Xstr(step)+0.5*(velx+velx1)*Dt;
 Ystr(step+1) = Ystr(step)+0.5*(vely+vely1)*Dt;

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

plot(Xstr,Ystr,'k')
plot(Xstr,-Ystr,'k')

end
%---

