close all
clear all

%------------
% FDLIB BEMLIB
%
% draw streamlines due to a periodic array of 
% rings of of point forces
%------------


Iopt=1; % compute only the Green's function

Np = 2;
Nsum = 4;

%---
Xring=0.0;
Yring=1.0;
%---

%---
% preferences
%---

dir=2;  % radial
dir=1;  % axial

%---
figure(1)
hold on
plot([-2, 2],[0 0],'k')
set(gca,'fontsize',14)
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
axis square
axis([-2 2 -2 2])
box
%---


%---
if(dir==1)
%---

 L = 1.0; % ring separation

 Xstream =[ 0.01 0.01 0.01 0.01 ...
            0.01 0.01 0.01 0.01 ...
            0.01 0.01 ...
            0.01 0.01 0.01 0.01 ...
            0.01 0.01 0.01 0.01 ...
            0.01 0.01 0.01 0.01]';
 Ystream =[ 0.01 0.05 0.10 0.20 ...
            0.30 0.40 0.50 0.60 ...
            0.70 0.80 ...
            0.90 1.00 1.10 1.20 ...
            1.30 1.40 1.50 1.60 ...
            1.70 1.80 1.90 1.95]';
 Nstr = size(Xstream);
 Nstep = 2*128;
 Ds = 0.02;

%---
else
%---

L = 2.0; % ring separation

 Xstream = 0.5*L*[ 0.01 0.05 0.1 0.2 0.3  ...
                  0.4 ...
                  0.5 0.6 0.7 0.8 0.9 ...
                  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
                 ]';
 Ystream =[ 1.99 1.99 1.99 1.99 1.99 ...
            1.99 ...
            1.99 1.99 1.99 1.99 1.99 ...
            Yring Yring Yring Yring Yring Yring Yring  Yring Yring 
            ]';
 Nstr =size(Xstream);
 Nstep =2*128;
 Ds =-0.02;
%---
end
%---


%---
for istr=1:Nstr
%---

 clear Xstr Ystr

 Xstr(1) = Xstream(istr);
 Ystr(1) = Ystream(istr);

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
  = sgf_ax_1p (Iopt ...
              ,Xvel,Yvel ...
              ,Xring,Yring ...
              ,L ...
              ,Nsum,Np);

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
  = sgf_ax_1p (Iopt ...
              ,Xvel,Yvel ...
              ,Xring,Yring ...
              ,L ...
              ,Nsum,Np);

 if(dir==1)
  velx1=SXX;
  vely1=SYX;
 else
  velx1=SXY;
  vely1=SYY;
 end

 Xstr(step+1) = Xstr(step)+0.5*(velx+velx1)*Dt;
 Ystr(step+1) = Ystr(step)+0.5*(vely+vely1)*Dt;

 if(Xstr(step+1)>L)
  break
 end
 if(Xstr(step+1)<-L)
  break
 end
 if(Ystr(step+1)>2)
  break
 end
 if(Ystr(step+1)<-2)
  break
 end

%---
end

for irepeat=-1:2
 Xstrr = Xstr+irepeat*L;
 plot( Xstrr, Ystr,'k')
 plot( Xstrr,-Ystr,'k')
 if(dir==2)
  Xstrr = Xstr+irepeat*L;
  plot(-Xstrr, Ystr,'k')
  plot(-Xstrr,-Ystr,'k')
 end
end

%---
end
%---
plot(Xring, Yring,'ko')
plot(Xring,-Yring,'ko')

