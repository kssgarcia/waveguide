    close all
    clear all

%------------------------------------------------------
% Hydrostatic shape of an axisymmetric
% sessile drop resting on a horizontal plate,
% or pendant drop hanging underneath a horizontal plate,
% for a specified volume and contact angle.
%------------------------------------------------------

  Jsp= 1; % sessile
  Jsp=-1; % pendant

  gac=1.0;   % acceleration of gravity
  gamma=1.8; % surface tension
  rhod=1.0;  % density of the drop
  rhoa=0.0;  % density of the ambient fluid
  volume=1.0;   % drop volume

  alpha=0.75*pi;   % contact angle
  npts =32;     %  number of interfacial markers
  epsilon=0.01; % for the shooting method

  maxiter=16;
  tol=0.000000001;

  [req,x,s] = drop_ax (Jsp ...
   ,gac,gamma,rhod,rhoa,volume ...
   ,alpha,npts,epsilon,maxiter,tol ...
   );

%---
% plot
%---

 hold on

 plot( s,Jsp*x,'-o')
 plot(-s,Jsp*x,'-o')

 xwall(1)=-1.5;
 ywall(1)= 0;
 xwall(2)= 1.5;
 ywall(2)= 0;

 dtheta=alpha/64;
 xc=req*cos(alpha);
 for i=1:64+1
  theta=pi-alpha+(i-1)*dtheta;
  xcrc(i)=xc+req*cos(theta);
  ycrc(i)=   req*sin(theta);
 end

 if(Jsp==1)
  patch([xwall xwall(2) xwall(1)] ,[ywall ywall(2)-0.2 ywall(1)-0.2],'g')
  plot( ycrc,-xcrc,'--r')
  plot(-ycrc,-xcrc,'--r')
 else
  patch([xwall xwall(2) xwall(1)] ,[ywall ywall(2)+0.2 ywall(1)+0.2],'g')
  plot( ycrc,xcrc,'--r')
  plot(-ycrc,xcrc,'--r')
 end
 plot(xwall,ywall,'k')

xlabel('y','fontsize',15)
ylabel('x','fontsize',15)
set(gca,'fontsize',15)
axis([-1.0 1.0 -0.5 1.0])
axis equal
box

%---
% done
%---
 


