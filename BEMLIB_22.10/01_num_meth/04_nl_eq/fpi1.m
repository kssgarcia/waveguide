clear all
close all

%-----
% graphical depiction of
% fixed-point iterations
% for one equation
%-----

figure(1)
hold on

%---
% choose the case
%---

icase=5;

if(icase==1)
 xmin=0.001; xmax=5.0;
 x=1.5;
 itmax=5;
 length= 0.075;
 axis([xmin,xmax,0,5])
elseif(icase==2)
 xmin=0.0; xmax=1.5;
 x=0.1;
 itmax=5;
 length= 0.025;
elseif(icase==3)
 xmin=-0.5; xmax=2.0;
 x=0.9;
 itmax=5;
 length= 0.075;
 axis([xmin,xmax,-2,4])
elseif(icase==4)
 xmin=0.0; xmax=2.5;
 x=1.0;
 itmax=7;
 length= 0.050;
 axis([0,2,0,2.5])
elseif(icase==5)
 xmin=0; xmax=1.5;
 x=0.4;
 itmax=3;
 length= 0.050;
 axis([xmin,1.5,-1,2.0])
 xmin=0.001; xmax=5.0;
end
%===

i=1;
xplot(i)=x; yplot(i)=0;

%----
for k=1:itmax
%----
 if(icase==1)
  g=log(3)+ 2.0*log(x);
 elseif(icase==2)
  g=exp(-x);
 elseif(icase==3)
  g=-1.8+exp(x);
 elseif(icase==4)
  g=2.0*cos(x);
 elseif(icase==5)
  g=x*(2-x);
 end
%----
 i=i+1;
 xplot(i)=x; yplot(i)=g;
 i=i+1;
 xplot(i)=g; yplot(i)=g;
 x=g;
end
plot(xplot,yplot)

%===
for i=1:2*itmax
%length = 0.10*sqrt((xplot(i)-xplot(i+1))^2+(yplot(i)-yplot(i+1))^2);
 drawarrow(xplot(i),xplot(i+1),yplot(i),yplot(i+1),length);
% plot_arrow(xplot(i),yplot(i),xplot(i+1),yplot(i+1),'linewidth',1);
end
%===

nplt=64;
for j=1:nplt+1
 xgplt(j) = (j-1)*(xmax-xmin)/64+xmin;
 if(icase==1)
 ygplt(j) = log(3)+ 2.0*log(xgplt(j));
 elseif(icase==2)
 ygplt(j) = exp(-xgplt(j));
 elseif(icase==3)
  ygplt(j) = -1.8+exp(xgplt(j));
 elseif(icase==4)
 ygplt(j) = 2.0*cos(xgplt(j));
 elseif(icase==5)
 ygplt(j) = xgplt(j)*(2-xgplt(j));
 end
end
plot(xgplt,ygplt,'--')

%===

xdiag(1)=xmin;ydiag(1)=xmin;
xdiag(2)=xmax;ydiag(2)=xmax;
plot(xdiag,ydiag)

xdiag(1)=xmin;ydiag(1)=0.0;
xdiag(2)=xmax;ydiag(2)=0.0;
plot(xdiag,ydiag)


%=====
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box on
print -deps foo.eps
%=====
