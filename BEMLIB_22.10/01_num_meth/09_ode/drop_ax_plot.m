%========
% plot the data contained in file: drop_ax.out
%========

close all
clear all
hold on

file2 = fopen('drop_ax.out');

Nwall = fscanf(file2,'%f',[1,1]);
 idle = fscanf(file2,'%f',[1,1]);
 xwall(1) = fscanf(file2,'%f',[1,1]);
 ywall(1) = fscanf(file2,'%f',[1,1]);
  idle = fscanf(file2,'%f',[1,1]);
 xwall(2) = fscanf(file2,'%f',[1,1]);
 ywall(2) = fscanf(file2,'%f',[1,1]);

more=1;
while (more>0)
 npts   = fscanf(file2,'%f',[1,1]);
 if(npts==0) break; end
 alpha = fscanf(file2,'%f',[1,1]);
 a     = fscanf(file2,'%f',[1,1]);
 alpha=alpha*pi;
 Jsp   = fscanf(file2,'%d',[1,1]);
 for i=1:npts
  idle = fscanf(file2,'%f',[1,1]);
  x(i) = fscanf(file2,'%f',[1,1]);
  y(i) = fscanf(file2,'%f',[1,1]);
 end
 patch(x,y,'y')
 plot(x,y,'.')

 xc=a*cos(alpha)
 dtheta=alpha/64;
 for i=1:64+1
  theta=pi-alpha+(i-1)*dtheta;
  xcrc(i)=xc+a*cos(theta);
  ycrc(i)=   a*sin(theta);
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

end

xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis([-1.0 1.0 -0.5 1.0])
axis equal
box
