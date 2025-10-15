%-------------------------------------------
% grid_2d
%
% Division of a planar contour consisting of 
% line segments and circular arcs in boundary
% elements with corresponding shapes
%
% SYMBOLS:
% --------
%
% Isym = 0 if the element distribution
%          on a segment is not symmetric
% Isym = 1 if the element distribution on a segment 
%          is symmetric with respect to the mid-point
%---------------------------------------------------

%-----------------------------------------
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis([-2 6 -2 6])
axis square
box

%--------
% prepare
%--------
 
  file1 = fopen('grid_2d_m.dat');
  xstart = fscanf(file1,'%f',[1,1]);
  ystart = fscanf(file1,'%f',[1,1]);
  Itype = fscanf(file1,'%f',[1,1]);

  sinit = 0.0;     % origin of arc length
  Iseg = 0;        % counts the number of segments

%-------------------
% begin the gridding
%-------------------

while(Itype~=0)

    Iseg = Iseg+1;

%--------------
   if(Itype==1) % line segment
%--------------

     xend = fscanf(file1,'%f',[1,1]);
     yend = fscanf(file1,'%f',[1,1]);
     N    = fscanf(file1,'%f',[1,1]);
     ratio= fscanf(file1,'%f',[1,1]);
     Isym = fscanf(file1,'%f',[1,1]);

     x1 = xstart;
     y1 = ystart;
     x2 = xend;
     y2 = yend;

     [xe,ye,se,xm,ym,sm]...
    ...
      = elm_line (N,ratio,x1,y1,x2,y2,sinit,Isym);

%-----------------------
    elseif(Itype==2) % circular arc
%----------------------

     radius = fscanf(file1,'%f',[1,1]);
     angle1 = fscanf(file1,'%f',[1,1]);
     angle2 = fscanf(file1,'%f',[1,1]);
     N      = fscanf(file1,'%f',[1,1]);
     ratio  = fscanf(file1,'%f',[1,1]);
     Isym   = fscanf(file1,'%f',[1,1]);

     angle1 = angle1*pi;
     angle2 = angle2*pi;

     xcnt = xstart-radius*cos(angle1);
     ycnt = ystart-radius*sin(angle1);

     [xe,ye,te,se,xm,ym,tm,sm]...
      ...
      = elm_arc (N,ratio,xcnt,ycnt,radius,angle1,angle2,sinit,Isym);

%---
 end
%---

     if(Iseg==1) plot(xe,ye,'-o'); end
     if(Iseg==2) plot(xe,ye,'-og'); end
     if(Iseg==3) plot(xe,ye,'-or'); end
     if(Iseg==4) plot(xe,ye,'-oc'); end
     if(Iseg==5) plot(xe,ye,'-om'); end
     if(Iseg==6) plot(xe,ye,'-ok'); end

%---
% reset starting point
%---

  xstart = xe(N+1);
  ystart = ye(N+1);
  sinit  = se(N+1);

  Itype = fscanf(file1,'%f',[1,1]);

    clear xe ye

%------------
 end      % while:  return for another segment
%------------
