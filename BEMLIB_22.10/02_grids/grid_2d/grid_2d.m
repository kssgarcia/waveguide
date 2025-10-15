%-------------------------------------------
% grid_2d
%
% Discretization of a planar contour consisting of 
% line segments and circular arcs in boundary
% elements with corresponding shapes
%
% SYMBOLS:
% --------
%
% Xeg(i,j) is the x coordinate of the ith point
%          of the jth segment
% Yeg(i,j) is the y coordinate of the ith point
%          of the jth segment
%
% Isym = 0 if the element distribution
%          on a segment is not symmetric
% Isym = 1 if the element distribution on a segment 
%          is symmetric with respect to the mid-point
%
% If the jth segment is an arc, then
%        Teg(i,j) is the polar angle of the ith point subtended
%        from the arc center
%
% seg(i,j) is the cummulative arc length  at the ith point  of the
%          jth segment
%---------------------------------------------------

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

hold on

%--------
% prepare
%--------
 
  file1 = fopen('grid_2d_m.dat');
  xstart = fscanf(file1,'%f',[1,1]);
  ystart = fscanf(file1,'%f',[1,1]);
  Itype = fscanf(file1,'%f',[1,1]);

  sinit = 0.0;     % origin of arc length
  Iseg = 1;        % count the number of segments

%---------------
% start gridding
%---------------

while(Itype~=0)

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

     clear xe
     clear ye

     [xe,ye,se,xm,ym,sm]...
    ...
      = elm_line (N,ratio,x1,y1,x2,y2,sinit,Isym);

%-----------------------
    elseif(Itype==2) % circular arc
%----------------------

     radius = fscanf(file1,'%f',[1,1]);
     angle1 = fscanf(file1,'%f',[1,1]);
     angle2 = fscanf(file1,'%f',[1,1]);
     N    = fscanf(file1,'%f',[1,1]);
     ratio= fscanf(file1,'%f',[1,1]);
     Isym = fscanf(file1,'%f',[1,1]);

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

     figure(1)
     if(Iseg==1) plot(xe,ye,'-o'); end
     if(Iseg==2) plot(xe,ye,'-og'); end
     if(Iseg==3) plot(xe,ye,'-or'); end
     if(Iseg==4) plot(xe,ye,'-oc'); end
     if(Iseg==5) plot(xe,ye,'-om'); end
     if(Iseg==6) plot(xe,ye,'-ok'); end

%------------------------
% save into global arrays
%------------------------

 for i=1:N+1
   xeg(i,Iseg) = xe(i);   % end points
   yeg(i,Iseg) = ye(i);
   seg(i,Iseg) = se(i);
   if(Itype==2) % circular arc
   teg(i,Iseg) = te(i);
   end
 end


 for i=1:N
   xmg(i,Iseg) = xm(i);   % mid points
   ymg(i,Iseg) = ym(i);
   smg(i,Iseg) = sm(i);
   if(Itype==2) % circular arc
   tmg(i,Iseg) = tm(i);
   end
 end

  xstart = xe(N+1);   % reset starting point
  ystart = ye(N+1);
  sinit  = se(N+1);

  Iseg = Iseg+1;

  Itype = fscanf(file1,'%f',[1,1]);

%------------
 end      % while:  return for another segment
%------------

  Iseg = Iseg-1;

%-----
% done
%-----

xlabel('x','fontsize',13);
ylabel('y','fontsize',13);
set(gca,'fontsize',15)
axis equal
axis([-2 6 -1 5])
box on

