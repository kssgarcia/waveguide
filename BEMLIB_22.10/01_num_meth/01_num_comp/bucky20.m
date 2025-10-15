%==========
% dodecahedron
%
% 20 vertices
% 12 pentagonal faces
% 30 edges
%
% each vertex has 3 closest neighbors
%==========

f=(1+sqrt(5))/2.0;  % golden ratio

Nv=20;  % number of vertices

%---
% define the vertices
%---

fc = 1/sqrt(3);

v(1,:) = fc*[ 1, 1, 1];
v(2,:) = fc*[-1, 1, 1];
v(3,:) = fc*[ 1,-1, 1];
v(4,:) = fc*[-1,-1, 1];
v(5,:) = fc*[ 1, 1,-1];
v(6,:) = fc*[-1, 1 -1];
v(7,:) = fc*[ 1,-1,-1];
v(8,:) = fc*[-1,-1 -1];

v(9,:)  = fc*[0, 1/f, f];
v(10,:) = fc*[0,-1/f, f];
v(11,:) = fc*[0, 1/f,-f];
v(12,:) = fc*[0,-1/f,-f];

v(13,:) = fc*[ f, 0, 1/f];
v(14,:) = fc*[-f, 0, 1/f];
v(15,:) = fc*[ f, 0,-1/f];
v(16,:) = fc*[-f, 0,-1/f];

v(17,:) = fc*[ 1/f, f, 0];
v(18,:) = fc*[-1/f, f, 0];
v(19,:) = fc*[ 1/f,-f, 0];
v(20,:) = fc*[-1/f,-f, 0];

% radii should be 1.0:

for i=1:Nv
 radius(i) = v(i,1)^2+v(i,2)^2+v(i,3)^2;
end

%======
% connectivity table
%
% c(i,j) = 1 if ith and jth vertices are neighbors
%
% neigh(i,1) is the first  neighbor of the ith vertex
% neigh(i,2) is the second neighbor of the ith vertex
% neigh(i,3) is the third  neighbor of the ith vertex
%======

c=zeros(Nv,Nv);

for i=1:Nv

 min = 100.0;
 imin1 = 1;
 for j=1:Nv
  if(j~=i)
   dist2 = (v(i,1)-v(j,1))^2+(v(i,2)-v(j,2))^2+(v(i,3)-v(j,3))^2;
   if(dist2<min)
    min = dist2;
    imin1=j;
   end
  end
 end
 c(i,imin1)=1;
 neigh(i,1)=imin1;

 min = 100.0;
 imin2 = 1;
 for j=1:Nv
  if(j~=i & j~=imin1)
   dist2 = (v(i,1)-v(j,1))^2+(v(i,2)-v(j,2))^2+(v(i,3)-v(j,3))^2;
   if(dist2<min)
    min = dist2;
    imin2=j;
   end
  end
 end
 c(i,imin2)=1;
 neigh(i,2)=imin2;

 min = 100.0;
 imin3 = 1;
 for j=1:Nv
  if(j~=i & j~=imin1 & j~=imin2)
   dist2 = (v(i,1)-v(j,1))^2+(v(i,2)-v(j,2))^2+(v(i,3)-v(j,3))^2;
   if(dist2<min)
    min = dist2;
    imin3=j;
   end
  end
 end
 c(i,imin3)=1;
 neigh(i,3)=imin3;

end

%---
% connectivity matrix is symmetric
%---

for i=1:Nv
 for j=1:i-1
  c(j,i)=c(i,j);
 end
end

%======
% identify the faces
%======

for i=1:Nv  % vertex tag
 vtag(i) = 0;
end

nface=0;  % count faces

for repeat=1:3

for i=1:Nv

 nface=nface+1;

%---
% find the face vertices
%---

if(repeat==1)
 X1= v(i,1)-v(neigh(i,1),1);
 Y1= v(i,2)-v(neigh(i,1),2);
 Z1= v(i,3)-v(neigh(i,1),3);
 X2= v(i,1)-v(neigh(i,2),1);
 Y2= v(i,2)-v(neigh(i,2),2);
 Z2= v(i,3)-v(neigh(i,2),3);
elseif(repeat==2)
 X1= v(i,1)-v(neigh(i,2),1);
 Y1= v(i,2)-v(neigh(i,2),2);
 Z1= v(i,3)-v(neigh(i,2),3);
 X2= v(i,1)-v(neigh(i,3),1);
 Y2= v(i,2)-v(neigh(i,3),2);
 Z2= v(i,3)-v(neigh(i,3),3);
else
 X1= v(i,1)-v(neigh(i,3),1);
 Y1= v(i,2)-v(neigh(i,3),2);
 Z1= v(i,3)-v(neigh(i,3),3);
 X2= v(i,1)-v(neigh(i,1),1);
 Y2= v(i,2)-v(neigh(i,1),2);
 Z2= v(i,3)-v(neigh(i,1),3);
end

 prx = Y1*Z2-Z1*Y2;
 pry = Z1*X2-X1*Z2;
 prz = X1*Y2-Y1*X2;

 Ic=0;

 for j=1:Nv
  x1 = v(j,1)-v(i,1);
  y1 = v(j,2)-v(i,2);
  z1 = v(j,3)-v(i,3);
  prj = prx*x1+pry*y1+prz*z1;
  if(abs(prj)<0.00001)
   Ic=Ic+1;
   face(nface,Ic)=j;
   vtag(j)=vtag(j)+1;
  end
 end

%---
% rearrange the vertices
% in cyclic order
%---

 for j=1:Ic-2
  if(c(face(nface,j),face(nface,j+1))==0)
    for k=j+2:Ic
      temp=face(nface,j+1);
      if(c(face(nface,j),face(nface,k))==1)
       face(nface,j+1)=face(nface,k);
       face(nface,k)=temp;
       break
      end
    end
  end
 end

end
end

%-------------
% unique faces
%-------------

nfaceu=1;
for k=1:5
  faceu(1,k)=face(1,k);
end

for i=2:nface

 count = 0;
 for k=1:5
  count = count+face(i,k)^2;
 end

 flag=0;
 for j=1:nfaceu
   count1 = 0;
   for k=1:5
    count1 = count1+faceu(j,k)^2;
   end
   if(abs(count-count1)<0.0001) flag=1; end
 end

 if(flag==0)
  nfaceu=nfaceu+1;
  for k=1:5
   faceu(nfaceu,k)=face(i,k);
  end
 end

end

%============ PLOTTING ======================

%h1 =  surf(x1,y1,z1,w1);
%h5 = plot3(xplot,yplot,zplot,'b','EraseMode','xor','Linewidth',lw);
%h6 = plot3(xplot,yplot,zplot,'r','EraseMode','xor','Linewidth',2);

%----
% print spheres at the vertices
%----

[x,y,z]=sphere(16);

rad=0.05;
x = rad*x; y=rad*y; z=rad*z;
u = ones(size(x));
for i=1:Nv
 x1 = x + v(i,1)*u;
 y1 = y + v(i,2)*u;
 z1 = z + v(i,3)*u;
 w1 =ones(size(x1));
 surf(x1,y1,z1,w1);
 colormap(jet)
 f = findobj('Type','surface');
 set(f,'FaceLighting','phong');
 material shiny
 shading interp
% light
 hold on
end
axis([-1 1 -1 1 -1 1])
axis('square')
grid off

%----
% plot struts between neighbors
%----

for i=1:Nv
for j=1:Nv
 if(c(i,j)==1)
  xplot(1)= v(i,1);
  yplot(1)= v(i,2);
  zplot(1)= v(i,3);
  xplot(2)= v(j,1);
  yplot(2)= v(j,2);
  zplot(2)= v(j,3);
%  set(h5,'XData',xplot,'YData',yplot,'ZData',zplot)
%  set(h6,'XData',xplot,'YData',yplot,'ZData',zplot)
%  drawnow  % drawnow is need on the under MS-Windows
%  pause(0.1)
   plot3(xplot,yplot,zplot,'b','Linewidth',5)
   hold on
   plot3(xplot,yplot,zplot,'r','Linewidth',2)
 end
end
end

%----------------------
% draw the unique faces
%----------------------

for i=1:nfaceu

 for k=1:5
  polyx(k) =v(faceu(i,k),1); 
  polyy(k) =v(faceu(i,k),2); 
  polyz(k) =v(faceu(i,k),3); 
  polyc(k) =10.0; 
 end
 patch(polyx,polyy,polyz,polyz)
 hold on
 clear polyx;clear polyy;clear polyz;clear polyc;

end

%---
% done
%---
