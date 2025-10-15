close all
clear all

%========================================
% Generate and display a buckyball
% (Buckminsterfullerene)
% 60 vertices
%
% 32 faces consisting of:
% 12 pentagons
% 20 hexagons
%
% Each vertex has three closest neighbors
%========================================

f = (1+sqrt(5))/2.0;  % golden ratio

Nv = 60;  % number of vertices

%---
% hard-code the vertices on the unit sphere
%---

fc = 1/sqrt(1+9*f^2);
f0 = f+2;
f1 = 1+2*f;
f2 = 2*f;
f3 = 3*f;

v( 1,:) = fc*[  0,  1, f3];
v( 2,:) = fc*[  0,  1,-f3];
v( 3,:) = fc*[  0, -1, f3];
v( 4,:) = fc*[  0, -1,-f3];
v( 5,:) = fc*[  1, f3,  0];
v( 6,:) = fc*[  1,-f3,  0];
v( 7,:) = fc*[ -1, f3,  0];
v( 8,:) = fc*[ -1,-f3,  0];
v( 9,:) = fc*[ f3,  0,  1];
v(10,:) = fc*[ f3,  0, -1];
v(11,:) = fc*[-f3,  0,  1];
v(12,:) = fc*[-f3,  0, -1];

v(13,:) = fc*[ 2, f1, f];
v(14,:) = fc*[ 2, f1,-f];
v(15,:) = fc*[ 2,-f1, f];
v(16,:) = fc*[ 2,-f1,-f];
v(17,:) = fc*[-2, f1, f];
v(18,:) = fc*[-2, f1,-f];
v(19,:) = fc*[-2,-f1, f];
v(20,:) = fc*[-2,-f1,-f];

v(21,:) = fc*[ f1, f, 2];
v(22,:) = fc*[ f1, f,-2];
v(23,:) = fc*[ f1,-f, 2];
v(24,:) = fc*[ f1,-f,-2];
v(25,:) = fc*[-f1, f, 2];
v(26,:) = fc*[-f1, f,-2];
v(27,:) = fc*[-f1,-f, 2];
v(28,:) = fc*[-f1,-f,-2];

v(29,:) = fc*[ f, 2, f1];
v(30,:) = fc*[ f, 2,-f1];
v(31,:) = fc*[ f,-2, f1];
v(32,:) = fc*[ f,-2,-f1];
v(33,:) = fc*[-f, 2, f1];
v(34,:) = fc*[-f, 2,-f1];
v(35,:) = fc*[-f,-2, f1];
v(36,:) = fc*[-f,-2,-f1];

v(37,:) = fc*[ 1, f0, f2];
v(38,:) = fc*[ 1, f0,-f2];
v(39,:) = fc*[ 1,-f0, f2];
v(40,:) = fc*[ 1,-f0,-f2];
v(41,:) = fc*[-1, f0, f2];
v(42,:) = fc*[-1, f0,-f2];
v(43,:) = fc*[-1,-f0, f2];
v(44,:) = fc*[-1,-f0,-f2];

v(45,:) = fc*[ f0, f2, 1];
v(46,:) = fc*[ f0, f2,-1];
v(47,:) = fc*[ f0,-f2, 1];
v(48,:) = fc*[ f0,-f2,-1];
v(49,:) = fc*[-f0, f2, 1];
v(50,:) = fc*[-f0, f2,-1];
v(51,:) = fc*[-f0,-f2, 1];
v(52,:) = fc*[-f0,-f2,-1];

v(53,:) = fc*[ f2, 1, f0];
v(54,:) = fc*[ f2, 1,-f0];
v(55,:) = fc*[ f2,-1, f0];
v(56,:) = fc*[ f2,-1,-f0];
v(57,:) = fc*[-f2, 1, f0];
v(58,:) = fc*[-f2, 1,-f0];
v(59,:) = fc*[-f2,-1, f0];
v(60,:) = fc*[-f2,-1,-f0];

% radii should be 1.0:

for i=1:Nv
 radius(i) = v(i,1)^2+v(i,2)^2+v(i,3)^2;
end

%======================================================
% connectivity table
%
% each vertex has three neighbors
%
% c(i,j) = 1 if ith and jth vertices are neighbors
%        = 0 if ith and jth vertices are not neighbors
%
% neigh(i,1) is the first  neighbor of the ith vertex
% neigh(i,2) is the second neighbor of the ith vertex
% neigh(i,3) is the third  neighbor of the ith vertex
%======================================================

c = zeros(Nv,Nv);

% run over all vertices:

for i=1:Nv

% first neighbor:

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

% second neighbor:

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

% third neighbor:

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

%---------------------------------
% connectivity matrix is symmetric
%---------------------------------

for i=1:Nv
 for j=1:i-1
  c(j,i)=c(i,j);
 end
end

%==========================================
% identify the vertices composing the faces
%
% facetype = 5 for a pentagon 
% facetype = 6 for a hexagon
%==========================================

for i=1:Nv  % vertex tag
 vtag(i) = 0;
end

nface=0;  % count faces

%---
% Each vertex belongs to three polygonal faces
% Scan the vertices three times to pick up all faces
%---

for repeat=1:3  

for i=1:Nv  % run over vertices

 nface=nface+1;

%---
% find the face vertices
%---

if(repeat==1)   % distances from 1,2 neighbors
 X1 = v(i,1)-v(neigh(i,1),1);
 Y1 = v(i,2)-v(neigh(i,1),2);
 Z1 = v(i,3)-v(neigh(i,1),3);
 X2 = v(i,1)-v(neigh(i,2),1);
 Y2 = v(i,2)-v(neigh(i,2),2);
 Z2 = v(i,3)-v(neigh(i,2),3);

elseif(repeat==2)  % distances from 2,3 neighbors
 X1 = v(i,1)-v(neigh(i,2),1);
 Y1 = v(i,2)-v(neigh(i,2),2);
 Z1 = v(i,3)-v(neigh(i,2),3);
 X2 = v(i,1)-v(neigh(i,3),1);
 Y2 = v(i,2)-v(neigh(i,3),2);
 Z2 = v(i,3)-v(neigh(i,3),3);

else  % distances from 3,1 neighbors
 X1 = v(i,1)-v(neigh(i,3),1);
 Y1 = v(i,2)-v(neigh(i,3),2);
 Z1 = v(i,3)-v(neigh(i,3),3);
 X2 = v(i,1)-v(neigh(i,1),1);
 Y2 = v(i,2)-v(neigh(i,1),2);
 Z2 = v(i,3)-v(neigh(i,1),3);
end

 prx = Y1*Z2-Z1*Y2;   % cross product
 pry = Z1*X2-X1*Z2;   % is normal to the face
 prz = X1*Y2-Y1*X2;

 Ic=0;

%---
% run over all vertices
% if inner product is zero,
% the vertex is accepted
%
% face(nface,j) contains the vertex labels
%---

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

 facetype(nface)=Ic;   % set the face type, 5 or 6

%---
% rearrange the face vertices in cyclic order
% with the help of the connectivity table
% so that all vertices are neighbors
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

end % of run over vertices
end % of repeat

%-------------
% unique faces
%
% nfaceu: number of unique vertices
%  faceu: vertices of unique faces
%-------------

% accept the first face:

nfaceu=1;
facetypeu(1) = facetype(1);
for k=1:facetype(1)
  faceu(1,k)=face(1,k);
end

% run over all other faces:

%------------
for i=2:nface
%------------

% comparison measure is count:

 count = 0;
 for k=1:facetype(i)
  count = count+face(i,k)^2;
 end

% compare with current unique faces:

 flag=0;
 for j=1:nfaceu
   count1 = 0;
   for k=1:facetypeu(j)
    count1 = count1+faceu(j,k)^2;
   end
   if(abs(count-count1)<0.0001) flag=1; end
 end

% accept a new face:

 if(flag==0)
  nfaceu=nfaceu+1;
  facetypeu(nfaceu) = facetype(i);
  for k=1:facetype(i)
   faceu(nfaceu,k)=face(i,k);
  end
 end

%------------
end
%------------

%============ GRPAPHICS ======================

%h1 =  surf(x1,y1,z1,w1);
%h5 = plot3(xplot,yplot,zplot,'b','EraseMode','xor','Linewidth',lw);
%h6 = plot3(xplot,yplot,zplot,'r','EraseMode','xor','Linewidth',2);

%----
% a sphere
%----

[x,y,z] = sphere(16);

rad=0.05;
x = rad*x; y=rad*y; z=rad*z;

%----
% draw spheres at the vertices
%----

for i=1:Nv
 x1 = x + v(i,1);
 y1 = y + v(i,2);
 z1 = z + v(i,3);
 c1 = ones(size(x1));
 surf(x1,y1,z1,c1);
 colormap(jet)
 f = findobj('Type','surface');
 set(f,'FaceLighting','phong');
 material shiny
 shading interp
% light
 hold on
end

%----
% draw struts between neighbors
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
   plot3(xplot,yplot,zplot,'b','Linewidth',5)
   plot3(xplot,yplot,zplot,'r','Linewidth',2)
  end
 end
end

%----------------------
% paint the unique faces
%----------------------

for i=1:nfaceu

 for k=1:facetypeu(i)
  polyx(k) = v(faceu(i,k),1); 
  polyy(k) = v(faceu(i,k),2); 
  polyz(k) = v(faceu(i,k),3); 
  polyc(k) = 10.0; 
  if(facetypeu(i)==5) polyc(k)=5.00; end
 end
 patch(polyx,polyy,polyz,polyc)
 hold on
 clear polyx;clear polyy;clear polyz;clear polyc;

end

%-------
% finish
%-------

% title('60-vertex fullerene')
axis([-1 1 -1 1 -1 1])
axis('square')
grid on
xlabel('x','fontsize',13)
ylabel('y','fontsize',13)
zlabel('z','fontsize',13)
set(gca,'fontsize',13)
box on

%---
% done
%---
