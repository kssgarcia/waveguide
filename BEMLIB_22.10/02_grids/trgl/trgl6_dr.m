close all
clear all

%==================================
% Driver for triangulating a sphere
% into a surface mesh of 6-node elements
% originating from an icosahedron
%==================================

%---
% refinement level
%---

ndiv = 1;
ndiv = 0;

%---
% triangulate
% to generate an icosahedron
%---

[Npts,Nelm,p,ne,n,nbe] = trgl6_icos(ndiv);

%----
% plot elements subdivided into four
% three-node triangles
%----

figure(1)
hold on
axis square
box on
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
hold on
axis square
box on
xlabel('x')
ylabel('y')
zlabel('z')

for i=1:Nelm
 j1 = n(i,1);
 j2 = n(i,2);
 j3 = n(i,3);
 j4 = n(i,4);
 j5 = n(i,5);
 j6 = n(i,6);
 px1 = p(j1,1); px2 = p(j2,1); px3 = p(j3,1);
 px4 = p(j4,1); px5 = p(j5,1); px6 = p(j6,1); 
 py1 = p(j1,2); py2 = p(j2,2); py3 = p(j3,2);
 py4 = p(j4,2); py5 = p(j5,2); py6 = p(j6,2); 
 pz1 = p(j1,3); pz2 = p(j2,3); pz3 = p(j3,3);
 pz4 = p(j4,3); pz5 = p(j5,3); pz6 = p(j6,3); 
 figure(1)
%if(pz1>=-0.01 & pz4>=-0.01 & pz6>=-0.01)
  patch([px1,px4,px6],[py1,py4,py6],[pz1,pz4,pz6],'y')
%end
%if(pz2>=-0.01 & pz5>=-0.01 & pz4>=-0.01)
  patch([px2,px5,px4],[py2,py5,py4],[pz2,pz5,pz4],'r')
%end
%if(pz3>=-0.01 & pz6>=-0.01 & pz5>=-0.01)
  patch([px3,px6,px5],[py3,py6,py5],[pz3,pz6,pz5],'w')
%end
%if(pz4>=-0.01 & pz5>=-0.01 & pz6>=-0.01)
  patch([px4,px5,px6],[py4,py5,py6],[pz4,pz5,pz6],'c')
%end
 figure(2)
 patch([py1,py4,py6],[pz1,pz4,pz6],[px1,px4,px6],'w')
 patch([py2,py5,py4],[pz2,pz5,pz4],[px2,px5,px4],'r')
 patch([py3,py6,py5],[pz3,pz6,pz5],[px3,px6,px5],'y')
 patch([py4,py5,py6],[pz4,pz5,pz6],[px4,px5,px6],'c')
 pause(0.1)
end
