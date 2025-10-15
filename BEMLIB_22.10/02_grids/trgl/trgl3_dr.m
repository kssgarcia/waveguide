close all
clear all

%============================================
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%
% Ishape = 1: octahedron
% Ishape = 2: icosahedron
%============================================

Ishape = 1;
Ishape = 2;

%----
% triangulate
%----

if(Ishape==1)
 Ndiv = 3;
 [Npts,Nelm,p,ne,n,nbe] = trgl3_octa (Ndiv);
else
 Ndiv = 2;
 [Npts,Nelm,p,ne,n,nbe] = trgl3_icos (Ndiv);
end

%----
% plot elements 
%----

figure(1)
hold on
axis equal
box on
if(Ishape==1) title('Descendant of the octahedron'); end
if(Ishape==2) title('Descendant of the icosahedron'); end

for i=1:Nelm
 j1= n(i,1);
 j2= n(i,2);
 j3= n(i,3);
 px1 = p(j1,1); px2 = p(j2,1); px3 = p(j3,1);
 py1 = p(j1,2); py2 = p(j2,2); py3 = p(j3,2);
 pz1 = p(j1,3); pz2 = p(j2,3); pz3 = p(j3,3);
 patch([px1,px2,px3],[py1,py2,py3],[pz1,pz2,pz3],'y');
end
