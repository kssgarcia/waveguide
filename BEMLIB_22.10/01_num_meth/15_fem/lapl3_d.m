close all
clear all

%========================================
% CODE lapl3_d
%
% Solution of Laplace's equation
% with the Dirichlet boundary condition
% in a disk-like domain
% using 3-node triangles
%========================================

%-----------
% input data
%-----------

ndiv = 3;  % discretization level

%------------
% triangulate
%------------

[ne,ng,p,c,efl,gfl] = trgl3_disk (ndiv);
% [ne,ng,p,c,efl,gfl] = trgl3_delaunay;

%-------
% deform
%-------

defx = 0.6;
defx = 0.0;

for i=1:ng
 p(i,1)=p(i,1)*(1.0-defx*p(i,2)^2 );
end

%-----------------------------------------
% specify the Dirichlet boundary condition
%-----------------------------------------

for i=1:ng
 if(gfl(i,1)==1)
   gfl(i,2) = sin(pi*p(i,2));    % example
   gfl(i,2) = p(i,1);            % another example
   gfl(i,2) = p(i,1)^2;            % another example
   gfl(i,2) = p(i,1)*sin(0.5*pi*p(i,2));   % another example
 end
end

%-------------------------------------
% assemble the global diffusion matrix
%-------------------------------------

gdm = zeros(ng,ng); % initialize

for l=1:ne          % loop over the elements

% compute the element diffusion matrix

j=c(l,1); x1=p(j,1); y1=p(j,2);
j=c(l,2); x2=p(j,1); y2=p(j,2);
j=c(l,3); x3=p(j,1); y3=p(j,2);

[edm_elm] = edm3 (x1,y1,x2,y2,x3,y3);

   for i=1:3
     i1 = c(l,i);
     for j=1:3
       j1 = c(l,j);
       gdm(i1,j1) = gdm(i1,j1) + edm_elm(i,j);
     end
   end
end

% disp (gdm);

%---------------------------------------------
% set the right-hand side of the linear system
% and implement the Dirichlet boundry condition
%----------------------------------------------

for i=1:ng
 b(i) = 0.0;
end

for j=1:ng
 if(gfl(j,1)==1) 
   for i=1:ng 
    b(i) = b(i) - gdm(i,j)*gfl(j,2);
    gdm(i,j)=0; gdm(j,i)=0;
   end
   gdm(j,j)=1.0;
   b(j)=gfl(j,2);
 end
end

%------------------------
% solve the linear system
%------------------------

f = b/gdm'; 

%-----
% plot
%-----

plot_3 (ne,ng,p,c,f);
trimesh (c,p(:,1),p(:,2),f);
%trisurf (c,p(:,1),p(:,2),f,f);

%-----
% done
%-----
