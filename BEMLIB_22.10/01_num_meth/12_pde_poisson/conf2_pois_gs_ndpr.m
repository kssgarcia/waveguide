function [f, iter, Iflag] = conf2_pois_gs_ndpr ...
...
   (Nxi,Neta,xi,eta,x,y,h,g,itermax,tol,relax,fbot,f);

%------------------------------------------
% Solution of Poisson's equation
% in orthogonal curvilinear coordinates
% with the homogeneous Neumann boundary condition
% at the top boundary,
% the Dirichlet boundary condition
% at the bottom boundary,
% and periodic boundary condition
% at the left and right boundaries
%
% The solution is found by point Gauss-Siedel
% iterations
%------------------------------------------

%--------
% compute the differentiation coefficients
%--------

[bl,bc,br, bb,bx,bt, cl,cc,cr,cb,ct] = conf2_diff  ...
...
   (Nxi,Neta,xi,eta,h);

%---
% prepare
%---

Iflag = 0;  % convergence flag; 0 is good

%---
% initialize
%---

%for i=1:Nxi+2
% for j=1:Neta+1
%  f(i,j) = 0.0;
% end
%end

%---
% implement the Dirichlet
% boundary condition at the bottom
%---

for i=1:Nxi+1
 f(i,1) = fbot(i);  % bottom
end

%------------------------
% Gauss-Siedel iterations
%------------------------

for iter=1:itermax

%---
% homogeneous Neumann condition on top
%---

for i=1:Nxi+2
  f(i,Neta+1)=f(i,Neta-1);
end

%------------------------
% update nodes row-by-row
%------------------------

%---
% save
%---

 fsave = f;

%---
% iterate
%---

 cormax=0.0;

 for j=2:Neta
  for i=2:Nxi+1
   res = -(cl(i,j)*f(i-1,j)+cr(i,j)*f(i+1,j) ...
          +cb(i,j)*f(i,j-1)+ct(i,j)*f(i,j+1) ...
          +g(i,j) )/cc(i,j) ...
          -f(i,j);

   f(i,j) = f(i,j) + relax*res;

   cor = abs(f(i,j)-fsave(i,j));
   if(cor > cormax)
     cormax = cor;
   end

  end
 end

% cormax

%------------------------------------------------
% wrap: implement the periodic boundary condition
%------------------------------------------------

 for j=1:Neta+1
   f(1,j)    = f(Nxi+1,j);
   f(Nxi+2,j) = f(2,j);
  end

 if(cormax<tol) 
   Iflag=1;
   break
 end

%---
end % of iterations
%---

%-----
% done
%-----

return
