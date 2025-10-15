function e = mg_gs(Niter,Nx,Ny,e,r)

%---
% perform Niter Gauss-Siedel iterations
% on the equation A e= r
%---

for i=1:Niter

 for j=2:Ny
  for i=2:Nx
   e(i,j)= 0.25*(e(i+1,j)+e(i-1,j)+e(i,j+1)+e(i,j-1)+r(i,j));
  end
 end

end

%-----
% done
%-----

return
