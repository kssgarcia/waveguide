  function [x,f,Iflag] = newton1_2 ...
     ...
    (menu ...
    ,Niter ...
    ,eps ...
    ,x ...
    ,italk ...
    )

%==========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licencing agreement
%==========================================

%------------------------------------------------
% This program accompanies the book:
%             C. Pozrikidis
% "Numerical Computation in Science and Engineering"
%        Oxford University Press
%------------------------------------------------

%-----------------------------------
% Solve one nonlinear equation by the
% second-order Newton method
%-----------------------------------

   tol = 0.0000001;
   relax = 1.0;

%---------------------
% start the iterations
%---------------------

  for i=1:Niter

   f = newton1_2_fun(menu,x);
   x1 = x + eps;     % derivative by finite differences
   f1 = newton1_2_fun(menu,x1);
   Df = (f1-f)/eps;
   Dx = -f/Df;
   x  = x + Dx;

     if(italk==1)  
       format long;
       disp([x,f]);
       format short;
     end

   iescape = 1;
   if(abs(Dx) > tol) iescape = 0; end

   if(iescape==1)
      Iflag = 0;
      f = newton1_2_fun(menu,x);
     return
   end

  end

%---
% done
%---

 return
