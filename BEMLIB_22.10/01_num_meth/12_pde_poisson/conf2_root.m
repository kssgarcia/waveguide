function [x, y] = conf2_root(alpha,k,xi,eta,x,y)

%-------------------------------------------------
% node position in the xy plane by Newton's method
%-------------------------------------------------

 eps=0.0001; %small number for computing the Jacobian
 tol=0.000000001; % tolerance
 Niter = 20;  % max number of iterations

%-----------
% iterations
%-----------

  for iter=1:Niter

    f = conf2_fnc (alpha,k,xi,eta,x,y);
    error = sqrt(f(1)^2+f(2)^2);
    if(error<tol) break; end;

%---------------------
% Jacobian by numerical differentiation
%---------------------

     x = x+eps;      % perturb
     f1 = conf2_fnc(alpha,k,xi,eta,x,y);
     x = x-eps;       % reset
     Jac(1,1) = (f1(1)-f(1))/eps;
     Jac(2,1) = (f1(2)-f(2))/eps;

     y = y + eps;      % perturb
     f1 = conf2_fnc (alpha,k,xi,eta,x,y);
     y = y - eps;       % reset
     Jac(1,2) = (f1(1)-f(1))/eps;
     Jac(2,2) = (f1(2)-f(2))/eps;

     b1  = -f(1);
     b2  = -f(2);
     Det = Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1);
     dx = (b1*Jac(2,2)-Jac(1,2)*b2)/Det;
     dy = (b2*Jac(1,1)-Jac(2,1)*b1)/Det;
     x = x+dx;
     y = y+dy;

 end

%-----
% done
%-----

return
