function [sigma] = prony (Np,f,Dt,K)

%====================================
% FDLIB
%
% Prony fitting of a signal consisting
% of a sum of complex exponentials
%
% Np; number of data
%====================================

%--------
% prepare
%--------

N = 2*K;     % number of complex exponentials
M = Np-N-1;  % number of equations 

%---------
% warning
%---------

if(M<N) 
  disp 'sorry, not enough data in the time series'
  return
end

%--------------------------------------
% form the over-determined Prony system
%--------------------------------------

for i=1:M
  for j=1:N
    am(i,j) = f(i+j-1);
  end
  rhs(i) = -f(i+N);
end

%----------------------------
% subtract adjacent equations
%----------------------------

for i=1:M-1
   for j=1:N
    amr(i,j) = am(i,j) -am(i+1,j);
   end
   rhsr(i) = rhs(i)-rhs(i+1);
end

%------------------------
% solve the normal system
%
% matlab will do this for you
%------------------------

sol = amr\rhsr';

%---
% compute the polynomial coefficients
%---

poly_cf(1) = 1.0;

for j=1:N
  poly_cf(j+1)= sol(N+1-j);
end

z = roots(poly_cf);

%---------------------
% recover sigmas from: z = exp(-i*sigma*Dt)
%---------------------

for j=1:N
  sigma(j,1) = -angle(z(j))/Dt;          % real part
  sigma(j,2) = 0.5*log(z(j)*z(j)')/Dt;   % imaginary part
end

%----------------------------------------------
%  Form the over-determined system for the
%  coefficients of the exponentials
%
%  Coefficients correspond to sigma(1), sigma(3), ...
%
%  The coefficients corresponding to sigma(2), sigma(4), ...
%  arise by switching the sign of the imaginary part
%----------------------------------------------

%      Do i=1,Np

%       time = Dt*(i-1.0D0)
%       Ic = 0

%       Do j=1,Nexp2,2
%        expp  = Dexp(sigma(j,2)*time)
%        angle =      sigma(j,1)*time
%        Ic = Ic+1
%        am(i,ic) = expp*Dcos(angle)
%        Ic = Ic+1
%        am(i,Ic) = expp*Dsin(angle)
%       End Do
%       rhs(i) = signal(i)

%      End Do

%-------------------------------------------
% Form the normal system of 2*Nexp equations
%-------------------------------------------

%      Do i=1,Nexp2

%        Do j=1,Nexp2
%         am_n(i,j) = 0.0D0
%         Do l=1,Npoints
%          am_n(i,j) = am_n(i,j)+am(l,i)*am(l,j)
%         End Do
%        End Do

%        rhs_n(i) = 0.0D0
%        Do l=1,Npoints
%         rhs_n(i) = rhs_n(i)+am(l,i)*rhs(l)
%        End Do

%      End Do

%--------------------------------------------------
% Solve the normal system by Cholesky decomposition
%--------------------------------------------------

%      Isee = 0   ! do not display

%      call chol_c
%     +
%     +  (Nexp2
%     +  ,am_n
%     +  ,rhs_n
%     +  ,Isee
%     +  ,sol
%     +  )

%--------------------
% Assign the solution
%--------------------

%     Ic = 0  ! counter

%     Do i=1,Nexp
%      j1 = 2*i-1
%      j2 = 2*i
%      ic = Ic+1
%      coeff(j1,1) = sol(ic)
%      Ic = Ic+1
%      coeff(j1,2) = sol(Ic)
%      coeff(j2,1) =  coeff(j1,1)
%      coeff(j2,2) = -coeff(j1,2)
%     End Do

%-----
% Done
%-----

return
