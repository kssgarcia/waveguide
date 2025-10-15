      function c = leverrier(N,A)

%======================================
% FDLIB
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%======================================

%------------------------------------------------
% This program accompanies the book:
%          C. Pozrikidis
% Numerical Computation in Science and Engineering
%        Oxford University Press
%------------------------------------------------

%-----------------------------------------
%  Coefficients of the 
%  characteristic polynomial of an nxn matrix 
%  by Leverrier' method
%
%  Algorithm (5.4.3)
%
%  SYMBOLS:
%  --------
%       
%  N .... size (rows/columns) of matrix A
%  A .... square matrix
%  c .... vector of polynomial coefficients
%
%  B .... iterated matrix for computation of c
%  D .... iterated matrix for computation of c
%
% P_N = (-1)^N * (lambda^N + c_1*lambda^(N-1)
%     + ... + c_N
%-----------------------------------------

%-----------------------------
% initial elements of matrix B
%-----------------------------

  for i=1:N
    for j=1:N
      B(i,j) = 0.0;
    end
    B(i,i) = 1.0;
  end

%---------------------------------------
% outer loop for polynomial coefficients
%---------------------------------------

 for i=1:N   
   
   c(i) = 0.0;

%---
% update the matrix D
%---

  for j=1:N
   for k=1:N
    D(j,k) = 0.0;
     for l=1:N
      D(j,k) = D(j,k) + A(j,l)*B(l,k);
     end
   end
  end

%---
% compute c(i)
%---
 
  trace = 0.0;
 
  for j=1:N
   trace = trace + D(j,j);
  end

  c(i) = -trace/i;

%---
% update the matrix B
%---

   for j=1:N
     for k=1:N
       B(j,k) = D(j,k);
     end
     B(j,j) = B(j,j) + c(i);
   end

  end    % of outer loop

%----------------------------------------------
% update coefficients for leading factor (-1)^N
%----------------------------------------------

  factor = (-1.0)^N;

  for i=1:N
     c(i) = factor*c(i);
  end

%-----
% done
%-----

  return
