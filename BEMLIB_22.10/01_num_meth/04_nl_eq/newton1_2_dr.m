  clear all
  close all

%===========
% FDLIB
%
% Driver for one nonlinear equation
% using the second-order
% Newton's method 
%===========

  global rlp
  global mexp alpha beta

%---
% paremeters
%---

  Niter = 10;
  eps = 0.0001;
  italk = 1;

  menu = 110;
  menu = 111;
  menu = 200;

%---
% prepare
%---

  if(menu==110)
   rlp = 0.0078974814;
   x = 0.3482780010;
   x = 0.3;
  elseif(menu==111)
   mexp=0.5;
   alpha = 1.0;
   beta = -1.0;
   x = beta/(2.0*alpha);
  elseif(menu==200)
   x = 0.4;
  else
   disp ('this choice is not available')
   break
  end

  [x,f,Iflag] = newton1_2 ...
     ...
    (menu ...
    ,Niter ...
    ,eps ...
    ,x ...
    ,italk ...
    );

