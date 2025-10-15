function [f] = bisectionfun(x)

%============================
% Function: bisection_fun
%
% Inputs:   x
%
% Outputs:  f - the function value
%               at x
%%%
% Note the dot operator is used
% here so the function can handle
% an array as an input. It will operate
% on each element of the array in turn
%%%
%
%
% Mark Blyth 2007
%============================

%f = x.^5-3*x.^4+2*x.^2-4;
 f = x.^2-2*x.+0.9;

return
