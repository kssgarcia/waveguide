clear all
close all

%===============================
% Driver for sum_fast
%
% Sum an infinite series whose
% terms decay like 1/i^s
% for i=1, 2, ...
% usine Aitken extrapolation
%===============================

%---
% parameters
%---

N = 1;
p = 2;

kmax = 5;
kmax = 15;

tol = 0.00000001;
tol = 0.00001;

menu = 2;

menu = 1;
s = 1.1;

%---
% sum
%---

[Sum,a,b,kmax] = sum_fast (menu,s,N,p,kmax,tol);

%===
% verbose
%===

dsp(1,1)=b(1);

 for k=2:kmax
   dsp(k,1)=b(k);
   dsp(k,2)=a(k);
 end

 dsp(kmax+1,1) = b(kmax+1);
 dsp(kmax+1,2) = a(kmax+1);

format long
dsp
