close all
clear all

%====================
% draw a fractal tree
% in the xy plane
%====================

%---
h = 1.0;  % height
r = 0.45;  % branch point
b = 0.70  % branch length
theta = 0.30*pi;  % branch angle
mgen = 5;  % generation index
%---

%---
% prepare
%---

%---
cs = cos(theta);
sn = sin(theta);
rc = 1-r;
rcb = b*rc;
%---

%---
A(1,1) =  cs;
A(1,2) = -sn;
A(2,1) = -A(1,2);
A(2,2) =  A(1,1);
%---

%---
% first generation
%
% k(i)=1 indicates a branch-end
%---

n=8;
x(1) =  0.0;     y(1) = 0.0;           k(1)=0;
x(2) =  0.0;     y(2) = r*h;           k(2)=0;
x(3) =  0.0;     y(3) = h ;            k(3)=1;
x(4) =  0.0;     y(4) = r*h;           k(4)=0;
x(5) =-rcb*h*sn; y(5) = h*(r+rcb*cs);  k(5)=1;
x(6) =  0.0;     y(6) = r*h;           k(6)=0;
x(7) = rcb*h*sn; y(7) = h*(r+rcb*cs);  k(7)=1;
x(8) =  0.0;     y(8) = r*h;           k(8)=0;

%---
% run over generations
%---

for gen=1:mgen

nnew = n;
xnew = x;
ynew = y;
knew = k;

 for i=n:-1:1

  if(k(i)==1)

   for j=nnew:-1:i+1
    xnew(j+6)=xnew(j);
    ynew(j+6)=ynew(j);
    knew(j+6)=knew(j);
   end

   x1=x(i-1);y1=y(i-1);
   x2=x(i)  ;y2=y(i);
   dl = sqrt((x2-x1)^2+(y2-y1)^2);
   xnew(i) = rc*x1+r*x2;
   ynew(i) = rc*y1+r*y2;
   knew(i) = 0;
   xnew(i+1) = x2;
   ynew(i+1) = y2;
   knew(i+1) = 1;
   xnew(i+2) = xnew(i);
   ynew(i+2) = ynew(i);
   knew(i+2) = 0;
   xnew(i+3) = xnew(i)+rcb*(A(1,1)*(x2-x1)+A(1,2)*(y2-y1));
   ynew(i+3) = ynew(i)+rcb*(A(2,1)*(x2-x1)+A(2,2)*(y2-y1));
   knew(i+3) =1;
   xnew(i+4) = xnew(i);
   ynew(i+4) = ynew(i);
   knew(i+4) = 0;
   xnew(i+5) = xnew(i)+rcb*(A(1,1)*(x2-x1)+A(2,1)*(y2-y1));
   ynew(i+5) = ynew(i)+rcb*(A(1,2)*(x2-x1)+A(2,2)*(y2-y1));
   knew(i+5) = 1;
   xnew(i+6) = xnew(i);
   ynew(i+6) = ynew(i);
   knew(i+6) = 0;
   nnew = nnew+6;
  end

 end

 n = nnew;
 x = xnew;
 y = ynew;
 k = knew;

end

%---
% end of run
%---


figure(1)
axis square;
axis equal;
axis off;
plot(x,y,'k')
