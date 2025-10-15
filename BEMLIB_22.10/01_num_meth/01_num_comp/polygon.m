clc;
clear all;
close all;

%==============
% does a point lie inside a polygon?
%==============

tol = 0.0001;  % tolerance of total sum

%------
vertices = ...
[10, 10;
 6, 8;
 7, 8;
 6, 6;
 8, 2;
 10, 4;
 9, 5;
 11, 7;
 10, 10];

x0=2; y0=6;
x0=8; y0=6;
%------

clear vertices

%------
vertices =  ...
[0  1  
 2  0 ;
 4  3 ;
 8  2 ;
 7  4 ;
 3  7 ;
 2  4 ;
 1  3 ;
 1  1 ;
 0  1 ;
];

x0=7.0; y0 = 6;
x0=2.1; y0 = 2;
%------

x = vertices(:,1);          % x-coordinates of polygon vertices.
y = vertices(:,2);          % y-coordinates of polygon vertices.

%---
% plot
%---

figure
hold on
grid on;
plot(x0, y0, 'r*', x, y);
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
axis equal
set(gca,'fontsize',15)
box
%---

theta = 0;

M = length(x)-1;

for i=1:M

 X1 = x(i)-x0;
 Y1 = y(i)-y0;
 X2 = x(i+1)-x0;
 Y2 = y(i+1)-y0;
 norm1 = sqrt(X1*X1+Y1*Y1);
 norm2 = sqrt(X2*X2+Y2*Y2);
 crossprd = X1*Y2-X2*Y1;
 den = norm1*norm2;
 dth = asin(crossprd/den);
 innerprd = X1*X2+Y1*Y2;

 if(innerprd<0)
  dth = dth*(pi/abs(dth)-1.0);
 end

 theta = theta + dth;

end

theta_over_2pi = theta/(2*pi)

if(abs(theta)<tol)
    fprintf('Point lies outside polygon \n\n');
else
    fprintf('Point lies inside polygon \n\n');
end
