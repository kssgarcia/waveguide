close all
clear

%===
% visualization of a Julia set
%===

xdiv = linspace(-1.25,1.25,100);   % grid
ydiv = linspace(-1.25,1.25,100);   % grid

[X,Y] = meshgrid(xdiv,ydiv);       % grid

c =-0.488679-0.5679*i;
l = 30;

Z = X+i*Y;

for k=1:l; 
  Z = Z.^2+c;
  W = exp(-abs(Z));
end

%===
% plot
%===

figure(1)
colormap copper(256)
pcolor(W);
shading flat;
axis('square','equal','off');

