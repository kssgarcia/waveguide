file1 = fopen('susp_3d.trj')
orbit  = fscanf(file1,'%f',[11,inf]);
                                                                                
si = size(orbit);
npts=si(2);
                                                                                
% skip = 27;
 skip = 5;
Icount = skip;

xmin=8.5; xmax=12.5; ymin=-2; ymax=2; zmin=-2; zmax=2;
xmin=0; xmax=30; ymin=-15; ymax=15; zmin=-15; zmax=15;
xmin=-1.5; xmax=2.5; ymin=-2; ymax=2; zmin=-1; zmax=3;

%---
% wall
%---

xp(1) = xmin;
yp(1) =  0;
zp(1) = -2;
cp(1) = 0.9;

xp(2) =  xmax;
yp(2) =  0;
zp(2) = -2;
cp(2) = 0.9;

xp(3) =  xmax;
yp(3) =  0;
zp(3) =  2;
cp(3) = 0.9;

xp(4) = xmin;
yp(4) =  0;
zp(4) =  2;
cp(4) = 0.9;

patch(xp,zp,yp,cp)

%---
% marker
%---

xmp = 11.25; ymp=2.0; zmp=-0.2;

%---
% particle
%---

for l=1:npts
 for nprtcl=1:1

xc=0;
yc=0;
zc=0;
a1=1.0;a2=1.0;
a3=0.25;

Ndiv=20;
[x,y,z]=ellipsoid(xc,yc,zc,a1,a2,a3,Ndiv);

theta = pi/2;   % around y
cs = cos(theta);
sn = sin(theta);
for i=1:Ndiv+1
 for j=1:Ndiv+1
  tmpx = cs*x(i,j)-sn*z(i,j);
  tmpz = sn*x(i,j)+cs*z(i,j);
  x(i,j) = tmpx;
  z(i,j) = tmpz;
end
end
  
if(nprtcl==1)
theta = orbit(10,l)*pi; % around z
else
theta = orbit(19,l)*pi; % around z
end

cs = cos(theta);
sn = sin(theta);
for i=1:Ndiv+1
 for j=1:Ndiv+1
  tmpx = cs*x(i,j)-sn*y(i,j);
  tmpy = sn*x(i,j)+cs*y(i,j);
  x(i,j) = tmpx;
  y(i,j) = tmpy;
end
end

if(nprtcl==1)
phi = orbit(11,l)*pi;  % around x
else
phi = orbit(20,l)*pi; % around x
end

cs = cos(phi);
sn = sin(phi);
for i=1:Ndiv+1
 for j=1:Ndiv+1
  tmpy = cs*y(i,j)-sn*z(i,j);
  tmpz = sn*y(i,j)+cs*z(i,j);
  y(i,j) = tmpy;
  z(i,j) = tmpz;
end
end

if(nprtcl==1)
cx = orbit(4,l);
cy = orbit(5,l);
cz = orbit(6,l);
else
cx = orbit(13,l);
cy = orbit(14,l);
cz = orbit(15,l);
end
x=x+cx; y=y+cy; z= z+cz;

hold on

 if(Icount==skip)
  mesh(x,z,y)
  surf(x,z,y)
  patch(xp,zp,yp,cp)
  axis([xmin xmax ymin ymax zmin zmax])
  hold on;
% plot3(xmp,zmp,ymp,'o'); hold off
%  view(0,90)
% figure
 end

end
                                                                                
 if(Icount==skip)
  Icount = 0;
  pause(0.2)
 end
                                                                                
 Icount=Icount+1;

end

 mesh(x,z,y)
 surf(x,z,y)
 patch(xp,zp,yp,cp)
axis([xmin xmax ymin ymax zmin zmax])
%  plot3(xmp,ymp,zmp)
%  view(0,90)

xlabel('x');
ylabel('z');
zlabel('y');
