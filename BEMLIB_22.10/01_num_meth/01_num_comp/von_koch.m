clear all
close all

%-------------------------------------
% von_koch
%
% construct the periodic von_koch line
%
% 1/3 < alpha < 1
%-------------------------------------

L=1.0;      % period
m=4;        % generation
alpha=0.36;
alpha=0.80;
alpha=0.70;

alphac=1.0-alpha;

%---------------
% basic line m=1
%
% The array "in" is a vertex index that is
% equal to 1 if the corresponding corner
% angle is obtuse, and 2 otherwise
%---------------

height = sqrt(3.0*alpha^2+2.0*alpha-1.0);

nv=4;

 x(4)=L;
 y(4)=L*height;
in(4)=1;

 x(3)=L*alpha;
 y(3)=0;
in(3)=2;

 x(2)=-x(3);
 y(2)= y(3);
in(2)=2;
 x(1)=-x(4);
 y(1)= y(4);
in(1)=1;

%---
% run over generations
%---

if(m>1)

%---
for i=2:m
%---

%---------------------------
% Save the vertices 
% of the previous generation
%---------------------------

  for j=1:nv
   xaux(j)=x(j); yaux(j)=y(j); iaux(j)=in(j);
  end

  kv=nv;
  nv=3*4^(i-1)+1;

  x(nv)=xaux(kv); y(nv)=yaux(kv); in(nv)=iaux(kv);
  j=nv;

  for l=1:kv-1

    sgn=1.0;
%    if(rand>0.5)
%     sgn=-1.0;
%    end

    k=kv-l;

    xa=xaux(k);       ya=yaux(k);
    xb=xaux(k+1);     yb=yaux(k+1);
    xc=0.5D0*(xa+xb); yc=0.5D0*(ya+yb);
    dx=xb-xa;         dy=yb-ya;
    ds = sqrt(dx^2+dy^2);
    height = sqrt(3.0*(4.0*alpha-1.0))/6.0;

    j=j-1;
    x(j)=xc+alphac*dx/3.0;
    y(j)=yc+alphac*dy/3.0;
    in(j)=2;
    j=j-1;
    x(j)=xc-sgn*dy*height;
    y(j)=yc+sgn*dx*height;
    in(j)=1;
    j=j-1;
    x(j)=xc-alphac*dx/3.0;
    y(j)=yc-alphac*dy/3.0;
    in(j)=2;
    j=j-1;
    x(j)=xa;
    y(j)=ya;
    in(j)=in(k);

  end

%---
end
%---

end


%---
% plotting
%---

hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis equal
axis off

plot(x,y)

for i=1:nv
 if(in(i)==1)
   plot(x(i),y(i),'.r')
 else
   plot(x(i),y(i),'.g')
 end
end

%------------------------
% rectangular projections
%-------

Ic=0;     % vertex counter

for i=1:nv
 if(in(i)==1)
  if(i==1)
   dx=x(2)-x(1);
   dy=0.0;
  elseif(i==nv)
   dx=x(nv)-x(nv-1);
   dy=0.0;
  else
   dx=0.5*(x(i+1)-x(i-1));
   dy=0.5*(y(i+1)-y(i-1));
  end
  Ic=Ic+1;
  xs(Ic)=x(i)-dx;
  ys(Ic)=y(i)-dy;
  Ic=Ic+1;
  xs(Ic)=x(i)+dx;
  ys(Ic)=y(i)+dy;
 else
  Ic=Ic+1;
  xs(Ic)=x(i);
  ys(Ic)=y(i);
 end
end

nvs=Ic;

figure
hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis equal
axis off

plot(xs,ys,'.-')
