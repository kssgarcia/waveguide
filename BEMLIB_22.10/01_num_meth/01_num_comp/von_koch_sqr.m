clear all
close all

%------------
% von_koch
%
% construct the periodic von_koch line
%-----------

L=1.0;      % period
m=1;        % generation
alpha=0.40;
alpha=0.80;

alphac=1.0-alpha;

%---------------
% basic line m=1
%---------------

height = sqrt(3.0*alpha^2+2.0*alpha-1.0);

nv=6;

 x(6)=L;
 y(6)=L*height;
in(6)=1;

 x(5)=L*alpha;
 y(5)=L*height;
in(5)=1;

 x(4)=L*alpha;
 y(4)=0;
in(4)=2;

 x(3)=-x(4);
 y(3)= y(4);
in(3)=2;

 x(2)=-x(5);
 y(2)= y(5);
in(2)=1;

 x(1)=-x(6);
 y(1)= y(6);
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
box

plot(x,y)

for i=1:nv
 if(in(i)==1)
%  plot(xfr(i),yfr(i),'.r')
 else
%  plot(xfr(i),yfr(i),'.g')
 end
end
