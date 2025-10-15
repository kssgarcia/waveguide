%=======
% stability graphs
% for the convection--diffusion equation
%=====

close all
clear all
hold on

im = sqrt(-1);
Nx=64;
Ny=64;

%-------------
opt= 4;  % DF
opt= 6;  % FRAC STEPS
opt= 5;  % MACCORMACK
opt= 2;  % UPWIND
opt= 1;  % FTCS
opt= 3;  % UPWIND3
%-------------

if(opt==1)
 ax=-0.2;bx= 1;
 ay=-1.2;by= 1.2;
else
 ax=-0.2;bx= 1;
 ay=-1.2;by= 1.2;
end
ax= 0.0;bx= 1;
ay=-1.2;by= 1.2;

xx = linspace(ax, bx, Nx+1);
yy = linspace(ay, by, Ny+1);

Ntheta=32;
Dtheta=2*pi/Ntheta;

%---
for i=1:Nx+1

 al=xx(i);

 for j=1:Ny+1
  c=yy(j);
 
  max = 0.0;

  for k=1:Ntheta+1
   theta=(k-1)*Dtheta;
%---
   if(opt==1)
     G = 1-2*al+2*al*cos(theta)+im*c*sin(theta);   % FTCS
%---
   elseif(opt==2)
    G = 1-(c+2*al)+(c+2*al)*cos(theta)+im*c*sin(theta); % UPWIND
%---
   elseif(opt==3)
    G = -c/6*exp(2*im*theta)+(c+al)*exp(im*theta) ...
         + 1-c/2-2*al +(-c/3+al)*exp(-im*theta);
%---
   elseif(opt==4)
     A = 1+2*al;
     B = -2*(2*al*cos(theta)+im*c*sin(theta));
     C = -1+2*al;
     disc = B*B-4*A*C;
     G1 = 0.5*(-B+sqrt(disc))/A;
     G2 = 0.5*(-B-sqrt(disc))/A;
     G = G1;
     if(abs(G2)>G)
      G = abs(G2);
     end
%---
   elseif(opt==5)
     eB = exp(-2*im*theta);
     eA = exp(-im*theta);
     e1 = exp( im*theta);
     e2 = exp(2*im*theta);
     fAs = (1+c)*e1-c    + al*(e2-2*e1+1);
     fs  = (1+c)   -c*eA + al*(e1-2   +eA);
     f1s = (1+c)*eA-c*eB + al*(1 -2*eA+eB);
     G = 0.5*(1+fs-c*(fs-fAs)+al*(fAs-2*fs+f1s));
%---
   elseif(opt==6)
     eB = exp(-2*im*theta);
     eA = exp(-im*theta);
     e1 = exp( im*theta);
     e2 = exp(2*im*theta);
     fAs = 0.5*c*e2+e1 -0.5*c;
     fs  = 0.5*c*e1+1.0-0.5*c*eA;
     f1s = 0.5*c   +eA -0.5*c*eB;
     G = al*fAs+(1-2*al)*fs+al*f1s;
   end
%---

   if(abs(G)>max)
    max = abs(G);
   end

  end

  zz(j,i)=max;
  if(max<1.001)
   plot(al,c,'kx')
  end

 end
end

%cn = contour(xx,yy,zz,[1 1]);
axis square
xlabel('\alpha','fontsize',15)
ylabel('c','fontsize',15)
axis([ax bx ay by])
set(gca,'fontsize',15)
plot([ax bx],[0 0])
plot([0 0],[ay by])
box

break

n = cn(2,1);
for k=1:n
 xc(k) = cn(1,k+1);
 yc(k) = cn(2,k+1);
end

if(opt==21)
 patch([ax, bx, bx, ax],[ay, ay, by, by],'y')
 patch(xc,yc,'w')
elseif(opt==22)
 patch([0, bx, bx, 0],[ay, ay, by, by],'y')
%elseif(opt==23|opt==24)
% xc(n+1)=bx;yc(n+1)=ay;
% xc(n+2)=bx;yc(n+2)=by;
% xc(n+3)=ax;yc(n+3)=by;
% patch(xc,yc,'y')
else
 patch(xc,yc,'y')
end
