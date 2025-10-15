clear all;
close all;

%===============================
% shear flow past periodic bumps
%===============================

%---
% parameters
%---

alpha=0.6;
alpha=0.2;
alpha=0.3;
alpha=0.4;

L=1.0;

Nxil=8;   % number of left xi divisions
Nxic=32;  % number of central xi divisions
Nxir=8;  % number of right xi divisions
Neta=32;  % number eta divisions

rho=1.0;

visc=0.1;
visc=0.05;
visc=0.02;
visc=0.01;
visc=1.0;

relax=0.25;
relax=0.10;
relax=0.5;
tol=0.0000001;

Niteri=20;       % number of inner iterations
Niteri=3;       % number of inner iterations
Niteri=5;       % number of inner iterations
Niterg=6000;   % number of global iterations

vort_init=0.0;

%---
% prepare
%---

Nxi=Nxil+Nxic+Nxir;

nu=visc/rho;


%================
% grid generation
%================

[xi,eta,x,y,h]= conf2_grid(alpha,L,Nxil,Nxic,Nxir,Neta);

 xi(Nxi+2)= xi(2)+L;  % wrap
eta(Nxi+2)=eta(2);

for j=1:Neta+1      % wrap
 x(Nxi+2,j)=x(2,j)+L;
 y(Nxi+2,j)=y(2,j);
end

%========
% compute the differentiation coefficients
%========

[bl,bc,br, bb,bx,bt, cl,cc,cr,cb,ct] = conf2_diff  ...
...
   (Nxi,Neta,xi,eta,h);

%====================================
% stream function boundary conditions
%====================================

for i=1:Nxi+1
 psib(i)=0.0;  % bottom
 psit(i)=1.0;  % top
end

%===========
% initialize
%===========

for i=1:Nxi+2
 for j=1:Neta+1
  psi(i,j) = 0.0;
  vort(i,j) = vort_init;
 end
end

%==================
% global iterations
%==================

for iterg=1:Niterg

 save = psi;

%===========================
% solve the Poisson equation
% for the stream function
%===========================

[psi, iter, Iflag] = conf2_pois_gs_dpr ...
...
   (Nxi,Neta,xi,eta,x,y,h,vort,Niteri,tol,relax,psib,psit,psi);

%=====
% compute the velocity
%======

%---
% zero lower wall velocity
%----

for i=1:Nxi+1
 uxi(i,1)  =0.0;
 ueta(i,1) =0.0;
 ux(i,1)   =0.0;
 uy(i,1)   =0.0;
end

for j=2:Neta
 for i=2:Nxi+1

  uxi(i,j) = bb(i,j)*psi(i,j-1)+ bx(i,j)*psi(i,j) ...
           + bt(i,j)*psi(i,j+1);
  teta_x(i,j) =  bb(i,j)*x(i,j-1)+ bx(i,j)*x(i,j) ...   % eta tangential vector
               + bt(i,j)*x(i,j+1);
  teta_y(i,j) =  bb(i,j)*y(i,j-1)+ bx(i,j)*y(i,j) ...
               + bt(i,j)*y(i,j+1);
%  norm = sqrt(teta_x(i,j)^2+teta_y(i,j)^2)
  ueta(i,j) = -(bl(i,j)*psi(i-1,j)+bc(i,j)*psi(i,j) ...
               +br(i,j)*psi(i+1,j) );
  txi_x(i,j) =  bl(i,j)*x(i-1,j)+ bc(i,j)*x(i,j) ...   % xi tangential vector
              + br(i,j)*x(i+1,j);
  txi_y(i,j) =  bl(i,j)*y(i-1,j)+ bc(i,j)*y(i,j) ...
              + br(i,j)*y(i+1,j);
%  norm = sqrt(txi_x(i,j)^2+txi_y(i,j)^2)
  ux(i,j) = uxi(i,j)*txi_x(i,j)+ueta(i,j)*teta_x(i,j);
  uy(i,j) = uxi(i,j)*txi_y(i,j)+ueta(i,j)*teta_y(i,j);
 end
end

%---
% wrap
%---

for j=2:Neta
 uxi(1,j) = uxi(Nxi+1,j); 
 ueta(1,j) = ueta(Nxi+1,j);
 txi_x(1,j) = txi_x(Nxi+1,j);
 txi_y(1,j) = txi_y(Nxi+1,j);
 teta_x(1,j) = teta_x(Nxi+1,j);
 teta_y(1,j) = teta_y(Nxi+1,j); 
 ux(1,j) = ux(Nxi+1,j); 
 uy(1,j) = uy(Nxi+1,j); 
end

%quiver(x(1:Nxi+1,1:Neta),y(1:Nxi+1,1:Neta),uxi,ueta);
%quiver(x(1:Nxi+1,1:Neta),y(1:Nxi+1,1:Neta),txi_x,txi_y);
%quiver(x(1:Nxi+1,2:Neta),y(1:Nxi+1,2:Neta) ...
%      ,ux(1:Nxi+1,2:Neta),uy(1:Nxi+1,2:Neta));
%quiver(x(1:Nxi+1,1:Neta),y(1:Nxi+1,1:Neta),teta_x,teta_y);

%figure
%hold on
% for i=1:Nxi
%  for j=2:Neta-1
%  A=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)];
%  B=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1)];
%  C=[ux(i,j),ux(i+1,j),ux(i+1,j+1),ux(i,j+1)];
%  patch(A,B,C,C);
%  end
% end

%=====================================
% source of Poisson equation for omega
%=====================================

for j=2:Neta
 for i=2:Nxi+1
   Det = eta(j+1)-eta(j-1);
   Dxi = xi(i+1)-xi(i-1);
   vort_src(i,j) = uxi(i,j)*(vort(i+1,j)-vort(i-1,j))/(Dxi*h(i,j)) ...
                 + ueta(i,j)*(vort(i,j+1)-vort(i,j-1))/(Det*h(i,j));
   vort_src(i,j) = -vort_src(i,j)/nu;
 end
end

%----
% boundary condition for the vorticity
% on the lower wall
%----

for i=1:Nxi+1
 ds = sqrt( (x(i,2)-x(i,1))^2+(y(i,2)-y(i,1))^2 );
 vort_bot(i)=-uxi(i,2)/ds;
end

%--------------
% solve for the vorticity
%--------------

[vort, iter, Iflag] = conf2_pois_gs_ndpr ...
...
   (Nxi,Neta,xi,eta,x,y,h,vort_src,Niteri,tol,relax,vort_bot,vort);

%------------------
% monitor the error
%------------------

cormax=0.0;

for i=1:Nxi
 for j=1:Neta
  res = abs(psi(i,j)-save(i,j));
  if(res>cormax)
    cormax = res;
  end
 end
end

cormax

if(cormax<tol) 
 break
end

%===
end
%===


iterg

%=========
% plotting
%=========

%---
% streamfunction
%---

Iskip=0;
if(Iskip==0)
 figure
 hold on

 for i=1:Nxi
  for j=1:Neta
  A=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)];
  B=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1)];
  C=[psi(i,j),psi(i+1,j),psi(i+1,j+1),psi(i,j+1)];
  patch(A,B,C,C);
  end
 end
 set(gca,'fontsize',15)
 xlabel('x','fontsize',15)
 ylabel('y','fontsize',15)
 zlabel('\psi','fontsize',15)
 %axis equal
 axis ([-0.5 0.5 0 1.5 0 1.0]);
 box
end

%---
% vorticity
%---

Iskip=0;
if(Iskip==0)
 figure
 hold on

 for i=1:Nxi
  for j=1:Neta
  A=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)];
  B=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1)];
  C=[vort(i,j),vort(i+1,j),vort(i+1,j+1),vort(i,j+1)];
  patch(A,B,C,C);
  end
 end
 set(gca,'fontsize',15)
 xlabel('x','fontsize',15)
 ylabel('y','fontsize',15)
 zlabel('\omega','fontsize',15)
 %axis equal
 %axis ([-0.5 0.5 0 1.5 0 1.0]);
 box
end

%---
% vorticity
%---

figure
contour(x,y,vort,64)
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)

%---
% velocity vector field
%---

figure
quiver(x(1:Nxi+1,1:Neta),y(1:Nxi+1,1:Neta) ...
      ,ux(1:Nxi+1,1:Neta),uy(1:Nxi+1,1:Neta));

%---
% velocity vector field
%---

figure
 for i=1:Nxi
  for j=1:Neta-1
  A=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)];
  B=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1)];
  u1=sqrt(ux(i,j)^2    +uy(i,j)^2);
  u2=sqrt(ux(i+1,j)^2  +uy(i+1,j)^2);
  u3=sqrt(ux(i+1,j+1)^2+uy(i+1,j+1)^2);
  u4=sqrt(ux(i,j+1)^2  +uy(i,j+1)^2);
  C=[ux(i,j),ux(i+1,j),ux(i+1,j+1),ux(i,j+1)];
  C=[u1,u2,u3,u4];
  C=[uy(i,j),uy(i+1,j),uy(i+1,j+1),uy(i,j+1)];
  patch(A,B,C,C);
  end
 end
 set(gca,'fontsize',15)
 xlabel('x','fontsize',15)
 ylabel('y','fontsize',15)
 zlabel('|u|','fontsize',15)
 zlabel('u_y','fontsize',15)
 %axis equal
 %axis ([-0.5 0.5 0 1.5 0 1.0]);
 box
