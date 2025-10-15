%=======
% stability graphs
%=====

close all
clear all

im = sqrt(-1);
Nx=128;
Ny=128;

opt= 3;  % RK3
opt= 4;  % RK4
opt= 2;  % RK2
opt= 1;  % RK1

opt= 12;  % AB2
opt= 13;  % AB3
opt= 14;  % AB4
opt= 15;  % AB5

opt= 21;  % AM1
opt= 22;  % AM2
opt= 23;  % AM3
opt= 24;  % AM4

if(opt==1)
 ax=-2;bx= 4;
 ay=-3;by= 3;
elseif(opt==2)
 ax=-2;bx= 4;
 ay=-3;by= 3;
elseif(opt==3)
 ax=-2;bx= 4;
 ay=-3;by= 3;
elseif(opt==4)
 ax=-2;bx= 4;
 ay=-3;by= 3;
elseif(opt==21)
 ax=-4;bx= 2;
 ay=-3;by= 3;
elseif(opt==23)
 ax=-1;bx= 9;
 ay=-5;by= 5;
elseif(opt==24)
 ax=-1;bx= 5;
 ay=-3;by= 3;
else
 ax=-0.5;bx= 1.5;
 ay=-1;by= 1;
end

xx = linspace(ax, bx, Nx+1);
yy = linspace(ay, by, Ny+1);

for i=1:Nx+1

  x=xx(i);
 for j=1:Ny+1
  y=yy(j);
  z = x+im*y;
  z = -z;

  if(opt==1)

   zz(j,i) = abs(1+z);

  elseif(opt==2)

   zz(j,i) = abs(1+z+0.5*z^2);

  elseif(opt==3)

   zz(j,i) = abs(1+z+0.5*z^2+z^3/6);

  elseif(opt==4)

   zz(j,i) = abs(1+z+0.5*z^2+z^3/6+z^4/24);

  elseif(opt==12)

   A = 1.0;
   B = -(1+1.5*z);
   C = 0.5*z;
   D = sqrt(B*B-4*A*C);
   root1 = (-B+D)/(2*A);
   root2 = (-B-D)/(2*A);
   cand1 = abs(root1);
   cand2 = abs(root2);
   zz(j,i) = cand1;
   if(cand2>cand1)
    zz(j,i) =cand2;
   end

  elseif(opt==13)

   MAT=zeros(3,3);
   MAT(1,1) = 1+23/12*z; MAT(1,2)=-4/3*z; MAT(1,3)=5/12*z;
   MAT(2,1) =1.0; MAT(3,2)=1.0;
   eigen=eig(MAT);
   cand1 = abs(eigen(1));
   cand2 = abs(eigen(2));
   cand3 = abs(eigen(3));
   zz(j,i) = cand1;

   if(cand2>zz(j,i))
    zz(j,i) =cand2;
   end
   if(cand3>zz(j,i))
    zz(j,i) =cand3;
   end

  elseif(opt==14)

   MAT=zeros(4,4);
   MAT(1,1)=1+55/24*z;MAT(1,2)=-59/24*z;
   MAT(1,3)=37/24*z; MAT(1,4)=-9/24*z;
   MAT(2,1)=1;
   MAT(3,2)=1;
   MAT(4,3)=1;
   eigen=eig(MAT);
   cand1 = abs(eigen(1));
   cand2 = abs(eigen(2));
   cand3 = abs(eigen(3));
   cand4 = abs(eigen(4));
   zz(j,i) = cand1;
   if(cand2>zz(j,i))
    zz(j,i) =cand2;
   end
   if(cand3>zz(j,i))
    zz(j,i) =cand3;
   end
   if(cand4>zz(j,i))
    zz(j,i) =cand4;
   end

  elseif(opt==15)

   MAT=zeros(5,5);
   MAT(1,1)=1+1901/720*z;MAT(1,2)=-2774/720*z;
   MAT(1,3)=2616/720*z; MAT(1,4)=-1274/720*z;
   MAT(1,5)=251/720*z;
   MAT(2,1)=1;
   MAT(3,2)=1;
   MAT(4,3)=1;
   MAT(5,4)=1;
   eigen=eig(MAT);
   cand1 = abs(eigen(1));
   cand2 = abs(eigen(2));
   cand3 = abs(eigen(3));
   cand4 = abs(eigen(4));
   cand5 = abs(eigen(5));
   zz(j,i) = cand1;
   if(cand2>zz(j,i))
    zz(j,i) =cand2;
   end
   if(cand3>zz(j,i))
    zz(j,i) =cand3;
   end
   if(cand4>zz(j,i))
    zz(j,i) =cand4;
   end
   if(cand5>zz(j,i))
    zz(j,i) =cand5;
   end

  elseif(opt==21)

   zz(j,i) = abs( 1/(1-z) );

  elseif(opt==22)

   zz(j,i) = abs( (1+0.5*z)/(1-0.5*z) );

  elseif(opt==23)

   MAT=zeros(2,2);
   MAT(1,1)=(1+8/12*z)/(1-5/12*z);
   MAT(1,2)=-1/12*z/(1-5/12*z);
   MAT(2,1)=1; 
   eigen=eig(MAT);
   cand1 = abs(eigen(1));
   cand2 = abs(eigen(2));
   zz(j,i) = cand1;
   if(cand2>cand1)
    zz(j,i) =cand2;
   end

  elseif(opt==24)

   MAT=zeros(3,3);
   fc = 1/(1-9/24*z);
   MAT(1,1)=(1+19/24*z)*fc;
   MAT(1,2)=-5*z/24*fc;
   MAT(1,3)= z/24*fc;
   MAT(2,1)=1;
   MAT(3,2)=1;
   eigen=eig(MAT);
   cand1 = abs(eigen(1));
   cand2 = abs(eigen(2));
   cand3 = abs(eigen(3));
   zz(j,i) = cand1;
   if(cand2>zz(j,i)) zz(j,i)=cand2; end
   if(cand3>zz(j,i)) zz(j,i)=cand3; end

   end

 end
end

cn = contour(xx,yy,zz,[1 1]);
axis equal
hold on

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


axis([ax bx ay by])
plot([ax bx],[0 0])
plot([0 0],[ay by])
xlabel('-\Delta t \lambda_{A,R}','fontsize',15)
ylabel('-\Delta t \lambda_{A,I}','fontsize',15)
set(gca,'fontsize',15)
