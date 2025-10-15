clear all
close all

%=========================
%
% FDLIB
%
% temperature distribution
% across a chimney wall
%=========================

a = 0.5;
Tinner = 300;
Touter = 10;
k = 0.15;      % thermal conductivity

N1=32; ratio1=10;  Isym1=1;   % boundary element preferences
N2=16; ratio2=5.0; Isym2=0;
N3=16; ratio3=2.0; Isym3=1;
N4=16; ratio4=0.2; Isym4=0;

NGL=6;  % Gauss-Legendre base points

%================================
% Boundary element discretization
%================================

%---
% diagonal segment 1
%---

x1=0;y1=a; x2=a; y2=0; sinit=0.0;

[xe1,ye1,se1,xm1,ym1,sm1]...
      = elm_line (N1,ratio1,x1,y1,x2,y2,sinit,Isym1);

Ic=0;    % collocation point counter
for i=1:N1
 Ic=Ic+1;
 xcl(Ic)=xm1(i);
 ycl(Ic)=ym1(i);
end

%---
% inner wall segment 2
%---

x1 = a;   y1=0.0;
x2 = 2*a; y2=0.0;
sinit=se1(N1+1);

[xe2,ye2,se2,xm2,ym2,sm2]...
      = elm_line (N2,ratio2,x1,y1,x2,y2,sinit,Isym2);

for i=1:N2
 Ic=Ic+1;
 xcl(Ic)=xm2(i);
 ycl(Ic)=ym2(i);
end

%---
% vertical segment 3
%---

x1=2*a; y1=0;
x2=2*a; y2=a;
sinit=se2(N2+1);

[xe3,ye3,se3,xm3,ym3,sm3]...
      = elm_line (N3,ratio3,x1,y1,x2,y2,sinit,Isym3);

for i=1:N3
 Ic=Ic+1;
 xcl(Ic)=xm3(i);
 ycl(Ic)=ym3(i);
end

%---
% horizontal segment 4
%---

x1=2*a; y1=a;
x2=0; y2=a;
sinit=se3(N3+1);

[xe4,ye4,se4,xm4,ym4,sm4]...
      = elm_line (N4,ratio4,x1,y1,x2,y2,sinit,Isym4);

for i=1:N4
 Ic=Ic+1;
 xcl(Ic)=xm4(i);
 ycl(Ic)=ym4(i);
end

Ncl =Ic;

%===============================================
% compute the slp and dlp influence coefficients
% at the collocation points
%===============================================

Itype=1;

%---------
for i=1:Ncl
%---------

 Jc=0;
 test=0.0;

 for j=1:N1
 Jc=Jc+1;
 Ising =0; if(Jc==i) Ising=1; end
 [SLP1(i,j), DLP1(i,j)] = sdlp ...
      ...
       (xcl(i),ycl(i),0.0 ...    % Evaluation point
       ,xe1(j),ye1(j),0.0 ...   % First element point
       ,xe1(j+1),ye1(j+1),0.0 ...   % Second element point
       ,NGL ...         % Gauss-Legendre quadrature order
       ,Ising ...       % Element singularity index
       ,Itype ...       % Element type
       ,0.0,0.0,0.0);   % Arc radius and center
 test = test+DLP1(i,j);
 end

 for j=1:N2
 Jc=Jc+1;
 Ising =0; if(Jc==i) Ising=1; end
 [SLP2(i,j), DLP2(i,j)] = sdlp ...
      ...
       (xcl(i),ycl(i),0.0 ...    % Evaluation point
       ,xe2(j),ye2(j),0.0 ...   % First element point
       ,xe2(j+1),ye2(j+1),0.0 ...   % Second element point
       ,NGL ...         % Gauss-Legendre quadrature order
       ,Ising ...       % Element singularity index
       ,Itype ...       % Element type
       ,0.0,0.0,0.0);   % Arc radius and center
 test = test+DLP2(i,j);
 end

 for j=1:N3
 Jc=Jc+1;
 Ising =0; if(Jc==i) Ising=1; end
 [SLP3(i,j), DLP3(i,j)] = sdlp ...
      ...
       (xcl(i),ycl(i),0.0 ...    % Evaluation point
       ,xe3(j),ye3(j),0.0 ...   % First element point
       ,xe3(j+1),ye3(j+1),0.0 ...   % Second element point
       ,NGL ...         % Gauss-Legendre quadrature order
       ,Ising ...       % Element singularity index
       ,Itype ...       % Element type
       ,0.0,0.0,0.0);   % Arc radius and center
 test = test+DLP3(i,j);
 end

 for j=1:N4
 Jc=Jc+1;
 Ising =0; if(Jc==i) Ising=1; end
 [SLP4(i,j), DLP4(i,j)] = sdlp ...
      ...
       (xcl(i),ycl(i),0.0 ...    % Evaluation point
       ,xe4(j),ye4(j),0.0 ...   % First element point
       ,xe4(j+1),ye4(j+1),0.0 ...   % Second element point
       ,NGL ...         % Gauss-Legendre quadrature order
       ,Ising ...       % Element singularity index
       ,Itype ...       % Element type
       ,0.0,0.0,0.0);   % Arc radius and center
 test = test+DLP4(i,j);
 end

% test

%--
end
%--

%==========================
% compile the linear system
%
% MAT x = rhs
%==========================

for i=1:Ncl

 rhs(i)=0.0;

 Jc=0;

 for j=1:N1
  Jc=Jc+1;
  MAT(i,Jc) = DLP1(i,j);
 end

 for j=1:N2
  Jc=Jc+1;
  MAT(i,Jc) =-SLP2(i,j);
  rhs(i) = rhs(i)-DLP2(i,j)*Tinner;
 end

 for j=1:N3
  Jc=Jc+1;
  MAT(i,Jc) = DLP3(i,j);
 end

 for j=1:N4
  Jc=Jc+1;
  MAT(i,Jc) =-SLP4(i,j);
  rhs(i) = rhs(i)-DLP4(i,j)*Touter;
 end
 
end

for i=1:N1
 MAT(i,i)=MAT(i,i)-0.5;
end
for i=1:N2
 l=N1+i;
 rhs(l)=rhs(l)+0.5*Tinner;
end
for i=1:N3
 l=N1+N2+i;
 MAT(l,l)=MAT(l,l)-0.5;
end
for i=1:N4
 l=N1+N2+N3+i;
 rhs(l)=rhs(l)+0.5*Touter;
end

%======
% solve
%======

solution = rhs/MAT';

%==============================
% boundary temperature and flux
%==============================

for i=1:N1
 temp(i)=solution(i);
 flux(i)=00;
end
for i=1:N2
 l=N1+i;
 temp(l)=Tinner;
 flux(l)=-k*solution(l);
end
for i=1:N3
 l=N1+N2+i;
 temp(l)=solution(l);
 flux(l)=0.0;
end
for i=1:N4
 l=N1+N2+N3+i;
 temp(l)=Touter;
 flux(l)=k*solution(l);
end

%==============================
% plot the boundary temperature
%==============================

hold on
plot(xe1,ye1,'-')
plot(xe2,ye2,'-')
plot(xe3,ye3,'-')
plot(xe4,ye4,'-')
plot(xcl,ycl,'o')
for i=1:Ncl
 plot3([xcl(i),xcl(i)],[ycl(i),ycl(i)],[0, temp(i)],'-o')
end
set(gca,'fontsize',15)
xlabel('x(m)','fontsize',15)
ylabel('y(m)','fontsize',15)
zlabel('T (C)','fontsize',15)

%=======================
% plot the boundary flux
%=======================

figure
hold on
plot(xe1,ye1,'-')
plot(xe2,ye2,'-')
plot(xe3,ye3,'-')
plot(xe4,ye4,'-')
plot(xcl,ycl,'o')
for i=1:Ncl
 plot3([xcl(i),xcl(i)],[ycl(i),ycl(i)],[0, flux(i)],'o-')
end

xlabel('x(m)','fontsize',15)
ylabel('y(m)','fontsize',15)
zlabel('flux (Watt)','fontsize',15)
set(gca,'fontsize',15)
