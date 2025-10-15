close all
clear all

%--------------------
% flow in a cavity in primary
% variables using the velocity/pressure
% formulation
%--------------------

%-----------------------
% settings and parameters
%-----------------------

Lx=1.0;  % cavity size
Ly=1.0;

Nx=32; % grid size
Ny=32;

Dt=0.001;

visc=0.01; % viscosity
dens=1.0; % density

lidU=1.0; % lid velocity

Nstep=2000;

%--------------------------------------------
% parameters for solving the Poisson equation
%--------------------------------------------

itermax=50000;
tol=0.000001;
relax=0.2;

qleft=0.0;  % Neumann boundary condition
qright=0.0;
qbot=0.0;
qtop=0.0;

Ishift=1;

%----------------------------------------
% parameters for slip velocity iterations
%----------------------------------------

slipN = 50;        % max number of inner iterations
sliptol = 0.00001; % tolerance
sliprel = 1.0;     % relaxation

%--------
% prepare
%--------

nu=visc/dens; % kinematic viscosity

Dx = Lx/Nx;
Dy = Ly/Ny;

Dx2 = 2.0*Dx;
Dy2 = 2.0*Dy;

Dxs = Dx*Dx;
Dys = Dy*Dy;

Dtor = Dt/dens;

for i=1:Nx+1
  Ulid(i) = lidU;
end

%-----------------------------------
% define grid lines
% initialize time, velocity  (U, V)
% and the projection function (T)
%-----------------------------------

time = 0.0;

for j=1:Ny+1
 for i=1:Nx+1
   X(i,j) = (i-1.0)*Dx;
   Y(i,j) = (j-1.0)*Dy;
   U(i,j) = 0.0;   % x velocity
   V(i,j) = 0.0;   % y velocity
   T(i,j) = 0.0;   % projection function
 end
end

for i=1:Nx+1
   U(i,Ny+1) = Ulid(i);
end

%-----------------------------------
% naive velocity boundary conditions
%-----------------------------------

for i=1:Nx+1
  BCxt(i) = Ulid(i);  % top wall
  BCxb(i) = 0.0;      % bottom wall
end

for j=1:Ny+1
  BCyl(j) = 0.0; % left wall
  BCyr(j) = 0.0; % right wall
end

%--------------
% time stepping
%--------------

%===============
for step=1:Nstep
%===============

%---
% animation
%---

if(step==1)
  Handle1 = quiver(X,Y,U,V);
  set(Handle1, 'erasemode', 'xor');
  axis ([0, Lx, 0, Ly])
  axis equal
  set(gca,'fontsize',15)
  xlabel('x','fontsize',15)
  ylabel('y','fontsize',15)
  hold on
  plot([0, Lx, Lx, 0, 0],[0, 0, Ly, Ly, 0]);
else
  set(Handle1,'UData',U,'VData',V);
  pause(0.01)
  drawnow
end

%-------------------------------------
% initialize the intermediate velocity
%-------------------------------------

for j=1:Ny+1
  for i=1:Nx+1
   UU(i,j) = U(i,j);
   VV(i,j) = V(i,j);
  end
end

%=============================
% will perform inner iterations
% for the projection function
% to satisfy the no-slip
% boundary condition
%=============================

    for inner=1:slipN

%---
% zero the entries of the tridiagonal matrix
%---

 for i=1:Nx-1
   atr(i) = 0.0;
   btr(i) = 0.0;
   ctr(i) = 0.0;
 end

%-------------------------------------------------
%  Integrate Conv-Diff equation in the x direction
%  using the Crank-Nicolson method
%  Advance the velocity from u^n to u*
%-------------------------------------------------

  Iskip=0;
  if(Iskip==0)

  AL = nu*Dt/Dxs;

%--------------
  for j=2:Ny          % run over rows
%----

    for i=2:Nx
      RC = U(i,j)*Dt/Dx;     % cell Reynolds number
      C1 =   RC + 2.0*AL;
      C2 =        4.0*(1.0-AL);
      C3 =  -RC + 2.0*AL;
      ctr(i-1) =  -RC - 2.0*AL;
      atr(i-1) =        4.0*(1.0+AL);
      btr(i-1) =   RC - 2.0*AL;
      Blu(i-1) = C1*U(i-1,j) + C2*U(i,j) + C3*U(i+1,j);  % right-hand side
      Blv(i-1) = C1*V(i-1,j) + C2*V(i,j) + C3*V(i+1,j);  % right-hand side
    end

    Blv(1)    = Blv(1)    - ctr(1)   *BCyl(j);
    Blv(Nx-1) = Blv(Nx-1) - btr(Nx-1)*BCyr(j);

    Xsu = thomas (Nx-1,atr,btr,ctr,Blu);  % x component
    Xsv = thomas (Nx-1,atr,btr,ctr,Blv);  % y component

    for k=1:Nx-1
     UU(k+1,j) = Xsu(k);
     VV(k+1,j) = Xsv(k);
    end

%---
   end  % End of running over rows
%-----------

  for i=1:Nx-1     % reset the tridiagonal matrix
    atr(i) = 0.0;
    btr(i) = 0.0;
    ctr(i) = 0.0;
  end

 end  % of skip

%-------------------------------------------------
%  Integrate Conv-Diff equation in the y direction
%  using the Crank-Nicolson method
%  Advance the velocity from u* to u**
%-------------------------------------------------

  Iskip=0;
  if(Iskip==0)

   AL = nu*Dt/Dys;

%---
   for i=2:Nx     % run over columns
%---
     for j=2:Ny     % from bottom to top
       RC = V(i,j)*Dt/Dy;       % cell Reynolds number
       C1 =  RC +2.0*AL;
       C2 =      4.0*(1.0-AL);
       C3 = -RC +2.0*AL;
       ctr(j-1) = -RC -2.0*AL;
       atr(j-1) =      4.0*(1.0+AL);
       btr(j-1) =  RC -2.0*AL;
       Blu(j-1) = C1*UU(i,j-1) + C2*UU(i,j) + C3*UU(i,j+1); % right-hand side
       Blv(j-1) = C1*VV(i,j-1) + C2*VV(i,j) + C3*VV(i,j+1); % right-hand side
     end

     Blu(1)    = Blu(1)    - ctr(1)   *BCxb(i); %    boundary conditions
     Blu(Ny-1) = Blu(Ny-1) - btr(Ny-1)*BCxt(i); %    boundary conditions

     Xsu = thomas (Ny-1,atr,btr,ctr,Blu);  % x component
     Xsv = thomas (Ny-1,atr,btr,ctr,Blv);  % y component

     for k=1:Ny-1
       UU(i,k+1) = Xsu(k);
       VV(i,k+1) = Xsv(k);
     end
%---
  end   % of run over columns
%---

 end  % of skip

%-------------------------------------
% Compute intermediate compressibility
% by centered differences
%
%       Divus = Div u**
%-------------------------------------

% initialize

  for j=1:Ny+1
    for i=1:Nx+1
      Divus(i,j) = 0.0;
    end
  end

% interior nodes

  for i=2:Nx
    for j=2:Ny
      DuDx = (UU(i+1,j)-UU(i-1,j))/Dx2;
      DvDy = (VV(i,j+1)-VV(i,j-1))/Dy2;
      Divus(i,j) = DuDx+DvDy;
     end
   end


% left wall

  for j=1:Ny+1
   DuDx = (-3.0*UU(1,j)+4.0*UU(2,j)-UU(3,j))/Dx2;
   DvDy = 0.0;
   Divus(1,j) = DuDx+DvDy;
  end

   save11l = Divus(1,j);     % save for corners
   save12l = Divus(1,Ny+1);   % save for corners

% bottom wall

   for i=1:Nx+1
     DuDx = 0.0;
     DvDy = (-3.0*VV(i,1)+4.0*VV(i,2)-VV(i,3))/Dy2;
     Divus(i,1) = DuDx+DvDy;
   end

   save11b = Divus(1,1);     % save for corners
   save21b = Divus(Nx+1,1);   % save for corners

% right wall

   for j=1:Ny+1
     DuDx = (3.0*UU(Nx+1,j)-4.0*UU(Nx,j)+UU(Nx-1,j))/Dx2;
     DvDy = 0.0;
     Divus(Nx+1,j) = DuDx+DvDy;
   end

   save21r = Divus(Nx+1,1);     % save for corners
   save22r = Divus(Nx+1,Ny+1);   % save for corners

% top wall

   for i=1:Nx+1
%     DuDx = (Ulid(i+1)-Ulid(i-1))/Dx2;
     DuDx = 0.0;
     DvDy = (3.0*VV(i,Ny+1)-4.0*VV(i,Ny)+VV(i,Ny-1))/Dy2;
     Divus(i,Ny+1) = DuDx+DvDy;
   end

   save12t = Divus(1,   Ny+1);      % save for corners
   save22t = Divus(Nx+1,Ny+1);      % save for corners

% corners by averaging

   Divus(1,      1) = 0.5*(save11l+save11b);
   Divus(1,   Ny+1) = 0.5*(save12l+save12t);
   Divus(Nx+1,   1) = 0.5*(save21b+save21r);
   Divus(Nx+1,Ny+1) = 0.5*(save22r+save22t);

%----------------------------------
% Solve for the projection function
% by Gauss-Siedel (GS) iterations
%----------------------------------

   Iskip=0;
   if(Iskip==0)

%---
% source term
%-----

   for i=1:Nx+1
    for j=1:Ny+1
      source(i,j) = -Divus(i,j)/Dtor;
     end
   end

  [T,iter,Iflag] = pois_gs_nnnn ...
 ...
    (Nx,Ny,Dx,Dy,source,itermax,tol,relax ...
 ...
    ,qleft,qright,qbot,qtop,T,Ishift);

  if(Iflag==0)
   disp "Poisson solver did not converge"
   break
  end

%----------------------------------
% Project the velocity at all nodes
% except at the corner nodes
%----------------------------------

%---
% interior nodes
%----

  for i=2:Nx
    for j=2:Ny
     DTDx = (T(i+1,j)-T(i-1,j))/Dx2;
     DTDy = (T(i,j+1)-T(i,j-1))/Dy2;
     UU(i,j) = UU(i,j) - Dtor*DTDx;
     VV(i,j) = VV(i,j) - Dtor*DTDy;
    end
  end

%--------------
% lower boundary
%
% Use forward differences with
% special treatment of the near-corner nodes
%--------------

   for i=2:Nx
      DTDy = 0.0;
      if(i==2)
         DTDx = (-3.0*T(2,1)+4.0*T(3,1)-T(4,1))/Dx2;
      elseif(i==Nx)
         DTDx =-(-3.0*T(Nx,1)+4.0*T(Nx-1,1)-T(Nx-2,1))/Dx2;
      else
         DTDx = (T(i+1,1)-T(i-1,1))/Dx2;
      end
      UU(i,1) = BCxb(i) - Dtor*DTDx;
      VV(i,1) =         - Dtor*DTDy;
    end

%--------------
% Upper boundary
%
% Use backward differences
% special treatment for the near-corner nodes
%--------------

   for i=2:Nx
      DTDy = 0.0;
      if(i==2)
        DTDx = (-3.0*T(2,Ny+1)+4.0*T(3,Ny+1)-T(4,Ny+1))/Dx2;
      elseif(i==Nx)
        DTDx =-(-3.0*T(Nx,Ny+1)+4.0*T(Nx-1,Ny+1)-T(Nx-2,Ny+1))/Dx2;
      else
        DTDx = (T(i+1,Ny+1)-T(i-1,Ny+1))/Dx2;
      end
      UU(i,Ny+1) = BCxt(i) - Dtor*DTDx;
      VV(i,Ny+1) =         - Dtor*DTDy;
   end

%-------------
% Left boundary
%
% Use forward differences
% special treatment for the near-corner nodes
%-------------

   for j=2:Ny
      DTDx = 0.0;
      if(j==2)
        DTDy = (-3.0*T(1,2)+4.0*T(1,3)-T(1,4))/Dy2;
      elseif(j==Ny)
        DTDy =-(-3.0*T(1,Ny)+4.0*T(1,Ny-1)-T(1,Ny-2))/Dy2;
      else
        DTDy = (T(1,j+1)-T(1,j-1))/Dy2;
      end
      UU(1,j) =         - Dtor*DTDx;
      VV(1,j) = BCyl(j) - Dtor*DTDy;
   end

%--------------
% Right boundary
%
% Use backward differences
% special treatment for the near-corner nodes
%--------------

    for j=2:Ny
      DTDx = 0.0;
      if(j==2)
        DTDy = (-3.0*T(Nx+1,2)+4.0*T(Nx+1,3)-T(Nx+1,4))/Dy2;
      elseif(j==Ny)
        DTDy =-(-3.0*T(Nx+1,Ny)+4.0*T(Nx+1,Ny-1)-T(Nx+1,Ny-2))/Dy2;
      else
        DTDy = (T(Nx+1,j+1)-T(Nx+1,j-1))/Dy2;
      end
      UU(Nx+1,j) =         - Dtor*DTDx;
      VV(Nx+1,j) = BCyr(j) - Dtor*DTDy;
    end

  end % of Iskip

%-------------------------------------- 
% Compute the wall slip velocity
%
% and modify the boundary conditions
% for the intermediate (star) velocities
%-------------------------------------- 

    slipmax = 0.0;

% top and bottom:

    for i=1:Nx+1
      cor = UU(i,Ny+1)-Ulid(i);
      if(abs(cor)>slipmax) slipmax = cor; end
      BCxt(i) = BCxt(i)-sliprel*cor;
        cor = UU(i,1);
        if(abs(cor)>slipmax) slipmax = cor; end
        BCxb(i) = BCxb(i)-sliprel*cor;
    end

% left and right:

    for j=1:Ny+1
      cor = VV(1,j);
      if(abs(cor)>slipmax) slipmax = cor; end
      BCyl(j) = BCyl(j)-sliprel*cor;
        corr = VV(Nx+1,j);
        if(abs(cor)>slipmax) slipmax = cor; end
        BCyr(j) = BCyr(j)-sliprel*cor;
    end

    slipmax

    if(slipmax<sliptol) break; end


  end % of inner iterations


%-----------------------------------
% Update velocity to the final value
%-----------------------------------

  for j=1:Ny+1
    for i=1:Nx+1
      U(i,j) = UU(i,j);
      V(i,j) = VV(i,j);
    end
  end

%-----------
% Reset time
%-----------

  time = time + Dt

%===========
 end  % of time stepping
%===========

figure
mesh(X,Y,U);
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('u','fontsize',15)
set(gca,'fontsize',15)
box

figure
mesh(X,Y,V);
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('V','fontsize',15)
set(gca,'fontsize',15)
box

figure
mesh(X,Y,T);
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('T','fontsize',15)
set(gca,'fontsize',15)
box
