clear all

%=======================================
% Solution of the boundary-value problem
% for the Blasius equation
%
% will solve the boundary-value problem
% with newton iterations
%
% will solve an extended system of six ODEs 
%=======================================

 tspan = [0, 10];

%---
% initial guess
%----

 y0 = [0.0 0.0 0.0 0.0 0.0 1.0];

%---
% iterate
%---

 for iter=1:10

  [t, y]= ode23(@blas_fnc,tspan,y0);

  N = size(y);
  residual = y(N(1),2)-1.0;

  y0 = [0.0; ...
        0.0; ...
        y0(3)-residual/y(N(1),5);...
        0.0; ...
        0.0; ...
        1.0];

  if(abs(residual)<0.00000001) break; end

 end

%---
% plot the solution
%---

plot(t,y(:,2))
xlabel('\eta','fontsize',15)
ylabel('f','fontsize',15)
set(gca,'fontsize',15)
box
