%==========================================
% CODE sdl
%
% Finite-element code for steady one-dimensional 
% diffusion with linear elements
%==========================================

%-----------
% input data
%-----------

L=1.0; k=1.0; q0=-1.0; fL=0.0;

% q0=-1; fL=exp(1.0);  % CASE B

ne=10; ratio=1.0;
ne=512; ratio=1.0;
ne=16; ratio=2.0;

%----------------
% grid generation
%----------------

xe = elm_line1 (0,L,ne,ratio);

%-------------------
% specify the source
%-------------------

for i=1:ne+1
% s(i) = - exp(xe(i));   % CASE B
  s(i) = 0.0;
  s(i) = 1.0;
  s(i) = 10.0*exp(-5.0*xe(i)^2/L^2);
end

%-----------------
% compact assembly
%-----------------

[at,bt,ct,b] = sdl_sys (ne,xe,q0,fL,k,s);

%--------------
% linear solver
%--------------

f = thomas (ne,at,bt,ct,b);

f(ne+1) = fL;

%-----
% plot
%-----

plot(xe, f,'-o');
xlabel('x','fontsize',15);
ylabel('f','fontsize',15);
set(gca,'fontsize',15)

%-----
% done
%-----
