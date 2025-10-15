clear all
close all

%=====================
% particle solution of 
% the convection equation
%=====================

N  = 2*2*64;
a  =-10;
b  = 10;
Dt = 0.0010;

%---
% prepare
%---

Dx = (b-a)/N;

%---
% initial condition
%---

for i=1:N+1
 x(i) = a+(i-1)*Dx;

 F(i) = tanh(x(i));
 F(i) = exp(-x(i)*x(i));

 v(i) = F(i)*F(i);
 v(i) = 1.0;
 v(i) = tanh(x(i));
 v(i) = F(i); % Burgers

end

time = 0.0;

%-----
for istep=1:30000
%-----

 x = x + Dt*v;

 time = time + Dt;

 if(istep==1)
  Handle1 = plot(x,F,'ko-');
  set(Handle1, 'erasemode','xor');
  axis([-2 2 -0.0 1.2])
% axis([-2 2 -1.2 1.2])
  xlabel('x','fontsize',15)
  ylabel('f','fontsize',15)
  set(gca,'fontsize',15)
 else
  set(Handle1,'XData',x,'YData',F);
  pause(0.01)
  drawnow
 end

%-----
end
%-----
