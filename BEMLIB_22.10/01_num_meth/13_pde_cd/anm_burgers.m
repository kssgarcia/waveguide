clear all
close all

%-----
% animation of the solution to burgers equation
% using different finite-difference methods
%------

r=0.50;

N=2*2*32;
a=-2;b= 2;
Dx=(b-a)/N;

Dt=r*Dx;

method=1;  % Lax
method=2;  % Lax_Wendroff

%---
% initial condition
%----

for i=1:N+1
 x(i)=a+(i-1)*Dx;
 f(i) = exp(-x(i)*x(i));
 q(i) = 0.5*f(i)^2;
end


%---------
for step=1:30000

for i=2:N
 if(method==1)%--- Lax
  fnew(i) = 0.5*(1+r*f(i-1))*f(i-1)+0.5*(1-r*f(i+1))*f(i+1);
 elseif(method==2)%--- Lax-Wendroff
   A = f(i-1)+f(i);
   B = f(i-1)+2*f(i)+f(i+1);
   C = f(i)+f(i+1);
   fnew(i) = f(i) + 0.5*r*(q(i-1)-q(i+1)) + 0.25*r*r*(A*q(i-1)-B*q(i)+C*q(i+1));
 end
%---
end

f(2:N)=fnew(2:N);

for i=1:N+1
 q(i) = 0.5*f(i)^2;
end

if(step==1)
  Handle1 = plot(x,f,'o-');
  set(Handle1, 'erasemode', 'xor');
  set(gca,'fontsize',15)
  axis([a b 0.0 1.2])
  xlabel('x','fontsize',15)
  ylabel('f','fontsize',15)
else
  set(Handle1,'XData',x,'YData',f);
  pause(0.1)
  drawnow
end

end
