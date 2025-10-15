clear all
close all

%-----
% animation of the solution of the convection
% with various methods
%------

c=1.50;

N=2*2*2*2*32;
a=-10;
b= 30;
L=b-a;
Dx=L/N;
U=1.0;
Dt=c*Dx/U;
nstep=30000;

%-----
% grid and initial condition
%------

for i=1:N+1
 x(i)=a+(i-1)*Dx;
 f(i) = exp(-x(i)*x(i));
 f(i) = 0.5D0*cos(2.0*pi*x(i)/L);
 f(i) = - cosh(x(i))/( sinh(x(i))+ 1);
end

fnew=f;

for step=1:nstep

t=t+Dt;

for i=2:N
 fnew(i) = 0.5*c*f(i-1)+f(i)-0.5*c*f(i+1);
 fnew(i) = 0.5*(1+c)*f(i-1)+0.5*(1-c)*f(i+1); % Lax
 fnew(i) = 0.5*c*(1+c)*f(i-1)+(1-c*c)*f(i) ...
                 +0.5*c*(c-1)*f(i+1);
 fnew(i) = c*f(i-1)+(1-c)*f(i); % FTBS
 f(i) = - cosh(x(i))/( sinh(x(i))+exp(-t));
end

f(2:N)=fnew(2:N);

if(step==1)
  Handle1 = plot(x,f,'o-');
  set(Handle1, 'erasemode', 'xor');
  set(gca,'fontsize',15)
  axis([a b -0.2 1.2])
  xlabel('x','fontsize',15)
  ylabel('f','fontsize',15)
else
  set(Handle1,'XData',x,'YData',f);
  pause(0.1)
  drawnow
end

end
