clear all
close all

%---
% animate an exact solution of Burger's equation
%---

c=0.50;

N=2*2*32;
a=-2;
b= 2;
Dx=(b-a)/N;
U=1.0;
Dt=c*Dx/U;

for i=1:N+1
 x(i) = a+(i-1)*Dx;
 f(i) = exp(-x(i)*x(i));
 f(i) = - cosh(x(i))/( sinh(x(i))+ 1);
end

fnew=f;
t=0.0;

for step=1:30000
t=t+Dt

for i=2:N
 fnew(i) = - cosh(x(i))/( sinh(x(i))+exp(-t));
end

f(2:N)=fnew(2:N);

if(step==1)
  Handle1 = plot(x,f,'-');
  set(Handle1, 'erasemode', 'xor');
  set(gca,'fontsize',15)
  axis([a b -10.0 10.0])
  xlabel('x','fontsize',15)
  ylabel('f','fontsize',15)
else
  set(Handle1,'XData',x,'YData',f);
  pause(0.1)
  drawnow
end

end
