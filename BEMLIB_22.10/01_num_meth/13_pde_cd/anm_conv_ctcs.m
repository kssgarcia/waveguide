clear all
close all

%---
% animation of the solution of the
% linear convection equation with ctcs
%-----

c=0.50;

N=2*2*2*2*32;
a=-10;
b= 30;
Dx=(b-a)/N;
U=1.0;
Dt=c*Dx/U;

for i=1:N+1
 x(i)=a+(i-1)*Dx;
 f(i) = exp(-x(i)*x(i));
end
fold=f;
fnew=f;

for i=2:N
 fnew(i) = 0.5*c*(1+c)*f(i-1)+(1-c*c)*f(i) ...
                 +0.5*c*(c-1)*f(i+1);
end
f=fnew;

time(1)=0;
ggg(1)=f(N/2);

for step=1:300

for i=2:N
 fctcs(i) = fold(i)+c*f(i-1)-c*f(i+1);
end

fold(2:N)=f(2:N);
f(2:N)=fctcs(2:N);

time(step+1)=time(step)+Dt;
ggg(step+1)=f(N/2);


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

figure
plot(time,ggg)
