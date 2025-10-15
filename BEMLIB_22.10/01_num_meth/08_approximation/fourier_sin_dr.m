close all
clear all

%========
% Driver for evaluating sine Fourier coefficients
%========

N = 16;
a = 1.0;
b = 2.0;

L = b-a;
h = L/N;

for i=1:N+1
 x(i) = a+(i-1)*h;
 f(i) = exp(x(i));
end

[bs] = fourier_sin(N,L,f);

k = 2*pi/L;
np = 32*N;
dx = 3.0*L/np;

xx = a-L;

for i=1:np+1
   xh  = xx-a;
   zcm = 0.0;
       for p = 1:N-1
         arg = 0.5*p*k*xh;
         sn  = sin(arg);
         zcm = zcm+bs(p+1)*sn;
       end
   xpl(i) = xx;
   ypl(i) = zcm;
   exa(i) = exp(xx);
   xx = xx+dx;
end
 
%==========
% plotting
%==========

figure(1)
hold on
axis([0 3 -15 15])
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
plot(x,f,'ko')
plot(xpl,ypl,'k--')
plot(xpl,exa,'k')

