close all
clear all

%========

N = 16;
a = 1.0;
b =  2.0;

L = b-a;
h = L/N;

for i=1:N+1
 x(i) = a+(i-1)*h;
 f(i) = exp(x(i));
end

%---
% cosine Fourier coefficients
%---

ac = fourier_cos(N,L,f);

%---
% Fourier reproduction
%---

k =2*pi/L;
np = 32*N;
dx = 3.0*L/np;

xx = a-L;

for i=1:np+1
   xh  = xx-a;
   zcm = 0.5*ac(1);
       for p = 1:N
         arg = 0.5*p*k*xh;
         cs  = cos(arg);
         zcm = zcm+ac(p+1)*cs;
       end
   xpl(i)=xx;
   ypl(i)=zcm;
   exa(i)=exp(xx);
   xx = xx+dx;
end

%---
% plotting
%---
 
figure(1)
hold on
axis([0 3 0 12])
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box on

plot(x,f,'ko')
plot(xpl,ypl,'k--')
plot(xpl,exa,'r')
