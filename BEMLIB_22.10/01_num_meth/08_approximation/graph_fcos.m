clear all
close all

%------
% prepare a graph of the cosine Fourier series
% of the sine function
%------

a = 0.1;
b = 2.45;

nplt =128;

%---
% prepare
%---

L = b-a;
Dx = L/nplt;

figure(1)
hold on
axis([-1, 2, 0, 1.2])
xlabel('(x-a)/L','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box

nf = 1; % number of Fourier terms
        % (will be doubled)

%---
for iloopf = 1:4
%---

nf = 2*nf;

%---
for ip=-1:1

for i=1:nplt+1
 x(i) = a +(i-1)*Dx+ip*L;
 xh = x(i)-a;

 ssm = 0.5;
 for p=2:2:nf
  ssm = ssm - cos(p*pi*xh/L)/(p*p-1);
 end

 smm = 0.0;
 for p=2:2:nf
  smm = smm - p^2*cos(p*pi*xh/L)/(p*p-1);
 end

 xplot(i) = (x(i)-a)/L;
 y(i) = 4.0*ssm/pi;
 w(i) = 4.0*smm/pi;
 z(i) = sin(pi*xh/L);
end

figure(1)
plot(xplot,y,'k--');
%plot(xplot,w,'k--');

if(ip==0)
 plot(xplot,z,'k','linewidth',2);
end

end
%---

%---
end
%---
