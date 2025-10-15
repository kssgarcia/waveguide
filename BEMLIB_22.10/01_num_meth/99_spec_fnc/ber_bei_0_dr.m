clear all
close all

%---
% plot the Kelvin functions ber_0 and bei_0
%---

Iopt=0;

for i=1:32
 x(i) = 0.01+8*(i-1)/32;
 X = x(i);
 [Y,Z,ber_0_p,bei_0_p] = ber_bei_0 (Iopt,X);
 y(i)=Y;
 z(i)=Z;
end

hold on
plot(x,y,'k')
plot(x,z,'k--')
xlabel('z','fontsize',15)
ylabel('ber_0(z),   bei_0(z)','fontsize',15)
set(gca,'fontsize',15)
axis([0 8 -10 10])
box

