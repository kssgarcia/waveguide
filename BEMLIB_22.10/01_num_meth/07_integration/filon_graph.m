clear all
close all

%------
% prepare graphs of a Filon function
%------

N=2*2*128;
a=0.0;
b=3.0;

k=100.0;

Dx=(b-a)/N;

hold on

%---
for repeat=1:1

for i=1:N+1
 x(i) = 0.001+ a +(i-1)*Dx;
 y(i) = exp(-x(i)) *cos(k*x(i));
end
if(repeat==1)
plot(x,y);
else
end

end
%---

%axis([0, 1, 0, 10])
axis([0, 3, -1, 1])
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box
