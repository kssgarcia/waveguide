clear all
close all

%------
% prepare graphs of the Fourier series
% of the cosine function
%------

N=2*128;
a=-pi;
b=2*pi;

Dx=(b-a)/N;

hold on

%---
for repeat=1:1

for i=1:N+1
 x(i) = a +(i-1)*Dx;
 y(i) = 0.0;
 for p=2:2:50
  y(i) = y(i) + p*sin(p*x(i))/(p*p-1);
 end
 y(i)=4.0*y(i)/pi;
 z(i)=cos(x(i));
end
if(repeat==1)
plot(x/pi,y,'--');
plot(x/pi,z);
else
end

end
%---

%axis([0, 1, 0, 10])
%axis([0, 3, 0, 1.2])
xlabel('x/\pi','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box
