%=======================
% plot a sine line in 3D
%=======================

N = 64;
dp = 2.0*pi/N;

for i=1:65
 phase = (i-1.0)*dp;
 x(i) = 0.8*phase;
 y(i) = 0.2*phase;
 z(i) = sin(phase);
end

figure(1)
plot3(x,y,z,'k');
xlabel('x','fontsize',13)
ylabel('y','fontsize',13)
zlabel('z','fontsize',13)
set(gca,'fontsize',13)
axis square
box on
