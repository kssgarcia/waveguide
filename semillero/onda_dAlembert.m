L = 5;
x = linspace(-L,3*L,500);

f = @(x) exp(-x.^2);
g = @(x) 0.5*(tanh(5*(x+2)) - tanh(5*(x-2)));

c = 1;
t = 0:10;
% onda que viaja a la derecha
figure(1)
for t = 0:0%10
    plot(x,g(x-c*t))
    pbaspect([4,1,1])
    ylim([0,1.1])
    title(sprintf('t=%f',t))
    drawnow
    pause(0.5)
end

% onda que viaja hacia la izquierda
figure(2)
x0 = 2.5*L;
for t = 0:0%10
    plot(x, f(x - x0 + c*t))
    pbaspect([4,1,1])
    ylim([0,1.1])
    title(sprintf('t=%f',t))
    drawnow
    pause(0.5)
end

% superposición de ondas
figure(3)
x0 = 2.5*L;
for t = 0:0%:0.5:10
    plot(x, f(x - c*t) + g(x - x0 + c*t))
    pbaspect([4,1,1])
    ylim([0,2.1])
    title(sprintf('t=%f',t))
    drawnow
    pause(0.5)
end

% Reflexión y transmisión
v1 = 2;
v2 = 1;
A = (v2-v1)/(v2+v1);
B = 2*v2/(v2+v1);

R1 = @(s) exp(-(s+5).^2)
L1 = @(s) A*R1(-s)
R2 = @(s) B*R1(v1/v2*s)

x1 = linspace(-6,0,500);
x2 = linspace(0,6,500);

figure(4); clf

for t=0:0.5:5
    plot([0,0],[-1,1],'k')
    hold on
    plot(x1, R1(x1 - v1*t) + L1(x1 + v1*t))
    plot(x2, R2(x2 - v2*t))
    hold off
    title(sprintf('t=%f',t))
    xlim([-6,6])
    ylim([-1.1,1.1])
    pbaspect([4,1,1])
    drawnow
    pause(0.5)
end




