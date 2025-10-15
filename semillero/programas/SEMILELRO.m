clear all;

x = linspace(-20,20,5000);
%funcion gaussiana
c=1;
f = @(s) exp(-s.^2);
g = @(x) 0.5*(tanh(5*(x+2)) - tanh(5*(x-2)));
% for t=0:10
% plot(x,g(x + c*t));
% title(sprintf('t=%f', t))
% ylim([0, 1.1])
% pbaspect([4,1,1])
% drawnow
% pause (0.5);
% end

figure(1)
for t=0:20
    plot(x,f(x + 8 - c*t) + g(x - 8 + c*t));
    title(sprintf('t=%f', t))
    ylim([0, 2])
    pbaspect([4,1,1])
    drawnow
    pause (0.5);
end

