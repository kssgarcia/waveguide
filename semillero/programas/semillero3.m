
k =0 ;
x = linspace(0,1,500) ;
u =0;
nk=50;
for k=0:nk 
    ter = (4/pi^2) * ((-1)^k/(2*k+1)^2) * sin((2*k+1)*pi*x); 
    u=u+ter ;
end

fx = zeros(length(x));
for i=1:length(x); fx(i)=f(x(i)); end

figure(1);clf
plot (x,fx);
hold on
plot (x,u);
hold off

% for t=0.1*(0:20)
% figure(2);clf
% plot(x,fx);
% hold on
% plot(x,fun(t,x));
% hold off
% title(sprintf('t=%f',t))
% ylim([-0.8, 0.8])
% drawnow
% pause(0.2)
% end

figure(1)
subplot(2,2,1)
t=0;
plot(x,fun(t,x), 'LineWidth', 2)
ylim([-0.8, 0.8])
ylabel('u')
t = text(0.05, 0.7, sprintf('t=%.1f', t), 'HorizontalAlignment','Left', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';

subplot(2,2,2)
t=0.2;
plot(x,fun(t,x), 'LineWidth', 2)
ylim([-0.8, 0.8])
% ylabel('u')
t = text(0.05, 0.7, sprintf('t=%.1f', t), 'HorizontalAlignment','Left', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';

subplot(2,2,3)
t=0.4;
plot(x,fun(t,x), 'LineWidth', 2)
ylim([-0.8, 0.8])
ylabel('u')
xlabel('x')
t = text(0.05, 0.7, sprintf('t=%.1f', t), 'HorizontalAlignment','Left', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';

subplot(2,2,4)
t=1;
plot(x,fun(t,x), 'LineWidth', 2)
ylim([-0.8, 0.8])
xlabel('x')
t = text(0.05, 0.7, sprintf('t=%.1f', t), 'HorizontalAlignment','Left', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';

outfile = '../poster_ondas/onda_triangulo';
print('-depsc', outfile)
[status,output]=system(['epstopdf ',outfile,'.eps']);

