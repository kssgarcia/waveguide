clear all;

x1 = linspace(-20,0,5000);
x2 = linspace(0,20,5000);
V1=2;
V2=1;
A=(V2-V1)/(V2+V1);
B=(2*V2)/(V2+V1);
% A=1;
% B=0;
%funcion gaussiana
R1 = @(s) exp(-(s+8).^2);
L1= @(s) A*R1(-s);
R2= @(s) B*R1((V1/V2)*s);
figure(1);
clf;
% for t=0:20
%     plot(x1,R1(x1 - V1*t) + L1(x1 + V1*t), 'LineWidth', 2);
%     hold on
%     plot(x2,R2(x2 - V2*t), 'LineWidth', 2)
%     plot([0,0],[-1,1.1],'k')
%     hold off 
%     title(sprintf('t=%f', t))
%     ylim([-1, 1.1])
%     xlabel('x')
%     pbaspect([4,1,1])
%     drawnow
%     pause (0.1);
% end

figure(1)
subplot(3,1,1)
t = 0;
plot(x1,R1(x1 - V1*t) + L1(x1 + V1*t), 'LineWidth', 2);
hold on
plot(x2,R2(x2 - V2*t), 'LineWidth', 2)
plot([0,0],[-1,1.1],'k')
hold off
% title(sprintf('t=%d', t))
t = text(18, 1, sprintf('t=%d', t), 'HorizontalAlignment','right', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';
ylim([-1, 1.1])
% xlabel('x')
pbaspect([4,1,1])

subplot(3,1,2)
t = 4;
plot(x1,R1(x1 - V1*t) + L1(x1 + V1*t), 'LineWidth', 2);
hold on
plot(x2,R2(x2 - V2*t), 'LineWidth', 2)
plot([0,0],[-1,1.1],'k')
hold off
% title(sprintf('t=%d', t))
t = text(18, 1, sprintf('t=%d', t), 'HorizontalAlignment','right', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';
ylim([-1, 1.1])
% xlabel('x')
pbaspect([4,1,1])

subplot(3,1,3)
t = 8;
plot(x1,R1(x1 - V1*t) + L1(x1 + V1*t), 'LineWidth', 2);
hold on
plot(x2,R2(x2 - V2*t), 'LineWidth', 2)
plot([0,0],[-1,1.1],'k')
hold off
% title(sprintf('t=%d', t))
t = text(18, 1, sprintf('t=%d', t), 'HorizontalAlignment','right', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';
ylim([-1, 1.1])
xlabel('x')
pbaspect([4,1,1])

outfile = '../poster_ondas/reflex_trans_1';
print('-depsc', outfile)
[status,output]=system(['epstopdf ',outfile,'.eps']);
