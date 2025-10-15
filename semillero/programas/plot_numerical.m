clear all
load u.mat ut x dt

[nx, nt] = size(ut);
nx = nx + 1;
vec = zeros(1,nt);
ut = [vec; ut; vec];

% figure(1); clf
% for i=1:nt
%     t = (i-1)*dt
%     plot(x, ut(:,i), 'LineWidth', 2)
%     ylim([-1.1, 1.1])
%     t=text(0.05, 1, sprintf('t=%.2f', t), 'HorizontalAlignment','Left', 'VerticalAlignment', 'Top');
%     t.FontWeight = 'bold';

%     drawnow
%     pause(0.5)
% end

figure(1);clf
x_text = 0.95;

subplot(2,2,1)
i=1;
t=(i-1)*dt;
plot(x, ut(:,i), 'LineWidth', 2)
ylim([-1.1, 1.1])
ylabel('u')
t=text(x_text, 1, sprintf('t=%.2f', t), 'HorizontalAlignment','right', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';

subplot(2,2,2)
i=11;
t=(i-1)*dt;
plot(x, ut(:,i), 'LineWidth', 2)
ylim([-1.1, 1.1])
t=text(x_text, 1, sprintf('t=%.2f', t), 'HorizontalAlignment','right', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';

subplot(2,2,3)
i=21;
t=(i-1)*dt;
plot(x, ut(:,i), 'LineWidth', 2)
ylim([-1.1, 1.1])
ylabel('u')
xlabel('x')
t=text(x_text, 1, sprintf('t=%.2f', t), 'HorizontalAlignment','right', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';

subplot(2,2,4)
i=33;
t=(i-1)*dt;
plot(x, ut(:,i), 'LineWidth', 2)
ylim([-1.1, 1.1])
xlabel('x')
t=text(x_text, 1, sprintf('t=%.2f', t), 'HorizontalAlignment','right', 'VerticalAlignment', 'Top');
t.FontWeight = 'bold';

outfile = '../poster_ondas/onda_numerico';
print('-depsc', outfile)
[status,output]=system(['epstopdf ',outfile,'.eps']);
