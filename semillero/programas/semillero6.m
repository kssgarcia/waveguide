n =0 ;
x = linspace(0,1,500) ;
u =0;
N=50;
L=1;
b=zeros(N,1);
for n=1:N 
    b(n)= integral(@(x) f(x).*sin(n*pi*x/L),0,L); 
    
end


% c=1; 
% u=0;
% t=0;
% for n=1:N
%     u=u+b(n)*cos(n*pi*c*t/L)*sin(n*pi*x/L);
% end

figure(1)
for i=1:100
    t=(i-1)*dt;
    figure(1)
    plot(x,fourier(t,x,b))
    ylim([-1.1, 1.1])
    title(sprintf('t=%f', t))
    drawnow
    pause(0.5)
    grid on
end
% fx = zeros(length(x));
% for i=1:length(x); fx(i)=f(x(i)); end
% 
% figure(1);clf
% plot (x,fx);
% hold on
% plot (x,u);
% hold off
% 
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
