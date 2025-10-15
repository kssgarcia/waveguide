% x = linspace(0,1,500) 
f=@(x)exp(-400*(x-0.3).^2) % condicion inicial 
N=50;
L=1;
b=zeros(N,1);
for n=1:N 
    b(n)= 2*integral(@(x) f(x).*sin(n*pi*x/L),0,L); 
    
end

n=90;
dx=L/n;
dt=0.01;
c=1;
sigma2=((c*dt)/dx)^2;
B=zeros(n-1);
for i=2: n-2
    B(i,i)=2*(1-sigma2);
    B(i,i+1)=sigma2;
    B(i,i-1)=sigma2;
end
B(1,1)=2*(1-sigma2);
B(n-1,n-1)=2*(1-sigma2);
B(1,2)=sigma2;
B(n-1,n-2)=sigma2;
x=linspace(0,L,n+1);
u0=f(x);
figure(1)
plot(x,u0)

u0 = u0(:);
u0 = u0(2:end-1);

u1=0.5*B*u0;

nt=100;
ut = zeros(n-1,nt+1);
ut(:,1) = u0;
for i= 1:nt
    u=B*u1-u0;
    ut(:, i+1) = u;
    
    u0=u1;
    u1=u;
    figure(1)
    t=(i+1)*dt;
    plot(x,[0;u;0])
    hold on 
    plot(x,fourier(t,x,b),"--")
    hold off
    ylim([-1.1, 1.1])
    title(sprintf('t=%f', t))
    drawnow
    pause(0.1)
end

save u.mat ut x dt

