function u = fourier(t,x,b)
N=length(b);
L=1;
c=1; 
u=0;
for n=1:N
    u=u+b(n)*cos(n*pi*c*t/L)*sin(n*pi*x/L);
end