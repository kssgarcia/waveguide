u=@(x)exp(x.^2);
x=[0,1];
u(x)
u2=@(x)(4*x.^2 + 2).*exp(x.^2);
u2exact=u2(1);
h_list=[0.1, 0.01, 1e-3, 1e-4];
x=1;
for h=h_list
    u2aprox=(u(x+h)-2*u(x)+u(x-h))/h^2;
    error=abs(u2exact-u2aprox);
   %% disp(sprintf('%f %.8f %.8f', h, u2aprox, error))
   disp(sprintf('%e %.8f %.8e', h, u2aprox, error))

end
