%----
% mathieu equation
%---

N=32;

%===========
for rep=1:16
%===========

q=rep-1.0;

A=zeros(N,N);

A(1,1)=0.0; A(1,2)=2*q;
for i=2:N-1
 A(i,i-1)=q;A(i,i)=4*(i-1)^2;A(i,i+1)=q;
end
A(N,N-1)=q;A(N,N)=4*(N-1)^2;

lam=eig(A);

a0(rep)=lam(1);
a2(rep)=lam(2);
a4(rep)=lam(3);
a6(rep)=lam(4);

%---

A=zeros(N,N);

A(1,1)=4.0; A(1,2)=q;
for i=2:N-1
 A(i,i-1)=q;A(i,i)=4*i^2;A(i,i+1)=q;
end
A(N,N-1)=q;A(N,N)=4*N^2;

lam=eig(A);
b2(rep)=lam(1);
b4(rep)=lam(2);
b6(rep)=lam(3);

%---

A=zeros(N,N);

A(1,1)=1+q; A(1,2)=q;
for i=2:N-1
 A(i,i-1)=q;A(i,i)=(2*i-1)^2;A(i,i+1)=q;
end
A(N,N-1)=q;A(N,N)=(2*N-1)^2;
lam=eig(A);
a1(rep)=lam(1);
a3(rep)=lam(2);
a5(rep)=lam(3);

%---

A=zeros(N,N);

A(1,1)=1-q; A(1,2)=q;
for i=2:N-1
 A(i,i-1)=q;A(i,i)=(2*i-1)^2;A(i,i+1)=q;
end
A(N,N-1)=q;A(N,N)=(2*N-1)^2;
lam=eig(A);
b1(rep)=lam(1);
b3(rep)=lam(2);
b5(rep)=lam(3);

xplot(rep)=q;

end

hold on
plot(xplot,a0)
plot(xplot,a1)
plot(xplot,a2)
plot(xplot,a3)
plot(xplot,a4)
plot(xplot,a5)
plot(xplot,a6)

plot(xplot,b1,'--')
plot(xplot,b2,'--')
plot(xplot,b3,'--')
plot(xplot,b4,'--')
plot(xplot,b5,'--')
plot(xplot,b6,'--')

xlabel('q','fontsize',15)
ylabel('a_i  or   b_i','fontsize',15)
set(gca,'fontsize',15)
box
