hold on

%=====
% plot the Phi function
%=====

a=-1.0;
b= 1.0;

%---
for K=1:4
%---

N=2^K;

step=(b-a)/N;
for i=1:N+1
 xx(i)=a+step*(i-1.0);
 xx(i) = cos((i-0.5)*pi/(N+1));
end

M=128;
step=(b-a)/M;

for i=1:M+1
 x(i)=a+step*(i-1.0);
 prod = 1.0;
 for j=1:N+1
  prod=prod*(x(i)-xx(j));
 end
 y(i)=2^N*prod;
end

if(K==1)
 plot(x,y,':');
elseif(K==2)
 plot(x,y,'--');
elseif(K==3)
 plot(x,y,'-.');
else
 plot(x,y);
end

end

hold on
plot(xx,0,'o','markersize',5,'color','red')
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis([-1 1 -4 4])
axis([-1 1 -1.5 1.5])
box
