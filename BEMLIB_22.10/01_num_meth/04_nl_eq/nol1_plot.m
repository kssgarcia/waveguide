close all
clear all
%----

file2 = fopen('nol1.out')

N   = fscanf(file2,'%d',[1,1])
for i=1:N
 x(i)  = fscanf(file2,'%f',[1,1]);
 y(i)  = fscanf(file2,'%f',[1,1]);
 x(i) = log10(x(i));
end

plot(x,y,'k')
xlabel('log(\mu L_p/a)','fontsize',15)
ylabel('ka','fontsize',15)
set(gca,'fontsize',15)
axis([-5 5, 0 2.5])



