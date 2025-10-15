clear all
%===
% integration of the lorenze system
%===

 tspan=[0, 10];
 y0=[0];
 [t, y]= ode23(@fstiff,tspan,y0);

%plot(t,y(:,1))
hold on
plot(t,y,'o')
xlabel('t','fontsize',15)
ylabel('x','fontsize',15)
set(gca,'fontsize',15)

