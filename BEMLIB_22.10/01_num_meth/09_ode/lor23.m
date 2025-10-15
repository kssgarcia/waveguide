clear all

%=================================
% integration of the lorenz system
% using a matlab solver
%=================================

tspan = [0, 100];

y0 = [0.2 0.2 0.1];

[t, y] = ode23(@lor_fnc,tspan,y0);

%plot(t,y(:,1))
hold on
plot3(y(1,1),y(1,2),y(1,3),'o')
plot3(y(:,1),y(:,2),y(:,3))
xlabel('x_1','fontsize',15)
ylabel('x_2','fontsize',15)
zlabel('x_3','fontsize',15)
set(gca,'fontsize',15)

