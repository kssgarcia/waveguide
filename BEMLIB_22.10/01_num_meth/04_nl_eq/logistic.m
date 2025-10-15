hold on

%---
% visualization of logistic mapping
%---

lambda= 2.8:0.01:4;

for i=1:length(lambda)
    x(1) = 0.8;
    for j = 1:200
        x(j+1) = lambda(i)*x(j)*(1.0-x(j));
    end
    plot (lambda(i),x(100:200),'k')
end

axis([2.8 4 0 1])
xlabel('\lambda','fontsize',15)
ylabel('x','fontsize',15)
set(gca,'fontsize',15)
box

