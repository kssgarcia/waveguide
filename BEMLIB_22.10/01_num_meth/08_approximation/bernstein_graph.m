%---
% plot the bernstein polynomials
%---
clear all
close all
hold on
%---

M=1;

Npl=128;
Dt = 1/Npl;

for k=1:Npl
 t= (k-1)*Dt;
 x(k) = t;

 bern = bernstein (M,t);
 for j=1:M+1
  y(k,j)=bern(M,j);
 end

end

%---

for j=1:M+1
 plot(x,y(:,j))
end

xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box
