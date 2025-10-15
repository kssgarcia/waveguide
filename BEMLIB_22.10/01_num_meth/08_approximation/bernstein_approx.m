%---
% approximation using the bernstein polynomials
%---
clear all
close all
hold on
%---

Npl=128;
Dt = 1/Npl;

%---
for M=2:2:20

  for j=1:M+1
  xdata(j)=(j-1)/M;
  h(j) = exp(xdata(j));
  xx = xdata(j);
  h(j) = 1/(1+25*xx*xx);
  end


 for k=1:Npl
 t= (k-1)*Dt;
 x(k) = t;
 bern = bernstein (M,t);
 y(k) = 0.0;
 for j=1:M+1
  y(k)=y(k)+h(j)*bern(M,j);
 end

 end

 if(M==2)
 plot(x,y,'--')
 else
 plot(x,y)
 end

end

plot(xdata,h,'o')

xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box
