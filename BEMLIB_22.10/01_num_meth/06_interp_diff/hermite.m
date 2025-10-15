clear all
close all
hold on

%---
% hermite intepolation functions
%---

N=4;
x=[0.0, 0.20, 0.5, 0.80, 1.0];
m=[3, 0, 2, 2, 1];
m=[3, 1, 2, 2, 1];

nplt=256;
step=(x(N+1)-x(1))/nplt;

%---
% compute the matrix A(i,k,s)
%---

A = hermite_A(N,x,m);

%---
% run over abscissas
%for i=1:N+1;
%for i=2:2; 
for i=1:1; 
%---

for j=1:nplt+1

  xint = x(1)+step*(j-1.0);
  xplot(j)=xint;

%---
  k=m(i)+1;
  q(i,k)= hermite_r(xint,N,x,m,i,k);

  yplot(j,k)=q(i,k);

  for p=1:m(i)
    k=m(i)+1-p;
    q(i,k) = hermite_r(xint,N,x,m,i,k);
    for s=k+1:m(i)+1;
     q(i,k) = q(i,k)-A(i,k,s)*q(i,s);
    end
    yplot(j,k)=q(i,k);
  end
%---

end

 for k=1:m(i)+1
  if(k==1)
   plot(xplot,yplot(:,k),'r');
  elseif(k==2)
%   plot(xplot,yplot(:,k),'r');
  elseif(k==3)
%   plot(xplot,yplot(:,k),'r');
  elseif(k==4)
%   plot(xplot,yplot(:,k),'r');
  end
 end

%---
end
%---

plot(x,zeros(N+1),'ro')

xlabel('x','fontsize',15)
ylabel('q(1,1)','fontsize',15)
set(gca,'fontsize',15)
box

