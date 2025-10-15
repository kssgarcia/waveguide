c=-8;
x=200.0;
for k=1:30
% xnew=x-1/3*(x*x*x-c)/x^2
 xnew=sqrt(1+x)
 x = xnew;
end
