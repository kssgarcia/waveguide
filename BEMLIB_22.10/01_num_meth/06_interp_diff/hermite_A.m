function A = hermite_A(N,x,m)

%---
% derivatives or r
%---

eps=0.0001;

for i=1:N+1
 for k=1:m(i)

  for s=1:m(i)+1

%---
   if(s==1)
    xint=x(i);
    r00 = hermite_r(xint,N,x,m,i,k);
    A(i,k,s)=r00;
%---
   elseif(s==2)
    xint=x(i)+eps;
    r01 = hermite_r(xint,N,x,m,i,k);
    xint=x(i)-eps;
    r11 = hermite_r(xint,N,x,m,i,k);
    A(i,k,s)=(r01-r11)/(2.0*eps);
%---
   elseif(s==3)
    A(i,k,s)=(r01-2.0*r00+r11)/eps^2;
%---
   elseif(s==4)
    xint=x(i)+2*eps;
    r02 = hermite_r(xint,N,x,m,i,k);
    xint=x(i)-2*eps;
    r12 = hermite_r(xint,N,x,m,i,k);
    A(i,k,s)=(r02-2.0*r01+2.0*r11-r12)/(2.0*eps^3);
   end
%---
  end
  A(i,k,k)=1.0;

 end
end
