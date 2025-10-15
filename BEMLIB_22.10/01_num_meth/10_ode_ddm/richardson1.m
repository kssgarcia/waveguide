%----
% graph of rho vs omega
% for several N
%---

hold on
N=16;

omopt=0.5/(sin(0.5*pi/(N+1))^2+sin(0.5*N*pi/(N+1))^2);
omax=0.5/sin(0.5*N*pi/(N+1))^2;

for times=1:7
 if(times==1)
  omega=omopt;
 elseif(times==2)
  omega=0.2*omax;
 elseif(times==3)
  omega=0.4*omax;
 elseif(times==4)
  omega=0.6*omax;
 elseif(times==5)
  omega=0.8*omax;
 elseif(times==6)
  omega=omax;
 elseif(times==7)
  omega=0.0;
 end

 for i=1:N
  l(i)=1-4*omega*sin(0.5*i*pi/(N+1))^2;
 end
 x=linspace(1,N,N)

  if(times==1)
   plot(x,l,'--d')
  else
   plot(x,l,':o')
  end

end

%---
set(gca,'fontsize',15)
xlabel('i','fontsize',15)
ylabel('\lambda','fontsize',15)
axis([0 N -1 1]);
box
%---

