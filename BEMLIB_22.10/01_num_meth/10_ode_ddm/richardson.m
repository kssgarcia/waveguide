%----
% graph of the spectral radius vs omega
% for several N
%---

hold on
Npl=128;

N=1;

for nn=1:5

omax=0.5/sin(0.5*N*pi/(N+1))^2;

for m=1:Npl+1

 omega(m)=(m-1)*omax/Npl;
 maxx=0.0;
 for i=1:N
  l=1-4*omega(m)*sin(0.5*i*pi/(N+1))^2;
  if(abs(l)>maxx)
   maxx=abs(l);
  end
  mmxx(m)=maxx;
 end
end

N=2*N
if(nn==1)
plot(omega/omax,mmxx,'--')
else
plot(omega/omax,mmxx,'-')
end

end

%---
set(gca,'fontsize',15)
xlabel('\omega/\omega_{max}','fontsize',15)
ylabel('\rho','fontsize',15)
box
%---

