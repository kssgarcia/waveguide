clear all
close all

%=============================
% analysis of the SOR method 
% for a model problem with N=2
%=============================

Nomega=128;
domega=2.0/Nomega;

%---
for aloop=1:8
%---

 if(aloop==1)
  a=0.0;
 elseif(aloop==2)
  a=0.2;
 elseif(aloop==3)
  a=0.5;
 elseif(aloop==4)
  a=0.8;
 elseif(aloop==5)
  a=0.9;
 elseif(aloop==6)
  a=0.98;
 elseif(aloop==7)
  a=1.0;
 elseif(aloop==8)
  a=1.02;
 end

 for i=1:Nomega+1

  omega = 0.0001+(i-1.0)*domega;
  Omega = a^2*omega^2;
  discr = Omega*(Omega+4*(1-omega));

  if(discr==0)
   labs1(i) = abs(1.0-omega);
   labs2(i) = labs1(i);
  elseif(discr>0)
   tmp = sqrt(discr);
   labs1(i) = abs(1.0-omega + 0.5*(Omega+tmp));
   labs2(i) = abs(1.0-omega + 0.5*(Omega-tmp));
  elseif(discr<0)
   tmp = sqrt(-discr);
   real =  1.0-omega + 0.5*Omega;
   imag =  0.5*tmp;
   labs1(i) = sqrt(real^2+imag^2);
   labs2(i) = labs1(i);
  end

  omg(i)=omega;

 end

 if(aloop==1)
  plot(omg,labs1,'--');
% plot(omg,labs2,'--');
 else
  plot(omg,labs1,'b');
% plot(omg,labs2,'b');
 end
 hold on

%---
end
%---

xlabel('\omega','fontsize',15)
ylabel('|\lambda_2|','fontsize',15)
ylabel('|\lambda_1|','fontsize',15)
set(gca,'fontsize',15)

axis([0 2 0 1.4])
