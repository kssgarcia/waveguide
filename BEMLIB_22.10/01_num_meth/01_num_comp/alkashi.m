clear all
close all

%===============
% compute the power b=a^p
% by the Al-kashi algorithm
%===============

a = 1.10;

%---
for iloop=1:7
%---

initime = cputime;

n = 128;
n = 2^iloop;

b = 1.0;

iops = 0;

while(n>0)

  if(mod(n,2)==0)
   n = n/2;
   iops = iops+1;
   a = a*a;
   iops = iops+1;
  else
   n = n-1;
   b = b*a;
   iops = iops+1;
  end

end

b
[iloop, iops]

fintime = cputime-initime
plotx(iloop) = iloop;
ploty(iloop) = iops;

end

figure(1)
hold on
box on
xlabel('ln_2 n','fontsize',15)
ylabel('operation count (n_c)','fontsize',15)
set(gca,'fontsize',15)
plot(plotx,ploty,'k-o')
axis([1 7 0 15])

