clear all
close all

%---
for p=1:10
%---

N = 2^p;

for i=1:N
 for j=1:N
  A(i,j) = rand;
 end
 b(i) = rand;
end

initime = cputime

sol = b/A';

exectime = cputime-initime

plotx(p) = log10(N);;
ploty(p) = log10(exectime);

%---
end
%---

figure(1)
hold on
hold on
box on
xlabel('log(N)','fontsize',15)
ylabel('log(CPU time)','fontsize',15)
set(gca,'fontsize',15)
plot(plotx,ploty,'k-o')
plot([2 3.5],[-3, 1.5],'r-')
axis([1.5 3.5 -3.0 1.5])


