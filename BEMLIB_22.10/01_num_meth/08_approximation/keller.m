%========
% plot the data contained in file: keller.dat
%========

file2 = fopen('keller.dat')

for i=1:2*2*512
  x(i) = fscanf(file2,'%f',[1,1]);
  u(i) = fscanf(file2,'%f',[1,1]);
  v(i) = fscanf(file2,'%f',[1,1]);
  T(i) = fscanf(file2,'%f',[1,1]);
end

plot(x,u,'.')
xlabel('time','fontsize',15)
ylabel('velocity','fontsize',15)
set(gca,'fontsize',15)
axis([0 2000 1.75 2.75])
