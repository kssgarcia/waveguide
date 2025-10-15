clear all

%======================
% bubble sort algorithm
%======================

%---
% random data
%---

N = 16;

 for i=1:N
  x(i)=rand;
 end

fprintf('%4.2f ', x);

%clear x;
%x = randperm(N);
%fprintf('%d ', x);

fprintf('\n');

%---

for k=N-1:-1:1  % number of comparisons

 Istop = 1;

 for i=1:k
   if(x(i)>x(i+1))
     save = x(i);
     x(i) = x(i+1);
     x(i+1) = save;
     Istop = 0;
   end
 end

 if(Istop==1) break; end

end

%---

%fprintf('%d ', x);
fprintf('%4.2f ', x);
fprintf('\n');
