format long

%---
% sum to pi
%---

N = 1;

for j=1:7
 sum = 0;
 for i=0:N
   sum = sum + 1.0/16.0^i * ( 4/(8*i + 1)  ...
             - 2.0/(8.0*i+4.0) - 1/(8.0*i+5.0) - 1/(8.0*i+6.0) );
 end
 pie(j) = sum;
 N = 2*N;
end

pie'

