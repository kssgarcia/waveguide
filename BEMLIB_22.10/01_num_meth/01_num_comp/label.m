clc; clear all; close all;

%======================
% sort and label a list
%======================

%x = [9, 8, 12, -2, 8, 3, 4, 6, 3, 0, 3];
x = [9, 12, -2, 8, 3];
n = length(x);

fprintf('initial vector:');
x'

%---
% launch
%---

m = [1:n];  % initialize

for k=n-1:-1:1
    Istop=1;
    for i=1:k
      if(x(i)>x(i+1))
          save = x(i);
          x(i) = x(i+1);
        x(i+1) = save;
         Istop = 0;
         isave = m(i);
          m(i) = m(i+1);
        m(i+1) = isave;
      end
    end
    if(Istop==1) break; end
end

%---
% done
%---

fprintf('sorted x:');
x'
fprintf('sorted indices:');
m'
