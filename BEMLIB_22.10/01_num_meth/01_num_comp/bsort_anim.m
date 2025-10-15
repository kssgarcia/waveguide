clear all
close all

%=========================
% Animation of bubble sort
%
% Roger Jang, 980703
%=========================

%---
% generate the data
%---

data_n = 64;
x = 1:data_n;
y = rand(1, data_n);
y = randperm(length(x));

%---
% plotting bar
%---

r=0.3;
rect = [-r -r r r] + sqrt(-1)*[0 1 1 0];

%---
% construct the patch array
%---

patchH = [];

for i=1:data_n,
  xx = real(rect)+x(i);
  yy = imag(rect)*y(i);
  h = patch(xx, yy, 'c');
  patchH = [patchH, h];  % array of patches is constructed recursively
end

axis([0 data_n+1 0 max(y)]);
title('Animation of bubble sort');
set(gca, 'box', 'on');
set(patchH, 'erasemode', 'xor');

%---
% bubble sort
%---

fprintf('%d ', y);
fprintf('\n');

for i=data_n:-1:2,
  for j=1:i-1,

    if (y(j)>y(j+1))

        %---
        % animation
        %---

          set(patchH(j),   'FaceColor', [1 0 1]);
	  set(patchH(j+1), 'FaceColor', [1 0 1]);
          bsort_anim_switch(patchH, j, j+1);
	  set(patchH(j),   'FaceColor', [0 1 1]);
          set(patchH(j+1), 'FaceColor', [0 1 1]);

       %---
       % switch elements
       %---

          temp = y(j);
          y(j) = y(j+1);
          y(j+1) = temp;
          temp = patchH(j);
          patchH(j) = patchH(j+1);
          patchH(j+1) = temp;
          fprintf('%d ', y);
          fprintf('\n');

    end

  end
end
