function bsort_anum_switch(patchH, bar1, bar2)

%==========
% Switch two bars, used in sorting animation.
%
%  Roger Jang, 980703
%==========

r = 0.3;
rect = [-r -r r r] + sqrt(-1)*[0 1 1 0];

for i=0:abs(bar2-bar1)*100;
   pos_x = bar1 + i*0.01*sign(bar2-bar1);
   xx = real(rect)+pos_x;
   set(patchH(bar1), 'xdata', xx);
   pos_x = bar2 - i*0.01*sign(bar2-bar1);
   xx = real(rect)+pos_x;
   set(patchH(bar2), 'xdata', xx);
   drawnow;
end
