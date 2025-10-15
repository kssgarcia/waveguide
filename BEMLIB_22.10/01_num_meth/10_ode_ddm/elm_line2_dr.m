%========================
% Code: elm_line2_dr
%
% Driver for elm_line2
%=======================

%-----------
% input data
%-----------

x1=0.0; x2=1.0;
n=20;        % must be even
ratio=10.0;

%-----------
% discretize
%-----------

xe = elm_line2 (x1,x2,n,ratio);

%-----
% plot
%-----

ye= zeros(n+1,1);
plot(xe, ye,'-o');
set(gca,'fontsize',15)

%-----
% done
%-----
