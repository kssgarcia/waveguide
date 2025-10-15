%========================
% Code: elm_line3_dr
%
% Element discretization
%=======================
                                                                                
%-----------
% input data
%-----------
                                                                                
x1=-0.1; x2=1.4;
n=21;        % must be odd
ratio=10.0;
                                                                                
%-----------
% discretize
%-----------

xe = elm_line3 (x1,x2,n,ratio);

%-----
% plot
%-----

ye= zeros(n+1,1);
plot(xe, ye,'-o');
set(gca,'fontsize',15)


%-----
% done
%-----
