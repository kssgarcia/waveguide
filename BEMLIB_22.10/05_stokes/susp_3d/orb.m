file1 = fopen('susp_3d.orb')
orbit  = fscanf(file1,'%f',[3,inf]);


%-- prolate:
 plot3(orbit(1,:),-orbit(3,:),orbit(2,:),'ok')
%---

%plot3(orbit(1,:),-orbit(3,:),-orbit(2,:),'k')
%plot3(-orbit(1,:),-orbit(3,:),-orbit(2,:),'k')
%--

hold on;

[z, y, x]= sphere;
mesh(0.99*x, 0.99*y, 0.99*z);

xlabel('x')
ylabel('z')
zlabel('y')
view(30,18)


