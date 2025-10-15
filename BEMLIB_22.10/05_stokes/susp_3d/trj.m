file1 = fopen('susp_3d.trj')
orbit  = fscanf(file1,'%f',[9,inf]);

%plot3(orbit(4,:),orbit(6,:),orbit(5,:),'ok')

si = size(orbit);
npts=si(2);

fc=0.05;
skip = 2;

Icount = skip

for i=1:npts
arx=fc*orbit(7,i);
ary=fc*orbit(8,i);
arz=fc*orbit(9,i);
px(1) = orbit(4,i)-arx;
py(1) = orbit(5,i)-arx;
pz(1) = orbit(6,i)-arx;
px(1) = orbit(4,i);
py(1) = orbit(5,i);
pz(1) = orbit(6,i);
px(1)=px(1)+i;
px(2) = px(1)+2.0*arx;
py(2) = py(1)+2.0*ary;
pz(2) = pz(1)+2.0*arz;
ppx=px(2);
ppy=py(2);
ppz=pz(2);
if(Icount==skip) 
hold on;
% plot3(px,-py,pz,'-')
% plot3(ppx,-ppy,ppz,'o')
 plot3(pz,px,py,'-')
 plot3(ppz,ppx,ppy,'o')
 Icount = 0;
end
Icount=Icount+1;
end

%axis([-1.1 30 -1 1 -1 1]);
%axis([0 400 -0.10 0.10 -0.10 0.10]);
%view(90,-90)
%axis off
%axis([-0.10 0.10 0 400 -0.10 0.10]);
view(17,-70)
view(74,-56)
view(66,-68)

xlabel('z')
ylabel('x')
zlabel('y')

