figure(2)

file2 = fopen('flow_1d_1p.vel')
N  = fscanf(file2,'%f',[1,1])
M  = fscanf(file2,'%f',[1,1])
Dx = fscanf(file2,'%f',[1,1])
Dy = fscanf(file2,'%f',[1,1])
L = fscanf(file2,'%f',[1,1])
ybottom = fscanf(file2,'%f',[1,1])
P = fscanf(file2,'%f',[N+1,M+1]);

for i=1:N+1             % set up plotting vectors
   X(i)=-L+Dx*(i-1.0);
end
for j=1:M+1
   Y(j)=ybottom + Dy*(j-1.0);
end
colormap(gray)
view(-74,10)
axis([-1 1 -1 1 -1 1])

mesh(X,Y,P')   % plotting, labelling, and formatting

xlabel('x/L')
ylabel('y/L')
zlabel('v')
hold on;

%%%%%%%%
%N  = fscanf(file2,'%f',[1,1])
%M  = fscanf(file2,'%f',[1,1])
%Dx = fscanf(file2,'%f',[1,1])
%Dy = fscanf(file2,'%f',[1,1])
%ybottom = fscanf(file2,'%f',[1,1])
%T = fscanf(file2,'%f',[N+1,M+1]);
%for i=1:N+1             % set up plotting vectors
%   X(i)=-L+Dx*(i-1.0);
%end
%for j=1:M+1
%   Y(j)=ybottom+Dy*(j-1.0);
%end
%colormap(gray)
%axis([-1 1 -1 1 -1 1])
%mesh(X,Y,T')   % plotting, labelling, and formatting
%xlabel('x')
%ylabel('y')
%zlabel('u')
%%%%%%%%5
recx = fscanf(file2,'%f',[1,1])
recy = fscanf(file2,'%f',[1,1])
%%%%%%%%5
for i=1:3
shift = (i-2.0)*L;
xp(1) = recx+shift; yp(1) = - 2.0*recy; zp(1) = 0;
xp(2) = recx+shift; yp(2) =      0; zp(2) = 0;
xp(3) =-recx+shift; yp(3) =      0; zp(3) = 0;
xp(4) =-recx+shift; yp(4) =  -2.0*recy; zp(4) = 0;
xp(5) = recx+shift; yp(5) =  -2.0*recy; zp(5) = 0;
plot3(xp,yp,zp,'k')
end
%%%%%%%%5

fclose(file2)

