function stairs (xs,ys,zs,steph,stepl,stepw,nsteps)

%=========================
% draw a straight staircase
%=========================

for i=1:nsteps

 % vertical sides of the steps

 length = (i-1)*stepl;
 height = (i-1)*steph;
 ylevel = ys+length;

 x(1) = 0.0;       y(1)=ylevel; z(1)=zs+height;
 x(2) = xs+stepw;  y(2)=ylevel; z(2)=z(1);
 x(3) = x(2);      y(3)=ylevel; z(3)=z(2)+steph;
 x(4) = 0.0;       y(4)=ylevel; z(4)=z(3);

 for j=1:4
  c(j)=z(j);
 end
 patch(x,y,z,c); hold on

 % horizontal sides of the steps

 zlevel = zs+height+steph;


 x(1)=0.0;       y(1)=ys+length;  z(1)=zlevel;
 x(2)=xs+stepw;  y(2)=y(1);       z(2)=zlevel;
 x(3)=x(2);      y(3)=y(2)+stepl; z(3)=zlevel;
 x(4)=0.0;       y(4)=y(3);       z(4)=zlevel;
 for j=1:4
  c(j)=z(j);
 end
 patch(x,y,z,c); hold on
end

% side walls

for k=1:2

xdisp = stepw*(k-1.0);

Ic=0; % counter

for i=1:nsteps
  length = (i-1)*stepl;
  height = (i-1)*steph;
  Ic=Ic+1;
  xx(Ic) = xs+xdisp;
  yy(Ic) = ys+length;
  zz(Ic) = zs+height;
  Ic=Ic+1;
  xx(Ic) = xx(Ic-1);
  yy(Ic) = yy(Ic-1);
  zz(Ic) = zz(Ic-1)+steph;
  Ic=Ic+1;
  xx(Ic) = xx(Ic-1);
  yy(Ic) = yy(Ic-1)+stepl;
  zz(Ic) = zz(Ic-1);
end

  Ic = Ic+1;
  xx(Ic) = xx(Ic-1);
  yy(Ic) = yy(Ic-1);
  zz(Ic) = zz(Ic-1)-steph;
  Ic = Ic+1;
  xx(Ic) = xx(3);
  yy(Ic) = yy(3);
  zz(Ic) = zz(3)-steph;
  Ic = Ic+1;
  xx(Ic) = xs+xdisp;
  yy(Ic) = ys;
  zz(Ic) = zs;

  for i=1:Ic
   cc(i) = zz(i);
  end

  patch(xx,yy,zz,cc);

end

%----
% done
%----

%---
return
%---
