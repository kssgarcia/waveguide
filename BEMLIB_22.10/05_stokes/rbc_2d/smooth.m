function fsmooth = smooth(N,f,Npass)

%============
% Smoothing using the five-point formula
% of Longuet-Higgins and Cokelet (1976)
%============

%---
if(Npass>0)
%---

for pass=1:Npass

%---
% wrap
%---

 f(N+1)=f(1);
 f(N+2)=f(2);
 f(N+3)=f(3);

 fsm(1) = (-f(N-1)+4*f(N)+10*f(1)+4*f(2)-f(3))/16.0;
 fsm(2) = (-f(N)  +4*f(1)+10*f(2)+4*f(3)-f(4))/16.0;

 for i=3:N+1
  fsm(i) = (-f(i-2)+4*f(i-1)+10*f(i)+4*f(i+1)-f(i+2))/16.0;
 end

 f=fsm;

end

%---
end
%---

fsmooth = f(1:N+1);

%---
% done
%---

return
