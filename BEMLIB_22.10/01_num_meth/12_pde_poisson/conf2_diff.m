function [bl,bc,br, bb,bx,bt, cl,cc,cr,cb,ct] = conf2_diff ...
...
   (Nxi,Neta,xi,eta,h)

%------------------------------------------
% Differentiation coefficients for the first
% and second derivative
%
% first derivative:
%
%                  bt
%  bl bc br        bx
%                  bb
%
% second derivative:
%
%     ct
%  cl cc cr
%     cb
%
%------------------------------------------

for j=2:Neta

 for i=2:Nxi+1

 Dx1 = xi(i)-xi(i-1);
 Dx2 = xi(i+1)-xi(i);
 Dxt = Dx1+Dx2;
 hs=h(i,j)^2;

 bl(i,j) = -Dx2/(Dx1*Dxt*h(i,j));
 br(i,j) =  Dx1/(Dx2*Dxt*h(i,j));
 bc(i,j) = -bl(i,j)-br(i,j);

 cl(i,j)  = 2.0/(Dx1*Dxt*hs);
 cr(i,j)  = 2.0/(Dx2*Dxt*hs);
 cc(i,j)  = -cl(i,j)-cr(i,j);

 Dy1 = eta(j)-eta(j-1);
 Dy2 = eta(j+1)-eta(j);
 Dyt = Dy1+Dy2;

 bb(i,j) = -Dy2/(Dy1*Dyt*h(i,j));
 bt(i,j) =  Dy1/(Dy2*Dyt*h(i,j));
 bx(i,j) = -bb(i,j)-bt(i,j);

 cb(i,j) = 2.0/(Dy1*Dyt*hs);
 ct(i,j) = 2.0/(Dy2*Dyt*hs);
 cc(i,j) = cc(i,j)-cb(i,j)-ct(i,j);

 end

 bl(1,j) = bl(Nxi+1,j);
 bc(1,j) = bc(Nxi+1,j);
 br(1,j) = br(Nxi+1,j);
 bb(1,j) = bb(Nxi+1,j);
 bx(1,j) = bx(Nxi+1,j);
 bt(1,j) = bt(Nxi+1,j);

 cl(1,j) = cl(Nxi+1,j);
 cc(1,j) = cc(Nxi+1,j);
 cr(1,j) = cr(Nxi+1,j);
 cb(1,j) = cb(Nxi+1,j);
 ct(1,j) = ct(Nxi+1,j);

end

%---
% done
%---

return
