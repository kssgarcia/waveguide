function f = twospheres_vel(x)

%=======================
% Velocity of the center of a sphere relative
% to another sphere placed at the origin
% in simple shear flow:
%
% u_x = y, u_y=0, u_z=0
%
% Use equations in Section 2.2
% of da Cunha & Hinch (1996)
% J. Fluid Mech. 309, 211-223
%
% SYMBOLS:
%
% x: particle separation
%=======================

rs = x(1)^2+x(2)^2+x(3)^2;
r  = sqrt(rs);
rc = r*rs;
rq = rs*rs;
rp = rc*rs;
rx = rc*rc;
re = rp*rc;
rn = re*r;
rt = re*rs;
rl = rt*r;
rv = rt*rs;

if(r > 2.5) 
 A = 5.0D0/rc - 8.0D0/rp + 25.0D0/rx - 35.0D0/re + 125.0D0/rn ...
     -102.0D0/rt + 12.5D0/rl + 430.0D0/rv;
 B = (16.0D0/rp + 10.0D0/re - 36.0D0/rt - 25.0D0/rl - 36.0D0/rv)/3.00;

elseif(r<2.5 & r> 2.01) 

 A = -4.3833D0 + 17.7176D0/r + 14.8204D0/rs - 92.4471D0/rc ...
     -46.3151D0/rq + 232.2304D0/rp;
 B = - 3.1918D0 + 12.3641D0/r + 11.4615D0/rs - 65.2926/rc ...
     -36.4909D0/rq + 154.8074D0/rp;

else

 C = -log(r-2.0D0);
 Cs = C^2;
 A = 16.3096D0/r-7.1548;
 B = 2.0D0*(0.4056D0*Cs+1.49681D0*C-1.9108D0)/ ...
   (r*(Cs+6.04250D0*C+6.32549D0));
end

e = x(1)*x(2)*(B-A)/rs;

f(1) = x(2) + e*x(1) - 0.50D0*B*x(2);
f(2) =        e*x(2) - 0.50D0*B*x(1);
f(3) =        e*x(3);

%---
% done
%---

return
